// ============================================================
// SVC_GP_vecchia_optimized_fixed.stan
//
// Optimizations vs original:
//   1. MEMORY: A[p,N,N] eliminated → X_svc[N,p] passed, product on-the-fly
//   2. MEMORY: dist[N,N] eliminated → nn_dist[N,M,M] (compact neighbor dists)
//   3. SPEED:  exp(-d/phi) factored out of SVC loop (shared phi)
//   4. SAMPLING: log-parameterization for sigma, phi, tau
//   5. SAMPLING: tighter, better-centered priors
// ============================================================

functions {
  // Vecchia conditional: returns [B_1..B_{l-1}, F] where
  //   B = regression coefficients, F = conditional variance
  vector get_B_F(matrix mini_cov) {
    int l = rows(mini_cov);
    vector[l] result;
    if (l == 1) {
      result[1] = mini_cov[1, 1];
    } else {
      row_vector[l-1] Bsi = mdivide_right_spd(
        mini_cov[1, 2:l], mini_cov[2:l, 2:l]);
      result[1:(l-1)] = to_vector(Bsi);
      result[l] = mini_cov[1, 1] - Bsi * mini_cov[2:l, 1];
    }
    return result;
  }

  real vecchia_lpdf(vector Z,
                    vector sigma, real phi, real sigma_e,
                    matrix X_svc,           // N × p  (replaces A[p,N,N])
                    array[,,] real nn_dist,  // N × M × M (replaces dist[N,N])
                    vector beta,
                    array[] int nb_NN, array[,] int NN, matrix X,
                    array[] int uncensored_idx) {

    int n = num_elements(Z);
    int p = cols(X_svc);
    real lp = 0;

    // Pre-compute X * beta once
    vector[n] Xbeta = X * beta;

    for (i in 1:n) {
      int m_i = nb_NN[i];

      // Build mini covariance matrix using COMPACT data
      matrix[m_i, m_i] minicov;
      for (j in 1:m_i) {
        for (k in j:m_i) {
          real d_jk = nn_dist[i, j, k];

          // OPTIMIZATION: factor exp() out of the p-loop (fixed phi)
          real exp_val = exp(-d_jk / phi);
          real cov_val = 0.0;
          for (l in 1:p) {
            // A[l, neigh[j], neigh[k]] = X_svc[neigh[j], l] * X_svc[neigh[k], l]
            // But we use nn_Xsvc via NN indices
            cov_val += sigma[l] * X_svc[NN[i,j], l] * X_svc[NN[i,k], l];
          }
          cov_val *= exp_val;

          minicov[j, k] = cov_val;
          minicov[k, j] = cov_val;
        }
      }
      // Add nugget to diagonal
      for (j in 1:m_i)
        minicov[j, j] += sigma_e;

      // Conditional regression + variance
      vector[m_i] bf = get_B_F(minicov);
      real cond_var = bf[m_i];

      // Numerical guard
      if (is_nan(cond_var) || cond_var <= 0)
        return negative_infinity();

      // Conditional mean
      real mu_i = Xbeta[NN[i, 1]];
      if (m_i > 1) {
        for (j in 2:m_i)
          mu_i += bf[j-1] * (Z[NN[i, j]] - Xbeta[NN[i, j]]);
      }

      if (is_nan(mu_i) || is_inf(mu_i))
        return negative_infinity();

      // Likelihood contribution
      real sd_i = sqrt(cond_var);
      if (uncensored_idx[i] == 0)
        lp += normal_lpdf(Z[i] | mu_i, sd_i);
      else
        lp += normal_lcdf(Z[i] | mu_i, sd_i);
    }
    return lp;
  }
}

data {
  int<lower=1> N;
  vector[N] Z;
  int<lower=2> M;                    // max neighbors (including self)
  int<lower=1> p;                    // number of SVC components
  int<lower=1> pX;                   // number of covariates
  matrix[N, pX] X;                   // design matrix
  matrix[N, p] X_svc;               // SVC columns (NEW — replaces A[p,N,N])
  array[N, M] int NN;               // neighbor indices
  array[N] int nb_NN;               // number of neighbors per point
  array[N, M, M] real nn_dist;      // compact distances (NEW — replaces dist[N,N])
  array[N] int uncensored_idx;      // 0 = observed, 1 = censored

  // Prior hyperparameters (on LOG scale for sigma, phi, tau)
  vector[pX] prior_mean_alpha;
  vector<lower=0>[pX] prior_sd_alpha;
  vector[p] prior_mean_log_sigma;         // log(MLE_sigma)
  real<lower=0> prior_sd_log_sigma;  // e.g. 1.0
  real prior_mean_log_phi;           // log(MLE_phi)
  real<lower=0> prior_sd_log_phi;    // e.g. 1.0
  real prior_mean_log_tau;           // log(MLE_tau)
  real<lower=0> prior_sd_log_tau;    // e.g. 1.0
}

parameters {
  vector[pX] alpha;
  vector[p] log_sigma;     // unconstrained
  real log_phi;             // unconstrained (fixed range)
  real log_tau;             // unconstrained
}

transformed parameters {
  vector<lower=0>[p] sigma = exp(log_sigma);
  real<lower=0> phi = exp(log_phi);
  real<lower=0> sigma_e = exp(log_tau);
}

model {
  // Priors — normal on log scale = lognormal on original scale
  alpha ~ normal(prior_mean_alpha, prior_sd_alpha);
  log_sigma ~ normal(prior_mean_log_sigma, prior_sd_log_sigma);
  log_phi ~ normal(prior_mean_log_phi, prior_sd_log_phi);
  log_tau ~ normal(prior_mean_log_tau, prior_sd_log_tau);

  // Vecchia log-likelihood
  Z ~ vecchia(sigma, phi, sigma_e, X_svc, nn_dist,
              alpha, nb_NN, NN, X, uncensored_idx);
}

generated quantities {
  // Per-observation log-likelihood for LOO-CV / WAIC (optional)
  // Uncomment if needed — adds ~50% runtime
  // vector[N] log_lik;
  // ... (similar loop as vecchia_lpdf)
}
