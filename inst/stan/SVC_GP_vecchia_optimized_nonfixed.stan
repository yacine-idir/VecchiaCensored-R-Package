// ============================================================
// SVC_GP_vecchia_optimized_nonfixed.stan
//
// Same optimizations as fixed version, but phi is a vector (one per SVC)
// ============================================================

functions {
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
                    vector sigma, vector phi, real sigma_e,
                    matrix X_svc,
                    array[,,] real nn_dist,
                    vector beta,
                    array[] int nb_NN, array[,] int NN, matrix X,
                    array[] int uncensored_idx) {

    int n = num_elements(Z);
    int p = cols(X_svc);
    real lp = 0;
    vector[n] Xbeta = X * beta;

    for (i in 1:n) {
      int m_i = nb_NN[i];
      matrix[m_i, m_i] minicov;

      for (j in 1:m_i) {
        for (k in j:m_i) {
          real d_jk = nn_dist[i, j, k];
          real cov_val = 0.0;
          // Non-fixed: each SVC component has its own phi
          for (l in 1:p) {
            cov_val += sigma[l] * exp(-d_jk / phi[l])
                       * X_svc[NN[i,j], l] * X_svc[NN[i,k], l];
          }
          minicov[j, k] = cov_val;
          minicov[k, j] = cov_val;
        }
      }
      for (j in 1:m_i)
        minicov[j, j] += sigma_e;

      vector[m_i] bf = get_B_F(minicov);
      real cond_var = bf[m_i];

      if (is_nan(cond_var) || cond_var <= 0)
        return negative_infinity();

      real mu_i = Xbeta[NN[i, 1]];
      if (m_i > 1) {
        for (j in 2:m_i)
          mu_i += bf[j-1] * (Z[NN[i, j]] - Xbeta[NN[i, j]]);
      }

      if (is_nan(mu_i) || is_inf(mu_i))
        return negative_infinity();

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
  int<lower=2> M;
  int<lower=1> p;
  int<lower=1> pX;
  matrix[N, pX] X;
  matrix[N, p] X_svc;
  array[N, M] int NN;
  array[N] int nb_NN;
  array[N, M, M] real nn_dist;
  array[N] int uncensored_idx;

  vector[pX] prior_mean_alpha;
  vector<lower=0>[pX] prior_sd_alpha;
  real prior_mean_log_sigma;
  real<lower=0> prior_sd_log_sigma;
  real prior_mean_log_phi;
  real<lower=0> prior_sd_log_phi;
  real prior_mean_log_tau;
  real<lower=0> prior_sd_log_tau;
}

parameters {
  vector[pX] alpha;
  vector[p] log_sigma;
  vector[p] log_phi;       // one phi per SVC component
  real log_tau;
}

transformed parameters {
  vector<lower=0>[p] sigma = exp(log_sigma);
  vector<lower=0>[p] phi = exp(log_phi);
  real<lower=0> sigma_e = exp(log_tau);
}

model {
  alpha ~ normal(prior_mean_alpha, prior_sd_alpha);
  log_sigma ~ normal(prior_mean_log_sigma, prior_sd_log_sigma);
  log_phi ~ normal(prior_mean_log_phi, prior_sd_log_phi);
  log_tau ~ normal(prior_mean_log_tau, prior_sd_log_tau);

  Z ~ vecchia(sigma, phi, sigma_e, X_svc, nn_dist,
              alpha, nb_NN, NN, X, uncensored_idx);
}
