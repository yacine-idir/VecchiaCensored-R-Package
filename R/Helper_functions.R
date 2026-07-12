# ============================================================================
# FONCTIONS BAYESIENNES
# ============================================================================

#' @importFrom GpGp order_maxmin find_ordered_nn fit_model
#' @importFrom Rfast Dist
#' @importFrom Matrix t crossprod diag
#' @importFrom cmdstanr cmdstan_model
fit_SVC_nocensored_bayesian_Vecchia <- function(Y_obs, locs_obs, X_obs, M,
                                                svc_indices = NULL,
                                                censored_indices,
                                                chains, iter_warmup,
                                                iter_sampling,
                                                parallel_chains, FIXED) {
  n <- length(Y_obs)
  X_obs <- as.matrix(X_obs)

  # ── 1. Order data + build NN ─────────────────────────────
  if (sum(censored_indices) > 0) {
    censored_indices_pos <- which(censored_indices == 1)
    non_censored_indices <- setdiff(1:n, censored_indices_pos)

    locs_cen   <- locs_obs[censored_indices_pos, , drop = FALSE]
    locs_n_cen <- locs_obs[non_censored_indices, , drop = FALSE]
    Y_cen      <- Y_obs[censored_indices_pos]
    Y_n_cen    <- Y_obs[non_censored_indices]
    X_cen      <- X_obs[censored_indices_pos, , drop = FALSE]
    X_n_cen    <- X_obs[non_censored_indices, , drop = FALSE]

    ord <- GpGp::order_maxmin(locs_n_cen)
    locs_n_cen <- locs_n_cen[ord, , drop = FALSE]
    Y_n_cen    <- Y_n_cen[ord]
    X_n_cen    <- X_n_cen[ord, , drop = FALSE]   # Bug 4 fix: drop=FALSE

    all_locs <- rbind(locs_n_cen, locs_cen)
    all_X    <- rbind(X_n_cen, X_cen)
    all_Y    <- c(Y_n_cen, Y_cen)

    censored_indices <- c(rep(0, length(non_censored_indices)),
                          rep(1, length(censored_indices_pos)))

    M <- min(M, n - 1)
    NN <- GpGp::find_ordered_nn(all_locs[, 1:2], M)

    # Fix NN for censored points: neighbors from non-censored only
    dmat_full <- Rfast::Dist(all_locs)
    n_ncen <- length(non_censored_indices)
    dmat_sub <- dmat_full[(n_ncen + 1):n, 1:n_ncen, drop = FALSE]
    for (i in 1:nrow(dmat_sub))
      NN[i + n_ncen, 2:(M + 1)] <- order(dmat_sub[i, ])[1:M]

    NN[which(is.na(NN))] <- 0
    nb_NN <- apply(NN, 1, function(x) sum(x != 0))

  } else {
    ord <- GpGp::order_maxmin(locs_obs)
    all_locs <- locs_obs[ord, , drop = FALSE]
    all_X    <- X_obs[ord, , drop = FALSE]        # Bug 4 fix: drop=FALSE
    all_Y    <- Y_obs[ord]

    M <- min(M, n - 1)
    NN <- GpGp::find_ordered_nn(all_locs, M)
    NN[which(is.na(NN))] <- 0
    nb_NN <- apply(NN, 1, function(x) sum(x != 0))

    dmat_full <- NULL
  }

  M_stan <- M + 1  # Stan M includes self

  # ── 2. GpGp frequentist init ─────────────────────────────
  gpgp_fit <- GpGp::fit_model(
    y = all_Y, locs = all_locs, X = all_X,
    covfun_name = "exponential_isotropic", silent = TRUE, m_seq = c(10))

  mle_sigma <- gpgp_fit$covparms[1]
  mle_phi   <- gpgp_fit$covparms[2]
  mle_tau   <- max(mle_sigma * gpgp_fit$covparms[3], 1e-4)
  mle_beta  <- as.vector(gpgp_fit$betahat)
  mle_se    <- as.vector(gpgp_fit$sebeta)

  cat(sprintf("  [GpGp MLE] sigma=%.2f, phi=%.0f, tau=%.4f\n",
              mle_sigma, mle_phi, mle_tau))

  # ── 3. SVC setup ─────────────────────────────────────────
  if (is.null(svc_indices)) {
    p_svc <- 1
    X_svc <- matrix(1, nrow = n, ncol = 1)
  } else {
    p_svc <- length(svc_indices)
    X_svc <- all_X[, svc_indices, drop = FALSE]
  }

  # ── 4. COMPACT DATA: nn_dist[N, M, M] ───────────────────
  cat("  [Compact data] Computing nn_dist...")
  t0 <- proc.time()

  nn_dist <- array(0.0, dim = c(n, M_stan, M_stan))

  if (!is.null(dmat_full)) {
    for (i in 1:n) {
      nneigh <- nb_NN[i]
      for (j in 1:nneigh) {
        nj <- NN[i, j]
        for (k in j:nneigh) {
          nk <- NN[i, k]
          d <- dmat_full[nj, nk]
          nn_dist[i, j, k] <- d
          nn_dist[i, k, j] <- d
        }
      }
    }
    rm(dmat_full)
  } else {
    for (i in 1:n) {
      nneigh <- nb_NN[i]
      for (j in 1:nneigh) {
        nj <- NN[i, j]
        for (k in j:nneigh) {
          nk <- NN[i, k]
          d <- sqrt(sum((all_locs[nj, ] - all_locs[nk, ])^2))
          nn_dist[i, j, k] <- d
          nn_dist[i, k, j] <- d
        }
      }
    }
  }

  t1 <- (proc.time() - t0)[3]
  cat(sprintf(" done (%.1fs)\n", t1))

  old_mem <- (n * n * 8 + p_svc * n * n * 8) / 1e6
  new_mem <- (n * M_stan * M_stan * 8 + n * p_svc * 8) / 1e6
  cat(sprintf("  [Memory] Old: %.0f MB (dist + A) -> New: %.0f MB (nn_dist + X_svc) = %.0fx\n",
              old_mem, new_mem, old_mem / new_mem))

  # ── 5. PRIORS ────────────────────────────────────────────
  PRIOR_SD_LOG <- 1.0
  K_PRIOR_BETA <- 10
  prior_sd_alpha <- K_PRIOR_BETA * pmax(mle_se, 0.1)

  int_idx <- NULL
  for (j in seq_len(ncol(all_X))) {
    if (all(abs(all_X[, j] - 1) < 1e-10) && j %in% svc_indices) {
      int_idx <- which(svc_indices == j)
      break
    }
  }

  init_sigma <- rep(0.1, p_svc)
  if (!is.null(int_idx)) {
    init_sigma[int_idx] <- mle_sigma
  } else {
    init_sigma <- rep(mle_sigma / p_svc, p_svc)
  }

  # ── 6. Stan data ─────────────────────────────────────────
  data_stan <- list(
    N   = n,
    Z   = all_Y,
    M   = M_stan,
    p   = p_svc,
    pX  = ncol(all_X),
    X   = all_X,
    X_svc = X_svc,
    NN  = NN,
    nb_NN = nb_NN,
    nn_dist = nn_dist,
    uncensored_idx = censored_indices,
    prior_mean_alpha     = mle_beta,
    prior_sd_alpha       = prior_sd_alpha,
    log_sigma_prior_means <- rep(log(0.1), p_svc)       # par défaut : petit
    if (!is.null(int_idx))
    log_sigma_prior_means[int_idx] <- log(mle_sigma)  # intercept : MLE
    prior_mean_log_sigma = log_sigma_prior_means         # vecteur    prior_sd_log_sigma   = PRIOR_SD_LOG,
    prior_mean_log_phi   = log(max(mle_phi, 0.01)),
    prior_sd_log_phi     = PRIOR_SD_LOG,
    prior_mean_log_tau   = log(max(mle_tau, 0.01)),
    prior_sd_log_tau     = PRIOR_SD_LOG
  )

  cat(sprintf("  [Priors (log scale, SD=%.1f)]\n", PRIOR_SD_LOG))
  cat(sprintf("    sigma: 95%% CI ~ [%.1f, %.1f]\n",
              mle_sigma * exp(-2 * PRIOR_SD_LOG), mle_sigma * exp(2 * PRIOR_SD_LOG)))
  cat(sprintf("    phi:   95%% CI ~ [%.0f, %.0f]\n",
              mle_phi * exp(-2 * PRIOR_SD_LOG), mle_phi * exp(2 * PRIOR_SD_LOG)))
  cat(sprintf("    tau:   95%% CI ~ [%.4f, %.4f]\n",
              mle_tau * exp(-2 * PRIOR_SD_LOG), mle_tau * exp(2 * PRIOR_SD_LOG)))

  # ── 7. Stan model ───────────────────────────────────────
  if (FIXED) {
    stan_file <- system.file("stan", "SVC_GP_vecchia_optimized_fixed.stan",
                             package = "VecchiaCensored")
    if (stan_file == "" || !file.exists(stan_file))
      stan_file <- "SVC_GP_vecchia_optimized_fixed.stan"
    cat("  [Stan] Fixed range model\n")
  } else {
    stan_file <- system.file("stan", "SVC_GP_vecchia_optimized_nonfixed.stan",
                             package = "VecchiaCensored")
    if (stan_file == "" || !file.exists(stan_file))
      stan_file <- "SVC_GP_vecchia_optimized_nonfixed.stan"
    cat("  [Stan] Non-fixed range model\n")
  }

  if (!file.exists(stan_file))
    stop("Optimized Stan file not found: ", stan_file,
         "\nPlace it in inst/stan/ or the working directory.")

  model <- cmdstanr::cmdstan_model(stan_file)

  # ── 8. Init at frequentist estimates ─────────────────────
  init_fn <- function() {
    init <- list(
      alpha     = mle_beta,
      log_sigma = log(pmax(init_sigma, 0.01)),
      log_tau   = log(max(mle_tau, 0.01))
    )
    if (FIXED) {
      init$log_phi <- log(max(mle_phi, 0.01))
    } else {
      init$log_phi <- rep(log(max(mle_phi, 0.01)), p_svc)
    }
    return(init)
  }

  cat("  [Init] log_sigma =", round(log(pmax(init_sigma, 0.01)), 2),
      "| log_phi =", round(log(mle_phi), 2),
      "| log_tau =", round(log(mle_tau), 2), "\n")

  # ── 9. Sample ───────────────────────────────────────────
  cat("  [Sampling] chains=", chains, ", warmup=", iter_warmup,
      ", sampling=", iter_sampling, "\n")

  MCMC <- model$sample(
    data            = data_stan,
    chains          = chains,
    parallel_chains = parallel_chains,
    iter_warmup     = iter_warmup,
    iter_sampling   = iter_sampling,
    init            = init_fn,
    adapt_delta     = 0.7
  )

  return(list(MCMC, data_stan))
}


# ============================================================================
# FONCTIONS FREQUENTISTES
# ============================================================================

# ────────────────────────────────────────────────────────────
# SVC non censuré : vraisemblance profilée + log-param
# ────────────────────────────────────────────────────────────

fit_SVC_nocensored_freq_Vecchia <- function(Y_obs, locs_obs, X_obs, M,
                                            svc_indices, fixed_range) {

  ord      <- GpGp::order_maxmin(locs_obs)
  Y_obs    <- Y_obs[ord]
  locs_obs <- locs_obs[ord, , drop = FALSE]
  X_obs    <- X_obs[ord, , drop = FALSE]

  n     <- length(Y_obs)
  p     <- ncol(X_obs)
  n_svc <- length(svc_indices)
  M     <- min(M, n - 1)
  NN    <- GpGp::find_ordered_nn(locs_obs, M)
  X_svc <- X_obs[, svc_indices, drop = FALSE]

  # ── GpGp init ───────────────────────────────────────────
  gpgp_init <- GpGp::fit_model(
    y = Y_obs, locs = locs_obs, X = X_obs,
    covfun_name = "exponential_isotropic", silent = TRUE, m_seq = c(10))

  sigma_total <- gpgp_init$covparms[1]
  phi0        <- gpgp_init$covparms[2]
  tau0        <- max(sigma_total * gpgp_init$covparms[3], 1e-4)

  cat(sprintf("  [GpGp] sigma=%.2f, phi=%.0f, tau=%.4f, loglik=%.1f\n",
              sigma_total, phi0, tau0, gpgp_init$loglik))

  # ── Identifier la colonne intercept ──────────────────────
  int_idx <- NULL
  for (j in seq_len(p)) {
    if (all(abs(X_obs[, j] - 1) < 1e-10) && (j %in% svc_indices)) {
      int_idx <- which(svc_indices == j)
      break
    }
  }

  # ── Helpers pour la paramétrisation ──────────────────────
  log_to_covlist <- function(log_cov) {
    theta <- exp(log_cov)
    cov_list <- vector("list", n_svc + 1)
    if (fixed_range && n_svc > 1) {
      phi_val <- theta[n_svc + 1]
      for (i in seq_len(n_svc))
        cov_list[[i]] <- c(sigma_sq = theta[i], phi = phi_val)
      cov_list$nugget <- theta[length(theta)]
    } else {
      for (i in seq_len(n_svc)) {
        idx <- (2 * (i - 1) + 1):(2 * i)
        cov_list[[i]] <- c(sigma_sq = theta[idx[1]], phi = theta[idx[2]])
      }
      cov_list$nugget <- theta[length(theta)]
    }
    return(cov_list)
  }

  # ── Bornes sur l'échelle log ─────────────────────────────
  LOG_LB_SIGMA <- log(1e-6);  LOG_UB_SIGMA <- log(1e4)
  LOG_LB_PHI   <- log(1e-1);  LOG_UB_PHI   <- log(1e8)
  LOG_LB_TAU   <- log(1e-8);  LOG_UB_TAU   <- log(1e4)

  if (fixed_range && n_svc > 1) {
    log_lower <- c(rep(LOG_LB_SIGMA, n_svc), LOG_LB_PHI, LOG_LB_TAU)
    log_upper <- c(rep(LOG_UB_SIGMA, n_svc), LOG_UB_PHI, LOG_UB_TAU)
  } else {
    log_lower <- numeric(2 * n_svc + 1)
    log_upper <- numeric(2 * n_svc + 1)
    for (i in seq_len(n_svc)) {
      log_lower[2 * (i - 1) + 1] <- LOG_LB_SIGMA
      log_lower[2 * (i - 1) + 2] <- LOG_LB_PHI
      log_upper[2 * (i - 1) + 1] <- LOG_UB_SIGMA
      log_upper[2 * (i - 1) + 2] <- LOG_UB_PHI
    }
    log_lower[2 * n_svc + 1] <- LOG_LB_TAU
    log_upper[2 * n_svc + 1] <- LOG_UB_TAU
  }

  # ── Vraisemblance profilée (β éliminé par GLS) ──────────
  profile_nll <- function(log_cov) {
    if (any(log_cov < LOG_LB_SIGMA - 1) || any(log_cov > LOG_UB_PHI + 1))
      return(1e10)

    tryCatch({
      cov_list <- log_to_covlist(log_cov)
      U <- precision_decomp_sparse(cov_list, X_svc, NN, locs_obs)

      UX   <- U %*% X_obs
      UY   <- U %*% Y_obs
      XtX  <- as.matrix(Matrix::crossprod(UX))
      XtY  <- as.matrix(Matrix::crossprod(UX, UY))

      if (rcond(XtX) < 1e-15) return(1e10)

      R    <- chol(XtX)
      beta <- backsolve(R, forwardsolve(t(R), XtY))

      resid <- drop(Y_obs - X_obs %*% beta)
      z     <- U %*% resid
      quad  <- sum(z * z)
      logdet <- 2 * sum(log(Matrix::diag(U)))

      nll <- -0.5 * logdet + 0.5 * (n * log(2 * pi) + quad)
      if (!is.finite(nll)) return(1e10)
      return(nll)
    }, error = function(e) 1e10)
  }

  # ── Point de départ ──────────────────────────────────────
  make_start <- function(log_sigmas, log_phi, log_tau) {
    if (fixed_range && n_svc > 1) {
      c(log_sigmas, log_phi, log_tau)
    } else {
      log_cov <- numeric(2 * n_svc + 1)
      for (i in seq_len(n_svc)) {
        log_cov[2 * (i - 1) + 1] <- log_sigmas[i]
        log_cov[2 * (i - 1) + 2] <- log_phi[min(i, length(log_phi))]
      }
      log_cov[2 * n_svc + 1] <- log_tau
      log_cov
    }
  }

  ls1 <- rep(log(1e-3), n_svc)
  if (!is.null(int_idx)) ls1[int_idx] <- log(sigma_total)
  else ls1 <- rep(log(sigma_total / n_svc), n_svc)
  start_par <- make_start(ls1, log(phi0), log(tau0))
  start_par <- pmax(pmin(start_par, log_upper), log_lower)

  # ── Optimization ─────────────────────────────────────────
  nll0 <- profile_nll(start_par)
  cat(sprintf("  [Init] loglik = %.1f -> ", -nll0))

  best_opt <- tryCatch(
    optim(start_par, profile_nll, method = "L-BFGS-B",
          lower = log_lower, upper = log_upper,
          control = list(maxit = 500, factr = 1e7)),
    error = function(e) {
      cat("(L-BFGS-B failed, fallback Nelder-Mead) ")
      optim(start_par, profile_nll, method = "Nelder-Mead",
            control = list(maxit = 10000))
    })

  best_nll <- best_opt$value
  cat(sprintf("final loglik = %.1f (conv=%d)\n", -best_nll, best_opt$convergence))

  # ── Récupérer β au point optimal ─────────────────────────
  cov_list_final <- log_to_covlist(best_opt$par)
  U_final <- precision_decomp_sparse(cov_list_final, X_svc, NN, locs_obs)
  UX <- U_final %*% X_obs
  UY <- U_final %*% Y_obs
  XtX <- as.matrix(Matrix::crossprod(UX))
  XtY <- as.matrix(Matrix::crossprod(UX, UY))
  beta_final <- drop(solve(XtX, XtY))

  # ── Reconstruire la sortie au format standard ────────────
  theta_final <- exp(best_opt$par)

  if (fixed_range && n_svc > 1) {
    phi_final   <- theta_final[n_svc + 1]
    nug_final   <- theta_final[length(theta_final)]
    out_par     <- numeric(p + 2 * n_svc + 1)
    out_par[1:p] <- beta_final
    for (i in seq_len(n_svc)) {
      idx <- p + (2 * (i - 1) + 1):(2 * i)
      out_par[idx[1]] <- theta_final[i]
      out_par[idx[2]] <- phi_final
    }
    out_par[length(out_par)] <- nug_final
  } else {
    out_par <- c(beta_final, theta_final)
  }

  return(list(
    par         = out_par,
    value       = best_nll,
    convergence = best_opt$convergence,
    message     = best_opt$message
  ))
}


# ────────────────────────────────────────────────────────────
# Vraisemblance SVC non censuré (inchangée)
# ────────────────────────────────────────────────────────────

Like_SVC_nocensored_freq_Vecchia <- function(mean_parameters,
                                             cov_parameters,
                                             Y_obs,
                                             locs_obs,
                                             X_obs,
                                             NN,
                                             svc_indices) {

  X_svc <- X_obs[, svc_indices, drop = FALSE]
  L <- precision_decomp_sparse(cov_parameters, X_svc, NN, locs_obs)

  mu    <- X_obs %*% mean_parameters
  resid <- drop(Y_obs - mu)

  z    <- L %*% resid
  quad <- sum(z * z)

  logdet_prec <- 2 * sum(log(Matrix::diag(L)))

  n  <- length(resid)
  c0 <- n * log(2 * pi)
  ll <- 0.5 * logdet_prec - 0.5 * (c0 + quad)

  return(ll)
}


# ============================================================================
# PREDICTION HELPER
# ============================================================================

prediction_freq <- function(beta, sigma, phi, tau, M, Y_obs, X_obs, locs_obs,
                            svc_indices, cen_indices, X_pred, locs_pred,
                            simulation, pred_coef = FALSE) {

  X_obs  <- as.matrix(X_obs)
  X_pred <- as.matrix(X_pred)

  # ── Traitement des observations censurées ────────────────
  if (sum(cen_indices) > 0) {

    cen_indx   <- which(cen_indices == 1)
    n_cen_indx <- which(cen_indices == 0)

    n_cen   <- length(cen_indx)
    n_n_cen <- length(n_cen_indx)
    n_obs   <- length(Y_obs)

    Y_obs_n_cen    <- Y_obs[n_cen_indx]
    Y_obs_cen      <- Y_obs[cen_indx]
    X_obs_n_cen    <- X_obs[n_cen_indx, , drop = FALSE]    # Bug 3 fix
    X_obs_cen      <- X_obs[cen_indx, , drop = FALSE]      # Bug 3 fix
    locs_obs_n_cen <- locs_obs[n_cen_indx, , drop = FALSE]
    locs_obs_cen   <- locs_obs[cen_indx, , drop = FALSE]

    ord <- GpGp::order_maxmin(locs_obs_n_cen)

    ALL_X_obs_ordered    <- rbind(X_obs_n_cen[ord, , drop = FALSE], X_obs_cen)
    ALL_Y_obs_ordered    <- c(Y_obs_n_cen[ord], Y_obs_cen)
    ALL_locs_obs_ordered <- rbind(locs_obs_n_cen[ord, , drop = FALSE], locs_obs_cen)

    if (is.null(svc_indices)) {
      X_svc <- matrix(1, nrow = n_obs, ncol = 1)
      n_svc <- 1
    } else {
      n_svc <- length(svc_indices)
      X_svc <- ALL_X_obs_ordered[, svc_indices, drop = FALSE]
    }

    cov_list <- vector("list", n_svc)
    for (i in seq_len(n_svc))
      cov_list[[i]] <- c(sigma_sq = sigma[i], phi = phi[i])
    cov_list$nugget <- as.numeric(tau)

    NN <- GpGp::find_ordered_nn(ALL_locs_obs_ordered, m = M)
    dmat_full <- Rfast::Dist(ALL_locs_obs_ordered)
    dmat_sub <- dmat_full[(n_n_cen + 1):n_obs, 1:n_n_cen, drop = FALSE]
    for (i in 1:nrow(dmat_sub))
      NN[i + n_n_cen, 2:(M + 1)] <- order(dmat_sub[i, ])[1:M]

    L <- Matrix::t(precision_decomp_sparse(
      cov_list, X_svc, NN, ALL_locs_obs_ordered, n_n_cen + 1))

    y_cens_pred <- X_obs_cen %*% beta - forwardsolve(
      Matrix::t(L[(n_n_cen + 1):n_obs, (n_n_cen + 1):n_obs]),
      Matrix::t(L[1:n_n_cen, (n_n_cen + 1):n_obs]) %*%
        (Y_obs_n_cen[ord] - X_obs_n_cen[ord, , drop = FALSE] %*% beta))

    # Bug 2 fix: cond_sd = 1/L[i,i] (was sqrt of that before)
    diag_indices <- (n_n_cen + 1):n_obs
    L_diag  <- L[cbind(diag_indices, diag_indices)]
    cond_sd <- 1 / L_diag

    truncated_values <- numeric(n_cen)
    for (j in 1:n_cen) {
      if (cond_sd[j] > 0) {
        z <- (Y_obs_cen[j] - y_cens_pred[j]) / cond_sd[j]
        if (is.finite(z)) {
          pnorm_z <- pnorm(z)
          if (pnorm_z > 1e-10) {
            mills_ratio <- dnorm(z) / pnorm_z
            truncated_values[j] <- y_cens_pred[j] - cond_sd[j] * mills_ratio
          } else {
            truncated_values[j] <- Y_obs_cen[j]
          }
        } else {
          truncated_values[j] <- Y_obs_cen[j]
        }
      } else {
        truncated_values[j] <- y_cens_pred[j]
      }
    }
    Y_obs[cen_indx] <- truncated_values
  }

  # ── Prediction at new locations ──────────────────────────
  n_pred <- nrow(X_pred)
  n_obs  <- length(Y_obs)
  n      <- n_pred + n_obs

  ord      <- GpGp::order_maxmin(locs_obs)
  ord_pred <- GpGp::order_maxmin(locs_pred)

  ALL_X_ordered    <- rbind(X_obs[ord, , drop = FALSE], X_pred[ord_pred, , drop = FALSE])
  Y_obs_ordered    <- Y_obs[ord]
  ALL_locs_ordered <- rbind(locs_obs[ord, , drop = FALSE], locs_pred[ord_pred, , drop = FALSE])

  if (is.null(svc_indices)) {
    X_svc <- matrix(1, nrow = n, ncol = 1)
    n_svc <- 1
  } else {
    n_svc <- length(svc_indices)
    X_svc <- ALL_X_ordered[, svc_indices, drop = FALSE]
  }

  cov_list <- vector("list", n_svc)
  for (i in seq_len(n_svc))
    cov_list[[i]] <- c(sigma_sq = sigma[i], phi = phi[i])
  cov_list$nugget <- as.numeric(tau)

  NN <- GpGp::find_ordered_nn(ALL_locs_ordered, m = M)
  dmat <- Rfast::Dist(ALL_locs_ordered)
  dmat <- dmat[(n_obs + 1):n, 1:n_obs, drop = FALSE]
  for (i in 1:nrow(dmat))
    NN[i + n_obs, 2:(M + 1)] <- order(dmat[i, ])[1:M]

  L <- Matrix::t(precision_decomp_sparse(
    cov_list, X_svc, NN, ALL_locs_ordered, n_obs + 1))

  prediction <- X_pred[ord_pred, , drop = FALSE] %*% beta - forwardsolve(
    Matrix::t(L[(n_obs + 1):n, (n_obs + 1):n]),
    Matrix::t(L[1:n_obs, (n_obs + 1):n]) %*%
      (Y_obs[ord] - X_obs[ord, , drop = FALSE] %*% beta))

  if (pred_coef) {
    NN_coef <- GpGp::find_ordered_nn(locs_obs[ord, , drop = FALSE], M)
    L_i <- Matrix::t(precision_decomp_sparse(
      cov_list, X_obs[ord, svc_indices, drop = FALSE],
      NN_coef, locs_obs[ord, , drop = FALSE], 1))
    coef_prediction <- matrix(NA, nrow = n_pred, ncol = length(svc_indices))
    for (i in 1:length(svc_indices)) {
      cov_Y_i <- cov_list[[i]][1] * exp(-dmat / cov_list[[i]][2])
      cov_Y_i <- sweep(cov_Y_i, 2, X_obs[ord, svc_indices, drop = FALSE][, i], `*`)
      coef_prediction[ord_pred, i] <- beta[svc_indices[i]] +
        as.vector((cov_Y_i %*% L_i %*% Matrix::t(L_i) %*%
                     (Y_obs[ord] - X_obs[ord, , drop = FALSE] %*% beta)))
    }
  }

  if (simulation) {
    iid <- rnorm(n_pred)
    sim <- prediction + forwardsolve(
      Matrix::t(L[(n_obs + 1):n, (n_obs + 1):n]), iid)
    sim[ord_pred] <- sim
    prediction[ord_pred] <- prediction
    final_list <- list(prediction = prediction, simulation = sim)
  } else {
    prediction[ord_pred] <- prediction
    final_list <- list(prediction = prediction)
  }

  if (pred_coef) {
    final_list <- append(final_list, list(coef_prediction))
  }

  return(final_list)
}


# ============================================================================
# ESTIMATION HELPER : CENSORED
# ============================================================================

fit_censored_freq_Vecchia <- function(Y_obs, locs_obs, X_obs, M,
                                      svc_indices, censored_indices,
                                      fixed_range) {

  n <- length(Y_obs)
  p <- ncol(X_obs)
  censored_indices_pos <- which(censored_indices == 1)
  non_censored_indices <- setdiff(1:n, censored_indices_pos)

  locs_cen   <- locs_obs[censored_indices_pos, , drop = FALSE]
  locs_n_cen <- locs_obs[non_censored_indices, , drop = FALSE]
  Y_cen      <- Y_obs[censored_indices_pos]
  Y_n_cen    <- Y_obs[non_censored_indices]
  X_cen      <- X_obs[censored_indices_pos, , drop = FALSE]
  X_n_cen    <- X_obs[non_censored_indices, , drop = FALSE]

  ord        <- GpGp::order_maxmin(locs_n_cen)
  locs_n_cen <- locs_n_cen[ord, , drop = FALSE]
  Y_n_cen    <- Y_n_cen[ord]
  X_n_cen    <- X_n_cen[ord, , drop = FALSE]

  all_locs <- rbind(locs_n_cen, locs_cen)
  all_X    <- rbind(X_n_cen, X_cen)
  all_Y    <- c(Y_n_cen, Y_cen)

  n_svc <- if (is.null(svc_indices)) 1 else length(svc_indices)
  M     <- min(M, n - 1)
  NN    <- GpGp::find_ordered_nn(all_locs, M)
  dmat  <- Rfast::Dist(all_locs)

  if (length(non_censored_indices) != n) {
    dmat_sub <- dmat[(length(Y_n_cen) + 1):n, 1:length(Y_n_cen), drop = FALSE]
    for (i in 1:nrow(dmat_sub))
      NN[i + length(non_censored_indices), 2:(M + 1)] <- order(dmat_sub[i, ])[1:M]
  }

  # GpGp init (Bug 5 fix)
  gpgp_init <- GpGp::fit_model(
    y = all_Y, locs = all_locs, X = all_X,
    covfun_name = "exponential_isotropic", silent = TRUE, m_seq = c(10))

  sigma_total <- gpgp_init$covparms[1]
  phi0        <- gpgp_init$covparms[2]
  tau0        <- max(sigma_total * gpgp_init$covparms[3], 1e-4)
  mean_start  <- as.vector(gpgp_init$betahat)

  # Intercept detection
  int_idx <- NULL
  for (j in seq_len(p)) {
    col_vals <- all_X[, j]
    if (all(abs(col_vals - 1) < 1e-10)) {
      if (!is.null(svc_indices) && j %in% svc_indices)
        int_idx <- which(svc_indices == j)
      else if (is.null(svc_indices))
        int_idx <- 1
      break
    }
  }

  n_mean  <- length(mean_start)
  n_n_cen <- length(non_censored_indices)

  # Bornes pour les cov params (log scale)
  LOG_LB_SIGMA <- log(1e-6);  LOG_UB_SIGMA <- log(1e4)
  LOG_LB_PHI   <- log(1e-1);  LOG_UB_PHI   <- log(1e8)
  LOG_LB_TAU   <- log(1e-8);  LOG_UB_TAU   <- log(1e4)

  if (fixed_range && n_svc > 1) {
    log_lower_cov <- c(rep(LOG_LB_SIGMA, n_svc), LOG_LB_PHI, LOG_LB_TAU)
    log_upper_cov <- c(rep(LOG_UB_SIGMA, n_svc), LOG_UB_PHI, LOG_UB_TAU)
  } else {
    log_lower_cov <- numeric(2 * n_svc + 1)
    log_upper_cov <- numeric(2 * n_svc + 1)
    for (i in seq_len(n_svc)) {
      log_lower_cov[2 * (i - 1) + 1] <- LOG_LB_SIGMA
      log_lower_cov[2 * (i - 1) + 2] <- LOG_LB_PHI
      log_upper_cov[2 * (i - 1) + 1] <- LOG_UB_SIGMA
      log_upper_cov[2 * (i - 1) + 2] <- LOG_UB_PHI
    }
    log_lower_cov[2 * n_svc + 1] <- LOG_LB_TAU
    log_upper_cov[2 * n_svc + 1] <- LOG_UB_TAU
  }

  full_lower <- c(rep(-Inf, n_mean), log_lower_cov)
  full_upper <- c(rep( Inf, n_mean), log_upper_cov)

  censored_nll <- function(par) {
    mean_par    <- par[1:n_mean]
    log_cov_par <- par[(n_mean + 1):length(par)]

    if (any(log_cov_par < LOG_LB_SIGMA - 1) || any(log_cov_par > LOG_UB_PHI + 1))
      return(1e10)

    cov_par <- exp(log_cov_par)

    cov_list <- vector("list", n_svc + 1)
    if (fixed_range && n_svc > 1) {
      phi_val <- cov_par[n_svc + 1]
      for (i in seq_len(n_svc))
        cov_list[[i]] <- c(sigma_sq = cov_par[i], phi = phi_val)
      cov_list$nugget <- cov_par[length(cov_par)]
    } else {
      for (i in seq_len(n_svc)) {
        idx <- (2 * (i - 1) + 1):(2 * i)
        cov_list[[i]] <- c(sigma_sq = cov_par[idx[1]], phi = cov_par[idx[2]])
      }
      cov_list$nugget <- cov_par[length(cov_par)]
    }

    tryCatch({
      ll <- Like_SVC_censored_freq_Vecchia(
        mean_par, cov_list, all_Y, all_locs, all_X,
        NN, svc_indices, n_no_cen = n_n_cen)
      if (!is.finite(ll)) return(1e10)
      return(-ll)
    }, error = function(e) 1e10)
  }

  make_cov_start <- function(log_sigmas, log_phi, log_tau) {
    if (fixed_range && n_svc > 1) {
      c(log_sigmas, log_phi, log_tau)
    } else {
      lc <- numeric(2 * n_svc + 1)
      for (i in seq_len(n_svc)) {
        lc[2 * (i - 1) + 1] <- log_sigmas[i]
        lc[2 * (i - 1) + 2] <- log_phi[min(i, length(log_phi))]
      }
      lc[2 * n_svc + 1] <- log_tau
      lc
    }
  }

  ls1 <- rep(log(1e-3), n_svc)
  if (!is.null(int_idx)) ls1[int_idx] <- log(sigma_total)
  else ls1 <- rep(log(sigma_total / n_svc), n_svc)
  start_par <- c(mean_start, make_cov_start(ls1, log(phi0), log(tau0)))

  cov_part <- start_par[(n_mean + 1):length(start_par)]
  cov_part <- pmax(pmin(cov_part, log_upper_cov), log_lower_cov)
  start_par[(n_mean + 1):length(start_par)] <- cov_part

  nll0 <- censored_nll(start_par)
  cat(sprintf("  [Censored init] loglik = %.1f -> ", -nll0))

  best_opt <- tryCatch(
    optim(start_par, censored_nll, method = "L-BFGS-B",
          lower = full_lower, upper = full_upper,
          control = list(maxit = 500, factr = 1e7)),
    error = function(e) {
      cat("(L-BFGS-B failed, fallback Nelder-Mead) ")
      optim(start_par, censored_nll, method = "Nelder-Mead",
            control = list(maxit = 10000))
    })

  best_nll <- best_opt$value
  cat(sprintf("final loglik = %.1f (conv=%d)\n", -best_nll, best_opt$convergence))

  # Reconstruire les paramètres
  mean_final  <- best_opt$par[1:n_mean]
  theta_final <- exp(best_opt$par[(n_mean + 1):length(best_opt$par)])

  if (fixed_range && n_svc > 1) {
    phi_final <- theta_final[n_svc + 1]
    nug_final <- theta_final[length(theta_final)]
    out_par   <- numeric(n_mean + 2 * n_svc + 1)
    out_par[1:n_mean] <- mean_final
    for (i in seq_len(n_svc)) {
      idx <- n_mean + (2 * (i - 1) + 1):(2 * i)
      out_par[idx[1]] <- theta_final[i]
      out_par[idx[2]] <- phi_final
    }
    out_par[length(out_par)] <- nug_final
  } else {
    out_par <- c(mean_final, theta_final)
  }

  return(list(
    par         = out_par,
    value       = best_nll,
    convergence = best_opt$convergence,
    message     = best_opt$message
  ))
}


# ────────────────────────────────────────────────────────────
# Vraisemblance SVC censuré (inchangée)
# ────────────────────────────────────────────────────────────

Like_SVC_censored_freq_Vecchia <- function(mean_parameters, cov_parameters, Y_obs,
                                           locs_obs, X_obs, NN, svc_indices, n_no_cen) {
  n <- length(Y_obs)

  if (is.null(svc_indices)) {
    X_svc <- matrix(1, nrow = n, ncol = 1)
  } else {
    X_svc <- X_obs[, svc_indices, drop = FALSE]
  }

  L <- precision_decomp_sparse(cov_parameters, X_svc, NN, locs_obs)
  mu <- X_obs %*% mean_parameters
  resid <- drop(Y_obs - mu)

  # Non-censored likelihood
  L_sub <- L[1:n_no_cen, 1:n_no_cen]
  z <- L_sub %*% resid[1:n_no_cen]
  quad <- sum(z * z)

  diag_indices <- 1:n_no_cen
  diag_vals <- L[cbind(diag_indices, diag_indices)]
  logdet_prec <- 2 * sum(log(diag_vals))

  c0 <- n_no_cen * log(2 * pi)
  ll_no_censored <- 0.5 * logdet_prec - 0.5 * (c0 + quad)

  # Censored likelihood
  L <- Matrix::t(L)

  cond_exp <- mu[(n_no_cen + 1):n] -
    forwardsolve(
      Matrix::t(L[(n_no_cen + 1):n, (n_no_cen + 1):n]),
      Matrix::t(L[1:n_no_cen, (n_no_cen + 1):n]) %*% (resid[1:n_no_cen])
    )

  diag_indices_cen <- (n_no_cen + 1):n
  cond_diag_vals <- L[cbind(diag_indices_cen, diag_indices_cen)]
  cond_sd <- 1 / cond_diag_vals

  ll_cen <- pnorm(q = Y_obs[(n_no_cen + 1):n], mean = cond_exp, sd = cond_sd, log.p = TRUE)

  ll <- ll_no_censored + sum(ll_cen)
  return(ll)
}
