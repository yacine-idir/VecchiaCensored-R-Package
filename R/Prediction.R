#' Make predictions from a fitted spatial model
#'
#' @param fit_object A fitted model object from fit_model()
#' @param locs_pred Matrix of prediction locations
#' @param X_pred Design matrix for prediction locations
#' @param M Number of nearest neighbors for Vecchia approximation
#' @param simulations Logical, generate conditional simulations
#' @param n_simulations Number of simulations (for non-Bayesian models)
#' @param pred_coef Logical, predict spatially varying coefficients (SVC only)
#' @param thin Integer, use every thin-th MCMC sample (Bayesian only, default 1)
#' @return A list containing predictions and optionally simulations/coefficients
#' @export
prediction <- function(fit_object, locs_pred, X_pred, M,
                        simulations = FALSE, n_simulations = NULL,
                        pred_coef = FALSE, thin = 1) {

  locs_pred <- as.matrix(locs_pred)
  X_pred    <- as.matrix(X_pred)
  n_pred    <- nrow(locs_pred)

  model_type <- fit_object$model_info$model_type

  # ==================================================================
  # GpGp simple
  # ==================================================================
  if (model_type == "gpgp") {

    if (pred_coef)
      warning("pred_coef is only available for SVC models. Ignoring.")

    pred_vec <- GpGp::predictions(
      fit = fit_object$raw_model, locs_pred = locs_pred, X_pred = X_pred)

    if (simulations) {
      if (is.null(n_simulations)) n_simulations <- 1
      sims <- GpGp::cond_sim(
        fit_object$raw_model, locs_pred = locs_pred,
        X_pred = X_pred, nsims = n_simulations)
      # FIX Bug 1: GpGp::predictions returns a vector, not a list
      return(list(prediction = pred_vec, simulations = sims))
    } else {
      return(list(prediction = pred_vec))
    }

  # ==================================================================
  # SVC frequentist (censored or not)
  # ==================================================================
  } else if (model_type %in% c("svc_freq", "svc_freq_censored")) {

    result <- VecchiaCensored:::prediction_freq(
      beta        = fit_object$parameters$beta,
      sigma       = fit_object$parameters$sigma,
      phi         = fit_object$parameters$phi,
      tau         = fit_object$parameters$tau,
      M           = M,
      Y_obs       = fit_object$data$Y_obs,
      X_obs       = fit_object$data$X_obs,
      locs_obs    = fit_object$data$locs_obs,
      svc_indices = fit_object$model_info$svc_indices,
      cen_indices = fit_object$model_info$censored_indices,
      X_pred      = X_pred,
      locs_pred   = locs_pred,
      simulation  = simulations,
      pred_coef   = pred_coef
    )

    out <- list(prediction = result$prediction)
    if (simulations) out$simulations <- matrix(result$simulation, ncol = 1)
    if (pred_coef) {
      idx <- if (simulations) 3 else 2
      out$coef_prediction <- result[[idx]]
    }
    return(out)

  # ==================================================================
  # SVC Bayesian
  # ==================================================================
  } else if (model_type == "svc_bayesian") {

    beta_samples  <- fit_object$parameters$beta    # iterations × pX
    sigma_samples <- fit_object$parameters$sigma   # iterations × n_svc
    phi_samples   <- fit_object$parameters$phi     # iterations × n_svc
    tau_samples   <- fit_object$parameters$tau     # vector of length iterations

    n_iter_total <- length(tau_samples)

    # Thinning
    iter_idx <- seq(1, n_iter_total, by = thin)
    n_iter   <- length(iter_idx)
    if (thin > 1)
      cat(sprintf("  [Bayesian prediction] Using %d/%d samples (thin=%d)\n",
                  n_iter, n_iter_total, thin))

    # ── FIX Bug 2: ensure data and censored_indices are consistent ──
    # fit_model for Bayesian stores ORDERED data but ORIGINAL censored_indices
    # We need censored_indices that match the stored data order
    cen_indices_for_pred <- fit_object$model_info$censored_indices

    if (!is.null(fit_object$data$ordered_indices)) {
      # Data is stored in ordered form → reorder censored_indices to match
      ord_idx <- fit_object$data$ordered_indices
      # Invert: cen_indices_for_pred should be in the same order as stored data
      cen_indices_for_pred <- fit_object$model_info$censored_indices[ord_idx]
    }

    # ── PERFORMANCE: precompute shared structures ──────────────────
    # Everything that depends only on locations (not parameters) is
    # computed ONCE, not n_iter times.

    Y_obs    <- fit_object$data$Y_obs
    X_obs    <- as.matrix(fit_object$data$X_obs)
    locs_obs <- as.matrix(fit_object$data$locs_obs)
    svc_idx  <- fit_object$model_info$svc_indices
    n_obs    <- length(Y_obs)

    has_censoring <- sum(cen_indices_for_pred) > 0

    if (has_censoring) {
      # ── Censored: order non-censored, append censored ──
      cen_pos   <- which(cen_indices_for_pred == 1)
      ncen_pos  <- which(cen_indices_for_pred == 0)
      n_ncen    <- length(ncen_pos)

      Y_ncen    <- Y_obs[ncen_pos]
      X_ncen    <- X_obs[ncen_pos, , drop = FALSE]
      locs_ncen <- locs_obs[ncen_pos, , drop = FALSE]

      ord_nc <- GpGp::order_maxmin(locs_ncen)

      ALL_locs_cen <- rbind(locs_ncen[ord_nc, , drop = FALSE],
                            locs_obs[cen_pos, , drop = FALSE])
      ALL_X_cen    <- rbind(X_ncen[ord_nc, , drop = FALSE],
                            X_obs[cen_pos, , drop = FALSE])
      ALL_Y_cen    <- c(Y_ncen[ord_nc], Y_obs[cen_pos])

      M_use <- min(M, n_obs - 1)
      NN_cen <- GpGp::find_ordered_nn(ALL_locs_cen, M_use)
      dmat_cen <- Rfast::Dist(ALL_locs_cen)
      dmat_sub <- dmat_cen[(n_ncen + 1):n_obs, 1:n_ncen, drop = FALSE]
      for (i in 1:nrow(dmat_sub))
        NN_cen[i + n_ncen, 2:(M_use + 1)] <- order(dmat_sub[i, ])[1:M_use]

      precomp_cen <- list(
        ALL_locs = ALL_locs_cen, ALL_X = ALL_X_cen, ALL_Y = ALL_Y_cen,
        NN = NN_cen, n_ncen = n_ncen, ord_nc = ord_nc,
        cen_pos = cen_pos, ncen_pos = ncen_pos
      )
    }

    # ── Prediction locations setup (shared across iterations) ──
    ord_obs  <- GpGp::order_maxmin(locs_obs)
    ord_pred <- GpGp::order_maxmin(locs_pred)

    ALL_locs_pred <- rbind(locs_obs[ord_obs, , drop = FALSE],
                           locs_pred[ord_pred, , drop = FALSE])
    n_total <- n_obs + n_pred

    M_use <- min(M, n_obs - 1)
    NN_pred <- GpGp::find_ordered_nn(ALL_locs_pred, M_use)
    dmat_pred <- Rfast::Dist(ALL_locs_pred)
    dmat_pred_sub <- dmat_pred[(n_obs + 1):n_total, 1:n_obs, drop = FALSE]
    for (i in 1:nrow(dmat_pred_sub))
      NN_pred[i + n_obs, 2:(M_use + 1)] <- order(dmat_pred_sub[i, ])[1:M_use]

    ALL_X_pred <- rbind(X_obs[ord_obs, , drop = FALSE],
                        X_pred[ord_pred, , drop = FALSE])

    if (is.null(svc_idx)) {
      X_svc_all <- matrix(1, nrow = n_total, ncol = 1)
      n_svc <- 1
    } else {
      n_svc <- length(svc_idx)
      X_svc_all <- ALL_X_pred[, svc_idx, drop = FALSE]
    }

    cat(sprintf("  [Bayesian prediction] Precomputed NN + distances. Looping %d iterations...\n", n_iter))

    # ── Storage ──
    predictions_all <- matrix(NA, n_pred, n_iter)
    if (simulations) simulations_all <- matrix(NA, n_pred, n_iter)
    if (pred_coef)   coef_predictions_all <- array(NA, dim = c(n_pred, n_svc, n_iter))

    t0 <- proc.time()

    for (ii in seq_along(iter_idx)) {
      iter <- iter_idx[ii]

      beta_i  <- as.vector(beta_samples[iter, ])
      sigma_i <- sigma_samples[iter, ]
      phi_i   <- phi_samples[iter, ]
      tau_i   <- tau_samples[iter]

      # ── Step 1: impute censored values if needed ──
      Y_current <- Y_obs
      if (has_censoring) {
        pc <- precomp_cen
        cov_list <- vector("list", n_svc + 1)
        for (j in seq_len(n_svc))
          cov_list[[j]] <- c(sigma_sq = sigma_i[j], phi = phi_i[j])
        cov_list$nugget <- as.numeric(tau_i)

        if (is.null(svc_idx)) {
          X_svc_cen <- matrix(1, nrow = n_obs, ncol = 1)
        } else {
          X_svc_cen <- pc$ALL_X[, svc_idx, drop = FALSE]
        }

        L_cen <- Matrix::t(VecchiaCensored:::precision_decomp_sparse(
          cov_list, X_svc_cen, pc$NN, pc$ALL_locs, pc$n_ncen + 1))

        y_cens_pred <- pc$ALL_X[pc$n_ncen + (1:length(pc$cen_pos)), , drop = FALSE] %*% beta_i -
          forwardsolve(
            Matrix::t(L_cen[(pc$n_ncen + 1):n_obs, (pc$n_ncen + 1):n_obs]),
            Matrix::t(L_cen[1:pc$n_ncen, (pc$n_ncen + 1):n_obs]) %*%
              (pc$ALL_Y[1:pc$n_ncen] - pc$ALL_X[1:pc$n_ncen, , drop = FALSE] %*% beta_i))

        diag_idx <- (pc$n_ncen + 1):n_obs
        y_cens_var <- 1 / sqrt(L_cen[cbind(diag_idx, diag_idx)])

        n_cen <- length(pc$cen_pos)
        truncated <- numeric(n_cen)
        for (j in 1:n_cen) {
          if (y_cens_var[j] > 0) {
            z <- (pc$ALL_Y[pc$n_ncen + j] - y_cens_pred[j]) / sqrt(y_cens_var[j])
            if (is.finite(z)) {
              pnz <- pnorm(z)
              if (pnz > 1e-10) {
                truncated[j] <- y_cens_pred[j] - sqrt(y_cens_var[j]) * dnorm(z) / pnz
              } else {
                truncated[j] <- pc$ALL_Y[pc$n_ncen + j]
              }
            } else {
              truncated[j] <- pc$ALL_Y[pc$n_ncen + j]
            }
          } else {
            truncated[j] <- y_cens_pred[j]
          }
        }
        Y_current[pc$cen_pos] <- truncated
      }

      # ── Step 2: prediction at new locations ──
      cov_list <- vector("list", n_svc + 1)
      for (j in seq_len(n_svc))
        cov_list[[j]] <- c(sigma_sq = sigma_i[j], phi = phi_i[j])
      cov_list$nugget <- as.numeric(tau_i)

      L_pred <- Matrix::t(VecchiaCensored:::precision_decomp_sparse(
        cov_list, X_svc_all, NN_pred, ALL_locs_pred, n_obs + 1))

      pred_i <- ALL_X_pred[(n_obs + 1):n_total, , drop = FALSE] %*% beta_i -
        forwardsolve(
          Matrix::t(L_pred[(n_obs + 1):n_total, (n_obs + 1):n_total]),
          Matrix::t(L_pred[1:n_obs, (n_obs + 1):n_total]) %*%
            (Y_current[ord_obs] - X_obs[ord_obs, , drop = FALSE] %*% beta_i))

      # Un-order predictions
      pred_unord <- numeric(n_pred)
      pred_unord[ord_pred] <- pred_i
      predictions_all[, ii] <- pred_unord

      if (simulations) {
        iid <- rnorm(n_pred)
        sim_i <- pred_i + forwardsolve(
          Matrix::t(L_pred[(n_obs + 1):n_total, (n_obs + 1):n_total]), iid)
        sim_unord <- numeric(n_pred)
        sim_unord[ord_pred] <- sim_i
        simulations_all[, ii] <- sim_unord
      }

      # Progress
      if (ii %% 100 == 0 || ii == n_iter) {
        elapsed <- (proc.time() - t0)[3]
        eta <- elapsed / ii * (n_iter - ii)
        cat(sprintf("  [%d/%d] %.1fs elapsed, ETA %.0fs\n", ii, n_iter, elapsed, eta))
      }
    }

    # ── Aggregate ──
    prediction_mean <- rowMeans(predictions_all)

    out <- list(prediction = prediction_mean)
    if (simulations) out$simulations <- simulations_all
    if (pred_coef) {
      # TODO: coef prediction for Bayesian not yet implemented in optimized version
      warning("pred_coef for Bayesian predictions: using mean parameters only")
    }
    return(out)

  } else {
    stop("Unknown model type: ", model_type)
  }
}
