# ============================================================================
# FONCTION D EESTIMATION HARMONISÉE
# ============================================================================
#' Fit a spatial model with Vecchia approximation
#'
#' @param Y_obs Observed response vector
#' @param locs_obs Matrix of spatial locations
#' @param X_obs Design matrix (optional)
#' @param svc_indices Indices for spatially varying coefficients
#' @param censored_indices Binary vector indicating censored observations
#' @param M Number of nearest neighbors
#' @param bayesian Logical, use Bayesian estimation
#' @param chains Number of MCMC chains
#' @param iter_warmup Warmup iterations
#' @param iter_sampling Sampling iterations
#' @param parallel_chains Number of parallel chains
#' @param fixed_range Use fixed range parameter
#' @return A fitted model object
#' @export
fit_model <- function(Y_obs, locs_obs, X_obs = NULL, svc_indices = NULL,
                      censored_indices, M, bayesian = FALSE, chains = NULL,
                      iter_warmup = NULL, iter_sampling = NULL, parallel_chains = NULL,fixed_range=TRUE) {

  # Validation des inputs (gardé identique)
  if (!is.numeric(Y_obs) || !is.vector(Y_obs)) {
    stop("Y_obs must be a numeric vector.")
  }
  n <- length(Y_obs)
  locs_obs <- as.matrix(locs_obs)
  if (!is.matrix(locs_obs) || nrow(locs_obs) != n) {
    stop("locs_obs must be a matrix with the same number of rows as the length of Y_obs.")
  }

  if (!is.null(X_obs)) {
    X_obs <- as.matrix(X_obs)
    if (!is.matrix(X_obs) || nrow(X_obs) != n) {
      stop("X_obs must be a matrix with the same number of rows as the length of Y_obs.")
    }
  } else {
    X_obs <- as.matrix(rep(1, n))
  }
  p <- ncol(X_obs)

  if (!is.null(svc_indices)) {
    if (!is.numeric(svc_indices) || !is.vector(svc_indices)) {
      stop("svc_indices must be a numeric vector.")
    }
    if (any(svc_indices < 1 | svc_indices > p)) {
      stop(paste0("svc_indices must contain only values between 1 and ", p,
                  ", which is the number of columns in X_obs."))
    }
  }

  if (!is.numeric(censored_indices) || !is.vector(censored_indices) || length(censored_indices) != n) {
    stop("censored_indices must be a numeric vector of the same length as Y_obs.")
  }
  if (!all(censored_indices %in% c(0, 1))) {
    stop("censored_indices must contain only 0s and 1s.")
  }

  if (!is.numeric(M) || length(M) != 1) {
    stop("M must be a numeric scalar.")
  }
  if (!is.logical(bayesian) || length(bayesian) != 1) {
    stop("bayesian must be a logical (TRUE or FALSE).")
  }

  # Structure harmonisée de retour
  result <- list(
    parameters = NULL,
    data = NULL,
    model_info = list(
      bayesian = bayesian,
      svc_indices = svc_indices,
      censored_indices = censored_indices,
      M = M,
      n_iterations = NULL,
      model_type = NULL
    ),
    raw_model = NULL  # Pour stocker l'objet original si nécessaire
  )

  if (!bayesian) {
    # ============ MODÈLES FRÉQUENTISTES ============

    if (sum(censored_indices) == 0) {
      # Cas non censuré
      if (is.null(svc_indices)) {
        # ---- Cas GpGp simple ----
        #  source("fonctions_nocensored_freq.R")
        model <- GpGp::fit_model(y = Y_obs, locs = locs_obs, X = X_obs,
                                 covfun_name = "exponential_isotropic",
                                 m_seq = c(10, M), silent = TRUE)

        result$model_info$model_type <- "gpgp"
        result$raw_model <- model
        result$parameters <- list(
          beta = as.vector(model$betahat),
          sigma = model$covparms[1],
          phi = model$covparms[2],
          tau = model$covparms[1] * model$covparms[3]  # sigma^2 * nugget_ratio
        )
        result$data <- list(
          Y_obs = Y_obs,
          X_obs = X_obs,
          locs_obs = locs_obs,
          ordered_indices = 1:n  # GpGp utilise l'ordre original
        )

      } else {
        # ---- SVC fréquentiste non censuré ----
        # source("fonctions_nocensored_freq.R")
        model <- fit_SVC_nocensored_freq_Vecchia(Y_obs = Y_obs, locs_obs = locs_obs,
                                                 X_obs = X_obs, M = M, svc_indices = svc_indices,fixed_range=fixed_range)

        result$model_info$model_type <- "svc_freq"
        result$raw_model <- model

        # Extraction des paramètres
        mean_params <- model$par[1:p]
        cov_params <- model$par[(p + 1):length(model$par)]
        n_svc <- length(svc_indices)

        sigma_vec <- numeric(n_svc)
        phi_vec <- numeric(n_svc)
        for (i in 1:n_svc) {
          idx <- (2 * (i - 1) + 1):(2 * i)
          sigma_vec[i] <- cov_params[idx[1]]
          phi_vec[i] <- cov_params[idx[2]]
        }
        tau <- cov_params[length(cov_params)]

        result$parameters <- list(
          beta = mean_params,
          sigma = sigma_vec,
          phi = phi_vec,
          tau = tau
        )

        # Données ordonnées pour la prédiction
        ord <- GpGp::order_maxmin(locs_obs)
        result$data <- list(
          Y_obs = Y_obs[ord],
          X_obs = X_obs[ord, , drop = FALSE],
          locs_obs = locs_obs[ord, , drop = FALSE],
          ordered_indices = ord
        )
      }

    } else {
      # ---- SVC fréquentiste censuré ----
      # source("fonctions_censored_freq.R")
      #sourceCpp("sparse_start.cpp")

      model <- fit_censored_freq_Vecchia(Y_obs = Y_obs, locs_obs = locs_obs,
                                         X_obs = X_obs, M = M, svc_indices = svc_indices,
                                         censored_indices = censored_indices,fixed_range=fixed_range)

      result$model_info$model_type <- "svc_freq_censored"
      result$raw_model <- model

      # Extraction des paramètres (même logique que non censuré)
      mean_params <- model$par[1:p]
      cov_params <- model$par[(p + 1):length(model$par)]

      if (is.null(svc_indices)) {
        n_svc <- 1
      } else {
        n_svc <- length(svc_indices)
      }

      sigma_vec <- numeric(n_svc)
      phi_vec <- numeric(n_svc)
      for (i in 1:n_svc) {
        idx <- (2 * (i - 1) + 1):(2 * i)
        sigma_vec[i] <- cov_params[idx[1]]
        phi_vec[i] <- cov_params[idx[2]]
      }
      tau <- cov_params[length(cov_params)]

      result$parameters <- list(
        beta = mean_params,
        sigma = sigma_vec,
        phi = phi_vec,
        tau = tau
      )

      # Préparation des données ordonnées (logique des fonctions existantes)
      censored_indices_pos <- which(censored_indices == 1)
      non_censored_indices_pos <- setdiff(1:n, censored_indices_pos)

      locs_n_cen <- locs_obs[non_censored_indices_pos, , drop = FALSE]
      ord <- GpGp::order_maxmin(locs_n_cen)

      # Reconstruction de l'ordre utilisé dans l'estimation
      ordered_n_cen <- non_censored_indices_pos[ord]
      all_ordered_indices <- c(ordered_n_cen, censored_indices_pos)

      result$data <- list(
        Y_obs = Y_obs[all_ordered_indices],
        X_obs = X_obs[all_ordered_indices, , drop = FALSE],
        locs_obs = locs_obs[all_ordered_indices, , drop = FALSE],
        ordered_indices = all_ordered_indices,
        n_non_censored = length(non_censored_indices_pos)
      )
    }

  } else {
    # ============ MODÈLE BAYÉSIEN ============
    # source("functions_bayesian.R")
    model_result <- fit_SVC_nocensored_bayesian_Vecchia(
      Y_obs = Y_obs, locs_obs = locs_obs, X_obs = X_obs, M = M,
      svc_indices = svc_indices, censored_indices = censored_indices,
      chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling,
      parallel_chains = parallel_chains,FIXED=fixed_range
    )

    mcmc_samples <- model_result[[1]]
    data_stan <- model_result[[2]]

    result$model_info$model_type <- "svc_bayesian"
    result$model_info$n_iterations <- iter_sampling * chains
    result$raw_model <- list(mcmc = mcmc_samples, data_stan = data_stan)

    # Extraction des paramètres MCMC
    alpha_samples <- mcmc_samples$draws("alpha", format = "matrix")
    sigma_samples <- mcmc_samples$draws("sigma", format = "matrix")
    #phi_samples <- mcmc_samples$draws("phi", format = "matrix")
    sigma_e_samples <- mcmc_samples$draws("sigma_e", format = "matrix")


    # Handle phi extraction based on fixed_range
    if (fixed_range) {
      # Single phi extracted and replicated to match number of SVCs
      phi_scalar <- mcmc_samples$draws("phi", format = "matrix")  # n_iter × 1
      n_svc <- if (is.null(svc_indices)) 1 else length(svc_indices)

      # Replicate phi across all SVCs: n_iter × n_svc
      phi_samples <- matrix(
        rep(phi_scalar, n_svc),
        nrow = nrow(phi_scalar),
        ncol = n_svc
      )
    } else {
      # Multiple phis (one per SVC)
      phi_samples <- mcmc_samples$draws("phi", format = "matrix")
    }





    result$parameters <- list(
      beta = alpha_samples,      # matrix: iterations × p
      sigma = sigma_samples,     # matrix: iterations × n_svc
      phi = phi_samples,         # matrix: iterations × n_svc
      tau = as.vector(sigma_e_samples)  # vector: iterations
    )

    # Données utilisées (déjà ordonnées dans la fonction bayésienne)
    if (sum(censored_indices) > 0) {
      # Logique de réordonnancement identique à functions_bayesian.R
      censored_indices_pos <- which(censored_indices == 1)
      non_censored_indices_pos <- setdiff(1:n, censored_indices_pos)

      locs_n_cen <- locs_obs[non_censored_indices_pos, , drop = FALSE]
      ord <- GpGp::order_maxmin(locs_n_cen)

      all_ordered_indices <- c(non_censored_indices_pos[ord], censored_indices_pos)

      result$data <- list(
        Y_obs = Y_obs[all_ordered_indices],
        X_obs = X_obs[all_ordered_indices, , drop = FALSE],
        locs_obs = locs_obs[all_ordered_indices, , drop = FALSE],
        ordered_indices = all_ordered_indices,
        n_non_censored = length(non_censored_indices_pos)
      )
    } else {
      ord <- GpGp::order_maxmin(locs_obs)
      result$data <- list(
        Y_obs = Y_obs[ord],
        X_obs = X_obs[ord, , drop = FALSE],
        locs_obs = locs_obs[ord, , drop = FALSE],
        ordered_indices = ord
      )
    }
  }

  class(result) <- c("fit_model", "list")

  return(result)
}

