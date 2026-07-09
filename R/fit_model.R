# ============================================================================
# FONCTION D'ESTIMATION HARMONISÉE
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
                      iter_warmup = NULL, iter_sampling = NULL,
                      parallel_chains = NULL, fixed_range = TRUE) {

  # Validation des inputs
  if (!is.numeric(Y_obs) || !is.vector(Y_obs))
    stop("Y_obs must be a numeric vector.")
  n <- length(Y_obs)
  locs_obs <- as.matrix(locs_obs)
  if (!is.matrix(locs_obs) || nrow(locs_obs) != n)
    stop("locs_obs must be a matrix with the same number of rows as the length of Y_obs.")
  if (!is.null(X_obs)) {
    X_obs <- as.matrix(X_obs)
    if (!is.matrix(X_obs) || nrow(X_obs) != n)
      stop("X_obs must be a matrix with the same number of rows as the length of Y_obs.")
  } else {
    X_obs <- as.matrix(rep(1, n))
  }
  p <- ncol(X_obs)
  if (!is.null(svc_indices)) {
    if (!is.numeric(svc_indices) || !is.vector(svc_indices))
      stop("svc_indices must be a numeric vector.")
    if (any(svc_indices < 1 | svc_indices > p))
      stop(paste0("svc_indices must contain only values between 1 and ", p, "."))
  }
  if (!is.numeric(censored_indices) || !is.vector(censored_indices) || length(censored_indices) != n)
    stop("censored_indices must be a numeric vector of the same length as Y_obs.")
  if (!all(censored_indices %in% c(0, 1)))
    stop("censored_indices must contain only 0s and 1s.")
  if (!is.numeric(M) || length(M) != 1)
    stop("M must be a numeric scalar.")
  if (!is.logical(bayesian) || length(bayesian) != 1)
    stop("bayesian must be a logical (TRUE or FALSE).")

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
    raw_model = NULL
  )

  if (!bayesian) {
    # ============ MODÈLES FRÉQUENTISTES ============

    if (sum(censored_indices) == 0) {
      if (is.null(svc_indices)) {
        # ---- Cas GpGp simple ----
        model <- GpGp::fit_model(y = Y_obs, locs = locs_obs, X = X_obs,
                                 covfun_name = "exponential_isotropic",
                                 m_seq = c(10, M), silent = TRUE)

        result$model_info$model_type <- "gpgp"
        result$raw_model <- model
        result$parameters <- list(
          beta = as.vector(model$betahat),
          sigma = model$covparms[1],
          phi = model$covparms[2],
          tau = model$covparms[1] * model$covparms[3]
        )
        # Bug 1 fix: store data in ORIGINAL order
        result$data <- list(
          Y_obs = Y_obs,
          X_obs = X_obs,
          locs_obs = locs_obs,
          ordered_indices = 1:n
        )

      } else {
        # ---- SVC fréquentiste non censuré ----
        model <- fit_SVC_nocensored_freq_Vecchia(
          Y_obs = Y_obs, locs_obs = locs_obs,
          X_obs = X_obs, M = M, svc_indices = svc_indices,
          fixed_range = fixed_range)

        result$model_info$model_type <- "svc_freq"
        result$raw_model <- model

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

        result$parameters <- list(
          beta = mean_params,
          sigma = sigma_vec,
          phi = phi_vec,
          tau = cov_params[length(cov_params)]
        )
        # Bug 1 fix: store data in ORIGINAL order
        result$data <- list(
          Y_obs = Y_obs,
          X_obs = X_obs,
          locs_obs = locs_obs
        )
      }

    } else {
      # ---- Fréquentiste censuré ----
      model <- fit_censored_freq_Vecchia(
        Y_obs = Y_obs, locs_obs = locs_obs,
        X_obs = X_obs, M = M, svc_indices = svc_indices,
        censored_indices = censored_indices,
        fixed_range = fixed_range)

      result$model_info$model_type <- "svc_freq_censored"
      result$raw_model <- model

      mean_params <- model$par[1:p]
      cov_params <- model$par[(p + 1):length(model$par)]
      n_svc <- if (is.null(svc_indices)) 1 else length(svc_indices)

      sigma_vec <- numeric(n_svc)
      phi_vec <- numeric(n_svc)
      for (i in 1:n_svc) {
        idx <- (2 * (i - 1) + 1):(2 * i)
        sigma_vec[i] <- cov_params[idx[1]]
        phi_vec[i] <- cov_params[idx[2]]
      }

      result$parameters <- list(
        beta = mean_params,
        sigma = sigma_vec,
        phi = phi_vec,
        tau = cov_params[length(cov_params)]
      )
      # Bug 1 fix: store data in ORIGINAL order
      result$data <- list(
        Y_obs = Y_obs,
        X_obs = X_obs,
        locs_obs = locs_obs
      )
    }

  } else {
    # ============ MODÈLE BAYÉSIEN ============
    model_result <- fit_SVC_nocensored_bayesian_Vecchia(
      Y_obs = Y_obs, locs_obs = locs_obs, X_obs = X_obs, M = M,
      svc_indices = svc_indices, censored_indices = censored_indices,
      chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling,
      parallel_chains = parallel_chains, FIXED = fixed_range
    )

    mcmc_samples <- model_result[[1]]
    data_stan <- model_result[[2]]

    result$model_info$model_type <- "svc_bayesian"
    result$model_info$n_iterations <- iter_sampling * chains
    result$raw_model <- list(mcmc = mcmc_samples, data_stan = data_stan)

    # Extract alpha (always named alpha in both old and new Stan)
    alpha_samples <- mcmc_samples$draws("alpha", format = "matrix")

    # Extract sigma: try direct first, then log_sigma (new Stan)
    sigma_samples <- tryCatch(
      mcmc_samples$draws("sigma", format = "matrix"),
      error = function(e) {
        log_sigma <- mcmc_samples$draws("log_sigma", format = "matrix")
        exp(log_sigma)
      })

    # Extract tau: try sigma_e first (old Stan), then log_tau (new Stan)
    tau_samples <- tryCatch(
      as.vector(mcmc_samples$draws("sigma_e", format = "matrix")),
      error = function(e) {
        tryCatch({
          as.vector(exp(mcmc_samples$draws("log_tau", format = "matrix")))
        }, error = function(e2) NULL)
      })

    # Extract phi: try phi first (old Stan), then log_phi (new Stan)
    n_svc <- if (is.null(svc_indices)) 1 else length(svc_indices)
    phi_samples <- tryCatch({
      phi_raw <- mcmc_samples$draws("phi", format = "matrix")
      if (fixed_range && ncol(phi_raw) == 1) {
        matrix(rep(phi_raw, n_svc), nrow = nrow(phi_raw), ncol = n_svc)
      } else {
        phi_raw
      }
    }, error = function(e) {
      tryCatch({
        log_phi_raw <- mcmc_samples$draws("log_phi", format = "matrix")
        phi_raw <- exp(log_phi_raw)
        if (fixed_range && ncol(phi_raw) == 1) {
          matrix(rep(phi_raw, n_svc), nrow = nrow(phi_raw), ncol = n_svc)
        } else {
          phi_raw
        }
      }, error = function(e2) NULL)
    })

    result$parameters <- list(
      beta = alpha_samples,
      sigma = sigma_samples,
      phi = phi_samples,
      tau = tau_samples
    )

    # Bug 1 fix: store data in ORIGINAL order
    result$data <- list(
      Y_obs = Y_obs,
      X_obs = X_obs,
      locs_obs = locs_obs
    )
  }

  class(result) <- c("fit_model", "list")
  return(result)
}
