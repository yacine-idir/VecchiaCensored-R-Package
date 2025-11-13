
#' Make predictions from a fitted spatial model
#'
#' @param fit_object A fitted model object from fit_model()
#' @param locs_pred Matrix of prediction locations
#' @param X_pred Design matrix for prediction locations
#' @param M Number of nearest neighbors for Vecchia approximation
#' @param simulations Logical, generate conditional simulations
#' @param n_simulations Number of simulations (for non-Bayesian models)
#' @param pred_coef Logical, predict spatially varying coefficients (SVC only)
#' @return A list containing predictions and optionally simulations/coefficients
#' @export
prediction <- function(fit_object, locs_pred, X_pred, M, simulations = FALSE, n_simulations = NULL, pred_coef = FALSE) {

  locs_pred <- as.matrix(locs_pred)
  X_pred <- as.matrix(X_pred)
  n_pred <- nrow(locs_pred)

  model_type <- fit_object$model_info$model_type

  if (model_type == "gpgp") {
    # ---- Cas GpGp simple ----
    if (pred_coef) {
      warning("pred_coef is only available for SVC models. Ignoring this option.")
    }

    result <- GpGp::predictions(fit = fit_object$raw_model, locs_pred = locs_pred, X_pred = X_pred)

    if (simulations) {
      if(is.null(n_simulations)) n_simulations <- 1
      sims <- GpGp::cond_sim(fit_object$raw_model, locs_pred = locs_pred, X_pred = X_pred, nsims = n_simulations)
      return(list(
        prediction = result$mean,
        simulations = sims
      ))
    } else {
      return(list(prediction = result))
    }

  } else if (model_type %in% c("svc_freq", "svc_freq_censored")) {
    # ---- Cas SVC fréquentiste ----
    #source("prediction_function.R")

    # Reconstruction de cov_parameters pour prediction_freq
    sigma <- fit_object$parameters$sigma
    phi <- fit_object$parameters$phi
    tau <- fit_object$parameters$tau
    beta <- fit_object$parameters$beta

    result <- prediction_freq(
      beta = beta,
      sigma = sigma,
      phi = phi,
      tau = tau,
      M = M,
      Y_obs = fit_object$data$Y_obs,
      X_obs = fit_object$data$X_obs,
      locs_obs = fit_object$data$locs_obs,
      svc_indices = fit_object$model_info$svc_indices,
      cen_indices = fit_object$model_info$censored_indices,
      X_pred = X_pred,
      locs_pred = locs_pred,
      simulation = simulations,
      pred_coef = pred_coef
    )

    if (simulations && pred_coef) {
      return(list(
        prediction = result$prediction,
        simulations = matrix(result$simulation, ncol = 1),
        coef_prediction = result[[3]]  # coef_prediction is 3rd element
      ))
    } else if (simulations) {
      return(list(
        prediction = result$prediction,
        simulations = matrix(result$simulation, ncol = 1)
      ))
    } else if (pred_coef) {
      return(list(
        prediction = result$prediction,
        coef_prediction = result[[2]]  # coef_prediction is 2nd element when no simulations
      ))
    } else {
      return(list(prediction = result$prediction))
    }

  } else if (model_type == "svc_bayesian") {
    # ---- Cas bayésien ----
    # source("prediction_function.R")

    beta_samples <- fit_object$parameters$beta      # iterations × p
    sigma_samples <- fit_object$parameters$sigma    # iterations × n_svc
    phi_samples <- fit_object$parameters$phi        # iterations × n_svc
    tau_samples <- fit_object$parameters$tau        # vector: iterations

    n_iter <- length(tau_samples)

    # Stockage des résultats
    predictions_all <- matrix(NA, n_pred, n_iter)
    if (simulations) {
      simulations_all <- matrix(NA, n_pred, n_iter)
    }
    if (pred_coef) {
      n_svc <- length(fit_object$model_info$svc_indices)
      coef_predictions_all <- array(NA, dim = c(n_pred, n_svc, n_iter))
    }

    # Prédiction pour chaque itération
    for (iter in 1:n_iter) {
      beta_i <- as.vector(beta_samples[iter, ])
      sigma_i <- sigma_samples[iter, ]
      phi_i <- phi_samples[iter, ]
      tau_i <- tau_samples[iter]

      result_i <- prediction_freq(
        beta = beta_i,
        sigma = sigma_i,
        phi = phi_i,
        tau = tau_i,
        M = M,
        Y_obs = fit_object$data$Y_obs,
        X_obs = fit_object$data$X_obs,
        locs_obs = fit_object$data$locs_obs,
        svc_indices = fit_object$model_info$svc_indices,
        cen_indices = fit_object$model_info$censored_indices,
        X_pred = X_pred,
        locs_pred = locs_pred,
        simulation = simulations,
        pred_coef = pred_coef
      )

      predictions_all[, iter] <- result_i$prediction

      if (simulations) {
        simulations_all[, iter] <- result_i$simulation
      }

      if (pred_coef) {
        # Extract coef_prediction from result_i
        if (simulations) {
          coef_predictions_all[, , iter] <- result_i[[3]]
        } else {
          coef_predictions_all[, , iter] <- result_i[[2]]
        }
      }
    }

    # Moyenne des prédictions
    prediction_mean <- rowMeans(predictions_all)

    if (simulations && pred_coef) {
      coef_prediction_mean <- apply(coef_predictions_all, c(1, 2), mean)
      return(list(
        prediction = prediction_mean,
        simulations = simulations_all,
        coef_prediction = coef_prediction_mean
      ))
    } else if (simulations) {
      return(list(
        prediction = prediction_mean,
        simulations = simulations_all
      ))
    } else if (pred_coef) {
      coef_prediction_mean <- apply(coef_predictions_all, c(1, 2), mean)
      return(list(
        prediction = prediction_mean,
        coef_prediction = coef_prediction_mean
      ))
    } else {
      return(list(prediction = prediction_mean))
    }

  } else {
    stop("Unknown model type: ", model_type)
  }
}
