#' Summary method for fit_model objects
#'
#' Provides customized summary output based on model type
#'
#' @param object A fit_model object returned by fit_model()
#' @param ... Additional arguments (not used)
#' @return Summary output appropriate for the model type
#' @export
summary.fit_model <- function(object, ...) {
  
  if (!inherits(object, "fit_model")) {
    stop("Object must be of class 'fit_model'")
  }
  
  model_type <- object$model_info$model_type
  
  # Check model type and dispatch accordingly
  if (model_type == "gpgp") {
    # GpGp model: return summary of raw_model
    cat("=== GpGp Model Summary ===\n\n")
    return(summary(object$raw_model))
    
  } else if (model_type %in% c("svc_freq", "svc_freq_censored")) {
    # Frequentist (non-GpGp): print parameters
    cat("=== Frequentist Model Summary ===\n")
    cat("Model type:", model_type, "\n\n")
    
    cat("Parameters:\n")
    print(object$parameters)
    
    cat("\n")
    cat("Model Information:\n")
    cat("  Bayesian:", object$model_info$bayesian, "\n")
    cat("  Number of neighbors (M):", object$model_info$M, "\n")
    cat("  Sample size:", length(object$data$Y_obs), "\n")
    if (!is.null(object$model_info$svc_indices)) {
      cat("  SVC indices:", paste(object$model_info$svc_indices, collapse = ", "), "\n")
    }
    if (sum(object$model_info$censored_indices) > 0) {
      cat("  Number of censored observations:", sum(object$model_info$censored_indices), "\n")
    }
    
    return(invisible(object$parameters))
    
  } else if (model_type == "svc_bayesian") {
    # Bayesian: print summary and plot densities
    cat("=== Bayesian Model Summary ===\n\n")
    
    # Print MCMC summary
    cat("MCMC Summary:\n")
    print(object$raw_model$mcmc$summary())
    
    cat("\n")
    cat("Model Information:\n")
    cat("  Number of iterations:", object$model_info$n_iterations, "\n")
    cat("  Number of neighbors (M):", object$model_info$M, "\n")
    cat("  Sample size:", length(object$data$Y_obs), "\n")
    if (!is.null(object$model_info$svc_indices)) {
      cat("  SVC indices:", paste(object$model_info$svc_indices, collapse = ", "), "\n")
    }
    if (sum(object$model_info$censored_indices) > 0) {
      cat("  Number of censored observations:", sum(object$model_info$censored_indices), "\n")
    }
    
    # Plot density plots for all parameters
    cat("\nGenerating density plots for parameters...\n")
    
    # Use bayesplot for density plots
    if (requireNamespace("bayesplot", quietly = TRUE)) {
      dens_plot <- bayesplot::mcmc_dens(object$raw_model$mcmc$draws())
      print(dens_plot)
    } else {
      warning("bayesplot package not available. Install it to see density plots.")
    }
    
    return(invisible(object$raw_model$mcmc$summary()))
    
  } else {
    # Unknown model type
    cat("Unknown model type:", model_type, "\n")
    return(invisible(object))
  }
}
