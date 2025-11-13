#' Print method for fit_model objects
#'
#' Provides concise display of fit_model objects
#'
#' @param x A fit_model object returned by fit_model()
#' @param ... Additional arguments (not used)
#' @return Invisibly returns the input object
#' @export
print.fit_model <- function(x, ...) {
  
  cat("Fitted Vecchia Model\n")
  cat("====================\n\n")
  
  cat("Model type:", x$model_info$model_type, "\n")
  cat("Bayesian:", x$model_info$bayesian, "\n")
  cat("Sample size:", length(x$data$Y_obs), "\n")
  cat("Number of neighbors (M):", x$model_info$M, "\n")
  
  if (!is.null(x$model_info$svc_indices)) {
    cat("SVC indices:", paste(x$model_info$svc_indices, collapse = ", "), "\n")
  }
  
  if (sum(x$model_info$censored_indices) > 0) {
    cat("Censored observations:", sum(x$model_info$censored_indices), "\n")
  }
  
  if (x$model_info$bayesian) {
    cat("MCMC iterations:", x$model_info$n_iterations, "\n")
  }
  
  cat("\nUse summary() for detailed results\n")
  
  invisible(x)
}
