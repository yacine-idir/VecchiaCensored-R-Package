#' Training dataset for Gaussian Process model
#'
#' A simulated spatial dataset with censored observations for testing
#' standard Gaussian Process models without spatially varying coefficients.
#'
#' @format A data frame with 160 rows and 7 columns:
#' \describe{
#'   \item{X}{X spatial coordinate (1 to 100)}
#'   \item{Y}{Y spatial coordinate (1 to 100)}
#'   \item{X1}{Transformed covariate: sqrt(|Y - 60|)}
#'   \item{X2}{Transformed covariate: sqrt(|X - 20|)}
#'   \item{Z}{True (uncensored) response variable}
#'   \item{censored}{Censored response variable (observations below 20th percentile are censored)}
#'   \item{is_censored}{Binary indicator: 1 if censored, 0 otherwise}
#' }
#'
#' @details
#' This dataset was generated with:
#' \itemize{
#'   \item Grid: 10 x 20 locations (200 total, 160 in training set)
#'   \item Covariance: Exponential with sigma² = 15, phi = 40
#'   \item Nugget: 0.1
#'   \item True model: Z = -5*X1 + 10*X2 + spatial error
#'   \item Censoring: ~20% of observations censored at lower tail
#' }
#'
#' @examples
#' # Load the data
#' data(grille_GP)
#'
#' # Visualize spatial locations
#' plot(grille_GP$X, grille_GP$Y, pch = 19,
#'      col = ifelse(grille_GP$is_censored == 1, "red", "blue"),
#'      main = "Training data (red = censored)")
#'
#' # Fit a simple GP model without censoring
#' fit_gp <- fit_model(
#'   Y_obs = grille_GP$censored,
#'   locs_obs = as.matrix(grille_GP[, c("X", "Y")]),
#'   X_obs = as.matrix(grille_GP[, c("X1", "X2")]),
#'   censored_indices = rep(0, nrow(grille_GP)),
#'   M = 30,
#'   bayesian = FALSE
#' )
#' summary(fit_gp)
#'
#' @source Simulated data for package demonstration
"grille_GP"

#' Test dataset for Gaussian Process model
#'
#' Hold-out test set corresponding to \code{\link{grille_GP}}, containing
#' 20% of the original simulated data for model validation.
#'
#' @format A data frame with 40 rows and 7 columns:
#' \describe{
#'   \item{X}{X spatial coordinate (1 to 100)}
#'   \item{Y}{Y spatial coordinate (1 to 100)}
#'   \item{X1}{Transformed covariate: sqrt(|Y - 60|)}
#'   \item{X2}{Transformed covariate: sqrt(|X - 20|)}
#'   \item{Z}{True (uncensored) response variable}
#'   \item{censored}{Censored response variable}
#'   \item{is_censored}{Binary indicator: 1 if censored, 0 otherwise}
#' }
#'
#' @details
#' This is the test split from the same simulation that produced \code{\link{grille_GP}}.
#' Use this dataset to evaluate prediction performance.
#'
#' @examples
#' # Use with grille_GP for prediction
#' data(grille_GP)
#' data(test_GP)
#'
#' # Fit model on training data
#' fit_gp <- fit_model(
#'   Y_obs = grille_GP$censored,
#'   locs_obs = as.matrix(grille_GP[, c("X", "Y")]),
#'   X_obs = as.matrix(grille_GP[, c("X1", "X2")]),
#'   censored_indices = rep(0, nrow(grille_GP)),
#'   M = 30,
#'   bayesian = FALSE
#' )
#'
#' # Make predictions on test set (when predict function is available)
#' # predictions <- predict(fit_gp,
#' #                        newlocs = as.matrix(test_GP[, c("X", "Y")]),
#' #                        newX = as.matrix(test_GP[, c("X1", "X2")]))
#'
#' @seealso \code{\link{grille_GP}}
#' @source Simulated data for package demonstration
"test_GP"

#' Training dataset for Spatially Varying Coefficients model
#'
#' A simulated spatial dataset with censored observations for testing
#' Spatially Varying Coefficients (SVC) models where the effect of
#' covariates varies smoothly across space.
#'
#' @format A data frame with 400 rows and 7 columns:
#' \describe{
#'   \item{X}{X spatial coordinate (1 to 100)}
#'   \item{Y}{Y spatial coordinate (1 to 100)}
#'   \item{X1}{Transformed covariate: sqrt(|Y - 60|) with spatially varying effect}
#'   \item{X2}{Transformed covariate: sqrt(|X - 20|) with spatially varying effect}
#'   \item{Z}{True (uncensored) response variable}
#'   \item{censored}{Censored response variable (observations below 20th percentile are censored)}
#'   \item{is_censored}{Binary indicator: 1 if censored, 0 otherwise}
#' }
#'
#' @details
#' This dataset was generated with spatially varying coefficients:
#' \itemize{
#'   \item Grid: 20 x 25 locations (500 total, 400 in training set)
#'   \item SVC for X1: Exponential covariance with sigma² = 15, phi = 80
#'   \item SVC for X2: Exponential covariance with sigma² = 15, phi = 15
#'   \item Nugget: 0.1
#'   \item True model: Z = -5*X1 + 10*X2 + SVC error
#'   \item Censoring: ~20% of observations censored at lower tail
#' }
#'
#' @examples
#' # Load the data
#' data(grille_SVC)
#'
#' # Visualize spatial locations
#' plot(grille_SVC$X, grille_SVC$Y, pch = 19,
#'      col = ifelse(grille_SVC$is_censored == 1, "red", "blue"),
#'      main = "SVC Training data (red = censored)")
#'
#' # Fit SVC model with censoring
#' fit_svc <- fit_model(
#'   Y_obs = grille_SVC$censored,
#'   locs_obs = as.matrix(grille_SVC[, c("X", "Y")]),
#'   X_obs = as.matrix(grille_SVC[, c("X1", "X2")]),
#'   svc_indices = c(1, 2),  # Both covariates have spatially varying effects
#'   censored_indices = grille_SVC$is_censored,
#'   M = 30,
#'   bayesian = FALSE
#' )
#' summary(fit_svc)
#'
#' # Fit SVC model with Bayesian approach
#' \dontrun{
#' fit_svc_bayes <- fit_model(
#'   Y_obs = grille_SVC$censored,
#'   locs_obs = as.matrix(grille_SVC[, c("X", "Y")]),
#'   X_obs = as.matrix(grille_SVC[, c("X1", "X2")]),
#'   svc_indices = c(1, 2),
#'   censored_indices = grille_SVC$is_censored,
#'   M = 30,
#'   bayesian = TRUE,
#'   chains = 2,
#'   iter_warmup = 500,
#'   iter_sampling = 500,
#'   parallel_chains = 2
#' )
#' summary(fit_svc_bayes)
#' }
#'
#' @source Simulated data for package demonstration
"grille_SVC"

#' Test dataset for Spatially Varying Coefficients model
#'
#' Hold-out test set corresponding to \code{\link{grille_SVC}}, containing
#' 20% of the original simulated data for model validation.
#'
#' @format A data frame with 100 rows and 7 columns:
#' \describe{
#'   \item{X}{X spatial coordinate (1 to 100)}
#'   \item{Y}{Y spatial coordinate (1 to 100)}
#'   \item{X1}{Transformed covariate: sqrt(|Y - 60|)}
#'   \item{X2}{Transformed covariate: sqrt(|X - 20|)}
#'   \item{Z}{True (uncensored) response variable}
#'   \item{censored}{Censored response variable}
#'   \item{is_censored}{Binary indicator: 1 if censored, 0 otherwise}
#' }
#'
#' @details
#' This is the test split from the same simulation that produced \code{\link{grille_SVC}}.
#' Use this dataset to evaluate SVC model prediction performance.
#'
#' @examples
#' # Use with grille_SVC for prediction
#' data(grille_SVC)
#' data(test_SVC)
#'
#' # Fit SVC model on training data
#' fit_svc <- fit_model(
#'   Y_obs = grille_SVC$censored,
#'   locs_obs = as.matrix(grille_SVC[, c("X", "Y")]),
#'   X_obs = as.matrix(grille_SVC[, c("X1", "X2")]),
#'   svc_indices = c(1, 2),
#'   censored_indices = grille_SVC$is_censored,
#'   M = 30,
#'   bayesian = FALSE
#' )
#'
#' # Make predictions on test set (when predict function is available)
#' # predictions <- predict(fit_svc,
#' #                        newlocs = as.matrix(test_SVC[, c("X", "Y")]),
#' #                        newX = as.matrix(test_SVC[, c("X1", "X2")]))
#'
#' @seealso \code{\link{grille_SVC}}
#' @source Simulated data for package demonstration
"test_SVC"
