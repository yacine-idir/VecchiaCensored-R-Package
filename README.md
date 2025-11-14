# VecchiaCensored

An R package for efficient spatial modeling with Vecchia approximation, supporting Gaussian Processes (GP) and Spatially Varying Coefficients (SVC) models with optional censoring and Bayesian inference.

## Features

- **Fast spatial modeling** using Vecchia approximation
- **Flexible model types**:
  - Standard Gaussian Process (GP)
  - Spatially Varying Coefficients (SVC)
- **Censoring support**: Handle left-censored observations
- **Estimation methods**: Both frequentist (MLE) and Bayesian (MCMC)
- **Unified interface**: Just two main functions - `fit_model()` for estimation and `prediction()` for prediction
- **Additional capabilities**:
  - Simulate from fitted models
  - Predict spatially varying coefficients
  - Flexible range parameters (shared or varying)

## Installation
After updating R, you can do it from "https://cran.r-project.org/", you will need Rtools to install from github, install it from "https://cran.r-project.org/bin/windows/Rtools/". Finally, you need the package "cmdstanr", in R :  
install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
 
Install 'VecchiaCensored" from GitHub using:

```r
# Install devtools if needed
install.packages("devtools")

# Install VecchiaCensored
devtools::install_github("yacine-idir/VecchiaCensored")
```

## Quick Start

```r
library(VecchiaCensored)

# Load example data
data(grille_GP)    # Training data for GP models
data(test_GP)      # Test data for GP models
data(grille_SVC)   # Training data for SVC models
data(test_SVC)     # Test data for SVC models
```

## Main Functions

### `fit_model()`
Estimate model parameters from observed data.

**Key arguments:**
- `Y_obs`: Response vector
- `locs_obs`: Matrix of spatial coordinates
- `X_obs`: Covariate matrix
- `censored_indices`: Binary vector (1 = censored, 0 = observed)
- `svc_indices`: Indices of covariates with spatially varying effects
- `M`: Number of nearest neighbors for Vecchia approximation
- `bayesian`: Logical, use Bayesian inference (default: FALSE)
- `fixed_range`: Logical, use same range for all SVCs (default: TRUE)
- `chains`: Number of chains in the MCMC sampling, only when bayesian=TRUE.
- `iter_warmup`: Number of warmup iterations in the MCMC sampling, only when bayesian=TRUE.
- `iter_sampling`: Number of iterations in the MCMC sampling, only when bayesian=TRUE.
- `parallel_chains`: Number of parallel chains in the MCMC sampling, only when bayesian=TRUE.

### `prediction()`
Predict at new locations using a fitted model.

**Key arguments:**
- `model`: Fitted model object from `fit_model()`
- `locs_pred`: Matrix of prediction locations
- `X_pred`: Covariate matrix for predictions
- `M`: Number of nearest neighbors
- `simulations`: Generate posterior samples (Bayesian only)
- `pred_coef`: Predict spatially varying coefficients (SVC only)

## Examples

### 1. Gaussian Process Models

#### Frequentist GP (No Censoring)

```r
# Fit the model
model_freq_no_cen <- fit_model(
  Y_obs = grille_GP$Z,
  locs_obs = grille_GP[, 1:2],
  X_obs = grille_GP[, c("X1", "X2")],
  M = 30,
  censored_indices = rep(0, nrow(grille_GP))
)

# View parameter estimates
summary.fit_model(model_freq_no_cen)

# Predict at new locations
prediction <- prediction(
  model_freq_no_cen,
  locs_pred = test_GP[, 1:2],
  X_pred = test_GP[, c("X1", "X2")],
  M = 30
)

# Evaluate predictions
plot(test_GP$Z, prediction$prediction)
cor(test_GP$Z, prediction$prediction)
```

#### Frequentist GP (With Censoring)

```r
# Fit model with censored observations
model_freq_cen <- fit_model(
  Y_obs = grille_GP$censored,
  locs_obs = grille_GP[, 1:2],
  X_obs = grille_GP[, c("X1", "X2")],
  M = 30,
  censored_indices = grille_GP$is_censored
)

summary.fit_model(model_freq_cen)

# Predict
prediction <- prediction(
  model_freq_cen,
  locs_pred = test_GP[, 1:2],
  X_pred = test_GP[, c("X1", "X2")],
  M = 30
)

# Same accuracy despite 20% censoring!
cor(test_GP$Z, prediction$prediction)
```

#### Bayesian GP (No Censoring)

```r
# Bayesian estimation
model_bay_no_cen <- fit_model(
  Y_obs = grille_GP$Z,
  locs_obs = grille_GP[, 1:2],
  X_obs = grille_GP[, c("X1", "X2")],
  M = 10,
  censored_indices = rep(0, nrow(grille_GP)),
  bayesian = TRUE,
  chains = 3,
  iter_warmup = 100,
  iter_sampling = 300,
  parallel_chains = 3
)

# View posterior summaries and density plots
summary.fit_model(model_bay_no_cen)

# Predict
prediction <- prediction(
  model_bay_no_cen,
  locs_pred = test_GP[, 1:2],
  X_pred = test_GP[, c("X1", "X2")],
  M = 30
)

plot(test_GP$Z, prediction$prediction)
```

#### Bayesian GP (With Censoring)

```r
model_bay_cen <- fit_model(
  Y_obs = grille_GP$censored,
  locs_obs = grille_GP[, 1:2],
  X_obs = grille_GP[, c("X1", "X2")],
  M = 10,
  censored_indices = grille_GP$is_censored,
  bayesian = TRUE,
  chains = 3,
  iter_warmup = 100,
  iter_sampling = 300,
  parallel_chains = 3
)

summary.fit_model(model_bay_cen)

prediction <- prediction(
  model_bay_cen,
  locs_pred = test_GP[, 1:2],
  X_pred = test_GP[, c("X1", "X2")],
  M = 30
)

# Accurate predictions even with censored data!
cor(test_GP$Z, prediction$prediction)
```

### 2. Spatially Varying Coefficients (SVC) Models

#### Frequentist SVC (No Censoring)

```r
# Fit SVC model where both covariates have spatially varying effects
model_freq_no_cen <- fit_model(
  Y_obs = grille_SVC$Z,
  locs_obs = grille_SVC[, 1:2],
  X_obs = grille_SVC[, c("X1", "X2")],
  M = 30,
  censored_indices = rep(0, nrow(grille_SVC)),
  svc_indices = 1:2  # Both covariates vary spatially
)

summary.fit_model(model_freq_no_cen)

# Predict
prediction <- prediction(
  model_freq_no_cen,
  locs_pred = test_SVC[, 1:2],
  X_pred = test_SVC[, c("X1", "X2")],
  M = 30
)

plot(test_SVC$Z, prediction$prediction)
cor(test_SVC$Z, prediction$prediction)
```

#### Bayesian SVC (With Censoring + Advanced Features)

```r
# Bayesian SVC with censoring
model_bay_cen <- fit_model(
  Y_obs = grille_SVC$censored,
  locs_obs = grille_SVC[, 1:2],
  X_obs = grille_SVC[, c("X1", "X2")],
  M = 10,
  censored_indices = grille_SVC$is_censored,
  bayesian = TRUE,
  chains = 3,
  iter_warmup = 100,
  iter_sampling = 300,
  parallel_chains = 3,
  svc_indices = c(1, 2)
)

summary.fit_model(model_bay_cen)

# Prediction with simulations and coefficient predictions
prediction <- prediction(
  model_bay_cen,
  locs_pred = test_SVC[, 1:2],
  X_pred = test_SVC[, c("X1", "X2")],
  M = 30,
  simulations = TRUE,    # Generate posterior samples
  pred_coef = TRUE       # Predict spatially varying coefficients
)

# Access results
prediction$prediction           # Point predictions
prediction$coef_prediction     # Predicted SVC coefficients
dim(prediction$simulations)    # n_test × n_iterations matrix
```

## Comparison with Other Packages

VecchiaCensored offers significant speed advantages over exact methods:

```r
# Compare with varycoef (exact, no approximation)
library(varycoef)
fit_vc <- SVC_mle(Z ~ X1 + X2 - 1, 
                  data = grille_SVC, 
                  locs = grille_SVC[, 1:2])

vc_pred <- predict(fit_vc, 
                   newlocs = as.matrix(test_SVC[, 1:2]),
                   newX = as.matrix(test_SVC[, c("X1", "X2")]), 
                   newW = as.matrix(test_SVC[, c("X1", "X2")]))

# VecchiaCensored achieves similar accuracy but faster!
cor(test_SVC$Z, vc_pred$y.pred)        # varycoef
cor(test_SVC$Z, prediction$prediction)  # VecchiaCensored (faster!)
```

## Key Advantages

1. **Speed**: Vecchia approximation enables efficient computation for large datasets
2. **Flexibility**: Handle censoring, SVC, and mixed models with a single interface
3. **Bayesian inference**: Full posterior distributions with MCMC
4. **User-friendly**: Only two main functions to learn
5. **Comprehensive**: Supports prediction, simulation, and coefficient estimation

## Package Data

The package includes four example datasets:

- `grille_GP`: Training data for GP models (160 observations, 20% censored)
- `test_GP`: Test data for GP models (40 observations)
- `grille_SVC`: Training data for SVC models (400 observations, 20% censored)
- `test_SVC`: Test data for SVC models (100 observations)

All datasets include spatial coordinates, covariates, true responses, and censoring indicators.

## Dependencies

Main dependencies:
- `GpGp`: Core Vecchia approximation functionality
- `Matrix`: Fast matrix operations
- `cmdstanr`: Bayesian inference 
- `bayesplot`: Visualization of Bayesian results 

## Citation

If you use this package in your research, please cite:

```
Paper refrence comming soon
```

## License

[Your license here, e.g., MIT, GPL-3]

## Contact

For questions, issues, or contributions, please [open an issue](https://github.com/yacine-idir/VecchiaCensored/issues) on GitHub.

## Acknowledgments

This package builds upon the excellent work of the `GpGp` package for Vecchia approximation.
