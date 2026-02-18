#' Single Unit Long-run Treatment Effect Estimation
#' 
#' Estimate time-varying treatment effects for a single unit using the two-step Bernstein expansion method.
#' 
#' @param y Numeric vector of outcome variable
#' @param D Numeric vector of treatment indicator (0/1)
#' @param z Numeric vector of exogenous variable (optional)
#' @param w Numeric vector of confounder variable (optional)
#' @param psi Numeric vector of user-defined variable (replaces time index)
#' @param t_post Time point in post-treatment period to estimate (default: 1)
#' @param previous_estimates Numeric vector of previously estimated treatment effects (for sequential estimation)
#' 
#' @return Estimated treatment effect at the specified time point
#' 
#' @details
#' This function implements the two-step Bernstein expansion method for estimating
#' time-varying treatment effects in the presence of endogeneity and long-run dependence.
#' The method involves:
#' 1. Bernstein polynomial expansion of treatment and outcome processes using psi
#' 2. Construction of common proximate variables
#' 3. Two-stage estimation with bias correction
#' 
#' @examples
#' # Generate synthetic data for a single unit
#' set.seed(123)
#' n <- 100
#' t <- 1:n  # time
#' psi <- 1:n  # user-defined common proximal variable (e.g., time)
#' D <- c(rep(0, 50), rep(1, 50))
#' w <- arima.sim(n = n, model = list(ar = 0.7))
#' z <- rnorm(n)
#' tr <- 10 + 0.5 * sin(2 * t * t / n)
#' y <- 5 + 2 * w + 0.5 * z + tr * D + rnorm(n)
#' 
#' # Estimate treatment effect at first post-treatment period
#' effect <- estimate_single_series(y = y, D = D, z = z, w = w, psi = psi, t_post = 1)
#' print(effect)
#' 
#' @export
estimate_single_series <- function(y, D, z = NULL, w = NULL, psi, t_post = 1, previous_estimates = NULL) {
  
  # Input validation
  if (missing(y) || missing(D) || missing(psi)) {
    stop("y, D, and psi must be specified")
  }
  
  if (length(y) != length(D) || length(y) != length(psi)) {
    stop("y, D, and psi must have the same length")
  }
  
  n <- length(y)
  
  # Handle optional variables
  if (is.null(z)) {
    z <- rep(0, n)
  }
  
  if (is.null(w)) {
    w <- rep(0, n)
  }
  
  if (length(z) != n || length(w) != n) {
    stop("z and w must have the same length as y")
  }
  
  # Standardize psi variable (shift to start from 1)
  t1 <- psi - min(psi) + 1
  t2 <- t1^2
  
  # Determine pre-treatment period
  T0 <- sum(D == 0)
  
  if (T0 == 0) {
    warning("No pre-treatment periods found")
    return(NA)
  }
  
  if (t_post > (n - T0)) {
    warning("t_post exceeds available post-treatment periods")
    t_post <- n - T0
  }
  
  # Extract relevant data
  sample_idx <- 1:(T0 + t_post)
  
  DD <- D[sample_idx]
  yy <- y[sample_idx]
  tt1 <- t1[sample_idx]
  tt2 <- t2[sample_idx]
  zz <- z[sample_idx]
  ww <- w[sample_idx]
  
  # Adjust for previously estimated effects if provided
  if (t_post > 1 && !is.null(previous_estimates)) {
    if (length(previous_estimates) >= (t_post - 1)) {
      treatment_indicator <- c(rep(0, T0), rep(1, t_post - 1), 0)
      estimated_effects <- c(rep(0, T0), previous_estimates[1:(t_post - 1)], 0)
      
      # Ensure length matches
      if (length(treatment_indicator) == length(DD) && 
          length(estimated_effects) == length(yy)) {
        DD <- DD - treatment_indicator
        yy <- yy - estimated_effects
      }
    }
  }
  
  # Step 1: Bernstein expansion for treatment variable
  lm_D <- try(stats::lm(DD ~ tt2 + tt1), silent = TRUE)
  if (inherits(lm_D, "try-error") || any(is.na(stats::coef(lm_D)))) {
    return(NA)
  }
  
  coef_D <- stats::coef(lm_D)
  a_x <- 2 * coef_D[2]  # Quadratic coefficient
  b_x <- coef_D[3]      # Linear coefficient
  
  # Step 2: Bernstein expansion for outcome variable
  lm_y <- try(stats::lm(yy ~ tt2 + tt1), silent = TRUE)
  if (inherits(lm_y, "try-error") || any(is.na(stats::coef(lm_y)))) {
    return(NA)
  }
  
  coef_y <- stats::coef(lm_y)
  a_y <- 2 * coef_y[2]  # Quadratic coefficient
  b_y <- coef_y[3]      # Linear coefficient
  
  # Check for numerical stability
  if (abs(a_x) < 1e-10) {
    return(NA)
  }
  
  # Construct common proximate variable Y_t
  x_t <- 2 * a_x * tt1 + b_x
  Y_t <- (a_y * x_t + a_x * b_y - a_y * b_x) / a_x
  
  # Calculate key quantities for estimation
  n_obs <- length(DD)
  T1 <- T0 + 1  # First treatment period
  
  # Mean differences in Y_t
  mean_Y_t_treated <- mean(Y_t[T1:min(T0 + t_post, n_obs)])
  mean_Y_t_control <- mean(Y_t[1:T0])
  
  # A1: Treatment variable coefficient
  numerator_A1 <- n_obs * sum(DD * Y_t) - sum(DD) * sum(Y_t)
  denominator_A1 <- n_obs * sum(Y_t^2) - sum(Y_t)^2
  
  if (abs(denominator_A1) < 1e-10) {
    return(NA)
  }
  
  A1 <- (mean_Y_t_treated - mean_Y_t_control) * (numerator_A1 / denominator_A1)
  
  # A3: Exogenous variable coefficient
  numerator_A3 <- n_obs * sum(zz * Y_t) - sum(zz) * sum(Y_t)
  A3 <- (mean_Y_t_treated - mean_Y_t_control) * (numerator_A3 / denominator_A1)
  
  # Auxiliary regression
  lm_aux <- try(stats::lm(yy ~ Y_t + zz + DD), silent = TRUE)
  if (inherits(lm_aux, "try-error") || any(is.na(stats::coef(lm_aux)))) {
    return(NA)
  }
  
  coef_aux <- stats::coef(lm_aux)
  beta_y <- coef_aux[2]  # Coefficient of Y_t
  gamma <- coef_aux[3]   # Coefficient of z
  
  # Calculate the new estimator
  lm_y_Yt <- try(stats::lm(yy ~ Y_t), silent = TRUE)
  if (inherits(lm_y_Yt, "try-error")) {
    return(NA)
  }
  
  y_pred <- stats::predict(lm_y_Yt)
  new_estimator <- mean(y_pred[T1:min(T0 + t_post, n_obs)]) - mean(y_pred[1:T0])
  
  # Final treatment effect estimator
  if (abs(A1) < 1e-10) {
    return(NA)
  }
  
  treatment_effect <- as.numeric(
    (new_estimator - beta_y * (mean_Y_t_treated - mean_Y_t_control) - A3 * gamma) / A1
  )
  
  return(treatment_effect)
}