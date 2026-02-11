#' Panel Data Long-run Treatment Effect Estimation
#' 
#' Main function for estimating heterogeneous treatment effects with endogeneity and long-run dependence in panel data.
#' 
#' @param data A data frame in long format containing panel data
#' @param unit_id Name of the unit identifier variable (character)
#' @param time_var Name of the time variable (character)
#' @param outcome Name of the outcome variable (character)
#' @param treatment Name of the treatment indicator variable (character)
#' @param covariates Optional vector of covariate names (character)
#' @param confounders Optional vector of confounder names (character)
#' @param true_effects Name of the true treatment effect variable, if available (character)
#' @param n_post_periods Number of post-treatment periods to estimate (default: 10)
#' @param bootstrap_iter Number of bootstrap iterations for inference (default: 100)
#' @param parallel Logical indicating whether to use parallel processing (default: FALSE)
#' @param n_cores Number of cores for parallel processing (default: 1)
#' @param seed Random seed for reproducibility
#' 
#' @return An object of class 'endoLTE_result' containing estimation results
#' 
#' @details
#' This function implements the extended method for panel data, which applies the single-unit
#' estimator to each treated unit and then aggregates the results using two pooling strategies:
#' 1. Unit-average pooling: Average within each unit first, then average across units
#' 2. Comprehensive pooling: Pool all estimates together and compute overall average
#' 
#' The function handles staggered adoption designs and provides bootstrap inference.
#' 
#' @examples
#' \dontrun{
#' # Generate synthetic panel data
#' panel_data <- generate_endoLTE_data(N = 30, T = 60, seed = 123)
#' 
#' # Estimate treatment effects
#' result <- estimate_panel_lte(
#'   data = panel_data,
#'   unit_id = "unit_id",
#'   time_var = "time",
#'   outcome = "y",
#'   treatment = "D",
#'   covariates = "z",
#'   confounders = "w",
#'   true_effects = "tr",
#'   n_post_periods = 8,
#'   bootstrap_iter = 50,
#'   parallel = FALSE,
#'   seed = 123
#' )
#' 
#' # View summary
#' summary(result)
#' 
#' # Plot results
#' plot_endoLTE(result)
#' }
#' 
#' @export
estimate_panel_lte <- function(data, unit_id, time_var, outcome, treatment,
                               covariates = NULL, confounders = NULL, 
                               true_effects = NULL,
                               n_post_periods = 10, bootstrap_iter = 100,
                               parallel = FALSE, n_cores = 1, seed = NULL) {
  
  # Set seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Input validation
  required_vars <- c(unit_id, time_var, outcome, treatment)
  if (!all(required_vars %in% names(data))) {
    missing_vars <- setdiff(required_vars, names(data))
    stop("Required variables not found in data: ", paste(missing_vars, collapse = ", "))
  }
  
  # Convert to data frame if not already
  data <- as.data.frame(data)
  
  # Order data by unit and time
  data <- data[order(data[[unit_id]], data[[time_var]]), ]
  
  # Identify treated units
  treated_units <- unique(data[[unit_id]][data[[treatment]] == 1])
  
  if (length(treated_units) == 0) {
    stop("No treated units found in the data")
  }
  
  message("Processing ", length(treated_units), " treated units...")
  
  # Initialize storage for unit results
  unit_results <- list()
  
  # Process each treated unit
  for (i in seq_along(treated_units)) {
    unit <- treated_units[i]
    message(sprintf("  Unit %d/%d: %s", i, length(treated_units), unit))
    
    # Extract unit data
    unit_data <- data[data[[unit_id]] == unit, ]
    
    # Prepare variables
    y <- unit_data[[outcome]]
    D <- unit_data[[treatment]]
    time <- unit_data[[time_var]]
    
    # Extract covariates if provided
    if (!is.null(covariates)) {
      z <- as.matrix(unit_data[, covariates, drop = FALSE])
    } else {
      z <- NULL
    }
    
    # Extract confounders if provided
    if (!is.null(confounders)) {
      w <- as.matrix(unit_data[, confounders, drop = FALSE])
    } else {
      w <- NULL
    }
    
    # Extract true effects if available
    if (!is.null(true_effects) && true_effects %in% names(unit_data)) {
      true_tr <- unit_data[[true_effects]]
    } else {
      true_tr <- NULL
    }
    
    # Estimate treatment effects for this unit
    unit_result <- estimate_unit_effects(
      y = y, D = D, z = z, w = w, time = time, 
      true_effects = true_tr,
      n_post_periods = n_post_periods,
      bootstrap_iter = bootstrap_iter
    )
    
    # Add unit identifier
    unit_result$unit_id <- unit
    
    # Store results
    unit_results[[as.character(unit)]] <- unit_result
  }
  
  # Calculate pooled effects
  pooled_results <- calculate_pooled_effects(unit_results)
  
  # Calculate performance metrics if true effects are available
  performance <- NULL
  if (!is.null(true_effects)) {
    performance <- calculate_performance(pooled_results)
  }
  
  # Create result object
  result <- list(
    unit_results = unit_results,
    pooled_results = pooled_results,
    performance = performance,
    parameters = list(
      unit_id = unit_id,
      time_var = time_var,
      outcome = outcome,
      treatment = treatment,
      covariates = covariates,
      confounders = confounders,
      true_effects = true_effects,
      n_post_periods = n_post_periods,
      bootstrap_iter = bootstrap_iter,
      n_units = length(treated_units)
    ),
    call = match.call()
  )
  
  class(result) <- "endoLTE_result"
  
  message("Estimation completed successfully!")
  
  return(result)
}

# Helper function: Estimate effects for a single unit
estimate_unit_effects <- function(y, D, z = NULL, w = NULL, time, true_effects = NULL,
                                  n_post_periods = 10, bootstrap_iter = 100) {
  
  n <- length(y)
  T0 <- sum(D == 0)
  
  # Determine how many post-treatment periods we can estimate
  max_post <- min(n_post_periods, n - T0)
  
  if (max_post <= 0) {
    return(list(
      estimates = numeric(0),
      true_effects = numeric(0),
      ci_lower = numeric(0),
      ci_upper = numeric(0),
      se = numeric(0),
      event_time = numeric(0)
    ))
  }
  
  # Initialize storage
  estimates <- numeric(max_post)
  true_vec <- if (!is.null(true_effects)) {
    true_effects[(T0 + 1):(T0 + max_post)]
  } else {
    rep(NA, max_post)
  }
  
  ci_lower <- numeric(max_post)
  ci_upper <- numeric(max_post)
  se <- numeric(max_post)
  
  previous_estimates <- NULL
  
  # Estimate for each post-treatment period
  for (t in 1:max_post) {
    # Handle multivariate covariates/confounders
    if (!is.null(z) && is.matrix(z)) {
      z_vec <- if (ncol(z) == 1) z[, 1] else rowMeans(z)
    } else if (!is.null(z)) {
      z_vec <- z
    } else {
      z_vec <- NULL
    }
    
    if (!is.null(w) && is.matrix(w)) {
      w_vec <- if (ncol(w) == 1) w[, 1] else rowMeans(w)
    } else if (!is.null(w)) {
      w_vec <- w
    } else {
      w_vec <- NULL
    }
    
    # Point estimate
    estimates[t] <- estimate_single_series(
      y = y, D = D, z = z_vec, w = w_vec, 
      time = time, t_post = t, 
      previous_estimates = previous_estimates
    )
    
    # Bootstrap inference if estimate is valid
    if (!is.na(estimates[t])) {
      bootstrap_result <- bootstrap_unit(
        y = y, D = D, z = z_vec, w = w_vec, time = time,
        t_post = t, previous_estimates = previous_estimates,
        B = bootstrap_iter
      )
      
      ci_lower[t] <- bootstrap_result$ci_lower
      ci_upper[t] <- bootstrap_result$ci_upper
      se[t] <- bootstrap_result$se
    } else {
      ci_lower[t] <- NA
      ci_upper[t] <- NA
      se[t] <- NA
    }
    
    # Update previous estimates for sequential estimation
    if (!is.na(estimates[t])) {
      previous_estimates <- estimates[1:t]
    }
  }
  
  return(list(
    estimates = estimates,
    true_effects = true_vec,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    se = se,
    event_time = 1:max_post,
    T0 = T0
  ))
}

# Helper function: Bootstrap inference for a single unit
bootstrap_unit <- function(y, D, z = NULL, w = NULL, time, t_post, 
                           previous_estimates = NULL, B = 100) {
  
  n <- length(y)
  T0 <- sum(D == 0)
  
  # Initialize bootstrap storage
  bootstrap_estimates <- numeric(B)
  
  for (b in 1:B) {
    tryCatch({
      # Generate bootstrap sample using residual bootstrap
      # Fit a simple model
      if (is.null(z) && is.null(w)) {
        model <- stats::lm(y ~ D + time)
      } else if (is.null(w)) {
        model <- stats::lm(y ~ D + time + z)
      } else if (is.null(z)) {
        model <- stats::lm(y ~ D + time + w)
      } else {
        model <- stats::lm(y ~ D + time + z + w)
      }
      
      # Generate new outcome with same residuals
      y_boot <- stats::fitted(model) + sample(stats::resid(model), replace = TRUE)
      
      # Re-estimate on bootstrap sample
      bootstrap_estimates[b] <- estimate_single_series(
        y = y_boot, D = D, z = z, w = w, time = time,
        t_post = t_post, previous_estimates = previous_estimates
      )
      
    }, error = function(e) {
      bootstrap_estimates[b] <- NA
    })
  }
  
  # Remove NA values
  valid_estimates <- bootstrap_estimates[!is.na(bootstrap_estimates)]
  
  if (length(valid_estimates) == 0) {
    return(list(
      ci_lower = NA,
      ci_upper = NA,
      se = NA
    ))
  }
  
  # Calculate bootstrap statistics
  se_est <- stats::sd(valid_estimates)
  ci_lower <- stats::quantile(valid_estimates, 0.025, na.rm = TRUE)
  ci_upper <- stats::quantile(valid_estimates, 0.975, na.rm = TRUE)
  
  return(list(
    ci_lower = as.numeric(ci_lower),
    ci_upper = as.numeric(ci_upper),
    se = se_est
  ))
}

# Helper function: Calculate pooled effects
calculate_pooled_effects <- function(unit_results) {
  
  if (length(unit_results) == 0) {
    return(NULL)
  }
  
  # Find maximum event time across all units
  max_event_time <- max(sapply(unit_results, function(x) length(x$estimates)))
  
  # Method 1: Unit-average pooling (average within each unit first)
  unit_averages_est <- numeric(0)
  unit_averages_true <- numeric(0)
  
  for (unit in unit_results) {
    valid_estimates <- unit$estimates[!is.na(unit$estimates)]
    valid_true <- unit$true_effects[!is.na(unit$true_effects)]
    
    if (length(valid_estimates) > 0) {
      unit_averages_est <- c(unit_averages_est, mean(valid_estimates))
    }
    
    if (length(valid_true) > 0) {
      unit_averages_true <- c(unit_averages_true, mean(valid_true))
    }
  }
  
  # Method 1 results
  if (length(unit_averages_est) > 0) {
    method1_estimate <- mean(unit_averages_est)
    method1_se <- stats::sd(unit_averages_est) / sqrt(length(unit_averages_est))
    method1_ci_lower <- method1_estimate - 1.96 * method1_se
    method1_ci_upper <- method1_estimate + 1.96 * method1_se
    method1_true <- if (length(unit_averages_true) > 0) mean(unit_averages_true) else NA
  } else {
    method1_estimate <- method1_se <- method1_ci_lower <- method1_ci_upper <- method1_true <- NA
  }
  
  # Method 2: Comprehensive pooling (pool all estimates together)
  all_estimates <- numeric(0)
  all_true <- numeric(0)
  
  for (unit in unit_results) {
    valid_idx <- !is.na(unit$estimates) & !is.na(unit$true_effects)
    if (any(valid_idx)) {
      all_estimates <- c(all_estimates, unit$estimates[valid_idx])
      all_true <- c(all_true, unit$true_effects[valid_idx])
    }
  }
  
  # Method 2 results
  if (length(all_estimates) > 0) {
    method2_estimate <- mean(all_estimates)
    method2_se <- stats::sd(all_estimates) / sqrt(length(all_estimates))
    method2_ci_lower <- method2_estimate - 1.96 * method2_se
    method2_ci_upper <- method2_estimate + 1.96 * method2_se
    method2_true <- if (length(all_true) > 0) mean(all_true) else NA
  } else {
    method2_estimate <- method2_se <- method2_ci_lower <- method2_ci_upper <- method2_true <- NA
  }
  
  # Event-time specific averages
  event_time_est <- numeric(max_event_time)
  event_time_true <- numeric(max_event_time)
  event_time_se <- numeric(max_event_time)
  event_time_ci_lower <- numeric(max_event_time)
  event_time_ci_upper <- numeric(max_event_time)
  n_units_per_time <- numeric(max_event_time)
  
  # Store all estimates by event time
  all_estimates_by_time <- vector("list", max_event_time)
  all_true_by_time <- vector("list", max_event_time)
  
  for (t in 1:max_event_time) {
    estimates_t <- numeric(0)
    true_t <- numeric(0)
    
    for (unit in unit_results) {
      if (t <= length(unit$estimates)) {
        if (!is.na(unit$estimates[t]) && !is.na(unit$true_effects[t])) {
          estimates_t <- c(estimates_t, unit$estimates[t])
          true_t <- c(true_t, unit$true_effects[t])
        }
      }
    }
    
    all_estimates_by_time[[t]] <- estimates_t
    all_true_by_time[[t]] <- true_t
    
    if (length(estimates_t) > 0) {
      event_time_est[t] <- mean(estimates_t)
      event_time_true[t] <- mean(true_t)
      event_time_se[t] <- stats::sd(estimates_t) / sqrt(length(estimates_t))
      event_time_ci_lower[t] <- event_time_est[t] - 1.96 * event_time_se[t]
      event_time_ci_upper[t] <- event_time_est[t] + 1.96 * event_time_se[t]
      n_units_per_time[t] <- length(estimates_t)
    } else {
      event_time_est[t] <- NA
      event_time_true[t] <- NA
      event_time_se[t] <- NA
      event_time_ci_lower[t] <- NA
      event_time_ci_upper[t] <- NA
      n_units_per_time[t] <- 0
    }
  }
  
  return(list(
    # Unit-average pooling (Method 1)
    unit_average = list(
      estimate = method1_estimate,
      true = method1_true,
      se = method1_se,
      ci_lower = method1_ci_lower,
      ci_upper = method1_ci_upper,
      n_units = length(unit_averages_est)
    ),
    
    # Comprehensive pooling (Method 2)
    comprehensive = list(
      estimate = method2_estimate,
      true = method2_true,
      se = method2_se,
      ci_lower = method2_ci_lower,
      ci_upper = method2_ci_upper,
      n_estimates = length(all_estimates)
    ),
    
    # Event-time specific results
    event_time = list(
      event_time = 1:max_event_time,
      estimates = event_time_est,
      true_effects = event_time_true,
      se = event_time_se,
      ci_lower = event_time_ci_lower,
      ci_upper = event_time_ci_upper,
      n_units_per_time = n_units_per_time
    ),
    
    # Raw data for further analysis
    all_estimates_by_time = all_estimates_by_time,
    all_true_by_time = all_true_by_time
  ))
}

# Helper function: Calculate performance metrics
calculate_performance <- function(pooled_results) {
  
  if (is.null(pooled_results)) {
    return(NULL)
  }
  
  # Extract all estimates and true values
  all_estimates <- unlist(pooled_results$all_estimates_by_time)
  all_true <- unlist(pooled_results$all_true_by_time)
  
  # Remove NA values
  valid_idx <- !is.na(all_estimates) & !is.na(all_true)
  all_estimates <- all_estimates[valid_idx]
  all_true <- all_true[valid_idx]
  
  if (length(all_estimates) == 0) {
    return(list(
      bias = NA,
      mse = NA,
      rmse = NA,
      coverage = NA,
      correlation = NA,
      n_estimates = 0
    ))
  }
  
  # Calculate metrics
  bias <- mean(all_estimates - all_true)
  mse <- mean((all_estimates - all_true)^2)
  rmse <- sqrt(mse)
  
  # Calculate coverage probability
  coverage_flags <- numeric(0)
  event_time_data <- pooled_results$event_time
  
  for (t in 1:length(event_time_data$event_time)) {
    if (!is.na(event_time_data$true_effects[t]) && 
        !is.na(event_time_data$ci_lower[t]) && 
        !is.na(event_time_data$ci_upper[t])) {
      is_covered <- as.numeric(
        event_time_data$true_effects[t] >= event_time_data$ci_lower[t] &&
        event_time_data$true_effects[t] <= event_time_data$ci_upper[t]
      )
      coverage_flags <- c(coverage_flags, is_covered)
    }
  }
  
  coverage <- if (length(coverage_flags) > 0) mean(coverage_flags) else NA
  
  # Calculate correlation
  correlation <- if (length(all_estimates) > 1) {
    stats::cor(all_estimates, all_true)
  } else {
    NA
  }
  
  return(list(
    bias = bias,
    mse = mse,
    rmse = rmse,
    coverage = coverage,
    correlation = correlation,
    n_estimates = length(all_estimates)
  ))
}