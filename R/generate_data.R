#' Generate Synthetic Panel Data with Endogenous Treatment Effects
#' 
#' Generate realistic panel data with endogenous treatment assignment and long-run dependence.
#' 
#' @param N Number of units (default: 30)
#' @param T Number of time periods per unit (default: 60)
#' @param H Hurst parameter for long-range dependence (default: 0.7)
#' @param treatment_frac Fraction of units that receive treatment (default: 0.3)
#' @param staggered Logical indicating whether treatment adoption is staggered (default: TRUE)
#' @param effect_size Average treatment effect size (default: 10)
#' @param effect_heterogeneity Degree of treatment effect heterogeneity (default: 2)
#' @param endogeneity_strength Strength of endogeneity in treatment assignment (default: 0.3)
#' @param seed Random seed for reproducibility
#' 
#' @return A data frame in long format with synthetic panel data
#' 
#' @details
#' Generates panel data with the following features:
#' 1. Unit-specific fixed effects
#' 2. Time fixed effects (common shocks)
#' 3. Long-run dependent confounders (fractional Brownian motion)
#' 4. Endogenous treatment assignment
#' 5. Time-varying heterogeneous treatment effects
#' 6. Staggered adoption design (optional)
#' 
#' @examples
#' # Generate a small dataset for testing
#' test_data <- generate_endoLTE_data(N = 10, T = 20, seed = 123)
#' head(test_data)
#' 
#' # Generate a larger dataset with staggered adoption
#' full_data <- generate_endoLTE_data(N = 50, T = 100, staggered = TRUE, seed = 456)
#' 
#' @export
generate_endoLTE_data <- function(N = 30, T = 60, H = 0.7, treatment_frac = 0.3,
                                  staggered = TRUE, effect_size = 10,
                                  effect_heterogeneity = 2, endogeneity_strength = 0.3,
                                  seed = NULL) {
  
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Validate inputs
  if (N <= 0 || T <= 0) {
    stop("N and T must be positive integers")
  }
  
  if (H <= 0 || H >= 1) {
    stop("Hurst parameter H must be between 0 and 1")
  }
  
  if (treatment_frac < 0 || treatment_frac > 1) {
    stop("treatment_frac must be between 0 and 1")
  }
  
  # Generate unit-specific fixed effects
  unit_fe <- stats::rnorm(N, mean = 0, sd = 1)
  
  # Generate time fixed effects (common shocks)
  time_fe <- 0.5 * sin(seq(0, 4 * pi, length.out = T)) + stats::rnorm(T, 0, 0.1)
  
  # Determine which units are treated
  n_treated <- max(1, round(N * treatment_frac))
  treated_units <- sample(1:N, n_treated)
  
  # Determine treatment timing
  if (staggered) {
    # Staggered adoption: different start times
    treatment_starts <- sample(round(T * 0.3):round(T * 0.7), n_treated, replace = TRUE)
  } else {
    # Simultaneous adoption: all start at the same time
    treatment_starts <- rep(round(T * 0.5), n_treated)
  }
  
  # Initialize data frame
  panel_data <- data.frame()
  
  # Generate data for each unit
  for (i in 1:N) {
    # Generate unit-specific long-run dependent process
    n_points <- 1024
    fbm_process <- fbm(H = H, n = n_points, T = 1, seed = ifelse(is.null(seed), NULL, seed + i))
    fbm_process <- fbm_process[1:T]
    
    # Unit fixed effect
    alpha_i <- unit_fe[i]
    
    # Generate confounders with long-run dependence
    w_it <- stats::runif(T) + fbm_process + stats::rnorm(T, 0, 0.1)
    
    # Generate exogenous variable
    z_it <- stats::runif(T, -0.01, 0.01) + stats::rnorm(T, 0, 0.05)
    
    # Determine treatment status
    if (i %in% treated_units) {
      idx <- which(treated_units == i)
      treatment_start <- treatment_starts[idx]
      
      # Treatment indicator (1 from treatment_start onward)
      D_it <- ifelse(1:T >= treatment_start, 1, 0)
      
      # Generate time-varying treatment effects
      base_effect <- stats::runif(T, effect_size * 0.5, effect_size * 1.5)
      unit_heterogeneity <- stats::rnorm(1, 0, effect_heterogeneity)
      
      # Add time pattern to treatment effects
      time_pattern <- 0.5 * sin(2 * pi * (1:T) / T)
      tr_it <- base_effect + unit_heterogeneity + time_pattern
    } else {
      # Control unit: never treated
      D_it <- rep(0, T)
      tr_it <- rep(0, T)
      treatment_start <- NA
    }
    
    # Generate outcome with endogeneity
    # Treatment assignment depends on confounders (endogeneity)
    propensity <- plogis(endogeneity_strength * scale(w_it) + stats::rnorm(T, 0, 0.5))
    D_it_endo <- ifelse(propensity > 0.5, 1, 0)
    
    # Ensure treatment starts at the right time for treated units
    if (i %in% treated_units) {
      D_it <- ifelse(1:T >= treatment_start, D_it_endo, 0)
    }
    
    # Generate error term with unit-specific variance
    sigma_i <- stats::runif(1, 0.5, 1.5)
    e_it <- sigma_i * stats::rnorm(T, 0, 0.5)
    
    # Generate outcome variable
    y_it <- alpha_i +                    # Unit fixed effect
      time_fe +                    # Time fixed effect
      5.0 * w_it +                 # Confounder effect
      0.1 * z_it +                 # Exogenous variable effect
      tr_it * D_it +               # Treatment effect (only when treated)
      e_it                         # Error term
    
    # Create unit data frame
    unit_df <- data.frame(
      unit_id = i,
      time = 1:T,
      is_treated = as.numeric(i %in% treated_units),
      treatment_start = ifelse(i %in% treated_units, treatment_start, NA),
      D = D_it,
      y = y_it,
      w = w_it,
      z = z_it,
      tr = tr_it * D_it,  # True treatment effect (0 when not treated)
      e = e_it
    )
    
    panel_data <- rbind(panel_data, unit_df)
  }
  
  # Add time fixed effects as a column
  time_fe_df <- data.frame(time = 1:T, time_fe = time_fe)
  panel_data <- merge(panel_data, time_fe_df, by = "time")
  
  # Reorder columns
  panel_data <- panel_data[, c("unit_id", "time", "is_treated", "treatment_start", 
                               "D", "y", "w", "z", "tr", "e", "time_fe")]
  
  # Sort by unit_id and time
  panel_data <- panel_data[order(panel_data$unit_id, panel_data$time), ]
  
  # Reset row names
  rownames(panel_data) <- NULL
  
  message(sprintf("Generated panel data with %d units, %d time periods", N, T))
  message(sprintf("%d treated units (%.1f%%)", n_treated, 100 * n_treated/N))
  
  return(panel_data)
}