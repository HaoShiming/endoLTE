#' Summarize Endogenous Long-run Treatment Effect Results
#' 
#' Print a comprehensive summary of estimation results from estimate_panel_lte.
#' 
#' @param object An object of class 'endoLTE_result' returned by estimate_panel_lte
#' @param digits Number of digits to display (default: 4)
#' @param ... Additional arguments passed to print method
#' 
#' @return Invisibly returns the input object
#' 
#' @method summary endoLTE_result
#' @export
summary.endoLTE_result <- function(object, digits = 4, ...) {
  
  cat("================================================================\n")
  cat("         Endogenous Long-run Treatment Effect Estimation        \n")
  cat("================================================================\n\n")
  
  # Display call information
  if (!is.null(object$call)) {
    cat("Call:\n")
    print(object$call)
    cat("\n")
  }
  
  # Display basic information
  cat("Data Summary:\n")
  cat("-------------\n")
  cat(sprintf("Number of treated units: %d\n", object$parameters$n_units))
  cat(sprintf("Number of post-treatment periods estimated: %d\n", object$parameters$n_post_periods))
  cat(sprintf("Bootstrap iterations: %d\n", object$parameters$bootstrap_iter))
  cat("\n")
  
  # Display pooled effects
  if (!is.null(object$pooled_results)) {
    cat("Pooled Average Treatment Effects:\n")
    cat("---------------------------------\n")
    
    # Method 1: Unit-average pooling
    method1 <- object$pooled_results$unit_average
    cat("Method 1: Unit-Average Pooling\n")
    cat("  Estimated ATE: ", format(method1$estimate, digits = digits), "\n", sep = "")
    cat("  95% CI: [", format(method1$ci_lower, digits = digits), ", ", 
        format(method1$ci_upper, digits = digits), "]\n", sep = "")
    cat("  Standard Error: ", format(method1$se, digits = digits), "\n", sep = "")
    
    if (!is.na(method1$true)) {
      cat("  True ATE: ", format(method1$true, digits = digits), "\n", sep = "")
      cat("  Bias: ", format(method1$estimate - method1$true, digits = digits), "\n", sep = "")
    }
    cat("  Number of units: ", method1$n_units, "\n\n", sep = "")
    
    # Method 2: Comprehensive pooling
    method2 <- object$pooled_results$comprehensive
    cat("Method 2: Comprehensive Pooling\n")
    cat("  Estimated ATE: ", format(method2$estimate, digits = digits), "\n", sep = "")
    cat("  95% CI: [", format(method2$ci_lower, digits = digits), ", ", 
        format(method2$ci_upper, digits = digits), "]\n", sep = "")
    cat("  Standard Error: ", format(method2$se, digits = digits), "\n", sep = "")
    
    if (!is.na(method2$true)) {
      cat("  True ATE: ", format(method2$true, digits = digits), "\n", sep = "")
      cat("  Bias: ", format(method2$estimate - method2$true, digits = digits), "\n", sep = "")
    }
    cat("  Number of estimates: ", method2$n_estimates, "\n\n", sep = "")
    
    # Comparison
    cat("Comparison of Pooling Methods:\n")
    cat("  Difference (Method1 - Method2): ", 
        format(method1$estimate - method2$estimate, digits = digits), "\n", sep = "")
    
    if (!is.na(method1$true) && !is.na(method2$true)) {
      cat("  True difference: ", format(method1$true - method2$true, digits = digits), "\n", sep = "")
    }
    cat("\n")
    
    # Event-time specific results (first 5 periods)
    event_time <- object$pooled_results$event_time
    cat("Event-Time Specific Estimates (first 5 periods):\n")
    cat("Event Time  Estimate  95% CI               True Effect  N Units\n")
    cat("----------  --------  -------------------  -----------  -------\n")
    
    for (t in 1:min(5, length(event_time$event_time))) {
      if (!is.na(event_time$estimates[t])) {
        ci_text <- sprintf("[%.3f, %.3f]", 
                          event_time$ci_lower[t], 
                          event_time$ci_upper[t])
        
        true_text <- if (!is.na(event_time$true_effects[t])) {
          format(event_time$true_effects[t], digits = 3)
        } else {
          "NA"
        }
        
        cat(sprintf("%10d  %8.3f  %19s  %11s  %7d\n",
                   t, 
                   event_time$estimates[t],
                   ci_text,
                   true_text,
                   event_time$n_units_per_time[t]))
      }
    }
    cat("\n")
  }
  
  # Display performance metrics
  if (!is.null(object$performance)) {
    perf <- object$performance
    cat("Performance Metrics:\n")
    cat("--------------------\n")
    cat(sprintf("Coverage Probability: %.4f (target: 0.95)\n", perf$coverage))
    cat(sprintf("Average Bias: %.4f\n", perf$bias))
    cat(sprintf("Mean Squared Error (MSE): %.4f\n", perf$mse))
    cat(sprintf("Root Mean Squared Error (RMSE): %.4f\n", perf$rmse))
    cat(sprintf("Correlation (True vs Estimated): %.4f\n", perf$correlation))
    cat(sprintf("Number of estimates: %d\n", perf$n_estimates))
    cat("\n")
  }
  
  cat("================================================================\n")
  
  invisible(object)
}

#' Print Endogenous Long-run Treatment Effect Results
#' 
#' Print a concise summary of estimation results.
#' 
#' @param x An object of class 'endoLTE_result' returned by estimate_panel_lte
#' @param ... Additional arguments passed to summary method
#' 
#' @method print endoLTE_result
#' @export
print.endoLTE_result <- function(x, ...) {
  
  cat("Endogenous Long-run Treatment Effect Results\n")
  cat("===========================================\n")
  
  cat(sprintf("Number of treated units: %d\n", x$parameters$n_units))
  
  if (!is.null(x$pooled_results)) {
    method1 <- x$pooled_results$unit_average
    method2 <- x$pooled_results$comprehensive
    
    cat("\nPooled Average Treatment Effects:\n")
    cat(sprintf("  Unit-Average Pooling: %.4f [%.4f, %.4f]\n", 
                method1$estimate, method1$ci_lower, method1$ci_upper))
    cat(sprintf("  Comprehensive Pooling: %.4f [%.4f, %.4f]\n", 
                method2$estimate, method2$ci_lower, method2$ci_upper))
  }
  
  if (!is.null(x$performance)) {
    cat(sprintf("\nCoverage Probability: %.4f\n", x$performance$coverage))
    cat(sprintf("RMSE: %.4f\n", x$performance$rmse))
  }
  
  invisible(x)
}