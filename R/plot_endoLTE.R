#' Visualize Endogenous Long-run Treatment Effect Results
#' 
#' Create comprehensive visualizations for endoLTE estimation results.
#' 
#' @param x An object of class 'endoLTE_result' returned by estimate_panel_lte
#' @param type Type of plot to create:
#'   "all_units" - Individual treatment effect trajectories for all units
#'   "event_time" - Event-time specific average effects
#'   "comparison" - Comparison of pooling methods
#'   "distribution" - Distribution of unit average effects
#'   "performance" - Performance metrics
#'   "scatter" - True vs estimated scatter plot
#'   "all" - Create all plots in a 2x3 grid (default)
#' @param ... Additional arguments passed to plotting functions
#' 
#' @return Invisibly returns the input object
#' 
#' @examples
#' \dontrun{
#' # After running estimate_panel_lte
#' plot_endoLTE(result, type = "all")
#' plot_endoLTE(result, type = "comparison")
#' }
#' 
#' @export
plot_endoLTE <- function(x, type = "all", ...) {
  
  if (!inherits(x, "endoLTE_result")) {
    stop("x must be an object of class 'endoLTE_result'")
  }
  
  # Save current par settings
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (type == "all") {
    # Create a 2x3 grid of plots
    graphics::par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
    
    plot_all_units(x, ...)
    plot_event_time(x, ...)
    plot_pooling_comparison(x, ...)
    plot_distribution(x, ...)
    plot_performance(x, ...)
    plot_scatter(x, ...)
    
  } else {
    switch(type,
           "all_units" = plot_all_units(x, ...),
           "event_time" = plot_event_time(x, ...),
           "comparison" = plot_pooling_comparison(x, ...),
           "distribution" = plot_distribution(x, ...),
           "performance" = plot_performance(x, ...),
           "scatter" = plot_scatter(x, ...),
           stop("Unknown plot type: ", type)
    )
  }
  
  invisible(x)
}

# Helper: Plot all individual unit trajectories
plot_all_units <- function(result, main = "Individual Unit Treatment Effects", ...) {
  
  unit_results <- result$unit_results
  
  if (length(unit_results) == 0) {
    message("No unit results to plot")
    return()
  }
  
  # Determine plot limits
  all_estimates <- unlist(lapply(unit_results, function(x) x$estimates))
  valid_estimates <- all_estimates[!is.na(all_estimates)]
  
  if (length(valid_estimates) == 0) {
    message("No valid estimates to plot")
    return()
  }
  
  y_lim <- range(valid_estimates)
  y_lim <- y_lim + c(-0.1, 0.1) * diff(y_lim)
  
  max_event_time <- max(sapply(unit_results, function(x) length(x$estimates)))
  x_lim <- c(1, max_event_time)
  
  # Create empty plot
  graphics::plot(NA, xlim = x_lim, ylim = y_lim,
       xlab = "Event Time (Periods Since Treatment)",
       ylab = "Treatment Effect",
       main = main,
       ...)
  
  # Add grid
  graphics::grid()
  
  # Plot each unit's trajectory
  colors <- grDevices::rainbow(length(unit_results), alpha = 0.6)
  
  for (i in seq_along(unit_results)) {
    unit <- unit_results[[i]]
    estimates <- unit$estimates
    valid_idx <- !is.na(estimates)
    
    if (any(valid_idx)) {
      graphics::lines(which(valid_idx), estimates[valid_idx], 
            col = colors[i], lwd = 1.5, type = "b", pch = 19, cex = 0.5)
    }
  }
  
  # Add legend
  graphics::legend("topright", legend = sprintf("%d units", length(unit_results)),
         bty = "n", cex = 0.8)
}

# Helper: Plot event-time specific averages
plot_event_time <- function(result, main = "Event-Time Specific Average Effects", ...) {
  
  pooled <- result$pooled_results
  
  if (is.null(pooled)) {
    message("No pooled results to plot")
    return()
  }
  
  event_time_data <- pooled$event_time
  valid_idx <- !is.na(event_time_data$estimates)
  
  if (sum(valid_idx) == 0) {
    message("No valid event-time estimates to plot")
    return()
  }
  
  # Determine y limits
  y_vals <- c(event_time_data$ci_lower[valid_idx], 
              event_time_data$ci_upper[valid_idx],
              event_time_data$true_effects[valid_idx])
  y_lim <- range(y_vals, na.rm = TRUE)
  y_lim <- y_lim + c(-0.1, 0.1) * diff(y_lim)
  
  # Create plot
  graphics::plot(event_time_data$event_time[valid_idx], 
       event_time_data$estimates[valid_idx],
       type = "b", pch = 19, col = "blue", lwd = 2,
       ylim = y_lim,
       xlab = "Event Time",
       ylab = "Treatment Effect",
       main = main,
       ...)
  
  # Add grid
  graphics::grid()
  
  # Add confidence intervals
  graphics::arrows(event_time_data$event_time[valid_idx], event_time_data$ci_lower[valid_idx],
         event_time_data$event_time[valid_idx], event_time_data$ci_upper[valid_idx],
         angle = 90, code = 3, length = 0.05, col = "gray")
  
  # Add true effects if available
  if (any(!is.na(event_time_data$true_effects[valid_idx]))) {
    graphics::lines(event_time_data$event_time[valid_idx], 
          event_time_data$true_effects[valid_idx],
          col = "red", lwd = 2, lty = 2)
  }
  
  # Add horizontal line at 0
  graphics::abline(h = 0, lty = 3, col = "darkgray")
  
  # Add legend
  legend_items <- c("Estimated", "95% CI")
  legend_cols <- c("blue", "gray")
  legend_lty <- c(1, 1)
  legend_pch <- c(19, NA)
  
  if (any(!is.na(event_time_data$true_effects[valid_idx]))) {
    legend_items <- c(legend_items, "True")
    legend_cols <- c(legend_cols, "red")
    legend_lty <- c(legend_lty, 2)
    legend_pch <- c(legend_pch, NA)
  }
  
  graphics::legend("topright", legend = legend_items,
         col = legend_cols, lty = legend_lty, pch = legend_pch,
         lwd = c(2, 1, 2), bg = "white")
}

# Helper: Plot comparison of pooling methods
plot_pooling_comparison <- function(result, main = "Comparison of Pooling Methods", ...) {
  
  pooled <- result$pooled_results
  
  if (is.null(pooled)) {
    message("No pooled results to plot")
    return()
  }
  
  method1 <- pooled$unit_average
  method2 <- pooled$comprehensive
  
  # Prepare data for bar plot
  estimates <- c(method1$estimate, method2$estimate)
  ci_lower <- c(method1$ci_lower, method2$ci_lower)
  ci_upper <- c(method1$ci_upper, method2$ci_upper)
  true_values <- c(method1$true, method2$true)
  
  # Determine y limits
  y_vals <- c(ci_lower, ci_upper, true_values, estimates)
  y_vals <- y_vals[!is.na(y_vals)]
  y_lim <- range(y_vals)
  y_lim <- y_lim + c(-0.2, 0.2) * diff(y_lim)
  
  # Create bar plot
  bar_names <- c("Unit-Average\nPooling", "Comprehensive\nPooling")
  bar_colors <- c("lightblue", "lightgreen")
  
  y_pos <- graphics::barplot(estimates, names.arg = bar_names, col = bar_colors,
                       ylim = y_lim,
                       main = main,
                       ylab = "Average Treatment Effect",
                       las = 1,
                       ...)
  
  # Add error bars (confidence intervals)
  graphics::arrows(y_pos, ci_lower, y_pos, ci_upper,
         angle = 90, code = 3, length = 0.1, lwd = 1.5)
  
  # Add true values as points
  if (any(!is.na(true_values))) {
    graphics::points(y_pos, true_values, pch = 19, col = "red", cex = 1.5)
  }
  
  # Add legend
  legend_items <- c("Estimated", "95% CI")
  legend_cols <- c("lightblue", "black")
  legend_pch <- c(15, NA)
  legend_lty <- c(NA, 1)
  
  if (any(!is.na(true_values))) {
    legend_items <- c(legend_items, "True")
    legend_cols <- c(legend_cols, "red")
    legend_pch <- c(legend_pch, 19)
    legend_lty <- c(legend_lty, NA)
  }
  
  graphics::legend("topright", legend = legend_items,
         col = legend_cols, pch = legend_pch, lty = legend_lty,
         pt.cex = c(2, NA, 1.5), bg = "white")
}

# Helper: Plot distribution of unit average effects
plot_distribution <- function(result, main = "Distribution of Unit Average Effects", ...) {
  
  pooled <- result$pooled_results
  
  if (is.null(pooled) || is.null(pooled$all_estimates_by_time)) {
    message("No data for distribution plot")
    return()
  }
  
  # Calculate unit averages
  unit_averages <- sapply(result$unit_results, function(x) {
    mean(x$estimates, na.rm = TRUE)
  })
  
  unit_averages <- unit_averages[!is.na(unit_averages)]
  
  if (length(unit_averages) == 0) {
    message("No unit averages to plot")
    return()
  }
  
  # Create histogram
  hist_breaks <- seq(min(unit_averages), max(unit_averages), length.out = 20)
  
  graphics::hist(unit_averages, breaks = hist_breaks, col = "lightgreen",
       main = main,
       xlab = "Average Treatment Effect per Unit",
       ylab = "Frequency",
       ...)
  
  # Add vertical lines for mean and true average
  mean_est <- mean(unit_averages)
  graphics::abline(v = mean_est, col = "red", lwd = 2, lty = 2)
  
  if (!is.null(pooled$unit_average$true) && !is.na(pooled$unit_average$true)) {
    graphics::abline(v = pooled$unit_average$true, col = "blue", lwd = 2, lty = 3)
  }
  
  # Add text with statistics
  stats_text <- sprintf("Mean: %.3f\nSD: %.3f\nN: %d",
                        mean_est,
                        stats::sd(unit_averages),
                        length(unit_averages))
  
  graphics::text(graphics::par("usr")[1] + diff(graphics::par("usr")[1:2]) * 0.7,
       graphics::par("usr")[3] + diff(graphics::par("usr")[3:4]) * 0.9,
       stats_text, adj = 0, cex = 0.8)
  
  # Add legend
  legend_items <- c("Sample Mean")
  legend_cols <- c("red")
  legend_lty <- c(2)
  
  if (!is.null(pooled$unit_average$true) && !is.na(pooled$unit_average$true)) {
    legend_items <- c(legend_items, "True Average")
    legend_cols <- c(legend_cols, "blue")
    legend_lty <- c(legend_lty, 3)
  }
  
  graphics::legend("topright", legend = legend_items,
         col = legend_cols, lty = legend_lty, lwd = 2, bg = "white")
}

# Helper: Plot performance metrics
plot_performance <- function(result, main = "Performance Metrics", ...) {
  
  performance <- result$performance
  
  if (is.null(performance)) {
    message("No performance metrics to plot")
    return()
  }
  
  # Prepare metrics for bar plot
  metrics <- c(performance$bias, performance$rmse, performance$coverage)
  names(metrics) <- c("Bias", "RMSE", "Coverage")
  
  # Determine bar colors based on performance
  bar_colors <- c(
    "Bias" = ifelse(abs(performance$bias) > 1, "red", "lightblue"),
    "RMSE" = ifelse(performance$rmse > 2, "orange", "lightgreen"),
    "Coverage" = ifelse(abs(performance$coverage - 0.95) > 0.05, "pink", "lightyellow")
  )
  
  # Determine y limits
  y_max <- max(c(metrics, 1.2), na.rm = TRUE)
  
  # Create bar plot
  graphics::barplot(metrics, col = bar_colors, ylim = c(0, y_max),
          main = main,
          ylab = "Value",
          las = 1,
          ...)
  
  # Add reference lines
  graphics::abline(h = 0, lty = 2, col = "gray")
  graphics::abline(h = 0.95, lty = 3, col = "blue", lwd = 1.5)
  
  # Add text with sample size
  if (!is.null(performance$n_estimates)) {
    graphics::text(2, y_max * 0.9, 
         sprintf("N = %d", performance$n_estimates),
         cex = 0.9, font = 2)
  }
}

# Helper: Plot true vs estimated scatter plot
plot_scatter <- function(result, main = "True vs Estimated Effects", ...) {
  
  pooled <- result$pooled_results
  
  if (is.null(pooled) || is.null(pooled$all_estimates_by_time)) {
    message("No data for scatter plot")
    return()
  }
  
  # Extract all estimates and true values
  all_estimates <- unlist(pooled$all_estimates_by_time)
  all_true <- unlist(pooled$all_true_by_time)
  
  # Match lengths and remove NAs
  min_len <- min(length(all_estimates), length(all_true))
  all_estimates <- all_estimates[1:min_len]
  all_true <- all_true[1:min_len]
  
  valid_idx <- !is.na(all_estimates) & !is.na(all_true)
  all_estimates <- all_estimates[valid_idx]
  all_true <- all_true[valid_idx]
  
  if (length(all_estimates) == 0) {
    message("No valid data for scatter plot")
    return()
  }
  
  # Calculate correlation
  corr_val <- stats::cor(all_estimates, all_true)
  
  # Create scatter plot
  graphics::plot(all_true, all_estimates,
       xlab = "True Treatment Effects",
       ylab = "Estimated Treatment Effects",
       main = main,
       pch = 19, col = grDevices::rgb(0, 0, 1, 0.5), cex = 0.8,
       ...)
  
  # Add 45-degree line (perfect fit)
  graphics::abline(0, 1, col = "red", lwd = 2, lty = 2)
  
  # Add regression line
  if (length(all_estimates) > 1) {
    lm_fit <- stats::lm(all_estimates ~ all_true)
    graphics::abline(lm_fit, col = "blue", lwd = 2)
    
    # Calculate R-squared
    r_squared <- summary(lm_fit)$r.squared
    corr_text <- sprintf("Correlation: %.3f\nR-squared: %.3f", corr_val, r_squared)
  } else {
    corr_text <- sprintf("Correlation: %.3f", corr_val)
  }
  
  # Add correlation text
  graphics::text(graphics::par("usr")[1] + diff(graphics::par("usr")[1:2]) * 0.1,
       graphics::par("usr")[3] + diff(graphics::par("usr")[3:4]) * 0.9,
       corr_text, adj = 0, cex = 0.9, font = 2)
  
  # Add legend
  graphics::legend("topleft", legend = c("45° line", "Regression line"),
         col = c("red", "blue"), lty = c(2, 1), lwd = 2, bg = "white")
}