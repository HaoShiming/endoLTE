# endoLTE
An R package for estimating time-varying treatment effects with endogeneity and long range dependency in untreated units absent time-series natural experiments. Implements robust methods for heterogeneous treatment effect estimation with staggered panel treatments adoption designs.

This package implements in R the social policy effects evaluation method by

Shiming Hao (2026): "Estimating Heterogeneous Treatments Effects with Endogeneity and Long Range Dependency in Untreated Units Absent Time-Series Natural Experiments" (Revise & Resubmit)

*This package is under active development...*

The current beta version can be installed from Github by:

```r
library(devtools)
devtools::install_github("HaoShiming/endoLTE", INSTALL_opts=c("--no-multiarch"))
library(endoLTE)
```

Examples:
```r
# endoLTE Package Usage Example
# ==============================

# Install and load package
# devtools::install_github("HaoShiming/endoLTE")
library(endoLTE)

# Example 1: Basic usage with example data
# -----------------------------------------
custom_data <- generate_endoLTE_data(
  N = 100,
  T = 80,
  H = 0.8,
  treatment_frac = 0.4,
  staggered = TRUE,
  effect_size = 15,
  effect_heterogeneity = 3,
  endogeneity_strength = 0.5,
  seed = 456
)

# Estimate with parallel processing
result2 <- estimate_panel_lte(
  data = custom_data,
  unit_id = "unit_id",
  time_var = "time",
  outcome = "y",
  treatment = "D",
  covariates = "z",
  confounders = "w",
  n_post_periods = 10,
  bootstrap_iter = 200,
  parallel = TRUE,
  n_cores = 4,
  seed = 789
)

# Extract specific results
unit_avg_ate <- result2$pooled_results$unit_average$estimate
comp_ate <- result2$pooled_results$comprehensive$estimate

cat(sprintf("Unit-average ATE: %.3f\n", unit_avg_ate))
cat(sprintf("Comprehensive ATE: %.3f\n", comp_ate))

# Performance metrics
if (!is.null(result2$performance)) {
  cat(sprintf("Coverage: %.3f\n", result2$performance$coverage))
  cat(sprintf("RMSE: %.3f\n", result2$performance$rmse))
  cat(sprintf("Correlation: %.3f\n", result2$performance$correlation))
}

# Example 3: Using your own data
# -------------------------------
# my_data <- read.csv("your_data.csv")
# 
# result3 <- estimate_panel_lte(
#   data = my_data,
#   unit_id = "id",
#   time_var = "period",
#   outcome = "outcome",
#   treatment = "treatment",
#   covariates = c("x1", "x2"),
#   confounders = "c1",
#   n_post_periods = 5,
#   bootstrap_iter = 1000
# )

# Save results
saveRDS(result1, file = "endoLTE_results.rds")

# Load results
loaded_results <- readRDS("endoLTE_results.rds")
summary(loaded_results)
```
