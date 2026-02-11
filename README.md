# endoLTE
An R package for estimating time-varying treatment effects with endogeneity and long range dependency in untreated units absent time-series natural experiments. Implements robust methods for heterogeneous treatment effect estimation with staggered panel treatments adoption designs.

This package implements in R the social policy effects evaluation method by

Shiming Hao (2026): "Estimating Heterogeneous Treatments Effects with Endogeneity and Long Range Dependency in Untreated Units Absent Time-Series Natural Experiments" (Revise & Resubmission)

This package is under active development...

The current beta version can be installed from Github by:

```r
library(devtools)
devtools::install_github("HaoShiming/endoLTE", INSTALL_opts=c("--no-multiarch"))
library(endoLTE)

