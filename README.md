# Extreme-BART: Bayesian Additive Regression Trees for High-Dimensional Extreme Value Analysis

This repository contains the R implementation and simulation scripts to reproduce the results presented in the paper:
" Extreme-BART: Bayesian Additive Regression Trees for High-Dimensional Extreme Value Analysis "


## Overview
`Extreme-BART` is a flexible Bayesian tree-based method for extreme value analysis under the Peaks-Over-Threshold (POT) framework. By combining Generalized Pareto Distribution (GPD) modeling with Bayesian Additive Regression Trees (BART) and Dirichlet Additive Regression Trees (DART), the method captures complex non-linear relationships between covariates and tail parameters while maintaining robustness against high-dimensional noise.

This repository provides the scripts to replicate the simulation studies of the paper:
*   `sim_m1.R`: Simulation code for Model 1 (dynamic GPD scale parameter $\sigma(\mathbf{x})$ and constant shape parameter $\xi$).
*   `sim_m2.R`: Simulation code for Model 2 (both GPD scale $\sigma(\mathbf{x})$ and shape $\xi(\mathbf{x})$ parameters are dynamic non-linear functions of covariates).


## Installation & Prerequisites
To run the simulation scripts, please make sure you have the following R packages installed.


```R

### 1. Install CRAN Packages
# You can install the standard R dependencies directly from CRAN:
install.packages(c("BART", "grf", "dplyr", "ggplot2", "evd", "FNN", "tidyr", "devtools"))

### 2. Install Development Packages from GitHub
# The competitor methods erf (Extremal Random Forests) and gbex (Gradient Boosting for Extremes) should be installed from their respective development repositories using devtools:
# In some R versions or system configurations, installation via
# install_github() may fail. If this occurs, please download
# the package source from the corresponding GitHub repository
# and install it locally.


# Install erf
devtools::install_github("nicolagnecco/erf")

# Install gbex
devtools::install_github("JVelthoen/gbex")
```

## Note on Execution Time

By default, the number of Monte Carlo replications M is set to 100 (to replicate the full results of the paper), which might take significant computational time to run. For a quick test run to verify the scripts, you can edit the script and set M to a smaller value (e.g., M <- 5):
