rm(list=ls(all=TRUE))

# Required sources
library(DBI)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(rstan)
library(statmod)
library(loo)
library(MASS)
library(emdbook)
library(pracma)
library(parallel)
library(cmdstanr)
library(crfsuite)
library(forcats)

# ====================================================================== #
# Loading data and applying inclusion/exclusion criteria
# - Combine Kappa/Lambda-FLC into a single biomarker
# - Cut-off: 99th percentile for biomarkers
# - Eliminate jumps between LoTs
source("functions/1 - Load and Clear Data.R")
# ====================================================================== #
# Standardising continuous variables (in log scale) and imputing
# their missing data via mean imputation (zero)
source("functions/2 - Standardisation and Imputation.R")
# ====================================================================== #
# Creating longitudinal and categorical data
source("functions/3 - Long and Short Formats.R")
# ====================================================================== #
# Building the joint model in Stan
source("functions/4 - Stan Model.R")
# ====================================================================== #
# Auxiliary functions
# - Calculate the maximum a posteriori (MAP)
# - Call model fitting functions
# - Calculate individual weighted residuals (IWRES)
# - Calculate multiclass classification metrics
# - Calculate the variable importance
source("functions/5 - Auxiliary Functions.R")
# ====================================================================== #


# ====================================================================== #
#     JOINT MODEL OF MULTIPLE LONGITUDINAL AND CATEGORICAL OUTCOMES      #
# ====================================================================== #
# Fitting the joint model
fit <- fit_jm(data = datalot2)

# Posterior summary
fit$fit$summary(variables = c("theta_M","theta_F","sigma2_M","sigma2_F","Omega_M","Omega_F","beta_raw","alpha_raw"))

# Generated quantities from the fitted joint model
gq <- gq_jm(fit = fit$fit, data = datalot2)

# Performance metrics for the joint model
metrics_jm(fit1 = fit$fit, fit2 = gq$fit, data = datalot2)

# Variable importance
vi_jm(fit1 = fit$fit, fit2 = gq$fit, data = datalot2)
# ====================================================================== #
