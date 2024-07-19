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


