## Title: Mosquito Thermal Trait Data processing ###############################
##
## Project: Global zoonoses - spillover of mosquito-borne pathogens
##
## Purpose: Conduct all parts of the analysis
##
## Contents: 0) Set-up, load in necessary packages and data-sets
##           1) Load empirical trait data
##           2) Fit trait thermal performance curves to trait data
##           3) Translate traits into model parameters
##           4) Build data set incorporating all axes of variation
##           5) Calculate model outputs
##           6) Illustrate model outputs (might be done separately)
##           7) Conduct sensitivity analysis
##
##
## Inputs:  
##          
##
##          
##          
##
## Outputs: 
##
## Written and maintained by: Kyle Dahlin, kydahlin@gmail.com
## Initialized March 2023


# 0) Set-up, load in necessary packages and data-sets ---------------------
library(tidyverse)

# 1) Load empirical trait data --------------------------------------------

# # Run this to produce dataset from raw data
# source("code/data-cleaning.R")
# data.in.TPC <- data.Reduced

# # Run this to load pre-processed dataset
# # load dataset for fitting trait TPCs
# data.in.TPC <- read_rds("data/clean/data_for_TPC_fitting.rds")

# 2) Fit trait thermal performance curves to trait data -------------------

# # Run this to generate samples of trait TPC parameters from informed posterior distributions
# source("code/get-thermal-trait-priors.R")

# Run this to load pre-processed data set
data.in.transform <- read_rds("data/clean/TPC_param_samples.rds")

# 3) Translate traits into model parameters -------------------------------

# Define temperature range of study
Temps <- seq(5, 50, length.out = 200)

# Thin samples
thin_size <- 100

source("code/trait-transform.R")

# 4) Build data set incorporating all axes of variation -------------------

# data.in.params <- read_rds("data/clean/parameter_TPCs.rds")
# 
source("code/trait-variation.R")


# 5) Calculate model outputs ----------------------------------------------

# data.in.analysis <- read_rds("data/clean/full_traitset.rds")

source("code/get-outputs.R")

# 6) Illustrate model outputs (might be done separately) ------------------

# !!! figures.R = find a better name

# 7) Conduct sensitivity analysis -----------------------------------------

# sensitivity-analysis.R = conduct sensitivity analyses as described in Shocket 2018 and 2020
