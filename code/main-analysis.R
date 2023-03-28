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

# # Run this to load pre-processed dataset
# # load dataset for fitting trait TPCs
# data.in.TPC <- read_rds("data/clean/data_for_TPC_fitting.rds")

# 2) Fit trait thermal performance curves to trait data -------------------

# Set parameters for MCMC
n.chains <-3 # 5
n.adapt <- 100 # 5000
n.samps <-100 # 5000

# # Run this to generate samples of trait TPC parameters from informed posterior distributions
# source("code/get-thermal-trait-priors.R")

# # Run this to load pre-processed data set
# data.in.transform <- read_rds("data/clean/TPC_param_samples.rds")

# 3) Translate traits into model parameters -------------------------------

# Define temperature range of study
Temps <- seq(5, 50, by = 0.2) #0.1)

# Thin samples
thin_size <- 100#1000

# source("code/trait-transform.R")

# 4) Build data set incorporating all axes of variation -------------------

## Set resolution for host trait variation
# Host density vector: Number of values to include to consider for vertebrate host density
KH_vec_length <- 100

# Biting tolerance vector: Number of values to consider for biting tolerance
sigmaH_vec_length <- 100

# data.in.params <- read_rds("data/clean/parameter_TPCs.rds")
# 
source("code/trait-variation.R")


# 5) Calculate model outputs ----------------------------------------------

# data.in.analysis <- read_rds("data/clean/full_traitset.rds")

source("code/get-outputs.R")

# 6) Illustrate model outputs (might be done separately) ------------------

source("code/create-figures.R")

# 7) Conduct sensitivity analysis -----------------------------------------

# sensitivity-analysis.R = conduct sensitivity analyses as described in Shocket 2018 and 2020
