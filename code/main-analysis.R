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
# # data-cleaning.R = produce trait dataset in readable form
# source("code/data-cleaning.R")
# data.in.TPC <- data.Reduced

# Run this to load pre-processed dataset
# load dataset for fitting trait TPCs
data.in.TPC <- read_csv("data/clean/data_for_TPC_fitting.csv")

# 2) Fit trait thermal performance curves to trait data -------------------

# get-thermal-trait-priors.R = calculate distributions of thermal trait parameters

# 3) Translate traits into model parameters -------------------------------

# trait-transform.R = perform necessary transformations to traits to get model parameters

# 4) Build data set incorporating all axes of variation -------------------

# get-analysis-dfs.R = produce data.frames incorporating all axes of variation

# 5) Calculate model outputs ----------------------------------------------

# !!! get-analysis-dfs.R = compute model outputs on data.frames from above. needs to be separate script

# 6) Illustrate model outputs (might be done separately) ------------------

# !!! figures.R = find a better name

# 7) Conduct sensitivity analysis -----------------------------------------

# sensitivity-analysis.R = conduct sensitivity analyses as described in Shocket 2018 and 2020
