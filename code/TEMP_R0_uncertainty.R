## Title: Prior distributions of mosquito thermal traits #######################
##
## Project: Global zoonoses - spillover of mosquito-borne pathogens
##
## Purpose: Compute prior distributions for thermal trait parameters and provide
##          functions for sampling from these distributions
##
## Contents: 1) Set-up, load in necessary packages and data-sets
##           2) Transform thermal trait parameters into model parameters
##           3) Calculate R0 across all thermal trait parameter samples
##           4) Create visualizations of R0 distributions
##           5) Functions to sample from thermal trait parameter distributions
##
##
## Inputs:  data - data/clean/ThermalTraitSamples.csv
##                 
##
##          code - code/Mordecai2017/mcmc_utils_all.R
##                 code/Mordecai2017/temp_functions_all.R
##
## Outputs: functions:
##          data:
##
## Written and maintained by: Kyle Dahlin, kydahlin@gmail.com
## Initialized February 2023
## _____________________________________________________________________________

# 1) Set-up,load packages, get data, etc. ----

### Load Libraries ----
library(tidyverse)

### Load in data and functions ----
# Load samples of thermal trait parameters
samples.All <- read.csv("data/clean/ThermalTraitSamples.csv")

# Load functions from Mordecai et al., 2017
# This file contains tools for analysis and visualization.
source("code/Mordecai2017/mcmc_utils_all.R")

# This file contains the thermal response functions and their derivatives.
source("code/Mordecai2017/temp_functions_all.R")


### Define names of mosquitoes and pathogens under consideration ----

# Define mosquito species names and codes:
sp_names <- c(
  "Aedes aegypti", "Aedes albopictus", "Culex quinquefasciatus",
  "Anopheles spp."
)
sp_code <- c("AE", "AL", "CQ", "AN")
names(sp_names) <- sp_code

# Define pathogen names
path_names <- c("DENV", "ZIKV", "WNV", "P.fal.")

# Define mosquito-pathogen pairs
## The generated parameter table should include a set of parameters for
## each pair defined here
MosqPathPairs <- tibble(
  Mosquito_species = c(
    "Aedes aegypti", "Aedes albopictus",
    "Aedes aegypti",
    "Culex quinquefasciatus",
    "Anopheles spp."
  ),
  Pathogen = c("DENV", "DENV", "ZIKV", "WNV", "P.fal.")
)

# Keep track of list of trait names
trait_names <- c(
  "GCD", # gonotrophic cycle duration
  "GCR", #
  "a", # biting rate
  "1/MDR", # mosquito development time
  "TFD", # total female fecundity
  "pEA", # the probability a mosquito will survive from hatching to maturation
  "e2a", #
  "p", # the survival probability of an adult mosquito
  "b", # probability of becoming infected
  "c", # probability of becoming infectious
  "EIP" # extrinsic incubation period
  # 'PDR' # parasite development rate = 1/EIP
)

# 2) Transform thermal trait parameters into model parameters ----
# !!! This should be a separate script / set of functions later

# Biting rate: GCR, GCD, a

# Fecundity: 0.5*EFD/a, 0.5*ER*pO/a

# Egg development probability: pEA, EV*pLA

# Larval development rate: MDR

# Adult mortality rate: 1/lf, -log(p)

# Extrinsic incubation rate: PDR

# Vector competence: bc


# 3) Calculate R0 across all thermal trait parameter samples ----

# 4) Create visualizations of R0 distributions ----

# 5) Functions to sample from thermal trait parameter distributions----








