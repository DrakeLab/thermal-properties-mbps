## Title: Prior distributions of mosquito thermal traits #######################
##
## Project: Global zoonoses - spillover of mosquito-borne pathogens
##
## Purpose: Compute prior distributions for thermal trait parameters and provide
##          functions for sampling from these distributions
##
## Contents: 1) Set-up, load in necessary packages and data-sets
##           2) Instantiate data frame incorporating all axes of variation
##           3) Transform thermal trait parameters into model parameters
##           4) Calculate R0 across all thermal trait parameter samples
##           5) Create visualizations of R0 distributions
##           6) Functions to sample from thermal trait parameter distributions
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

# This file contains functions to compute model outputs like R0 and sensitivities
source("code/output-functions.R")

res_reduce <- function(df, new_length) {
  old_length <- length(unique(df))
  if (old_length < new_length) {
    warning("new_length is larger than old length. Vector will be unchanged")
    new_length <- old_length
  }
  ret <- unique(df)[seq(1, length(unique(df)), length.out = new_length)]
}


# Keep track of list of trait names
output_trait_names <- c("a", "TFD", "EFD", "MDR", "e2a", "b", "c", "PDR", "lf")

# 2) Instantiate data frame incorporating all axes of variation ----
# !!! Move this to separate "initialization" script later

# System: mosquito species x pathogen
MosqPathPairs <- tibble(
  Mosquito_species = c(
    "Aedes aegypti", "Aedes albopictus",
    "Aedes aegypti",
    "Culex quinquefasciatus",
    "Anopheles spp."
  ),
  Pathogen = c("DENV", "DENV", "ZIKV", "WNV", "P.fal.")
)


# Sample number of thermal trait prior
sample_vec <- unique(samples.All$sample_num)

# Temperature
min_Temp <- 10
max_Temp <- 40
res_Temp <- 0.5
Temp_vec <- seq(min_Temp, max_Temp, by = res_Temp)

# Model type (Ross-Macdonald or Chitnis)
Models <- c("RM", "Chitnis")


# Create data frame incorporating all variations for the host
# *Host recruitment rate
lambdaH <- .005
# upper estimate from range for Primate traits: (0.001150685, 0.009624300)

# *Host mortality rate
muH <- 1 / (365 * 20)
# estimate from range for Primate traits: lifespan (8.6, 60) years

# *Host maximum biting tolerance / annoyance threshold (mosquitoes bites per day)
sigmaH_vec <- sort(c(10^seq(-0.25,3.25, by = .0625), 20, 50, Inf))

# *Host carrying capacity
KH_vec_length <- 10
KH_vec <- 10^seq(-2, 5, length.out = KH_vec_length + 6)
KH_vec <-  KH_vec %>%
  c(10^(seq(5,5+log10(2), by = diff(log10(KH_vec), lag = 1)[1]))) %>%
  unique()

## Get host-related pathogen parameters
# *Probability of becoming infected after bite from infectious mosquito
betaH <- 1
# operates as a scaling parameter

# *Host recovery rate
gammaH <- 1 / 5
# plausible estimate for infectious period = (2-14 days)

Host_df <- expand_grid(
  # Host parameters
  lambdaH = lambdaH,
  muH = muH,
  KH = KH_vec,
  sigmaH = sigmaH_vec,
  # Pathogen-specific host parameters
  gammaH = gammaH,
  betaH = betaH
) %>% 
  expand_grid(Models)
  

# 3) Calculate input thermal traits as functions of temperature ----
# Briere function
Briere <- function(q, Tmin, Tmax) {
  function(t) {
    pmax(q * t * (t - Tmin) * (Tmax - t)^(1 / 2), 0, na.rm = TRUE)
  }
}

# Quadratic function
Quadratic <- function(q, Tmin, Tmax) {
  function(t) {
    pmax(-q * (t - Tmin) * (t - Tmax), 0, na.rm = TRUE)
  }
}

# Linear
Linear <- function(q, z) {
  function(t) {
    pmax(-q * t + z, 0, na.RM = FALSE)
  }
}

# Function: designate proper thermal response function
# - output is a function of temperature
get.thermal.response <- function(data_in, Temperature) {
  parms <- dplyr::select(data_in, c, T0, Tm)
  function_type <- dplyr::select(data_in, func)
  
  temp_function <- case_when(
    function_type == "Briere" ~ Briere(parms$c, parms$T0, parms$Tm),
    function_type == "Quadratic" ~ Quadratic(parms$c, parms$T0, parms$Tm),
    function_type == "Linear" ~ Linear(parms$c, parms$T0)
  )
  
  out <- temp_function(Temperature)
}

# For each mosquito species, trait, and sample, get a thermal response curve
TPC_df <- samples.All %>%
  filter(sample_num %in% res_reduce(sample_num, 100)) %>% # to thin samples
  full_join(list(Temperature = Temp_vec), by = character(), copy = TRUE) %>%
  mutate(value = case_when(
    func == "Briere" ~ Briere(c, T0, Tm)(Temperature),
    func == "Quadratic" ~ Quadratic(c, T0, Tm)(Temperature),
    func == "Linear" ~ Linear(c, T0)(Temperature)
  )) %>% 
  dplyr::select(-c("c", "T0", "Tm"))

# 4) Transform thermal trait parameters into model parameters ----

# small constant to avoid INFs
eps <- .Machine$double.eps

# !!! This should be a separate script / set of functions later

# Collect mosquito parameters as functions of temperature (keep track of sample number)
MosqParam_df <- TPC_df %>% 
  dplyr::select(-"func") %>% 
  pivot_wider(id_cols = c("Species", "sample_num", "Temperature"), 
              names_from = trait, values_from = value) %>% 
  # Biting rate (sigmaV): GCR, GCD, a
  mutate(sigmaV = a) %>% 
  # Fecundity (f): 0.5*EFD/a, 0.5*TFD*muV/a, 0.5*ER*pO/a
  mutate(fecundity = 0.5*EFD/(a+eps)) %>% 
  # Egg development probability (deltaL): pEA, EV*pLA, e2a
  mutate(deltaL = e2a) %>% 
  # Larval development rate (rhoL): MDR
  mutate(rhoL = MDR) %>% 
  # Adult mortality rate (muV): 1/lf, -log(p)
  mutate(muV = 1/(lf+eps)) %>% 
  # Extrinsic incubation rate (etaV): PDR
  mutate(etaV = PDR) %>% 
  # Vector competence (betaV): bc
  mutate(betaV = b * c) %>% 
  # Select the relevant model parameters
  dplyr::select(Species, sample_num, Temperature, sigmaV, fecundity, deltaL, 
                rhoL, muV, lf, etaV, betaV)

KL_num <- 3000
MosqParam_df$KL <- KL_num

# 5) Combine data frames
full_df <- full_join(Host_df, MosqParam_df, by = character()) 


# 4) Calculate R0 across all thermal trait parameter samples ----
R0_df <- full_df %>% 
  mutate(V0 = compute.V0(.)) %>% 
  mutate(R0 = compute.R0(.))

# Remove un-needed dataframes from memory
rm("full_df", "Host_df", "MosqParam_df", "TPC_df", "samples.All")

# 5) Create visualizations of R0 distributions ----


# Biggest new development here will be to include 89% highest density interval
# lines

# 6) Functions to sample from thermal trait parameter distributions----








