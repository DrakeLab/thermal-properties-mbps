## For questions, contact Kyle Dahlin, kydahlin@gmail.com
## Originated September 2021
##
## Title: Analysis code for "Thermal optima" manuscript ########################
##
## Project: Global zoonoses - spillover of mosquito-borne pathogens
##
## Purpose: Build data frames for: 1) vector traits, 2) all traits and model
##          outputs, and 3) thermal characteristics of transmission
##
## Contents: 0) Load in necessary packages, functions, settings
##           1) Build vector trait data frame
##           2) Build all traits and model outputs data frame
##           3) Build thermal characteristics data frame
##
##
## Settings:  Options for reducing the resolution of variables for memory
##            allocation and figure plotting purposes
##
## Inputs:  data - data/clean/trait_transforms.rds
##
## Outputs: (in ./results)
##          1) VectorTraits.csv - vector trait data frame
##          2) AllOutputs.csv - all traits and model outputs data frame
##          3) ThermalCharacteristics.csv - thermal characteristics data frame
##
## Written and maintained by: Kyle Dahlin, kydahlin@gmail.com
## Initialized March 2023
# _______________________________________________________________________________

# 0) Load in necessary packages, functions, settings ###########################
# Packages
require(tidyverse)

# Functions for computing transmission measures
source("code/output-functions.R")

## Settings---------------------------------------------------------------------

# Resolution settings: (set these to Inf to keep original resolution)
# Temperature: Factor by which to reduce resolution
Temp_vec_length <- 500

# Host density vector: Number of values to include to consider for vertebrate host density
KH_vec_length <- 100

# Biting tolerance vector: Number of values to consider for biting tolerance
sigmaH_vec_length <- 100

# 1) Build vector trait data frame #############################################
## Assign names to mosquitoes and pathogens
# Define mosquito species names and codes:
mosquito_names <- c(
  "Aedes aegypti", "Aedes albopictus",
  "Culex quinquefasciatus", "Anopheles spp."
)

# Define pathogen names and codes
pathogen_names <- c("DENV", "ZIKV", "WNV", "P.fal.")

## Build table of mosquito parameters-------------------------------------------
## Get vector (mosquito) parameters, which are functions of temperature

## NB: Carrying capacity for larval mosquitoes
# In the absence of good estimates for each species or temperature-dependence of
# this trait, we assume that this parameter is constant. It can be used as a 
# scaling parameter for overall mosquito abundance 
# (it could alternately be used to fix the maximum adult mosquito density across species)
larval_mosquito_carrying_capacity <- 300

# Get mosquito, pathogen trait thermal performance curve (TPC) data from:
#   "results/MosquitoThermalResponse.csv"
VectorTraits_df <- read_csv("results/MosquitoThermalResponse.csv") %>%
  # Lower the resolution of temperature to length Temp_vec_length
  filter(Temperature %in% res_reduce(Temperature, Temp_vec_length)) %>%
  # Larval mosquito carrying capacity
  mutate(KL = larval_mosquito_carrying_capacity) %>%
  # Average adult lifespan
  mutate(lf = 1 / muV) %>%
  # Compute equilibrium mosquito abundance (V0) from other traits
  mutate(V0 = compute.V0(.)) %>%
  # Combine mosquito species and pathogen name into a single column
  unite(
    col = "system_ID",
    c("Mosquito_species", "Pathogen"),
    sep = "+",
    remove = TRUE
  )

# Save data frame---------------------------------------------------------------
write_csv(VectorTraits_df, "results/VectorTraits.csv")

# 2) Build all traits and model outputs data frame #############################

## Define vertebrate host parameters---------------------------------------------

# *Host recruitment rate
lambdaH_baseline <- .005
# upper estimate from range for Primate traits: (0.001150685, 0.009624300)

# *Host mortality rate
muH_vec <- 1 / (365 * c(1, 25))
muH_baseline <- 1 / (365 * 20)
# estimate from range for Primate traits: lifespan (8.6, 60) years

# *Host maximum biting tolerance / annoyance threshold (mosquitoes bites per day)
sigmaH_vec <- sort(c(10^seq(-0.25,3.25, by = .0625), 20, 50, Inf))
sigmaH_baseline <- 100

# *Host carrying capacity
KH_vec <- 10^seq(-2, 5, length.out = KH_vec_length + 6)
KH_vec <-  KH_vec %>%
  c(10^(seq(5,5+log10(2), by = diff(log10(KH_vec), lag = 1)[1]))) %>%
  unique()
            
## Get host-related pathogen parameters
# *Probability of becoming infected after bite from infectious mosquito
betaH_vec <- c(.25, .75)
betaH_baseline <- 1
# operates as a scaling parameter

# *Host recovery rate
gammaH_vec <- 1 / c(5, 14)
gammaH_baseline <- 1 / 5
# plausible estimate for infectious period = (2-14 days)

Host_df <- expand_grid(
  # Host parameters
  lambdaH = lambdaH_baseline,
  muH = muH_baseline,
  KH = KH_vec,
  sigmaH = sigmaH_vec,
  # Pathogen-specific host parameters
  gammaH = gammaH_baseline,
  betaH = betaH_baseline
)

## Combine all parameter combinations into a large table
AllOutputs_df <- expand_grid(VectorTraits_df, Host_df) %>%
  # Lower the resolution of temperature to length Temp_vec_length
  filter(Temperature %in% res_reduce(Temperature, Temp_vec_length)) %>%
  ## Compute outputs ##
  # Name the model (just in case this is handier than referring to sigmaH)
  mutate(Model = ifelse(is.infinite(sigmaH), "Ross-Macdonald model", "Chitnis model")) %>%
  # Basic reproduction number
  mutate(R0 = sqrt(compute.RH(.) * compute.RV(.))) %>%
  # If any R0 = NA, replace it with R0 = 0
  # (we only get R0=NA when V0<0, i.e. when mosquito recruitment is less than mortality)
  replace_na(list(R0 = 0)) %>%
  # Lower host density limit
  mutate(CHmin = compute.CHmin(.)) %>%
  # Upper host density limit
  mutate(CHmax = compute.CHmax(.)) %>%
  distinct()

AllOutputs_df <- AllOutputs_df %>%
  ## Select only variables used for visualization
  select(
    # Characteristics
    Model, system_ID, Temperature,
    # Host traits
    KH, sigmaH,
    # Vector abundance
    sigmaV, V0,
    # Basic reproduction number
    R0,
    # Critical minimum host density (only applicable when sigmaH is finite)
    CHmin,
    # Critical maximum host density
    CHmax
  ) %>%
  distinct()

# Save data frame---------------------------------------------------------------
write_csv(AllOutputs_df, "results/AllOutputs.csv")

# 3) Build thermal characteristics data frame ##################################

# Dataframes collecting thermal characteristics of R0:
# Topt  = temperature at which R0 is maximized
# CTmin = lowest temperature at which R0>1
# CTmax = highest temperature at which R0>1
# Tbr   = temperature range at half of max performance

# Estimate the thermal optimum for transmission
Topt_df <- AllOutputs_df %>%
  group_by(system_ID, sigmaH, KH) %>%
  filter(R0 == max(R0)) %>%
  # Get temperature at which R0 is maximized
  mutate(Topt = Temperature) %>%
  # Get the largest value of R0
  mutate(R0opt = R0) %>%
  # Keep track of whether R0 exceeds one
  mutate(threshold_bool = R0 > 1) %>%
  # remove duplicate rows
  distinct() %>%
  select(
    Model, system_ID, sigmaH, KH, Topt, R0opt, threshold_bool, CHmin, CHmax
  )

# Get the thermal range of parasite as a function of host traits
TempRange_df <- AllOutputs_df %>%
  # Restrict to rows where R0 exceeds one
  filter(R0 > 1) %>%
  group_by(system_ID, sigmaH, KH) %>%
  # Get lowest temperature at which R0 exceeds one
  mutate(CTmin = min(Temperature)) %>%
  # Get highest temperature at which R0 exceeds one
  mutate(CTmax = max(Temperature)) %>%
  select(Model, system_ID, sigmaH, KH, CTmin, CTmax) %>%
  # remove duplicate rows
  distinct()

ThermalCharacteristics_df <- left_join(Topt_df, TempRange_df,
  by = c("Model", "system_ID", "sigmaH", "KH")
) 

# Save data frame---------------------------------------------------------------
write_csv(ThermalCharacteristics_df, "results/ThermalCharacteristics.csv")
