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


# 1) Calculate outputs ----------------------------------------------------

## Combine all parameter combinations into a large table
AllOutputs_df <- data.in.analysis %>%
  # Lower the resolution of temperature to length Temp_vec_length
  filter(Temperature %in% res_reduce(Temperature, Temp_vec_length)) %>%
  # Name the model (just in case this is handier than referring to sigmaH)
  mutate(Model = ifelse(is.infinite(sigmaH), "Ross-Macdonald model", "Chitnis model")) %>%
  ## Compute outputs ##
  # Vector lifespan !!! do this in an earlier dataset instead
  mutate(lf = 1/muV) %>% 
  # Vector abundance
  mutate(V0 = compute.V0(.)) %>%
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
  dplyr::select(
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
write_rds(AllOutputs_df, "data/clean/AllOutputs.rds",
          compress = "gz")


# 3) Build thermal characteristics data frames ----------------------------

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
  dplyr::select(
    Model, system_ID, sigmaH, KH, Topt, R0opt, threshold_bool, CHmin, CHmax
  )

write_rds(Topt_df, "results/Topt_vals.rds",
          compress = "gz")

# Get the thermal range of parasite as a function of host traits
TempRange_df <- AllOutputs_df %>%
  # Restrict to rows where R0 exceeds one
  filter(R0 > 1) %>%
  group_by(system_ID, sigmaH, KH) %>%
  # Get lowest temperature at which R0 exceeds one
  mutate(CTmin = min(Temperature)) %>%
  # Get highest temperature at which R0 exceeds one
  mutate(CTmax = max(Temperature)) %>%
  dplyr::select(Model, system_ID, sigmaH, KH, CTmin, CTmax) %>%
  # remove duplicate rows
  distinct()

write_rds(TempRange_df, "results/TempRange_vals.rds",
          compress = "gz")

ThermalCharacteristics_df <- left_join(Topt_df, TempRange_df,
  by = c("Model", "system_ID", "sigmaH", "KH")
) 

# Save data frame---------------------------------------------------------------
write_rds(ThermalCharacteristics_df, "results/AllThermChar_vals.rds",
          compress = "gz")
