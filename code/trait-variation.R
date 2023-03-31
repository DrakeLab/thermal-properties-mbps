## Title: Transform TPC traits into model parameters #######################
##
## Project: Global zoonoses - spillover of mosquito-borne pathogens
##
## Purpose: Transform trait values, which we have fit TPCs to from data, in 
##          order to use in model analyses
##
## Contents: 0) Set-up, load in necessary packages and data-sets
##           1) Define accessory functions
##           2) Set *host* parameters
##           3) Set non-TPC *mosquito* parameters
##           4) Create host trait data frame
##           5) Combine and save data frames
##           *) Diagnostics and visualizations
##
##
## Inputs:  data - data/clean/trait_transforms.rds
##
## Outputs: data - data/clean/parameter_TPCs.rds
##
## Written and maintained by: Kyle Dahlin, kydahlin@gmail.com
## Initialized March 2023
# _______________________________________________________________________________

# 0) Set-up, load in necessary packages and data-sets ---------------------

# Load Libraries
library(tidyverse)
library(reshape2)
library(fst)

# 1) Define accessory functions -------------------------------------------



# 2) Set *host* parameters --------------------------------

## Host life history & behavioral traits ----
# Host recruitment rate:
# upper estimate from range for Primate traits: (0.001150685, 0.009624300)
lambdaH_baseline <- .005

# Host mortality rate:
# estimate from range for Primate traits: lifespan (8.6, 60) years
# muH_vec <- 1 / (365 * c(1, 25))
muH_baseline <- 1 / (365 * 20)

# Host maximum biting tolerance (mosquitoes bites per day)
sigmaH_vec <- 10^seq(-0.25,3.25, length.out = sigmaH_vec_length) %>% 
  c(1, 20, 50, Inf) %>% 
  unique() %>% sort()
sigmaH_baseline <- 100

# Host carrying capacity
KH_vec <- 10^seq(-2, 5, length.out = KH_vec_length) %>% 
  c(10^seq(-2,5)) %>% 
  unique() %>% sort()
      
# KH_vec2 <-  KH_vec %>%
#   c(10^(seq(5,5+log10(2), by = diff(log10(KH_vec), lag = 1)[1]))) %>%
#   unique()

## Host-related pathogen parameters
# Probability of becoming infected after bite from infectious mosquito
# operates as a scaling parameter
# betaH_vec <- c(.25, .75)
betaH_baseline <- 1

# Host recovery rate
# plausible estimate for infectious period = (2-14 days)
gammaH_vec <- 1 / c(5, 14)
gammaH_baseline <- 1 / 5

# 3) Create host trait data frame -----------------------------------------

data.Host <- expand_grid(
  # Life-history parameters
  lambdaH = lambdaH_baseline,
  muH = muH_baseline,
  KH = KH_vec,
  sigmaH = sigmaH_vec,
  # Infection-related parameters
  gammaH = gammaH_baseline,
  betaH = betaH_baseline
) %>% as.data.frame()

# 3) Set non-TPC *mosquito* parameters ------------------------------------

## Carrying capacity for larval mosquitoes
# NB: In the absence of good estimates for each species or temperature-dependence of this trait, we assume that this parameter is constant. It can be used as a  scaling parameter for overall mosquito abundance 
# (it could alternately be used to fix the maximum adult mosquito density across species)
larval_mosquito_carrying_capacity <- 300

data.Vec <- data.in.params %>% 
  mutate(KL = larval_mosquito_carrying_capacity)


# 4) Combine and save data frames -----------------------------------------
Vec_dim <- dim(data.Vec)[1]
Host_dim <- dim(data.Host)[1]

# These are still >100 MB even after compression.

# filter(data.Vec, system_ID == "Aedes aegypti / DENV") %>% 
#   expand_grid(data.Host) %>% 
#   write_rds("data/clean/AeDENV_data.rds", compress = "gz")
#   # write_fst("data/clean/AeDENV_data.fst", compress = 100)
# 
# filter(data.Vec, system_ID == "Aedes aegypti / ZIKV") %>% 
#   expand_grid(data.Host) %>% 
#   write_rds("data/clean/AeZIKV_data.rds", compress = "gz")
# 
# 
# filter(data.Vec, system_ID == "Aedes albopictus / DENV") %>% 
#   expand_grid(data.Host) %>% 
#   write_rds("data/clean/AlDENV_data.rds", compress = "gz")
# 
# filter(data.Vec, system_ID == "Culex quinquefasciatus / WNV") %>% 
#   expand_grid(data.Host) %>% 
#   write_rds("data/clean/CxWNV_data.rds", compress = "gz")
# 
# data.AnPlas <- filter(data.Vec, system_ID == "Anopheles gambiae / Plasmodium") %>% # "Anopheles gambiae / Plasmodium falciparum") %>% 
#   expand_grid(data.Host) %>% 
#   write_rds("data/clean/AnPlas_data.rds", compress = "gz")


# data.in.analysis <- expand_grid(data.Vec, data.Host)

# write_rds(data.in.analysis, "data/clean/full_traitset.rds", compress = "gz")


# *) Diagnostics & visualizations -----------------------------------------

plot_bool = FALSE
if (plot_bool) {
 
  
  # # Save figure
  # ggsave("figures/param_TPC_plot.svg",
  #        plot = TPC_plot,
  #        device = "svg",
  #        width = 16, height = 9, units = "in")
}
