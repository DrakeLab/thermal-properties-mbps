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
sigmaH_vec <- 10^seq(-0.25,2.25, length.out = sigmaH_vec_length - 6) %>% 
  c(1, 10, 20, 50, 100, Inf) %>% 
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
rm(data.in.params)

# 4) Combine and save data frames -----------------------------------------
Vec_dim <- dim(data.Vec)[1]
Host_dim <- dim(data.Host)[1]

# These are still >100 MB even after compression.
source("code/output-functions.R")


get.outputs <- function(in_df) {
  out_df <- in_df %>%
    # Name the model (just in case this is handier than referring to sigmaH)
    mutate(Model = ifelse(is.infinite(sigmaH), "Ross-Macdonald model", "Chitnis model")) %>%
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
    distinct() %>% 
    ## Select only variables used for visualization
    dplyr::select(
      # Characteristics
      Model, system_ID, Temperature,
      # Uncertainty (sample number)
      sample_num,
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
}

summarize.TempRange <- function(in_df) {
  out_df <- in_df %>%
    # Restrict to rows where R0 exceeds one
    filter(R0 > 1) %>%
    group_by(system_ID, sample_num, sigmaH, KH) %>%
    # Get lowest temperature at which R0 exceeds one
    mutate(CTmin = min(Temperature)) %>%
    # Get highest temperature at which R0 exceeds one
    mutate(CTmax = max(Temperature)) %>%
    dplyr::select(Model, system_ID, sample_num, sigmaH, KH, CTmin, CTmax) %>%
    # remove duplicate rows
    distinct() %>%
    arrange(system_ID, sample_num, sigmaH)%>%
    group_by(system_ID, Model, sigmaH, KH) %>%
    summarise(
      lowHCI_val = quantile(CTmin, 0.055),
      highHCI_val = quantile(CTmin, 0.945),
      mean_val = mean(norm_R0),
      median_val = median(norm_R0),
      # mode_val = mlv(norm_R0, method = 'mfv'),
      .groups = "keep"
    ) %>%
    arrange(system_ID, sigmaH, KH, mean_val, median_val)
}

summarize.R0 <- function(in_df) {
  out_df <- in_df %>%
    dplyr::select(system_ID, sample_num, Temperature, Model, sigmaH, KH, R0) %>%
    ungroup() %>% # !!!
    # Get mean and CIs of R0 across samples
    group_by(system_ID, Temperature, Model, sigmaH, KH) %>%
    summarise(
      lowHCI_val = quantile(R0, 0.055),
      highHCI_val = quantile(R0, 0.945),
      mean_val = mean(R0),
      median_val = median(R0),
      # mode_val = mlv(norm_R0, method = 'mfv'),
      .groups = "keep"
    ) %>%
    # Normalize R0 means across temperature
    group_by(system_ID, Model, sigmaH, KH) %>%
    mutate(norm_mean_val = mean_val / max(mean_val, na.rm = TRUE),
           norm_median_val = median_val / max(median_val, na.rm = TRUE),
           norm_lowHCI_val = lowHCI_val / max(lowHCI_val, na.rm = TRUE),
           norm_highHCI_val = highHCI_val / max(highHCI_val, na.rm = TRUE)) %>%
    # in the case that ALL lower HCI values are zero (resulting in NaN from above), replace with
    mutate(norm_lowHCI_val = ifelse(is.nan(norm_lowHCI_val), 0, norm_lowHCI_val)) %>%
    arrange(system_ID, sigmaH, KH, Temperature, mean_val, median_val) %>% # , mode_val) %>%
    distinct()
}

# 
# 
summarize.Topt <- function(in_df) {
  out_df <- in_df %>%
    dplyr::select(system_ID, sample_num, Model, sigmaH, KH, Topt) %>%
    ungroup() %>% # !!!
    # Get mean and CIs of variable across samples
    group_by(system_ID, Model, sigmaH, KH) %>%
    summarise(
      lowHCI_val = quantile(Topt, 0.055),
      highHCI_val = quantile(Topt, 0.945),
      mean_val = mean(Topt),
      median_val = median(Topt),
      # mode_val = mlv(norm_R0, method = 'mfv'),
      .groups = "keep"
    ) 
}


# Get and save data needed for figures ------------------------------------


# # Mean and quantiles of R0 TPCs (across KH and sigmaH) (use for Figure 4)
# data.Host %>%
#   # To set num. of curves, change "length.out" to be the number of curves you want
#   filter(KH %in% unique(KH)[seq(1, length(unique(KH)), length.out = 21)]) %>%
#   filter(sigmaH %in% c(100, Inf)) %>%
#   expand_grid(data.Vec) %>%
#   get.outputs(.) %>%
#   summarize.R0(.) %>%
#   write_rds("results/R0_TPC_data.rds")


# # Mean and quantiles of Topt (across KH and sigmaH) (use for Figures S2 and S3)
# data.Host %>% 
#   # To set num. of curves, change "length.out" to be the number of curves you want
#   filter(KH %in% unique(KH)[seq(1, length(unique(KH)), length.out = 51)]) %>%
#   filter(sigmaH %in% c(1, 10, 100, Inf)) %>%
#   expand_grid(data.Vec %>% 
#                 filter(between(Temperature, 20, 32))) %>% # known range of Topt
#   # filter(Temperature %in% unique(Temperature)[seq(1,
#   #                                                 length(unique(Temperature)),
#   #                                                 length.out = 51)])) %>%  # 51 works
#   # Name the model (just in case this is handier than referring to sigmaH)
#   mutate(Model = ifelse(is.infinite(sigmaH), "Ross-Macdonald model", "Chitnis model")) %>%
#   # Vector abundance
#   mutate(V0 = compute.V0(.)) %>%
#   # Basic reproduction number
#   mutate(R0 = sqrt(compute.RH(.) * compute.RV(.))) %>%
#   group_by(system_ID, sample_num, sigmaH, KH) %>%
#   filter(R0 == max(R0)) %>%
#   # Get temperature at which R0 is maximized
#   rename(Topt = Temperature) %>%
#   # If any R0 = NA, replace it with R0 = 0
#   # (we only get R0=NA when V0<0, i.e. when mosquito recruitment is less than mortality)
#   replace_na(list(R0 = 0)) %>%
#   ## Select only variables used for visualization
#   dplyr::select(Model, system_ID, sample_num, KH, sigmaH, Topt) %>%
#   distinct() %>% 
#   summarize.Topt(.) %>% 
#   write_rds("results/Topt_quantile_data.rds")

# Mean of Topt across KH and sigmaH
 

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
