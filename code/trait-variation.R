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
library(multidplyr)

source("code/output-functions.R") # needed for compute.variable functions

# Set up parallel
cluster <- new_cluster(parallel::detectCores() - 1)
cluster_library(cluster, c("dplyr", "tidyr"))

# 1) Define accessory functions -------------------------------------------
# 
# get.outputs <- function(in_df) {
#   out_df <- in_df %>%
#     # Name the model (just in case this is handier than referring to sigmaH)
#     mutate(Model = ifelse(is.infinite(sigmaH), "Ross-Macdonald model", "Chitnis model")) %>%
#     # Vector abundance
#     mutate(V0 = compute.V0(.)) %>%
#     # Basic reproduction number
#     mutate(R0 = sqrt(compute.RH(.) * compute.RV(.))) %>%
#     # If any R0 = NA, replace it with R0 = 0
#     # (we only get R0=NA when V0<0, i.e. when mosquito recruitment is less than mortality)
#     replace_na(list(R0 = 0)) %>%
#     # Lower host density limit
#     mutate(CHmin = compute.CHmin(.)) %>%
#     # Upper host density limit
#     mutate(CHmax = compute.CHmax(.)) %>%
#     distinct() %>% 
#     ## Select only variables used for visualization
#     dplyr::select(
#       # Characteristics
#       Model, system_ID, Temperature,
#       # Uncertainty (sample number)
#       sample_num,
#       # Host traits
#       KH, sigmaH,
#       # Vector abundance
#       sigmaV, V0,
#       # Basic reproduction number
#       R0,
#       # Critical minimum host density (only applicable when sigmaH is finite)
#       CHmin,
#       # Critical maximum host density
#       CHmax
#     ) %>%
#     distinct()
# }
# 
# summarize.R0 <- function(in_df) {
#   out_df <- in_df %>%
#     dplyr::select(system_ID, sample_num, Temperature, Model, sigmaH, KH, R0) %>%
#     ungroup() %>% # !!!
#     # Get mean and CIs of R0 across samples
#     group_by(system_ID, Temperature, Model, sigmaH, KH) %>%
#     summarise(
#       lowHCI_val = quantile(R0, 0.055),
#       highHCI_val = quantile(R0, 0.945),
#       mean_val = mean(R0),
#       median_val = median(R0),
#       # mode_val = mlv(norm_R0, method = 'mfv'),
#       .groups = "keep"
#     ) %>%
#     # Normalize R0 means across temperature
#     group_by(system_ID, Model, sigmaH, KH) %>%
#     mutate(norm_mean_val = mean_val / max(mean_val, na.rm = TRUE),
#            norm_median_val = median_val / max(median_val, na.rm = TRUE),
#            norm_lowHCI_val = lowHCI_val / max(lowHCI_val, na.rm = TRUE),
#            norm_highHCI_val = highHCI_val / max(highHCI_val, na.rm = TRUE)) %>%
#     # in the case that ALL lower HCI values are zero (resulting in NaN from above), replace with
#     mutate(norm_lowHCI_val = ifelse(is.nan(norm_lowHCI_val), 0, norm_lowHCI_val)) %>%
#     arrange(system_ID, sigmaH, KH, Temperature, mean_val, median_val) %>% # , mode_val) %>%
#     distinct()
# }

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

# thin_size <- 300
# num_samples <- length(unique(data.in.params$sample_num))
# sample_inds <- sample(unique(data.in.params$sample_num), thin_size, replace = FALSE)

data.Vec <- data.in.params %>% 
  # filter(sample_num %in% sample_inds) %>%
  mutate(KL = larval_mosquito_carrying_capacity)

# rm(data.in.params)

# 4) Combine and save data frames -----------------------------------------
Vec_dim <- dim(data.Vec)[1]
Host_dim <- dim(data.Host)[1]


# Get and save data needed for figures ------------------------------------

# Each of these proceeds by by producing a data frame of outputs, restricting
# to the appropriate case for the output of interest, then summarizing the
# outputs by calculating the mean, median, and the endpoints of the 89% HPIs


get.summary <- function(in_df, summary_vars, group_vars) {
  out_df <- in_df %>%
    pivot_longer(cols = {{summary_vars}}, names_to = "variable", values_to = "value") %>% 
    group_by(c(!!!syms(group_vars), variable)) %>%
    partition(cluster) %>%
    summarise(
      lowHCI = quantile(value, 0.055),
      highHCI = quantile(value, 0.945),
      mean = mean(value),
      median = median(value)
    ) %>%
    arrange(system_ID, sigmaH, KH) %>% 
    collect() %>%
    distinct()
}
eps <- .Machine$double.eps

R0_TPC_func <- function(in_df, system_name) {
  out_df <- in_df %>%
    expand_grid(data.Vec %>% filter(system_ID == system_name)) %>%
    mutate(Model = ifelse(is.infinite(sigmaH), "Ross-Macdonald model", "Chitnis model")) %>%
    mutate(V0 = ifelse( sigmaV_f * deltaL < (1 / lf),
                        0,
                        KL * rhoL * lf * (1 - 1 / (lf * sigmaV_f * deltaL)))) %>%
    # Basic reproduction number
    mutate(R0 = ifelse(V0 <= 0, 0,
                       ifelse(is.infinite(sigmaH),
                              sigmaV * sigmaV * betaH / (1 / (lf)), # Ross-Macdonald model
                              (sigmaV * sigmaH * KH / (sigmaH * KH + sigmaV * V0)) * betaV * V0 * exp(-1 / (lf * etaV)) / (KH * (gammaH + muH)) * sigmaH * sigmaV * betaH * KH / ((1 / (lf)) * (sigmaH * KH + sigmaV * V0))))) %>% 
    dplyr::select(system_ID, sample_num, Temperature, Model, sigmaH, KH, V0, R0) %>%
    # Normalize R0 across temperature
    group_by(system_ID, Model, sigmaH, KH, sample_num) %>%
    mutate(norm_R0 = R0 / max(R0, eps)) %>%
    ungroup() %>%
    dplyr::select(system_ID, Temperature, Model, sigmaH, KH, norm_R0) %>%
    pivot_longer(cols = norm_R0, names_to = "variable", values_to = "value") %>% 
    group_by(system_ID, Temperature, Model, sigmaH, KH, variable) %>%
    partition(cluster) %>%
    summarise(
      lowHCI = quantile(value, 0.055),
      highHCI = quantile(value, 0.945),
      mean = mean(value),
      median = median(value)
    ) %>%
    collect() %>%
    arrange(system_ID, sigmaH, KH) %>% 
    distinct()
}

Topt_heat_func <- function(in_df, system_name) {
  out_df <- in_df %>%
    expand_grid(data.Vec %>%
                  filter(system_ID == system_name) %>%
                  filter(between(Temperature, 20, 29))) %>% # known range of Topt
    data.table::data.table() %>%
    # Name the model (just in case this is handier than referring to sigmaH)
    mutate(Model = ifelse(is.infinite(sigmaH), "Ross-Macdonald model", "Chitnis model")) %>%
    group_by(sigmaH) %>%
    # partition(cluster) %>%
    mutate(V0 = ifelse( sigmaV_f * deltaL < (1 / lf),
                        0,
                        KL * rhoL * lf * (1 - 1 / (lf * sigmaV_f * deltaL)))) %>%
    # Basic reproduction number
    mutate(R0 = ifelse(V0 <= 0, 0,
                       ifelse(is.infinite(sigmaH),
                              sigmaV * sigmaV * betaH / (1 / (lf)), # Ross-Macdonald model
                              sigmaV * sigmaH * KH / (sigmaH * KH + sigmaV * V0)) * betaV * V0 * exp(-1 / (lf * etaV)) / (KH * (gammaH + muH)) * sigmaH * sigmaV * betaH * KH / ((1 / (lf)) * (sigmaH * KH + sigmaV * V0)))) %>%
    group_by(sample_num, sigmaH, KH) %>%
    # Filter to maximum value of R0
    filter(R0 == max(R0)) %>%
    # Get temperature at which R0 is maximized
    rename(Topt = Temperature) %>%
    dplyr::select(system_ID, sample_num, Model, sigmaH, KH, Topt) %>%
    # If any R0 = NA, replace it with R0 = 0
    # (we only get R0=NA when V0<0, i.e. when mosquito recruitment is less than mortality)
    replace_na(list(R0 = 0)) %>%
    ungroup() %>%
    pivot_longer(cols = Topt, names_to = "variable", values_to = "value") %>% 
    group_by(system_ID, Model, sigmaH, KH, variable) %>%
    partition(cluster) %>%
    summarise(
      lowHCI = quantile(value, 0.055),
      highHCI = quantile(value, 0.945),
      mean = mean(value),
      median = median(value)
    ) %>%
    collect() %>%
    arrange(system_ID, sigmaH, KH) %>% 
    distinct()
}

CT_heat_func <- function(in_df, system_name) {
  out_df <- in_df %>%
    expand_grid(data.Vec %>%
                  filter(between(Temperature, 13, 34)) %>% # known range of CT
                  filter(system_ID == system_name)) %>%
    # Name the model (just in case this is handier than referring to sigmaH)
    mutate(Model = ifelse(is.infinite(sigmaH), "Ross-Macdonald model", "Chitnis model")) %>%
    mutate(V0 = ifelse( sigmaV_f * deltaL < (1 / lf),
                        0,
                        KL * rhoL * lf * (1 - 1 / (lf * sigmaV_f * deltaL)))) %>%
    # Basic reproduction number
    mutate(R0 = ifelse(V0 <= 0, 0,
                       ifelse(is.infinite(sigmaH),
                              sigmaV * sigmaV * betaH / (1 / (lf)), # Ross-Macdonald model
                              sigmaV * sigmaH * KH / (sigmaH * KH + sigmaV * V0)) * betaV * V0 * exp(-1 / (lf * etaV)) / (KH * (gammaH + muH)) * sigmaH * sigmaV * betaH * KH / ((1 / (lf)) * (sigmaH * KH + sigmaV * V0)))) %>% 
    # Filter to maximum value of R0
    group_by(system_ID, sample_num, sigmaH, KH) %>%
    filter(R0 > 1) %>%
    # Get lowest temperature at which R0 exceeds one
    mutate(CTmin = min(Temperature)) %>%
    # Get highest temperature at which R0 exceeds one
    mutate(CTmax = max(Temperature)) %>%
    # Get width of critical thermal interval
    mutate(CTwidth = CTmax - CTmin) %>% 
    ungroup() %>%
    dplyr::select(system_ID, Model, sigmaH, KH, CTmin, CTmax, CTwidth) %>%
    pivot_longer(cols = c(CTwidth, CTmin, CTmax), names_to = "variable", values_to = "value") %>%
    group_by(system_ID, Model, sigmaH, KH, variable) %>%
    partition(cluster) %>%
    summarise(
      lowHCI = quantile(value, 0.055),
      highHCI = quantile(value, 0.945),
      mean = mean(value),
      median = median(value)
    ) %>%
    arrange(system_ID, sigmaH, KH) %>%
    collect() %>%
    distinct()
}


# # System 1: Aedes aegypti / DENV ------------------------------------------
# # Code = AeDENV
# system_name = "Aedes aegypti / DENV"
# 
# # Topt TPCs
# gc()
# system.time(
#   data.Host %>%
#     filter(sigmaH %in% c(100, Inf))  %>%
#     R0_TPC_func(., system_name) %>%
#     write_rds("results/AeDENV/R0_TPC_data.rds", compress = "gz")
# )
# 
# # Topt vs. sigmaH and KH
# gc()
# system.time(
#   data.Host %>%
#     filter(KH %in% unique(KH)[seq(1, length(unique(KH)), length.out = 61)]) %>%
#     filter(sigmaH %in% unique(sigmaH)[seq(1, length(unique(sigmaH)), length.out = 61)]) %>%
#     Topt_heat_func(., system_name) %>%
#     write_rds("results/AeDENV/Topt_vals.rds", compress = "gz")
# )
# gc()
# 
# # CTmin, CTmax, CTwidth vs. sigmaH and KH
# gc()
# system.time(
#   data.Host %>%
#     filter(KH %in% unique(KH)[seq(1, length(unique(KH)), length.out = 61)]) %>%
#     filter(sigmaH %in% unique(sigmaH)[seq(1, length(unique(sigmaH)), length.out = 61)]) %>%
#     CT_heat_func(., system_name) %>%
#     write_rds("results/AeDENV/CT_vals.rds", compress = "gz")
# )
# gc()
# 
# # System 2: Aedes aegypti / ZIKV ------------------------------------------
# # Code = AeZIKV
# system_name = "Aedes aegypti / ZIKV"
# 
# # Topt TPCs
# gc()
# system.time(
#   data.Host %>%
#     filter(sigmaH %in% c(100, Inf))  %>%
#     R0_TPC_func(., system_name) %>%
#     write_rds("results/AeZIKV/R0_TPC_data.rds", compress = "gz")
# )
# 
# # Topt vs. sigmaH and KH
# gc()
# system.time(
#   data.Host %>%
#     filter(KH %in% unique(KH)[seq(1, length(unique(KH)), length.out = 61)]) %>%
#     filter(sigmaH %in% unique(sigmaH)[seq(1, length(unique(sigmaH)), length.out = 61)]) %>%
#     Topt_heat_func(., system_name) %>%
#     write_rds("results/AeZIKV/Topt_vals.rds", compress = "gz")
# )
# gc()
# 
# # CTmin, CTmax, CTwidth vs. sigmaH and KH
# gc()
# system.time(
#   data.Host %>%
#     filter(KH %in% unique(KH)[seq(1, length(unique(KH)), length.out = 61)]) %>%
#     filter(sigmaH %in% unique(sigmaH)[seq(1, length(unique(sigmaH)), length.out = 61)]) %>%
#     CT_heat_func(., system_name) %>%
#     write_rds("results/AeZIKV/CT_vals.rds", compress = "gz")
# )
# gc()
# 
# # System 3: Aedes albopictus / DENV ---------------------------------------
# # Code = AlDENV
# system_name = "Aedes albopictus / DENV"
# 
# # Topt TPCs
# gc()
# system.time(
#   data.Host %>%
#     filter(sigmaH %in% c(100, Inf))  %>%
#     R0_TPC_func(., system_name) %>%
#     write_rds("results/AlDENV/R0_TPC_data.rds", compress = "gz")
# )
# 
# # Topt vs. sigmaH and KH
# gc()
# system.time(
#   data.Host %>%
#     filter(KH %in% unique(KH)[seq(1, length(unique(KH)), length.out = 61)]) %>%
#     filter(sigmaH %in% unique(sigmaH)[seq(1, length(unique(sigmaH)), length.out = 61)]) %>%
#     Topt_heat_func(., system_name) %>%
#     write_rds("results/AlDENV/Topt_vals.rds", compress = "gz")
# )
# gc()
# 
# # CTmin, CTmax, CTwidth vs. sigmaH and KH
# gc()
# system.time(
#   data.Host %>%
#     filter(KH %in% unique(KH)[seq(1, length(unique(KH)), length.out = 61)]) %>%
#     filter(sigmaH %in% unique(sigmaH)[seq(1, length(unique(sigmaH)), length.out = 61)]) %>%
#     CT_heat_func(., system_name) %>%
#     write_rds("results/AlDENV/CT_vals.rds", compress = "gz")
# )
# gc()

# System 4: Anopheles gambiae / Plasmodium relictum -----------------------
# Code = AgPlfa
system_name = "Anopheles gambiae / Plasmodium relictum"

# Topt TPCs
gc()
system.time(
  data.Host %>%
    filter(sigmaH %in% c(100, Inf))  %>%
    R0_TPC_func(., system_name) %>%
    write_rds("results/AgPlfa/R0_TPC_data.rds", compress = "gz")
)

# Topt vs. sigmaH and KH
gc()
system.time(
  data.Host %>%
    filter(KH %in% unique(KH)[seq(1, length(unique(KH)), length.out = 61)]) %>%
    filter(sigmaH %in% unique(sigmaH)[seq(1, length(unique(sigmaH)), length.out = 61)]) %>%
    Topt_heat_func(., system_name) %>%
    write_rds("results/AgPlfa/Topt_vals.rds", compress = "gz")
)
gc()

# CTmin, CTmax, CTwidth vs. sigmaH and KH
gc()
system.time(
  data.Host %>%
    filter(KH %in% unique(KH)[seq(1, length(unique(KH)), length.out = 61)]) %>%
    filter(sigmaH %in% unique(sigmaH)[seq(1, length(unique(sigmaH)), length.out = 61)]) %>%
    CT_heat_func(., system_name) %>%
    write_rds("results/AgPlfa/CT_vals.rds", compress = "gz")
)
gc()


# System 5: Culex quinquefasciatus / WNV ----------------------------------
# Code = CxWNV
system_name = "Culex quinquefasciatus / WNV"

# Topt TPCs
gc()
system.time(
  data.Host %>%
    filter(sigmaH %in% c(100, Inf))  %>%
    R0_TPC_func(., system_name) %>%
    write_rds("results/CxWNV/R0_TPC_data.rds", compress = "gz")
)

# Topt vs. sigmaH and KH
gc()
system.time(
  data.Host %>%
    filter(KH %in% unique(KH)[seq(1, length(unique(KH)), length.out = 61)]) %>%
    filter(sigmaH %in% unique(sigmaH)[seq(1, length(unique(sigmaH)), length.out = 61)]) %>%
    Topt_heat_func(., system_name) %>%
    write_rds("results/CxWNV/Topt_vals.rds", compress = "gz")
)
gc()

# CTmin, CTmax, CTwidth vs. sigmaH and KH
gc()
system.time(
  data.Host %>%
    filter(KH %in% unique(KH)[seq(1, length(unique(KH)), length.out = 61)]) %>%
    filter(sigmaH %in% unique(sigmaH)[seq(1, length(unique(sigmaH)), length.out = 61)]) %>%
    CT_heat_func(., system_name) %>%
    write_rds("results/CxWNV/CT_vals.rds", compress = "gz")
)
gc()



# *) Diagnostics & visualizations -----------------------------------------