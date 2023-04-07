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
library(foreach)

source("code/output-functions.R") # needed for compute.variable functions

# Set up parallel
# my.cluster <- parallel::makeCluster(
#   parallel::detectCores() - 1, 
#   type = "PSOCK"
# )

if (!exists("cluster")) {
  cluster <- new_cluster(parallel::detectCores() - 1)
  cluster_library(cluster, c("dplyr", "tidyr"))
}

# 1) Define accessory functions -------------------------------------------
slice<-function(x,n) {
  N<-length(x);
  lapply(seq(1,N,n),function(i) x[i:min(i+n-1,N)])
}

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
) %>% as.data.frame() %>% 
  mutate(Model = ifelse(is.infinite(sigmaH), "Ross-Macdonald model", "Chitnis model"))

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
  mutate(KL = larval_mosquito_carrying_capacity) %>% 
  mutate(V0 = ifelse( sigmaV_f * deltaL < (1 / lf),
                      0,
                      KL * rhoL * lf * (1 - 1 / (lf * sigmaV_f * deltaL))))

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
    mutate(RV = ifelse(is.infinite(sigmaH),
                       sigmaV * betaH / (1 / (lf + eps)), # Ross-Macdonald
                       sigmaH * sigmaV * betaH * KH / ((1 / (lf + eps)) * (sigmaH * KH + sigmaV * V0)))) %>% 
    mutate(bV = ifelse(is.infinite(sigmaH),
                       sigmaV, # Ross-Macdonald model
                       sigmaV * sigmaH * KH / (sigmaH * KH + sigmaV * V0 + eps))) %>% 
    mutate(RH = ifelse(V0 == 0, 
                       0, 
                       bV * betaV * V0 * exp(-1 / (lf * etaV)) / (KH * (gammaH + muH) + eps))) %>% 
    # Basic reproduction number
    mutate(R0 = sqrt(RV*RH))%>% 
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
                  filter(system_ID == system_name) #%>%
                # filter(between(Temperature, 20, 29)) # w/ host trait resolution at 40, HCImin of 11, HCImax of 39
    ) %>% # known range of Topt
    data.table::data.table() %>%
    # Name the model (just in case this is handier than referring to sigmaH)
    mutate(RV = ifelse(is.infinite(sigmaH),
                       sigmaV * betaH / (1 / (lf + eps)), # Ross-Macdonald
                       sigmaH * sigmaV * betaH * KH / ((1 / (lf + eps)) * (sigmaH * KH + sigmaV * V0)))) %>% 
    mutate(bV = ifelse(is.infinite(sigmaH),
                       sigmaV, # Ross-Macdonald model
                       sigmaV * sigmaH * KH / (sigmaH * KH + sigmaV * V0 + eps))) %>% 
    mutate(RH = ifelse(V0 == 0, 
                       0, 
                       bV * betaV * V0 * exp(-1 / (lf * etaV)) / (KH * (gammaH + muH) + eps))) %>% 
    # Basic reproduction number
    mutate(R0 = sqrt(RV*RH)) %>% 
    # Filter to maximum value of R0
    group_by(system_ID, sample_num, sigmaH, KH) %>%
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
      median = median(value)#,
      # .groups = "keep"
    ) %>%
    collect() %>%
    arrange(system_ID, sigmaH, KH) %>% 
    distinct()
}

CT_heat_func <- function(in_df, system_name) {
  out_df <- in_df %>%
    expand_grid(data.Vec %>%
                  # need to reduce sample size for memory purposes
                  # filter(between(Temperature, 13, 34)) %>% # known range of CT
                  filter(system_ID == system_name)) %>%
    mutate(RV = ifelse(is.infinite(sigmaH),
                       sigmaV * betaH / (1 / (lf + eps)), # Ross-Macdonald
                       sigmaH * sigmaV * betaH * KH / ((1 / (lf + eps)) * (sigmaH * KH + sigmaV * V0)))) %>% 
    mutate(bV = ifelse(is.infinite(sigmaH),
                       sigmaV, # Ross-Macdonald model
                       sigmaV * sigmaH * KH / (sigmaH * KH + sigmaV * V0 + eps))) %>% 
    mutate(RH = ifelse(V0 == 0, 
                       0, 
                       bV * betaV * V0 * exp(-1 / (lf * etaV)) / (KH * (gammaH + muH) + eps))) %>% 
    # Basic reproduction number
    mutate(R0 = sqrt(RV*RH)) %>%
    # Filter to maximum value of R0
    group_by(system_ID, sample_num, sigmaH, KH) %>%
    filter(R0 > 1) %>%
    # Get lowest temperature at which R0 exceeds one
    mutate(CTmin = min(Temperature)) %>%
    # Get highest temperature at which R0 exceeds one
    mutate(CTmax = max(Temperature)) %>%
    # Get width of critical thermal interval
    mutate(CTwidth = CTmax - CTmin)%>% 
    ungroup() %>%
    dplyr::select(system_ID, Model, sigmaH, KH, CTmin, CTmax, CTwidth) %>%
    pivot_longer(cols = c(CTwidth, CTmin, CTmax), names_to = "variable", values_to = "value")
  
  # If R0 < 1 across all temperatures, report the following:
  #  critical thermal minimum = Inf
  #  critical thermal maximum = -Inf
  #  critical thermal width = 0
  if (dim(out_df)[1] == 0) {
    out_df <- expand_grid(system_ID = system_name, sigmaH = in_df$sigmaH, KH = in_df$KH,
                          variable = c("CTmax", "CTmin", "CTwidth")) %>% 
      distinct() %>% 
      mutate(Model = ifelse(is.infinite(sigmaH), "Ross-Macdonald model", "Chitnis model")) %>%
      mutate(mean = case_when(
        variable == "CTmax" ~ -Inf,
        variable == "CTmin" ~ Inf,
        variable == "CTwidth" ~ 0
      )) %>% 
      mutate(median = mean) %>% 
      mutate(highHCI = mean) %>% 
      mutate(lowHCI = mean)
  } else {
    out_df <- out_df %>% 
      group_by(system_ID, Model, sigmaH, KH, variable) %>%
      # partition(cluster) %>%
      summarise(
        lowHCI = quantile(value, 0.055),
        highHCI = quantile(value, 0.945),
        mean = mean(value),
        median = median(value),
        .groups = "keep"
      ) %>%
      # collect() %>%
      arrange(system_ID, sigmaH, KH) %>%
      distinct()
  }
  return(out_df)
}

# Initialize a data frame to save our analyses
init.df <- tibble(system_ID = c(), Temperature = c(), Model = c(),
                  sigmaH = c(), KH = c(), variable = c(), 
                  lowHCI = c(), highHCI = c(), mean = c(), median = c())

# Remove the parameter data frame to free up memory
rm(data.in.params)

# R0 TPC data -------------------------------------------------------------

data.Host.R0_TPC <- data.Host %>%
  filter(KH %in% unique(KH)[seq(1, length(unique(KH)), length.out = 21)]) %>%
  filter(sigmaH %in% c(100, Inf))

R0_TPC.df <- init.df
for (system_name in unique(data.Vec$system_ID)) {
  print(paste0("R0 TPCs: ", system_name))
  system.time(
    R0_TPC.df <- data.Host.R0_TPC %>%
      R0_TPC_func(., system_name) %>%
      rbind(R0_TPC.df)
  )
  gc()
}

write_rds(R0_TPC.df, "results/R0_TPC_data.rds", compress = "gz")


# Topt vs. sigmaH data ----------------------------------------------------

data.Topt <- data.Host #%>%
# filter(KH %in% c(1, 10, 100, 1000, 10000)) # %>%
# filter(sigmaH %in% unique(sigmaH)[seq(1, length(unique(sigmaH)), length.out = 50)])

KH_slices <- slice(unique(data.Topt$KH), 11)

Topt.df <- init.df

gc()
for (system_name in unique(data.Vec$system_ID)) {
  print(paste0("CT vals: ", system_name))
  slice_num <- 1
  for(index_KH in KH_slices) {
    print(paste0("Slice number ", slice_num, " out of ", length(KH_slices)))
    for (index_sigmaH in unique(data.Topt$sigmaH)) {
      print(paste0(which(unique(data.Topt$sigmaH)== index_sigmaH)/length(unique(data.Topt$sigmaH))*100, "% complete"))
      # parallel compute CT values over host trait values
      # for (index_KH in unique(data.Host.CT$KH)) {
      system.time(Topt.df <- data.Topt %>% 
                    filter(sigmaH == index_sigmaH,
                           KH %in% index_KH) %>% 
                    Topt_heat_func(., system_name) %>% 
                    rbind(Topt.df))
      # }
    }
    slice_num <- slice_num +1
    gc()
  }
}

write_rds(Topt.df, "results/Topt_vals.rds", compress = "gz")


# Topt vs. KH -------------------------------------------------------------

data.Topt.KH <- data.Host %>%
  # filter(KH %in% unique(KH)[seq(1, length(unique(KH)), length.out = 10)]) %>%
  filter(sigmaH %in% c(1, 10, 100, Inf))

Topt.df <- init.df

gc()
for (system_name in unique(data.Vec$system_ID)) {
  print(paste0("Topt vs. KH: ", system_name))
  system.time(
    Topt.df <- data.Topt.KH %>%
      Topt_heat_func(., system_name) %>%
      rbind(Topt.df)
  )
  gc()
}

write_rds(Topt.df, "results/Topt_KH.rds", compress = "gz")


# CTmin, max, and width vs. sigma and KH ----------------------------------

data.Host.CT <- data.Host #%>%
# filter(KH %in% unique(KH)[seq(1, length(unique(KH)), length.out = 4)]) %>%
# filter(sigmaH %in% unique(sigmaH)[seq(1, length(unique(sigmaH)), length.out = 4)])

CT.df <- init.df

gc()
system.time(
  for (system_name in unique(data.Vec$system_ID)) {
    print(paste0("CT vals: ", system_name))
    for (index_sigmaH in unique(data.Host.CT$sigmaH)) {
      print(paste0(which(unique(data.Host.CT$sigmaH)== index_sigmaH)/length(unique(data.Host.CT$sigmaH))*100, "% complete"))
      # parallel compute CT values over host trait values
      # for (index_KH in unique(data.Host.CT$KH)) {
      CT.df <- data.Host.CT %>% 
        filter(sigmaH == index_sigmaH) %>% 
        CT_heat_func(., system_name) %>% 
        rbind(CT.df)
      # }
    }
    gc()
  }
)

# smallest CTmin ~13.85
# largest CTmax ~34.45

write_rds(CT.df, "results/CT_vals.rds", compress = "gz")



# *) Diagnostics & visualizations -----------------------------------------