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
library(doParallel)
library(progress)

# 1) Define accessory functions -------------------------------------------

eps <- .Machine$double.eps

# Function: Slice data to optimize usage of parallel processing and memory
slice<-function(x,n) {
  N<-length(x);
  lapply(seq(1,N,n),function(i) x[i:min(i+n-1,N)])
}

# Function: evaluate R0 (as a function of temperature), its mean, median, and 
#           highest prob. intervals across samples
R0_TPC_func <- function(in_df, system_name) {
  out_df <- in_df %>%
    expand_grid(filter(data.Vec, system_ID == system_name)) %>%
    data.table::data.table() %>%
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
    dplyr::select(system_ID, sample_num, Temperature, Model, sigmaH, KH, R0) %>%
    # Normalize R0 across temperature
    group_by(system_ID, Model, sigmaH, KH) %>%
    ungroup() %>%
    pivot_longer(cols = c(R0), names_to = "variable", values_to = "value") %>%
    group_by(system_ID, Temperature, Model, sigmaH, KH, variable) %>%
    partition(cluster) %>%
    summarise(
      lowHCI = quantile(value, 0.055),
      highHCI = quantile(value, 0.945),
      mean = mean(value),
      median = median(value)
    ) %>%
    collect()
}

# Function: evaluate Topt and its mean, median, and highest prob. intervals across
#           samples
Topt_heat_func <- function(in_df, system_name) {
  out_df <- in_df %>%
    expand_grid(filter(data.Vec, system_ID == system_name)) %>% 
    data.table::data.table() %>%
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
    ungroup() %>%
    pivot_longer(cols = Topt, names_to = "variable", values_to = "value") %>%
    group_by(system_ID, Model, sigmaH, KH, variable) %>%
    # partition(cluster) %>%
    summarise(
      lowHCI = quantile(value, 0.055),
      highHCI = quantile(value, 0.945),
      mean = mean(value),
      median = median(value),
      .groups = "keep"
    ) #%>%
    # collect()
}

# Function: evaluate CTmin, CTmax, CTwidth and their means, medians, and highest
#           prob. intervals across samples
CT_heat_func <- function(in_df, system_name) {
  out_df <- in_df %>%
    expand_grid(filter(data.Vec, system_ID == system_name)) %>%
    data.table::data.table() %>%
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
    filter(R0 > 1) 
  
  # If R0 < 1 across all temperatures, report the following:
  #  CTmin = Inf, CTmax = -Inf, CTwidth = 0
  if (dim(out_df)[1] == 0) {
    out_df <- expand_grid(select(in_df, sigmaH, KH, Model),
                          system_ID = system_name,
                          variable = c("CTmax", "CTmin", "CTwidth")) %>%
      distinct() %>%
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
      # Get lowest temperature at which R0 exceeds one
      mutate(CTmin = min(Temperature)) %>%
      # Get highest temperature at which R0 exceeds one
      mutate(CTmax = max(Temperature)) %>%
      # Get width of critical thermal interval
      mutate(CTwidth = ifelse(is.finite(CTmax), 
                              CTmax - CTmin,
                              0))%>%
      ungroup() %>%
      dplyr::select(system_ID, Model, sigmaH, KH, CTmin, CTmax, CTwidth) %>%
      pivot_longer(cols = c(CTwidth, CTmin, CTmax), names_to = "variable", values_to = "value") %>% 
      group_by(system_ID, Model, sigmaH, KH, variable) %>%
      # partition(cluster) %>%
      summarise(
        lowHCI = quantile(value, 0.055),
        highHCI = quantile(value, 0.945),
        mean = mean(value),
        median = median(value),
        .groups = "keep"
      ) #%>%
    # collect()
  }
  return(out_df)
}
# 2) Set *host* parameters --------------------------------

## Host life history & behavioral traits
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

## Host-related pathogen parameters
# Probability of becoming infected after bite from infectious mosquito
# operates as a scaling parameter
# betaH_vec <- c(.25, .75)
betaH_baseline <- 1

# Host recovery rate
# plausible estimate for infectious period = (2-14 days)
gammaH_vec <- 1 / c(5, 14)
gammaH_baseline <- 1 / 5

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

# 4) Combine and save data frames -----------------------------------------

# Initialize a data frame to save our analyses
init.df <- tibble(system_ID = c(), Temperature = c(), Model = c(),
                  sigmaH = c(), KH = c(), variable = c(),
                  lowHCI = c(), highHCI = c(), mean = c(), median = c())

# Remove the parameter data frame to free up memory
rm(data.in.params)

### R0 TPCs ----

# Set up parallel
if (!exists("cluster")) {
  cluster_size <- parallel::detectCores() - 1
  cluster <- new_cluster(cluster_size)
  cluster_library(cluster, c("dplyr", "tidyr"))
}

# Set up host trait data frame (for future visualization)
data.R0 <- data.Host %>%
  filter(sigmaH %in% c(100, Inf)) %>% 
  filter(KH %in% unique(KH)[seq(1, length(unique(KH)), length.out = 21)])

# Slice host trait data
sigmaH_slices <- slice(unique(data.R0$sigmaH), 1)
KH_slices <- slice(unique(data.R0$KH), 1)

# Initialize R0 data frame
R0.df <- init.df

gc()

# Collect R0 TPC data across systems and host trait values
for (system_name in unique(data.Vec$system_ID)) {
  print(paste0("Topt vals: ", system_name))
  KHslice_num <- 1
  for(index_KH in KH_slices) {
    print(paste0("KH slice number ", KHslice_num, " out of ", length(KH_slices)))
    sigmaHslice_num <- 1
    for (index_sigmaH in sigmaH_slices) {
      print(paste0("sigmaH slice number ", sigmaHslice_num, " out of ", length(sigmaH_slices)))
      R0.df <- data.R0 %>%
        filter(sigmaH %in% index_sigmaH,
               KH %in% index_KH) %>%
        R0_TPC_func(., system_name) %>%
        rbind(R0.df)
      sigmaHslice_num <- sigmaHslice_num  + 1
    }
    KHslice_num <- KHslice_num +1
    gc()
  }
}

# Save data
proper_dim <- 2 * dim(data.R0)[1] * length(unique(data.Vec$system_ID)) * length(unique(data.Vec$Temperature))

if (exists("R0.df") & dim(R0.df)[1] == proper_dim)
{write_rds(R0.df, "results/R0_vals.rds", compress = "gz")
} else {
  warning("No file written. R0.df either empty or not complete.")
}


### Topt ----
# Set up new cluster
# Close old cluster connections
rm(cluster)
gc()
# Start new cluster for doParallel
cluster_size <- parallel::detectCores() - 2

my.cluster <- parallel::makeCluster(
  cluster_size, 
  type = "PSOCK"
)
# Register cluster for doParallel
doSNOW::registerDoSNOW(cl = my.cluster)

# Set up iteration grid
iter_grid <- expand_grid(system_ID = unique(data.Vec$system_ID),
                         tibble(KH = data.Host$KH,
                                sigmaH = data.Host$sigmaH))
# Set up progress bar
iterations <- dim(iter_grid)[1]

pb <- progress_bar$new(
  format = ":spin :system progress = :percent [:bar] :elapsed | eta: :eta",
  total = iterations,
  width = 120)                                                                                                         

progress <- function(n){
  pb$tick(tokens = list(system = iter_grid$system_ID[n]))
}
opts <- list(progress = progress)

Topt.df <- foreach(
  system_name = iter_grid$system_ID,
  index_KH = iter_grid$KH,
  index_sigmaH = iter_grid$sigmaH,
  .packages = "tidyverse",
  .combine = rbind,
  .options.snow = opts) %dopar% {
    data.Host %>%
      filter(sigmaH %in% index_sigmaH,
             KH %in% index_KH) %>%
      Topt_heat_func(., system_name)
  }

close(pb)


# Save Topt data
proper_dim <- (dim(data.Topt)[1] * length(unique(data.Vec$system_ID)))
dim(Topt.df)[1] == proper_dim

if (exists("Topt.df") & dim(Topt.df)[1] == proper_dim)
{write_rds(Topt.df, "results/Topt_vals.rds", compress = "gz")
} else {
  warning("No file written. Topt.df either empty or not complete.")
}


### CTmin, CTmax, and CTwidth -------------------------------------------

CT.df <- foreach(
  system_name = iter_grid$system_ID,
  index_KH = iter_grid$KH,
  index_sigmaH = iter_grid$sigmaH,
  .packages = "tidyverse",
  .combine = rbind,
  .options.snow = opts) %dopar% {
    data.CT %>%
      filter(sigmaH %in% index_sigmaH,
             KH %in% index_KH) %>%
      CT_heat_func(., system_name)
  }

close(pb)

# Save CT data
(proper_dim <- 3 * dim(data.CT)[1] * length(unique(data.Vec$system_ID)))

if (exists("CT.df") & dim(CT.df)[1] == proper_dim)
{write_rds(CT.df, "results/CT_vals.rds", compress = "gz")
} else {
  warning("No file written. CT.df either empty or not complete.")
}

stopCluster(cl) 

# 5) Diagnostics & visualizations -----------------------------------------
# 
# # Figure out what's missing from CT.df
# CT.df <- read_rds("results/CT_vals.rds")
# 
# # Make a full list of all the combinations we should have
# full_combo.df <- expand_grid(system_ID = unique(data.Vec$system_ID),
#                              sigmaH = unique(data.Host$sigmaH),
#                              KH = unique(data.Host$KH),
#                              variable = unique(CT.df$variable)
# ) %>%
#   arrange(system_ID, sigmaH, KH, variable) %>% 
#   # combine all columns into a single one for easier comparison
#   unite("all", sep = "_")
# 
# reduce_CT.df <- select(CT.df, system_ID, sigmaH, KH, variable) %>% 
#   arrange(system_ID, sigmaH, KH, variable) %>% 
#   # combine all columns into a single one for easier comparison
#   unite("all", sep = "_")
# 
# # Number of missing entries
# (num_missing_entries = dim(full_combo.df)[1] - dim(reduce_CT.df)[1])
# 
# # Find rows of full_combo.df that don't appear in CT.df
# missing_entries <- full_combo.df$all[which(!(full_combo.df$all %in% reduce_CT.df$all))]
# 
# # Check we caught all the missing entries
# (num_missing_entries == length(missing_entries))
# 
# # Get the parameters of the missing entries
# missing_CT.df <- as.data.frame(missing_entries) %>%
#   rename(united_vals = missing_entries) %>% 
#   separate_wider_delim(united_vals, delim = "_", 
#                        names = c("system_ID", "sigmaH", "KH", "variable")) %>% 
#   select(-variable) %>% 
#   distinct() %>% 
#   mutate(sigmaH = as.double(sigmaH),
#          KH = as.double(KH))
# 
# # Go back and compute CT vals for these entries
# # Set up host trait data frame
# data.CT <- data.Host %>% 
#   select(-c(sigmaH,KH, Model)) %>% 
#   distinct() %>% 
#   cross_join(select(missing_CT.df, sigmaH,KH) %>% distinct()) %>%
#   mutate(Model = ifelse(is.infinite(sigmaH), "Ross-Macdonald model", "Chitnis model"))
# 
# # Slice host trait data 
# sigmaH_slices <- slice(unique(data.CT$sigmaH), 10)
# KH_slices <- slice(unique(data.CT$KH), cluster_size)
# 
# # Initialize CT data frame
# newCT.df <- init.df
# 
# gc()
# # Collect CTmin/max/width data across systems and host trait values
# for (system_name in unique(data.Vec$system_ID)) {
#   print(paste0("CTmin/max/width vals: ", system_name))
#   KHslice_num <- 1
#   for(index_KH in KH_slices) {
#     print(paste0("KH slice number ", KHslice_num, " out of ", length(KH_slices)))
#     sigmaHslice_num <- 1
#     for (index_sigmaH in sigmaH_slices) {
#       print(paste0("sigmaH slice number ", sigmaHslice_num, " out of ", length(sigmaH_slices)))
#       newCT.df <- data.CT %>%
#         filter(sigmaH %in% index_sigmaH,
#                KH %in% index_KH) %>%    
#         CT_heat_func(., system_name) %>%
#         rbind(newCT.df)
#       
#       # Check outputs as we go along
#       test_df <- filter(newCT.df,
#                         system_ID == system_name,
#                         sigmaH %in% index_sigmaH,
#                         KH %in% index_KH) %>% 
#         ungroup() %>% 
#         select(-c(system_ID, Model, lowHCI, highHCI, median)) %>% 
#         distinct() %>% 
#         pivot_wider(names_from = "variable", values_from = "mean") %>% 
#         select(-c(sigmaH, KH)) %>% 
#         distinct()
#       
#       print(test_df)
#       
#       sigmaHslice_num <- sigmaHslice_num  + 1
#     }
#     
#     KHslice_num <- KHslice_num +1
#     gc()
#   }
# }
# # smallest CTmin ~13.85
# # largest CTmax ~34.45
