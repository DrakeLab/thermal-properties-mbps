## Title: Calculate Topt and critical thermal extrema #######################
##
## Project: Global zoonoses - spillover of mosquito-borne pathogens
##
## Purpose: Transform trait values, which we have fit TPCs to from data, in
##          order to use in model analyses
##
## Contents: 0) Set-up, load in necessary packages and data-sets
##           1) Define accessory functions
##           2) Calculate R0 TPCs
##           3) Calculate Topt
##           4) Calculate Topt alternate: restrict to R0 > 1
##           5) Calculate CTmin, CTmax, and CTwidth
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

eps <- .Machine$double.eps

# Initialize a data frame to save our analyses
init.df <- tibble(system_ID = c(), Temperature = c(), Model = c(),
                  sigmaH = c(), KH = c(), variable = c(),
                  lowHCI = c(), highHCI = c(), mean = c(), median = c())

# Thin out the vector data set as necessary

thin_size <- 300


data.Vec <- read_rds("results/VecTPC_vals.rds")
samples <- unique(data.Vec$sample_num)
num_samples <- length(samples)
sample_inds <- sample(samples, min(num_samples, thin_size), replace = FALSE)

data.Vec <- data.Vec%>% 
  filter(sample_num %in% sample_inds)


# 1) Define accessory functions -------------------------------------------

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

# ALTERNATE Function: evaluate Topt ONLY WHEN R0>1 and its mean, median, and highest prob. intervals across
#           samples
Topt_heat_func_restricted <- function(in_df, system_name) {
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
    filter(R0 > 1) %>% 
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

# 1) Set *host* parameters --------------------------------

## Set resolution for host trait variation
# Host density vector: Number of values to include to consider for vertebrate host density
KH_vec_length <- 300 # full = 300, thin = 20

# Biting tolerance vector: Number of values to consider for biting tolerance
sigmaH_vec_length <- 300 # full = 300, thin = 20


## Host life history & behavioral traits
# Host recruitment rate:
# upper estimate from range for Primate traits: (0.001150685, 0.009624300)
lambdaH_baseline <- .005

# Host mortality rate:
# estimate from range for Primate traits: lifespan (8.6, 60) years
# muH_vec <- 1 / (365 * c(1, 25))
muH_baseline <- 1 / (365 * 20)

# Host maximum biting tolerance (mosquitoes bites per day)
sigmaH_vec <- 10^seq(-0.25, 3.25, length.out = sigmaH_vec_length - 7) %>%
  c(1, 10, 20, 50, 100, 1000, Inf) %>%
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
  betaH = betaH_baseline) %>%
  as.data.frame() %>%
  mutate(Model = ifelse(is.infinite(sigmaH), "Ross-Macdonald model", "Chitnis model"))




# 2) Calculate R0 TPCs -----------------------------------------

# Set up parallel
if (!exists("cluster")) {
  cluster_size <- min(21, parallel::detectCores() - 2)
  cluster <- new_cluster(cluster_size)
  cluster_library(cluster, c("dplyr", "tidyr"))
}

# Set up host trait data frame (for future visualization)
data.R0 <- as.data.frame(data.Host) %>%
  # filter(sigmaH %in% c(1e-1, 1, 10, 100, Inf)) %>%
  filter(sigmaH %in% c(10^seq(-1,2), Inf,
                       unique(sigmaH)[seq(1, length(unique(sigmaH)), length.out = 51)])) %>%
  filter(KH %in% c(10^seq(-2,5) ,
                   unique(KH)[seq(1, length(unique(KH)), length.out = 51)]))

# Slice host trait data
sigmaH_slices <- slice(unique(data.R0$sigmaH), 2)
KH_slices <- slice(unique(data.R0$KH), cluster_size)

# Initialize R0 data frame
R0.df <- init.df

gc()

# Collect R0 TPC data across systems and host trait values
for (system_name in unique(data.Vec$system_ID)) {
  print(paste0("R0 vals: ", system_name))
  KHslice_num <- 1
  pb <- progress_bar$new(
    format = ":spin :system progress = :percent [:bar] :elapsedfull | eta: :eta",
    total = length(KH_slices) * length(sigmaH_slices),
    width = 120)  
  for(index_KH in KH_slices) {
    pb$tick()
    # print(paste0("KH slice number ", KHslice_num, " out of ", length(KH_slices)))
    sigmaHslice_num <- 1
    for (index_sigmaH in sigmaH_slices) {
      pb$tick()
      # print(paste0("sigmaH slice number ", sigmaHslice_num, " out of ", length(sigmaH_slices)))
      temp_df <- data.R0 %>%
        filter(sigmaH %in% index_sigmaH,
               KH %in% index_KH) %>%
        R0_TPC_func(., system_name)
      
      R0.df <- rbind(temp_df, R0.df)
      sigmaHslice_num <- sigmaHslice_num  + 1
    }
    KHslice_num <- KHslice_num +1
    gc()
  }
}
pb$terminate()

# Save data
proper_dim <- dim(data.R0)[1] * length(unique(data.Vec$system_ID)) * length(unique(data.Vec$Temperature))

(dim(R0.df)[1] == proper_dim)

R0_TPCs.df <- R0.df %>% 
  ungroup() %>% 
  filter(sigmaH %in% c(100, Inf))
write_rds(R0_TPCs.df, "results/R0_TPC_vals.rds", compress = "gz")

R0_heat.df <- R0.df %>% 
  ungroup() %>% 
  filter(KH %in% c(0.1, 1, 10, 100))
write_rds(R0_heat.df, "results/R0_vals.rds", compress = "gz")


# 3) Calculate R0 derivatives ---------------------------------------------

# Set up parallel
if (!exists("cluster")) {
  cluster_size <- min(21, parallel::detectCores() - 2)
  cluster <- new_cluster(cluster_size)
  cluster_library(cluster, c("dplyr", "tidyr"))
}

# Set up host trait data frame (for future visualization)
data.dR0dk <- data.Host %>%
  filter(sigmaH %in% c(100, Inf)) #%>%
  # filter(sigmaH %in% c(10^seq(-1,2), Inf,
  #                      unique(sigmaH)[seq(1, length(unique(sigmaH)), length.out = 51)])) %>%
  # filter(KH %in% c(10^seq(-2,5)))
# filter(KH %in% c(10^seq(-2,5) ,
#                  unique(KH)[seq(1, length(unique(KH)), length.out = 51)]))

# Slice host trait data
sigmaH_slices <- slice(unique(data.dR0dk$sigmaH), 2)
KH_slices <- slice(unique(data.dR0dk$KH), cluster_size)

# Initialize R0 data frame
dR0dvar.df <- init.df

gc()

temp_vars <- data.Vec %>% 
  select(-c(system_ID, mosquito_species, pathogen, Temperature, sample_num,
            etaL, muL, KL, V0,)) %>% 
  colnames()

R0_deriv_func <- function(in_df, var_name) {
  out_df <- in_df %>% # Pr(surviving the EIP)
    mutate(thetaV = exp(-1 / (lf * etaV + eps))) %>% 
    # Host constant
    mutate(C = betaH / (KH * (gammaH + muH))) %>% 
    # Biting rate term
    mutate(K = ifelse(is.infinite(sigmaH), 1, sigmaH * KH / (sigmaH * KH + sigmaV * V0))) %>% 
    # Basic reproduction number
    mutate(R0 = K * sqrt(C * sigmaV^2 * V0 * betaV * thetaV * lf)) %>% 
    # Derivative of biting rate term wrt vector abundance
    mutate(dKdV0 = -sigmaV * K^2 / (sigmaH * KH)) %>% 
    # Derivative of vector abundance wrt variable
    mutate(dV0dvar = case_when(
      var_name == "lf" ~ rhoL * KL,
      var_name == "sigmaV_f" ~ (rhoL * KL)/ ((sigmaV_f^2) * deltaL + eps),
      var_name == "deltaL" ~ (rhoL * KL)/ (sigmaV_f * (deltaL^2) + eps),
      var_name == "rhoL" ~ (KL * (1 - deltaL))/ (sigmaV_f * deltaL + eps),
      TRUE ~ 0
    )) %>% 
    # Derivative of R0 wrt variable
    mutate(dR0dk = (C / (2 * R0 + eps)) * case_when(
      var_name == "lf" ~ sigmaV^2 * betaV * K * thetaV * (2 * dKdV0 * dV0dvar * V0 * lf + dV0dvar * K * lf + K * V0 * lf / (etaV * (lf^2) + eps) + K * V0),
      var_name == "sigmaV" ~ 2 * V0 * betaV * thetaV * lf * K * sigmaV * (K + sigmaV * dKdV0 * dV0dvar),
      var_name == "sigmaV_f" ~ sigmaV^2 * betaV * thetaV * lf * K * dV0dvar * (2 * dKdV0 * V0 + K),
      var_name == "deltaL" ~ sigmaV^2 * betaV * thetaV * lf * K * dV0dvar * (2 * dKdV0 * V0 + K),
      var_name == "rhoL" ~ sigmaV^2 * betaV * thetaV * lf * K * dV0dvar * (2 * dKdV0 * V0 + K),
      var_name == "etaV" ~ K^2 *sigmaV^2 * V0 * betaV * thetaV / (etaV^2 + eps),
      var_name == "betaV" ~ K^2 *sigmaV^2 * V0 * thetaV * lf,
      TRUE ~ NA
    ))
}

# Collect R0 TPC data across systems and host trait values
print_index = 1
for (system_name in unique(data.Vec$system_ID)) {
  for (var_name in temp_vars) {
    print(paste0("(",print_index, "/", length(temp_vars) * length(unique(data.Vec$system_ID)),") dR0 / d", var_name,": ", system_name))
    KHslice_num <- 1
    pb <- progress_bar$new(
      format = ":spin :system progress = :percent [:bar] :elapsedfull | eta: :eta",
      total = length(KH_slices) * length(sigmaH_slices),
      width = 120)  
    for(index_KH in KH_slices) {
      for (index_sigmaH in sigmaH_slices) {
        pb$tick()
        
        temp_df <- data.dR0dk %>%
          filter(sigmaH %in% index_sigmaH,
                 KH %in% index_KH) %>%
          expand_grid(filter(data.Vec, system_ID == system_name)) %>%
          # Calculate derivatives
          R0_deriv_func(., var_name) %>% 
          mutate(focal_var = var_name) %>% 
          dplyr::select(system_ID, sample_num, Temperature, Model, sigmaH, KH, focal_var, dR0dk) %>%
          group_by(system_ID, Temperature, Model, sigmaH, KH) %>%
          partition(cluster) %>%
          summarise(
            lowHCI = quantile(dR0dk, 0.055),
            highHCI = quantile(dR0dk, 0.945),
            mean = mean(dR0dk),
            median = median(dR0dk)
          ) %>%
          collect()
        
        dR0dvar.df <- rbind(temp_df, dR0dvar.df)
      }
    }
    gc()
    print_index <- print_index + 1
  }
}
pb$terminate()

# Save data
proper_dim <- dim(data.dR0dk)[1] * length(unique(data.Vec$system_ID)) * length(unique(data.Vec$Temperature)) * length(temp_vars)

(dim(dR0dvar.df)[1] == proper_dim)

write_rds(dR0dvar.df, "results/dR0dk_vals.rds", compress = "gz")


# 4) Calculate Topt -------------------------------------------------------

# Set up new cluster
# Close old cluster connections
rm(cluster)
gc()
# Start new cluster for doParallel
cluster_size <- parallel::detectCores()-1

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
proper_dim <- (dim(iter_grid)[1])
dim(Topt.df)[1] == proper_dim

if (exists("Topt.df") & dim(Topt.df)[1] == proper_dim)
{write_rds(Topt.df, "results/Topt_vals.rds", compress = "gz")
} else {
  warning("No file written. Topt.df either empty or not complete.")
}


# 4) Calculate Topt alternate: restrict to R0 > 1 -------------------------


Topt_alt.df <- foreach(
  system_name = iter_grid$system_ID,
  index_KH = iter_grid$KH,
  index_sigmaH = iter_grid$sigmaH,
  .packages = "tidyverse",
  .combine = rbind,
  .options.snow = opts) %dopar% {
    data.Host %>%
      filter(sigmaH %in% index_sigmaH,
             KH %in% index_KH) %>%
      Topt_heat_func_restricted(., system_name)
  }

close(pb)


# Save Topt data
proper_dim <- (dim(iter_grid)[1])
dim(Topt.df)[1] == proper_dim

if (exists("Topt_alt.df") & dim(Topt_alt.df)[1] == proper_dim)
{write_rds(Topt_alt.df, "results/Topt_alt_vals.rds", compress = "gz")
} else {
  warning("No file written. Topt.df either empty or not complete.")
}


# 5) Calculate CTmin, CTmax, and CTwidth ----------------------------------
# Set up new cluster
# Close old cluster connections
rm(cluster)
gc()
# Start new cluster for doParallel
cluster_size <- parallel::detectCores()-1

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

CT.df <- foreach(
  system_name = iter_grid$system_ID,
  index_KH = iter_grid$KH,
  index_sigmaH = iter_grid$sigmaH,
  .packages = "tidyverse",
  .combine = rbind,
  .options.snow = opts) %dopar% {
    data.Host %>%
      filter(sigmaH %in% index_sigmaH,
             KH %in% index_KH) %>%
      CT_heat_func(., system_name)
  }

close(pb)

# Save CT data
(proper_dim <- 3 * dim(iter_grid)[1])

(proper_dim == dim(CT.df)[1])

if (exists("CT.df") & dim(CT.df)[1] == proper_dim)
{write_rds(CT.df, "results/CT_vals.rds", compress = "gz")
} else {
  warning("No file written. CT.df either empty or not complete.")
}

stopCluster(cl) 

# *) Diagnostics & visualizations -----------------------------------------
