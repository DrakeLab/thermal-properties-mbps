## Title: Visualization of thermal traits ######################################
##
## Project: Global zoonoses - spillover of mosquito-borne pathogens
##
## Purpose: Produce visualizations of thermal trait distributions for a sanity check
##
## Contents: 1) Set-up, load in necessary packages and data-sets
##           2) Histograms of parameter distributions of thermal traits
##           3) Plots of thermal trait functions with 89% HCI
##           4)
##
##
## Inputs:  data - data/clean/gamma_fits.csv
##
##          code - code/Mordecai2017/mcmc_utils_all.R
##                 code/Mordecai2017/temp_functions_all.R
##
## Outputs: functions:
##          data - data/clean/ThermalTraitSamples.csv
##
## Written and maintained by: Kyle Dahlin, kydahlin@gmail.com
## Initialized February 2023

# 1) Set-up, load in necessary packages and data-sets ----
### Load Libraries ----
library(tidyverse)
library(reshape2)
library(cowplot)

### Load in data and functions ----
# Load samples of thermal trait parameters
samples.All <- read.csv("data/clean/ThermalTraitSamples.csv")

# Load functions from Mordecai et al., 2017
# This file contains tools for analysis and visualization.
source("code/Mordecai2017/mcmc_utils_all.R")

# This file contains the thermal response functions and their derivatives.
source("code/Mordecai2017/temp_functions_all.R")

# 2) Histograms of parameter distributions of thermal traits ----
plot_df <- samples.All %>%
  dplyr::select(-func) %>% 
  melt(id = c("Species", "trait", "sample_num"))

###* Figure: thermal trait parameter posterior distributions ----
parm_hists <- plot_df %>%
  ggplot(aes(value, color = Species, fill = Species)) +
  geom_histogram(aes(), bins = 100) +
  # geom_density() +
  facet_wrap(trait ~ variable, scales = "free") +
  theme_minimal_grid(16)

# 3) Plots of thermal trait functions with 89% HCI ----

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

### Temperature vector used for visualizations ----
Temps <- seq(10, 45, length.out = 100)

# Thinning intervals for samples
res_reduce <- function(df, new_length) {
  old_length <- length(unique(df))
  if (old_length < new_length) {
    warning("new_length is larger than old length. Vector will be unchanged")
    new_length <- old_length
  }
  ret <- unique(df)[seq(1, length(unique(df)), length.out = new_length)]
}

thin_size <- 100


# For each mosquito species, trait, and sample, get a thermal response curve
TPC_df <- samples.All %>%
  filter(sample_num %in% res_reduce(sample_num, thin_size)) %>%
  full_join(list(Temperature = Temps), by = character(), copy = TRUE) %>%
  # group_by(sample_num) %>%
  # try sapply or lapply
  # mutate(Trait_val = get.thermal.response(.,Temperature))Linear(parms$c, parms$T0)
  mutate(Trait_val = case_when(
    func == "Briere" ~ Briere(c, T0, Tm)(Temperature),
    func == "Quadratic" ~ Quadratic(c, T0, Tm)(Temperature),
    func == "Linear" ~ Linear(c, T0)(Temperature)
  )) %>% 
  dplyr::select(-c("c", "T0", "Tm"))

# get mean TPC from samples
meanTPC_df <- TPC_df %>% 
  group_by(Species, trait, Temperature) %>% 
  summarise(mean_val = mean(Trait_val), .groups = "keep")


# get edges of 89% HCI of samples
quantsTPC_df <- TPC_df %>% 
  group_by(Species, trait, Temperature) %>% 
  mutate(lowHCI_val = quantile(Trait_val, 0.055)) %>% 
  mutate(highHCI_val = quantile(Trait_val, 0.945)) %>% 
  dplyr::select(-c("sample_num", "Trait_val", "func"))


###* Figure: TPC curves with 89% high confidence intervals ---- 
TPC_plot <- TPC_df %>%
  group_by(sample_num) %>% 
  arrange(Temperature) %>% 
  filter(Species == "Aedes albopictus") %>% 
  # group_by()
  ggplot(aes(x = Temperature, y = Trait_val, color = Species)) +
  geom_line(data = meanTPC_df, aes(x = Temperature, y = mean_val)) +
  geom_line(data = quantsTPC_df, aes(x = Temperature, y = lowHCI_val),
            linewidth = 1, linetype = "dashed") +
  geom_line(data = quantsTPC_df, aes(x = Temperature, y = highHCI_val),
            linewidth = 1, linetype = "dashed") +
  # geom_point() +
  facet_grid(rows = vars(trait), scales = "free") +
  theme_minimal_grid(16)


# 4)  ----

