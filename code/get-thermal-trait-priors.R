## Title: Prior distributions of mosquito thermal traits #######################
##
## Project: Global zoonoses - spillover of mosquito-borne pathogens
##
## Purpose: Compute prior distributions for thermal trait parameters and provide
##          functions for sampling from these distributions
##
## Contents: 0) Set-up, load in necessary packages and data-sets
##           1) Define accessory functions
##           2) Calculate prior distributions of thermal trait parameters from data
##           3) Save thermal trait parameter distributions
##           4) Functions to sample from thermal trait parameter distributions
##
##
## Inputs:  data - data/clean/data_for_TPC_fitting.csv
##                 data/clean/trait_table.csv
##                 data/clean/gamma_fits.csv
##
##          code - code/Mordecai2017/mcmc_utils_all.R
##                 code/Mordecai2017/temp_functions_all.R
##
## Outputs: functions:
##          data - data/clean/ThermalTraitSamples.csv
##
## Written and maintained by: Kyle Dahlin, kydahlin@gmail.com
## Initialized February 2023
## Modified from original code provided in the following articles:
# * Mordecai, E. A., J. M. Caldwell, M. K. Grossman, C. A. Lippi, L. R. Johnson, M. Neira, J. R. Rohr, S. J. Ryan, V. Savage, M. S. Shocket, R. Sippy, A. M. Stewart Ibarra, M. B. Thomas, and O. Villena. 2019. Thermal biology of mosquito-borne disease. Ecology letters 22:1690–1708.
# * Mordecai, E. A., J. M. Cohen, M. V. Evans, P. Gudapati, L. R. Johnson, C. A. Lippi, K. Miazgowicz, C. C. Murdock, J. R. Rohr, S. J. Ryan, V. Savage, M. S. Shocket, A. Stewart Ibarra, M. B. Thomas, and D. P. Weikel. 2017. Detecting the impact of temperature on transmission of Zika, dengue, and chikungunya using mechanistic models. PLoS neglected tropical diseases 11:e0005568.
# * Mordecai, E. A., K. P. Paaijmans, L. R. Johnson, C. Balzer, T. Ben-Horin, E. de Moor, A. McNally, S. Pawar, S. J. Ryan, T. C. Smith, and K. D. Lafferty. 2013. Optimal temperature for malaria transmission is dramatically lower than previously predicted. Ecology letters 16:22–30.
# * Shocket, M. S., S. J. Ryan, and E. A. Mordecai. 2018. Temperature explains broad patterns of Ross River virus transmission. eLife 7.
# * Shocket, M. S., A. B. Verwillow, M. G. Numazu, H. Slamani, J. M. Cohen, F. El Moustaid, J. Rohr, L. R. Johnson, and E. A. Mordecai. 2020. Transmission of West Nile and five other temperate mosquito-borne viruses peaks at temperatures between 23°C and 26°C. eLife 9.
# * Tesla, B., L. R. Demakovsky, E. A. Mordecai, S. J. Ryan, M. H. Bonds, C. N. Ngonghala, M. A. Brindley, and C. C. Murdock. 2018. Temperature drives Zika virus transmission: evidence from empirical and mathematical models. Proceedings. Biological sciences / The Royal Society 285:20180795.
# _______________________________________________________________________________

# 0) Set-up, load in necessary packages and data-sets ---------------------

# Load Libraries
library(tidyverse)
library(IDPmisc)
library("rjags") # Make sure you have installed JAGS-4.x.y.exe (for any x >=0, y>=0) from http://www.sourceforge.net/projects/mcmc-jags/files
library("MASS")
library("here")

# Load utility functions (from Mordecai et al., 2017)
# This file contains tools for analysis and visualization.
source("code/Mordecai_2017/mcmc_utils_all.R")

# This file contains the thermal response functions and their derivatives.
source("code/Mordecai_2017/temp_functions_all.R")

### * Load in data ----
data.in <- read.csv("data/clean/data_for_TPC_fitting.csv") %>% dplyr::select(-X)


traits <- tibble(trait.name = unique(data.in$trait.name))
systems 

trait_table <- distinct(data.in, trait.name, mosquito_species, pathogen) %>% 
  arrange(trait.name, mosquito_species)
write_csv(trait_table, "data/clean/temp_trait_table.csv")

# Set up data frame of traits and mosquito pathogen pairs
full_df <- expand_grid(Mosquito = MosqPathPairs, trait = as_tibble(trait_names))

# This table shows which functional form to use for each trait
trait_lookup_table <- read.csv("data/clean/trait_table.csv", header = TRUE)

### Set parameters for Monte Carlo sampling ----

# Specify the parameters that control the MCMC (these will be used throughout the code).
n.chains <-5 # 2 # 5
n.adapt <- 5000 # 100 # 5000
n.samps <- 5000 # 100 # 5000

# 2) Define accessory functions ----

# Function: get samples of thermal trait parameters
thermtrait.prior.sample <- function(data_in, trait_name,
                                    n.chains = 5, n.adapt = 5000, n.samps = 5000) {
  with(as.list(data_in), {
    # data_in should be specified to the species and pathogen level

    # filter to just the trait of interest
    # might have to bind multiple traits together, use look-up table
    lookup_table <- read.csv("data/clean/trait_table.csv", header = TRUE) %>%
      filter(trait.name == trait_name | computed.from == trait_name)

    # trait_name <- unique(lookup_table$trait.name)

    from_names <- c(trait_name, unique(lookup_table$computed.from))

    data <- filter(data_in, trait.name %in% from_names)

    # determine appropriate model using lookup table
    thermal_function <- unique(lookup_table$thermal.function)

    jags_choice <- case_when(
      thermal_function == "Briere" ~ "code/jags-models/jags-briere-informative.bug",
      thermal_function == "Quadratic" ~ "code/jags-models/jags-quad-neg-informative.bug"
    )
    # initial values
    inits_list <- if (thermal_function == "Briere") {
      list(Tm = 31, T0 = 5, c = 0.00007)
    } else if (thermal_function == "Quadratic") {
      list(T0 = 5, Tm = 33, n.qd = 0.005)
    }
    # names of variables being fit
    variable.names <- case_when(
      thermal_function == "Briere" ~ c("c", "Tm", "T0", "sigma"),
      thermal_function == "Quadratic" ~ c("n.qd", "Tm", "T0", "sigma")
    )

    # load in hyperparameters from prior fits
    hypers_table <- read.csv("data/clean/gamma_fits.csv", header = TRUE) %>%
      filter(trait %in% from_names) %>%
      unique() %>% 
      # adjust by the appropriate multiplier
      mutate(value = multiplier * value) %>%
      dplyr::select(-multiplier) %>%
      # all the rest of this is to get the data back in the form that "jags" wants
      dplyr::select(-trait) %>%
      pivot_wider(names_from = Var2, values_from = value) %>%
      as.data.frame() %>%
      `rownames<-`(.[, 1]) %>%
      dplyr::select(-Var1)
    
    print(hypers_table)

    # create MCMC samples from the model with default priors
    jags <- jags.model(jags_choice,
      data = list("Y" = data$trait, "T" = data$T, "N" = length(data$T), "hypers" = hypers_table),
      n.chains = n.chains, inits = inits_list,
      n.adapt = n.adapt
    )
    # The coda.samples() function takes n.samps new samples, and saves
    # them in the coda format, which we use for visualization and
    # analysis.
    coda.samps <- coda.samples(jags, variable.names, n.samps)

    # This command combines the samples from the n.chains into a format
    # that will be used for further analyses.
    if (thermal_function == "Briere") {
      samps <- make.briere.samps(coda.samps, nchains = n.chains, samp.lims = c(1, n.samps))
    } else if (thermal_function == "Quadratic") {
      samps <- make.quad.samps(coda.samps, nchains = n.chains, samp.lims = c(1, n.samps), sig = TRUE)
    } else {
      samps <- make.linear.samps(coda.samps, nchains = n.chains, samp.lims = c(1, n.samps), sig = TRUE)
    }

    samps$tau <- 1 / samps$sigma
    if (is.null(samps$c)) {
      samps$c <- samps$qd
    }
    samps$func <- thermal_function
    out.samps <- dplyr::select(samps, T0, Tm, c, func)

    return(out.samps)
  })
}


# 3) Calculate prior distributions of thermal trait parameters from data ----

# !!! I will come back and generalize this to take any choice of mosquito species

# Aedes albopictus
temp.Albopictus <- data.Albopictus %>%
  thermtrait.transform()

samples.Albopictus <- tibble(
  trait = as.character(),
  T0 = as.double(),
  Tm = as.double(),
  c = as.double() # !!! note that we're using c as a generic parameter for Briere or Quadratic
)

for (trait_ID in unique(temp.Albopictus$trait.name)) {
  samples <- thermtrait.prior.sample(
    temp.Albopictus, trait_ID,
    n.chains, n.adapt, n.samps
  ) %>%
    mutate(trait = trait_ID) %>%
    mutate(sample_num = row_number())

  samples.Albopictus <- rbind(samples.Albopictus, samples)
}

# Aedes aegypti
temp.Aegypti <- data.Aegypti %>%
  thermtrait.transform()

samples.Aegypti <- tibble(
  trait = as.character(),
  T0 = as.double(),
  Tm = as.double(),
  c = as.double() # !!! note that we're using c as a generic parameter for Briere or Quadratic
)

for (trait_ID in unique(temp.Aegypti$trait.name)) {
  samples <- thermtrait.prior.sample(
    temp.Aegypti, trait_ID,
    n.chains, n.adapt, n.samps
  ) %>%
    mutate(trait = trait_ID) %>%
    mutate(sample_num = row_number())

  samples.Aegypti <- rbind(samples.Aegypti, samples)
}

# Combine samples
samples.All <- tibble(Species = "Aedes albopictus", samples.Albopictus) %>%
  rbind(tibble(Species = "Aedes aegypti", samples.Aegypti)) %>%
  # reorder columns to match old parameter table
  relocate(Species, trait, c, T0, Tm)

write_csv(samples.All, "data/clean/ThermalTraitSamples.csv")
