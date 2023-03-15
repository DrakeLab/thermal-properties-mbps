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
##          code - code/Mordecai_2017/mcmc_utils_all.R
##                 code/Mordecai_2017/temp_functions_all.R
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
# Trait data
data.in <- read_csv("data/clean/data_for_TPC_fitting.csv") # !!! remove this later so this can interact with the main-analysis.R script

# Load table of TPC forms to use for traits
TPC_forms <- read_csv("data/clean/trait_TPC_forms.csv") %>%
  # when TPC function can't be found from literature (a total of 10)
  # use Briere (less restrictive than quadratic)
  mutate(TPC.function = ifelse(is.na(TPC.function), "Briere", TPC.function))

# Set up data frame of traits and mosquito pathogen pairs
# full_df <- expand_grid(Mosquito = MosqPathPairs, trait = as_tibble(trait_names))

# This table shows which functional form to use for each trait
trait_lookup_table <- read_csv("data/clean/trait_transforms.csv")


# 1) Define accessory functions -------------------------------------------

# Function: get samples of thermal trait parameters
thermtrait.prior.sample <- function(data_in, trait_in, mosquito_in, pathogen_in,
                                    n.chains = 5, n.adapt = 5000, n.samps = 5000,
                                    old_informative = FALSE) {
  with(as.list(data_in), {
    # restrict to the trait we care about
    data <- filter(data_in, trait.name == trait_in)
    
    # filter to just the trait of interest
    # might have to bind multiple traits together, use look-up table
    # !!! still necessary?
    # transform_table <- read_csv("data/clean/trait_transforms.csv") %>%
    #   filter(trait.to == trait_in)
    
    # Get the proper TPC function
    TPC_function <- TPC_forms %>%
      # For cases with multiple possible TPC functions, choose the one used in later studies
      mutate(TPC.function = str_extract(TPC.function, "[^/]+$")) %>% 
      filter(trait.name == trait_in) %>%
      filter(mosquito_species == mosquito_in) %>%
      filter(pathogen == pathogen_in) %>%
      dplyr::select(TPC.function)
    
    # trait_name <- unique(transform_table$trait.name)
    
    # from_names <- c(trait_in, unique(transform_table$trait.from)) %>% unique()
    
    # initial values
    inits_list <- if (TPC_function == "Briere") {
      list(Tm = 31, T0 = 5, c = 0.00007)
    } else if (TPC_function == "Quadratic") {
      list(T0 = 5, Tm = 33, n.qd = 0.005)
    } else if (TPC_function == "Linear") {
      list(T0 = 5, Tm = 40, c = 0.005)
    }
    
    # names of TPC parameters being fit
    variable_names <- if (TPC_function == "Briere") {
      c("c", "Tm", "T0", "sigma")
    } else if (TPC_function == "Quadratic") {
      c("n.qd", "Tm", "T0", "sigma")
    } else if (TPC_function == "Linear") {
      c("c", "Tm", "T0", "sigma")
    }
    
    # If you want to use the hyperparameters from Mordecai et al., 2017, load them in
    if (old_informative == TRUE) {
      
      prev_hypers <- read_csv("data/clean/gamma_fits.csv") %>%
        filter(trait %in% trait_in) %>%
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
      
      if (dim(prev_hypers)[1] == 0) {stop("No saved TPC hyperparameter data. Switch old_informative to false")}
      
      jags_choice <- case_when(
        TPC_function == "Briere" ~ "code/jags-models/jags-briere-informative.bug",
        TPC_function == "Quadratic"  ~ "code/jags-models/jags-quad-neg-informative.bug",
        TPC_function == "Linear" ~ "code/jags-models/jags-linear-informative.bug"
      )
      
      jags_data <-  list("Y" = data$trait, "T" = data$T, 
                         "N" = length(data$T), "hypers" = prev_hypers)
      
      samps <- run.jags(jags_data, TPC_function, variable_names,
                        jags_choice, inits_list,
                        n.chains, n.adapt, n.samps)
      
      # Otherwise:
      # 1) use default priors
      # 2) update these using any related species
      # 3) sequentially update with more recent datasets
    } else {
      
      # Gather data from other species of the same genus
      other_species <- data %>%
        filter(stringr::word(mosquito_species, 1, 1) == stringr::word(mosquito_in, 1, 1)) %>%
        filter(mosquito_species != mosquito_in)
      
      # if such data is unavailable, use data from species outside of the genus
      if (dim(other_species)[1] == 0) {
        other_species <- data %>%
          filter(mosquito_species == "other spp.")
      }
      
      # First, inform prior distribution of TPC hyperparameters using data from other species
      if (dim(other_species)[1] > 0) { 
        
        # # !!! placeholder, to deal with inability to create informed priors
        # if (mosquito_in == "Culex quinquefasciatus" & pathogen_in == "none" & trait_in == "EPR") {
        #   assumed_data <- tibble(T = c(10,50), trait = c(0,0))
        #   other_species <- other_species %>% 
        #     dplyr::select(trait, T) %>% 
        #     rbind(assumed_data)
        # }
        # 
        # # Add zeroes at extreme temperatures to help with convergence of fitdistr below
        # assumed_data <- tibble(T = c(0,60), trait = c(0,0))
        # other_species <- other_species %>% 
        #   dplyr::select(trait, T) %>% 
        #   rbind(assumed_data)
        
        other_data <-  list("Y" = other_species$trait, "T" = other_species$T, 
                            "N" = length(other_species$T))
        
        # Select the appropriate bugs model
        jags_other <- case_when(
          TPC_function == "Briere" ~ "code/jags-models/jags-briere.bug",
          TPC_function == "Quadratic" ~ "code/jags-models/jags-quad-neg.bug",
          TPC_function == "Linear" ~  "code/jags-models/jags-linear.bug"
        )
        
        other_samps <- run.jags(other_data, TPC_function, variable_names,
                                jags_other, inits_list,
                                n.chains, n.adapt, 1000)
        prev_hypers = apply(other_samps, 2, function(df) fitdistr(df, "gamma")$estimate)
        
      } else {prev_hypers = c()} # if no trait data is available for any other species, start with uninformed priors
      
      # List of studies with data reported for the particular trait and mosquito species
      other_studies <- data %>%
        filter(mosquito_species == mosquito_in) %>%
        arrange(year, lead_author) %>% 
        distinct(lead_author, year)
      
      # Select the appropriate bug model
      jags_choice <- case_when(
        TPC_function == "Briere" & is.null(prev_hypers) ~ "code/jags-models/jags-briere.bug",
        TPC_function == "Briere" & !is.null(prev_hypers) ~ "code/jags-models/jags-briere-informative.bug",
        TPC_function == "Quadratic" & is.null(prev_hypers) ~ "code/jags-models/jags-quad-neg.bug",
        TPC_function == "Quadratic" & !is.null(prev_hypers) ~ "code/jags-models/jags-quad-neg-informative.bug",
        TPC_function == "Linear" & is.null(prev_hypers) ~ "code/jags-models/jags-linear.bug",
        TPC_function == "Linear" & !is.null(prev_hypers) ~ "code/jags-models/jags-linear-informative.bug"
      )
      
      # Initialize the previous hyperparameters
      
      for (ii in dim(other_studies)[1]) {
        # Identify study
        index_author <- other_studies$lead_author[ii]
        index_year <- other_studies$year[ii]
        
        data_temp <- filter(data, lead_author == index_author, year == index_year)
        
        jags_data <-  if (is.null(prev_hypers)) {
          list("Y" = data_temp$trait, "T" = data_temp$T, 
               "N" = length(data_temp$T))
        } else {
          list("Y" = data_temp$trait, "T" = data_temp$T, 
               "N" = length(data_temp$T), "hypers" = prev_hypers)
        }
        
        samps <- run.jags(jags_data, TPC_function, variable_names,
                          jags_choice, inits_list,
                          n.chains, n.adapt, n.samps)
        
        if (ii < dim(other_studies)[1]) {
          # for older entries, generate posterior distributions to estimate 
          # hyperparameters for future fitting
          
          prev_hypers = apply(samps, 2,
                              function(df) fitdistr(df, "gamma")$estimate)
          
          rm(out.samps) # !!! remove after debugging
          # for the final entry, just generate samples from the jags.model
        } else {
          samps <- mutate(samps, func = as.character(TPC_function))
          
          if (is.null(samps$T0)) {samps$T0 <- NA}
          
        }
        
      }
    }
    
    out.samps <- dplyr::select(samps, -c(sigma, tau, func))
    
    return(out.samps)
  })
}

# Function: get samples of trait TPC hyperparameters from running jags
run.jags <- function(jags_data, TPC_function, variable_names,
                     jags_choice, inits_list,
                     n.chains, n.adapt, n.samps) {
  # create MCMC samples from the model with default priors
  jags <- jags.model(jags_choice,
                     data = jags_data,
                     n.chains = n.chains, inits = inits_list,
                     n.adapt = n.adapt,
                     quiet = TRUE # switch to FALSE to show messages and progress bars
  )
  # The coda.samples() function takes n.samps new samples, and saves
  # them in the coda format, which we use for visualization and
  # analysis.
  coda.samps <- coda.samples(jags, variable_names, n.samps)
  
  # This command combines the samples from the n.chains into a format
  # that will be used for further analyses.
  if (TPC_function == "Briere") {
    samps <- make.briere.samps(coda.samps, nchains = n.chains, samp.lims = c(1, n.samps))
  } else if (TPC_function == "Quadratic") {
    samps <- make.quad.samps(coda.samps, nchains = n.chains, samp.lims = c(1, n.samps), sig = TRUE)
  } else {
    samps <- make.linear.samps(coda.samps, nchains = n.chains, samp.lims = c(1, n.samps), sig = TRUE) %>% 
      mutate(Tm = n.inter/slope) %>% 
      mutate(c = slope) %>% 
      dplyr::select(c, Tm, sigma)
  }
  
  samps$tau <- 1 / samps$sigma
  if (is.null(samps$c)) {
    samps <- mutate(samps, c = qd, .keep = "unused")
  }
  return(samps)
}

# 2) Calculate prior distributions of thermal trait parameters from data ----

### Set parameters for Monte Carlo sampling ----

# Specify the parameters that control the MCMC (these will be used throughout the code).
n.chains <- 2 # 5
n.adapt <- 100 # 5000
n.samps <- 100 # 5000

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

# 3) Save thermal trait parameter distributions ---------------------------

# Combine samples
samples.All <- tibble(Species = "Aedes albopictus", samples.Albopictus) %>%
  rbind(tibble(Species = "Aedes aegypti", samples.Aegypti)) %>%
  # reorder columns to match old parameter table
  relocate(Species, trait, c, T0, Tm)

write_csv(samples.All, "data/clean/ThermalTraitSamples.csv")

# 4) Functions to sample from thermal trait parameter distributions -------
