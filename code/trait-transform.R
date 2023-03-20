## Title: Transform TPC traits into model parameters #######################
##
## Project: Global zoonoses - spillover of mosquito-borne pathogens
##
## Purpose: Transform trait values, which we have fit TPCs to from data, in 
##          order to use in model analyses
##
## Contents: 0) Set-up, load in necessary packages and data-sets
##           1) Define accessory functions
##           2) 
##           3) 
##           4) 
##
##
## Inputs:  data - data/clean/TPC_param_samples.csv
##                 data/clean/trait_transforms.csv
##
## Outputs: data - data/clean/ThermalTraitSamples.csv
##
## Written and maintained by: Kyle Dahlin, kydahlin@gmail.com
## Initialized March 2023
# _______________________________________________________________________________

# 0) Set-up, load in necessary packages and data-sets ---------------------

# Load Libraries
library(tidyverse)

### * Load in data ----
# Trait TPC parameter samples 
data.in <- data.in.transform

# 1) Define accessory functions -------------------------------------------
# Define TPC functions
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
Linear <- function(q, Tmax) {
  function(t) {
    pmax(q * (Tmax - t), 0, na.RM = FALSE)
  }
}

# Function: designate proper thermal response function # !!! This should go somewhere else
# - output is a function of temperature
get.thermal.response <- function(data_in, Temperature) {
  parms <- dplyr::select(data_in, c, T0, Tm)
  function_type <- dplyr::select(data_in, func)
  
  temp_function <- case_when(
    function_type == "Briere" ~ Briere(parms$c, parms$T0, parms$Tm),
    function_type == "Quadratic" ~ Quadratic(parms$c, parms$T0, parms$Tm),
    function_type == "Linear" ~ Linear(parms$c, parms$Tm)
  )
  
  out <- temp_function(Temperature)
}


# 2)  ----

# Function: restrict to focal systems
focal_filter_func <- function(df) {
  df <- filter(df, system_ID %in% c(
  "Aedes aegypti / DENV", "Aedes aegypti / none",
  "Aedes aegypti / ZIKV", "Aedes aegypti / none",
  "Aedes albopictus / DENV", "Aedes albopictus / none",
  "Culex quinquefasciatus / WNV", "Culex quinquefasciatus / none",
  "Anopheles spp. / Plasmodium spp.",
  "Anopheles spp. / none"
))
}

# Define temperature range of study
Temps <- seq(0, 50, length.out = 20)

# Thin samples
thin_size <-  20

# Create data frame of TPCs
TPC_df <- data.in %>%
  filter(sample_num %in% seq(1, thin_size)) %>%
  full_join(list(Temperature = Temps), by = character(), copy = TRUE) %>%
  mutate(Trait_val = case_when(
    func == "Briere" ~ Briere(c, T0, Tm)(Temperature),
    func == "Quadratic" ~ Quadratic(c, T0, Tm)(Temperature),
    func == "Linear" ~ Linear(c, Tm)(Temperature)
  )) %>%
  dplyr::select(-c("c", "T0", "Tm"))

temp_df <- TPC_df %>% 
  pivot_wider(id_cols = c("system_ID", "sample_num", "Temperature"), 
              names_from = "trait",
              values_from = "Trait_val")
# List of model parameters: sigmaV, f, deltaL, rhoL, muV, etaV, betaV

# List of traits with TPCs fit from data: a, bc, PDR, e2a, MDR, EFD, lf, c, b, 
#                                         TFD, pLA, pRH, nLR, EFOC, EPR, pO, EV


# Combine traits into intermediate parameters as necessary
# i.e. putting together reproductive traits to estimate eggs per female per day (EFD)
intermediate_df <- temp_df %>% 
  mutate(e2a = case_when(
    !(is.na(e2a)) ~ e2a,
    !(is.na(pRH * nLR * pLA)) ~ pRH * nLR * pLA,
    !(is.na(EV * pLA)) ~ EV * pLA
  )) %>% 
  mutate(EFD = case_when(
    !(is.na(EFD)) ~ EFD,
    !(is.na(EFOC * a)) ~ EFOC * a,
    !(is.na(EPR * pO)) ~ EPR * pO
  )) %>% 
  mutate(bc = ifelse(is.na(bc), b * c, bc)) %>% 
  # throw out traits that are no longer needed
  dplyr::select(system_ID, sample_num, Temperature,
                a, bc, PDR, e2a, EFD, lf, MDR)

intermediate_df <- intermediate_df %>% 
  # separate out mosquito species and pathogen names
  mutate(mosquito_species = stringr::word(system_ID, 1, 2)) %>% 
  mutate(pathogen = stringr::word(system_ID, 4, 4)) %>% 
  # "forget" system_ID for now
  dplyr::select(-system_ID)

# Combine parasite relevant data with mosquito life history
noInfection_df <- intermediate_df %>% 
  filter(pathogen == "none") %>% 
  dplyr::select(-pathogen) %>% 
  dplyr::select(-c("bc", "PDR"))

Infection_df <- intermediate_df %>% 
  filter(pathogen != "none") %>% 
  dplyr::select(-c("a", "e2a", "EFD", "lf", "MDR"))

combined_df <- left_join(Infection_df, noInfection_df)%>% 
  unite(
    col = "system_ID",
    c("mosquito_species", "pathogen"),
    sep = " / ",
    remove = FALSE
  ) %>% 
  relocate(system_ID, mosquito_species, pathogen, Temperature, sample_num)


# Check what data we're missing (we'll go back and use substitutes for these)
missing_traits_df <- combined_df %>% 
  pivot_longer(cols = bc:MDR) %>% 
  dplyr::filter(is.na(value)) %>% 
  dplyr::select(-c(sample_num, Temperature)) %>% 
  unique()

# Following Shocket 2020: use Culex univittatus / WNV / bc for Culex quinquefasciatus / WNV / bc


# Deal with any duplicates: What do we do if we have two estimates for the same intermediate parameter?
# i.e. We have e2a for Culex but also pO and EV



# 3)  --------------


# *) Diagnostics & visualizations -----------------------------------------
