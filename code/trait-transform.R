## Title: Transform TPC traits into model parameters #######################
##
## Project: Global zoonoses - spillover of mosquito-borne pathogens
##
## Purpose: Transform trait values, which we have fit TPCs to from data, in 
##          order to use in model analyses
##
## Contents: 0) Set-up, load in necessary packages and data-sets
##           1) Define accessory functions
##           2) Create TPCs and deal with missing data
##           3) Create parameters data frame
##           4) Save parameters data frame
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


# 2) Create TPCs and deal with missing data ----

# choose a random set of samples from the TPC parameters
num_samples <- length(unique(data.in.transform$sample_num))
sample_inds <- sample(1:num_samples, thin_size, replace = FALSE)

# Create data frame of TPCs
TPC_df <- data.in.transform %>%
  filter(sample_num %in% sample_inds) %>%
  cross_join(list(Temperature = Temps), copy = TRUE) %>%
  mutate(Trait_val = case_when(
    func == "Briere" ~ Briere(c, T0, Tm)(Temperature),
    func == "Quadratic" ~ Quadratic(c, T0, Tm)(Temperature),
    func == "Linear" ~ Linear(c, Tm)(Temperature)
  )) %>%
  dplyr::select(-c("c", "T0", "Tm")) %>% 
  pivot_wider(id_cols = c("system_ID", "sample_num", "Temperature"), 
              names_from = "trait",
              values_from = "Trait_val")  %>% 
  # Combine traits into intermediate parameters as necessary
  # i.e. putting together reproductive traits to estimate eggs per female per day (EFD)
  mutate(e2a = case_when(
    !(is.na(e2a)) ~ e2a,
    # !(is.na(pRH * nLR * pLA)) ~ pRH * nLR * pLA,
    !(is.na(EV * pLA)) ~ EV * pLA
  )) %>% 
  mutate(EFD = case_when(
    !(is.na(EFD)) ~ EFD,
    !(is.na(EFOC * a)) ~ EFOC * a,
    !(is.na(EPR * pO)) ~ EV * EPR * pO * a
  )) %>% 
  mutate(bc = ifelse(is.na(bc), b * c, bc)) %>% 
  # throw out traits that are no longer needed
  dplyr::select(system_ID, sample_num, Temperature,
                a, bc, PDR, e2a, EFD, lf, MDR) %>% 
  # separate out mosquito species and pathogen names
  separate_wider_delim(system_ID, delim = " / ", names = c("mosquito_species","pathogen"))

# Combine parasite relevant traits with mosquito life history traits according to system
noInfection_df <- TPC_df %>% 
  filter(pathogen == "none") %>% 
  dplyr::select(-pathogen) %>% 
  dplyr::select(-c("bc", "PDR"))

Infection_df <- TPC_df %>% 
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
  filter(system_ID %in% c(
    "Aedes aegypti / DENV", "Aedes aegypti / none",
    "Aedes aegypti / ZIKV", "Aedes aegypti / none",
    "Aedes albopictus / DENV", "Aedes albopictus / none",
    "Culex quinquefasciatus / WNV", "Culex quinquefasciatus / none",
    "Anopheles gambiae / Plasmodium falciparum",
    "Anopheles gambiae / none"
  )) %>% 
  unique()
# Should just show: Culex quinquefasciatus / WNV / bc

# Similar to Shocket 2020: use Culex spp. / WNV / bc data for Culex quinquefasciatus / WNV / bc
# This combines data for Culex univittatus, tarsalis, and pipiens
combined_df <- combined_df %>% 
  # temporarily remove Cx. quinquefasciatus / WNV rows
  filter(system_ID != "Culex quinquefasciatus / WNV") %>% 
  # add rows back in after switching in "bc" values from Cx. univittatus / WNV
  rbind(filter(combined_df, system_ID == "Culex quinquefasciatus / WNV") %>% 
          # remove original bc values (all NA)
          dplyr::select(-bc) %>% 
          # join with Cx. univittatus data
          right_join(filter(combined_df, system_ID == "Culex spp. / WNV") %>% 
                       dplyr::select(Temperature:bc))) %>% 
  # remove Cx. univittatus data
  filter(system_ID != "Culex spp. / WNV") 

# Deal with any duplicates: What do we do if we have two estimates for the same intermediate parameter?
# i.e. We have e2a for Culex but also pO and EV


# 3) Make parameter dataframe ---------------------------------------------
eps <- .Machine$double.eps

data.in.params <- combined_df %>%
  # Biting rate.
  mutate(sigmaV = a) %>% 
  # Fecundity. Simplified from a * 0.5 * EFD / a
  mutate(sigmaV_f = 0.5 * EFD) %>% # = female eggs per female per day
  # Aquatic-stage mosquito survival probability.
  mutate(deltaL = e2a) %>% 
  # Aquatic-stage mosquito development rate.
  mutate(rhoL = MDR) %>% 
  mutate(etaL = rhoL / (deltaL + eps)) %>%
  # Aquatic-stage mosquito mortality rate. This should be ~infinite if deltaL = 0
  mutate(muL = etaL - rhoL) %>%   
  # Adult mosquito average lifespan
  mutate(lf = lf) %>% 
  # Pathogen development rate.
  mutate(etaV = PDR) %>%
  # Mosquito infection probability.
  mutate(betaV = bc) %>% 
  dplyr::select(system_ID:sample_num, lf, sigmaV:betaV)


# 4) Save parameter data frame --------------------------------------------


# *) Diagnostics & visualizations -----------------------------------------

if (plot_bool) {
library(cowplot)
# For each mosquito species, trait, and sample, get a thermal response curve
TPC_df <- data.in.params %>% 
  ungroup() %>% 
  dplyr::select(-c(muL, etaL)) %>%
  melt(id = c("system_ID", "mosquito_species","pathogen", "Temperature", "sample_num"),
       variable.name = "trait",
       value.name = "Trait_val") %>% 
  filter(system_ID %in% c("Aedes aegypti / DENV", "Aedes aegypti / ZIKV"))

# get mean TPC from samples
meanTPC_df <- TPC_df %>%
  group_by(system_ID, trait, Temperature) %>%
  summarise(mean_val = mean(Trait_val), .groups = "keep") %>%
  unique() %>% ungroup()

# get edges of 89% HCI of samples # !!! not calculated correctly?
quantsTPC_df <- TPC_df %>%
  group_by(system_ID, trait, Temperature) %>%
  mutate(lowHCI_val = quantile(Trait_val, 0.055, na.rm = TRUE)) %>%
  mutate(highHCI_val = quantile(Trait_val, 0.945, na.rm = TRUE)) %>%
  dplyr::select(-c("sample_num", "Trait_val")) %>%
  unique() %>% ungroup()

TPC_plot <- TPC_df %>%
  # group_by(sample_num) %>%
  arrange(Temperature) %>%
  # group_by()
  ggplot() +
  # means of TPC curves
  geom_line(
    data = meanTPC_df,
    aes(x = Temperature, y = mean_val, color = system_ID)
  ) +
  # 89% HCI of TPC curves
  geom_ribbon(
    data = quantsTPC_df,
    aes(x = Temperature, ymin = lowHCI_val, ymax = highHCI_val, fill = system_ID),
    alpha = 0.1
  ) +
  ylab("") +
  facet_wrap(~trait, scales = "free", ncol = 2) +
  theme_minimal_grid(12)

# Save figure
ggsave("figures/param_TPC_plot.svg",
       plot = TPC_plot,
       device = "svg",
       width = 16, height = 9, units = "in")
}
