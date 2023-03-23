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


# 2) Create TPCs and deal with missing data ----

# Create data frame of TPCs
TPC_df <- data.in %>%
  filter(sample_num %in% seq(1, thin_size)) %>%
  cross_join(list(Temperature = Temps), copy = TRUE) %>%
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
    # !(is.na(pRH * nLR * pLA)) ~ pRH * nLR * pLA,
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
  filter(system_ID %in% c(
    "Aedes aegypti / DENV", "Aedes aegypti / none",
    "Aedes aegypti / ZIKV", "Aedes aegypti / none",
    "Aedes albopictus / DENV", "Aedes albopictus / none",
    "Culex quinquefasciatus / WNV", "Culex quinquefasciatus / none",
    "Anopheles spp. / Plasmodium spp.",
    "Anopheles spp. / none"
  )) %>% 
  unique()

# Following Shocket 2020: use Culex univittatus / WNV / bc for Culex quinquefasciatus / WNV / bc
combined_df <- combined_df %>% 
  # temporarily remove Cx. quinquefasciatus / WNV rows
  filter(system_ID != "Culex quinquefasciatus / WNV") %>% 
  # add rows back in after switching out bc for Cx. univittatus / WNV rows
  rbind(filter(combined_df, system_ID == "Culex quinquefasciatus / WNV") %>% 
          # remove original bc values (all NA)
          dplyr::select(-bc) %>% 
          # join with Cx. univittatus data
          right_join(filter(combined_df, system_ID == "Culex univittatus / WNV") %>% 
                       dplyr::select(Temperature:bc))) %>% 
  # remove Cx. univittatus data
  filter(system_ID != "Culex univittatus / WNV") 

# Deal with any duplicates: What do we do if we have two estimates for the same intermediate parameter?
# i.e. We have e2a for Culex but also pO and EV


# 3) Make parameter dataframe ---------------------------------------------
eps <- .Machine$double.eps

parameter_df <- combined_df %>%
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
  # mutate(muL = ifelse(deltaL != 0, rhoL * (1 - deltaL) / (deltaL + eps), Inf)) %>%  
  # Adult mosquito mortality rate. This should be ~infinite if lf = 0
  mutate(muV = ifelse(lf != 0, 1/(lf + eps), Inf)) %>% 
  # Pathogen development rate.
  mutate(etaV = PDR) %>%
  # Mosquito infection probability.
  mutate(betaV = bc) %>% 
  dplyr::select(system_ID:sample_num, sigmaV:betaV)



# 4) Save parameter data frame --------------------------------------------

write_rds(parameter_df, "data/clean/parameter_TPCs.rds")



# *) Diagnostics & visualizations -----------------------------------------

plot_bool = FALSE
if (plot_bool) {
library(cowplot)
# For each mosquito species, trait, and sample, get a thermal response curve
TPC_df <- parameter_df %>% 
  ungroup() %>% 
  mutate(lf = 1/muV, .keep = "unused") %>%
  dplyr::select(-c(muL, etaL)) %>% 
  # mutate(lfL = 1/etaL, .keep = "unused") %>%
  # select(-c())
  melt(id = c("system_ID", "mosquito_species","pathogen", "Temperature", "sample_num"),
       variable.name = "trait",
       value.name = "Trait_val") #%>% 
  # mutate(Trait_val = ifelse(is.infinite(Trait_val), 0, Trait_val)) %>%
  #remove unreasonably long lifespan values
  # mutate(Trait_val = ifelse(trait %in% c("lfL", "lf"), ifelse(Trait_val>200, NA, Trait_val), Trait_val))

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
