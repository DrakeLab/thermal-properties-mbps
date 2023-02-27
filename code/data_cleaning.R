## Title: Mosquito Thermal Trait Data processing ###############################
##
## Project: Global zoonoses - spillover of mosquito-borne pathogens
##
## Purpose: Translate data from multiple sources into the format used in our
##          analyses
##
## Contents: 1) Set-up, load in necessary packages and data-sets
##           2) Define accessory functions
##           ) Exploratory analysis of the data
##           3) Translate traits and trait names to standard form
##           ) Apply any transformations needed to get the correct trait
##           ) Assign the proper thermal response function to the trait
##           ) Assign mosquito species and pathogen names to data
##
##
## Inputs:  data is found in data/raw/ in folders labeled by the first author
##          and year of publication of the article we obtained the data from
##
##          A table explaining how traits are transformed into the ones used in
##          analyses
##
## Outputs: data - data/clean/ThermalTraitData.csv
##
## Written and maintained by: Kyle Dahlin, kydahlin@gmail.com
## Initialized February 2023

# 1) Set-up, load in necessary packages and data-sets ----
library(tidyverse)

###* Data sets ----
# Mordecai 2013
# -- Mordecai_2013_supp_data.csv
# -- survival_data.csv
data.Mordecai2013 <- read.csv("data/raw/Mordecai_2013/Mordecai_2013_supp_data.csv", header = TRUE) %>%
  mutate(lead_author = "Mordecai") %>%
  mutate(year = "2013") %>%
  select(trait.name, T, trait, mosquito_species, pathogen, lead_author, year)

# Treat survival data separately.
# - Use data from the first day of the experiment, one day before reaching 0.01
#   threshold to three days after the threshold, giving 6 data points
# - Use the exponent(?) of the constant mortality function as an independent 
#   observation of mortality at each temperature
survival.Mordecai2013 <- read.csv("data/raw/Mordecai_2013/survival_data.csv", header = TRUE) %>%
  pivot_longer(cols = X5.C:X40.C, names_to = "T", values_to = "prop.alive") %>% 
  rename(day = D) %>% 
  mutate(T = readr::parse_number(T)) %>% 
  mutate(prop.alive = ifelse(day > 5 & is.na(prop.alive), 0, prop.alive))

survival_fit.df <- tibble(T = numeric(), y0 = numeric(), ymax = numeric(), 
                          k = numeric(), lag = numeric())
for (TT in c(15, 20, 25, 30)) {# unique(survival.Mordecai2013$T)) {
  df <- filter(survival.Mordecai2013, T == TT)
  
  Gomp_fit <- nls(1-prop.alive ~ Gompertz(day, y0, ymax, k, lag),
                  data = df,
                  start = list(y0 = 0.9, ymax = -0.03, k = -0.1, lag = 10))
  Gomp_coefs <- coef(Gomp_fit)
  
  survival_fit.df <- add_row(survival_fit.df, T = TT, y0 = Gomp_coefs["y0"], 
                             ymax = Gomp_coefs["ymax"], k = Gomp_coefs["k"], 
                             lag = Gomp_coefs["lag"])
}

# Gompertz function
Gompertz <- function(x, y0, ymax, k, lag){
  result <- y0 + (ymax -y0)*exp(-exp(k*(lag-x)/(ymax-y0) + 1) )
  return(result)
}

Gomp_df <- tibble(day = unique(survival.Mordecai2013$day)) %>% 
  full_join(survival_fit.df, by = character()) %>% 
  mutate(pred.p = Gompertz(x = day, y0 = y0, ymax = ymax, k = k, lag = lag))

survival_plot <- survival.Mordecai2013 %>% 
  arrange(T, prop.alive) %>% 
  # filter(prop.alive != 0) %>%
  ggplot(aes(x = day, y = 1-prop.alive, color = as.factor(T), group = T)) +
  # ggplot(aes(x = day, y = ((log((-log(1-prop.alive) )/ day))/day)/day, color = as.factor(T), group = T)) +
  # geom_path() +
  geom_point() +
  geom_path(data = Gomp_df, aes(x = day, y = pred.p))

survival_test <- survival.Mordecai2013 %>% 
  # This is roughly linear
  mutate(c = (-(log((-log(1-prop.alive) )/ day)/day))) %>% 
  filter(is.finite(c)) %>% 
  ggplot(aes(x = day, y = c, group = T, color = as.factor(T))) +
  geom_point() +
  # ylim(-0.5, 0.5)


# Mordecai 2017
# -- aegyptiDENVmodelTempData_2016-03-30.csv
# -- albopictusCHIKVmodelTempData_2016-03-26.csv
data.Mordecai2017.Aegypti <- read.csv("data/raw/Mordecai_2017/aegyptiDENVmodelTempData_2016-03-30.csv", header = TRUE) %>%
  # Exclude the Focks & Barrera 2006 data because they're from a model
  filter(ref != "Focks_Barrera_2006_Research&TrainingTropicalDis_Geneva_Paper") %>%
  # Exclude the Rohani et al 2009 data because it has "unrealistically long lifespans"
  filter(ref != "Rohani_et_al_2009_SEJTropMedPH") %>%
  mutate(mosquito_species = "Aedes aegypti") %>%
  mutate(pathogen = "DENV")

data.Mordecai2017.Albopictus <- read.csv("data/raw/Mordecai_2017/albopictusCHIKVmodelTempData_2016-03-26.csv", header = TRUE) %>%
  # The starved mosquitoes had much shorter survival than all other data, so remove them
  filter(trait2 %in% c("sugar-fed", NA)) %>%
  mutate(mosquito_species = "Aedes albopictus") %>%
  mutate(pathogen = "CHIKV")

data.Mordecai2017 <- rbind(data.Mordecai2017.Aegypti, data.Mordecai2017.Albopictus) %>%
  mutate(lead_author = "Mordecai") %>%
  mutate(year = "2017") %>%
  select(trait.name, T, trait, mosquito_species, pathogen, lead_author, year)

# Shocket 2018
# -- RRVPriorData.csv
# -- RRVTraitData.csv
data.Shocket2018 <- read.csv("data/raw/Shocket_2018/RRVTraitData.csv", header = TRUE) %>%
  mutate(mosquito_species = case_when(
    host.code == "Ovig" ~ "Ochlerotatus vigilax",
    host.code == "Cann" ~ "Culex annulirostris",
    host.code == "Anot" ~ "Aedes notoscriptus",
    host.code == "Ocam" ~ "Ochlerotatus camptorhynchus"
  )) %>%
  mutate(lead_author = "Shocket") %>%
  mutate(year = "2018") %>%
  rename(pathogen = paras.code) %>%
  select(trait.name, T, trait, mosquito_species, pathogen, lead_author, year)

# Tesla 2018
# -- InfectionData.csv
# -- SurvivalData.csv
# -- ZIKV_trait_data_fits.csv
# -- zikv_traits.csv
data.Tesla2018 <- read.csv("data/raw/Tesla_2018/zikv_traits.csv", header = TRUE) %>%
  mutate(mosquito_species = "Aedes aegypti") %>%
  mutate(pathogen = "ZIKV") %>%
  mutate(lead_author = "Tesla") %>%
  mutate(year = "2018") %>%
  select(trait.name, T, trait, mosquito_species, pathogen, lead_author, year)

# Mordecai 2019
# -- raw data not provided

# Shocket 2020
# - This article has separate tables for each trait
# -- TraitData_a.csv
# -- TraitData_bc.csv
# -- TraitData_EFD.csv
# -- TraitData_EV.csv
# -- TraitData_lf.csv
# -- TraitData_MDR.csv
# -- TraitData_PDR.csv
# -- TraitData_pLA.csv

data.Shocket2020 <- read.csv("data/raw/Shocket_2020/TraitData_a.csv", header = TRUE) %>%
  rename(trait.name = Trait.Name) %>%
  select(trait.name, T, trait, host.code, paras.code) %>%
  rbind(read.csv("data/raw/Shocket_2020/TraitData_bc.csv", header = TRUE) %>%
          select(trait.name, T, trait, host.code, paras.code)) %>%
  rbind(read.csv("data/raw/Shocket_2020/TraitData_EFD.csv", header = TRUE) %>%
          rename(trait.name = Trait.Name) %>%
          select(trait.name, T, trait, host.code, paras.code)) %>%
  rbind(read.csv("data/raw/Shocket_2020/TraitData_EV.csv", header = TRUE) %>%
          select(trait.name, T, trait, host.code, paras.code)) %>%
  rbind(read.csv("data/raw/Shocket_2020/TraitData_lf.csv", header = TRUE) %>%
          select(trait.name, T, trait, host.code, paras.code)) %>%
  rbind(read.csv("data/raw/Shocket_2020/TraitData_MDR.csv", header = TRUE) %>%
          select(trait.name, T, trait, host.code, paras.code)) %>%
  rbind(read.csv("data/raw/Shocket_2020/TraitData_PDR.csv", header = TRUE) %>%
          rename(trait.name = Trait.Name) %>%
          select(trait.name, T, trait, host.code, paras.code)) %>%
  rbind(read.csv("data/raw/Shocket_2020/TraitData_pLA.csv", header = TRUE) %>%
          select(trait.name, T, trait, host.code, paras.code)) %>%
  rename(pathogen = paras.code) %>%
  mutate(lead_author = "Shocket") %>%
  mutate(year = "2020") %>%
  mutate(mosquito_species = case_when(
    host.code == "Cpip" ~ "Culex pipiens",
    host.code == "Cqui" ~ "Culex quinquefasciatus",
    host.code == "Ctar" ~ "Culex tarsalis",
    host.code == "Cuni" ~ "Culex univittatus",
    host.code == "Cthe" ~ "Culex theileri",
    host.code == "Atri" ~ "Aedes triseriatus",
    host.code == "Atae" ~ "Aedes taeniorhynchus",
    host.code == "Avex" ~ "Aedes vexans",
    host.code == "Cmel" ~ "Culiseta melanura",
    host.code == "Cmol" ~ "Culex pipiens molestus",
    host.code == "Cpal" ~ "Culex pipiens pallens",
    host.code == "Cres" ~ "Culex restuans",
    host.code == "Ador" ~ "Aedes dorsalis",
    host.code == "Anig" ~ "Aedes nigromaculis",
    host.code == "Asol" ~ "Aedes sollicitans",
    host.code == "Asal" ~ "Aedes salinarius",
    TRUE ~ "missing code"
  )) %>%
  filter(mosquito_species != "missing code") %>%
  select(trait.name, T, trait, mosquito_species, pathogen, lead_author, year)

# Combine all data frames
data.All <- rbind(data.Mordecai2013, data.Mordecai2017, data.Shocket2018, data.Tesla2018, data.Shocket2020)

# write.csv(data.All, "data/clean/data.All.csv")

# 2) Define accessory functions ----

# List of traits from data
in_trait_list <- unique(data.All$trait.name)

# List of intermediate traits (ones that TPCs are fit to)
inter_trait_list <- c("a", "TFD", "EFD", "MDR", "e2a", "b", "c", "PDR", "lf")

# Table explaining how to convert from initial to intermediate traits


# List of parameters we need for the model
out_params_list <- c(
  "sigmaV", "fecundity", "deltaL", "lf",
  "KL", "rhoL", "muV", "betaV", "etaV"
)

# Functions for converting from intermediate traits to model parameters


# 3) Exploratory analysis of the data ----

# Check if data is repeated across studies

# Get all unique entries (ignoring the study it comes from)
data.Anon <- select(data.All, -c("lead_author", "year")) %>% unique()

dim(data.All)[1] - dim(data.Anon)[1] # 344 entries are repeated

###* Modify data.All to only keep unique entries ----
# assigning to them the oldest associated reference
data.All <- data.All %>%
  arrange(lead_author, desc(year)) %>%
  distinct(across(trait.name:pathogen), .keep_all = TRUE)

###* Filter to chosen mosquito species and pathogens ----

# Chosen species:
# -- Aedes aegypti, Aedes albopictus, Culex quinquefasciatus, Anopheles spp.
data.Temp <- data.All %>%
  dplyr::mutate(Species = stringr::word(mosquito_species, 2, 2, sep = " ")) %>%
  dplyr::mutate(Genus = stringr::word(mosquito_species, 1, 1, sep = " ")) %>%
  filter(Genus %in% c("Aedes", "Anopheles", "Culex")) %>%
  # Temporary: if data isn't from the exact species, combine it at the genus level
  # Anopheles data is always kept at the genus level (?)
  mutate(species_label = case_when(
    (Genus == "Aedes" & !(Species %in% c("aegypti", "albopictus"))) ~ "spp.",
    (Genus == "Culex" & Species != "quinquefasciatus") ~ "spp.",
    Genus == "Anopheles" ~ "spp.",
    TRUE ~ Species
  )) %>%
  mutate(mosquito_species_new = paste(Genus, species_label, sep = " "))

# To get data that exactly matches our chosen species
# data.Only <- filter(data.Temp, mosquito_species == mosquito_species_new)

# 3) Translate traits and trait names to standard form ----


# ) Apply any transformations needed to get the correct trait ----


# ) Assign the proper thermal response function to the trait ----


# ) Assign mosquito species and pathogen names to data ----
