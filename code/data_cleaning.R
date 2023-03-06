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
# -- Bayoh2001_mortality.csv (survival data for Anopheles gambiae from Bayoh 2001,
#    taken from Supplementary Data of Mordecai 2013. See Mordecai_2013_lifespan.R for details)
data.Mordecai2013 <- read.csv("data/raw/Mordecai_2013/Mordecai_2013_supp_data.csv", header = TRUE) %>%
  # Add in lifespan data
  rbind(read_csv("data/clean/Bayoh2001_mortality.csv")) %>% 
  mutate(lead_author = "Mordecai") %>%
  mutate(year = "2013") %>%
  select(trait.name, T, trait, mosquito_species, pathogen, lead_author, year)%>%
  unique()

# Mordecai 2017
# -- aegyptiDENVmodelTempData_2016-03-30.csv
# -- albopictusCHIKVmodelTempData_2016-03-26.csv
data.Mordecai2017.Aegypti <- read.csv("data/raw/Mordecai_2017/aegyptiDENVmodelTempData_2016-03-30.csv", header = TRUE) %>%
  # Exclude the Focks & Barrera 2006 data because they're from a model
  filter(ref != "Focks_Barrera_2006_Research&TrainingTropicalDis_Geneva_Paper") %>%
  # Exclude the Rohani et al 2009 data because it has "unrealistically long lifespans"
  filter(ref != "Rohani_et_al_2009_SEJTropMedPH") %>%
  mutate(mosquito_species = "Aedes aegypti") %>%
  mutate(pathogen = ifelse(trait.name %in% c("b", "c", "PDR"), "DENV", NA)) %>%
  select(-c("trait2", "trait2.name"))

data.Mordecai2017.Albopictus <- read.csv("data/raw/Mordecai_2017/albopictusCHIKVmodelTempData_2016-03-26.csv", header = TRUE) %>%
  # The starved mosquitoes had much shorter survival than all other data, so remove them
  filter(trait2 %in% c("sugar-fed", NA)) %>%
  mutate(mosquito_species = "Aedes albopictus") %>%
  mutate(pathogen = ifelse(trait.name %in% c("b", "c", "PDR"), "DENV", NA)) %>%
  select(-c("trait2", "trait2.name"))

data.Mordecai2017.PDRaddl <- read.csv("data/raw/Mordecai_2017/EIP_priors_2015-12-04.csv", header = TRUE) %>%
  rename(pathogen = virus, mosquito_species = mosquito)
  

data.Mordecai2017 <- rbind(data.Mordecai2017.Aegypti, 
                           data.Mordecai2017.Albopictus,
                           data.Mordecai2017.PDRaddl,
                           data.Yang2008) %>%
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
  # Rename pEA traits with the actual trait measured
  mutate(trait.name = ifelse(trait.name == "pEA", notes, trait.name)) %>% 
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
  # NB: MDR, pEA, EFD, and a are same as Mordecai 2017
  filter(!trait.name %in% c("MDR", "pEA", "EFD", "a")) %>% 
  mutate(mosquito_species = "Aedes aegypti") %>%
  mutate(infection_status = stringr::word(rep, 2, 2, sep = "-")) %>%
  mutate(pathogen = case_when(
    trait.name %in% c("bc", "EIR") ~ "ZIKV",
    infection_status == "inf" ~ "ZIKV")
    ) %>%
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

###* Visualize traits as functions of temperature

# Set up data frame for visualization
# Focal systems: 
# - Aedes aegypti / DENV
# - Aedes aegypti / ZIKV
# - Aedes albopictus / DENV
# - Culex quinquefasciatus / WNV
# - Anopheles spp. / Plasmodium falciparum
data.Viz <- data.All %>%
  dplyr::mutate(Species = stringr::word(mosquito_species, 2, 2, sep = " ")) %>%
  dplyr::mutate(Genus = stringr::word(mosquito_species, 1, 1, sep = " ")) %>%
  mutate(species_label = case_when(
    (Genus == "Aedes" & !(Species %in% c("aegypti", "albopictus"))) ~ "spp.",
    (Genus == "Culex" & Species != "quinquefasciatus") ~ "spp.",
    (Genus == "Anopheles" & Species != "gambiae") ~ "spp.",
    TRUE ~ Species
  )) %>%
  mutate(mosquito_species = paste(Genus, species_label, sep = " ")) %>% 
  mutate(pathogen_1 = stringr::word(pathogen, 1, 1, sep = " ")) %>% 
  mutate(pathogen_2 = stringr::word(pathogen, 2, 2, sep = " ")) %>% 
  mutate(pathogen_new = case_when(
    # pathogen %in% c("WNV-SA", "WNV-NY99", "WNV") ~ "WNV",
    pathogen == "DENV" ~ "DENV",
    pathogen == "ZIKV" ~ "ZIKV",
    pathogen %in% c("WNV-NY99", "WNV-SA", "WNV") ~ "WNV",
    (pathogen_1 == "Plasmodium"  & pathogen_2 != "falciparum") ~ "Plasmodium spp.",
    (pathogen_1 == "Plasmodium"  & pathogen_2 == "falciparum") ~ "Plasmodium falciparum",
    pathogen %in% c("SLEV", "MVE") ~ "flavivirus",
    pathogen %in% c("WEEV", "SINV", "EEEV", "RRV") ~ "togavirus",
    pathogen == "RVFV" ~ "togavirus",
    is.na(pathogen) ~ "none",
    TRUE ~ pathogen
    # TRUE ~ NA
  )) %>% 
  unite(
    col = "system_ID",
    c("mosquito_species", "pathogen_new"),
    sep = " / ",
    remove = TRUE
  ) 

# show thermal response of all traits across all systems
trait_plots <- data.Viz %>%
  ggplot(aes(x = T, y = trait, color = as.factor(system_ID), group = system_ID)) +
  geom_point() +
  facet_wrap(~ trait.name, scales = "free")

# show thermal response of selected traits across selected systems
select_trait_plots <- data.Viz %>% 
  filter(system_ID %in% c("Aedes aegypti / DENV", "Aedes aegypti / none", 
                          "Aedes aegypti / ZIKV", "Aedes aegypti / none",
                          "Aedes albopictus / DENV", "Aedes albopictus / none", 
                          "Culex quinquefasciatus / WNV", "Culex quinquefasciatus / none",
                          "Anopheles spp. / Plasmodium falciparum",
                          "Anopheles spp. / none"
                          )) %>% 
  ggplot(aes(x = T, y = trait, color = as.factor(system_ID), group = system_ID)) +
  geom_point() +
  facet_wrap(~ trait.name, scales = "free")

###* Filter to chosen mosquito species and pathogens ----

# Chosen species:
# -- Aedes aegypti, Aedes albopictus, Culex quinquefasciatus, Anopheles spp.
data.Temp <- data.All %>%
  dplyr::mutate(Species = stringr::word(mosquito_species, 2, 2, sep = " ")) %>%
  dplyr::mutate(Genus = stringr::word(mosquito_species, 1, 1, sep = " ")) %>%
  # filter(Genus %in% c("Aedes", "Anopheles", "Culex")) %>%
  # !!! Temporary: if data isn't from the exact species, combine it at the genus level
  # Can use data from these other species to create informative priors for the focal species
  # Anopheles data is always kept at the genus level (?)
  mutate(species_label = case_when(
    (Genus == "Aedes" & !(Species %in% c("aegypti", "albopictus"))) ~ "spp.",
    (Genus == "Culex" & Species != "quinquefasciatus") ~ "spp.",
    (Genus == "Anopheles" & Species != "gambiae") ~ "spp.",
    # Genus == "Anopheles" ~ "spp.",
    TRUE ~ Species
  )) %>%
  mutate(mosquito_species_new = paste(Genus, species_label, sep = " "))

# To get data that exactly matches our chosen species
# data.Only <- filter(data.Temp, mosquito_species == mosquito_species_new)

# 3) Translate traits and trait names to standard form ----


# ) Apply any transformations needed to get the correct trait ----


# ) Assign the proper thermal response function to the trait ----


# ) Assign mosquito species and pathogen names to data ----
