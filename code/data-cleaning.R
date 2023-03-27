## Title: Mosquito Thermal Trait Data processing ###############################
##
## Project: Global zoonoses - spillover of mosquito-borne pathogens
##
## Purpose: Translate data from multiple sources into the format used in our
##          analyses
##
## Contents: 0) Set-up, load in necessary packages
##           1) Load in all data sets
##           2) Combine data sets
##           3) Rename non-focal focal systems and synonymous traits
##           4) Data visualizations / diagnostics
##           5) Export data set
##
##
## Inputs:  data is found in data/raw/ in folders labeled by the first author
##          and year of publication of the article we obtained the data from
##
##          A table explaining how traits are transformed into the ones used in
##          analyses
##
## Outputs: data - data/clean/ThermalTraitData.rds
##
## Written and maintained by: Kyle Dahlin, kydahlin@gmail.com
## Initialized February 2023


# 0) Set-up, load in necessary packages and data-sets ---------------------
library(tidyverse)


# 1) Load in all data sets ------------------------------------------------

### * Mordecai 2013 ----
# -- Mordecai_2013_supp_data.csv
# -- Bayoh2001_mortality.csv (survival data for Anopheles gambiae from Bayoh 2001,
#    taken from Supplementary Data of Mordecai 2013. See Bayoh2001-lifespan-analysis.R for details)
data.Mordecai2013 <- read.csv("data/raw/Mordecai_2013/Mordecai_2013_supp_data.csv", header = TRUE) %>%
  # Add in lifespan data
  rbind(read_csv("data/clean/Bayoh2001_mortality.csv")) %>% 
  mutate(lead_author = "Mordecai") %>%
  # Make bc a probability
  mutate(trait = ifelse(trait.name == "bc", trait / 100, trait)) %>% 
  # Fix EIP entry errors (EIP = 0 really means parasite never developed)
  mutate(trait = ifelse(trait.name == "EIP" & trait == 0, Inf, trait)) %>% 
  mutate(year = "2013") %>%
  dplyr::select(trait.name, T, trait, mosquito_species, pathogen, lead_author, year) %>%
  unique()

### * Mordecai 2017 ----
# -- aegyptiDENVmodelTempData_2016-03-30.csv
# -- albopictusCHIKVmodelTempData_2016-03-26.csv
data.Mordecai2017.Aegypti <- read.csv("data/raw/Mordecai_2017/aegyptiDENVmodelTempData_2016-03-30.csv", header = TRUE) %>%
  # Exclude the Focks & Barrera 2006 data because they're from a model
  filter(ref != "Focks_Barrera_2006_Research&TrainingTropicalDis_Geneva_Paper") %>%
  # Exclude the Rohani et al 2009 data because it has "unrealistically long lifespans"
  filter(ref != "Rohani_et_al_2009_SEJTropMedPH") %>%
  mutate(mosquito_species = "Aedes aegypti") %>%
  mutate(pathogen = ifelse(trait.name %in% c("b", "c", "PDR", "EIP"), case_when(
    ref == "Alto&Bettinardi_2013_AJTMH" ~ "DENV",
    ref == "Carrington_et_al_2013_PNTD" ~ "DENV",
    ref == "Davis_1932_AmJEpidemiology" ~ "YFV",
    ref == "Focks_et_al_1995_AJTMH" ~ "DENV",
    ref == "Lambrects_et_al_2011_PNAS" ~ "DENV",
    ref == "McLeah_et_al_1975_MosquitoNews" ~ "DENV",
    ref == "McLean_et_al_1974_CanJMicobiol" ~ "DENV",
    ref == "Watts_et_al_1987_AJTMH" ~ "DENV"
    ), NA)) %>% 
  # mutate(pathogen = ifelse(trait.name %in% c("b", "c", "PDR"), "DENV", NA)) %>%
  dplyr::select(-c("trait2", "trait2.name"))

data.Mordecai2017.Albopictus <- read.csv("data/raw/Mordecai_2017/albopictusCHIKVmodelTempData_2016-03-26.csv", header = TRUE) %>%
  # The starved mosquitoes had much shorter survival than all other data, so remove them
  filter(trait2 %in% c("sugar-fed", NA)) %>%
  # This study only considered Aedes aegypti
  filter(ref != "Lambrects_et_al_2011_PNAS") %>% 
  mutate(mosquito_species = "Aedes albopictus") %>%
  mutate(pathogen = ifelse(trait.name %in% c("b", "c", "PDR", "EIP"), case_when(
    ref == "Alto&Bettinardi_2013_AJTMH" ~ "DENV",
    ref == "Westbrook_et_al_2010_VecBorn&ZooDis" ~ "CHIKV",
    ref == "Westbrook_Thesis_2010" ~ "CHIKV",
    ref == "Xiao_et_al_2014_Arch Virol" ~ "DENV"
  ), NA)) %>% 
  # mutate(pathogen = ifelse(trait.name %in% c("b", "c", "PDR"), "CHIKV", NA)) %>%
  dplyr::select(-c("trait2", "trait2.name"))

data.Mordecai2017.PDRaddl <- read.csv("data/raw/Mordecai_2017/EIP_priors_2015-12-04.csv", header = TRUE) %>%
  rename(pathogen = virus, mosquito_species = mosquito)
  
data.Mordecai2017 <- rbind(data.Mordecai2017.Aegypti, 
                           data.Mordecai2017.Albopictus,
                           data.Mordecai2017.PDRaddl) %>%
  mutate(lead_author = "Mordecai") %>%
  mutate(year = "2017") %>%
  dplyr::select(trait.name, T, trait, mosquito_species, pathogen, lead_author, year)

### *  Shocket 2018 ----
# -- RRVPriorData.csv
# -- RRVTraitData.csv

# Process EFD data separately
data.EFD.proc <- read.csv("data/raw/Shocket_2018/RRVTraitData.csv", header = TRUE) %>%
  # remove McDonald_1980_AusJofEco EFD data which will be dealt with separately
  filter(ref == "McDonald_1980_AusJofEco" & trait.name == "EFD") %>% 
  mutate(trait2 = as.double(trait2)) %>% 
  group_by(T) %>%
  mutate(trait = sum(trait)/max(trait2)) %>%
  dplyr::select(trait.name, T, trait, ref, trait2, trait2.name, host.code, paras.code, notes) %>% 
  unique() %>% 
  ungroup()

data.Shocket2018 <- read.csv("data/raw/Shocket_2018/RRVTraitData.csv", header = TRUE) %>%
  # remove EFD data from Williams_2011_AusJEntrom which did not count actual eggs laid 
  filter(!(ref == "Williams_2011_AusJEntrom" & trait.name == "EFD")) %>% 
  # remove McDonald_1980_AusJofEco EFD data which will be dealt with separately
  filter(!(ref == "McDonald_1980_AusJofEco" & trait.name == "EFD")) %>% 
  # add in processed McDonald_1980_AusJofEco EFD data
  rbind(data.EFD.proc) %>% 
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
  dplyr::select(trait.name, T, trait, mosquito_species, pathogen, lead_author, year)

### * Tesla 2018 ----
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
  dplyr::select(trait.name, T, trait, mosquito_species, pathogen, lead_author, year)

### * Mordecai 2019 ----
# -- raw data not provided

### * Shocket 2020 ----
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
  dplyr::select(trait.name, T, trait, host.code, paras.code) %>%
  rbind(read.csv("data/raw/Shocket_2020/TraitData_bc.csv", header = TRUE) %>%
          dplyr::select(trait.name, T, trait, host.code, paras.code)) %>%
  rbind(read.csv("data/raw/Shocket_2020/TraitData_EFD.csv", header = TRUE) %>%
          rename(trait.name = Trait.Name) %>%
          dplyr::select(trait.name, T, trait, host.code, paras.code)) %>%
  rbind(read.csv("data/raw/Shocket_2020/TraitData_EV.csv", header = TRUE) %>%
          dplyr::select(trait.name, T, trait, host.code, paras.code)) %>%
  rbind(read.csv("data/raw/Shocket_2020/TraitData_lf.csv", header = TRUE) %>%
          dplyr::select(trait.name, T, trait, host.code, paras.code)) %>%
  rbind(read.csv("data/raw/Shocket_2020/TraitData_MDR.csv", header = TRUE) %>%
          dplyr::select(trait.name, T, trait, host.code, paras.code)) %>%
  rbind(read.csv("data/raw/Shocket_2020/TraitData_PDR.csv", header = TRUE) %>%
          rename(trait.name = Trait.Name) %>%
          dplyr::select(trait.name, T, trait, host.code, paras.code)) %>%
  rbind(read.csv("data/raw/Shocket_2020/TraitData_pLA.csv", header = TRUE) %>%
          dplyr::select(trait.name, T, trait, host.code, paras.code)) %>%
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
  dplyr::select(trait.name, T, trait, mosquito_species, pathogen, lead_author, year)

### * Anopheles fecundity data ----
data.Anga_fecundity <- read.csv("data/raw/Anopheles_fecundity.csv", header = TRUE) %>%
  mutate(lead_author = "Dahlin") %>%
  mutate(year = "2023") %>%
  mutate(mosquito_species = case_when(
    host.code == "Anga" ~ "Anopheles gambiae",
    host.code == "Ansu" ~ "Anopheles superpictus"
  )) %>%
  mutate(pathogen = NA) %>% 
  dplyr::select(trait.name, T, trait, mosquito_species, pathogen, lead_author, year)

# 2) Combine data sets ----------------------------------------------------
# Combine all data frames
data.All <- rbind(data.Mordecai2013, data.Mordecai2017, data.Shocket2018, 
                  data.Tesla2018, data.Shocket2020, data.Anga_fecundity)

# 3) Rename non-focal systems and combine synonymous traits -----------------

# List of focal disease systems: 
# - Aedes aegypti / DENV
# - Aedes aegypti / ZIKV
# - Aedes albopictus / DENV
# - Culex quinquefasciatus / WNV
# - Anopheles spp. / Plasmodium falciparum

# Reduce all non-focal mosquito and parasite species to the genus level
data.Reduced <- data.All %>%
  # Separate mosquito binomial
  mutate(Species = stringr::word(mosquito_species, 2, 2, sep = " ")) %>%
  mutate(Genus = stringr::word(mosquito_species, 1, 1, sep = " ")) %>%
  mutate(Genus = ifelse(Genus %in% c("Aedes", "Culex", "Anopheles"), Genus, "Other")) %>% 
  # Reduce non-focal MOSQUITO species to genus level
  mutate(species_label = case_when(
    (Genus == "Aedes" & !(Species %in% c("aegypti", "albopictus"))) ~ "spp.",
    (Genus == "Culex" & !(Species %in% "quinquefasciatus")) ~ "spp.",
    Genus == "Anopheles" ~ "spp.",
    # (Genus == "Anopheles" & Species != "gambiae") ~ "spp.", # Not enough data for Anopheles gambiae alone to fit model, expand to genus level
    Genus == "Other" ~ "spp.",
    TRUE ~ Species
  )) %>%
  mutate(mosquito_species = paste(Genus, species_label, sep = " ")) %>% 
  # Separate parasite binomial
  mutate(pathogen_1 = stringr::word(pathogen, 1, 1, sep = " ")) %>% 
  mutate(pathogen_2 = stringr::word(pathogen, 2, 2, sep = " ")) %>% 
  # Reduce non-focal PARASITES to genus level
  mutate(pathogen = case_when(
    pathogen == "DENV" ~ "DENV",
    pathogen == "ZIKV" ~ "ZIKV",
    pathogen %in% c("WNV-NY99", "WNV-SA", "WNV") ~ "WNV",
    pathogen_1 == "Plasmodium" ~ "Plasmodium spp.",
    # (pathogen_1 == "Plasmodium"  & pathogen_2 != "falciparum") ~ "Plasmodium spp.", # Not enough data for Plasmodium falciparum alone to fit model.
    # (pathogen_1 == "Plasmodium"  & pathogen_2 == "falciparum") ~ "Plasmodium falciparum",
    pathogen %in% c("SLEV", "MVE", "YFV") ~ "other flavivirus",
    pathogen %in% c("RVFV", "WEEV", "SINV", "EEEV", "RRV", "CHIKV") ~ "other togavirus", # might want to separate out CHIKV later
    is.na(pathogen) ~ "none",
    TRUE ~ pathogen
    # TRUE ~ NA
  )) %>% 
  unite(
    col = "system_ID",
    c("mosquito_species", "pathogen"),
    sep = " / ",
    remove = FALSE
  ) %>% 
  dplyr::select(-c(Species, Genus, species_label, pathogen_1, pathogen_2))

# Load the table explaining how to convert from initial to intermediate traits
transform_table <- read_csv("data/clean/trait_transforms.csv") %>% 
  dplyr::select(-notes) %>% 
  rename(trait.name = trait.from) %>% 
  arrange(final.trait, trait.to, trait.name)

# Define small value for inverting traits
eps <- .Machine$double.eps

# Rename synonymous traits
data.in.TPC <- data.Reduced %>% 
  right_join(transform_table) %>% 
  mutate(trait = case_when(
    transform == "Identity" ~ trait,
    transform == "Inverse" ~ 1/(trait + eps),
    transform == "NegativeLogDifference" ~ 1/(-log(1-trait)), # following Mordecai et al., 2017 for prop.dead data from Alto et al., 2001b
    NA ~ trait,
    TRUE ~ trait
    # Others will be transformed after fitting thermal performance curves
  )) %>% 
  mutate(trait.from = trait.name) %>% 
  mutate(trait.name = trait.to, .keep = "unused") %>% 
  filter(!is.na(trait))

# Combine lifespan data for Aedes aegypti with and without ZIKV since we don't account for virulence in this model
data.in.TPC <- data.in.TPC %>% 
  mutate(pathogen = ifelse(mosquito_species == "Aedes aegypti" & final.trait == "lf",
    "none",
    pathogen
  ))

# 4) Save and export data set ---------------------------------------------

# Save data.frame to file
write_rds(data.in.TPC, "data/clean/data_for_TPC_fitting.rds")



# 5) Data visualizations / diagnostics ------------------------------------

###* Visualize traits as functions of temperature
plot_bool <- FALSE # decide whether you'd like to generate a diagnostic plot

if (plot_bool) {
library(cowplot)
# Set up data frame for visualization
data.Viz <- data.Out

# show thermal response of all traits across all systems
trait_plots <- data.Viz %>%
  ggplot(aes(x = T, y = trait, color = as.factor(mosquito_species),
             shape = as.factor(pathogen), group = system_ID)) +
  geom_point() +
  scale_shape_manual(name = "Pathogen",
                     values = 1:length(unique(data.Viz$pathogen))) +
  scale_color_discrete(name = "Mosquito species") +
  labs(x = "Temperature",
       y = "Trait value") +
  facet_wrap(~ trait.name, scales = "free") +
  theme_minimal(16)

# show thermal response data for focal systems only
select_trait_plots <- data.Viz %>% 
  filter(system_ID %in% c("Aedes aegypti / DENV", "Aedes aegypti / none", 
                          "Aedes aegypti / ZIKV", "Aedes aegypti / none",
                          "Aedes albopictus / DENV", "Aedes albopictus / none", 
                          "Culex quinquefasciatus / WNV", "Culex quinquefasciatus / none",
                          "Anopheles spp. / Plasmodium spp.",
                          "Anopheles spp. / none"
  )) %>%
  ggplot(aes(x = T, y = trait, color = as.factor(mosquito_species),
             shape = as.factor(pathogen), group = system_ID)) +
  geom_point() +
  scale_shape_manual(name = "Pathogen",
                     values = 1:length(unique(data.Viz$pathogen))) +
  scale_color_discrete(name = "Mosquito species") +
  labs(x = "Temperature",
       y = "Trait value") +
  facet_wrap(~ trait.name, scales = "free") +
  theme_minimal(16)
}