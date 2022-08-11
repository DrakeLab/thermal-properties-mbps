# Code for performing analyses
output-functions.R
  * Contains functions for calculating model outputs (e.g. R0, equilibrium prevalence)

get-mosquito-traits.R
  * Functions for determining mosquito trait values for the transmission model.
  * Traits vary by mosquito species, pathogen of interest, and temperature.
  * Creates MosquitoThermalResponse.csv in '/results' using data from '/data/clean/MosquitoThermalParameters.csv'

get-analysis-dfs.R
  * Builds data frames for: 1) vector traits, 2) all traits and model outputs, and 3) thermal characteristics of transmission
  * Creates 1) VectorTraits.csv, 2) AllOutputs.csv, and 3) ThermalCharacteristics.csv in '/results'
