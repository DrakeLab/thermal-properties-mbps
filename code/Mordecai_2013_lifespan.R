## Title: Fitting lifespan data from Mordecai 2013 ###############################
##
## Project: Global zoonoses - spillover of mosquito-borne pathogens
##
## Purpose: Estimate mortality rate from data provided in Bayoh 2001
##
## Contents: 1) Set-up, load in necessary packages and data-sets
##           2) Load in and clean data as necessary
##           3) Filter data to values prescribed in Mordecai 2013 supplement
##           4) Fit exponential survival curves
##           5) Get estimates of mu and translate to lifespan
##           6) Output data in correct format for further processing (to data_cleaning.R)
##
##
## Inputs:  data/raw/Mordecai_2013/Mordecai_2013_supp_data.csv
##
##
## Outputs: data - Mordecai_2013_lifespan.Rdata
##
## Written and maintained by: Kyle Dahlin, kydahlin@gmail.com
## Initialized March 2023


# 1) Set-up load in necessary packages ----
require(tidyverse)

# 2) Load in and clean data as necessary ----
###* Using data from (Bayoh 2001). 
mortality.data <- read.csv("data/raw/Mordecai_2013/survival_data.csv", header = TRUE) %>% 
  # change columns to be temperature entries
  pivot_longer()

  
# 3) Filter data to values prescribed in Mordecai 2013 supplement ----

###* They determined the final day where the population proportion exceeded 1%.


###* They then filtered the data to only include the following days: 
###* the first day of the experiment, one day preceding reaching the 1% 
###* threshold, the day where the 1% threshold is met, then the three days 
###* following. 

# 4) Fit exponential survival curves ----


# 5) Get estimates of mu and translate to lifespan ----


# 6) Output data in correct format for further processing (to data_cleaning.R) ----







###* 
###* The parameter of the exponential survival function was used for 𝜇, an 
###* independent observation of mortality at each temperature. 
###* This captured temperature-associated trends observed in the data, but the 
###* shape of these fitted curves were substantially different from those 
###* derived from a (more complex) Gompertz curve fit, especially between days 
###* 20 and 40. 
###* 
###* Future work could use a partial differential equation model of mosquito 
###* population dynamics with age as an additional parameter to overcome this 
###* misalignment between the data and the conventional assumption of constant 
###* mortality with age.
