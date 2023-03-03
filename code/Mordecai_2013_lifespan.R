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