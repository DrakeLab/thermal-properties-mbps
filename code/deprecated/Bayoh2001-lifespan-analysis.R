## Title: Fitting lifespan data from Mordecai 2013 ###############################
##
## Project: Global zoonoses - spillover of mosquito-borne pathogens
##
## Purpose: Estimate mortality rate from data provided in Bayoh 2001, following 
##          process outlined in Mordecai et al., 2013
##
## Contents: 1) Set-up, load in necessary packages and data-sets
##           2) Load in and clean data as necessary
##           3) Filter data to values prescribed in Mordecai 2013 supplement
##           4) Fit exponential survival curves
##           5) Get estimates of mu and translate to lifespan
##           6) Output data in correct format for further processing (to data_cleaning.R)
##
##
## Inputs:  data/raw/Mordecai_2013/survival_data.csv
##
##
## Outputs: data - data/clean/Bayoh2001_mortality.csv
##
## Written and maintained by: Kyle Dahlin, kydahlin@gmail.com
## Initialized March 2023



# 1) Set-up load in necessary packages ----
require(tidyverse)
require(stringr)

# 2) Load in and clean data as necessary ----
###* Using data from (Bayoh 2001).
mortality.data <- read.csv("data/raw/Mordecai_2013/survival_data.csv", header = TRUE) %>%
  # change columns to be: day, temperature (T), and proportion alive (prop.alive)
  pivot_longer(cols = X5.C:X40.C, names_to = "T", values_to = "prop.alive") %>% 
  rename(day = D) %>% 
  mutate(T = readr::parse_number(T)) %>% 
  relocate(T, .before = day) %>%
  arrange(T, day)


# 3) Filter data to values prescribed in Mordecai 2013 supplement ----

###* Find the final day where the population proportion exceeded 1%.
# Get the first day at which zero living mosquitoes were reported
first_zero <- mortality.data %>%
  group_by(T) %>%
  filter(prop.alive == 0) %>%
  summarize(first.zero = min(day)) %>%
  # population was persisted through the entire experiment at 20C
  add_row(T = 20, first.zero = Inf)

# Get the last day the proportion alive exceeded 0.1%
cutoff_days <- mortality.data %>%
  # Ignore days after the first zero was reported
  right_join(first_zero) %>%
  filter(day < first.zero) %>%
  # for each temperature, find the last day at which prop.alive > 0.001
  group_by(T) %>%
  filter(prop.alive > 0.001) %>%
  summarize(last.exceed.001 = max(day))

# Get the first day experimental data was reported
first_days <- mortality.data %>%
  # for each temperature, find the first day after day zero at which prop.alive
  # is recorded
  group_by(T) %>%
  filter(day > 0 & !is.na(prop.alive)) %>%
  summarize(first.day = min(day))

###* They then filtered the data to only include the following days:
###* - the first day of the experiment,
###* - one day preceding reaching the 1% threshold,
###* - the day where the 1% threshold is met,
###* - then the three days following.
reduced.mortality.data <- mortality.data %>%
  right_join(cutoff_days) %>%
  right_join(first_days) %>%
  rowwise() %>%
  filter(day %in% c(0, first.day, last.exceed.001 + seq(-1, 3))) %>%
  select(-c(last.exceed.001, first.day))

# # Visualize reduced mortality data
# test_plot <- reduced.mortality.data %>%
#   arrange(day) %>%
#   ggplot(aes(x = day, y = prop.alive, group = T, color = as.factor(T))) +
#   geom_path(lwd = 1) +
#   geom_point(data = mortality.data)

# 4) Fit exponential survival curves and get coefficient (mu) ----
###* The parameter of the exponential survival function was used for mu, an
###* independent observation of mortality at each temperature.
###* This captured temperature-associated trends observed in the data, but the
###* shape of these fitted curves were substantially different from those
###* derived from a (more complex) Gompertz curve fit, especially between days
###* 20 and 40.

# Define exponential function
Exp_func <- function(x, mu){
  result <- exp(-mu * x)
  return(result)
}

# Fit reduced data to negative exponential
fit.df <- tibble(T = numeric(), mu = numeric(), lf = numeric())
for (TT in unique(reduced.mortality.data$T)) {
  df <- filter(reduced.mortality.data, T == TT)
  # Fit Gompertz curves to the data
  Exp_fit <- nls(prop.alive ~ Exp_func(day, mu),
                  data = df,
                  start = list(mu = 0.01))
  mu <- coef(Exp_fit)
  
  fit.df <- add_row(fit.df, T = TT, mu = mu, lf = 1/mu)
}

# # Visualize fit
# Exp_df <- tibble(day = unique(reduced.mortality.data$day)) %>% 
#   full_join(fit.df, by = character()) %>% 
#   mutate(pred.p = Exp_func(x = day, mu = mu)) %>% 
#   arrange(day)
# 
# survival_plot <- reduced.mortality.data %>% 
#   arrange(day) %>% 
#   # filter(prop.alive != 0) %>%
#   ggplot(aes(x = day, y = prop.alive, color = as.factor(T), group = T)) +
#   geom_point() +
#   geom_path(data = Exp_df, aes(x = day, y = pred.p), lwd = 1)

# 6) Output data in correct format for further processing (to data_cleaning.R) ----
data.Mordecai2013.mortality <- fit.df %>% 
  select(-mu) %>% 
  pivot_longer(cols = lf, names_to = "trait.name", values_to = "trait") %>% 
  mutate(ref = "Bayoh2001", trait2 = NA, trait2.name = NA, 
         mosquito_species = "Anopheles gambiae", pathogen = NA)

write_csv(data.Mordecai2013.mortality,"data/clean/Bayoh2001_mortality.csv")





