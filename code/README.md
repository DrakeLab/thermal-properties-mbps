# Code for performing analyses
host-traits.R
  * Sets up vertebrate host parameters

vector-traits.R
  * Sets up vector trait parameters, these are functions of temperature and are sampled from posterior distributions constructed in get-thermal-trait-priors.R

  - data-cleaning.R
    * Processes data from data/raw/ folder into a format used to derive mosquito trait TPC posterior distributions in get-thermal-trait-priors.R

  - get-thermal-trait-priors.R
    * Samples from data-informed mosquito trait TPC hyperparameters using models in code/jags-models folder and functions from code/Mordecai_2017

  - trait-transform.R
    * Transform traits into parameters to be input into the model

get-outputs.R
  * Calculate the basic reproduction number (R0), transmission thermal optimum (Topt), critical thermal minimum (CTmin), critical thermal maximum (CTmax), and width of the parasite population thermal niche (CTwidth)

sensitivity.R
  * Sensitivity and uncertainty analysis of R0, Topt, CTmin, CTmax, and CTwidth
