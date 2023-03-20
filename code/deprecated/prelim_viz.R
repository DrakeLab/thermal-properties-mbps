## Title: Visualization of thermal traits ######################################
##
## Project: Global zoonoses - spillover of mosquito-borne pathogens
##
## Purpose: Produce visualizations of thermal trait distributions for a sanity check
##
## Contents: 1) Set-up, load in necessary packages and data-sets
##           2) Histograms of parameter distributions of thermal traits
##           3) Plots of thermal trait functions with 89% HCI
##           4)
##
##
## Inputs:  data - data/clean/gamma_fits.csv
##
##          code - code/Mordecai2017/mcmc_utils_all.R
##                 code/Mordecai2017/temp_functions_all.R
##
## Outputs: functions:
##          data - data/clean/ThermalTraitSamples.csv
##
## Written and maintained by: Kyle Dahlin, kydahlin@gmail.com
## Initialized February 2023

# 1) Set-up, load in necessary packages and data-sets ----
### Load Libraries ----
library(tidyverse)
library(reshape2)
library(cowplot)
library(latex2exp) # nec
library(viridis) # nec

### Load in data and functions ----
# Load samples of thermal trait parameters
samples.All <- read.csv("data/clean/ThermalTraitSamples.csv")

# Load functions from Mordecai et al., 2017
# This file contains tools for analysis and visualization.
source("code/Mordecai2017/mcmc_utils_all.R")

# This file contains the thermal response functions and their derivatives.
source("code/Mordecai2017/temp_functions_all.R")

# This file contains functions to compute model outputs like R0 and sensitivities
source("code/output-functions.R")

# 2) Histograms of parameter distributions of thermal traits ----
plot_df <- samples.All %>%
  dplyr::select(-func) %>% 
  melt(id = c("Species", "trait", "sample_num"))

###* Figure: thermal trait parameter posterior distributions ----
parm_hists <- plot_df %>%
  ggplot(aes(value, color = Species, fill = Species)) +
  geom_histogram(aes(), bins = 100) +
  # geom_density() +
  facet_wrap(trait ~ variable, scales = "free") +
  theme_minimal_grid(16)

# 3) Plots of thermal trait functions with 89% HCI ----

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
Linear <- function(q, z) {
  function(t) {
    pmax(-q * t + z, 0, na.RM = FALSE)
  }
}

# Function: designate proper thermal response function
# - output is a function of temperature
get.thermal.response <- function(data_in, Temperature) {
  parms <- dplyr::select(data_in, c, T0, Tm)
  function_type <- dplyr::select(data_in, func)
  
  temp_function <- case_when(
    function_type == "Briere" ~ Briere(parms$c, parms$T0, parms$Tm),
    function_type == "Quadratic" ~ Quadratic(parms$c, parms$T0, parms$Tm),
    function_type == "Linear" ~ Linear(parms$c, parms$T0)
  )
  
  out <- temp_function(Temperature)
}

### Temperature vector used for visualizations ----
Temps <- seq(10, 45, length.out = 200)

# Thinning intervals for samples
res_reduce <- function(df, new_length) {
  old_length <- length(unique(df))
  if (old_length < new_length) {
    warning("new_length is larger than old length. Vector will be unchanged")
    new_length <- old_length
  }
  ret <- unique(df)[seq(1, length(unique(df)), length.out = new_length)]
}

thin_size <- 100


# For each mosquito species, trait, and sample, get a thermal response curve
TPC_df <- samples.All %>%
  # filter(sample_num %in% seq(1, thin_size)) %>%
  full_join(list(Temperature = Temps), by = character(), copy = TRUE) %>%
  # group_by(sample_num) %>%
  # try sapply or lapply
  # mutate(Trait_val = get.thermal.response(.,Temperature))Linear(parms$c, parms$T0)
  mutate(Trait_val = case_when(
    func == "Briere" ~ Briere(c, T0, Tm)(Temperature),
    func == "Quadratic" ~ Quadratic(c, T0, Tm)(Temperature),
    func == "Linear" ~ Linear(c, T0)(Temperature)
  )) %>% 
  dplyr::select(-c("c", "T0", "Tm"))

# get mean TPC from samples
meanTPC_df <- TPC_df %>% 
  group_by(Species, trait, Temperature) %>% 
  summarise(mean_val = mean(Trait_val), .groups = "keep")


# get edges of 89% HCI of samples
quantsTPC_df <- TPC_df %>% 
  group_by(Species, trait, Temperature) %>% 
  mutate(lowHCI_val = quantile(Trait_val, 0.055)) %>% 
  mutate(highHCI_val = quantile(Trait_val, 0.945)) %>% 
  dplyr::select(-c("sample_num", "Trait_val", "func"))


###* Figure: TPC curves with 89% high confidence intervals ---- 
TPC_plot <- TPC_df %>%
  group_by(sample_num) %>% 
  arrange(Temperature) %>% 
  filter(Species == "Aedes albopictus") %>% 
  # group_by()
  ggplot(aes(x = Temperature)) +
  geom_line(data = meanTPC_df, aes(y = mean_val, color = Species)) +
  # 89% HCI of R0 TPC curves
  geom_ribbon(
    data = quantsTPC_df,
    aes(ymin = lowHCI_val, ymax = highHCI_val, fill = Species),
    alpha = 0.1
  ) +
  facet_wrap(~ trait, scales = "free", ncol = 2) +
  theme_minimal_grid(16)


# 4) Plots of R0 and HCIs ----

### Build full data frame ----
# System: mosquito species x pathogen
MosqPathPairs <- tibble(
  Mosquito_species = c(
    "Aedes aegypti", "Aedes albopictus",
    "Aedes aegypti",
    "Culex quinquefasciatus",
    "Anopheles spp."
  ),
  Pathogen = c("DENV", "DENV", "ZIKV", "WNV", "P.fal.")
)

# Temperature
min_Temp <- 10
max_Temp <- 40
res_Temp <- 0.0025
Temp_vec <- seq(min_Temp, max_Temp, by = res_Temp)

### Create host parameter data frame ----
# *Host recruitment rate
lambdaH <- .005
# *Host mortality rate
muH <- 1 / (365 * 20)
# *Host maximum biting tolerance / annoyance threshold (mosquitoes bites per day)
sigmaH_vec <- sort(c(10^seq(-0.25, 3.25, by = .25), 20, 50, Inf))
# *Host carrying capacity
KH_vec <- 10^seq(-2,5, length.out = 21)

## Get host-related pathogen parameters
# *Probability of becoming infected after bite from infectious mosquito
betaH <- 1
# *Host recovery rate
gammaH <- 1 / 5

Host_df <- expand_grid(
  # Host parameters
  lambdaH = lambdaH,
  muH = muH,
  KH = KH_vec,
  sigmaH = sigmaH_vec,
  # Pathogen-specific host parameters
  gammaH = gammaH,
  betaH = betaH
) %>%
  mutate(Model = case_when(
    is.finite(sigmaH) ~ "Chitnis",
    TRUE ~ "RM"
  ))

#### Create mosquito parameter data frame ----
# small constant to avoid INFs
eps <- .Machine$double.eps

# Get a smaller dataframe for TPCs, for visualization
TPC_df.reduced <- TPC_df %>% 
  filter(sample_num %in% seq(1, 200))

MosqParam_df <- TPC_df.reduced %>%
  dplyr::select(-"func") %>%
  pivot_wider(
    id_cols = c("Species", "sample_num", "Temperature"),
    names_from = trait, values_from = Trait_val
  ) %>%
  # Biting rate (sigmaV): GCR, GCD, a
  mutate(sigmaV = a) %>%
  # Fecundity (f): 0.5*EFD/a, 0.5*TFD*muV/a, 0.5*ER*pO/a
  mutate(fecundity = 0.5 * EFD / (a + eps)) %>%
  # Egg development probability (deltaL): pEA, EV*pLA, e2a
  mutate(deltaL = e2a) %>%
  # Larval development rate (rhoL): MDR
  mutate(rhoL = MDR) %>%
  # Adult mortality rate (muV): 1/lf, -log(p)
  mutate(muV = 1 / (lf + eps)) %>%
  # Extrinsic incubation rate (etaV): PDR
  mutate(etaV = PDR) %>%
  # Vector competence (betaV): bc
  mutate(betaV = b * c) %>%
  # Select the relevant model parameters
  dplyr::select(
    Species, sample_num, Temperature, sigmaV, fecundity, deltaL,
    rhoL, muV, lf, etaV, betaV
  )

KL_num <- 3000
MosqParam_df$KL <- KL_num

### Combine data frames and calculate R0 ----
R0_df <- full_join(Host_df, MosqParam_df, by = character())%>%
  mutate(V0 = compute.V0(.)) %>%
  mutate(R0 = compute.R0(.))

###* Plot sample R0 TPCs for a given species, host density and biting tolerance

simpleR0_df <-  R0_df %>%
  filter(sample_num %in% seq(1,20)) %>% 
  # To set num. of curves, change "length.out" to be the number of curves you want
  filter(KH %in% c(1, 10, 100)) %>%
  # Just consider one finite value of sigmaH and the Ross-Macdonald model
  filter(sigmaH == 100 | sigmaH == Inf) %>%
  filter(Species == "Aedes aegypti") %>% 
  filter(Temperature > 15 & Temperature < 35) %>% 
  group_by(sample_num) %>%
  # Normalize R0 for each curve_ID x system_ID combination, so that the maximum is always at one
  mutate(R0 = ifelse(is.nan(R0), 0, R0)) %>%
  mutate(norm_R0 = R0 / max(R0)) %>%
  mutate(norm_R0 = ifelse(is.nan(norm_R0), 0, norm_R0)) %>%
  # Restrict temperatures to where R0 is positive
  # filter(norm_R0 > 0) %>%
  ungroup() %>%
  # dplyr::select(Temperature, R0, norm_R0, sample_num) %>%
  arrange(Temperature, R0)

simpleR0plot <- simpleR0_df %>% 
  ggplot(aes(x = Temperature, y = R0, group = sample_num, color = sample_num)) +
  geom_path(alpha = 0.5, lwd = 1) +
  facet_wrap(sigmaH ~ KH, scales = "free",
             labeller = label_both)

### Example R0 TPCs with HCIs ----
newR0_df <-  R0_df %>%
  # filter(sample_num %in% seq(1,20)) %>% 
  # To set num. of curves, change "length.out" to be the number of curves you want
  filter(KH %in% c(1, 10, 100)) %>%
  # Just consider one finite value of sigmaH and the Ross-Macdonald model
  filter(sigmaH == 100 | sigmaH == Inf) %>%
  filter(Species == "Aedes aegypti") %>% 
  filter(Temperature > 15 & Temperature < 35) %>% 
  group_by(sample_num) %>%
  # Normalize R0 for each curve_ID x system_ID combination, so that the maximum is always at one
  mutate(R0 = ifelse(is.nan(R0), 0, R0)) %>%
  mutate(norm_R0 = R0 / max(R0)) %>%
  mutate(norm_R0 = ifelse(is.nan(norm_R0), 0, norm_R0)) %>%
  # Restrict temperatures to where R0 is positive
  # filter(norm_R0 > 0) %>%
  ungroup() %>%
  # dplyr::select(Temperature, R0, norm_R0, sample_num) %>%
  arrange(Temperature, R0)

meanR0TPC_df <- newR0_df %>%
  group_by(Species, Temperature, Model, sigmaH, KH) %>%
  summarise(
    mean_R0 = mean(R0),
    median_val = median(R0),
    # mode_val = mlv(norm_R0, method = 'mfv'),
    .groups = "keep"
  ) %>%
  arrange(Species, sigmaH, KH, Temperature, mean_R0, median_val)# , mode_val)

quantsR0TPC_df <- newR0_df %>%
  group_by(Species, Temperature, Model, sigmaH, KH) %>%
  mutate(lowHCI_val = quantile(R0, 0.055)) %>%
  mutate(highHCI_val = quantile(R0, 0.945)) %>%
  arrange(Species, sigmaH, KH, Temperature, lowHCI_val, highHCI_val) %>%
  dplyr::select(-c("sample_num"))

newR0plot <- simpleR0_df %>% 
  ggplot(aes(x = Temperature)) +
  geom_path(data = meanR0TPC_df, aes(y = mean_R0)) +
  geom_ribbon(
    data = filter(quantsR0TPC_df),
    aes(ymin = lowHCI_val, ymax = highHCI_val),
    alpha = 0.1
  ) +
  facet_wrap(sigmaH ~ KH, scales = "free",
             labeller = label_both) 

# 5) Plots of Topt ----
# Calculate Topt
Topt_df <- R0_df %>%
  group_by(Species, sample_num, sigmaH, KH) %>%
  filter(R0 == max(R0)) %>%
  # Get temperature at which R0 is maximized
  mutate(Topt = Temperature) %>%
  # Get the largest value of R0
  mutate(R0opt = R0) %>%
  # Keep track of whether R0 exceeds one
  mutate(threshold_bool = R0 > 1) %>%
  # Keep only relevant column labels
  dplyr::select(
    Model, Species, sample_num, sigmaH, KH, Topt, R0opt, threshold_bool # , CHmin, CHmax
  ) %>%
  # remove duplicate rows
  distinct()

# Calculate CTmin and CTmax
# Get the thermal range of parasite as a function of host traits
TempRange_df <- R0_df %>%
  # Restrict to rows where R0 exceeds one
  filter(R0 > 1) %>%
  group_by(Species, sample_num, sigmaH, KH) %>%
  # Get lowest temperature at which R0 exceeds one
  mutate(CTmin = min(Temperature)) %>%
  # Get highest temperature at which R0 exceeds one
  mutate(CTmax = max(Temperature)) %>%
  dplyr::select(Model, Species, sample_num, sigmaH, KH, CTmin, CTmax) %>%
  # remove duplicate rows
  distinct()

# Plot of Topt for each sample
Topt_plot <- Topt_df %>% 
  # filter(sample_num %in% seq(1,20)) %>% 
  # filter(KH %in% c(1,10,100)) %>% 
  filter(sigmaH %in% c(100, Inf)) %>% 
  filter(Species == "Aedes aegypti") %>% 
  ggplot(aes(x = log10(KH), y = Topt, color = sample_num, group = sample_num)) +
  geom_path() +
  facet_wrap(~sigmaH, labeller = label_both) +
  theme_minimal_grid(16)

# Plot Topt divergence for each sample
Topt_diverge_df <- Topt_df %>%
  ungroup() %>% 
  # filter(sample_num %in% seq(1,20)) %>% 
  # filter(KH %in% c(1,10,100)) %>% 
  filter(sigmaH %in% c(100, Inf)) %>% 
  filter(Species == "Aedes aegypti") %>% 
  select(-c(Model, Species, R0opt, threshold_bool)) %>% 
  group_by(sigmaH) %>% 
  pivot_wider(id_cols = c("sample_num", "KH"),
              names_from = c(sigmaH), names_prefix = "sigmaH=",
              values_from = "Topt") %>% 
  mutate(Topt_div = `sigmaH=100` - `sigmaH=Inf`) %>%
  select(-c(`sigmaH=100`, `sigmaH=Inf`)) %>% 
  arrange(KH, sample_num)

meanTopt_diverge_df <- Topt_diverge_df %>% 
  group_by(KH) %>%
  # mutate(mean_Topt_div = mean(Topt_div)) %>% 
  summarise(
    mean_val = mean(Topt_div),
    median_val = median(Topt_div),
    # mode_val = mlv(norm_R0, method = 'mfv'),
    .groups = "keep"
  )

Topt_diverge_plot <- Topt_diverge_df %>% 
  ggplot(aes(x = log10(KH))) +
  geom_path(aes(y = Topt_div, group = sample_num), alpha = 0.1) +
  geom_path(data = meanTopt_diverge_df,
            aes(y = mean_val), lwd = 1) +
  theme_minimal_grid(16)


# Plot of CTmin for each sample
CTmin_plot <- TempRange_df %>% 
  # filter(sample_num %in% seq(1,20)) %>% 
  # filter(KH %in% c(1,10,100)) %>% 
  filter(sigmaH %in% c(100, Inf)) %>% 
  filter(Species == "Aedes aegypti") %>% 
  ggplot(aes(x = log10(KH), y = CTmin, color = sample_num, group = sample_num)) +
  geom_path() +
  facet_wrap(~sigmaH, labeller = label_both) +
  theme_minimal_grid(16)

# Plot of CTmax for each sample
CTmax_plot <- TempRange_df %>% 
  # filter(sample_num %in% seq(1,20)) %>% 
  # filter(KH %in% c(1,10,100)) %>% 
  filter(sigmaH %in% c(100, Inf)) %>% 
  filter(Species == "Aedes aegypti") %>% 
  ggplot(aes(x = log10(KH), y = CTmax, color = sample_num, group = sample_num)) +
  geom_path() +
  facet_wrap(~sigmaH, labeller = label_both) +
  theme_minimal_grid(16)
