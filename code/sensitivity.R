## Title: Sensitivity and uncertainty analysis of Topt, CTmin, CTmax, and CTwidth
##
## Project: Global zoonoses - spillover of mosquito-borne pathogens
##
## Purpose: Determine which traits most underlie the variation in key outputs
##
## Contents: 0) Set up, load in necessary packages and data-sets
##           1) Define accessory functions
##           2) Set up data frames
##           3) Topt local, global sensitivity, uncertainty
##           4) CTmin local, global sensitivity, uncertainty
##           5) CTmax local, global sensitivity, uncertainty
##           6) CTwidth local, global sensitivity, uncertainty
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
library(multidplyr)
library(foreach)
library(doParallel)
library(progress)
library(latex2exp) # nec
library(viridis) # nec
library(cowplot) # nec
library(MetBrewer) #nec
library(svglite) #nec
library(HDInterval)

# Set up parallel
if (!exists("cluster")) {
  cluster_size <- min(21, parallel::detectCores() - 2)
  cluster <- makeCluster(cluster_size) #new_cluster(cluster_size)
  registerDoParallel(cluster)
  # cluster_library(cluster, c("dplyr", "tidyr"))
}

eps <- .Machine$double.eps

# Initialize a data frame to save our analyses
init.df <- tibble(system_ID = c(), Temperature = c(), Model = c(),
                  sigmaH = c(), KH = c(), variable = c(),
                  lowHCI = c(), highHCI = c(), mean = c(), median = c())

# Load needed datasets
# R0 data
data.R0 <- read_rds("results/R0_vals.rds") # !!! needs to be updated to include 1e-1, 1e0, 1e1
# Topt data
data.Topt <- read_rds("results/Topt_vals.rds")
# Critical thermal extrema data
data.CT <- read_rds("results/CT_vals.rds")


# KH_vec_length <- 300 # full = 300, thin = 20
# 
# # Biting tolerance vector: Number of values to consider for biting tolerance
# sigmaH_vec_length <- 300 # full = 300, thin = 20
# 
# source("code/trait-variation.R")

# 1) Define accessory functions -------------------------------------------

# Function: Slice data to optimize usage of parallel processing and memory
slice<-function(x,n) {
  N<-length(x);
  lapply(seq(1,N,n),function(i) x[i:min(i+n-1,N)])
}

# Label functions: 
appender_sigmaH <- function(string) {
  unname(TeX(paste("$\\sigma_H = $", string)))}
appender_KH <- function(string) {
  unname(TeX(paste("$K_H = $", string)))}


# 2) Set up data frames ---------------------------------------------------

###* Mean traits data frame ----

# Get mean values of uncertain traits

data.Vec <- read_rds("results/VecTPC_vals.rds")

mean.Vec <- data.Vec %>% 
  select(-c(KL, mosquito_species, pathogen, muL, etaL)) %>% 
  pivot_longer(cols = lf:V0, names_to = "variable", values_to = "value") %>%
  group_by(system_ID, Temperature, variable) %>% 
  summarise(mean = mean(value)) %>% 
  pivot_wider(names_from = "variable", values_from = "mean") %>% 
  mutate(etaL = rhoL / (deltaL + eps)) %>%
  # Aquatic-stage mosquito mortality rate. This should be ~infinite if deltaL = 0
  mutate(muL = etaL - rhoL)

# # Diagnostic: plot means
# plot.mean <- mean.Vec %>%
#   pivot_longer(cols = V0:sigmaV_f, names_to = "variable", values_to = "value") %>%
#   ggplot(mapping = aes(x = Temperature, y = value, color = system_ID)) +
#   geom_path() +
#   facet_wrap(~variable, scales = "free")
# plot.mean


# 3) R0 uncertainty -------------------------------------------------------

# Calculate the width of the 95% HPD for the full posterior of R0 across temperatures

sigmaH_vec <- c(50, Inf)#unique(data.Host$sigmaH)
KH_vec <- 10^seq(-2,5)

data.R0HPD <- dplyr::filter(data.Host,
                            sigmaH %in% sigmaH_vec,
                            KH %in% KH_vec)

full.R0.HPD <- expand_grid(data.R0HPD, 
                           data.Vec) %>% 
  data.table::data.table() %>%
  mutate(RV = ifelse(is.infinite(sigmaH),
                     sigmaV * betaH / (1 / (lf + eps)), # Ross-Macdonald
                     sigmaH * sigmaV * betaH * KH / ((1 / (lf + eps)) * (sigmaH * KH + sigmaV * V0)))) %>%
  mutate(bV = ifelse(is.infinite(sigmaH),
                     sigmaV, # Ross-Macdonald model
                     sigmaV * sigmaH * KH / (sigmaH * KH + sigmaV * V0 + eps))) %>%
  mutate(RH = ifelse(V0 == 0,
                     0,
                     bV * betaV * V0 * exp(-1 / (lf * etaV)) / (KH * (gammaH + muH) + eps))) %>%
  # Basic reproduction number
  mutate(R0 = sqrt(RV*RH)) %>%
  select(system_ID, sample_num, sigmaH, KH, Temperature, R0) %>% 
  group_by(system_ID, sigmaH, KH, Temperature) %>%
  summarise(
    HPD_low = hdi(R0, credMass = 0.95)[1],
    HPD_high = hdi(R0, credMass = 0.95)[2],
    HPD_width = max(eps,HPD_high-HPD_low),
    .groups = "keep"
  ) %>% 
  group_by(system_ID, sigmaH, KH) %>%
  mutate(norm_HPD_width = HPD_width / max(HPD_width))

# Save R0 highest posterior density data
# write_rds(full.R0.HPD, "results/full_R0_HPD.rds")
full.R0.HPD <- read_rds("results/full_R0_HPD.rds")

# R0 HPD diagnostic plot
R0.HPD.plot <- full.R0.HPD %>%
  # dplyr::filter(KH %in% c(1,100)) %>%
  arrange(system_ID, sigmaH) %>% 
  # dplyr::filter(is.finite(sigmaH)) %>%
  ggplot(aes(x = Temperature, y = norm_HPD_width, color = as.factor(KH))) +
  geom_path(lwd = 1) +
  scale_color_viridis_d(
    name = "Vertebrate host population\ndensity (ind/ha)",
    breaks = 10^seq(-2, 4),
    labels = unname(c(0.01, 0.1, 1, 10, 100, TeX("$10^3$"), TeX("$10^4$"))),
    limits = 10^seq(-2, 4)
  ) +
  scale_y_continuous(name = "Normalized HPD width") +
  facet_grid(rows = vars(system_ID), cols =vars(sigmaH),
             labeller = labeller(sigmaH = as_labeller(appender_sigmaH,
                                                      default = label_parsed),
                                 KH = as_labeller(appender_KH, 
                                                  default = label_parsed)),
             scales = "free") +
  theme_minimal(16)
R0.HPD.plot

# Save full R0 highest posterior density data
ggsave("figures/results/R0_HPD.svg", 
       plot = R0.HPD.plot,
       width = 16, height = 9)


# Get focal parameter names
temp_vars <- data.Vec %>% 
  select(-c(system_ID, mosquito_species, pathogen, Temperature, sample_num,
            etaL, muL, KL, V0, lf)) %>% 
  colnames() %>% c("muV")

# Initialize data frame
R0.HPD <- tibble(system_ID = c(), Temperature = c(), focal_var = c(),
                 sigmaH = c(), KH = c(), 
                 HPD_width = c(), full_HPD_width = c(), rel_HPD_width = c())


# # For each focal parameter, calculate relative HPD width
# R0.HPD <- foreach(var_name = temp_vars,
#                   .packages = c("tidyverse", "HDInterval"),
#                   combine = rbind) %dopar% {
#   # a) Set all but focal parameter to its posterior mean
#   data.HPD.Vec <- data.Vec %>% 
#     mutate(muV = 1/lf) %>% 
#     select(system_ID, Temperature, sample_num, all_of(var_name), KL) %>% 
#     full_join(mean.Vec %>% 
#                 mutate(muV = 1/lf) %>% 
#                 select(-all_of(var_name))%>% 
#                 distinct()) %>% 
#     mutate(lf = 1/muV) %>%
#     mutate(V0 = ifelse(sigmaV_f * deltaL < (1 / lf),
#                        0,
#                        KL * rhoL * lf * (1 - 1 / (lf * sigmaV_f * deltaL)))) %>%
#     select(-lf)
#   # b) Get posterior samples of R0 (as a function of temperature)
#   R0.HPD <- expand_grid(data.R0HPD, data.HPD.Vec) %>%
#     mutate(lf = 1/muV) %>% 
#     data.table::data.table() %>%
#     mutate(RV = ifelse(is.infinite(sigmaH),
#                        sigmaV * betaH / (1 / (lf + eps)), # Ross-Macdonald
#                        sigmaH * sigmaV * betaH * KH / ((1 / (lf + eps)) * (sigmaH * KH + sigmaV * V0)))) %>%
#     mutate(bV = ifelse(is.infinite(sigmaH),
#                        sigmaV, # Ross-Macdonald model
#                        sigmaV * sigmaH * KH / (sigmaH * KH + sigmaV * V0 + eps))) %>%
#     mutate(RH = ifelse(V0 == 0,
#                        0,
#                        bV * betaV * V0 * exp(-1 / (lf * etaV)) / (KH * (gammaH + muH) + eps))) %>%
#     # Basic reproduction number
#     mutate(R0 = sqrt(RV*RH)) %>%
#     select(system_ID, sample_num, sigmaH, KH, Temperature, R0) %>% 
#     # c) Calculate the 95% HPD at each temperature
#     group_by(system_ID, sigmaH, KH, Temperature) %>% 
#     summarise(
#       HPD_low = hdi(R0, credMass = 0.95)[1],
#       HPD_high = hdi(R0, credMass = 0.95)[2],
#       HPD_width = max(eps, HPD_high-HPD_low),
#       .groups = "keep"
#     ) %>% 
#     select(system_ID, sigmaH, KH, Temperature, HPD_width) %>% 
#     right_join(full.R0.HPD %>% 
#                  select(-c(HPD_low, HPD_high)) %>% 
#                  rename(full_HPD_width = HPD_width)) %>% 
#     mutate(focal_var = var_name) %>% 
#     group_by(system_ID, sigmaH, KH, focal_var, Temperature) %>% 
#     # d) Normalize HPD width by the full HPD width when all parameters are allowed to vary
#     mutate(rel_HPD_width = ifelse(full_HPD_width %in% c(0,eps), 0, HPD_width / full_HPD_width))
#     
#     R0.HPD
# }

# For each focal parameter, calculate relative HPD width
for (var_name in temp_vars) {
  # a) Set all but focal parameter to its posterior mean
  data.HPD.Vec <- data.Vec %>% 
    mutate(muV = 1/lf) %>% 
    select(system_ID, Temperature, sample_num, all_of(var_name), KL) %>% 
    full_join(mean.Vec %>% 
                mutate(muV = 1/lf) %>% 
                select(-all_of(var_name))%>% 
                distinct()) %>% 
    mutate(lf = 1/muV) %>%
    mutate(V0 = ifelse(sigmaV_f * deltaL < (1 / lf),
                       0,
                       KL * rhoL * lf * (1 - 1 / (lf * sigmaV_f * deltaL)))) %>%
    select(-lf)
  # b) Get posterior samples of R0 (as a function of temperature)
  R0.HPD <- expand_grid(data.R0HPD, data.HPD.Vec) %>%
    mutate(lf = 1/muV) %>% 
    data.table::data.table() %>%
    mutate(RV = ifelse(is.infinite(sigmaH),
                       sigmaV * betaH / (1 / (lf + eps)), # Ross-Macdonald
                       sigmaH * sigmaV * betaH * KH / ((1 / (lf + eps)) * (sigmaH * KH + sigmaV * V0)))) %>%
    mutate(bV = ifelse(is.infinite(sigmaH),
                       sigmaV, # Ross-Macdonald model
                       sigmaV * sigmaH * KH / (sigmaH * KH + sigmaV * V0 + eps))) %>%
    mutate(RH = ifelse(V0 == 0,
                       0,
                       bV * betaV * V0 * exp(-1 / (lf * etaV)) / (KH * (gammaH + muH) + eps))) %>%
    # Basic reproduction number
    mutate(R0 = sqrt(RV*RH)) %>%
    select(system_ID, sample_num, sigmaH, KH, Temperature, R0) %>% 
    # c) Calculate the 95% HPD at each temperature
    group_by(system_ID, sigmaH, KH, Temperature) %>% 
    summarise(
      HPD_low = hdi(R0, credMass = 0.95)[1],
      HPD_high = hdi(R0, credMass = 0.95)[2],
      HPD_width = max(eps, HPD_high-HPD_low),
      .groups = "keep"
    ) %>% 
    select(system_ID, sigmaH, KH, Temperature, HPD_width) %>% 
    right_join(full.R0.HPD %>% 
                 select(-c(HPD_low, HPD_high)) %>% 
                 rename(full_HPD_width = HPD_width)) %>% 
    mutate(focal_var = var_name) %>% 
    group_by(system_ID, sigmaH, KH, focal_var, Temperature) %>% 
    # d) Normalize HPD width by the full HPD width when all parameters are allowed to vary
    mutate(rel_HPD_width = ifelse(full_HPD_width %in% c(0,eps), 0, HPD_width / full_HPD_width)) %>% 
    rbind(R0.HPD)
}

# Save R0 relative highest posterior density data
# write_rds(R0.HPD, "results/R0_HPD_sens.rds")
R0.HPD <- read_rds("results/R0_HPD_sens.rds")

## Plot R0 uncertainty 
Temp_range_for_plot <- R0.HPD %>% 
  ungroup() %>%
  dplyr::filter(rel_HPD_width > eps) %>%
  # dplyr::filter(KH == KH_select) %>%
  select(Temperature) %>% range()

var_name_table <- list(
  betaV = "Vector competence",
  deltaL = "Pr(larval survival)",
  etaV = "Parasite development rate",
  muV = "Mortality rate",
  rhoL = "Immature development rate",
  sigmaV = "Max. biting rate",
  sigmaV_f = "Eggs per female per day"
)

# !!! replace this with "relative contributions to sensitivity" plot
R0.uncertainty.plot <- R0.HPD %>% 
  ungroup() %>% 
  # dplyr::filter(focal_var %in% c("sigmaV_f")) %>%
  mutate(system_ID = case_when(
    system_ID == "Anopheles gambiae / Plasmodium falciparum" ~ "An. gamb. / P. falciparum",
    system_ID == "Aedes aegypti / DENV" ~ "Ae. aegypti / DENV",
    system_ID == "Aedes albopictus / DENV" ~ "Ae. albopictus / DENV",
    system_ID == "Aedes aegypti / ZIKV" ~ "Ae. aegypti / ZIKV",
    system_ID == "Culex quinquefasciatus / WNV" ~ "Cx. quin. / WNV"
  )) %>% 
  arrange(system_ID, sigmaH) %>% 
  dplyr::filter(KH %in% c(1, 100)) %>%
  dplyr::filter(rel_HPD_width < 1.01) %>%
  # arrange(system_ID, sigmaH, focal_var, Temperature) %>%
  ggplot(aes(x = Temperature, y = rel_HPD_width, color = focal_var, fill = focal_var)) +
  geom_path(linewidth = 1) +
  scale_color_discrete(
    name = "Focal parameter",
    breaks = c("betaV", "deltaL", "etaV", "muV", "rhoL", "sigmaV", "sigmaV_f"),
    labels =  c("Vector competence", "Pr(larval survival)", 
                "Parasite development rate", "Mortality rate", 
                "Immature development rate", "Max. biting rate", 
                "Eggs per female per day")
  ) +
  scale_x_continuous(
    limits = Temp_range_for_plot
  ) +
  scale_y_continuous(
    name = "Relative HPD width",
    breaks = seq(0, 1, by = 0.2)
  ) +
  facet_grid(rows = vars(system_ID), cols = vars(KH, sigmaH), 
             labeller = labeller(sigmaH = as_labeller(appender_sigmaH, default = label_parsed),
                                 KH = as_labeller(appender_KH, default = label_parsed)),
             scales = "free") +
  ggtitle("Uncertainty analysis of R0") +
  theme_minimal_grid(12)
R0.uncertainty.plot

ggsave("figures/results/R0_uncertainty.svg", 
       plot = R0.uncertainty.plot,
       width = 16, height = 9)

# "Relative sensitivity" plot
plot.R0.rel.sens <- R0.HPD %>% 
  ungroup() %>% 
  mutate(system_label = case_when(
    system_ID == "Aedes aegypti / DENV" ~ "(a) *Ae. aegypti* and DENV",
    system_ID == "Aedes aegypti / ZIKV" ~ "(b) *Ae. aegypti* and ZIKV",
    system_ID == "Aedes albopictus / DENV" ~ "(c\u200b) *Ae. albopictus* and DENV",
    system_ID == "Anopheles gambiae / Plasmodium falciparum" ~ "(d) *An. gambiae* and *P. falciparum*",
    system_ID == "Culex quinquefasciatus / WNV" ~ "(e) *Cx. quinquefasciatus* and WNV"
  )) %>%
  filter(KH %in% c(1,100) & sigmaH %in% c(50, Inf)) %>% 
  filter(!(KH == 100 & sigmaH == Inf)) %>% 
  filter(rel_HPD_width != 0) %>% 
  # For each KH value, sum up all sensitivity values to get total sensitivity across parameters
  group_by(system_ID, sigmaH, KH, Temperature) %>% 
  mutate(sens_total = sum(rel_HPD_width)) %>% 
  # Divide each sensitivity value by the total sensitivity
  mutate(sens_rel = rel_HPD_width / sens_total) %>% 
  arrange(sens_rel) %>% 
  # Make a stacked plot, colored by parameter
  ggplot(aes(x = Temperature, y = sens_rel, fill = focal_var)) +
  geom_area() +
  scale_fill_discrete(
    name = "Focal parameter",
    breaks = c("betaV", "deltaL", "etaV", "muV", "rhoL", "sigmaV", "sigmaV_f"),
    labels =  c("Vector competence", "Pr(larval survival)", 
                "Parasite development rate", "Mortality rate", 
                "Immature development rate", "Max. biting rate", 
                "Eggs per female per day")
  ) +
  scale_x_continuous(
    limits = Temp_range_for_plot,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = unname(TeX("Relative contribution to $R_0$ HPD width"))
  ) +
  facet_grid(cols = vars(system_label),
             rows = vars(KH, sigmaH), 
             labeller = labeller(sigmaH = as_labeller(appender_sigmaH, default = label_parsed),
                                 KH = as_labeller(appender_KH, default = label_parsed)),
             scales = "free") +
  theme_minimal_hgrid(16) +
  theme(
    strip.text.x = element_text(size = 10)
  )
plot.R0.rel.sens


# 4) ddT R0 uncertainty -------------------------------------------------------

# Calculate the width of the 95% HPD for the full posterior of ddT R0 across temperatures

sigmaH_vec <- c(50, Inf)#unique(data.Host$sigmaH)
KH_vec <- 10^seq(-2,5)

data.ddTR0HPD <- dplyr::filter(data.Host,
                               sigmaH %in% sigmaH_vec,
                               KH %in% KH_vec)

full.ddTR0.HPD <- tibble(system_ID = c(), Temperature = c(), Model = c(),
                         sigmaH = c(), KH = c(), 
                         HPD_low = c(), HPD_high = c(), HPD_width = c())

full.ddTR0.HPD <- expand_grid(data.ddTR0HPD, 
                              data.Vec) %>% 
  data.table::data.table() %>%
  mutate(RV = ifelse(is.infinite(sigmaH),
                     sigmaV * betaH / (1 / (lf + eps)), # Ross-Macdonald
                     sigmaH * sigmaV * betaH * KH / ((1 / (lf + eps)) * (sigmaH * KH + sigmaV * V0)))) %>%
  mutate(bV = ifelse(is.infinite(sigmaH),
                     sigmaV, # Ross-Macdonald model
                     sigmaV * sigmaH * KH / (sigmaH * KH + sigmaV * V0 + eps))) %>%
  mutate(RH = ifelse(V0 == 0,
                     0,
                     bV * betaV * V0 * exp(-1 / (lf * etaV)) / (KH * (gammaH + muH) + eps))) %>%
  # Basic reproduction number
  mutate(R0 = sqrt(RV*RH)) %>%
  arrange(system_ID, sample_num, sigmaH, KH, Temperature) %>% 
  mutate(ddTR0 = (R0 - lag(R0)) / (Temperature - lag(Temperature))) %>% 
  select(system_ID, sample_num, sigmaH, KH, Temperature, ddTR0) %>% 
  group_by(system_ID, sigmaH, KH, Temperature) %>% 
  summarise(
    HPD_low = hdi(ddTR0, credMass = 0.95)[1],
    HPD_high = hdi(ddTR0, credMass = 0.95)[2],
    HPD_width = max(eps,HPD_high-HPD_low),
    .groups = "keep"
  )

# Save R0 highest posterior density data
# write_rds(full.ddTR0.HPD, "results/full_ddTR0_HPD.rds")
full.R0.HPD <- read_rds("results/full_R0_HPD.rds")

# # Diagnostic plot
# test.plot <- full.ddTR0.HPD %>%
#   dplyr::filter(KH == 1) %>%
#   ggplot(aes(x = Temperature, y = HPD_width)) +
#   geom_path() +
#   facet_grid(rows = vars(system_ID), cols = vars(sigmaH), scales = "free")
# test.plot

# Get focal parameter names
temp_vars <- data.Vec %>% 
  select(-c(system_ID, mosquito_species, pathogen, Temperature, sample_num,
            etaL, muL, KL, V0, lf)) %>% 
  colnames() %>% c("muV")

# Initialize data frame
ddTR0.HPD <- tibble(system_ID = c(), Temperature = c(), focal_var = c(),
                    sigmaH = c(), KH = c(), 
                    HPD_width = c(), full_HPD_width = c(), rel_HPD_width = c())

# For each focal parameter, calculate relative HPD width
for (var_name in temp_vars) {
  # a) Set all but focal parameter to its posterior mean
  data.HPD.Vec <- data.Vec %>% 
    mutate(muV = 1/lf) %>% 
    select(system_ID, Temperature, sample_num, all_of(var_name), KL) %>% 
    full_join(mean.Vec %>% 
                mutate(muV = 1/lf) %>% 
                select(-all_of(var_name))%>% 
                distinct()) %>% 
    mutate(lf = 1/muV) %>%
    mutate(V0 = ifelse(sigmaV_f * deltaL < (1 / lf),
                       0,
                       KL * rhoL * lf * (1 - 1 / (lf * sigmaV_f * deltaL)))) %>%
    select(-lf)
  # b) Get posterior samples of R0 (as a function of temperature)
  ddTR0.HPD <- expand_grid(data.R0HPD, data.HPD.Vec) %>%
    mutate(lf = 1/muV) %>% 
    data.table::data.table() %>%
    mutate(RV = ifelse(is.infinite(sigmaH),
                       sigmaV * betaH / (1 / (lf + eps)), # Ross-Macdonald
                       sigmaH * sigmaV * betaH * KH / ((1 / (lf + eps)) * (sigmaH * KH + sigmaV * V0)))) %>%
    mutate(bV = ifelse(is.infinite(sigmaH),
                       sigmaV, # Ross-Macdonald model
                       sigmaV * sigmaH * KH / (sigmaH * KH + sigmaV * V0 + eps))) %>%
    mutate(RH = ifelse(V0 == 0,
                       0,
                       bV * betaV * V0 * exp(-1 / (lf * etaV)) / (KH * (gammaH + muH) + eps))) %>%
    # Basic reproduction number
    mutate(R0 = sqrt(RV*RH)) %>%
    arrange(system_ID, sample_num, sigmaH, KH, Temperature) %>% 
    mutate(ddTR0 = (R0 - lag(R0)) / (Temperature - lag(Temperature))) %>% 
    select(system_ID, sample_num, sigmaH, KH, Temperature, ddTR0) %>% 
    # c) Calculate the 95% HPD at each temperature
    group_by(system_ID, sigmaH, KH, Temperature) %>% 
    summarise(
      HPD_low = hdi(ddTR0, credMass = 0.95)[1],
      HPD_high = hdi(ddTR0, credMass = 0.95)[2],
      HPD_width = max(eps, HPD_high-HPD_low),
      .groups = "keep"
    ) %>% 
    select(system_ID, sigmaH, KH, Temperature, HPD_width) %>% 
    right_join(full.ddTR0.HPD %>% 
                 select(-c(HPD_low, HPD_high)) %>% 
                 rename(full_HPD_width = HPD_width)) %>% 
    mutate(focal_var = var_name) %>% 
    group_by(system_ID, sigmaH, KH, focal_var, Temperature) %>% 
    # d) Normalize HPD width by the full HPD width when all parameters are allowed to vary
    mutate(rel_HPD_width = ifelse(full_HPD_width %in% c(0,eps), 0, HPD_width / full_HPD_width)) %>% 
    rbind(ddTR0.HPD)
}

# Save R0 relative highest posterior density data
# write_rds(ddTR0.HPD, "results/ddTR0_HPD_sens.rds")
ddTR0.HPD <- read_rds("results/ddTR0_HPD_sens.rds")

## Plot R0 uncertainty 
Temp_range_for_plot <- ddTR0.HPD %>% ungroup() %>% 
  dplyr::filter(HPD_width > eps) %>% 
  dplyr::filter(KH == KH_select) %>% 
  select(Temperature) %>%  range()

var_name_table <- list(
  betaV = c(TeX("$\\beta_V$")),
  deltaL = c(TeX("$\\delta_L$")),
  etaV = c(TeX("$\\eta_V$")),
  muV = c(TeX("$\\mu_V$")),
  rhoL = c(TeX("$\\rho_L$")),
  sigmaV = c(TeX("$\\sigma_V$")),
  sigmaV_f = c(TeX("$\\sigma_v f$"))
)

# !!! [x] create relative sensitivity plot, like for R0 HPDs
ddTR0.uncertainty.plot <- ddTR0.HPD %>% 
  ungroup() %>% 
  # dplyr::filter(focal_var %in% c("sigmaV_f")) %>%
  mutate(system_ID = case_when(
    system_ID == "Anopheles gambiae / Plasmodium falciparum" ~ "An. gamb. / P. falciparum",
    system_ID == "Aedes aegypti / DENV" ~ "Ae. aegypti / DENV",
    system_ID == "Aedes albopictus / DENV" ~ "Ae. albopictus / DENV",
    system_ID == "Aedes aegypti / ZIKV" ~ "Ae. aegypti / ZIKV",
    system_ID == "Culex quinquefasciatus / WNV" ~ "Cx. quin. / WNV"
  )) %>% 
  arrange(system_ID, sigmaH) %>% 
  dplyr::filter(KH %in% c(1, 100)) %>%
  dplyr::filter(rel_HPD_width < 1.05) %>%
  # arrange(system_ID, sigmaH, focal_var, Temperature) %>%
  ggplot(aes(x = Temperature, y = rel_HPD_width, color = focal_var)) +
  geom_path(linewidth = 1) +
  scale_color_discrete(
    name = "Focal parameter",
    breaks = c("betaV", "deltaL", "etaV", "muV", "rhoL", "sigmaV", "sigmaV_f"),
    labels =  unname(TeX(c("$\\beta_V$", "$\\delta_L$", "$\\eta_V$", 
                           "$\\mu_V$", "$\\rho_L$", "$\\sigma_V$",
                           "$\\sigma_v f$")))
  ) +
  scale_x_continuous(
    limits = Temp_range_for_plot
  ) +
  scale_y_continuous(
    name = "Relative HPD width",
    breaks = seq(0, 1, by = 0.2)
  ) +
  facet_grid(rows = vars(system_ID), cols = vars(KH, sigmaH), 
             labeller = labeller(sigmaH = as_labeller(appender_sigmaH, default = label_parsed),
                                 KH = as_labeller(appender_KH, default = label_parsed)),
             scales = "free") +
  ggtitle("Uncertainty analysis of R0") +
  theme_minimal_grid(12)
ddTR0.uncertainty.plot

ggsave("figures/results/ddTR0_uncertainty.svg", 
       plot = ddTR0.uncertainty.plot,
       width = 16, height = 9)

# "Relative sensitivity" plot
plot.ddTR0.rel.sens <- ddTR0.HPD %>% 
  ungroup() %>% 
  mutate(system_label = case_when(
    system_ID == "Aedes aegypti / DENV" ~ "(a) *Ae. aegypti* and DENV",
    system_ID == "Aedes aegypti / ZIKV" ~ "(b) *Ae. aegypti* and ZIKV",
    system_ID == "Aedes albopictus / DENV" ~ "(c\u200b) *Ae. albopictus* and DENV",
    system_ID == "Anopheles gambiae / Plasmodium falciparum" ~ "(d) *An. gambiae* and *P. falciparum*",
    system_ID == "Culex quinquefasciatus / WNV" ~ "(e) *Cx. quinquefasciatus* and WNV"
  )) %>%
  filter(KH %in% c(1,100) & sigmaH %in% c(50, Inf)) %>% 
  filter(!(KH == 100 & sigmaH == Inf)) %>% 
  filter(rel_HPD_width != 0) %>% 
  # For each KH value, sum up all sensitivity values to get total sensitivity across parameters
  group_by(system_ID, sigmaH, KH, Temperature) %>% 
  mutate(sens_total = sum(rel_HPD_width)) %>% 
  # Divide each sensitivity value by the total sensitivity
  mutate(sens_rel = rel_HPD_width / sens_total) %>% 
  arrange(sens_rel) %>% 
  # Make a stacked plot, colored by parameter
  ggplot(aes(x = Temperature, y = sens_rel, fill = focal_var)) +
  geom_area() +
  scale_fill_discrete(
    name = "Focal parameter",
    breaks = c("betaV", "deltaL", "etaV", "muV", "rhoL", "sigmaV", "sigmaV_f"),
    labels =  c("Vector competence", "Pr(larval survival)", 
                "Parasite development rate", "Mortality rate", 
                "Immature development rate", "Max. biting rate", 
                "Eggs per female per day")
  ) +
  scale_x_continuous(
    limits = Temp_range_for_plot,
    expand = c(0,0)
  ) +
  scale_y_continuous(
    name = unname(TeX("Relative contribution to $\\frac{dR_0}{dT}$ HPD width"))
  ) +
  facet_grid(cols = vars(system_label),
             rows = vars(KH, sigmaH), 
             labeller = labeller(sigmaH = as_labeller(appender_sigmaH, default = label_parsed),
                                 KH = as_labeller(appender_KH, default = label_parsed)),
             scales = "free") +
  theme_minimal_hgrid(16) +
  theme(
    strip.text.x = element_text(size = 10)
  )
plot.ddTR0.rel.sens

# 5) Topt uncertainty -------------------------------------------------------

# Calculate the width of the 95% HPD for the full posterior of Topt across vertebrate host abundance

sigmaH_vec <- c(50, Inf)
KH_vec <- 10^seq(-2,5)

vec <- unique(data.Host$KH)

data.ToptHPD <- dplyr::filter(data.Host, 
                              sigmaH == 50,
                              KH < 1e3) #%>% 
# dplyr::filter(KH %in% unique(KH)[seq(1, length(unique(KH)), length.out = 31)])

full.Topt.HPD <- tibble(system_ID = c(), Temperature = c(), Model = c(),
                        sigmaH = c(), KH = c(), 
                        HPD_low = c(), HPD_high = c(), HPD_width = c())

for (index_KH in unique(data.ToptHPD$KH)) {
  full.Topt.HPD <- expand_grid(dplyr::filter(data.ToptHPD, KH == index_KH), 
                               data.Vec) %>%
    data.table::data.table() %>%
    mutate(RV = ifelse(is.infinite(sigmaH),
                       sigmaV * betaH / (1 / (lf + eps)), # Ross-Macdonald
                       sigmaH * sigmaV * betaH * KH / ((1 / (lf + eps)) * (sigmaH * KH + sigmaV * V0)))) %>%
    mutate(bV = ifelse(is.infinite(sigmaH),
                       sigmaV, # Ross-Macdonald model
                       sigmaV * sigmaH * KH / (sigmaH * KH + sigmaV * V0 + eps))) %>%
    mutate(RH = ifelse(V0 == 0,
                       0,
                       bV * betaV * V0 * exp(-1 / (lf * etaV)) / (KH * (gammaH + muH) + eps))) %>%
    # Basic reproduction number
    mutate(R0 = sqrt(RV*RH)) %>%
    # dplyr::filter to maximum value of R0
    dplyr::filter(R0>0) %>% 
    select(system_ID, sample_num, sigmaH, KH, Temperature, R0) %>% 
    group_by(system_ID, sample_num, sigmaH, KH) %>%
    dplyr::filter(R0 == max(R0)) %>%
    distinct() %>% 
    # Get temperature at which R0 is maximized
    rename(Topt = Temperature) %>% 
    select(system_ID, sigmaH, KH, sample_num, Topt) %>% 
    group_by(system_ID, sigmaH, KH) %>% 
    summarise(
      HPD_low = hdi(Topt, credMass = 0.95)[1],
      HPD_high = hdi(Topt, credMass = 0.95)[2],
      HPD_width = max(eps,HPD_high-HPD_low),
      .groups = "keep"
    ) %>% 
    rbind(full.Topt.HPD)
}



# Save Topt highest posterior density data
# write_rds(full.Topt.HPD, "results/full_Topt_HPD.rds")
full.Topt.HPD <- read_rds("results/full_Topt_HPD.rds")

# # Diagnostic plot
# test.plot <- full.Topt.HPD %>%
#   # dplyr::filter(KH == 1e-2) %>%
#   ggplot(aes(x = KH, y = HPD_width)) +
#   geom_path() +
#   scale_x_continuous(
#     trans = 'log10'
#   ) +
#   facet_grid(rows = vars(system_ID), cols = vars(sigmaH), scales = "free")
# test.plot

# Get focal parameter names
temp_vars <- data.Vec %>% 
  select(-c(system_ID, mosquito_species, pathogen, Temperature, sample_num,
            etaL, muL, KL, V0, lf)) %>% 
  colnames() %>% c("muV")

# Initialize data frame
Topt.HPD <- tibble(system_ID = c(), Temperature = c(), focal_var = c(),
                   sigmaH = c(), KH = c(), 
                   HPD_width = c(), full_HPD_width = c(), rel_HPD_width = c())

# For each focal parameter, calculate relative HPD width
for (var_name in temp_vars) {
  # a) Set all but focal parameter to its posterior mean
  data.HPD.Vec <- data.Vec %>% 
    mutate(muV = 1/lf) %>% 
    select(system_ID, Temperature, sample_num, all_of(var_name), KL) %>% 
    full_join(mean.Vec %>% 
                mutate(muV = 1/lf) %>% 
                select(-all_of(var_name))%>% 
                distinct()) %>% 
    mutate(lf = 1/muV) %>%
    mutate(V0 = ifelse(sigmaV_f * deltaL < (1 / lf),
                       0,
                       KL * rhoL * lf * (1 - 1 / (lf * sigmaV_f * deltaL)))) %>%
    select(-lf)
  
  for (index_KH in unique(data.ToptHPD$KH)) {
    # b) Get posterior samples of R0 (as a function of temperature)
    Topt.HPD <- expand_grid(dplyr::filter(data.ToptHPD, KH == index_KH), 
                            data.HPD.Vec) %>%
      mutate(lf = 1/muV) %>% 
      data.table::data.table() %>%
      mutate(RV = ifelse(is.infinite(sigmaH),
                         sigmaV * betaH / (1 / (lf + eps)), # Ross-Macdonald
                         sigmaH * sigmaV * betaH * KH / ((1 / (lf + eps)) * (sigmaH * KH + sigmaV * V0)))) %>%
      mutate(bV = ifelse(is.infinite(sigmaH),
                         sigmaV, # Ross-Macdonald model
                         sigmaV * sigmaH * KH / (sigmaH * KH + sigmaV * V0 + eps))) %>%
      mutate(RH = ifelse(V0 == 0,
                         0,
                         bV * betaV * V0 * exp(-1 / (lf * etaV)) / (KH * (gammaH + muH) + eps))) %>%
      # Basic reproduction number
      mutate(R0 = sqrt(RV*RH)) %>%
      # dplyr::filter to maximum value of R0
      dplyr::filter(R0>0) %>% 
      select(system_ID, sample_num, sigmaH, KH, Temperature, R0) %>% 
      group_by(system_ID, sample_num, sigmaH, KH) %>%
      dplyr::filter(R0 == max(R0)) %>%
      distinct() %>% 
      # Get temperature at which R0 is maximized
      rename(Topt = Temperature) %>% 
      select(system_ID, sigmaH, KH, sample_num, Topt) %>% 
      group_by(system_ID, sigmaH, KH) %>% 
      summarise(
        HPD_low = hdi(Topt, credMass = 0.95)[1],
        HPD_high = hdi(Topt, credMass = 0.95)[2],
        HPD_width = max(eps, HPD_high-HPD_low),
        .groups = "keep"
      ) %>% 
      select(system_ID, sigmaH, KH, HPD_width) %>% 
      right_join(full.Topt.HPD %>% 
                   select(-c(HPD_low, HPD_high)) %>% 
                   rename(full_HPD_width = HPD_width)) %>% 
      mutate(focal_var = var_name) %>% 
      group_by(system_ID, sigmaH, KH, focal_var) %>% 
      # d) Normalize HPD width by the full HPD width when all parameters are allowed to vary
      mutate(rel_HPD_width = ifelse(full_HPD_width %in% c(0,eps), 0, HPD_width / full_HPD_width)) %>% 
      rbind(Topt.HPD)
  }
}

# Save Topt relative highest posterior density data
# write_rds(Topt.HPD, "results/Topt_HPD_sens.rds")
Topt.HPD <- read_rds("results/Topt_HPD_sens.rds")

## Plot Topt uncertainty 
var_name_table <- list(
  betaV = "Vector competence",
  deltaL = "Pr(larval survival)",
  etaV = "Parasite development rate",
  muV = "Mortality rate",
  rhoL = "Immature development rate",
  sigmaV = "Max. biting rate",
  sigmaV_f = "Eggs per female per day"
)

Topt.uncertainty.plot <- Topt.HPD %>% 
  ungroup() %>% 
  # dplyr::filter(focal_var %in% c("sigmaV_f")) %>%
  mutate(system_ID = case_when(
    system_ID == "Anopheles gambiae / Plasmodium falciparum" ~ "An. gamb. / P. falciparum",
    system_ID == "Aedes aegypti / DENV" ~ "Ae. aegypti / DENV",
    system_ID == "Aedes albopictus / DENV" ~ "Ae. albopictus / DENV",
    system_ID == "Aedes aegypti / ZIKV" ~ "Ae. aegypti / ZIKV",
    system_ID == "Culex quinquefasciatus / WNV" ~ "Cx. quin. / WNV"
  )) %>% 
  arrange(system_ID, sigmaH, KH) %>% 
  # dplyr::filter(KH %in% c(1, 100)) %>%
  dplyr::filter(rel_HPD_width < 1.01) %>%
  # arrange(system_ID, sigmaH, focal_var, Temperature) %>%
  ggplot(aes(x = KH, y = rel_HPD_width, color = focal_var)) +
  geom_path(linewidth = 1) +
  scale_color_discrete(
    name = "Focal parameter",
    breaks = c("betaV", "deltaL", "etaV", "muV", "rhoL", "sigmaV", "sigmaV_f"),
    labels =  c("Vector competence", "Pr(larval survival)", 
                "Parasite development rate", "Mortality rate", 
                "Immature development rate", "Max. biting rate", 
                "Eggs per female per day")
  )  +
  scale_x_continuous(
    name = "Vertebrate host population density",
    trans = 'log10'
  ) +
  scale_y_continuous(
    name = "Relative HPD width",
    breaks = seq(0, 2, by = 0.2)
  ) +
  facet_grid(cols = vars(system_ID), 
             scales = "free") +
  ggtitle("Uncertainty analysis of Topt") +
  theme_minimal_grid(12)
Topt.uncertainty.plot

ggsave("figures/results/Topt_uncertainty.svg", 
       plot = Topt.uncertainty.plot,
       width = 16, height = 9)

# !!! [x] use facet_wrap instead
# !!! [x] make plot "vertical"
# !!! [x] make x and y axis labels descriptive
# !!! [x] make x ticks nice
# !!! [x] make system labels nice
# !!! [x] make focal variable names nice
# "Relative sensitivity" plot
plot.Topt.rel.sens <- Topt.HPD %>%
  ungroup() %>% 
  mutate(system_label = case_when(
    system_ID == "Aedes aegypti / DENV" ~ "(a) *Ae. aegypti* and DENV",
    system_ID == "Aedes aegypti / ZIKV" ~ "(b) *Ae. aegypti* and ZIKV",
    system_ID == "Aedes albopictus / DENV" ~ "(c\u200b) *Ae. albopictus* and DENV",
    system_ID == "Anopheles gambiae / Plasmodium falciparum" ~ "(d) *An. gambiae* and *P. falciparum*",
    system_ID == "Culex quinquefasciatus / WNV" ~ "(e) *Cx. quinquefasciatus* and WNV"
  )) %>%
  filter(!is.na(rel_HPD_width)) %>% 
  # For each KH value, sum up all sensitivity values to get total sensitivity
  group_by(system_ID, sigmaH, KH) %>% 
  mutate(sens_total = sum(rel_HPD_width)) %>% 
  # Divide each sensitivity value by the total sensitivity
  mutate(sens_rel = rel_HPD_width / sens_total) %>% 
  arrange(sens_rel) %>% 
  # Make a stacked plot, colored by parameter
  ggplot(aes(x = KH, y = sens_rel, fill = focal_var)) +
  geom_area() +
  scale_fill_discrete(
    name = "Focal parameter",
    breaks = c("betaV", "deltaL", "etaV", "muV", "rhoL", "sigmaV", "sigmaV_f"),
    labels =  c("Vector competence", "Pr(larval survival)", 
                "Parasite development rate", "Mortality rate", 
                "Immature development rate", "Max. biting rate", 
                "Eggs per female per day")
  )  +
  scale_x_log10(
    name = TeX("Vertebrate host population density (ind/ha)"),
    expand = c(0, 0),
    # limits = c(0.05, 2E3),
    breaks = 10^seq(-1, 3, 1),
    labels = TeX(c("$0.1$", "$1$", "$10", "$100$", "$10^3$"))
  ) +
  scale_y_continuous(
    name = unname(TeX("Relative contribution to $T_{opt}$ HPD width"))
  ) +
  # scale_y_continuous(
  #   limits = c(0,1)
  # ) +
  facet_grid(cols = vars(system_label), scales = "free") +
  theme_minimal_hgrid(16) +
  theme(
    strip.text.x = element_text(size = 14)
  )
plot.Topt.rel.sens

# 6) CTmin/max/width uncertainty ------------------------------------------

# Calculate the width of the 95% HPD for the full posterior of Topt across vertebrate host abundance

sigmaH_vec <- c(50, Inf)
KH_vec <- 10^seq(-2,5)

data.CTHPD <- dplyr::filter(data.Host, 
                            sigmaH %in% sigmaH_vec,
                            KH < 1e3) #%>% 
#dplyr::filter(KH %in% unique(KH)[seq(1, length(unique(KH)), length.out = 11)])

full.CT.HPD <- tibble(system_ID = c(), Temperature = c(), Model = c(),
                      sigmaH = c(), KH = c(), variable = c(),
                      HPD_low = c(), HPD_high = c(), HPD_width = c())


for (index_KH in unique(data.CTHPD$KH)) {
  full.CT.HPD <- expand_grid(dplyr::filter(data.CTHPD, KH == index_KH), 
                             data.Vec) %>%
    data.table::data.table() %>%
    mutate(RV = ifelse(is.infinite(sigmaH),
                       sigmaV * betaH / (1 / (lf + eps)), # Ross-Macdonald
                       sigmaH * sigmaV * betaH * KH / ((1 / (lf + eps)) * (sigmaH * KH + sigmaV * V0)))) %>%
    mutate(bV = ifelse(is.infinite(sigmaH),
                       sigmaV, # Ross-Macdonald model
                       sigmaV * sigmaH * KH / (sigmaH * KH + sigmaV * V0 + eps))) %>%
    mutate(RH = ifelse(V0 == 0,
                       0,
                       bV * betaV * V0 * exp(-1 / (lf * etaV)) / (KH * (gammaH + muH) + eps))) %>%
    # Basic reproduction number
    mutate(R0 = sqrt(RV*RH)) %>%
    # dplyr::filter to maximum value of R0
    group_by(system_ID, sample_num, sigmaH, KH) %>%
    dplyr::filter(R0 > 1) %>%
    # Get lowest temperature at which R0 exceeds one
    mutate(CTmin = min(Temperature)) %>%
    # Get highest temperature at which R0 exceeds one
    mutate(CTmax = max(Temperature)) %>%
    # Get width of critical thermal interval
    mutate(CTwidth = CTmax - CTmin) %>% 
    # Add back values removed from dplyr::filtering R0>1 above
    # and assign the right NA values
    full_join(expand_grid(dplyr::filter(data.CTHPD, KH == index_KH), 
                          data.Vec) %>% 
                select(c(Model, system_ID, sample_num, sigmaH,KH)) %>% 
                distinct()) %>% 
    replace_na(list(CTmin = Inf,
                    CTmax = -Inf,
                    CTwidth = 0)) %>% 
    pivot_longer(cols = c(CTmin, CTmax, CTwidth), names_to = "variable", values_to = "value") %>% 
    select(system_ID, sample_num, sigmaH, KH, variable, value) %>%
    group_by(system_ID, sigmaH, KH, variable) %>% 
    summarise(
      HPD_low = hdi(value, credMass = 0.95)[1],
      HPD_high = hdi(value, credMass = 0.95)[2],
      HPD_width = max(eps, HPD_high-HPD_low),
      .groups = "keep"
    ) %>% 
    rbind(full.CT.HPD)
}

# Save Topt highest posterior density data
# write_rds(full.CT.HPD, "results/full_CT_HPD.rds")
full.CT.HPD <- read_rds("results/full_CT_HPD.rds")

# # Diagnostic plot
# test.plot <- full.CT.HPD %>%
#   filter(is.finite(HPD_width)) %>%
#   # filter(variable == "CTmin") %>%
#   # dplyr::filter(KH == 1e-2) %>%
#   ggplot(aes(x = KH, y = HPD_width, color = as.factor(sigmaH))) +
#   geom_path() +
#   scale_x_continuous(
#     trans = 'log10'
#   ) +
#   facet_grid(rows = vars(system_ID), cols = vars(variable), scales = "free")
# test.plot
# ggsave("figures/results/test_plot.svg", 
#        plot = test.plot,
#        width = 16, height = 9)


# Get focal parameter names
temp_vars <- data.Vec %>% 
  select(-c(system_ID, mosquito_species, pathogen, Temperature, sample_num,
            etaL, muL, KL, V0, lf)) %>% 
  colnames() %>% c("muV")

# Initialize data frame
CT.HPD <- tibble(system_ID = c(), Temperature = c(), focal_var = c(),
                 sigmaH = c(), KH = c(), variable = c(),
                 HPD_width = c(), full_HPD_width = c(), rel_HPD_width = c())



# For each focal parameter, calculate relative HPD width
for (var_name in temp_vars) {
  print(paste0("Focal variable: ", var_name))
  # a) Set all but focal parameter to its posterior mean
  data.HPD.Vec <- data.Vec %>% 
    mutate(muV = 1/lf) %>% 
    select(system_ID, Temperature, sample_num, all_of(var_name), KL) %>% 
    full_join(mean.Vec %>% 
                mutate(muV = 1/lf) %>% 
                select(-all_of(var_name))%>% 
                distinct(),
              by = c("system_ID", "Temperature")) %>% 
    mutate(lf = 1/muV) %>%
    mutate(V0 = ifelse(sigmaV_f * deltaL < (1 / lf),
                       0,
                       KL * rhoL * lf * (1 - 1 / (lf * sigmaV_f * deltaL)))) %>%
    select(-lf)
  pb <- progress_bar$new(
    format = ":spin :system progress = :percent [:bar] :elapsed | eta: :eta",
    total = length(unique(data.CTHPD$KH)),
    width = 120)    
  for (index_KH in unique(data.CTHPD$KH)) {
    pb$tick()
    # b) Get posterior samples of R0 (as a function of temperature)
    CT.HPD <- expand_grid(dplyr::filter(data.CTHPD, KH == index_KH), 
                          data.HPD.Vec) %>%
      mutate(lf = 1/muV) %>% 
      data.table::data.table() %>%
      mutate(RV = ifelse(is.infinite(sigmaH),
                         sigmaV * betaH / (1 / (lf + eps)), # Ross-Macdonald
                         sigmaH * sigmaV * betaH * KH / ((1 / (lf + eps)) * (sigmaH * KH + sigmaV * V0)))) %>%
      mutate(bV = ifelse(is.infinite(sigmaH),
                         sigmaV, # Ross-Macdonald model
                         sigmaV * sigmaH * KH / (sigmaH * KH + sigmaV * V0 + eps))) %>%
      mutate(RH = ifelse(V0 == 0,
                         0,
                         bV * betaV * V0 * exp(-1 / (lf * etaV)) / (KH * (gammaH + muH) + eps))) %>%
      # Basic reproduction number
      mutate(R0 = sqrt(RV*RH)) %>%
      # dplyr::filter to maximum value of R0
      group_by(system_ID, sample_num, sigmaH, KH) %>%
      dplyr::filter(R0 > 1) %>%
      # Get lowest temperature at which R0 exceeds one
      mutate(CTmin = min(Temperature)) %>%
      # Get highest temperature at which R0 exceeds one
      mutate(CTmax = max(Temperature)) %>%
      # Get width of critical thermal interval
      mutate(CTwidth = CTmax - CTmin) %>% 
      select(system_ID, sample_num, sigmaH, KH, CTmin, CTmax, CTwidth) %>% 
      distinct() %>% 
      # Add back values removed from dplyr::filtering R0>1 above
      # and assign the right NA values
      # full_join(expand_grid(dplyr::filter(data.CTHPD, KH == index_KH), 
      #                       data.Vec) %>% 
      #             select(c(Model, system_ID, sample_num, sigmaH,KH)) %>% 
      #             distinct()) %>% 
      replace_na(list(CTmin = Inf,
                      CTmax = -Inf,
                      CTwidth = 0)) %>% 
      pivot_longer(cols = c(CTmin, CTmax, CTwidth), names_to = "variable", values_to = "value") %>% 
      group_by(system_ID, sigmaH, KH, variable) %>% 
      summarise(
        HPD_low = hdi(value, credMass = 0.95)[1],
        HPD_high = hdi(value, credMass = 0.95)[2],
        HPD_width = max(eps, HPD_high-HPD_low),
        .groups = "keep"
      ) %>% 
      select(system_ID, sigmaH, KH, variable, HPD_width) %>% 
      right_join(full.CT.HPD %>% 
                   select(-c(HPD_low, HPD_high)) %>% 
                   rename(full_HPD_width = HPD_width) %>% 
                   filter(KH == index_KH),
                 by = c("system_ID", "sigmaH", "KH", "variable")) %>% 
      mutate(focal_var = var_name) %>% 
      group_by(system_ID, sigmaH, KH, focal_var) %>% 
      # d) Normalize HPD width by the full HPD width when all parameters are allowed to vary
      mutate(rel_HPD_width = ifelse(full_HPD_width %in% c(0,eps), 0, HPD_width / full_HPD_width)) %>% 
      rbind(CT.HPD)
  }
  pb$terminate()
}


# Save CT relative highest posterior density data
# write_rds(CT.HPD, "results/CT_HPD_sens.rds", compress = 'gz')
CT.HPD <- read_rds("results/CT_HPD_sens.rds")

## Plot Topt uncertainty 
var_name_table <- list(
  betaV = c(TeX("$\\beta_V$")),
  deltaL = c(TeX("$\\delta_L$")),
  etaV = c(TeX("$\\eta_V$")),
  muV = c(TeX("$\\mu_V$")),
  rhoL = c(TeX("$\\rho_L$")),
  sigmaV = c(TeX("$\\sigma_V$")),
  sigmaV_f = c(TeX("$\\sigma_v f$"))
)

# !!! [x] re-run analysis since CTwidth is coming up as NA
# !!! [x] make relative sensitivity plot like was done for R0
# CTwidth
CTwidth.uncertainty.plot <- CT.HPD %>%   
  ungroup() %>% 
  filter(variable == "CTwidth") %>% 
  filter(!is.na(HPD_width)) %>% 
  # dplyr::filter(focal_var %in% c("sigmaV_f")) %>%
  mutate(system_ID = case_when(
    system_ID == "Anopheles gambiae / Plasmodium falciparum" ~ "An. gamb. / P. falciparum",
    system_ID == "Aedes aegypti / DENV" ~ "Ae. aegypti / DENV",
    system_ID == "Aedes albopictus / DENV" ~ "Ae. albopictus / DENV",
    system_ID == "Aedes aegypti / ZIKV" ~ "Ae. aegypti / ZIKV",
    system_ID == "Culex quinquefasciatus / WNV" ~ "Cx. quin. / WNV"
  )) %>% 
  arrange(system_ID, sigmaH) %>% 
  # dplyr::filter(KH %in% c(1, 100)) %>%
  dplyr::filter(rel_HPD_width < 1.2) %>%
  arrange(system_ID, sigmaH, KH) %>%
  ggplot(aes(x = KH, y = rel_HPD_width, color = focal_var, linetype = sigmaH)) +
  geom_path(linewidth = 1) +
  scale_linetype_binned() +
  scale_color_discrete(
    name = "Focal parameter",
    breaks = c("betaV", "deltaL", "etaV", "muV", "rhoL", "sigmaV", "sigmaV_f"),
    labels =  c("Vector competence", "Pr(larval survival)", 
                "Parasite development rate", "Mortality rate", 
                "Immature development rate", "Max. biting rate", 
                "Eggs per female per day")
  )  +
  scale_x_continuous(
    name = "Vertebrate host population density",
    trans = 'log10'
  ) +
  scale_y_continuous(
    name = "Relative HPD width",
    breaks = seq(0, 2, by = 0.2)
  ) +
  facet_grid(cols = vars(sigmaH), rows = vars(system_ID), scales = "free",
             labeller = labeller(sigmaH = as_labeller(appender_sigmaH,
                                                      default = label_parsed),
                                 KH = as_labeller(appender_KH, 
                                                  default = label_parsed)),) +
  ggtitle("Uncertainty analysis of CTwidth") +
  theme_minimal_grid(12)
CTwidth.uncertainty.plot

ggsave("figures/results/CTwidth_uncertainty.svg", 
       plot = CTwidth.uncertainty.plot,
       width = 16, height = 9)

plot.CTwidth.rel.sens <- CT.HPD %>%
  ungroup() %>% 
  filter(variable == "CTwidth") %>% 
  mutate(system_label = case_when(
    system_ID == "Aedes aegypti / DENV" ~ "(a) *Ae. aegypti* and DENV",
    system_ID == "Aedes aegypti / ZIKV" ~ "(b) *Ae. aegypti* and ZIKV",
    system_ID == "Aedes albopictus / DENV" ~ "(c\u200b) *Ae. albopictus* and DENV",
    system_ID == "Anopheles gambiae / Plasmodium falciparum" ~ "(d) *An. gambiae* and *P. falciparum*",
    system_ID == "Culex quinquefasciatus / WNV" ~ "(e) *Cx. quinquefasciatus* and WNV"
  )) %>%
  filter(!is.na(rel_HPD_width)) %>% 
  # For each KH value, sum up all sensitivity values to get total sensitivity
  group_by(system_ID, sigmaH, KH) %>% 
  mutate(sens_total = sum(rel_HPD_width)) %>% 
  # Divide each sensitivity value by the total sensitivity
  mutate(sens_rel = rel_HPD_width / sens_total) %>% 
  arrange(sens_rel) %>% 
  # Make a stacked plot, colored by parameter
  ggplot(aes(x = KH, y = sens_rel, fill = focal_var)) +
  geom_area() +
  scale_fill_discrete(
    name = "Focal parameter",
    breaks = c("betaV", "deltaL", "etaV", "muV", "rhoL", "sigmaV", "sigmaV_f"),
    labels =  c("Vector competence", "Pr(larval survival)", 
                "Parasite development rate", "Mortality rate", 
                "Immature development rate", "Max. biting rate", 
                "Eggs per female per day")
  )  +
  scale_x_log10(
    name = TeX("Vertebrate host population density (ind/ha)"),
    expand = c(0, 0),
    # limits = c(0.9, 1.1E3),
    breaks = 10^seq(-2, 3, 1),
    labels = TeX(c("$0.01","$0.1$", "$1$", "$10", "$100$", "$10^3$"))
  ) +
  scale_y_continuous(
    name = unname(TeX("Relative contribution to $CT_{width}$ HPD width"))
  ) +
  # scale_y_continuous(
  #   limits = c(0,1)
  # ) +
  facet_grid(cols = vars(system_label),
             rows = vars(sigmaH), 
             labeller = labeller(sigmaH = as_labeller(appender_sigmaH,
                                                      default = label_parsed)),
             scales = "free") +
  theme_minimal_hgrid(16) +
  theme(
    strip.text.x = element_text(size = 10)
  )
plot.CTwidth.rel.sens

# !!! [x] make relative sensitivity plot like was done for R0
# CTmin
CTmin.uncertainty.plot <- CT.HPD %>% 
  filter(variable == "CTmin") %>% 
  filter(!is.na(HPD_width)) %>% 
  ungroup() %>% 
  # dplyr::filter(focal_var %in% c("sigmaV_f")) %>%
  mutate(system_ID = case_when(
    system_ID == "Anopheles gambiae / Plasmodium falciparum" ~ "An. gamb. / P. falciparum",
    system_ID == "Aedes aegypti / DENV" ~ "Ae. aegypti / DENV",
    system_ID == "Aedes albopictus / DENV" ~ "Ae. albopictus / DENV",
    system_ID == "Aedes aegypti / ZIKV" ~ "Ae. aegypti / ZIKV",
    system_ID == "Culex quinquefasciatus / WNV" ~ "Cx. quin. / WNV"
  )) %>% 
  arrange(system_ID, sigmaH) %>% 
  # dplyr::filter(KH %in% c(1, 100)) %>%
  dplyr::filter(rel_HPD_width < 1.2) %>%
  arrange(system_ID, sigmaH, KH) %>%
  ggplot(aes(x = KH, y = rel_HPD_width, color = focal_var, linetype = sigmaH)) +
  geom_path(linewidth = 1) +
  scale_linetype_binned() +
  scale_color_discrete(
    name = "Focal parameter",
    breaks = c("betaV", "deltaL", "etaV", "muV", "rhoL", "sigmaV", "sigmaV_f"),
    labels =  c("Vector competence", "Pr(larval survival)", 
                "Parasite development rate", "Mortality rate", 
                "Immature development rate", "Max. biting rate", 
                "Eggs per female per day")
  )  +
  scale_x_continuous(
    name = "Vertebrate host population density",
    trans = 'log10'
  ) +
  scale_y_continuous(
    name = "Relative HPD width",
    breaks = seq(0, 2, by = 0.2)
  ) +
  facet_grid(cols = vars(sigmaH), rows = vars(system_ID), scales = "free",
             labeller = labeller(sigmaH = as_labeller(appender_sigmaH,
                                                      default = label_parsed),
                                 KH = as_labeller(appender_KH, 
                                                  default = label_parsed)),) +
  ggtitle("Uncertainty analysis of CTmin") +
  theme_minimal_grid(12)
CTmin.uncertainty.plot

ggsave("figures/results/CTmin_uncertainty.svg", 
       plot = CTmin.uncertainty.plot,
       width = 16, height = 9)

plot.CTmin.rel.sens <- CT.HPD %>%
  ungroup() %>% 
  filter(variable == "CTmin") %>% 
  mutate(system_label = case_when(
    system_ID == "Aedes aegypti / DENV" ~ "(a) *Ae. aegypti* and DENV",
    system_ID == "Aedes aegypti / ZIKV" ~ "(b) *Ae. aegypti* and ZIKV",
    system_ID == "Aedes albopictus / DENV" ~ "(c\u200b) *Ae. albopictus* and DENV",
    system_ID == "Anopheles gambiae / Plasmodium falciparum" ~ "(d) *An. gambiae* and *P. falciparum*",
    system_ID == "Culex quinquefasciatus / WNV" ~ "(e) *Cx. quinquefasciatus* and WNV"
  )) %>%
  filter(!is.na(rel_HPD_width)) %>% 
  # For each KH value, sum up all sensitivity values to get total sensitivity
  group_by(system_ID, sigmaH, KH) %>% 
  mutate(sens_total = sum(rel_HPD_width)) %>% 
  # Divide each sensitivity value by the total sensitivity
  mutate(sens_rel = rel_HPD_width / sens_total) %>% 
  arrange(sens_rel) %>% 
  # Make a stacked plot, colored by parameter
  ggplot(aes(x = KH, y = sens_rel, fill = focal_var)) +
  geom_area() +
  scale_fill_discrete(
    name = "Focal parameter",
    breaks = c("betaV", "deltaL", "etaV", "muV", "rhoL", "sigmaV", "sigmaV_f"),
    labels =  c("Vector competence", "Pr(larval survival)", 
                "Parasite development rate", "Mortality rate", 
                "Immature development rate", "Max. biting rate", 
                "Eggs per female per day")
  )  +
  scale_x_log10(
    name = TeX("Vertebrate host population density (ind/ha)"),
    expand = c(0, 0),
    # limits = c(0.9, 1.1E3),
    breaks = 10^seq(-2, 3, 1),
    labels = TeX(c("$0.01","$0.1$", "$1$", "$10", "$100$", "$10^3$"))
  ) +
  scale_y_continuous(
    name = unname(TeX("Relative contribution to $CT_{min}$ HPD width"))
  ) +
  # scale_y_continuous(
  #   limits = c(0,1)
  # ) +
  facet_grid(cols = vars(system_label),
             rows = vars(sigmaH), 
             labeller = labeller(sigmaH = as_labeller(appender_sigmaH,
                                                      default = label_parsed)),
             scales = "free") +
  theme_minimal_hgrid(16) +
  theme(
    strip.text.x = element_text(size = 10)
  )
plot.CTmin.rel.sens


# !!! [x] make relative sensitivity plot like was done for R0
# CTmax
CTmax.uncertainty.plot <- CT.HPD %>% 
  filter(variable == "CTmax") %>% 
  filter(!is.na(HPD_width)) %>% 
  ungroup() %>% 
  # dplyr::filter(focal_var %in% c("sigmaV_f")) %>%
  mutate(system_ID = case_when(
    system_ID == "Anopheles gambiae / Plasmodium falciparum" ~ "An. gamb. / P. falciparum",
    system_ID == "Aedes aegypti / DENV" ~ "Ae. aegypti / DENV",
    system_ID == "Aedes albopictus / DENV" ~ "Ae. albopictus / DENV",
    system_ID == "Aedes aegypti / ZIKV" ~ "Ae. aegypti / ZIKV",
    system_ID == "Culex quinquefasciatus / WNV" ~ "Cx. quin. / WNV"
  )) %>% 
  arrange(system_ID, sigmaH) %>% 
  # dplyr::filter(KH %in% c(1, 100)) %>%
  dplyr::filter(rel_HPD_width < 1.2) %>%
  arrange(system_ID, sigmaH, KH) %>%
  ggplot(aes(x = KH, y = rel_HPD_width, color = focal_var, linetype = sigmaH)) +
  geom_path(linewidth = 1) +
  scale_linetype_binned() +
  scale_color_discrete(
    name = "Focal parameter",
    breaks = c("betaV", "deltaL", "etaV", "muV", "rhoL", "sigmaV", "sigmaV_f"),
    labels =  c("Vector competence", "Pr(larval survival)", 
                "Parasite development rate", "Mortality rate", 
                "Immature development rate", "Max. biting rate", 
                "Eggs per female per day")
  )  +
  scale_x_continuous(
    name = "Vertebrate host population density",
    trans = 'log10'
  ) +
  scale_y_continuous(
    name = "Relative HPD width",
    breaks = seq(0, 2, by = 0.2)
  ) +
  facet_grid(cols = vars(sigmaH), rows = vars(system_ID), scales = "free",
             labeller = labeller(sigmaH = as_labeller(appender_sigmaH,
                                                      default = label_parsed),
                                 KH = as_labeller(appender_KH, 
                                                  default = label_parsed)),) +
  ggtitle("Uncertainty analysis of CTmax") +
  theme_minimal_grid(12)
CTmax.uncertainty.plot

ggsave("figures/results/CTmax_uncertainty.svg", 
       plot = CTmax.uncertainty.plot,
       width = 16, height = 9)

plot.CTmax.rel.sens <- CT.HPD %>%
  ungroup() %>% 
  filter(variable == "CTmax") %>% 
  mutate(system_label = case_when(
    system_ID == "Aedes aegypti / DENV" ~ "(a) *Ae. aegypti* and DENV",
    system_ID == "Aedes aegypti / ZIKV" ~ "(b) *Ae. aegypti* and ZIKV",
    system_ID == "Aedes albopictus / DENV" ~ "(c\u200b) *Ae. albopictus* and DENV",
    system_ID == "Anopheles gambiae / Plasmodium falciparum" ~ "(d) *An. gambiae* and *P. falciparum*",
    system_ID == "Culex quinquefasciatus / WNV" ~ "(e) *Cx. quinquefasciatus* and WNV"
  )) %>%
  filter(!is.na(rel_HPD_width)) %>% 
  # For each KH value, sum up all sensitivity values to get total sensitivity
  group_by(system_ID, sigmaH, KH) %>% 
  mutate(sens_total = sum(rel_HPD_width)) %>% 
  # Divide each sensitivity value by the total sensitivity
  mutate(sens_rel = rel_HPD_width / sens_total) %>% 
  arrange(sens_rel) %>% 
  # Make a stacked plot, colored by parameter
  ggplot(aes(x = KH, y = sens_rel, fill = focal_var)) +
  geom_area() +
  scale_fill_discrete(
    name = "Focal parameter",
    breaks = c("betaV", "deltaL", "etaV", "muV", "rhoL", "sigmaV", "sigmaV_f"),
    labels =  c("Vector competence", "Pr(larval survival)", 
                "Parasite development rate", "Mortality rate", 
                "Immature development rate", "Max. biting rate", 
                "Eggs per female per day")
  )  +
  scale_x_log10(
    name = TeX("Vertebrate host population density (ind/ha)"),
    expand = c(0, 0),
    # limits = c(0.9, 1.1E3),
    breaks = 10^seq(-2, 3, 1),
    labels = TeX(c("$0.01","$0.1$", "$1$", "$10", "$100$", "$10^3$"))
  ) +
  scale_y_continuous(
    name = unname(TeX("Relative contribution to $CT_{max}$ HPD width"))
  ) +
  # scale_y_continuous(
  #   limits = c(0,1)
  # ) +
  facet_grid(cols = vars(system_label),
             rows = vars(sigmaH), 
             labeller = labeller(sigmaH = as_labeller(appender_sigmaH,
                                                      default = label_parsed)),
             scales = "free") +
  theme_minimal_hgrid(16) +
  theme(
    strip.text.x = element_text(size = 10)
  )
plot.CTmax.rel.sens


# 4) Densities of Topt, CT min, max, and width ----------------------------

###* Density of Topt ----
data.Topt <- data.Host %>%
  dplyr::filter(sigmaH %in% c(100, Inf)) %>%
  dplyr::filter(KH %in% 10^seq(-2,4))

Topt.density.df <- expand_grid(data.Topt, data.Vec) %>% 
  data.table::data.table() %>%
  mutate(RV = ifelse(is.infinite(sigmaH),
                     sigmaV * betaH / (1 / (lf + eps)), # Ross-Macdonald
                     sigmaH * sigmaV * betaH * KH / ((1 / (lf + eps)) * (sigmaH * KH + sigmaV * V0)))) %>%
  mutate(bV = ifelse(is.infinite(sigmaH),
                     sigmaV, # Ross-Macdonald model
                     sigmaV * sigmaH * KH / (sigmaH * KH + sigmaV * V0 + eps))) %>%
  mutate(RH = ifelse(V0 == 0,
                     0,
                     bV * betaV * V0 * exp(-1 / (lf * etaV)) / (KH * (gammaH + muH) + eps))) %>%
  # Basic reproduction number
  mutate(R0 = sqrt(RV*RH)) %>%
  # dplyr::filter to maximum value of R0
  dplyr::filter(R0>0) %>% 
  select(system_ID, sample_num, Model, sigmaH, KH, Temperature, R0) %>% 
  group_by(system_ID, sample_num, sigmaH, KH) %>%
  dplyr::filter(R0 == max(R0)) %>%
  distinct() %>% 
  # Get temperature at which R0 is maximized
  rename(Topt = Temperature) %>%
  dplyr::select(system_ID, sample_num, Model, sigmaH, KH, Topt)

Topt.density.df <- full_join(Topt.density.df,
                             expand_grid(data.Topt, data.Vec) %>%
                               select(c(Model, system_ID, sample_num, sigmaH,KH)) %>%
                               distinct())

plot.Topt.density <- Topt.density.df %>%
  filter(is.infinite(sigmaH)) %>% 
  ggplot(aes(x = Topt)) +
  geom_density(aes(color = as.factor(KH), linetype = as.factor(sigmaH)),
               lwd = 1, adjust = 1) +
  scale_color_viridis_d(
    name = "Vertebrate host population\ndensity (ind/ha)",
    breaks = 10^seq(-2, 5),
    labels = unname(c(0.01, 0.1, 1, 10, 100, TeX("$10^3$"), TeX("$10^4$"), TeX("$10^5$"))),
    option = "plasma"
  ) +
  facet_wrap(~system_ID, ncol = 1, scales = "free") +
  theme_minimal_grid()

###* Densities of CTmin,max,width ----

# Initialize data frame
CT.density.df <- init.df

data.CT <- data.Host %>%
  dplyr::filter(sigmaH %in% c(100, Inf)) %>%
  dplyr::filter(KH %in% 10^seq(-2,4))

# Slice host trait data
sigmaH_slices <- slice(unique(data.CT$sigmaH), 2)
KH_slices <- slice(unique(data.CT$KH), cluster_size)

CT.density.df <- expand_grid(data.CT, data.Vec) %>% 
  data.table::data.table() %>%
  mutate(RV = ifelse(is.infinite(sigmaH),
                     sigmaV * betaH / (1 / (lf + eps)), # Ross-Macdonald
                     sigmaH * sigmaV * betaH * KH / ((1 / (lf + eps)) * (sigmaH * KH + sigmaV * V0)))) %>%
  mutate(bV = ifelse(is.infinite(sigmaH),
                     sigmaV, # Ross-Macdonald model
                     sigmaV * sigmaH * KH / (sigmaH * KH + sigmaV * V0 + eps))) %>%
  mutate(RH = ifelse(V0 == 0,
                     0,
                     bV * betaV * V0 * exp(-1 / (lf * etaV)) / (KH * (gammaH + muH) + eps))) %>%
  # Basic reproduction number
  mutate(R0 = sqrt(RV*RH)) %>%
  # filter to temperatures where R0 exceeds one
  group_by(system_ID, sample_num, sigmaH, KH) %>%
  dplyr::filter(R0 > 1) %>%
  # Get lowest temperature at which R0 exceeds one
  mutate(CTmin = min(Temperature)) %>%
  # Get highest temperature at which R0 exceeds one
  mutate(CTmax = max(Temperature)) %>%
  # Get width of critical thermal interval
  mutate(CTwidth = CTmax - CTmin) %>%
  dplyr::select(system_ID, sample_num, Model, sigmaH, KH, CTmin, CTmax, CTwidth) %>% 
  distinct()

# Add back values removed from dplyr::filtering R0>1 above
CT.density.df <- full_join(CT.density.df,
                           expand_grid(data.CT, data.Vec) %>%
                             select(c(Model, system_ID, sample_num, sigmaH,KH)) %>%
                             distinct()) %>%
  replace_na(list(CTmin = Inf,
                  CTmax = -Inf,
                  CTwidth = 0))

# Density plots of CTmin, Topt, CTmax

# Combine data frames
all.density.df <- right_join(Topt.density.df, CT.density.df) %>% 
  select(-CTwidth) %>% 
  pivot_longer(cols = c(CTmin, Topt, CTmax), names_to = "variable", values_to = "value") %>% 
  mutate(system_label = case_when(
    system_ID == "Aedes aegypti / DENV" ~ "(a) *Ae. aegypti* and DENV",
    system_ID == "Aedes aegypti / ZIKV" ~ "(b) *Ae. aegypti* and ZIKV",
    system_ID == "Aedes albopictus / DENV" ~ "(c\u200b) *Ae. albopictus* and DENV",
    system_ID == "Anopheles gambiae / Plasmodium falciparum" ~ "(d) *An. gambiae* and *P. falciparum*",
    system_ID == "Culex quinquefasciatus / WNV" ~ "(e) *Cx. quinquefasciatus* and WNV"
  )) %>%
  mutate(variable = case_when(
    variable == "CTmin" ~ "Critical thermal minimum",
    variable == "Topt" ~ "Thermal optimum",
    variable == "CTmax" ~ "Critical thermal maximum"
  )) %>% 
  filter(KH < 1e3) %>% 
  mutate(KH_factor = factor(KH, levels = c("0.01", "0.1", "1", "10", "100")))

all.density.df$variable <- factor(all.density.df$variable, 
                                  levels = c("Critical thermal minimum", "Thermal optimum","Critical thermal maximum"))
all.density.df$sigmaH <- factor(all.density.df$sigmaH, 
                                levels = c(100, Inf))

write_rds(all.density.df, "results/CTminmaxTopt_densities.rds")

# !!! [] clean up code
# !!! [] make y label nice
# !!! [] make system labels nice (decrease font size)
# !!! [x] make single box and whisker plot? x = Temperature, y = pop. density, show distributions with box&whisker

### Option 1: inline violin plots ----
plot.densities_violin_inline <- ggplot() +
  # geom_boxplot(data = all.density.df,
  #             aes(x = value, y = KH_factor, fill = variable, linetype = sigmaH)) +
  geom_violin(data = filter(all.density.df, variable == "Thermal optimum"),
              aes(x = value, y = KH_factor, linetype = sigmaH, fill = "Topt")
  ) +
  geom_violin(data = filter(all.density.df, variable == "Critical thermal minimum"),
              aes(x = value, y = KH_factor, linetype = sigmaH, fill = "CTmin"),
  ) +
  geom_violin(data = filter(all.density.df, variable == "Critical thermal maximum"),
              aes(x = value, y = KH_factor, linetype = sigmaH, fill = "CTmax"),
  ) +
  scale_x_continuous(
    name = "Temperature (°C)",
    expand = c(0,0)
  ) +
  scale_y_discrete(
    name = "Vertebrate host population density (ind/ha)"
  ) +
  scale_linetype(
    name = "Biting tolerance\n(bites per host per day)"
  ) +
  scale_fill_discrete(
    name = "",
    breaks = c("CTmin", "Topt", "CTmax"),
    labels = c("Critical thermal minimum", "Thermal optimum", "Critical thermal maximum")
  ) +
  facet_grid(cols = vars(system_label)) +
  theme_minimal_grid()
plot.densities_violin_inline

### Option 2: inline box and whisker plots ----
plot.densities_boxwhisker_inline <- ggplot() +
  # geom_boxplot(data = all.density.df,
  #             aes(x = value, y = KH_factor, fill = variable, linetype = sigmaH)) +
  geom_boxplot(data = filter(all.density.df, variable == "Thermal optimum"),
               aes(x = value, y = KH_factor, linetype = sigmaH, fill = "Topt")
  ) +
  geom_boxplot(data = filter(all.density.df, variable == "Critical thermal minimum"),
               aes(x = value, y = KH_factor, linetype = sigmaH, fill = "CTmin"),
  ) +
  geom_boxplot(data = filter(all.density.df, variable == "Critical thermal maximum"),
               aes(x = value, y = KH_factor, linetype = sigmaH, fill = "CTmax"),
  ) +
  scale_x_continuous(
    name = "Temperature (°C)",
    expand = c(0,0)
  ) +
  scale_y_discrete(
    name = "Vertebrate host population density (ind/ha)"
  ) +
  scale_linetype(
    name = "Biting tolerance\n(bites per host per day)"
  ) +
  scale_fill_discrete(
    name = "",
    breaks = c("CTmin", "Topt", "CTmax"),
    labels = c("Critical thermal minimum", "Thermal optimum", "Critical thermal maximum")
  ) +
  facet_grid(cols = vars(system_label)) +
  theme_minimal_grid()
plot.densities_boxwhisker_inline

### Option 3: stacked violin plots ----
plot.densities_violin_stacked <- ggplot() +
  geom_violin(data = all.density.df,
              aes(x = value, y = KH_factor, fill = variable, linetype = sigmaH)) +
  scale_x_continuous(
    name = "Temperature (°C)",
    expand = c(0,0)
  ) +
  scale_y_discrete(
    name = "Vertebrate host population density (ind/ha)"
  ) +
  scale_linetype(
    name = "Biting tolerance\n(bites per host per day)"
  ) +
  scale_fill_discrete(
    name = ""
  ) +
  facet_grid(cols = vars(system_label)) +
  theme_minimal_grid()
plot.densities_violin_stacked

### Option 4: stacked box and whisker plots ----
plot.densities_boxwhisker_stacked <- ggplot() +
  geom_boxplot(data = all.density.df,
              aes(x = value, y = KH_factor, fill = variable, linetype = sigmaH)) +
  scale_x_continuous(
    name = "Temperature (°C)",
    expand = c(0,0)
  ) +
  scale_y_discrete(
    name = "Vertebrate host population density (ind/ha)"
  ) +
  scale_linetype(
    name = "Biting tolerance\n(bites per host per day)"
  ) +
  scale_fill_discrete(
    name = ""
  ) +
  facet_grid(cols = vars(system_label)) +
  theme_minimal_grid()
plot.densities_boxwhisker_stacked

# Dynamic: sigmaH = 100
plot.densities_dynamic <- all.density.df  %>% 
  dplyr::filter(sigmaH %in% 100) %>%
  ggplot(aes(x = value)) +
  geom_density(aes(color = as.factor(KH)),
               lwd = 1.5) +
  scale_color_viridis_d(
    name = "Vertebrate host population\ndensity (ind/ha)",
    breaks = 10^seq(-2, 4),
    labels = unname(c(0.01, 0.1, 1, 10, 100, TeX("$10^3$"), TeX("$10^4$"))),
    limits = 10^seq(-2, 4),
    option = "plasma"
  ) +
  scale_x_continuous(name = "Temperature") +
  facet_grid(rows = vars(system_ID), cols = vars(variable), scales = "free") +
  theme_minimal_grid() 

plot.densities_dynamic

ggsave("figures/results/thermal_densities_dynamic.svg",
       plot = plot.densities_dynamic,
       width = 16, height = 9)

# Static: sigmaH = Inf
plot.densities_RM <- all.density.df  %>% 
  dplyr::filter(sigmaH %in% Inf) %>%
  # dplyr::filter(KH < 1E3) %>%
  ggplot(aes(x = value)) +
  geom_density(aes(color = as.factor(KH)),
               lwd = 1.5) +
  scale_color_viridis_d(
    name = "Vertebrate host population\ndensity (ind/ha)",
    breaks = 10^seq(-2, 4),
    labels = unname(c(0.01, 0.1, 1, 10, 100, TeX("$10^3$"), TeX("$10^4$"))),
    limits = 10^seq(-2, 4),
    option = "plasma"
  ) +
  scale_x_continuous(name = "Temperature") +
  facet_grid(rows = vars(system_ID), cols = vars(variable), scales = "free") +
  theme_minimal_grid() 

plot.densities_RM

ggsave("figures/results/thermal_densities_RM.svg",
       plot = plot.densities_RM,
       width = 16, height = 9)

# !!! Consider adding an additional plot showing how the extrema vary more with sigmaH when KH is low

###* Density plot of CTwidth ----
plot.CTwidth.density <- CT.density.df %>% 
  dplyr::filter(CTwidth > 0) %>% 
  ggplot(aes(x = CTwidth)) +
  geom_density(aes(color = as.factor(KH), linetype = as.factor(sigmaH)),
               lwd = 1, adjust = 1) +
  # geom_density(data = dplyr::filter(Topt.density.df, sigmaH == Inf),
  #              linetype = 1,
  #              lwd = 2) +
  # color: scaled log10, color-blind friendly
  scale_color_viridis_d(
    name = "Vertebrate host population\ndensity (ind/ha)",
    breaks = 10^seq(-2, 5),
    labels = unname(c(0.01, 0.1, 1, 10, 100, TeX("$10^3$"), TeX("$10^4$"), TeX("$10^5$"))),
    option = "plasma"
  ) +
  facet_wrap(~system_ID, scales = "free", ncol = 1) +
  theme_minimal_grid()
plot.CTwidth.density

plot.R0.density <- data.R0 %>% 
  dplyr::filter(Temperature == 10) %>% 
  ggplot(aes(x = median)) +
  geom_density(aes(color = as.factor(KH), linetype = as.factor(sigmaH)),
               lwd = 1, adjust = 1) +
  # geom_density(data = dplyr::filter(Topt.density.df, sigmaH == Inf),
  #              linetype = 1,
  #              lwd = 2) +
  # color: scaled log10, color-blind friendly
  scale_color_viridis_d(
    name = "Vertebrate host population\ndensity (ind/ha)",
    breaks = 10^seq(-2, 5),
    labels = unname(c(0.01, 0.1, 1, 10, 100, TeX("$10^3$"), TeX("$10^4$"), TeX("$10^5$"))),
    option = "plasma"
  ) +
  facet_wrap(~system_ID, scales = "free", ncol = 1) +
  theme_minimal_grid()
plot.R0.density












