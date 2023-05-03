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

# Set up parallel
if (!exists("cluster")) {
  cluster_size <- min(21, parallel::detectCores() - 2)
  cluster <- new_cluster(cluster_size)
  cluster_library(cluster, c("dplyr", "tidyr"))
}

eps <- .Machine$double.eps

# Initialize a data frame to save our analyses
init.df <- tibble(system_ID = c(), Temperature = c(), Model = c(),
                  sigmaH = c(), KH = c(), variable = c(),
                  lowHCI = c(), highHCI = c(), mean = c(), median = c())


# 1) Define accessory functions -------------------------------------------

# Function: Slice data to optimize usage of parallel processing and memory
slice<-function(x,n) {
  N<-length(x);
  lapply(seq(1,N,n),function(i) x[i:min(i+n-1,N)])
}

# Function: evaluate R0 (as a function of temperature), its mean, median, and 
#           highest prob. intervals across samples
R0_TPC_func <- function(in_df, system_name) {
  out_df <- in_df %>%
    expand_grid(filter(data.Vec, system_ID == system_name)) %>%
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
    mutate(R0 = sqrt(RV*RH))%>%
    dplyr::select(system_ID, sample_num, Temperature, Model, sigmaH, KH, R0) %>%
    # Normalize R0 across temperature
    group_by(system_ID, Model, sigmaH, KH) %>%
    ungroup() %>%
    pivot_longer(cols = c(R0), names_to = "variable", values_to = "value") %>%
    group_by(system_ID, Temperature, Model, sigmaH, KH, variable) %>%
    partition(cluster) %>%
    summarise(
      lowHCI = quantile(value, 0.055),
      highHCI = quantile(value, 0.945),
      mean = mean(value),
      median = median(value)
    ) %>%
    collect()
}

# Function: evaluate Topt and its mean, median, and highest prob. intervals across
#           samples
Topt_heat_func <- function(in_df, system_name) {
  out_df <- in_df %>%
    expand_grid(filter(data.Vec, system_ID == system_name)) %>% 
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
    # Filter to maximum value of R0
    group_by(system_ID, sample_num, sigmaH, KH) %>%
    filter(R0 == max(R0)) %>%
    # Get temperature at which R0 is maximized
    rename(Topt = Temperature) %>%
    dplyr::select(system_ID, sample_num, Model, sigmaH, KH, Topt) %>%
    ungroup() %>%
    pivot_longer(cols = Topt, names_to = "variable", values_to = "value") %>%
    group_by(system_ID, Model, sigmaH, KH, variable) %>%
    # partition(cluster) %>%
    summarise(
      lowHCI = quantile(value, 0.055),
      highHCI = quantile(value, 0.945),
      mean = mean(value),
      median = median(value),
      .groups = "keep"
    ) #%>%
  # collect()
}

# ALTERNATE Function: evaluate Topt ONLY WHEN R0>1 and its mean, median, and highest prob. intervals across
#           samples
Topt_heat_func_restricted <- function(in_df, system_name) {
  out_df <- in_df %>%
    expand_grid(filter(data.Vec, system_ID == system_name)) %>% 
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
    # Filter to maximum value of R0
    group_by(system_ID, sample_num, sigmaH, KH) %>%
    filter(R0 == max(R0)) %>%
    filter(R0 > 1) %>% 
    # Get temperature at which R0 is maximized
    rename(Topt = Temperature) %>%
    dplyr::select(system_ID, sample_num, Model, sigmaH, KH, Topt) %>%
    ungroup() %>%
    pivot_longer(cols = Topt, names_to = "variable", values_to = "value") %>%
    group_by(system_ID, Model, sigmaH, KH, variable) %>%
    # partition(cluster) %>%
    summarise(
      lowHCI = quantile(value, 0.055),
      highHCI = quantile(value, 0.945),
      mean = mean(value),
      median = median(value),
      .groups = "keep"
    ) #%>%
  # collect()
}

# Function: evaluate CTmin, CTmax, CTwidth and their means, medians, and highest
#           prob. intervals across samples
CT_heat_func <- function(in_df, system_name) {
  out_df <- in_df %>%
    expand_grid(filter(data.Vec, system_ID == system_name)) %>%
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
    # Filter to maximum value of R0
    group_by(system_ID, sample_num, sigmaH, KH) %>%
    filter(R0 > 1) 
  
  # If R0 < 1 across all temperatures, report the following:
  #  CTmin = Inf, CTmax = -Inf, CTwidth = 0
  if (dim(out_df)[1] == 0) {
    out_df <- expand_grid(select(in_df, sigmaH, KH, Model),
                          system_ID = system_name,
                          variable = c("CTmax", "CTmin", "CTwidth")) %>%
      distinct() %>%
      mutate(mean = case_when(
        variable == "CTmax" ~ -Inf,
        variable == "CTmin" ~ Inf,
        variable == "CTwidth" ~ 0
      )) %>%
      mutate(median = mean) %>%
      mutate(highHCI = mean) %>%
      mutate(lowHCI = mean)
  } else {
    out_df <- out_df %>%
      # Get lowest temperature at which R0 exceeds one
      mutate(CTmin = min(Temperature)) %>%
      # Get highest temperature at which R0 exceeds one
      mutate(CTmax = max(Temperature)) %>%
      # Get width of critical thermal interval
      mutate(CTwidth = ifelse(is.finite(CTmax), 
                              CTmax - CTmin,
                              0))%>%
      ungroup() %>%
      dplyr::select(system_ID, Model, sigmaH, KH, CTmin, CTmax, CTwidth) %>%
      pivot_longer(cols = c(CTwidth, CTmin, CTmax), names_to = "variable", values_to = "value") %>% 
      group_by(system_ID, Model, sigmaH, KH, variable) %>%
      # partition(cluster) %>%
      summarise(
        lowHCI = quantile(value, 0.055),
        highHCI = quantile(value, 0.945),
        mean = mean(value),
        median = median(value),
        .groups = "keep"
      ) #%>%
    # collect()
  }
  return(out_df)
}

# 2) Set up data frames ---------------------------------------------------

###* Mean traits data frame ----

# Get mean values of uncertain traits

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

# 3) Topt local, global sensitivity, uncertainty --------------------------

###* Density of Topt ----

data.Topt <- data.Host %>%
  filter(sigmaH %in% c(100, Inf)) %>%
  filter(KH %in% 10^seq(-2,4))

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
  # Filter to maximum value of R0
  filter(R0>0) %>% 
  select(system_ID, sample_num, Model, sigmaH, KH, Temperature, R0) %>% 
  group_by(system_ID, sample_num, sigmaH, KH) %>%
  filter(R0 == max(R0)) %>%
  distinct() %>% 
  # Get temperature at which R0 is maximized
  rename(Topt = Temperature) %>%
  dplyr::select(system_ID, sample_num, Model, sigmaH, KH, Topt, R0) #%>%
  #distinct() 

Topt.density.df <- full_join(Topt.density.df,
                           expand_grid(data.CT, data.Vec) %>% 
                             select(c(Model, system_ID, sample_num, sigmaH,KH)) %>% 
                             distinct())


plot.Topt.density <- Topt.density.df %>% 
  # filter(sigmaH %in% Inf) %>%
  # filter(KH < 1E3) %>%
  ggplot(aes(x = Topt)) +
  geom_density(aes(color = as.factor(KH), linetype = as.factor(sigmaH)),
               lwd = 1, adjust = 1) +
  # geom_density(data = filter(Topt.density.df, sigmaH == Inf),
  #              linetype = 1,
  #              lwd = 2) +
  # color: scaled log10, color-blind friendly
  scale_color_viridis_d(
    name = "Vertebrate host population\ndensity (ind/ha)",
    breaks = 10^seq(-2, 5),
    labels = unname(c(0.01, 0.1, 1, 10, 100, TeX("$10^3$"), TeX("$10^4$"), TeX("$10^5$"))),
    option = "plasma"
  ) +
  facet_wrap(~system_ID, scales = "free") +
  theme_minimal_grid()

plot.Topt.density

###* Global sensitivity (following Johnson 2015) ----
# 


# 4) CTmin local, global sensitivity, uncertainty -------------------------

###* Density of CTmin,max,width ----

# Initialize R0 data frame
CT.density.df <- init.df

gc()

data.CT <- data.Host %>%
  filter(sigmaH %in% c(100, Inf)) %>%
  filter(KH %in% 10^seq(-2,4))

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
  # Filter to maximum value of R0
  group_by(system_ID, sample_num, sigmaH, KH) %>%
  filter(R0 > 1) %>%
  # Get lowest temperature at which R0 exceeds one
  mutate(CTmin = min(Temperature)) %>%
  # Get highest temperature at which R0 exceeds one
  mutate(CTmax = max(Temperature)) %>%
  # Get width of critical thermal interval
  mutate(CTwidth = CTmax - CTmin) %>%
  dplyr::select(system_ID, sample_num, Model, sigmaH, KH, CTmin, CTmax, CTwidth) %>% 
  distinct()

# Add back values removed from filtering R0>1 above
CT.density.df <- full_join(CT.density.df,
                  expand_grid(data.CT, data.Vec) %>% 
                    select(c(Model, system_ID, sample_num, sigmaH,KH)) %>% 
                    distinct()) %>% 
  replace_na(list(CTmin = Inf,
             CTmax = -Inf,
             CTwidth = 0))


# 5) CTmax local, global sensitivity, uncertainty -------------------------


# 6) CTwidth local, global sensitivity, uncertainty -----------------------



# *) Diagnostics & visualizations -----------------------------------------


###* Density plots of CTmin, Topt, CTmax ----

# Combine data frames
all.density.df <- right_join(Topt.density.df, CT.density.df) %>% 
  select(-CTwidth) %>% 
  pivot_longer(cols = c(CTmin, Topt, CTmax), names_to = "variable", values_to = "value") %>% 
  mutate(system_ID = case_when(
    system_ID == "Anopheles gambiae / Plasmodium falciparum" ~ "An. gamb. / P. falciparum",
    system_ID == "Aedes aegypti / DENV" ~ "Ae. aegypti / DENV",
    system_ID == "Aedes albopictus / DENV" ~ "Ae. albopictus / DENV",
    system_ID == "Aedes aegypti / ZIKV" ~ "Ae. aegypti / ZIKV",
    system_ID == "Culex quinquefasciatus / WNV" ~ "Cx. quin. / WNV"
  )) %>% 
  mutate(variable = case_when(
    variable == "CTmin" ~ "Critical thermal minimum",
    variable == "Topt" ~ "Thermal optimum",
    variable == "CTmax" ~ "Critical thermal maximum"
  ))

all.density.df$variable <- factor(all.density.df$variable, 
                                  levels = c("Critical thermal minimum", "Thermal optimum","Critical thermal maximum"))

# Dynamic: sigmaH = 100
plot.densities_dynamic <- all.density.df  %>% 
  filter(sigmaH %in% 100) %>%
  # filter(KH < 1E3) %>%
  ggplot(aes(x = value)) +
  geom_density(aes(color = as.factor(KH)),
               lwd = 1.5,
               adjust = 2) +
  scale_color_viridis_d(
    name = "Vertebrate host population\ndensity (ind/ha)",
    breaks = 10^seq(-2, 4),
    labels = unname(c(0.01, 0.1, 1, 10, 100, TeX("$10^3$"), TeX("$10^4$"))),
    limits = 10^seq(-2, 4),
    option = "plasma"
  ) +
  facet_grid(rows = vars(system_ID), cols = vars(variable), scales = "free") +
  theme_minimal_grid() 

plot.densities_dynamic

ggsave("figures/results/thermal_densities_dynamic.svg",
       plot = plot.densities_dynamic,
       width = 16, height = 9)
# Static: sigmaH = Inf
plot.densities_RM <- all.density.df  %>% 
  filter(sigmaH %in% Inf) %>%
  # filter(KH < 1E3) %>%
  ggplot(aes(x = value)) +
  geom_density(aes(color = as.factor(KH)),
               lwd = 1.5,
               adjust = 2) +
  scale_color_viridis_d(
    name = "Vertebrate host population\ndensity (ind/ha)",
    breaks = 10^seq(-2, 4),
    labels = unname(c(0.01, 0.1, 1, 10, 100, TeX("$10^3$"), TeX("$10^4$"))),
    limits = 10^seq(-2, 4),
    option = "plasma"
  ) +
  facet_grid(rows = vars(system_ID), cols = vars(variable), scales = "free") +
  theme_minimal_grid() 

plot.densities_RM

ggsave("figures/results/thermal_densities_RM.svg",
       plot = plot.densities_RM,
       width = 16, height = 9)

###* Density plot of CTwidth ----
plot.CTwidth.density <- CT.density.df %>% 
filter(CTwidth > 0) %>% 
  ggplot(aes(x = CTwidth)) +
  geom_density(aes(color = as.factor(KH), linetype = as.factor(sigmaH)),
               lwd = 1, adjust = 1) +
  # geom_density(data = filter(Topt.density.df, sigmaH == Inf),
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















