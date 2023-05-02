## Title: Generate figures for "Thermal optima" manuscript ########################
##
## Project: Global zoonoses - spillover of mosquito-borne pathogens
##
## Purpose: Create figures used in the manuscript
##
## Contents: 0) Load in necessary packages, functions, settings
##           1) Build vector trait data frame
##           2) Build all traits and model outputs data frame
##           3) Build thermal characteristics data frame
##
##
## Settings:  Options for reducing the resolution of variables for memory
##            allocation and figure plotting purposes
##
## Inputs:  data - results/AllThermChar_thin.rds
##                 data/clean/AllOutputs_thin.rds
##
## Outputs: (in ./figures)
##
## Written and maintained by: Kyle Dahlin, kydahlin@gmail.com
## Initialized March 2023
# _______________________________________________________________________________


# 0) Load libraries and data sets. Define helper functions ----------------

# Load libraries
library(tidyverse)
library(latex2exp) # nec
library(viridis) # nec
library(cowplot) # nec
library(modeest)
library(MetBrewer) #nec
library(svglite) #nec

# Load in data sets
data.CT <- read_rds("results/CT_vals.rds")
data.R0 <- read_rds("results/R0_vals.rds")
data.Topt <- read_rds("results/Topt_vals.rds")
data.Topt.restrict <- read_rds("results/Topt_alt_vals.rds")

# Helper function to place legends in empty facets of plot grids
# Code by: Artem Sokolov, found here: https://stackoverflow.com/questions/54438495/shift-legend-into-empty-facets-of-a-faceted-plot-in-ggplot2
shift_legend <- function(p) {
  pnls <- cowplot::plot_to_gtable(p) %>%
    gtable::gtable_filter("panel") %>%
    with(setNames(grobs, layout$name)) %>%
    purrr::keep(~ identical(.x, zeroGrob()))
  
  if (length(pnls) == 0) {return(p)} #stop("No empty facets in the plot")
  
  lemon::reposition_legend(p, "center", panel = names(pnls))
}

# Figure 3: Trait TPCs ----------------------------------------------------



# Figure 4: R0 TPCs -------------------------------------------------------
# R0 as a function of temperature, seperate curves across values of KH and sigmaH

# Plot R0 vs. temperature curves. One plot for each system
# One curve for Ross-Macdonald, a few for Chitnis

# x-axes = temperature
# y-axes = R0 (normalized so that max(R0) = 1)
# plots = one for each mosquito+pathogen combination, "system_ID"
# curve 1 = Ross-Macdonald R0 (shape doesn't depend on KH)
# curves  = Chitnis dynamic R0,
#   - colour = value of KH. place it on a color scale in the bottom right corner
newKH_vals <- unique(data.R0$KH)[seq(1, length(unique(data.R0$KH)), length.out = 21)]

data.R0_norm <- data.R0 %>% 
  filter(KH %in% newKH_vals) %>% 
  filter(variable == "R0") %>%
  select(-variable) %>%
  # normalize variables to have a maximum of 1
  group_by(system_ID, Model, sigmaH, KH) %>% 
  mutate(mean_norm = mean / max(mean)) %>% 
  mutate(median_norm = median / max(median)) %>%
  mutate(highHCI_norm = highHCI / max(median)) %>% 
  mutate(lowHCI_norm = lowHCI / max(median)) %>% 
  arrange(system_ID, sigmaH, KH, Temperature, mean_norm)

R0_plot <- ggplot(mapping = aes(x = Temperature, group = KH)) +
  # path of mean, normalized R0 as a function of temperature (finite sigmaH):
  geom_path(
    data = filter(data.R0_norm, sigmaH == 100),
    aes(y = median_norm, colour = KH, linetype = "finite"),
    lwd = 1
  ) +
  # path of mean, normalized R0 as a function of temperature (infinite sigmaH):
  geom_path(
    data = filter(data.R0_norm, sigmaH == Inf),
    aes(y = median_norm, linetype = "infinite", colour = "black"),
    colour = "black",
    lwd = 1.5
  ) +
  # # 89% HCI of R0 TPC curves
  # geom_ribbon(
  #   data = filter(data.R0_norm, sigmaH == 100),
  #   aes(ymin = lowHCI_norm, ymax = highHCI_norm, fill = KH),
  #   alpha = 0.1
  # ) +
  # x-axis:
  scale_x_continuous(
    name = "Temperature (C)",
    expand = c(0, 0),
    limits = c(14, 34)
  ) +
  # y-axis:
  scale_y_continuous(
    name = TeX("Normalized $R_0$"),
    expand = c(0, .05)
  ) +
  # color: scaled log10, color-blind friendly
  scale_color_viridis(
    name = "Vertebrate host population\ndensity (ind/ha)",
    trans = "log10",
    breaks = 10^seq(-2, 4),
    labels = unname(c(0.01, 0.1, 1, 10, 100, TeX("$10^3$"), TeX("$10^4$"))),
    option = "plasma"
  ) +
  # fill: scaled log10, color-blind friendly
  scale_fill_viridis(
    name = "Vertebrate host population\ndensity (ind/ha)",
    trans = "log10",
    breaks = 10^seq(-2, 4),
    labels = unname(c(0.01, 0.1, 1, 10, 100, TeX("$10^3$"), TeX("$10^4$"))),
    option = "plasma"
  ) +
  # linetype:
  scale_linetype_manual(
    name = "Biting tolerance\n(bites per host per day)",
    breaks = c("finite", "infinite"),
    labels = unname(c(TeX("100 (Chitnis)"), TeX("$\\infty$ (Ross-Macdonald)"))),
    values = c("solid", "22")
  ) +
  # faceting:
  facet_wrap(~system_ID, scales = "free", nrow = 2, labeller = labeller()) +
  # legend
  guides(
    color = guide_colourbar(
      title.position = "top",
      nrow = 1,
      barwidth = 12,
      show.limits = TRUE,
      direction = "horizontal"
    ),
    linetype = guide_legend(
      title.position = "top",
      nrow = 2,
      keywidth = 4,
      direction = "vertical"
    )
  ) +
  # theme options:
  theme_minimal_vgrid(16) +
  theme(
    legend.box = "vertical",
    legend.margin = margin(),
    legend.title = element_text(size = 10),
    legend.text = element_text(
      size = 10,
      hjust = 0
    ),
    legend.justification = "left",
    axis.title.x = element_text(hjust = 0.5),
    strip.text = ggtext::element_markdown(
      size = 11,
      hjust = 0,
      padding = margin(0, 0, 0, 0),
      margin = margin(1, 1, 1, 1)
    )
  )

# Plot
R0_plot <- shift_legend(R0_plot)

ggsave("figures/results/R0_TPCs_median.svg", R0_plot,
       width = 16, height = 9)

# Figure 5: Topt heatmap --------------------------------------------------
# Mean value of Topt as a function of KH and sigmaH
Topt_heat_df <- read_rds("results/Topt_vals.rds")

Toptalt_df <- Topt_heat_df %>%
  ungroup() %>%
  distinct() %>%
  # Arrange along the plotting variables
  arrange(system_ID, sigmaH, KH, mean) %>% 
  mutate(sigmaH_proxy = sigmaH)

# add many proxy sigmaH points at infinity to make the values easier to see there
## Add a "point at infinity" for sigmaH (NB: to be added to get-analysis-dfs.R)
# max finite value of sigmaH in the original data set
max_sigmaH <- max(filter(Toptalt_df, is.finite(sigmaH))$sigmaH)
# Value of sigmaH to replace 'Inf' with for plotting
infinite_skew <- 10^log10(max_sigmaH*1) # 1.1*max_sigmaH #10^2 + 80
# Used to create a gap dividing finite and infinite sigmaH
infinite_divide <- 10^log10(max_sigmaH*0.75)

temp <- filter(Toptalt_df, is.infinite(sigmaH)) %>%
  mutate(sigmaH_proxy = ifelse(is.finite(sigmaH), sigmaH,
                               mean(c(1.25 * infinite_divide, infinite_skew))
  ))
temp2 <- filter(Toptalt_df, is.infinite(sigmaH)) %>%
  mutate(sigmaH_proxy = ifelse(is.finite(sigmaH), sigmaH,
                               mean(c(1.15 * infinite_divide, infinite_skew))
  ))
temp3 <- filter(Toptalt_df, is.infinite(sigmaH)) %>%
  mutate(sigmaH_proxy = ifelse(is.finite(sigmaH), sigmaH,
                               mean(c(1.05 * infinite_divide, infinite_skew))
  ))
temp4 <- filter(Toptalt_df, is.infinite(sigmaH)) %>%
  mutate(sigmaH_proxy = ifelse(is.finite(sigmaH), sigmaH,
                               mean(c(1.3 * infinite_divide, infinite_skew))
  ))
temp5 <- filter(Toptalt_df, is.infinite(sigmaH)) %>%
  mutate(sigmaH_proxy = ifelse(is.finite(sigmaH), sigmaH,
                               mean(c(1.2 * infinite_divide, infinite_skew))
  ))
temp6 <- filter(Toptalt_df, is.infinite(sigmaH)) %>%
  mutate(sigmaH_proxy = ifelse(is.finite(sigmaH), sigmaH,
                               mean(c(1.1 * infinite_divide, infinite_skew))
  ))
temp7 <- filter(Toptalt_df, is.infinite(sigmaH)) %>%
  mutate(sigmaH_proxy = ifelse(is.finite(sigmaH), sigmaH,
                               mean(c(infinite_divide, infinite_skew))
  ))

Toptalt_df <- add_row(Toptalt_df, temp)
Toptalt_df <- add_row(Toptalt_df, temp2)
Toptalt_df <- add_row(Toptalt_df, temp3)
Toptalt_df <- add_row(Toptalt_df, temp4)
Toptalt_df <- add_row(Toptalt_df, temp5)
Toptalt_df <- add_row(Toptalt_df, temp6)
Toptalt_df <- add_row(Toptalt_df, temp7) %>%
  filter(is.finite(sigmaH_proxy)) %>% 
  arrange(sigmaH_proxy) %>%
  distinct()


# Define sigmaH tick breaks and labels
sigmaH_breaks <- c(10^c(seq(-1, 2, by = 1)), infinite_divide, 10^log10(mean(c(infinite_divide,infinite_skew))))
sigmaH_labels <- c(paste(c(10^c(seq(-1, 2, by = 1)))), TeX("$\\ldots$"), TeX("$\\infty$"))

min_Topt <- floor(min(min(Toptalt_df$mean),min(Toptalt_df$median)))
max_Topt <- ceiling(max(max(Toptalt_df$mean),max(Toptalt_df$median)))
Topt_breaks <- seq(min_Topt, max_Topt, by = 1)

# Plot mean of Topt heat map
Toptalt_plots <- Toptalt_df %>%
  ## Set up plot ##
  # x = biting tolerance (with a point at infinity),
  # y = vertebrate host population density,
  # color = transmission thermal optimum
  ggplot(aes(x = KH, y = sigmaH_proxy, z = median)) +
  # Topt for FINITE sigmaH:
  geom_contour_filled(
    data = filter(Toptalt_df, sigmaH_proxy < infinite_divide**0.99),
    breaks = Topt_breaks
  ) +
  # Topt for INFINITE sigmaH:
  # take mean of mean to hide variation from artifacts of uncertainty 
  geom_contour_filled(
    data = filter(Toptalt_df, sigmaH_proxy > infinite_divide**1.01) %>% 
      group_by(system_ID) %>% 
      mutate(mean = mean(mean), .groups = "keep"),
    size = 1,
    breaks = Topt_breaks
  ) +
  # # R0 = 1 contour: cheat and use where CTwidth = 0 (i.e. where R0<1 at any temperature)
  geom_contour(data = filter(data.CT, variable == "CTwidth"),
               aes(y = sigmaH, z = mean, colour = "R0 = 1"), 
               breaks = min(filter(data.CT, variable == "CTwidth", mean > 0)$mean),
               size = 1) +
  # geom_contour(data = data.R0_norm,
  #              aes(y = sigmaH, z = mean, colour = "R0=1"), breaks = c(1), size = 1) +
  # x-axis: log10-scale with no buffer space
  scale_x_log10(
    name = TeX("Vertebrate host population density (ind/ha)"),
    expand = c(0, 0),
    limits = c(0.05, 1E4),
    breaks = 10^seq(-1, 4, 1),
    labels = TeX(c("$0.1$", "$1$", "$10", "$100$", "$10^3$", "$10^4$"))
  ) +
  # y-axis: log10-scale y axis, with no buffer space
  scale_y_log10(
    name = TeX("Biting tolerance (bites per host per day)"),
    expand = c(0, 0),
    breaks = sigmaH_breaks,
    labels = sigmaH_labels
  ) +
  
  # fill:
  scale_fill_manual(
    name = "Transmission\nthermal optimum (°C)",
    values = met.brewer("Johnson",
                        n = length(Topt_breaks),
                        direction = -1
    )
  ) +
  # color:
  scale_colour_manual(
    name = "",
    labels = c("R0=1" = unname(TeX("$R_0=1$"))),
    values = c("black"),
    guide = guide_legend(
      order = 1,
      keyheight = unit(1, "cm"),
      override.aes = list(size = 6)
    )
  ) +
  # faceting:
  facet_wrap(~system_ID, scales = "free", nrow = 2, labeller = labeller()) +
  # legend:
  guides(fill = guide_coloursteps(
    title.position = "top", title.hjust = 0,
    barheight = unit(5, "cm"),
    show.limits = TRUE
  )) +
  # theme options:
  theme_minimal(16) +
  theme(
    legend.box = "vertical",
    legend.position = "right",
    legend.margin = margin(0, 0, 0, 0),
    legend.title = element_text(size = 10),
    legend.direction = "vertical",
    legend.justification = "left",
    axis.title.x = element_text(hjust = 0.5),
    strip.text = ggtext::element_markdown(
      size = 11,
      hjust = 0,
      padding = margin(0, 0, 0, 0),
      margin = margin(1, 1, 1, 1)
    )
  )

# Plot
Toptalt_plots <- shift_legend(Toptalt_plots)

ggsave("figures/results/Topt_median.svg", Toptalt_plots,
       width = 16, height = 9)


# Figure 6: CTwidth heatmap -----------------------------------------------
# Shows how the parasite thermal niche (temperatures between CT_min and CT_max) 
# depends on vertebrate host availability (density and biting tolerance) in each
# of the focal systems.

#   x    = biting tolerance (with a point at infinity)
#   y    = vertebrate host population density
# colour = width of parasite thermal niche (CTmax - CTmin)
# Facets = systems

CTwidth_df <- read_rds("results/CT_vals.rds") %>%
  filter(variable == "CTwidth") %>% 
  mutate(sigmaH_proxy = sigmaH) %>% 
  # Arrange along the plotting variables
  arrange(system_ID, KH, sigmaH_proxy, mean) %>% 
  ungroup()

# add many third proxy sigmaH points at infinity to make the values easier to see there
temp <- filter(CTwidth_df, is.infinite(sigmaH)) %>%
  mutate(sigmaH_proxy = ifelse(is.finite(sigmaH), sigmaH,
                               mean(c(1.25 * infinite_divide, infinite_skew))
  ))
temp2 <- filter(CTwidth_df, is.infinite(sigmaH)) %>%
  mutate(sigmaH_proxy = ifelse(is.finite(sigmaH), sigmaH,
                               mean(c(1.15 * infinite_divide, infinite_skew))
  ))
temp3 <- filter(CTwidth_df, is.infinite(sigmaH)) %>%
  mutate(sigmaH_proxy = ifelse(is.finite(sigmaH), sigmaH,
                               mean(c(1.05 * infinite_divide, infinite_skew))
  ))
temp4 <- filter(CTwidth_df, is.infinite(sigmaH)) %>%
  mutate(sigmaH_proxy = ifelse(is.finite(sigmaH), sigmaH,
                               mean(c(1.3 * infinite_divide, infinite_skew))
  ))
temp5 <- filter(CTwidth_df, is.infinite(sigmaH)) %>%
  mutate(sigmaH_proxy = ifelse(is.finite(sigmaH), sigmaH,
                               mean(c(1.2 * infinite_divide, infinite_skew))
  ))
temp6 <- filter(CTwidth_df, is.infinite(sigmaH)) %>%
  mutate(sigmaH_proxy = ifelse(is.finite(sigmaH), sigmaH,
                               mean(c(1.1 * infinite_divide, infinite_skew))
  ))
temp7 <- filter(CTwidth_df, is.infinite(sigmaH)) %>%
  mutate(sigmaH_proxy = ifelse(is.finite(sigmaH), sigmaH,
                               mean(c(infinite_divide, infinite_skew))
  ))

CTwidth_df <- add_row(CTwidth_df, temp)
CTwidth_df <- add_row(CTwidth_df, temp2)
CTwidth_df <- add_row(CTwidth_df, temp3)
CTwidth_df <- add_row(CTwidth_df, temp4)
CTwidth_df <- add_row(CTwidth_df, temp5)
CTwidth_df <- add_row(CTwidth_df, temp6)
CTwidth_df <- add_row(CTwidth_df, temp7) %>%
  arrange(sigmaH_proxy) %>%
  distinct()

CTWidth_heat_plots <- CTwidth_df %>%
  ## Set up plot ##
  # x = biting tolerance (with a point at infinity),
  # y = vertebrate host population density,
  # color = width of parasite thermal mean_CTwidth
  ggplot(aes(x = KH, y = sigmaH_proxy, z = mean)) +
  # CTwidth for FINITE sigmaH:
  geom_contour_filled(data = filter(CTwidth_df, sigmaH_proxy < infinite_divide**0.99)) +
  # CTwidth for INFINITE sigmaH:
  geom_contour_filled(data = filter(CTwidth_df, sigmaH_proxy > infinite_divide**1.01)) +
  # x-axis: log10-scale with no buffer space
  scale_x_log10(
    name = TeX("Vertebrate host population density (ind/ha)"),
    expand = c(0, 0),
    limits = c(0.05, 1E4),
    breaks = 10^seq(-1, 4, 1),
    labels = TeX(c("$0.1$", "$1$", "$10", "$100$", "$10^3$", "$10^4$"))
  ) +
  # y-axis: log10-scale y axis, with no buffer space
  scale_y_log10(
    name = TeX("Biting tolerance (bites per host per day)"),
    expand = c(0, 0),
    breaks = sigmaH_breaks,
    labels = sigmaH_labels
  ) +
  # color:
  scale_fill_viridis_d(
    name = "Mean parasite\nthermal tolerance\nrange width (°C)",
    option = "plasma"
  ) +
  # faceting:
  facet_wrap(~system_ID, scales = "free", nrow = 2, labeller = labeller()) +
  # legend:
  guides(fill = guide_coloursteps(
    title.position = "top", title.hjust = 0.5,
    barwidth = unit(5, "cm"),
    show.limits = TRUE
  )) +
  # theme options:
  theme_minimal(16) +
  theme(
    legend.box = "horizontal",
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0),
    legend.title = element_text(size = 10),
    legend.direction = "horizontal",
    legend.justification = "left",
    axis.title.x = element_text(hjust = 0.5),
    strip.text = ggtext::element_markdown(
      size = 11,
      hjust = 0,
      padding = margin(0, 0, 0, 0),
      margin = margin(1, 1, 1, 1)
    )
  )

# Plot
CTWidth_heat_plots <- shift_legend(CTWidth_heat_plots)

ggsave("figures/results/CTwidth_mean.svg", CTWidth_heat_plots,
       width = 16, height = 9)

# *Supplementary Figures* ----


# Figure S1: mean R0 as function of temperature and biting tolerance --------


# Figure S2: Topt as a function of KH -------------------------------------
Topt_df <- read_rds("results/Topt_vals.rds") %>% 
  pivot_wider(id_cols = c("system_ID", "Model", "sigmaH", "KH"), 
              names_from = variable,
              values_from = c("mean", "median", "lowHCI", "highHCI"))  %>% 
  right_join(read_rds("results/CT_vals.rds") %>% 
               filter(variable == "CTwidth") %>% 
               pivot_wider(id_cols = c("system_ID", "Model", "sigmaH", "KH"), 
                           names_from = "variable",
                           values_from = c("mean", "median", "lowHCI", "highHCI"))) %>% 
  filter(sigmaH %in% c(1, 10, 100, Inf)) %>% 
  filter(mean_CTwidth>0) %>% 
  select(-c(mean_CTwidth, median_CTwidth, lowHCI_CTwidth, highHCI_CTwidth)) %>%
  rename(mean = mean_Topt, median = median_Topt, lowHCI = lowHCI_Topt, 
         highHCI = highHCI_Topt) %>% 
  # filter(KH > 0.01) %>%
  # filter(sigmaH %in% 10^seq(0, 2) | is.infinite(sigmaH)) %>%
  group_by(KH) %>%
  # Arrange along the plotting variables
  arrange(system_ID, KH, mean) %>%
  arrange(sigmaH) %>%
  ungroup()

# PLOTTING
Topt_plot <- Topt_df %>%
  ## Set up plot ##
  # color = sigmaH
  ggplot(aes(x = KH, colour = as.factor(sigmaH))) +
  # Topt curves:
  geom_path(aes(y = median), lwd = 1) +
  # Add dotted lines showing 89% HCI limits
  geom_path(aes(y = lowHCI), linetype = "dashed") +
  geom_path(aes(y = highHCI), linetype = "dashed") +
  # x-axis: log10 scale, no buffer space
  scale_x_log10(
    name = TeX("Vertebrate host population density (ind/ha)"),
    expand = c(0, 0),
    limits = c(0.01, 2E4),
    breaks = 10^seq(-2, 4),
    labels = TeX(c(
      "0.01", "0.1", "$1$", "$10$", "$100$", "$10^3$", "$10^4$"
    ))
  ) +
  # y-axis
  scale_y_continuous(
    name = expression("Thermal optimum "(degree * C)),
    expand = c(0.05, 0.05),
    # limits = c(21,32), # mean: 21-31, median: 21-32
    # breaks = seq(21,32, by = 1)
  ) +
  # color:
  scale_colour_manual(
    name = "Vertebrate host biting tolerance\n(bites per host per day)",
    values = c(met.brewer("VanGogh3", 4, direction = 1), "black"),
    # values = c(brewer.pal(9, "YlOrRd")[c(3, 5, 7)], "black"),
    breaks = c(1, 10, 100, Inf),
    labels = unname(c("1 (Chitnis)", "10 (Chitnis)", "100 (Chitnis)", TeX("$\\infty$ (Ross-Macdonald)")))
  ) +
  scale_fill_manual(
    name = "Vertebrate host biting tolerance\n(bites per host per day)",
    values = c(met.brewer("VanGogh3", 4, direction = 1), "black"),
    # values = c(brewer.pal(9, "YlOrRd")[c(3, 5, 7)], "black"),
    breaks = c(1, 10, 100, Inf),
    labels = unname(c("1 (Chitnis)", "10 (Chitnis)", "100 (Chitnis)", TeX("$\\infty$ (Ross-Macdonald)")))
  ) +
  # faceting:
  facet_wrap(~system_ID, scales = "free", nrow = 2, labeller = labeller()) +
  # theme options:
  theme_minimal(16) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(),
    legend.direction = "vertical",
    legend.title = element_text(size = 10),
    legend.text = element_text(
      size = 10,
      hjust = 0
    ),
    legend.justification = "left",
    axis.title.x = element_text(hjust = 0.5),
    strip.text = ggtext::element_markdown(
      size = 11,
      hjust = 0,
      padding = margin(0, 0, 0, 0),
      margin = margin(1, 1, 1, 1)
    )
  )


# Plot
Topt_plot <- shift_legend(Topt_plot)

ggsave("figures/results/Topt_KH_median.svg", Topt_plot,
       width = 16, height = 9)

# Figure S3: Topt as a function of sigmaH ---------------------------------
# Mean Topt as a function of sigmaH, across discrete values of KH

Topt_plot <- read_rds("results/Topt_vals.rds") %>% 
  # filter(sigmaH %in% c(1, 10, 100, Inf)) %>% 
  # Choose appropriate pop. dens. values and group
  filter(KH %in% 10^seq(0, 3)) %>%
  arrange(system_ID, KH, sigmaH, median) %>% 
  ## Set up plot ##
  # color = sigmaH
  ggplot(aes(x = sigmaH, colour = as.factor(KH))) +
  # Topt curves:
  geom_path(aes(y = mean), lwd = 1) +
  # Dotted lines showing limits of 89% HCI
  geom_path(aes(y = lowHCI),linetype = "dashed") +
  geom_path(aes(y = highHCI),linetype = "dashed") +
  # white lines breaking up finite/infinite sigmaH values:
  geom_vline(
    aes(xintercept = 0.95 * infinite_divide),
    colour = "white",
    lwd = .9,
    show.legend = FALSE
  ) +
  geom_vline(
    aes(xintercept = 1.05 * infinite_divide),
    colour = "white",
    lwd = .9,
    show.legend = FALSE
  ) +
  geom_vline(
    aes(xintercept = 1.15 * infinite_divide),
    colour = "white",
    lwd = .9,
    show.legend = FALSE
  ) +
  # x-axis: Use a log10-scale with no buffer space
  scale_x_log10(
    name = TeX("Vertebrate host biting tolerance (bites per host per day)"),
    expand = c(0, 0),
    limits = c(0.5, 150),
    breaks = sigmaH_breaks,
    labels = sigmaH_labels
  ) +
  # y-axis
  scale_y_continuous(
    name = expression("Thermal optimum "(degree * C)),
    expand = c(0.05, 0.05),
    limits = c(21,32), # mean: 21-31, median: 21-32
    breaks = seq(21,32, by = 1)
  ) +
  # color:
  scale_colour_manual(
    name = "Vertebrate host population\ndensity (ind/ha)",
    values = met.brewer("Hokusai3", 4,
                        # option = "plasma"
    )) +
  scale_fill_manual(
    name = "Vertebrate host population\ndensity (ind/ha)",
    values = met.brewer("Hokusai3", 4,
                        # option = "plasma"
    )) +
  # faceting:
  facet_wrap(~system_ID,
             scales = "free",
             nrow = 2,
             labeller = labeller()  ) +
  # theme options:
  theme_minimal(16) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(),
    legend.direction = "vertical",
    legend.title = element_text(size = 10),
    legend.text = element_text(
      size = 10,
      hjust = 0
    ),
    legend.justification = "left",
    axis.title.x = element_text(hjust = 0.5),
    strip.text = ggtext::element_markdown(
      size = 11,
      hjust = 0,
      padding = margin(0, 0, 0, 0),
      margin = margin(1, 1, 1, 1)
    )
  )


# Plot
Topt_plot <- shift_legend(Topt_plot)

ggsave("figures/results/Topt_sigmaH_mean.svg", Topt_plot,
       width = 16, height = 9)

# Figure S4: CTmin heatmap ------------------------------------------------
# Shows how the parasite thermal niche (temperatures between CT_min and CT_max) 
# depends on vertebrate host availability (density and biting tolerance) in each
# of the focal systems.

#   x    = biting tolerance (with a point at infinity)
#   y    = vertebrate host population density
# colour = CTmin
# Facets = systems

CTmin_df <- read_rds("results/CT_vals.rds") %>%
  filter(variable == "CTmin") %>% 
  filter(is.finite(mean)) %>%
  mutate(sigmaH_proxy = sigmaH) %>% 
  # Arrange along the plotting variables
  arrange(system_ID, KH, sigmaH_proxy, mean) %>% 
  ungroup()

# add many third proxy sigmaH points at infinity to make the values easier to see there
temp <- filter(CTmin_df, is.infinite(sigmaH)) %>%
  mutate(sigmaH_proxy = ifelse(is.finite(sigmaH), sigmaH,
                               mean(c(1.25 * infinite_divide, infinite_skew))
  ))
temp2 <- filter(CTmin_df, is.infinite(sigmaH)) %>%
  mutate(sigmaH_proxy = ifelse(is.finite(sigmaH), sigmaH,
                               mean(c(1.15 * infinite_divide, infinite_skew))
  ))
temp3 <- filter(CTmin_df, is.infinite(sigmaH)) %>%
  mutate(sigmaH_proxy = ifelse(is.finite(sigmaH), sigmaH,
                               mean(c(1.05 * infinite_divide, infinite_skew))
  ))
temp4 <- filter(CTmin_df, is.infinite(sigmaH)) %>%
  mutate(sigmaH_proxy = ifelse(is.finite(sigmaH), sigmaH,
                               mean(c(1.3 * infinite_divide, infinite_skew))
  ))
temp5 <- filter(CTmin_df, is.infinite(sigmaH)) %>%
  mutate(sigmaH_proxy = ifelse(is.finite(sigmaH), sigmaH,
                               mean(c(1.2 * infinite_divide, infinite_skew))
  ))
temp6 <- filter(CTmin_df, is.infinite(sigmaH)) %>%
  mutate(sigmaH_proxy = ifelse(is.finite(sigmaH), sigmaH,
                               mean(c(1.1 * infinite_divide, infinite_skew))
  ))
temp7 <- filter(CTmin_df, is.infinite(sigmaH)) %>%
  mutate(sigmaH_proxy = ifelse(is.finite(sigmaH), sigmaH,
                               mean(c(infinite_divide, infinite_skew))
  ))

CTmin_df <- add_row(CTmin_df, temp)
CTmin_df <- add_row(CTmin_df, temp2)
CTmin_df <- add_row(CTmin_df, temp3)
CTmin_df <- add_row(CTmin_df, temp4)
CTmin_df <- add_row(CTmin_df, temp5)
CTmin_df <- add_row(CTmin_df, temp6)
CTmin_df <- add_row(CTmin_df, temp7) %>%
  arrange(sigmaH_proxy) %>%
  distinct()

CTmin_heat_plots <- CTmin_df %>%
  ## Set up plot ##
  # x = biting tolerance (with a point at infinity),
  # y = vertebrate host population density,
  # color = width of parasite thermal mean_CTwidth
  ggplot(aes(x = KH, y = sigmaH_proxy, z = mean)) +
  
  # CTmin for FINITE sigmaH:
  geom_contour_filled(
    data = filter(CTmin_df, sigmaH_proxy < 100)
  ) +
  
  # CTmin for INFINITE sigmaH:
  geom_contour_filled(
    data = filter(CTmin_df, is.infinite(sigmaH))
  ) +
  
  # x-axis: log10-scale with no buffer space
  scale_x_log10(
    name = TeX("Vertebrate host population density (ind/ha)"),
    expand = c(0, 0),
    limits = c(0.05, 1E4),
    breaks = 10^seq(-1, 4, 1),
    labels = TeX(c("$0.1$", "$1$", "$10", "$100$", "$10^3$", "$10^4$"))
  ) +
  # y-axis: log10-scale y axis, with no buffer space
  scale_y_log10(
    name = TeX("Biting tolerance (bites per host per day)"),
    expand = c(0, 0),
    breaks = sigmaH_breaks,
    labels = sigmaH_labels
  ) +
  # color:
  scale_fill_viridis_d(
    name = "Mean critical thermal minimum (°C)",
    option = "plasma"
  ) +
  # faceting:
  facet_wrap(~system_ID, scales = "free", nrow = 2, labeller = labeller()) +
  # legend:
  guides(fill = guide_coloursteps(
    title.position = "top", title.hjust = 0.5,
    barwidth = unit(10, "cm"),
    show.limits = TRUE
  )) +
  # theme options:
  theme_minimal(16) +
  theme(
    legend.box = "horizontal",
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0),
    legend.title = element_text(size = 10),
    legend.direction = "horizontal",
    legend.justification = "left",
    axis.title.x = element_text(hjust = 0.5),
    strip.text = ggtext::element_markdown(
      size = 11,
      hjust = 0,
      padding = margin(0, 0, 0, 0),
      margin = margin(1, 1, 1, 1)
    )
  )

# Plot
CTmin_heat_plots <- shift_legend(CTmin_heat_plots)

ggsave("figures/results/CTmin_mean.svg", CTmin_heat_plots,
       width = 16, height = 9)

# Figure S5: CTmax heatmap ------------------------------------------------
CTmax_df <- read_rds("results/CT_vals.rds") %>%
  filter(variable == "CTmax") %>% 
  filter(is.finite(mean)) %>%
  mutate(sigmaH_proxy = sigmaH) %>% 
  # Arrange along the plotting variables
  arrange(system_ID, KH, sigmaH_proxy, mean) %>% 
  ungroup()

# add many third proxy sigmaH points at infinity to make the values easier to see there
temp <- filter(CTmax_df, is.infinite(sigmaH)) %>%
  mutate(sigmaH_proxy = ifelse(is.finite(sigmaH), sigmaH,
                               mean(c(1.25 * infinite_divide, infinite_skew))
  ))
temp2 <- filter(CTmax_df, is.infinite(sigmaH)) %>%
  mutate(sigmaH_proxy = ifelse(is.finite(sigmaH), sigmaH,
                               mean(c(1.15 * infinite_divide, infinite_skew))
  ))
temp3 <- filter(CTmax_df, is.infinite(sigmaH)) %>%
  mutate(sigmaH_proxy = ifelse(is.finite(sigmaH), sigmaH,
                               mean(c(1.05 * infinite_divide, infinite_skew))
  ))
temp4 <- filter(CTmax_df, is.infinite(sigmaH)) %>%
  mutate(sigmaH_proxy = ifelse(is.finite(sigmaH), sigmaH,
                               mean(c(1.3 * infinite_divide, infinite_skew))
  ))
temp5 <- filter(CTmax_df, is.infinite(sigmaH)) %>%
  mutate(sigmaH_proxy = ifelse(is.finite(sigmaH), sigmaH,
                               mean(c(1.2 * infinite_divide, infinite_skew))
  ))
temp6 <- filter(CTmax_df, is.infinite(sigmaH)) %>%
  mutate(sigmaH_proxy = ifelse(is.finite(sigmaH), sigmaH,
                               mean(c(1.1 * infinite_divide, infinite_skew))
  ))
temp7 <- filter(CTmax_df, is.infinite(sigmaH)) %>%
  mutate(sigmaH_proxy = ifelse(is.finite(sigmaH), sigmaH,
                               mean(c(infinite_divide, infinite_skew))
  ))

CTmax_df <- add_row(CTmax_df, temp)
CTmax_df <- add_row(CTmax_df, temp2)
CTmax_df <- add_row(CTmax_df, temp3)
CTmax_df <- add_row(CTmax_df, temp4)
CTmax_df <- add_row(CTmax_df, temp5)
CTmax_df <- add_row(CTmax_df, temp6)
CTmax_df <- add_row(CTmax_df, temp7) %>%
  arrange(sigmaH_proxy) %>%
  distinct()

min_CTmax <- floor(min(min(CTmax_df$mean),min(CTmax_df$median)))
max_CTmax <- ceiling(max(max(CTmax_df$mean),max(CTmax_df$median)))
CTmax_breaks <- seq(min_CTmax, max_CTmax, by = 1)

CTmax_heat_plots <- CTmax_df %>%
  ## Set up plot ##
  # x = biting tolerance (with a point at infinity),
  # y = vertebrate host population density,
  # color = CTmax
  ggplot(aes(x = KH, y = sigmaH_proxy, z = mean)) +
  # CTmax for FINITE sigmaH:
  geom_contour_filled(
    data = filter(CTmax_df, sigmaH_proxy < 100)
  ) +
  
  # CTmax for INFINITE sigmaH:
  geom_contour_filled(
    data = filter(CTmax_df, is.infinite(sigmaH))
  ) +
  # x-axis: log10-scale with no buffer space
  scale_x_log10(
    name = TeX("Vertebrate host population density (ind/ha)"),
    expand = c(0, 0),
    limits = c(0.05, 1E4),
    breaks = 10^seq(-1, 4, 1),
    labels = TeX(c("$0.1$", "$1$", "$10", "$100$", "$10^3$", "$10^4$"))
  ) +
  # y-axis: log10-scale y axis, with no buffer space
  scale_y_log10(
    name = TeX("Biting tolerance (bites per host per day)"),
    expand = c(0, 0),
    breaks = sigmaH_breaks,
    labels = sigmaH_labels
  ) +
  # color:
  scale_fill_viridis_d(
    name = "Mean critical thermal maximum (°C)",
    option = "plasma"
  ) +
  # faceting:
  facet_wrap(~system_ID, scales = "free", nrow = 2, labeller = labeller()) +
  # legend:
  guides(fill = guide_coloursteps(
    title.position = "top", title.hjust = 0.5,
    barwidth = unit(10, "cm"),
    show.limits = TRUE
  )) +
  # theme options:
  theme_minimal(16) +
  theme(
    legend.box = "horizontal",
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0),
    legend.title = element_text(size = 10),
    legend.direction = "horizontal",
    legend.justification = "left",
    axis.title.x = element_text(hjust = 0.5),
    strip.text = ggtext::element_markdown(
      size = 11,
      hjust = 0,
      padding = margin(0, 0, 0, 0),
      margin = margin(1, 1, 1, 1)
    )
  )

# Plot
CTmax_heat_plots <- shift_legend(CTmax_heat_plots)

ggsave("figures/results/CTmax_mean.svg", CTmax_heat_plots,
       width = 16, height = 9)

# Figure S6: CTmin/max as functions of KH ---------------------------------


# Figure S7: CHmin/max as functions of KH ---------------------------------


###* Figure: Topt variance as function of KH and sigmaH ----
Toptstdev_df <- Toptalt_df %>% 
  select(system_ID, sigmaH_proxy, KH, std_Topt, mean_Topt) %>% 
  distinct()

# Plot standard deviation of Topt heat map
Toptstdev_plots <- Toptalt_df %>% 
  select(system_ID, sigmaH_proxy, KH, std_Topt, mean_Topt) %>% 
  distinct() %>% 
  ## Set up plot ##
  # x = biting tolerance (with a point at infinity),
  # y = vertebrate host population density,
  # color = transmission thermal optimum
  ggplot(aes(x = KH, y = sigmaH_proxy, z = std_Topt)) +
  # Topt for FINITE sigmaH:
  geom_contour_filled(
    data = filter(Toptstdev_df, sigmaH_proxy < infinite_divide**0.98)
  ) +
  # Topt for INFINITE sigmaH:
  geom_contour_filled(
    data = filter(Toptstdev_df, sigmaH_proxy > infinite_divide**1.02),
    size = 1
  ) +
  # x-axis: log10-scale with no buffer space
  scale_x_log10(
    name = TeX("Vertebrate host population density (ind/ha)"),
    expand = c(0, 0),
    # limits = c(0.05, 1E4),
    # breaks = 10^seq(-1, 4, 1),
    # labels = TeX(c("$0.1$", "$1$", "$10", "$100$", "$10^3$", "$10^4$"))
  ) +
  # y-axis: log10-scale y axis, with no buffer space
  scale_y_log10(
    name = TeX("Biting tolerance (bites per host per day)"),
    expand = c(0, 0),
    breaks = sigmaH_breaks,
    labels = sigmaH_labels
  ) +
  # fill:
  scale_fill_manual(
    name = "Standard deviation\nof transmission thermal\noptimum (°C)",
    values = met.brewer("Johnson",
                        n = 13,
                        direction = -1)) +
  # color:
  # faceting:
  facet_wrap(~system_ID, scales = "free", nrow = 2, labeller = labeller()) +
  # legend:
  guides(fill = guide_coloursteps(
    title.position = "top", title.hjust = 0,
    barheight = unit(5, "cm"),
    show.limits = TRUE
  )) +
  # theme options:
  theme_minimal(16) +
  theme(
    legend.box = "vertical",
    legend.position = "right",
    legend.margin = margin(0, 0, 0, 0),
    legend.title = element_text(size = 10),
    legend.direction = "vertical",
    legend.justification = "left",
    axis.title.x = element_text(hjust = 0.5),
    strip.text = ggtext::element_markdown(
      size = 11,
      hjust = 0,
      padding = margin(0, 0, 0, 0),
      margin = margin(1, 1, 1, 1)
    )
  )

# Plot
Toptstdev_plots <- shift_legend(Toptstdev_plots)
ggsave("figures/results/Topt_stdev.svg", Toptstdev_plots,
       width = 16, height = 9)

