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
## Inputs:  data - data/clean/trait_transforms.rds
##
## Outputs: (in ./results)
##          1) VectorTraits.csv - vector trait data frame
##          2) AllOutputs.csv - all traits and model outputs data frame
##          3) ThermalCharacteristics.csv - thermal characteristics data frame
##
## Written and maintained by: Kyle Dahlin, kydahlin@gmail.com
## Initialized March 2023
# _______________________________________________________________________________

# 0) ----

library(tidyverse)
library(latex2exp) # nec
library(viridis) # nec
library(cowplot) # nec
library(modeest)

data.in.thermchars <- read_rds("results/AllThermChar_thin.rds")
data.in.outputs <- read_rds("data/clean/AllOutputs_thin.rds")

# 1) ----
# Helper function to place legends in empty facets of plot grids
# Code by: Artem Sokolov, found here: https://stackoverflow.com/questions/54438495/shift-legend-into-empty-facets-of-a-faceted-plot-in-ggplot2
shift_legend <- function(p) {
  pnls <- cowplot::plot_to_gtable(p) %>%
    gtable::gtable_filter("panel") %>%
    with(setNames(grobs, layout$name)) %>%
    purrr::keep(~ identical(.x, zeroGrob()))
  
  if (length(pnls) == 0) stop("No empty facets in the plot")
  
  lemon::reposition_legend(p, "center", panel = names(pnls))
}

# 2) -----

R0_df <- data.in.outputs %>%
  # To set num. of curves, change "length.out" to be the number of curves you want
  filter(KH %in% unique(KH)[seq(1, length(unique(KH)), length.out = 21)]) %>%
  # Just consider one finite value of sigmaH and the Ross-Macdonald model
  filter(sigmaH == 100 | is.infinite(sigmaH)) %>%
  # Group by curve_ID and system_ID
  group_by(system_ID, Model, KH, sample_num) %>%
  # Normalize R0 for each curve_ID x system_ID combination, so that the maximum is always at one
  mutate(R0 = ifelse(is.nan(R0), 0, R0)) %>%
  mutate(norm_R0 = R0 / max(R0)) %>%
  mutate(norm_R0 = ifelse(is.nan(norm_R0), 0, norm_R0)) %>%
  # Restrict temperatures to where R0 is positive
  # filter(norm_R0 > 0) %>%
  ungroup() %>%
  dplyr::select(system_ID, sample_num, Temperature, Model, sigmaH, KH, norm_R0) %>%
  arrange(system_ID, sample_num, sigmaH, Temperature)

meanR0TPC_df <- R0_df %>%
  group_by(system_ID, Temperature, Model, sigmaH, KH) %>%
  summarise(
    mean_val = mean(norm_R0),
    median_val = median(norm_R0),
    # mode_val = mlv(norm_R0, method = 'mfv'),
    .groups = "keep"
  ) %>%
  arrange(system_ID, sigmaH, KH, Temperature, mean_val, median_val) %>% # , mode_val) %>% 
  distinct()

quantsR0TPC_df <- R0_df %>%
  group_by(system_ID, Temperature, Model, sigmaH, KH) %>%
  mutate(lowHCI_val = quantile(norm_R0, 0.055)) %>%
  mutate(highHCI_val = quantile(norm_R0, 0.945)) %>%
  arrange(system_ID, sigmaH, KH, Temperature, lowHCI_val, highHCI_val) %>%
  dplyr::select(-c("sample_num")) %>% 
  distinct()

###* Figure: R0 vs. temperature with HCIs ----
# Plot R0 vs. temperature curves. One plot for each mosquito+pathogen pair.
# One curve for Ross-Macdonald, a few for Chitnis

# x-axes = temperature
# y-axes = R0 (normalized so that max(R0) = 1)
# plots = one for each mosquito+pathogen combination, "system_ID"
# curve 1 = Ross-Macdonald R0 (shape doesn't depend on KH)
# curves  = Chitnis dynamic R0,
#   - colour = value of KH. place it on a color scale in the bottom right corner

## PLOTTING ##
R0_plot <- ggplot(mapping = aes(x = Temperature,group = KH)) +
  # path of mean, normalized R0 as a function of temperature (finite sigmaH):
  geom_path(
    data = filter(meanR0TPC_df, sigmaH == 100),
    aes(y = mean_val, colour = KH, linetype = "finite")
  ) +
  # path of mean, normalized R0 as a function of temperature (infinite sigmaH):
  geom_path(
    data = filter(meanR0TPC_df, sigmaH == Inf),
    aes(y = mean_val, linetype = "infinite", colour = "black"),
    colour = "black",
    lwd = 1.5
  ) +
  # # 89% HCI of R0 TPC curves
  # geom_ribbon(
  #   data = filter(quantsR0TPC_df, sigmaH == 100),
  #   aes(ymin = lowHCI_val, ymax = highHCI_val, fill = KH),
  #   alpha = 0.1
  # ) +
  # x-axis:
  scale_x_continuous(
    name = "Temperature (C)",
    expand = c(0, 0),
    limits = c(15, 35)
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

ggsave("results/R0_tpcs.png", R0_plot,
       width = 16, height = 9)

###* Figure: Topt curves with HCIs ----

plotTopt_df <- Topt_df %>%
  filter(KH > 0.01) %>%
  filter(sigmaH %in% 10^seq(0, 2) | is.infinite(sigmaH)) %>%
  group_by(KH) %>%
  # Arrange along the plotting variables
  arrange(system_ID, KH, Topt) %>%
  arrange(sigmaH) %>%
  ungroup()

meanTopt_df <- plotTopt_df %>%
  group_by(system_ID, Model, sigmaH, KH) %>%
  summarise(
    mean_val = mean(Topt),
    median_val = median(Topt),
    # mode_val = mlv(norm_R0, method = 'mfv'),
    .groups = "keep"
  ) %>%
  arrange(system_ID, sigmaH, KH, mean_val, median_val)# , mode_val)

quantsTopt_df <- plotTopt_df %>%
  group_by(system_ID, Model, sigmaH, KH) %>%
  mutate(lowHCI_val = quantile(Topt, 0.055)) %>%
  mutate(highHCI_val = quantile(Topt, 0.945)) %>%
  arrange(system_ID, sigmaH, KH, lowHCI_val, highHCI_val) %>%
  dplyr::select(-c("sample_num"))


library(MetBrewer) #nec

Topt_plots <- meanTopt_df %>%
  ## Set up plot ##
  # color = sigmaH
  ggplot(aes(
    x = KH
  )) +
  # Topt curves:
  geom_path(aes(y = mean_val, colour = as.factor(sigmaH)), lwd = 1) +
  # 89% HCI of R0 TPC curves
  geom_ribbon(
    data = quantsTopt_df,
    aes(ymin = lowHCI_val, ymax = highHCI_val, 
        fill = as.factor(sigmaH)),
    alpha = 0.05
  ) +
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
    expand = c(0.05, 0.05)
  ) +
  # color:
  scale_colour_manual(
    name = "Vertebrate host biting tolerance\n(bites per host per day)",
    values = c(met.brewer("VanGogh3", 3, direction = 1), "black"),
    # values = c(brewer.pal(9, "YlOrRd")[c(3, 5, 7)], "black"),
    breaks = c(1, 10, 100, Inf),
    labels = unname(c("1 (Chitnis)", "10 (Chitnis)", "100 (Chitnis)", TeX("$\\infty$ (Ross-Macdonald)")))
  ) +
  scale_fill_manual(
    name = "Vertebrate host biting tolerance\n(bites per host per day)",
    values = c(met.brewer("VanGogh3", 3, direction = 1), "black"),
    # values = c(brewer.pal(9, "YlOrRd")[c(3, 5, 7)], "black"),
    breaks = c(1, 10, 100, Inf),
    labels = unname(c("1 (Chitnis)", "10 (Chitnis)", "100 (Chitnis)", TeX("$\\infty$ (Ross-Macdonald)")))
  ) +
  
  # faceting:
  facet_wrap(~system_ID,
             scales = "free",
             nrow = 2,
             labeller = labeller()
  ) +
  
  # theme options:
  theme_minimal_hgrid(11) +
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
ggsave("results/Topt_plots.png", Topt_plots,
       width = 16, height = 9)

###* Figure: Topt mean as function of KH and sigmaH ----

Topt_heat_df <- Topt_df %>%
  group_by(system_ID, Model, sigmaH, KH) %>%
  summarise(
    mean_val = mean(Topt),
    median_val = median(Topt),
    std_val = sd(Topt),
    # mode_val = mlv(norm_R0, method = 'mfv'),
    .groups = "keep",
    across()
  ) %>%
  arrange(system_ID, sigmaH, KH, mean_val, median_val)# , mode_val)


Toptalt_df <- Topt_heat_df %>%
  ungroup() %>%
  distinct() %>%
  # Arrange along the plotting variables
  arrange(system_ID, sigmaH, KH, mean_val) %>% 
  mutate(sigmaH_proxy = sigmaH)

# add many proxy sigmaH points at infinity to make the values easier to see there
## Add a "point at infinity" for sigmaH (NB: to be added to get-analysis-dfs.R)
# Value of sigmaH to replace 'Inf' with for plotting
infinite_skew <- 10^3 + 2500
# Used to create a gap dividing finite and infinite sigmaH
infinite_divide <- 10^mean(c(3, log10(infinite_skew)))

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
  arrange(sigmaH_proxy) %>%
  distinct()

Toptalt_df <- Toptalt_df %>% 
  select(-sample_num) %>% 
  distinct()


# Define sigmaH tick breaks and labels
sigmaH_breaks <- c(10^c(seq(-1, 3, by = 1)), infinite_divide, infinite_skew)
sigmaH_labels <- c(paste(c(10^c(seq(-1, 3, by = 1)))), TeX("$\\ldots$"), TeX("$\\infty$"))
Topt_breaks <- seq(21, 31, by = 1)

# Plot mean of Topt heat map
Toptalt_plots <- Toptalt_df %>%
  ## Set up plot ##
  # x = biting tolerance (with a point at infinity),
  # y = vertebrate host population density,
  # color = transmission thermal optimum
  ggplot(aes(x = KH, y = sigmaH_proxy, z = mean_val)) +
  
  # Topt for FINITE sigmaH:
  geom_contour_filled(
    data = filter(Toptalt_df, sigmaH_proxy < infinite_divide**0.98),
    breaks = Topt_breaks
  ) +
  
  # Topt for INFINITE sigmaH:
  geom_contour_filled(
    data = filter(Toptalt_df, sigmaH_proxy > infinite_divide**1.02),
    size = 1,
    breaks = Topt_breaks
  ) +
  
  # R0 = 1 contour:
  geom_contour(data = Toptalt_df,
               aes(z = R0opt, colour = "R0=1"), breaks = c(1), size = 1) +
  
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
    name = "Transmission thermal optimum (°C)",
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
    title.position = "top", title.hjust = 0.5,
    barwidth = unit(5, "cm"),
    show.limits = TRUE
  )) +
  
  # theme options:
  theme_minimal(11) +
  theme(
    legend.box = "vertical",
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
ggsave("results/Topt_mean.png", Toptalt_plots,
       width = 16, height = 9)


###* Figure: Topt variance as function of KH and sigmaH ----
test_df <- Toptalt_df %>% 
  select(system_ID, sigmaH_proxy, KH, std_val, mean_val) %>% 
  distinct()

# Plot mean of Topt heat map
Toptstdev_plots <- Toptalt_df %>% 
  select(system_ID, sigmaH_proxy, KH, std_val, mean_val) %>% 
  distinct() %>% 
  ## Set up plot ##
  # x = biting tolerance (with a point at infinity),
  # y = vertebrate host population density,
  # color = transmission thermal optimum
  ggplot(aes(x = KH, y = sigmaH_proxy, z = std_val)) +
  
  # Topt for FINITE sigmaH:
  geom_contour_filled(
    data = filter(test_df, sigmaH_proxy < infinite_divide**0.98)
  ) +
  
  # Topt for INFINITE sigmaH:
  geom_contour_filled(
    data = filter(test_df, sigmaH_proxy > infinite_divide**1.02),
    size = 1
  ) +
  
  # R0 = 1 contour:
  # geom_contour(data = Toptalt_df,
  #              aes(z = R0opt, colour = "R0=1"), breaks = c(1), size = 1) +
  
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
    name = "Standard deviation of\ntransmission thermal optimum (°C)",
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
    title.position = "top", title.hjust = 0.5,
    barwidth = unit(5, "cm"),
    show.limits = TRUE
  )) +
  
  # theme options:
  theme_minimal(11) +
  theme(
    legend.box = "vertical",
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
ggsave("results/Topt_stdev.png", Toptstdev_plots,
       width = 16, height = 9)


# 6) Functions to sample from thermal trait parameter distributions----
