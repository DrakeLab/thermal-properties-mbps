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
data.in.thermchars <- read_rds("results/AllThermChar_thin.rds") #%>% 
  # Focus on two systems for now
  # filter(system_ID %in% c("Aedes albopictus / DENV", "Anopheles spp. / Plasmodium"))
data.in.outputs <- read_rds("data/clean/AllOutputs_thin.rds") #%>% 
  # Focus on two systems for now
  # filter(system_ID %in% c("Aedes albopictus / DENV", "Anopheles spp. / Plasmodium"))

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


R0_df <- read_rds("results/R0_TPC_data.rds") %>% 
  filter(KH %in% 10^seq(-2,5)) %>% 
  arrange(system_ID, KH, sigmaH, Temperature)

R0_plot <- ggplot(mapping = aes(x = Temperature, group = KH)) +
  # path of mean, normalized R0 as a function of temperature (finite sigmaH):
  geom_path(
    data = filter(R0_df, sigmaH == 100),
    aes(y = mean, colour = KH, linetype = "finite")
  ) +
  # path of mean, normalized R0 as a function of temperature (infinite sigmaH):
  geom_path(
    data = filter(R0_df, sigmaH == Inf),
    aes(y = mean, linetype = "infinite", colour = "black"),
    colour = "black",
    lwd = 1.5
  ) +
  # 89% HCI of R0 TPC curves
  geom_ribbon(
    data = filter(R0_df, sigmaH == 100),
    aes(ymin = lowHCI, ymax = highHCI, fill = KH),
    alpha = 0.1
  ) +
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

# Plot
R0_plot <- shift_legend(R0_plot)

ggsave("figures/R0_TPCs.svg", R0_plot,
       width = 16, height = 9)

# Figure 5: Topt heatmap --------------------------------------------------
# Mean value of Topt as a function of KH and sigmaH
Topt_heat_df <- read_rds("results/Topt_heat_data.rds")

Toptalt_df <- Topt_heat_df %>%
  ungroup() %>%
  distinct() %>%
  # Arrange along the plotting variables
  rename(mean_Topt = mean_val) %>% 
  arrange(system_ID, sigmaH, KH, mean_Topt) %>% 
  mutate(sigmaH_proxy = sigmaH)

# add many proxy sigmaH points at infinity to make the values easier to see there
## Add a "point at infinity" for sigmaH (NB: to be added to get-analysis-dfs.R)
# Value of sigmaH to replace 'Inf' with for plotting
infinite_skew <- 10^2 + 80
# Used to create a gap dividing finite and infinite sigmaH
infinite_divide <- 10^mean(c(2, log10(infinite_skew)))

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


# Define sigmaH tick breaks and labels
sigmaH_breaks <- c(10^c(seq(-1, 2, by = 1)), infinite_divide, infinite_skew)
sigmaH_labels <- c(paste(c(10^c(seq(-1, 2, by = 1)))), TeX("$\\ldots$"), TeX("$\\infty$"))
Topt_breaks <- seq(floor(min(Toptalt_df$mean_Topt)), 
                   ceiling(max(Toptalt_df$mean_Topt)), by = 1)

# Plot mean of Topt heat map
Toptalt_plots <- Toptalt_df %>%
  ## Set up plot ##
  # x = biting tolerance (with a point at infinity),
  # y = vertebrate host population density,
  # color = transmission thermal optimum
  ggplot(aes(x = KH, y = sigmaH_proxy, z = mean_Topt)) +
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
  # # R0 = 1 contour: # !!! figure this out
  # geom_contour(data = Toptalt_df,
  #              aes(z = mean_R0opt, colour = "R0=1"), breaks = c(1), size = 1) +
  
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
  theme_minimal(11) +
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
ggsave("figures/Topt_mean.svg", Toptalt_plots,
       width = 16, height = 9)


# Figure 6: CTwidth heatmap -----------------------------------------------
# Shows how the parasite thermal niche (temperatures between CT_min and CT_max) 
# depends on vertebrate host availability (density and biting tolerance) in each
# of the focal systems.

#   x    = biting tolerance (with a point at infinity)
#   y    = vertebrate host population density
# colour = width of parasite thermal niche (CTmax - CTmin)
# Facets = systems

CTwidth_df <- data.in.thermchars %>%
  # remove irrelevant columns
  select(-c(threshold_bool, CHmin, CHmax, R0opt, Topt)) %>%
  ungroup() %>%
  # Calculate thermal niche width
  mutate(CTwidth = CTmax - CTmin, .keep = "unused") %>% 
  # Replace NAs with zeroes for CTwidth 
  # NAs correspond to settings where R0 is always less than one, so CTmin/max are undefined
  mutate(CTwidth = replace(CTwidth, is.na(CTwidth), 0)) %>%
  group_by(system_ID, Model, sigmaH, KH) %>%
  mutate(
    mean_CTwidth = mean(CTwidth),
    median_CTwidth = median(CTwidth),
    std_CTwidth = sd(CTwidth),
    # mode_val = mlv(norm_R0, method = 'mfv'),
    # .groups = "keep"
  ) %>% 
  # collapse groups
  select(-c(sample_num, CTwidth)) %>% 
  unique() %>% 
  mutate(sigmaH_proxy = sigmaH) %>% 
  # Arrange along the plotting variables
  arrange(system_ID, KH, sigmaH_proxy, mean_CTwidth) %>% 
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
  ggplot(aes(x = KH, y = sigmaH_proxy, z = mean_CTwidth)) +
  
  # CTwidth for FINITE sigmaH:
  geom_contour_filled(
    data = filter(CTwidth_df, sigmaH_proxy < infinite_divide**0.98)
  ) +
  
  # CTwidth for INFINITE sigmaH:
  geom_contour_filled(
    data = filter(CTwidth_df, sigmaH_proxy > infinite_divide**1.02)
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
  theme_minimal(11) +
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
ggsave("figures/CTwidth_mean.svg", CTWidth_heat_plots,
       width = 16, height = 9)

# *Supplementary Figures* ----


# Figure S1: mean R0 as function of temperature and biting tolerance --------


# Figure S2: Topt as a function of KH -------------------------------------
Topt_df <- read_rds("results/Topt_KH.rds") %>% 
  select(-c(variable)) %>%
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
  ggplot(aes(
    x = KH
  )) +
  # Topt curves:
  geom_path(aes(y = mean, colour = as.factor(sigmaH)), lwd = 1) +
  # # 89% HCI of R0 TPC curves
  # geom_ribbon(
  #   
  #   aes(ymin = lowHCI, ymax = highHCI,
  #       fill = as.factor(sigmaH)),
  #   alpha = 0.05
  # ) +
  # Add dotted lines showing limits of ribbons
  geom_path(
    aes(y = lowHCI,
        color = as.factor(sigmaH)),
    linetype = "dashed"
  ) +
  geom_path(
    
    aes(y = highHCI,
        color = as.factor(sigmaH)),
    linetype = "dashed"
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
    expand = c(0.05, 0.05),
    limits = c(11,39)
  ) +
  # # color:
  # scale_colour_manual(
  #   name = "Vertebrate host biting tolerance\n(bites per host per day)",
  #   values = c(met.brewer("VanGogh3", 4, direction = 1), "black"),
  #   # values = c(brewer.pal(9, "YlOrRd")[c(3, 5, 7)], "black"),
  #   breaks = c(1, 10, 100, Inf),
  #   labels = unname(c("1 (Chitnis)", "10 (Chitnis)", "100 (Chitnis)", TeX("$\\infty$ (Ross-Macdonald)")))
  # ) +
  # scale_fill_manual(
  #   name = "Vertebrate host biting tolerance\n(bites per host per day)",
  #   values = c(met.brewer("VanGogh3", 4, direction = 1), "black"),
  #   # values = c(brewer.pal(9, "YlOrRd")[c(3, 5, 7)], "black"),
  #   breaks = c(1, 10, 100, Inf),
  #   labels = unname(c("1 (Chitnis)", "10 (Chitnis)", "100 (Chitnis)", TeX("$\\infty$ (Ross-Macdonald)")))
  # ) +
  # faceting:
  facet_wrap(~system_ID,
             scales = "free",
             nrow = 2,
             labeller = labeller()) +
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
Topt_plot <- shift_legend(Topt_plot)

ggsave("figures/ToptKH_plot.svg", Topt_plot,
       width = 16, height = 9)

# Figure S3: Topt as a function of sigmaH ---------------------------------
# Mean Topt as a function of sigmaH, across discrete values of KH

Topt_plot <- read_rds("results/Topt_KH.rds") %>% 
  # Choose appropriate pop. dens. values and group
  # filter(KH %in% 10^seq(0, 4)) %>%
  group_by(KH) %>%
  # Arrange along the plotting variables
  arrange(system_ID, sigmaH, mean_Topt) %>%
  arrange(KH) %>%
  ungroup()%>%
  ## Set up plot ##
  # color = sigmaH
  ggplot(aes(
    x = sigmaH
  )) +
  # Topt curves:
  geom_path(aes(y = mean, colour = as.factor(KH)), lwd = 1) +
  # 89% HCI of R0 TPC curves
  geom_ribbon(
    
    aes(ymin = lowHCI, ymax = highHCI, 
        fill = as.factor(KH)),
    alpha = 0.05
  ) +
  # Add dotted lines showing limits of ribbons
  geom_path(
    
    aes(y = lowHCI,
        color = as.factor(KH)),
    linetype = "dashed"
  ) +
  geom_path(
    
    aes(y = highHCI,
        color = as.factor(KH)),
    linetype = "dashed"
  ) +
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
    limits = c(5, 3500),
    breaks = sigmaH_breaks,
    labels = sigmaH_labels
  ) +
  # y-axis
  scale_y_continuous(
    name = expression("Thermal optimum "(degree * C)),
    expand = c(0.05, 0.05)
  ) +
  
  # color:
  scale_colour_manual(
    name = "Vertebrate host population\ndensity (ind/ha)",
    values = met.brewer("Hokusai3", 11,
                        # option = "plasma"
    )) +
  scale_fill_manual(
    name = "Vertebrate host population\ndensity (ind/ha)",
    values = met.brewer("Hokusai3", 11,
                        # option = "plasma"
    )) +
  # faceting:
  facet_wrap(~system_ID,
             scales = "free",
             nrow = 2,
             labeller = labeller()  ) +
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
Topt_plot <- shift_legend(Topt_plot)

ggsave("figures/ToptsigmaH_plot.svg", Topt_plot,
       width = 16, height = 9)

# Figure S4: CTmin heatmap ------------------------------------------------


# Figure S5: CTmax heatmap ------------------------------------------------


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
  theme_minimal(11) +
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
ggsave("figures/Topt_stdev.svg", Toptstdev_plots,
       width = 16, height = 9)

