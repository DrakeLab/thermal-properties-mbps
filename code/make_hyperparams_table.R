load("code/Mordecai2017/aedes_prior_gamma_fits.Rsave")

library(reshape2)
library(tidyverse)

trait_names <- c('a', 'b', 'c', 'e2a', 'EFD', 'lf', 'MDR', 'PDR', 'TFD')

hyperparams_df <- data.frame(trait = "a", melt(gamma.fits.a), multiplier = 0.1) %>% 
  rbind(data.frame(trait = "b", melt(gamma.fits.b), multiplier = 0.5)) %>%  
  rbind(data.frame(trait = "c", melt(gamma.fits.c), multiplier = 0.5)) %>%  
  rbind(data.frame(trait = "e2a", melt(gamma.fits.e2a), multiplier = 0.1)) %>%  
  rbind(data.frame(trait = "EFD", melt(gamma.fits.EFD), multiplier = 0.1)) %>%  
  rbind(data.frame(trait = "lf", melt(gamma.fits.lf), multiplier = 0.01) )%>%  
  rbind(data.frame(trait = "MDR", melt(gamma.fits.MDR), multiplier = 0.1)) %>%  
  rbind(data.frame(trait = "PDR", melt(gamma.fits.PDR), multiplier = 0.5) )%>% 
  rbind(data.frame(trait = "TFD", melt(gamma.fits.TFD), multiplier = 0.1))

write_csv(hyperparams_df, "data/clean/gamma_fits.csv")
