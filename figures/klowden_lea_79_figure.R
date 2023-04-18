library(tidyverse)
library(cowplot)

defense_df <- tibble(
  Species = c("Aedes triseriatus", "Aedes aegypti", "Aedes sollicitans",
              "Aedes taeniorhynchus", "Anopheles crucians","Aedes triseriatus", "Aedes aegypti", "Aedes sollicitans",
              "Aedes taeniorhynchus", "Anopheles crucians"),
  restrained = c(1,1,1,1,1,0,0,0,0,0),
  fed_percent = c(.86,1,.67,1,NA,.1,.01,.02,.04,.06),
  n = c(83,117,164,20,0,240,122,1252,796,32)
) %>%
  mutate(estimate_fed = floor(n*fed_percent)) %>% 
  group_by(Species, restrained) %>% 
  rowwise() %>% 
  mutate(low_CI = ifelse(!is.na(fed_percent), 
                         binom.test(estimate_fed, n, p = 0.5, conf.level = 0.95)$conf.int[1],
                         NA)) %>% 
  mutate(high_CI = ifelse(!is.na(fed_percent), 
                         binom.test(estimate_fed, n, p = 0.5, conf.level = 0.95)$conf.int[2],
                         NA)) %>% 
  mutate(restrained = ifelse(restrained, "Restrained", "Unrestrained"))


defense_plot <- defense_df %>% 
  ggplot(aes(x = restrained, y = fed_percent, fill = Species)) +
  geom_boxplot(aes(ymin = low_CI, lower = low_CI, middle = fed_percent, upper = high_CI, ymax = high_CI),
               stat = "identity",
               outlier.shape = NA) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Percent of mosquitoes fed") +
  ggtitle("Across species, mosquito feeding success is significantly higher when hosts are restrained (Klowden & Lea 1979)") +
  theme_minimal_hgrid(16) 