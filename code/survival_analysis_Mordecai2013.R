library(tidyverse)

# Treat survival data separately.
# - Use data from the first day of the experiment, one day before reaching 0.01
#   threshold to three days after the threshold, giving 6 data points
# - Use the exponent(?) of the constant mortality function as an independent 
#   observation of mortality at each temperature
survival.Mordecai2013 <- read.csv("data/raw/Mordecai_2013/survival_data.csv", header = TRUE) %>%
  pivot_longer(cols = X5.C:X40.C, names_to = "T", values_to = "prop.alive") %>% 
  rename(day = D) %>% 
  mutate(T = readr::parse_number(T)) %>% 
  # replace NAs (beyond the first five days) with 0's
  mutate(prop.alive = ifelse(day > 5 & is.na(prop.alive), 0, prop.alive))

# !!! add in days prior to experiment, so that Gompertz fit is non-singular. This is not justified.
fake_data <- expand.grid(day = seq(-20,0), T = unique(survival.Mordecai2013$T), prop.alive = 1)

survival.Mordecai2013 <- rbind(survival.Mordecai2013, fake_data) %>% 
  arrange(day)

# Gompertz function
Gompertz <- function(x, k, lag){
  ymax <- 1
  y0 <- 0
  result <- y0 + (ymax -y0)*exp(-exp(k*(lag-x)/(ymax-y0) + 1) )
  # result <- pmin(1, result)
  return(result)
}

survival_fit.df <- tibble(T = numeric(), y0 = numeric(), ymax = numeric(), 
                          k = numeric(), lag = numeric())
# had to remove extreme high and low temperatures because fit is singular 
for (TT in unique(survival.Mordecai2013$T)) {
  df <- filter(survival.Mordecai2013, T == TT)
  # Fit Gompertz curves to the data
  Gomp_fit <- nls(1-prop.alive ~ Gompertz(day, k, lag),
                  data = df,
                  start = list(k = 0.1, lag = 0))
  Gomp_coefs <- coef(Gomp_fit)
  
  survival_fit.df <- add_row(survival_fit.df, T = TT, k = Gomp_coefs["k"], 
                             lag = Gomp_coefs["lag"])
}

Gomp_df <- tibble(day = unique(survival.Mordecai2013$day)) %>% 
  full_join(survival_fit.df, by = character()) %>% 
  mutate(pred.p = Gompertz(x = day, k = k, lag = lag))

survival_plot <- survival.Mordecai2013 %>% 
  arrange(day) %>% 
  # filter(prop.alive != 0) %>%
  ggplot(aes(x = day, y = 1-prop.alive, color = as.factor(T), group = T)) +
  geom_point() +
  geom_path(data = Gomp_df, aes(x = day, y = pred.p))

survival_test <- survival.Mordecai2013 %>% 
  # This is roughly linear
  mutate(k = (1-log(-log(1-prop.alive)))/(day - 8)) %>% 
  filter(is.finite(k)) %>% 
  # mutate(c = (-(log((-log(1-prop.alive) )/ day)/day))) %>% 
  # filter(is.finite(c)) %>% 
  ggplot(aes(x = day, y = k, group = T, color = as.factor(T))) +
  geom_point() +
  ylim(-1, 1)
  