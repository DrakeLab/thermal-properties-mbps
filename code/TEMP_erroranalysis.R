library(tidyverse)
library(reshape2)

temp_vec <- seq(10,40)

briere_func <- function(c,T0,Tm) {
  function(t) {pmax(c*t*(t-T0)*(Tm-t)^(1/2),0,na.rm = TRUE)}
}

mean_curve <- briere_func(2.02E-4, 13.35, 40.08)(temp_vec)
lower_curve <- briere_func(1.20E-4, 8.27, 40.00)(temp_vec)
upper_curve <- briere_func(2.80E-4, 17.41, 40.28)(temp_vec)

plot_df <- data.frame(temp = temp_vec, mean = mean_curve, lower = lower_curve, upper = upper_curve)

plot_df %>% 
  melt(id = "temp") %>% 
  ggplot(aes(x = temp, y = value, color = variable)) +
  geom_line()