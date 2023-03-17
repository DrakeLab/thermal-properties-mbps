data_in <- data.in

distinct_combos <- distinct(data_in, trait.name, system_ID)
sample_num <- sample(1:dim(distinct_combos)[1], 1)


system_sample <- distinct_combos$system_ID[sample_num]
trait_in <- distinct_combos$trait.name[sample_num]

mosquito_in <-  filter(data.in, system_ID == system_sample) %>% 
  dplyr::select(mosquito_species) %>% unique() %>% as.character()
pathogen_in <-  filter(data.in, system_ID == system_sample) %>% 
  dplyr::select(pathogen) %>% unique() %>% as.character()


# trait_in = "EFD"


# mosquito_in = "Culex quinquefasciatus"

# pathogen_in = "none"


n.chains <- 2 # 5
n.adapt <- 100 # 5000
n.samps <- 100 # 5000

data_in <- data.in

samples <- tibble(
  trait = as.character(),
  system_ID = as.character(),
  T0 = as.double(),
  Tm = as.double(),
  c = as.double(), # !!! note that we're using c as a generic parameter for Briere or Quadratic
  sample_num = as.double(),
  func = as.character()
)

distinct_combos <- distinct(data_in, trait.name, system_ID)


for (sample_num in 1:dim(distinct_combos)[1]) {
# for (sample_num in problem_set) {
  # if (sample_num %in% c(27,28, 29, 62, 63)) {next} # skip the problem samples
  print(paste0("System # ", sample_num, " --------------------------------------------------------------------------"))
  system_sample <- distinct_combos$system_ID[sample_num]
  trait_in <- distinct_combos$trait.name[sample_num]

  mosquito_in <-  filter(data.in, system_ID == system_sample) %>%
    dplyr::select(mosquito_species) %>% unique() %>% as.character()
  pathogen_in <-  filter(data.in, system_ID == system_sample) %>%
    dplyr::select(pathogen) %>% unique() %>% as.character()


  temp_sample <- thermtrait.prior.sample(data_in, trait_in, mosquito_in, pathogen_in,
                                     n.chains, n.adapt, n.samps,
                                     old_informative = FALSE)
  
  temp_sample <- temp_sample %>%
    mutate(trait = trait_in,
           system_ID = system_sample)
  samples <- rbind(samples, temp_sample)

}

# Save samples for now
write_csv(samples, "data/clean/temp_samples.csv")

# -------------------------------------------------------------------------


# samples <- read_csv("data/clean/temp_samples.csv")

library(reshape2)
library(cowplot)

# 2) Histograms of parameter distributions of thermal traits ----
plot_df <- samples %>%
  dplyr::select(-func) %>% 
  mutate(temp_diff = Tm - T0) %>% 
  mutate(logc = log(c)) %>% 
  melt(id = c("system_ID", "trait", "sample_num")) %>% 
  filter(system_ID %in% c("Aedes aegypti / DENV", "Aedes aegypti / none",
                          "Aedes aegypti / ZIKV", "Aedes aegypti / none",
                          "Aedes albopictus / DENV", "Aedes albopictus / none",
                          "Culex quinquefasciatus / WNV", "Culex quinquefasciatus / none",
                          "Anopheles spp. / Plasmodium spp.",
                          "Anopheles spp. / none"
  )) 

###* Figure: thermal trait parameter posterior distributions ----
parm_hists <- plot_df %>%
  filter(variable != "c") %>% 
  ggplot(aes(value, color = system_ID, fill = system_ID)) +
  geom_histogram(aes(), bins = 100) +
  # geom_density() +
  facet_grid(trait ~ variable, scales = "free") +
  theme_minimal_grid(12)

# distributions should be clumped (except for c)
# logc should be clumped
# temp_diff should be positive

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
Linear <- function(q, Tmax) {
  function(t) {
    pmax(q*(Tmax - t), 0, na.RM = FALSE)
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
    function_type == "Linear" ~ Linear(parms$c, parms$Tm)
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
TPC_df <- samples %>%
  # filter(sample_num %in% seq(1, thin_size)) %>%
  full_join(list(Temperature = Temps), by = character(), copy = TRUE) %>%
  # group_by(sample_num) %>%
  # try sapply or lapply
  # mutate(Trait_val = get.thermal.response(.,Temperature))Linear(parms$c, parms$T0)
  mutate(Trait_val = case_when(
    func == "Briere" ~ Briere(c, T0, Tm)(Temperature),
    func == "Quadratic" ~ Quadratic(c, T0, Tm)(Temperature),
    func == "Linear" ~ Linear(c, Tm)(Temperature)
  )) %>% 
  dplyr::select(-c("c", "T0", "Tm")) #%>% 
  # filter(system_ID %in% c("Aedes aegypti / DENV", "Aedes aegypti / none",
  #                         "Aedes aegypti / ZIKV", "Aedes aegypti / none",
  #                         "Aedes albopictus / DENV", "Aedes albopictus / none",
  #                         "Culex quinquefasciatus / WNV", "Culex quinquefasciatus / none",
  #                         "Anopheles spp. / Plasmodium spp.",
  #                         "Anopheles spp. / none"
  # ))
  # 

# get mean TPC from samples
meanTPC_df <- TPC_df %>% 
  group_by(system_ID, trait, Temperature) %>% 
  summarise(mean_val = mean(Trait_val), .groups = "keep")


# get edges of 89% HCI of samples
quantsTPC_df <- TPC_df %>% 
  group_by(system_ID, trait, Temperature) %>% 
  mutate(lowHCI_val = quantile(Trait_val, 0.055)) %>% 
  mutate(highHCI_val = quantile(Trait_val, 0.945)) %>% 
  dplyr::select(-c("sample_num", "Trait_val", "func"))


###* Figure: TPC curves with 89% high confidence intervals ---- 
TPC_plot <- TPC_df %>%
  group_by(sample_num) %>%
  arrange(Temperature)  %>%
  # group_by()
  ggplot() +
  geom_line(data = meanTPC_df,
            aes(x = Temperature, y = mean_val, color = system_ID)) +
  # 89% HCI of R0 TPC curves
  geom_ribbon(
    data = quantsTPC_df,
    aes(x = Temperature, ymin = lowHCI_val, ymax = highHCI_val, fill = system_ID),
    alpha = 0.1
  ) +
  facet_wrap(~ trait, scales = "free", ncol = 2) +
  theme_minimal_grid(12) #+
  # geom_point(test_df, mapping = aes(x = T, y = trait))
