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
n.samps <- 10 # 5000

test_df <- thermtrait.prior.sample(data_in, trait_in, mosquito_in, pathogen_in,
                                   n.chains = 2, n.adapt = 100, n.samps = 100,
                                   old_informative = FALSE)
samples <- tibble(
  trait = as.character(),
  system_ID = as.character(),
  T0 = as.double(),
  Tm = as.double(),
  c = as.double() # !!! note that we're using c as a generic parameter for Briere or Quadratic
)
distinct_combos <- distinct(data_in, trait.name, system_ID)
for (sample_num in 1:dim(distinct_combos)[1]) {
  
  system_sample <- distinct_combos$system_ID[sample_num]
  trait_in <- distinct_combos$trait.name[sample_num]
  
  mosquito_in <-  filter(data.in, system_ID == system_sample) %>% 
    dplyr::select(mosquito_species) %>% unique() %>% as.character()
  pathogen_in <-  filter(data.in, system_ID == system_sample) %>% 
    dplyr::select(pathogen) %>% unique() %>% as.character()
  
  
  temp_sample <- thermtrait.prior.sample(data_in, trait_in, mosquito_in, pathogen_in,
                                     n.chains = 2, n.adapt = 100, n.samps = 100,
                                     old_informative = FALSE)
  temp_sample <- temp_sample %>% 
    mutate(trait = trait_in,
           system_ID = system_sample)
  
  samples <- rbind(samples, temp_sample)
}












