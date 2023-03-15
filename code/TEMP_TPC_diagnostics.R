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

data_in <- data.in

distinct_combos <- distinct(data_in, trait.name, system_ID)

for (sample_num in 69:dim(distinct_combos)[1]) {
  if (sample_num %in% c(27,28, 29, 62, 63)) {next} # skip the problem samples
  # if (sample_num %in% c(27, 28, 29, 62, 63)) {next} # skip the problem samples
  system_sample <- distinct_combos$system_ID[sample_num]
  trait_in <- distinct_combos$trait.name[sample_num]
  
  mosquito_in <-  filter(data.in, system_ID == system_sample) %>% 
    dplyr::select(mosquito_species) %>% unique() %>% as.character()
  pathogen_in <-  filter(data.in, system_ID == system_sample) %>% 
    dplyr::select(pathogen) %>% unique() %>% as.character()
  
  
  temp_sample <- thermtrait.prior.sample(data_in, trait_in, mosquito_in, pathogen_in,
                                     n.chains = 2, n.adapt = 100, n.samps = 10,
                                     old_informative = FALSE)
  temp_sample <- temp_sample %>% 
    mutate(trait = trait_in,
           system_ID = system_sample)
  print(sample_num)
  samples <- rbind(samples, temp_sample)
  
}

###* trouble systems: 
###* 27 [Error in stats::optim(x = c(238.546983753754, 390.491816190212, 208.251158143909,  : non-finite finite-difference value [1]]
###* 28, [Error in stats::optim(x = c(451.773190350689, 425.454518866436, 363.892114765623,  : non-finite finite-difference value [1]]
###* 29, [Error in stats::optim(x = c(862.729906398741, 808.199692220595, 632.023043163558,  : non-finite finite-difference value [1]]
###* 62, [Error in stats::optim(x = c(1502.33981220254, 1555.74156690337, 3083.40032794262,  : non-finite finite-difference value [1]]
###* 63, [Error in stats::optim(x = c(3275.04483451616, 6311.67010262804, 6066.91078603536,  : non-finite finite-difference value [1]]











