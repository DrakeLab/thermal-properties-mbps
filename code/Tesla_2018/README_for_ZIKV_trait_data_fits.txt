# Metadata for Tesla et al. "Temperature drives Zika virus transmission: evidence from empirical and mathematical models"
# code files

# Contents:
# R scripts:
# 1. Informative_ZIKV_IndividualParameterFitCode.R
# 2. ZIKV_IndividualParameterFitCode.R
# 3. ZIKV_R0_plotting.R
# 4. mcmc_utils_all.R
# 5. temp_functions_all.R

# Rsave data files:
# 1. LifespanFits_2016-03-30-informative.Rsave
# 2. Informative_Aegypti_DENV_ParameterFits_2016-03-30.Rsave
# 3. zikv_trait_fits_informative.Rsave
# 4. zikv_trait_fits_uninformative.Rsave
# 5. previous_DENV_fits.Rsave
# 6. ZIKV_model_outputs-uninformative.Rsave

# CSV data files:
# 1. zikv_traits.csv
# 2. ZIKV_trait_data_fits.csv

# Data in CSV file 1, is used to fit trait thermal responses with either less-informative priors (R script 2) or data-informed priors (R script 1 with Rsave files 1-2). These fitting procedures draw on functions described in R scripts 4-5

# Fitted thermal response functions are plotted and incorporated into R0 in R script 3. The fitted R0 is then compared with previous R0 fits (Rsave file 5).

# All outputs are saved for future use (Rsave files 3, 4, 6 and CSV file 2).