## For questions, contact Kyle Dahlin, kydahlin@gmail.com
## Initialized Oct 2020
##
## Title: Functions for getting model outputs like R0 #############
##
## Project: Global zoonoses - spillover of mosquito-borne pathogens
##
## Purpose: Provide functions for evaluating several model properties and equilibrium values
##
## Contents: A set of functions listed below
##
##          compute.RH: Evaluates the host type reproduction number
##
##          compute.RV: Evaluates the mosquito type reproduction number
##
##          compute.R0: Evaluates the basic reproduction number (actually just sqrt(RH*RV))
##
##          compute.CHmin: Evaluates the critical minimum host density for pathogen persistence
##
##          compute.CHmax: Evaluates the critical maximum host density for pathogen persistence


# Helper function: Reduces resolution of a vector by sub-sampling to length new_length
res_reduce <- function(df, new_length) {
  old_length <- length(unique(df))
  if (old_length < new_length) {
    warning("new_length is larger than old length. Vector will be unchanged")
    new_length <- old_length
  }
  ret <- unique(df)[seq(1, length(unique(df)), length.out = new_length)]
}

## 1) Functions for computing outputs ##########################################

### Equilibrium total vector population-----------------------------------------
compute.V0 <- function(input) {
  with(input, {
    eps <- .Machine$double.eps
    if(is.null(lf)) {lf = 1 / (muV + eps)}
    temp <- sigmaV * fecundity * deltaL
    # check if net reproduction rate exceeds mortality rate
    temp_bool <- temp > (1 / lf)
    ifelse(temp_bool,
           V0 <- KL * rhoL * lf * (1 - 1 / (lf * temp + eps)) ,
           V0 <- 0 # if mortality exceeds reproduction, set to zero
    )
  })
}

### Host-to-mosquito reproduction number----------------------------------------
compute.RH <- function(input) {
  with(input, {
    eps <- .Machine$double.eps
    # Derive intermediate quantities
    bV <- ifelse(is.infinite(sigmaH),
                 sigmaV, # Ross-Macdonald model
                 sigmaV * sigmaH * KH / (sigmaH * KH + sigmaV * V0 + eps)
    )
    RH <- ifelse(V0 == 0, 
                 0, 
                 bV * betaV * V0 * exp(-1 / (lf * etaV)) / (KH * (gammaH + muH) + eps))
  })
}

### Mosquito-to-host reproduction number----------------------------------------
compute.RV <- function(input) {
  with(input, {
    eps <- .Machine$double.eps
    RV <- ifelse(is.infinite(sigmaH),
                 sigmaV * betaH / (1 / (lf + eps)), # Ross-Macdonald
                 sigmaH * sigmaV * betaH * KH / ((1 / (lf + eps)) * (sigmaH * KH + sigmaV * V0))
    )
  })
}

### Basic reproduction number---------------------------------------------------
# Computed via the type reproduction numbers (see above)
compute.R0 <- function(input) {
  RH <- compute.RH(input)
  RV <- compute.RV(input)
  # Return R0
  return(sqrt(RH * RV))
}

### Host density limits---------------------------------------------------------
# Critical minimum host density for pathogen persistence (R0>1)
# Only applicable for the Chitnis dynamic model
compute.CHmin <- function(input) {
  with(input, {
    eps <- .Machine$double.eps
    # intermediate quantities
    chi <- 0.5 * sigmaV * betaV * exp(-muV / etaV) * sigmaH * betaH *
      (1 / (gammaH + muH)) * (1 / (muV + eps))
    prefactor <- sigmaV * V0 / (sigmaH + eps)
    Hmin <- ifelse(is.infinite(sigmaH), 
                   0,
                   prefactor * (chi - 1 - sqrt(chi^2 - 2 * chi))
    )
  })
}

# Critical maximum host density for pathogen persistence (R0>1)
compute.CHmax <- function(input) {
  with(input, {
    eps <- .Machine$double.eps
    # intermediate quantities
    chi <- 0.5 * sigmaV * betaV * exp(-muV / etaV) * sigmaH * betaH *
      (1 / (gammaH + muH)) * (1 / (muV + eps))
    prefactor <- sigmaV * V0 / (sigmaH + eps)
    
    Hmax <- ifelse(is.infinite(sigmaH),
                   (sigmaV^2) * V0 * betaV * exp(-muV / etaV) * betaH *
                     (1 / (gammaH + muH)) * (1 / muV),
                   prefactor * (chi - 1 + sqrt(chi^2 - 2 * chi)))
  })
}