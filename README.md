# `maldetectr`: an R package that implements site occupancy models for estimating prevelance of malaria parasites

## Installation

Install `devtools`:

    install.packages( "devtools" )

Install package:

    devtools::install_github( "andersgs/maldetectr" )

## Current version

Implements a zero-inflated binomial model.

## Usage

Maximum likelihood estimate using a zero-inflated binomial model.

    library( maldetectr )
    #simulate some data
    n_sites <- 30 # number of individuals tested
    n_samples <- 5 # number of times each indivdual was tested
    p_occupancy <- 0.2 # prevalence of malaria in population
    p_detection <- 0.8 # probability of test detecting malaria given that individual
                       # has malaria
    data <- maldetectr::sim_zib_data( n_sites = n_sites, 
                                      n_samples = n_tests, 
                                      p_oc = p_occupancy, 
                                      p_dt = p_detection, 
                                      seed = 42 )
    ml_est <- maldetectr::get_zib_ml( data = data, 
                                              ci = 0.89 )

## Issues

Please use the GitHub issues page
