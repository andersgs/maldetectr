# `maldetectr`: an R package that implements site occupancy models for estimating prevelance of malaria parasites

## Installation

Install `devtools`:

    install.packages( "devtools" )

Install package:

    devtools::install_git( "andersgs/maldetectr" )

## Current version

Implements a zero-inflated binomial model.

## Usage

Maximum likelihood estimate using a zero-inflated binomial model.

    library( maldetectr )
    #simulate some data
    n_sites <- 30
    n_tests <- 5
    p_occupancy <- 0.2
    p_detection <- 0.8
    data <- maldetectr::sim_zib_data( n_sites = n_sites, 
                                      n_tests = n_tests, 
                                      p_oc = p_occupancy, 
                                      p_dt = p_detection, 
                                      seed = 42 )
    ml_est <- maldetectr::get_zib_likelihood( data = data, 
                                              ci = 0.89 )

## Issues

Please use the GitHub issues page