# Forecasting of cohort fertility under a hierarchical Bayesian approach
R code to generate cohort fertility forecasts using the hierarchical Bayesian model (Model hB) proposed in: "Forecasting of cohort fertility under a hierarchical Bayesian approach" by Joanne Ellison, Erengul Dodd and Jonathan J. Forster.

## Set-up
To run this code you will need R version 3.4.0 or above and be able to install packages, in particular rstan, which is required to implement the Hamiltonian Monte Carlo methodology to fit the model. This will probably require you to also install the package Rtools, explained in the rstan installation instructions that can be found at https://github.com/stan-dev/rstan/wiki.

## Data
This code downloads data from the Human Fertility Database (HFD) directly. You will need to have set up an HFD account (go to https://www.humanfertility.org/Home/Index), the username and password for which will be requested in the code.

## Generating forecasts
Follow the instructions in the source file (source_file.r) to complete the set-up, generate and process the Model hB forecasts.
