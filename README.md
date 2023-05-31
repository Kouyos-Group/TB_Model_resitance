# TB__MODEL_RESISTANCE 

In this directory, there are all the scripts used during the internship. 

Warning: you have to change the pathways to load scripts/data. 

## Likelihood

likelihood_gaussian_Georgia.R

likelihood_gaussian_Georgia_no_normalize.R

likelihood_gaussian_Georgia_sd_scale.R

see_result.R

### scenarios

1. scenario_2000-2015.R 
2. scenario_2015-2021.R

### data

1. MDR_RR_burden_estimates_2023-02-22.csdv
  
2. TB_burden_countries_2023-02-22.csv

## Result/Georgia/rho_low

The results for the different likelihood method techniques are shown here.

fit.RData

### sd_scale

fit.RData

### no_normalize

fit.RData

## Model

This folder contains everything needed to simulate our model and answer our question "Is the TRUNCATE strategy a solution if we look at the population level?" It contains the file execution.R as well as the folder scenarios. 

**execution.R**

This file allows to simulate the model and to see the difference between two scenarios in the number of infected, the number resistant to the first line of treatment, the number resistant to the second line of treatment, the proportion of infected resistant to the first line of treatment and finally the proportion of infected resistant to the second line of treatment. 

The two scenarios tested are with the TRUNCATE strategy and with the standard strategy. 

### scenarios

1. standard_strategy.R

2. Truncate_strategy.R
