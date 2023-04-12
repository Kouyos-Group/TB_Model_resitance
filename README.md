# TB__MODEL_RESISTANCE 

In this directory, there are all the scripts used during the internship. 

Warning: you have to change the pathways to load scripts/data. 

## Script 

The script folder contains the Model/deterministic folder, which is used to simulate the model. It also contains the Likelihood folder, which we use to determine some key parameters. 

### model/deterministic

This folder contains everything needed to simulate our model and answer our question "Is the TRUNCATE strategy a solution if we look at the population level?" It contains the file execution.R as well as the folder scenarios. 

**execution.R**

This file allows to simulate the model and to see the difference between two scenarios in the number of infected, the number resistant to the first line of treatment, the number resistant to the second line of treatment, the proportion of infected resistant to the first line of treatment and finally the proportion of infected resistant to the second line of treatment. 

The two scenarios tested are with the TRUNCATE strategy and with the standard strategy. 

The initial conditions and the parameters are fixed through the literature or through the likelihood method.  
In the first stage, we simulate from 1971 (introduction of the first line of treatment with rifampicin) to 2015 (introduction of the second line of treatment with bedaquine) with only two treatments, the first line and another treatment in case of failure. In a second phase, from 2015 to 2060, a second line of treatment is introduced in addition to the two previous treatments. 

#### scenarios

1. intro_rifampicin_44_year.R 
2. standard_strategy.R
3. TRUNCATE_strategy.R

#### data

fit.RData

### Parametrization

In this file, there is all the content used to estimate some parameters. The likelihood_gaussian.R script allows to launch the estimation. In the data file, there are the different data used to fit the model. These data are available on the WHO website. In the scenarios folder, there is the scenario_1971-2015.csv with the introduction of the first line of treatment (rifampicin) and scenario_2015-2021.R with the introduction of the second line of treatment (bedaqualine).  

likelihood_gaussian.R

##### data

1. MDR_RR_burden_estimates_2023-02-22.csdv
  
2. TB_burden_countries_2023-02-22.csv

##### scenarios

1. scenario_1971-2015.R

2. scenario_2015-2022.R
