# Remove current working environment
rm(list = ls())
# Necessary libraries
if(!requireNamespace("pacman", quietly = T))
  install.packages("pacman")
pacman::p_load("deSolve", "chron", "tidyverse", "deSolve", "tidyverse", "stringr", "dplyr", "ggplot2", "ggpubr")

##### 2000-2022 ####
fit<-readRDS("~/Bureau/stage/script/T_T/result/Georgia/rho_low/fit.RData")


fixed <- c( 
  #birth
  pi = 0.0008991667, 
  
  #rate for infection (alpha)
  infection = exp(fit$par["infection"]), 
  resistance_cost= exp(fit$par["resistance_cost"]),
  
  #progression to active disease          
  delta1 =0.03333 , delta2 = 0.00017, 
  
  #progression rate from early latent to late latent
  epsilon = 0.304 , 
  
  ####treated  rate
  beta= 0.62,
  
  ####cured with different treatment (t = treatment)
  cured_t1 = 0.86,
  cured_t2 = 0.59,
  cured_t3 = 0.59,
  
  #time treamtent
  time_1=1/6, 
  time_2=1/6, 
  time_3=1/6,
  
  #self cured
  omicron=0.01667, 
  
  #acquisition resistance
  rho= exp(fit$par["rho"]), 
  
  #death due to tuberculosis
  theta= exp(fit$par["theta"]),  
  
  #lambda
  l = exp(fit$par["l"]),
  #l=exp(fit$par["l"]),
  
  #natural death
  mu= 0.0008966667, 
  
  #relapse
  phi = 0.000, 
  
  #number in differents compartments
  e = exp(fit$par["e"]), 
  late = exp(fit$par["late"]),
  i = exp(fit$par["i"]), 
  t = exp(fit$par["t"]), 
  r = exp(fit$par["r"]), 
  
  #part of resistance in the pop
  resist = exp(fit$par["resist"]), 
  treat = exp(fit$par["treat"])
)

names(fixed) <- c("pi", "infection", "resistance_cost", "delta1", "delta2", "epsilon", "beta", "cured_t1", "cured_t2", "cured_t3", "time_1", "time_2", 
                  "time_3", "omicron", "rho", "theta", "l", "mu","phi","e", "late", "i", "t", "r", "resist", "treat")

init <- c(S = 70000, 
          
          Es00 = 30000 * fixed["e"] * (1 - fixed["resist"]) / (fixed["e"] + fixed["late"] + fixed["i"] + fixed["t"] + fixed["r"])  , 
          Es10 = 30000 * fixed["e"] * fixed["resist"] / (fixed["e"] + fixed["late"] + fixed["i"] + fixed["t"] + fixed["r"]), 
          Es01 = 0, 
          Es11 = 0,
          
          
          Ls00 = 30000 * fixed["late"] * (1 - fixed["resist"]) / (fixed["e"] + fixed["late"] + fixed["i"] + fixed["t"] + fixed["r"]), 
          Ls10 = 30000 * fixed["late"] * fixed["resist"] / (fixed["e"] + fixed["late"] + fixed["i"] + fixed["t"] + fixed["r"]),
          Ls01 = 0, 
          Ls11 = 0, 
          
          
          Is00 = 30000 * fixed["i"] * (1 - fixed["resist"]) / (fixed["e"] + fixed["late"] + fixed["i"] + fixed["t"] + fixed["r"]), 
          Is10 = 30000 * fixed["i"] * fixed["resist"] / (fixed["e"] + fixed["late"] + fixed["i"] + fixed["t"] + fixed["r"]), 
          Is01 = 0, 
          Is11 = 0, 
          
          Ts00t1 = 30000 * fixed["t"] * (1 - fixed["resist"]) * (1 - fixed["treat"]) / (fixed["e"] + fixed["late"] + fixed["i"] + fixed["t"] + fixed["r"]), 
          Ts10t1 = 30000 * fixed["t"] * fixed["resist"] * (1 - fixed["treat"]) / (fixed["e"] + fixed["late"] + fixed["i"] + fixed["t"] + fixed["r"]) , 
          Ts01t1 = 0, 
          Ts11t1 = 0, 
          Ts00t2 = 0, 
          Ts10t2 = 0, 
          Ts01t2 = 0, 
          Ts11t2 = 0, 
          Ts00t3 = 30000 * fixed["t"] * (1 - fixed["resist"]) * fixed["treat"] / (fixed["e"] + fixed["late"] + fixed["i"] + fixed["t"] + fixed["r"]), 
          Ts10t3 = 30000 * fixed["t"] * fixed["resist"] * fixed["treat"] / (fixed["e"] + fixed["late"] + fixed["i"] + fixed["t"] + fixed["r"]), 
          Ts01t3 = 0, 
          Ts11t3 = 0, 
          
          Rs00 = 30000 * fixed["r"] * (1 - fixed["resist"]) / (fixed["e"] + fixed["late"] + fixed["i"] + fixed["t"] + fixed["r"]), 
          Rs10 = 30000 * fixed["r"] * fixed["resist"] / (fixed["e"] + fixed["late"] + fixed["i"] + fixed["t"] + fixed["r"]), 
          Rs01 = 0, 
          Rs11 = 0, 
          
          incidence = 0 , 
          M=0, 
          rR=0
)

names(init)<-c("S", 
               
               "Es00", 
               "Es10", 
               "Es01", 
               "Es11",
               
               
               "Ls00", 
               "Ls10",
               "Ls01", 
               "Ls11", 
               
               
               "Is00", 
               "Is10", 
               "Is01", 
               "Is11", 
               
               "Ts00t1", 
               "Ts10t1", 
               "Ts01t1", 
               "Ts11t1", 
               "Ts00t2", 
               "Ts10t2", 
               "Ts01t2", 
               "Ts11t2", 
               "Ts00t3", 
               "Ts10t3", 
               "Ts01t3", 
               "Ts11t3", 
               
               "Rs00", 
               "Rs10", 
               "Rs01", 
               "Rs11", 
               
               "incidence", 
               "M", 
               "rR")


source("/home/louis/Bureau/stage/script/T_T/likelihood/scenarios/scenario_2000-2015.R")
times <- seq(0,15 * 12, by=1)
simulation <- as.data.frame(ode(init,times,sir_equations,fixed))
simulation <- simulation[-1,]

#second period (introduction bedaqualine)
init<-tail(simulation, n=1)[2:33]
init2<-unlist(init)
names(init2)<-names(init)
init<-init2
times <- seq(15*12, 22*12 , by=1)
source("/home/louis/Bureau/stage/script/T_T/likelihood/scenarios/scenario_2015-2021.R")
simulation2 <-  as.data.frame(ode(init,times,sir_equations,fixed))
simulation2 <- simulation2[-1,]

#keep only the data that correspond to the data of WHO
simulation<- rbind(simulation, simulation2)


sir_values_0<-data.frame(simulation)

data0<-data.frame("time" = sir_values_0[,1], 
                   "S" = sir_values_0[,2], 
                   "E" = rowSums(sir_values_0[,3:6]), 
                   "L"= rowSums(sir_values_0[,7:10]),
                   "I"=rowSums(sir_values_0[,11:14]), 
                   "T"= rowSums(sir_values_0[,15:26]), 
                   "R" = rowSums(sir_values_0[,27:30]),
                   "Infected" = rowSums(sir_values_0[,11:14]),
                   "number_sensible" = rowSums(sir_values_0[,str_sub(names(sir_values_0), 1 , 1) == "I" & 
                                                              (str_sub(names(sir_values_0), 3 , 3) == "0" | str_sub(names(sir_values_0), 4 , 4) == "0")]),
                   "number_1st_line_resistant" = rowSums(sir_values_0[,str_sub(names(sir_values_0), 1 , 1) == "I" & 
                                                                        (str_sub(names(sir_values_0), 3 , 3) == "1")]),
                   "number_2nd_line_resistant" = rowSums(sir_values_0[,str_sub(names(sir_values_0), 1 , 1) == "I" & str_sub(names(sir_values_0), 4 , 4) == "1"]),
                  "freq_1st_line_resistant" = sir_values_0[,12]/rowSums(sir_values_0[,11:14]),
                  "freq_2nd_line_resistant" = sir_values_0[,13]/rowSums(sir_values_0[,11:14]), 
                  "freq_both_line_resistant" = sir_values_0[,14]/rowSums(sir_values_0[,11:14]),
                  "freq_resistant" = (sir_values_0[,13] + sir_values_0[,12] + sir_values_0[,14])/rowSums(sir_values_0[,11:14]),
                   "Incidence" = sir_values_0[,31],
                   "scenario" = rep("2000 - 2015", each = length(sir_values_0[,1])))

init<-tail(simulation, n=1)[2:33]
init2<-unlist(init)
names(init2)<-names(init)
init<-c(init2, nT = 0)
times <- seq(22 * 12, 40 *12 , by=1)

##### Keep standard scenario ####
source("/home/louis/Bureau/stage/script/T_T/model/scenarios/standard_strategy.R")
simulation <-  as.data.frame(ode(init,times,sir_equations,fixed))

sir_values_1<-data.frame(simulation)

data1<-data.frame("time" = sir_values_1[,1], 
                  "S" = sir_values_1[,2], 
                  "E" = rowSums(sir_values_1[,3:6]), 
                  "L"= rowSums(sir_values_1[,7:10]),
                  "I"=rowSums(sir_values_1[,11:14]), 
                  "T"= rowSums(sir_values_1[,15:26]), 
                  "R" = rowSums(sir_values_1[,27:30]),
                  "Infected" = rowSums(sir_values_1[,11:14]),
                  "Prevalence" = rowSums(sir_values_1[,11:14])/rowSums(sir_values_1[,2:30]),
                  "number_sensible" = rowSums(sir_values_1[,str_sub(names(sir_values_1), 1 , 1) == "I" & 
                                                             (str_sub(names(sir_values_1), 3 , 3) == "0" | str_sub(names(sir_values_1), 4 , 4) == "0")]),
                  "number_1st_line_resistant" = rowSums(sir_values_1[,str_sub(names(sir_values_1), 1 , 1) == "I" & 
                                                                       (str_sub(names(sir_values_1), 3 , 3) == "1")]),
                  "number_2nd_line_resistant" = rowSums(sir_values_1[,str_sub(names(sir_values_1), 1 , 1) == "I" & str_sub(names(sir_values_1), 4 , 4) == "1"]),
                  "freq_1st_line_resistant" = (sir_values_1[,12] + sir_values_1[,14])/rowSums(sir_values_1[,11:14]),
                  "freq_2nd_line_resistant" = (sir_values_1[,13] + sir_values_1[,14])/rowSums(sir_values_1[,11:14]), 
                  "freq_resistant" = (sir_values_1[,13] + sir_values_1[,12] + sir_values_1[,14])/rowSums(sir_values_1[,11:14]),
                  "freq_both_line_resistant" = sir_values_1[,14]/rowSums(sir_values_1[,11:14]),
                  "Incidence" = sir_values_1[,31], 
                  "scenario" = rep("Standard Strategy", each = length(sir_values_1[,1])))

##### TRUNCATE strategy ####
init_T <- c(init, 
            Ts00t4 = 0 , Ts10t4 = 0, Ts01t4 = 0, Ts11t4 = 0,
            Ts00t5 = 0 , Ts10t5 = 0, Ts01t5 = 0, Ts11t5 = 0)

fixed_T2 <- c(fixed, cured_t4 = 0.575, cured_t5 = 0.575, time_4 = 1/2, time_5 = 1, rho_T = 0.005) 
source("~/Bureau/stage/script/T_T/model/scenarios/TRUNCATE_strategy.R")
simulation <-  as.data.frame(ode(init_T,times,sir_equations,fixed_T2))


sir_values_3<-data.frame(simulation)

data3<-data.frame("time" = sir_values_3[,1], 
                  "S" = sir_values_3[,2], 
                  "E" = rowSums(sir_values_3[,3:6]), 
                  "L"= rowSums(sir_values_3[,7:10]),
                  "I"=rowSums(sir_values_3[,11:14]), 
                  "T"= rowSums(sir_values_3[,15:26]) + rowSums(sir_values_3[,34:41]), 
                  "R" = rowSums(sir_values_3[,27:30]),
                  "Infected" = rowSums(sir_values_3[,11:14]),
                  "Prevalence" = rowSums(sir_values_3[,11:14])/(rowSums(sir_values_3[,2:30])+rowSums(sir_values_3[35:42])),
                  "number_sensible" = rowSums(sir_values_3[,str_sub(names(sir_values_3), 1 , 1) == "I" & 
                                                             (str_sub(names(sir_values_3), 3 , 3) == "0" | str_sub(names(sir_values_3), 4 , 4) == "0")]),
                  "number_1st_line_resistant" = rowSums(sir_values_3[,str_sub(names(sir_values_3), 1 , 1) == "I" & 
                                                                       (str_sub(names(sir_values_3), 3 , 3) == "1")]),
                  "number_2nd_line_resistant" = rowSums(sir_values_3[,str_sub(names(sir_values_3), 1 , 1) == "I" & str_sub(names(sir_values_3), 4 , 4) == "1"]),
                  "freq_1st_line_resistant" = (sir_values_3[,12] + sir_values_3[,14])/rowSums(sir_values_3[,11:14]),
                  "freq_2nd_line_resistant" = (sir_values_3[,13] + sir_values_3[,14])/rowSums(sir_values_3[,11:14]), 
                  "freq_resistant" = (sir_values_3[,13] + sir_values_3[,12] + sir_values_3[,14])/rowSums(sir_values_3[,11:14]),
                  "freq_both_line_resistant" = sir_values_3[,14]/rowSums(sir_values_3[,11:14]),
                  "Incidence" = sir_values_3[,31], 
                  "scenario" = rep("TS - treatment success 0.575", each = length(sir_values_3[,1])))


##### TRUNCATE strategy 3####
fixed_T3 <- c(fixed, cured_t4 = 0.647, cured_t5 = 0.647, time_4 = 1/2, time_5 = 1, rho_T = 0.005) 
source("~/Bureau/stage/script/T_T/model/scenarios/TRUNCATE_strategy.R")
simulation <-  as.data.frame(ode(init_T,times,sir_equations,fixed_T3))


sir_values_4<-data.frame(simulation)

data4<-data.frame("time" = sir_values_4[,1], 
                  "S" = sir_values_4[,2], 
                  "E" = rowSums(sir_values_4[,3:6]), 
                  "L"= rowSums(sir_values_4[,7:10]),
                  "I"=rowSums(sir_values_4[,11:14]), 
                  "T"= rowSums(sir_values_4[,15:26]) + rowSums(sir_values_4[,34:41]), 
                  "R" = rowSums(sir_values_4[,27:30]),
                  "Infected" = rowSums(sir_values_4[,11:14]),
                  "Prevalence" = rowSums(sir_values_4[,11:14])/(rowSums(sir_values_4[,2:30])+rowSums(sir_values_4[35:42])),
                  "number_sensible" = rowSums(sir_values_4[,str_sub(names(sir_values_4), 1 , 1) == "I" & 
                                                             (str_sub(names(sir_values_4), 3 , 3) == "0" | str_sub(names(sir_values_4), 4 , 4) == "0")]),
                  "number_1st_line_resistant" = rowSums(sir_values_4[,str_sub(names(sir_values_4), 1 , 1) == "I" & 
                                                                       (str_sub(names(sir_values_4), 3 , 3) == "1")]),
                  "number_2nd_line_resistant" = rowSums(sir_values_4[,str_sub(names(sir_values_4), 1 , 1) == "I" & str_sub(names(sir_values_4), 4 , 4) == "1"]),
                  "freq_1st_line_resistant" = (sir_values_4[,12] + sir_values_4[,14])/rowSums(sir_values_4[,11:14]),
                  "freq_2nd_line_resistant" = (sir_values_4[,13] + sir_values_4[,14])/rowSums(sir_values_4[,11:14]), 
                  "freq_resistant" = (sir_values_4[,13] + sir_values_4[,12] + sir_values_4[,14])/rowSums(sir_values_4[,11:14]),
                  "freq_both_line_resistant" = sir_values_4[,14]/rowSums(sir_values_4[,11:14]),
                  "Incidence" = sir_values_4[,31],
                  "scenario" = rep("TS - treatment success 0.647", each = length(sir_values_4[,1])))


##### TRUNCATE strategy 3####
fixed_T5 <- c(fixed, cured_t4 = 0.862, cured_t5 = 0.862, time_4 = 1/2, time_5 = 1, rho_T = 0.005) 
source("~/Bureau/stage/script/T_T/model/scenarios/TRUNCATE_strategy.R")
simulation <-  as.data.frame(ode(init_T,times,sir_equations,fixed_T5))


sir_values_5<-data.frame(simulation)

data5<-data.frame("time" = sir_values_5[,1], 
                  "S" = sir_values_5[,2], 
                  "E" = rowSums(sir_values_5[,3:6]), 
                  "L"= rowSums(sir_values_5[,7:10]),
                  "I"=rowSums(sir_values_5[,11:14]), 
                  "T"= rowSums(sir_values_5[,15:26]) + rowSums(sir_values_5[,34:41]), 
                  "R" = rowSums(sir_values_5[,27:30]),
                  "Infected" = rowSums(sir_values_5[,11:14]),
                  "Prevalence" = rowSums(sir_values_5[,11:14])/(rowSums(sir_values_5[,2:30])+rowSums(sir_values_5[35:42])),
                  "number_sensible" = rowSums(sir_values_5[,str_sub(names(sir_values_5), 1 , 1) == "I" & 
                                                             (str_sub(names(sir_values_5), 3 , 3) == "0" | str_sub(names(sir_values_5), 4 , 4) == "0")]),
                  "number_1st_line_resistant" = rowSums(sir_values_5[,str_sub(names(sir_values_5), 1 , 1) == "I" & 
                                                                       (str_sub(names(sir_values_5), 3 , 3) == "1")]),
                  "number_2nd_line_resistant" = rowSums(sir_values_5[,str_sub(names(sir_values_5), 1 , 1) == "I" & str_sub(names(sir_values_5), 4 , 4) == "1"]),
                  "freq_1st_line_resistant" = (sir_values_5[,12] + sir_values_5[,14])/rowSums(sir_values_5[,11:14]),
                  "freq_2nd_line_resistant" = (sir_values_5[,13] + sir_values_5[,14])/rowSums(sir_values_5[,11:14]), 
                  "freq_resistant" = (sir_values_5[,13] + sir_values_5[,12] + sir_values_5[,14])/rowSums(sir_values_5[,11:14]),
                  "freq_both_line_resistant" = sir_values_5[,14]/rowSums(sir_values_5[,11:14]),
                  "Incidence" = sir_values_5[,31],
                  "scenario" = rep("TS - treatment success 0.862", each = length(sir_values_5[,1])))


##### TRUNCATE strategy 7####
fixed_T7 <- c(fixed, cured_t4 = 0.575, cured_t5 = 0.575, time_4 = 1/2, time_5 = 1, rho_T = 0.01) 
source("~/Bureau/stage/script/T_T/model/scenarios/TRUNCATE_strategy.R")
simulation <-  as.data.frame(ode(init_T,times,sir_equations,fixed_T7))


sir_values_7<-data.frame(simulation)

data7<-data.frame("time" = sir_values_7[,1], 
                  "S" = sir_values_7[,2], 
                  "E" = rowSums(sir_values_7[,3:6]), 
                  "L"= rowSums(sir_values_7[,7:10]),
                  "I"=rowSums(sir_values_7[,11:14]), 
                  "T"= rowSums(sir_values_7[,15:26]) + rowSums(sir_values_7[,34:41]), 
                  "R" = rowSums(sir_values_7[,27:30]),
                  "Infected" = rowSums(sir_values_7[,11:14]),
                  "Prevalence" = rowSums(sir_values_7[,11:14])/(rowSums(sir_values_7[,2:30])+rowSums(sir_values_7[35:42])),
                  "number_sensible" = rowSums(sir_values_7[,str_sub(names(sir_values_7), 1 , 1) == "I" & 
                                                             (str_sub(names(sir_values_7), 3 , 3) == "0" | str_sub(names(sir_values_7), 4 , 4) == "0")]),
                  "number_1st_line_resistant" = rowSums(sir_values_7[,str_sub(names(sir_values_7), 1 , 1) == "I" & 
                                                                       (str_sub(names(sir_values_7), 3 , 3) == "1")]),
                  "number_2nd_line_resistant" = rowSums(sir_values_7[,str_sub(names(sir_values_7), 1 , 1) == "I" & str_sub(names(sir_values_7), 4 , 4) == "1"]),
                  "freq_1st_line_resistant" = (sir_values_7[,12] + sir_values_7[,14])/rowSums(sir_values_7[,11:14]),
                  "freq_2nd_line_resistant" = (sir_values_7[,13] + sir_values_7[,14])/rowSums(sir_values_7[,11:14]), 
                  "freq_resistant" = (sir_values_7[,13] + sir_values_7[,12] + sir_values_7[,14])/rowSums(sir_values_7[,11:14]),
                  "freq_both_line_resistant" = sir_values_7[,14]/rowSums(sir_values_7[,11:14]),
                  "Incidence" = sir_values_7[,31], 
                  "scenario" = rep("TS - treatment success 0.575", each = length(sir_values_7[,1])))


##### TRUNCATE strategy 8####
fixed_T8 <- c(fixed, cured_t4 = 0.647 , cured_t5 = 0.647, time_4 = 1/2, time_5 = 1, rho_T = 0.01) 
source("~/Bureau/stage/script/T_T/model/scenarios/TRUNCATE_strategy.R")
simulation <-  as.data.frame(ode(init_T,times,sir_equations,fixed_T8))


sir_values_8<-data.frame(simulation)

data8<-data.frame("time" = sir_values_8[,1], 
                  "S" = sir_values_8[,2], 
                  "E" = rowSums(sir_values_8[,3:6]), 
                  "L"= rowSums(sir_values_8[,7:10]),
                  "I"=rowSums(sir_values_8[,11:14]), 
                  "T"= rowSums(sir_values_8[,15:26]) + rowSums(sir_values_8[,34:41]), 
                  "R" = rowSums(sir_values_8[,27:30]),
                  "Infected" = rowSums(sir_values_8[,11:14]),
                  "Prevalence" = rowSums(sir_values_8[,11:14])/(rowSums(sir_values_8[,2:30])+rowSums(sir_values_8[35:42])),
                  "number_sensible" = rowSums(sir_values_8[,str_sub(names(sir_values_8), 1 , 1) == "I" & 
                                                             (str_sub(names(sir_values_8), 3 , 3) == "0" | str_sub(names(sir_values_8), 4 , 4) == "0")]),
                  "number_1st_line_resistant" = rowSums(sir_values_8[,str_sub(names(sir_values_8), 1 , 1) == "I" & 
                                                                       (str_sub(names(sir_values_8), 3 , 3) == "1")]),
                  "number_2nd_line_resistant" = rowSums(sir_values_8[,str_sub(names(sir_values_8), 1 , 1) == "I" & str_sub(names(sir_values_8), 4 , 4) == "1"]),
                  "freq_1st_line_resistant" = (sir_values_8[,12] + sir_values_8[,14])/rowSums(sir_values_8[,11:14]),
                  "freq_2nd_line_resistant" = (sir_values_8[,13] + sir_values_8[,14])/rowSums(sir_values_8[,11:14]), 
                  "freq_resistant" = (sir_values_8[,13] + sir_values_8[,12] + sir_values_8[,14])/rowSums(sir_values_8[,11:14]),
                  "freq_both_line_resistant" = sir_values_8[,14]/rowSums(sir_values_8[,11:14]),
                  "Incidence" = sir_values_8[,31],
                  "scenario" = rep("TS - treatment success 0.647", each = length(sir_values_8[,1])))


##### TRUNCATE strategy 9####
fixed_T9 <- c(fixed, cured_t4 = 0.862 , cured_t5 = 0.862, time_4 = 1/2, time_5 = 1, rho_T = 0.01) 
source("~/Bureau/stage/script/T_T/model/scenarios/TRUNCATE_strategy.R")
simulation <-  as.data.frame(ode(init_T,times,sir_equations,fixed_T9))


sir_values_9<-data.frame(simulation)

data9<-data.frame("time" = sir_values_9[,1], 
                  "S" = sir_values_9[,2], 
                  "E" = rowSums(sir_values_9[,3:6]), 
                  "L"= rowSums(sir_values_9[,7:10]),
                  "I"=rowSums(sir_values_9[,11:14]), 
                  "T"= rowSums(sir_values_9[,15:26]) + rowSums(sir_values_9[,34:41]), 
                  "R" = rowSums(sir_values_9[,27:30]),
                  "Infected" = rowSums(sir_values_9[,11:14]),
                  "Prevalence" = rowSums(sir_values_9[,11:14])/(rowSums(sir_values_9[,2:30])+rowSums(sir_values_9[35:42])),
                  "number_sensible" = rowSums(sir_values_9[,str_sub(names(sir_values_9), 1 , 1) == "I" & 
                                                             (str_sub(names(sir_values_9), 3 , 3) == "0" | str_sub(names(sir_values_9), 4 , 4) == "0")]),
                  "number_1st_line_resistant" = rowSums(sir_values_9[,str_sub(names(sir_values_9), 1 , 1) == "I" & 
                                                                       (str_sub(names(sir_values_9), 3 , 3) == "1")]),
                  "number_2nd_line_resistant" = rowSums(sir_values_9[,str_sub(names(sir_values_9), 1 , 1) == "I" & str_sub(names(sir_values_9), 4 , 4) == "1"]),
                  "freq_1st_line_resistant" = (sir_values_9[,12] + sir_values_9[,14])/rowSums(sir_values_9[,11:14]),
                  "freq_2nd_line_resistant" = (sir_values_9[,13] + sir_values_9[,14])/rowSums(sir_values_9[,11:14]), 
                  "freq_resistant" = (sir_values_9[,13] + sir_values_9[,12] + sir_values_9[,14])/rowSums(sir_values_9[,11:14]),
                  "freq_both_line_resistant" = sir_values_9[,14]/rowSums(sir_values_9[,11:14]),
                  "Incidence" = sir_values_9[,31] ,
                  "scenario" = rep("TS - treatment success 0.862", each = length(sir_values_9[,1])))

##### TRUNCATE strategy 7####
fixed_T11 <- c(fixed, cured_t4 = 0.575, cured_t5 = 0.575, time_4 = 1/2, time_5 = 1, rho_T = 0.015) 
source("~/Bureau/stage/script/T_T/model/scenarios/TRUNCATE_strategy.R")
simulation <-  as.data.frame(ode(init_T,times,sir_equations,fixed_T11))


sir_values_11<-data.frame(simulation)

data11<-data.frame("time" = sir_values_11[,1], 
                  "S" = sir_values_11[,2], 
                  "E" = rowSums(sir_values_11[,3:6]), 
                  "L"= rowSums(sir_values_11[,7:10]),
                  "I"=rowSums(sir_values_11[,11:14]), 
                  "T"= rowSums(sir_values_11[,15:26]) + rowSums(sir_values_11[,34:41]), 
                  "R" = rowSums(sir_values_11[,27:30]),
                  "Infected" = rowSums(sir_values_11[,11:14]),
                  "Prevalence" = rowSums(sir_values_11[,11:14])/(rowSums(sir_values_11[,2:30])+rowSums(sir_values_11[35:42])),
                  "number_sensible" = rowSums(sir_values_11[,str_sub(names(sir_values_11), 1 , 1) == "I" & 
                                                             (str_sub(names(sir_values_11), 3 , 3) == "0" | str_sub(names(sir_values_11), 4 , 4) == "0")]),
                  "number_1st_line_resistant" = rowSums(sir_values_11[,str_sub(names(sir_values_11), 1 , 1) == "I" & 
                                                                       (str_sub(names(sir_values_11), 3 , 3) == "1")]),
                  "number_2nd_line_resistant" = rowSums(sir_values_11[,str_sub(names(sir_values_11), 1 , 1) == "I" & str_sub(names(sir_values_11), 4 , 4) == "1"]),
                  "freq_1st_line_resistant" = (sir_values_11[,12] + sir_values_11[,14])/rowSums(sir_values_11[,11:14]),
                  "freq_2nd_line_resistant" = (sir_values_11[,13] + sir_values_11[,14])/rowSums(sir_values_11[,11:14]), 
                  "freq_resistant" = (sir_values_11[,13] + sir_values_11[,12] + sir_values_11[,14])/rowSums(sir_values_11[,11:14]),
                  "freq_both_line_resistant" = sir_values_11[,14]/rowSums(sir_values_11[,11:14]),
                  "Incidence" = sir_values_11[,31], 
                  "scenario" = rep("TS - treatment success 0.575", each = length(sir_values_11[,1])))


##### TRUNCATE strategy 8####
fixed_T12 <- c(fixed, cured_t4 = 0.647 , cured_t5 = 0.647, time_4 = 1/2, time_5 = 1, rho_T = 0.015) 
source("~/Bureau/stage/script/T_T/model/scenarios/TRUNCATE_strategy.R")
simulation <-  as.data.frame(ode(init_T,times,sir_equations,fixed_T12))


sir_values_12<-data.frame(simulation)
data12<-data.frame("time" = sir_values_12[,1], 
                  "S" = sir_values_12[,2], 
                  "E" = rowSums(sir_values_12[,3:6]), 
                  "L"= rowSums(sir_values_12[,7:10]),
                  "I"=rowSums(sir_values_12[,11:14]), 
                  "T"= rowSums(sir_values_12[,15:26]) + rowSums(sir_values_12[,34:41]), 
                  "R" = rowSums(sir_values_12[,27:30]),
                  "Infected" = rowSums(sir_values_12[,11:14]),
                  "Prevalence" = rowSums(sir_values_12[,11:14])/(rowSums(sir_values_12[,2:30])+rowSums(sir_values_12[35:42])),
                  "number_sensible" = rowSums(sir_values_12[,str_sub(names(sir_values_12), 1 , 1) == "I" & 
                                                             (str_sub(names(sir_values_12), 3 , 3) == "0" | str_sub(names(sir_values_12), 4 , 4) == "0")]),
                  "number_1st_line_resistant" = rowSums(sir_values_12[,str_sub(names(sir_values_12), 1 , 1) == "I" & 
                                                                       (str_sub(names(sir_values_12), 3 , 3) == "1")]),
                  "number_2nd_line_resistant" = rowSums(sir_values_12[,str_sub(names(sir_values_12), 1 , 1) == "I" & str_sub(names(sir_values_12), 4 , 4) == "1"]),
                  "freq_1st_line_resistant" = (sir_values_12[,12] + sir_values_12[,14])/rowSums(sir_values_12[,11:14]),
                  "freq_2nd_line_resistant" = (sir_values_12[,13] + sir_values_12[,14])/rowSums(sir_values_12[,11:14]), 
                  "freq_resistant" = (sir_values_12[,13] + sir_values_12[,12] + sir_values_12[,14])/rowSums(sir_values_12[,11:14]),
                  "freq_both_line_resistant" = sir_values_12[,14]/rowSums(sir_values_12[,11:14]),
                  "Incidence" = sir_values_12[,31],
                  "scenario" = rep("TS - treatment success 0.647", each = length(sir_values_12[,1])))


##### TRUNCATE strategy 9####
fixed_T13 <- c(fixed, cured_t4 = 0.862, cured_t5 = 0.862, time_4 = 1/2, time_5 = 1, rho_T = 0.015) 
source("~/Bureau/stage/script/T_T/model/scenarios/TRUNCATE_strategy.R")
simulation <-  as.data.frame(ode(init_T,times,sir_equations,fixed_T13))


sir_values_13<-data.frame(simulation)

data13<-data.frame("time" = sir_values_13[,1], 
                  "S" = sir_values_13[,2], 
                  "E" = rowSums(sir_values_13[,3:6]), 
                  "L"= rowSums(sir_values_13[,7:10]),
                  "I"=rowSums(sir_values_13[,11:14]), 
                  "T"= rowSums(sir_values_13[,15:26]) + rowSums(sir_values_13[,34:41]), 
                  "R" = rowSums(sir_values_13[,27:30]),
                  "Infected" = rowSums(sir_values_13[,11:14]),
                  "Prevalence" = rowSums(sir_values_13[,11:14])/(rowSums(sir_values_13[,2:30])+rowSums(sir_values_13[35:42])),
                  "number_sensible" = rowSums(sir_values_13[,str_sub(names(sir_values_13), 1 , 1) == "I" & 
                                                             (str_sub(names(sir_values_13), 3 , 3) == "0" | str_sub(names(sir_values_13), 4 , 4) == "0")]),
                  "number_1st_line_resistant" = rowSums(sir_values_13[,str_sub(names(sir_values_13), 1 , 1) == "I" & 
                                                                       (str_sub(names(sir_values_13), 3 , 3) == "1")]),
                  "number_2nd_line_resistant" = rowSums(sir_values_13[,str_sub(names(sir_values_13), 1 , 1) == "I" & str_sub(names(sir_values_13), 4 , 4) == "1"]),
                  "freq_1st_line_resistant" = (sir_values_13[,12] + sir_values_13[,14])/rowSums(sir_values_13[,11:14]),
                  "freq_2nd_line_resistant" = (sir_values_13[,13] + sir_values_13[,14])/rowSums(sir_values_13[,11:14]), 
                  "freq_both_line_resistant" = sir_values_13[,14]/rowSums(sir_values_13[,11:14]),
                  "freq_resistant" = (sir_values_13[,13] + sir_values_13[,12] + sir_values_13[,14])/rowSums(sir_values_13[,11:14]),
                  "Incidence" = sir_values_13[,31] ,
                  "scenario" = rep("TS - treatment success 0.862", each = length(sir_values_13[,1])))



##########################################################
dataA<-rbind(data1,data3, data4, data5)
dataA$time <- 2000 + (dataA$time/12) 
dataA$time <- as.factor(dataA$time)

dataB <- rbind(data1, data7, data8, data9)
dataB$time <- 2000 + (dataB$time/12) 
dataB$time <- as.factor(dataB$time)

dataC <- rbind(data1, data11, data12, data13)
dataC$time <- 2000 + (dataC$time/12) 
dataC$time <- as.factor(dataC$time)

a <- ggplot(dataA, aes(x=time, y= Prevalence*100000, group=scenario), size = 0.8) +
  geom_line(aes(color=scenario),size = 1) + 
  theme_linedraw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) + 
  scale_x_discrete(breaks = seq(2020, 2040, by = 5)) +
  ylim(23,26) +
  labs(color = "Strategy") +
  theme(axis.title= element_text(size=15),
        legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=15), #change legend text font sizeaxis.title= element_text(size=14),
        axis.text.x = element_text(size=15, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15))


b <- ggplot(dataB, aes(x=time, y= Prevalence*100000, group=scenario), size = 0.8) +
  geom_line(aes(color=scenario),size = 1) + 
  theme_linedraw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) + 
  scale_x_discrete(breaks = seq(2020, 2040, by = 5)) +
  ylim(23,26) +
  theme(axis.title= element_text(size=15),
        legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=15), #change legend text font sizeaxis.title= element_text(size=14),
        axis.text.x = element_text(size=15, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15))


c <- ggplot(dataC, aes(x=time, y= Prevalence*100000, group=scenario), size = 0.8) +
  geom_line(aes(color=scenario),size = 1) + 
  theme_linedraw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) + 
  scale_x_discrete(breaks = seq(2020, 2040, by = 5 )) +
  ylim(23,26) +
  theme(axis.title= element_text(size=15),
        legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=15), #change legend text font sizeaxis.title= element_text(size=14),
        axis.text.x = element_text(size=15, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15))
  
  
annotate_figure(ggarrange(a, b, c, 
                          labels = c("A - TS mutation rate 0.005", "B - TS mutation rate 0.01", " C - TS mutation rate 0.015"),
                          font.label = list(size = 14, face = "bold"),
                          common.legend = T,
                          vjust=0,
                          hjust=-0.15,
                          legend = "right",
                          ncol = 3, nrow = 1),
                bottom = text_grob("Year",size = 20, hjust = 3),
                top = text_grob("", size = 20),
                left = text_grob("Prevalence per 100,000 individuals", rot=90, size = 20),
  
)

g <- ggplot(dataA, aes(x=time, y= freq_1st_line_resistant * Prevalence*100000, group=scenario)) +
  geom_line(aes(color=scenario),size = 1) + 
  theme_linedraw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) + 
  scale_x_discrete(breaks = seq(2020, 2040, by = 5)) +
 # ylim(min(c(dataA$freq_1st_line_resistant,dataB$freq_1st_line_resistant, dataC$freq_1st_line_resistant)),
  #     max(c(dataA$freq_1st_line_resistant,dataB$freq_1st_line_resistant, dataC$freq_1st_line_resistant))) +
  theme(axis.title= element_text(size=15),
        legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=15), #change legend text font sizeaxis.title= element_text(size=14),
        axis.text.x = element_text(size=15, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15))

h <- ggplot(dataB, aes(x=time, y= freq_1st_line_resistant * Prevalence*100000 , group=scenario)) +
  geom_line(aes(color=scenario),size = 1) + 
  theme_linedraw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) + 
  #ylim(min(c(dataA$freq_1st_line_resistant,dataB$freq_1st_line_resistant, dataC$freq_1st_line_resistant)),
   #    max(c(dataA$freq_1st_line_resistant,dataB$freq_1st_line_resistant, dataC$freq_1st_line_resistant))) +
  scale_x_discrete(breaks = seq(2020, 2040, by = 5)) +
  theme(axis.title= element_text(size=15),
        legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=15), #change legend text font sizeaxis.title= element_text(size=14),
        axis.text.x = element_text(size=15, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15))

i <- ggplot(dataC, aes(x=time, y= freq_1st_line_resistant * Prevalence*100000, group=scenario)) +
  geom_line(aes(color=scenario),size = 1) + 
  theme_linedraw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) + 
  scale_x_discrete(breaks = seq(2020, 2040, by = 5)) +
  #ylim(min(c(dataA$freq_1st_line_resistant,dataB$freq_1st_line_resistant, dataC$freq_1st_line_resistant)),
   #    max(c(dataA$freq_1st_line_resistant,dataB$freq_1st_line_resistant, dataC$freq_1st_line_resistant))) + 
  theme(axis.title= element_text(size=15),
        legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=15), #change legend text font sizeaxis.title= element_text(size=14),
        axis.text.x = element_text(size=15, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15))

annotate_figure(ggarrange(g, h, i, 
                          labels = c("acquisition resistance 0.005", "acquisition resistance 0.01", "acquisition resistance 0.015"),
                          font.label = list(size = 14, face = "bold"),
                          common.legend = T,
                          vjust=0,
                          hjust=-0.1,
                          legend = "right",
                          ncol = 3, nrow = 1),
                bottom = text_grob("Year",size = 20, hjust = 2),
                top = text_grob("", size = 20),
                left = text_grob("Prevalence Rifampicin resistant for 100,000 individuals", rot=90, size = 15),
                
)

j <- ggplot(dataA, aes(x=time, y= freq_2nd_line_resistant  * Prevalence*100000, group=scenario)) +
  geom_line(aes(color=scenario),size = 1) + 
  theme_linedraw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) + 
  scale_x_discrete(breaks = seq(2020, 2040, by = 5)) +
  ylim(min(c(dataA$freq_2nd_line_resistant * dataA$Prevalence *100000,
             dataB$freq_2nd_line_resistant * dataB$Prevalence *100000, 
             dataC$freq_2nd_line_resistant * dataC$Prevalence *100000)),
       max(c(dataA$freq_2nd_line_resistant * dataA$Prevalence *100000,
             dataB$freq_2nd_line_resistant * dataB$Prevalence *100000, 
             dataC$freq_2nd_line_resistant * dataC$Prevalence *100000))) + 
  theme(axis.title= element_text(size=15),
        legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=15), #change legend text font sizeaxis.title= element_text(size=14),
        axis.text.x = element_text(size=15, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15))

k <- ggplot(dataB, aes(x=time, y= freq_2nd_line_resistant  * Prevalence*100000, group=scenario)) +
  geom_line(aes(color=scenario),size = 1) + 
  theme_linedraw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) + 
  ylim(min(c(dataA$freq_2nd_line_resistant * dataA$Prevalence *100000,
             dataB$freq_2nd_line_resistant * dataB$Prevalence *100000, 
             dataC$freq_2nd_line_resistant * dataC$Prevalence *100000)),
       max(c(dataA$freq_2nd_line_resistant * dataA$Prevalence *100000,
             dataB$freq_2nd_line_resistant * dataB$Prevalence *100000, 
             dataC$freq_2nd_line_resistant * dataC$Prevalence *100000))) + 
  scale_x_discrete(breaks = seq(2020, 2040, by = 5)) +
  theme(axis.title= element_text(size=15),
        legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=15), #change legend text font sizeaxis.title= element_text(size=14),
        axis.text.x = element_text(size=15, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15))

l <- ggplot(dataC, aes(x=time, y= freq_2nd_line_resistant  * Prevalence*100000, group=scenario)) +
  geom_line(aes(color=scenario),size = 1) + 
  theme_linedraw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) + 
  scale_x_discrete(breaks = seq(2020, 2040, by = 5)) +
  ylim(min(c(dataA$freq_2nd_line_resistant * dataA$Prevalence *100000,
             dataB$freq_2nd_line_resistant * dataB$Prevalence *100000, 
             dataC$freq_2nd_line_resistant * dataC$Prevalence *100000)),
       max(c(dataA$freq_2nd_line_resistant * dataA$Prevalence *100000,
             dataB$freq_2nd_line_resistant * dataB$Prevalence *100000, 
             dataC$freq_2nd_line_resistant * dataC$Prevalence *100000))) + 
  theme(axis.title= element_text(size=15),
        legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=15), #change legend text font sizeaxis.title= element_text(size=14),
        axis.text.x = element_text(size=15, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15))

annotate_figure(ggarrange(j,k,l, 
                          labels = c("acquisition resistance 0.005", "acquisition resistance 0.01", "acquisition resistance 0.015"),
                          font.label = list(size = 14, face = "bold"),
                          common.legend = T,
                          vjust=0,
                          hjust=-0.1,
                          legend = "right",
                          ncol = 3, nrow = 1),
                bottom = text_grob("Year",size = 20, hjust = 3),
                top = text_grob("", size = 20),
                left = text_grob("Prevalence of Bedaqualine resistant per 100,000", rot=90, size = 15),
                
)

m <- ggplot(dataA, aes(x=time, y= freq_both_line_resistant  * Prevalence * 100000, group=scenario)) +
  geom_line(aes(color=scenario),size = 1) + 
  theme_linedraw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) + 
  scale_x_discrete(breaks = seq(2020, 2040, by = 5)) +
  ylim(min(c(dataA$freq_both_line_resistant * dataA$Prevalence *100000,
             dataB$freq_both_line_resistant * dataB$Prevalence *100000, 
             dataC$freq_both_line_resistant * dataC$Prevalence *100000)),
       max(c(dataA$freq_both_line_resistant * dataA$Prevalence *100000,
             dataB$freq_both_line_resistant * dataB$Prevalence *100000, 
             dataC$freq_both_line_resistant * dataC$Prevalence *100000))) + 
  theme(axis.title= element_text(size=15),
        legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=15), #change legend text font sizeaxis.title= element_text(size=14),
        axis.text.x = element_text(size=15, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15))

n <- ggplot(dataB, aes(x=time, y= freq_both_line_resistant  * Prevalence *100000, group=scenario)) +
  geom_line(aes(color=scenario),size = 1) + 
  theme_linedraw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) + 
  ylim(min(c(dataA$freq_both_line_resistant * dataA$Prevalence *100000,
             dataB$freq_both_line_resistant * dataB$Prevalence *100000, 
             dataC$freq_both_line_resistant * dataC$Prevalence *100000)),
       max(c(dataA$freq_both_line_resistant * dataA$Prevalence *100000,
             dataB$freq_both_line_resistant * dataB$Prevalence *100000, 
             dataC$freq_both_line_resistant * dataC$Prevalence *100000))) + 
  scale_x_discrete(breaks = seq(2020, 2040, by = 5)) +
  theme(axis.title= element_text(size=15),
        legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=15), #change legend text font sizeaxis.title= element_text(size=14),
        axis.text.x = element_text(size=15, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15))

o <- ggplot(dataC, aes(x=time, y= freq_both_line_resistant  * Prevalence *100000 , group=scenario)) +
  geom_line(aes(color=scenario),size = 1) + 
  theme_linedraw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) + 
  scale_x_discrete(breaks = seq(2020, 2040, by = 5)) +
  ylim(min(c(dataA$freq_both_line_resistant * dataA$Prevalence *100000,
             dataB$freq_both_line_resistant * dataB$Prevalence *100000, 
             dataC$freq_both_line_resistant * dataC$Prevalence *100000)),
       max(c(dataA$freq_both_line_resistant * dataA$Prevalence *100000,
             dataB$freq_both_line_resistant * dataB$Prevalence *100000, 
             dataC$freq_both_line_resistant * dataC$Prevalence *100000))) + 
  theme(axis.title= element_text(size=15),
        legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=15), #change legend text font sizeaxis.title= element_text(size=14),
        axis.text.x = element_text(size=15, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15))

annotate_figure(ggarrange(m,n,o, 
                          labels = c("A - mutation rate 0.005", "B - mutation rate 0.01", "C - mutation rate 0.015"),
                          font.label = list(size = 14, face = "bold"),
                          common.legend = T,
                          vjust=0,
                          hjust=-0.3,
                          legend = "right",
                          ncol = 3, nrow = 1),
                bottom = text_grob("Year",size = 20, hjust = 3),
                top = text_grob("", size = 20),
                left = text_grob("Frequency resistant both Antibiotics", rot=90, size = 20),
                
)

p <- ggplot(dataA, aes(x=time, y= freq_resistant, group=scenario)) +
  geom_line(aes(color=scenario)) + 
  theme_linedraw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) + 
  scale_x_discrete(breaks = seq(2022, 2040, by = 5)) +
  ylim(min(c(dataA$freq_resistant,dataB$freq_resistant, dataC$freq_resistant)),
       max(c(dataA$freq_resistant,dataB$freq_resistant, dataC$freq_resistant))) +
  theme(axis.title= element_text(size=15),
        legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=15), #change legend text font sizeaxis.title= element_text(size=14),
        axis.text.x = element_text(size=15, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15))

q <- ggplot(dataB, aes(x=time, y= freq_resistant, group=scenario)) +
  geom_line(aes(color=scenario)) + 
  theme_linedraw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) + 
  scale_x_discrete(breaks = seq(2022, 2040, by = 5)) +
  ylim(min(c(dataA$freq_resistant,dataB$freq_resistant, dataC$freq_resistant)),
       max(c(dataA$freq_resistant,dataB$freq_resistant, dataC$freq_resistant))) +
  theme(axis.title= element_text(size=15),
        legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=15), #change legend text font sizeaxis.title= element_text(size=14),
        axis.text.x = element_text(size=15, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15))


r <- ggplot(dataC, aes(x=time, y= freq_resistant, group=scenario)) +
 geom_line(aes(color=scenario)) + 
  theme_linedraw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) + 
  scale_x_discrete(breaks = seq(2022, 2040, by = 5)) +
  ylim(min(c(dataA$freq_resistant,dataB$freq_resistant, dataC$freq_resistant)),
       max(c(dataA$freq_resistant,dataB$freq_resistant, dataC$freq_resistant))) +
  theme(axis.title= element_text(size=15),
        legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=15), #change legend text font sizeaxis.title= element_text(size=14),
        axis.text.x = element_text(size=15, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=15))


annotate_figure(ggarrange(p,q,r,
                          labels = c("A", "B", "C"),
                          font.label = list(size = 14, face = "bold"),
                          common.legend = T,
                          vjust=45,
                          legend = "right",
                          ncol = 3, nrow = 1),
                bottom = text_grob("Time",size = 20),
                left = text_grob("freq resistant", rot=90, size = 20),
                
)
