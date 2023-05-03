# Remove current working environment
rm(list = ls())
# Necessary libraries
if(!requireNamespace("pacman", quietly = T))
  install.packages("pacman")
pacman::p_load("deSolve", "chron", "tidyverse", "deSolve", "tidyverse", "stringr", "dplyr", "ggplot2")

##### 2000-2022 ####
fit<-readRDS("~/Bureau/stage/script/T_T/result/fit.RData")


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
  l = 0.4,
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
                   "freq_1st_line_resistant" = (rowSums(sir_values_0[,str_sub(names(sir_values_0), 1 , 1) == "I" & 
                                                                       (str_sub(names(sir_values_0), 3 , 3) == "1")]) /
                                                  rowSums(sir_values_0[,11:14])),
                   "freq_2nd_line_resistant" = (rowSums(sir_values_0[,str_sub(names(sir_values_0), 1 , 1) == "I" & str_sub(names(sir_values_0), 4 , 4) == "1"]) / 
                                                  rowSums(sir_values_0[,11:14])), 
                   "Incidence" = sir_values_0[,31],
                   "scenario" = rep("2000 - 2015", each = length(sir_values_0[,1])))

init<-tail(simulation, n=1)[2:33]
init2<-unlist(init)
names(init2)<-names(init)
init<-c(init2, nT = 0)
times <- seq(22 * 12, 40 *12 , by=1/3)

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
                  "number_sensible" = rowSums(sir_values_1[,str_sub(names(sir_values_1), 1 , 1) == "I" & 
                                                             (str_sub(names(sir_values_1), 3 , 3) == "0" | str_sub(names(sir_values_1), 4 , 4) == "0")]),
                  "number_1st_line_resistant" = rowSums(sir_values_1[,str_sub(names(sir_values_1), 1 , 1) == "I" & 
                                                                       (str_sub(names(sir_values_1), 3 , 3) == "1")]),
                  "number_2nd_line_resistant" = rowSums(sir_values_1[,str_sub(names(sir_values_1), 1 , 1) == "I" & str_sub(names(sir_values_1), 4 , 4) == "1"]),
                  "freq_1st_line_resistant" = (rowSums(sir_values_1[,str_sub(names(sir_values_1), 1 , 1) == "I" & 
                                                                      (str_sub(names(sir_values_1), 3 , 3) == "1")]) /
                                                 rowSums(sir_values_1[,11:14])),
                  "freq_2nd_line_resistant" = (rowSums(sir_values_1[,str_sub(names(sir_values_1), 1 , 1) == "I" & str_sub(names(sir_values_1), 4 , 4) == "1"]) / 
                                                 rowSums(sir_values_1[,11:14])), 
                  "Incidence" = sir_values_1[,31],
                  "scenario" = rep("standard strategy", each = length(sir_values_1[,1])))

##### TRUNCATE strategy ####
init_T <- c(init, 
            Ts00t4 = 0 , Ts10t4 = 0, Ts01t4 = 0, Ts11t4 = 0,
            Ts00t5 = 0 , Ts10t5 = 0, Ts01t5 = 0, Ts11t5 = 0)

fixed_T1 <- c(fixed, cured_t4 = 0.10 , cured_t5 = 0.10, time_4 = 1/2, time_5 = 1, rho_T = 1) 
source("~/Bureau/stage/script/T_T/model/scenarios/TRUNCATE_strategy.R")
simulation <-  as.data.frame(ode(init_T,times,sir_equations,fixed_T1))


sir_values_2<-data.frame(simulation)

data2<-data.frame("time" = sir_values_2[,1], 
                  "S" = sir_values_2[,2], 
                  "E" = rowSums(sir_values_2[,3:6]), 
                  "L"= rowSums(sir_values_2[,7:10]),
                  "I"=rowSums(sir_values_2[,11:14]), 
                  "T"= rowSums(sir_values_2[,15:26]) + rowSums(sir_values_2[,34:41]), 
                  "R" = rowSums(sir_values_2[,27:30]),
                  "Infected" = rowSums(sir_values_2[,11:14]),
                  "number_sensible" = rowSums(sir_values_2[,str_sub(names(sir_values_2), 1 , 1) == "I" & 
                                                             (str_sub(names(sir_values_2), 3 , 3) == "0" | str_sub(names(sir_values_2), 4 , 4) == "0")]),
                  "number_1st_line_resistant" = rowSums(sir_values_2[,str_sub(names(sir_values_2), 1 , 1) == "I" & 
                                                                       (str_sub(names(sir_values_2), 3 , 3) == "1")]),
                  "number_2nd_line_resistant" = rowSums(sir_values_2[,str_sub(names(sir_values_2), 1 , 1) == "I" & str_sub(names(sir_values_2), 4 , 4) == "1"]),
                  "freq_1st_line_resistant" = (rowSums(sir_values_2[,str_sub(names(sir_values_2), 1 , 1) == "I" & 
                                                                      (str_sub(names(sir_values_2), 3 , 3) == "1")]) /
                                                 rowSums(sir_values_2[,11:14])),
                  "freq_2nd_line_resistant" = (rowSums(sir_values_2[,str_sub(names(sir_values_2), 1 , 1) == "I" & str_sub(names(sir_values_2), 4 , 4) == "1"]) / 
                                                 rowSums(sir_values_2[,11:14])), 
                  "Incidence" = sir_values_2[,31],
                  "scenario" = rep("Truncate strategy 1", each = length(sir_values_2[,1])))


##### TRUNCATE strategy 2####
fixed_T2 <- c(fixed, cured_t4 = 0.10, cured_t5 = 0.10, time_4 = 1/2, time_5 = 1, rho_T = 3) 
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
                  "number_sensible" = rowSums(sir_values_3[,str_sub(names(sir_values_3), 1 , 1) == "I" & 
                                                             (str_sub(names(sir_values_3), 3 , 3) == "0" | str_sub(names(sir_values_3), 4 , 4) == "0")]),
                  "number_1st_line_resistant" = rowSums(sir_values_3[,str_sub(names(sir_values_3), 1 , 1) == "I" & 
                                                                       (str_sub(names(sir_values_3), 3 , 3) == "1")]),
                  "number_2nd_line_resistant" = rowSums(sir_values_3[,str_sub(names(sir_values_3), 1 , 1) == "I" & str_sub(names(sir_values_3), 4 , 4) == "1"]),
                  "freq_1st_line_resistant" = (rowSums(sir_values_3[,str_sub(names(sir_values_3), 1 , 1) == "I" & 
                                                                      (str_sub(names(sir_values_3), 3 , 3) == "1")]) /
                                                 rowSums(sir_values_3[,11:14])),
                  "freq_2nd_line_resistant" = (rowSums(sir_values_3[,str_sub(names(sir_values_3), 1 , 1) == "I" & str_sub(names(sir_values_3), 4 , 4) == "1"]) / 
                                                 rowSums(sir_values_3[,11:14])), 
                  "Incidence" = sir_values_3[,31],
                  "scenario" = rep("Truncate strategy 2", each = length(sir_values_3[,1])))


##### TRUNCATE strategy 3####
fixed_T3 <- c(fixed, cured_t4 = 0.10 , cured_t5 = 0.10, time_4 = 1/2, time_5 = 1, rho_T = 5) 
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
                  "number_sensible" = rowSums(sir_values_4[,str_sub(names(sir_values_4), 1 , 1) == "I" & 
                                                             (str_sub(names(sir_values_4), 3 , 3) == "0" | str_sub(names(sir_values_4), 4 , 4) == "0")]),
                  "number_1st_line_resistant" = rowSums(sir_values_4[,str_sub(names(sir_values_4), 1 , 1) == "I" & 
                                                                       (str_sub(names(sir_values_4), 3 , 3) == "1")]),
                  "number_2nd_line_resistant" = rowSums(sir_values_4[,str_sub(names(sir_values_4), 1 , 1) == "I" & str_sub(names(sir_values_4), 4 , 4) == "1"]),
                  "freq_1st_line_resistant" = (rowSums(sir_values_4[,str_sub(names(sir_values_4), 1 , 1) == "I" & 
                                                                      (str_sub(names(sir_values_4), 3 , 3) == "1")]) /
                                                 rowSums(sir_values_4[,11:14])),
                  "freq_2nd_line_resistant" = (rowSums(sir_values_4[,str_sub(names(sir_values_4), 1 , 1) == "I" & str_sub(names(sir_values_4), 4 , 4) == "1"]) / 
                                                 rowSums(sir_values_4[,11:14])), 
                  "Incidence" = sir_values_4[,31],
                  "scenario" = rep("Truncate strategy 3", each = length(sir_values_4[,1])))


##### TRUNCATE strategy 4####
fixed_T4 <- c(fixed, cured_t4 = 0.10 , cured_t5 = 0.10, time_4 = 1/2, time_5 = 1, rho_T = 10) 
source("~/Bureau/stage/script/T_T/model/scenarios/TRUNCATE_strategy.R")
simulation <-  as.data.frame(ode(init_T,times,sir_equations,fixed_T4))


sir_values_5<-data.frame(simulation)

data5<-data.frame("time" = sir_values_5[,1], 
                  "S" = sir_values_5[,2], 
                  "E" = rowSums(sir_values_5[,3:6]), 
                  "L"= rowSums(sir_values_5[,7:10]),
                  "I"=rowSums(sir_values_5[,11:14]), 
                  "T"= rowSums(sir_values_5[,15:26]) + rowSums(sir_values_5[,34:41]), 
                  "R" = rowSums(sir_values_5[,27:30]),
                  "Infected" = rowSums(sir_values_5[,11:14]),
                  "number_sensible" = rowSums(sir_values_5[,str_sub(names(sir_values_5), 1 , 1) == "I" & 
                                                             (str_sub(names(sir_values_5), 3 , 3) == "0" | str_sub(names(sir_values_5), 4 , 4) == "0")]),
                  "number_1st_line_resistant" = rowSums(sir_values_5[,str_sub(names(sir_values_5), 1 , 1) == "I" & 
                                                                       (str_sub(names(sir_values_5), 3 , 3) == "1")]),
                  "number_2nd_line_resistant" = rowSums(sir_values_5[,str_sub(names(sir_values_5), 1 , 1) == "I" & str_sub(names(sir_values_5), 4 , 4) == "1"]),
                  "freq_1st_line_resistant" = (rowSums(sir_values_5[,str_sub(names(sir_values_5), 1 , 1) == "I" & 
                                                                      (str_sub(names(sir_values_5), 3 , 3) == "1")]) /
                                                 rowSums(sir_values_5[,11:14])),
                  "freq_2nd_line_resistant" = (rowSums(sir_values_5[,str_sub(names(sir_values_5), 1 , 1) == "I" & str_sub(names(sir_values_5), 4 , 4) == "1"]) / 
                                                 rowSums(sir_values_5[,11:14])), 
                  "Incidence" = sir_values_5[,31],
                  "scenario" = rep("Truncate strategy 4", each = length(sir_values_5[,1])))



##########################################################

data<-rbind(data1, data2, data3, data4, data5)

ggplot(data, aes(x=2000 + (time/12), y= Infected, group=scenario)) +
  geom_line(aes(color=scenario)) + 
  xlab("time (years)") + ylab("Infected")

ggplot(data, aes(x=2000 + (time/12), y= number_sensible, group=scenario)) +
  geom_line(aes(color=scenario)) + 
  xlab("time (years)") + ylab("number sensible")

ggplot(data, aes(x=2000 + (time/12), y= number_1st_line_resistant, group=scenario)) +
  geom_line(aes(color=scenario)) + 
  xlab("time (years)") + ylab("number_1st_line_resistant")

ggplot(data, aes(x=2000 + (time/12), y= number_2nd_line_resistant, group=scenario)) +
  geom_line(aes(color=scenario)) + 
  xlab("time (years)") + ylab("number_2nd_line_resistant")

for (i in 1:length(fixed)){ assign(names(fixed[i]),fixed[i])}
state <- c()
for (i in 1:length(init)){ state <- c(state, assign(names(init[i]),1))}
names(state) <- names(init)

for (i in 1:length(fixed_T1)){ assign(names(fixed_T1[i]),fixed_T1[i])}
state <- c()
for (i in 1:length(init_T)){ state <- c(state, assign(names(init_T[i]),1))}
names(state) <- names(init_T)


plot(sir_values_1$time,sir_values_1$Ts00t1, col = "red", type = "l", ylim = c(0, max(sir_values_1$Ts00t1)))
lines(sir_values_2$time,sir_values_2$Ts00t1 + sir_values_2$Ts00t4, col = "blue", type = "l")
lines(sir_values_2$time,sir_values_2$Ts00t4, col = "blue", type = "l")

plot(sir_values_1$time,sir_values_1$Ts10t1, col = "red", type = "l", ylim = c(0, max(sir_values_1$Ts00t1)))
lines(sir_values_2$time,sir_values_2$Ts01t1, col = "blue", type = "l")
lines(sir_values_2$time,sir_values_2$Ts01t4, col = "blue", type = "l")
lines(sir_values_2$time,sir_values_2$Ts01t5, col = "blue", type = "l")

sir_values_1 <- sir_values_2

plot(sir_values_1$time,sir_values_1$S, type="l", ylim = c(0,100000), col="green")
lines(sir_values_1$time, sir_values_1$E, type = "l", lty = 1, col="red")
lines(sir_values_1$time, sir_values_1$L, type = "l", lty = 1, col="blue")
lines(sir_values_1$time, sir_values_1$I, type = "l", lty = 1, col="orange")
lines(sir_values_1$time, sir_values_1$T, type = "l", lty = 1, col="pink")
lines(sir_values_1$time, sir_values_1$R, type = "l", lty = 1, col="purple")

ggplot(data2)

sir_values_1 <- sir_values_2
plot(sir_values_1$time, sir_values_1$Ts11t2, type = "l", lty = 1, col="orange")
