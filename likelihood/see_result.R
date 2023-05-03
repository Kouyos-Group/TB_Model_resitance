# Remove current working environment
rm(list = ls())
# Necessary libraries
if(!requireNamespace("pacman", quietly = T))
  install.packages("pacman")
pacman::p_load("deSolve", "chron", "tidyverse", "deSolve", "tidyverse", "stringr", "dplyr", "ggplot2")

# Read the data (and set the working directory using setwd() )
data1 <- read.csv2("/home/louis/Bureau/stage/script/likelihood/data/TB_burden_countries_2023-02-22.csv", sep=",")
data1 <-subset(data1, data1[,1]== "Georgia")
data2<-read.csv2("/home/louis/Bureau/stage/script/likelihood/data/MDR_RR_TB_burden_estimates_2023-02-22.csv", sep=",")
data2<-subset(data2, data2[,1]== "Georgia")


data<-cbind(data1[,c(6,8, 35)], c(rep(NA, each=nrow(data1)-nrow(data2)), data2[,8]))
colnames(data)<-c("time","Incidence", "death", "Rresist")
data[,2]<-as.integer(data[,2])
data[,3]<-as.integer(data[,3])
data[,4]<-as.numeric(data[,4])

#### PRIOR ####
#All parameters
#All parameters
fixed <- c( 
  #birth
  pi = 0.0008991667, 
  
  #rate for infection (alpha)
  infection = 7*10^-6, 
  resistance_cost=0.22,
  
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
  rho= 0.036, 
  
  #death due to tuberculosis
  theta= 0.01667,  
  
  #lambda
  l=0.0175,
  
  #natural death
  mu= 0.0008966667, 
  
  #relapse
  phi = 0.0008, 
  
  #number in differents compartments
  e = 0.1, 
  late = 0.6,
  i = 0.1, 
  t = 0.1, 
  r = 0.1, 
  
  #part of resistance in the pop
  resist = 0.02, 
  treat = 0.019
)

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

#first period (introduction rifampicin)
source("/home/louis/Bureau/stage/script/T_T/likelihood/scenarios/scenario_2000-2015.R")
times <- seq(0, 15*12, by=1)
simulation <- as.data.frame(ode(init,times,sir_equations,fixed))
simulation <- simulation[-1,]

#second period (introduction bedaqualine)
init<-tail(simulation, n=1)[2:33]
init2<-unlist(init)
names(init2)<-names(init)
init<-init2
times <- seq(15*12, 22*12, by=1)
source("/home/louis/Bureau/stage/script/T_T/likelihood/scenarios/scenario_2015-2021.R")
simulation2 <-  as.data.frame(ode(init,times,sir_equations,fixed))
simulation2 <- simulation2[-1,]

#keep only the data that correspond to the data of WHO
simulation<- rbind(simulation, simulation2)

#compare simulate data and real data
Incidence_prior<-as.vector(rowsum(simulation[,31], rep(1:22, each=12)))/as.vector(rowsum(rowSums(simulation[,2:30]), rep(1:22, each=12)))
dead_prior <-as.vector(rowsum(simulation[,32], rep(1:22, each=12)))/as.vector(rowsum(rowSums(simulation[,2:30]), rep(1:22, each=12)))
rR_prior<-as.vector(rowsum(simulation[,33], rep(1:22, each=12)))/as.vector(rowsum(simulation[,31], rep(1:22, each=12)))



#### sd ####
#conditions
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
  l=exp(fit$par["l"]),
  
  #natural death
  mu= 0.0008966667, 
  
  #relapse
  phi = 0.0008, 
  
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


#first period (introduction rifampicin)
source("/home/louis/Bureau/stage/script/T_T/likelihood/scenarios/scenario_2000-2015.R")
times <- seq(0, 15*12, by=1)
simulation <- as.data.frame(ode(init,times,sir_equations,fixed))
simulation <- simulation[-1,]

#second period (introduction bedaqualine)
init<-tail(simulation, n=1)[2:33]
init2<-unlist(init)
names(init2)<-names(init)
init<-init2
times <- seq(15*12, 22*12, by=1)
source("/home/louis/Bureau/stage/script/T_T/likelihood/scenarios/scenario_2015-2021.R")
simulation2 <-  as.data.frame(ode(init,times,sir_equations,fixed))
simulation2 <- simulation2[-1,]

#keep only the data that correspond to the data of WHO
simulation<- rbind(simulation, simulation2)


#compare simulate data and real data
Incidence_sd<-as.vector(rowsum(simulation[,31], rep(1:22, each=12)))/as.vector(rowsum(rowSums(simulation[,2:30]), rep(1:22, each=12)))
dead_sd <-as.vector(rowsum(simulation[,32], rep(1:22, each=12)))/as.vector(rowsum(rowSums(simulation[,2:30]), rep(1:22, each=12)))
rR_sd<-as.vector(rowsum(simulation[,33], rep(1:22, each=12)))/as.vector(rowsum(simulation[,31], rep(1:22, each=12)))

#### sd scale ####
#conditions
fit<-readRDS("~/Bureau/stage/script/T_T/result/sd_scale/fit.RData")

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
  l=exp(fit$par["l"]),
  
  #natural death
  mu= 0.0008966667, 
  
  #relapse
  phi = 0.0008, 
  
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
#first period (introduction rifampicin)
source("/home/louis/Bureau/stage/script/T_T/likelihood/scenarios/scenario_2000-2015.R")
times <- seq(0, 15*12, by=1)
simulation <- as.data.frame(ode(init,times,sir_equations,fixed))
simulation <- simulation[-1,]

#second period (introduction bedaqualine)
init<-tail(simulation, n=1)[2:33]
init2<-unlist(init)
names(init2)<-names(init)
init<-init2
times <- seq(15*12, 22*12, by=1)
source("/home/louis/Bureau/stage/script/T_T/likelihood/scenarios/scenario_2015-2021.R")
simulation2 <-  as.data.frame(ode(init,times,sir_equations,fixed))
simulation2 <- simulation2[-1,]

#keep only the data that correspond to the data of WHO
simulation<- rbind(simulation, simulation2)


#compare simulate data and real data
Incidence_sd_scale<-as.vector(rowsum(simulation[,31], rep(1:22, each=12)))/as.vector(rowsum(rowSums(simulation[,2:30]), rep(1:22, each=12)))
dead_sd_scale <-as.vector(rowsum(simulation[,32], rep(1:22, each=12)))/as.vector(rowsum(rowSums(simulation[,2:30]), rep(1:22, each=12)))
rR_sd_scale<-as.vector(rowsum(simulation[,33], rep(1:22, each=12)))/as.vector(rowsum(simulation[,31], rep(1:22, each=12)))

#plot
#Incidence
Incidence<-data.frame(rbind(cbind(time = data$time, incidence = data$Incidence/100000, type = "real data"), 
                            cbind(time = data$time, incidence = round(Incidence_prior, 10), type = "simulation prior"),
                            cbind(time = data$time, incidence = round(Incidence_sd, 10), type = "simulation sd"),
                            cbind(time = data$time, incidence = round(Incidence_sd_scale, 4), type = "simulation sd scale")
                            )
                      )

Incidence$incidence<-as.numeric(Incidence$incidence)


#Death
death<-data.frame(rbind(cbind(time = data$time, death = data$death/100000, type = "real data"), 
                        cbind(time = data$time, death = round(dead_prior, 10), type = "simulation prior"),
                        cbind(time = data$time, death = round(dead_sd, 10), type = "simulation sd"),
                        cbind(time = data$time, death = round(dead_sd_scale, 4), type = "simulation sd scale")
)
)

death$death<-as.numeric(death$death)

#Resitance to the rifampcin
rifampicin_resistance <-data.frame(rbind(cbind(time = data$time, "proportion rifampicin resistance" = data$Rresist/100, type = "real data"), 
                        cbind(time = data$time, "proportion rifampicin resistance" = rR_prior, type = "simulation prior"),
                        cbind(time = data$time, "proportion rifampicin resistance" = rR_sd, type = "simulation sd"),
                        cbind(time = data$time, "proportion rifampicin resistance" = rR_sd_scale, type = "simulation sd scale")
)
)

rifampicin_resistance$proportion.rifampicin.resistance<-as.numeric(rifampicin_resistance$proportion.rifampicin.resistance)

########################################################
ggplot(rifampicin_resistance, aes(x=time, y=proportion.rifampicin.resistance, group=type)) +
  geom_line(aes(color=type)) + 
  xlab("time (years)") + ylab("proportion rifampicin resistance")

ggplot(Incidence, aes(x=time, y=incidence, group=type)) +
  geom_line(aes(color=type)) + 
  xlab("time (years)")

ggplot(death, aes(x=time, y=death, group=type)) +
  geom_line(aes(color=type)) + 
  xlab("time (years)")


