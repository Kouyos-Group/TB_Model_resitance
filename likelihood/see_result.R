# Remove current working environment
rm(list = ls())
# Necessary libraries
if(!requireNamespace("pacman", quietly = T))
  install.packages("pacman")
pacman::p_load("deSolve", "chron", "tidyverse", "deSolve", "tidyverse", "stringr", "dplyr", "ggplot2", "ggpubr", "gtExtras")

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
fixed <- c(
  #birth
  pi = 0.0008991667,
  
  #rate for infection (alpha)
  infection = 7*10^-5,
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
  l=0.3,
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
fit<-readRDS("~/Bureau/stage/script/T_T/result/Georgia/rho_low/sd_scale/fit.RData")

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

### no normalize ####
#conditions
fit<-readRDS("~/Bureau/stage/script/T_T/result/Georgia/rho_low/no_normalize/fit.RData")

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
Incidence_no_normalize <- as.vector(rowsum(simulation[,31], rep(1:22, each=12)))/as.vector(rowsum(rowSums(simulation[,2:30]), rep(1:22, each=12)))
dead_no_normalize <- as.vector(rowsum(simulation[,32], rep(1:22, each=12)))/as.vector(rowsum(rowSums(simulation[,2:30]), rep(1:22, each=12)))
rR_no_normalize <- as.vector(rowsum(simulation[,33], rep(1:22, each=12)))/as.vector(rowsum(simulation[,31], rep(1:22, each=12)))

#plot
#Incidence
Incidence <-data.frame(time = data$time, 
                   "real data" = data$Incidence/100000,
                   "simulation sd" = Incidence_sd)#,
                   #"simulation sd scale"= Incidence_sd_scale,
                   #"simulation no sd" = Incidence_no_normalize)


#Death
Death <-data.frame(time = data$time, 
                  "real data" = data$death/100000,
                   "simulation sd" = dead_sd)#,
                   #"simulation sd scale"= dead_sd_scale,
                   #"simulation no sd" = dead_no_normalize)


#Resitance to the rifampcin
rifampicin_resistance <-data.frame(time = data$time, 
                                   "real data" = data$Rresist/100,
                                   "simulation sd" = rR_sd)#,
                                   #"simulation sd scale"= rR_sd_scale,
                                   #"simulation no sd" = rR_no_normalize)
                        
rifampicin_resistance$time <- as.factor(rifampicin_resistance$time)
Incidence$time <- as.factor(Incidence$time)
Death$time <- as.factor(Death$time)

########################################################


# Créer le graphique avec ggplot
a <- ggplot() +
  # Ajouter les lignes y_lignes
  geom_line(data = rifampicin_resistance, aes(x = time, y = simulation.sd * 100,color = "simulation", group = "simulation"),
            size = 1.5) +
  #geom_line(data = rifampicin_resistance, aes(x = time, y = simulation.sd.scale, color = "simulation sd scale")) +
  #geom_line(data = rifampicin_resistance, aes(x = time, y = simulation.no.sd, color = "simulation no sd")) +
  # Ajouter les points y_points
  geom_point(data = rifampicin_resistance, aes(x = time, y = real.data * 100,color = "real data", group = "real data")) + 
  xlab("Year") + 
  theme(legend.position = "none") + 
  ylab("Proportion Rifampicin resistant per 100 infected") +
  scale_x_discrete(breaks = seq(2000, 2022, by = 10)) +
  theme_linedraw() + 
  theme(axis.title= element_text(size=15),
        legend.title = element_text(size=20), #change legend title font size
        legend.text = element_text(size=15), #change legend text font sizeaxis.title= element_text(size=14),
        axis.text.x = element_text(size=20, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=20))


                            
# Créer le graphique avec ggplot
b <- ggplot() +
  # Ajouter les lignes y_lignes
  geom_line(data = Incidence, aes(x = time, y = simulation.sd * 1000000, group = "simulation", color = "simulation"), size = 1.5) +
  #geom_line(data = Incidence, aes(x = time, y = simulation.sd.scale, color = "simulation sd scale")) +
  #geom_line(data = Incidence, aes(x = time, y = simulation.no.sd, color = "simulation no sd")) +
  # Ajouter les points y_points
  geom_point(data = Incidence, aes(x = time, y = real.data * 1000000, group = "real data", color = "real data")) + 
  xlab("Year") + 
  theme(legend.position = "none") + 
  ylab("Incidence per 100,000 individuals") +
  scale_x_discrete(breaks = seq(2000, 2022, by = 10)) +
  theme_linedraw() + 
  theme(axis.title= element_text(size=15),
        legend.title = element_text(size=20), #change legend title font size
        legend.text = element_text(size=15), #change legend text font sizeaxis.title= element_text(size=14),
        axis.text.x = element_text(size=20, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=20))

# Créer le graphique avec ggplot
c <- ggplot() +
  # Ajouter les lignes y_lignes
  #
  geom_line(data = Death, aes(x = time, y = simulation.sd * 100000, group = "simulation", color = "simulation"), size = 1.5) +
  #geom_line(data = Death, aes(x = time, y = simulation.sd.scale, color = "simulation sd scale")) +
  #geom_line(data = Death, aes(x = time, y = simulation.no.sd, color = "simulation no sd")) +
  # Ajouter les points y_points
  geom_point(data = Death, aes(x = time, y = real.data*100000,group = "real data", color = "real data")) + 
  xlab("Year") + 
  theme(legend.position = "none") + 
  ylab("Death due to TB per 100,000 indiviuduals") +
  scale_x_discrete(breaks = seq(2000, 2022, by = 10)) +
  theme_linedraw()+ 
  theme(axis.title= element_text(size=15),
        legend.title = element_text(size=20), #change legend title font size
        legend.text = element_text(size=15), #change legend text font sizeaxis.title= element_text(size=14),
        axis.text.x = element_text(size=20, angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=20))


ggarrange(b, c, a, 
          labels = c("A", "B", "C"),
          legend = "right",
          common.legend = T,
          font.label = list(size = 15),
          ncol = 3, nrow = 1)


cost_no_normalize <- sum(Incidence_no_normalize- data$Incidence)^2 + sum(dead_no_normalize-data$death)^2 + sum(rR_no_normalize[16:22] - data$Rresist[16:22])^2


error_incidence <- (Incidence_sd-data$incidence/100000) / sd(data$incidence/100000)
error_death <- (dead_sd-data$death/100000) / sd(data$death/100000)
error_Rresist <- (rR_sd[16:22]-data$Rresist[16:22]/100) / sd(data$death[16:22]/100)

cost_sd <- sum(error_incidence^2) + sum(error_death^2) + sum(error_Rresist^2)

#####Calcul of sd due to scale #####
data$Incidence
### Incidence ###
total_var_incidence <- var(data$Incidence/100000)

# Régression linéaire pour calculer la variance expliquée
fit_incidence <- lm(data$Incidence/100000 ~ data$time)
exp_var_incidence <- summary(fit_incidence)$r.squared * total_var_incidence

# Calcul de la variance résiduelle
res_var_incidence <- total_var_incidence - exp_var_incidence

# Calcul de la part de la variance totale due au bruit et à l'échelle
scale_part_incidence <- exp_var_incidence / total_var_incidence

sd_scale_incidence <- sd(data$Incidence/100000) * sqrt(scale_part_incidence)

### Death ###
total_var_death <- var(data$death/100000)

# Régression linéaire pour calculer la variance expliquée
fit_death <- lm(data$death/100000 ~ data$time)
exp_var_death <- summary(fit_death)$r.squared * total_var_death

# Calcul de la variance résiduelle
res_var_death <- total_var_death - exp_var_death

# Calcul de la part de la variance totale due au bruit et à l'échelle
scale_part_death <- exp_var_death / total_var_death

sd_scale_death <- sd(data$death/100000) * sqrt(scale_part_death)

### Incidence resistance ###
total_var_Rresist <- var(data$Rresist[16:22]/100)

# Régression linéaire pour calculer la variance expliquée
fit_Rresist <- lm(data$Rresist[16:22]/100 ~ data$time[16:22])
exp_var_Rresist <- summary(fit_Rresist)$r.squared * total_var_Rresist

# Calcul de la variance résiduelle
res_var_Rresist <- total_var_Rresist - exp_var_Rresist

# Calcul de la part de la variance totale due au bruit et à l'échelle
scale_part_Rresist <- exp_var_Rresist / total_var_Rresist

sd_scale_Rresist <- sd(data$Rresist[16:22]/100) * sqrt(scale_part_Rresist)


error_incidence <- (Incidence_sd_scale-data$Incidence/100000) / sd_scale_incidence
error_death <- (dead_sd_scale-data$death/100000) / sd_scale_death
error_Rresist <- (rR_sd_scale[16:22]-data$Rresist[16:22]/100) / sd_scale_Rresist

cost_sd_scale <- sum(error_incidence^2) + sum(error_death^2) + sum(error_Rresist^2)

data <- data.frame(
  x = c("not standardized", "standardized", "standardized by the scale"),
  y = c(cost_no_normalize, cost_sd, cost_sd_scale)
)

# Horizontal version
ggplot(data, aes(x=x, y=log(y))) +
  geom_segment( aes(x=x, xend=x, y=0, yend=log(y)), color="skyblue") +
  geom_point( color="blue", size=4, alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank()
  )

