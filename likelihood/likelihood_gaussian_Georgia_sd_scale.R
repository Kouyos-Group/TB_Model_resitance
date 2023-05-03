# Remove current working environment
rm(list = ls())
# Necessary libraries
if(!requireNamespace("pacman", quietly = T))
  install.packages("pacman")
pacman::p_load("deSolve", "chron", "tidyverse", "deSolve", "tidyverse", "stringr", "dplyr", "ggplot2")


#### Load data ####
#load
data1 <- read.csv2("/home/locoll/data/data/TB_burden_countries_2023-02-22.csv", sep=",")
data1 <- subset(data1, data1[,1] == "Georgia")
data2 <- read.csv2("/home/locoll/data/data/MDR_RR_TB_burden_estimates_2023-02-22.csv", sep=",")
data2 <- subset(data2, data2[,1] == "Georgia")

#structure the data (here we keep incidence (incidence), death due to Tb (death) and proportion of infection with rifampicin strain (Rresist))
data <- cbind(data1[,c(6,8, 35)], c(rep(NA, each=nrow(data1)-nrow(data2)), data2[,8]))
colnames(data) <- c("times","incidence", "death", "Rresist")
data[,2] <- as.integer(data[,2])
data[,3] <- as.integer(data[,3])
data[,4] <- as.numeric(data[,4])

######Calcul of sd due to scale #####
### Incidence ###
total_var_incidence <- var(data$incidence/100000)

# Régression linéaire pour calculer la variance expliquée
fit_incidence <- lm(data$incidence/100000 ~ data$times)
exp_var_incidence <- summary(fit_incidence)$r.squared * total_var_incidence

# Calcul de la variance résiduelle
res_var_incidence <- total_var_incidence - exp_var_incidence

# Calcul de la part de la variance totale due au bruit et à l'échelle
scale_part_incidence <- exp_var_incidence / total_var_incidence

sd_scale_incidence <- sd(data$incidence/100000) * sqrt(scale_part_incidence)

### Death ###
total_var_death <- var(data$death/100000)

# Régression linéaire pour calculer la variance expliquée
fit_death <- lm(data$death/100000 ~ data$times)
exp_var_death <- summary(fit_death)$r.squared * total_var_death

# Calcul de la variance résiduelle
res_var_death <- total_var_death - exp_var_death

# Calcul de la part de la variance totale due au bruit et à l'échelle
scale_part_death <- exp_var_death / total_var_death

sd_scale_death <- sd(data$death/100000) * sqrt(scale_part_death)

### Incidence resistance ###
total_var_Rresist <- var(data$Rresist[16:22]/100)

# Régression linéaire pour calculer la variance expliquée
fit_Rresist <- lm(data$Rresist[16:22]/100 ~ data$times[16:22])
exp_var_Rresist <- summary(fit_Rresist)$r.squared * total_var_Rresist

# Calcul de la variance résiduelle
res_var_Rresist <- total_var_Rresist - exp_var_Rresist

# Calcul de la part de la variance totale due au bruit et à l'échelle
scale_part_Rresist <- exp_var_Rresist / total_var_Rresist

sd_scale_Rresist <- sd(data$Rresist[16:22]/100) * sqrt(scale_part_Rresist)



#### Build Function cost for optim ####
# Cost function to calculate the sum of squared residuals (SSR) of predicted vs observed data
cost <- function (free,fixed,data) {
  
  
  #conditions
  pars <- c(free,fixed)
  pars <- trans(pars)
  init <- c(S = 70000, 
            
            Es00 = 30000 * pars["e"] * (1 - pars["resist"]) / (pars["e"] + pars["late"] + pars["i"] + pars["t"] + pars["r"])  , 
            Es10 = 30000 * pars["e"] * pars["resist"] / (pars["e"] + pars["late"] + pars["i"] + pars["t"] + pars["r"]), 
            Es01 = 0, 
            Es11 = 0,
            
            
            Ls00 = 30000 * pars["late"] * (1 - pars["resist"]) / (pars["e"] + pars["late"] + pars["i"] + pars["t"] + pars["r"]), 
            Ls10 = 30000 * pars["late"] * pars["resist"] / (pars["e"] + pars["late"] + pars["i"] + pars["t"] + pars["r"]),
            Ls01 = 0, 
            Ls11 = 0, 
            
            
            Is00 = 30000 * pars["i"] * (1 - pars["resist"]) / (pars["e"] + pars["late"] + pars["i"] + pars["t"] + pars["r"]), 
            Is10 = 30000 * pars["i"] * pars["resist"] / (pars["e"] + pars["late"] + pars["i"] + pars["t"] + pars["r"]), 
            Is01 = 0, 
            Is11 = 0, 
            
            Ts00t1 = 30000 * pars["t"] * (1 - pars["resist"]) * (1 - pars["treat"]) / (pars["e"] + pars["late"] + pars["i"] + pars["t"] + pars["r"]), 
            Ts10t1 = 30000 * pars["t"] * pars["resist"] * (1 - pars["treat"]) / (pars["e"] + pars["late"] + pars["i"] + pars["t"] + pars["r"]) , 
            Ts01t1 = 0, 
            Ts11t1 = 0, 
            Ts00t2 = 0, 
            Ts10t2 = 0, 
            Ts01t2 = 0, 
            Ts11t2 = 0, 
            Ts00t3 = 30000 * pars["t"] * (1 - pars["resist"]) * pars["treat"] / (pars["e"] + pars["late"] + pars["i"] + pars["t"] + pars["r"]), 
            Ts10t3 = 30000 * pars["t"] * pars["resist"] * pars["treat"] / (pars["e"] + pars["late"] + pars["i"] + pars["t"] + pars["r"]), 
            Ts01t3 = 0, 
            Ts11t3 = 0, 
            
            Rs00 = 30000 * pars["r"] * (1 - pars["resist"]) / (pars["e"] + pars["late"] + pars["i"] + pars["t"] + pars["r"]), 
            Rs10 = 30000 * pars["r"] * pars["resist"] / (pars["e"] + pars["late"] + pars["i"] + pars["t"] + pars["r"]), 
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
  
  #Check that there are 100,000 individuals at the beginning of the simulation
  if (any(init<0)){
    cost <- Inf
    return(cost)
  } else {  
    #first period (introduction rifampicin)
    source("/home/locoll/data/script/R/scenario_1971-2015.R")
    times <- seq(0, 15*12, by=1)
    simulation <- as.data.frame(ode(init,times,sir_equations,parms=pars))
    simulation <- simulation[-1,]
    
    #second period (introduction bedaqualine)
    init<-tail(simulation, n=1)[2:33]
    init2<-unlist(init)
    names(init2)<-names(init)
    init<-init2
    times <- seq(15*12, 22*12, by=1)
    source("/home/locoll/data/script/R/scenario_2015-2021.R")
    simulation2 <-  as.data.frame(ode(init,times,sir_equations,parms=pars))
    simulation2 <- simulation2[-1,]
    
    #keep only the data that correspond to the data of WHO
    simulation<- rbind(simulation, simulation2)
    #simulation <- tail(simulation, n = 22*12)
    
  if (nrow(simulation) == 22*12) {
    #compare simulate data and real data
    Incidence<-as.vector(rowsum(simulation[,31], rep(1:22, each=12)))/as.vector(rowsum(rowSums(simulation[,2:30]), rep(1:22, each=12)))
    dead <-as.vector(rowsum(simulation[,32], rep(1:22, each=12)))/as.vector(rowsum(rowSums(simulation[,2:30]), rep(1:22, each=12)))
    rR<-as.vector(rowsum(simulation[,33], rep(1:22, each=12)))/as.vector(rowsum(simulation[,31], rep(1:22, each=12)))
    rR<-rR[16:22]
    
    error_incidence <- (Incidence-data$incidence/100000) / sd_scale_incidence
    error_death <- (dead-data$death/100000) / sd_scale_death
    error_Rresist <- (rR-data$Rresist[16:22]/100) / sd_scale_Rresist
    
    cost <- sum(error_incidence^2) + sum(error_death^2) + sum(error_Rresist^2)
  
  } else {
    cost <- Inf
    return(cost)
  }
  }  
}

#Parameters we want to test
trans <- function(pars) {
  pars["infection"] <- exp(pars["infection"])
  pars["rho"] <- exp(pars["rho"])
  pars["resistance_cost"] <- exp(pars["resistance_cost"])
  pars["l"] <- exp(pars["l"])
  pars["theta"] <- exp(pars["theta"])
  pars["e"] <- exp(pars["e"])
  pars["late"] <- exp(pars["late"])
  pars["i"] <- exp(pars["i"])
  pars["t"] <- exp(pars["t"])
  pars["r"] <- exp(pars["r"])
  pars["resist"] <- exp(pars["resist"])
  pars["treat"] <- exp(pars["treat"])
  return(pars)
}

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
  t = 0.05, 
  r = 0.05, 
  
  #part of resistance in the pop
  resist = 0.02, 
  treat = 0.019
)


#parameters we want to test
free <- c(log(fixed["infection"]), 
          log(fixed["resistance_cost"]),
          log(fixed["l"]),
          log(fixed["rho"]),
          log(fixed["theta"]),
          log(fixed["e"]), 
          log(fixed["late"]), 
          log(fixed["i"]), 
          log(fixed["t"]),
          log(fixed["r"]), 
          log(fixed["resist"]),
          log(fixed["treat"]) 
)

#remove to fixed
fixed<-fixed[ ! names(fixed) %in%  names(free) ]

#### Apply function optim ####
#count the simulation time (1)
start.time <- Sys.time()
# Fit the model to the data
fit <- optim(free,cost,gr=NULL,fixed,data,method="Nelder-Mead",hessian=TRUE)
#count the simulation time (2)
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken
#print fit output
fit 
saveRDS(fit, file="fit.RData")

#### Apply function hessian ####
#start.time <- Sys.time()
#hess <- hessian(fit, func=cost, free,fixed,init, data , method="Richardson")
#end.time <- Sys.time()
#time.taken <- round(end.time - start.time,2)
#time.taken
#hess
#saveRDS(hess, file="hessian.RData")


#### creates a table with the estimated values at 95% confidence ####
# create table that summarizes the results:
cols <- c("Estimate","2.5%","97.5%")
rows <- names(fit$par)
est <- as.data.frame(matrix(NA,nrow=length(rows),ncol=length(cols)))
rownames(est) <- rows
colnames(est) <- cols

# CAVE model estimates are on log scale, so we have to exp them (and round)
# start with first column:
est[1] <- round(exp(fit$par),8)
# define 95%CI
up <- fit$par + 1.96*sqrt(diag(fit$hessian))
lo <- fit$par - 1.96*sqrt(diag(fit$hessian))
# second and third columns:
est[2] <- round(exp(lo), 8)
est[3] <- round(exp(up),8)

est

write.table(est, file = "estimate_95.txt", sep = ",", quote = FALSE, row.names = F)
