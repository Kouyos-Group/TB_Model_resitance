# Model by Christian Althaus
# Remove current working environment
rm(list = ls())# Necessary libraries
if(!requireNamespace("pacman", quietly = T))
  install.packages("pacman")
pacman::p_load("deSolve", "chron", "tidyverse")

# Read the data (and set the working directory using setwd() )
data1 <- read.csv2("/home/louis/Bureau/stage/script/likelihood/data/TB_burden_countries_2023-02-22.csv", sep=",")
data1 <-subset(data1, data1[,1]== "South Africa")
data2<-read.csv2("/home/louis/Bureau/stage/script/likelihood/data/MDR_RR_TB_burden_estimates_2023-02-22.csv", sep=",")
data2<-subset(data2, data2[,1]== "South Africa")



data<-cbind(data1[,c(6,8, 35)], c(rep(NA, each=nrow(data1)-nrow(data2)), data2[,8]))
colnames(data)<-c("times","incidence", "death", "Rresist")
data[,2]<-as.integer(data[,2])
data[,3]<-as.integer(data[,3])
data[,4]<-as.numeric(data[,4])


# Cost function to calculate the sum of squared residuals (SSR) of predicted vs observed data
cost <- function(free,fixed,data) {
  
  #conditions
  pars <- c(free,fixed)
  pars <- trans(pars)
  init <- c("S" = 100000 - pars["E"] - pars["L"] - pars["I"], 
            Es0000 = pars["E"], Es1000 = 0, Es0100 = 0, Es1100 = 0, Es0010 = 0, Es1010 = 0,Es0110 = 0,Es1110 = 0,
            Es0001 = 0, Es1001 = 0, Es0101 = 0, Es1101 = 0, Es0011 = 0, Es1011 = 0,Es0111 = 0,Es1111 = 0,
            
            
            Ls0000 =  pars["L"], Ls1000 = 0, Ls0100 = 0, Ls1100 = 0, Ls0010 = 0, Ls1010 = 0,Ls0110 = 0,Ls1110 = 0,
            Ls0001 = 0, Ls1001 = 0, Ls0101 = 0, Ls1101 = 0, Ls0011 = 0, Ls1011 = 0,Ls0111 = 0,Ls1111 = 0,
            
            
            Is0000 = pars["I"], Is1000 = 0, Is0100 = 0, Is1100 = 0, Is0010 = 0, Is1010 = 0,Is0110 = 0,Is1110 = 0,
            Is0001 = 0, Is1001 = 0, Is0101 = 0, Is1101 = 0, Is0011 = 0, Is1011 = 0,Is0111 = 0,Is1111 = 0,
            
            Ts0000t1 = 0, Ts1000t1 = 0, Ts0100t1 = 0, Ts1100t1 = 0, Ts0010t1 = 0, Ts1010t1 = 0,Ts0110t1 = 0,Ts1110t1 = 0,
            Ts0001t1 = 0, Ts1001t1 = 0, Ts0101t1 = 0, Ts1101t1 = 0, Ts0011t1 = 0, Ts1011t1 = 0,Ts0111t1 = 0,Ts1111t1 = 0,
            Ts0000t2 = 0, Ts1000t2 = 0, Ts0100t2 = 0, Ts1100t2 = 0, Ts0010t2 = 0, Ts1010t2 = 0,Ts0110t2 = 0,Ts1110t2 = 0,
            Ts0001t2 = 0, Ts1001t2 = 0, Ts0101t2 = 0, Ts1101t2 = 0, Ts0011t2 = 0, Ts1011t2 = 0,Ts0111t2 = 0,Ts1111t2 = 0,
            Ts0000t3 = 0, Ts1000t3 = 0, Ts0100t3 = 0, Ts1100t3 = 0, Ts0010t3 = 0, Ts1010t3 = 0,Ts0110t3 = 0,Ts1110t3 = 0,
            Ts0001t3 = 0, Ts1001t3 = 0, Ts0101t3 = 0, Ts1101t3 = 0, Ts0011t3 = 0, Ts1011t3 = 0,Ts0111t3 = 0,Ts1111t3 = 0,
            Ts0000t4 = 0, Ts1000t4 = 0, Ts0100t4 = 0, Ts1100t4 = 0, Ts0010t4 = 0, Ts1010t4 = 0,Ts0110t4 = 0,Ts1110t4 = 0,
            Ts0001t4 = 0, Ts1001t4 = 0, Ts0101t4 = 0, Ts1101t4 = 0, Ts0011t4 = 0, Ts1011t4 = 0,Ts0111t4 = 0,Ts1111t4 = 0,
            
            Rs0000 = 0, Rs1000 = 0, Rs0100 = 0, Rs1100 = 0, Rs0010 = 0, Rs1010 = 0,Rs0110 = 0,Rs1110 = 0,
            Rs0001 = 0, Rs1001 = 0, Rs0101 = 0, Rs1101 = 0, Rs0011 = 0, Rs1011 = 0,Rs0111 = 0,Rs1111 = 0,
            
            incidence = 0 , 
            M=0, 
            rR=0
  )
  
  names(init)[names(init) == "S.E"] <- "S"
  names(init)[names(init) == "Es0000.E"] <- "Es0000"
  names(init)[names(init) == "Ls0000.L"] <- "Ls0000"
  names(init)[names(init) == "Is0000.I"] <- "Is0000"
  
  #first period (introduction rifampicin)
  source("/home/louis/Bureau/stage/script/likelihood/scenarios/scenario_1971-2015.R")
  times <- seq(0, 44*12, by=1)
  simulation <- as.data.frame(ode(init,times,sir_equations,parms=pars))
  simulation <- simulation[-1,]
  
  #second period (introduction bedaqualine)
  init<-tail(simulation, n=1)[2:133]
  init2<-unlist(init)
  names(init2)<-names(init)
  init<-init2
  times <- seq(44*12, 51*12, by=1)
  source("/home/louis/Bureau/stage/script/likelihood/scenarios/scenario_2015-2021.R")
  simulation2 <-  as.data.frame(ode(init,times,sir_equations,parms=pars))
  simulation2 <- simulation2[-1,]
  
  #keep only the data that correspond to the data of WHO
  simulation<- rbind(simulation, simulation2)
  simulation <- tail(simulation, n = 22*12)
  
  if (nrow(simulation) == 22*12){
  #compare simulate data and real data
  Incidence<-as.vector(rowsum(simulation[,131], rep(1:22, each=12)))/as.vector(rowsum(rowSums(simulation[,2:130]), rep(1:22, each=12)))
  dead <-as.vector(rowsum(simulation[,132], rep(1:22, each=12)))/as.vector(rowsum(rowSums(simulation[,2:130]), rep(1:22, each=12)))
  rR<-as.vector(rowsum(simulation[,133], rep(1:22, each=12)))/as.vector(rowsum(simulation[,131], rep(1:22, each=12)))
  rR<-rR[16:22]
  cost <- sum((Incidence-data$incidence/100000)^2) + sum((dead-data$death/100000)^2) + sum((rR-data$Rresist[16:22]/100)^2)
  }else{ cost <- Inf}
  
  return(cost)
}

trans <- function(pars) {
  pars["infection"] <- exp(pars["infection"])
  pars["beta"] <- exp(pars["beta"])
  pars["s"] <- exp(pars["s"])
  pars["t"] <- exp(pars["t"])
  pars["I"] <- exp(pars["I"])
  pars["E"] <- exp(pars["E"])
  pars["L"] <- exp(pars["L"])
  return(pars)
}

fixed <- c( pi = 0.00155, #birth
            
            #infection
            infection = exp(-11.76509894), 
            s=exp(-1.50691782),
            
            delta1 =0.03333 , delta2 = 0.00017, #progression to active disease
            
            epsilon = 0.304 , #progression rate from early latent to late latent
            
            ####treated  rate
            beta=exp(-3.17473866),
            
            ###test rate
            sigma=1,#sigma = 0.121, 
            
            ### good result for a test rate 
            omega=1,#omega=0.63,
            
            ####cured with different treatment 
            t = exp(-0.03092109),
            
            omicron=0.01667, #self cured
            
            rho= 0.0036, #acquisition resistance
            
            theta= 0.01667, #death due to tuberculosis 
            #(1. FigurRs of the dead: a decade of tuberculosis morta1sity registrations in South Africa. https://journa1ss.co.za/doi/epdf/10.7196/SAMJ.2019.v109i10.14073 doi:10.7196/SAMJ.2019.v109i10.14073.)
            
            phi=0.00008, #re-activation from recovered to infected active stage
            
            l1=0.175,
            l2=0.175,
            
            
            mu= 0.00077, 
            
            I = 1000, 
            
            E = 1000, 
            
            L = 10000
) #naturals death


free <- c(log(fixed["infection"]), 
          log(fixed["s"]),
          log(fixed["t"]),
          log(fixed["beta"]),
          log(fixed["I"]),
          log(fixed["E"]),
          log(fixed["L"])
)

#remove to fixed
fixed<-fixed[ ! names(fixed) %in%  names(free) ]

start.time <- Sys.time()

# Fit the model to the data
fit <- optim(free,cost,gr=NULL,fixed,data,method="Nelder-Mead",hessian=TRUE)

end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken

fit 

saveRDS(fit, file="fit.RData")


start.time <- Sys.time()

hess <- hessian(fit, func=cost, free,fixed,init, data , method="Richardson")

end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken

hess

saveRDS(hess, file="hessian.RData")


# create table that summarizes the results:
cols <- c("Estimate","2.5%","97.5%")
rows <- names(fit$par)
est <- as.data.frame(matrix(NA,nrow=length(rows),ncol=length(cols)))
rownames(est) <- rows
colnames(est) <- cols

# hint: R0=beta0/gamma
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
