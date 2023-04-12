rm(list=ls())
#### Load libraries ####
library(deSolve)
library(tidyverse)
library(stringr)
library(dplyr)
library(ggplot2)

###Initial States 
fit<-readRDS("~/Bureau/stage/script/likelihood/results/gaussian/fit.RData")

parameters_values <- c( pi = 0.00155, #birth
            
            #infection
            infection = exp(fit$par["infection"]), 
            s=exp(fit$par["s"]),
            
            delta1 =0.03333 , delta2 = 0.00017, #progression to active disease
            
            epsilon = 0.304 , #progression rate from early latent to late latent
            
            ####treated  rate
            beta=exp(fit$par["beta"]),
            
            ###test rate
            sigma=1,#sigma = 0.121, 
            
            ### good result for a test rate 
            omega=1,#omega=0.63,
            
            ####cured with different treatment 
            t = exp(fit$par["t"]),
            
            omicron=0.01667, #self cured
            
            rho= 0.0036, #acquisition resistance
            
            theta= 0.01667, #death due to tuberculosis 
            #(1. FigurRs of the dead: a decade of tuberculosis morta1sity registrations in South Africa. https://journa1ss.co.za/doi/epdf/10.7196/SAMJ.2019.v109i10.14073 doi:10.7196/SAMJ.2019.v109i10.14073.)
            
            phi=0.00008, #re-activation from recovered to infected active stage
            
            l1=0.175,
            l2=0.175,
            
            
            mu= 0.00077, 
            
            I = 2000, 
            
            E = 2000, 
            
            L = 30000
) #naturals death
names(parameters_values)[names(parameters_values) == "infection.infection"] <- "infection"
names(parameters_values)[names(parameters_values) == "s.s"] <- "s"
names(parameters_values)[names(parameters_values) == "t.t"] <- "t"
names(parameters_values)[names(parameters_values) == "beta.beta"] <- "beta"


state <- c("S" = 100000 - exp(fit$par["E"]) - exp(fit$par["L"]) - exp(fit$par["I"]), 
          Es0000 = exp(fit$par["E"]), Es1000 = 0, Es0100 = 0, Es1100 = 0, Es0010 = 0, Es1010 = 0,Es0110 = 0,Es1110 = 0,
          Es0001 = 0, Es1001 = 0, Es0101 = 0, Es1101 = 0, Es0011 = 0, Es1011 = 0,Es0111 = 0,Es1111 = 0,
          
          
          Ls0000 =  exp(fit$par["L"]), Ls1000 = 0, Ls0100 = 0, Ls1100 = 0, Ls0010 = 0, Ls1010 = 0,Ls0110 = 0,Ls1110 = 0,
          Ls0001 = 0, Ls1001 = 0, Ls0101 = 0, Ls1101 = 0, Ls0011 = 0, Ls1011 = 0,Ls0111 = 0,Ls1111 = 0,
          
          
          Is0000 = exp(fit$par["I"]), Is1000 = 0, Is0100 = 0, Is1100 = 0, Is0010 = 0, Is1010 = 0,Is0110 = 0,Is1110 = 0,
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

names(state)[names(state) == "S.E"] <- "S"
names(state)[names(state) == "Es0000.E"] <- "Es0000"
names(state)[names(state) == "Ls0000.L"] <- "Ls0000"
names(state)[names(state) == "Is0000.I"] <- "Is0000"

#function
#source("/home/louis/Bureau/stage/script/model/deterministic/generalized_model.R")

#scenario
source("/home/louis/Bureau/stage/script/model/deterministic/scenarios/bedaqualine/intro_rifampicin_44_year.R")
####

#  numerically solving the SIR model 
# Run
i=0

assign(paste0("sir_values_", i), ode(
  y = state,
  times = time_values,
  func = sir_equations,
  parms = parameters_values 
)
)

sir_values_0<-data.frame(sir_values_0)

data0<-data.frame("time" = sir_values_0[,1], 
                  "S" = sir_values_0[,2], 
                  "E" = rowSums(sir_values_0[,3:18]), 
                  "L"= rowSums(sir_values_0[,19:34]),
                  "I"=rowSums(sir_values_0[,35:50]), 
                  "T"= rowSums(sir_values_0[,51:114]), 
                  "R" = rowSums(sir_values_0[,115:130]),
                  "Infected" = rowSums(sir_values_0[,35:50]),
                  "number_sensible" = rowSums(sir_values_0[,str_sub(names(sir_values_0), 1 , 1) == "I" & 
                                                                       (str_sub(names(sir_values_0), 3 , 3) == "0" | str_sub(names(sir_values_0), 4 , 4) == "0")]),
                  "number_1st_line_resistant" = rowSums(sir_values_0[,str_sub(names(sir_values_0), 1 , 1) == "I" & 
                                                                       (str_sub(names(sir_values_0), 3 , 3) == "1" | str_sub(names(sir_values_0), 4 , 4) == "1")]),
                  "number_2nd_line_resistant" = rowSums(sir_values_0[,str_sub(names(sir_values_0), 1 , 1) == "I" & str_sub(names(sir_values_0), 6 , 6) == "1"]),
                  "freq_1st_line_resistant" = (rowSums(sir_values_0[,str_sub(names(sir_values_0), 1 , 1) == "I" & 
                                                                      (str_sub(names(sir_values_0), 3 , 3) == "1" | str_sub(names(sir_values_0), 4 , 4) == "1")]) /
                                                 rowSums(sir_values_0[,35:50])),
                  "freq_2nd_line_resistant" = (rowSums(sir_values_0[,str_sub(names(sir_values_0), 1 , 1) == "I" & str_sub(names(sir_values_0), 6 , 6) == "1"]) / 
                                                 rowSums(sir_values_0[,35:50])), 
                  "Incidence" = sir_values_0[,131],
                  "scenario" = rep("intro rifampcin", each = length(sir_values_0[,1])))


state<-tail(sir_values_0[2:133], n=1)
state2<-unlist(state)
names(state2)<-names(state)
state<-state2


#scenario
source("/home/louis/Bureau/stage/script/model/deterministic/scenarios/bedaqualine/scenario_1st_line.R")
####

#  numerically solving the SIR model 
# Run
i=1

assign(paste0("sir_values_", i), ode(
  y = state,
  times = time_values,
  func = sir_equations,
  parms = parameters_values 
)
)
sir_values_1<-data.frame(sir_values_1)

data1<-data.frame("time" = sir_values_1[,1], 
                  "S" = sir_values_1[,2], 
                  "E" = rowSums(sir_values_1[,3:18]), 
                  "L"= rowSums(sir_values_1[,19:34]),
                  "I"=rowSums(sir_values_1[,35:50]), 
                  "T"= rowSums(sir_values_1[,51:114]), 
                  "R" = rowSums(sir_values_1[,115:130]),
                 "Infected" = rowSums(sir_values_1[,35:50]),
                 "number_sensible" = rowSums(sir_values_1[,str_sub(names(sir_values_1), 1 , 1) == "I" & 
                                                                      (str_sub(names(sir_values_1), 3 , 3) == "0" | str_sub(names(sir_values_1), 4 , 4) == "0")]),
                 "number_1st_line_resistant" = rowSums(sir_values_1[,str_sub(names(sir_values_1), 1 , 1) == "I" & 
                                                     (str_sub(names(sir_values_1), 3 , 3) == "1" | str_sub(names(sir_values_1), 4 , 4) == "1")]),
                 "number_2nd_line_resistant" = rowSums(sir_values_1[,str_sub(names(sir_values_1), 1 , 1) == "I" & str_sub(names(sir_values_1), 6 , 6) == "1"]),
                 "freq_1st_line_resistant" = (rowSums(sir_values_1[,str_sub(names(sir_values_1), 1 , 1) == "I" & 
                                                                  (str_sub(names(sir_values_1), 3 , 3) == "1" | str_sub(names(sir_values_1), 4 , 4) == "1")]) /
                                                rowSums(sir_values_1[,35:50])),
                 "freq_2nd_line_resistant" = (rowSums(sir_values_1[,str_sub(names(sir_values_1), 1 , 1) == "I" & str_sub(names(sir_values_1), 6 , 6) == "1"]) / 
                                                rowSums(sir_values_1[,35:50])), 
                 "Incidence" = sir_values_1[,131],
                 "scenario" = rep("standard strategy", each = length(sir_values_1[,1])))

#scenario
source("/home/louis/Bureau/stage/script/model/deterministic/scenarios/bedaqualine/scenario_2nd_line.R")
####

#  numerically solving the SIR model 
# Run
i=2

assign(paste0("sir_values_", i), ode(
  y = state,
  times = time_values,
  func = sir_equations,
  parms = parameters_values 
)
)

sir_values_2<-data.frame(sir_values_2)
data2<-data.frame("time" = sir_values_2[,1], 
                  "S" = sir_values_2[,2], 
                  "E" = rowSums(sir_values_2[,3:18]), 
                  "L"= rowSums(sir_values_2[,19:34]),
                  "I"=rowSums(sir_values_2[,35:50]), 
                  "T"= rowSums(sir_values_2[,51:114]), 
                  "R" = rowSums(sir_values_2[,115:130]),
                  "Infected" = rowSums(sir_values_2[,35:50]),
                  "number_sensible" = rowSums(sir_values_2[,str_sub(names(sir_values_2), 1 , 1) == "I" & 
                                                             (str_sub(names(sir_values_2), 3 , 3) == "0" | str_sub(names(sir_values_2), 4 , 4) == "0")]),
                  "number_1st_line_resistant" = rowSums(sir_values_2[,str_sub(names(sir_values_2), 1 , 1) == "I" & 
                                                                       (str_sub(names(sir_values_2), 3 , 3) == "1" | str_sub(names(sir_values_2), 4 , 4) == "1")]),
                  "number_2nd_line_resistant" = rowSums(sir_values_2[,str_sub(names(sir_values_2), 1 , 1) == "I" & str_sub(names(sir_values_2), 6 , 6) == "1"]),
                  "freq_1st_line_resistant" = (rowSums(sir_values_2[,str_sub(names(sir_values_2), 1 , 1) == "I" & 
                                                                      (str_sub(names(sir_values_2), 3 , 3) == "1" | str_sub(names(sir_values_2), 4 , 4) == "1")]) /
                                                 rowSums(sir_values_2[,35:50])),
                  "freq_2nd_line_resistant" = (rowSums(sir_values_2[,str_sub(names(sir_values_2), 1 , 1) == "I" & str_sub(names(sir_values_2), 6 , 6) == "1"]) / 
                                                 rowSums(sir_values_2[,35:50])), 
                  "Incidence" = sir_values_2[,131],
                  "scenario" = rep("TRUNCATE strategy", each = length(sir_values_2[,1])))


data<-rbind(data0,data1, data2)

ggplot(data, aes(x= time, y=Infected, group=scenario)) +
    geom_line(aes(color=scenario)) + 
    xlab("time (years)")

ggplot(data, aes(x= time, y=number_sensible, group=scenario)) +
  geom_line(aes(color=scenario)) + 
  xlab("time (years)")

ggplot(data, aes(x=time, y=number_1st_line_resistant, group=scenario)) +
  geom_line(aes(color=scenario)) + 
  xlab("time (years)")

ggplot(data, aes(x=time, y=number_2nd_line_resistant, group=scenario)) +
  geom_line(aes(color=scenario)) + 
  xlab("time (years)")


ggplot(data, aes(x=time, y=freq_1st_line_resistant, group=scenario)) +
  geom_line(aes(color=scenario)) + 
  xlab("time (years)")

ggplot(data, aes(x=time, y=freq_2nd_line_resistant, group=scenario)) +
  geom_line(aes(color=scenario)) + 
  xlab("time (years)")



compare1<-rbind(data0, data1)
compare1[,1]<-compare1[,1]+1971
#####
c<-data.frame(time = compare1[,1], 
                 S = compare1[,2], 
                 E = compare1[,3], 
                 L= compare1[,4],
                 I=compare1[,5], 
                 T= compare1[,6], 
                 R = compare1[,7])
c<-gather(c , "Compartment", "Number_of_Individuals", -time)

ggplot(c, aes(x=time, y=Number_of_Individuals, group=Compartment)) +
  geom_line(aes(color=Compartment)) + 
  xlab("time (years)")

#LOAD THE REAL DATA
#data_ref1<-read.csv2("/home/louis/Bureau/stage/script/likelihood/data/TB_burden_countries_2023-02-22.csv", sep=",")
#data_ref1<-subset(data_ref1, data_ref1[,1]== "South Africa")
#data_ref2<-read.csv2("/home/louis/Bureau/stage/script/likelihood/data/MDR_RR_TB_burden_estimates_2023-02-22.csv", sep=",")
#data_ref2<-subset(data_ref2, data_ref2[,1]== "South Africa")

#data_ref<-data_ref1[,c(6,7,8, 35)]
#colnames(data)<-c("times", "pop", "e_inc_100k", "e_mort_100k")
#data_ref[,2]<-as.integer(data_ref[,2])
#data_ref[,3]<-as.integer(data_ref[,3])
#data_ref[,4]<-as.integer(data_ref[,4])
#pourcent_rif_resis<-data_ref2[,c(6,8)]
#colnames(pourcent_rif_resis)<-c("times", "pourcent_resist")
#pourcent_rif_resis[,2]<-as.numeric(pourcent_rif_resis[,2])

#compare1<-compare1[compare1[,1] > 2000 & compare1[,1] < 2022,]
#Incidence<-as.vector(rowsum(compare1[,"Incidence"], rep(1:22, each=12)))/as.vector(rowsum(rowSums(compare1[,2:7]), rep(1:22, each=12)))
#data_ref$e_inc_100k
#compare1$Incidence


#incidence<-data.frame(time = compare1[,1], 
#                      Incidence = compare1[,"Incidence"], 
#                      type = "simulation")
#incidence_ref<-data.frame(time = data_ref[,1], 
#                          Incidence = data_ref[,3], 
#                          type = "real data")#


#incidence<-rbind(incidence, incidence_ref)

#ggplot(incidence, aes(x= time, y=Incidence, group=type)) +
#  geom_line(aes(color=type)) + 
#  xlab("time (years)")

#Infected<-data.frame(time = sir_values[,1]/12, sir_values[,35:50], all = rowSums(sir_values[,35:50]))
#Infected<-gather(Infected , "Strain", "Number_of_Individuals", -time)

#ggplot(Infected, aes(x=time, y=Number_of_Individuals, group=Strain)) +
#  geom_line(aes(color=Strain)) + 
#  xlab("time (years)")


#resistance <- data.frame(sir_values)
#Resistance<-data.frame(time = resistance[,1]/12, 
#                       "sensible" = resistance[,str_sub(names(resistance), 1 , 1) == "I" & str_count(names(resistance), "1") == 0], 
#                       "t1" = rowSums(resistance[,str_sub(names(resistance), 1 , 1) == "I" & 
#                                                  (str_sub(names(resistance), 3 , 3) == "1" | str_sub(names(resistance), 4 , 4) == "1")]),
#                       "t2" = rowSums(resistance[,str_sub(names(resistance), 1 , 1) == "I" & str_sub(names(resistance), 5 , 5) == "1"]),
#                       "t3" = rowSums(resistance[,str_sub(names(resistance), 1 , 1) == "I" & str_sub(names(resistance), 6 , 6) == "1"])
#                       )
                                  
#Resistance<-gather(Resistance , "Strain", "Number_of_Individuals", -time)
#ggplot(Resistance, aes(x=time, y=Number_of_Individuals, group=Strain)) +
#  geom_line(aes(color=Strain)) + 
#  xlab("time (years)")

#write.table(data.frame(tail(sir_values_1[,2:132], n=1)), "/home/louis/Bureau/stage/script/model/deterministic/initial_conditions.txt", row.names=FALSE, col.names = TRUE)
#for (i in 1:length(parameters_values)){ assign(names(parameters_values[i]),parameters_values[i])}


##########################################################


