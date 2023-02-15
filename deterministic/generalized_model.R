#### Load libraries ####
library(deSolve)
library(tidyverse)
library(stringr)
library(dplyr)

####Initial States ####
state <-c(S =10000, 
          
          Es0000 = 4000, Es1000 = 0, Es0100 = 0, Es1100 = 0, Es0010 = 0, Es1010 = 0,Es0110 = 0,Es1110 = 0,
          Es0001 = 0, Es1001 = 0, Es0101 = 0, Es1101 = 0, Es0011 = 0, Es1011 = 0,Es0111 = 0,Es1111 = 0,
          
          
          Ls0000 = 1200, Ls1000 = 0, Ls0100 = 0, Ls1100 = 0, Ls0010 = 0, Ls1010 = 0,Ls0110 = 0,Ls1110 = 0,
          Ls0001 = 0, Ls1001 = 0, Ls0101 = 0, Ls1101 = 0, Ls0011 = 0, Ls1011 = 0,Ls0111 = 0,Ls1111 = 0,
          
          
          Is0000 = 1, Is1000 = 0, Is0100 = 0, Is1100 = 0, Is0010 = 0, Is1010 = 0,Is0110 = 0,Is1110 = 0,
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
          Rs0001 = 0, Rs1001 = 0, Rs0101 = 0, Rs1101 = 0, Rs0011 = 0, Rs1011 = 0,Rs0111 = 0,Rs1111 = 0
          
)

#### Parameters ####
parameters_values <- c( pi = 0.01856/12, #birth
                        
                        #infection 
                        alphas0000 = (615/100000)/12, alphas1000 = (615/100000)/12, alphas0100 = (615/100000)/12, alphas0010 = (615/100000)/12, 
                        alphas1100 = (615/100000)/12, alphas1010 = (615/100000)/12, alphas0110 = (615/100000)/12, alphas1110 = (615/100000)/12,
                        alphas0001 = (615/100000)/12, alphas1001 = (615/100000)/12, alphas0101 = (615/100000)/12, alphas0011 = (615/100000)/12, 
                        alphas1101 = (615/100000)/12, alphas1011 = (615/100000)/12, alphas0111 = (615/100000)/12, alphas1111 = (615/100000)/12,
                     
                        ####reinfection 
                        #L -> E 
                        Lambda1s0000=0.21/12,Lambda1s1000=0.21/12,Lambda1s0100=0.21/12,Lambda1s1100=0.21/12,
                        Lambda1s0010=0.21/12,Lambda1s1010=0.21/12,Lambda1s0110=0.21/12,Lambda1s1110=0.21/12,
                        Lambda1s0001=0.21/12,Lambda1s1001=0.21/12,Lambda1s0101=0.21/12,Lambda1s1101=0.21/12,
                        Lambda1s0011=0.21/12,Lambda1s1011=0.21/12,Lambda1s0111=0.21/12,Lambda1s1111=0.21/12,
                        
                        #R -> E
                        Lambda2s0000=0.21/12,Lambda2s1000=0.21/12,Lambda2s0100=0.21/12,Lambda2s1100=0.21/12,
                        Lambda2s0010=0.21/12,Lambda2s1010=0.21/12,Lambda2s0110=0.21/12,Lambda2s1110=0.21/12,
                        Lambda2s0001=0.21/12,Lambda2s1001=0.21/12,Lambda2s0101=0.21/12,Lambda2s1101=0.21/12,
                        Lambda2s0011=0.21/12,Lambda2s1011=0.21/12,Lambda2s0111=0.21/12,Lambda2s1111=0.21/12,
                      
                       
                        epsilon = 8/12, #progression rate from early latent to late latent
                       
                       
                        delta1 = 0.4/12, delta2 = 0.0020/12, #progression to active disease
                      
                        ####treated  rate
                        beta=0.0,
                        
                        ###test rate
                        sigma = 0.9, 
                        
                        ### good result for a test rate 
                        omega=0.9,
                
                        ####cured with different treatment 
                        taus0000t1=0.008,taus1000t1=0,taus0100t1=0,taus1100t1=0,taus0010t1=0.008,taus1010t1=0,taus0110t1=0,taus1110t1=0,
                        taus0001t1=0.008,taus1001t1=0,taus0101t1=0,taus1101t1=0,taus0011t1=0.008,taus1011t1=0,taus0111t1=0,taus1111t1=0,
                        
                        taus0000t2=0.004,taus1000t2=0.004,taus0100t2=0.004,taus1100t2=0.004,taus0010t2=0,taus1010t2=0,taus0110t2=0,taus1110t2=0,
                        taus0001t2=0.004,taus1001t2=0.004,taus0101t2=0.004,taus1101t2=0.004,taus0011t2=0,taus1011t2=0,taus0111t2=0,taus1111t2=0,
                        
                        taus0000t3=0.004,taus1000t3=0.004,taus0100t3=0.004,taus1100t3=0.004,taus0010t3=0.004,taus1010t3=0.004,taus0110t3=0.004,taus1110t3=0.004,
                        taus0001t3=0,taus1001t3=0,taus0101t3=0,taus1101t3=0,taus0011t3=0,taus1011t3=0,taus0111t3=0,taus1111t3=0,
                        
                        taus0000t4=0.002,taus1000t4=0.002,taus0100t4=0.002,taus1100t4=0.002,taus0010t4=0.002,taus1010t4=0.002,taus0110t4=0.002,taus1110t4=0.002,
                        taus0001t4=0.002,taus1001t4=0.002,taus0101t4=0.002,taus1101t4=0.002,taus0011t4=0.002,taus1011t4=0.002,taus0111t4=0.002,taus1111t4=0.002,
                      
                       
                        omicron=0.2/12, #self cured
                       
                        rho= 0.01/12, #acquisition resistance
                       
                        ita=0.02/12, #ita_old= 0, ita_new=0,  #treatment shift due to rRsistance
                       
                        theta= 0.026/12, #death due to tuberculosis 
                        #(1. FigurRs of the dead: a decade of tuberculosis morta1sity registrations in South Africa. https://journa1ss.co.za/doi/epdf/10.7196/SAMJ.2019.v109i10.14073 doi:10.7196/SAMJ.2019.v109i10.14073.)
                       
                        phi=0.01/12, #re-activation from recovered to infected active stage
                      
                    
                        mu= 0.00926/12 ) #natura1s death

####
#### MATRICES #####
##row and col name 
strains<-paste0(c("0","1"), rep(c("0","1"), each=2),rep(c("0","1"), each=4), rep(c("0","1"), each=8))
treatments<-c("t1","t2","t3", "t4")

### Resistance matrix ###
#resistance acquisition with 1st treatment
resistance_matrix_t1 <- matrix(c(0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 
                                 0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), ncol=16, byrow=T)

colnames(resistance_matrix_t1)<-strains
rownames(resistance_matrix_t1)<-strains

#resistance acquisition with 2nd treatment
resistance_matrix_t2 <- matrix(c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0, 
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1, 
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), ncol=16, byrow=T)

colnames(resistance_matrix_t2)<-strains
rownames(resistance_matrix_t2)<-strains

#resistance acquisition with third treatment
resistance_matrix_t3 <- matrix(c(0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0, 
                                 0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), ncol=16, byrow=T)

colnames(resistance_matrix_t3)<-strains
rownames(resistance_matrix_t3)<-strains

#resistance acquisition with fourth treatment
resistance_matrix_t4 <- matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), ncol=16, byrow=T)

colnames(resistance_matrix_t4)<-strains
rownames(resistance_matrix_t4)<-strains

#list with all resistance matrices
resistance_matrix<-list(resistance_matrix_t1, resistance_matrix_t2, resistance_matrix_t3, resistance_matrix_t4)


### New treatment if not adapted against the strain ###
#Here, we assume the shift it's always good

#the adapted treatment is the first one
new_treatment_t1 <- matrix(c(0,1,1,1,
                             0,0,0,0,
                             0,0,0,0,
                             0,0,0,0,
                             0,1,1,1,
                             0,0,0,0,
                             0,0,0,0,
                             0,0,0,0,
                             0,1,1,1,
                             0,0,0,0,
                             0,0,0,0,
                             0,0,0,0,
                             0,1,1,1,
                             0,0,0,0,
                             0,0,0,0,
                             0,0,0,0), ncol=4, byrow=T)
colnames(new_treatment_t1)<-treatments
rownames(new_treatment_t1)<-strains

#the adapted treatment is the second one
new_treatment_t2 <- matrix(c(0,0,0,0,
                             1,0,1,1,
                             1,0,1,1,
                             1,0,1,1,
                             0,0,0,0,
                             0,0,0,0,
                             0,0,0,0,
                             0,0,0,0,
                             0,0,0,0,
                             1,0,1,1,
                             1,0,1,1,
                             1,0,1,1,
                             0,0,0,0,
                             0,0,0,0,
                             0,0,0,0,
                             0,0,0,0), ncol=4, byrow=T)
colnames(new_treatment_t2)<-treatments
rownames(new_treatment_t2)<-strains

##the adapted treatment is the third one
new_treatment_t3 <- matrix(c(0,0,0,0,
                             0,0,0,0,
                             0,0,0,0,
                             0,0,0,0,
                             0,0,0,0,
                             0,0,0,0,
                             0,0,0,0,
                             0,0,0,0,
                             0,0,0,0,
                             0,0,0,0,
                             0,0,0,0,
                             0,0,0,0,
                             0,0,0,0,
                             0,0,0,0,
                             0,0,0,0,
                             0,0,0,0), ncol=4, byrow=T)

colnames(new_treatment_t3)<-treatments
rownames(new_treatment_t3)<-strains

#the adapted treatment is the fourth one
new_treatment_t4 <- matrix(c(0,0,0,0,
                             0,0,0,0,
                             0,0,0,0,
                             0,0,0,0,
                             0,0,0,0,
                             1,1,1,0,
                             1,1,1,0,
                             1,1,1,0,
                             0,0,0,0,
                             0,0,0,0,
                             0,0,0,0,
                             0,0,0,0,
                             0,0,0,0,
                             1,1,1,0,
                             1,1,1,0,
                             1,1,1,0), ncol=4, byrow=T)

colnames(new_treatment_t4)<-treatments
rownames(new_treatment_t4)<-strains

#list with all new treattment matrices
new_treatment<-list(new_treatment_t1, new_treatment_t2, new_treatment_t3, new_treatment_t4)


#### Treatment choice ###
#If there is no test
no_test<-matrix(c(1,0,0,0,
                  1,0,0,0,
                  1,0,0,0,
                  1,0,0,0,
                  1,0,0,0,
                  1,0,0,0,
                  1,0,0,0,
                  1,0,0,0,
                  1,0,0,0,
                  1,0,0,0,
                  1,0,0,0,
                  1,0,0,0,
                  1,0,0,0,
                  1,0,0,0,
                  1,0,0,0,
                  1,0,0,0), ncol=4, byrow=T)

colnames(no_test)<-treatments
rownames(no_test)<-strains

#if there is a test and the reslt are good
good_test_result<-matrix(c(1,0,0,0,
                                    0,1,0,0,
                                    0,1,0,0,
                                    0,1,0,0,
                                    1,0,0,0,
                                    0,0,0,1,
                                    0,0,0,1,
                                    0,0,0,1,
                                    1,0,0,0,
                                    0,1,0,0,
                                    0,1,0,0,
                                    0,1,0,0,
                                    1,0,0,0,
                                    0,0,0,1,
                                    0,0,0,1,
                                    0,0,0,1), ncol=4, byrow=T)
colnames(good_test_result)<-treatments
rownames(good_test_result)<-strains


#if there is a test and the reslt are false
false_result<-matrix(c(1,1,1,1,
                       1,1,1,1,
                       1,1,1,1,
                       1,1,1,1,
                       1,1,1,1,
                       1,1,1,1,
                       1,1,1,1,
                       1,1,1,1,
                       1,1,1,1,
                       1,1,1,1,
                       1,1,1,1,
                       1,1,1,1,
                       1,1,1,1,
                       1,1,1,1,
                       1,1,1,1,
                       1,1,1,1), ncol=4, byrow=T)

colnames(false_result)<-treatments
rownames(false_result)<-strains
nb_treatment<-1/length(treatments)

####place of specific character####
#strains
s_first<-3
s_last<-nchar(strains[1]) + 2
#treatments
t_last<-nchar(names(subset(state, str_sub(names(state), 1, 1) == "T")[1]))
t_first<-nchar(names(subset(state, str_sub(names(state), 1, 1) == "T")[1])) - 1


#### Vector with parameters ####
alpha<-subset(parameters_values, str_sub(names(parameters_values), 1, 5) == "alpha")
lambda_1<- subset(parameters_values, str_sub(names(parameters_values), 1, 7) == "Lambda1")
lambda_2<-subset(parameters_values, str_sub(names(parameters_values), 1, 7) == "Lambda2")

tau<-subset(parameters_values, str_sub(names(parameters_values), 1, 3) == "tau")


#### write differentials equations ####

sir_equations <- function(time, state, parameters_values) {
  with(as.list(c(state, parameters_values)), {
    
    ### Groups ###
    S <- sum(subset(state, str_sub(names(state), 1, 1) == "S"))
    E <- sum(subset(state, str_sub(names(state), 1, 1) == "E"))
    L <- sum(subset(state, str_sub(names(state), 1, 1) == "L"))
    I <- sum(subset(state, str_sub(names(state), 1, 1) == "I"))
    T <- sum(subset(state, str_sub(names(state), 1, 1) == "T"))
    R <-sum(subset(state, str_sub(names(state), 1, 1) == "R"))
    N=S+E+L+I+T+R
    
    
    
    #vector
    vE  <-subset(state, str_sub(names(state), 1, 1) == "E") #vector with all the Early Latent types
    vL  <-subset(state, str_sub(names(state), 1, 1) == "L") #vector with all the Late latent types
    vI  <-subset(state, str_sub(names(state), 1, 1) == "I") #vector with all the Infected types
    vTs <-tapply(subset(state, str_sub(names(state), 1, 1) == "T"), str_sub(names(subset(state, str_sub(names(state), 1, 1) == "T")), s_first, s_last), sum) 
    #vector with all strains of Treated (don't take account about treatment here)
    vT  <-subset(state, str_sub(names(state), 1, 1) == "T") #vector with all the Recovered types
    vR  <-subset(state, str_sub(names(state), 1, 1) == "R") #vector with all the Recovered types
    
    
    
    
    
    
    ### Equation ###
    
    #ODE for Susceptible compartment
    assign(paste0("dS"),
            pi  * N 
           - mu * S 
           - sum((vTs + vI) * alpha * S) 
    )
           #Check if after the next step, individual in the compartment is < 0, in the other case, individual in the compartment = 0
  #         if (S+ eval(parse(text=paste0("d", names(subset(state, str_sub(names(state), 1, 1) == "S"))))) < 0){
  #           state[names(subset(state, str_sub(names(state), 1, 1) == "S"))] <- 0
  #         }
           
             
    #ODE for Early Latent stage compartments
    for (i in 1: length(vE)){
      assign(paste0("d", names(vE)[i]), 
           (vTs[i] + vI[i]) * alpha[i] * S +
             (vTs[i] + vI[i]) * lambda_1[i] * L +
             (vTs[i] + vI[i]) * lambda_2[i] * R - 
             vE[i] * (epsilon + mu + delta1) 
              )
      #Check if after the next step, individual in the compartment is < 0, in the other case, individual in the compartment = 0
   #   if (vE[i] + eval(parse(text=paste0("d", names(subset(state, str_sub(names(state), 1, 1) == "E")[i])))) < 0){
  #      state[names(subset(state, str_sub(names(state), 1, 1) == "E")[i])] <- 0
  #    }
     
   }
    #ODE for Late Latent stage compartments
    for (i in 1:length(vL)){
      assign(paste0("d", names(vL[i])),
              epsilon * vE[i] -
               (mu + delta2)  * vL[i] -
              sum((vTs + vI) * lambda_1 * vL[i]) 

         )
      #Check if after the next step, individual in the compartment is < 0, in the other case, individual in the compartment = 0
   #   if (vL[i] + eval(parse(text=paste0("d", names(subset(state, str_sub(names(state), 1, 1) == "L")[i])))) < 0){
  #      state[names(subset(state, str_sub(names(state), 1, 1) == "L")[i])] <- 0
  #    }
    }
    
    #ODE for Infected stage compartments
    for (i in 1:length(vI)){
      assign(paste0("d", names(vI[i])), 
            delta1 * vE[i] + 
            delta2 * vL[i] + 
            phi * vR[i] -
            (mu + theta + omicron + beta) * vI[i]
  )
      #Check if after the next step, individual in the compartment is < 0, in the other case, individual in the compartment = 0
   #   if (vI[i] + eval(parse(text=paste0("d", names(subset(state, str_sub(names(state), 1, 1) == "I")[i])))) < 0){
  #      state[names(subset(state, str_sub(names(state), 1, 1) == "I")[i])] <- 0
  #    }
 
        }
   
     #ODE for Treated stage compartments
    for (i in 1:length(vT)){
      assign(paste0("d", names(vT)[i]), 
             - ( mu + theta + omicron + tau[i]) * vT[i] -
            
              sum( rho * (resistance_matrix[[as.numeric(str_sub(names(vT[i]),-1))]]
                            [str_sub(names(vT), s_first, s_last)[i],]) * vT[i]) +
             
               sum( rho * (resistance_matrix[[as.numeric(str_sub(names(vT[i]),-1))]]
                    [,str_sub(names(vT), s_first, s_last)[i]]) * 
                      subset(vT, str_sub(names(vT), -1) == str_sub(names(vT[i]), -1))) -
             
              sum( ita *  (new_treatment_t1[str_sub(names(vT), s_first, s_last)[i], str_sub(names(vT), t_first, t_last)[i]]
                    +       new_treatment_t2[str_sub(names(vT), s_first, s_last)[i], str_sub(names(vT), t_first, t_last)[i]]
                    +       new_treatment_t3[str_sub(names(vT), s_first, s_last)[i], str_sub(names(vT), t_first, t_last)[i]]
                    +       new_treatment_t4[str_sub(names(vT), s_first, s_last)[i], str_sub(names(vT), t_first, t_last)[i]]) *
               vT[i]) +
  
              sum(ita * ((eval(parse(text=paste0("new_treatment_",str_sub(names(vT),t_first, t_last)))[i])
                    [str_sub(names(vT[i]), s_first, s_last),]) 
                   *  subset(vT, str_sub(names(vT), s_first, s_last) == str_sub(names(vT[i]), s_first, s_last)))) +
                   
              beta * subset(vI, str_sub(names(vI), s_first, s_last) ==  str_sub(names(vT[i]), s_first, s_last)) * 
               ((1-sigma) * no_test[str_sub(names(vT)[i], s_first, s_last), 
                                    str_sub(names(vT)[i], t_first, t_last)] +
                 sigma * omega * good_test_result[str_sub(names(vT)[i], s_first, s_last), 
                                                   str_sub(names(vT)[i], t_first, t_last)] +
                 nb_treatment * sigma * (1-omega) * false_result[str_sub(names(vT)[i], s_first, s_last), 
                                                        str_sub(names(vT)[i], t_first, t_last)])
      )
      
      #Check if after the next step, individual in the compartment is < 0, in the other case, individual in the compartment = 0
   #   if (vT[i] + eval(parse(text=paste0("d", names(subset(state, str_sub(names(state), 1, 1) == "T")[i])))) < 0){
  #      state[names(subset(state, str_sub(names(state), 1, 1) == "T")[i])] <- 0
  #    }
    
    }
    
    #ODE for Recovered stage compartments
    for (i in 1:length(vR)){
      assign(paste0("d", names(vR)[i]), 
            - ( phi + mu)  * vR[i] -
               sum((vTs + vI)* lambda_2 * 
                       vR[i]) +
              omicron * vI[i] +
             
             sum(subset(tau, str_sub(names(tau), 5, 8) == str_sub(names(vR)[i], s_first, s_last))   * 
                     subset(vT, str_sub(names(vT), s_first, s_last) == str_sub(names(vR)[i], s_first, s_last))) +
            omicron *  vTs[i]
             )
      #Check if after the next step, individual in the compartment is < 0, in the other case, individual in the compartment = 0
   #   if (vR[i] + eval(parse(text=paste0("d", names(subset(state, str_sub(names(state), 1, 1) == "R")[i])))) < 0){
  #      state[names(subset(state, str_sub(names(state), 1, 1) == "R")[i])] <- 0
  #    }

    }
    
    return(list(c(dS, 
                  
                  dEs0000, dEs1000, dEs0100, dEs1100, dEs0010, dEs1010, dEs0110, dEs1110,
                  dEs0001, dEs1001, dEs0101, dEs1101, dEs0011, dEs1011, dEs0111, dEs1111,
                  
                  dLs0000, dLs1000, dLs0100, dLs1100, dLs0010, dLs1010, dLs0110, dLs1110,
                  dLs0001, dLs1001, dLs0101, dLs1101, dLs0011, dLs1011, dLs0111, dLs1111,
                  
                  dIs0000, dIs1000, dIs0100, dIs1100, dIs0010, dIs1010, dIs0110, dIs1110,
                  dIs0001, dIs1001, dIs0101, dIs1101, dIs0011, dIs1011, dIs0111, dIs1111,
                  
                  dTs0000t1, dTs1000t1, dTs0100t1, dTs1100t1, dTs0010t1, dTs1010t1, dTs0110t1, dTs1110t1,
                  dTs0001t1, dTs1001t1, dTs0101t1, dTs1101t1, dTs0011t1, dTs1011t1, dTs0111t1, dTs1111t1,
                  dTs0000t2, dTs1000t2, dTs0100t2, dTs1100t2, dTs0010t2, dTs1010t2, dTs0110t2, dTs1110t2,
                  dTs0001t2, dTs1001t2, dTs0101t2, dTs1101t2, dTs0011t2, dTs1011t2, dTs0111t2, dTs1111t2,
                  dTs0000t3, dTs1000t3, dTs0100t3, dTs1100t3, dTs0010t3, dTs1010t3, dTs0110t3, dTs1110t3,
                  dTs0001t3, dTs1001t3, dTs0101t3, dTs1101t3, dTs0011t3, dTs1011t3, dTs0111t3, dTs1111t3,
                  dTs0000t4, dTs1000t4, dTs0100t4, dTs1100t4, dTs0010t4, dTs1010t4, dTs0110t4, dTs1110t4,
                  dTs0001t4, dTs1001t4, dTs0101t4, dTs1101t4, dTs0011t4, dTs1011t4, dTs0111t4, dTs1111t4,
                  
                  dRs0000, dRs1000, dRs0100, dRs1100, dRs0010, dRs1010, dRs0110, dRs1110,
                  dRs0001, dRs1001, dRs0101, dRs1101, dRs0011, dRs1011, dRs0111, dRs1111
                  
                  )))
  })
}

####  the points in time where to calculate variables values ####
time_values <- seq(0, 5*12, by = 1) 
####

####  numerically solving the SIR model ####

sir_values_1 <- ode(
  y = state,
  times = time_values,
  func = sir_equations,
  parms = parameters_values 
)
  
#####

time<-sir_values_1[,1]
S<-sir_values_1[,2]
E<-rowSums(sir_values_1[,3:18])
L <-rowSums(sir_values_1[,19:34])
I <-rowSums(sir_values_1[,35:50])
T <-rowSums(sir_values_1[,51:114])
R <-rowSums(sir_values_1[,115:130])

a<-plot(time, rowSums(sir_values_1[,2:130]))

# plotting the time seriRs of susceptiblRs:
plot(time, S, type = "l", col = "blue",
     xlab = "time (months)", ylab = "number of people", ylim = c(0,13000))
# adding the time seriRs of early latent:
lines(time, E, col = "darkorchid")
# adding the time seriRs of late latent:
lines(time, L, col = "darkmagenta")
# adding the time seriRs of infectious:
lines(time, I, col = "red")
# adding the time seriRs of treated:
lines(time, T, col = "darkgoldenrod2")
# adding the time seriRs of recovered:
lines(time, R, col = "green")

# adding a legend:
legend("topright", c("Susceptibles","Early Latent", "Late latent", "Infectious", "Treated","Recovered"),
       col = c("blue","darkorchid","darkmagenta", "red", "darkgoldenrod2" ,"green"), lty = 1, bty = "n")

S


#plot(time, sir_values_1[,35], type = "l", col = "blue",
#     xlab = "time (days)", ylab = "number of people", ylim = c(0,1000))
# adding the time seriRs of early latent:
#lines(time, sir_values_1[,36], col = "darkorchid")
# adding the time seriRs of late latent:
#lines(time, sir_values_1[,37], col = "darkmagenta")
# adding the time seriRs of late latent:
#lines(time, sir_values_1[,38], col = "red")



#Verification
#####
for ( i in 1:length(parameters_values)){
  
  assign(paste0(names(parameters_values[i])),parameters_values[i] )
}

state<-sir_values_1[2,2:130]
