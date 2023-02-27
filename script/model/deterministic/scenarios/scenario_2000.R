##### 70 first years #####
###Initial States 

state <-c(S = 97500, 
          
          Es0000 = 700, Es1000 = 0, Es0100 = 0, Es1100 = 0, Es0010 = 0, Es1010 = 0,Es0110 = 0,Es1110 = 0,
          Es0001 = 0, Es1001 = 0, Es0101 = 0, Es1101 = 0, Es0011 = 0, Es1011 = 0,Es0111 = 0,Es1111 = 0,
          
          
          Ls0000 = 1300, Ls1000 = 0, Ls0100 = 0, Ls1100 = 0, Ls0010 = 0, Ls1010 = 0,Ls0110 = 0,Ls1110 = 0,
          Ls0001 = 0, Ls1001 = 0, Ls0101 = 0, Ls1101 = 0, Ls0011 = 0, Ls1011 = 0,Ls0111 = 0,Ls1111 = 0,
          
          
          Is0000 = 500, Is1000 = 0, Is0100 = 0, Is1100 = 0, Is0010 = 0, Is1010 = 0,Is0110 = 0,Is1110 = 0,
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

####place of specific character
#strains
s_first<-3
s_last<-nchar(strains[1]) + 2
#treatments
t_last<-nchar(names(subset(state, str_sub(names(state), 1, 1) == "T")[1]))
t_first<-nchar(names(subset(state, str_sub(names(state), 1, 1) == "T")[1])) - 1

### Parameters 
infection = 1/120 * (615/100000)/12
infection = 1/120 * 10/12

parameters_values <- c( pi = 0.01856/12, #birth
                        
                        #infection 
                        alphas0000 = infection, alphas1000 = infection, alphas0100 = infection, alphas0010 = infection, 
                        alphas1100 = infection, alphas1010 = infection, alphas0110 = infection, alphas1110 = infection,
                        alphas0001 = infection, alphas1001 = infection, alphas0101 = infection, alphas0011 = infection, 
                        alphas1101 = infection, alphas1011 = infection, alphas0111 = infection, alphas1111 = infection,
                        
                        ####reinfection 
                        #L -> E 
                        Lambda1s0000=infection * 0.1,Lambda1s1000=infection * 0.1,Lambda1s0100=infection * 0.1,Lambda1s1100=infection * 0.1,
                        Lambda1s0010=infection * 0.1,Lambda1s1010=infection * 0.1,Lambda1s0110=infection * 0.1,Lambda1s1110=infection * 0.1,
                        Lambda1s0001=infection * 0.1,Lambda1s1001=infection * 0.1,Lambda1s0101=infection * 0.1,Lambda1s1101=infection * 0.1,
                        Lambda1s0011=infection * 0.1,Lambda1s1011=infection * 0.1,Lambda1s0111=infection * 0.1,Lambda1s1111=infection * 0.1,
                        
                        #R -> E
                        Lambda2s0000=infection * 0.1,Lambda2s1000=infection * 0.1,Lambda2s0100=infection * 0.1,Lambda2s1100=infection * 0.1,
                        Lambda2s0010=infection * 0.1,Lambda2s1010=infection * 0.1,Lambda2s0110=infection * 0.1,Lambda2s1110=infection * 0.1,
                        Lambda2s0001=infection * 0.1,Lambda2s1001=infection * 0.1,Lambda2s0101=infection * 0.1,Lambda2s1101=infection * 0.1,
                        Lambda2s0011=infection * 0.1,Lambda2s1011=infection * 0.1,Lambda2s0111=infection * 0.1,Lambda2s1111=infection * 0.1,
                        
                        
                        
                        delta1 = (10.25/(12*1000)) , delta2 = 2 / ( 12 * 1000), #progression to active disease
                        
                        epsilon = (1/30) * (1 - (10.25/(12*1000))) , #progression rate from early latent to late latent
                        
                        ####treated  rate
                        beta=0,
                        
                        ###test rate
                        sigma = 1, 
                        
                        ### good result for a test rate 
                        omega=0.5,
                        
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
                        
                        theta= 0.016/12, #death due to tuberculosis 
                        #(1. FigurRs of the dead: a decade of tuberculosis morta1sity registrations in South Africa. https://journa1ss.co.za/doi/epdf/10.7196/SAMJ.2019.v109i10.14073 doi:10.7196/SAMJ.2019.v109i10.14073.)
                        
                        phi=0.01/12, #re-activation from recovered to infected active stage
                        
                        
                        mu= 0.00926/12 ) #natura1s death


####
#### MATRICES #####
##row and col name# 
strains<-paste0(c("0","1"), rep(c("0","1"), each=2),rep(c("0","1"), each=4), rep(c("0","1"), each=8))
treatments<-c("t1","t2","t3", "t4")

### Resistance matrix ### (These matrices do not change over time )
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
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), ncol=16, byrow=TRUE)

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
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), ncol=16, byrow=TRUE)

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
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), ncol=16, byrow=TRUE)

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
                                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), ncol=16, byrow=TRUE)

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
                             0,0,0,0), ncol=4, byrow=TRUE)
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
                             0,0,0,0), ncol=4, byrow=TRUE)
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
                             0,0,0,0), ncol=4, byrow=TRUE)

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
                             1,1,1,0), ncol=4, byrow=TRUE)

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
                  1,0,0,0), ncol=4, byrow=TRUE)

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
                           0,0,0,1), ncol=4, byrow=TRUE)
colnames(good_test_result)<-treatments
rownames(good_test_result)<-strains


#if there is a test and the reslt are false
false_result<-matrix(c(1,1,0,1,
                       1,1,0,1,
                       1,1,0,1,
                       1,1,0,1,
                       1,1,0,1,
                       1,1,0,1,
                       1,1,0,1,
                       1,1,0,1,
                       1,1,0,1,
                       1,1,0,1,
                       1,1,0,1,
                       1,1,0,1,
                       1,1,0,1,
                       1,1,0,1,
                       1,1,0,1,
                       1,1,0,1), ncol=4, byrow=TRUE)

colnames(false_result)<-treatments
rownames(false_result)<-strains
nb_treatment<-1/length(treatments)

# Vector with parameters 
alpha<-subset(parameters_values, str_sub(names(parameters_values), 1, 5) == "alpha") #vector infection rate
lambda_1<- subset(parameters_values, str_sub(names(parameters_values), 1, 7) == "Lambda1") #vector re-infection rate for Late Latent
lambda_2<-subset(parameters_values, str_sub(names(parameters_values), 1, 7) == "Lambda2") #vector re-infection rate for Recovered
tau<-subset(parameters_values, str_sub(names(parameters_values), 1, 3) == "tau") #vector of treatment success rate

time_values <- seq(0, 1*12, by = 1) #time for simulation


