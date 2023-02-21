##### 70-120 years #####
###Initial States 

state2 <-as.vector(tail(sir_values_1[,2:130], n=1))
names(state2)<-names(state)


state2[str_sub(names(state2), 1, 1) == "T" & str_sub(names(state2), -1, -1) == "3"] <- state2[str_sub(names(state2), 1, 1) == "T" & str_sub(names(state2), -1, -1) == "2"]
state2[str_sub(names(state2), 1, 1) == "T" & str_sub(names(state2), -1, -1) == "2"] <- 0


### Parameters 
infection = 0.05 * (615/100000)/12

parameters_values <- c( pi = 0.01856/12, #birth
                        
                        #infection 
                        alphas0000 = infection, alphas1000 = infection, alphas0100 = infection, alphas0010 = infection, 
                        alphas1100 = infection, alphas1010 = infection, alphas0110 = infection, alphas1110 = infection,
                        alphas0001 = infection, alphas1001 = infection, alphas0101 = infection, alphas0011 = infection, 
                        alphas1101 = infection, alphas1011 = infection, alphas0111 = infection, alphas1111 = infection,
                        
                        ####reinfection 
                        #L -> E 
                        Lambda1s0000=infection,Lambda1s1000=infection,Lambda1s0100=infection,Lambda1s1100=infection,
                        Lambda1s0010=infection,Lambda1s1010=infection,Lambda1s0110=infection,Lambda1s1110=infection,
                        Lambda1s0001=infection,Lambda1s1001=infection,Lambda1s0101=infection,Lambda1s1101=infection,
                        Lambda1s0011=infection,Lambda1s1011=infection,Lambda1s0111=infection,Lambda1s1111=infection,
                        
                        #R -> E
                        Lambda2s0000=infection,Lambda2s1000=infection,Lambda2s0100=infection,Lambda2s1100=infection,
                        Lambda2s0010=infection,Lambda2s1010=infection,Lambda2s0110=infection,Lambda2s1110=infection,
                        Lambda2s0001=infection,Lambda2s1001=infection,Lambda2s0101=infection,Lambda2s1101=infection,
                        Lambda2s0011=infection,Lambda2s1011=infection,Lambda2s0111=infection,Lambda2s1111=infection,
                        
                        
                        epsilon = 8/12, #progression rate from early latent to late latent
                        
                        
                        delta1 = 0.4/12, delta2 = 0.0020/12, #progression to active disease
                        
                        ####treated  rate
                        beta=0.001,
                        
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
                        
                        theta= 0.026/12, #death due to tuberculosis 
                        #(1. FigurRs of the dead: a decade of tuberculosis morta1sity registrations in South Africa. https://journa1ss.co.za/doi/epdf/10.7196/SAMJ.2019.v109i10.14073 doi:10.7196/SAMJ.2019.v109i10.14073.)
                        
                        phi=0.01/12, #re-activation from recovered to infected active stage
                        
                        
                        mu= 0.00926/12 ) #natura1s death


### Matrices specifiy to 70 first years 
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

##the adapted treatment is the third one
new_treatment_t2 <- matrix(c(0,0,0,0,
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

colnames(new_treatment_t2)<-treatments
rownames(new_treatment_t2)<-strains

#the adapted treatment is the second one
new_treatment_t3 <- matrix(c(0,0,0,0,
                             1,1,0,1,
                             1,1,0,1,
                             1,1,0,1,
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
                           0,0,1,0,
                           0,0,1,0,
                           0,0,1,0,
                           1,0,0,0,
                           0,0,0,1,
                           0,0,0,1,
                           0,0,0,1,
                           1,0,0,0,
                           0,0,1,0,
                           0,0,1,0,
                           0,0,1,0,
                           1,0,0,0,
                           0,0,0,1,
                           0,0,0,1,
                           0,0,0,1), ncol=4, byrow=TRUE)
colnames(good_test_result)<-treatments
rownames(good_test_result)<-strains


#if there is a test and the reslt are false
false_result<-matrix(c(1,0,1,1,
                       1,0,1,1,
                       1,0,1,1,
                       1,0,1,1,
                       1,0,1,1,
                       1,0,1,1,
                       1,0,1,1,
                       1,0,1,1,
                       1,0,1,1,
                       1,0,1,1,
                       1,0,1,1,
                       1,0,1,1,
                       1,0,1,1,
                       1,0,1,1,
                       1,0,1,1,
                       1,1,1,1), ncol=4, byrow=TRUE)

colnames(false_result)<-treatments
rownames(false_result)<-strains
nb_treatment<-1/length(treatments)


# Vector with parameters 
alpha<-subset(parameters_values, str_sub(names(parameters_values), 1, 5) == "alpha")
lambda_1<- subset(parameters_values, str_sub(names(parameters_values), 1, 7) == "Lambda1")
lambda_2<-subset(parameters_values, str_sub(names(parameters_values), 1, 7) == "Lambda2")
tau<-subset(parameters_values, str_sub(names(parameters_values), 1, 3) == "tau")



#  the points in time where to calculate variables values 
time_values <- seq(24, 48, by = 1) 
####