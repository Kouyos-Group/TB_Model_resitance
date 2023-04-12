#### write differentials equations ####
sir_equations <- function(time, state, parameters_values) {
  with(as.list(c(state, parameters_values)), {
    
### Parameters 
strains<-paste0(c("0","1"), rep(c("0","1"), each=2),rep(c("0","1"), each=4), rep(c("0","1"), each=8))
treatments<-c("t1","t2","t3", "t4")
####place of specific character
#strains
s_first<-3
s_last<-nchar(strains[1]) + 2
#treatments
t_last<-nchar(names(subset(state, str_sub(names(state), 1, 1) == "T")[1]))
t_first<-nchar(names(subset(state, str_sub(names(state), 1, 1) == "T")[1])) - 1


#matrix success treatment
treatments_success <- matrix(c(1,1,1,1,
                               0,1,1,1,
                               0,1,1,1,
                               0,1,1,1,
                               1,0,1,1,
                               0,0,1,1,
                               0,0,1,1,
                               0,0,1,1,
                               1,1,0,1,
                               0,1,0,1,
                               0,1,0,1,
                               0,1,0,1,
                               1,0,0,1,
                               0,0,0,1,
                               0,0,0,1,
                               0,0,0,1), ncol=4, byrow=TRUE)
colnames(treatments_success)<-treatments
rownames(treatments_success)<-strains

#matrix fail treatment
treatments_fail <- matrix(c(0,0,0,0,
                            1,0,0,0,
                            1,0,0,0,
                            1,0,0,0,
                            0,1,0,0,
                            1,1,0,0,
                            1,1,0,0,
                            1,1,0,0,
                            0,0,1,0,
                            1,0,1,0,
                            1,0,1,0,
                            1,0,1,0,
                            0,1,1,0,
                            1,1,1,0,
                            1,1,1,0,
                            1,1,1,0), ncol=4, byrow=TRUE)
colnames(treatments_fail)<-treatments
rownames(treatments_fail)<-strains

success_treatment_rate<- c(t1=1/6, t2=1/9, t3=1/2, t4=1/6)
tau_names<-paste0("tau", strains, rep(treatments, each=16))
tau<-sapply(tau_names, function(x) assign(x, t * treatments_success[substr(x, 4,7), substr(x, 8,9)] * success_treatment_rate[substr(x, 8, 9)]))
names(tau)<-tau_names
#pars <- c(pars, tau)

ita_rate <- c(t1=1/5.5, t2=1/9, t3=1/1.8, t4=(1/5.5))
ita_names<-paste0("ita", strains, rep(treatments, each=16))
ita<-sapply(tau_names, function(x) assign(x, treatments_fail[substr(x, 4,7), substr(x, 8,9)] * ita_rate[substr(x, 8, 9)]))
names(ita)<-ita_names
#pars <- c(pars, ita)

alpha_names<-paste0("alpha", strains)
alpha<-sapply(alpha_names, function(x) infection * (1- s)^ str_count(x, "1"))
names(alpha)<-alpha_names
#pars <- c(pars, alpha)

lambda_1_names<-paste0("lambda_1", strains)
lambda_1<-sapply(lambda_1_names, function(x) infection * (1- l1) * (1- s)^ str_count(x, "1"))
names(lambda_1)<-lambda_1_names
#pars <- c(pars, lambda_1)

lambda_2_names<-paste0("lambda_2", strains)
lambda_2<-sapply(lambda_2_names, function(x) infection * (1 - l2) * (1- s)^ str_count(x, "1"))
names(lambda_2)<-lambda_2_names


### New treatment if not adapted against the strain ###
#Here, we assume the shift it's always good

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
no_test<-matrix(c(0,0,1,0,
                  0,0,1,0,
                  0,0,1,0,
                  0,0,1,0,
                  0,0,1,0,
                  0,0,1,0,
                  0,0,1,0,
                  0,0,1,0,
                  0,0,1,0,
                  0,0,1,0,
                  0,0,1,0,
                  0,0,1,0,
                  0,0,1,0,
                  0,0,1,0,
                  0,0,1,0,
                  0,0,1,0), ncol=4, byrow=TRUE)

colnames(no_test)<-treatments
rownames(no_test)<-strains

#if there is a test and the reslt are good
good_test_result<-matrix(c(0,0,1,0,
                           0,0,1,0,
                           0,0,1,0,
                           0,0,1,0,
                           0,0,1,0,
                           0,0,1,0,
                           0,0,1,0,
                           0,0,1,0,
                           0,0,0,1,
                           0,0,0,1,
                           0,0,0,1,
                           0,0,0,1,
                           0,0,0,1,
                           0,0,0,1,
                           0,0,0,1,
                           0,0,0,1), ncol=4, byrow=TRUE)
colnames(good_test_result)<-treatments
rownames(good_test_result)<-strains


#if there is a test and the reslt are false
false_result<-matrix(c(sum(good_test_result[,1])/sum(good_test_result),sum(good_test_result[,2])/sum(good_test_result),sum(good_test_result[,3])/sum(good_test_result),sum(good_test_result[,4])/sum(good_test_result),
                       sum(good_test_result[,1])/sum(good_test_result),sum(good_test_result[,2])/sum(good_test_result),sum(good_test_result[,3])/sum(good_test_result),sum(good_test_result[,4])/sum(good_test_result),
                       sum(good_test_result[,1])/sum(good_test_result),sum(good_test_result[,2])/sum(good_test_result),sum(good_test_result[,3])/sum(good_test_result),sum(good_test_result[,4])/sum(good_test_result),
                       sum(good_test_result[,1])/sum(good_test_result),sum(good_test_result[,2])/sum(good_test_result),sum(good_test_result[,3])/sum(good_test_result),sum(good_test_result[,4])/sum(good_test_result),
                       sum(good_test_result[,1])/sum(good_test_result),sum(good_test_result[,2])/sum(good_test_result),sum(good_test_result[,3])/sum(good_test_result),sum(good_test_result[,4])/sum(good_test_result),
                       sum(good_test_result[,1])/sum(good_test_result),sum(good_test_result[,2])/sum(good_test_result),sum(good_test_result[,3])/sum(good_test_result),sum(good_test_result[,4])/sum(good_test_result),
                       sum(good_test_result[,1])/sum(good_test_result),sum(good_test_result[,2])/sum(good_test_result),sum(good_test_result[,3])/sum(good_test_result),sum(good_test_result[,4])/sum(good_test_result),
                       sum(good_test_result[,1])/sum(good_test_result),sum(good_test_result[,2])/sum(good_test_result),sum(good_test_result[,3])/sum(good_test_result),sum(good_test_result[,4])/sum(good_test_result),
                       sum(good_test_result[,1])/sum(good_test_result),sum(good_test_result[,2])/sum(good_test_result),sum(good_test_result[,3])/sum(good_test_result),sum(good_test_result[,4])/sum(good_test_result),
                       sum(good_test_result[,1])/sum(good_test_result),sum(good_test_result[,2])/sum(good_test_result),sum(good_test_result[,3])/sum(good_test_result),sum(good_test_result[,4])/sum(good_test_result),
                       sum(good_test_result[,1])/sum(good_test_result),sum(good_test_result[,2])/sum(good_test_result),sum(good_test_result[,3])/sum(good_test_result),sum(good_test_result[,4])/sum(good_test_result),
                       sum(good_test_result[,1])/sum(good_test_result),sum(good_test_result[,2])/sum(good_test_result),sum(good_test_result[,3])/sum(good_test_result),sum(good_test_result[,4])/sum(good_test_result),
                       sum(good_test_result[,1])/sum(good_test_result),sum(good_test_result[,2])/sum(good_test_result),sum(good_test_result[,3])/sum(good_test_result),sum(good_test_result[,4])/sum(good_test_result),
                       sum(good_test_result[,1])/sum(good_test_result),sum(good_test_result[,2])/sum(good_test_result),sum(good_test_result[,3])/sum(good_test_result),sum(good_test_result[,4])/sum(good_test_result),
                       sum(good_test_result[,1])/sum(good_test_result),sum(good_test_result[,2])/sum(good_test_result),sum(good_test_result[,3])/sum(good_test_result),sum(good_test_result[,4])/sum(good_test_result),
                       sum(good_test_result[,1])/sum(good_test_result),sum(good_test_result[,2])/sum(good_test_result),sum(good_test_result[,3])/sum(good_test_result),sum(good_test_result[,4])/sum(good_test_result)), ncol=4, byrow=TRUE)


colnames(false_result)<-treatments
rownames(false_result)<-strains

#vector of each state
vE  <-state[str_sub(names(state), 1, 1) == "E"] #vector with all the Early Latent types
vL  <-state[str_sub(names(state), 1, 1) == "L"] #vector with all the Late latent types
vI  <-state[str_sub(names(state), 1, 1) == "I"] #vector with all the Infected types
vT  <-state[str_sub(names(state), 1, 1) == "T"] #vector with all the Recovered types
vTs <-vT[str_sub(names(vT),-1,-1)=="1"] + vT[str_sub(names(vT),-1,-1)=="2"] + vT[str_sub(names(vT),-1,-1)=="3"] + 
  vT[str_sub(names(vT),-1,-1)=="4"] #vector with all strains of Treated (don't take account about treatment here)
vT_inf<-vT * c(treatments_fail[,1], treatments_fail[,2],treatments_fail[,3],treatments_fail[,4]) #vector with the treated individual allow to infected 
vTs_inf<-vT_inf[str_sub(names(vT_inf),-1,-1)=="1"] + vT_inf[str_sub(names(vT_inf),-1,-1)=="2"] + vT_inf[str_sub(names(vT_inf),-1,-1)=="3"] + 
  vT_inf[str_sub(names(vT_inf),-1,-1)=="4"] #vector with the treated individual allow to infected according to their treatments
vR  <-state[str_sub(names(state), 1, 1) == "R"] #vector with all the Recovered types

### Groups ###
S <-sum(state[str_sub(names(state), 1, 1) == "S"])
E <- sum(vE)
L <- sum(vL)
I <- sum(vI)
T <- sum(vT)
R <-sum(vR)

N=S+E+L+I+T+R

#ODE for Susceptible compartment
assign("dS",
       pi  * N #birth
       - mu * S #natural death 
       - sum((vTs_inf + vI) * alpha * S) #infection
)

#ODE for Early Latent stage compartments
for (i in 1: length(vE)){
  assign(paste0("d", names(vE)[i]), 
         (vTs_inf[i] + vI[i]) * alpha[i] * S + #infection from Suscepitble
           (vTs_inf[i] + vI[i]) * lambda_1[i] * L + #infection from Late latent
           (vTs_inf[i] + vI[i]) * lambda_2[i] * R - #infection from Recovered
           vE[i] * (epsilon + mu + delta1)  #transition to late latent, natural death, transition to active stage
  )
}

#ODE for Late Latent stage compartments
for (i in 1:length(vL)){
  assign(paste0("d", names(vL[i])),
         epsilon * vE[i] - #transition from early latent
           (mu + delta2)  * vL[i] - #death, and transition to active stage
           sum((vTs_inf + vI) * lambda_1 * vL[i]) #re-infection
  )
}

#ODE for Infected stage compartments
for (i in 1:length(vI)){
  assign(paste0("d", names(vI[i])), 
         delta1 * vE[i] +  #transition from Early latent
           delta2 * vL[i] + #transition from Late latent
           phi * vR[i] - #relaspe from Recovered
           (mu + theta + omicron + beta) * vI[i] #natural death, death due to tb, natural cured, go to treatment
  )
}

#ODE for Treated stage compartments
for (i in 1:length(vT)){
  assign(paste0("d", names(vT)[i]), 
         #natural death, death due to tb, natural cured, treatment cured
         - ( mu + theta + omicron + tau[i] + ita[i]) * vT[i] - 
           
           #transition to other strain due to resistance acquisition
           sum( rho * (resistance_matrix[[as.numeric(str_sub(names(vT[i]),-1))]] 
                       [str_sub(names(vT), s_first, s_last)[i],]) * vT[i]) +
           
           
           #transition from other strain due to resistance acquisition
           sum( rho * (resistance_matrix[[as.numeric(str_sub(names(vT[i]),-1))]]
                       [,str_sub(names(vT[i]), s_first, s_last)]) * 
                  subset(vT, str_sub(names(vT), -1) == str_sub(names(vT[i]), -1))) + 
           
           
           #transition from other treatment
           sum(subset(ita, str_sub(names(ita), 4, 7) == str_sub(names(vT[i]), s_first, s_last)) * 
                 ((eval(parse(text=paste0("new_treatment_",str_sub(names(vT),t_first, t_last)))[i])
                   [str_sub(names(vT[i]), s_first, s_last),]) 
                  *  subset(vT, str_sub(names(vT), s_first, s_last) == str_sub(names(vT[i]), s_first, s_last)))) +
           
           
           #start to treatment from infected stage
           beta * subset(vI, str_sub(names(vI), s_first, s_last) ==  str_sub(names(vT[i]), s_first, s_last)) * 
           ((1-sigma) * no_test[str_sub(names(vT)[i], s_first, s_last), 
                                str_sub(names(vT)[i], t_first, t_last)] +
              sigma * omega * good_test_result[str_sub(names(vT)[i], s_first, s_last), 
                                               str_sub(names(vT)[i], t_first, t_last)] +
              sigma * (1-omega) * false_result[str_sub(names(vT)[i], s_first, s_last), 
                                               str_sub(names(vT)[i], t_first, t_last)])
         
  )
}

#ODE for Recovered stage compartments
for (i in 1:length(vR)){
  assign(paste0("d", names(vR)[i]), 
         - ( phi + mu)  * vR[i]- #relapse, natural death
           #re-infection 
           sum((vTs_inf + vI) * lambda_2 * 
                 vR[i]) + 
           omicron * (vI[i] + vTs[i]) + #natural cured 
           
           #cured due to the treatment
           sum(subset(tau, str_sub(names(tau), 4, 7) == str_sub(names(vR)[i], s_first, s_last))   * 
                 subset(vT, str_sub(names(vT), s_first, s_last) == str_sub(names(vR)[i], s_first, s_last)))   
         
         
  )
  
  #ODE for death people due to tb (for likelihood)
  assign("dincidence", sum((vTs_inf + vI) * alpha * S) #infection 
         + sum((vTs_inf + vI) * lambda_1 * L) + #infection from Late latent
           sum((vTs_inf + vI) * lambda_2 * R)
         - state[str_sub(names(state), 1, 3) == "inc"]
  )
  
  #ODE for death people due to tb (for likelihood)
  assign("dMtb", theta * (I + T) - state[str_sub(names(state), 1, 1) == "M"])
  
  #ODE for perople infected with a rifampcin resistant strain (for likelihood)
  assign("drR", (delta1 * sum(vE[(str_sub(names(vE), 3 , 3) == "1" | str_sub(names(vE), 4 , 4) == "1")]) +  #transition from Early latent
                   delta2 * sum(vI[(str_sub(names(vI), 3 , 3) == "1" | str_sub(names(vI), 4 , 4) == "1")]))#transition from Late latent
         - state[str_sub(names(state), 1, 2) == "rR"])
  
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
              dRs0001, dRs1001, dRs0101, dRs1101, dRs0011, dRs1011, dRs0111, dRs1111, 
              
              dincidence,
              dMtb, 
              drR
              
)))
  })
}

time_values<-seq(2012, 2060, by = 1/12)
