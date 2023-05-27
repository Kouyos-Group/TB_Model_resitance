#### write differentials equations ####
sir_equations <- function(time, state, parameters_values) {
  with(as.list(c(state, parameters_values)), {
    
    #### strains, treatment #### 
    strains <- paste0(c("0","1"), rep(c("0","1"), each=2))
    treatments <- c("t1","t2","t3")
    ##place of specific character
    #strains
    s_first <- 3
    s_last <- nchar(strains[1]) + 2
    
    #treatments
    t_first <- s_last+1
    t_last <- s_last+2
    
    #### Matrices #####
    #matrix success treatment
    treatments_success <- matrix(c(1,1,1,
                                   0,1,1,
                                   1,0,1,
                                   0,0,1), ncol=3, byrow=TRUE)
    
    colnames(treatments_success) <- treatments
    rownames(treatments_success) <- strains
    
    #matrix fail treatment
    treatments_fail <- ifelse(treatments_success == 0, 1, 0)
    
    colnames(treatments_fail)<-treatments
    
    rownames(treatments_fail)<-strains
    
    
    #treatment success vector
    t <- c(cured_t1,cured_t2,cured_t3)
    names(t) <- treatments
    time_treatment <- c(time_1, time_2, time_3)
    names(time_treatment) <- c("t1", "t2", "t3")
    tau_names<-paste0("tau", strains, rep(treatments, each=4))
    tau<-sapply(tau_names, function(x) assign(x,treatments_success[substr(x, 4,5), substr(x, 6,7)] *  time_treatment[substr(x, 6, 7)]))
    names(tau)<-tau_names
    
    #vector treatment shift due to resistance (at the end of the treatment) 
    ita_names<-paste0("ita", strains, rep(treatments, each=4))
    ita<-sapply(tau_names, function(x) assign(x, treatments_fail[substr(x, 4,5), substr(x, 6,7)] *  time_treatment[substr(x, 6, 7)]))
    names(ita)<-ita_names
    
    #vector Susceptible infection
    alpha_names<-paste0("alpha", strains)
    alpha<-sapply(alpha_names, function(x) infection * (1 - resistance_cost)^ str_count(x, "1"))
    names(alpha)<-alpha_names
    
    #vector Late latent infection
    lambda_names<-paste0("lambda", strains)
    lambda<-sapply(lambda_names, function(x) infection * (1 - l) * (1- resistance_cost)^ str_count(x, "1"))
    names(lambda)<-lambda_names
    
    
    ### New treatment if not adapted against the strain ###
    #Here, we assume the shift it's always good
    
    ### Resistance matrix ### (These matrices do not change over time )
    #resistance acquisition with 1st treatment
    resistance_matrix_t1 <- matrix(c(0,1,0,0,
                                     0,0,0,0,
                                     0,0,0,1,
                                     0,0,0,0), ncol=4, byrow=TRUE)
    
    colnames(resistance_matrix_t1)<-strains
    rownames(resistance_matrix_t1)<-strains
    
    #resistance acquisition with 2nd treatment
    resistance_matrix_t2 <- matrix(c(0,0,1,0,
                                     0,0,0,1,
                                     0,0,0,0,
                                     0,0,0,0), ncol=4, byrow=TRUE)
    
    colnames(resistance_matrix_t2)<-strains
    rownames(resistance_matrix_t2)<-strains
    
    #resistance acquisition with third treatment
    resistance_matrix_t3 <- matrix(c(0,0,0,0,
                                     0,0,0,0,
                                     0,0,0,0,
                                     0,0,0,0), ncol=4, byrow=TRUE)
    
    colnames(resistance_matrix_t3)<-strains
    rownames(resistance_matrix_t3)<-strains
    
    #list with all resistance matrices
    resistance_matrix<-list(resistance_matrix_t1, resistance_matrix_t2, resistance_matrix_t3)
    
    
    #start_treatment (according to the treatment)
    #change according to the scenario
    start_treatment <- matrix(c(1,0,0,
                                0,0,1,
                                1,0,0,
                                0,0,1), ncol=3, byrow=TRUE)
    
    colnames(start_treatment)<-treatments
    
    rownames(start_treatment)<-strains
    
    #continuity of treatment
    c_treat_t1 <- matrix(c(1,1,1,
                           0,0,0,
                           1,1,1,
                           0,0,0), ncol=3, byrow=TRUE)
    
    colnames(c_treat_t1)<-treatments
    
    rownames(c_treat_t1)<-strains
    
    c_treat_t2 <- matrix(c(0,0,0,
                           0,0,0,
                           0,0,0,
                           0,0,0), ncol=3, byrow=TRUE)
    
    colnames(c_treat_t2)<-treatments
    
    rownames(c_treat_t2)<-strains
    
    c_treat_t3 <- matrix(c(0,0,0,
                           1,1,1,
                           0,0,0,
                           1,1,1), ncol=3, byrow=TRUE)
    
    colnames(c_treat_t3)<-treatments
    
    rownames(c_treat_t3)<-strains

    continuous_treatment<-list(c_treat_t1, c_treat_t2, c_treat_t3)
    #list with all new treattment matrices
    #new_treatment<-list(new_treatment_t1, new_treatment_t2, new_treatment_t3)
    
    #vector of each state
    vE  <-state[str_sub(names(state), 1, 1) == "E"] #vector with all the Early Latent types
    vL  <-state[str_sub(names(state), 1, 1) == "L"] #vector with all the Late latent types
    vI  <-state[str_sub(names(state), 1, 1) == "I"] #vector with all the Infected types
    vT  <-state[str_sub(names(state), 1, 1) == "T"] #vector with all the Recovered types
    vTs <-vT[str_sub(names(vT),-1,-1)=="1"] + vT[str_sub(names(vT),-1,-1)=="2"] + vT[str_sub(names(vT),-1,-1)=="3"] #vector with all strains of Treated (don't take account about treatment here)
    vT_inf<-vT * c(treatments_fail[,1], treatments_fail[,2],treatments_fail[,3]) #vector with the treated individual allow to infected 
    vTs_inf<-vT_inf[str_sub(names(vT_inf),-1,-1)=="1"] + vT_inf[str_sub(names(vT_inf),-1,-1)=="2"] + vT_inf[str_sub(names(vT_inf),-1,-1)=="3"]  #vector with the treated individual allow to infected according to their treatments
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
    for (i in 1:length(vE)){
      assign(paste0("d", names(vE)[i]), 
             (vTs_inf[i] + vI[i]) * alpha[i] * S  #infection from Suscepitble
             + (vTs_inf[i] + vI[i]) * lambda[i] * L  #infection from Late latent
             + (vTs_inf[i] + vI[i]) * lambda[i] * R  #infection from Recovered
             -  vE[i] * (epsilon+ mu + delta1)  #transition to late latent, natural death, transition to active stage
      )
    }
    
    #ODE for Late Latent stage compartments
    for (i in 1:length(vL)){
      assign(paste0("d", names(vL[i])),
             epsilon * vE[i] #transition from early latent
             - (mu + delta2)  * vL[i] #death, and transition to active stage
             - sum((vTs_inf + vI) * lambda * vL[i]) #re-infection
      )
    }
    
    #ODE for Infected stage compartments
    for (i in 1:length(vI)){
      assign(paste0("d", names(vI[i])),  
             delta1 * vE[i]  #transition from Early latent
             + delta2 * vL[i]  #transition from Late latent
             - (mu + theta + omicron + beta) * vI[i] #natural death, death due to tb, natural cured, go to treatment
            # + phi * vR[i] #relapse
      )
    }
    
    #ODE for Treated stage compartments
    for (i in 1:length(vT)){
      assign(paste0("d", names(vT)[i]), 
             #natural death, death due to tb, natural cured, shift tratment due ot resistance, end of treatment and no resistance
             - (mu + theta + omicron + ita[i] + tau[i]) * vT[i] 
             
             #start to treatment from infected stage
             +  (beta * subset(vI, str_sub(names(vI), s_first, s_last) ==  str_sub(names(vT[i]), s_first, s_last)) * 
                   start_treatment[str_sub(names(vT)[i], s_first, s_last), str_sub(names(vT)[i], t_first, t_last)])
             
             #transition to other strain due to resistance acquisition
             - sum( rho * (resistance_matrix[[as.numeric(str_sub(names(vT[i]),-1))]] 
                           [str_sub(names(vT), s_first, s_last)[i],]) * vT[i]) 
             
             
             #transition from other strain due to resistance acquisition
             + sum( rho * (resistance_matrix[[as.numeric(str_sub(names(vT[i]),-1))]]
                           [,str_sub(names(vT[i]), s_first, s_last)]) * 
                      subset(vT, str_sub(names(vT), -1) == str_sub(names(vT[i]), -1))) 
             
             #fail treatment
             #treatment appropriate
             + vT[i] * (1-t[str_sub(names(vT[i]), 5,6)]) * tau[i]
             
             #treatment not appropriate
            + sum(subset(ita, str_sub(names(ita), s_first+1, s_last+1) == str_sub(names(vT[i]), s_first, s_last)) *
                     (continuous_treatment[[as.numeric(str_sub(names(vT[i]),-1))]][str_sub(names(vT[i]), s_first, s_last),]) *
                     subset(vT, str_sub(names(vT), s_first, s_last) == str_sub(names(vT[i]), s_first, s_last)))
             
      )
    }
    
    #ODE for Recovered stage compartments
    for (i in 1:length(vR)){
      assign(paste0("d", names(vR)[i]), 0
             - mu  * vR[i]  #natural death
             #reinfection  
             - sum((vTs_inf + vI) * lambda * 
                     vR[i])  
             
             + omicron * (vI[i] + vTs[i]) #natural cured 
             
             #cured due to the treatment
             + sum(subset(tau, str_sub(names(tau), s_first+1, s_last+1) == str_sub(names(vR)[i], s_first, s_last)) * 
                     t * subset(vT,str_sub(names(vT), s_first, s_last) == str_sub(names(vR)[i], s_first, s_last)))
            
      )
    }
    
    
    #ODE for death people due to tb (for likelihood)
    assign("dincidence", 
           sum((vTs_inf + vI) * alpha * S) #infection 
           + sum((vTs_inf + vI) * lambda * L) + #infection from Late latent
             sum((vTs_inf + vI) * lambda * R)
           - incidence
    )
    
    #ODE for death people due to tb (for likelihood)
    assign("dM", 
           theta * (I + T) - M 
    )
    
    #ODE for people infected with a rifampcin resistant strain (for likelihood)
    assign("drR", 
           (sum((vTs_inf[(str_sub(names(alpha), s_first+3 , s_first+3) == "1")] + 
                   vI[(str_sub(names(alpha), s_first+3 , s_first+3) == "1")]) * 
                  S * alpha[(str_sub(names(alpha), s_first+3 , s_first+3) == "1")])
            
            #infection from Late latent
            +  (sum((vTs_inf[(str_sub(names(alpha), s_first+3 , s_first+3) == "1")] + 
                       vI[(str_sub(names(alpha), s_first+3 , s_first+3) == "1")]) * 
                      lambda[(str_sub(names(lambda), s_first+4 , s_first+4) == "1")] * L)) 
            
            #infection from Recovered
            + (sum((vTs_inf[(str_sub(names(alpha), s_first+3 , s_first+3) == "1")] + 
                      vI[(str_sub(names(alpha), s_first+3 , s_first+3) == "1")]) * 
                     lambda[(str_sub(names(lambda), s_first+4 , s_first+4) == "1")] * R)))
           - rR 
    )
    
    
    
    return(list(c(dS, 
                  
                  dEs00, dEs10, dEs01, dEs11,
                  
                  dLs00, dLs10, dLs01, dLs11,
                  
                  dIs00, dIs10, dIs01, dIs11, 
                  
                  dTs00t1, dTs10t1, dTs01t1, dTs11t1, 
                  dTs00t2, dTs10t2, dTs01t2, dTs11t2, 
                  dTs00t3, dTs10t3, dTs01t3, dTs11t3, 
                  
                  dRs00, dRs10, dRs01, dRs11, 
                  
                  dincidence,
                  dM, 
                  drR
                  
    )))
  })
}
