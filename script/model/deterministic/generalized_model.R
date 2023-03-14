#### Load libraries ####
library(deSolve)
library(tidyverse)
library(stringr)
library(dplyr)

#### write differentials equations ####
sir_equations <- function(time, state, parameters_values) {
  with(as.list(c(state, parameters_values)), {
    
    #vector of each state
    vE  <-state[str_sub(names(state), 1, 1) == "E"] #vector with all the Early Latent types
    vL  <-state[str_sub(names(state), 1, 1) == "L"] #vector with all the Late latent types
    vI  <-state[str_sub(names(state), 1, 1) == "I"] #vector with all the Infected types
    vT  <-state[str_sub(names(state), 1, 1) == "T"] #vector with all the Recovered types
    vTs <-vT[str_sub(names(vT),-1,-1)=="1"] + vT[str_sub(names(vT),-1,-1)=="2"] + vT[str_sub(names(vT),-1,-1)=="3"] + 
      vT[str_sub(names(vT),-1,-1)=="4"] #vector with all strains of Treated (don't take account about treatment here)
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
           - mu * S #naural death 
           - sum((vTs + vI) * alpha * S) #infection
    )
  
    #ODE for Early Latent stage compartments
    for (i in 1: length(vE)){
      assign(paste0("d", names(vE)[i]), 
             (vTs[i] + vI[i]) * alpha[i] * S + #infection from Suscepitble
               (vTs[i] + vI[i]) * lambda_1[i] * L + #infection from Late latent
               (vTs[i] + vI[i]) * lambda_2[i] * R - #infection from Suscepitble
               vE[i] * (epsilon + mu + delta1)  #transition to late latent, death, transition to active stage
      )
    }
    
    #ODE for Late Latent stage compartments
    for (i in 1:length(vL)){
      assign(paste0("d", names(vL[i])),
             epsilon * vE[i] - #transition from early latent
               (mu + delta2)  * vL[i] - #death, and transition to active stage
               sum((vTs + vI) * lambda_1 * vL[i]) #infection
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
             - ( mu + theta + omicron ) * vT[i] - 
               #natural death, death due to tb, natural cured
               
               sum( rho * (resistance_matrix[[as.numeric(str_sub(names(vT[i]),-1))]] 
                           [str_sub(names(vT), s_first, s_last)[i],]) * vT[i]) +
               #transition to other strain due to resistance acquisition
               
               sum( rho * (resistance_matrix[[as.numeric(str_sub(names(vT[i]),-1))]]
                           [,str_sub(names(vT[i]), s_first, s_last)]) * 
                      subset(vT, str_sub(names(vT), -1) == str_sub(names(vT[i]), -1))) -
               #transition from other strain due to resistance acquisition
               
               sum( ita *  sapply(new_treatment, function(x) sum(x[str_sub(names(vT), s_first, s_last)[i], str_sub(names(vT), t_first, t_last)[i]])) * vT[i]) +
               
               #transition to other treatment
            
               sum(ita * ((eval(parse(text=paste0("new_treatment_",str_sub(names(vT),t_first, t_last)))[i])
                           [str_sub(names(vT[i]), s_first, s_last),]) 
                          *  subset(vT, str_sub(names(vT), s_first, s_last) == str_sub(names(vT[i]), s_first, s_last)))) +
               #transition from other treatment
               
               beta * subset(vI, str_sub(names(vI), s_first, s_last) ==  str_sub(names(vT[i]), s_first, s_last)) * 
               ((1-sigma) * no_test[str_sub(names(vT)[i], s_first, s_last), 
                                 str_sub(names(vT)[i], t_first, t_last)] +
                  sigma * omega * good_test_result[str_sub(names(vT)[i], s_first, s_last), 
                                                   str_sub(names(vT)[i], t_first, t_last)] +
                  nb_treatment * sigma * (1-omega) * false_result[str_sub(names(vT)[i], s_first, s_last), 
                                                                str_sub(names(vT)[i], t_first, t_last)])
             #transition from infected stage
      )
    }
    
    #ODE for Recovered stage compartments
    for (i in 1:length(vR)){
      assign(paste0("d", names(vR)[i]), 
            - ( phi + mu)  * vR[i] - #relapse, natural death
              sum((vTs + vI)* lambda_2 * 
                     vR[i]) + #re-infection 
               omicron * (vI[i] + vTs[i]) + #natural cured 
               
               sum(subset(tau, str_sub(names(tau), 5, 8) == str_sub(names(vR)[i], s_first, s_last))   * 
                     subset(vT, str_sub(names(vT), s_first, s_last) == str_sub(names(vR)[i], s_first, s_last))) 
            #cured due to the treatment 
      )
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
