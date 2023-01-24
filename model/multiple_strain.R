###Load packages
library(TiPS) #Tuto to use TiPS: https://cran.r-project.org/web/packages/TiPS/vignettes/TiPS.html


###Write Reactions
#####
#INFECTED STRAIN (0: not resistant, 1:resistant)
###first X--: rifampicin resitant
###second -X-: isoniazid resistant
##third --X: belaquilin resistant
#paste0("s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4))

#TREATMENT
#t0=>no treatment
#t1=>old treatment
#t2=>short-course
#paste0("t",c("0","1","2"))

#Latent compartments
#paste0("L","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4))
#Infected compartments
#paste0("I","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4),"t",rep(c("0","1","2"), each=8))

#####
reactions <- c(
  #####birth
  paste0("0 [ pi * ", c("S",paste0("L","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4)),paste0("E","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4)),
                        paste0("I","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4),"t",rep(c("0","1","2"), each=8)),
                        paste0("R","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4),"t",rep(c("0","1","2"), each=8))), " ] -> S"),
  
  
  ####infection
  paste0("S [ S * alpha","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), " * I", "s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4),
         "t",rep(c("0","1","2"), each=8), " ] -> E","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4)),
  
  
  ####Progression from early latent to late latent
  paste0("E","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), " [ E","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), " * epsilon ] -> L",
         "s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4)),
  
  
  ####disease progression from early latent to active disease
  paste0("E","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4)," [ ", "E","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), " * (1-psi) * delta1 ] -> I",
         "s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), "t0"), #not detected
  
  paste0("E","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4)," [ ", "E","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), " * psi * delta1 * (1-sigma) ] -> I",
         "s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), "t1"), #detected but not test to know with which strain
  
  paste0("E","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4)," [ ", "E","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), " * psi * delta1 * psi * (1-omega) * iota", 
         rep(c("0","1","2"),each=8)," ] -> I","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), "t", rep(c("0","1","2"),each=8)), #detected, tested but wrong resuEt
  
  paste0("E","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4)," [ ", "E","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), " * psi * delta1 * sigma * omega ] -> I",
         "s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), "t", c("1","2","2","2","1","0","0","0")), #detected, tested, good result 
  
  
  ####re-infection from late latent to early latent stage
  paste0("L","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), " [ L","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), " * lambda1 ] -> E",
         "s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4)),
  
  
  ####disease progression from late latent to active disease
  paste0("L","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4)," [ ", "L","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), " * (1-psi) * delta2 ] -> I",
         "s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), "t0"), #not detected
  
  paste0("L","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4)," [ ", "L","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), " * psi * delta2 * (1-sigma) ] -> I",
         "s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), "t1"), #detected but not test to know with which strain
  
  paste0("L","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4)," [ ", "L","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), " * psi * delta2 * psi * (1-omega) * iota", 
         rep(c("0","1","2"),each=8)," ] -> I","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), "t", rep(c("0","1","2"),each=8)), #detected, tested but wrong result
  
  paste0("L","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4)," [ ", "L","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), " * psi * delta2 * sigma * omega ] -> I",
         "s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4), "t", c("1","2","2","2","1","0","0","0")), #detected, tested, good result
  
  
  ####successful treatment
  paste0("I","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4),"t",rep(c("0","1","2"), each=8),
         " [ tau","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4),"t",rep(c("0","1","2"), each=8),
         " * I","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4),"t",rep(c("0","1","2"), each=8), " ] -> R",
         "s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4),"t",rep(c("0","1","2"), each=8)),
  
  
  ####cured by itself
  paste0("I","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4),"t",rep(c("0","1","2"), each=8),
         " [ omicron * I","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4),"t",rep(c("0","1","2"), each=8), " ] -> R",
         "s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4),"t",rep(c("0","1","2"), each=8)),
  
  
  ####treatment fail, strains become resistant
  #acquisition rifampicin resistnace
  paste0("I","s0",c("0","1"),rep(c("0","1"),each=2), "t1 [ ","I","s0",c("0","1"),rep(c("0","1"),each=2),
         "t1", " * rho","100 ] -> Is1",c("0","1"),rep(c("0","1"),each=2), "t1" ),
    
  #acquisistion isoniazid
  paste0("I","s",c("0","1"),"0",rep(c("0","1"),each=2), "t1", " [ ","I","s",c("0","1"),"0",rep(c("0","1"),each=2),
         "t1 * rho010 ] -> I","s",c("0","1"),"1",rep(c("0","1"),each=2), "t1"),
    
  #acquisistion short course resistance
  paste0("I","s",c("0","1"),rep(c("0","1"),each=2),"0", "t1", " [ ","I","s",c("0","1"),rep(c("0","1"),each=2),"0",
         "t1 * rho001 ] -> Is",c("0","1"),rep(c("0","1"),each=2),"1","t1"),
  
  
  ####death due to TB
  paste0("I","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4),"t",rep(c("0","1","2"), each=8),
         " [ theta", " * I","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4),"t",rep(c("0","1","2"), each=8),
         " ] -> S"),
  
  
  ####Reactivation from cured to active stage
  paste0("R","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4),"t",rep(c("0","1","2"), each=8),
         " [ phi * I","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4),"t",rep(c("0","1","2"), each=8), " ] -> I",
         "s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4),"t",rep(c("0","1","2"), each=8)),
  
  
  ####Reinfection from cured to early latent stage
  paste0("R","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4),"t",rep(c("0","1","2"), each=8),
         " [ R","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4),"t",rep(c("0","1","2"), each=8), " * lambda2 ] -> E",
         "s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4)),

  
  ####natural death
  paste0(c("S",paste0("E","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4)),
           paste0("L","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4)),
           paste0("I","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4),"t",rep(c("0","1","2"), each=8)),
           paste0("R","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4),"t",rep(c("0","1","2"), each=8))),
         " [  mu * ", c("S",paste0("E","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4)),
                        paste0("L","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4)),
                        paste0("I","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4),"t",rep(c("0","1","2"), each=8)),
                        paste0("R","s",c("0","1"),rep(c("0","1"),each=2),rep(c("0","1"),each=4),"t",rep(c("0","1","2"), each=8)))," ] -> S")
  
  
  
)
#####

###Build the simulator that will allow to run multiple trajectories
sir_simu <- build_simulator(reactions)

#####
###Initial States


initialStates <- c(S =88513, 
                   
                   Es000 = 57892, Es100 = 0, Es010 = 0, Es110 = 0, Es001 = 0, Es101 = 0,Es011 = 0,Es111 = 0,
                   
                   Ls000 = 3259, Ls100 = 0, Ls010 = 0, Ls110 = 0, Ls001 = 0, Ls101 = 0,Ls011 = 0,Ls111 = 0,
                   
                   Is000t0 = 176, Is100t0 = 0, Is010t0 = 0, Is110t0 = 0, Is001t0 = 0, Is101t0 = 0,Is011t0 = 0,Is111t0 = 0,
                   Is000t1 = 0, Is100t1 = 0, Is010t1 = 0, Is110t1 = 0, Is001t1 = 0, Is101t1 = 0,Is011t1 = 0,Is111t1 = 0,
                   Is000t2 = 0, Is100t2 = 0, Is010t2 = 0, Is110t2 = 0, Is001t2 = 0, Is101t2 = 0,Is011t2 = 0,Is111t2 = 0,
                   
                   Rs000t0 = 0, Rs100t0 = 0, Rs010t0 = 0, Rs110t0 = 0, Rs001t0 = 0, Rs101t0 = 0,Rs011t0 = 0,Rs111t0 = 0,
                   Rs000t1 = 0, Rs100t1 = 0, Rs010t1 = 0, Rs110t1 = 0, Rs001t1 = 0, Rs101t1 = 0,Rs011t1 = 0,Rs111t1 = 0,
                   Rs000t2 = 0, Rs100t2 = 0, Rs010t2 = 0, Rs110t2 = 0, Rs001t2 = 0, Rs101t2 = 0,Rs011t2 = 0,Rs111t2 = 0
                   
)
#####

###Duration of the simuI_isorn = 0, I_isoriso=0, I_isorr=0, I_isore=0, I_isorp=0, I_isorb=0, I_isord=0ation
time <- c(0,70, 800)
dT <- round(1/365, 4) # use a daily time step

#####
###Parameters value
theta <- list(pi = 0.01856/365.25, #birth
             
              alphas000 = 3.779603e-05, alphas100 = 1.779603e-05, alphas010 = 1.779603e-05, alphas001 = 1.779603e-05, alphas110 = 1.779603e-05,
              alphas101 = 1.779603e-05, alphas011 = 1.779603e-05, alphas111 = 1.779603e-05, #infection 
              
              #alphas000 = 1.779603e-05, alphas100 = 1.779603e-05, alphas010 = 1.779603e-05, alphas001 = 1.779603e-05, alphas110 = 1.779603e-05,
              #alphas101 = 1.779603e-05, alphas011 = 1.779603e-05, alphas111 = 1.779603e-05, #infection 
              
              epsilon = 0.2/365.25, #progression rate from early latent to late latent
             
              delta1 = 0.01/365.25, delta2 = 0.0008/365.25, #progression to active disease
              
              lambda1=0.21/365.25, lambda2=0.21/365.25, #reinfection to early latent stage
              
              psi=1, #detected rate
              
              sigma=c(0,1), #tested to know which strain is it
              
              omega=1, #test is right
              
              iota0=1/3, iota1=1/3, iota2=1/3, #rate wrong result for the 0,1,2 treatment
              
              taus000t0 = 0, taus100t0 = 0,taus010t0 = 0,taus001t0 = 0,taus110t0 = 0,taus101t0 = 0,taus011t0 = 0,taus111t0 = 0, #cured with no treatment
              taus000t1 = 0.005479452, taus100t1 = 0.0, taus010t1 = 0, taus001t1 = 0.005479452,taus110t1 = 0,taus101t1 = 0,taus011t1 = 0,taus111t1 = 0, #cured with old treatment
              taus000t2 = 0.005479452,taus100t2 = 0.005479452,taus010t2 = 0.005479452,taus001t2 = 0,taus110t2 = 0.005479452,taus101t2 = 0,taus011t2 = 0,taus111t2 = 0, #cured with short course treatment
              
              omicron=0,   #0.2/365.25, #self cured
              
              phi=0.01/365.25, #re-activation from recovered to infected active stage
              
              rho100= 30/365.25, rho010=30/365.25, rho001=30/365.25, #acquisition of a new resistance
              #rho100=0.8/365.25, rho010=0.8/365.25, rho001=0.08/365.25, #acquisition of a new resistance
              
              theta= 0.06/365.25, #death due to tuberculosis 
              #(1. Figures of the dead: a decade of tuberculosis mortality registrations in South Africa. https://journals.co.za/doi/epdf/10.7196/SAMJ.2019.v109i10.14073 doi:10.7196/SAMJ.2019.v109i10.14073.)
              
              mu= 0.00926) #natural death

#####
###In some simulations, the population size of a deme compartment may be zero before the upper time limit is reached, because of stochasticity or parameter values. 
#In this case, the simulation is considered to have failed and is halted.
#To bypass these failures, we can define the following wrapper:
safe_run <- function(f, ...) {
  out <- list()
  while(! length(out)) {out <- f(...)}
  out
}

#A safe version of our simulator sir_simu() is then:
safe_sir_simu <- function(...) safe_run(sir_simu, ...)



###Direct method
traj_dm <- safe_sir_simu(
  paramValues = theta,
  initialStates = initialStates,
  times = time,
  tau = dT,
  method = "mixed",
  nTrials=5)
#####
colnames(traj_dm$traj)
tail(traj_dm$traj)

plot(traj_dm$traj[,1],log(traj_dm$traj[,4]), type = "l", ylim=c(0,12) , xlab="time (years)", ylab="number of individuals")
lines(traj_dm$traj[,1], log(traj_dm$traj[,21]), type= "l", col = "blue")
lines(traj_dm$traj[,1], log(traj_dm$traj[,13]), type= "l", col = "purple")
lines(traj_dm$traj[,1], log(traj_dm$traj[,5]), type= "l", col = "red")
lines(traj_dm$traj[,1], log(traj_dm$traj[,53]), type= "l", col = "green")

plot(traj_dm$traj[,1], traj_dm$traj[,49], type= "l", col = "green")


layout(matrix(1:4,2,2))
plot(traj_dm$traj[,1],traj_dm$traj[,4], type = "l")
plot(traj_dm$traj[,1], traj_dm$traj[,21], type= "l", col = "blue")
plot(traj_dm$traj[,1], traj_dm$traj[,13], type= "l", col = "purple")
plot(traj_dm$traj[,1], traj_dm$traj[,5], type= "l", col = "red")



Is000<-as.vector(rowSums(traj_dm$traj[,c(21,29,37)]))
Is100<-as.vector(rowSums(traj_dm$traj[,c(22,30,38)]))
Is010<-as.vector(rowSums(traj_dm$traj[,c(23,31,39)]))
Is110<-as.vector(rowSums(traj_dm$traj[,c(24,32,40)]))
Is001<-as.vector(rowSums(traj_dm$traj[,c(25,33,41)]))
Is101<-as.vector(rowSums(traj_dm$traj[,c(26,34,42)]))
Is011<-as.vector(rowSums(traj_dm$traj[,c(27,35,43)]))
Is111<-as.vector(rowSums(traj_dm$traj[,c(28,36,44)]))

plot(traj_dm$traj[,1], Is000, type= "l", col = "blue")
lines(traj_dm$traj[,1], Is100, type= "l", col = "red")
lines(traj_dm$traj[,1], Is010, type= "l", col = "green")
lines(traj_dm$traj[,1], Is110, type= "l", col = "black")
lines(traj_dm$traj[,1], Is001, type= "l", col = "purple")
lines(traj_dm$traj[,1], Is101, type= "l", col = "grey")
lines(traj_dm$traj[,1], Is011, type= "l", col = "orange")
lines(traj_dm$traj[,1], Is111, type= "l", col = "yellow")


ncol(traj_dm$traj)


Susceptible<-as.vector(traj_dm$traj[,4])
EarlyLatent<-as.vector(rowSums(traj_dm$traj[,5:12]))
LateLatent<-as.vector(rowSums(traj_dm$traj[,13:20]))
Infected<-as.vector(rowSums(traj_dm$traj[,21:44]))
Recovered<-as.vector(rowSums(traj_dm$traj[,44:68]))
Infected_Resistant<-as.vector(rowSums(traj_dm$traj[,c(22,23,24,25,26,27,28,30,31,32,33,34,35,36,38,39,40,41,42,43,44)]))

layout(matrix(1:1,1,1))
plot(traj_dm$traj[,1], log(Susceptible),ylim=c(0,12), type= "l", col = "black", ylab="log number of individuals", xlab="time (years)")
lines(traj_dm$traj[,1], log(EarlyLatent), type= "l", col = "red")
lines(traj_dm$traj[,1], log(LateLatent), type= "l", col = "purple")
lines(traj_dm$traj[,1], log(Infected), type= "l", col = "blue")
lines(traj_dm$traj[,1], log(Recovered), type= "l", col = "green")



layout(matrix(1:4,2,2))
plot(traj_dm$traj[,1],Susceptible, type = "l", xlab="time (years)")
plot(traj_dm$traj[,1], Recovered, type= "l", col = "green",xlab="time (years)")
plot(traj_dm$traj[,1], Infected, type= "l", col = "blue", xlab="time (years)")
plot(traj_dm$traj[,1], Infected_Resistant, type="l", col= "darkblue", xlab="time (years)")
