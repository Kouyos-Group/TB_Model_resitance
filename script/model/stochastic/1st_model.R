###Load packages
library(TiPS)

###Write Reactions
reactions <- c(
  #birth
  paste0("0 [ pi * ", c("S","L_s","L_r","I_s","I_r","R"), " ] -> S"),
  
  #infection
  paste0("S [ alpha_",c("s","r")," * S * I_", c("s","r")," ] -> L_", c("s","r")),
  
  #disease progression
  paste0("L_", c("s","r"), " [ delta * ", "L_", c("s","r"), " ] -> I_", c("s","r")),
  
  #successful treatment
  paste0("I_", c("s","r"), " [ tau_", c("s", "r"), " * I_", c("s","r"), " ] -> R"),
  
  #treatment fail, strains become resistant
  "I_s [ I_s * rho ] -> I_r",
  
  #death due to TB
  paste0("I_", c("s","r"), " [ theta_", c("s", "r"), " * I_", c("s","r"), " ] -> S"),
  
  #natural death
  paste0(c("S","L_s","L_r","I_s","I_r","R"), " [  mu * ", c("S","L_s","L_r","I_s","I_r","R"), " ] -> S")
)



###Build the simulator that will allow to run multiple trajectories
sir_simu <- build_simulator(reactions)

###Initial States
initialStates <- c(S = 7700, L_s = 1500, L_r = 0, I_s = 1225, I_r= 0, R = 0)

###Duration of the simulation
time <- c(0, 200)
dT <- round(1/365, 4) # use a daily time step

###Parameters value
theta <- list(pi = 0, #birth
              alpha_s = 2*10^-5, alpha_r= 2*10^-5, #infection 
              delta= 0.116, #progression
              tau_s= 0.2873, tau_r = 0, #cured 
              rho=0.2, #become resistant
              theta_s= 0.123, theta_r = 0.123, #death due to tuberculosis
              mu= 1/70) #natural death

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
  method = "approximate")

###Visualisation
# If trajectory simulated
#count of each compartments 
Susceptible = traj_dm$traj[,4]
Latent_sensible = traj_dm$traj[,5]
Latent_resistant = traj_dm$traj[,6]
Infected_sensible = traj_dm$traj[,7]
Infected_resistant = traj_dm$traj[,8]
Recovered = traj_dm$traj[,9]

Times = traj_dm$traj[,1]

par(mar=c(5, 5, 3, 13),xpd=TRUE)

plot(Times, Susceptible, type='l', col='green', ylim=c(0,10000) ,lty=1, lwd = 2,xlab="Time (in year)", 
     ylab="Number of Indivuals",cex.lab=2.5,cex.axis=1.5 )
lines(Times, Latent_sensible, col="orange", lty=1, lwd=2)
lines(Times, Latent_resistant, col="yellow", lty=1, lwd=2)
lines(Times, Infected_sensible, col="red", lty=1, lwd=2)
lines(Times, Infected_resistant, col="purple", lty=1, lwd=2)
lines(Times, Recovered, col="blue", lty=1, lwd=2)

legend(x='bottomright', legend=c("Susceptible", "Latent sensible", "Latent resistant", "Infected sensible", "Infected resistant", "Recovered"), 
       fill = c("green", 'orange', 'yellow', 'red', 'purple', 'blue'),cex=1.2, bg="transparent", 
       box.lty=0,inset=c(-0.5,0))

tail(traj_dm$traj)




