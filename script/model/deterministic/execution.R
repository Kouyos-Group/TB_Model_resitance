library(here)
source(here("./model/deterministic/generalized_model.R"))
#before 70
source(here("./model/deterministic/1st_scenario_before_70.R"))
#  the points in time where to calculate variables values 
####

#  numerically solving the SIR model 
# Run before 70
i=1

assign(paste0("sir_values_", i), ode(
  y = state,
  times = time_values,
  func = sir_equations,
  parms = parameters_values 
)
)

# Run after 70
source(here("./model/deterministic/2nd_scenario_after_70.R"))
i=2

assign(paste0("sir_values_", i), ode(
  y = state,
  times = time_values,
  func = sir_equations,
  parms = parameters_values 
)
)

sir_values <-rbind(sir_values_1, sir_values_2[-1,])



#####

time<-sir_values[,1]
S<-sir_values[,2]
E<-rowSums(sir_values[,3:18])
L <-rowSums(sir_values[,19:34])
I <-rowSums(sir_values[,35:50])
T <-rowSums(sir_values[,51:114])
R <-rowSums(sir_values[,115:130])

a<-plot(time, rowSums(sir_values_1[,2:130]))

# plotting the time seriRs of susceptiblRs:
plot(time, S, type = "l", col = "blue",
     xlab = "time (months)", ylab = "number of people", ylim = c(0,16000))
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

