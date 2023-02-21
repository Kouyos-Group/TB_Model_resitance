##### See Graph ####

time<-sir_values[,1]
S<-sir_values[,2]
E<-rowSums(sir_values[,3:18])
L <-rowSums(sir_values[,19:34])
I <-rowSums(sir_values[,35:50])
T <-rowSums(sir_values[,51:114])
R <-rowSums(sir_values[,115:130])

a<-plot(time, rowSums(sir_values[,2:130]))

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

#Verification
#####
for ( i in 1:length(parameters_values)){
  
  assign(paste0(names(parameters_values[i])),parameters_values[i] )
}
