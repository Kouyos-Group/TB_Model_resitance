#LOAD THE PACKAGES
rm(list = ls())
# Necessary libraries
if(!requireNamespace("pacman", quietly = T))
  install.packages("pacman")
pacman::p_load("deSolve", "chron", "tidyverse")
library(here)

#LOAD THE REAL DATA

data<-data.frame(cbind(times=c(0:12), infected=c(1.00,1.35, 1.79, 2.33, 3.04, 3.98, 5.19, 6.77, 8.37, 11.53,15.03, 19.58, 25.51)))

#LOAD THE MODEL
source(here("./model/deterministic/generalized_model.R"))
source(here("./model/deterministic/1st_scenario_before_70.R"))

#####
# Cost function to calculate the sum of squared residuals (SSR) of predicted vs observed data
cost <- function(free,fixed,init,data) {
  pars <- c(free,fixed) 
  pars <- trans(pars)
  times <- c(0,data$times)
  simulation <- as.data.frame(ode(init,times,sir_equations,parms=pars))
  simulation <- simulation[-1,]
  cost <- sum((rowSums(simulation[,35:50])-data$infected)^2)
  return(cost)
}

# Parameter transformation
trans <- function(pars) {
  pars["alphas0000"] <- exp(pars["alphas0000"])
  pars["omicron"] <- exp(pars["omicron"])
  return(pars)
}

# Fit the model to the data
init<-state
data$times <- data$times - data$times[1]

#what we want to test
free <- c(alphas0000 = log(infection), omicron= log(0.2/12))
#remove to fixed
fixed<-parameters_values[ ! names(parameters_values) %in%  names(free) ]


# Fit the model to the data
fit <- optim(free,cost,gr=NULL,fixed,init,data,method="Nelder-Mead",hessian=TRUE)



# create table that summarizes the results:
cols <- c("Estimate","2.5%","97.5%")
rows <- c(names(fit$par))
est <- as.data.frame(matrix(NA,nrow=length(rows),ncol=length(cols)))
rownames(est) <- rows
colnames(est) <- cols

# hint: R0=beta0/gamma
# CAVE model estimates are on log scale, so we have to exp them (and round)
# start with first column:
est[1] <- round(c(exp(fit$par[1]),exp(fit$par[2])),4)
# define 95%CI
up <- fit$par + 1.96*sqrt(diag(solve(fit$hessian)))
lo <- fit$par - 1.96*sqrt(diag(solve(fit$hessian)))
# second and third columns:
est[2] <- round(c(exp(lo[1]),exp(lo[2])),5)
est[3] <- round(c(exp(up[1]),exp(up[2])),5)
est

# Save results for individual country
fitG <- fit
estG <- est

# Plot the best-fit model to the data
pars <- trans(c(fit$par,fixed))
timepoints <- seq(0,12,1)
simulation <- as.data.frame(ode(init,timepoints,sir_equations,parms=pars))

# calculate incident cases using diff: remove last timepoint entry
# If only diff(simulation$C) is entered, then x and y not the same length
# Hence: remove last timepoint entry
plot(timepoints[1:(length(timepoints)-1)], diff(rowSums(simulation[,35:50])), type="l", col="red", ylim=c(0,10), xlab="times",
     ylab="number of infected",frame=FALSE,axes=FALSE,main="Guinea")
axis(1,timepoints)
axis(2)
points(data$times,data$infected,col="red")
legend("topleft",inset=0.05,legend=c("New Cases"), col=c("red"), lty=c(1,2))

