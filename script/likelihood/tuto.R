# Model by Christian Althaus
# Remove current working environment
rm(list = ls())
# Necessary libraries
if(!requireNamespace("pacman", quietly = T))
  install.packages("pacman")
pacman::p_load("deSolve", "chron", "tidyverse")
library(here)


# Read the data (and set the working directory using setwd() )
ebola <- read.csv("Ebola_outbreak_West_Africa_data.csv")
str(ebola)

#Column date: format is problematic --> convert to date format
lct <- Sys.getlocale("LC_TIME"); Sys.setlocale("LC_TIME", "C")
ebola$Date<-as.Date(ebola$Date, "%d %b %Y")

# plot the cumulative number of cases versus time for each of the three countries
# add the number of deaths (using "points")
# Guinea:
plot(ebola$Date, ebola$Guinea_Cases)
points(ebola$Date, ebola$Guinea_Death, pch='x')

# Sierra Leone:
plot(ebola$Date, ebola$SierraLeone_Cases)
points(ebola$Date, ebola$SierraLeone_Death, pch='x')

# Liberia:
plot(ebola$Date[is.na(ebola$Liberia_Cases)==F], ebola$Liberia_Cases[is.na(ebola$Liberia_Cases)==F])
points(ebola$Date[is.na(ebola$Liberia_Cases)==F], ebola$SierraLeone_Death[is.na(ebola$Liberia_Cases)==F])

# Estimating R0 using linear regression
# Restrict to the first month
Guinea_first_month <- ebola[which(ebola$Date < min(ebola$Date)+30),
                            c("Date","Guinea_Cases", "Date") ]
plot(x=Guinea_first_month$Date, Guinea_first_month$Guinea_Cases, col="blue")

# Add a new column with a log transformed number of cases
Guinea_first_month<-cbind(Guinea_first_month,
                          log_GC=log(Guinea_first_month$Guinea_Cases))
# Outbreak day
Guinea_first_month$outbreak_day<-as.numeric(Guinea_first_month$Date)-
  min(as.numeric(Guinea_first_month$Date))+1
linear_model<- lm(formula = log_GC~outbreak_day, Guinea_first_month)
sum_linear_model<- summary(linear_model)
plot(x=Guinea_first_month$outbreak_day, y=Guinea_first_month$log_GC,
     col="blue", main="Association between calendar time and logCases")
abline(linear_model, col="red", lwd=3)

# Next, calculate the R0 using: R0 = (1 + delta*incubation period (days)) * (1 + delta * infectious peri
# You can use incubation and infectious periods that are stated in the paper. R0 Guinea = xxx
# Compare your estimate to Table 2.
delta<-as.numeric(linear_model$coefficients["outbreak_day"])
# apply to formula:
# R0
# incubation period 5.3
# infectious period 5.61
R_null<-((1+delta*5.3)*(1+delta*5.61))
# Get confidence intervals
# Note: basic formula for CI of means: 95%CI = mean estimate +/- 1.96*SD
UL_delta <- delta+1.96*sum_linear_model$coefficients["outbreak_day","Std. Error"]
LL_delta <- delta-1.96*sum_linear_model$coefficients["outbreak_day","Std. Error"]
# Upper limit R0
R_null_UL <- ((1+UL_delta*5.3)*(1+UL_delta*5.61))
# Lower limit R0
R_null_LL <- ((1+LL_delta*5.3)*(1+LL_delta*5.61))

(Guinea_R0_CI <- c(R0 = R_null,LL = R_null_LL,
                   UL = R_null_UL))

# Definition of the SEIR model
SEIR <- function(t, x, parms) {
  with(as.list(c(parms,x)),{
    if(t < tau1) beta <- beta0
    else beta <- beta0*exp(-k*(t-tau1))
    # beta0 = per contact infection rate
    #k = decay induced by the control measures
    # tau1 = time of introduction of control measures
    N <- S + E + I + R
    # N = Total, S = Susceptible, E = Exposed, I = infected, R = Recovered
    dS <- - beta*S*I/N
    dE <- beta*S*I/N - sigma*E
    dI <- sigma*E - gamma*I
    # 1/sigma = average duration of incubation
    # 1/gamma = average duration of infectiousness
    dR <- (1-f)*gamma*I
    #f = fatality (death) rate
    dD <- f*gamma*I
    #D = deaths
    dC <- sigma*E
    #C = Cumulative cases of infected
    der <- c(dS,dE,dI,dR,dD,dC)
    list(der)
  })
}

# Cost function to calculate the sum of squared residuals (SSR) of predicted vs observed data
cost <- function(free,fixed,init,data) {
  pars <- c(free,fixed)
  pars <- trans(pars)
  times <- c(0,data$times+pars["tau0"])
  simulation <- as.data.frame(ode(init,times,SEIR,parms=pars))
  simulation <- simulation[-1,]
  cost <- sum((simulation$C-data$cases)^2)
  return(cost)
}

trans <- function(pars) {
  pars["beta0"] <- exp(pars["beta0"])
  pars["k"] <- exp(pars["k"])
  pars["f"] <- plogis(pars["f"])
  pars["tau0"] <- exp(pars["tau0"])
  pars["tau1"] <- exp(pars["tau1"])
  return(pars)
}
# Fit the model to the data
#####################################################
# GUINEA: Prepare the data and set the initial values
#####################################################
data <- na.omit(ebola[c("Date","Guinea_Cases","Guinea_Death")])
names(data) <- c("times","cases","deaths")
# data$times <- chron(as.character(data$times), format=c(dates = "day mon year"))
data$times <- as.Date(as.character(data$times), format="%Y-%m-%d")
begin <- chron("2 Dec 2013", format=c(dates = "day mon year"))
delay <- as.numeric(data$times[1] - as.Date(begin))
## Warning: Incompatible methods ("-.Date", "Ops.dates") for "-"
data$times <- data$times - data$times[1]
N <- 1e6
init <- c(S = N - 1, E = 0, I = 1, R = 0, D = 0, C = 1)
fixed <- c(tau0 = log(delay), tau1 = -Inf, sigma = 1/5.3, gamma = 1/5.61)
free <- c(beta0 = log(0.2), k = log(0.001), f = 0)
# Fit the model to the data
fit <- optim(free,cost,gr=NULL,fixed,init,data,method="Nelder-Mead",hessian=TRUE)
# create table that summarizes the results:
cols <- c("Estimate","2.5%","97.5%")
rows <- c("R0",names(fit$par))
est <- as.data.frame(matrix(NA,nrow=length(rows),ncol=length(cols)))
rownames(est) <- rows
colnames(est) <- cols
# hint: R0=beta0/gamma
# CAVE model estimates are on log scale, so we have to exp them (and round)
# start with first column:
est[1] <- round(c(exp(fit$par[1])/fixed["gamma"],exp(fit$par[1]),
                  exp(fit$par[2]),plogis(fit$par[3])),4)
# define 95%CI
up <- fit$par + 1.96*sqrt(diag(solve(fit$hessian)))
lo <- fit$par - 1.96*sqrt(diag(solve(fit$hessian)))
# second and third columns:
est[2] <- round(c(exp(lo[1])/fixed["gamma"],exp(lo[1]),exp(lo[2]),plogis(lo[3])),4)
est[3] <- round(c(exp(up[1])/fixed["gamma"],exp(up[1]),exp(up[2]),plogis(up[3])),4)
est

# Save results for individual country
fitG <- fit
estG <- est

# Plot the best-fit model to the data
pars <- trans(c(fit$par,fixed))
end <- chron("1 Sep 2014", format=c(dates = "day mon year"))
months <- chron(c("1 Dec 2013","1 Mar 2014","1 Jun 2014","1 Sep 2014"),
                format=c(dates = "day mon year"))
timepoints <- seq(0,end-begin,1)
simulation <- as.data.frame(ode(init,timepoints,SEIR,parms=pars))
# calculate incident cases using diff: remove last timepoint entry
# If only diff(simulation$C) is entered, then x and y not the same length
# Hence: remove last timepoint entry
plot(timepoints[1:(length(timepoints)-1)], diff(simulation$C), type="l", col="red", ylim=c(0,10), xlab=N,
     ylab="Incidence number of cases/deaths",frame=FALSE,axes=FALSE,main="Guinea")
axis(1,months-begin,months)
axis(2)
lines(timepoints[1:(length(timepoints)-1)], diff(simulation$D), lty=2)
points(data$times+pars["tau0"],data$cases,col="red")
points(data$times+pars["tau0"],data$deaths,pch=0)
legend("topleft",inset=0.05,legend=c("New Cases","New Deaths"), col=c("red","black"), lty=c(1,2))

