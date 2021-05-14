source("HMM_functions_from_appendix.R")
#Get earthquake data
dat <- read.table( "http://www.hmms-for-time-series.de/second/data/earthquakes.txt")

x <-dat[,2] #observation of number of earthquakes per year
d <-dat[,1] #years
n <-length(x) #number of observations


######### Fit 2-state HMM ########
summary(x) #this can tell us a good lamda0... mean=19.36 so 20+/-5 are good guesses for lambdas
m       <- 2 #number of states
lambda0 <- c(15,25) #inital guess for means of poisson distributions
gamma0  <- matrix(c(0.9,0.1,0.1,0.9),m,m,byrow=TRUE) #initial guess for transition probability matrix (off-diagonals equal)

mod2s   <- pois.HMM.mle(x,m,lambda0,gamma0,stationary=TRUE) #find HMM that fits data (assume stationary distribution) 

delta0  <- c(1,1)/2 #use initial distribution with equal probabilities for each state
mod2h   <- pois.HMM.mle(x,m,lambda0,gamma0,delta=delta0,stationary=FALSE) #find HMM that fits data (non-stationary) 

mod2s; mod2h
#mod2s$mllk = 342.3183 ; mod2h$mllk = 341.8787 ; mod2s$mllk > mod2h$mllk
#mod2s$AIC = 692.6365 ; mod2h$AIC = 693.7574 ; mod2s$AIC < mod2h$AIC
#mod2s$BIC = 703.3278 ; mod2s$BIC = 707.1216 ; mod2s$BIC < mod2h$BIC
#the results suggest that mod2s is a better fit than mod2h


######### Fit 3-state HMM ########
summary(x) #mean~20 so guess lambda0 from here
m       <-3 
lambda0 <-c(10,20,30)
gamma0  <-matrix(c(0.8 ,0.1 ,0.1,
                   0.1 ,0.8 ,0.1,
                   0.1 ,0.1 ,0.8),m,m,byrow=TRUE) 

mod3s   <-pois.HMM.mle(x,m,lambda0,gamma0,stationary=TRUE)

delta0  <- c(1,1,1)/3 
mod3h   <-pois.HMM.mle(x,m,lambda0,gamma0,delta=delta0,stationary=FALSE) 

mod3s; mod3h
#mod3s$mllk = 329.4603 ; mod3h$mllk = 328.5275 ; mod3s$mllk > mod3h$mllk
#mod3s$AIC = 676.9206 ; mod3h$AIC = 679.055 ; mod3s$AIC < mod3h$AIC
#mod3s$BIC = 700.976 ; mod3h$BIC = 708.4561 ; mod3s$BIC < mod3h$BIC
#the results suggest that mod3s is a better fit than mod3h


######### Fit 4-state HMM ########
m       <- 4
lambda0 <- c(10,15,20,30)
gamma0  <- matrix(c(0.85, 0.05 ,0.05 ,0.05, 
                    0.05 ,0.85 ,0.05 ,0.05, 
                    0.05 ,0.05 ,0.85 ,0.05,  
                    0.05 ,0.05 ,0.05 ,0.85), m, m, byrow=TRUE)

mod4s   <- pois.HMM.mle(x,m,lambda0,gamma0,stationary=TRUE)

delta0  <- c(1,1,1,1)/4
mod4h   <- pois.HMM.mle(x,m,lambda0,gamma0,delta=delta0,stationary=FALSE)

mod4s; mod4h
#mod4s$mllk = 327.8316 ; mod4h$mllk = 326.6749 ; mod4s$mllk > mod4h$mllk
#mod4s$AIC = 687.6632 ; mod4h$AIC = 691.3499 ; mod4s$AIC < mod4h$AIC
#mod4s$BIC = 730.4284 ; mod4h$BIC = 742.1336 ; mod4s$BIC < mod4h$BIC
#the results suggest that mod4s is a better fit than mod4h


######### Global Decoding using 3-state HMM ########
state_seq <- pois.HMM.viterbi(x, mod3s)
Earthquake3 <- data.frame(Time = d,
                          Observation = x,
                          State = factor(state_seq))

#function to plot observation sequence and underlying state sequence
pois.HMM.plot_seq <- function(dat, lambda, m,
                              title = 'Timeseries of Observations and Underlying States',
                              xlabel = 'Time',
                              ylabel = 'Observation') {
  p <- ggplot(dat, aes(x = Time, y = Observation)) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = xlabel, y = ylabel) +
    geom_line()+
    geom_hline(yintercept = lambda, alpha = 0.3) + 
    geom_point(aes(x = Time, y = lambda[State], colour = State)) +
    labs(colour = "State")
  return(p)
}

pois.HMM.plot_seq(Earthquake3, mod3s$lambda, 3, 
                  title = 'Number of Earthquakes Per Year and the Underlying State Sequence',
                  xlabel = 'Year',
                  ylabel = 'Number of Earthquakes')
