source("HMM functions from appendix.R")
source("normal modified functions appendix.R")

########## Generate a poisson distributed sample from an HMM ########
T=1000
pois3s <- list(m = 3,
               lambda = c(12, 19, 29),
               gamma = matrix(c(1/3, 1/3, 1/3, 2/3, 0, 1/3, 1/2, 1/2, 0), nrow=3, ncol=3, byrow = TRUE),
               delta = 1/32 * c(15, 9, 8))

sample <- pois.HMM.generate_sample(T, pois3s)
pois3sdata = data.frame(Time = c(1:T), Observation = sample$observ, State = factor(sample$state))

########## Fit 3-state HMM to generated data ########
summary(pois3sdata$Observation) 
m       <-3 
lambda0 <-c(12,17,23)
gamma0  <-matrix(c(0.6 ,0.2 ,0.2,
                   0.2 ,0.6 ,0.2,
                   0.2 ,0.2 ,0.6),m,m,byrow=TRUE) 

modpois3s   <-pois.HMM.mle(pois3sdata$Observation, m, lambda0, gamma0, stationary=TRUE)

delta0  <- c(1,1,1)/3 
modpois3h   <-pois.HMM.mle(pois3sdata$Observation, m, lambda0, gamma0, delta=delta0, stationary=FALSE) 

modpois3s; modpois3h
#modpois3s$mllk > modpois3h$mllk so we will go with the stationary one

abs(pois3s$lambda - modpois3s$lambda)
abs(pois3s$gamma - modpois3s$gamma)


#Find underlying state sequence 
state_seq3 <- pois.HMM.viterbi(pois3sdata$Observation, modpois3s)
percent_accuracy = sum(pois3sdata$State==state_seq3)/T

sate_seq3sm <- state_seq3[1:100]

pois3sdatasm <-head(pois3sdata,100)
pois3sdatasm$GuessState <- factor(sate_seq3sm)


#plot first 100 points
p <- ggplot(pois3sdatasm, aes(x = Time, y = Observation)) +
  ggtitle('Timeseries of Observations and Underlying States') +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = 'Time', y = 'Observation') +
  geom_line() +
  geom_hline(yintercept = pois3s$lambda, alpha = 0.3) + 
  geom_point(aes(x = Time, y = pois3s$lambda[State], colour = State)) +
  labs(colour = "State")
#overlay decoded state sequence
p + geom_hline(yintercept = modpois3s$lambda, alpha = 0.4, colour = 'purple') + 
  geom_point(aes(Time, y = modpois3s$lambda[GuessState], colour = State)) +
  labs(colour = "State")


#plot last 100 points
sate_seq3tail <- state_seq3[(T-99):T]

pois3sdatatail <-tail(pois3sdata,100)
pois3sdatatail$GuessState <- factor(sate_seq3tail)

p2 <- ggplot(pois3sdatatail, aes(x = Time, y = Observation)) +
  ggtitle('Timeseries of Observations and Underlying States') +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = 'Time', y = 'Observation') +
  geom_line() +
  geom_hline(yintercept = pois3s$lambda, alpha = 0.3) + 
  geom_point(aes(x = Time, y = pois3s$lambda[State], colour = State)) +
  labs(colour = "State")
#overlay decoded state sequence
p2 + geom_hline(yintercept = modpois3s$lambda, alpha = 0.4, colour = 'purple') + 
  geom_point(aes(Time, y = modpois3s$lambda[GuessState], colour = State)) +
  labs(colour = "State")


############ Try fitting 2-state HMM to generated data ########
summary(pois3sdata$Observation) 
m       <-2
lambda0 <-c(13, 23)
gamma0  <-matrix(c(0.8, 0.2,
                   0.2, 0.8),m,m,byrow=TRUE) 

modpois2s <-pois.HMM.mle(pois3sdata$Observation, m, lambda0, gamma0, stationary=TRUE)

#for T=1000, modpois2s$mllk > modpois3s$mllk and modpois2s$AIC < modpois3s$AIC???
#for T=10000, modpois2s$mllk > modpois3s$mllk and modpois2s$AIC > modpois3s$AIC
#for T=100000, modpois2s$mllk > modpois3s$mllk and modpois2s$AIC > modpois3s$AIC



############ Generate a normal distribution HMM sample ##############
T=1000
norm3s <- list(m = 3,
               mu = c(-3, 10, 18),
               sigma = c(1, 2, 1.5),
               gamma = matrix(c(0.955, 0.024, 0.021,
                                0.050, 0.899, 0.051,
                                0.000, 0.197, 0.803), nrow=3, ncol=3, byrow = TRUE),
               delta = c(0.4436, 0.4045, 0.1519))

sample <- norm.HMM.generate_sample(T, norm3s)
normdata = data.frame(Time = c(1:T), Observation = sample$observ, State = factor(sample$state))

########## Fit 3-state HMM to generated data ########
summary(normdata$Observation) 
m       <- 3 
mu0     <- c(-3, 9, 11)
sigma0  <- c(1, 2, 1.5) #figure out how to make inital guess for sd
gamma0  <- matrix(c(0.9 ,0.05 ,0.05,
                    0.05 ,0.9 ,0.05,
                    0.05 ,0.05 ,0.9),m,m,byrow=TRUE) 

modnorm3s   <-norm.HMM.mle(normdata$Observation, m, mu0, sigma0, gamma0, stationary=TRUE)

delta0  <- c(1,1,1)/3 
modnorm3h   <-norm.HMM.mle(normdata$Observation, m, mu0, sigma0, gamma0, delta=delta0, stationary=FALSE) 

