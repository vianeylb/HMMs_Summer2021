library(ggplot2)

#function to remove +0i from real numbers
rid_Im <- function(x) {
  if (all(Im(z <- zapsmall(x))==0)) as.numeric(z) else x
}

#function to find stationary distribution (initial state distribution)
compute_delta <- function(A){
  e = eigen(t(A)) #find eigenvalues and left eigenvectors of A
  vec = rid_Im(e$vectors[,rid_Im(e$values)==1.0])/sum(rid_Im(e$vectors[,rid_Im(e$values)==1.0])) #eigenvector with eigenvalue 1 with components summing to 1
  return(vec)
}


#function to simulate data from an HMM with NORMAL state dependent distributions
sim_HMM_norm <- function(T, N, trans.prob, delta, state.dep.norm) {
  output = data.frame("Index" = c(1:T),
                      "State" = NA,
                      "Observation" = NA)
  states = c(1:N)
  
  #get inital state from initial state disribution delta
  output$State[1] = sample(x=states, size=1, prob=delta)
  #get initial observation using inital state found and the state-dependent normal distributions state.dep.norm
  output$Observation[1] = rnorm(n=1, mean=state.dep.norm[output$State[1],1], sd=state.dep.norm[output$State[1],2])
  
  for (i in 2:T) {
    #using the previous state and the transition probability, get the following state and observation
    transition = trans.prob[output$State[i-1],]
    output$State[i] = sample(x=states, size=1, prob=transition)
    output$Observation[i] = rnorm(n=1, mean=state.dep.norm[output$State[i],1], sd=state.dep.norm[output$State[i],2])
  }

  output$State <- as.factor(output$State)
  return(output)
}


#Example inputs for sim_HMM_norm of an HMM with two states
#(note that state.dep.norm must have the format column1=mean, column2=sd and each row is a different state in order)
A = matrix(c(0.1, 0.9, 0.4, 0.6), nrow=2, ncol=2, byrow = TRUE) #transition probablity matrix example
d = compute_delta(A) #inital state distribution example
normdist = data.frame("mean" = c(4, 10), "sd" = c(1, 2)) #state dependent normal distributions example

#output of sim_HMM_norm with example HMM
norm.hmm = sim_HMM_norm(T=10000, N=2, trans.prob = A, delta = d, state.dep.norm = normdist)

#function to plot a histogram of the simulated data from sim_HMM_norm overlayed with the state dependent normal distributions
hist_results_norm <- function(sim.data, 
                              state.dep.norm, 
                              title = 'Histogram of HMM Simulation',
                              xlabel = 'Observation', 
                              ylabel = 'Frequency', 
                              numbreaks=15){
  observ = sim.data$Observation
  h <- hist(sim.data$Observation, main=title, xlab=xlabel, ylab=ylabel, breaks=numbreaks)
  xfit <- seq(min(observ), max(observ), length=150)
  #more colours and names can be added if number of states exceeds 9
  colours = c('blue', 'red', 'green', 'darkgreen', 'gold', 'deeppink', 'purple', 'black', 'brown')
  name = c('State 1', 'State 2', 'State 3', 'State 4', 'State 5', 'State 6', 'State 7', 'State 8', 'State 9')
  
  for (i in 1:nlevels(sim.data$State)) {
    yfit <- dnorm(xfit, mean = state.dep.norm$mean[i], sd = state.dep.norm$sd[i]) 
    yfit <- yfit * diff(h$mids[1:2]) * sum(sim.data$State==i) 
    lines(xfit, yfit, lwd = 2, col=colours[i])
  }
  legend("topright", name[1:nlevels(sim.data$State)], col=colours[1:nlevels(sim.data$State)], lwd=2)
}


#plot histogram for the simulated data from our example HMM
hist_results_norm(norm.hmm, normdist)

#function to plot the time series of the simulated data from sim_HMM_norm
timeseries <- function(sim.data,
                       title = 'Time Series of Simulated HMM',
                       xlabel = 'Time',
                       ylabel = 'Observation') {
  ggplot(sim.data, aes(x = Index, y = Observation)) +
    theme_light() +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_point(aes(color = State))+
    geom_line(colour='grey', alpha=0.8, lwd=0.3) +
    labs(x = xlabel, y = ylabel)
}

#plot the time series for simulated HMM data example
small.norm.hmm = sim_HMM_norm(T=500, N=2, trans.prob = A, delta = d, state.dep.norm = normdist)
timeseries(small.norm.hmm)


####### sim_HMM for state dependent distributions of type normal, gamma, beta, poisson, uniform, logistic #######

sim_HMM <- function(T, N, trans.prob, delta, state.dist) {
  output = data.frame("Index" = c(1:T),
                      "State" = NA,
                      "Observation" = NA)
  states = c(1:N)
  
  output$State[1] = sample(x=states, size=1, prob=delta)
  
  if (state.dist$type[output$State[1]] == 'normal') {
    output$Observation[1] = rnorm(n=1,
                                  mean=state.dist$norm.mean[output$State[1]], 
                                  sd=state.dist$norm.sd[output$State[1]])
    
  } else if (state.dist$type[output$State[1]] == 'gamma') {
    output$Observation[1] = rgamma(n=1,
                                   shape=state.dist$gamma.shape[output$State[1]],
                                   rate=state.dist$gamma.rate[output$State[1]])
    
  } else if (state.dist$type[output$State[1]] == 'beta') {
    output$Observation[1] = rbeta(n=1,
                                  shape1=state.dist$beta.shape1[output$State[1]],
                                  shape2=state.dist$beta.shape2[output$State[1]])
    
  } else if (state.dist$type[output$State[1]] == 'poisson') {
    output$Observation[1] = rpois(n=1,
                                  lambda=state.dist$pois.lambda[output$State[1]])
    
  } else if (state.dist$type[output$State[1]] == 'uniform') {
    output$Observation[1] = runif(n=1,
                                  min=state.dist$unif.min[output$State[1]],
                                  max=state.dist$unif.max[output$State[1]])
    
  } else if (state.dist$type[output$State[1]] == 'logistic') {
    output$Observation[1] = rlogis(n=1,
                                  location=state.dist$log.location[output$State[1]],
                                  scale=state.dist$log.scale[output$State[1]])
  }
  
  for (i in 2:T) {
    transition = trans.prob[output$State[i-1],]
    output$State[i] = sample(x=states, size=1, prob=transition)
    if (state.dist$type[output$State[i]] == 'normal') {
      output$Observation[i] = rnorm(n=1,
                                    mean=state.dist$norm.mean[output$State[i]], 
                                    sd=state.dist$norm.sd[output$State[i]])
      
    } else if (state.dist$type[output$State[i]] == 'gamma') {
      output$Observation[i] = rgamma(n=1,
                                     shape=state.dist$gamma.shape[output$State[i]],
                                     rate=state.dist$gamma.rate[output$State[i]])
      
    } else if (state.dist$type[output$State[i]] == 'beta') {
      output$Observation[i] = rbeta(n=1,
                                    shape1=state.dist$beta.shape1[output$State[i]],
                                    shape2=state.dist$beta.shape2[output$State[i]])
      
    } else if (state.dist$type[output$State[i]] == 'poisson') {
      output$Observation[i] = rpois(n=1,
                                    lambda=state.dist$pois.lambda[output$State[i]])
      
    } else if (state.dist$type[output$State[i]] == 'uniform') {
      output$Observation[i] = runif(n=1,
                                    min=state.dist$unif.min[output$State[i]],
                                    max=state.dist$unif.max[output$State[i]])
      
    } else if (state.dist$type[output$State[i]] == 'logistic') {
      output$Observation[i] = rlogis(n=1,
                                     location=state.dist$log.location[output$State[i]],
                                     scale=state.dist$log.scale[output$State[i]])
    }
  }
  output$State <- as.factor(output$State)
  return(output)
}


#Example inputs for sim_HMM for an HMM with three states
B = matrix(c(1/3, 1/3, 1/3, 2/3, 0, 1/3, 1/2, 1/2, 0), nrow=3, ncol=3, byrow = TRUE) #transition probablity matrix example
delt = compute_delta(B) #initial state distribution example
#state dependent distributions (note that each row is a different state distribution and each column is a different 
#parameter for the distributions that sim_HMM allows)
#sim_HMM only takes state.dist written with these column names (you can omit columns that are not needed for desired distributions)
dist <- data.frame("type" = as.character(c('normal', 'gamma', 'poisson')),
                    "norm.mean" = as.numeric(c(10, NA, NA)), 
                    "norm.sd" = as.numeric(c(1, NA, NA)),
                    "gamma.shape" = as.numeric(c(NA, 5, NA)),
                    "gamma.rate" = as.numeric(c(NA, 1, NA)),
                    "beta.shape1" = as.numeric(c(NA, NA, NA)),
                    "beta.shape2" = as.numeric(c(NA, NA, NA)),
                    "pois.lambda" = as.numeric(c(NA, NA, 10)),
                    "unif.min" = as.numeric(c(NA, NA, NA)),
                    "unif.max" = as.numeric(c(NA, NA, NA)),
                    "log.location" = as.numeric(c(NA, NA, NA)),
                    "log.scale" = as.numeric(c(NA, NA, NA)),
                    stringsAsFactors = FALSE)

#run sim_HMM with our example inputs
hmm.example = sim_HMM(1000, 3, B, delt, dist)

#function to plot a histogram of the simulated data from sim_HMM overlayed with the state dependent distributions
hist_results <- function(sim.data, 
                         state.dist, 
                         title = 'Histogram of HMM Simulation',
                         xlabel = 'Observation', 
                         ylabel = 'Frequency', 
                         numbreaks=15){
  observ = sim.data$Observation
  h <- hist(sim.data$Observation, main=title, xlab=xlabel, ylab=ylabel, breaks=numbreaks)
  xfit <- seq(min(round(sim.data$Observation)), max(round(sim.data$Observation)))
  #more colours and names can be added if number of states exceeds 9
  colours = c('blue', 'red', 'green', 'darkgreen', 'gold', 'deeppink', 'purple', 'black', 'brown')
  name = c('State 1', 'State 2', 'State 3', 'State 4', 'State 5', 'State 6', 'State 7', 'State 8', 'State 9')
  
  for (i in 1:nlevels(sim.data$State)) {
    if (state.dist$type[i] == 'normal'){
      yfit <- dnorm(xfit, mean = state.dist$norm.mean[i], sd = state.dist$norm.sd[i]) 
      yfit <- yfit * diff(h$mids[1:2]) * sum(sim.data$State==i) 
      lines(xfit, yfit, lwd = 2, col=colours[i])
      
    } else if (state.dist$type[i] == 'gamma') {
      yfit <- dgamma(xfit, shape=state.dist$gamma.shape[i], rate=state.dist$gamma.rate[i]) 
      yfit <- yfit * diff(h$mids[1:2]) * sum(sim.data$State==i) 
      lines(xfit, yfit, lwd = 2, col=colours[i])
      
    } else if (state.dist$type[i] == 'beta') {
      yfit <- dbeta(xfit, shape1=state.dist$beta.shape1[i], shape2=state.dist$beta.shape2[i]) 
      yfit <- yfit * diff(h$mids[1:2]) * sum(sim.data$State==i) 
      lines(xfit, yfit, lwd = 2, col=colours[i])
      
    } else if (state.dist$type[i] == 'poisson') {
      yfit <- dpois(xfit, lambda=state.dist$pois.lambda[i])
      yfit <- yfit * diff(h$mids[1:2]) * sum(sim.data$State==i) 
      lines(xfit, yfit, lwd = 2, col=colours[i])
      
    } else if (state.dist$type[i] == 'uniform') {
      yfit <- dunif(xfit, min=state.dist$unif.min[i], max=state.dist$unif.max[i])
      yfit <- yfit * diff(h$mids[1:2]) * sum(sim.data$State==i) 
      lines(xfit, yfit, lwd = 2, col=colours[i])
      
    } else if (state.dist$type[i] == 'logistic') {
      yfit <- dlogis(xfit, location=state.dist$log.location[i], scale=state.dist$log.scale[i])
      yfit <- yfit * diff(h$mids[1:2]) * sum(sim.data$State==i) 
      lines(xfit, yfit, lwd = 2, col=colours[i])
    }
  }
  legend("topright", name[1:nlevels(sim.data$State)], col=colours[1:nlevels(sim.data$State)], lwd=2)
}

#plot histogram for example HMM
hist_results(hmm.example, dist, numbreaks = 25)

#plot timeseries for HMM example
timeseries(hmm.example)


