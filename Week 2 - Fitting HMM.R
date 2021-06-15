library(tidyverse)
library(ggplot2)
library(plyr)

#Transform Poisson natural parameters to working parameters
#Discussed in section 3.3, starting page 50 
#m: number of states
#lambda: vector of parameters for Poisson dist
#gamma: matrix of transition probabilities
#delta: vector of initial distribution
pois.HMM.pn2pw <- function(m,lambda,gamma,delta=NULL,stationary=TRUE)
{
  #Transform lambda by taking log
  tlambda <- log(lambda)
  #Transform gamma by log(gamma_ij/gamma_ii)
  foo <- log(gamma/diag(gamma))
  #Take only off-diagonal elements of resulting matrix
  tgamma <- as.vector(foo[!diag(m)])
  #Transform initial distribution parameter if the MC is not stationary 
  if(stationary) 
    {tdelta <-NULL} 
  else 
    #Take log of each element except the first divided by the first element
    {tdelta<-log(delta[-1]/delta[1])}
  parvect <- c(tlambda,tgamma,tdelta)
  return(parvect)
}

# Transform Poisson working parameters to natural parameters
pois.HMM.pw2pn <- function(m,parvect,stationary=TRUE)
{
  #Transform by taking exponent of tlambda
  lambda <- exp(parvect[1:m])
  #Create mxm identity matrix 
  gamma <- diag(m)
  #Set off-diagonal elements of matrix, starting from 1st col and moving down
  gamma[!gamma] <- exp(parvect[(m+1):(m*m)])
  #Divide each element by the row sum, so the rows will sum to 1
  gamma <- gamma/apply(gamma,1,sum)
  if(stationary) 
    #Obtain stationary dist by solving system of equations
    {delta<-solve(t(diag(m)-gamma+1),rep(1,m))} 
  else
    #Create vector of 1 and the elements of tdelta
    {foo<-c(1,exp(parvect[(m*m+1):(m*m+m-1)]))
    #Divide each element of vector by the sum of the elements
    delta<-foo/sum(foo)}
  return(list(lambda=lambda,gamma=gamma,delta=delta))
}

#Computing minus the log-likelihood from the working parameters
#Algorithm is given on page 49 
#parvect: working parameters
#x: observations
#m: number of states
pois.HMM.mllk <- function(parvect, x, m, stationary=TRUE, ...)
  {
    if(m==1) 
      #If only one state, parvect should only contain log(lambda) (?) 
      #Then log-lik is sum of log(p(x))
      return (-sum(dpois(x, exp(parvect), log=TRUE)))
    n <- length(x)
    pn <- pois.HMM.pw2pn(m, parvect, stationary=stationary)
    #foo is delta*P(x1)
    foo <- pn$delta*dpois(x[1], pn$lambda)
    #sumfoo is w1 = delta*P(x1)*1'
    sumfoo <- sum(foo)
    lscale <- log(sumfoo)
    #foo is phi1 
    foo <- foo/sumfoo
    for (i in 2:n)
      {
        #If data is not missing, define P(xi)
        if (!is.na(x[i]))
          {P<-dpois(x[i], pn$lambda)}
        #Else P is identity matrix
        else
          {P<- rep(1,m)}
        #%*% is matrix multiplication
        #Here foo is v
        foo <- foo%*%pn$gamma*P
        sumfoo <- sum(foo)
        lscale <- lscale + log(sumfoo)
        #Here foo is phi_t 
        foo <- foo/sumfoo
      }
    mllk <- -lscale
    return(mllk)
}

#Computing MLE from natural parameters
#lambda0, gamma0, delta0 are starting values for estimation 
pois.HMM.mle <- function(x, m, lambda0, gamma0, delta0=NULL, stationary=TRUE,...)
  {
    #Get working parameters
    parvect0 <- pois.HMM.pn2pw(m, lambda0, gamma0, delta0, stationary=stationary)
    #nlm is an unconstrained minimizer
    mod <- nlm(pois.HMM.mllk, parvect0, x=x,m=m, stationary=stationary)
    #Transform estimates of working parameters into estimates of natural parameters
    pn <- pois.HMM.pw2pn (m=m, mod$estimate, stationary=stationary)
    #Get the negative log likelihood value
    mllk <- mod$minimum
    #np is the number of estimated parameters
    np <- length(parvect0)
    AIC <- 2*(mllk+np)
    #n is the number of estimated parameters
    n <- sum(!is.na(x))
    BIC <- 2*mllk+np*log(n)
    list(m=m, lambda=pn$lambda, gamma=pn$gamma, delta=pn$delta,
               code=mod$code, mllk=mllk, AIC=AIC, BIC=BIC)
  }

#Generating a sample
#ns: length of sample
#mod: list with elements m, lambda, gamma, delta
pois.HMM.generate_sample <- function(ns, mod)
  {
    #mvect is list of states to sample from
    mvect <- 1:mod$m
    state <- numeric(ns)
    #Get initial state 
    state[1] <- sample(mvect, 1, prob=mod$delta)
    #Get remaining states
    for (i in 2:ns) state[i] <- sample(mvect, 1, prob=mod$gamma[state[i-1],])
    #Generate observations
    x <- rpois(ns, lambda = mod$lambda[state])
    return(data.frame(index=c(1:ns), state=state, obs=x)) 
  }

#Global decoding by the Viterbi algorithm
#Discussed in Section 5.4.2 (page 88)
pois.HMM.viterbi <- function (x, mod)
  {
    n <- length(x)
    #xi is a nxm matrix of zeros 
    xi <- matrix(0, n, mod$m)
    #foo is delta*P(x1)
    foo <- mod$delta*dpois(x[1], mod$lambda)
    #Set first row of xi
    xi[1,] <- foo/sum(foo)
    for (t in 2:n)
      {
        #foo is xi_t, this is equation 5.9
        foo <- apply(xi[t-1,]*mod$gamma, 2, max)*dpois(x[t], mod$lambda)
        xi[t,] <- foo/sum(foo)
    }
    #iv is the maximizing sequence of states 
    iv <- numeric(n)
    iv[n] <- which.max(xi[n,])
    for (t in (n-1):1)
      #This is equation 5.11
      iv[t] <- which.max(mod$gamma[,iv[t+1]]*xi[t,])
    return(iv)
  }

graph_hmm_output <- function(output){
  ggplot(data = output, mapping = aes(x=output$index, y=output$obs, 
                                      color = output$state))+
    geom_line()+
    theme_minimal()+
    scale_colour_continuous(type="viridis")
}

#Testing with simulated Poisson HMM
m <- 3
lambda <- c(2, 5, 8)
gamma <- matrix(c(0.1, 0.2, 0.7,
                  0.3, 0.4, 0.3,
                  0.6, 0.3, 0.1), nrow = m, byrow = TRUE)
delta <- c(0.4, 0.2, 0.4)
mod <- list(m=m, lambda=lambda, gamma=gamma, delta=delta)
#Generate sample
pois.sample <- pois.HMM.generate_sample(1000, mod)
head(pois.sample)
graph_hmm_output(pois.sample)
#Get MLE
lambda0 <- c(1, 4, 7)
gamma0 <- matrix(rep(1/m, m*m), nrow = m, byrow = TRUE)
delta0 <- rep(1/m, m)
pois.mle <- pois.HMM.mle(pois.sample$obs, m, lambda0, gamma0, delta0, stationary=TRUE)
pois.mle
pois.decoding <- pois.HMM.viterbi(pois.sample$obs, mod)
count(pois.decoding-pois.sample$state)

 
#Transform normal natural parameters to working parameters
#mu does not need to be transformed
#Transform sigma by taking log
norm.HMM.pn2pw <- function(m,mu,sigma,gamma,delta=NULL,stationary=TRUE)
{
  tsigma <- log(sigma)
  foo <- log(gamma/diag(gamma))
  tgamma <- as.vector(foo[!diag(m)])
  if(stationary) 
  {tdelta <-NULL} 
  else 
  {tdelta<-log(delta[-1]/delta[1])}
  parvect <- c(mu,tsigma,tgamma,tdelta)
  return(parvect)
}

# Transform normal working parameters to natural parameters
norm.HMM.pw2pn <- function(m,parvect,stationary=TRUE)
{
  #mu is untransformed
  mu <- parvect[1:m]
  #Take exponent of tsigma
  sigma <- exp(parvect[(m+1):(2*m)])
  gamma <- diag(m)
  gamma[!gamma] <- exp(parvect[(2*m+1):(m+m*m)])
  gamma <- gamma/apply(gamma,1,sum)
  if(stationary) 
  {delta<-solve(t(diag(m)-gamma+1),rep(1,m))} 
  else
  {foo<-c(1,exp(parvect[(m+m*m+1):(m*m+2*m-1)]))
  delta<-foo/sum(foo)}
  return(list(mu=mu,sigma=sigma,gamma=gamma,delta=delta))
}

#Computing minus the log-likelihood from the working parameters
norm.HMM.mllk <- function(parvect, x, m, stationary=TRUE, ...)
{
  if(m==1) 
    return (-sum(dnorm(x, mean=parvect[1], sd=exp(parvect[2]), log=TRUE)))
  n <- length(x)
  pn <- norm.HMM.pw2pn(m, parvect, stationary=stationary)
  foo <- pn$delta*dnorm(x[1], pn$mu, pn$sigma)
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo/sumfoo
  for (i in 2:n)
  {
    if (!is.na(x[i]))
    {P<-dnorm(x[i], pn$mu, pn$sigma)}
    else
    {P<- rep(1,m)}
    foo <- foo%*%pn$gamma*P
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo/sumfoo
  }
  mllk <- -lscale
  return(mllk)
}

#Computing MLE from natural parameters
norm.HMM.mle <- function(x, m, mu0, sigma0, gamma0, delta0=NULL, stationary=TRUE,...)
{
  parvect0 <- norm.HMM.pn2pw(m, mu0, sigma0, gamma0, delta0, stationary=stationary)
  mod <- nlm(norm.HMM.mllk, parvect0, x=x,m=m, stationary=stationary)
  pn <- norm.HMM.pw2pn (m=m, mod$estimate, stationary=stationary)
  mllk <- mod$minimum
  np <- length(parvect0)
  AIC <- 2*(mllk+np)
  n <- sum(!is.na(x))
  BIC <- 2*mllk+np*log(n)
  list(m=m, mu=pn$mu, sigma=pn$sigma, gamma=pn$gamma, delta=pn$delta,
       code=mod$code, mllk=mllk, AIC=AIC, BIC=BIC)
}

#Generating a sample from normal distributions
norm.HMM.generate_sample <- function(ns, mod)
{
  mvect <- 1:mod$m
  state <- numeric(ns)
  state[1] <- sample(mvect, 1, prob=mod$delta)
  for (i in 2:ns) state[i] <- sample(mvect, 1, prob=mod$gamma[state[i-1],])
  x <- rnorm(ns, mean = mod$mu[state], sd = mod$sigma[state])
  return(data.frame(index=c(1:ns), state=state, obs=x)) 
}

#Global decoding by the Viterbi algorithm
pois.HMM.viterbi <- function (x, mod)
{
  n <- length(x)
  xi <- matrix(0, n, mod$m)
  foo <- mod$delta*dnorm(x[1], pn$mu, pn$sigma)
  xi[1,] <- foo/sum(foo)
  for (t in 2:n)
  {
    foo <- apply(xi[t-1,]*mod$gamma, 2, max)*dnorm(x[t], pn$mu, pn$sigma)
    xi[t,] <- foo/sum(foo)
  }
  iv <- numeric(n)
  iv[n] <- which.max(xi[n,])
  for (t in (n-1):1)
    iv[t] <- which.max(mod$gamma[,iv[t+1]]*xi[t,])
  return(iv)
}

#Testing functions for normal distributions
m <- 3
mu <- c(2, 5, 8)
sigma <- c(2, 4, 6)
gamma <- matrix(c(0.1, 0.2, 0.7,
                  0.3, 0.4, 0.3,
                  0.6, 0.3, 0.1), nrow = m, byrow = TRUE)
delta <- c(0.4, 0.2, 0.4)
mod2 <- list(m=m, mu=mu, sigma=sigma, gamma=gamma, delta=delta)
norm.sample <- norm.HMM.generate_sample(1000, mod2)
head(norm.sample)
graph_hmm_output(norm.sample)
mu0 <- c(5, 5, 5)
sigma0 <- c(5, 5, 5)
gamma0 <- matrix(rep(1/m, m*m), nrow = m, byrow = TRUE)
delta0 <- rep(1/m, m)
norm.mle <- norm.HMM.mle(norm.sample$obs, m, mu0, sigma0, gamma0, delta0, stationary=TRUE)
norm.mle
norm.decoding <- norm.HMM.viterbi(norm.sample$obs, mod2)
count(norm.decoding-norm.sample$state)
