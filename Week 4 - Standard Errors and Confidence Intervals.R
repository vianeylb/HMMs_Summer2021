library(tidyverse)
library(ggplot2)
library(ggfortify)
library(plyr)

source("Week 2 - Fitting HMM.R")
source("Week 3 - Functions.R")

#Using Hessian
norm.HMM.mle <- function(x, m, mu0, sigma0, gamma0, delta0=NULL, stationary=TRUE,...)
{
  parvect0 <- norm.HMM.pn2pw(m, mu, sigma0, gamma0, delta0, stationary=stationary)
  mod <- nlm(norm.HMM.mllk, parvect0, x=x,m=m, stationary=stationary, hessian=TRUE)
  pn <- norm.HMM.pw2pn(m=m, mod$estimate, stationary=stationary)
  mllk <- mod$minimum
  if (det(mod$hessian) != 0){
    h <- solve(mod$hessian)
    jacobian <- norm.jacobian(m, mod$estimate, stationary=stationary) 
    h <- t(jacobian)%*%h%*%jacobian
  }
  np <- length(parvect0)
  AIC <- 2*(mllk+np)
  n <- sum(!is.na(x))
  BIC <- 2*mllk+np*log(n)
  if (det(mod$hessian) != 0){
    list(m=m, mu=pn$mu, sigma=pn$sigma, gamma=pn$gamma, delta=pn$delta,
         code=mod$code, mllk=mllk, AIC=AIC, BIC=BIC, invhessian=h)
  }
  else{
    list(m=m, mu=pn$mu, sigma=pn$sigma, gamma=pn$gamma, delta=pn$delta,
         code=mod$code, mllk=mllk, AIC=AIC, BIC=BIC)
  }
}

norm.jacobian <- function(m, parvect, stationary=TRUE){
  jacobian <- matrix(0, nrow=length(parvect), ncol=length(parvect))
  jacobian[1:m, 1:m] <- diag(m)
  jacobian[(m+1):(2*m), (m+1):(2*m)] <- diag(exp(parvect[(m+1):(2*m)]))
  gamma0 <- diag(m)
  gamma0[!gamma0] <- exp(parvect[(2*m+1):(m+m*m)])
  count <- 0
  for (i in 1:m){
    for (j in 1:m){
      if (j != i){
        count <- count + 1
        k <- (sum(gamma0[i,]))
        foo <- -(gamma0[i,]*gamma0[i,j])/(k^2)
        foo[j] <- (gamma0[i,j]*k - gamma0[i,j]^2)/(k^2)
        foo <- foo[-i]
        jacobian[2*m+count, (2*m+(i-1)*(m-1)+1):(2*m+i*(m-1))] <- foo
      }
    }
  }
  if(!stationary)
  {
    delta0 <- c(1, exp(parvect[(m+m*m+1):(m*m+2*m-1)]))
    k <- (sum(gamma0[i,]))
    for (i in 2:m){
      foo <- -(delta0*delta0[i])/(k^2)
      foo[i] <- (delta0[i]*k - delta0[i]^2)/(k^2)
      foo <- foo[-1]
      jacobian[(m+1)*m, ((m+1)*m):((m+2)*m-1)]<- foo
    }
  }
  return(jacobian)
}

#Using bootstrapping
norm.bootstrap.estimates <- function(mod, n, len, stationary){
  m <- mod$m
  mu_estimate <- numeric(n*m)
  sigma_estimate <- numeric(n*m)
  gamma_estimate <- numeric(n*m*m)
  delta_estimate <- numeric(n*m)
  for (i in 1:n){
    sample <- norm.HMM.generate_sample(len, mod)
    mod2 <- norm.HMM.mle(sample$obs, m, mod$mu, mod$sigma, mod$gamma, mod$delta, stationary=stationary)
    mu_estimate[((i-1)*m+1):(i*m)] <- mod2$mu
    sigma_estimate[((i-1)*m+1):(i*m)] <- mod2$sigma
    gamma_estimate[((i-1)*m*m+1):(i*m*m)] <- mod2$gamma
    delta_estimate[((i-1)*m+1):(i*m)] <- mod2$delta
  }
  return(list(mu=mu_estimate, sigma=sigma_estimate, gamma=gamma_estimate, delta=delta_estimate))
}

norm.bootstrap.covariance <- function(bootstrap, m, n){
  size <- (m+3)*m
  cov <- matrix(rep(0, size*size), size)
  foo <- rep(0, size)
  for (i in 1:n){
    estimates <- c(bootstrap$mu[((i-1)*m+1):(i*m)], bootstrap$sigma[((i-1)*m+1):(i*m)], bootstrap$gamma[((i-1)*m*m+1):(i*m*m)], bootstrap$delta[((i-1)*m+1):(i*m)])
    foo <- foo + estimates
  }
  foo <- foo/n
  for (i in 1:n){
    estimates <- c(bootstrap$mu[((i-1)*m+1):(i*m)], bootstrap$sigma[((i-1)*m+1):(i*m)], bootstrap$gamma[((i-1)*m*m+1):(i*m*m)], bootstrap$delta[((i-1)*m+1):(i*m)])
    cov <- cov + ((estimates - foo) %o% (estimates - foo))
  }
  cov <- cov/(n-1)
  return(cov)
}

norm.bootstrap.ci <- function(mod, bootstrap, alpha, m, n){
  mu_lower <- rep(NA, m)
  mu_upper <- rep(NA, m)
  sigma_lower <- rep(NA, m)
  sigma_upper <- rep(NA, m)
  gamma_lower <- rep(NA, m*m)
  gamma_upper <- rep(NA, m*m)
  delta_lower <- rep(NA, m)
  delta_upper <- rep(NA, m)
  bootstrap1 <- data.frame(mu=bootstrap$mu, sigma=bootstrap$sigma, delta=bootstrap$delta)
  bootstrap2 <- data.frame(gamma=bootstrap$gamma)
  for (i in 0:(m-1)){
    foo <- bootstrap1 %>% filter((row_number() %% m) == i)
    mu_lower[m-i] <- 2*mod$mu[m-i] - quantile(foo$mu, 1-(alpha/2), names=FALSE)
    mu_upper[m-i] <- 2*mod$mu[m-i] - quantile(foo$mu, alpha/2, names=FALSE)
    sigma_lower[m-i] <- 2*mod$sigma[m-i] - quantile(foo$sigma, 1-(alpha/2), names=FALSE)
    sigma_upper[m-i] <- 2*mod$sigma[m-i] - quantile(foo$sigma, alpha/2, names=FALSE)
    delta_lower[m-i] <- 2*mod$delta[m-i] - quantile(foo$delta, 1-(alpha/2), names=FALSE)
    delta_upper[m-i] <- 2*mod$delta[m-i] - quantile(foo$delta, alpha/2, names=FALSE)
  }
  for (i in 0:(m*m-1)){
    foo <- bootstrap2 %>% filter((row_number() %% (m*m)) == i)                                                                                                                                                                                                                                                                                                                                                 
    gamma_lower[m-i] <- 2*mod$gamma[m-i] - quantile(foo$gamma, 1-(alpha/2), names=FALSE)
    gamma_upper[m-i] <- 2*mod$gamma[m-i] - quantile(foo$gamma, alpha/2, names=FALSE)
  }
  return(list(mu_lower=mu_lower, mu_upper=mu_upper, 
              sigma_lower=sigma_lower, sigma_upper=sigma_upper, 
              gamma_lower=gamma_lower, gamma_upper=gamma_upper, 
              delta_lower=delta_lower, delta_upper=delta_upper))
}

#Testing 
m <- 3
mu <- c(2, 5, 8)
sigma <- c(2, 4, 6)
gamma <- matrix(c(0.1, 0.2, 0.7,
                  0.3, 0.4, 0.3,
                  0.6, 0.3, 0.1), nrow = m, byrow = TRUE)
delta <- c(0.4, 0.2, 0.4)
mod <- list(m=m, mu=mu, sigma=sigma, gamma=gamma, delta=delta)
norm.sample <- norm.HMM.generate_sample(1000, mod)
mu0 <- c(5, 5, 5)
sigma0 <- c(5, 5, 5)
gamma0 <- matrix(rep(1/m, m*m), nrow = m, byrow = TRUE)
delta0 <- rep(1/m, m)
norm.mle <- norm.HMM.mle(norm.sample$obs, m, mu0, sigma0, gamma0, delta0, stationary=TRUE)
norm.mle
sd <- sqrt(diag(norm.mle$invhessian))
theta <- c(norm.mle$mu, norm.mle$sigma, 
           norm.mle$gamma[1,2:3],
           norm.mle$gamma[2,-2],
           norm.mle$gamma[3, 1:2])
norm.cih <- list(lower=theta-1.96*sd, upper=theta+1.96*sd)
norm.cih
#Testing bootstrapping 
n <- 100
norm.bootstrap <- norm.bootstrap.estimates(norm.mle, n, TRUE)
norm.cov <- norm.bootstrap.covariance(norm.bootstrap, m, n)
norm.cov
norm.cib <- norm.bootstrap.ci(norm.mle, norm.bootstrap, 0.05, m, n)
norm.cib
