library(tidyverse)
library(ggplot2)
library(ggfortify)
library(plyr)

source("Week 2 - Fitting HMM.R")
source("Week 3 - Functions.R")

#Using bootstrapping
norm.bootstrap.estimates <- function(mod, n, stationary){
  mu_estimate <- numeric()
  sigma_estimate <- numeric()
  gamma_estimate <- numeric()
  delta_estimate <- numeric()
  for (i in 1:n){
    sample <- norm.HMM.generate_sample(1000, mod)
    mod2 <- norm.HMM.mle(sample$obs, m, mod$mu, mod$sigma, mod$gamma, mod$delta, stationary=stationary)
    mu_estimate <- c(mu_estimate, mod2$mu) 
    sigma_estimate <- c(sigma_estimate, mod2$sigma) 
    gamma_estimate <- c(gamma_estimate, mod2$gamma) 
    delta_estimate <- c(delta_estimate, mod2$delta) 
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
norm.mle <- norm.HMM.mle(norm.sample$obs, m, mu0, sigma0, gamma0, delta0, stationary=FALSE)
norm.mle
#Testing bootstrapping 
n <- 10
norm.bootstrap <- norm.bootstrap.estimates(norm.mle, n, FALSE)
norm.bootstrap$mu
norm.cov <- norm.bootstrap.covariance(norm.bootstrap, m, n)
norm.cov
norm.ci <- norm.bootstrap.ci(mod, norm.bootstrap, 0.05, m, n)
norm.ci

