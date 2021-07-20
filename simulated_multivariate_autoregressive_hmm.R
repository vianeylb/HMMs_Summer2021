library(tidyverse)
library(ggplot2)
library(ggfortify)
library(plyr)
library(bayesforecast)
library(gridExtra)
library(mvtnorm)

source("multivariate_autoregressive_hmm_functions.R")

m <- 3
k <- 2
q <- 3
mu1 <- c(0, 2)
mu2 <- c(5, 8)
mu3 <- c(8, 10)
mu <- list(mu1 = mu1, mu2 = mu2, mu3 = mu3)
sigma1 <- diag(2)
sigma2 <- diag(2, 2)
sigma3 <- diag(2) + 3
sigma <- list(sigma1 = sigma1, sigma2 = sigma2, sigma3 = sigma3)
gamma <- matrix(c(
  0.1, 0.2, 0.7,
  0.3, 0.4, 0.3,
  0.6, 0.3, 0.1
), nrow = m, byrow = TRUE)
delta <- c(0.4, 0.2, 0.4)
phi1 <- matrix(nrow = k, ncol = k * q)
phi1[, 1:2] <- diag(0.4, 2)
phi1[, 3:4] <- diag(0.2, 2)
phi1[, 5:6] <- diag(0.1, 2)
phi2 <- matrix(nrow = k, ncol = k * q)
phi2[, 1:2] <- diag(0.4, 2) + 0.1
phi2[, 3:4] <- diag(0.2, 2) + 0.1
phi2[, 5:6] <- diag(0.1, 2) + 0.1
phi3 <- matrix(nrow = k, ncol = k * q)
phi3[, 1:2] <- diag(0.5, 2) 
phi3[, 3:4] <- diag(0.3, 2) 
phi3[, 5:6] <- diag(0.1, 2) 
phi <- list(phi1 = phi1, phi2 = phi2, phi3 = phi3)
mod <- list(m = m, k = k, mu = mu, sigma = sigma, gamma = gamma, phi = phi, delta = delta, q = q)
mar_sample <- mar_hmm_generate_sample(1000, mod)
mu0 <- list(c(5, 5), c(5, 5), c(5, 5))
sigma0 <- list(diag(2, 2), diag(2, 2), diag(2, 2))
gamma0 <- matrix(rep(1 / m, m * m), nrow = m, byrow = TRUE)
delta0 <- rep(1 / m, m)
phi0 <- list(matrix(rep(0.2, k * k * q), nrow = k, byrow = TRUE),
             matrix(rep(0.2, k * k * q), nrow = k, byrow = TRUE),
             matrix(rep(0.2, k * k * q), nrow = k, byrow = TRUE))
mar_mle <- mar_hmm_mle(mar_sample$obs, m, q, k,
                       mu0, sigma0, gamma0, phi, delta0,
                       stationary = FALSE, hessian = TRUE)
mar_mle
mar_decoding <- mar_hmm_viterbi(mar_sample$obs, mod)
mar_decoding
count(mar_decoding$state)
count(mar_sample$state)
count(mar_decoding$state - mar_sample$state)
# Get pseudo-residuals
mar_pr <- mar_hmm_pseudo_residuals(mar_sample$obs, mar_mle,
                                   "forecast", stationary = FALSE)
view(mar_pr)
# Index plot of pseudo-residuals
ggplot(mar_pr) +
  geom_point(aes(x = index, y = npsr), size = 0.5, colour = "black") +
  theme_minimal()
# Histogram of pseudo-residuals
ggplot(mar_pr, aes(npsr)) +
  geom_histogram(aes(y = ..density..), colour = "navy", fill = "light blue") +
  stat_function(fun = dnorm, colour = "red") +
  theme_minimal()
# QQ plot of pseudo-residuals
ggplot(mar_pr, aes(sample = npsr)) +
  stat_qq() +
  stat_qq_line() +
  theme_minimal()
# ACF of pseudo-residuals
ggacf(mar_pr$npsr) +
  theme_minimal()

