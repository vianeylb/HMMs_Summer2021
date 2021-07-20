library(tidyverse)
library(ggplot2)
library(ggfortify)
library(plyr)
library(bayesforecast)
library(gridExtra)

source("univariate_autoregressive_hmm_functions.R")

m <- 3
mu <- c(2, 5, 8)
sigma <- c(2, 4, 6)
gamma <- matrix(c(
  0.1, 0.2, 0.7,
  0.3, 0.4, 0.3,
  0.6, 0.3, 0.1
), nrow = m, byrow = TRUE)
delta <- c(0.4, 0.2, 0.4)
q <- 3
phi <- list(c(0.5, 0.3, 0.1),
            c(0.4, 0.2, 0.1),
            c(0.1, 0.1, 0.1))
mod <- list(m = m, mu = mu, sigma = sigma, gamma = gamma, phi = phi, delta = delta, q = q)
ar_sample <- ar_hmm_generate_sample(1000, mod)
# Plots 
graph_hmm_output(ar_sample)
graph_hmm_hist(ar_sample)
ggacf(ar_sample$obs) +
  theme_minimal()
ggpacf(ar_sample$obs) +
  theme_minimal()
# Get MLE
mu0 <- c(5, 5, 5)
sigma0 <- c(5, 5, 5)
gamma0 <- matrix(rep(1 / m, m * m), nrow = m, byrow = TRUE)
delta0 <- rep(1 / m, m)
phi0 <- list(rep(0.3, 3),
             rep(0.3, 3),
             rep(0.3, 3))
ar_mle <- ar_hmm_mle(ar_sample$obs, m, q, mu0, sigma0, gamma0, phi0,
                     delta0, stationary = FALSE, hessian = TRUE)
ar_mle
ar_decoding <- ar_hmm_viterbi(ar_sample$obs, ar_mle)
ar_decoding
count(ar_decoding$state)
count(ar_sample$state)
count(ar_decoding$state - ar_sample$state)
# Get pseudo-residuals
ar_pr <- ar_hmm_pseudo_residuals(ar_sample$obs, ar_mle,
                                 "forecast", stationary = FALSE)
# Index plot of pseudo-residuals
ggplot(ar_pr) +
  geom_point(aes(x = index, y = npsr), size = 0.5, colour = "black") +
  theme_minimal()
# Histogram of pseudo-residuals
ggplot(ar_pr, aes(npsr)) +
  geom_histogram(aes(y = ..density..), colour = "navy", fill = "light blue") +
  stat_function(fun = dnorm, colour = "red") +
  theme_minimal()
# QQ plot of pseudo-residuals
ggplot(ar_pr, aes(sample = npsr)) +
  stat_qq() +
  stat_qq_line() +
  theme_minimal()
# ACF of pseudo-residuals
ggacf(ar_pr$npsr) +
  theme_minimal()
# Getting confidence intervals via Hessian method
sd <- sqrt(diag(ar_mle$invhessian))
theta <- c(
  ar_mle$mu, ar_mle$sigma,
  as.vector(ar_mle$gamma[!diag(m)]),
  unlist(ar_mle$phi, use.names = FALSE)
)
norm_cih <- list(lower = theta - 1.96 * sd, upper = theta + 1.96 * sd)
norm_cih
# Getting confidence intervals via bootstrapping
n <- 3
len <- length(ar_sample$obs)
ar_bootstrap <- ar_bootstrap_estimates(ar_mle, n, len, stationary = FALSE)
ar_bootstrap
ar_cib <- ar_bootstrap_ci(ar_mle, ar_bootstrap, 0.05)
ar_cib