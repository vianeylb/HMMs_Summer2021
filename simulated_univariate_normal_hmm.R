library(tidyverse)
library(ggplot2)
library(ggfortify)
library(plyr)
library(bayesforecast)
library(gridExtra)

source("univariate_normal_hmm_functions.R")

# Testing functions for normal distributions
m <- 3
mu <- c(2, 5, 8)
sigma <- c(2, 4, 6)
gamma <- matrix(c(
  0.1, 0.2, 0.7,
  0.3, 0.4, 0.3,
  0.6, 0.3, 0.1
), nrow = m, byrow = TRUE)
delta <- c(0.4, 0.2, 0.4)
mod <- list(m = m, mu = mu, sigma = sigma, gamma = gamma, delta = delta)
# Generate sample
norm_sample <- norm_hmm_generate_sample(1000, mod)
# Plots of sample
graph_hmm_output(norm_sample)
graph_hmm_hist(norm_sample)
ggacf(norm_sample$obs) +
  theme_minimal()
# Get MLE
mu0 <- c(5, 5, 5)
sigma0 <- c(5, 5, 5)
gamma0 <- matrix(rep(1 / m, m * m), nrow = m, byrow = TRUE)
delta0 <- rep(1 / m, m)
norm_mle <- norm_hmm_mle(norm_sample$obs, m, mu0, sigma0, gamma0,
                         delta0, stationary = TRUE, hessian = TRUE)
norm_mle
# Decode states
norm_decoding <- norm_hmm_viterbi(norm_sample$obs, norm_mle)
count(norm_decoding$state - norm_sample$state)
# Get marginal distribution
start <- min(norm_sample$obs)
end <- max(norm_sample$obs)
mnorm <- norm_marginal(start, end, 1000, norm_mle)
ggplot() +
  geom_histogram(data = norm_sample, aes(x = obs, y = ..density..),
                 binwidth = 1, colour = "navy", fill = "light blue") +
  geom_line(data = mnorm, aes(x = x, y = mnorm),
            colour = "red", size = 1) +
  theme_minimal()
# Get pseudo-residuals
norm_pr <- norm_hmm_pseudo_residuals(norm_sample$obs, norm_mle,
                                     "forecast", stationary = TRUE)
# Index plot of pseudo-residuals
ggplot(norm_pr) +
  geom_point(aes(x = index, y = npsr), size = 0.5, colour = "black") +
  theme_minimal()
# Histogram of pseudo-residuals
ggplot(norm_pr, aes(npsr)) +
  geom_histogram(aes(y = ..density..), colour = "navy", fill = "light blue") +
  stat_function(fun = dnorm, colour = "red") +
  theme_minimal()
# QQ plot of pseudo-residuals
ggplot(norm_pr, aes(sample = npsr)) +
  stat_qq() +
  stat_qq_line() +
  theme_minimal()
# ACF of pseudo-residuals
ggacf(norm_pr$npsr) +
  theme_minimal()
# Getting confidence intervals via Hessian method
sd <- sqrt(diag(norm_mle$invhessian))
theta <- c(
  norm_mle$mu, norm_mle$sigma,
  as.vector(norm_mle$gamma[!diag(m)])
)
norm_cih <- list(lower = theta - 1.96 * sd, upper = theta + 1.96 * sd)
norm_cih
# Getting confidence intervals via bootstrapping
n <- 10
len <- length(norm_sample$obs)
norm_bootstrap <- norm_bootstrap_estimates(norm_mle, n, len, FALSE)
norm_cov <- norm_bootstrap_covariance(norm_bootstrap, m, n)
norm_cov
norm_cib <- norm_bootstrap_ci(norm_mle, norm_bootstrap, 0.05, m)
norm_cib

parvect <- norm_hmm_pn2pw(m, mu, sigma, gamma,
                        delta, stationary = FALSE)
pn <- norm_hmm_pw2pn(m, parvect, stationary = FALSE)
pn
