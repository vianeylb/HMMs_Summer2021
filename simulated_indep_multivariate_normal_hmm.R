library(tidyverse)
library(ggplot2)
library(ggfortify)
library(plyr)
library(bayesforecast)
library(gridExtra)
library(mvtnorm)

setwd("C:/Jessica/UofT Y4/Research/Coding")
source("indep_multivariate_normal_hmm_functions.R")

# Fitting simulated 3-state bivariate normal HMM
m <- 3
k <- 2
mu1 <- c(0, 2)
mu2 <- c(5, 8)
mu3 <- c(8, 10)
mu <- list(mu1 = mu1, mu2 = mu2, mu3 = mu3)
sigma1 <- c(1, 2)
sigma2 <- c(2, 3)
sigma3 <- c(2, 1)
sigma <- list(sigma1 = sigma1, sigma2 = sigma2, sigma3 = sigma3)
gamma <- matrix(c(
  0.1, 0.2, 0.7,
  0.3, 0.4, 0.3,
  0.6, 0.3, 0.1
), nrow = m, byrow = TRUE)
delta <- c(0.4, 0.2, 0.4)
mod <- list(m = m, k = k, mu = mu, sigma = sigma, gamma = gamma, delta = delta)
inmvnorm_sample <- inmvnorm_hmm_generate_sample(1000, mod)
mu0 <- list(c(5, 5), c(5, 5), c(5, 5))
sigma0 <- list(c(2, 2), c(2, 2), c(2, 2))
gamma0 <- matrix(rep(1 / m, m * m), nrow = m, byrow = TRUE)
delta0 <- rep(1 / m, m)
inmvnorm_mle <- inmvnorm_hmm_mle(inmvnorm_sample$obs, m, k,
                             mu0, sigma0, gamma, delta0,
                             stationary = FALSE, hessian = TRUE)
inmvnorm_mle
inmvnorm_decoding <- inmvnorm_hmm_viterbi(inmvnorm_sample$obs, inmvnorm_mle)
count(inmvnorm_decoding$state)
count(inmvnorm_sample$state)
count(inmvnorm_decoding$state - inmvnorm_sample$state)

# Get psuedo-residuals
inmvnorm_pr <- inmvnorm_hmm_pseudo_residuals(inmvnorm_sample$obs,
                                         inmvnorm_mle, "forecast", stationary = FALSE)
# Index plot of pseudo-residuals
ggplot(inmvnorm_pr) +
  geom_point(aes(x = index, y = npsr), size = 0.5, colour = "black") +
  theme_minimal()
# Histogram of pseudo-residuals
ggplot(inmvnorm_pr, aes(npsr)) +
  geom_histogram(aes(y = ..density..), colour = "navy", fill = "light blue") +
  stat_function(fun = dnorm, colour = "red") +
  theme_minimal()
# QQ plot of pseudo-residuals
ggplot(inmvnorm_pr, aes(sample = npsr)) +
  stat_qq() +
  stat_qq_line() +
  theme_minimal()
# ACF of pseudo-residuals
ggacf(inmvnorm_pr$npsr) +
  theme_minimal()

# Confidence intervals using Hessian
invhessian <- inmvnorm_inv_hessian(inmvnorm_mle, stationary = FALSE)
sd <- sqrt(diag(invhessian))
mu_estimate <- unlist(inmvnorm_mle$mu, use.names = FALSE)
sigma_estimate <- unlist(inmvnorm_mle$sigma, use.names = FALSE)
gamma_estimate <- as.vector(inmvnorm_mle$gamma[!diag(m)])
theta <- c(mu_estimate, sigma_estimate, gamma_estimate)
inmvnorm_cih <- list(lower = theta - 1.96 * sd, upper = theta + 1.96 * sd)
inmvnorm_cih

# Confidence intervals using bootstrapping
n <- 10
len <- length(inmvnorm_sample$obs)
inmvnorm_bootstrap <- inmvnorm_bootstrap_estimates(inmvnorm_mle, n, k, len, FALSE)
inmvnorm_cib <- inmvnorm_bootstrap_ci(inmvnorm_mle, inmvnorm_bootstrap, 0.05, m, k)
inmvnorm_cib

x <- inmvnorm_sample$obs
inmvnorm_sample$obs
parvect <- inmvnorm_hmm_pn2pw(m, mu, sigma, gamma, delta, stationary = FALSE)
parvect
pn <- inmvnorm_hmm_pw2pn(m, k, parvect, stationary = FALSE)
pn
pmat <- inmvnorm_densities(x, mod, m, k, 1000)
inmvnorm_hmm_mllk(parvect, x, m, k, stationary = FALSE)
inmvnorm_hmm_lforward(x, mod)[, 1:6]
mvnorm_hmm_lforward(x, mod2)[, 1:6]
inmvnorm_dist_mat(x, mod, 1000)
inmvnorm_jacobian(mod, 18)
