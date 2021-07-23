library(tidyverse)
library(ggplot2)
library(ggfortify)
library(plyr)
library(bayesforecast)
library(gridExtra)
library(mvtnorm)

source("indep_multivariate_autoregressive_hmm_functions.R")

m <- 3
k <- 2
q <- 3
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
phi1 <- matrix(c(0.4, 0.2, 0.1, 0.4, 0.2, 0.1), k, q, byrow = TRUE)
phi2 <- matrix(c(-0.4, -0.2, -0.1, -0.4, -0.2, -0.1), k, q, byrow = TRUE)
phi3 <- matrix(c(0.3, 0.2, 0.3, 0.1, 0.2, 0.1), k, q, byrow = TRUE)
phi <- list(phi1 = phi1, phi2 = phi2, phi3 = phi3)
mod <- list(m = m, k = k, mu = mu,
            sigma = sigma, gamma = gamma, phi = phi, delta = delta, q = q)
inmar_sample <- inmar_hmm_generate_sample(1000, mod)
mu0 <- list(c(5, 5), c(5, 5), c(5, 5))
sigma0 <- list(c(2, 2), c(2, 2), c(2, 2))
gamma0 <- matrix(rep(1 / m, m * m), nrow = m, byrow = TRUE)
delta0 <- rep(1 / m, m)
phi0 <- list(
  matrix(rep(0.2, k * q), nrow = k, byrow = TRUE),
  matrix(rep(0.2, k * q), nrow = k, byrow = TRUE),
  matrix(rep(0.2, k * q), nrow = k, byrow = TRUE)
)
start_time <- Sys.time()
inmar_mle <- inmar_hmm_mle(inmar_sample$obs, m, q, k,
                       mu0, sigma0, gamma0, phi0, delta = delta0,
                       stationary = FALSE, hessian = TRUE
)
end_time <- Sys.time()
inmar_mle
inmar_decoding <- inmar_hmm_viterbi(inmar_sample$obs, inmar_mle)
inmar_decoding
count(inmar_decoding$state)
count(inmar_sample$state)
count(inmar_decoding$state - inmar_sample$state)
# Get pseudo-residuals
inmar_pr <- inmar_hmm_pseudo_residuals(inmar_sample$obs, inmar_mle,
                                   "ordinary",
                                   stationary = FALSE
)
# Index plot of pseudo-residuals
ggplot(inmar_pr) +
  geom_point(aes(x = index, y = npsr), size = 0.5, colour = "black") +
  theme_minimal()
# Histogram of pseudo-residuals
ggplot(inmar_pr, aes(npsr)) +
  geom_histogram(aes(y = ..density..), colour = "navy", fill = "light blue") +
  stat_function(fun = dnorm, colour = "red") +
  theme_minimal()
# QQ plot of pseudo-residuals
ggplot(inmar_pr, aes(sample = npsr)) +
  stat_qq() +
  stat_qq_line() +
  theme_minimal()
# ACF of pseudo-residuals
ggacf(inmar_pr$npsr) +
  theme_minimal()

# Confidence intervals using Hessian
invhessian <- inmar_inv_hessian(inmar_mle, stationary = FALSE)
sd <- sqrt(diag(invhessian))
mu_estimate <- unlist(inmar_mle$mu, use.names = FALSE)
sigma_estimate <- unlist(inmar_mle$sigma, use.names = FALSE)
gamma_estimate <- as.vector(inmar_mle$gamma[!diag(m)])
phi_estimate <- unlist(inmar_mle$phi, use.names = FALSE)
theta <- c(mu_estimate, sigma_estimate, gamma_estimate, phi_estimate)
inmar_cih <- list(lower = theta - 1.96 * sd, upper = theta + 1.96 * sd)
inmar_cih

# Confidence intervals using bootstrapping
n <- 2
len <- length(inmar_sample$obs)
inmar_bootstrap <- inmar_bootstrap_estimates(inmar_mle, n, len, FALSE)
inmar_bootstrap
inmar_cib <- inmar_bootstrap_ci(mod, inmar_bootstrap, 0.05)
inmar_cib