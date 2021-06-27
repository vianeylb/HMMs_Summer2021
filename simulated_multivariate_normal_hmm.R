library(tidyverse)
library(ggplot2)
library(ggfortify)
library(plyr)
library(bayesforecast)
library(gridExtra)
library(mvtnorm)

source("multivariate_normal_hmm_functions.R")

# Fitting simulated 3-state bivariate normal HMM
m <- 3
k <- 2
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
mod <- list(m = m, k = k, mu = mu, sigma = sigma, gamma = gamma, delta = delta)
mvnorm_sample <- mvnorm_hmm_generate_sample(1000, mod)
mu0 <- list(c(5, 5), c(5, 5), c(5, 5))
sigma0 <- list(diag(2, 2), diag(2, 2), diag(2, 2))
gamma0 <- matrix(rep(1 / m, m * m), nrow = m, byrow = TRUE)
delta0 <- rep(1 / m, m)
mvnorm_mle <- mvnorm_hmm_mle(mvnorm_sample$obs, m, k,
                             mu0, sigma0, gamma0, delta0,
                             stationary = FALSE, hessian = TRUE)
mvnorm_mle
mvnorm_decoding <- mvnorm_hmm_viterbi(mvnorm_sample$obs, mvnorm_mle)
count(mvnorm_decoding$state)

# Get psuedo-residuals
mvnorm_pr <- mvnorm_hmm_pseudo_residuals(mvnorm_sample$obs,
                                         mvnorm_mle, k, "ordinary")
# Index plot of pseudo-residuals
ggplot(mvnorm_pr) +
  geom_point(aes(x = index, y = npsr), size = 0.5, colour = "black") +
  theme_minimal()
# Histogram of pseudo-residuals
ggplot(mvnorm_pr, aes(npsr)) +
  geom_histogram(aes(y = ..density..), colour = "navy", fill = "light blue") +
  stat_function(fun = dnorm, colour = "red") +
  theme_minimal()
# QQ plot of pseudo-residuals
ggplot(mvnorm_pr, aes(sample = npsr)) +
  stat_qq() +
  stat_qq_line() +
  theme_minimal()
# ACF of pseudo-residuals
ggacf(mvnorm_pr$npsr) +
  theme_minimal()

# Confidence intervals using Hessian
sd <- sqrt(diag(mvnorm_mle$invhessian))
mu_estimate <- unlist(mvnorm_mle$mu, use.names = FALSE)
sigma_estimate <- mvnorm_mle$sigma
for (i in 1:m) {
  foo <- sigma_estimate[[i]]
  sigma_estimate[[i]] <- foo[lower.tri(foo, diag = TRUE)]
}
sigma_estimate <- unlist(sigma_estimate, use.names = FALSE)
gamma_estimate <- as.vector(mvnorm_mle$gamma[!diag(m)])
theta <- c(mu_estimate, sigma_estimate, gamma_estimate)
mvnorm_cih <- list(lower = theta - 1.96 * sd, upper = theta + 1.96 * sd)
mvnorm_cih

# Confidence intervals using bootstrapping
n <- 2
len <- length(mvnorm_sample$obs)
mvnorm_bootstrap <- mvnorm_bootstrap_estimates(mvnorm_mle, n, k, len, FALSE)
mvnorm_cib <- mvnorm_bootstrap_ci(mvnorm_mle, mvnorm_bootstrap, 0.05, m, k)
mvnorm_cib
