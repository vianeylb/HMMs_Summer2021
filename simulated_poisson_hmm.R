library(tidyverse)
library(ggplot2)
library(ggfortify)
library(plyr)
library(bayesforecast)
library(gridExtra)

source("univariate_poisson_hmm_functions.R")

# Testing functions for Poisson distributions
m <- 3
lambda <- c(2, 5, 8)
gamma <- matrix(c(
  0.1, 0.2, 0.7,
  0.3, 0.4, 0.3,
  0.6, 0.3, 0.1
), nrow = m, byrow = TRUE)
delta <- c(0.4, 0.2, 0.4)
mod <- list(m = m, lambda = lambda, gamma = gamma, delta = delta)
# Generate sample
pois_sample <- pois_hmm_generate_sample(1000, mod)
# Plots of sample
graph_hmm_output(pois_sample)
graph_hmm_hist(pois_sample)
ggacf(pois_sample$obs) +
  theme_minimal()
# Get MLE
lambda0 <- c(1, 4, 7)
gamma0 <- matrix(rep(1 / m, m * m), nrow = m, byrow = TRUE)
delta0 <- rep(1 / m, m)
pois_mle <- pois_hmm_mle(pois_sample$obs, m, lambda0, gamma0, delta0,
                         stationary = TRUE)
pois_mle
# Decode states
pois_decoding <- pois_hmm_viterbi(pois_sample$obs, pois_mle)
count(pois_decoding$state - pois_sample$state)
# Get marginal distribution
nx <- max(pois_sample$obs)
mpois <- pois_marginal(nx, pois_mle)
ggplot() +
  geom_histogram(data = pois_sample, aes(x = obs, y = ..density..),
                 binwidth = 1, colour = "navy", fill = "light blue") +
  geom_line(data = mpois, aes(x = x, y = mpois), colour = "red", size = 1) +
  theme_minimal()
# Get pseudo-residuals
pois_pr <- pois_hmm_pseudo_residuals(pois_sample$obs, pois_mle, "forecast")
pois_pr <- pois_pr %>% filter(lo != -Inf)
# Index plot of pseudo-residuals
ggplot(pois_pr) +
  geom_point(aes(x = seq_len(length(lo)), y = lo),
             size = 0.5, colour = "black") +
  geom_point(aes(x = seq_len(length(hi)), y = hi),
             size = 0.5, colour = "red") +
  theme_minimal()
# Histograms of pseudo-residuals
histlo <- ggplot(pois_pr, aes(lo)) +
  geom_histogram(aes(y = ..density..), colour = "navy", fill = "light blue") +
  stat_function(fun = dnorm, colour = "red") +
  theme_minimal()
histmi <- ggplot(pois_pr, aes(mi)) +
  geom_histogram(aes(y = ..density..), colour = "navy", fill = "light blue") +
  stat_function(fun = dnorm, colour = "red") +
  theme_minimal()
histhi <- ggplot(pois_pr, aes(hi)) +
  geom_histogram(aes(y = ..density..), colour = "navy", fill = "light blue") +
  stat_function(fun = dnorm, colour = "red") +
  theme_minimal()
grid.arrange(histlo, histmi, histhi, ncol = 3)
# QQ plot of pseudo-residuals
ggplot(pois_pr, aes(sample = mi)) +
  stat_qq() +
  stat_qq_line() +
  theme_minimal()
# ACF of pseudo-residuals
acflo <- ggacf(pois_pr$lo) +
  theme_minimal()
acfmi <- ggacf(pois_pr$mi) +
  theme_minimal()
acfhi <- ggacf(pois_pr$hi) +
  theme_minimal()
grid.arrange(acflo, acfmi, acfhi, ncol = 3)
