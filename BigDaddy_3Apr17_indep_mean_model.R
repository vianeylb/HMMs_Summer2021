library(tidyverse)
library(ggplot2)
library(ggfortify)
library(plyr)
library(bayesforecast)
library(gridExtra)
library(mvtnorm)
library(parallel)
library(Rcpp)

setwd("C:/Jessica/UofT Y4/Research/Coding")

source("indep_multivariate_autoregressive_hmm_functions.R")
sourceCpp("foralg.cpp")

setwd("C:/Jessica/UofT Y4/Research/Coding/Lab Data")

pseudo_residual_plot <- function(data) {
  # Index plot of pseudo-residuals
  plot_index <- ggplot(data) +
    geom_point(aes(x = index, y = npsr), size = 0.5, colour = "black") +
    theme_minimal()
  # Histogram of pseudo-residuals
  plot_hist <- ggplot(data, aes(npsr)) +
    geom_histogram(aes(y = ..density..), colour = "navy", fill = "light blue") +
    stat_function(fun = dnorm, colour = "red") +
    theme_minimal()
  # QQ plot of pseudo-residuals
  plot_qq <- ggplot(data, aes(sample = npsr)) +
    stat_qq() +
    stat_qq_line() +
    theme_minimal()
  # ACF of pseudo-residuals
  plot_acf <- ggacf(data$npsr) +
    theme_minimal()
  plot <- grid.arrange(plot_index, plot_hist, plot_qq, plot_acf,
                       nrow = 2, ncol = 2)
  return(plot)
}

#BigDaddy_3Apr17 mean dynamic
filename <- "Mean_BigDaddy_3Apr17_dynamic.csv"
data <- read.csv(filename)
head(data)

subdata <- data %>% filter(Behavior == "Swimming")
x <- subdata$X_dynamic[-1]
x_lag <- subdata$X_dynamic[-length(subdata$X_dynamic)]
regx <- lm(x ~ x_lag)
suminmary(regx)
y <- subdata$Y_dynamic[-1]
y_lag <- subdata$Y_dynamic[-length(subdata$Y_dynamic)]
regy <- lm(y ~ y_lag)
suminmary(regy)
z <- subdata$Z_dynamic[-1]
z_lag <- subdata$Z_dynamic[-length(subdata$Z_dynamic)]
regz <- lm(z ~ z_lag)
suminmary(regz)
mean(subdata$X_dynamic)
mean(subdata$Y_dynamic)
mean(subdata$Z_dynamic)
sd(subdata$X_dynamic)
sd(subdata$Y_dynamic)
sd(subdata$Z_dynamic)

n <- length(data$X_dynamic)
obs <- matrix(NA, 3, n)
obs[1, ] <- data$X_dynamic
obs[2, ] <- data$Y_dynamic
obs[3, ] <- data$Z_dynamic

mu4 <- list(c(0, 0, 0), c(0, 0, 0), c(0, 0, 0), c(0, 0, 0))
sigma4 <- list(c(0.03, 0.02, 0.033),
               c(0.028, 0.018, 0.025),
               c(0.004, 0.003, 0.004),
               c(0.015, 0.013, 0.018))
gamma4 <- matrix(c(
  0.85, 0.05, 0.05, 0.05,
  0.05, 0.85, 0.05, 0.05,
  0.01, 0.01, 0.97, 0.01,
  0.01, 0.01, 0.01, 0.97
), nrow = 4, byrow = TRUE)
phi4 <- list(matrix(c(-0.37, -0.46, -0.47), nrow = 3),
             matrix(c(-0.33, -0.34, -0.36), nrow = 3),
             matrix(c(-0.14, -0.14, -0.12), nrow = 3),
             matrix(c(-0.21, -0.23, -0.25), nrow = 3))
delta4 <- c(1/4, 1/4, 1/4, 1/4)

start_time <- Sys.time()
mean_dynamic_mle4 <- ininmar_hmm_mle(obs, m = 4, q = 1, k = 3,
                                 mu4, sigma4, gamma4, phi4, delta4,
                                 stationary = FALSE, hessian = TRUE
)
end_time <- Sys.time()
end_time - start_time

mean_dynamic_mle4
save(mean_dynamic_mle4, file = "BigDaddy_3Apr17_indep_mean_dynamic_mle4.RData")
load("BigDaddy_3Apr17_indep_mean_dynamic_mle4.RData")
mean_dynamic_pr4 <- inmar_hmm_pseudo_residuals(obs, mean_dynamic_mle4,
                                             "forecast",
                                             stationary = FALSE)
mean_dynamic_pr_plot4 <- pseudo_residual_plot(mean_dynamic_pr4)
ggsave("BigDaddy_3Apr17_indep_mean_dynamic_pr_plot4.png", mean_dynamic_pr_plot4)
mean_dynamic_decoding4 <- inmar_hmm_viterbi(obs, mean_dynamic_mle4)
count(mean_dynamic_decoding4$state)

#BigDaddy_3Apr17 mean static
filename <- "Mean_BigDaddy_3Apr17_static.csv"
data <- read.csv(filename)

subdata <- data %>% filter(Behavior == "Swimming")
x <- subdata$X_static[-1]
x_lag <- subdata$X_static[-length(subdata$X_static)]
regx <- lm(x ~ x_lag)
suminmary(regx)
y <- subdata$Y_static[-1]
y_lag <- subdata$Y_static[-length(subdata$Y_static)]
regy <- lm(y ~ y_lag)
suminmary(regy)
z <- subdata$Z_static[-1]
z_lag <- subdata$Z_static[-length(subdata$Z_static)]
regz <- lm(z ~ z_lag)
suminmary(regz)
mean(subdata$X_static)
mean(subdata$Y_static)
mean(subdata$Z_static)
sd(subdata$X_static)
sd(subdata$Y_static)
sd(subdata$Z_static)

n <- length(data$X_static)
obs <- matrix(NA, 3, n)
obs[1, ] <- data$X_static
obs[2, ] <- data$Y_static
obs[3, ] <- data$Z_static

# 4-state model
# States: feeding, non-directed motion, resting, swimming 
mu4 <- list(c(-0.36, -0.89, 0.08),
            c(-0.4, -0.87, 0.13),
            c(-0.42, -0.88, -0.03),
            c(-0.41, -0.87, -0.02))
sigma4 <- list(c(0.03, 0.02, 0.033),
               c(0.13, 0.07, 0.16),
               c(0.07, 0.04, 0.1),
               c(0.08, 0.13, 0.18))
gamma4 <- matrix(c(
  0.85, 0.05, 0.05, 0.05,
  0.05, 0.85, 0.05, 0.05,
  0.01, 0.01, 0.97, 0.01,
  0.01, 0.01, 0.01, 0.97
), nrow = 4, byrow = TRUE)
phi4 <- list(matrix(c(0.96, 0.93, 0.98), nrow = 3),
             matrix(c(0.94, 0.96, 0.97), nrow = 3),
             matrix(c(1, 0.96, 1), nrow = 3),
             matrix(c(0.97, 0.9, 0.99), nrow = 3))
delta4 <- c(1/4, 1/4, 1/4, 1/4)

start_time <- Sys.time()
mean_static_mle4 <- inmar_hmm_mle(obs, m = 4, q = 1, k = 3,
                                mu4, sigma4, gamma4, phi4, delta4,
                                stationary = FALSE, hessian = TRUE
)
end_time <- Sys.time()
end_time - start_time

mean_static_mle4
save(mean_static_mle4, file = "BigDaddy_3Apr17_mean_static_mle4.RData")
load("BigDaddy_3Apr17_indep_mean_static_mle4.RData")
mean_static_pr4 <- inmar_hmm_pseudo_residuals(obs, mean_static_mle4,
                                            "forecast",
                                            stationary = FALSE)
mean_static_pr_plot4 <- pseudo_residual_plot(mean_static_pr4)
ggsave("BigDaddy_3Apr17_indep_mean_static_pr_plot4.png", mean_static_pr_plot4)
mean_static_decoding4 <- inmar_hmm_viterbi(obs, mean_static_mle4)
count(mean_static_decoding4$state)