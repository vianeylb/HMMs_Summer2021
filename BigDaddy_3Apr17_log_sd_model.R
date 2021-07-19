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

source("univariate_autoregressive_hmm_functions.R")
source("multivariate_autoregressive_hmm_functions.R")
sourceCpp("foralg.cpp")
sourceCpp("dmvnrm_arma.cpp")

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

#BigDaddy_3Apr17 log sd dynamic
filename <- "Log_SD_BigDaddy_3Apr17_dynamic.csv"
data <- read.csv(filename)
head(data)

subdata <- data %>% filter(Behavior == "Swimming")
x <- subdata$X_dynamic[-1]
x_lag <- subdata$X_dynamic[-length(subdata$X_dynamic)]
regx <- lm(x ~ x_lag)
summary(regx)
y <- subdata$Y_dynamic[-1]
y_lag <- subdata$Y_dynamic[-length(subdata$Y_dynamic)]
regy <- lm(y ~ y_lag)
summary(regy)
z <- subdata$Z_dynamic[-1]
z_lag <- subdata$Z_dynamic[-length(subdata$Z_dynamic)]
regz <- lm(z ~ z_lag)
summary(regz)
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

# 4-state model
# States: feeding, non-directed motion, resting, swimming, 
mu4 <- list(c(-4.05, -4.51, -3.93),
            c(-4.5, -4.86, -4.28),
            c(-6, -5.99, -5.64),
            c(-5.59, -5.7, -5.37))
sigma4 <- list(matrix(c(0.9, 0, 0, 0, 0.92, 0, 0, 0, 1.09), nrow = 3),
               matrix(c(1.07, 0, 0, 0, 0.97, 0, 0, 0, 1.14), nrow = 3),
               matrix(c(0.36, 0, 0, 0, 0.3, 0, 0, 0, 0.36), nrow = 3),
               matrix(c(0.99, 0, 0, 0, 0.67, 0, 0, 0, 0.94), nrow = 3))
gamma4 <- matrix(c(
  0.85, 0.05, 0.05, 0.05,
  0.05, 0.85, 0.05, 0.05,
  0.01, 0.01, 0.97, 0.01,
  0.01, 0.01, 0.01, 0.97
), nrow = 4, byrow = TRUE)
phi4 <- list(matrix(c(0.69, 0, 0, 0, 0.72, 0, 0, 0, 0.67), nrow = 3),
             matrix(c(0.75, 0, 0, 0, 0.72, 0, 0, 0, 0.73), nrow = 3),
             matrix(c(0.63, 0, 0, 0, 0.58, 0, 0, 0, 0.63), nrow = 3),
             matrix(c(0.89, 0, 0, 0, 0.88, 0, 0, 0, 0.86), nrow = 3))
delta4 <- c(1/4, 1/4, 1/4, 1/4)

start_time <- Sys.time()
log_sd_dynamic_mle4 <- mar_hmm_mle(obs, m = 4, q = 1, k = 3,
                                 mu4, sigma4, gamma4, phi4, delta4,
                                 stationary = FALSE, hessian = TRUE
)
end_time <- Sys.time()
runningtime2 <- end_time - start_time

log_sd_dynamic_mle4
save(log_sd_dynamic_mle4, file = "BigDaddy_3Apr17_log_sd_dynamic_mle4.RData")
load("BigDaddy_3Apr17_log_sd_dynamic_mle4.RData")
log_sd_dynamic_pr4 <- mar_hmm_pseudo_residuals(obs, log_sd_dynamic_mle4,
                                             "forecast", stationary = FALSE)
log_sd_dynamic_pr_plot4 <- pseudo_residual_plot(log_sd_dynamic_pr4)
log_sd_dynamic_decoding4 <- mar_hmm_viterbi(obs, log_sd_dynamic_mle4)
count(log_sd_dynamic_decoding4$state)


#BigDaddy_3Apr17 log sd static
filename <- "Log_SD_BigDaddy_3Apr17_static.csv"
data <- read.csv(filename)

subdata <- data %>% filter(Behavior == "Feeding")
x <- subdata$X_static[-1]
x_lag <- subdata$X_static[-length(subdata$X_static)]
regx <- lm(x ~ x_lag)
summary(regx)
y <- subdata$Y_static[-1]
y_lag <- subdata$Y_static[-length(subdata$Y_static)]
regy <- lm(y ~ y_lag)
summary(regy)
z <- subdata$Z_static[-1]
z_lag <- subdata$Z_static[-length(subdata$Z_static)]
regz <- lm(z ~ z_lag)
summary(regz)
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
mu4 <- list(c(-5.63, -6.28, -5.62),
            c(-5.67, -6.37, -5.56),
            c(-8.78, -8.84, -8.14),
            c(-7.9, -8.16, -7.59))
sigma4 <- list(matrix(c(0.95, 0, 0, 0, 1.11, 0, 0, 0, 1.16), nrow = 3),
               matrix(c(1.37, 0, 0, 0, 1.48, 0, 0, 0, 1.32), nrow = 3),
               matrix(c(0.87, 0, 0, 0, 0.75, 0, 0, 0, 0.83), nrow = 3),
               matrix(c(1.69, 0, 0, 0, 1.5, 0, 0, 0, 1.59), nrow = 3))
gamma4 <- matrix(c(
  0.85, 0.05, 0.05, 0.05,
  0.05, 0.85, 0.05, 0.05,
  0.01, 0.01, 0.97, 0.01,
  0.01, 0.01, 0.01, 0.97
), nrow = 4, byrow = TRUE)
phi4 <- list(matrix(c(0.96, 0, 0, 0, 0.93, 0, 0, 0, 0.98), nrow = 3),
             matrix(c(0.76, 0, 0, 0, 0.75, 0, 0, 0, 0.75), nrow = 3),
             matrix(c(0.68, 0, 0, 0, 0.59, 0, 0, 0, 0.66), nrow = 3),
             matrix(c(0.89, 0, 0, 0, 0.85, 0, 0, 0, 0.87), nrow = 3))
delta4 <- c(1/4, 1/4, 1/4, 1/4)

start_time <- Sys.time()
log_sd_static_mle4 <- mar_hmm_mle(obs, m = 4, q = 1, k = 3,
                                mu4, sigma4, gamma4, phi4, delta4,
                                stationary = FALSE, hessian = TRUE
)
end_time <- Sys.time()
end_time - start_time

log_sd_static_mle4
save(log_sd_static_mle4, file = "BigDaddy_3Apr17_log_sd_static_mle4.RData")
load("BigDaddy_3Apr17_log_sd_static_mle4.RData")
log_sd_static_pr4 <- mar_hmm_pseudo_residuals(obs, log_sd_static_mle4,
                                            "forecast",
                                            stationary = FALSE)
log_sd_static_pr_plot4 <- pseudo_residual_plot(log_sd_static_pr4)
log_sd_static_decoding4 <- mar_hmm_viterbi(obs, log_sd_static_mle4)
count(log_sd_static_decoding4$state)
