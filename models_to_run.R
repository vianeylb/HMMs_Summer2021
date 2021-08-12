library(tidyverse)
library(ggplot2)
library(ggfortify)
library(plyr)
library(bayesforecast)
library(gridExtra)
library(mvtnorm)
library(parallel)
library(Rcpp)

source("univariate_autoregressive_hmm_functions.R")
source("indep_multivariate_autoregressive_hmm_functions.R")
sourceCpp("foralg.cpp")

# Returns dataframe with column that has behavior labels as positive integers
# labels is the integers assigned to Feeding, Non directed motion, Resting
# Swimming, Tag Roll respectively 
get_behavior_states <- function(data, labels){
  data$State[data$Behavior == "Feeding"] <- labels[1]
  data$State[data$Behavior == "Eating"] <- labels[1]
  data$State[data$Behavior == "Non directed motion"] <- labels[2]
  data$State[data$Behavior == "Resting"] <- labels[3]
  data$State[data$Behavior == "Swimming"] <- labels[4]
  data$State[data$Behavior == "Tag Roll"] <- labels[5]
  data$State[data$Behavior == "Unconfirmed tag roll"] <- labels[5]
  return(data)
}

get_confusion_matrix <- function(data, decoding){
  
  con_mat <- data_frame(Predicted = rep(0, 5), State_1 = rep(0, 5), State_2 = rep(0, 5),
                        State_3 = rep(0, 5), State_4 = rep(0, 5), Total = rep(0, 5))
  
  for (i in 1:4){
    for (j in 1:4){
      con_mat[i, j + 1] <- length(which(decoding$state == i & data$State  == j))
    }
  }
  
  con_mat$Total <- rowSums(con_mat)
  con_mat[5, ] <- as.list(colSums(con_mat))
  con_mat$Predicted <- c("State 1", "State 2", "State 3", "State 4", "Total")
  
  return(con_mat)
}

# 4 variable model
filename_dynamic <- "Mean_BigDaddy_3Apr17_dynamic.csv"
data_dynamic <- read.csv(filename_dynamic)
filename_static <- "Mean_BigDaddy_3Apr17_static.csv"
data_static <- read.csv(filename_static)

n <- length(data_dynamic$X_dynamic)
obs <- matrix(NA, 4, n)
obs[1, ] <- data_dynamic$X_dynamic
obs[2, ] <- data_static$X_static
obs[3, ] <- data_static$Y_static
obs[4, ] <- data_static$Z_static

mu4 <- list(c(0, -0.36, -0.89, 0.08),
            c(0, -0.4, -0.87, 0.13),
            c(0, -0.42, -0.88, -0.03),
            c(0, -0.41, -0.87, -0.02))
sigma4 <- list(c(0.03,  0.03, 0.02, 0.033),
               c(0.028, 0.13, 0.07, 0.16),
               c(0.004, 0.07, 0.04, 0.1),
               c(0.015,  0.08, 0.13, 0.18))
gamma4 <- matrix(c(
  0.85, 0.05, 0.05, 0.05,
  0.05, 0.85, 0.05, 0.05,
  0.01, 0.01, 0.97, 0.01,
  0.01, 0.01, 0.01, 0.97
), nrow = 4, byrow = TRUE)
phi4 <- list(matrix(c(-0.37, 0.96, 0.93, 0.98), nrow = 4),
             matrix(c(-0.33, 0.94, 0.96, 0.97), nrow = 4),
             matrix(c(-0.14, 1, 0.96, 1), nrow = 4),
             matrix(c(-0.21, 0.97, 0.9, 0.99), nrow = 4))
delta4 <- c(1/4, 1/4, 1/4, 1/4)

start_time <- Sys.time()
mean4_mle4 <- inmar_hmm_mle(obs, m = 4, q = 1, k = 4,
                            mu4, sigma4, gamma4, phi4, delta4,
                            stationary = FALSE, hessian = TRUE,
                            steptol = 0.001, iterlim = 400, stepmax = 20
)
end_time <- Sys.time()
end_time - start_time

mean4_mle4
save(mean4_mle4, file = "BigDaddy_3Apr17_indep_mean4_mle4.RData")
load("BigDaddy_3Apr17_indep_mean4_mle4.RData")

mean4_decoding4 <- inmar_hmm_viterbi(obs, mean4_mle4)
count(mean4_decoding4$state)
labels <- c(1, 2, 3, 4, 5)
data <- get_behavior_states(data, labels)
count(data$State)
count(mean4_decoding4$state - data$State)
con_mat <- get_confusion_matrix(data, mean4_decoding4) 
con_mat

# Dynamic + static model
filename_dynamic <- "Mean_BigDaddy_3Apr17_dynamic.csv"
data_dynamic <- read.csv(filename_dynamic)
filename_static <- "Mean_BigDaddy_3Apr17_static.csv"
data_static <- read.csv(filename_static)

n <- length(data_dynamic$X_dynamic)
obs <- matrix(NA, 6, n)
obs[1, ] <- data_dynamic$X_dynamic
obs[2, ] <- data_dynamic$Y_dynamic
obs[3, ] <- data_dynamic$Z_dynamic
obs[4, ] <- data_static$X_static
obs[5, ] <- data_static$Y_static
obs[6, ] <- data_static$Z_static

mu4 <- list(c(0, 0, 0, -0.36, -0.89, 0.08),
            c(0, 0, 0, -0.4, -0.87, 0.13),
            c(0, 0, 0, -0.42, -0.88, -0.03),
            c(0, 0, 0, -0.41, -0.87, -0.02))
sigma4 <- list(c(0.03, 0.02, 0.033, 0.03, 0.02, 0.033),
               c(0.028, 0.018, 0.025, 0.13, 0.07, 0.16),
               c(0.004, 0.003, 0.004, 0.07, 0.04, 0.1),
               c(0.015, 0.013, 0.018, 0.08, 0.13, 0.18))
gamma4 <- matrix(c(
  0.85, 0.05, 0.05, 0.05,
  0.05, 0.85, 0.05, 0.05,
  0.01, 0.01, 0.97, 0.01,
  0.01, 0.01, 0.01, 0.97
), nrow = 4, byrow = TRUE)
phi4 <- list(matrix(c(-0.37, -0.46, -0.47, 0.96, 0.93, 0.98), nrow = 6),
             matrix(c(-0.33, -0.34, -0.36, 0.94, 0.96, 0.97), nrow = 6),
             matrix(c(-0.14, -0.14, -0.12, 1, 0.96, 1), nrow = 6),
             matrix(c(-0.21, -0.23, -0.25, 0.97, 0.9, 0.99), nrow = 6))
delta4 <- c(1/4, 1/4, 1/4, 1/4)

start_time <- Sys.time()
mean_mle4 <- inmar_hmm_mle(obs, m = 4, q = 1, k = 6,
                           mu4, sigma4, gamma4, phi4, delta4,
                           stationary = FALSE, hessian = TRUE
)
end_time <- Sys.time()
end_time - start_time

mean_mle4
save(mean_mle4, file = "BigDaddy_3Apr17_indep_mean_mle4.RData")
load("BigDaddy_3Apr17_indep_mean_mle4.RData")


mean_decoding4 <- inmar_hmm_viterbi(obs, mean_mle4)
count(mean_decoding4$state)
labels <- c(1, 2, 3, 4, 5)
data <- get_behavior_states(data, labels)
count(data$State)
count(mean_decoding4$state - data$State)
con_mat <- get_confusion_matrix(data, mean_decoding4) 
con_mat