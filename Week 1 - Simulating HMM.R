library(tidyverse)
library(ggplot2)

sim_hmm_norm <- function(t, N, transition, initial, mu, sigma){
  df = data.frame(index = c(1:t), state = numeric(N), observation = numeric(N))
  
  state = sample(x = c(1:N), size = 1, prob = initial)
  observation = rnorm(1, mu[state], sigma[state])
  df$state[1] <- state
  df$observation[1] <- observation
  
  for (i in 2:t){
    probabilities = transition[state,]
    state = sample(x = c(1:N), size = 1, prob = probabilities)
    observation = rnorm(1, mu[state], sigma[state])
    df$state[i] <- state
    df$observation[i] <- observation
  }
  return(df)
}

sim_hmm_pois <- function(t, N, transition, initial, lambda){
  df = data.frame(index = c(1:t), state = numeric(N), observation = numeric(N))
  
  state = sample(x = c(1:N), size = 1, prob = initial)
  observation = rpois(1, lambda[state])
  df$state[1] <- state
  df$observation[1] <- observation
  
  for (i in 2:t){
    probabilities = transition[state,]
    state = sample(x = c(1:N), size = 1, prob = probabilities)
    observation = rnorm(1, lambda[state])
    df$state[i] <- state
    df$observation[i] <- observation
  }
  return(df)
}

graph_hmm_output <- function(output){
  ggplot(data = output, mapping = aes(x=output$index, y=output$observation, 
                                      color = output$state))+
    geom_line()+
    theme_minimal()+
    scale_colour_continuous(type="viridis")
}

transition <- matrix(c(0.8, 0.05, 0.05, 0.05, 0.05,
                       0.05, 0.8, 0.05, 0.05, 0.05,
                       0.05, 0.05, 0.8, 0.05, 0.05,
                       0.05, 0.05, 0.05, 0.8, 0.05,
                       0.05, 0.05, 0.05, 0.05, 0.8), nrow = 5, byrow = TRUE)
initial <- c(0.25, 0.25, 0.25, 0.25, 0.25)
mu <- c(1, 2, 3, 4, 5)
sigma <- c(1, 1, 1, 1, 1)
output <- sim_hmm_norm(1000, 5, transition, initial, mu, sigma) 
head(output)
graph_hmm_output(output)

lambda <- c(1, 2, 3, 4, 5)
output2 <- sim_hmm_pois(1000, 5, transition, initial, lambda) 
head(output2)
graph_hmm_output(output2)
