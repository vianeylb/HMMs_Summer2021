library(tidyverse)
library(ggplot2)
library(ggfortify)
library(plyr)
library(bayesforecast)
library(gridExtra)

source("Week 2 - Fitting HMM.R")
source("Week 3 - Functions.R")

#Testing functions for Poisson distributions
m <- 3
lambda <- c(2, 5, 8)
gamma <- matrix(c(0.1, 0.2, 0.7,
                  0.3, 0.4, 0.3,
                  0.6, 0.3, 0.1), nrow = m, byrow = TRUE)
delta <- c(0.4, 0.2, 0.4)
mod <- list(m=m, lambda=lambda, gamma=gamma, delta=delta)
#Generate sample
pois.sample <- pois.HMM.generate_sample(1000, mod)
head(pois.sample)
#Plots of sample
graph_hmm_output(pois.sample)
graph_hmm_hist(pois.sample)
ggacf(pois.sample$obs)+
  theme_minimal()
#Get MLE
lambda0 <- c(1, 4, 7)
gamma0 <- matrix(rep(1/m, m*m), nrow = m, byrow = TRUE)
delta0 <- rep(1/m, m)
pois.mle <- pois.HMM.mle(pois.sample$obs, m, lambda0, gamma0, delta0, stationary=TRUE)
pois.mle
#Decode states
pois.decoding <- pois.HMM.viterbi(pois.sample$obs, pois.mle)
count(pois.decoding$state-pois.sample$state)
#Get marginal distribution
nx <- max(pois.sample$obs)
mpois <- pois.marginal(nx, pois.mle)
ggplot()+
  geom_histogram(data=pois.sample, aes(x=obs, y=..density..), binwidth=1, colour="navy", fill="light blue")+
  geom_line(data=mpois, aes(x=x,y=mpois), colour="red", size=1)+
  theme_minimal()
#Get pseudo-residuals
pois.pr <- pois.HMM.pseudo_residuals(pois.sample$obs, pois.mle, "forecast") 
pois.pr <- pois.pr %>% filter(lo != -Inf)
view(pois.opr)
#Index plot of pseudo-residuals
ggplot(pois.pr)+
  geom_point(aes(x=c(1:length(lo)), y=lo), size=0.5, colour="black")+
  geom_point(aes(x=c(1:length(hi)), y=hi), size=0.5, colour="red")+
  theme_minimal()
#Histograms of pseudo-residuals
histlo <- ggplot(pois.pr, aes(lo))+
  geom_histogram(aes(y=..density..),  colour="navy", fill="light blue")+
  stat_function(fun=dnorm, colour="red")+
  theme_minimal()
histmi <- ggplot(pois.pr, aes(mi))+
  geom_histogram(aes(y=..density..), colour="navy", fill="light blue")+
  stat_function(fun=dnorm, colour="red")+
  theme_minimal()
histhi <- ggplot(pois.pr, aes(hi))+
  geom_histogram(aes(y=..density..), colour="navy", fill="light blue")+
  stat_function(fun=dnorm, colour="red")+
  theme_minimal()
grid.arrange(histlo, histmi, histhi, ncol=3)
#QQ plot of pseudo-residuals
ggplot(pois.pr, aes(sample=mi))+
  stat_qq()+
  stat_qq_line()+
  theme_minimal()
#ACF of pseudo-residuals
acflo <- ggacf(pois.pr$lo)+
  theme_minimal()
acfmi <- ggacf(pois.pr$mi)+
  theme_minimal()
acfhi <- ggacf(pois.pr$hi)+
  theme_minimal()
grid.arrange(acflo, acfmi, acfhi, ncol=3)

#Testing functions for normal distributions
m <- 3
mu <- c(2, 5, 8)
sigma <- c(2, 4, 6)
gamma <- matrix(c(0.1, 0.2, 0.7,
                  0.3, 0.4, 0.3,
                  0.6, 0.3, 0.1), nrow = m, byrow = TRUE)
delta <- c(0.4, 0.2, 0.4)
mod2 <- list(m=m, mu=mu, sigma=sigma, gamma=gamma, delta=delta)
#Generate sample
norm.sample <- norm.HMM.generate_sample(1000, mod2)
head(norm.sample)
#Plots of sample
graph_hmm_output(norm.sample)
graph_hmm_hist(norm.sample)
ggacf(norm.sample$obs)+
  theme_minimal()
#Get MLE
mu0 <- c(5, 5, 5)
sigma0 <- c(5, 5, 5)
gamma0 <- matrix(rep(1/m, m*m), nrow = m, byrow = TRUE)
delta0 <- rep(1/m, m)
norm.mle <- norm.HMM.mle(norm.sample$obs, m, mu0, sigma0, gamma0, delta0, stationary=TRUE)
norm.mle
#Decode states
norm.decoding <- norm.HMM.viterbi(norm.sample$obs, norm.mle)
count(norm.decoding$states-norm.sample$state)
#Get marginal distribution
start <- min(norm.sample$obs)
end <- max(norm.sample$obs)
mnorm <- norm.marginal(start, end, 1000, norm.mle)
ggplot()+
  geom_histogram(data=norm.sample, aes(x=obs, y=..density..), binwidth=1, colour="navy", fill="light blue")+
  geom_line(data=mnorm, aes(x=x,y=mnorm), colour="red", size=1)+
  theme_minimal()
#Get pseudo-residuals
norm.pr <- norm.HMM.pseudo_residuals(norm.sample$obs, norm.mle, "forecast") 
#Index plot of pseudo-residuals
ggplot(norm.pr)+
  geom_point(aes(x=index), y=npsr), size=0.5, colour="black")+
  theme_minimal()
#Histogram of pseudo-residuals
ggplot(norm.pr, aes(npsr))+
  geom_histogram(aes(y=..density..),  colour="navy", fill="light blue")+
  stat_function(fun=dnorm, colour="red")+
  theme_minimal()
#QQ plot of pseudo-residuals
ggplot(norm.pr, aes(sample=npsr))+
  stat_qq()+
  stat_qq_line()+
  theme_minimal()
#ACF of pseudo-residuals
ggacf(norm.pr$npsr)+
  theme_minimal()
