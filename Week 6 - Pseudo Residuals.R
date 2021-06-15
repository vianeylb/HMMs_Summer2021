library(tidyverse)
library(ggplot2)
library(ggfortify)
library(plyr)

#Computing log(forward probabilities) for normal distribution
norm.HMM.lforward <- function(x, mod)
{
  n <- length(x)
  lalpha <- matrix(NA, mod$m, n)
  foo <- mod$delta*dnorm(x[1], mod$mu, mod$sigma)
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo/sumfoo
  lalpha [ ,1] <- lscale + log(foo)
  for (i in 2:n)
  {
    foo <- foo%*%mod$gamma*dnorm(x[i], mod$mu, mod$sigma)
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo/sumfoo
    lalpha[,i] <- log(foo) + lscale
  }
  return(lalpha)
}

#Computing log(backward probabilities) for normal distribution
norm.HMM.lbackward <- function(x, mod)
{
  n <- length(x)
  m <- mod$m
  lbeta <- matrix(NA,m,n)
  lbeta[,n] <- rep(0,m)
  foo <- rep(1/m,m)
  lscale <- log(m)
  for (i in (n-1):1)
  {
    foo <- mod$gamma%*%(dnorm(x[i+1], mod$mu, mod$sigma)*foo)
    lbeta [,i] <- log(foo) + lscale
    sumfoo <- sum(foo)
    foo <- foo/sumfoo
    lscale <- lscale + log(sumfoo)
  }
  return(lbeta)
}


#Normal pseudo-residuals for Normal HMM
#Type can be "ordinary" or "forecast"
norm.HMM.pseudo_residuals <- function(x, mod, type)
{
  if (type=="ordinary"){
    n <- length(x)
    la <- norm.HMM.lforward(x, mod)
    lb <- norm.HMM.lbackward(x, mod)
    lafact <- apply(la, 2, max)
    lbfact <- apply(lb, 2, max)
    
    P <- matrix(NA, n, mod$m)
    for (i in 1:n){
      P[i,] <- pnorm(x[i], mean=mod$mu, sd=mod$sigma)
    }
    
    npsr <- rep(NA, n)
    npsr[1] <- qnorm(mod$delta %*% P[1,])
    for (i in 2:n) {
      a <- exp(la[,i-1]-lafact[i])
      b <- exp(lb[,i]-lbfact[i])
      foo <- (a%*%mod$gamma)*b
      foo <- foo/sum(foo)
      npsr[i] <- qnorm(foo%*%P[i,])
    }
    
    return(data.frame(npsr, x, index=c(1:n)))
  }
  else if (type=="forecast"){
    n <- length(x)
    la <- norm.HMM.lforward(x, mod)
    
    P <- matrix(NA, n, mod$m)
    for (i in 1:n){
      P[i,] <- pnorm(x[i], mean=mod$mu, sd=mod$sigma)
    }
    
    npsr <- rep(NA, n)
    npsr[1] <- qnorm(mod$delta %*% P[1,])
    for (i in 2:n) {
      la_max <- max(la[,i-1])
      a <- exp(la[,i-1]-la_max)
      npsr[i] <- qnorm(t(a)%*%(gamma/sum(a))%*%P[i,])
    }
    
    return(data.frame(npsr, x, index=c(1:n)))
  }
}

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
norm.pr <- norm.HMM.pseudo_residuals(norm.sample$obs, norm.mle, "ordinary") 
#Index plot of pseudo-residuals
ggplot(norm.pr)+
  geom_point(aes(x=index, y=npsr), size=0.5, colour="black")+
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
