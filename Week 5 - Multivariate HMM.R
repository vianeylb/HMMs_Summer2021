library(tidyverse)
library(ggplot2)
library(ggfortify)
library(plyr)
library(mvtnorm)

setwd("C:/Jessica/UofT Y4/Research/Coding")

#Changes for multivariate normal
#mu is a list of vectors
#sigma is a list of matrices 
#parvect is a list
#k is the number of variables

#Transform normal natural parameters to working parameters
mvnorm.HMM.pn2pw <- function(m,mu,sigma,gamma,delta=NULL,stationary=TRUE)
{
  #Put all means into one vector
  mu <- unlist(mu, use.names = FALSE)
  #Only need to transform diagonal elements of sigma
  #Include only lower triangle & diag of matrix, since covariance matrix must be symmetric
  tsigma <- lapply(sigma, diag_log_lower)
  tsigma <- unlist(tsigma, use.names = FALSE)
  foo <- log(gamma/diag(gamma))
  tgamma <- as.vector(foo[!diag(m)])
  if(stationary) 
  {tdelta <-NULL} 
  else 
  {tdelta<-log(delta[-1]/delta[1])}
  #parvect is one vector of values
  parvect <- c(mu,tsigma,tgamma,tdelta)
  return(parvect)
}

#Applies log to only diagonal elements of matrix,
#then returns only lower triangular entries of matrix as vector
#going down column by column
diag_log_lower <- function(mat){
  diag(mat) <- log(diag(mat))
  vect <- mat[lower.tri(mat, diag=TRUE)]
  return(vect)
}

#Applies exp to only diagonal elements of matrix
diag_exp <- function(mat){
  diag(mat) <- exp(diag(mat))
  return(mat)
}

#Get nth triangular number (0 is 0th num)
triangular_num <- function(n){
  nums <- choose(seq(n+1),2)
  return(nums[n+1])
}

# Transform normal working parameters to natural parameters
mvnorm.HMM.pw2pn <- function(m,k,parvect,stationary=TRUE)
{
  #Change mu to list of vectors format
  mu <- list()
  count <- 1
  for (i in 1:m){
    mu[[i]] <- parvect[count:(i*k)]
    count <- count + k
  }

  #Change sigma to list of matrices format
  tsigma <- list()
  #Get number of elements in lower triangle (including diag) of matrix
  t <- triangular_num(k)
  for (i in 1:m){
    tsigma_vals <- parvect[count:(count + t - 1)]
    foo <- diag(k)
    foo[lower.tri(foo, diag=TRUE)] <- tsigma_vals
    foo <- t(foo)
    foo[lower.tri(foo, diag=TRUE)] <- tsigma_vals
    tsigma[[i]] <- foo
    count <- count + t
  }
  sigma <- lapply(tsigma, diag_exp)
  
  tgamma <- parvect[count:(count + m*(m-1) - 1)]
  count <- count + m*(m-1)
  gamma <- diag(m)
  gamma[!gamma] <- exp(tgamma)
  gamma <- gamma/apply(gamma,1,sum)
  
  if(stationary) 
  {delta<-solve(t(diag(m)-gamma+1),rep(1,m))} 
  else {
    tdelta <- parvect[count:(count + m - 2)]
    foo<-c(1,exp(tdelta))
    delta<-foo/sum(foo)
  }
  return(list(mu=mu,sigma=sigma,gamma=gamma,delta=delta))
}

#Computing minus the log-likelihood from the working parameters
mvnorm.HMM.mllk <- function(parvect, x, m, k, stationary=TRUE, ...)
{
  n <- length(x)
  pn <- mvnorm.HMM.pw2pn(m, k, parvect, stationary=stationary)
  P <- get_densities(x[[1]], pn, m)
  foo <- pn$delta*P
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo/sumfoo
  for (i in 2:n)
  {
    P <- get_densities(x[[i]], pn, m)
    foo <- foo%*%pn$gamma*P
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo/sumfoo
  }
  mllk <- -lscale
  return(mllk)
}

#Get state dependent probability densities given x and pn
get_densities <- function(x, pn, m){
  pvect <- numeric(m)
  for (i in 1:m){
    pvect[i] <- dmvnorm(x, pn$mu[[i]], pn$sigma[[i]])
  }
  return(pvect)
}

#Computing MLE from natural parameters
mvnorm.HMM.mle <- function(x, m, k, mu0, sigma0, gamma0, delta0=NULL, stationary=TRUE, hessian=TRUE,...)
{
  parvect0 <- mvnorm.HMM.pn2pw(m=m, mu=mu0, sigma=sigma0, gamma=gamma0, delta=delta0, stationary=stationary)
  mod <- nlm(mvnorm.HMM.mllk, parvect0, x=x,m=m, k=k,stationary=stationary)
  pn <- mvnorm.HMM.pw2pn(m=m, k=k, parvect=mod$estimate, stationary=stationary)
  mllk <- mod$minimum
  np <- length(parvect0)
  AIC <- 2*(mllk+np)
  n2 <- sum(!is.na(x))
  BIC <- 2*mllk+np*log(n2)
  list(m=m, mu=pn$mu, sigma=pn$sigma, gamma=pn$gamma, delta=pn$delta,
       code=mod$code, mllk=mllk, AIC=AIC, BIC=BIC)
}

#Generating a sample from normal distributions
mvnorm.HMM.generate_sample <- function(ns, mod)
{
  mvect <- 1:mod$m
  state <- numeric(ns)
  state[1] <- sample(mvect, 1, prob=mod$delta)
  if (ns>1){
    for (i in 2:ns) state[i] <- sample(mvect, 1, prob=mod$gamma[state[i-1],])
  }
  x <- lapply(state, mvnorm.HMM.sample_one, mod=mod)
  return(list(index=c(1:ns), state=state, obs=x)) 
}

mvnorm.HMM.sample_one <- function(state, mod){
  x <- rmvnorm(1, mean=mod$mu[[state]], sigma=mod$sigma[[state]])
  return(x)
}

#Global decoding by the Viterbi algorithm
mvnorm.HMM.viterbi <- function(x, mod)
{
  n <- length(x)
  xi <- matrix(0, n, mod$m)
  P <- get_densities(x[[1]], mod, m)
  foo <- mod$delta*P
  xi[1,] <- foo/sum(foo)
  for (t in 2:n)
  {
    P <- get_densities(x[[t]], mod, m)
    foo <- apply(xi[t-1,]*mod$gamma, 2, max)*P
    xi[t,] <- foo/sum(foo)
  }
  iv <- numeric(n)
  iv[n] <- which.max(xi[n,])
  for (t in (n-1):1)
    iv[t] <- which.max(mod$gamma[,iv[t+1]]*xi[t,])
  return(data.frame(index=1:n, states=iv))
}

#Computing log(forward probabilities) for normal distribution
#Here x is a list of matricies 
mvnorm.HMM.lforward <- function(x, mod)
{
  n <- length(x)
  lalpha <- matrix(NA, mod$m, n)
  P <- get_densities(x[[1]], mod, m)
  foo <- mod$delta*P
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo/sumfoo
  lalpha [ ,1] <- lscale + log(foo)
  for (i in 2:n)
  {
    P <- get_densities(x[[i]], mod, m)
    foo <- foo%*%mod$gamma*P
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo/sumfoo
    lalpha[,i] <- log(foo) + lscale
  }
  return(lalpha)
}

#Computing log(backward probabilities) for normal distribution
mvnorm.HMM.lbackward <- function(x, mod)
{
  n <- length(x)
  m <- mod$m
  lbeta <- matrix(NA,m,n)
  lbeta[,n] <- rep(0,m)
  foo <- rep(1/m,m)
  lscale <- log(m)
  for (i in (n-1):1)
  {
    P <- get_densities(x[[i+1]], mod, m)
    foo <- mod$gamma%*%(P*foo)
    lbeta [,i] <- log(foo) + lscale
    sumfoo <- sum(foo)
    foo <- foo/sumfoo
    lscale <- lscale + log(sumfoo)
  }
  return(lbeta)
}

#Conditional probabilities for normal distribution
#Conditional probability that observation at time t equals
#xc , given all observations other than that at time t.
#Note: xc is a vector and the result (dxc) is a matrix .
mvnorm.HMM.conditional <- function(xc, x, mod)
{
  n <- length(x)
  m <- mod$m
  nxc <- length(xc)
  dxc <- matrix(NA, nrow=nxc, ncol=n)
  Px <- matrix (NA, nrow=m, ncol=nxc)
  for (j in 1:nxc) Px[,j]<-get_densities(xc[[j]], mod, m)
  la <- mvnorm.HMM.lforward(x, mod)
  lb <- mvnorm.HMM.lbackward(x, mod)
  la <- cbind(log(mod$delta), la)
  lafact <- apply(la, 2, max)
  lbfact <- apply(lb, 2, max)
  for (i in 1:n)
  {
    foo <- (exp(la[,i]-lafact[i])%*%mod$gamma)*exp(lb[,i]-lbfact[i])
    foo <- foo/sum(foo)
    dxc[,i] <- foo%*%Px
  }
  return(dxc)
}

#Conditional probability that observation at time t equals
#xc, given observations at time 1 to t-1. 
#Note: xc is a vector and the result (dxc) is a matrix .
mvnorm.HMM.conditional2 <- function(xc, x, mod)
{
  n <- length(x)
  m <- mod$m
  nxc <- length(xc)
  dxc <- matrix(NA, nrow=nxc, ncol=n)
  Px <- matrix(NA, nrow=m, ncol=nxc)
  for (j in 1:nxc) Px[,j]<-get_densities(xc[[j]], mod, m)
  la <- mvnorm.HMM.lforward(x, mod)
  la <- cbind(log(mod$delta), la)
  lafact <- apply(la, 2, max)
  for (i in 1:n)
  {
    foo <- exp(la[,i]-lafact[i])
    foo <- foo/sum(foo)
    dxc[,i] <- foo%*%mod$gamma%*%Px
  }
  return(dxc)
}

m <- 3
k <- 2
mu1 <- c(0, 2)
mu2 <- c(5, 8)
mu3 <- c(8, 10)
mu <- list(mu1=mu1, mu2=mu2, mu3=mu3)
sigma1 <- diag(2)
sigma2 <- diag(2,2)
sigma3 <- diag(2) + 3
sigma <- list(sigma1=sigma1, sigma2=sigma2, sigma3=sigma3)
gamma <- matrix(c(0.1, 0.2, 0.7,
                  0.3, 0.4, 0.3,
                  0.6, 0.3, 0.1), nrow = m, byrow = TRUE)
delta <- c(0.4, 0.2, 0.4)
mod <- list(m=m, k=k, mu=mu, sigma=sigma, gamma=gamma, delta=delta)
mvnorm.sample <- mvnorm.HMM.generate_sample(1000, mod)
head(mvnorm.sample)
mu0 <- list(c(5, 5), c(5, 5), c(5, 5))
sigma0 <- list(diag(2, 2), diag(2, 2), diag(2, 2))
gamma0 <- matrix(rep(1/m, m*m), nrow = m, byrow = TRUE)
delta0 <- rep(1/m, m)
mvnorm.mle <- mvnorm.HMM.mle(mvnorm.sample$obs, m, k, mu0, sigma0, gamma0, delta0, stationary=FALSE)
mvnorm.mle 
mvnorm.decoding <- mvnorm.HMM.viterbi(mvnorm.sample$obs, mvnorm.mle)
head(mvnorm.decoding)
count(mvnorm.decoding$states)
