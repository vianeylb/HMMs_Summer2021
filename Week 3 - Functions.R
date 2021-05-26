library(tidyverse)
library(ggplot2)
library(ggfortify)
library(plyr)

source("Week 2 - Fitting HMM.R")

#Get Poisson marginal distribution 
pois.marginal <- function(n, mod){
  #Get stationary distribution
  delta <-solve(t(diag(mod$m)-mod$gamma+1),rep(1,mod$m))
  mpois <- delta[1]*dpois(x=0:n, lambda=mod$lambda[1])
  for (i in 2:mod$m){
    mpois <- mpois + delta[i]*dpois(x=0:n, lambda=mod$lambda[i])
  }
  return(data.frame(x=0:n, mpois=mpois))
}

#Computing log(forward probabilities) for Poisson distribution
#See Section 4.1.1
pois.HMM.lforward <- function(x, mod)
  {
    n <- length(x)
    #Empty matrix of size mxn
    lalpha <- matrix(NA, mod$m, n)
    foo <- mod$delta*dpois(x[1], mod$lambda)
    sumfoo <- sum(foo)
    lscale <- log(sumfoo)
    foo <- foo/sumfoo
    #Set first column of the matrix lalpha
    lalpha [ ,1] <- lscale + log(foo)
    for (i in 2:n)
      {
        foo <- foo%*%mod$gamma*dpois(x[i], mod$lambda)
        sumfoo <- sum(foo)
        lscale <- lscale + log(sumfoo)
        foo <- foo/sumfoo
        #Set ith column of the matrix lalpha to alpha_i
        lalpha[,i] <- log(foo) + lscale
      }
    return(lalpha)
}

#Computing log(backward probabilities) for Poisson distribution
#See Section 4.1.2
pois.HMM.lbackward <- function(x, mod)
  {
    n <- length(x)
    m <- mod$m
    lbeta <- matrix(NA,m,n)
    #Set last column of matrix lbeta to 0
    lbeta[,n] <- rep(0,m)
    foo <- rep(1/m,m)
    lscale <- log(m)
    for (i in (n-1):1)
      {
        foo <- mod$gamma %*%(dpois (x[i+1] , mod$lambda)*foo)
        #Set ith column of the matrix lbeta to beta_i
        lbeta [,i] <- log(foo) + lscale
        sumfoo <- sum(foo)
        foo <- foo/sumfoo
        lscale <- lscale + log(sumfoo)
      }
    return(lbeta)
}

#Conditional probabilities for Poisson distribution
#Conditional probability that observation at time t equals
#xc , given all observations other than that at time t.
#Note: xc is a vector and the result (dxc) is a matrix .
pois.HMM.conditional <- function(xc, x, mod)
  {
    n <- length(x)
    m <- mod$m
    nxc <- length(xc)
    dxc <- matrix(NA, nrow=nxc, ncol=n)
    Px <- matrix (NA, nrow=m, ncol=nxc)
    #Set elements of Px to p_i(xc_j)
    for (j in 1:nxc) Px[,j]<-dpois(xc[j], mod$lambda)
    la <- pois.HMM.lforward(x, mod)
    lb <- pois.HMM.lbackward(x, mod)
    la <- cbind(log(mod$delta), la)
    #Get maximum value in each column 
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
pois.HMM.conditional2 <- function(xc, x, mod)
{
  n <- length(x)
  m <- mod$m
  nxc <- length(xc)
  dxc <- matrix(NA, nrow=nxc, ncol=n)
  Px <- matrix (NA, nrow=m, ncol=nxc)
  for (j in 1:nxc) Px[,j]<-dpois(xc[j], mod$lambda)
  la <- pois.HMM.lforward(x, mod)
  la <- cbind(log(mod$delta), la)
  lafact <- apply(la, 2, max)
  for (i in 2:n)
  {
    foo <- exp(la[,i]-lafact[i])
    foo <- foo/sum(foo)
    dxc[,i] <- foo%*%mod$gamma%*%Px
  }
  return(dxc)
}

#Normal pseudo-residuals for Poisson HMM
#Type can be "ordinary" or "forecast"
#Output is lower, upper and middle pseudo-residuals
pois.HMM.pseudo_residuals <- function(x, mod, type)
{
  n <- length(x)
  if (type=="ordinary"){
    cdists <- pois.HMM.conditional(xc=0:max(x), x, mod)
  }
  else if (type=="forecast"){
    cdists <- pois.HMM.conditional2(xc=0:max(x), x, mod)
  }
  #Each column of cumdists starts with 0, remaining elements i are cumulative sum of 
  #elements in the column of cdists up to that cell 
  cumdists <- rbind(rep(0,n), apply(cdists, 2, cumsum))
  ulo <- uhi <- rep(NA, n)
  for (i in 1:n)
  {
    ulo[i] <- cumdists[x[i]+1 ,i]
    uhi[i] <- cumdists[x[i]+2 ,i]
  }
  umi <- 0.5*(ulo+uhi)
  npsr <- qnorm(rbind(ulo, umi, uhi))
  return(data.frame(lo=npsr[1,], mi=npsr[2,], hi=npsr[3,]))
}

#Get normal marginal distribution 
norm.marginal <- function(start, end, n, mod){
  #Get stationary distribution
  delta <-solve(t(diag(mod$m)-mod$gamma+1),rep(1,mod$m))
  x <- seq(start, end, length.out = n)
  mnorm <- delta[1]*dnorm(x, mean=mod$mu[1], sd=mod$sigma[1])
  for (i in 2:mod$m){
    mnorm <- mnorm + delta[i]*dnorm(x, mean=mod$mu[i], sd=mod$sigma[i])
  }
  return(data.frame(x=x, mnorm=mnorm))
}

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

#Conditional probabilities for normal distribution
#Conditional probability that observation at time t equals
#xc , given all observations other than that at time t.
#Note: xc is a vector and the result (dxc) is a matrix .
norm.HMM.conditional <- function(xc, x, mod)
{
  n <- length(x)
  m <- mod$m
  nxc <- length(xc)
  dxc <- matrix(NA, nrow=nxc, ncol=n)
  Px <- matrix (NA, nrow=m, ncol=nxc)
  for (j in 1:nxc) Px[,j]<-dnorm(xc[j], mod$mu, mod$sigma)
  la <- norm.HMM.lforward(x, mod)
  lb <- norm.HMM.lbackward(x, mod)
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
norm.HMM.conditional2 <- function(xc, x, mod)
{
  n <- length(x)
  m <- mod$m
  nxc <- length(xc)
  dxc <- matrix(NA, nrow=nxc, ncol=n)
  Px <- matrix (NA, nrow=m, ncol=nxc)
  for (j in 1:nxc) Px[,j]<-dnorm(xc[j], mod$mu, mod$sigma)
  la <- norm.HMM.lforward(x, mod)
  la <- cbind(log(mod$delta), la)
  lafact <- apply(la, 2, max)
  for (i in 2:n)
  {
    foo <- exp(la[,i]-lafact[i])
    foo <- foo/sum(foo)
    dxc[,i] <- foo%*%mod$gamma%*%Px
  }
  return(dxc)
}

#Normal pseudo-residuals for Normal HMM
#Type can be "ordinary" or "forecast"
norm.HMM.pseudo_residuals <- function(x, mod, type)
{
  n <- length(x)
  xord <- sort(x)
  if (type=="ordinary"){
    cdists <- norm.HMM.conditional(xc=xord, x, mod)
  }
  else if (type=="forecast"){
    cdists <- norm.HMM.conditional2(xc=xord, x, mod)
  }
  xdiff <- c(0, xord[2:n] - xord[1:n-1])
  cdists2 <- xdiff*cdists
  cdf <- apply(cdists2, 2, cumsum)
  cdf <- cdf + matrix(rep(1-cdf[n,],n),nrow=n, byrow=TRUE)
  upr <- rep(NA, n)
  for (i in 1:n)
  {
    index <- which(xord==x[i])
    upr[i] <- cdf[index,i]
  }
  npsr <- qnorm(upr)
  return(data.frame(npsr, x))
}
