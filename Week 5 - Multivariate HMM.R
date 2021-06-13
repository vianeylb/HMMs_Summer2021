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
  P <- mvnorm.densities(x[[1]], pn, m)
  foo <- pn$delta*P
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo/sumfoo
  for (i in 2:n)
  {
    P <- mvnorm.densities(x[[i]], pn, m)
    foo <- foo%*%pn$gamma*P
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo/sumfoo
  }
  mllk <- -lscale
  return(mllk)
}

#Get state dependent probability densities given x and mod
mvnorm.densities <- function(x, mod, m){
  pvect <- numeric(m)
  for (i in 1:m){
    pvect[i] <- dmvnorm(x, mod$mu[[i]], mod$sigma[[i]])
  }
  return(pvect)
}

#Computing MLE from natural parameters
mvnorm.HMM.mle <- function(x, m, k, mu0, sigma0, gamma0, delta0=NULL, stationary=TRUE, hessian=FALSE,...)
{
  parvect0 <- mvnorm.HMM.pn2pw(m=m, mu=mu0, sigma=sigma0, gamma=gamma0, delta=delta0, stationary=stationary)
  mod <- nlm(mvnorm.HMM.mllk, parvect0, x=x,m=m, k=k,stationary=stationary, hessian=hessian)
  pn <- mvnorm.HMM.pw2pn(m=m, k=k, parvect=mod$estimate, stationary=stationary)
  mllk <- mod$minimum
  
  np <- length(parvect0)
  AIC <- 2*(mllk+np)
  n2 <- sum(!is.na(x))
  BIC <- 2*mllk+np*log(n2)
  
  if (hessian){
    if (!stationary){
      np2 <- np - m + 1
      h <- mod$hessian[1:np2,1:np2]
    }
    else{
      np2 <- np
      h <- mod$hessian
    }
    if (det(h) != 0){
      h <- solve(h)
      jacobian <- mvnorm.jacobian(m, k, n=np2, mod$estimate, stationary=stationary) 
      h <- t(jacobian)%*%h%*%jacobian
      return(list(m=m, mu=pn$mu, sigma=pn$sigma, gamma=pn$gamma, delta=pn$delta,
                  code=mod$code, mllk=mllk, AIC=AIC, BIC=BIC, invhessian=h))
    }
    else{
      return(list(m=m, mu=pn$mu, sigma=pn$sigma, gamma=pn$gamma, delta=pn$delta,
                  code=mod$code, mllk=mllk, AIC=AIC, BIC=BIC, hessian=mod$hessian))
    }
  }
  else{
    return(list(m=m, mu=pn$mu, sigma=pn$sigma, gamma=pn$gamma, delta=pn$delta,
                code=mod$code, mllk=mllk, AIC=AIC, BIC=BIC))
  }
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
  P <- mvnorm.densities(x[[1]], mod, mod$m)
  foo <- mod$delta*P
  xi[1,] <- foo/sum(foo)
  for (t in 2:n)
  {
    P <- mvnorm.densities(x[[t]], mod, mod$m)
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
mvnorm.HMM.lforward <- function(x, mod)
{
  n <- length(x)
  lalpha <- matrix(NA, mod$m, n)
  foo <- mod$delta*mvnorm.densities(x[[1]], mod, mod$m)
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo/sumfoo
  lalpha [ ,1] <- lscale + log(foo)
  for (i in 2:n)
  {
    foo <- foo%*%mod$gamma*mvnorm.densities(x[[i]], mod, mod$m)
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
    foo <- mod$gamma%*%(mvnorm.densities(x[[i+1]], mod, mod$m)*foo)
    lbeta [,i] <- log(foo) + lscale
    sumfoo <- sum(foo)
    foo <- foo/sumfoo
    
    lscale <- lscale + log(sumfoo)
  }
  return(lbeta)
}


#Normal pseudo-residuals for Normal HMM
#Type can be "ordinary" or "forecast"
mvnorm.HMM.pseudo_residuals <- function(x, mod, k, type)
{
  if (type=="ordinary"){
    n <- length(x)
    la <- mvnorm.HMM.lforward(x, mod)
    lb <- mvnorm.HMM.lbackward(x, mod)
    lafact <- apply(la, 2, max)
    lbfact <- apply(lb, 2, max)
    
    P <- matrix(NA, n, mod$m)
    for (i in 1:n){
      P[i,] <- mvnorm.dist(x[[i]][1,], mod, mod$m, k)
    }
    
    npsr <- rep(NA, n)
    npsr[1] <- qnorm(rbind(c(1, 0, 0)) %*% P[1,])
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
    la <- mvnorm.HMM.lforward(x, mod)
    
    P <- matrix(NA, n, mod$m)
    for (i in 1:n){
      P[i,] <- mvnorm.dist(x[[i]][1,], mod, mod$m, k)
    }
    
    npsr <- rep(NA, n)
    npsr[1] <- qnorm(rbind(c(1, 0, 0)) %*% P[1,])
    for (i in 2:n) {
      la_max <- max(la[,i-1])
      a <- exp(la[,i-1]-la_max)
      npsr[i] <- qnorm(t(a)%*%(gamma/sum(a))%*%P[i,])
    }
    
    return(data.frame(npsr, x, index=c(1:n)))
  }
}

#Get multivariate normal distribution given mod and x 
mvnorm.dist <- function(x, mod, m, k){
  pvect <- numeric(m)
  for (i in 1:m){
    pvect[i] <- pmvnorm(lower=rep(-Inf, k), upper=x, mean=mod$mu[[i]], sigma=mod$sigma[[i]])
  }
  return(pvect)
}

#Jacobian matrix for parameters
#n should be total number of parameters estimated, excluding delta
mvnorm.jacobian <- function(m, k, n, parvect, stationary=TRUE){
  pn <- mvnorm.HMM.pw2pn(m, k, parvect, stationary)
  jacobian <- matrix(0, nrow=n, ncol=n)
  #Jacobian for mu only is a m*k identity matrix
  jacobian[1:(m*k), 1:(m*k)] <- diag(m*k)
  #count is the row at which the current parameter's derivatives are placed
  rowcount <- m*k + 1
  #There are t*m sigma parameters
  t <- triangular_num(k)
  for (i in 1:m){
    sigma <- pn$sigma[[i]]
    sigma[lower.tri(sigma, diag=FALSE)] <- rep(1, length(sigma[lower.tri(sigma, diag=FALSE)]))
    sigma <- sigma[lower.tri(sigma, diag=TRUE)]
    jacobian[rowcount:(rowcount+t-1), rowcount:(rowcount+t-1)] <- diag(sigma)
    rowcount <- rowcount + t
  }
  colcount <- rowcount
  for (i in 1:m){
    for (j in 1:m){
      if (j != i){
        foo <- -pn$gamma[i,j]*pn$gamma[i,]
        foo[j] <- pn$gamma[i,j]*(1-pn$gamma[i,j])
        foo <- foo[-i]
        jacobian[rowcount, colcount:(colcount+m-2)] <- foo
        rowcount <- rowcount + 1
      }
    }
    colcount <- colcount + m - 1
  }
  return(jacobian)
}

#Bootstrapping estimates
mvnorm.bootstrap.estimates <- function(mod, n, k, len, stationary){
  m <- mod$m
  mu_estimate <- numeric(n*m*k)
  sigma_estimate <- numeric(n*m*k*k)
  gamma_estimate <- numeric(n*m*m)
  delta_estimate <- numeric(n*m)
  for (i in 1:n){
    sample <- mvnorm.HMM.generate_sample(len, mod)
    mod2 <- mvnorm.HMM.mle(sample$obs, m, k, mod$mu, mod$sigma, mod$gamma, mod$delta, stationary=stationary)
    mu_estimate[((i-1)*m*k+1):(i*m*k)] <- unlist(mod2$mu, use.names = FALSE)
    sigma_estimate[((i-1)*m*k*k+1):(i*m*k*k)] <- unlist(mod2$sigma, use.names = FALSE)
    gamma_estimate[((i-1)*m*m+1):(i*m*m)] <- mod2$gamma
    delta_estimate[((i-1)*m+1):(i*m)] <- mod2$delta
  }
  return(list(mu=mu_estimate, sigma=sigma_estimate, gamma=gamma_estimate, delta=delta_estimate))
}

mvnorm.bootstrap.ci <- function(mod, bootstrap, alpha, m, k){
  #Rows of matrix are states, columns are variables
  mu_lower <- matrix(NA, m, k)
  mu_upper <- matrix(NA, m, k)
  bootstrap_mu <- data.frame(mu=bootstrap$mu)
  mu <- unlist(mvnorm.mle$mu, use.names = FALSE)
  for (i in 1:m){
    for (j in 1:k){
      if (i==m & j==k){
        foo <- bootstrap_mu %>% filter((row_number() %% (m*k)) == 0)
      }
      else{
        foo <- bootstrap_mu %>% filter((row_number() %% (m*k)) == (i-1)*k+j)
      }
      mu_lower[i, j] <- 2*mu[(i-1)*k+j] - quantile(foo$mu, 1-(alpha/2), names=FALSE)
      mu_upper[i, j] <- 2*mu[(i-1)*k+j] - quantile(foo$mu, alpha/2, names=FALSE)
    }
  }
  
  #Only want lower triangle of each sigma matrix, since is symmetric 
  t <- triangular_num(k)
  mat <- matrix(c(1:(k*k)), k)
  tvect <- mat[lower.tri(mat, diag=TRUE)]
  sigma_lower = matrix(NA, 3, t)
  sigma_upper = matrix(NA, 3, t)
  bootstrap_sigma<- data.frame(sigma=bootstrap$sigma)
  sigma <- unlist(mvnorm.mle$sigma, use.names = FALSE)
  for (i in 1:m){
    for (j in 1:t){
      tj <- tvect[j]
      if (i==m & j==t){
        foo <- bootstrap_sigma %>% filter((row_number() %% (m*k*k)) == 0)
      }
      else{
        foo <- bootstrap_sigma %>% filter((row_number() %% (m*k*k)) == (i-1)*k*k+tj)
      }
      sigma_lower[i, j] <- 2*sigma[(i-1)*k*k+tj] - quantile(foo$sigma, 1-(alpha/2), names=FALSE)
      sigma_upper[i, j] <- 2*sigma[(i-1)*k*k+tj] - quantile(foo$sigma, alpha/2, names=FALSE)
    }
  }
  
  gamma_lower <- rep(NA, m*m)
  gamma_upper <- rep(NA, m*m)
  bootstrap_gamma <- data.frame(gamma=bootstrap$gamma)
  gamma <- mod$gamma
  for (i in 1:(m*m)){
    if (i==(m*m)){
      foo <- bootstrap_gamma %>% filter((row_number() %% (m*m)) == 0)
    }
    else{
      foo <- bootstrap_gamma %>% filter((row_number() %% (m*m)) == i)
    }
    gamma_lower[i] <- 2*gamma[i] - quantile(foo$gamma, 1-(alpha/2), names=FALSE)
    gamma_upper[i] <- 2*gamma[i] - quantile(foo$gamma, alpha/2, names=FALSE)
  }
  
  delta_lower <- rep(NA, m)
  delta_upper <- rep(NA, m)
  bootstrap_delta <- data.frame(delta=bootstrap$delta)
  delta <- mod$delta
  for (i in 1:m){
    if (i==m){
      foo <- bootstrap_delta  %>% filter((row_number() %% m) == 0)
    }
    else{
      foo <- bootstrap_delta  %>% filter((row_number() %% m) == i)
    }
    delta_lower[i] <- 2*delta[i] - quantile(foo$delta, 1-(alpha/2), names=FALSE)
    delta_upper[i] <- 2*delta[i] - quantile(foo$delta, alpha/2, names=FALSE)
  }
  
  return(list(mu_lower=mu_lower, mu_upper=mu_upper, 
              sigma_lower=sigma_lower, sigma_upper=sigma_upper,
              gamma_lower=gamma_lower, gamma_upper=gamma_upper, 
              delta_lower=delta_lower, delta_upper=delta_upper))
}

#Fitting 3-state bivariate normal 
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
mu0 <- list(c(5, 5), c(5, 5), c(5, 5))
sigma0 <- list(diag(2, 2), diag(2, 2), diag(2, 2))
gamma0 <- matrix(rep(1/m, m*m), nrow = m, byrow = TRUE)
delta0 <- rep(1/m, m)
mvnorm.mle <- mvnorm.HMM.mle(mvnorm.sample$obs, m, k, mu0, sigma0, gamma0, delta0, stationary=FALSE, hessian=TRUE)
mvnorm.mle 
mvnorm.decoding <- mvnorm.HMM.viterbi(mvnorm.sample$obs, mvnorm.mle)
count(mvnorm.decoding$states)

#Get psuedo-residuals
mvnorm.pr <- mvnorm.HMM.pseudo_residuals(mvnorm.sample$obs, mvnorm.mle, k, "ordinary")
#Index plot of pseudo-residuals
ggplot(mvnorm.pr)+
  geom_point(aes(x=index, y=npsr), size=0.5, colour="black")+
  theme_minimal()
#Histogram of pseudo-residuals
ggplot(mvnorm.pr, aes(npsr))+
  geom_histogram(aes(y=..density..),  colour="navy", fill="light blue")+
  stat_function(fun=dnorm, colour="red")+
  theme_minimal()
#QQ plot of pseudo-residuals
ggplot(mvnorm.pr, aes(sample=npsr))+
  stat_qq()+
  stat_qq_line()+
  theme_minimal()
#ACF of pseudo-residuals
ggacf(mvnorm.pr$npsr)+
  theme_minimal()

#Confidence intervals using Hessian
sd <- sqrt(diag(mvnorm.mle$invhessian))
mu_estimate <- unlist(mvnorm.mle$mu, use.names = FALSE)
sigma_estimate <- mvnorm.mle$sigma
for (i in 1:m){
  foo <- sigma_estimate[[i]]
  sigma_estimate[[i]] <- foo[lower.tri(foo, diag=TRUE)]
  
}
sigma_estimate <- unlist(sigma_estimate, use.names = FALSE)
gamma_estimate <- as.vector(mvnorm.mle$gamma[!diag(m)])
theta <- c(mu_estimate, sigma_estimate, gamma_estimate)
mvnorm.cih <- list(lower=theta-1.96*sd, upper=theta+1.96*sd)
mvnorm.cih


#Bootstrapping confidence intervals
n <- 2
len <- length(mvnorm.sample$obs)
mvnorm.bootstrap <- mvnorm.bootstrap.estimates(mvnorm.mle, n, k, len, FALSE)
mvnorm.bootstrap 
mvnorm.cib <- mvnorm.bootstrap.ci(mvnorm.mle, mvnorm.bootstrap, 0.05, m, k)
mvnorm.cib
