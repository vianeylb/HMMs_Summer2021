library(ggpubr)
library(ggplot2)
library(forecast)
library(gridExtra)
library(plyr)
# Single Variable and Single Subject -------------------
#A.1.1 (modified for normal)
norm.HMM.pn2pw <- function(m, mu, sigma, gamma, delta=NULL, stationary=TRUE) {
  #' Transform natural parameters to working
  #'
  #' This function is for normal distributions.
  #' 
  #' m = number of states,
  #' mu = vector of means for each state dependent normal distribution
  #' sigma = vector of standard deviations for each state dependent normal distribution
  #' gamma = transition probability matrix
  #' delta = inital state distribution
  
  tmu <- mu 
  tsigma <- log(sigma)
  
  if(m==1) {
    return(tmu, tsigma)}
  
  foo <- log(gamma/diag(gamma)) 
  tgamma <- as.vector(foo[!diag(m)]) 
  
  if(stationary) {
    tdelta <- NULL}
  else {
    tdelta <- log(delta[-1]/delta[1])} 
  
  parvect <- c(tmu, tsigma, tgamma, tdelta) 
  return(parvect)
}


#A.1.2 (modified for normal)
norm.HMM.pw2pn <- function(m, parvect, stationary=TRUE) {
  #' Transform working parameters to natural
  #'
  #' This function is for normal distributions.
  #' 
  #' m = number of states,
  #' parvect = (working means, working sd, working trans prob matrix entries, working initial dist) 
  
  mu <- parvect[1:m]
  sigma <- exp(parvect[(m + 1):(2*m)]) 
  gamma <- diag(m) 
  
  if (m==1) {
    return(list(mu=mu, sigma=sigma, gamma=gamma, delta=1))}
  
  gamma[!gamma] <- exp(parvect[(2*m + 1):(2*m*m)]) 
  gamma <- gamma/apply(gamma, 1, sum) 
  
  if(stationary) {
    delta<-solve(t(diag(m)-gamma + 1),rep(1,m))}
  else {
    foo<-c(1,exp(parvect[(2*m*m + 1):(2*m*m + m-1)])) 
    delta <-foo/sum(foo)}
  
  return(list(mu=mu, sigma=sigma, gamma=gamma, delta=delta))
}



#A.1.3 (modified for normal)
norm.HMM.mllk <- function(parvect, x, m, stationary=TRUE, ...) {
  #' Compute -log-likelihood from working parameters
  #'
  #' This function is for normal distributions.
  #' 
  #' parvect = (working means, working sds, working trans prob matrix entries, working initial dist),
  #' x = observations,
  #' m = number of states,
  
  if(m==1) {return(-sum(dnorm(x, parvect[1], exp(parvect[2]), log=TRUE)))} 
  
  n       <- length(x) 
  pn      <- norm.HMM.pw2pn(m,parvect, stationary=stationary) 
  foo     <- pn$delta*dnorm(x[1], pn$mu, pn$sigma) 
  sumfoo  <- sum(foo) 
  lscale  <- log(sumfoo)
  foo     <- foo/sumfoo
  
  for (i in 2:n) {
    if (!is.na(x[i])) {P <- dnorm(x[i], pn$mu, pn$sigma)}  
    else {P <- rep(1,m)} 
    
    foo     <- foo %*% pn$gamma*P 
    sumfoo  <- sum(foo) 
    lscale  <- lscale+log(sumfoo) 
    foo     <- foo/sumfoo
  }
  mllk <- -lscale
  return(mllk)
}


#' #A.1.4 (modified for normal)
#' norm.HMM.mle <- function(x, m, mu0, sigma0, gamma0, delta0=NULL, stationary=TRUE,...) {
#'   #' Compute Maximum Likelihood Estimate
#'   #'
#'   #' This function is for normal distributions starting with natural parameters
#'   #' 
#'   #' x        = observations,
#'   #' m        = number of states,
#'   #' mu0      = inital guess for natural means
#'   #' sigma0   = initial guess for natural standard deviations
#'   #' gamma0   = initial guess for natural transition probability matrix
#'   #' delta0   = initial guess for initial state distribution
#'   
#'   parvect0 <- norm.HMM.pn2pw(m, mu0, sigma0, gamma0, delta0, stationary=stationary)
#'   mod <- nlm(norm.HMM.mllk, parvect0, x=x, m=m, stationary=stationary) 
#'   
#'   pn    <- norm.HMM.pw2pn(m=m, mod$estimate, stationary=TRUE) 
#'   mllk  <- mod$minimum 
#'   np    <- length(parvect0) 
#'   AIC   <- 2*(mllk+np) 
#'   n     <- sum(!is.na(x)) 
#'   BIC   <- 2*mllk+np*log(n) 
#'   
#'   list(m=m, 
#'        mu=pn$mu, 
#'        sigma=pn$sigma,
#'        gamma=pn$gamma, 
#'        delta=pn$delta, 
#'        code=mod$code, 
#'        mllk=mllk,
#'        AIC=AIC,
#'        BIC=BIC)
#' }


#A.1.4 (modified for normal)
norm.HMM.fit <- function(x, m, mu0, sigma0, gamma0, delta0=NULL, stationary=TRUE,...) {
  #' Compute Maximum Likelihood Estimate
  #'
  #' This function is for normal distributions starting with natural parameters
  #' 
  #' x        = observations,
  #' m        = number of states,
  #' mu0      = inital guess for natural means
  #' sigma0   = initial guess for natural standard deviations
  #' gamma0   = initial guess for natural transition probability matrix
  #' delta0   = initial guess for initial state distribution
  
  parvect0 <- norm.HMM.pn2pw(m, mu0, sigma0, gamma0, delta0, stationary=stationary)
  mod <- nlm(norm.HMM.mllk, parvect0, x=x, m=m, stationary=stationary) 
  
  pn    <- norm.HMM.pw2pn(m=m, mod$estimate, stationary=TRUE) 
  
  list(m=m, 
       mu=pn$mu, 
       sigma=pn$sigma,
       gamma=pn$gamma, 
       delta=pn$delta)
}


#A.1.5 (modified for normal)
norm.HMM.generate_sample <- function(ns, mod) {
  #' Generate a sample realization of an HMM
  #'
  #' This function is for normal distributions.
  #' 
  #' ns = length of realization 
  #' mod = HMM
  
  mvect <- 1:mod$m 
  state <- numeric(ns)
  state[1]<- sample(mvect, 1, prob=mod$delta)
  
  for (i in 2:ns) {
    state[i] <- sample(mvect, 1, prob=mod$gamma[state[i-1],])} 
  
  x <- rnorm(ns, mean=mod$mu[state], sd=mod$sigma[state]) 
  return(list(state = state, observ = x)) 
}


#A.1.6 (modified for normal)
norm.HMM.viterbi <-function(x, mod){
  #' Global decoding by the Viterbi algorithm
  #'
  #' This function is for normal distributions.
  #' 
  #' x = sequence of observations
  #' mod = HMM
  
  n       <- length(x) 
  xi      <- matrix(0,n,mod$m)  
  foo     <- mod$delta*dnorm(x[1], mod$mu, mod$sigma) 
  xi[1,]  <- foo/sum(foo) 
  
  for (i in 2:n){
    foo<-apply(xi[i-1,]*mod$gamma, 2,max)*dnorm(x[i], mod$mu, mod$sigma)
    xi[i,] <- foo/sum(foo) 
  }
  iv<-numeric(n) 
  iv[n] <-which.max(xi[n,])
  
  for (i in (n-1):1){
    iv[i] <- which.max(mod$gamma[,iv[i+1]]*xi[i,])
  }
  return(iv)
}


#A.1.7 Computing log(forward probabilities)
norm.HMM.lforward <-function(x,mod){
  n           <- length(x) #number of observations
  lalpha      <- matrix(NA,mod$m,n) #mxn matrix
  foo         <- mod$delta*dnorm(x[1],mod$mu, mod$sigma) #probability vector that observ1 came from each state
  sumfoo      <- sum(foo) #sum of prob vector
  lscale      <- log(sumfoo) #log of sum of prob vector
  foo         <- foo/sumfoo #divide entries of prop vector by log of sum
  lalpha[,1]  <- lscale+log(foo) #set first column of matrix
  
  for (i in 2:n){
    foo         <- foo%*%mod$gamma*dnorm(x[i], mod$mu, mod$sigma) 
    sumfoo      <- sum(foo)
    lscale      <- lscale+log(sumfoo)
    foo         <- foo/sumfoo
    lalpha[,i]  <- log(foo)+lscale
  }
  return(lalpha)
}


#A.1.8 Computing log(backward probabilities)
norm.HMM.lbackward <-function(x,mod){
  n           <- length(x) #number of observations
  m           <- mod$m #number of states
  lbeta       <- matrix(NA,m,n) #mxn matrix
  lbeta[,n]   <- rep(0,m) #fill last column with zeros
  foo         <- rep(1/m,m)
  lscale      <- log(m)
  
  for (i in (n-1):1){
    foo         <- mod$gamma%*%(dnorm(x[i+1], mod$mu, mod$sigma)*foo) 
    lbeta[,i]   <- log(foo)+lscale
    sumfoo      <- sum(foo)
    foo         <- foo/sumfoo
    lscale      <- lscale+log(sumfoo)
  }
  return(lbeta)
}


#A.1.9 Conditional probabilities
norm.HMM.conditional <- function(xc,x,mod){
  n     <- length(x)
  m     <- mod$m
  nxc   <- length(xc)
  dxc   <- matrix(NA,nrow=nxc,ncol=n)
  Px    <- matrix(NA,nrow=m,ncol=nxc)
  
  for (j in 1:nxc) {
    Px[,j] <-dnorm(xc[j], mod$mu, mod$sigma)}
  
  la      <- norm.HMM.lforward(x,mod)
  lb      <- norm.HMM.lbackward(x,mod)
  la      <- cbind(log(mod$delta),la) 
  lafact  <- apply(la,2,max) #get max of each column of la (max prob that observ arose from given state)
  lbfact  <- apply(lb,2,max)
  
  for (i in 1:n){
    foo     <- (exp(la[,i]-lafact[i])%*%mod$gamma)*exp(lb[,i]-lbfact[i])
    foo     <- foo/sum(foo)
    dxc[,i] <- foo%*%Px
  }
  return(dxc)
}


#A.1.10 Pseudo-residuals
norm.HMM.pseudo_residuals <- function(x,mod) {
  n <- length(x)
  xc <- sort(x)
  width <- diff(xc)
  cdists <- norm.HMM.conditional(xc,x,mod) 
  cdistno1 <- cdists[1:(n-1),]
  mult <- apply(cdistno1, 2, "*", width)
  mult <- rbind(mult, width[n-1]*cdists[n,])
  multsum <- apply(mult, 2,cumsum)
  
  for (i in 1:(n-1)){
    if (max(multsum[,i]) > 1){
      multsum[,i] <- multsum[,i]/max(multsum[,i])
    }
  }
  df <- rep(0,n)
  for (i in 1:n){
    df[i] = multsum[which(xc==x[i]),i]
  }
  npr <- qnorm(df)
  return(npr)
}


#Get standard error of all fitted parameters using parametric bootstrap method
norm.HMM.params_SE <- function(x, n, modfit, stationary=TRUE){
  ns <- length(x)
  m <- modfit$m
  mus = matrix(numeric(m*n), nrow = n, ncol=m) #matrix to be filled with fitted mus
  sigmas = matrix(numeric(m*n), nrow = n, ncol=m) #matrix to be filled with fitted mus
  gammas = matrix(numeric(m*m*n), nrow = n, ncol = m*m) #matrix to be filled with entries of gamma
  deltas = matrix(numeric(m*n), nrow = n, ncol=m) #matrix to be filled with fitted deltas
  
  for (i in 1:n){
    sample <- norm.HMM.generate_sample(ns, modfit) #generate observations based on modfit
    x <- sample$observ #get observations from generated sample
    mod <- norm.HMM.fit(x,m,modfit$mu,modfit$sigma, modfit$gamma, modfit$delta, stationary=stationary) #fit model to generated observations
    mus[i,] <- mod$mu #add fitted mu to mus matrix
    sigmas[i,] <- mod$sigma
    gammas[i,] <- as.vector(t(mod$gamma)) #add fitted gamma as row
    deltas[i,] = mod$delta #add fitted delta to deltas matrix
  }
  
  mu.cov = cov(mus) #get var-covar matrix of mus
  mu.SE = sqrt(diag(mu.cov))
  #mu.upper = modfit$mu + (1.96 * mu.SE) #calculate upper 95% CI from var
  #mu.lower = modfit$mu - (1.96 * mu.SE) #calculate lower 95% CI from var
  
  sigma.cov = cov(sigmas) #get var-covar matrix of sigmas
  sigma.SE = sqrt(diag(sigma.cov))
  #sigma.upper = modfit$sigma + (1.96 * sigma.SE) #calculate upper 95% CI from var
  #sigma.lower = pmax(modfit$sigma - (1.96 * sigma.SE),0) #calculate lower 95% CI from var
  
  delta.cov = cov(deltas) #get var-covar matrix of lambdas
  delta.SE = sqrt(diag(delta.cov))
  #delta.upper = modfit$delta + (1.96 * delta.SE) #calculate upper 95% CI from var
  #delta.lower = pmax(modfit$delta - (1.96 * delta.SE),0) #calculate lower 95% CI from var
  
  gammafit = as.vector(t(modfit$gamma))
  gamma.cov = cov(gammas)
  gamma.SE = sqrt(diag(gamma.cov))
  gamma.SE = matrix(gamma.SE, m,m, byrow=TRUE)
  #gamma.upper = gammafit + (1.96 * sqrt(diag(gamma.cov))) #calculate upper 95% CI from var
  #gamma.upper = matrix(gamma.upper, m,m, byrow=TRUE)
  #gamma.lower = pmax(gammafit - (1.96 * sqrt(diag(gamma.cov))),0) #calculate lower 95% CI from var
  #gamma.lower = matrix(gamma.lower, m,m, byrow=TRUE)
  
  result = list("mu" = modfit$mu, 
                "mu.SE" = mu.SE,
                #"mu.upper.conf" = mu.upper,
                #"mu.lower.conf" = mu.lower,
                "sigma" = modfit$sigma, 
                "sigma.SE" = sigma.SE,
                #"sigma.upper.conf" = sigma.upper,
                #"sigma.lower.conf" = sigma.lower,
                "gamma" = modfit$gamma,
                "gamma.SE" = gamma.SE,
                #"gamma.upper.conf" = gamma.upper,
                #"gamma.lower.conf" = gamma.lower,
                "delta" = modfit$delta, 
                "delta.SE" = delta.SE)
                #"delta.upper.conf" = delta.upper,
                #"delta.lower.conf" = delta.lower
  return(result)
}


#Compute CIs of mu and sigma using Monte Carlo approach
norm.HMM.CI_MonteCarlo <- function(range, m, n=100, params_SE, level=0.975){
  xc = length(range)
  mu = params_SE$mu
  mu.SE = params_SE$mu.SE
  sigma = params_SE$sigma
  sigma.SE = params_SE$sigma.SE
  
  density.lst <- list(matrix(numeric(xc*n), ncol = xc, nrow = n))
  
  for (k in 1:m){
    densities <- matrix(numeric(xc*n), ncol = xc, nrow = n)
    for (i in 1:n){
      sample.mu <- rnorm(1, mu[k], mu.SE[k])
      sample.sigma <- rnorm(1, sigma[k], sigma.SE[k])
      densities[i,] <- dnorm(range, sample.mu, sample.sigma)
    }
    density.lst[[k]] <- densities
  }
  
  upper <- matrix(numeric(xc*m), ncol=xc, nrow=m)
  lower <- matrix(numeric(xc*m), ncol=xc, nrow=m)

  for (k in 1:m){
    densities <- density.lst[[k]]
    for (j in 1:xc){
      upper[k,j] <- quantile(densities[,j], probs = level, na.rm= TRUE)
      lower[k,j] <- quantile(densities[,j], probs = 1-level, na.rm= TRUE)
    }
  }
  return(list(range=range, upper=upper, lower=lower))
}

#MLE that also conputed the inverse hessian (variance-covariance matrix)
# norm.HMM.mle <- function(x, m, mu0, sigma0, gamma0, delta0=NULL, stationary=TRUE, hessian=TRUE,...){
#   parvect0 <- norm.HMM.pn2pw(m, mu0, sigma0, gamma0, delta0, stationary=stationary)
#   mod <- nlm(norm.HMM.mllk, parvect0, x=x,m=m, stationary=stationary, hessian=hessian)
#   pn <- norm.HMM.pw2pn(m=m, mod$estimate, stationary=stationary)
#   mllk <- mod$minimum
#   if (hessian){
#     h <- mod$hessian
#     if (det(h) != 0){
#       h <- solve(h)
#       jacobian <- norm.jacobian(m, mod$estimate, stationary=stationary)
#       h <- t(jacobian)%*%h%*%jacobian
#     }
#   }
#   np <- length(parvect0)
#   AIC <- 2*(mllk+np)
#   n <- sum(!is.na(x))
#   BIC <- 2*mllk+np*log(n)
#   if (hessian){
#     return(list(m=m, mu=pn$mu, sigma=pn$sigma, gamma=pn$gamma, 
#                 delta=pn$delta, code=mod$code, mllk=mllk, 
#                 AIC=AIC, BIC=BIC, hessian = mod$hessian, invhessian = h))
#   }
#   else{
#     return(list(m=m, mu=pn$mu, sigma=pn$sigma, gamma=pn$gamma, 
#                 delta=pn$delta, code=mod$code, mllk=mllk, 
#                 AIC=AIC, BIC=BIC))
#   }
# }


#compute the jacobian for the conversion of the hessian
# norm.jacobian <- function(m, parvect, stationary=TRUE){
#   n <- 2*m+m*(m-1)
#   pn <- norm.HMM.pw2pn(m=m, parvect=parvect, stationary=stationary)
#   jacobian <- matrix(0, nrow=n, ncol=n)
#   jacobian[1:m, 1:m] <- diag(m)
#   jacobian[(m+1):(2*m), (m+1):(2*m)] <- diag(pn$sigma)
#   mat <- pn$gamma
#   diag(mat) <- NA
#   mat <-matrix(mat[which(!is.na(mat))],nrow=m,ncol=m-1, byrow=TRUE)
#   for (k in 1:m){
#     M <- matrix(numeric((m-1)*(m-1)),nrow=m-1,ncol=m-1)
#     for (i in 1:(m-1)){
#       for (j in 1:(m-1)){
#         if (i==j){
#           M[i,j] <- mat[k,i]*(1 - mat[k,i])}
#         else {
#           M[i,j] <- -mat[k,i]*mat[k,j]}
#       }
#     }
#     jacobian[(2*m+1+(k-1)*(m-1)):(2*m+k*(m-1)),(2*m+1+(k-1)*(m-1)):(2*m+k*(m-1))] <- M
#   }
#   return(jacobian)
# }


#gets the standard errors from the inverse hessian
norm.HMM.SE_hessian <- function(m, mod){
  h <- sqrt(diag(mod$invhessian))
  mu.SE <- h[1:m]
  sigma.SE <- h[(m+1):(2*m)]
  gamma.SE <- diag(m)
  gamma.SE[!gamma.SE] <- h[(2*m+1):(2*m+m*(m-1))]
  diag(gamma.SE) <- NA
  return(list(mu=mod$mu, mu.SE=mu.SE, sigma=mod$sigma, sigma.SE=sigma.SE,
              gamma=mod$gamma, gamma.SE=gamma.SE))
}


# Multivariate, Single Subject ---------------------
#Transform natural parameters to working
mvnorm.HMM.pn2pw <- function(m, mu, sigma, gamma, delta=NULL, stationary=TRUE){
  tmu <- unlist(mu)
  tsigma <- log(unlist(sigma))
  
  if(m==1) {
    return(tmu, tsigma)}
  
  foo <- log(gamma/diag(gamma)) 
  tgamma <- as.vector(foo[!diag(m)]) 
  
  if(stationary) {
    tdelta <- NULL}
  else {
    tdelta <- log(delta[-1]/delta[1])} 
  
  parvect <- c(tmu, tsigma, tgamma, tdelta) 
  return(parvect)
}


#Transform working parameters to natural
#v is the number of response variables
mvnorm.HMM.pw2pn <- function(m, v, parvect, stationary=TRUE){
  mu <- split(parvect[1:(m*v)], ceiling(seq_along(parvect[1:(m*v)])/m))
  sigma <- split(exp(parvect[(m*v+1):(2*m*v)]), ceiling(seq_along(exp(parvect[(m*v+1):(2*v*m)]))/m))
  gamma <- diag(m) 
  
  if (m==1) {
    return(list(mu=mu, sigma=sigma, gamma=gamma, delta=1))}
  
  gamma[!gamma] <- exp(parvect[(2*m*v+1):(2*m*v + m*m-m)]) 
  gamma <- gamma/apply(gamma, 1, sum) 
  
  if(stationary){
    delta<-solve(t(diag(m)-gamma+1),rep(1,m))}
  else {
    foo<-c(1,exp(parvect[(2*m*v + m*m-m+1):(length(parvect)-1)])) 
    delta <-foo/sum(foo)}
  
  return(list(mu=mu, sigma=sigma, gamma=gamma, delta=delta))
}


#Compute -log-likelihood from working parameters
mvnorm.HMM.mllk <- function(parvect, x, m, v, stationary=TRUE, ...){
  n       <- nrow(x)
  pn      <- mvnorm.HMM.pw2pn(m, v, parvect, stationary=stationary) 
  P <- rep(1,m)
  for (k in 1:v){
    P <- P*dnorm(x[1,k], pn$mu[[k]], pn$sigma[[k]])}
  foo     <- pn$delta*P
  sumfoo  <- sum(foo) 
  lscale  <- log(sumfoo)
  foo     <- foo/sumfoo
  
  for (i in 2:n){
    P <- rep(1,m)
    for (k in 1:v){
      if (!is.na(x[i,k])){
        P <- P*dnorm(x[i,k], pn$mu[[k]], pn$sigma[[k]])}
    }
    foo     <- foo %*% pn$gamma*P 
    sumfoo  <- sum(foo) 
    lscale  <- lscale+log(sumfoo) 
    foo     <- foo/sumfoo
  }
  mllk <- -lscale
  return(mllk)
}


#Compute Maximum Likelihood Estimate
mvnorm.HMM.mle <- function(x, m, v, mu0, sigma0, gamma0, delta0=NULL, stationary=TRUE,...) {
  
  parvect0 <- mvnorm.HMM.pn2pw(m, mu0, sigma0, gamma0, delta0, stationary=stationary)
  mod <- nlm(mvnorm.HMM.mllk, parvect0, x=x, m=m, v=v, stationary=stationary) 
  
  pn    <- mvnorm.HMM.pw2pn(m=m, v=v, mod$estimate, stationary=TRUE) 
  mllk  <- mod$minimum 
  #np    <- length(parvect0) 
  #AIC   <- 2*(mllk+np) 
  #n     <- sum(!is.na(x)) 
  #BIC   <- 2*mllk+np*log(n) 
  
  list(m=m, 
       v=v,
       mu=pn$mu, 
       sigma=pn$sigma,
       gamma=pn$gamma, 
       delta=pn$delta, 
       code=mod$code, 
       mllk=mllk)
}


#Generate a sample realization of an HMM
mvnorm.HMM.generate_sample <- function(ns, mod){
  mvect <- 1:mod$m 
  v <- mod$v
  state <- numeric(ns)
  state[1] <- sample(mvect, 1, prob=mod$delta)
  
  for (i in 2:ns) {
    state[i] <- sample(mvect, 1, prob=mod$gamma[state[i-1],])
  } 
  x <- matrix(numeric(ns*v), ncol=v, nrow=ns)
  for (k in 1:v){
    x[,k] <- rnorm(ns, mean=mod$mu[[k]][state], sd=mod$sigma[[k]][state]) 
  }
  return(list(state = state, observ = x)) 
}

# Multivariate, Multiple Subjects (Longitudinal) --------------------
#example observation data for mvnorm multiple subjects (3D array)
# T=10
# v=2
# s=3
# Time = c('t1','t2','t3','t4','t5','t6','t7','t8', 't9', 't10')
# Variable = c('Var1', 'Var2')
# Subject = c('Sub1', 'Sub2', 'Sub3')
# x <- array(dim=c(T, v, s), dimnames = list(Time, Variable, Subject))

#mu and sigma -> list of matrices each matrix is a different variable, row=subject, column=state
mvnorm.longitudinal.HMM.pn2pw <- function(m, s, mu, sigma, gamma, delta=NULL, stationary=TRUE){
  tmu <- numeric()
  tsigma <- numeric()
  for (i in 1:s){
    tmu <- c(tmu, as.vector(t(mu[[i]]))) #gives in order of variable, subject, state
    tsigma <- c(tsigma, log(as.vector(t(sigma[[i]]))))
  }
  if(m==1) {
    return(tmu, tsigma)}
  
  tgamma <- numeric()
  for (i in 1:s){
    foo <- log(gamma[[i]]/diag(gamma[[i]])) 
    tgamma <- c(tgamma, as.vector(foo[!diag(m)]))
  }
  if(stationary) {
    tdelta <- NULL}
  else {
    tdelta <- numeric()
    for (i in 1:s){
      tdelta <- c(tdelta,log(delta[[i]][-1]/delta[[i]][1]))}}
  parvect <- c(tmu, tsigma, tgamma, tdelta) 
  return(parvect)
}

# m   <-3 
# v   <- 2
# s   <- 2
# mu <- list(matrix(c(11, 19, 23, 
#                      12, 21, 26), ncol=m, nrow=s, byrow=TRUE), 
#             matrix(c(-10, -5, 1,
#                      -9, -5, 0), ncol=m, nrow=s, byrow=TRUE))
# sigma <- list(matrix(c(3, 3, 3, 
#                         3, 3, 3), ncol=m, nrow=s, byrow=TRUE), 
#                matrix(c(3, 3, 3, 
#                         3, 3, 3), ncol=m, nrow=s, byrow=TRUE))
# gamma  <- list(matrix(c(0.9, 0.05, 0.05,
#                          0.05, 0.9, 0.05,
#                          0.05, 0.05, 0.9),m,m,byrow=TRUE),
#                 matrix(c(0.9, 0.05, 0.05,
#                          0.05, 0.9, 0.05,
#                          0.05, 0.05, 0.9),m,m,byrow=TRUE))
# delta <- list(c(1/3,1/3,1/3), c(1/3,1/3,1/3))



#Transform working parameters to natural
mvnorm.longitudinal.HMM.pw2pn <- function(m, v, s, parvect, stationary=TRUE){
  mu <- split(parvect[1:(m*v*s)], ceiling(seq_along(parvect[1:(m*v*s)])/(m*s)))
  sigma <- split(exp(parvect[(m*v*s+1):(2*m*v*s)]), ceiling(seq_along(exp(parvect[(m*v*s+1):(2*v*m*s)]))/(m*s)))
  for (i in 1:v){
    mu[[i]] <- matrix(mu[[i]],ncol=m, nrow=s, byrow=TRUE)
    sigma[[i]] <- matrix(sigma[[i]],ncol=m, nrow=s, byrow=TRUE)
  }
  gamma <- list()
  delta <- list()
  for (i in 1:s){
    gamma[[i]] <- diag(m)
  }
  if (m==1) {
    for (i in 1:s){
      delta[[i]] = 1
    }
    return(list(mu=mu, sigma=sigma, gamma=gamma, delta=delta))}
  g <- split(parvect[(2*m*v*s+1):(2*m*v*s + s*m*m-s*m)], 
             ceiling(seq_along(parvect[(2*m*v*s+1):(2*m*v*s + s*m*m-s*m)])/(m*(m-1))))
  for (i in 1:s){
    gamma[[i]][!gamma[[i]]] <- exp(g[[i]]) 
    gamma[[i]] <- gamma[[i]]/apply(gamma[[i]], 1, sum) 
  }
  if(stationary){
    for (i in 1:s){
    delta[[i]]<-solve(t(diag(m)-gamma[[i]]+1),rep(1,m))}}
  else {
    d <- split(parvect[(2*m*v*s + s*m*m-s*m+1):length(parvect)], 
               ceiling(seq_along(parvect[(2*m*v*s + s*m*m-s*m+1):length(parvect)])/(m-1)))
    for (i in 1:s){
    foo<-c(1,exp(d[[i]])) 
    delta[[i]] <-foo/sum(foo)}}
  
  return(list(mu=mu, sigma=sigma, gamma=gamma, delta=delta))
}


#Compute -log-likelihood from working parameters
mvnorm.longitudinal.HMM.mllk <- function(parvect, x, m, v, s, stationary=TRUE, ...){
  n       <- nrow(x)
  pn      <- mvnorm.longitudinal.HMM.pw2pn(m, v, s, parvect, stationary=stationary) 
  
  longlscale <- 1
  for (j in 1:s){
    P <- rep(1,m)
    for (k in 1:v){
      P <- P*dnorm(x[1,k,j], pn$mu[[k]][j,], pn$sigma[[k]][j,])}
    foo     <- pn$delta[[j]]*P
    sumfoo  <- sum(foo) 
    lscale  <- log(sumfoo)
    foo     <- foo/sumfoo
    
    for (i in 2:n){
      P <- rep(1,m)
      for (k in 1:v){
        if (!is.na(x[i,k,j])){
          P <- P*dnorm(x[i,k,j], pn$mu[[k]][j,], pn$sigma[[k]][j,])}
      }
      foo     <- foo %*% pn$gamma[[j]]*P 
      sumfoo  <- sum(foo) 
      lscale  <- lscale+log(sumfoo) 
      foo     <- foo/sumfoo
    }
    longlscale <- longlscale + (-lscale) #adding because using loglikelihood
  }
  mllk <- longlscale
  return(mllk)
}


#Compute Maximum Likelihood Estimate
mvnorm.longitudinal.HMM.mle <- function(x, m, v, s, mu0, sigma0, gamma0, delta0=NULL, stationary=TRUE,...) {
  
  parvect0 <- mvnorm.longitudinal.HMM.pn2pw(m, s, mu0, sigma0, gamma0, delta0, stationary=stationary)
  mod <- nlm(mvnorm.longitudinal.HMM.mllk, parvect0, x=x, m=m, v=v, s=s, stationary=stationary) 
  
  pn    <- mvnorm.longitudinal.HMM.pw2pn(m=m, v=v, s=s, mod$estimate, stationary=TRUE) 
  mllk  <- mod$minimum 
  #np    <- length(parvect0) 
  #AIC   <- 2*(mllk+np) 
  #n     <- sum(!is.na(x)) 
  #BIC   <- 2*mllk+np*log(n) 
  
  list(m=m, 
       v=v,
       mu=pn$mu, 
       sigma=pn$sigma,
       gamma=pn$gamma, 
       delta=pn$delta, 
       code=mod$code, 
       mllk=mllk)
}


#Generate a sample realization of an HMM
mvnorm.longitudinal.HMM.generate_sample <- function(ns, mod){
  mvect <- 1:mod$m 
  v <- mod$v
  s <- mod$s
  state <- matrix(numeric(ns*s), nrow=ns, ncol=s)
  for (j in 1:s){
    state[1,j] <- sample(mvect, 1, prob=mod$delta[[j]])
    for (i in 2:ns) {
      state[i,j] <- sample(mvect, 1, prob=mod$gamma[[j]][state[(i-1),j],])
    } }
  x <- array(dim=c(ns,v,s))
  for (k in 1:v){
    for (j in 1:s){
      x[,k,j] <- rnorm(ns, mean=mod$mu[[k]][j,state[,j]], sd=mod$sigma[[k]][j,state[,j]]) 
    }
  }
  return(list(state = state, observ = x)) 
}


# HMM Multivaraite, Multisubject, pooling ------------------------
# example observation data for mvnorm multiple subjects (3D array)
# T=10
# v=2
# s=3
# Time = c('t1', 't2', 't3', 't4', 't5', 't6', 't7', 't8',  't9',  't10')
# Variable = c('Var1', 'Var2')
# Subject = c('Sub1', 'Sub2', 'Sub3')
# x <- array(dim = c(T, v, s), dimnames = list(Time, Variable, Subject))

norm_working_params <- function(num_states, num_subjects, 
                                mu, sigma, gamma, delta = NULL, 
                                stationary = TRUE, gamma_pooled = FALSE) {
  #' Transform normal HMM parameters from natural to working
  #' 
  #' This function transforms the natural normal HMM parameters that have 
  #' additional constraints into working parameters that incorporate the 
  #' constraints. The output is a single vector which includes all the working
  #' parameters.
  #' 
  #' @param num_states The number of states in the HMM.
  #' @param num_subjects The number of subjects that generated the data being
  #'   fit with the HMM.
  #' @param mu The means for the normally distributed state dependent 
  #'   distributions of the HMM. mu is a list of matrices with each matrix 
  #'   corresponding to a different variable in the data being fit. The columns
  #'   of the matrices correspond to the state number and the rows correspond to
  #'   the subject number.
  #' @param sigma The standard deviations for the normally distributed state 
  #'   dependent distributions of the HMM. sigma is a list of matrices with 
  #'   each matrix corresponding to a different variable in the data being fit. 
  #'   The columns of the matrices correspond to the state number and the rows
  #'   correspond to the subject number.
  #' @param gamma A list with each element being the transition probability 
  #'   matrix of the HMM for the subject corresponding to that index.
  #' @param delta A list with each element being the initial state distribution
  #'   vector of the HMM for the subject corresponding to that index.
  #' @param stationary A logical variable denoting whether to treat the HMM as 
  #'   one with a stationary distribution or without.
  #' @param gamma_pooled A logical variable indiacting whether the transition
  #'   probability matrix `gamma` and the resulting stationary distribution 
  #'   `delta` should be the same for all subjects.
  
  if (gamma_pooled) {
    ng <- 1
  } else {
    ng <- num_subjects
  }
  tmu     <- numeric()
  tsigma  <- numeric()
  for (j in 1:num_variables) {
    tmu    <- c(tmu, as.vector(t(mu[[j]]))) 
    tsigma <- c(tsigma, log(as.vector(t(sigma[[j]]))))
  }
  if (num_states == 1) {
    return(tmu, tsigma)
  }
  tgamma   <- numeric()
  for (i in 1:ng) {
      foo    <- log(gamma[[i]]/diag(gamma[[i]])) 
      tgamma <- c(tgamma, as.vector(foo[!diag(num_states)]))
  }
  if (stationary) {
    tdelta <- NULL
  } else {
    tdelta <- numeric()
    for (i in 1:ng) {
      tdelta <- c(tdelta, log(delta[[i]][-1]/delta[[i]][1]))
    }
  } 
  c(tmu, tsigma, tgamma, tdelta) 
}

split_vec <- function(vector, start, end, length, exp = FALSE) {
  #' Split vector into list
  #' 
  #' This function splits a vector from the index start to the index end into
  #' segments of length length and complies them into a list. This function can
  #' also apply the exponential function to the vector elements if desired.
  #' 
  #' @param vector The vector to be split into a list.
  #' @param start The first index of the vector vector that will be included in 
  #'   the list.
  #' @param end The last index of the vector vector that will be included in the
  #'   list.
  #' @param length A value indicating the length of the segments that vector
  #'   will be split into.
  #' @param exp A logical variable indicating whether the elements of vector
  #'   should be transformed with `exp()`.
  
  if (exp) {
    return(split(exp(vector[start:end]), 
                 ceiling(seq_along(exp(vector[start:end]))/length)))
  }
  split(vector[start:end], ceiling(seq_along(vector[start:end])/length))
}

norm_natural_params <- function(num_states, num_variables, num_subjects, 
                                working_params, stationary = TRUE,
                                state_dep_dist_pooled = FALSE,
                                gamma_pooled = FALSE) {
  #' Transform normal HMM parameters from working to natural
  #' 
  #' This function transforms the working normal HMM parameters back into the 
  #' original format of the natural parameters and outputs them as a list. This 
  #' function is the reverse of `norm_working_params()`.
  #' 
  #' @inheritParams norm_working_params
  #' @param num_variables The number of independent variables in the data being
  #'   fit with an HMM.
  #' @param working_params A vector of the working normal parameters for the 
  #'   HMM.
  #' @param state_dep_dist_pooled A logical variable indiacting whether the
  #'   state dependent distribution parameters `mu` and `sigma` should be 
  #'   treated as equal for all subjects.
  
  ns          <- num_subjects
  mu_start    <- 1
  mu_end      <- num_states*num_variables*num_subjects
  sigma_end   <- mu_end + num_states*num_variables*num_subjects
  mu_len      <- sigma_len <- num_subjects*num_states
  if (state_dep_dist_pooled) {
    mu_end    <- num_states*num_variables
    sigma_end <- mu_end + num_states*num_variables
    mu_len    <- sigma_len <- num_states
  }
  sigma_start <- mu_end + 1
  gamma_start <- sigma_end + 1
  gamma_end   <- sigma_end + num_subjects*num_states^2 - num_subjects*num_states
  if (gamma_pooled) {
    gamma_end <- sigma_end + num_states^2 - num_states
    ns        <- 1
  }
  delta_start <- gamma_end + 1
  delta_end   <- length(working_params)
  
  mu    <- split_vec(working_params, mu_start, mu_end, mu_len)       
  sigma <- split_vec(working_params, sigma_start, sigma_end, sigma_len, 
                     exp = TRUE)
  for (j in 1:num_variables) {
    mu[[j]]    <- matrix(mu[[j]], ncol = num_states, byrow = TRUE)
    sigma[[j]] <- matrix(sigma[[j]], ncol = num_states, byrow = TRUE)
  }
  gamma <- list()
  delta <- list()
  for (i in 1:ns) {
    gamma[[i]] <- diag(num_states)
  }
  if (num_states == 1) {
    for (i in 1:ns) {
      delta[[i]] = 1
    } 
    return(list(mu = mu, sigma = sigma, gamma = gamma, delta = delta))
  }
  g <- split_vec(working_params, gamma_start, gamma_end, 
                 num_states*(num_states - 1))
  for (i in 1:ns) {
    gamma[[i]][!gamma[[i]]] <- exp(g[[i]]) 
    gamma[[i]]              <- gamma[[i]]/apply(gamma[[i]], 1, sum) 
  }
  if (stationary) {
    for (i in 1:ns) {
      delta[[i]] <- solve(t(diag(num_states) - gamma[[i]] + 1), 
                          rep(1, num_states))
    }
  } else {
    d <- split_vec(working_params, delta_start, delta_end, num_states - 1)
    for (i in 1:ns) {
      foo        <- c(1, exp(d[[i]])) 
      delta[[i]] <- foo/sum(foo)
    }
  } 
  list(mu = mu, sigma = sigma, gamma = gamma, delta = delta)
}

norm_loglikelihood <- function(working_params, x, 
                               num_states, num_variables, num_subjects, 
                               stationary = TRUE, 
                               state_dep_dist_pooled = FALSE,
                               gamma_pooled = FALSE) {
  #' Compute negative log-likelihood of normal HMM parameters
  #' 
  #' This function computes the negative log-likelihood that the given normal 
  #' HMM parameters could have generated the data being fit.
  #' 
  #' @inheritParams norm_natural_params
  #' @param x The data to be fit with an HMM in the form of a 3D array. The 
  #'   first index (row) corresponds to time, the second (column) to the 
  #'   variable number, and the third (matrix number) to the subject number.
  
  num_time  <- nrow(x)
  pn        <- norm_natural_params(num_states, num_variables, num_subjects, 
                                   working_params, stationary = TRUE) 
  cum_loglikelihood <- 0
  for (i in 1:num_subjects) {
    s_ind   <- i
    g_ind   <- i
    if (state_dep_dist_pooled) {
      s_ind <- 1
    }
    if (gamma_pooled) {
      g_ind <- 1
    }
    P   <- rep(1, num_states)
    for (j in 1:num_variables) {
      P <- P*dnorm(x[1, j, i], pn$mu[[j]][s_ind, ], pn$sigma[[j]][s_ind, ])
    }
    forward_probs     <- pn$delta[[g_ind]]*P
    sum_forward_probs <- sum(forward_probs) 
    loglikelihood     <- log(sum_forward_probs)
    forward_probs     <- forward_probs/sum_forward_probs
    
    for (t in 2:num_time) {
      P     <- rep(1, num_states)
      for (j in 1:num_variables) {
        if (!is.na(x[t, j, i])) {
          P <- P*dnorm(x[t, j, i], pn$mu[[j]][s_ind, ], pn$sigma[[j]][s_ind, ])
        }
      }
      forward_probs     <- forward_probs %*% pn$gamma[[g_ind]]*P 
      sum_forward_probs <- sum(forward_probs) 
      loglikelihood     <- loglikelihood + log(sum_forward_probs) 
      forward_probs     <- forward_probs/sum_forward_probs
    }
    cum_loglikelihood <- cum_loglikelihood + loglikelihood
  }
  - cum_loglikelihood
}

norm_jacobian <- function(num_states, num_variables, num_subjects, 
                          working_params, natural_params,
                          state_dep_dist_pooled = FALSE,
                          gamma_pooled = FALSE) {
  #' Compute modified jacobian of natural parameters with respect to working
  #' 
  #' The function computes a modified version of the jacobian of the normal 
  #' natural HMM parameters with respect ot the normal working HMM parameters.
  #' The resulting matrix is useful for converting the hessian of the working
  #' parameters into the hessian of the natural parameters.
  #' 
  #' @inheritParams norm_natural_params
  #' @param natural_params A list of the normal natural HMM parameters `mu`, 
  #'   `sigma`, `gamma`, and `delta`.
  
  mu_len   <- num_subjects*num_variables*num_states
  ns       <- ng <- num_subjects
  if (state_dep_dist_pooled) {
    mu_len <- num_states*num_variables
    ns     <- 1
  }
  gamma_len   <- (num_states^2 - num_states)*num_subjects
  if (gamma_pooled) {
    gamma_len <- num_states^2 - num_states
    ng        <- 1
  }
  n <- 2*mu_len + gamma_len
  l <- num_states^2 - num_states
  
  sigma  <- numeric()
  for (i in 1:num_variables) {
    sigma <- c(sigma, as.vector(t(natural_params$sigma[[i]])))
  }
  jacobian <- matrix(0, nrow = n, ncol = n)
  jacobian[1:mu_len, 1:mu_len] <- diag(mu_len)
  jacobian[(mu_len + 1):(2*mu_len), (mu_len + 1):(2*mu_len)] <- diag(sigma)
  
  for (s in 1:ng) {
    gamma <- natural_params$gamma[[s]]
    gamma_ordered <- gamma
    diag(gamma_ordered) <- NA
    gamma_ordered <- gamma_ordered[which(!is.na(gamma_ordered))]
    w <- length(gamma_ordered)
    for (i in 1:length(gamma_ordered)) {
      for (j in 1:length(gamma_ordered)) {
        if (i == j) {
          jacobian[2*mu_len + (s - 1)*w + i, 2*mu_len + (s - 1)*w + i] <- 
            gamma_ordered[i]*(1 - gamma_ordered[j])
        } else if (which(gamma == gamma_ordered[i], arr.ind = TRUE)[1] == 
                   which(gamma == gamma_ordered[j], arr.ind = TRUE)[1]) {
          jacobian[2*mu_len + (s - 1)*w + i, 2*mu_len + (s - 1)*w + j] <- 
            -gamma_ordered[i]*gamma_ordered[j]
        }
      } 
    }
  } 
  jacobian
}

norm_fit_hmm <- function(x, num_states, num_variables, num_subjects,
                         mu0, sigma0, gamma0, delta0 = NULL, 
                         stationary = TRUE, 
                         state_dep_dist_pooled = FALSE, gamma_pooled = FALSE,
                         hessian = FALSE) {
  #' Fit an HMM
  #' 
  #' This function fits data with an HMM by maximizing the likelihood estimate
  #' given initial normal parameters. 
  #' 
  #' @inheritParams norm_loglikelihood
  #' @param mu0 The starting values for the means of the normally distributed 
  #'   state dependent distributions of the HMM. mu0 is a list of matrices with
  #'   each matrix corresponding to a different variable in the data being fit.
  #'   The columns of the matrices correspond to the state number and the rows
  #'   correspond to the subject number.
  #' @param sigma0 The starting values for the standard deviations of the 
  #'   normally distributed state dependent distributions of the HMM. sigma0 is 
  #'   a list of matrices with each matrix corresponding to a different variable
  #'   in the data being fit. The columns of the matrices correspond to the 
  #'   state number and the rows correspond to the subject number.
  #' @param gamma0 A list with each element being the starting transition 
  #'   probability matrix of the HMM for the subject corresponding to that 
  #'   index.
  #' @param delta0 A list with each element being the starting initial state 
  #'   distribution vector of the HMM for the subject corresponding to that 
  #'   index.
  #' @param hessian A logical variable denoting whether to compute the hessian 
  #'   at the minimum or not.
  
  working_params <- norm_working_params(num_states, num_subjects, 
                                        mu0, sigma0, gamma0, delta0, 
                                        stationary = stationary,
                                        gamma_pooled = gamma_pooled)
  hmm <- nlm(norm_loglikelihood, 
             working_params, 
             x = x, 
             num_states = num_states, 
             num_variables = num_variables, 
             num_subjects = num_subjects, 
             stationary = stationary,
             state_dep_dist_pooled = state_dep_dist_pooled,
             gamma_pooled = gamma_pooled,
             hessian = hessian) 
  
  pn    <- norm_natural_params(num_states = num_states, 
                               num_variables = num_variables, 
                               num_subjects = num_subjects,
                               working_params = hmm$estimate, 
                               stationary = TRUE,
                               state_dep_dist_pooled = state_dep_dist_pooled,
                               gamma_pooled = gamma_pooled) 
  if (hessian) {
    h <- hmm$hessian
    if (det(h) != 0) {
      h <- solve(h)
      jacobian <- norm_jacobian(num_states, num_variables, num_subjects, 
                                hmm$estimate, pn, 
                                state_dep_dist_pooled = state_dep_dist_pooled,
                                gamma_pooled = gamma_pooled)
      h <- t(jacobian)%*%h%*%jacobian
    }
  }
  if ((hessian) & (det(hmm$hessian) != 0)) {
    list(num_states = num_states, 
         num_variables = num_variables, 
         num_subjects = num_subjects, 
         mu = pn$mu, 
         sigma = pn$sigma,
         gamma = pn$gamma, 
         delta = pn$delta, 
         code = hmm$code, 
         max_loglikelihood = hmm$minimum,
         inverse_hessian = h)
  } else if (hessian) {
    list(num_states = num_states, 
         num_variables = num_variables, 
         num_subjects = num_subjects, 
         mu = pn$mu, 
         sigma = pn$sigma,
         gamma = pn$gamma, 
         delta = pn$delta, 
         code = hmm$code, 
         max_loglikelihood = hmm$minimum,
         hessian = hmm$hessian, 
         determinant = det(hmm$hessian))
  } else {
    list(num_states = num_states, 
         num_variables = num_variables, 
         num_subjects = num_subjects, 
         mu = pn$mu, 
         sigma = pn$sigma,
         gamma = pn$gamma, 
         delta = pn$delta, 
         code = hmm$code, 
         max_loglikelihood = hmm$minimum)
  }
}

norm_generate_sample <- function(num_sample, hmm, 
                                 state_dep_dist_pooled = FALSE, 
                                 gamma_pooled = FALSE) {
  #' Generate data from a normal HMM
  #' 
  #' This function generates data from a normal HMM with specified parameters.
  #' 
  #' @param num_sample The size of the desired sample (number of timesteps).
  #' @param hmm A list of parameters that specify the normal HMM, including 
  #'   `num_states`, `num_variables`, `num_subjects`, `mu`, `sigma`, `gamma`,
  #'   `delta`.
  
  state_vec     <- 1:hmm$num_states 
  num_variables <- hmm$num_variables
  num_subjects  <- hmm$num_subjects
  state         <- matrix(numeric(num_sample*num_subjects), 
                          nrow = num_sample, ncol = num_subjects)
  for (i in 1:num_subjects) {
    g_ind   <- i
    if (gamma_pooled) {
      g_ind <- 1
    }
    state[1, i]   <- sample(state_vec, 1, prob = hmm$delta[[g_ind]])
    for (t in 2:num_sample) {
      state[t, i] <- sample(state_vec, 1, 
                            prob = hmm$gamma[[g_ind]][state[(t - 1), i], ])
    }
  } 
  x <- array(dim = c(num_sample, num_variables, num_subjects))
  for (j in 1:num_variables) {
    for (i in 1:num_subjects) {
      s_ind   <- i
      if (state_dep_dist_pooled) {
        s_ind <- 1
      }
      x[, j, i] <- rnorm(num_sample, 
                         mean = hmm$mu[[j]][s_ind, state[, i]], 
                         sd = hmm$sigma[[j]][s_ind, state[, i]]) 
    }
  } 
  list(state = state, observ = x)
}

norm_viterbi <- function(x, hmm, 
                         state_dep_dist_pooled = FALSE, 
                         gamma_pooled = FALSE) {
  #' Global decoding by the Viterbi algorithm
  #'
  #' This function takes in data x assumed to be generated by the HMM hmm and 
  #' outputs the most likely sequence of states that could have generated the
  #' data using global decoding by the Viterbi algorithm.
  #' 
  #' @inheritParams norm_fit_hmm
  #' @inheritParams norm_viterbi
  
  n             <- nrow(x) 
  num_states    <- hmm$num_states
  num_variables <- hmm$num_variables
  num_subjects  <- hmm$num_subjects
  state_probs   <- list()
  sequence      <- matrix(0, nrow = n, ncol = num_subjects)
  
  for (i in 1:num_subjects) {
    s_ind   <- g_ind <- i
    if (state_dep_dist_pooled) {
      s_ind <- 1
    }
    if (gamma_pooled) {
      g_ind <- 1
    }
    state_probs[[i]] <- matrix(0, nrow = n, ncol = num_states)
    P                <- rep(1, num_states)
    for (j in 1:num_variables) {
      P <- P*dnorm(x[1, j, i], hmm$mu[[j]][s_ind, ], hmm$sigma[[j]][s_ind, ])
    } 
    forward_probs         <- hmm$delta[[g_ind]]*P
    state_probs[[i]][1, ] <- forward_probs/sum(forward_probs)
    
    for (t in 2:n) {
      P <- rep(1, num_states)
      for (j in 1:num_variables) {
        P <- P*dnorm(x[t, j, i], hmm$mu[[j]][s_ind, ], hmm$sigma[[j]][s_ind, ])
      }
      forward_probs <- apply(state_probs[[i]][t - 1, ]*
                               hmm$gamma[[g_ind]], 2, max)*P
      state_probs[[i]][t, ] <- forward_probs/sum(forward_probs)
    }
    sequence[n, i] <- which.max(state_probs[[i]][n, ])
    for (t in (n - 1):1){
      sequence[t, i] <- which.max(hmm$gamma[[g_ind]][, sequence[t + 1]]*
                                    state_probs[[i]][t, ])
    }
  }
  sequence
}

norm_logforward <- function(x, hmm, 
                            state_dep_dist_pooled = FALSE, 
                            gamma_pooled = FALSE) {
  #' Compute log forward probabilities
  #' 
  #' This function computes the log forward probabilities of the data based on
  #' the HMM hmm.
  #' 
  #' @inheritParams norm_viterbi
  
  n             <- nrow(x)
  num_states    <- hmm$num_states
  num_variables <- hmm$num_variables
  num_subjects  <- hmm$num_subjects
  lalpha        <- list() 
  for (i in 1:num_subjects) {
    s_ind   <- g_ind <- i
    if (state_dep_dist_pooled) {
      s_ind <- 1
    }
    if (gamma_pooled) {
      g_ind <- 1
    }
    lalpha[[i]] <- matrix(NA, nrow = num_states, ncol = n)
    P           <- rep(1, num_states)
    for (j in 1:num_variables) {
      P <- P*dnorm(x[1, j, i], hmm$mu[[j]][s_ind, ], hmm$sigma[[j]][s_ind, ])
    }
    forward_probs     <- hmm$delta[[g_ind]]*P
    sum_forward_probs <- sum(forward_probs) 
    loglikelihood     <- log(sum_forward_probs)
    forward_probs     <- forward_probs/sum_forward_probs
    lalpha[[i]][, 1]  <- loglikelihood + log(forward_probs)
    
    for (t in 2:n) {
      P     <- rep(1, num_states)
      for (j in 1:num_variables) {
        if (!is.na(x[t, j, i])) {
          P <- P*dnorm(x[t, j, i], hmm$mu[[j]][s_ind, ], 
                       hmm$sigma[[j]][s_ind, ])
        }
      }
      forward_probs     <- forward_probs %*% hmm$gamma[[g_ind]]*P 
      sum_forward_probs <- sum(forward_probs) 
      loglikelihood     <- loglikelihood + log(sum_forward_probs) 
      forward_probs     <- forward_probs/sum_forward_probs
      lalpha[[i]][, t]  <- loglikelihood + log(forward_probs)
    }
  } 
  lalpha
}


norm_forcast_psr <- function(x, hmm, state_dep_dist_pooled = FALSE, 
                             gamma_pooled = FALSE) {
  n             <- nrow(x)
  num_states    <- hmm$num_states
  num_variables <- hmm$num_variables
  num_subjects  <- hmm$num_subjects
  lalpha        <- norm_logforward(x, hmm, 
                                   state_dep_dist_pooled, 
                                   gamma_pooled)
  forcast_psr   <- list()
  for (i in 1:num_subjects) {
    s_ind   <- g_ind <- i
    if (state_dep_dist_pooled) {
      s_ind <- 1
    }
    if (gamma_pooled) {
      g_ind <- 1
    }
    pstepmat         <- matrix(NA, n, num_states)
    forcast_psr[[i]] <- rep(NA, n)
    ind.step         <- 1:n
    for (j in 1:num_variables) {
      ind.step <- sort(intersect(ind.step, which(!is.na(x[, j, i]))))
    }
    for (k in ind.step) {
      for (m in 1:num_states) {
        P   <- 1
        for (j in 1:num_variables) {
          P <- P*pnorm(x[k, j, i], mean = hmm$mu[[j]][s_ind, m], 
                       sd = hmm$sigma[[j]][s_ind, m])
        }
        pstepmat[k, m] <- P
      }
    }
    if (1 %in% ind.step) {
        forcast_psr[[i]][1] <- qnorm(hmm$delta[[g_ind]] %*% pstepmat[1, ])
    }
    for (t in 2:n) {
      c <- max(lalpha[[i]][, t - 1])
      a <- exp(lalpha[[i]][, t - 1] - c)
      if (t %in% ind.step) {
        forcast_psr[[i]][t] <- qnorm(t(a) %*% (hmm$gamma[[g_ind]]/sum(a)) %*% 
                                       pstepmat[t, ])
      }
    }
  }
  forcast_psr
}

norm_se <- function(hmm, state_dep_dist_pooled = FALSE, gamma_pooled = FALSE) {
  #' Compute standard errors of the fit HMM parameters
  #' 
  #' This function computes the standard errors of the fitted HMM parameters
  #' using the hessian computed from `norm_fit_hmm()`.
  #' 
  #' @inheritParams norm_generate_sample
  
  num_states    <- hmm$num_states
  num_variables <- hmm$num_variables
  num_subjects  <- hmm$num_subjects
  ns          <- num_subjects
  
  mu_start    <- 1
  mu_end      <- num_states*num_variables*num_subjects
  sigma_end   <- mu_end + num_states*num_variables*num_subjects
  mu_len      <- sigma_len <- num_subjects*num_states
  if (state_dep_dist_pooled) {
    mu_end    <- num_states*num_variables
    sigma_end <- mu_end + num_states*num_variables
    mu_len    <- sigma_len <- num_states
  }
  sigma_start <- mu_end + 1
  gamma_start <- sigma_end + 1
  gamma_end   <- sigma_end + num_subjects*num_states^2 - num_subjects*num_states
  if (gamma_pooled) {
    gamma_end <- sigma_end + num_states^2 - num_states
    ns        <- 1
  }
  h        <- sqrt(diag(hmm$inverse_hessian))
  mu.se    <- split_vec(h, mu_start, mu_end, mu_len)       
  sigma.se <- split_vec(h, sigma_start, sigma_end, sigma_len)
  
  for (j in 1:num_variables) {
    mu.se[[j]]    <- matrix(mu.se[[j]], ncol = num_states, byrow = TRUE)
    sigma.se[[j]] <- matrix(sigma.se[[j]], ncol = num_states, byrow = TRUE)
  }
  gamma.se <- list()
  g <- split_vec(h, gamma_start, gamma_end, num_states*(num_states - 1))
  for (i in 1:ns) {
    gamma.se[[i]] <- diag(num_states)
    gamma.se[[i]][!gamma.se[[i]]] <- g[[i]]
    diag(gamma.se[[i]]) <- NA
  }
  list(mu = hmm$mu,
       mu.se = mu.se,
       sigma = hmm$sigma,
       sigma.se = sigma.se,
       gamma = hmm$gamma,
       gamma.se = gamma.se)
}


norm_ci <- function(x, num_states, num_variables, num_subjects, params_se, 
                    state_dep_dist_pooled = FALSE,
                    x_step = 0.2, n = 100, level= 0.975) {
  #' Compute the confidence intervals of fitted normal distributions
  #' 
  #' This function computes the confidence intervals for the fitted normal state
  #' dependent distributions using the Monte Carlo approach.
  #' 
  #' @inheritParams norm_fit_hmm
  #' @param params_se A list of normal HMM parameters and thier standard errors
  #'   the output of `norm_se()`.
  #' @param x_step A value indicating the step length for the range of 
  #'   observation values.
  #' @param n A value indicating the number of samples used in the Monte Carlo
  #'   estimation of the confidence intervals.
  #' @param level The value indicating the desired confidence level.

  conf_intervals        <- list()
  for (i in 1:num_subjects) {
    conf_intervals[[i]] <- list()
    s_ind   <- i
    if (state_dep_dist_pooled) {
      s_ind <- 1
    }
    for (j in 1:num_variables) {
      range       <- seq(min(x[, j, i]), max(x[, j, i]), length.out = 100)
      xc          <- length(range)
      mu          <- params_se$mu[[j]][s_ind, ]
      mu.se       <- params_se$mu.se[[j]][s_ind, ]
      sigma       <- params_se$sigma[[j]][s_ind, ]
      sigma.se    <- params_se$sigma.se[[j]][s_ind, ]
      density.lst <- list()
      for (k in 1:num_states) {
        densities <- matrix(numeric(xc*n), ncol = xc, nrow = n)
        for (l in 1:n) {
          sample.mu      <- rnorm(1, mu[k], mu.se[k])
          sample.sigma   <- rnorm(1, sigma[k], sigma.se[k])
          densities[l, ] <- dnorm(range, sample.mu, sample.sigma)
        }
        density.lst[[k]] <- densities
      }
      upper <- matrix(numeric(xc*num_states), ncol = xc, nrow = num_states)
      lower <- matrix(numeric(xc*num_states), ncol = xc, nrow = num_states)
      
      for (k in 1:num_states) {
        densities <- density.lst[[k]]
        for (t in 1:xc) {
          upper[k, t] <- quantile(densities[, t], probs = level, na.rm = TRUE)
          lower[k, t] <- quantile(densities[, t], probs = 1 - level, 
                                  na.rm = TRUE)
        }
      }
      conf_intervals[[i]][[j]] <- list(range = range, 
                                       upper = upper, 
                                       lower = lower)
    }
  }
  conf_intervals
}

