#A.1.1 (modified for normal)
norm.HMM.pn2pw <- function(m, mu, sigma, gamma, delta=NULL, stationary=TRUE){
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
norm.HMM.pw2pn <- function(m, parvect, stationary=TRUE){
  #' Transform working parameters to natural
  #'
  #' This function is for normal distributions.
  #' 
  #' m = number of states,
  #' parvect = (working means, working sd, working trans prob matrix entries, working initial dist) 
  
  mu <- parvect[1:m]
  sigma <- exp(parvect[(m+1):(2*m)]) 
  gamma <- diag(m) 
  
  if (m==1) {
    return(list(mu=mu, sigma=sigma, gamma=gamma, delta=1))}
  
  gamma[!gamma] <- exp(parvect[(2*m+1):(2*m*m)]) 
  gamma <- gamma/apply(gamma, 1, sum) 
  
  if(stationary){
    delta<-solve(t(diag(m)-gamma+1),rep(1,m))}
  else {
    foo<-c(1,exp(parvect[(2*m*m+1):(2*m*m+m-1)])) 
    delta <-foo/sum(foo)}
  
  return(list(mu=mu, sigma=sigma, gamma=gamma, delta=delta))
}



#A.1.3 (modified for normal)
norm.HMM.mllk <- function(parvect, x, m, stationary=TRUE ,...){
  #' Compute -log-likelihood from working parameters
  #'
  #' This function is for normal distributions.
  #' 
  #' parvect = (working means, working sds, working trans prob matrix entries, working initial dist),
  #' x = observations,
  #' m = number of states,
  
  if(m==1) {return(-sum(dnorm(x, parvect[1], exp(parvect[2]), log=TRUE)))} 
  
  n       <- length(x) 
  pn      <- norm.HMM.pw2pn(m,parvect ,stationary=stationary) 
  foo     <- pn$delta*dnorm(x[1], pn$mu, pn$sigma) 
  sumfoo  <- sum(foo) 
  lscale  <- log(sumfoo)
  foo     <- foo/sumfoo
  
  for (i in 2:n){
    if (!is.na(x[i]))  {P <- dnorm(x[i], pn$mu, pn$sigma)}  
    else {P <- rep(1,m)} 
    
    foo     <- foo %*% pn$gamma*P 
    sumfoo  <- sum(foo) 
    lscale  <- lscale+log(sumfoo) 
    foo     <- foo/sumfoo
  }
  mllk <- -lscale
  return(mllk)
}


#A.1.4 (modified for normal)
norm.HMM.mle <- function(x, m, mu0, sigma0, gamma0, delta0=NULL, stationary=TRUE,...) {
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
  mllk  <- mod$minimum 
  np    <- length(parvect0) 
  AIC   <- 2*(mllk+np) 
  n     <- sum(!is.na(x)) 
  BIC   <- 2*mllk+np*log(n) 
  
  list(m=m, 
       mu=pn$mu, 
       sigma=pn$sigma,
       gamma=pn$gamma, 
       delta=pn$delta, 
       code=mod$code, 
       mllk=mllk,
       AIC=AIC,
       BIC=BIC)
}


#A.1.5 (modified for normal)
norm.HMM.generate_sample <- function(ns, mod){
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
    foo<-apply(xi[i-1,]*mod$gamma ,2,max)*dnorm(x[i], mod$mu, mod$sigma)
    xi[i,] <- foo/sum(foo) 
  }
  iv<-numeric(n) 
  iv[n] <-which.max(xi[n,])
  
  for (i in (n-1):1){
    iv[i] <- which.max(mod$gamma[,iv[i+1]]*xi[i,])
  }
  return(iv)
}


