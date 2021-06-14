library(mvtnorm)
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
#q is the number of response variables
mvnorm.HMM.pw2pn <- function(m, q, parvect, stationary=TRUE){
  mu <- split(parvect[1:(m*q)], ceiling(seq_along(parvect[1:(m*q)])/m))
  sigma <- split(exp(parvect[(m*q+1):(2*m*q)]), ceiling(seq_along(exp(parvect[(m*q+1):(2*q*m)]))/m))
  gamma <- diag(m) 
  
  if (m==1) {
    return(list(mu=mu, sigma=sigma, gamma=gamma, delta=1))}
  
  gamma[!gamma] <- exp(parvect[(2*m*q+1):(2*m*q + m*m-m)]) 
  gamma <- gamma/apply(gamma, 1, sum) 
  
  if(stationary){
    delta<-solve(t(diag(m)-gamma+1),rep(1,m))}
  else {
    foo<-c(1,exp(parvect[(2*m*q + m*m-m+1):(length(parvect)-1)])) 
    delta <-foo/sum(foo)}
  
  return(list(mu=mu, sigma=sigma, gamma=gamma, delta=delta))
}


#Compute -log-likelihood from working parameters
mvnorm.HMM.mllk <- function(parvect, x, m, q, stationary=TRUE ,...){
  n       <- nrow(x)
  pn      <- mvnorm.HMM.pw2pn(m, q, parvect ,stationary=stationary) 
  P <- rep(1,m)
  for (k in 1:q){
    P <- P*dnorm(x[1,k], pn$mu[[k]], pn$sigma[[k]])}
  foo     <- pn$delta*P
  sumfoo  <- sum(foo) 
  lscale  <- log(sumfoo)
  foo     <- foo/sumfoo
  
  for (i in 2:n){
    P <- rep(1,m)
    for (k in 1:q){
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
mvnorm.HMM.mle <- function(x, m, q, mu0, sigma0, gamma0, delta0=NULL, stationary=TRUE,...) {
  
  parvect0 <- mvnorm.HMM.pn2pw(m, mu0, sigma0, gamma0, delta0, stationary=stationary)
  mod <- nlm(mvnorm.HMM.mllk, parvect0, x=x, m=m, q=q, stationary=stationary) 
  
  pn    <- mvnorm.HMM.pw2pn(m=m, q=q, mod$estimate, stationary=TRUE) 
  mllk  <- mod$minimum 
  #np    <- length(parvect0) 
  #AIC   <- 2*(mllk+np) 
  #n     <- sum(!is.na(x)) 
  #BIC   <- 2*mllk+np*log(n) 
  
  list(m=m, 
       q=q,
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
  q <- mod$q
  state <- numeric(ns)
  state[1] <- sample(mvect, 1, prob=mod$delta)
  
  for (i in 2:ns) {
    state[i] <- sample(mvect, 1, prob=mod$gamma[state[i-1],])
  } 
  x <- matrix(numeric(ns*q), ncol=q, nrow=ns)
  for (k in 1:q){
    x[,k] <- rnorm(ns, mean=mod$mu[[k]][state], sd=mod$sigma[[k]][state]) 
  }
  return(list(state = state, observ = x)) 
}


