#the contraints for the shape parameter alpha and the scale parameter theta of the gamma distribution 
#are alpha, theta > 0 so we can use the same method as with lambda of the poisson distribution
#A.1.1 (modified for gamma)
gam.HMM.pn2pw <- function(m, alpha, theta, gamma,delta=NULL, stationary=TRUE){
  #' Transform natural parameters to working
  #'
  #' This function is for gamma distributions.
  #' 
  #' m = number of states,
  #' alpha = shape parameter(s)
  #' theta = scale parameter(s)
  #' gamma = transition probability matrix
  #' delta = inital state distribution
  
  talpha <- log(alpha) 
  ttheta <- log(theta) 
  
  if(m==1) {
    return(talpha, ttheta)}
  
  foo <- log(gamma/diag(gamma)) 
  tgamma <- as.vector(foo[!diag(m)]) 
  
  if(stationary) {
    tdelta <- NULL}
  else {
    tdelta <- log(delta[-1]/delta[1])} 
  
  parvect <- c(talpha, ttheta, tgamma, tdelta) 
  return(parvect)
}


#A.1.2 (modified for gamma)
gam.HMM.pw2pn <- function(m, parvect, stationary=TRUE){
  #' Transform working parameters to natural
  #'
  #' This function is for gamma distributions.
  #' 
  #' m = number of states,
  #' parvect = (working shapes, working scales, working trans prob matrix entries, working initial dist) 
  
  alpha <- exp(parvect[1:m])
  theta <- exp(parvect[(m+1):(2*m)])
  gamma <- diag(m) 
  
  if (m==1) {
    return(list(alpha=alpha, theta=theta, gamma=gamma, delta=1))}
  
  gamma[!gamma] <- exp(parvect[(2*m+1):(2*m*m)]) 
  gamma <- gamma/apply(gamma ,1,sum) 
  
  if(stationary){
    delta<-solve(t(diag(m)-gamma+1),rep(1,m))} 
  else {
    foo<-c(1,exp(parvect[(2*m*m+1):(2*m*m+m-1)])) 
    delta <-foo/sum(foo)}
  
  return(list(alpha=alpha, theta=theta, gamma=gamma, delta=delta))
}


#A.1.3 (modified for gamma)
gam.HMM.mllk <- function(parvect, x, m, stationary=TRUE ,...){
  #' Compute -log-likelihood from working parameters
  #'
  #' This function is for gamma distributions.
  #' 
  #' parvect = (working shapes, working scales, working trans prob matrix entries, working initial dist),
  #' x = observations,
  #' m = number of states,
  
  if(m==1) {return(-sum(dgamma(x, shape=exp(parvect[1]), scale=exp(parvect[2]), log=TRUE)))} 
  
  n       <- length(x) 
  pn      <- gam.HMM.pw2pn(m,parvect ,stationary=stationary)
  foo     <- pn$delta*dgamma(x[1], shape=pn$alpha, scale=pn$theta) 
  sumfoo  <- sum(foo) 
  lscale  <- log(sumfoo)
  foo     <- foo/sumfoo
  
  for (i in 2:n){
    if (!is.na(x[i]))  {P <- dgamma(x[i], shape=pn$alpha, scale=pn$theta)}
    else {P <- rep(1,m)} 
    
    foo     <- foo %*% pn$gamma*P 
    sumfoo  <- sum(foo) 
    lscale  <- lscale+log(sumfoo) 
    foo     <- foo/sumfoo
  }
  mllk <- -lscale
  return(mllk)
}


#A.1.4 (modified for gamma)
gam.HMM.mle <- function(x,m,alpha0,theta0,gamma0,delta0=NULL,stationary=TRUE,...) {
  #' Compute Maximum Likelihood Estimate
  #'
  #' This function is for poission distributions starting with natural parameters
  #' 
  #' x        = observations,
  #' m        = number of states,
  #' alpha0   = inital guess for natural shapes
  #' theta0   = inital guess for natural scales
  #' gamma0   = initial guess for natural transition probability matrix
  #' delta0   = initial guess for initial state distribution
  
  parvect0 <- gam.HMM.pn2pw(m,alpha0, theta0, gamma0, delta0 , stationary=stationary)
  mod <- nlm(gam.HMM.mllk ,parvect0 ,x=x,m=m, stationary=stationary) 
  
  pn    <- gam.HMM.pw2pn(m=m, mod$estimate, stationary=TRUE) 
  mllk  <- mod$minimum 
  np    <- length(parvect0)
  AIC   <- 2*(mllk+np)
  n     <- sum(!is.na(x))
  BIC   <- 2*mllk+np*log(n)
  
  list(m=m, 
       alpha=pn$alpha, 
       theta=pn$theta, 
       gamma=pn$gamma, 
       delta=pn$delta, 
       code=mod$code, 
       mllk=mllk,
       AIC=AIC,
       BIC=BIC)
}


#A.1.5 (modified for gamma)
gam.HMM.generate_sample <- function(ns, mod){
  #' Generate a sample realization of an HMM
  #'
  #' This function is for poission distributions.
  #' 
  #' ns = length of realization 
  #' mod = HMM
  
  mvect <- 1:mod$m 
  state <- numeric(ns) 
  state[1]<- sample(mvect, 1, prob=mod$delta) 
  
  for (i in 2:ns) {
    state[i] <- sample(mvect, 1, prob=mod$gamma[state[i-1],])} 
  
  x <- rgamma(ns, shape=mod$alpha[state], scale=mod$theta[state])
  return(list(state = state, observ = x)) #I changed this from return(x) to hold on to the underlying states
}


#A.1.6 (modified for gamma)
gam.HMM.viterbi <-function(x, mod){
  #' Global decoding by the Viterbi algorithm
  #'
  #' This function is for poission distributions.
  #' 
  #' x = sequence of observations
  #' mod = HMM
  
  n       <- length(x) 
  xi      <- matrix(0,n,mod$m)
  foo     <- mod$delta*dgamma(x[1], shape=mod$alpha, scale=mod$theta) 
  xi[1,]  <- foo/sum(foo) 
  
  for (i in 2:n){
    foo<-apply(xi[i-1,]*mod$gamma ,2,max)*dgamma(x[i], shape=mod$alpha, scale=mod$theta)
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
gam.HMM.lforward <-function(x,mod){
  n           <- length(x) #number of observations
  lalpha      <- matrix(NA,mod$m,n) #mxn matrix
  foo         <- mod$delta*dgamma(x[1], shape=mod$alpha, scale=mod$theta) #probability vector that observ1 came from each state
  sumfoo      <- sum(foo) #sum of prob vector
  lscale      <- log(sumfoo) #log of sum of prob vector
  foo         <- foo/sumfoo #divide entries of prop vector by log of sum
  lalpha[,1]  <- lscale+log(foo) #set first column of matrix
  
  for (i in 2:n){
    foo         <- foo%*%mod$gamma*dgamma(x[i], shape=mod$alpha, scale=mod$theta)
    sumfoo      <- sum(foo)
    lscale      <- lscale+log(sumfoo)
    foo         <- foo/sumfoo
    lalpha[,i]  <- log(foo)+lscale
  }
  return(lalpha)
}


#A.1.8 Computing log(backward probabilities)
gam.HMM.lbackward <-function(x,mod){
  n           <- length(x) #number of observations
  m           <- mod$m #number of states
  lbeta       <- matrix(NA,m,n) #mxn matrix
  lbeta[,n]   <- rep(0,m) #fill last column with zeros
  foo         <- rep(1/m,m)
  lscale      <- log(m)
  
  for (i in (n-1):1){
    foo         <- mod$gamma%*%(dgamma(x[i+1], shape=mod$alpha, scale=mod$theta)*foo) 
    lbeta[,i]   <- log(foo)+lscale
    sumfoo      <- sum(foo)
    foo         <- foo/sumfoo
    lscale      <- lscale+log(sumfoo)
  }
  return(lbeta)
}


#A.1.9 Conditional probabilities
gam.HMM.conditional <- function(xc,x,mod){
  n     <- length(x)
  m     <- mod$m
  nxc   <- length(xc)
  dxc   <- matrix(NA,nrow=nxc,ncol=n)
  Px    <- matrix(NA,nrow=m,ncol=nxc)
  
  for (j in 1:nxc) {
    Px[,j] <-dgamma(xc[j], shape=mod$alpha, scale=mod$theta)}
  
  la      <- gam.HMM.lforward(x,mod)
  lb      <- gam.HMM.lbackward(x,mod)
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
gam.HMM.pseudo_residuals <- function(x,mod) {
  n <- length(x)
  xc <- sort(x)
  width <- diff(xc)
  cdists <- gam.HMM.conditional(xc,x,mod) 
  cdistno1 <- cdists[1:(n-1),]
  mult <- apply(cdistno1, 2, "*", width)
  mult <- rbind(mult, width[n-1]*cdists[n,])
  multsum <- apply(mult ,2,cumsum)
  
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



