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

