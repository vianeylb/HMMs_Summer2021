#function to remove +0i from real numbers
rid_Im <- function(x) {
  if (all(Im(z <- zapsmall(x))==0)) as.numeric(z) else x
}

#function to find stationary distribution (initial state distribution)
compute_delta <- function(A){
  e = eigen(t(A)) #find eigenvalues and left eigenvectors of A
  vec = rid_Im(e$vectors[,rid_Im(e$values)==1.0])/sum(rid_Im(e$vectors[,rid_Im(e$values)==1.0])) #eigenvector with eigenvalue 1 with components summing to 1
  return(vec)
}


#A.1.1
pois.HMM.pn2pw <- function(m,lambda,gamma,delta=NULL, stationary=TRUE){
  #' Transform natural parameters to working
  #'
  #' This function is for poission distributions.
  #' 
  #' m = number of states,
  #' lambda = vector of means for each state dependent poission distribution
  #' gamma = transition probability matrix
  #' delta = inital state distribution

  tlambda <- log(lambda) #transform means to working params
  
  if(m==1) {
    return(tlambda)}
  
  foo <- log(gamma/diag(gamma)) #define new matrix foo where the entries are the working params
  tgamma <- as.vector(foo[!diag(m)]) #get off-diagonal entries of foo -> (a_21, a_31, a_12, a_32, a_13, a_23) for m=3
  
  if(stationary) {
    tdelta <- NULL}
    else {
      tdelta <- log(delta[-1]/delta[1])} #transform delta to working params (not sure how)
  
  parvect <- c(tlambda ,tgamma ,tdelta) 
  return(parvect)
}


#A.1.2
pois.HMM.pw2pn <- function(m, parvect, stationary=TRUE){
  #' Transform working parameters to natural
  #'
  #' This function is for poission distributions.
  #' 
  #' m = number of states,
  #' parvect = (working means, working trans prob matrix entries, working initial dist) 

  lambda <- exp(parvect[1:m]) #transform working means to natural means
  gamma <- diag(m) #gamma = identity matrix size m
  
  if (m==1) {
    return(list(lambda=lambda,gamma=gamma,delta=1))}
  
  gamma[!gamma] <- exp(parvect[(m+1):(m*m)]) #set off-diagonals of gamma to the tranformed natural params
  gamma <- gamma/apply(gamma ,1,sum) #scale rows of gamma so that they sum to 1
  
  if(stationary){
    delta<-solve(t(diag(m)-gamma+1),rep(1,m))} #get stationary distribution (eigenvector with eigenvalue 1)
    else {
      foo<-c(1,exp(parvect[(m*m+1):(m*m+m-1)])) #transform delta to natural param
      delta <-foo/sum(foo)}
  
  return(list(lambda=lambda ,gamma=gamma ,delta=delta))
}


#A.1.3
pois.HMM.mllk <- function(parvect, x, m, stationary=TRUE ,...){
  #' Compute -log-likelihood from working parameters
  #'
  #' This function is for poission distributions.
  #' 
  #' parvect = (working means, working trans prob matrix entries, working initial dist),
  #' x = observations,
  #' m = number of states,

  if(m==1) {return(-sum(dpois(x,exp(parvect),log=TRUE)))} 
  
  n       <- length(x) #number of observations
  pn      <- pois.HMM.pw2pn(m,parvect ,stationary=stationary) #transform working params to natural
  foo     <- pn$delta*dpois(x[1],pn$lambda) #get vector likelihood of initial observation from HMM
  sumfoo  <- sum(foo) #sum to get likelihood of initial observation from HMM
  lscale  <- log(sumfoo)
  foo     <- foo/sumfoo
  
  for (i in 2:n){
    if (!is.na(x[i]))  {P <- dpois(x[i],pn$lambda)}  #probability for observation resulting from each state
      else {P <- rep(1,m)} #if there is an observation missing, probalility for
                           #that observation resulting from each state is 1
    
    foo     <- foo %*% pn$gamma*P #get vector likelihood of observations from HMM up to ith
    sumfoo  <- sum(foo) #sum to get likelihood of observations resulting from HMM up to ith observation
    lscale  <- lscale+log(sumfoo) 
    foo     <- foo/sumfoo
  }
  mllk <- -lscale
  return(mllk)
}


#A.1.4
pois.HMM.mle <- function(x,m,lambda0,gamma0,delta0=NULL,stationary=TRUE,...) {
  #' Compute Maximum Likelihood Estimate
  #'
  #' This function is for poission distributions starting with natural parameters
  #' 
  #' x        = observations,
  #' m        = number of states,
  #' lambda0  = inital guess for natural means
  #' gamma0   = initial guess for natural transition probability matrix
  #' dalta0   = initial guess for initial state distribution

  parvect0 <- pois.HMM.pn2pw(m,lambda0 ,gamma0 ,delta0 , stationary=stationary)
  mod <- nlm(pois.HMM.mllk ,parvect0 ,x=x,m=m, stationary=stationary) #non-linear minimization of -log-likelihood
                                                                      #with initial working params parvect0
  
  pn    <- pois.HMM.pw2pn(m=m, mod$estimate, stationary=TRUE) #get natural params from parms with min mllk
  mllk  <- mod$minimum #the minimum -log-likelihood
  np    <- length(parvect0) #number of initial working params
  AIC   <- 2*(mllk+np) #Akaike Information Criterion
  n     <- sum(!is.na(x)) #number of observations
  BIC   <- 2*mllk+np*log(n) #Bayesian Information Criterion
  
  list(m=m, 
       lambda=pn$lambda, 
       gamma=pn$gamma, 
       delta=pn$delta, 
       code=mod$code, #integer indicating why the optimization process nlm of mllk terminated 
       mllk=mllk,
       AIC=AIC,
       BIC=BIC)
}


pois.HMM.fit <- function(x,m,lambda0,gamma0,delta0=NULL,stationary=TRUE,...) {

  parvect0 <- pois.HMM.pn2pw(m,lambda0 ,gamma0 ,delta0 , stationary=TRUE)
  mod <- nlm(pois.HMM.mllk ,parvect0 ,x=x,m=m, stationary=stationary)
  pn    <- pois.HMM.pw2pn(m=m, mod$estimate, stationary=stationary)
  list(m=m, 
       lambda=pn$lambda, 
       gamma=pn$gamma, 
       delta=pn$delta)
}

####### NOTE for A.1.5, A.1.6 ######
#I'm guessing mod is the data for an HMM and it should be of the following form with correct names:
#mod <- list(m,       #number of states
#            gamma,   #transition probability matrix
#            lambda,  #vector of means of the poisson distributions
#            delta)   #inital state distribution
####################################


#A.1.5
pois.HMM.generate_sample <- function(ns, mod){
  #' Generate a sample realization of an HMM
  #'
  #' This function is for poission distributions.
  #' 
  #' ns = length of realization 
  #' mod = HMM
  
  mvect <- 1:mod$m #vector of states 1 through m
  state <- numeric(ns) #create vector of length ns to be filled with state realization sequence
  state[1]<- sample(mvect, 1, prob=mod$delta) #get initial state from initial distribution
  
  for (i in 2:ns) {
    state[i] <- sample(mvect, 1, prob=mod$gamma[state[i-1],])} #get state sequence using trans prob matrix
  
  x <- rpois(ns,lambda=mod$lambda[state]) #generate observations based on state sequence and pois dists
  return(list(state = state, observ = x)) #I changed this from return(x) to hold on to the underlying states
}


#A.1.6
pois.HMM.viterbi <-function(x, mod){
  #' Global decoding by the Viterbi algorithm
  #'
  #' This function is for poission distributions.
  #' 
  #' x = sequence of observations
  #' mod = HMM

  n       <- length(x) #number of observations
  xi      <- matrix(0,n,mod$m) #nxm matrix of zeros 
  foo     <- mod$delta*dpois(x[1],mod$lambda) #vector of probability that initial observation
                                              #resulted from each state
  xi[1,]  <- foo/sum(foo) #ensure probability vector sums to 1, set as row 1 of xi
  
  for (i in 2:n){
    foo<-apply(xi[i-1,]*mod$gamma ,2,max)*dpois(x[i],mod$lambda) #get max prob of row x[i-1,] * prob that
                                                                 #observation i resulted from each state
    xi[i,] <- foo/sum(foo) 
  }
  iv<-numeric(n) #create vector of lengeth n
  iv[n] <-which.max(xi[n,]) #last element of iv is state that maximizes last row of xi
  
  for (i in (n-1):1){
    iv[i] <- which.max(mod$gamma[,iv[i+1]]*xi[i,]) #elements of vi are states that maximize xi *
                                                   #prob of transitioning to the next state
  }
  return(iv)
}


#A.1.7 Computing log(forward probabilities)
#this function gives a matrix of probabilities in logarithmic form rows=states, columns=observation
#probability that 
pois.HMM.lforward <-function(x,mod){
  n           <- length(x) #number of observations
  lalpha      <- matrix(NA,mod$m,n) #mxn matrix
  foo         <- mod$delta*dpois(x[1],mod$lambda) #probability vector of observ1 given state dist
  sumfoo      <- sum(foo) #sum of prob vector
  lscale      <- log(sumfoo) #log of sum of prob vector
  foo         <- foo/sumfoo #divide entries of prop vector by log of sum
  lalpha[,1]  <- lscale+log(foo) #set first column of matrix
  
  for (i in 2:n){
    foo         <- foo%*%mod$gamma*dpois(x[i],mod$lambda) #prob vector of observi given state dist and trans prob
    sumfoo      <- sum(foo)
    lscale      <- lscale+log(sumfoo)
    foo         <- foo/sumfoo
    lalpha[,i]  <- log(foo)+lscale
  }
  return(lalpha)
}


#A.1.8 Computing log(backward probabilities)
pois.HMM.lbackward <-function(x,mod){
  n           <- length(x) #number of observations
  m           <- mod$m #number of states
  lbeta       <- matrix(NA,m,n) #mxn matrix
  lbeta[,n]   <- rep(0,m) #fill last column with zeros
  foo         <- rep(1/m,m)
  lscale      <- log(m)

  for (i in (n-1):1){
    foo         <- mod$gamma%*%(dpois(x[i+1],mod$lambda)*foo) 
    lbeta[,i]   <- log(foo)+lscale
    sumfoo      <- sum(foo)
    foo         <- foo/sumfoo
    lscale      <- lscale+log(sumfoo)
  }
  return(lbeta)
}


#A.1.9 Conditional probabilities
pois.HMM.conditional <- function(xc,x,mod){
  n     <- length(x)
  m     <- mod$m
  nxc   <- length(xc)
  dxc   <- matrix(NA,nrow=nxc,ncol=n)
  Px    <- matrix(NA,nrow=m,ncol=nxc)
  
  for (j in 1:nxc) {
    Px[,j] <-dpois(xc[j],mod$lambda)}
  
  la      <- pois.HMM.lforward(x,mod)
  lb      <- pois.HMM.lbackward(x,mod)
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
pois.HMM.pseudo_residuals <- function(x,mod) {
  n         <- length(x)
  #matrix where rows->xc, columns->observation num
  #entry i,j is prob that x[j]=xc[i] given x[1:j] and x[(j+1):n]
  cdists    <- pois.HMM.conditional(xc=0:max(x),x,mod) 
  #cumulatively sum the columns of cdists (ex: cumsum(c(1,2,3))=c(1,3,6))
  #add a row of zeros to the top to account for x[1]=0
  cumdists  <- rbind(rep(0,n),apply(cdists ,2,cumsum)) 
  ulo <- uhi <- rep(NA,n) 
    
  for (i in 1:n){
      ulo[i] <- cumdists[x[i]+1,i] #prob that observ[i]<x[i]
      uhi[i] <- cumdists[x[i]+2,i] #prob that observ[i]<x[i+1]
  }
  umi   <- 0.5*(ulo+uhi) #get middle value
  #get value from N(0,1) with corresponding percentile of ulo, umi, and uhi
  npsr  <- qnorm(rbind(ulo,umi,uhi)) 
  return(data.frame(lo=npsr[1,], mi=npsr[2,], hi=npsr[3,]))
}


pois.HMM.initial_params <- function(x, m){
  sum = summary(x)
  mean= as.numeric(sum[4])
  min = as.numeric(sum[1])
  max = as.numeric(sum[6])
  Q1 = as.numeric(sum[2])
  Q3 = as.numeric(sum[5])
  med = as.numeric(sum[3])
  diff = max-min
  
  if (m==1){
    lambda0 = mean(x)
    gamma0 = matrix(c(1))
  }
  else if (m==2){
    lambda0 = c(Q1, Q3)
    gamma0 = matrix(c(0.9, 0.1, 
                      0.1, 0.9), 2, 2, byrow=TRUE)
  }
  else if (m==3){
    lambda0 = c(Q1, med, Q3)
    gamma0 = matrix(c(0.9, 0.05, 0.05,
                      0.05, 0.9, 0.05,
                      0.05, 0.05, 0.9), 3, 3, byrow=TRUE)
  }
  else if (m==4){
    d = diff/5
    lambda0 = c(min+d, min+2*d, min+3*d, min+4*d)
    gamma0 = matrix(c(0.91, 0.03, 0.03, 0.03,
                      0.03, 0.91, 0.03, 0.03,
                      0.03, 0.03, 0.91, 0.03,
                      0.03, 0.03, 0.03, 0.91), 4, 4, byrow=TRUE)
  }
  else if (m==5){
    d = diff/6
    lambda0 = c(min+d, min+2*d, min+3*d, min+4*d, min+5*d)
    gamma0 = matrix(c(0.90, 0.025, 0.025, 0.025, 0.025,
                      0.025, 0.90, 0.025, 0.025, 0.025,
                      0.025, 0.025, 0.90, 0.025, 0.025,
                      0.025, 0.025, 0.025, 0.90, 0.025,
                      0.025, 0.025, 0.025, 0.025, 0.90), 5, 5, byrow=TRUE)
  }
  return(list(lambda0=lambda0, gamma0=gamma0))
}


pois.HMM.CI <- function(x, n, modfit, level=0.975, stationary=TRUE){
  ns <- length(x)
  m <- modfit$m
  lambdas = matrix(numeric(m*n), nrow = n, ncol=m) #matrix to be filled with fitted lambdas
  gammas = matrix(numeric(m*m*n), nrow = n, ncol = m*m) #matrix to be filled with entries of gamma
  deltas = matrix(numeric(m*n), nrow = n, ncol=m) #matrix to be filled with fitted deltas
  
  for (i in 1:n){
    sample <- pois.HMM.generate_sample(ns, modfit) #generate observations based on modfit
    x <- sample$observ #get observations from generated sample
    mod <- pois.HMM.fit(x,m,modfit$lambda, modfit$gamma, modfit$delta, stationary = stationary) #fit model to generated observations
    lambdas[i,] <- mod$lambda #add fitted lambda to lambdas matrix
    gammas[i,] = as.vector(t(mod$gamma)) #add fitted gamma as row
    deltas[i,] = mod$delta #add fitted delta to deltas matrix
  }
  
  lambda.cov = cov(lambdas) #get var-covar matrix of lambdas
  lambda.upper = modfit$lambda + (qnorm(level) * sqrt(diag(lambda.cov))) #calculate upper 95% CI from var
  lambda.lower = pmax(modfit$lambda - (qnorm(level) * sqrt(diag(lambda.cov))),0) #calculate lower 95% CI from var
  
  delta.cov = cov(deltas) #get var-covar matrix of lambdas
  delta.upper = modfit$delta + (qnorm(level) * sqrt(diag(delta.cov))) #calculate upper 95% CI from var
  delta.lower = pmax(modfit$delta - (qnorm(level) * sqrt(diag(delta.cov))),0) #calculate lower 95% CI from var
  
  gammafit = as.vector(t(modfit$gamma))
  gamma.cov = cov(gammas)
  gamma.upper = gammafit + (qnorm(level) * sqrt(diag(gamma.cov))) #calculate upper 95% CI from var
  gamma.upper = matrix(gamma.upper, m,m, byrow=TRUE)
  gamma.lower = pmax(gammafit - (qnorm(level) * sqrt(diag(gamma.cov))),0) #calculate lower 95% CI from var
  gamma.lower = matrix(gamma.lower, m,m, byrow=TRUE)
  
  result = list("lambda" = modfit$lambda, 
                "lambda.upper.conf" = lambda.upper,
                "lambda.lower.conf" = lambda.lower,
                "gamma" = modfit$gamma,
                "gamma.upper.conf" = gamma.upper,
                "gamma.lower.conf" = gamma.lower,
                "delta" = modfit$delta, 
                "delta.upper.conf" = delta.upper,
                "delta.lower.conf" = delta.lower)
  return(result)
}



