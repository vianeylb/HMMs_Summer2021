library(tidyverse)
library(ggplot2)
library(ggfortify)
library(plyr)
library(bayesforecast)
library(gridExtra)

#q: order of autoregressive model
#phi: autoregressive model parameters, as list of vectors

ar_hmm_generate_sample <- function(ns, mod) {
  mvect <- 1:mod$m
  state <- numeric(ns)
  state[1] <- sample(mvect, 1, prob = mod$delta)
  for (i in 2:ns) state[i] <- sample(mvect, 1, prob = mod$gamma[state[i - 1], ])
  q <- mod$q
  x <- numeric(ns)
  x[1] <- rnorm(1, mean = mod$mu[state[1]], sd = mod$sigma[state[1]])
  for (i in 2:q){
    phi <- mod$phi[[state[i]]][1:(i-1)]
    mean <- mod$mu[state[i]] + phi%*%x[(i-1):1]
    x[i] <- rnorm(1, mean = mean, sd = mod$sigma[state[i]])
  }
  for (i in (q + 1):ns){
    mean <- mod$mu[state[i]]+ mod$phi[[state[i]]]%*%x[(i - 1):(i - q)]
    x[i] <- rnorm(1, mean = mean, sd = mod$sigma[state[i]])
  }
  return(data_frame(index = c(1:ns), state = state, obs = x))
}

ar_hmm_pn2pw <- function(m, mu, sigma, gamma, phi,
                           delta = NULL, stationary = TRUE) {
  tsigma <- log(sigma)
  foo <- log(gamma / diag(gamma))
  tgamma <- as.vector(foo[!diag(m)])
  tphi <- unlist(phi, use.names = FALSE)
  if (stationary) {
    tdelta <- NULL
  }
  else {
    tdelta <- log(delta[-1] / delta[1])
  }
  parvect <- c(mu, tsigma, tgamma, tphi, tdelta)
  return(parvect)
}

ar_hmm_pw2pn <- function(m, q, parvect, stationary = TRUE) {
  mu <- parvect[1:m]
  sigma <- exp(parvect[(m + 1):(2 * m)])
  gamma <- diag(m)
  gamma[!gamma] <- exp(parvect[(2 * m + 1):(m + m * m)])
  gamma <- gamma / apply(gamma, 1, sum)
  
  count <- m + m * m + 1
  phi <- list()
  for (i in 1:m) {
    phi[[i]] <- parvect[count:(count + q - 1)]
    count <- count + q
  }
  
  if (stationary) {
    delta <- solve(t(diag(m) - gamma + 1), rep(1, m))
  }
  else {
    foo <- c(1, exp(parvect[count:(count + m - 2)]))
    delta <- foo / sum(foo)
  }
  return(list(mu = mu, sigma = sigma, gamma = gamma, phi = phi, delta = delta))
}

ar_hmm_mllk <- function(parvect, x, m, q, stationary = TRUE, ...) {
  n <- length(x)
  pn <- ar_hmm_pw2pn(m, q, parvect, stationary = stationary)
  foo <- pn$delta * dnorm(x[1], pn$mu, pn$sigma)
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo / sumfoo
  
  phi <- matrix(unlist(mod$phi, use.names = FALSE), nrow=q, byrow=TRUE)
  for (i in 2:q){
    phi2 <- phi[, 1:(i - 1)]
    x_lag <-  x[(i - 1):1]
    if (i == 2){
      mean <- mod$mu + phi2 * x_lag
    }
    else {
      mean <- mod$mu + as.vector(phi2 %*% x_lag)
    }
    p <- dnorm(x[i], mean = mean, sd = mod$sigma)
    foo <- foo %*% pn$gamma * p
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
  }

  for (i in (q + 1):n) {
    x_lag <- x[(i - 1):(i - q)]
    mean <- mod$mu + as.vector(phi%*%x_lag)
    p <- dnorm(x[i], mean = mean, sd = mod$sigma)
    foo <- foo %*% pn$gamma * p
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
  }
  mllk <- -lscale
  return(mllk)
}

ar_hmm_mle <- function(x, m, q, mu0, sigma0, gamma0, phi0,
                         delta0 = NULL, stationary = TRUE,
                         hessian = FALSE, ...) {
  parvect0 <- ar_hmm_pn2pw(m, mu0, sigma0, gamma0, phi0, delta0,
                             stationary = stationary)
  mod <- nlm(ar_hmm_mllk, parvect0, x = x, m = m, q = q,
             stationary = stationary,
             hessian = hessian)
  pn <- ar_hmm_pw2pn(m, q, mod$estimate, 
                       stationary = stationary)
  mllk <- mod$minimum
  
  np <- length(parvect0)
  aic <- 2 * (mllk + np)
  n <- sum(!is.na(x))
  bic <- 2 * mllk + np * log(n)
  
  return(list(
    m = m, q = q, mu = pn$mu, sigma = pn$sigma, gamma = pn$gamma, phi = pn$phi, delta = pn$delta,
    code = mod$code, mllk = mllk, aic = aic, bic = bic)
  )
}

ar_hmm_viterbi <- function(x, mod) {
  n <- length(x)
  xi <- matrix(0, n, mod$m)
  foo <- mod$delta * dnorm(x[1], mod$mu, mod$sigma)
  xi[1, ] <- foo / sum(foo)
  
  q <- mod$q
  phi <- matrix(unlist(mod$phi, use.names = FALSE), nrow = q, byrow = TRUE)
  for (t in 2:q){
    phi2 <- phi[, 1:(t - 1)]
    x_lag <-  x[(t - 1):1]
    if (t == 2){
      mean <- mod$mu + phi2 * x_lag
    }
    else {
      mean <- mod$mu + as.vector(phi2 %*% x_lag)
    }
    p <- dnorm(x[t], mean = mean, sd = mod$sigma)
    foo <- apply(xi[t - 1, ] * mod$gamma, 2, max) * p
    xi[t, ] <- foo / sum(foo)
  }
  
  for (t in (q + 1):n) {
    x_lag <- x[(t - 1):(t - q)]
    mean <- mod$mu + as.vector(phi%*%x_lag)
    p <- dnorm(x[t], mean = mean, sd = mod$sigma)
    foo <- apply(xi[t - 1, ] * mod$gamma, 2, max) * p
    xi[t, ] <- foo / sum(foo)
  }
  
  
  iv <- numeric(n)
  iv[n] <- which.max(xi[n, ])
  for (t in (n - 1):1) {
    iv[t] <- which.max(mod$gamma[, iv[t + 1]] * xi[t, ])
  }
  return(data_frame(index = 1:n, states = iv))
}


ar_hmm_lforward <- function(x, mod) {
  n <- length(x)
  lalpha <- matrix(NA, mod$m, n)
  foo <- mod$delta * dnorm(x[1], mod$mu, mod$sigma)
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo / sumfoo
  lalpha[, 1] <- lscale + log(foo)
  
  q <- mod$q
  phi <- matrix(unlist(mod$phi, use.names = FALSE), nrow=q, byrow=TRUE)
  for (i in 2:q){
    phi2 <- phi[, 1:(i - 1)]
    x_lag <-  x[(i - 1):1]
    if (i == 2){
      mean <- mod$mu + phi2 * x_lag
    }
    else {
      mean <- mod$mu + as.vector(phi2 %*% x_lag)
    }
    p <- dnorm(x[i], mean = mean, sd = mod$sigma)
    foo <- foo %*% mod$gamma * p
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
    lalpha[, i] <- log(foo) + lscale
  }
  
  for (i in (q + 1):n) {
    x_lag <- x[(i - 1):(i - q)]
    mean <- mod$mu + as.vector(phi %*% x_lag)
    p <- dnorm(x[i], mean = mean, sd = mod$sigma)
    foo <- foo %*% mod$gamma * p
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
    lalpha[, i] <- log(foo) + lscale
  }
  
  return(lalpha)
}

ar_hmm_lbackward <- function(x, mod) {
  n <- length(x)
  m <- mod$m
  lbeta <- matrix(NA, m, n)
  lbeta[, n] <- rep(0, m)
  foo <- rep(1 / m, m)
  lscale <- log(m)
  
  q <- mod$q
  phi <- matrix(unlist(mod$phi, use.names = FALSE), nrow=q, byrow=TRUE)
  for (i in (n - 1):q) {
    x_lag <- x[i:(i - q + 1)]
    mean <- mod$mu + as.vector(phi %*% x_lag)
    p <- dnorm(x[i + 1], mean = mean, sd = mod$sigma)
    foo <- mod$gamma %*% p * foo
    lbeta[, i] <- log(foo) + lscale
    sumfoo <- sum(foo)
    foo <- foo / sumfoo
    lscale <- lscale + log(sumfoo)
  }
  
  for (i in (q - 1):1){
    phi2 <- phi[, 1:i]
    x_lag <-  x[i:1]
    if (i == 1){
      mean <- mod$mu + phi2 * x_lag
    }
    else {
      mean <- mod$mu + as.vector(phi2 %*% x_lag)
    }
    p <- dnorm(x[i + 1], mean = mean, sd = mod$sigma)
    foo <- mod$gamma %*% p * foo
    lbeta[, i] <- log(foo) + lscale
    sumfoo <- sum(foo)
    foo <- foo / sumfoo
    lscale <- lscale + log(sumfoo)
  }
  return(lbeta)
}

ar_hmm_pseudo_residuals <- function(x, mod, type, stationary) {
  if (stationary) {
    delta <- solve(t(diag(mod$m) - mod$gamma + 1), rep(1, mod$m))
  }
  else {
    delta <- mod$delta
  }
  if (type == "ordinary") {
    n <- length(x)
    la <- ar_hmm_lforward(x, mod)
    lb <- ar_hmm_lbackward(x, mod)
    lafact <- apply(la, 2, max)
    lbfact <- apply(lb, 2, max)
    
    q <- mod$q
    phi <- matrix(unlist(mod$phi, use.names = FALSE), nrow=q, byrow=TRUE)
    p <- matrix(NA, n, mod$m)
    p[1, ] <- pnorm(x[1], mean = mod$mu, sd = mod$sigma)
    for (i in 2:q){
      phi2 <- phi[, 1:(i - 1)]
      x_lag <-  x[(i - 1):1]
      if (i == 2){
        mean <- mod$mu + phi2 * x_lag
      }
      else {
        mean <- mod$mu + as.vector(phi2 %*% x_lag)
      }
      p[i, ] <- pnorm(x[i], mean = mean, sd = mod$sigma)
    }
    for (i in (q + 1):n) {
      x_lag <- x[(i - 1):(i - q)]
      mean <- mod$mu + as.vector(phi %*% x_lag)
      p[i, ] <- pnorm(x[i], mean = mean, sd = mod$sigma)
    }
    
    npsr <- rep(NA, n)
    npsr[1] <- qnorm(delta %*% p[1, ])
    for (i in 2:n) {
      a <- exp(la[, i - 1] - lafact[i])
      b <- exp(lb[, i] - lbfact[i])
      foo <- (a %*% mod$gamma) * b
      foo <- foo / sum(foo)
      npsr[i] <- qnorm(foo %*% p[i, ])
    }
    
    return(data_frame(npsr, x, index = c(1:n)))
  }
  else if (type == "forecast") {
    n <- length(x)
    la <- ar_hmm_lforward(x, mod)
    
    q <- mod$q
    phi <- matrix(unlist(mod$phi, use.names = FALSE), nrow=q, byrow=TRUE)
    p <- matrix(NA, n, mod$m)
    p[1, ] <- pnorm(x[1], mean = mod$mu, sd = mod$sigma)
    for (i in 2:q){
      phi2 <- phi[, 1:(i - 1)]
      x_lag <-  x[(i - 1):1]
      if (i == 2){
        mean <- mod$mu + phi2 * x_lag
      }
      else {
        mean <- mod$mu + as.vector(phi2 %*% x_lag)
      }
      p[i, ] <- pnorm(x[i], mean = mean, sd = mod$sigma)
    }
    for (i in (q + 1):n) {
      x_lag <- x[(i - 1):(i - q)]
      mean <- mod$mu + as.vector(phi %*% x_lag)
      p[i, ] <- pnorm(x[i], mean = mean, sd = mod$sigma)
    }
    
    npsr <- rep(NA, n)
    npsr[1] <- qnorm(delta %*% p[1, ])
    for (i in 2:n) {
      la_max <- max(la[, i - 1])
      a <- exp(la[, i - 1] - la_max)
      npsr[i] <- qnorm(t(a) %*% (mod$gamma / sum(a)) %*% p[i, ])
    }
    
    return(data_frame(npsr, x, index = c(1:n)))
  }
}


m <- 3
mu <- c(2, 5, 8)
sigma <- c(2, 4, 6)
gamma <- matrix(c(
  0.1, 0.2, 0.7,
  0.3, 0.4, 0.3,
  0.6, 0.3, 0.1
), nrow = m, byrow = TRUE)
delta <- c(0.4, 0.2, 0.4)
q <- 3
phi <- list(c(0.5, 0.3, 0.1),
            c(0.4, 0.2, 0.1),
            c(0.1, 0.1, 0.1))
mod <- list(m = m, mu = mu, sigma = sigma, gamma = gamma, phi = phi, delta = delta, q = q)
ar_sample <- ar_hmm_generate_sample(1000, mod)
# Plots 
graph_hmm_output(ar_sample)
graph_hmm_hist(ar_sample)
ggacf(ar_sample$obs) +
  theme_minimal()
ggpacf(ar_sample$obs) +
  theme_minimal()
# Get MLE
mu0 <- c(5, 5, 5)
sigma0 <- c(5, 5, 5)
gamma0 <- matrix(rep(1 / m, m * m), nrow = m, byrow = TRUE)
delta0 <- rep(1 / m, m)
phi0 <- list(rep(0.3, 3),
             rep(0.3, 3),
             rep(0.3, 3))
ar_mle <- ar_hmm_mle(ar_sample$obs, m, q, mu0, sigma0, gamma0, phi0,
                         delta0, stationary = TRUE, hessian = TRUE)
ar_mle
ar_decoding <- ar_hmm_viterbi(ar_sample$obs, ar_mle)
ar_decoding
count(ar_decoding$states)
count(ar_decoding$states - ar_sample$state)
# Get pseudo-residuals
ar_pr <- ar_hmm_pseudo_residuals(ar_sample$obs, ar_mle,
                                     "forecast", stationary = FALSE)
# Index plot of pseudo-residuals
ggplot(ar_pr) +
  geom_point(aes(x = index, y = npsr), size = 0.5, colour = "black") +
  theme_minimal()
# Histogram of pseudo-residuals
ggplot(ar_pr, aes(npsr)) +
  geom_histogram(aes(y = ..density..), colour = "navy", fill = "light blue") +
  stat_function(fun = dnorm, colour = "red") +
  theme_minimal()
# QQ plot of pseudo-residuals
ggplot(ar_pr, aes(sample = npsr)) +
  stat_qq() +
  stat_qq_line() +
  theme_minimal()
# ACF of pseudo-residuals
ggacf(ar_pr$npsr) +
  theme_minimal()

parvect <- ar_hmm_pn2pw(m, mu, sigma, gamma, phi,
                        delta, stationary = FALSE)
parvect
pn <- ar_hmm_pw2pn(m, q, parvect, stationary = FALSE)
pn
ar_hmm_mllk(parvect, ar_sample$obs, m, q, stationary = FALSE)
ar_hmm_lforward(ar_sample$obs, ar_mle)
ar_hmm_lbackward(ar_sample$obs, ar_mle)
