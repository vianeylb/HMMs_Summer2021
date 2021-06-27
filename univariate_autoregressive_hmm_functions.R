library(tidyverse)
library(ggplot2)
library(ggfortify)
library(plyr)
library(bayesforecast)
library(gridExtra)

#q: order of autoregressive model
#phi: autoregressive model parameters, as list of vectors

# Get mean for given x in autoregressive series
# phi is a qxq matrix 
# i is the time index
get_ar_mean <- function(mu, phi, x, q, i){
  if (i == 1) {
    mean <- mu
  }
  else if (i <= q){
    phi <- phi[, 1:(i - 1)]
    x_lag <-  x[(i - 1):1]
    if (i == 2)
      mean <- mu + phi * x_lag
    else{
      mean <- mu + as.vector(phi %*% x_lag)
    }
  } 
  else {
    x_lag <- x[(i - 1):(i - q)]
    mean <- mu + as.vector(phi %*% x_lag)
  }
  return(mean)
}

ar_hmm_generate_sample <- function(ns, mod) {
  mvect <- 1:mod$m
  state <- numeric(ns)
  state[1] <- sample(mvect, 1, prob = mod$delta)
  for (i in 2:ns) state[i] <- sample(mvect, 1, prob = mod$gamma[state[i - 1], ])
  x <- numeric(ns)
  phi <- matrix(unlist(mod$phi, use.names = FALSE), nrow = mod$q, byrow = TRUE)
  for (i in 1:ns){
    mean <- get_ar_mean(mod$mu, phi, x, mod$q, i)
    x[i] <- rnorm(1, mean = mean[state[i]], sd = mod$sigma[state[i]])
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
  
  phi <- matrix(unlist(pn$phi, use.names = FALSE), nrow=q, byrow=TRUE)
  for (i in 2:n) {
    if (!is.na(x[i])) {
      mean <- get_ar_mean(pn$mu, phi, x, q, i)
      p <- dnorm(x[i], mean, pn$sigma)
    }
    else {
      p <- rep(1, m)
    }
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
  
  if (hessian) {
    if (!stationary) {
      np2 <- np - m + 1
      h <- mod$hessian[1:np2, 1:np2]
    }
    else {
      np2 <- np
      h <- mod$hessian
    }
    if (det(h) != 0) {
      h <- solve(h)
      jacobian <- ar_jacobian(m, q, np2, mod$estimate, 
                                stationary = stationary)
      h <- t(jacobian) %*% h %*% jacobian
      return(list(
        m = m, q = q, mu = pn$mu, sigma = pn$sigma,
        gamma = pn$gamma, delta = pn$delta, phi = pn$phi,
        code = mod$code, mllk = mllk,
        aic = aic, bic = bic, invhessian = h
      ))
    }
    else {
      return(list(
        m = m, q = q, mu = pn$mu, sigma = pn$sigma,
        gamma = pn$gamma, delta = pn$delta, phi = pn$phi,
        code = mod$code, mllk = mllk,
        aic = aic, bic = bic, hessian = mod$hessian
      ))
    }
  }
  else {
    return(list(
      m = m, q = q, mu = pn$mu, sigma = pn$sigma,
      gamma = pn$gamma, delta = pn$delta, phi = pn$phi,
      code = mod$code, mllk = mllk, aic = aic, bic = bic
    ))
  }
}

ar_hmm_viterbi <- function(x, mod) {
  n <- length(x)
  xi <- matrix(0, n, mod$m)
  foo <- mod$delta * dnorm(x[1], mod$mu, mod$sigma)
  xi[1, ] <- foo / sum(foo)
  
  phi <- matrix(unlist(mod$phi, use.names = FALSE), nrow = mod$q, byrow = TRUE)
  for (t in 2:n) {
    mean <- get_ar_mean(mod$mu, phi, x, mod$q, t)
    p <- dnorm(x[t], mean = mean, sd = mod$sigma)
    foo <- apply(xi[t - 1, ] * mod$gamma, 2, max) * p 
    xi[t, ] <- foo / sum(foo)
  }
  
  iv <- numeric(n)
  iv[n] <- which.max(xi[n, ])
  for (t in (n - 1):1) {
    iv[t] <- which.max(mod$gamma[, iv[t + 1]] * xi[t, ])
  }
  return(data_frame(index = 1:n, state = iv))
}

ar_hmm_lforward <- function(x, mod) {
  n <- length(x)
  lalpha <- matrix(NA, mod$m, n)
  foo <- mod$delta * dnorm(x[1], mod$mu, mod$sigma)
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo / sumfoo
  lalpha[, 1] <- lscale + log(foo)
  
  phi <- matrix(unlist(mod$phi, use.names = FALSE), nrow = mod$q, byrow = TRUE)
  for (i in 2:n) {
    mean <- get_ar_mean (mod$mu, phi, x, mod$q, i)
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
  
  phi <- matrix(unlist(mod$phi, use.names = FALSE), nrow = mod$q, byrow = TRUE)
  for (i in (n - 1):1) {
    mean <- get_ar_mean (mod$mu, phi, x, mod$q, i + 1)
    p <- dnorm(x[i + 1], mean = mean, sd = mod$sigma)
    foo <- mod$gamma %*% (p * foo)
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
    
    phi <- matrix(unlist(mod$phi, use.names = FALSE), nrow = mod$q, byrow = TRUE)
    p <- matrix(NA, n, mod$m)
    for (i in 1:n) {
      mean <- get_ar_mean (mod$mu, phi, x, mod$q, i)
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
    
    phi <- matrix(unlist(mod$phi, use.names = FALSE), nrow = mod$q, byrow = TRUE)
    p <- matrix(NA, n, mod$m)
    for (i in 1:n) {
      mean <- get_ar_mean (mod$mu, phi, x, mod$q, i)
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

ar_jacobian <- function(m, q, n, parvect, stationary = TRUE) {
  pn <- ar_hmm_pw2pn(m, q, parvect, stationary)
  jacobian <- matrix(0, nrow = n, ncol = n)
  jacobian[1:m, 1:m] <- diag(m)
  jacobian[(m + 1):(2 * m), (m + 1):(2 * m)] <- diag(pn$sigma)
  count <- 0
  for (i in 1:m) {
    for (j in 1:m) {
      if (j != i) {
        count <- count + 1
        foo <- -pn$gamma[i, j] * pn$gamma[i, ]
        foo[j] <- pn$gamma[i, j] * (1 - pn$gamma[i, j])
        foo <- foo[-i]
        jacobian[2 * m + count,
                 (2 * m + (i - 1) * (m - 1) + 1):(2 * m + i * (m - 1))] <- foo
      }
    }
  }
  count <- 2 * m + count + 1
  phi <- unlist(pn$phi, use.names = FALSE)
  jacobian[count:n, count:n] <- diag(phi)
  return(jacobian)
}

ar_bootstrap_estimates <- function(mod, n, len, stationary) {
  m <- mod$m
  q <- mod$q
  mu_estimate <- numeric(n * m)
  sigma_estimate <- numeric(n * m)
  gamma_estimate <- numeric(n * m * m)
  delta_estimate <- numeric(n * m)
  phi_estimate <- numeric(n * m * q)
  for (i in 1:n) {
    sample <- ar_hmm_generate_sample(len, mod)
    mod2 <- ar_hmm_mle(sample$obs, m, q, mod$mu, mod$sigma, mod$gamma, mod$phi,
                         mod$delta, stationary = stationary, hessian = FALSE)
    mu_estimate[((i - 1) * m + 1):(i * m)] <- mod2$mu
    sigma_estimate[((i - 1) * m + 1):(i * m)] <- mod2$sigma
    gamma_estimate[((i - 1) * m * m + 1):(i * m * m)] <- mod2$gamma
    phi_estimate[((i - 1) * m * q + 1):(i * m * q)] <- 
      unlist(mod2$phi, use.names = FALSE)
    delta_estimate[((i - 1) * m + 1):(i * m)] <- mod2$delta
  }
  return(list(mu = mu_estimate,
              sigma = sigma_estimate,
              gamma = gamma_estimate, 
              phi = phi_estimate,
              delta = delta_estimate))
}

ar_bootstrap_ci <- function(mod, bootstrap, alpha) {
  m <- mod$m
  mu_lower <- rep(NA, m)
  mu_upper <- rep(NA, m)
  sigma_lower <- rep(NA, m)
  sigma_upper <- rep(NA, m)
  gamma_lower <- rep(NA, m * m)
  gamma_upper <- rep(NA, m * m)
  phi_lower <- rep(NA, m * q)
  phi_upper <- rep(NA, m * q)
  delta_lower <- rep(NA, m)
  delta_upper <- rep(NA, m)
  bootstrap1 <- data_frame(mu = bootstrap$mu,
                           sigma = bootstrap$sigma,
                           delta = bootstrap$delta)
  bootstrap2 <- data_frame(gamma = bootstrap$gamma)
  bootstrap3 <- data_frame(phi = bootstrap$phi)
  for (i in 1:m) {
    if (i == m) {
      foo <- bootstrap1 %>% filter((row_number() %% m) == 0)
    }
    else {
      foo <- bootstrap1 %>% filter((row_number() %% m) == i)
    }
    mu_lower[i] <- 2 * mod$mu[i] -
      quantile(foo$mu, 1 - (alpha / 2), names = FALSE)
    mu_upper[i] <- 2 * mod$mu[i] -
      quantile(foo$mu, alpha / 2, names = FALSE)
    sigma_lower[i] <- 2 * mod$sigma[i] -
      quantile(foo$sigma, 1 - (alpha / 2), names = FALSE)
    sigma_upper[i] <- 2 * mod$sigma[i] -
      quantile(foo$sigma, alpha / 2, names = FALSE)
    delta_lower[i] <- 2 * mod$delta[i] -
      quantile(foo$delta, 1 - (alpha / 2), names = FALSE)
    delta_upper[i] <- 2 * mod$delta[i] -
      quantile(foo$delta, alpha / 2, names = FALSE)
  }
  phi <- unlist(mod$phi, use.names = FALSE)
  for (i in 1:(m * m)) {
    if (i == (m * m)) {
      foo <- bootstrap2 %>% filter((row_number() %% (m * m)) == 0)
    }
    else {
      foo <- bootstrap2 %>% filter((row_number() %% (m * m)) == i)
    }
    
    gamma_lower[i] <- 2 * mod$gamma[i] -
      quantile(foo$gamma, 1 - (alpha / 2), names = FALSE)
    gamma_upper[i] <- 2 * mod$gamma[i] -
      quantile(foo$gamma, alpha / 2, names = FALSE)
  }
  for (i in 1:(m * q)) {
    if (i == (m * q)) {
      foo <- bootstrap3 %>% filter((row_number() %% (m * q)) == 0)
    }
    else {
      foo <- bootstrap3 %>% filter((row_number() %% (m * q)) == i)
    }
    
    phi_lower[i] <- 2 * phi[i] -
      quantile(foo$phi, 1 - (alpha / 2), names = FALSE)
    phi_upper[i] <- 2 * phi[i] -
      quantile(foo$phi, alpha / 2, names = FALSE)
  }
  return(list(
    mu_lower = mu_lower, mu_upper = mu_upper,
    sigma_lower = sigma_lower, sigma_upper = sigma_upper,
    gamma_lower = gamma_lower, gamma_upper = gamma_upper,
    phi_lower = phi_lower, phi_upper = phi_upper,
    delta_lower = delta_lower, delta_upper = delta_upper
  ))
}
