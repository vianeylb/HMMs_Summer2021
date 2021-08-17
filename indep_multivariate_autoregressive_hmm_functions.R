library(tidyverse)
library(ggplot2)
library(ggfortify)
library(plyr)
library(bayesforecast)
library(gridExtra)
library(mvtnorm)
library(parallel)
library(Rcpp)

sourceCpp("foralg.cpp")
sourceCpp("dmvnrm_arma.cpp")
#sourceCpp("dmvnrm_arma_v2.cpp")


# q: order of autoregressive model
# phi: autoregressive model parameters, as list of qxq matrices

# Get mean for given x in multivariate autoregressive series
# phi is a list of matrices
# i is the time index
get_inmar_mean <- function(mu, phi, x, m, q, k, i) {
  if (i == 1) {
    mean <- mu
  }
  else {
    mean <- list()
    if (i <= q) {
      x_lag <- x[, (i - 1):1]
      for (j in 1:m) {
        mean[[j]] <- mu[[j]] + diag(phi[[j]][, 1:(i - 1)] %*% t(x_lag))
      }
    }
    else {
      x_lag <- x[, (i - 1):(i - q)]
      for (j in 1:m) {
        mean[[j]] <- mu[[j]] + diag(phi[[j]] %*% t(x_lag))
      }
    }
  }
  return(mean)
}

get_all_inmar_means <- function(x, mod, m, q, k) {
  n <- ncol(x)
  means <- list()
  for(i in 1:m){
    mu_matrix <- matrix(mod$mu[[i]], nrow = n, ncol = k, byrow = TRUE)
    for (j in 1:k){
      x_lags <- matrix(nrow = n, ncol = q)
      for (h in 1:q){
        x_lag <- lag(x[j, ], h)
        x_lag[1:h] <- rep(0, h)
        x_lags[, h] <- x_lag
      }
      phi <- mod$phi[[i]][j, ]
      #mu_matrix[, j] <- mu_matrix[, j] + x_lags%*% phi
      mu_matrix[, j] <- mu_matrix[, j] + (x_lags - mu_matrix[,j])%*% phi
    }
    means[[i]] <- mu_matrix
  }
  return(means)
}

inmar_hmm_generate_sample <- function(ns, mod) {
  mvect <- 1:mod$m
  state <- numeric(ns)
  state[1] <- sample(mvect, 1, prob = mod$delta)
  if (ns > 1) {
    for (i in 2:ns) {
      state[i] <- sample(mvect, 1, prob = mod$gamma[state[i - 1], ])
    }
  }
  x <- matrix(nrow = mod$k, ncol = ns)
  for (i in 1:ns) {
    mean <- get_inmar_mean(mod$mu, mod$phi, x, mod$m, mod$q, mod$k, i)
    x[, i] <- rmvnorm(1, mean = mean[[state[i]]], sigma = diag(mod$sigma[[state[i]]])^2)
  }
  return(list(index = c(1:ns), state = state, obs = x))
}

inmar_hmm_pn2pw <- function(m, mu, sigma, gamma, phi,
                          delta = NULL, stationary = TRUE) {
  mu <- unlist(mu, use.names = FALSE)
  tsigma <- log(unlist(sigma, use.names = FALSE))
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

inmar_hmm_pw2pn <- function(m, q, k, parvect, stationary = TRUE) {
  mu <- list()
  count <- 1
  for (i in 1:m) {
    mu[[i]] <- parvect[count:(i * k)]
    count <- count + k
  }
  
  sigma <- list()
  for (i in 1:m) {
    sigma[[i]] <- exp(parvect[count:(count + k - 1)])
    count <- count + k
  }
  
  tgamma <- parvect[count:(count + m * (m - 1) - 1)]
  count <- count + m * (m - 1)
  gamma <- diag(m)
  gamma[!gamma] <- exp(tgamma)
  gamma <- gamma / apply(gamma, 1, sum)
  
  tphi <- parvect[(count:(count + k * q * m))]
  count <- count + k * q * m
  phi <- list()
  for (i in 1:m) {
    foo <- tphi[((i - 1) * k * q + 1):(i * k * q)]
    phi[[i]] <- matrix(foo, nrow = k)
  }
  
  if (stationary) {
    delta <- solve(t(diag(m) - gamma + 1), rep(1, m))
  }
  else {
    tdelta <- parvect[count:(count + m - 2)]
    foo <- c(1, exp(tdelta))
    delta <- foo / sum(foo)
  }
  return(list(mu = mu, sigma = sigma, gamma = gamma, phi = phi, delta = delta))
}

inmar_hmm_mllk <- function(parvect, x, m, q, k, stationary = TRUE, state = NULL) {
  n <- ncol(x)
  pn <- inmar_hmm_pw2pn(m, q, k, parvect, stationary = stationary)
  
  if (is.null(state)){
    p <- inmar_densities(x, pn, m, q, k, n)
  } else {
    p <- inmar_densities_labelled(x, pn, m, q, k, n, state)
  }
  
  foo <- matrix(pn$delta, ncol = m)
  lscale <- foralg(n, m, foo, pn$gamma, p)
  mllk <- -lscale
  return(mllk)
}

# Returns n * m matrix of state dependent probability densities
inmar_densities <- function(x, mod, m, q, k, n) {
  p <- matrix(1, nrow = n, ncol = m)
  means <- get_all_inmar_means(x, mod, m, q, k)
  for (i in 1:n) {
    for (j in 1:m) {
      for (l in 1:k){
        p[i, j] <- p[i, j] * dnorm(x[l, i], means[[j]][i, l], mod$sigma[[j]][l])
      }
    }
  }
  return(p)
}

inmar_densities_labelled <- function(x, mod, m, q, k, n, state) {
  p <- matrix(1, nrow = n, ncol = m)
  means <- get_all_inmar_means(x, mod, m, q, k)
  for (i in 1:n) {
    for (j in 1:m) {
      ## check if state is known
      if(state[i] == 0){
        for (l in 1:k){
          p[i, j] <- p[i, j] * dnorm(x[l, i], means[[j]][i, l], mod$sigma[[j]][l])
        }
      } else {
          if (j == state[i]){
            for (l in 1:k){
              p[i, j] <- p[i, j] * dnorm(x[l, i], means[[j]][i, l], mod$sigma[[j]][l])
            }
          } else{
            p[i, j] <- 0
          }
      }
      
    }
  }
  return(p)
}

# Computing MLE from natural parameters
inmar_hmm_mle <- function(x, m, q, k, mu0, sigma0, gamma0, phi0, delta0 = NULL,
                        stationary = TRUE, hessian = FALSE,
                        steptol = 1e-6, iterlim = 100,
                        stepmax = 100, state = NULL) {
  parvect0 <- inmar_hmm_pn2pw(m, mu0, sigma0, gamma0, phi0, delta0,
                            stationary = stationary
  )
  mod <- nlm(inmar_hmm_mllk, parvect0,
             x = x, m = m, q = q, k = k,
             stationary = stationary, hessian = hessian, state = state, 
             steptol = steptol, stepmax = stepmax, iterlim = iterlim, 
             print.level = 2
  )
  pn <- inmar_hmm_pw2pn(
    m = m, q = q, k = k, parvect = mod$estimate,
    stationary = stationary
  )
  mllk <- mod$minimum
  
  np <- length(parvect0)
  aic <- 2 * (mllk + np)
  n <- sum(!is.na(x))
  bic <- 2 * mllk + np * log(n)
  
  if (hessian) {
    return(list(
      m = m, q = q, k = k, mu = pn$mu, sigma = pn$sigma,
      gamma = pn$gamma, phi = pn$phi, delta = pn$delta,
      code = mod$code, mllk = mllk,
      aic = aic, bic = bic, hessian = mod$hessian, np = np
    ))
  }
  else {
    return(list(
      m = m, q = q, k = k, mu = pn$mu, sigma = pn$sigma,
      gamma = pn$gamma, phi = pn$phi, delta = pn$delta,
      code = mod$code, mllk = mllk, aic = aic, bic = bic
    ))
  }
}

inmar_hmm_viterbi <- function(x, mod) {
  n <- ncol(x)
  xi <- matrix(0, n, mod$m)
  p <- inmar_densities(x, mod, mod$m, mod$q, mod$k, n)
  foo <- mod$delta * p[1, ]
  xi[1, ] <- foo / sum(foo)
  for (t in 2:n) {
    foo <- apply(xi[t - 1, ] * mod$gamma, 2, max) * p[t, ]
    xi[t, ] <- foo / sum(foo)
  }
  iv <- numeric(n)
  iv[n] <- which.max(xi[n, ])
  for (t in (n - 1):1) {
    iv[t] <- which.max(mod$gamma[, iv[t + 1]] * xi[t, ])
  }
  return(data_frame(index = 1:n, state = iv))
}

inmar_hmm_lforward <- function(x, mod) {
  n <- ncol(x)
  lalpha <- matrix(NA, mod$m, n)
  p <- inmar_densities(x, mod, mod$m, mod$q, mod$k, n)
  foo <- mod$delta * p[1, ]
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo / sumfoo
  lalpha[, 1] <- lscale + log(foo)
  for (i in 2:n) {
    foo <- foo %*% mod$gamma * p[i, ]
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
    lalpha[, i] <- log(foo) + lscale
  }
  return(lalpha)
}

inmar_hmm_lbackward <- function(x, mod) {
  n <- ncol(x)
  m <- mod$m
  p <- inmar_densities(x, mod, mod$m, mod$q, mod$k, n)
  lbeta <- matrix(NA, m, n)
  lbeta[, n] <- rep(0, m)
  foo <- rep(1 / m, m)
  lscale <- log(m)
  for (i in (n - 1):1) {
    foo <- mod$gamma %*% (p[i + 1, ] * foo)
    lbeta[, i] <- log(foo) + lscale
    sumfoo <- sum(foo)
    foo <- foo / sumfoo
    lscale <- lscale + log(sumfoo)
  }
  return(lbeta)
}

inmar_hmm_pseudo_residuals <- function(x, mod, type, stationary) {
  if (stationary) {
    delta <- solve(t(diag(mod$m) - mod$gamma + 1), rep(1, mod$m))
  }
  else {
    delta <- mod$delta
  }
  if (type == "ordinary") {
    n <- ncol(x)
    la <- inmar_hmm_lforward(x, mod)
    lb <- inmar_hmm_lbackward(x, mod)
    lafact <- apply(la, 2, max)
    lbfact <- apply(lb, 2, max)
    p <- inmar_dist_mat(x, mod, n)
    npsr <- rep(NA, n)
    npsr[1] <- qnorm(delta %*% p[1, ])
    for (i in 2:n) {
      a <- exp(la[, i - 1] - lafact[i])
      b <- exp(lb[, i] - lbfact[i])
      foo <- (a %*% mod$gamma) * b
      foo <- foo / sum(foo)
      npsr[i] <- qnorm(foo %*% p[i, ])
    }
    return(data_frame(npsr, index = c(1:n)))
  }
  else if (type == "forecast") {
    n <- ncol(x)
    la <- inmar_hmm_lforward(x, mod)
    p <- inmar_dist_mat(x, mod, n)
    npsr <- rep(NA, n)
    npsr[1] <- qnorm(delta %*% p[1, ])
    for (i in 2:n) {
      la_max <- max(la[, i - 1])
      a <- exp(la[, i - 1] - la_max)
      npsr[i] <- qnorm(t(a) %*% (mod$gamma / sum(a)) %*% p[i, ])
    }
    return(data_frame(npsr, index = c(1:n)))
  }
}

# Get multivariate normal distribution given mod and x
# Returns n * m matrix
inmar_dist_mat <- function(x, mod, n) {
  p <- matrix(NA, n, mod$m)
  means <- get_all_inmar_means(x, mod, mod$m, mod$q, mod$k)
  for (i in 1:n) {
    for (j in 1:mod$m) {
      p[i, j] <- pmvnorm(
        lower = rep(-Inf, mod$k), upper = x[, i],
        mean = means[[j]][i, ], sigma = diag(mod$sigma[[j]]^2)
      )
    }
  }
  return(p)
}

inmar_inv_hessian <- function(mod, stationary = TRUE){
  if (!stationary) {
    np2 <- mod$np - mod$m + 1
    h <- mod$hessian[1:np2, 1:np2]
  }
  else {
    np2 <- mod$np
    h <- mod$hessian
  }
  h <- solve(h)
  jacobian <- inmar_jacobian(mod, np2)
  h <- t(jacobian) %*% h %*% jacobian
  return(h)
}

inmar_jacobian <- function(mod, n) {
  m <- mod$m
  q <- mod$q
  k <- mod$k
  
  jacobian <- matrix(0, nrow = n, ncol = n)
  jacobian[1:(m * k), 1:(m * k)] <- diag(m * k)
  
  sigma <- unlist(mod$sigma, use.names = FALSE)
  jacobian[(m * k + 1):(2 * m * k), (m * k + 1):(2 * m * k)] <- diag(sigma)
  rowcount <- 2 * m * k + 1  
  colcount <- rowcount

  for (i in 1:m) {
    for (j in 1:m) {
      if (j != i) {
        foo <- -mod$gamma[i, j] * mod$gamma[i, ]
        foo[j] <- mod$gamma[i, j] * (1 - mod$gamma[i, j])
        foo <- foo[-i]
        jacobian[rowcount, colcount:(colcount + m - 2)] <- foo
        rowcount <- rowcount + 1
      }
    }
    colcount <- colcount + m - 1
  }
  
  phi <- unlist(mod$phi, use.names = FALSE)
  jacobian[rowcount:n, colcount:n] <- diag(phi)
  
  return(jacobian)
}

inmar_bootstrap_estimates <- function(mod, n, len, stationary) {
  m <- mod$m
  k <- mod$k
  q <- mod$q
  mu_estimate <- numeric(n * m * k)
  sigma_estimate <- numeric(n * m * k)
  gamma_estimate <- numeric(n * m * m)
  phi_estimate <- numeric(n * m * k * k * q)
  delta_estimate <- numeric(n * m)
  for (i in 1:n) {
    sample <- inmar_hmm_generate_sample(len, mod)
    mod2 <- inmar_hmm_mle(sample$obs, m, q, k, mod$mu, mod$sigma,
                        mod$gamma, mod$phi, mod$delta,
                        stationary = stationary
    )
    mu_estimate[((i - 1) * m * k + 1):(i * m * k)] <-
      unlist(mod2$mu, use.names = FALSE)
    sigma_estimate[((i - 1) * m * k + 1):(i * m * k)] <-
      unlist(mod2$sigma, use.names = FALSE)
    gamma_estimate[((i - 1) * m * m + 1):(i * m * m)] <- mod2$gamma
    phi_estimate[((i - 1) * m * k * q + 1):(i * m * k * q)] <-
      unlist(mod2$phi, use.names = FALSE)
    delta_estimate[((i - 1) * m + 1):(i * m)] <- mod2$delta
  }
  return(list(
    mu = mu_estimate,
    sigma = sigma_estimate,
    gamma = gamma_estimate,
    phi = phi_estimate,
    delta = delta_estimate
  ))
}

inmar_bootstrap_ci <- function(mod, bootstrap, alpha) {
  m <- mod$m
  k <- mod$k
  
  mu_lower <- matrix(NA, m, k)
  mu_upper <- matrix(NA, m, k)
  bootstrap_mu <- data_frame(mu = bootstrap$mu)
  mu <- unlist(mod$mu, use.names = FALSE)
  for (i in 1:m) {
    for (j in 1:k) {
      if (i == m & j == k) {
        foo <- bootstrap_mu %>%
          filter((row_number() %% (m * k)) == 0)
      }
      else {
        foo <- bootstrap_mu %>%
          filter((row_number() %% (m * k)) == (i - 1) * k + j)
      }
      mu_lower[i, j] <- 2 * mu[(i - 1) * k + j] -
        quantile(foo$mu, 1 - (alpha / 2), names = FALSE)
      mu_upper[i, j] <- 2 * mu[(i - 1) * k + j] -
        quantile(foo$mu, alpha / 2, names = FALSE)
    }
  }
  
  sigma_lower <- matrix(NA, m, k)
  sigma_upper <- matrix(NA, m, k)
  bootstrap_sigma <- data_frame(sigma = bootstrap$sigma)
  sigma <- unlist(mod$sigma, use.names = FALSE)
  for (i in 1:m) {
    for (j in 1:k) {
      if (i == m & j == k) {
        foo <- bootstrap_sigma %>%
          filter((row_number() %% (m * k)) == 0)
      }
      else {
        foo <- bootstrap_sigma %>%
          filter((row_number() %% (m * k)) == (i - 1) * k + j)
      }
      sigma_lower[i, j] <- 2 * sigma[(i - 1) * k + j] -
        quantile(foo$sigma, 1 - (alpha / 2), names = FALSE)
      sigma_upper[i, j] <- 2 * sigma[(i - 1) * k + j] -
        quantile(foo$sigma, alpha / 2, names = FALSE)
    }
  }
  
  gamma_lower <- rep(NA, m * m)
  gamma_upper <- rep(NA, m * m)
  bootstrap_gamma <- data_frame(gamma = bootstrap$gamma)
  gamma <- mod$gamma
  for (i in 1:(m * m)) {
    if (i == (m * m)) {
      foo <- bootstrap_gamma %>%
        filter((row_number() %% (m * m)) == 0)
    }
    else {
      foo <- bootstrap_gamma %>%
        filter((row_number() %% (m * m)) == i)
    }
    gamma_lower[i] <- 2 * gamma[i] -
      quantile(foo$gamma, 1 - (alpha / 2), names = FALSE)
    gamma_upper[i] <- 2 * gamma[i] -
      quantile(foo$gamma, alpha / 2, names = FALSE)
  }
  
  phi_lower <- matrix(NA, nrow = m * k, ncol = q)
  phi_upper <- matrix(NA, nrow = m * k, ncol = q)
  bootstrap_phi <- data_frame(phi = bootstrap$phi)
  phi <- unlist(mod$phi, use.names = FALSE)
  for (i in 1:(m * k * q)) {
    if (i == (m * k * q)) {
      foo <- bootstrap_phi %>%
        filter((row_number() %% (m * k * q)) == 0)
    }
    else {
      foo <- bootstrap_phi %>%
        filter((row_number() %% (m * k * q)) == i)
    }
    phi_lower[i] <- 2 * phi[i] -
      quantile(foo$phi, 1 - (alpha / 2), names = FALSE)
    phi_upper[i] <- 2 * phi[i] -
      quantile(foo$phi, alpha / 2, names = FALSE)
  }
  
  delta_lower <- rep(NA, m)
  delta_upper <- rep(NA, m)
  bootstrap_delta <- data_frame(delta = bootstrap$delta)
  delta <- mod$delta
  for (i in 1:m) {
    if (i == m) {
      foo <- bootstrap_delta %>% filter((row_number() %% m) == 0)
    }
    else {
      foo <- bootstrap_delta %>% filter((row_number() %% m) == i)
    }
    delta_lower[i] <- 2 * delta[i] -
      quantile(foo$delta, 1 - (alpha / 2), names = FALSE)
    delta_upper[i] <- 2 * delta[i] -
      quantile(foo$delta, alpha / 2, names = FALSE)
  }
  
  return(list(
    mu_lower = mu_lower, mu_upper = mu_upper,
    sigma_lower = sigma_lower, sigma_upper = sigma_upper,
    gamma_lower = gamma_lower, gamma_upper = gamma_upper,
    phi_lower = phi_lower, phi_upper = phi_upper,
    delta_lower = delta_lower, delta_upper = delta_upper
  ))
}
