library(tidyverse)
library(ggplot2)
library(ggfortify)
library(plyr)
library(mvtnorm)
library(parallel)
library(Rcpp)

setwd("C:/Jessica/UofT Y4/Research/Coding")
sourceCpp("foralg.cpp")
sourceCpp("dmvnrm_arma.cpp")


# Changes for multivariate normal
# mu is a list of vectors
# sigma is a list of matrices
# parvect is a list
# k is the number of variables

# Transform normal natural parameters to working parameters
inmvnorm_hmm_pn2pw <- function(m, mu, sigma, gamma,
                             delta = NULL, stationary = TRUE) {
  tmu <- unlist(mu, use.names = FALSE)
  tsigma <- log(unlist(sigma, use.names = FALSE))
  foo <- log(gamma / diag(gamma))
  tgamma <- as.vector(foo[!diag(m)])
  if (stationary) {
    tdelta <- NULL
  }
  else {
    tdelta <- log(delta[-1] / delta[1])
  }
  parvect <- c(tmu, tsigma, tgamma, tdelta)
  return(parvect)
}


# Transform normal working parameters to natural parameters
inmvnorm_hmm_pw2pn <- function(m, k, parvect, stationary = TRUE) {
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
  
  if (stationary) {
    delta <- solve(t(diag(m) - gamma + 1), rep(1, m))
  }
  else {
    tdelta <- parvect[count:(count + m - 2)]
    foo <- c(1, exp(tdelta))
    delta <- foo / sum(foo)
  }
  return(list(mu = mu, sigma = sigma, gamma = gamma, delta = delta))
}

# Computing minus the log-likelihood from the working parameters
inmvnorm_hmm_mllk <- function(parvect, x, m, k, stationary = TRUE) {
  n <- ncol(x)
  pn <- inmvnorm_hmm_pw2pn(m, k, parvect, stationary = stationary)
  p <- inmvnorm_densities(x, pn, m, k, n)
  foo <- matrix(pn$delta, ncol = m)
  lscale <- foralg(n, m, foo, pn$gamma, p)
  mllk <- -lscale
  return(mllk)
}

# Returns n * m matrix of state dependent probability densities
inmvnorm_densities <- function(x, mod, m, k, n) {
  p <- matrix(1, nrow = n, ncol = m)
  for (i in 1:n) {
    for (j in 1:m) {
      for (l in 1:k){
        p[i, j] <- p[i, j] * dnorm(x[l, i], mod$mu[[j]][l], mod$sigma[[j]][l])
      }
    }
  }
  return(p)
}

# Computing MLE from natural parameters
inmvnorm_hmm_mle <- function(x, m, k, mu0, sigma0, gamma0, delta0 = NULL,
                           stationary = TRUE, hessian = FALSE) {
  parvect0 <- inmvnorm_hmm_pn2pw(
    m = m, mu = mu0, sigma = sigma0,
    gamma = gamma0, delta = delta0,
    stationary = stationary
  )
  mod <- nlm(inmvnorm_hmm_mllk, parvect0,
             x = x, m = m, k = k,
             stationary = stationary, hessian = hessian
  )
  pn <- inmvnorm_hmm_pw2pn(
    m = m, k = k, parvect = mod$estimate,
    stationary = stationary
  )
  mllk <- mod$minimum
  
  np <- length(parvect0)
  aic <- 2 * (mllk + np)
  n <- sum(!is.na(x))
  bic <- 2 * mllk + np * log(n)
  
  if (hessian) {
    return(list(
      m = m, k = k, mu = pn$mu, sigma = pn$sigma,
      gamma = pn$gamma, delta = pn$delta,
      code = mod$code, mllk = mllk,
      aic = aic, bic = bic, hessian = mod$hessian, np = np
    ))
  }
  else {
    return(list(
      m = m, k = k, mu = pn$mu, sigma = pn$sigma,
      gamma = pn$gamma, delta = pn$delta,
      code = mod$code, mllk = mllk, aic = aic, bic = bic
    ))
  }
}

# Generating a sample from normal distributions
# x is a m x ns matrix
inmvnorm_hmm_generate_sample <- function(ns, mod) {
  mvect <- 1:mod$m
  state <- numeric(ns)
  state[1] <- sample(mvect, 1, prob = mod$delta)
  if (ns > 1) {
    for (i in 2:ns) {
      state[i] <- sample(mvect, 1, prob = mod$gamma[state[i - 1], ])
    }
  }
  x <- sapply(state, inmvnorm_hmm_sample_one, mod = mod)
  return(list(index = c(1:ns), state = state, obs = x))
}

inmvnorm_hmm_sample_one <- function(state, mod) {
  x <- rmvnorm(1, mean = mod$mu[[state]], sigma = diag(mod$sigma[[state]]^2))
  return(x)
}

# Global decoding by the Viterbi algorithm
inmvnorm_hmm_viterbi <- function(x, mod) {
  n <- ncol(x)
  xi <- matrix(0, n, mod$m)
  p <- inmvnorm_densities(x, mod, mod$m, mod$k, n)
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

# Computing log(forward probabilities) for normal distribution
inmvnorm_hmm_lforward <- function(x, mod) {
  n <- ncol(x)
  lalpha <- matrix(NA, mod$m, n)
  p <- inmvnorm_densities(x, mod, mod$m, mod$k, n)
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

# Computing log(backward probabilities) for normal distribution
inmvnorm_hmm_lbackward <- function(x, mod) {
  n <- ncol(x)
  m <- mod$m
  p <- inmvnorm_densities(x, mod, mod$m, mod$k, n)
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

# Normal pseudo-residuals for Normal HMM
# Type can be "ordinary" or "forecast"
inmvnorm_hmm_pseudo_residuals <- function(x, mod, type, stationary = TRUE) {
  if (stationary) {
    delta <- solve(t(diag(mod$m) - mod$gamma + 1), rep(1, mod$m))
  }
  else {
    delta <- mod$delta
  }
  if (type == "ordinary") {
    n <- ncol(x)
    la <- inmvnorm_hmm_lforward(x, mod)
    lb <- inmvnorm_hmm_lbackward(x, mod)
    lafact <- apply(la, 2, max)
    lbfact <- apply(lb, 2, max)
    p <- inmvnorm_dist_mat(x, mod, n)
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
    la <- inmvnorm_hmm_lforward(x, mod)
    p <- inmvnorm_dist_mat(x, mod, n)
    npsr <- rep(NA, n)
    npsr[1] <- qnorm(delta %*% p[1, ])
    for (i in 2:n) {
      la_max <- max(la[, i - 1])
      a <- exp(la[, i - 1] - la_max)
      npsr[i] <- qnorm(t(a) %*% (gamma / sum(a)) %*% p[i, ])
    }
    return(data_frame(npsr, index = c(1:n)))
  }
}

# Get multivariate normal distribution given mod and x
inmvnorm_dist_mat <- function(x, mod, n) {
  p <- matrix(NA, n, mod$m)
  for (i in 1:n) {
    for (j in 1:m) {
      p[i, j] <- pmvnorm(
        lower = rep(-Inf, mod$k), upper = x[, i],
        mean = mod$mu[[j]], sigma = diag(mod$sigma[[j]]^2)
      )
    }
  }
  return(p)
}

inmvnorm_inv_hessian <- function(mod, stationary = TRUE){
  if (!stationary) {
    np2 <- mod$np - mod$m + 1
    h <- mod$hessian[1:np2, 1:np2]
  }
  else {
    np2 <- mod$np
    h <- mod$hessian
  }
  h <- solve(h)
  jacobian <- inmvnorm_jacobian(mod, np2)
  h <- t(jacobian) %*% h %*% jacobian
  return(h)
}


# Jacobian matrix for parameters
# n should be total number of parameters estimated, excluding delta
inmvnorm_jacobian <- function(mod, n) {
  m <- mod$m
  k <- mod$k
  jacobian <- matrix(0, nrow = n, ncol = n)
  # Jacobian for mu only is a m*k identity matrix
  jacobian[1:(m * k), 1:(m * k)] <- diag(m * k)
  sigma <- unlist(mod$sigma, use.names = FALSE)
  jacobian[(m * k + 1):(2 * m * k), (m * k + 1):(2 * m * k)] <- diag(sigma)
  # count is the row at which the current parameter's derivatives are placed
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
  return(jacobian)
}

# Bootstrapping estimates
inmvnorm_bootstrap_estimates <- function(mod, n, k, len, stationary) {
  m <- mod$m
  mu_estimate <- numeric(n * m * k)
  sigma_estimate <- numeric(n * m * k)
  gamma_estimate <- numeric(n * m * m)
  delta_estimate <- numeric(n * m)
  for (i in 1:n) {
    sample <- inmvnorm_hmm_generate_sample(len, mod)
    mod2 <- inmvnorm_hmm_mle(sample$obs, m, k, mod$mu, mod$sigma,
                           mod$gamma, mod$delta,
                           stationary = stationary
    )
    mu_estimate[((i - 1) * m * k + 1):(i * m * k)] <-
      unlist(mod2$mu, use.names = FALSE)
    sigma_estimate[((i - 1) * m * k + 1):(i * m * k)] <-
      unlist(mod2$sigma, use.names = FALSE)
    gamma_estimate[((i - 1) * m * m + 1):(i * m * m)] <- mod2$gamma
    delta_estimate[((i - 1) * m + 1):(i * m)] <- mod2$delta
  }
  return(list(
    mu = mu_estimate, sigma = sigma_estimate,
    gamma = gamma_estimate, delta = delta_estimate
  ))
}

inmvnorm_bootstrap_ci <- function(mod, bootstrap, alpha, m, k) {
  # Rows of matrix are states, columns are variables
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
    delta_lower = delta_lower, delta_upper = delta_upper
  ))
}
