library(tidyverse)
library(ggplot2)
library(ggfortify)
library(plyr)
library(mvtnorm)

setwd("C:/Jessica/UofT Y4/Research/Coding")

# Changes for multivariate normal
# mu is a list of vectors
# sigma is a list of matrices
# parvect is a list
# k is the number of variables

# Transform normal natural parameters to working parameters
mvnorm_hmm_pn2pw <- function(m, mu, sigma, gamma,
                             delta = NULL, stationary = TRUE) {
  # Put all means into one vector
  mu <- unlist(mu, use.names = FALSE)
  # Only need to transform diagonal elements of sigma
  # Include only lower triangle & diag of matrix,
  # since covariance matrix must be symmetric
  tsigma <- lapply(sigma, diag_log_lower)
  tsigma <- unlist(tsigma, use.names = FALSE)
  foo <- log(gamma / diag(gamma))
  tgamma <- as.vector(foo[!diag(m)])
  if (stationary) {
    tdelta <- NULL
  }
  else {
    tdelta <- log(delta[-1] / delta[1])
  }
  # parvect is one vector of values
  parvect <- c(mu, tsigma, tgamma, tdelta)
  return(parvect)
}

# Applies log to only diagonal elements of matrix,
# then returns only lower triangular entries of matrix as vector
# going down column by column
diag_log_lower <- function(mat) {
  diag(mat) <- log(diag(mat))
  vect <- mat[lower.tri(mat, diag = TRUE)]
  return(vect)
}

  # Applies exp to only diagonal elements of matrix
diag_exp <- function(mat) {
  diag(mat) <- exp(diag(mat))
  return(mat)
}

# Get nth triangular number (0 is 0th num)
triangular_num <- function(n) {
  nums <- choose(seq(n + 1), 2)
  return(nums[n + 1])
}

# Transform normal working parameters to natural parameters
mvnorm_hmm_pw2pn <- function(m, k, parvect, stationary = TRUE) {
  # Change mu to list of vectors format
  mu <- list()
  count <- 1
  for (i in 1:m) {
    mu[[i]] <- parvect[count:(i * k)]
    count <- count + k
  }

  # Change sigma to list of matrices format
  tsigma <- list()
  # Get number of elements in lower triangle (including diag) of matrix
  t <- triangular_num(k)
  for (i in 1:m) {
    tsigma_vals <- parvect[count:(count + t - 1)]
    foo <- diag(k)
    foo[lower.tri(foo, diag = TRUE)] <- tsigma_vals
    foo <- t(foo)
    foo[lower.tri(foo, diag = TRUE)] <- tsigma_vals
    tsigma[[i]] <- foo
    count <- count + t
  }
  sigma <- lapply(tsigma, diag_exp)

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
mvnorm_hmm_mllk <- function(parvect, x, m, k, stationary = TRUE, ...) {
  n <- ncol(x)
  pn <- mvnorm_hmm_pw2pn(m, k, parvect, stationary = stationary)
  p <- mvnorm_densities(x[, 1], pn, m)
  foo <- pn$delta * p
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo / sumfoo
  for (i in 2:n) {
    p <- mvnorm_densities(x[, i], pn, m)
    foo <- foo %*% pn$gamma * p
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
  }
  mllk <- -lscale
  return(mllk)
}

# Get state dependent probability densities given x and mod
mvnorm_densities <- function(x, mod, m) {
  pvect <- numeric(m)
  for (i in 1:m) {
    pvect[i] <- dmvnorm(x, mod$mu[[i]], mod$sigma[[i]])
  }
  return(pvect)
}

# Computing MLE from natural parameters
mvnorm_hmm_mle <- function(x, m, k, mu0, sigma0, gamma0, delta0 = NULL,
                           stationary = TRUE, hessian = FALSE, ...) {
  parvect0 <- mvnorm_hmm_pn2pw(m = m, mu = mu0, sigma = sigma0,
                               gamma = gamma0, delta = delta0,
                               stationary = stationary)
  mod <- nlm(mvnorm_hmm_mllk, parvect0, x = x, m = m, k = k,
             stationary = stationary, hessian = hessian)
  pn <- mvnorm_hmm_pw2pn(m = m, k = k, parvect = mod$estimate,
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
      jacobian <- mvnorm_jacobian(m, k, n = np2, mod$estimate,
                                  stationary = stationary)
      h <- t(jacobian) %*% h %*% jacobian
      return(list(
        m = m, k = k, mu = pn$mu, sigma = pn$sigma,
        gamma = pn$gamma, delta = pn$delta,
        code = mod$code, mllk = mllk,
        aic = aic, bic = bic, invhessian = h
      ))
    }
    else {
      return(list(
        m = m, k = k, mu = pn$mu, sigma = pn$sigma,
        gamma = pn$gamma, delta = pn$delta,
        code = mod$code, mllk = mllk,
        aic = aic, bic = bic, hessian = mod$hessian
      ))
    }
  }
  else {
    return(list(
      m = m, k = k, mu = pn$mu, sigma = pn$sigma, gamma = pn$gamma, delta = pn$delta,
      code = mod$code, mllk = mllk, aic = aic, bic = bic
    ))
  }
}

# Generating a sample from normal distributions
# x is a m x ns matrix
mvnorm_hmm_generate_sample <- function(ns, mod) {
  mvect <- 1:mod$m
  state <- numeric(ns)
  state[1] <- sample(mvect, 1, prob = mod$delta)
  if (ns > 1) {
    for (i in 2:ns) {
      state[i] <- sample(mvect, 1, prob = mod$gamma[state[i - 1], ])
    }
  }
  x <- sapply(state, mvnorm_hmm_sample_one, mod = mod)
  return(list(index = c(1:ns), state = state, obs = x))
}

mvnorm_hmm_sample_one <- function(state, mod) {
  x <- rmvnorm(1, mean = mod$mu[[state]], sigma = mod$sigma[[state]])
  return(x)
}

# Global decoding by the Viterbi algorithm
mvnorm_hmm_viterbi <- function(x, mod) {
  n <- ncol(x)
  xi <- matrix(0, n, mod$m)
  p <- mvnorm_densities(x[, 1], mod, mod$m)
  foo <- mod$delta * p
  xi[1, ] <- foo / sum(foo)
  for (t in 2:n) {
    p <- mvnorm_densities(x[, t], mod, mod$m)
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

# Computing log(forward probabilities) for normal distribution
mvnorm_hmm_lforward <- function(x, mod) {
  n <- ncol(x)
  lalpha <- matrix(NA, mod$m, n)
  foo <- mod$delta * mvnorm_densities(x[, 1], mod, mod$m)
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo / sumfoo
  lalpha[, 1] <- lscale + log(foo)
  for (i in 2:n) {
    foo <- foo %*% mod$gamma * mvnorm_densities(x[, i], mod, mod$m)
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
    lalpha[, i] <- log(foo) + lscale
  }
  return(lalpha)
}

# Computing log(backward probabilities) for normal distribution
mvnorm_hmm_lbackward <- function(x, mod) {
  n <- ncol(x)
  m <- mod$m
  lbeta <- matrix(NA, m, n)
  lbeta[, n] <- rep(0, m)
  foo <- rep(1 / m, m)
  lscale <- log(m)
  for (i in (n - 1):1) {
    foo <- mod$gamma %*% (mvnorm_densities(x[, i + 1], mod, mod$m) * foo)
    lbeta[, i] <- log(foo) + lscale
    sumfoo <- sum(foo)
    foo <- foo / sumfoo
    lscale <- lscale + log(sumfoo)
  }
  return(lbeta)
}


# Normal pseudo-residuals for Normal HMM
# Type can be "ordinary" or "forecast"
mvnorm_hmm_pseudo_residuals <- function(x, mod, type, stationary = TRUE) {
  if (stationary) {
    delta <- solve(t(diag(mod$m) - mod$gamma + 1), rep(1, mod$m))
  }
  else {
    delta <- mod$delta
  }
  if (type == "ordinary") {
    n <- ncol(x)
    la <- mvnorm_hmm_lforward(x, mod)
    lb <- mvnorm_hmm_lbackward(x, mod)
    lafact <- apply(la, 2, max)
    lbfact <- apply(lb, 2, max)
    p <- mvnorm_dist_mat(x, mod)
    npsr <- rep(NA, n)
    npsr[1] <- qnorm(mod$delta %*% p[1, ])
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
    la <- mvnorm_hmm_lforward(x, mod)
    p <- mvnorm_dist_mat(x, mod)
    npsr <- rep(NA, n)
    npsr[1] <- qnorm(mod$delta %*% p[1, ])
    for (i in 2:n) {
      la_max <- max(la[, i - 1])
      a <- exp(la[, i - 1] - la_max)
      npsr[i] <- qnorm(t(a) %*% (gamma / sum(a)) %*% p[i, ])
    }
    return(data_frame(npsr, index = c(1:n)))
  }
}

# Get multivariate normal distribution given mod and x
mvnorm_dist_mat <- function(x, mod) {
  p <- matrix(NA, n, mod$m)
  for (i in 1:n) {
    for (j in 1:m) {
      p[i, j] <- pmvnorm(lower = rep(-Inf, mod$k), upper = x[, i],
                        mean = mod$mu[[j]], sigma = mod$sigma[[j]])
    }
  }
  return(p)
}

# Jacobian matrix for parameters
# n should be total number of parameters estimated, excluding delta
mvnorm_jacobian <- function(m, k, n, parvect, stationary = TRUE) {
  pn <- mvnorm_hmm_pw2pn(m, k, parvect, stationary)
  jacobian <- matrix(0, nrow = n, ncol = n)
  # Jacobian for mu only is a m*k identity matrix
  jacobian[1:(m * k), 1:(m * k)] <- diag(m * k)
  # count is the row at which the current parameter's derivatives are placed
  rowcount <- m * k + 1
  # There are t*m sigma parameters
  t <- triangular_num(k)
  for (i in 1:m) {
    sigma <- pn$sigma[[i]]
    sigma[lower.tri(sigma, diag = FALSE)] <-
      rep(1, length(sigma[lower.tri(sigma, diag = FALSE)]))
    sigma <- sigma[lower.tri(sigma, diag = TRUE)]
    jacobian[rowcount:(rowcount + t - 1),
             rowcount:(rowcount + t - 1)] <- diag(sigma)
    rowcount <- rowcount + t
  }
  colcount <- rowcount
  for (i in 1:m) {
    for (j in 1:m) {
      if (j != i) {
        foo <- -pn$gamma[i, j] * pn$gamma[i, ]
        foo[j] <- pn$gamma[i, j] * (1 - pn$gamma[i, j])
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
mvnorm_bootstrap_estimates <- function(mod, n, k, len, stationary) {
  m <- mod$m
  mu_estimate <- numeric(n * m * k)
  sigma_estimate <- numeric(n * m * k * k)
  gamma_estimate <- numeric(n * m * m)
  delta_estimate <- numeric(n * m)
  for (i in 1:n) {
    sample <- mvnorm_hmm_generate_sample(len, mod)
    mod2 <- mvnorm_hmm_mle(sample$obs, m, k, mod$mu, mod$sigma,
                           mod$gamma, mod$delta, stationary = stationary)
    mu_estimate[((i - 1) * m * k + 1):(i * m * k)] <-
      unlist(mod2$mu, use.names = FALSE)
    sigma_estimate[((i - 1) * m * k * k + 1):(i * m * k * k)] <-
      unlist(mod2$sigma, use.names = FALSE)
    gamma_estimate[((i - 1) * m * m + 1):(i * m * m)] <- mod2$gamma
    delta_estimate[((i - 1) * m + 1):(i * m)] <- mod2$delta
  }
  return(list(mu = mu_estimate, sigma = sigma_estimate,
              gamma = gamma_estimate, delta = delta_estimate))
}

mvnorm_bootstrap_ci <- function(mod, bootstrap, alpha, m, k) {
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

  # Only want lower triangle of each sigma matrix, since is symmetric
  t <- triangular_num(k)
  mat <- matrix(c(1:(k * k)), k)
  tvect <- mat[lower.tri(mat, diag = TRUE)]
  sigma_lower <- matrix(NA, 3, t)
  sigma_upper <- matrix(NA, 3, t)
  bootstrap_sigma <- data_frame(sigma = bootstrap$sigma)
  sigma <- unlist(mod$sigma, use.names = FALSE)
  for (i in 1:m) {
    for (j in 1:t) {
      tj <- tvect[j]
      if (i == m & j == t) {
        foo <- bootstrap_sigma %>%
          filter((row_number() %% (m * k * k)) == 0)
      }
      else {
        foo <- bootstrap_sigma %>%
          filter((row_number() %% (m * k * k)) == (i - 1) * k * k + tj)
      }
      sigma_lower[i, j] <- 2 * sigma[(i - 1) * k * k + tj] -
        quantile(foo$sigma, 1 - (alpha / 2), names = FALSE)
      sigma_upper[i, j] <- 2 * sigma[(i - 1) * k * k + tj] -
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
