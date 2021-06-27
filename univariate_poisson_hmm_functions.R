library(tidyverse)
library(ggplot2)
library(ggfortify)
library(plyr)

graph_hmm_output <- function(output) {
  ggplot(output, aes(x = index, y = obs, color = state)) +
    geom_line() +
    theme_minimal() +
    scale_colour_continuous(type = "viridis")
}

graph_hmm_hist <- function(output) {
  ggplot(output, aes(obs)) +
    geom_histogram(binwidth = 1, colour = "navy", fill = "light blue") +
    theme_minimal()
}

# Transform Poisson natural parameters to working parameters
# Discussed in section 3.3, starting page 50
# m: number of states
# lambda: vector of parameters for Poisson dist
# gamma: matrix of transition probabilities
# delta: vector of initial distribution
pois_hmm_pn2pw <- function(m, lambda, gamma, delta = NULL, stationary = TRUE) {
  # Transform lambda by taking log
  tlambda <- log(lambda)
  # Transform gamma by log(gamma_ij/gamma_ii)
  foo <- log(gamma / diag(gamma))
  # Take only off-diagonal elements of resulting matrix
  tgamma <- as.vector(foo[!diag(m)])
  # Transform initial distribution parameter if the MC is not stationary
  if (stationary) {
    tdelta <- NULL
  }
  # Take log of each element except the first divided by the first element
  else {
    tdelta <- log(delta[-1] / delta[1])
  }
  parvect <- c(tlambda, tgamma, tdelta)
  return(parvect)
}

# Transform Poisson working parameters to natural parameters
pois_hmm_pw2pn <- function(m, parvect, stationary = TRUE) {
  # Transform by taking exponent of tlambda
  lambda <- exp(parvect[1:m])
  # Create mxm identity matrix
  gamma <- diag(m)
  # Set off-diagonal elements of matrix, starting from 1st col and moving down
  gamma[!gamma] <- exp(parvect[(m + 1):(m * m)])
  # Divide each element by the row sum, so the rows will sum to 1
  gamma <- gamma / apply(gamma, 1, sum)
  # Obtain stationary dist by solving system of equations
  if (stationary) {
    delta <- solve(t(diag(m) - gamma + 1), rep(1, m))
  }
  # Create vector of 1 and the elements of tdelta
  else {
    foo <- c(1, exp(parvect[(m * m + 1):(m * m + m - 1)]))
    # Divide each element of vector by the sum of the elements
    delta <- foo / sum(foo)
  }
  return(list(lambda = lambda, gamma = gamma, delta = delta))
}

# Computing minus the log-likelihood from the working parameters
# Algorithm is given on page 49
# parvect: working parameters
# x: observations
# m: number of states
pois_hmm_mllk <- function(parvect, x, m, stationary = TRUE, ...) {
  if (m == 1) {
    # If only one state, parvect should only contain log(lambda) (?)
    # Then log-lik is sum of log(p(x))
    return(-sum(dpois(x, exp(parvect), log = TRUE)))
  }
  n <- length(x)
  pn <- pois_hmm_pw2pn(m, parvect, stationary)
  # foo is delta*P(x1)
  foo <- pn$delta * dpois(x[1], pn$lambda)
  # sumfoo is w1 = delta*P(x1)*1'
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  # foo is phi1
  foo <- foo / sumfoo
  for (i in 2:n) {
    # If data is not missing, define p(xi)
    if (!is.na(x[i])) {
      p <- dpois(x[i], pn$lambda)
    }
    # Else p is identity matrix
    else {
      p <- rep(1, m)
    }
    # %*% is matrix multiplication
    # Here foo is v
    foo <- foo %*% pn$gamma * p
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    # Here foo is phi_t
    foo <- foo / sumfoo
  }
  mllk <- -lscale
  return(mllk)
}

# Computing MLE from natural parameters
# lambda0, gamma0, delta0 are starting values for estimation
pois_hmm_mle <- function(x, m, lambda0, gamma0,
                         delta0 = NULL, stationary = TRUE, ...) {
  # Get working parameters
  parvect0 <- pois_hmm_pn2pw(m, lambda0, gamma0, delta0, 
                             stationary = stationary)
  # nlm is an unconstrained minimizer
  mod <- nlm(pois_hmm_mllk, parvect0, x = x, m = m, 
             stationary = stationary)
  # Transform estimates of working parameters
  # into estimates of natural parameters
  pn <- pois_hmm_pw2pn(m, mod$estimate,
                       stationary = stationary)
  # Get the negative log likelihood value
  mllk <- mod$minimum
  # np is the number of estimated parameters
  np <- length(parvect0)
  aic <- 2 * (mllk + np)
  # n is the number of estimated parameters
  n <- sum(!is.na(x))
  bic <- 2 * mllk + np * log(n)
  list(
    m = m, lambda = pn$lambda, gamma = pn$gamma, delta = pn$delta,
    code = mod$code, mllk = mllk, aic = aic, bic = bic
  )
}

# Generating a sample
# ns: length of sample
# mod: list with elements m, lambda, gamma, delta
pois_hmm_generate_sample <- function(ns, mod) {
  # mvect is list of states to sample from
  mvect <- 1:mod$m
  state <- numeric(ns)
  # Get initial state
  state[1] <- sample(mvect, 1, prob = mod$delta)
  # Get remaining states
  for (i in 2:ns) state[i] <- sample(mvect, 1, prob = mod$gamma[state[i - 1], ])
  # Generate observations
  x <- rpois(ns, lambda = mod$lambda[state])
  return(data.frame(index = c(1:ns), state = state, obs = x))
}

# Global decoding by the Viterbi algorithm
# Discussed in Section 5.4.2 (page 88)
pois_hmm_viterbi <- function(x, mod) {
  n <- length(x)
  # xi is a nxm matrix of zeros
  xi <- matrix(0, n, mod$m)
  # foo is delta*P(x1)
  foo <- mod$delta * dpois(x[1], mod$lambda)
  # Set first row of xi
  xi[1, ] <- foo / sum(foo)
  for (t in 2:n) {
    # foo is xi_t, this is equation 5.9
    foo <- apply(xi[t - 1, ] * mod$gamma, 2, max) * dpois(x[t], mod$lambda)
    xi[t, ] <- foo / sum(foo)
  }
  # iv is the maximizing sequence of states
  iv <- numeric(n)
  iv[n] <- which.max(xi[n, ])
  for (t in (n - 1):1) {
    # This is equation 5.11
    iv[t] <- which.max(mod$gamma[, iv[t + 1]] * xi[t, ])
  }
  return(data.frame(index = c(1:n), state = iv))
}

# Get Poisson marginal distribution
pois_marginal <- function(n, mod) {
  # Get stationary distribution
  delta <- solve(t(diag(mod$m) - mod$gamma + 1), rep(1, mod$m))
  mpois <- delta[1] * dpois(x = 0:n, lambda = mod$lambda[1])
  for (i in 2:mod$m) {
    mpois <- mpois + delta[i] * dpois(x = 0:n, lambda = mod$lambda[i])
  }
  return(data.frame(x = 0:n, mpois = mpois))
}

# Computing log(forward probabilities) for Poisson distribution
# See Section 4.1.1
pois_hmm_lforward <- function(x, mod) {
  n <- length(x)
  # Empty matrix of size mxn
  lalpha <- matrix(NA, mod$m, n)
  foo <- mod$delta * dpois(x[1], mod$lambda)
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo / sumfoo
  # Set first column of the matrix lalpha
  lalpha[, 1] <- lscale + log(foo)
  for (i in 2:n) {
    foo <- foo %*% mod$gamma * dpois(x[i], mod$lambda)
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
    # Set ith column of the matrix lalpha to alpha_i
    lalpha[, i] <- log(foo) + lscale
  }
  return(lalpha)
}

# Computing log(backward probabilities) for Poisson distribution
# See Section 4.1.2
pois_hmm_lbackward <- function(x, mod) {
  n <- length(x)
  m <- mod$m
  lbeta <- matrix(NA, m, n)
  # Set last column of matrix lbeta to 0
  lbeta[, n] <- rep(0, m)
  foo <- rep(1 / m, m)
  lscale <- log(m)
  for (i in (n - 1):1) {
    foo <- mod$gamma %*% (dpois(x[i + 1], mod$lambda) * foo)
    # Set ith column of the matrix lbeta to beta_i
    lbeta[, i] <- log(foo) + lscale
    sumfoo <- sum(foo)
    foo <- foo / sumfoo
    lscale <- lscale + log(sumfoo)
  }
  return(lbeta)
}

# Conditional probabilities for Poisson distribution
# Conditional probability that observation at time t equals
# xc , given all observations other than that at time t.
# Note: xc is a vector and the result (dxc) is a matrix .
pois_hmm_conditional <- function(xc, x, mod) {
  n <- length(x)
  m <- mod$m
  nxc <- length(xc)
  dxc <- matrix(NA, nrow = nxc, ncol = n)
  px <- matrix(NA, nrow = m, ncol = nxc)
  # Set elements of px to p_i(xc_j)
  for (j in 1:nxc) px[, j] <- dpois(xc[j], mod$lambda)
  la <- pois_hmm_lforward(x, mod)
  lb <- pois_hmm_lbackward(x, mod)
  la <- cbind(log(mod$delta), la)
  # Get maximum value in each column
  lafact <- apply(la, 2, max)
  lbfact <- apply(lb, 2, max)
  for (i in 1:n) {
    foo <- (exp(la[, i] - lafact[i]) %*% mod$gamma) * exp(lb[, i] - lbfact[i])
    foo <- foo / sum(foo)
    dxc[, i] <- foo %*% px
  }
  return(dxc)
}

# Conditional probability that observation at time t equals
# xc, given observations at time 1 to t-1.
# Note: xc is a vector and the result (dxc) is a matrix .
pois_hmm_conditional2 <- function(xc, x, mod) {
  n <- length(x)
  m <- mod$m
  nxc <- length(xc)
  dxc <- matrix(NA, nrow = nxc, ncol = n)
  px <- matrix(NA, nrow = m, ncol = nxc)
  for (j in 1:nxc) px[, j] <- dpois(xc[j], mod$lambda)
  la <- pois_hmm_lforward(x, mod)
  la <- cbind(log(mod$delta), la)
  lafact <- apply(la, 2, max)
  for (i in 2:n) {
    foo <- exp(la[, i] - lafact[i])
    foo <- foo / sum(foo)
    dxc[, i] <- foo %*% mod$gamma %*% px
  }
  return(dxc)
}

# Normal pseudo-residuals for Poisson HMM
# Type can be "ordinary" or "forecast"
# Output is lower, upper and middle pseudo-residuals
pois_hmm_pseudo_residuals <- function(x, mod, type) {
  n <- length(x)
  if (type == "ordinary") {
    cdists <- pois_hmm_conditional(xc = 0:max(x), x, mod)
  }
  else if (type == "forecast") {
    cdists <- pois_hmm_conditional2(xc = 0:max(x), x, mod)
  }
  # Each column of cumdists starts with 0,
  # remaining elements i are cumulative sum of
  # elements in the column of cdists up to that cell
  cumdists <- rbind(rep(0, n), apply(cdists, 2, cumsum))
  ulo <- uhi <- rep(NA, n)
  for (i in 1:n) {
    ulo[i] <- cumdists[x[i] + 1, i]
    uhi[i] <- cumdists[x[i] + 2, i]
  }
  umi <- 0.5 * (ulo + uhi)
  npsr <- qnorm(rbind(ulo, umi, uhi))
  return(data.frame(lo = npsr[1, ], mi = npsr[2, ], hi = npsr[3, ]))
}
