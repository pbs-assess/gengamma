rgengamma <- function(n, mean, sigma, Q) {
  if (Q == 0) {
    message('Q == 0, using rlnorm')
    y <- rlnorm(n, mean, sigma)
  } else {
    # Get mu from mean
    k <- Q^-2
    beta <- Q / sigma
    log_theta <- log(mean) - lgamma( (k*beta+1)/beta ) + lgamma( k )
    mu <- log_theta + log(k) / beta
    w <- log(Q^2 * rgamma(n, 1 / Q^(2), 1)) / Q
    y <- exp(mu + (sigma * w))
  }
  return(y)
}

get_dgmu <- function(Q, sigma, mean) {
  assertthat::assert_that(Q != 0, msg = "Q == 0, cannot compute mu")
  k <- Q^-2
  beta <- Q / sigma
  log_theta <- log(mean) - lgamma( (k*beta+1)/beta ) + lgamma( k )
  mu <- log_theta + log(k) / beta
  mu
}

# This is modified from the sdmtmb cpp file
sdm_rgengamma <- function(n, mean, sigma, Q) {
 lambda = Q
 assertthat::assert_that(lambda != 0, msg = "Q == 0, use rlnorm")
 k = lambda^(-2)
 beta = sigma^(-1) * lambda
 log_theta = log(mean) - lgamma((k * beta + 1) / beta) + lgamma(k)
 w = log(rgamma(n, k, 1))
 y = w / beta + log_theta
  exp(y)
}

# check rng consistency with what is expected using what is in flexsurv documentation
# set.seed(1)
# sdm_rgengamma(n = 1, mean = 0.2, sigma = 0.5, Q = 0)

sigma <- 0.5
Q <- 0.5
set.seed(1)
rgengamma(n = 1, mean = 0.2, sigma = sigma, Q = Q)
set.seed(1)
rgamma(n = 1, shape = 1 / 0.5^2, rate = exp(-get_dgmu(mean = 0.2, sigma = sigma, Q = Q)) / 0.5^2)
# Test the case using flexsurv::rgengamma()
local({ set.seed(1)
  mu <- get_dgmu(mean = 0.2, sigma = sigma, Q = Q)
  w <- log(Q^2 * rgamma(n = 1, 1/Q^2, 1))/Q
  exp(mu + sigma * w)
})

