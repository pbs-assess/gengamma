rgengamma <- function(n, mean, sigma, Q) {
  lambda <- Q
  ind0 <- lambda == 0
  y <- numeric(n)

  if (any(ind0)) {
    message('Q == 0, using rlnorm')
    y[ind0] <- rlnorm(n, mean, sigma)
  } 
  
  if (any(!ind0)) {
    # Get mu from mean
    k <- lambda^-2
    beta <- lambda / sigma
    log_theta <- log(mean) - lgamma( (k*beta+1)/beta ) + lgamma( k )
    mu <- log_theta + log(k) / beta
    w <- log(lambda^2 * rgamma(n, 1 / lambda^(2), 1)) / lambda
    y[!ind0] <- exp(mu + (sigma * w))
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

get_phi <- function(cv, family, mu, p, Q) {
  # Numeric solution to get sigma for desired_cv
  f <- function(sigma, q = Q, desired_cv = cv) {
    print(sigma)
    set.seed(1)
    x <- rgengamma(3e6, mean = mu, sigma = sigma, Q = q)
    cv <- sd(x)/mean(x)
    (cv - desired_cv)^2
  }

  f2 <- function(sigma, desired_cv = cv) {
    print(sigma)
    set.seed(1)
    x <- rgengamma(3e6, mean = mu, sigma = sigma, Q = sigma)
    cv <- sd(x)/mean(x)
    (cv - desired_cv)^2
  }
  
  switch(family, 
    "gamma" = 1 / cv^2, 
    "lognormal" = log(1 + cv^2)^0.5, # simplification of `sdconv(): cv is independent of mu
    "tweedie" = ((mu * cv)^2) / (mu^p), 
    "gengamma" = optimize(f, c(0.0001, 3), tol = 0.00001, desired_cv = cv, q = Q)$minimum,
    "gengamma-gamma-case" = optimize(f2, c(0.0001, 3), tol = 0.00001, desired_cv = cv)$minimum
    )
}

# Fit models across families
fit_cross <- function(.data, .mesh, sim_fam, .Q = NA, fit_fam) {
  if (sim_fam == 'delta-gengamma') {
    dat <- .data |> filter(family == 'delta-gengamma', Q == .Q)
  } else {
    dat <- .data |> filter(family == sim_fam, is.na(.Q))
  }

  .family <- switch(fit_fam,
      "delta-gamma" = sdmTMB::delta_gamma(),
      "delta-lognormal" = sdmTMB::delta_lognormal(),
      "delta-gengamma" = sdmTMB::delta_gengamma(),
      "tweedie" = sdmTMB::tweedie(link = "log"),
      # Add more cases for other families if needed
      stop("Invalid family name")
  )

  message("\tFitting data simulated from: ", sim_fam, " with - ", fit_fam)
  sdmTMB(formula = observed ~ 1,
    data = dat,
    mesh = .mesh,
    family = .family, 
    spatial = "on", 
    spatiotemporal = 'off'
  )
}

get_index_summary <- function(predict_obj) {
  index <- get_index(predict_obj, bias_correct = TRUE)
  mutate(index, 
    sim_family = unique(predict_obj$fit_obj$data$family),
    sim_link = unique(predict_obj$fit_obj$data$link),
    fit_family1 = family(predict_obj$fit_obj)[[1]][[1]],
    fit_family2 = ifelse(fit_family1 == 'tweedie', 'tweedie', family(predict_obj$fit_obj)[[2]][[1]]),
    phi = unique(predict_obj$fit_obj$data$phi),
    Q = unique(predict_obj$fit_obj$data$Q),
    sigma_O = unique(predict_obj$fit_obj$data$sigma_O)
  ) |>
  bind_rows()
}
