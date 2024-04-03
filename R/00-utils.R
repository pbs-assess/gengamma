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

get_gamma_cv <- function(phi) 1 / sqrt(phi)

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

choose_family <- function(fit_fam) {
  switch(fit_fam,
      "lognormal" = sdmTMB::lognormal(link = "log"),
      "gamma" = Gamma(link = "log"),
      "gengamma" = sdmTMB::gengamma(link = "log"),
      "delta-gamma" = sdmTMB::delta_gamma(link2 = "log", type = "standard"),
      "delta-lognormal" = sdmTMB::delta_lognormal(link2 = "log", type = "standard"),
      "delta-gengamma" = sdmTMB::delta_gengamma(link2 = "log", type = "standard"),
      "tweedie" = sdmTMB::tweedie(link = "log"),
      # Add more cases for other families if needed
      stop("Invalid family name")
  )
}

# Fit models across families
fit_cross <- function(.data, .mesh, sim_fam, .Q = NA, fit_fam,
  sp = "on", st = "off") {
  if (grepl("gengamma", sim_fam)) {
    dat <- .data |> filter(grepl("gengamma", family), Q == .Q)
  } else {
    dat <- .data |> filter(family == sim_fam, is.na(.Q))
  }

  .family <- choose_family(fit_fam)

  message("\tFitting data simulated from: ", sim_fam, " with - ", fit_fam)
  sdmTMB(formula = observed ~ 1,
    data = dat,
    mesh = .mesh,
    family = .family,
    spatial = sp,
    spatiotemporal = st
  )
}

get_sanity_df <- function(fit_obj, silent = FALSE) {
  tibble(
    sim_family = unique(fit_obj$data$family),
    fit_family = ifelse(family(fit_obj)[[1]][[1]] == 'tweedie', 'tweedie', family(fit_obj)[[2]][[1]]),
    Q = unique(fit_obj$data$Q),
    sanity_allok = sanity(fit_obj, silent = silent)$all_ok
    )
}

get_fitted_estimates <- function(fit_obj) {
  sd_rep_est <- as.list(fit_obj$sd_report, what = "Estimate")
  sd_rep_se <- as.list(fit_obj$sd_report, what = "Std")

  gg_Q <- sd_rep_est$gengamma_Q
  gg_Q_se <- sd_rep_se$gengamma_Q

  family1 <- family(fit_obj)[[1]][[1]]
  family2 <- ifelse(family1 == 'tweedie', 'tweedie', family(fit_obj)[[2]][[1]])
  mod <- ifelse(family1 == 'tweedie', 1, 2)

  tidy_ran <- tidy(fit_obj, effects = 'ran_pars', model = mod) |>
    tidyr::pivot_wider(names_from = term, values_from = c(estimate, std.error))

  out <- tibble(
    sim_family = unique(fit_obj$data$family),
    fit_family = family2,
    Q = unique(fit_obj$data$Q),
    est_Q = gg_Q,
    est_Qse = gg_Q_se,
  )
  bind_cols(out, tidy_ran)
}

get_index_summary <- function(predict_obj) {
  fit_obj <- predict_obj$fit_obj
  sd_rep_est <- as.list(predict_obj$fit_obj$sd_report, what = "Estimate")
  sd_rep_se <- as.list(predict_obj$fit_obj$sd_report, what = "Std")

  ran_pars <- tidy(predict_obj$fit_obj, 'ran_pars')

  gg_Q <- sd_rep_est$gengamma_Q
  gg_Q_se <- sd_rep_se$gengamma_Q

  index <- get_index(predict_obj, bias_correct = TRUE)
  mutate(index,
    sim_family = unique(fit_obj$data$family),
    sim_link = unique(fit_obj$data$link),
    #fit_family1 = family(fit_obj)[[1]][[1]],
    fit_family = ifelse(family(fit_obj)[[1]][[1]] == 'tweedie', 'tweedie', family(fit_obj)[[2]][[1]]),
    phi = unique(fit_obj$data$phi),
    cv = unique(fit_obj$data$cv),
    Q = unique(fit_obj$data$Q),
    sigma_O = unique(fit_obj$data$sigma_O),
    aic = AIC(fit_obj),
    est_Q = gg_Q,
    est_Qse = gg_Q_se,
    est_sigmaO = gg_Q,
    est_sigmaOse = gg_Q_se
  ) |>
  bind_rows()
}

beep <- function() {
  if (Sys.info()[["user"]] == "jilliandunic") beepr::beep()
}