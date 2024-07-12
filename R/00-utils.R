rgengamma <- function(n, mean, sigma, Q) {
  lambda <- Q
  ind0 <- lambda == 0
  y <- numeric(n)

  if (any(ind0)) {
    message("Q == 0, using rlnorm")
    y[ind0] <- rlnorm(n, mean, sigma)
  }

  if (any(!ind0)) {
    # Get mu from mean
    k <- lambda^-2
    beta <- lambda / sigma
    log_theta <- log(mean) - lgamma((k * beta + 1) / beta) + lgamma(k)
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
    cv <- sd(x) / mean(x)
    (cv - desired_cv)^2
  }

  f2 <- function(sigma, desired_cv = cv) {
    print(sigma)
    set.seed(1)
    x <- rgengamma(3e6, mean = mu, sigma = sigma, Q = sigma)
    cv <- sd(x) / mean(x)
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
    "delta-gamma-poisson-link" = sdmTMB::delta_gamma(link2 = "log", type = "poisson-link"),
    "delta-lognormal-poisson-link" = sdmTMB::delta_lognormal(link2 = "log", type = "poisson-link"),
    "delta-gengamma-poisson-link" = sdmTMB::delta_gengamma(link2 = "log", type = "poisson-link"),
    "tweedie" = sdmTMB::tweedie(link = "log"),
    # Add more cases for other families if needed
    stop("Invalid family name")
  )
}

# Fit real data
get_fit <- function(survey_dat, formula, region, family, species = NULL,
                   cutoff = 8, time = "year",
                   sp = "on", st = "iid", offset = "offset",
                   use_priors = FALSE, ...) {
  survey_dat <- filter(survey_dat, survey_abbrev %in% region)

  if (is.null(species)) {
    species <- unique(survey_dat$species)
  } else {
    survey_dat <- filter(survey_dat, species == {{species}})
  }
  mesh <- sdmTMB::make_mesh(survey_dat, xy_cols = c("X", "Y"), cutoff = cutoff)

  ny <- length(unique(survey_dat$year))

  b_priors <- sdmTMB::normal(NA, NA)

  if (use_priors) {
    b_priors <- sdmTMB::normal(rep(0, ny), rep(30, ny))
  }

  fit <- tryCatch(
    sdmTMB::sdmTMB(
      formula = formula,
      data = survey_dat,
      mesh = mesh,
      time = "year",
      spatial = sp,
      spatiotemporal = st,
      offset = "offset",
      family = choose_family(family),
      priors = sdmTMB::sdmTMBpriors(b = b_priors),
      ...
    ),
    error = function(e) paste(species, region, family, "\n\tError:", e, sep = " - ")
  )

  if (inherits(fit, 'sdmTMB')) {
    sanity_check <- all(unlist(sdmTMB::sanity(fit, gradient_thresh = 0.005)))

    if (!sanity_check) {
      sigma_check <- check_sigma_collapse(fit)
      if ("sigma_O" %in% sigma_check) sp = "off"
      if ("sigma_E" %in% sigma_check) st = "off"

      # Deal with non-positive definite hessian - how to decide what random field to turn off
      # for now let's just turn off the spatiotemporal to start and see how it goes
      # because we aren't including any covariates I think it might make more
      # sense to turn off the st field?
      if (isFALSE(fit$pos_def_hessian)) {
        st <- "off"
      }
    }
  } else {
    sanity_check <- FALSE
  }

  # Turn off spatial and/or spatiotemporal rf that collapse to zero
  if (!sanity_check) {
    #message("\tFitting: sp = ", sp, ", st = ", st, " for ", species, "-", region, "-", family, " failed")
    message("\tUpdating with sp = ", sp, "st = ", st)
    fit <- tryCatch(
      update(fit, spatial = sp, spatiotemporal = st),
      error = function(e) paste(species, region, family, "\n\tError:", e, sep = " - ")
    )
  }
  fit
}

is_delta <- function(fit_obj) ifelse(family(fit_obj)[[1]][[1]] == "tweedie", FALSE, TRUE)

check_sigma_collapse <- function(fit) {
  rp <- sdmTMB::tidy(fit, "ran_pars", 1) |>
    mutate(m = 1)

  if (is_delta(fit)) {
    rp2 <- sdmTMB::tidy(fit, "ran_pars", 2) |>
      mutate(m = 2)
    rp <- bind_rows(rp, rp2)
  }

  # Find indices where term contains "sigma" and estimate < 0.01
  id_s <- grep("sigma", rp$term)
  id_cs <- id_s[rp$estimate[id_s] < 0.01]
  collapsed_sigma <- rp[id_cs, ] # Does it mean anything if it is the sigma in the encounter vs catch that collapses?
  unique(collapsed_sigma$term)
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
  sdmTMB::sdmTMB(formula = observed ~ 1,
    data = dat,
    mesh = .mesh,
    family = .family,
    spatial = sp,
    spatiotemporal = st
  )
}

get_sanity_df <- function(fit_obj, real_data = FALSE, silent = FALSE, .gradient_thresh = 0.001) {
  type <- family(fit_obj)$type
  if (!real_data) {
    out <- tibble(
      sim_family = unique(fit_obj$data$family),
      fit_family = ifelse(family(fit_obj)[[1]][[1]] == 'tweedie', 'tweedie', family(fit_obj)[[2]][[1]]),
      Q = unique(fit_obj$data$Q),
      sanity_allok = sanity(fit_obj, silent = silent, gradient_thresh = .gradient_thresh)$all_ok,
      gradient_thresh = .gradient_thresh,
      type = type
    )
  } else {
    out <- tibble(
      species = unique(fit_obj$data$species),
      region = unique(fit_obj$data$survey_abbrev),
      fit_family = ifelse(family(fit_obj)[[1]][[1]] == 'tweedie', 'tweedie', family(fit_obj)[[2]][[1]]),
      gradient_thresh = .gradient_thresh,
      type = type,
      spatial = unique(fit_obj$spatial),
      spatiotemporal = unique(fit_obj$spatiotemporal)
    ) |>
    bind_cols(tibble::tibble(!!!sanity(fit_obj, silent = silent, gradient_thresh = .gradient_thresh)))
  }
}

get_fitted_estimates <- function(fit_obj, real_data = FALSE) {
  if (inherits(fit_obj, "sdmTMB")) {
    type <- family(fit_obj)$type

    sd_rep_est <- as.list(fit_obj$sd_report, what = "Estimate")
    sd_rep_se <- as.list(fit_obj$sd_report, what = "Std")

    gg_Q <- sd_rep_est$gengamma_Q
    gg_Q_se <- sd_rep_se$gengamma_Q

    family1 <- family(fit_obj)[[1]][[1]]
    family2 <- ifelse(family1 == "tweedie", "tweedie", family(fit_obj)[[2]][[1]])
    mod <- ifelse(family1 == "tweedie", 1, 2)

    tidy_ran <- tidy(fit_obj, effects = "ran_pars", model = mod) |>
      tidyr::pivot_wider(names_from = term, values_from = c(estimate, std.error, conf.low, conf.high))

    if (!real_data) {
      out <- tibble(
        sim_family = unique(fit_obj$data$family),
        fit_family = family2,
        Q = unique(fit_obj$data$Q),
        est_Q = gg_Q,
        est_Qse = gg_Q_se,
        spatial = unique(fit_obj$spatial),
        spatiotemporal = unique(fit_obj$spatiotemporal),
        type = type
      )
    } else {
      region <- unique(fit_obj$data$survey_abbrev)
      species <- unique(fit_obj$data$species)
      out <- tibble(
        species = species,
        region = region,
        fit_family = family2,
        type = type,
        est_Q = gg_Q,
        est_Qse = gg_Q_se,
        spatial = unique(fit_obj$spatial),
        spatiotemporal = unique(fit_obj$spatiotemporal),
        aic = AIC(fit_obj)
      )
    }
    bind_cols(out, tidy_ran)
  } else {
    fit_obj
  }
}

get_index_summary <- function(predict_obj) {
  fit_obj <- predict_obj$fit_obj

  type <- family(fit_obj)$type
  sd_rep_est <- as.list(predict_obj$fit_obj$sd_report, what = "Estimate")
  sd_rep_se <- as.list(predict_obj$fit_obj$sd_report, what = "Std")

  # ran_pars <- tidy(predict_obj$fit_obj, 'ran_pars')
  # est_sigmaO <- ran_pars$estimate[ran_pars$term == 'sigma_O']
  # est_sigmaOse <- ran_pars$std.error[ran_pars$term == 'sigma_O']

  gg_Q <- sd_rep_est$gengamma_Q
  gg_Q_se <- sd_rep_se$gengamma_Q

  index <- get_index(predict_obj, bias_correct = TRUE)
  mutate(index,
    sim_family = unique(fit_obj$data$family),
    sim_link = unique(fit_obj$data$link),
    #fit_family1 = family(fit_obj)[[1]][[1]],
    fit_family = ifelse(family(fit_obj)[[1]][[1]] == 'tweedie', 'tweedie', family(fit_obj)[[2]][[1]]),
    type = type,
    phi = unique(fit_obj$data$phi),
    cv = unique(fit_obj$data$cv),
    Q = unique(fit_obj$data$Q),
    sigma_O = unique(fit_obj$data$sigma_O),
    aic = AIC(fit_obj),
    est_Q = gg_Q,
    est_Qse = gg_Q_se#,
    # est_sigmaO = est_sigmaO,
    # est_sigmaOse = est_sigmaOse
  ) |>
    bind_rows()
}

clean_name <- function(x) gsub("/", "-", gsub(" ", "-", x))

beep <- function() {
  if (Sys.info()[["user"]] == "jilliandunic") beepr::beep()
}

choose_multi_type <- function(cores = cores) {
  is_rstudio <- !is.na(Sys.getenv("RSTUDIO", unset = NA))
  is_unix <- .Platform$OS.type == "unix"
  if (!is_rstudio && is_unix) {
    future::plan(future::multicore, workers = cores)
  } else {
    future::plan(future::multisession, workers = cores)
  }
}

# ---------------
# Residuals
# ---------------
get_rqr <- function(fit_obj, id, m = NULL) {
  if (is.null(m)) {
    m <- if (family(fit_obj)[[1]][1] == "tweedie") 1 else 2
  } else {
    m <- m
  }
  tibble(r = residuals(fit_obj, type = "mle-mvn", model = m), id = id)
}

get_dr <- function(fit_obj, fit_id, nsim = 200, seed = sample.int(1e6, 1), type = "mle-mvn", ...) {
  set.seed(seed)
  simulate(fit_obj, nsim = nsim, type = type, ...) |>
    dharma_residuals(fit_obj, plot = FALSE, ...) |>
    mutate(id = fit_id, seed = seed)
}

plot_resid <- function(resids) {
  ggplot(resids, aes(x = expected, y = observed)) +
    geom_point(shape = 21) +
    geom_abline(intercept = 0, slope = 1) +
    facet_grid(fit_family ~ title) +
    guides(colour = "none")
}

# https://stackoverflow.com/questions/13690184/update-inside-a-function-only-searches-the-global-environment
my_update <- function(mod, formula = NULL, data = NULL) {
  call <- getCall(mod)
  if (is.null(call)) {
    stop("Model object does not support updating (no call)", call. = FALSE)
  }
  term <- terms(mod)
  if (is.null(term)) {
    stop("Model object does not support updating (no terms)", call. = FALSE)
  }

  if (!is.null(data)) call$data <- data
  if (!is.null(formula)) call$formula <- update.formula(call$formula, formula)
  env <- attr(term, ".Environment")

  eval(call, env, parent.frame())
}
