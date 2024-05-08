library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(sdmTMB)

source(here::here('R', '00-utils.R'))

out_dir <- here::here('data-outputs', 'cross-sim')
error_dir <- file.path(out_dir, 'errors')
fit_dir <- file.path(out_dir, 'fits')
fit_sum_dir <- file.path(out_dir, 'fit-summaries')
ind_dir <- file.path(out_dir, 'index')
fig_dir <- here::here('figures', 'cross-sim')

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(error_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(ind_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fit_sum_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fit_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Steps:
# 1. Simulate data from Lognormal, Gamma, Tweedie, and Gengamma
#   - use starting CV of 0.8, calculate respective phi for each family
#   - double check that gengamma can recover gengamma_Q
# 2. Cross fit across families

# Simulate observed data from lognormal, gamma, gengamma
# ------------------------------------------------------------------------------
sim_fit <- function(predictor_dat, mesh_sim,
                   cv, b0, sigma_O, tweedie_p,
                   Q_values, gengamma_phi,
                   type = "", # provide option to toggle poisson-link on/off
                   sp = "on", st = "off",
                   sample_size = 500,
                   rep = NULL, get_simulation_output = FALSE,
                   save_fits = FALSE, fits_only = FALSE) {
  # QUESTION: Does the seed need to be the same for the binom_sim component and the positive component?
  # Simulate from binomial
  binom_sim <- sdmTMB_simulate(
    formula = ~ 1,
    data = predictor_dat,
    time = "year",
    mesh = mesh_sim,
    family = binomial(),
    range = 0.5, # Parameter that controls the decay of spatial correlation. If length 2, the spatial and spatiotemporal ranges will be unique.
    sigma_E = 0, # #sigma_E = 0.1,; SD of spatiotemporal process (Epsilon).
    #seed = 42,
    sigma_O = sigma_O, # #sigma_O = 0.2,; SD of spatial process (Omega)
    B = c(b0) # B0 = intercept, B1 = a1 slope;  A vector of beta values (fixed-effect coefficient values).
    #B = c(0.2, -0.4) # B0 = intercept, B1 = a1 slope
    ) |>
    as_tibble()

  # Generate dataframe of true mu's for 2nd lp
  sim_dat <- sdmTMB_simulate(
      formula = ~ 1,
      data = predictor_dat,
      time = "year",
      mesh = mesh_sim,
      family = Gamma(link = 'log'), # Delta families are not supported. Instead, simulate the two component models separately and combine.
      range = 0.5, # Parameter that controls the decay of spatial correlation. If length 2, the spatial and spatiotemporal ranges will be unique.
      phi = get_phi(family = 'lognormal', cv = cv), # Observation error scale parameter (e.g., SD in Gaussian).
      sigma_E = 0, # #sigma_E = 0.1,; SD of spatiotemporal process (Epsilon).
      sigma_O = sigma_O, # #sigma_O = 0.2,; SD of spatial process (Omega)
      #seed = 42,
      B = c(b0) # B0 = intercept, B1 = a1 slope;  A vector of beta values (fixed-effect coefficient values).
      #B = c(0.2, -0.4) # B0 = intercept, B1 = a1 slope
    ) |>
    mutate(family = 'gamma', link = 'log', cv = cv,
           sigma_O = sigma_O, b0 = b0, Q = NA) |>
    tibble::as_tibble() |>
    mutate(encounter_mu = binom_sim$mu,
           encounter_observed = binom_sim$observed)

  # Get true index for later
  true_index <- sim_dat |>
    group_by(year) |>
    summarise(biomass = sum(encounter_mu * mu))

  sampled <- sample_n(sim_dat, size = sample_size)
  # QUESTION: How many samples should be drawn?
  # I was looking at coverage from the real surveys which is ~5% of the survey grid / year (I think??)
  # But I think there are convergence issues with 4000 / 10000 as done here.
  sampled_mesh <- make_mesh(sampled, c("X", "Y"), cutoff = 0.1)
  #plot(sampled_mesh[[1]])

  # Delta-lognormal
  # ----------------------------------
  # SIMULATE{y_i(i,m) = exp(rnorm(log(mu_i(i,m)) - pow(phi(m), Type(2)) / Type(2), phi(m)));}
  dl_sim <- sampled |>
    mutate(family = paste0('delta-lognormal', type),
          phi = get_phi(family = 'lognormal', cv = cv),
          catch_observed = exp(rnorm(n(), log(mu) - 0.5 * phi^2, phi))
    ) |>
    mutate(observed = encounter_observed * catch_observed)

  # Delta-gamma
  # ----------------------------------
  # s1 = exp(ln_phi(m));        // shape
  # s2 = mu_i(i,m) / s1;        // scale
  # SIMULATE{y_i(i,m) = rgamma(s1, s2);}
  dg_sim <- sampled |>
    mutate(family = paste0('delta-gamma', type),
          phi = get_phi(family = 'gamma', cv = cv),
          catch_observed = rgamma(n(), shape = phi, scale = mu / phi)
    ) |>
    mutate(observed = encounter_observed * catch_observed)

  # Tweedie
  # ----------------------------------
  tw_sim <- sampled |>
    mutate(family = 'tweedie',
          phi = get_phi(family = 'tweedie', cv = cv, p = tweedie_p, mu = exp(b0)), # Fix of tweedie phi for given mu
          tweedie_p = tweedie_p,
          catch_observed = fishMod::rTweedie(n(), mu = mu, phi = phi, p = tweedie_p)
    ) |>
    mutate(observed = encounter_observed * catch_observed)

  # Gengamma
  # ----------------------------------
  # tmp_ll = sdmTMB::dgengamma(y_i(i,m), mu_i(i,m), phi(m), gengamma_Q, true);
  # SIMULATE{y_i(i,m) = sdmTMB::rgengamma(mu_i(i,m), phi(m), gengamma_Q);}
  dgg_sim <- purrr::map2_dfr(Q_values, gengamma_phi, \(.Q, .phi) {
    sampled |>
      mutate(Q = .Q,
            phi = .phi,
            family = paste0('delta-gengamma', type),
            catch_observed = rgengamma(n(), mean = mu, sigma = .phi, Q = .Q))
      }
    ) |>
    mutate(observed = encounter_observed * catch_observed)

  # ------------------------------------------------------------------------------
  # Fit cross-simulations
  # ------------------------------------------------------------------------------
  .families <- c('delta-lognormal', 'delta-gamma', 'tweedie', 'delta-gengamma')
  # if (type == '-poisson-link') {
  #   .families <- c('delta-lognormal-poisson-link', 'delta-gamma-poisson-link', 'tweedie', 'delta-gengamma-poisson-link')
  # }
  sim_df <- bind_rows(dl_sim, dg_sim, tw_sim, dgg_sim)

  if (get_simulation_output) {
    sim_list <- list(sim_df = sim_df,
      sampled = sampled,
      sampled_mesh = sampled_mesh,
      true_index = true_index)
    return(sim_list)
  }

  cross_combos <- distinct(sim_df, family, Q) |>
    tidyr::crossing(.fit_fam = .families) |>
    rename(.sim_fam = 'family', .Q = "Q")

  fits <- cross_combos |>
    pmap(\(.sim_fam, .Q, .fit_fam) {
      tryCatch({
        fit_cross(
          .data = sim_df, # data filtering built into `fit_cross()`
          .mesh = sampled_mesh,
          sim_fam = .sim_fam,
          fit_fam = .fit_fam,
          .Q = .Q,
          sp = sp,
          st = st)
        },
        error = function(e) {
          error_out <- sim_df |>
            mutate(fit_family = .fit_fam, rep = rep)
          error_filename <- paste0(.fit_fam, '-', rep, '-cv', cv, '-sigmao', sigma_O, '.rds')
          saveRDS(error_out, file.path(error_dir, error_filename))
        }
      )
    }
  )

  fits <- keep(fits, ~inherits(.x, "sdmTMB"))

  if (save_fits) {
    fits_filename <- paste0(rep, '-cv', cv, '-sigmao', sigma_O, '-b', b0, '-n', sample_size, '.rds')
    message("\tSaving: ", file.path(fit_dir, fits_filename))
    saveRDS(fits, file.path(fit_dir, fits_filename))
    if (fits_only) {
      return(fits)
    }
  }

  fit_summary <- map_dfr(fits, get_fitted_estimates) |>
    mutate(rep = rep)

  fit_sanity <- map_dfr(fits, \(x)
    tibble(
      sim_family = unique(x$data$family),
      #fit_family1 = family(x)[[1]][[1]],
      fit_family = ifelse(family(x)[[1]][[1]] == 'tweedie', 'tweedie', family(x)[[2]][[1]]),
      Q = unique(x$data$Q),
      sanity_allok = sanity(x)$all_ok)
    )

  # Only get predictions for cases where sanity checks passed
  sanity_pass <- fit_sanity$sanity_allok
  pred_list <- fits[sanity_pass] |>
    map(predict, newdata = predictor_dat, return_tmb_object = TRUE)

  index_df <- map_dfr(pred_list, get_index_summary) |>
    mutate(rep = rep)

  index_df <- index_df |>
    mutate(true = pull(true_index, 'biomass')) |>
    as_tibble()

  list(index_df = index_df, fit_summary = fit_summary)
}

n_year = 1
grid_width = 100
set.seed = 42
predictor_dat <- tidyr::expand_grid(
    X = seq(0, 1, length.out = grid_width), Y = seq(0, 1, length.out = grid_width),
    year = seq_len(n_year)
  )

mesh_sim <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.1)

# Observation error scale parameter (e.g., SD in Gaussian).
#cv_values <- c(0.2, 0.5, 0.8)
cv <- c(0.8, 0.95)[1] # lingcod wcvi cv ~0.83; dogfish wcvi gamma cv ~0.95
b0 <- 0
sigma_O <- c(0.2, 0.6, 1.0, 1.75, 3.2)[2] #dogfish wcvi gamma ~1.0
tweedie_p <- 1.5

# Define the values of Q
# .gamma_q <- get_phi(family = 'gengamma-gamma-case', cv = cv, mu = 1)
# Q_values <- c(-5, -2, -1, -0.5, -0.001, 0.001, 0.5, round(.gamma_q, digits = 3), 1, 2, 5)

if (!file.exists(file.path(out_dir, 'Q-values.txt'))) {
  Q_values <- c(-5, -2, -1, -0.5, -0.001, 0.001, 0.5, 0.8, 1, 2, 5) # Hard code .gamma_q for now
  dput(Q_values, file.path(out_dir, 'Q-values.txt'))
} else {
  Q_values <- dget(file.path(out_dir, 'Q-values.txt'))
}

if (!file.exists(file.path(out_dir, 'gengamma-phi.txt'))) {
  gengamma_phi <- map_dbl(Q_values, ~ get_phi(family = 'gengamma', cv = cv, mu = 1, Q = .x))
  dput(round(gengamma_phi, 5), file.path(out_dir, 'gengamma-phi.txt'))
} else {
  gengamma_phi <- dget(file.path(out_dir, 'gengamma-phi.txt'))
}
# Save one run of simulated data objects for benchmarking
# set.seed(42)
# sim_list <- sim_fit(rep = 1,
  # predictor_dat = predictor_dat, mesh_sim = mesh_sim,
  # cv = cv, b0 = b0, sigma_O = sigma_O, tweedie_p = tweedie_p,
  # Q_values = Q_values, gengamma_phi = gengamma_phi,
  # get_simulation_output = TRUE)
# saveRDS(sim_list, file.path(out_dir, 'sim-list.rds'))

set.seed(42)
fit <- sim_fit(rep = 42,
  predictor_dat = predictor_dat, mesh_sim = mesh_sim,
  cv = 0.95, b0 = 0, sigma_O = 1, tweedie_p = tweedie_p,
  Q_values = Q_values, gengamma_phi = gengamma_phi,
  sp = "off", sample_size = 1000,
  save_fits = TRUE, fits_only = TRUE)
beep()
# ------------------------------------------------------------------------------
# Cross-simulation
# ------------------------------------------------------------------------------
# Set up info for parallel runs of replicates
is_rstudio <- !is.na(Sys.getenv("RSTUDIO", unset = NA))
is_unix <- .Platform$OS.type == "unix"
cores <- parallel::detectCores()

# test <- sim_fit(predictor_dat = predictor_dat, mesh_sim = mesh_sim,
#     cv = cv, b0 = b0, sigma_O = sigma_O, tweedie_p = tweedie_p,
#     Q_values = Q_values, gengamma_phi = gengamma_phi)

# Non progress-tracking version:
# -------------------------------
# index_df <- furrr::future_map_dfr(1:n_reps, ~sim_fit(rep = .x))
# beep()
# future::plan(future::sequential)

# -------------------------------
# Use version with progress bar
# -------------------------------
cv <- c(0.8, 0.95)[1] # lingcod wcvi cv ~0.83; dogfish wcvi gamma cv ~0.95
b0 <- 0
sigma_O <- c(0.2, 0.6, 1.0, 1.75)[2] #dogfish wcvi gamma ~1.0
tweedie_p <- 1.5
type <- ""
#type <- "-poisson-link"
n_reps <- 200
tag <- paste0('cv', cv, '-sigmao', sigma_O, '-nreps', n_reps, type)

progressr::handlers(global = TRUE)
progressr::handlers("progress")
prog_fxn <- function(xs) {
  p <- progressr::progressor(along = xs)
  furrr::future_map(xs, function(x) {
    p(sprintf("x=%g", x))
    sim_fit(rep = x,
           predictor_dat = predictor_dat, mesh_sim = mesh_sim,
           cv = cv, b0 = b0, sigma_O = sigma_O, tweedie_p = tweedie_p,
           Q_values = Q_values, gengamma_phi = gengamma_phi)
  })
}

message("\tChecking for index file tag: \n\t   ", tag)
if (!file.exists(file.path(ind_dir, paste0(tag, '.rds')))) {
  options(future.globals.maxSize = 800 * 1024 ^ 2) # 800 mb
  # Multicore crashes for me after updating to R 4.4 :(
  # if (!is_rstudio && is_unix) {
  #   future::plan(future::multicore, workers = cores)
  # } else {
    future::plan(future::multisession, workers = cores)
  # }
  out <- prog_fxn(1:n_reps)
  future::plan(future::sequential)

  index_df <- map_dfr(out, 'index_df')
  fit_summary <- map_dfr(out, 'fit_summary')

  saveRDS(index_df, file.path(ind_dir, paste0(tag, '.rds')))
  saveRDS(fit_summary, file.path(fit_sum_dir, paste0(tag, '.rds')))
} else {
  index_df <- readRDS(file.path(ind_dir, paste0(tag, '.rds'))) |> as_tibble()
  fit_summary <- readRDS(file.path(fit_sum_dir, paste0(tag, '.rds'))) |> as_tibble()
}
beep()
# -------------------------------

# Self check on gg Q estimation
# ------------------------------------------------------------------------------
# Slow; also note that sanity check fails if Q is very negative, e.g., I tested with Q = -5
run_self_check <- FALSE
if (run_self_check) {
  local({
    .gamma_q <- get_phi(family = 'gengamma-gamma-case', cv = cv, mu = 1)
    Q_values <- c(-5, -2, -1, -0.5, -0.001, -0.0001, -0.00001, 0.00001, 0.0001, 0.001, 0.5, .gamma_q, 1, 2, 5)
    gengamma_phi <- map_dbl(Q_values, ~ get_phi(family = 'gengamma', cv = cv, mu = 1, Q = .x))

    gg_sim <- purrr::map2_dfr(Q_values, gengamma_phi, \(.Q, .phi) {
      sampled |>
        mutate(Q = .Q,
              phi = .phi,
              family = 'gengamma',
              error = rgengamma(n(), mean = mu, sigma = .phi, Q = .Q))
        }
      ) |>
      mutate(observed = exp(eta + log(error)))

    gg_fits <- tibble(sim_fam = rep('gengamma', length(Q_values)), .Q = Q_values) |>
      purrr::pmap(fit_cross, .data = gg_sim, .mesh = sampled_mesh, fit_fam = 'gengamma')

    Q_ests <- map_df(gg_fits, \(x) tibble(gengamma_Q = as.list(x$sd_report, 'Est')$gengamma_Q,
      Q_se = as.list(x$sd_report, 'Std')$gengamma_Q,
      sanity_allok = sanity(x)$all_ok)) |>
      mutate(.Q = Q_values)

    Q_ests
    }
  )
  beep()
}
# ------------------------------------------------------------------------------
