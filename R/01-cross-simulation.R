library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(sdmTMB)

source(here::here('R', '00-utils.R'))

dir.create(here::here('data-outputs'), showWarnings = FALSE, recursive = TRUE)
dir.create(here::here('data-outputs', 'errors'), showWarnings = FALSE, recursive = TRUE)
dir.create(here::here('figures'), showWarnings = FALSE, recursive = TRUE)
out_dir <- here::here('data-outputs')
error_dir <- here::here('data-outputs', 'errors')
fig_dir <- here::here('figures')

# Steps:
# 1. Simulate data from Lognormal, Gamma, Tweedie, and Gengamma
#   - use starting CV of 0.8, calculate respective phi for each family
#   - double check that gengamma can recover gengamma_Q
# 2. Cross fit across families

# Simulate observed data from lognormal, gamma, gengamma
# ------------------------------------------------------------------------------
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
cv <- 0.8 # lingcod wcvi cv ~0.83
b0 <- 0
sigma_O <- 0.6
tweedie_p <- 1.5

# Define the values of Q
# .gamma_q <- get_phi(family = 'gengamma-gamma-case', cv = cv, mu = 1)
# Q_values <- c(-5, -2, -1, -0.5, -0.001, 0.001, 0.5, round(.gamma_q, digits = 3), 1, 2, 5)
Q_values <- c(-5, -2, -1, -0.5, -0.001, 0.001, 0.5, 0.8, 1, 2, 5) # Hard code .gamma_q for now
if (!file.exists(file.path(out_dir, 'gengamma-phi.txt'))) {
  gengamma_phi <- map_dbl(Q_values, ~ get_phi(family = 'gengamma', cv = cv, mu = 1, Q = .x))
  dput(round(gengamma_phi, 5), file.path(out_dir, 'gengamma-phi.txt'))
} else {
  gengamma_phi <- dget(file.path(out_dir, 'gengamma-phi.txt'))
}

sim_fit <- function(rep = NA, get_simulation_output = FALSE) {
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
      family = lognormal(link = 'log'), # Delta families are not supported. Instead, simulate the two component models separately and combine.
      range = 0.5, # Parameter that controls the decay of spatial correlation. If length 2, the spatial and spatiotemporal ranges will be unique.
      phi = get_phi(family = 'lognormal', cv = cv), # Observation error scale parameter (e.g., SD in Gaussian).
      sigma_E = 0, # #sigma_E = 0.1,; SD of spatiotemporal process (Epsilon).
      sigma_O = sigma_O, # #sigma_O = 0.2,; SD of spatial process (Omega)
      #seed = 42,
      B = c(b0) # B0 = intercept, B1 = a1 slope;  A vector of beta values (fixed-effect coefficient values).
      #B = c(0.2, -0.4) # B0 = intercept, B1 = a1 slope
    ) |>
    mutate(family = 'lognormal', link = 'log', cv = cv,
          sigma_O = sigma_O, b0 = b0, Q = NA) |>
    tibble::as_tibble() |>
    mutate(encounter_mu = binom_sim$mu,
          encounter_observed = binom_sim$observed)

  # Get true index for later
  true_index <- sim_dat |>
    group_by(year) |>
    summarise(biomass = sum(encounter_mu * mu))

  sampled <- sample_n(sim_dat, size = 500)
  # QUESTION: How many samples should be drawn?
  # I was looking at coverage from the real surveys which is ~5% of the survey grid / year (I think??)
  # But I think there are convergence issues with 4000 / 10000 as done here.
  sampled_mesh <- make_mesh(sampled, c("X", "Y"), cutoff = 0.1)
  #plot(sampled_mesh[[1]])

  # Delta-lognormal
  # ----------------------------------
  # SIMULATE{y_i(i,m) = exp(rnorm(log(mu_i(i,m)) - pow(phi(m), Type(2)) / Type(2), phi(m)));}
  dl_sim <- sampled |>
    mutate(family = 'delta-lognormal',
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
    mutate(family = 'delta-gamma',
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
            family = 'delta-gengamma',
            catch_observed = rgengamma(n(), mean = mu, sigma = .phi, Q = .Q))
      }
    ) |>
    mutate(observed = encounter_observed * catch_observed)

  # ------------------------------------------------------------------------------
  # Fit cross-simulations
  # ------------------------------------------------------------------------------
  .families <- c('delta-lognormal', 'delta-gamma', 'tweedie', 'delta-gengamma')
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
          .Q = .Q)
        }, 
        error = function(e) {
          error_out <- sim_df |>
            mutate(fit_family = .fit_fam, rep = rep)
          error_filename <- paste0(.fit_fam, '-', rep, '.rds')
          saveRDS(error_out, file.path(error_dir, error_filename))
        }
      )
    }
  )
  fits <- keep(fits, ~inherits(.x, "sdmTMB"))
  # saveRDS(fits, file.path(fit_dir, 'fits.rds'))
  # fits <- readRDS(file.path(fit_dir, 'fits.rds'))

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

  index_df <- map_dfr(pred_list, get_index_summary)
  
  if (!is.na(rep)) {
    index_df$rep = rep
  }

  index_df |>
    mutate(true = pull(true_index, 'biomass'))
}

# Save one run of simulated data objects for benchmarking
set.seed(42)
sim_list <- sim_fit(1, get_simulation_output = TRUE)
saveRDS(sim_list, file.path(out_dir, 'sim-list.rds'))

# ------------------------------------------------------------------------------
# Cross-simulation
# ------------------------------------------------------------------------------
# Set up info for parallel runs of replicates
is_rstudio <- !is.na(Sys.getenv("RSTUDIO", unset = NA))
is_unix <- .Platform$OS.type == "unix"
cores <- parallel::detectCores() - 2
options(future.globals.maxSize = 800 * 1024 ^ 2) # 800 mb
if (!is_rstudio && is_unix) {
  future::plan(future::multicore, workers = cores)
} else {
  future::plan(future::multisession, workers = cores)
}

n_reps <- 100

# Non progress-tracking version:
# -------------------------------
# index_df <- furrr::future_map_dfr(1:n_reps, ~sim_fit(rep = .x))
# beep()
# future::plan(future::sequential)

# -------------------------------
# Use version with progress bar
# -------------------------------
progressr::handlers(global = TRUE)
progressr::handlers("progress")
prog_fxn <- function(xs) {
  p <- progressr::progressor(along = xs)
  furrr::future_map_dfr(xs, function(x) {
    p(sprintf("x=%g", x))
    sim_fit(rep = x)
  })
}
future::plan(future::multicore, workers = cores)
index_df <- prog_fxn(1:n_reps)
beep()
future::plan(future::sequential)
# -------------------------------
#saveRDS(index_df, file.path(out_dir, 'index_df.rds'))
#index_df <- readRDS(file.path(out_dir, 'index_df.rds')) |> as_tibble()

cross_combos <- bind_rows(
  tibble(Q = NA, sim_family = c('delta-lognormal', 'delta-gamma', 'tweedie')),
  tidyr::crossing(Q = Q_values, sim_family = 'delta-gengamma')
  ) |>
  tidyr::crossing(fit_family = c('lognormal', 'Gamma', 'tweedie', 'gengamma'))

sanity_tally <- left_join(
    cross_combos |> tidyr::unite(col = 'sim_combo_key', sim_family, Q, fit_family, sep = ":", remove = TRUE), 
    index_df |> tidyr::unite(col = 'sim_combo_key', sim_family, Q, fit_family, sep = ":", remove = FALSE), 
    by = c('sim_combo_key'))

plot_df <- index_df |>
  group_by(rep, sim_family, fit_family, Q, `_sdmTMB_time`) |>
  mutate(RMSE = sqrt(mean((log(est) - log(true))^2)),
         MRE = mean((est - true) / true),
         covered = lwr < true & upr > true,
         ci_width = upr - lwr,
         title = ifelse(is.na(Q), sim_family, paste0(sim_family, ": Q=", signif(Q, digits = 2)))
  ) |>
  ungroup(fit_family) |>
  mutate(min_aic = min(aic),
         d_aic = aic - min_aic) |>
  ungroup()
title_levels <- unique(plot_df$title)
title_levels <- c(title_levels[-1], title_levels[1])
plot_df <- plot_df |>
  mutate(title = factor(title, levels = title_levels))


# - Compare: RMSE, MRE, coverage, AIC
# - look at the consequence of selecting the wrong family (what does AIC choose)
# - examine bias in estimates and how the above relates to the Q value
theme_set(theme_light(base_size = 12))

plot_violin <- function(.data, .x, .ncol = NULL) {
  ggplot(data = .data, aes(x = {{.x}}, y = fit_family)) +
    stat_summary(fun = mean, geom = "point") +
    geom_violin(aes(col = fit_family), alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    scale_color_brewer(palette = "Dark2") +
    facet_wrap(~ title, ncol = .ncol) + 
    guides(colour = 'none') + 
}

plot_linedot <- function(.data, .x, .ncol = NULL) {
  ggplot(data = .data, aes(x = {{.x}}, y = fit_family, colour = fit_family)) +
  geom_linerange(xmin = 0, mapping = aes(xmax = {{.x}})) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0.95, linetype = 'dashed') +
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(~ title, ncol = .ncol) + 
  guides(colour = 'none')
}

plot_violin(plot_df, .x = RMSE, .ncol = 5) + 
  ggtitle('RMSE')
ggsave(filename = file.path(fig_dir, 'rmse.png'), width = 11, height = 6.5)


plot_violin(plot_df, .x = MRE, .ncol = 5) +
  ggtitle('MRE')
ggsave(filename = file.path(fig_dir, 'mre.png'), width = 11, height = 6.5)

plot_violin(plot_df, .x = (d_aic + 1), .ncol = 5) + 
  scale_x_continuous(trans = 'log10') + 
  geom_vline(xintercept = 1, linetype = 'dashed') +
  ggtitle("Delta AIC")
ggsave(filename = file.path(fig_dir, 'daic.png'), width = 11, height = 6.5)

plot_df |>
  group_by(title, fit_family, Q, sim_family) |>
  summarise(n_sanity_pass = n(), 
            prop_covered = sum(covered) / n_sanity_pass
         ) |>
plot_linedot(.x = prop_covered, .ncol = 5) + 
  ggtitle("95% CI Coverage")
ggsave(filename = file.path(fig_dir, 'ci-coverage.png'), width = 11, height = 6.5)

plot_violin(plot_df, .x = ci_width, .ncol = 5) +
  scale_x_continuous(trans = 'log10') +
  ggtitle("95% CI Width")
ggsave(filename = file.path(fig_dir, 'ci-width.png'), width = 11, height = 6.5)

sanity_tally |>
  mutate(sanity_pass = ifelse(is.na(est), 0, 1)) |>
  group_by(sim_family, Q, fit_family) |>
  summarise(n_pass = sum(sanity_pass), pass_prop = n_pass / n_reps) |>
  mutate(title = ifelse(is.na(Q), sim_family, paste0(sim_family, ": Q=", signif(Q, digits = 2)))) |>
  mutate(title = factor(title, levels = title_levels)) |>
plot_linedot(.x = pass_prop, .ncol = 5) +
  ggtitle("Proportion of fits that passed sanity check")
ggsave(filename = file.path(fig_dir, 'sanity-pass.png'), width = 11, height = 6.5)

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
