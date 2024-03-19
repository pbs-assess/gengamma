library(dplyr)
library(tidyr)
library(purrr)
devtools::load_all('../sdmTMB')

source(here::here('R', '00-utils.R'))

dir.create(here::here('data-outputs'), showWarnings = FALSE, recursive = TRUE)
out_dir <- here::here('data-outputs')

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
#phi <- get_phi(family = 'lognormal', cv = cv)
b0 <- 0
sigma_O <- 0.2
tweedie_p <- 1.5

# Define the values of Q
#Q_values <- c(-2, -1, -0.5, -0.001, -0.0001, -0.00001, 0.00001, 0.0001, 0.001, 0.4, 0.5, 1, 2)
.gamma_q <- get_phi(family = 'gengamma-gamma-case', cv = cv, mu = 1)
Q_values <- c(-5, -2, -1, -0.5, -0.001, 0.001, 0.5, round(.gamma_q, digits = 3), 1, 2, 5)
gengamma_phi <- map_dbl(Q_values, ~ get_phi(family = 'gengamma', cv = cv, mu = 1, Q = .x))
# |--> SLOW; consider saving gengamma_phi values
#saveRDS(gengamma_phi, file.path(out_dir, 'gengamma-phi.rds'))
#gengamma_phi <- readRDS(file.path(out_dir, 'gengamma-phi.rds'))

# Simulate from binomial 
binom_sim <- sdmTMB_simulate(
  formula = ~ 1,
  data = predictor_dat,
  time = "year",
  mesh = mesh_sim,
  family = binomial(),
  range = 0.5, # Parameter that controls the decay of spatial correlation. If length 2, the spatial and spatiotemporal ranges will be unique.
  sigma_E = 0, # #sigma_E = 0.1,; SD of spatiotemporal process (Epsilon).
  seed = 42,
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
    seed = 42,
    B = c(b0) # B0 = intercept, B1 = a1 slope;  A vector of beta values (fixed-effect coefficient values).
    #B = c(0.2, -0.4) # B0 = intercept, B1 = a1 slope
  ) |>
  mutate(family = 'lognormal', link = 'log', cv = cv,
         sigma_O = sigma_O, b0 = b0, Q = NA) |>
  tibble::as_tibble() |>
  mutate(binom_mu = binom_sim$mu, 
         binom_obs = binom_sim$observed) |>
  rename(true = 'mu')
  #select(-observed, -family, -link) |>

# simulate a binomial example once 
# then multiply the probabilities from that by the mean from the positive part from each of the others
# i.e., the 'mu' from the binomial is the probability
# the CV is basically fixed for a binomial model with 1 trial
# all that matters is what the mean probability of observation is, 
# i.e. the factor year coefficient (which is in logit space)

# Get true index for later
true_index <- sim_dat |>
  group_by(year) |>
  summarise(biomass = sum(true * binom_mu))

sampled <- sample_n(sim_dat, size = 6000)
sampled_mesh <- make_mesh(sampled, c("X", "Y"), cutoff = 0.1)
#plot(sampled_mesh[[1]])

# Lognormal ----
# SIMULATE{y_i(i,m) = exp(rnorm(log(mu_i(i,m)) - pow(phi(m), Type(2)) / Type(2), phi(m)));}
ln_sim <- sampled |>
  mutate(family = 'lognormal', 
         phi = get_phi(family = 'lognormal', cv = cv),
         error = exp(rnorm(n(), log(true) - 0.5 * phi^2, phi))
  ) |>
  mutate(observed = exp(log(true) + log(error)))
# SIMULATE{y_i(i,m) = exp(rnorm(log(mu_i(i,m)) - pow(phi(m), Type(2)) / Type(2), phi(m)));}

# ----------------------------------
# trying to understand how to put together the delta model
# ----------------------------------
dln_sim <- ln_sim |>
  mutate(family = 'delta-lognormal', 
         observed = binom_obs * observed, 
         delta_true = binom_mu * true)

test <- sdmTMB(formula = observed ~ 1,
    data = dln_sim,
    mesh = sampled_mesh,
    family = sdmTMB::delta_lognormal(), 
    spatial = "on", 
    spatiotemporal = 'off'
  )

test_pred <- predict(test, newdata = predictor_dat, return_tmb_object = TRUE)

test_index <- get_index(test_pred, bias_correct = TRUE)
# ----------------------------------

# Gamma ----
ga_sim <- sampled |>
  mutate(family = 'gamma', 
         phi = get_phi(family = 'gamma', cv = cv),
         error = rgamma(n(), shape = phi, scale = true / phi)
  ) |>
  mutate(observed = exp(eta + log(error)))
# s1 = exp(ln_phi(m));        // shape
# s2 = mu_i(i,m) / s1;        // scale
# SIMULATE{y_i(i,m) = rgamma(s1, s2);}
# // s1 = Type(1) / (pow(phi, Type(2)));  // s1=shape, ln_phi=CV,shape=1/CV^2

# Tweedie ----
tw_sim <- sampled |>
  mutate(family = 'tweedie', 
         phi = get_phi(family = 'tweedie', cv = cv, p = tweedie_p, mu = exp(b0)), # QUESTION: correct fixing of mu???
         tweedie_p = tweedie_p,
         error = fishMod::rTweedie(n(), mu = true, phi = phi, p = tweedie_p)
  ) |>
  mutate(observed = exp(eta + log(error)))
# looks ok??
#sd(tw_sim$observed) / mean(tw_sim$observed)

# Gengamma ----

gg_sim <- purrr::map2_dfr(Q_values, gengamma_phi, \(.Q, .phi) {
  sampled |> 
    mutate(Q = .Q, 
           phi = .phi,
           family = 'gengamma', 
           error = rgengamma(n(), mean = true, sigma = .phi, Q = .Q))
    }
  ) |>
  mutate(observed = exp(eta + log(error)))
# tmp_ll = sdmTMB::dgengamma(y_i(i,m), mu_i(i,m), phi(m), gengamma_Q, true);
# SIMULATE{y_i(i,m) = sdmTMB::rgengamma(mu_i(i,m), phi(m), gengamma_Q);}



save(predictor_dat, true_index, sampled, sampled_mesh, ln_sim, ga_sim, gg_sim, Q_values, 
  file = file.path(out_dir, "sim-families.RData"))
# save(predictor_dat, true_index, sampled, sampled_mesh, ln_sim, ga_sim, tw_sim, gg_sim, Q_values, 
#   file = file.path(out_dir, "sim-families.RData"))

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
              error = rgengamma(n(), mean = true, sigma = .phi, Q = .Q))
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
  beepr::beep()
}
# ------------------------------------------------------------------------------
