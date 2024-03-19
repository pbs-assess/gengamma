library(dplyr)
library(ggplot2)
library(purrr)
devtools::load_all('../sdmTMB')
theme_set(theme_light())

source(here::here('R', '00-utils.R'))
load(file.path(here::here("data-outputs", "sim-families.RData")))
fit_dir <- file.path(here::here("data-outputs", "fits"))

# Deal with tweedie separately. Need to simulate model with binomial component?
#.families <- c('lognormal', 'gamma', 'tweedie', 'gengamma')
.families <- c('lognormal', 'gamma', 'gengamma')

# Compare duration of sdmTMB calls with different fit families
# ------------------------------------------------------------------------------
benchmark_expr <- function(sim_dat, fit_family) {
  sdmTMB(formula = observed ~ 1,
      data = sim_dat,
      mesh = sampled_mesh,
      family = fit_family, 
      spatial = "on", 
      spatiotemporal = 'off'
    )
}

test <- rbenchmark::benchmark(
  "lognormal fit to log" = benchmark_expr(sim_dat = ln_sim, fit_family = sdmTMB::lognormal(link = 'log')),
  "gamma fit to log" = benchmark_expr(sim_dat = ln_sim, fit_family = Gamma(link = 'log')),
  #benchmark_expr(sim_dat = ln_sim, fit_family = sdmTMB::tweedie(link = 'log')),
  "gengamma fit to log" = benchmark_expr(sim_dat = ln_sim, fit_family = sdmTMB::gengamma(link = 'log')),  
  "lognormal fit to ga sim" = benchmark_expr(sim_dat = ga_sim, fit_family = sdmTMB::lognormal(link = 'log')),
  "gamma fit to ga sim" = benchmark_expr(sim_dat = ga_sim, fit_family = Gamma(link = 'log')),
  #benchmark_expr(sim_dat = ga_sim, fit_family = sdmTMB::tweedie(link = 'log')),
  "gengamma fit to ga sim" = benchmark_expr(sim_dat = ga_sim, fit_family = sdmTMB::gengamma(link = 'log')),  
  "lognormal fit to gg sim" = benchmark_expr(sim_dat = gg_sim, fit_family = sdmTMB::lognormal(link = 'log')),
  "gamma fit to gg sim" = benchmark_expr(sim_dat = gg_sim, fit_family = Gamma(link = 'log')),
  #benchmark_expr(sim_dat = gg_sim, fit_family = sdmTMB::tweedie(link = 'log')),
  "gengamma fit to gg sim" = benchmark_expr(sim_dat = gg_sim, fit_family = sdmTMB::gengamma(link = 'log')),  
  replications = 10, 
  columns = c("test", "replications", "elapsed",
                      "relative", "user.self", "sys.self")
)
beepr::beep()

# ------------------------------------------------------------------------------
# Fit cross-simulation for data from all but gengamma
ln_fits <- .families |>
  purrr::map(\(x) fit_cross(fit_fam = x, .data = ln_sim, .mesh = sampled_mesh, sim_fam = 'lognormal'))

ga_fits <- .families |>
  purrr::map(\(x) fit_cross(fit_fam = x, .data = ga_sim, .mesh = sampled_mesh, sim_fam = 'gamma'))

# tw_fits <- .families |>
#   purrr::map(\(x) fit_cross(fit_fam = x, .data = tw_sim, .mesh = sampled_mesh, sim_fam = 'tweedie'))

# Fit observed from gengamma distribution across different Q values ------------
Q_fam_combos <- tidyr::crossing(.Q = Q_values, .fit_fam = .families)

gg_cross_fits <- Q_fam_combos |>
  pmap(\(.Q, .fit_fam) fit_cross(.data = gg_sim, .mesh = sampled_mesh, sim_fam = 'gengamma', 
    .Q = .Q, fit_fam = .fit_fam))

Q_sim_fam_combos <- tidyr::crossing(.Q = Q_values, 
  .fit_fam = .families, 
  .sim_fam = .families
)

fit_list <- c(ln_fits, ga_fits, gg_cross_fits)
#saveRDS(fit_list, file.path(fit_dir, 'fits.R'))
#fit_list <- readRDS(file.path(fit_dir, 'fits.R'))

fit_sanity <- map_dfr(fit_list, \(x) tibble(
  sim_family = unique(x$data$family),
  fit_family = x$family[[1]],
  Q = unique(x$data$Q),
  sanity_allok = sanity(x)$all_ok)
  )

sanity_pass <- fit_sanity$sanity_allok

pred_list <- fit_list[sanity_pass] |>
  map(predict, newdata = predictor_dat, return_tmb_object = TRUE)

index_df <- map_dfr(pred_list, get_index_summary)
#beepr::beep()
#saveRDS(index_df, file.path(here::here("data-outputs"), 'cross-fit-index-df.rds')
#index_df <- readRDS(file.path(here::here("data-outputs"), 'cross-fit-index-df.rds'))

AIC_df <- map_dfr(fit_list, \(fit) tibble(
  aic = AIC(fit),
  sim_family = unique(fit$data$family),
  fit_family = fit$family[[1]], 
  Q = unique(fit$data$Q)
  )
)

plot_df <- index_df |>
  left_join(AIC_df) |>
  mutate(true = pull(true_index, 'biomass')) |>
  # Add 'rep' to group_by once additional sims are added
  group_by(sim_family, sim_link, fit_family, fit_link, Q, `_sdmTMB_time`) |>
  mutate(RMSE = sqrt(mean((log(est) - log(true))^2)),
         MRE = mean((est - true) / true), 
         covered = lwr < true & upr > true,
         title = ifelse(is.na(Q), sim_family, paste0(sim_family, ": Q=", signif(Q, digits = 2)))
  ) |>
  ungroup(fit_family) |>
  mutate(min_aic = min(aic),
         d_aic = aic - min_aic) |>
  ungroup()
title_levels <- unique(plot_df$title)
plot_df <- plot_df |>
  mutate(title = factor(title, levels = title_levels))

# - Compare: RMSE, MRE, coverage, AIC
# - look at the consequence of selecting the wrong family (what does AIC choose)
# - examine bias in estimates and how the above relates to the Q value

plot_df |>
  filter(Q != -5 | is.na(Q)) |>
ggplot(data = _, aes(x = RMSE, y = fit_family)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  facet_wrap(~ title)

plot_df |>
  filter(Q != -5 | is.na(Q)) |>
ggplot(data = _, aes(x = MRE, y = fit_family)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  facet_wrap(~ title)

plot_df |>
  filter(Q != -5 | is.na(Q)) |>
ggplot(data = _, aes(x = covered, y = fit_family)) +
  geom_point() +
  geom_vline(xintercept = 0.96, linetype = 'dashed') +
  facet_wrap(~ title)

plot_df |>
  filter(Q != -5 | is.na(Q)) |>
ggplot(data = _, aes(x = d_aic + 1, y = fit_family)) +
  geom_point() +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  scale_x_continuous(trans = 'log10') +
  facet_wrap(~ title, scales = 'free_x')
