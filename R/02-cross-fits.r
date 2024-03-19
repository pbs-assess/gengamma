library(dplyr)
library(ggplot2)
library(purrr)
library(sdmTMB)
theme_set(theme_light())

source(here::here('R', '00-utils.R'))
fit_dir <- file.path(here::here("data-outputs", "fits"))

.families <- c('delta-lognormal', 'delta-gamma', 'tweedie', 'delta-gengamma')

load(file.path(here::here("data-outputs", "sim-families.RData")))
sim_df <- bind_rows(dl_sim, dg_sim, tw_sim, dgg_sim)
cross_combos <- distinct(sim_df, family, Q) |>
  tidyr::crossing(.fit_fam = .families) |>
  rename(.sim_fam = 'family', .Q = "Q")

# ------------------------------------------------------------------------------
# Fit cross-simulations
fits <- cross_combos |>
  pmap(\(.sim_fam, .Q, .fit_fam) {
    fit_cross(
    .data = sim_df, # data filtering built into `fit_cross()`
    .mesh = sampled_mesh, 
    sim_fam = .sim_fam,
    fit_fam = .fit_fam, 
    .Q = .Q)
  }
)
# saveRDS(fits, file.path(fit_dir, 'fits.rds'))
# fits <- readRDS(file.path(fit_dir, 'fits.rds'))

fit_sanity <- map_dfr(fits, \(x) 
  tibble(
    sim_family = unique(x$data$family),
    fit_family1 = family(x)[[1]][[1]],
    fit_family2 = ifelse(family(x)[[1]][[1]] == 'tweedie', 'tweedie', family(x)[[2]][[1]]),
    Q = unique(x$data$Q),
    sanity_allok = sanity(x)$all_ok)
  )

# Only get predictions for cases where sanity checks passed
sanity_pass <- fit_sanity$sanity_allok
pred_list <- fits[sanity_pass] |>
  map(predict, newdata = predictor_dat, return_tmb_object = TRUE)

index_df <- map_dfr(pred_list, get_index_summary)
index_df
#beepr::beep()
saveRDS(index_df, file.path(here::here("data-outputs"), 'cross-fit-index-df.rds')
index_df <- readRDS(file.path(here::here("data-outputs"), 'cross-fit-index-df.rds'))

AIC_df <- map_dfr(fits, \(fit) tibble(
  aic = AIC(fit),
  sim_family = unique(fit$data$family),
  fit_family = ifelse(family(x)[[1]][[1]] == 'tweedie', 'tweedie', family(x)[[2]][[1]]), 
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

# Speed test
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

speed_test <- rbenchmark::benchmark(
  "ln; dl" = benchmark_expr(sim_dat = dl_sim, fit_family = sdmTMB::delta_lognormal()),
  "ln; dg" = benchmark_expr(sim_dat = dl_sim, fit_family = sdmTMB::delta_gamma()),
  "ln; tw" = benchmark_expr(sim_dat = dl_sim, fit_family = sdmTMB::tweedie()),
  "ln; dgg" = benchmark_expr(sim_dat = dl_sim, fit_family = sdmTMB::delta_gengamma()),
  "dg; dl" = benchmark_expr(sim_dat = dg_sim, fit_family = sdmTMB::delta_lognormal()),
  "dg; dg" = benchmark_expr(sim_dat = dg_sim, fit_family = sdmTMB::delta_gamma()),
  "dg; tw" = benchmark_expr(sim_dat = dg_sim, fit_family = sdmTMB::tweedie()),
  "dg; dgg" = benchmark_expr(sim_dat = dg_sim, fit_family = sdmTMB::delta_gengamma()),
  "tw; dl" = benchmark_expr(sim_dat = tw_sim, fit_family = sdmTMB::delta_lognormal()),
  "tw; dg" = benchmark_expr(sim_dat = tw_sim, fit_family = sdmTMB::delta_gamma()),
  "tw; tw" = benchmark_expr(sim_dat = tw_sim, fit_family = sdmTMB::tweedie()),
  "tw; gg" = benchmark_expr(sim_dat = tw_sim, fit_family = sdmTMB::delta_gengamma()),
  "dgg Q = 0.5; dl" = benchmark_expr(sim_dat = dgg_sim |> filter(Q == 0.5), fit_family = sdmTMB::delta_lognormal()),
  "dgg Q = 0.5; dg" = benchmark_expr(sim_dat = dgg_sim |> filter(Q == 0.5), fit_family = sdmTMB::delta_gamma()),
  "dgg Q = 0.5; tw" = benchmark_expr(sim_dat = dgg_sim |> filter(Q == 0.5), fit_family = sdmTMB::tweedie()),
  "dgg Q = 0.5; dgg" = benchmark_expr(sim_dat = dgg_sim |> filter(Q == 0.5), fit_family = sdmTMB::delta_gengamma()),
  replications = 10, 
  columns = c("test", "replications", "elapsed",
                      "relative", "user.self", "sys.self")
) |>
  tidyr::separate(test, into = c('sim_family', 'fit_family'), sep = '; ')
beepr::beep()