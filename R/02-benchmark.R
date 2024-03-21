library(dplyr)
library(ggplot2)

theme_set(theme_light(base_size = 12))

sim_list <- readRDS(file.path(here::here('data-outputs'), 'sim-list.rds'))
sim_df <- sim_list$sim_df
sampled_mesh <- sim_list$sampled_mesh

dl_sim <- filter(sim_df, family == 'delta-lognormal')
dg_sim <- filter(sim_df, family == 'delta-gamma')
tw_sim <- filter(sim_df, family == 'tweedie')
dgg_sim <- filter(sim_df, family == 'delta-gengamma')

family_label <- tibble(
  fit_family = c('dl', 'dg', 'tw', 'dgg'), 
  family = c('lognormal', 'gamma', 'tweedie', 'gengamma'))

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
  "tw; dgg" = benchmark_expr(sim_dat = tw_sim, fit_family = sdmTMB::delta_gengamma()),
  "dgg Q = 0.001; dl" = benchmark_expr(sim_dat = dgg_sim |> filter(Q == 0.001), fit_family = sdmTMB::delta_lognormal()),
  "dgg Q = 0.001; dg" = benchmark_expr(sim_dat = dgg_sim |> filter(Q == 0.001), fit_family = sdmTMB::delta_gamma()),
  "dgg Q = 0.001; tw" = benchmark_expr(sim_dat = dgg_sim |> filter(Q == 0.001), fit_family = sdmTMB::tweedie()),
  "dgg Q = 0.001; dgg" = benchmark_expr(sim_dat = dgg_sim |> filter(Q == 0.001), fit_family = sdmTMB::delta_gengamma()),
  "dgg Q = 0.5; dl" = benchmark_expr(sim_dat = dgg_sim |> filter(Q == 0.5), fit_family = sdmTMB::delta_lognormal()),
  "dgg Q = 0.5; dg" = benchmark_expr(sim_dat = dgg_sim |> filter(Q == 0.5), fit_family = sdmTMB::delta_gamma()),
  "dgg Q = 0.5; tw" = benchmark_expr(sim_dat = dgg_sim |> filter(Q == 0.5), fit_family = sdmTMB::tweedie()),
  "dgg Q = 0.5; dgg" = benchmark_expr(sim_dat = dgg_sim |> filter(Q == 0.5), fit_family = sdmTMB::delta_gengamma()),
  replications = 10, 
  columns = c("test", "replications", "elapsed",
                      "relative", "user.self", "sys.self")
) |>
  tidyr::separate(test, into = c('sim_family', 'fit_family'), sep = '; ')
beep()

speed_test |>
  group_by(fit_family) |>
  mutate(max_elapsed = max(elapsed)) |>
  left_join(family_label) |>
ggplot(aes(x = forcats::fct_reorder(family, max_elapsed), y = elapsed)) + 
  geom_point() +
  labs(x = "Fit Family", y = "Time Elapsed (s)")
ggsave(filename = file.path(here::here('figures'), 'benchmark.png'), width = 4.5, height = 4)
