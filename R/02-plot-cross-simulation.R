library(dplyr)
library(ggplot2)

theme_set(theme_light(base_size = 12))

plot_violin <- function(.data, .x, .ncol = NULL) {
  ggplot(data = .data, aes(x = {{.x}}, y = fit_family)) +
    stat_summary(fun = mean, geom = "point") +
    geom_violin(aes(col = fit_family), alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    scale_color_brewer(palette = "Dark2") +
    facet_wrap(~ title, ncol = .ncol) +
    guides(colour = 'none')
}

plot_linedot <- function(.data, .x, .ncol = NULL, xint = 0.95) {
  ggplot(data = .data, aes(x = {{.x}}, y = fit_family, colour = fit_family)) +
  geom_linerange(xmin = 0, mapping = aes(xmax = {{.x}})) +
  geom_point(size = 3) +
  geom_vline(xintercept = xint, linetype = 'dashed') +
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(~ title, ncol = .ncol) +
  guides(colour = 'none')
}

fig_dir <- here::here('figures')
fit_dir <- here::here('data-outputs', 'fit-summaries')
ind_dir <- here::here('data-outputs', 'index')

# Load data
# ------------------------------------------------------------------------------
Q_values <- dget(file.path('data-outputs', 'Q-values.txt'))

relevel_fit_family <- function(df, fam_levels = c('gengamma', 'gamma', 'tweedie', 'lognormal' )) {
  df |>
  mutate(fit_family = tolower(fit_family)) |>
  mutate(fit_family = factor(fit_family, fam_levels))
}

# Look at real data to figure out what these should be.
Qs_to_plot <- c(-1, 0.001, 0.5, 0.8, 1, 2)

index_df <- file.path(ind_dir, list.files(ind_dir)) |>
  purrr::map_dfr(readRDS) |>

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
  group_by(rep, cv, sim_family, fit_family, Q, `_sdmTMB_time`) |>
  mutate(RMSE = sqrt(mean((log(est) - log(true))^2)),
         MRE = mean((est - true) / true),
         covered = lwr < true & upr > true,
         ci_width = upr - lwr,
         title = ifelse(is.na(Q), sim_family, paste0(sim_family, ": Q=", signif(Q, digits = 2)))
  ) |>
  ungroup(fit_family) |>
  mutate(min_aic = min(aic),
         d_aic = aic - min_aic) |>
  ungroup() |>
  arrange(sim_family, Q)
title_levels <- unique(plot_df$title)
title_levels <- c(title_levels[-1], title_levels[1])
plot_df <- plot_df |>
  mutate(title = factor(title, levels = title_levels)) |>
  filter(!(Q %in% c(-5, -0.001, 5))) |>
  mutate(xtitle = paste0('cv = ', cv, "\nsigma_O = ", sigma_O)) |>
  relevel_fit_family()


# - Compare: RMSE, MRE, coverage, AIC
# - look at the consequence of selecting the wrong family (what does AIC choose)
# - examine bias in estimates and how the above relates to the Q value

tag <- ""
plot_violin(plot_df, .x = RMSE, .ncol = 5) +
  ggtitle("RMSE") +
  facet_grid(xtitle ~ title)
ggsave(filename = file.path(fig_dir, paste0('rmse', tag, '.png')), width = 11, height = 6.5)


plot_violin(plot_df, .x = MRE, .ncol = 5) +
  ggtitle("MRE") +
  facet_grid(xtitle ~ title)
ggsave(filename = file.path(fig_dir, paste0('mre', tag, '.png')), width = 11, height = 6.5)

plot_violin(plot_df, .x = (d_aic + 1), .ncol = 5) +
  scale_x_continuous(trans = 'log10') +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  facet_grid(xtitle ~ title) +
  ggtitle("Delta AIC")

ggplot(plot_df, aes(x = (d_aic + 1), y = fit_family)) +
  geom_point(alpha = 0.5, colour = 'grey50') +
  stat_summary(fun = mean, geom = "point") +
  scale_x_continuous(trans = 'log10') +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  facet_grid(xtitle ~ title) +
  ggtitle("Delta AIC")

ggsave(filename = file.path(fig_dir, paste0('daic', tag, '.png')), width = 11, height = 6.5)

plot_df |>
  group_by(title, fit_family, Q, cv, sigma_O, sim_family) |>
  summarise(n_sanity_pass = n(),
            prop_covered = sum(covered) / n_sanity_pass
         ) |>
  mutate(xtitle = paste0('cv = ', cv, "\nsigma_O = ", sigma_O)) |>
plot_linedot(.x = prop_covered, .ncol = 5) +
  ggtitle("95% CI Coverage") +
  facet_grid(xtitle ~ title)
ggsave(filename = file.path(fig_dir, paste0('ci-coverage', tag, '.png')), width = 11, height = 6.5)

plot_violin(plot_df, .x = ci_width, .ncol = 5) +
  scale_x_continuous(trans = 'log10') +
  facet_grid(xtitle ~ title) +
  ggtitle("95% CI Width")
ggsave(filename = file.path(fig_dir, paste0('ci-width', tag, '.png')), width = 11, height = 6.5)

nreps <- 50
sanity_tally |>
  mutate(sanity_pass = ifelse(is.na(est), 0, 1)) |>
  group_by(sim_family, Q, cv, sigma_O, fit_family) |>
  summarise(n_pass = sum(sanity_pass), pass_prop = n_pass / nreps) |>
  mutate(title = ifelse(is.na(Q), sim_family, paste0(sim_family, ": Q=", signif(Q, digits = 2)))) |>
  mutate(title = factor(title, levels = title_levels)) |>
  mutate(xtitle = paste0('cv = ', cv, "\nsigma_O = ", sigma_O)) |>
plot_linedot(.x = pass_prop, .ncol = 5, xint = 1.0) +
  facet_grid(xtitle ~ title) +
  ggtitle("Passed sanity check")
ggsave(filename = file.path(fig_dir, paste0('sanity-pass', tag, '.png')), width = 11, height = 6.5)

# ------------------------------------------------------------------------------
# Get back ran pars?
# ------------------------------------------------------------------------------
fit_df <- file.path(fit_dir, list.files(fit_dir)) |>
  purrr::map_dfr(readRDS)

fit_df |>
  #filter(sim_family == "delta-gengamma", fit_family == "gengamma") |>
  group_by(sim_family, Q, fit_family) |>
  summarise()

mean_ests <- fit_df |>
  group_by(sim_family, Q, fit_family) |>
  summarise_at(vars(est_Q:std.error_tweedie_p), mean, na.rm = TRUE, groups = ".drop")

#fit_df |>
mean_ests |>
  filter(sim_family == 'delta-gengamma', fit_family == 'gengamma') |>
  select(sim_family:est_Qse) |>
ggplot(aes(x = Q, y = est_Q)) +
  geom_point(aes(y = est_Q + 1.96 * est_Qse), stroke = 3, size = 6, shape = "-", colour = '#2c7bb6', alpha = 1) +
  geom_point(aes(y = est_Q - 1.96 * est_Qse), stroke = 3, size = 6, shape = "-", colour = '#D7191C', alpha = 1) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Simulation Q", y = "Estimated Q")
ggsave(filename = file.path(fig_dir, paste0('Q-estimation', tag, '.png')), width = 4.5, height = 4)
# Notes:
# - cannot recover Q < -1 very well (I'm assuming because there is just too little
# information at the tails)
# - similarly cannot seem to recover Q at Q = 5