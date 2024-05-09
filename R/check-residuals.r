# library(sdmTMB)
devtools::load_all("../sdmTMB") # add-gengamma-rqr branch
library(dplyr)
library(ggplot2)
theme_set(theme_light())

source(here::here("R", "00-utils.R"))

# ------------
# fits <- readRDS(here::here('data-outputs', 'fits', '42-cv0.8-sigmao0-b0-n1000.rds'))
fits <- readRDS(here::here("data-outputs", "cross-sim", "fits", "42-cv0.95-sigmao1-b0-n1000-poisson-link.rds"))

sanity_df <- purrr::map_dfr(fits, ~ get_sanity_df(.x, silent = TRUE))

fit_df <- purrr::map_dfr(fits, get_fitted_estimates) |>
  mutate(id = row_number()) |>
  mutate(title = ifelse(is.na(Q), sim_family, paste0(sim_family, ": Q=", signif(Q, digits = 2)))) |>
  arrange(sim_family, Q) |>
  left_join(sanity_df)
title_levels <- c(unique(fit_df$title)[-1], "delta-gamma")

idx <- fit_df |>
  # filter(sim_family == "delta-gengamma", Q == -2, fit_family == "gengamma") |>
  filter(is.na(Q) | !(Q %in% c(-5, 5))) |>
  pull(id)

# RQR
# ------------------
# rqr_df <- purrr::map_dfr(1:length(fits), ~ get_rqr(fits[[.x]], id = .x)) |>
#   left_join(fit_df)
# beep()

rqr_df <- purrr::map_dfr(idx, ~ get_rqr(fits[[.x]], id = .x)) |>
  left_join(fit_df)
beep()


# Check gengamma rqr
# idx <- fit_df |>
#   filter(sim_family == "delta-gengamma", fit_family == "gengamma", sanity_allok == TRUE) |>
#   pull(id)

rqr_df |>
  # filter(id %in% idx) |>
  filter(sanity_allok == TRUE) |>
  mutate(title = factor(title, levels = title_levels)) |>
  filter(!Q %in% c(-5, 5)) |>
  ggplot(aes(sample = r)) +
  geom_qq() +
  geom_abline(intercept = 0, slope = 1) +
  facet_grid(fit_family ~ title)

# DHARMa residuals
# ------------------
dr_df <- purrr::map_dfr(fit_df$id[fit_df$sanity_allok], ~ get_dr(fits[[.x]], .x, seed = 100, nsim = 200)) |>
  tibble::as_tibble()
beep()

left_join(dr_df, fit_df) |>
  mutate(title = factor(title, levels = title_levels)) |>
  plot_resid()

left_join(dr_df, fit_df) |>
  filter(title == "delta-gengamma: Q=0.8") |>
  plot_resid()

# dr_df_mleeb <- purrr::map_dfr(fit_df$id[fit_df$sanity_allok], ~ get_dr(fits[[.x]], .x, seed = 100, nsim = 200, type = 'mle-eb')) |>
#   tibble::as_tibble()
# beep()

# left_join(dr_df_mleeb, fit_df) |>
#   mutate(title = factor(title, levels = title_levels)) |>
#   plot_resid()
