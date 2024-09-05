library(dplyr)
library(ggplot2)
library(sdmTMB)
library(patchwork)
library(ggtext)
library(ggh4x)

theme_set(
  theme_light(base_size = 12) +
    theme(strip.background = element_rect(fill = "grey95"),
          strip.text = element_text(colour = "black"))
)
#theme_set(gfplot::theme_pbs())
source(here::here("R", "00-utils.R"))

plot_violin <- function(.data, .x, .ncol = NULL,
  .scale = "area", .adjust = 1,
  .summary_fun = "mean") {
  ggplot(data = .data, aes(x = {{ .x }}, y = fit_family)) +
    geom_violin(aes(colour = fit_family), scale = .scale, adjust = .adjust, bw = "SJ") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    stat_summary(fun = .summary_fun, geom = "point", colour = "black") +
    scale_color_manual(values = family_colours_no_delta) +
    facet_wrap(~title, ncol = .ncol) +
    guides(colour = "none") +
    theme(#axis.title.y = element_blank(),
          strip.text.x = ggtext::element_markdown())
}

plot_linedot <- function(.data, .x, .ncol = NULL, xint = 0.95) {
  ggplot(data = .data, aes(x = {{ .x }}, y = fit_family, colour = fit_family)) +
    geom_linerange(xmin = 0, mapping = aes(xmax = {{ .x }})) +
    geom_point(size = 3) +
    geom_vline(xintercept = xint, linetype = "dashed") +
    scale_color_manual(values = family_colours_no_delta) +
    facet_wrap(~title, ncol = .ncol) +
    guides(colour = "none")
}

relevel_fit_family <- function(df, fam_levels = c("gengamma", "tweedie", "gamma", "lognormal")) {
  df |>
    mutate(fit_family = tolower(fit_family)) |>
    mutate(fit_family = factor(fit_family, fam_levels)) |>
    mutate(fit_family = forcats::fct_recode(fit_family, "Tweedie" = "tweedie"))
}

filter_plot_df <- function(x, apply_filter = TRUE) {
  if (apply_filter) {
    x <- x |>
      filter(is.na(Q) | (Q %in% c(-1, 2)))
      # filter(is.na(Q) | (Q %in% c(NA, -1, -0.5, 0.001, 0.8, 2))) |>
      # filter(cv == 0.8, sigma_O == 1)
      # filter(is.na(Q) | (Q %in% c(-1, -0.5, 0.001, 0.5, 2))) |>
      # filter(cv == 0.95, sigma_O == 1)
  }
  x
}

out_dir <- here::here("data-outputs", "cross-sim")
fig_dir <- here::here("figures", "cross-sim")
fit_dir <- here::here(out_dir, "fit-summaries")
ind_dir <- here::here(out_dir, "index")

# Load data
# ------------------------------------------------------------------------------
Q_values <- dget(file.path(out_dir, "Q-values.txt"))

f <- file.path(ind_dir, list.files(ind_dir))

index_df <- readRDS(file.path(ind_dir, "cv0.95-sigmao1-nreps1000.rds"))
#index_df <- readRDS(file.path(ind_dir, "cv0.95-sigmao1-nreps50.rds"))
tag <- ""

# index_df <- readRDS(file.path(ind_dir, "cv0.8-sigmao1-nreps200-poisson-link.rds"))
# tag <- "-poisson-link"

apply_filter <- TRUE

title_lu <- tibble(
  title = c("gengamma (Q = -1)", "lognormal (Q = 0)", "gamma (Q = &sigma;)", "tweedie", "gengamma (Q = 2)"),
  sim_family = c("delta-gengamma", "delta-lognormal", "delta-gamma", "tweedie", "delta-gengamma"),
  Q = c(-1, NA, NA, NA, 2)
)

rqr_title_levels <- c(
  "gg: Q=-2", "gg: Q=-1", "gg: Q=-0.5",
  "gg: Q=-0.001", "lognormal", "gg: Q=0.001",
  "gg: Q=0.5", "gg: Q=0.95", "gamma", "gg: Q=1", "gg: Q=2",
  "tweedie"
)

# Prepare full plotting df
plot_df <- index_df |>
  group_by(rep, cv, sim_family, fit_family, Q, `_sdmTMB_time`) |> # keeping in case we add more time-steps
  mutate(
    re = (est - true) / true, # this is just relative error since we only have one time step
    lwr_ci_50 = exp(log_est + se * qnorm(0.25)),
    upr_ci_50 = exp(log_est + se * qnorm(0.75)),
    covered_50 = lwr_ci_50 < true & upr_ci_50 > true,
    covered = lwr < true & upr > true, # coverage is going to reflect the bias observed in the relative error
    ci_width = upr - lwr,
    ci_width_50 = upr_ci_50 - lwr_ci_50,
    title = ifelse(is.na(Q), sim_family, paste0(sim_family, ": Q=", signif(Q, digits = 2)))
  ) |>
  ungroup(fit_family) |>
  mutate(
    min_aic = min(aic),
    daic = aic - min_aic,
    rel_lik = exp(-0.5 * daic)
  ) |>
  mutate(aic_w = rel_lik / sum(rel_lik)) |> # see: https://atsa-es.github.io/atsa-labs/sec-uss-comparing-models-with-aic-and-model-weights.html
  ungroup() |>
  arrange(sim_family, Q) |>
  mutate(title = gsub("delta-", "", title)) |>
  mutate(title = gsub("gengamma", "gg", title)) |>
  mutate(title = gsub("-poisson-link", "", title))
title_levels <- unique(plot_df$title)
title_levels <- c(title_levels[-1], title_levels[1])
plot_df <- plot_df |>
  mutate(title = factor(title, levels = title_levels)) |>
  mutate(xtitle = paste0("cv = ", cv, "\nsigma_O = ", sigma_O)) |>
  relevel_fit_family()

# ------------------
# Main text figures
# ------------------

# Residuals
# ----------
rqr_df <- readRDS(file.path(here::here("data-outputs", "cross-sim", "rqr-101-cv0.95-sigmao1-b0-n2000.rds")))
rqr_df <- rqr_df |>
  filter(sanity_allok == TRUE) |>
  mutate(title = gsub("delta-", "", title)) |>
  mutate(title = gsub("gengamma", "gg", title)) |>
  mutate(title = gsub("-poisson-link", "", title)) |>
  relevel_fit_family()

# QQ plot
rqr_plot <- rqr_df |>
  select(-title) |>
  left_join(title_lu) |>
  filter(!is.na(title)) |>
  mutate(title = factor(title, levels = title_lu$title)) |>
  mutate(title = forcats::fct_recode(title, "Tweedie" = "tweedie")) |>
  ggplot(aes(sample = r)) +
  geom_qq() +
  geom_abline(intercept = 0, slope = 1) +
  facet_grid(fit_family ~ title) +
  coord_fixed() +
  scale_y_continuous(limits = c(-6.2, 6.2), breaks = seq(-6, 6, by = 2), labels = seq(-6, 6, by = 2), position = "left") +
  scale_x_continuous(limits = c(-4.6, 4.6), breaks = seq(-6, 6, by = 2), labels = seq(-6, 6, by = 2)) +
  labs(x = "Theoretical", y = "Sample") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        plot.background = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        tagger.panel.tag.background = element_blank(),
        strip.background = element_blank(),
        strip.text.x = ggtext::element_markdown()) +
  tagger::tag_facets()

wrap_elements(plot = grid::textGrob("Simulated distribution", vjust = 0), clip = FALSE) +
  wrap_elements(full = rqr_plot, clip = FALSE) +
  wrap_elements(plot = grid::textGrob("Fitted distribution", rot = 270, vjust = 1), clip = FALSE) +
  plot_layout(widths = c(50, 0.5), heights = c(0.5, 50),
    design = "A#
    BC")

ggsave(width = 7.3, height = 7.5, filename = file.path(fig_dir, "../figure-2-cross-sim-qq.png"))

# Relative error & AIC weight
# ---------------------------
main_text_df <-
  plot_df |>
  select(-title) |>
  left_join(title_lu) |>
  filter(!is.na(title)) |>
  mutate(title = factor(title, levels = title_lu$title)) |>
  mutate(title = forcats::fct_recode(title, "Tweedie" = "tweedie"))

re <- main_text_df |>
  plot_violin(.x = re, .ncol = 5, .summary_fun = "median") +
  scale_x_continuous(limits = c(-0.25, 1.25), breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  guides(x = guide_axis(title = "Relative error"),
    x.sec = guide_none(title = "Simulated"),
    y = guide_axis(title = "Fitted")
  ) +
  theme(axis.title.x = element_text(vjust = -1, size = 14, margin = margin(0, 0, 10, 2)),
        axis.title.x.top = element_text(vjust = 1, size = 12))

# AIC weights
aic_weight <- main_text_df |>
  filter_plot_df(apply_filter = apply_filter) |>
  plot_violin(.x = aic_w, .ncol = 5,
    .scale = "width",
    .summary_fun = "median") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1)) +
  guides(x = guide_axis(title = "AIC weight (%)"),
    x.sec = guide_none(title = "Simulated"),
    y = guide_axis(title = "Fitted")
  ) +
  theme(axis.title.x = element_text(vjust = -1, size = 14),
        axis.title.x.top = element_text(vjust = 1, size = 12))

((re / aic_weight) + plot_annotation(tag_levels = 'a', tag_suffix = ")")) &
  theme(plot.tag.position  = c(0, .9))
ggsave(width = 8, height = 5, filename = file.path(fig_dir, "../figure-3-cross-sim-RE-AIC-weight.png"))

# ------

# Full data / old
# - Compare: RMSE, MRE, coverage, AIC
# - look at the consequence of selecting the wrong family (what does AIC choose)
# - examine bias in estimates and how the above relates to the Q value
# - histogram is an option like Figure 5 in Thorson et al 2021 - surprising scale)

# rmse <- plot_df |>
#   filter_plot_df(apply_filter = apply_filter) |>
#   plot_violin(.x = rmse, .ncol = 10) +
#   scale_x_continuous(trans = "log10") +
#   labs(x = "RMSE", y = "Fit family") +
#   ggtitle(paste0("Root mean square error ", tag))
# ggsave(rmse, filename = file.path(fig_dir, paste0("rmse", tag, ".png")), width = 11, height = 2.5)

re <- plot_df |>
  filter_plot_df(apply_filter = apply_filter) |>
  plot_violin(.x = re, .ncol = 10) +
  labs(x = "Relative error", y = "Fit family") +
  lims(x = c(-0.255, 0.7)) +
  ggtitle(paste0("Relative error ", tag))
re
ggsave(re, filename = file.path(fig_dir, paste0("re", tag, ".png")), width = 11, height = 2.5)

aic_weight_hist <-
  plot_df |>
  filter_plot_df(apply_filter = apply_filter) |>
  ggplot() +
  geom_histogram(aes(x = aic_w)) +
  scale_x_continuous(labels = scales::label_percent(suffix = "")) +
  facet_grid(fit_family ~ title) +
  xlab("AIC Weight (%)")
aic_weight_hist

aic_weight <-
  plot_df |>
  filter_plot_df(apply_filter = apply_filter) |>
  ggplot(aes(x = aic_w, y = fit_family)) +
  geom_violin(aes(colour = fit_family), scale = "width") +
  stat_summary(fun = mean, geom = "point", colour = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_brewer(values = family_colours) +
  scale_x_continuous(labels = scales::label_percent(suffix = "")) +
  facet_wrap(~title, ncol = 10) +
  guides(colour = "none") +
  labs(x = "AIC Weight (%)", y = "Fit family") +
  ggtitle(paste0("AIC weight ", tag))
aic_weight

ggsave(filename = file.path(fig_dir, paste0("aic_weight", tag, ".png")), width = 11, height = 2.5)

ci_coverage <-
  plot_df |>
  filter_plot_df(apply_filter = apply_filter) |>
  group_by(title, fit_family, Q, cv, sigma_O, sim_family) |>
  summarise(
    n_sanity_pass = n(),
    prop_covered = sum(covered_50) / n_sanity_pass
  ) |>
  mutate(xtitle = paste0("cv = ", cv, "\nsigma_O = ", sigma_O)) |>
  plot_linedot(.x = prop_covered, .ncol = 10, xint = 0.5) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
  # facet_grid(xtitle ~ title) +
  labs(x = "Coverage", y = "Fit family") +
  ggtitle(paste0("50% Confidence interval coverage ", tag))
ci_coverage
ggsave(filename = file.path(fig_dir, paste0("ci-coverage", tag, ".png")), width = 11, height = 2.5)

ci_width <-
  plot_df |>
  filter_plot_df(apply_filter = apply_filter) |>
  plot_violin(.x = ci_width_50, .ncol = 10) +
  scale_x_continuous(trans = "log10") +
  #facet_grid(xtitle ~ title) +
  labs(x = "CI Width", y = "Fit family") +
  ggtitle(paste0("50% Confidence interval width ", tag))
ci_width
ggsave(filename = file.path(fig_dir, paste0("ci-width", tag, ".png")), width = 11, height = 2.5)

if (tag != "-poisson-link") {
  cross_combos <- bind_rows(
    tibble(Q = NA, sim_family = c("delta-lognormal", "delta-gamma", "tweedie")),
    tidyr::crossing(Q = Q_values, sim_family = "delta-gengamma")
  ) |>
    tidyr::crossing(fit_family = c("lognormal", "Gamma", "tweedie", "gengamma"))
} else {
  cross_combos <- bind_rows(
    tibble(Q = NA, sim_family = c("delta-lognormal-poisson-link", "delta-gamma-poisson-link", "tweedie")),
    tidyr::crossing(Q = Q_values, sim_family = "delta-gengamma-poisson-link")
  ) |>
    tidyr::crossing(fit_family = c("lognormal", "Gamma", "tweedie", "gengamma"))
}

sanity_tally <- left_join(
  cross_combos |> tidyr::unite(col = "sim_combo_key", sim_family, Q, fit_family, sep = ":", remove = TRUE),
  index_df |> tidyr::unite(col = "sim_combo_key", sim_family, Q, fit_family, sep = ":", remove = FALSE),
  by = c("sim_combo_key")
)

sanity_plot <- sanity_tally |>
  filter_plot_df() |>
  mutate(sanity_pass = ifelse(is.na(est), 0, 1)) |>
  group_by(sim_family, Q, cv, sigma_O, fit_family) |>
  summarise(n_pass = sum(sanity_pass), pass_prop = n_pass / max(rep), nreps = max(rep)) |>
  mutate(title = ifelse(is.na(Q), sim_family, paste0(sim_family, ": Q=", signif(Q, digits = 2)))) |>
  # mutate(title = factor(title, levels = title_levels)) |>
  mutate(title = gsub("delta-", "", title)) |>
  mutate(title = gsub("gengamma", "gg", title)) |>
  mutate(title = gsub("-poisson-link", "", title)) |>
  mutate(title = factor(title, levels = title_levels)) |>
  mutate(xtitle = paste0("cv = ", cv, "\nsigma_O = ", sigma_O)) |>
  relevel_fit_family() |>
  plot_linedot(.x = pass_prop, .ncol = 5, xint = 1.0) +
  # facet_grid(xtitle ~ title) +
  facet_wrap(~title, ncol = 10) +
  theme(axis.text.x = element_text(size = 9)) +
  scale_x_continuous(labels = scales::label_percent(suffix = "")) +
  labs(x = "Convergence rate (%)", y = "Fit family") +
  ggtitle(paste0("Model convergence rate ", tag))

sanity_plot
ggsave(filename = file.path(fig_dir, paste0("sanity-pass", tag, ".png")), width = 11, height = 2.5)

# ONLY RUN WHEN WANT SINGLE.
rmse / re
ggsave(filename = file.path(fig_dir, paste0("rmse-re", tag, ".png")), width = 11, height = 5)

ci_coverage + ci_width + plot_layout(nrow = 2)
ggsave(filename = file.path(fig_dir, paste0("ci-summ", tag, ".png")), width = 11, height = 5)

aic_weight / sanity_plot
ggsave(filename = file.path(fig_dir, paste0("aic_weight-sanity", tag, ".png")), width = 11, height = 5)


# ------------------------------------------------------------------------------
# Residuals
# ------------------------------------------------------------------------------
# rqr_df <- readRDS(file.path(here::here("data-outputs", "cross-sim", "rqr-101-cv0.95-sigmao1-b0-n2000.rds")))

# rqr_df <- rqr_df |>
#   filter(sanity_allok == TRUE) |>
#   mutate(title = gsub("delta-", "", title)) |>
#   mutate(title = gsub("gengamma", "gg", title)) |>
#   mutate(title = gsub("-poisson-link", "", title))
# rqr_title_levels <- unique(rqr_df$title)
# rqr_title_levels <- c(rqr_title_levels[-1], rqr_title_levels[1])

# Check gengamma rqr
# idx <- fit_df |>
#   filter(sim_family == "delta-gengamma", fit_family == "gengamma", sanity_allok == TRUE) |>
#   pull(id)
rqr_plot <- rqr_df |>
  mutate(title = factor(gsub("gg: ", "gengamma\n", title), levels = gsub("gg: ", "gengamma\n", rqr_title_levels))) |>
  mutate(title = forcats::fct_recode(title, "lognormal\n" = "lognormal", "gamma\n" = "gamma", "Tweedie\n" = "tweedie")) |>
  relevel_fit_family() |>
  ggplot(aes(sample = r)) +
  geom_qq() +
  geom_abline(intercept = 0, slope = 1) +
  facet_grid(fit_family ~ title)
rqr_plot

ggsave(rqr_plot, filename = file.path(fig_dir, paste0("rqr-", tag, ".png")),
  width = 15, height = 7
)

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


# ------------------------------------------------------------------------------
# Get back ran pars?
# ------------------------------------------------------------------------------
# fit_df <- file.path(fit_dir, list.files(fit_dir)) |>
#   purrr::map_dfr(readRDS)
fit_df <- readRDS(file.path(fit_dir, "cv0.8-sigmao1-nreps200.rds"))

fit_df |>
  # filter(sim_family == "delta-gengamma", fit_family == "gengamma") |>
  group_by(sim_family, Q, fit_family) |>
  summarise()

mean_ests <- fit_df |>
  group_by(sim_family, Q, fit_family) |>
  summarise_at(vars(est_Q:std.error_tweedie_p), mean, na.rm = TRUE, groups = ".drop")

# fit_df |>
mean_ests |>
  filter(sim_family == "delta-gengamma", fit_family == "gengamma") |>
  filter(!(Q %in% c(-5, 5))) |>
  select(sim_family:est_Qse) |>
  ggplot(aes(x = Q, y = est_Q)) +
  geom_point(aes(y = est_Q + 1.96 * est_Qse), stroke = 3, size = 6, shape = "-", colour = "#2c7bb6", alpha = 1) +
  geom_point(aes(y = est_Q - 1.96 * est_Qse), stroke = 3, size = 6, shape = "-", colour = "#D7191C", alpha = 1) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Simulation Q", y = "Estimated Q")
ggsave(filename = file.path(fig_dir, paste0("Q-estimation", tag, ".png")), width = 4.5, height = 4)
# Notes:
# - cannot recover Q < -1 very well (I'm assuming because there is just too little
# information at the tails)
# - similarly cannot seem to recover Q at Q = 5
