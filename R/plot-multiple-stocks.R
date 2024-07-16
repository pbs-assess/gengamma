library(dplyr)
library(ggplot2)
library(janitor)
library(sdmTMB)
library(patchwork)
theme_set(gfplot::theme_pbs())

source("R/00-utils.R")

dir.create(here::here("figures", "bc-gf-data"), showWarnings = FALSE, recursive = TRUE)
fig_dir <- here::here("figures", "bc-gf-data")
df_dir <- here::here("data-outputs", "bcgf-outputs")

family_levels <- c(
  "delta-gengamma", "delta-lognormal", "delta-gamma", "tweedie",
  "delta-lognormal-poisson-link", "delta-gamma-poisson-link", "delta-gengamma-poisson-link"
)

# https://mk.bcgsc.ca/colorblind/palettes/8.color.blindness.palette.txt
family_colours <- c(
  "tweedie" = "#2271B2", #' honolulu blue'
  "delta-lognormal" = "#359B73", #' ocean green'
  "delta-gamma" = "#d55e00", # 'bamboo'
  "delta-gengamma" = "#F748A5"
) # 'barbie pink'
# family_shapes <- c('tweedie-not-converged' = 0,
#                    'delta-lognormal-not-converged' = 1,
#                    'delta-gamma-not-converged' = 2,
#                    'delta-gengamma-not-converged' = 5,
#                    'tweedie' = 15,
#                    'delta-lognormal' = 16,
#                    'delta-gamma' = 17,
#                    'delta-gengamma' = 18)

spp_list <- gfsynopsis::get_spp_names() |>
  select(species = species_common_name, species_code, species_science_name)

pos_sets <- readRDS(here::here("data", "clean-survey-data.rds")) |>
  distinct(species, region = survey_abbrev, mean_pos, mean_sets, mean_pos_sets, prop_pos) |>
  mutate(mean_pos = round(mean_pos), mean_sets = round(mean_sets))

lu_df <- readRDS(file.path(df_dir, "lu-df.rds")) |>
  janitor::clean_names() |>
  left_join(spp_list) |>
  left_join(pos_sets)
lu_df <- lu_df |>
  filter(type != "poisson-link") |>
  filter(n_regions == 3)

# fit estimates of models that converged
fit_ests <- readRDS(file.path(df_dir, "fitted-estimates-all.rds")) |>
  janitor::clean_names() |>
  mutate(type = gsub("poisson_link_delta", "poisson-link", type)) |>
  mutate(type = ifelse(fit_family == "tweedie", "standard", type)) |>
  right_join(lu_df)

# sanity summary for all models fit (not necessarily converged)
sanity_df <- readRDS(file.path(df_dir, "sanity-df-all.rds")) |>
  janitor::clean_names() |>
  right_join(lu_df)

summary_df <- left_join(fit_ests, sanity_df)

rqr_df <- readRDS(file.path(df_dir, "rqr-catch-df-all.rds")) |>
  right_join(lu_df) |>
  mutate(family = factor(family, levels = family_levels))

# indices of models that converged
index_df <- readRDS(file.path(df_dir, "index-df-all.rds")) |>
  janitor::clean_names() |>
  select(-type) |>
  right_join(lu_df) |>
  as_tibble()

# design index
design_df <- readRDS(here::here("data", "survey-design-index.rds"))
design_df <- design_df |>
  select(year, est = biomass, lwr = lowerci, upr = upperci, se = re,
         species = species_common_name, region = survey_abbrev) |>
  right_join(distinct(lu_df, species, region)) |>
  mutate(family = "design")

# Number of models that converged
sanity_df |>
  tabyl(family, all_ok)

# Frequency of estimated Q
q_ests <- summary_df |>
  filter(family == "delta-gengamma") |>
  select(species, region, est_q, est_qse) |>
  arrange(est_q) |>
  mutate(q_rank = row_number())

hist(q_ests$est_q)
hist(q_ests$est_q, breaks = "FD")

fd_breaks <- pretty(range(q_ests$est_q), n = nclass.FD(q_ests$est_q), min.n = 1)
summary_df |>
  filter(family == "delta-gengamma") |>
ggplot(aes(x = est_q)) +
  geom_histogram(bins = 16) +
  scale_x_continuous(breaks = seq(-1.5, 1.5, 0.25), labels = seq(-1.5, 1.5, 0.25)) +
  # geom_histogram(breaks = fd_breaks, labels = fd_breaks) +
  labs(x = "Q estimates", y = "Count")

# ------
aic_df <- summary_df |>
  group_by(species, region) |>
  mutate(
    min_aic = min(aic),
    daic = aic - min_aic,
    rel_lik = exp(-0.5 * daic)
  ) |>
  mutate(aic_w = rel_lik / sum(rel_lik)) |> # see: https://atsa-es.github.io/atsa-labs/sec-uss-comparing-models-with-aic-and-model-weights.html
  ungroup() |>
  left_join(q_ests) |>
  #arrange(species)
  arrange(q_rank)

aic_plot_df <- aic_df |>
  mutate(species = gsub("north", "", species)) |>
  mutate(species = stringr::str_to_title(species)) |>
  mutate(
    species_id = as.numeric(factor(species)),
    odd_species = ifelse(species_id %% 2 == 0, 1, 0),
    converged = ifelse(all_ok, 1, 0),
    x = case_when(
      !all_ok ~ 10^-0.6,
      # !all_ok & family == "tweedie" ~ 10^-1,
      # !all_ok & family == "delta-lognormal" ~ 10^-0.867,
      # !all_ok & family == "delta-gamma" ~ 10^-0.733,
      # !all_ok & family == "delta-gengamma" ~ 10^-0.6,
      TRUE ~ daic + 1
    )
  ) |>
  mutate(fspecies = forcats::fct_rev(species)) |>
  # mutate(fspecies2 = forcats::fct_reorder(fspecies, as.numeric(species_code))) |>
  mutate(fspecies2 = forcats::fct_reorder(fspecies, as.numeric(q_rank))) |>
  mutate(odd_species2 = ifelse(as.numeric(fspecies2) %% 2 == 0, 1, 0)) |>
  mutate(family = factor(family, levels = family_levels))

aic_plot <-
  ggplot() +
  gfplot::theme_pbs(base_size = 12) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(vjust = 0),
        axis.text = element_text(size = 9),
        panel.spacing.x = unit(0.2, "lines"),
        panel.border = element_rect(fill = NA, linewidth = 0.5),
        legend.title = element_blank(),
        legend.position = "top",
        legend.margin = margin(0, 0, -0.3, 0, "cm")) +
  coord_cartesian(clip = "off", expand = FALSE) +
  facet_wrap(~region, ncol = 3, drop = FALSE) +
  geom_tile(
    data = aic_plot_df |> distinct(species, region, .keep_all = TRUE),
    mapping = aes(
      x = 1, y = fspecies2,
      width = Inf, height = 1, fill = factor(odd_species2)
    )
  ) +
  scale_fill_manual(values = c("grey95", "white")) +
  annotate(geom = "rect", xmin = 0.6, xmax = 2, ymin = -Inf, ymax = Inf, fill = "white") +
  geom_point(
    data = aic_plot_df |> filter(all_ok), # models that converged
    aes(x = x, y = species, colour = family, shape = family),
    stroke = 0.5, size = 2, position = ggstance::position_dodgev(height = 0)
  ) +
  scale_colour_manual(values = family_colours) +
  labs(x = "Delta AIC + 1", y = "Species",
       #shape = "Family",
       colour = "Family") +
  guides(fill = "none") +
  scale_x_continuous(
    trans = "log10",
    labels = scales::label_number(trim = TRUE, accuracy = 1),
    limits = c(0.6, 3000)
  ) +
  annotate(geom = "text", x = 1600, y = length(unique(aic_plot_df$species)),
           vjust = -2, hjust = 0.75, label = "Q", colour = "grey30") +
  geom_text(data = filter(aic_plot_df, family == "delta-gengamma"), aes(
    x = 1600, y = fspecies2,
    label = round(est_q, digits = 1)
  ), hjust = 0.75, size = 3, colour = "grey50")
  # geom_point(
  #   data = aic_plot_df |> filter(!all_ok), # models that did not converge
  #   aes(x = x, y = species, colour = family),
  #   shape = 21, stroke = 0.7, size = 2,
  #   position = ggstance::position_dodgev(height = 0.5)
  # ) +
  # geom_text(
  #   data = aic_plot_df |> distinct(species, region, .keep_all = TRUE),
  #   aes(x = 0.5, y = fspecies2, label = round(mean_pos)), # paste0(round(mean_pos), "(", prop_pos, ")")),
  #   hjust = 1, colour = "grey50",
  #   size = 3
  # ) +
  # geom_text(data = distinct(aic_plot_df, region, mean_sets), aes(
  #   x = 0.5, y = length(unique(aic_plot_df$species)),
  #   label = mean_sets
  # ), vjust = -3.5, hjust = 1, size = 3) +
  # geom_text(data = aic_plot_df |> filter(prop_pos >= 0.05), aes(x = 0.04, y = fspecies2, label = signif(est_q, 2)),
  #           size = 2.5, hjust = 0) +
aic_plot
#aic_plot_filename <- "aic-plot-spp-code-order.png"
aic_plot_filename <- "aic-plot-qest-order.png"
ggsave(aic_plot, width = 7.5, height = 6, filename = file.path(fig_dir, aic_plot_filename))


# Plot indices
# -----------------
#p1 <-
pind <- index_df |>
  left_join(select(aic_df, fname, min_aic:q_rank)) |>
  filter(species == "lingcod")

ggplot(data = pind, aes(x = year, y = est, group = family)) +
    geom_line(aes(colour = family)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = family), alpha = 0.3) +
    geom_text()
    scale_y_continuous(trans = "log10") +
    facet_grid(species ~ region, scales = "free_y") +
    scale_colour_manual(values = family_colours) +
    scale_fill_manual(values = family_colours)

p2 <- rqr_df |>
  filter(species == "lingcod") |>
  mutate(fit_family = factor(tolower(fit_family), gsub("delta-", "", family_levels))) |>
  ggplot(aes(sample = r)) +
  geom_qq() +
  geom_abline(intercept = 0, slope = 1) +
  facet_grid(fit_family ~ region, scales = "free", switch = "y") +
  scale_y_continuous(breaks = seq(-6, 6, by = 2), labels = seq(-6, 6, by = 2), position = "right") +
  labs(x = "Theoretical", y = "Sample") +
  theme(strip.text.x = element_blank())

p1 / p2

# ------------------------------------------------------------------------------
# Will anyone want to see the estimated CV of the different fits to real data?
get_ln_phi <- function(model) {
  m <- ifelse(family(model)[[1]][[1]] == 'tweedie', 1, 2)
  tibble(
    ln_phi = model$parlist$ln_phi[[m]],
    species = unique(model$data$species),
    region = unique(model$data$survey_abbrev),
    spatial = unique(model$spatial),
    spatiotemporal = unique(model$spatiotemporal),
    fit_family = family(model)[[m]][[1]],
    type = family(model)$type
  )
}

f <- list.files(file.path(df_dir, c("fits_sp-off-st-iid", "fits_sp-off-st-off", "fits_sp-on-st-iid", "fits_sp-on-st-off")),
  full.names = TRUE)

ln_phi_df <- f |>
  purrr::map(~ {
    model <- readRDS(.x)
    try(get_ln_phi(model))
  }) |>
  purrr::keep(is.data.frame) |>
  bind_rows()

filter(ln_phi_df, fit_family == "Gamma", type == "standard") |>
  left_join(sanity_df) |>
  filter(all_ok) |>
  mutate(cv = get_gamma_cv(phi = exp(ln_phi))) |>
  pull('cv') |>
  hist()

get_gamma_cv <- function(phi) 1 / sqrt(phi)



# ggplot(data = _, aes(x = daic + 1, y = forcats::fct_rev(species))) +
#   geom_tile(data = filter(failed_df, species_id %% 2 == 0), aes(width = Inf, height = 1, fill = 'grey80'), colour = NA, alpha = 0.3) +
#   scale_fill_manual(values = c("grey", "white")) +
#   geom_jitter(aes(colour = family), width = 0, height = 0.5) +
#   geom_point(data = failed_df |> filter(is.na(daic), region != "SYN WCHG"),
#     aes(x = family_id, y = forcats::fct_rev(species), colour = family), shape = 21) +
#   scale_x_continuous(trans = "log10") +
#   facet_wrap(~ region, ncol = 3)
# # Question: do we ever fit delta models where spatial or spatiotemporal setting is different?
