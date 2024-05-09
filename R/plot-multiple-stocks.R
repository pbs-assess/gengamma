library(dplyr)
library(ggplot2)
library(janitor)
theme_set(gfplot::theme_pbs())

source("R/00-utils.R")

dir.create(here::here("figures", "bc-gf-data"), showWarnings = FALSE, recursive = TRUE)
fig_dir <- here::here("figures", "bc-gf-data")

dat <- readRDS(file.path("data-outputs", "clean-survey-data.rds")) |>
  rename(region = "survey_abbrev")
spp_lu <- distinct(dat, species, species_code)

mean_pos_sets <- dat |>
  group_by(species, region, year) |>
  summarise(pos = sum(present), n_sets = n()) |>
  summarise(mean_pos = mean(pos), mean_sets = mean(n_sets), .groups = "drop") |>
  mutate(
    mean_pos_sets = paste0(round(mean_pos), "/", round(mean_sets)),
    prop_pos = round(mean_pos / mean_sets, digits = 2)
  )

spp_region_keep <- mean_pos_sets |>
  filter(prop_pos >= 0.05) |>
  select(species, region)

df_dir <- here::here("data-outputs", "bcgf-outputs")

family_levels <- c("tweedie", "delta-lognormal", "delta-gamma", "delta-gengamma")

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


lu_df <- readRDS(file.path(df_dir, "lu-df.rds")) |>
  janitor::clean_names() |>
  filter(species != "longspine thornyhead") # never has data for the 3 main regions
# mutate(family_id = as.numeric(factor(family, family_levels))) # all models fit

# indices of models that converged
indices <- readRDS(file.path(df_dir, "index-df.rds")) |>
  tibble::as_tibble() |>
  janitor::clean_names() |>
  filter(species != "longspine thornyhead") # never has data for the 3 main regions
# fit estimates of models that converged
fit_ests <- readRDS(file.path(df_dir, "fitted-estimates-df.rds")) |>
  janitor::clean_names() |>
  filter(species != "longspine thornyhead") # never has data for the 3 main regions
# sanity summary for all models fit (not necessarily converged)
sanity_df <- readRDS(file.path(df_dir, "sanity-df.rds")) |>
  janitor::clean_names() |>
  filter(species != "longspine thornyhead") # never has data for the 3 main regions


gg_fit_ests <- readRDS(file.path(df_dir, "gengamma-all-fitted-estimates-df.rds")) |>
  janitor::clean_names() |>
  filter(species != "longspine thornyhead")
gg_Q <- gg_fit_ests |>
  distinct(species, region, fit_family, est_q, est_qse) |>
  mutate(pretty_q = ifelse(abs(est_q) > 1, round(est_q), round(est_q, 2)))


# Number of species:
length(unique(lu_df$species))
unique(lu_df$species)

# Number of models tried:
nrow(lu_df)

# Number of models that converged
# FIXME: NUMBERS LOOK INCORRECT
sanity_df |>
  filter(all_ok) |>
  tabyl(fit_family)

fit_ests |>
  filter(fit_family == "gengamma") |>
  select(species, region, est_q, est_qse) |>
  pull("est_q") |>
  hist()


# ------
aic_df <- fit_ests |>
  select(species, region, family, aic) |>
  group_by(species, region) |>
  mutate(
    min_aic = min(aic),
    daic = aic - min_aic
  ) |>
  ungroup() |>
  arrange(species)

aic_plot_df <- left_join(lu_df, aic_df) |>
  left_join(gg_Q) |>
  mutate(
    species_id = as.numeric(factor(species)),
    odd_species = ifelse(species_id %% 2 == 0, 1, 0),
    converged = ifelse(is.na(aic), 0, 1),
    x = case_when(
      is.na(daic) & family == "tweedie" ~ 10^-1,
      is.na(daic) & family == "delta-lognormal" ~ 10^-0.867,
      is.na(daic) & family == "delta-gamma" ~ 10^-0.733,
      is.na(daic) & family == "delta-gengamma" ~ 10^-0.6,
      TRUE ~ daic + 1
    )
  ) |>
  filter(region != "SYN WCHG") |>
  left_join(spp_lu) |>
  left_join(mean_pos_sets, by = c("species", "region")) |>
  mutate(fspecies = forcats::fct_rev(species)) |>
  mutate(fspecies2 = forcats::fct_reorder(fspecies, as.numeric(species_code))) |>
  mutate(odd_species2 = ifelse(as.numeric(fspecies2) %% 2 == 0, 1, 0))

aic_plot_df <- left_join(lu_df, aic_df) |>
  filter(!(species %in% c("butter sole", "copper rockfish", "shortraker rockfish"))) |> # only in one region
  left_join(gg_Q) |>
  mutate(
    species_id = as.numeric(factor(species)),
    odd_species = ifelse(species_id %% 2 == 0, 1, 0),
    converged = ifelse(is.na(aic), 0, 1),
    x = case_when(
      is.na(daic) & family == "tweedie" ~ 10^-1,
      is.na(daic) & family == "delta-lognormal" ~ 10^-0.867,
      is.na(daic) & family == "delta-gamma" ~ 10^-0.733,
      is.na(daic) & family == "delta-gengamma" ~ 10^-0.6,
      TRUE ~ daic + 1
    )
  ) |>
  filter(region != "SYN WCHG") |>
  left_join(spp_lu) |>
  left_join(mean_pos_sets, by = c("species", "region")) |>
  mutate(fspecies = forcats::fct_rev(species)) |>
  mutate(fspecies2 = forcats::fct_reorder(fspecies, as.numeric(species_code))) |>
  mutate(odd_species2 = ifelse(as.numeric(fspecies2) %% 2 == 0, 1, 0))

mean_sets <- aic_plot_df |>
  distinct(mean_sets, .keep_all = TRUE) |>
  mutate(mean_sets = round(mean_sets))

aic_plot <-
  ggplot() +
  geom_tile(
    data = aic_plot_df |> distinct(species, region, .keep_all = TRUE),
    mapping = aes(
      x = 1, y = fspecies2,
      width = Inf, height = 1, fill = factor(odd_species2)
    )
  ) +
  scale_fill_manual(values = c("grey95", "white")) +
  scale_x_continuous(
    trans = "log10", labels = scales::label_number(trim = TRUE),
    limits = c(0.04, 1100)
  ) +
  gfplot::theme_pbs(base_size = 13) +
  facet_wrap(~region, ncol = 3, drop = FALSE) +
  geom_point(
    data = aic_plot_df |> filter(converged == 1, prop_pos >= 0.05),
    aes(x = x, y = species, colour = family),
    stroke = 0.5, size = 2, position = ggstance::position_dodgev(height = 0.5)
  ) +
  geom_point(
    data = aic_plot_df |> filter(converged == 0, prop_pos >= 0.05),
    aes(x = x, y = species, colour = family),
    shape = 21, stroke = 0.7, size = 2
  ) +
  scale_colour_manual(values = family_colours) +
  labs(x = "Delta AIC + 1", y = "Species", colour = "Family") +
  guides(fill = "none") +
  geom_text(
    data = aic_plot_df |>
      distinct(species, region, .keep_all = TRUE) |>
      filter(prop_pos >= 0.05),
    aes(
      x = 0.08, y = fspecies2,
      label = round(mean_pos)
    ), # paste0(round(mean_pos), "(", prop_pos, ")")),
    hjust = 1, colour = "grey50",
    size = 3
  ) +
  geom_text(data = mean_sets, aes(
    x = 0.08, y = length(unique(aic_plot_df$species)),
    label = mean_sets
  ), vjust = -2, hjust = 1, size = 3) +
  # geom_text(data = aic_plot_df |> filter(prop_pos >= 0.05), aes(x = 0.04, y = fspecies2, label = signif(est_q, 2)),
  #           size = 2.5, hjust = 0) +
  coord_cartesian(clip = "off") +
  theme(legend.margin = margin(0, 0, 0, -0.3, "cm"))
aic_plot
ggsave(aic_plot, width = 13, height = 9, filename = file.path(fig_dir, "aic-plot-spp-code-order.png"))




# ------------------------------------------------------------------------------





# ggplot(data = _, aes(x = daic + 1, y = forcats::fct_rev(species))) +
#   geom_tile(data = filter(failed_df, species_id %% 2 == 0), aes(width = Inf, height = 1, fill = 'grey80'), colour = NA, alpha = 0.3) +
#   scale_fill_manual(values = c("grey", "white")) +
#   geom_jitter(aes(colour = family), width = 0, height = 0.5) +
#   geom_point(data = failed_df |> filter(is.na(daic), region != "SYN WCHG"),
#     aes(x = family_id, y = forcats::fct_rev(species), colour = family), shape = 21) +
#   scale_x_continuous(trans = "log10") +
#   facet_wrap(~ region, ncol = 3)
# # Question: do we ever fit delta models where spatial or spatiotemporal setting is different?
