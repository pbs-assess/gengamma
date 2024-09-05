library(dplyr)
library(ggplot2)
library(janitor)
library(sdmTMB)
library(patchwork)
theme_set(gfplot::theme_pbs())

source("R/00-utils.R")

dir.create(here::here("figures", "multi-species"), showWarnings = FALSE, recursive = TRUE)
df_dir <- here::here("data-outputs", "multi-species")

family_levels <- c(
  "delta-lognormal", "delta-gamma", "tweedie", "delta-gengamma"#,
  #"delta-lognormal-poisson-link", "delta-gamma-poisson-link", "delta-gengamma-poisson-link"
)

family_shapes <- c(
  "tweedie" =  3, # plus
  "Tweedie" =  3, # plus
  "delta-lognormal" = 17, # triangle
  "delta-gamma" = 15, # square
  "delta-gengamma" = 19 # point
)


# https://mk.bcgsc.ca/colorblind/palettes/8.color.blindness.palette.txt
# family_colours <- c(
#   "tweedie" = "#F748A5", # 'barbie pink'
#   "delta-lognormal" = "#359B73", #' ocean green'
#   "delta-gamma" = "#d55e00", # 'bamboo'
#   "delta-gengamma" = "#2271B2" #' honolulu blue'
# )
# family_colours <- c(
#   "tweedie" = "#2271B2", #' honolulu blue'
#   "delta-lognormal" = "#009F81", #'jeepers creepers'
#   "delta-gamma" = "#d55e00", # 'bamboo'
#   "delta-gengamma" = "#9f0162" # 'jazzberry jam'
# )
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
  distinct(species, region, mean_pos, mean_sets, mean_pos_sets, prop_pos) |>
  mutate(mean_pos = round(mean_pos), mean_sets = round(mean_sets))

lu_df <- readRDS(file.path(df_dir, "lu-df.rds")) |>
  janitor::clean_names() |>
  left_join(spp_list) |>
  left_join(pos_sets)
lu_df <- lu_df |>
  filter(type != "poisson-link") |>
  filter(n_regions == 3)

# fit estimates of models that converged
fit_ests <- readRDS(file.path(df_dir, "fitted-estimates.rds")) |>
  janitor::clean_names() |>
  mutate(type = gsub("poisson_link_delta", "poisson-link", type)) |>
  mutate(type = ifelse(fit_family == "tweedie", "standard", type)) |>
  right_join(lu_df)

# sanity summary for all models fit (not necessarily converged)
sanity_df <- readRDS(file.path(df_dir, "sanity-df.rds")) |>
  janitor::clean_names() |>
  right_join(lu_df)

rqr_df <- readRDS(file.path(df_dir, "rqr-catch-df.rds")) |>
  right_join(lu_df) |>
  mutate(family = factor(family, levels = family_levels))

# indices of models that converged
index_df <- readRDS(file.path(df_dir, "index-df.rds")) |>
  janitor::clean_names() |>
  select(-type) |>
  right_join(lu_df) |>
  as_tibble()

# additional exclusion criterion is if index CV is > 1
index_sanity <- index_df |>
  group_by(fname, id) |>
  summarise(mean_ind_cv = mean(sqrt(exp(se^2) - 1))) |>
  filter(mean_ind_cv < 1)

summary_df <-
  left_join(fit_ests, sanity_df) |>
  right_join(index_sanity)

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
q_ests <- fit_ests |>
  filter(family == "delta-gengamma") |>
  filter(region == "GOA") |>
  #filter(region == "SYN WCVI") |>
  select(species, region, est_q, est_qse) |>
  arrange(est_q) |>
  mutate(q_rank = row_number())

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
  left_join(q_ests |> select(species, q_rank)) |>
  ungroup()

aic_plot_df <- aic_df |>
  mutate(species = gsub("north", "", species)) |>
  mutate(species = stringr::str_to_title(species)) |>
  mutate(
    species_id = as.numeric(factor(species)),
    odd_species = ifelse(species_id %% 2 == 0, 1, 0),
    converged = ifelse(all_ok, 1, 0),
    x = case_when(
      !all_ok ~ 10^-0.6,
      TRUE ~ daic + 1
    )
  ) |>
  mutate(fspecies = forcats::fct_rev(species)) |>
  mutate(fspecies2 = forcats::fct_reorder(fspecies, as.numeric(q_rank))) |>
  mutate(odd_species2 = ifelse(as.numeric(fspecies2) %% 2 == 0, 1, 0)) |>
  mutate(family = factor(family, levels = family_levels)) |>
  mutate(family = forcats::fct_recode(family, Tweedie = "tweedie")) |>
  mutate(fregion = factor(gsub("SYN ", "", region), levels = c("GOA", "HS", "QCS", "HS-QCS", "WCVI"))) |>
  filter(fregion %in% c("GOA", "HS-QCS", "WCVI"))

aic_w_plot <- ggplot() +
  gfplot::theme_pbs(base_size = 12) +
  facet_wrap(~fregion, ncol = 3, drop = TRUE) +
  geom_tile(
    data = aic_plot_df,
    mapping = aes(
      x = 1, y = fspecies2,
      width = Inf, height = 1, fill = factor(odd_species2)
    )
  ) +
  scale_fill_manual(values = c("grey95", "white")) +
  geom_point(
    data = aic_plot_df |> filter(all_ok), # models that converged
    aes(x = aic_w, y = species, colour = family, shape = family),
    stroke = 0.75, size = 2, position = ggstance::position_dodgev(height = 0.8)
  )  +
  scale_colour_manual(values = family_colours) +
  scale_shape_manual(values = family_shapes) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = scales::label_percent(suffix = ""),
    limits = c(0, 1.27)) +
  annotate(geom = "text", x = 1.25, y = length(unique(aic_plot_df$species)),
           vjust = -3, hjust = 0.75, label = "Q", colour = "grey30") +
  geom_text(data = filter(aic_plot_df, family == "delta-gengamma"),
    aes(x = 1.25, y = fspecies2, label = round(est_q, digits = 1)),
    hjust = 0.75, size = 3, colour = "grey50") +
  labs(x = "AIC weight (%)", y = "Species",
       shape = "Family",
       colour = "Family") +
  guides(fill = "none") +
  coord_cartesian(clip = "off") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(vjust = 0),
        axis.text = element_text(size = 9),
        panel.spacing.x = unit(0.2, "lines"),
        panel.border = element_rect(fill = NA, linewidth = 0.5),
        legend.title = element_blank(),
        legend.position = "top",
        legend.margin = margin(0, 0, -0.3, 0, "cm"))
aic_w_plot
ggsave(aic_w_plot, width = 7.5, height = 7.5, filename = file.path("figures", "figure-4-aic-weight.png"))

# Used this to look at how close Q is to sigma:
aic_df |>
select(species, region, fit_family, estimate_phi, est_q, est_qse, std_error_phi, estimate_range) |>
filter(#species %in% c("arrowtooth flounder", "shortespine thornyhead", "rex sole", "dover sole"),
      #fit_family %in% c("Gamma", "gengamma"),
      fit_family == "gengamma",
      region == "GOA"
      #region %in% c("HS-QCS", "SYN WCVI")
      ) |>
mutate(diff = estimate_phi - est_q) |>
arrange(abs(diff))
# when difference between sigma and estimated Q was smallest in WCVI and HS-QCS


# Plot indices
# -----------------
# index_region <- "SYN QCS"
# index_region <- "SYN HS"
# index_region <- "SYN WCVI"
index_region <- "HS-QCS"
index_region <- "GOA"
# index_region <- unique(summary_df$region)

# index_species <- c("silvergray rockfish")
#index_species <- c("rex sole")
index_species <- c("north pacific spiny dogfish", "sablefish", "shortspine thornyhead")
# index_species <- c("north pacific spiny dogfish", "lingcod", "pacific ocean perch", "arrowtooth flounder")
# index_species <- c("pacific ocean perch")
# index_species <- unique(summary_df$species)
# index_species <- index_species[!grepl('silvergray', index_species)]


index_q_rank <- summary_df |>
  filter(family == "delta-gengamma") |>
  filter(region %in% index_region) |>
  select(species, region, est_q, est_qse) |>
  arrange(est_q) |>
  mutate(q_rank = row_number()) |>
  mutate(fspecies = forcats::fct_inorder(species))

scaled_design <- design_df |>
  filter(species %in% index_species, region %in% index_region) |>
  group_by(species, region) |>
  mutate(design_geo_mean = exp(mean(log(est)))) |>
  distinct(species, region, design_geo_mean)


pind <- index_df |>
  right_join(index_sanity) |>
  left_join(select(aic_df, fname, estimate_phi, min_aic:aic_w)) |>
  # filter(family != "tweedie") |>
  filter(species %in% index_species, region %in% index_region) |>
  mutate(family = factor(family, levels = family_levels)) |>
  mutate(aic_w_text = paste0(round(aic_w * 100), "%")) |>
  group_by(species, region, family) |>
  mutate(est = est * log(1e5), upr = upr * log(1e5), lwr = lwr * log(1e5)) |>
  mutate(geo_mean = exp(mean(log(est)))) |>
  mutate(
    est_scaled = est,
    lwr_scaled = lwr,
    upr_scaled = upr
  ) |>
  ungroup() |>
  left_join(scaled_design) |>
  mutate(
    est_scaled = est_scaled / design_geo_mean,
    lwr_scaled = lwr_scaled / design_geo_mean,
    upr_scaled = upr_scaled / design_geo_mean
  ) |>
  select(-c(id, type:mean_pos_sets))

pind <- pind |>
  left_join(pind |>
    filter(!(family == "tweedie" & species == "north pacific spiny dogfish")) |>
    group_by(species) |>
    #summarise(ymax = max(upr_scaled))
    summarise(ymax = max(upr))
  ) |>
  left_join(pind |>
    distinct(species, region, family) |>
    arrange(species, region, family) |>
    group_by(species, region) |>
    # mutate(x = rep(c(2005, 2010, 2015, 2020), length.out = n()))
    mutate(x = rep(seq(1995, 2020, 8), length.out = n()))
  ) |>
  left_join(index_q_rank |> select(species, est_q, q_rank)) |>
  mutate(species = gsub("north ", "", species)) |>
  mutate(species = stringr::str_to_title(species)) |>
  arrange(desc(q_rank)) |>
  mutate(species = forcats::fct_inorder(species))

p1 <- ggplot(data = pind |> filter(!(family == "tweedie" & species == "Pacific Spiny Dogfish")),
    aes(x = year, y = est, group = family)) +
  geom_line(data = pind |> filter(!(family == "tweedie" & species == "Pacific Spiny Dogfish")), aes(colour = family)) +
  geom_line(data = pind |> filter((family == "tweedie" & species == "Pacific Spiny Dogfish")),
    aes(colour = family), alpha = 0.4) +
  geom_ribbon(data = pind |> filter(!(family == "tweedie" & species == "Pacific Spiny Dogfish")),
    aes(ymin = lwr, ymax = upr, fill = family), alpha = 0.1) +
  geom_ribbon(data = pind |> filter(family == "tweedie" & species == "Pacific Spiny Dogfish"),
    aes(ymin = lwr, ymax = upr, fill = family), alpha = 0.05) +
  geom_rect(data = pind |> filter((family == "tweedie" & species == "Pacific Spiny Dogfish")) |> distinct(),
    aes(xmin = 1990, xmax = 2023, ymin = 115e5, ymax = 140e5), fill = "white") +
  geom_text(data = distinct(pind, species, family, region, aic_w_text, ymax, .keep_all = TRUE),
            aes(x = x, y = ymax, label = aic_w_text, colour = family)) +
  theme(plot.margin = margin(c(0, 0, 0, 0)),
    strip.text.x.top = element_text(hjust = 0, size = 11, margin = margin(c(0.1, 0.1, 0.14, 0.1), unit = "cm")),
  ) +
  scale_colour_manual(values = family_colours) +
  scale_fill_manual(values = family_colours) +
  geom_point(data = tibble(year = 2005, est = -1, family = 'gengamma'), alpha = 0) + # weird hack needed to get 0 on y-axis
  scale_y_continuous(labels = scales::label_number(scale = 1 / 1e5),
                     expand = expansion(mult = c(0, 0.1), add = c(-1, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0), add = c(0.2, 0.2))) +
  labs(x = "Year", y = "Biomass index (x 100 000)") +
  # geom_pointrange(
  #   data = design_df |>
  #     filter(species %in% index_species, region %in% index_region) |>
  #     mutate(species = gsub("north ", "", species)) |>
  #     mutate(species = stringr::str_to_title(species)) |>
  #     mutate(species = factor(species, levels = levels(pind$species))),
  #   mapping = aes(x = year, y = est, ymin = lwr, ymax = upr),
  #   size = 0.5) +
  guides(colour = "none", fill = "none") +
  # facet_wrap(species ~ ., scales = "free_y", ncol = 1) +
  # coord_cartesian(clip = 'off') +
  facet_wrap_custom(species ~ ., scales = "free_y", ncol = 1,
    scale_overrides = list(scale_override(3,
      scale_y_continuous(limits = c(-1, 125 * 1e5),
      labels = scales::label_number(scale = 1 / 1e5),
      expand = expansion(mult = c(0, 0.1), add = c(-1, 0)),
      oob = scales::squish_infinite
      #na.value = 200 * 1e5 + (200 * 1e5 * 0.1)
      ))))
p1

p2 <- rqr_df |>
  filter(species %in% index_species, region %in% index_region) |>
  left_join(select(aic_df, species, family, fname, q_rank)) |>
  mutate(fit_family = factor(tolower(fit_family), gsub("delta-", "", family_levels))) |>
  mutate(fit_family = forcats::fct_recode(fit_family, Tweedie = "tweedie")) |>
  mutate(species = factor(species)) |>
  mutate(species = forcats::fct_reorder(species, order(q_rank, decreasing = TRUE))) |>
  ggplot(aes(sample = r)) +
  geom_qq() +
  geom_abline(intercept = 0, slope = 1) +
  facet_grid(species ~ fit_family, switch = "y") +
  #facet_wrap(species ~ fit_family, labeller = labeller(.multi_line = FALSE)) +
  coord_fixed(xlim = c(-4, 4), ylim = c(-5.3, 5.3)) +
  #lims(x = c(-6, 6), y = c(-6, 6)) +
  scale_x_continuous(breaks = seq(-6, 6, by = 2), labels = seq(-6, 6, by = 2)) +
  scale_y_continuous(breaks = seq(-6, 6, by = 2), labels = seq(-6, 6, by = 2), position = "right") +
  labs(x = "Theoretical", y = "Sample") +
  theme(plot.margin = margin(c(0.5, 0, 0.5, 0)),
        strip.text.y = element_blank(),
        strip.text.x = element_text(),
        panel.spacing.y = unit(1.5, "lines"),
        panel.spacing.x = unit(0, "lines")
  )
#p2
# see: https://github.com/tidyverse/ggplot2/issues/2096#issuecomment-389825118
# dev.new(width = 9.372340, height = 7.031579)
{
  g <- ggplot_gtable(ggplot_build(p2))
  stript <- which(grepl('strip-t', g$layout$name))
  fills <- rep(scales::alpha(family_colours[family_levels], alpha = 0.3), 3)
  k <- 1
  for (i in stript) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  wrap_elements(full = p1) + wrap_elements(full = g) + plot_layout(widths = c(0.46, 0.8))
}

ggsave(width = 9.180851, height = 6.631579, filename = file.path("figures", "figure-5-index-rqr.png"))


# All indices - species & regions
# --------------------------------
plot_all_spp_index <- function(.region, exclusion_spp, .ncol = 3, guide_rows = 1) {
  # Custom labeller function
  custom_labeller <- function(df) {
    # Create a named vector with species and corresponding est_q
    labels <- paste0(df$species, " | Q = ", round(df$est_q, digits = 2), "")
    names(labels) <- df$species
    return(labels)
  }

  index_q_rank <- summary_df |>
    filter(family == "delta-gengamma") |>
    filter(region %in% .region) |>
    select(species, region, est_q, est_qse) |>
    arrange(est_q) |>
    mutate(q_rank = row_number()) |>
    mutate(fspecies = forcats::fct_inorder(species))

  scaled_design <- design_df |>
    filter(region %in% .region) |>
    group_by(species, region) |>
    mutate(design_geo_mean = exp(mean(log(est)))) |>
    distinct(species, region, design_geo_mean)

  pind <- index_df |>
    left_join(index_sanity) |>
    left_join(select(aic_df, fname, estimate_phi, min_aic:aic_w)) |>
    filter(region %in% .region) |>
    mutate(family = factor(family, levels = family_levels)) |>
    mutate(aic_w_text = ifelse(is.na(aic_w), "\u2013 ", paste0(round(aic_w * 100), "%"))) |>
    group_by(species, region, family) |>
    mutate(est = est * log(1e5), upr = upr * log(1e5), lwr = lwr * log(1e5)) |>
    mutate(geo_mean = exp(mean(log(est)))) |>
    mutate(
      est_scaled = est,
      lwr_scaled = lwr,
      upr_scaled = upr
    ) |>
    ungroup() |>
    left_join(scaled_design) |>
    mutate(
      est_scaled = est_scaled / design_geo_mean,
      lwr_scaled = lwr_scaled / design_geo_mean,
      upr_scaled = upr_scaled / design_geo_mean
    ) |>
    select(-c(id, type:mean_pos_sets))
  pind <- pind |>
    left_join(pind |>
      filter(!(family == "tweedie" & aic_w < 0.1)) |>
      group_by(species, region) |>
      summarise(ymax = max(upr, na.rm = TRUE))
    ) |>
    left_join(pind |>
      distinct(species, region, family) |>
      arrange(species, region, family) |>
      group_by(species, region) |>
      mutate(x = ifelse(region == "GOA",
        rep(seq(1995, 2020, 8), length.out = n()),
        rep(seq(2007, 2022, 4), length.out = n()))
      )
    ) |>
    left_join(index_q_rank |> select(species, est_q, q_rank)) |>
    mutate(species = gsub("north ", "", species)) |>
    mutate(species = stringr::str_to_title(species)) |>
    arrange(desc(q_rank)) |>
    mutate(species = forcats::fct_inorder(species))

  pind <- pind |>
    filter(region == .region) |>
    filter(!(species %in% exclusion_spp)) |>
    mutate(species = forcats::fct_reorder(species, rev(est_q))) |>
    mutate(family = forcats::fct_recode(family, "Tweedie" = "tweedie"))

  ggplot(data = pind, aes(x = year, y = est, group = family, colour = family, fill = family)) +
    geom_line(data = pind |> filter(!(family == "Tweedie" & aic_w < 0.1)), aes(colour = family)) +
    geom_ribbon(data = pind |> filter(!(family == "Tweedie" & aic_w < 0.1)),
      aes(ymin = lwr, ymax = upr, fill = family), colour = NA, alpha = 0.1) +
    geom_line(data = pind |> filter(family == "Tweedie") |> distinct(species, family, x, ymax),
      aes(x = x, y = ymax), alpha = 0) +
    geom_ribbon(data = pind |> filter(family == "Tweedie") |> distinct(species, family, x, ymax),
      aes(x = x, y = ymax, ymin = ymax, ymax = ymax), colour = NA, alpha = 0) +
    # geom_line(data = pind, aes(colour = family)) +
    # geom_ribbon(data = pind,
    #   aes(ymin = lwr, ymax = upr, fill = family), alpha = 0.1) +
    geom_text(data = distinct(pind, species, family, region, aic_w_text, ymax, .keep_all = TRUE),
              aes(x = x, y = ymax, label = aic_w_text, colour = family
              ), show.legend = FALSE) +
    theme(
      strip.text.x.top = element_text(hjust = 0, size = 11, margin = margin(c(0.1, 0.1, 0.14, 0.1), unit = "cm")),
    ) +
    scale_colour_manual(values = family_colours) +
    scale_fill_manual(values = family_colours) +
    geom_point(data = tibble(year = 2003, est = -1, family = 'gengamma'), alpha = 0) + # weird hack needed to get 0 on y-axis
    scale_y_continuous(labels = scales::label_number(scale = 1 / 1e5),
                      expand = expansion(mult = c(0, 0.1), add = c(-1, 0))) +
    scale_x_continuous(expand = expansion(mult = c(0, 0), add = c(0.2, 0.2))) +
    labs(x = "Year", y = "Biomass index (x 100 000)") +
    guides(colour = guide_legend(override.aes = list(alpha = 1, nrow = guide_rows))
      ) +
    facet_wrap(species ~ ., labeller = labeller(species = custom_labeller(pind)), scales = "free_y", ncol = .ncol) +
    ggtitle(gsub("SYN ", "", .region)) +
    theme(plot.title = element_text(size = 12, face = "bold"),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(size = 10))
  }

plot_all_spp_index(
  exclusion_spp =  c("Petrale Sole", "Walleye Pollock"),
  .region = "GOA",
  .ncol = 3
)
ggsave(width = 8.5, height = 9.3, filename = here::here("figures", "supp", "index-goa.png"))

plot_all_spp_index(
  exclusion_spp =  c("Walleye Pollock"),
  .region = "HS-QCS",
  .ncol = 3
) +
  theme(legend.position = c(0.66, 0.08),
        legend.direction = "horizontal")
ggsave(width = 8.5, height = 9.3, filename = here::here("figures", "supp", "index-hs-qcs.png"))

plot_all_spp_index(
  exclusion_spp = NULL,
  .region = "SYN WCVI",
  .ncol = 3,
  guide_rows = 2
) +
theme(legend.position = c(0.66, 0.08))
ggsave(width = 8.5, height = 9.3, filename = here::here("figures", "supp", "index-wcvi.png"))


# Compare scale of design and model
# -----------------------------------
# Could scale them all to the design based index - and then you would see bias based how far
# off from 1
p3 <- ggplot(data = pind, aes(x = year, y = est, group = family)) +
  geom_line(aes(colour = family)) +
  geom_pointrange(data = design_df |> filter(species %in% index_species, region %in% index_region) |>
    mutate(species = stringr::str_to_title(species)),
    mapping = aes(x = year, ymin = lwr / 1e5, ymax = upr)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = family), alpha = 0.1) +
  geom_text(data = distinct(pind, species, family, region, aic_w_text, ymax, .keep_all = TRUE),
    aes(x = x, y = ymax, label = aic_w_text, colour = family)) +
  scale_y_continuous() +
  facet_grid(species ~ region, scales = "free_y") +
  scale_colour_manual(values = family_colours) +
  scale_fill_manual(values = family_colours) +
  theme(legend.position = c(2010, 1e5))

p4 <- rqr_df |>
  filter(species %in% index_species, region %in% index_region) |>
  left_join(select(aic_df, species, family, fname, q_rank)) |>
  mutate(fit_family = factor(tolower(fit_family), gsub("delta-", "", family_levels))) |>
  mutate(fit_family = forcats::fct_recode(fit_family, Tweedie = "tweedie")) |>
  mutate(species = factor(species)) |>
  mutate(species = forcats::fct_reorder(species, order(q_rank, decreasing = TRUE))) |>
  ggplot(aes(sample = r)) +
  geom_qq() +
  geom_abline(intercept = 0, slope = 1) +
  facet_grid(species ~ fit_family, switch = "y") +
  #facet_wrap(species ~ fit_family, labeller = labeller(.multi_line = FALSE)) +
  coord_fixed(xlim = c(-4, 4), ylim = c(-5.3, 5.3)) +
  #lims(x = c(-6, 6), y = c(-6, 6)) +
  scale_x_continuous(breaks = seq(-6, 6, by = 2), labels = seq(-6, 6, by = 2)) +
  scale_y_continuous(breaks = seq(-6, 6, by = 2), labels = seq(-6, 6, by = 2), position = "right",) +
  labs(x = "Theoretical", y = "Sample") +
  theme(plot.margin = margin(c(0.5, 0, 0.5, 0.5)),
        strip.text.y = element_blank(),
        strip.text.x = element_text(),
        panel.spacing.y = unit(1.6, "lines"),
        panel.spacing.x = unit(0, "lines")
  )

# see: https://github.com/tidyverse/ggplot2/issues/2096#issuecomment-389825118
g <- ggplot_gtable(ggplot_build(p4))
stript <- which(grepl('strip-t', g$layout$name))
fills <- rep(scales::alpha(family_colours[family_levels], alpha = 0.3), 3)
k <- 1
for (i in stript) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
p3 + wrap_elements(full = g) + plot_layout(widths = c(0.4, 0.8))


# Get ratio of model-based index / design-based index - like surprising scale paper
# ------------------------------
# doesn't use the geometric mean?
dm_ratio <- left_join(
  design_df |>
    group_by(species, region) |>
    mutate(B = mean(est) / n()) |>
    distinct(species, region, B),
  index_df |>
    group_by(species, fit_family, region) |>
    mutate(I = mean(est * log(1e5)) / n()) |>
    distinct(species, region, I)
) |>
  mutate(ratio = I / B)

dm_ratio |> filter(species %in% index_species, region %in% index_region)

ggplot(data = dm_ratio, aes(x = log(ratio))) +
  geom_histogram() +
  geom_vline(xintercept = 0) +
  facet_grid(. ~ fit_family)



# ------------------------------------------------------------------------------
# Will anyone want to see the estimated CV of the different fits to real data?
get_ln_phi <- function(model) {
  m <- ifelse(family(model)[[1]][[1]] == 'tweedie', 1, 2)
  tibble(
    ln_phi = model$parlist$ln_phi[[m]],
    species = unique(model$data$species),
    region = unique(model$data$region),
    spatial = unique(model$spatial),
    spatiotemporal = unique(model$spatiotemporal),
    fit_family = family(model)[[m]][[1]],
    type = family(model)$type
  )
}

get_gamma_cv <- function(phi) 1 / sqrt(phi)

f <- list.files(file.path(df_dir, "fits"), full.names = TRUE)
f <- f[grepl(".rds", f)]

ln_phi_df <- f |>
  purrr::map(~ {
    message(stringr::str_extract(.x, ".*\\/(.*).rds$", group = 1))
    model <- readRDS(.x)
    try(get_ln_phi(model))
  }) |>
  purrr::keep(is.data.frame) |>
  bind_rows()

filter(ln_phi_df, fit_family == "Gamma", type == "standard") |>
  filter(!(region %in% c("SYN HS", "SYN QCS"))) |>
  #right_join(lu_df |> filter(fit_family == "Gamma", type == "standard")) |>
  mutate(cv = get_gamma_cv(phi = exp(ln_phi))) |>
  pull('cv') |>
  # hist()
  #modeest::mlv()
  mean()

filter(ln_phi_df, fit_family == "Gamma", type == "standard", region == "SYN WCVI") |>
  filter(species %in% c('pacific cod', 'longnose skate', 'arrowtooth flounder', 'dover sole', 'shortspine thornyhead')) |>
  mutate(cv = get_gamma_cv(phi = exp(ln_phi))) |>
  left_join(q_ests)




# ggplot(data = _, aes(x = daic + 1, y = forcats::fct_rev(species))) +
#   geom_tile(data = filter(failed_df, species_id %% 2 == 0), aes(width = Inf, height = 1, fill = 'grey80'), colour = NA, alpha = 0.3) +
#   scale_fill_manual(values = c("grey", "white")) +
#   geom_jitter(aes(colour = family), width = 0, height = 0.5) +
#   geom_point(data = failed_df |> filter(is.na(daic), region != "SYN WCHG"),
#     aes(x = family_id, y = forcats::fct_rev(species), colour = family), shape = 21) +
#   scale_x_continuous(trans = "log10") +
#   facet_wrap(~ region, ncol = 3)
# # Question: do we ever fit delta models where spatial or spatiotemporal setting is different?

# Compare aic results of the models where we used a prior
tr_df_dir <- file.path("data-outputs/multi-species/", 'priors')
tr_aic_df <- readRDS(file.path(tr_df_dir, 'aic-df.rds'))
tr_q <- tr_aic_df |>
  filter(family == "delta-gengamma") |>
  arrange(est_q) |>
  mutate(q_rank = row_number()) |>
  select(species, region, q = est_q, q_rank)

test <- aic_df |>
  select(-q_rank) |>
  mutate(priors = FALSE) |>
  filter((fname %in% tr_aic_df$fname)) |>
  bind_rows(tr_aic_df) |>
  mutate(species2 = ifelse(priors, paste0(species, "*"), species)) |>
  mutate(
    species_id = as.numeric(factor(species)),
    odd_species = ifelse(species_id %% 2 == 0, 1, 0)
    ) |>
  left_join(tr_q) |>
  mutate(species = paste0(species, "\n(", gsub("SYN ", "", region), ")")) |>
  mutate(fspecies = forcats::fct_rev(species)) |>
  mutate(fspecies2 = forcats::fct_reorder(fspecies, as.numeric(q_rank))) |>
  mutate(odd_species2 = ifelse(as.numeric(fspecies2) %% 2 == 0, 1, 0)) |>
  mutate(family = factor(family, levels = family_levels))

test |>
  mutate(priors = ifelse(priors, "Prior", "No prior")) |>
ggplot() +
  gfplot::theme_pbs(base_size = 12) +
  facet_wrap(~priors, ncol = 3, drop = FALSE) +
  geom_tile(mapping = aes(x = 1, y = fspecies2,
    width = Inf, height = 1, fill = factor(odd_species2)
  )) +
  scale_fill_manual(values = c("grey95", "white")) +
  geom_point(
    aes(x = aic_w, y = species, colour = family, shape = family),
    stroke = 0.5, size = 2, position = ggstance::position_dodgev(height = 1)
  )  +
  guides(fill = "none", shape = "none", colour = "none") +
  scale_colour_manual(values = family_colours) +
  scale_shape_manual(values = family_shapes) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = scales::label_percent(suffix = ""),
    limits = c(0, 1.27)) +
  annotate(geom = "text", x = 1.25, y = length(unique(test$species)),
        vjust = -2, hjust = 0.75, label = "Q", colour = "grey30") +
  geom_text(aes(x = 1.25, y = fspecies2, label = round(q, digits = 1)),
    hjust = 0.75, size = 3, colour = "grey50") +
  coord_cartesian(clip = "off")

tr_ind <- readRDS(file.path(tr_df_dir, 'index-df.rds')) |>
  filter(!stringr::str_detect(fname, 'poisson-link')) |>
  select(-type) |>
  left_join(lu_df) |>
  as_tibble() |>
  left_join(tr_q)

tr_ind_cv <- tr_ind |>
  group_by(fname) |>
  summarise(mean_ind_cv = mean(sqrt(exp(se^2) - 1))) |>
  left_join(lu_df) |>
  select(-c(fname, id, fit_family, n_regions, sp_hyp:mean_pos_sets)) |>
  mutate(f_family = as.numeric(factor(family))* 1.2 * 1e6,
         species = factor(species),
         label = paste0(gsub("delta-", "", family), ": ", round(mean_ind_cv, digits = 2))) |>
  left_join(tr_q) |>
  mutate(species = forcats::fct_reorder(species, as.numeric(rev(q_rank))))

tr_ind |>
  #filter(species == "silvergray rockfish") |>
  filter(species %in% unique(tr_ind$species)) |>
  ggplot(aes(x = year, y = est * log(1e5), group = family)) +
    # geom_ribbon(data = sg_i, aes(x = year, ymin = lwr * 10, ymax = upr * 10, group = NULL)) +
    # geom_line(data = sg_i, aes(x = year, y = est * 10, group = NULL), colour = 'red') +
    geom_line(aes(colour = family)) +
    geom_ribbon(aes(ymin = lwr * log(1e5), ymax = upr * log(1e5), fill = family), alpha = 0.1) +
    scale_colour_manual(values = family_colours) +
    scale_fill_manual(values = family_colours) +
    facet_wrap(forcats::fct_reorder(species, as.numeric(rev(q_rank))) ~ ., nrow = 3, scales = "free_y") +
    geom_text(data = (tr_ind_cv), aes(x = 2025, y = f_family, label = label), hjust = 1) +
    coord_cartesian(clip = "off") +
    geom_pointrange(
      data = design_df |>
        filter(species == 'silvergray rockfish', region == "SYN WCVI"),
      mapping = aes(x = year, y = est, ymin = lwr, ymax = upr))


ggplot(data = pind, ) +
  geom_line(aes(colour = family)) +
  geom_pointrange(data = design_df |> filter(species %in% index_species, region %in% index_region),
    mapping = aes(x = year, ymin = lwr / 1e5, ymax = upr)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = family), alpha = 0.1) +
  geom_text(data = distinct(pind, species, family, region, aic_w_text, ymax, .keep_all = TRUE),
    aes(x = x, y = ymax, label = aic_w_text, colour = family)) +
  scale_y_continuous() +
  facet_grid(species ~ region, scales = "free_y") +
  scale_colour_manual(values = family_colours) +
  scale_fill_manual(values = family_colours) +
  theme(legend.position = c(2010, 1e5))
