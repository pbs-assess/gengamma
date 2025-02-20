library(dplyr)
library(ggplot2)
library(janitor)
library(sdmTMB)
library(patchwork)
theme_set(gfplot::theme_pbs())

source("R/00-utils.R")

df_dir <- here::here("data-outputs", "multi-species")
df_dir_priors <- here::here(df_dir, "priors")
df_dir_scaled <- here::here(df_dir, "scaled")

family_levels <- c("delta-lognormal", "delta-gamma", "tweedie", "delta-gengamma")

family_shapes <- c(
  "tweedie" =  3, # plus
  "Tweedie" =  3, # plus
  "delta-lognormal" = 17, # triangle
  "lognormal" = 17, # triangle
  "delta-gamma" = 15, # square
  "gamma" = 15, # square
  "delta-gengamma" = 19, # point
  "gengamma" = 19 # point
)

lu_df <- readRDS(file.path(here::here(df_dir), "lu-df.rds")) |>
  janitor::clean_names()

# fit estimates of models that converged
fit_ests <- bind_rows(
  readRDS(file.path(df_dir_priors, "fitted-estimates.rds")),
  readRDS(file.path(df_dir_scaled, "fitted-estimates.rds")),
  readRDS(file.path(df_dir, "fitted-estimates.rds")) |>
    filter((species == "walleye pollock" & fit_family != "gengamma" & region == "GOA") |
          (species == "petrale sole" & fit_family != "gengamma" & region == "GOA")) |>
    mutate(scale_factor = 1)
) |>
  janitor::clean_names() |>
  left_join(lu_df)

# sanity summary for all models fit (not necessarily converged)
sanity_df <- bind_rows(
  readRDS(file.path(df_dir_priors, "sanity-df.rds")),
  readRDS(file.path(df_dir_scaled, "sanity-df.rds")),
  readRDS(file.path(df_dir, "sanity-df.rds")) |>
    filter((species == "walleye pollock" & fit_family != "gengamma" & region == "GOA") |
          (species == "petrale sole" & fit_family != "gengamma" & region == "GOA")) |>
    mutate(scale_factor = 1)
) |>
  janitor::clean_names() |>
  left_join(lu_df)

sanity_df |>
  distinct(species, region, fit_family, priors, scale_factor, all_ok) |>
  filter(all_ok == FALSE)

rqr_df <- bind_rows(
  readRDS(file.path(df_dir_priors, "rqr-catch-df.rds")),
  readRDS(file.path(df_dir_scaled, "rqr-catch-df.rds")),
  readRDS(file.path(df_dir, "rqr-catch-df.rds")) |>
    filter((stringr::str_detect(fname, "walleye-pollock-(GOA)-") &
         !stringr::str_detect(fname, "delta-gengamma")) |
         (stringr::str_detect(fname, "petrale-sole-(GOA)-") &
         !stringr::str_detect(fname, "delta-gengamma"))) |>
    mutate(scale_factor = 1)
) |>
  left_join(lu_df) |>
  mutate(family = factor(family, levels = family_levels)) |>
  mutate(family = forcats::fct_recode(family, "Tweedie" = "tweedie")) |>
  mutate(fix_type = case_when(priors == TRUE ~ "Penalty on year effect",
                             is.numeric(scale_factor) ~ "Scaled catch"))

# indices of models that converged
index_df <- bind_rows(
  readRDS(file.path(df_dir_priors, "index-df.rds")),
  readRDS(file.path(df_dir_scaled, "index-df.rds")),
  readRDS(file.path(df_dir, "index-df.rds")) |>
    filter((stringr::str_detect(fname, "walleye-pollock-(GOA)-") &
         !stringr::str_detect(fname, "delta-gengamma")) |
         (stringr::str_detect(fname, "petrale-sole-(GOA)-") &
         !stringr::str_detect(fname, "delta-gengamma"))) |>
    mutate(scale_factor = 1),
) |>
  janitor::clean_names() |>
  select(-type) |>
  left_join(lu_df) |>
  as_tibble() |>
  group_by(fname, id) |>
  mutate(mean_ind_cv = mean(sqrt(exp(se^2) - 1))) |>
  filter(mean_ind_cv < 1)

summary_df <-
  left_join(fit_ests, sanity_df) |>
  mutate(fix_type = case_when(priors == TRUE ~ "Penalty on year effect",
                             is.numeric(scale_factor) ~ "Scaled catch")) |>
  group_by(species, region, fix_type) |>
  mutate(
    min_aic = min(aic),
    daic = aic - min_aic,
    rel_lik = exp(-0.5 * daic)
  ) |>
  mutate(aic_w = rel_lik / sum(rel_lik)) |>
  ungroup() |>
  right_join(index_df) |>
  mutate(family = factor(family, levels = rev(family_levels))) |>
  mutate(family = forcats::fct_recode(family, "Tweedie" = "tweedie",
    "lognormal" = "delta-lognormal", "gengamma" = "delta-gengamma", "gamma" = "delta-gamma"))

aic_w <- ggplot(data = summary_df |> filter(species == "walleye pollock")) +
  gfplot::theme_pbs(base_size = 12) +
  geom_point(
    aes(x = aic_w, y = family, colour = family, shape = family),
    stroke = 0.75, size = 2, position = ggstance::position_dodgev(height = 0.8)
  ) +
  scale_colour_manual(values = family_colours_no_delta) +
  scale_shape_manual(values = family_shapes) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = scales::label_percent(suffix = "")) +
  guides(colour = "none", shape = "none") +
  labs(x = "AIC weight (%)", y = "Fit family") +
  ggh4x::facet_nested_wrap(~ fix_type + region, ncol = 3, drop = TRUE)
aic_w
ggsave(width = 5.3, height = 1.9, filename = here::here("figures", "supp", "pollock-prior-scaled-aic-w.png"))

# Index plot
index_species <- c("walleye pollock")

pind <- summary_df |>
  select(species, region, fit_family, year, est, lwr, upr, est_q, aic_w, fname, family, mean_ind_cv, fix_type, priors, scale_factor) |>
  mutate(est = ifelse(fix_type == "Scaled catch", est * scale_factor, est),
         lwr = ifelse(fix_type == "Scaled catch", lwr * scale_factor, lwr),
         upr = ifelse(fix_type == "Scaled catch", upr * scale_factor, upr)) |>
  #mutate(family = factor(family, levels = family_levels)) |>
  bind_rows()
  mutate(aic_w_text = ifelse(is.na(aic_w), "\u2013 ", paste0(round(aic_w * 100), "%"))) |>
  group_by(species, region, family, fix_type) |>
  mutate(est = est * log(1e5), upr = upr * log(1e5), lwr = lwr * log(1e5)) |>
  ungroup() |>
  group_by(region) |>
  mutate(ymax = max(upr, na.rm = TRUE),
         min_year = ifelse(region == "HS-QCS", min(year) + 2, min(year) + 3),
         max_year = ifelse(region == "HS-QCS", max(year) - 2, max(year) - 3)) |>
  rowwise() |>
  mutate(x = case_when(fit_family == "lognormal" ~ seq(min_year, max_year, length.out = 4)[1],
                       fit_family == "Gamma" ~ seq(min_year, max_year, length.out = 4)[2],
                       fit_family == "tweedie" ~ seq(min_year, max_year, length.out = 4)[3],
                       fit_family == "gengamma" ~ seq(min_year, max_year, length.out = 4)[4]
                       )) |>
  left_join(summary_df |>
    filter(fit_family == "gengamma") |>
    mutate(q_text = round(est_q, digits = 2)) |>
    distinct(region, fix_type, q_text)) |>
  mutate(strip_label = paste0(region, ": ", fix_type, " | Q = ", q_text))

pind1 <- pind |> filter(species %in% index_species)

p1 <- ggplot(data = pind1, aes(x = year, y = est, group = family)) +
  geom_line(aes(colour = family)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = family), alpha = 0.1) +
  geom_text(data = distinct(pind1, species, family, region, aic_w_text, ymax, strip_label, .keep_all = TRUE),
            aes(x = x, y = ymax, label = aic_w_text, colour = family)) +
  theme(plot.margin = margin(c(0, 0, 0, 0)),
        strip.text.x.top = element_text(hjust = 0, size = 11, margin = margin(c(0.1, 0.1, 0.14, 0.1), unit = "cm")),
  ) +
  scale_colour_manual(values = family_colours_no_delta) +
  scale_fill_manual(values = family_colours_no_delta) +
  geom_point(data = tibble(year = 2005, est = -1, family = 'gengamma'), alpha = 0) + # weird hack needed to get 0 on y-axis
  scale_y_continuous(labels = scales::label_number(scale = 1 / 1e5),
                     expand = expansion(mult = c(0, 0.1), add = c(-1, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0), add = c(0.2, 0.2))) +
  labs(x = "Year", y = "Biomass index (x 10<sup>5</sup>)") +
  guides(colour = "none", fill = "none") +
  facet_wrap(strip_label ~ ., scales = "free", ncol = 1) +
  theme(plot.background = element_blank(),
        axis.title.y = ggtext::element_markdown())
p1

p2 <- rqr_df |>
  filter(species %in% index_species) |>
  left_join(pind1 |> distinct(fname, species, region, priors, scale_factor, fix_type, strip_label)) |>
  mutate(fit_family = factor(tolower(fit_family), gsub("delta-", "", family_levels))) |>
  mutate(fit_family = forcats::fct_recode(fit_family, Tweedie = "tweedie")) |>
  ggplot(aes(sample = r)) +
  geom_qq() +
  geom_abline(intercept = 0, slope = 1) +
  facet_grid(strip_label ~ fit_family, switch = "y") +
  #facet_wrap(species ~ fit_family, labeller = labeller(.multi_line = FALSE)) +
  coord_fixed(xlim = c(-4, 4), ylim = c(-5.3, 5.3)) +
  #lims(x = c(-6, 6), y = c(-6, 6)) +
  scale_x_continuous(breaks = seq(-6, 6, by = 2), labels = seq(-6, 6, by = 2)) +
  scale_y_continuous(breaks = seq(-6, 6, by = 2), labels = seq(-6, 6, by = 2), position = "right") +
  labs(x = "Theoretical", y = "Sample") +
  theme(plot.margin = margin(c(0.5, 0, 0.5, 0)),
        strip.text.y = element_blank(),
        strip.text.x = element_text(),
        panel.spacing.y = unit(2.41, "lines"),
        panel.spacing.x = unit(0, "lines")
  )
# p2

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

}

d <- "
ADE
BDE
CDE"
wrap_elements(plot = grid::textGrob("a)", rot = 0, vjust = -8.25, hjust = 0), clip = FALSE) +
wrap_elements(plot = grid::textGrob("b)", rot = 0, vjust = -8.5, hjust = 0), clip = FALSE) +
wrap_elements(plot = grid::textGrob("c)", rot = 0, vjust = -8.8, hjust = 0), clip = FALSE) +
wrap_elements(full = p1) + wrap_elements(full = g) + plot_layout(widths = c(0, 0.45, 0.8), design = d)

ggsave(width = 9.372340, height = 7.031579, filename = file.path("figures", "supp", "pollock-index-rqr.png"))

# ---
# look at GOA petrale sole residuals
# ---
index_species <- c("petrale sole")

pind2 <- pind |>
  filter(species %in% index_species) |>
  filter(scale_factor == 1) |>
  group_by(species, region) |>
  mutate(ymax = max(upr, na.rm = TRUE),
         strip_label = paste0(stringr::str_to_title(species), " - ", region)) |>
  ungroup()

p1 <- ggplot(data = pind2, aes(x = year, y = est, group = family)) +
  geom_line(aes(colour = family)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = family), alpha = 0.1) +
  geom_text(data = distinct(pind2, species, family, region, aic_w_text, ymax, strip_label, .keep_all = TRUE),
            aes(x = x, y = ymax, label = aic_w_text, colour = family)) +
  theme(plot.margin = margin(c(0, 0, 0, 0)),
        strip.text.x.top = element_text(hjust = 0, size = 11, margin = margin(c(0.1, 0.1, 0.14, 0.1), unit = "cm")),
  ) +
  scale_colour_manual(values = family_colours_no_delta) +
  scale_fill_manual(values = family_colours_no_delta) +
  geom_point(data = tibble(year = 2005, est = -1, family = 'gengamma'), alpha = 0) + # weird hack needed to get 0 on y-axis
  scale_y_continuous(labels = scales::label_number(scale = 1 / 1e5),
                     expand = expansion(mult = c(0, 0.1), add = c(-1, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0), add = c(0.2, 0.2))) +
  labs(x = "Year", y = "Biomass index (x 10<sup>5</sup>)") +
  guides(colour = "none", fill = "none") +
  facet_wrap(strip_label ~ ., scales = "free_y", ncol = 1) +
  theme(plot.background = element_blank(),
        axis.title.y = ggtext::element_markdown())
p1

p2 <- rqr_df |>
  filter(species %in% index_species) |>
  filter(scale_factor == 1) |>
  left_join(pind2 |> distinct(fname, species, region, priors, scale_factor, fix_type, strip_label)) |>
  mutate(fit_family = factor(tolower(fit_family), gsub("delta-", "", family_levels))) |>
  mutate(fit_family = forcats::fct_recode(fit_family, Tweedie = "tweedie")) |>
  ggplot(aes(sample = r)) +
  geom_qq() +
  geom_abline(intercept = 0, slope = 1) +
  facet_grid(species ~ fit_family, switch = "y") +
  coord_fixed(xlim = c(-4, 4), ylim = c(-5.3, 5.3)) +
  scale_x_continuous(breaks = seq(-6, 6, by = 2), labels = seq(-6, 6, by = 2)) +
  scale_y_continuous(breaks = seq(-6, 6, by = 2), labels = seq(-6, 6, by = 2), position = "right") +
  labs(x = "Theoretical", y = "Sample") +
  theme(plot.margin = margin(c(0.5, 0, 0.5, 0)),
        strip.text.y = element_blank(),
        strip.text.x = element_text(),
        panel.spacing.y = unit(1.4, "lines"),
        panel.spacing.x = unit(0, "lines")
  )
p2

# dev.new(width = 9.372340, height = 3.063158)
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

}

wrap_elements(full = p1) + wrap_elements(full = g) + plot_layout(widths = c(0.45, 0.8))
ggsave(width = 9.372340, height = 3.063158, filename = file.path("figures", "supp", "petrale-index-rqr.png"))

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
  scale_colour_manual(values = family_colours_no_delta) +
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
  mutate(priors = ifelse(priors, "Penalty", "No prior")) |>
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
