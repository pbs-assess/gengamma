library(dplyr)
library(sdmTMB) # pak::pkg_install("pbs-assess/sdmTMB")
library(future)

source("R/00-utils.R")

# Make things faster by only fitting the 16 species that occur in all three surveys
spp <- c(
  "shortspine thornyhead",
  "rex sole",
  "longnose skate",
  # "silvergray rockfish", # not in the GoA data
  "arrowtooth flounder",
  "dover sole",
  "petrale sole",
  "english sole",
  "spotted ratfish",
  "flathead sole",
  "pacific ocean perch",
  "lingcod",
  "pacific cod",
  "sablefish",
  "walleye pollock",
  "north pacific spiny dogfish"
)

# Load data
# ---------------
if (file.exists(here::here("data", "clean-survey-data.rds")) &&
      file.exists(here::here("data", "clean-afsc-data.rds"))) {
  bc_dat <- readRDS(here::here("data", "clean-survey-data.rds")) |>
    filter(region != "SYN WCHG") |>
    filter(species %in% spp)

  goa_dat <- readRDS(here::here("data", "clean-afsc-data.rds")) |>
    mutate(longitude = lon_start, latitude = lat_start)

  dat <- bind_rows(
    bc_dat |> filter(prop_pos >= 0.2),
    goa_dat
  ) |>
    filter(!(region %in% c("SYN HS", "SYN QCS"))) |>
    select(-itis, -stratum, -lat_start, -lon_start, -effort)
  saveRDS(dat, here::here("data-outputs", "data-used.rds"))
} else {
  dat <- readRDS(here::here("data-outputs", "data-used.rds"))
}

spp_regions <- distinct(dat, species, region) |>
  rename(.species = "species", .region = "region") |>
  group_by(.species) |>
  mutate(n_regions = n()) |>
  ungroup()

# Setup files/dirs
# ---------------
tag <- ""

dc <- here::here("data", "raw")
out_dir <- here::here("data-outputs", "multi-species")
fit_dir <- here::here("data-outputs", "multi-species", paste0("fits", tag))
dir.create(fit_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Prep modelling inputs
# ----------------------
overwrite_cache <- FALSE
run_trouble_fits <- FALSE
.spatial <- "on"
.spatiotemporal <- "iid"
survey_dat <- dat
cutoff <- 8
regions <- c("SYN WCVI", "GOA", "HS-QCS")
families <- c(
  "tweedie",
  "delta-gamma", "delta-lognormal", "delta-gengamma"#,
  #"delta-gamma-poisson-link", "delta-lognormal-poisson-link", "delta-gengamma-poisson-link"
)
# "delta-plink-lognormal", "delta-plink-gamma", "delta-plink-gengamma")
family_lu <- c(
  "tweedie" = "tweedie",
  "delta-gamma" = "Gamma",
  "delta-lognormal" = "lognormal",
  "delta-gengamma" = "gengamma"#,
  # "delta-gamma-poisson-link" = "Gamma",
  # "delta-lognormal-poisson-link" = "lognormal",
  # "delta-gengamma-poisson-link" = "gengamma"
) %>%
  tibble(.family = names(x = .), fit_family = .) |>
  mutate(type = ifelse(grepl('poisson-link', .family), 'poisson-link', 'standard'))
tofit <- tidyr::expand_grid(
  # .region = regions,
  # .species = unique(spp_regions$.species),
  spp_regions |> select(-n_regions),
  .family = families
) |>
  mutate(.id = row_number())

lu_df <- tofit |>
  left_join(family_lu) |>
  left_join(spp_regions) |>
  mutate(sp_hyp = clean_name(.species)) |>
  tidyr::unite(col = "fname", sp_hyp, .region, .family, sep = "-", remove = FALSE)
saveRDS(lu_df, file.path(out_dir, "lu-df.rds"))

# ----------------------
# Fit models
# ----------------------
tofit <- tofit  |> filter(.region %in% regions) # only fit the big three

# overwrite_cache <- TRUE
future::plan(future::multisession, workers = 8)
progressr::with_progress({
  handler <- progressr::progressor(along = tofit$.id)
  tofit |>
    # purrr::pmap(\(.region, .family, .id, .species) {
    furrr::future_pmap(\(.region, .family, .id, .species) {
      handler() # Signal progress
      sp_file <- file.path(fit_dir, paste0(clean_name(.species), "-", .region, "-", .family, ".rds"))
      if (overwrite_cache | !file.exists(sp_file)) {
        message("Fitting: ", sp_file)
        f <- get_fit(
          survey_dat = survey_dat,
          formula = as.formula(catch_weight ~ 0 + as.factor(year)),
          offset = "offset",
          region = .region, family = .family, species = .species,
          cutoff = cutoff, sp = .spatial, st = .spatiotemporal,
          use_priors = FALSE
        )
        message("Saving fits: ", sp_file)
        saveRDS(f, sp_file)
      }
    })
  beep()
  future::plan(future::sequential())
})

# ------------------------------------------------------------------------------
# Look at trouble fits
# ------------------------------------------------------------------------------
if (run_trouble_fits) {
library(ggplot2)

trouble_fits <- tofit |>
  filter((.species == "walleye pollock" & .region == "GOA") |
        (.species == "petrale sole" & .region == "GOA") |
        (.species == "walleye pollock" & .region == "HS-QCS"))

prior_fit_dir <- file.path(out_dir, "priors", "fits")
dir.create(prior_fit_dir, showWarnings = FALSE, recursive = TRUE)

scale_fit_dir <- file.path(out_dir, "scale", "fits")
dir.create(scale_fit_dir, showWarnings = FALSE, recursive = TRUE)

hist_df <- semi_join(survey_dat, trouble_fits, by = c("species" = ".species", "region" = ".region")) |>
  filter(catch_weight > 0) |>
  mutate(species = stringr::str_to_title(species))

set_label <- hist_df |>
  group_by(species, region) |>
  mutate(bins = cut_interval(catch_weight, n = 30),
         max_catch = max(catch_weight)) |>
  add_count(bins) |>
  summarise(max_count = max(n), max_catch = max(max_catch)) |>
  left_join(hist_df |> distinct(species, region, mean_pos, mean_sets, mean_pos_sets, prop_pos)) |>
  mutate(percent_pos = round(prop_pos * 100)) |>
  # mutate(label = paste0("Mean~positive~sets:~",
  #                     "frac(", round(mean_pos), ",", round(mean_sets), ")",
  #                     "~`=`~", percent_pos, "*\"%\""))
  mutate(label = paste0("Mean positive sets: ", percent_pos, "%"))

ggplot(data = hist_df, mapping = aes(x = catch_weight)) +
  gfplot::theme_pbs(base_size = 12) +
  geom_histogram(bins = 30, closed = "right", boundary = 0.001) +
  geom_text(data = set_label, aes(x = max_catch, y = max_count,
    label = label),# parse = TRUE,
    hjust = 1, vjust = 3, size = 3.5) +
  scale_y_continuous(trans = "sqrt") +
  ggh4x::facet_nested_wrap(~ species + region, scale = "free") +
  theme(strip.text = element_text(size = 11)) +
  labs(x = "Catch weight (kg)", y = "Count")

ggsave(width = 7.8, height = 2.8, filename = here::here("figures", "supp", "petrale-pollock-dists.png"))

# Using priors
# -------------
trouble_fits |>
  purrr::pmap(\(.region, .family, .id, .species) {
    prior_sd = 50
    sp_file <- file.path(prior_fit_dir, paste0(clean_name(.species), "-", .region, "-", .family, ".rds"))
    message(clean_name(.species), "-", .region, "-", .family)
    f <- get_fit(
      survey_dat = survey_dat,
      formula = as.formula(catch_weight ~ 0 + as.factor(year)),
      offset = "offset",
      region = .region, family = .family, species = .species,
      cutoff = cutoff, sp = .spatial, st = .spatiotemporal,
      use_priors = TRUE, silent = FALSE
    )
    saveRDS(f, sp_file)
  })
beep()
# ---

# Scaling data
# -------------
trouble_fits |>
  filter(.species == "walleye pollock", .family == "delta-gengamma") |>
  purrr::pmap(\(.region, .family, .id, .species) {
    scale_factor <- 10
    sp_file <- file.path(scaled_fit_dir, paste0(clean_name(.species), "-", .region, "-", .family, '-', scale_factor, ".rds"))
    f <- get_fit(
      survey_dat = survey_dat,
      formula = as.formula(catch_weight / scale_factor ~ 0 + as.factor(year)),
      offset = "offset",
      region = .region, family = .family, species = .species,
      cutoff = cutoff, sp = .spatial, st = .spatiotemporal,
      use_priors = FALSE
    )
    saveRDS(f, sp_file)
  })
beep()
}