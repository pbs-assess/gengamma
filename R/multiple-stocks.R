library(dplyr)
#library(sdmTMB)
devtools::load_all('../sdmTMB')
library(ggplot2)
theme_set(gfplot::theme_pbs())

source("R/00-utils.R")

# Load data
dat <- readRDS(file.path("data-outputs", "clean-survey-data.rds"))
spp_list <- unique(dat$species)

# Simplified synoptic grid
syn_grid <- as_tibble(gfplot::synoptic_grid) |>
  select(survey, X, Y, depth)

dc <- here::here("data", "raw")
out_dir <- here::here("data-outputs", "bcgf-outputs")
fit_dir <- here::here("data-outputs", "bcgf-outputs", "fits")
dir.create(fit_dir, showWarnings = FALSE, recursive = TRUE)

get_fit <- function(survey_dat, formula, region, family, species = NULL, cutoff = 8, time = "year",
                   sp = "on", st = "iid", offset = "offset") {
  survey_dat <- filter(survey_dat, survey_abbrev %in% region)

  if (is.null(species)) {
    species <- unique(survey_dat$species)
  } else {
    survey_dat <- filter(survey_dat, species == {{species}})
  }
  mesh <- make_mesh(survey_dat, xy_cols = c("X", "Y"), cutoff = cutoff)
  fit <- tryCatch(
    sdmTMB(
      formula = formula,
      data = survey_dat,
      mesh = mesh,
      time = "year",
      spatial = sp,
      spatiotemporal = st,
      offset = "offset",
      family = choose_family(family)
    ),
    error = function(e) paste(species, region, family, "\n\tError:", e, sep = " - ")
  )

  if (inherits(fit, 'sdmTMB')) {
    sanity_check <- all(unlist(sdmTMB::sanity(fit, gradient_thresh = 0.005)))
  }
  # Turn off spatial field if model does not fit and spatiotemporal == "off"
  if ((!inherits(fit, 'sdmTMB') | !sanity_check) & (sp == "on" & st == "off")) {
    message("\tFitting: sp = ", sp, ", st = ", st, " for ", species, "-", region, "-", family, " failed")
    message("\tUpdating with sp = off")
    fit <- tryCatch(
      update(fit, spatial = "off"),
      error = function(e) paste(species, region, family, "\n\tError:", e, sep = " - ")
    )
  }
  fit
}


regions <- c("SYN WCVI", "SYN QCS", "SYN HS", "SYN WCHG")
families <- c("tweedie",
  "delta-gamma", "delta-lognormal", "delta-gengamma")
  #"delta-plink-lognormal", "delta-plink-gamma", "delta-plink-gengamma")
tofit <- tidyr::expand_grid(.region = regions, .family = families,
  .species = spp_list) |>
  mutate(.id = row_number())

#species <- spp_list[[4]]
species <- 'yellowtail rockfish'
sp_file <- paste0(clean_name(species), ".rds")
cutoff <- 8
save_fits <- TRUE

fits <-
  tofit |>
    filter(.species == species) |>
    filter(.region == 'SYN WCHG') |>
    purrr::pmap(\(.region, .family, .id, .species) {
        get_fit(survey_dat = dat,
            formula = as.formula(catch_weight ~ 0 + as.factor(year)),
            offset = "offset",
            region = .region, family = .family, species = .species,
            cutoff = cutoff, sp = "on", st = "off"
            )
    })
beep()
# if (save_fits) {
#   message("Saving fits: ", sp_file)
#   saveRDS(fits, file.path(fit_dir, sp_file))
# }
# beep()

# Individual model fitting for checks
# test <- local({
#   survey_dat <- dat
#   .region <- "SYN QCS"
#   .family <- "tweedie"
#   .cutoff <- 8
#   .species <- "butter sole"
#   get_fit(survey_dat = survey_dat,
#          formula = as.formula(catch_weight ~ 0 + as.factor(year)),
#          offset = "offset",
#          region = .region, family = .family, species = .species,
#          cutoff = cutoff, sp = "on", st = "iid"
#   )
# })
# test

survey_dat <- dat
cutoff <- 8
regions <- c("SYN WCVI", "SYN QCS", "SYN HS", "SYN WCHG")
families <- c("tweedie",
  "delta-gamma", "delta-lognormal", "delta-gengamma")
  #"delta-plink-lognormal", "delta-plink-gamma", "delta-plink-gengamma")
family_lu <- c(
  "tweedie" = "tweedie",
  "delta-gamma" = "Gamma",
  "delta-lognormal" = "lognormal",
  "delta-gengamma" = "gengamma") %>%
  tibble(.family = names(x = .), fit_family = .)
tofit <- tidyr::expand_grid(.region = regions, .family = families,
  .species = spp_list) |>
  mutate(.id = row_number())

lu_df <- tofit |>
  left_join(family_lu) |>
  mutate(sp_hyp = clean_name(.species)) |>
  tidyr::unite(col = "fname", sp_hyp, .region, .family, sep = "-", remove = FALSE)
saveRDS(lu_df, file.path(out_dir, 'lu-df.rds'))
# START WITH st = "off" for now.... in part because I want to better understand:
# https://github.com/pbs-assess/sdmTMB/issues/263
# it is also a lot faster
# Fit models
choose_multi_type(cores = 8)
progressr::with_progress({
  handler <- progressr::progressor(along = tofit$.id)
  fits <- tofit |>
    #purrr::pmap(\(.region, .family, .id, .species) {
    furrr::future_pmap(\(.region, .family, .id, .species) {
      handler() # Signal progress
      sp_file <- file.path(fit_dir, paste0(clean_name(.species), "-", .region, "-", .family, ".rds"))
      if (overwrite_cache | !file.exists(sp_file)) {
        message("Fitting: ", sp_file)
        f <- get_fit(survey_dat = survey_dat,
                formula = as.formula(catch_weight ~ 0 + as.factor(year)),
                offset = "offset",
                region = .region, family = .family, species = .species,
                cutoff = cutoff, sp = "on", st = "off")
        message("Saving fits: ", sp_file)
        saveRDS(f, sp_file)
      }
  })
  beep()
  future::plan(future::sequential())
})

# ------------------------------------------------------------------------------
# Examine outputs
# ------------------------------------------------------------------------------
f_name <- list.files(fit_dir)
f <- file.path(fit_dir, f_name)
# Just arrowtooth for testing code
fa <- lapply(f[grepl("arrowtooth-flounder", f)], readRDS)

fits <-lapply(f, readRDS) |>
  setNames(stringr::str_extract(f_name, ("(.*)\\.rds"), group = 1))

if (!file.exists(file.path(out_dir, "sanity-df.rds"))) {
sanity_df <- purrr::keep(fits, ~inherits(.x, "sdmTMB")) |>
  purrr::map_dfr(get_sanity_df, real_data = TRUE, silent = TRUE) |>
  left_join(lu_df, by = c("species" = ".species", "region" = ".region", "fit_family"))
saveRDS(sanity_df, file.path(out_dir, "sanity-df.rds"))
} else {message("\t Loading cached sanity-df")
  sanity_df <- readRDS(file.path(out_dir, "sanity-df.rds"))
}

ok_sanity <- sanity_df |> filter(all_ok == TRUE)

# Only predict and get index for fits that passed a sanity check
ok_fits <- fits[ok_sanity$fname]

predictions <- purrr::map(ok_fits, \(x) {
  predict(x, newdata = x$data, return_tmb_object = TRUE)
})

inds <- purrr::map(predictions, get_index, bias_correct = TRUE) |>
  purrr::map2(.x = _, .y = names(inds), ~ mutate(.x, fname = .y)) |>
  bind_rows() |>
  left_join(select(lu_df, .species, .region, .family)) |>
  select(-fname)
saveRDS(inds, file.path(out_dir, "index-df.rds"))

fit_ests <- purrr::map(ok_fits, get_fitted_estimates, real_data = TRUE) |>
  purrr::map2(.x = _, .y = names(ok_fits), ~ mutate(.x, fname = .y)) |>
  bind_rows() |>
  left_join(lu_df) |>
  select(-fname)
beep()
saveRDS(fit_ests, file.path(out_dir, "fitted-estimates-df.rds"))

rqr_residuals1 <- purrr::map2_dfr(ok_fits[1:30], names(ok_fits[1:30]), ~get_rqr(.x, id = .y))
rqr_residuals2 <- purrr::map2_dfr(ok_fits[31:100], names(ok_fits[31:100]), ~get_rqr(.x, id = .y))

rqr_residuals3 <- purrr::map2_dfr(ok_fits[101:length(ok_fits)], names(ok_fits[101:length(ok_fits)]),
  ~ {
  message(.y)
  get_rqr(.x, id = .y)
  }
  )


rqr_residuals <- bind_rows(rqr_residuals1, rqr_residuals2, rqr_residuals3)

test <- residuals(ok_fits[[5]], type = 'mle-mvn', model = 2)
str(test)

# RQR
# ------------------
rqr_df <- purrr::map_dfr(1:length(fits), ~ get_rqr(fits[[.x]], id = .x)) |>
  left_join(fit_df)
beep()