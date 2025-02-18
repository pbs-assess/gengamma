# ------------------------------------------------------------------------------
# Cross-validation
# ------------------------------------------------------------------------------
source(here::here("R", "00-utils.R"))

library(sdmTMB)
library(dplyr)
library(future)

cv_dir <- here::here("data-outputs", "cv")
dir.create(cv_dir, showWarnings = FALSE)

lu_df <- readRDS(file.path(here::here("data-outputs", "multi-species", "lu-df.rds"))) |>
  janitor::clean_names()

sanity_df <- readRDS(here::here("data-outputs", "multi-species", "sanity-df.rds"))

docv <- sanity_df |>
  right_join(lu_df) |>
  filter(all_ok)
docv$sp_list <- lapply(docv$spatial, function(x) as.list(strsplit(x, "-")[[1]]))
docv$st_list <- lapply(docv$spatiotemporal, function(x) as.list(strsplit(x, "-")[[1]]))

if (!file.exists(here::here("data-outputs", "cv-lu.rds"))) {
  docv <- docv |>
    select(species, region, family, spatial, spatiotemporal, sp_list, st_list) %>%
    setNames(paste0(".", names(.)))
  saveRDS(docv, here::here("data-outputs", "cv-lu.rds"))
} else {
  docv <- readRDS(here::here("data-outputs", "cv-lu.rds"))
}

survey_dat <- readRDS(here::here("data-outputs", "data-used.rds"))
# survey_dat <- readRDS("~/Downloads/data-used.rds")
# survey_dat <- filter(survey_dat, !(region == "GOA" & species == "petrale sole"))


fit_all <- function(.dat, .species, .region, .family, .spatial, .spatiotemporal) {
  set.seed(11199328)
  .dat <- filter(.dat, species == .species, region == .region)

  if (.region == "GOA") {
    mesh <- sdmTMB::make_mesh(.dat, xy_cols = c("X", "Y"), type = "cutoff_search", n_knots = 500)
  } else {
    mesh <- sdmTMB::make_mesh(.dat, xy_cols = c("X", "Y"), cutoff = 8)
  }
  fit <- tryCatch(
    {
      sdmTMB_cv(
        formula = as.formula(catch_weight ~ 0 + as.factor(year)),
        data = .dat,
        mesh = mesh,
        time = "year",
        spatial = .spatial,
        spatiotemporal = .spatiotemporal,
        offset = "offset",
        family = choose_family(.family),
        control = sdmTMBcontrol(multiphase = FALSE, profile = TRUE),
        k_folds = 50
      )
    },
    error = function(e) NA
  )
  if (length(fit) > 1) {
    data.frame(
      family = .family,
      region = unique(.dat$region),
      species = unique(.dat$species),
      sumloglik = fit$sum_loglik,
      foldloglik = paste(round(fit$fold_loglik, 5), collapse = ";")
    )
  } else {
    data.frame(
      family = .family,
      region = unique(.dat$region),
      species = unique(.dat$species),
      sumloglik = NA_real_,
      foldloglik = ""
    )
  }
}

cores <- availableCores()
plan(multicore, workers = min(cores, 50L))

out <- docv |>
  # filter(.region == "GOA") |>
  # filter(.species == "petrale sole") |>
  select(.species, .region, .family, .sp_list, .st_list) |>
  purrr::pmap(\(.species, .region, .family, .sp_list, .st_list) {
    fname <- file.path(cv_dir, paste0(.species, "-", .region, "-", .family, ".rds"))
    if (!file.exists(fname)) {
      message("Running: ", fname)
      out <- fit_all(.dat = survey_dat, .species, .region, .family, .sp_list, .st_list)
      saveRDS(out, fname)
    }
  })
future::plan(future::sequential)
# beepr::beep()

cv <- lapply(list.files(cv_dir, full.names = TRUE), readRDS) |>
  bind_rows() #|>
# filter(!(species %in% c("petrale sole", "walleye pollock")))

saveRDS(cv, here::here("data-outputs", "cv-summary-df.rds"))
