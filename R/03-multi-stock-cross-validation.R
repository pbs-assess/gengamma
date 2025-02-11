# ------------------------------------------------------------------------------
# Cross-validation
# ------------------------------------------------------------------------------
# library(purrr)
# source(here::here("R", "01-multiple-stocks.R")) # check run_trouble_fits <- FALSE
source(here::here("R", "00-utils.R"))

library(sdmTMB)
library(dplyr)
library(future)

# run_cv <- function(survey_dat, formula, region, family, species,
#   cutoff = cutoff, sp, st, kfolds = 10) {

#   survey_dat <- filter(survey_dat, region %in% {{region}})

#   if (is.null(species)) {
#     species <- unique(survey_dat$species)
#   } else {
#     survey_dat <- filter(survey_dat, species == {{species}})
#   }
#   if (region != "GOA"){
#     mesh <- sdmTMB::make_mesh(survey_dat, xy_cols = c("X", "Y"), cutoff = cutoff)
#   } else {
#     mesh <- sdmTMB::make_mesh(survey_dat, xy_cols = c("X", "Y"), type = "cutoff_search", n_knots = 500)
#   }

#   ny <- length(unique(survey_dat$year))
#   .family <- choose_family(family)
#   fit <- tryCatch(
#     sdmTMB::sdmTMB_cv(
#       formula = formula,
#       data = survey_dat,
#       mesh = mesh,
#       time = "year",
#       spatial = sp,
#       spatiotemporal = st,
#       offset = "offset",
#       family = .family,
#       k_folds = kfolds
#     ),
#     error = function(e) paste(species, region, family, "\n\tError:", e, sep = " - ")
#   )
# }

# lu_df <- readRDS(file.path(here::here("data-outputs", "multi-species", "lu-df.rds"))) |>
#   janitor::clean_names()

# sanity_df <- readRDS(here::here("data-outputs", "multi-species", "sanity-df.rds"))
# # tr_species <- sanity_df |> filter(all_ok == FALSE) |> select(species, region, fit_family)

# docv <- sanity_df |>
#   right_join(lu_df) |>
#   filter(all_ok)
# docv$sp_list <- lapply(docv$spatial, function(x) strsplit(x, "-")[[1]])
# docv$st_list <- lapply(docv$spatiotemporal, function(x) strsplit(x, "-")[[1]])

# docv <- docv |>
#   select(species, region, family, spatial, spatiotemporal, sp_list, st_list) %>%
#   setNames(paste0(".", names(.)))

# pluck(docv, "spatial_list", 1) |> as.list()

survey_dat <- readRDS(here::here("data-outputs", "data-used.rds"))
# survey_dat <- readRDS("~/Downloads/data-used.rds")
survey_dat <- filter(survey_dat, !(region == "GOA" & species == "petrale sole"))
unique(survey_dat$species)
unique(survey_dat$region)

select(survey_dat, species, region) |>
  distinct() |> as.data.frame()

# dat <- survey_dat |> filter(species == "arrowtooth flounder", region == "GOA")
# plot(mesh)
# .spatial <-
# .spatiotemporal <- docv |> slice(1) |> purrr::pull(".spatiotemporal")

fit_all <- function(.dat, .family) {
  set.seed(1199328)
  if (unique(.dat$region) == "GOA") {
    mesh <- sdmTMB::make_mesh(dat, xy_cols = c("X", "Y"), type = "cutoff_search", n_knots = 500)
  } else {
    mesh <- sdmTMB::make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 8)
  }
  fit <- tryCatch({sdmTMB_cv(
    formula = as.formula(catch_weight ~ 0 + as.factor(year)),
    data = .dat,
    mesh = mesh,
    time = "year",
    spatial = "on",
    spatiotemporal = "iid",
    offset = "offset",
    family = choose_family(.family),
    control = sdmTMBcontrol(multiphase = FALSE, profile = TRUE),
    k_folds = 25
  )}, error = function(e) NA)
  if (length(fit) > 1) {
    data.frame(family = .family, region = unique(.dat$region), species = unique(.dat$species), loglik = fit$sum_loglik)
  } else {
    data.frame(family = .family, region = unique(.dat$region), species = unique(.dat$species), loglik = NA_real_)
  }
}

dsplit <- group_by(survey_dat, species, region) |> group_split()
plan(multisession, workers = 6)
df_lg <- purrr::map_dfr(dsplit, fit_all, .family = "delta-lognormal")
df_dg <- purrr::map_dfr(dsplit, fit_all, .family = "delta-gamma")
df_tw <- purrr::map_dfr(dsplit, fit_all, .family = "tweedie")
df_dgg <- purrr::map_dfr(dsplit, fit_all, .family = "delta-gengamma")
future::plan(future::sequential)

saveRDS(df_lg, file = here::here("data-outputs", "cross-val-lognormal.rds"))
saveRDS(df_dg, file = here::here("data-outputs", "cross-val-gamma.rds"))
saveRDS(df_tw, file = here::here("data-outputs", "cross-val-tweedie.rds"))
saveRDS(df_dgg, file = here::here("data-outputs", "cross-val-gengamma.rds"))

