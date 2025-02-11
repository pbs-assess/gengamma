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

dat <- survey_dat |> filter(species == "arrowtooth flounder", region == "GOA")
mesh <- sdmTMB::make_mesh(dat, xy_cols = c("X", "Y"), type = "cutoff_search", n_knots = 500)
# .spatial <-
# .spatiotemporal <- docv |> slice(1) |> purrr::pull(".spatiotemporal")
.family <- "delta-gamma"

tictoc::tic()
plan(multisession, workers = 6)
test_gamma <- sdmTMB_cv(
  formula = as.formula(catch_weight ~ 0 + as.factor(year)),
  data = dat,
  mesh = mesh,
  time = "year",
  spatial = list("on", "on"),
  spatiotemporal = list("iid", "iid"),
  offset = "offset",
  family = choose_family(.family),
  k_folds = 12
)
tictoc::toc()
# beepr::beep()
future::plan(future::sequential)


test_gamma$models[[1]] |> sanity()
test_gamma$models[[2]] |> sanity()
test_gamma$models[[3]] |> sanity()
test_gamma$models[[4]] |> sanity()
test_gamma$models[[5]] |> sanity()
test_gamma$models[[6]] |> sanity()
test_gamma$models[[7]] |> sanity()
test_gamma$models[[8]] |> sanity()
test_gamma$models[[9]] |> sanity()
test_gamma$models[[10]] |> sanity()
test_gamma$models[[11]] |> sanity()
