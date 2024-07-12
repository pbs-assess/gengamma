library(dplyr)
library(sdmTMB) # pak::pkg_install("pbs-assess/sdmTMB")

source("R/00-utils.R")

# Load data
dat <- readRDS(file.path("data", "clean-survey-data.rds")) |>
  filter(survey_abbrev != "SYN WCHG") |>
  filter(prop_pos >= 0.2) |>
  group_by(species) |>
  mutate(n_regions = n()) |>
  ungroup()

spp_regions <- distinct(dat, species, survey_abbrev) |>
  rename(.species = "species", .region = "survey_abbrev") |>
    group_by(.species) |>
    mutate(n_regions = n()) |>
    ungroup()

# Simplified synoptic grid
syn_grid <- as_tibble(gfplot::synoptic_grid) |>
  select(survey, X, Y, depth)

# tag <- "_sp-on-st-iid"
# tag <- "_sp-on-st-off"
.spatial <- "on"
.spatiotemporal <- "iid"
#tag <- paste0("_sp-", .spatial, "-st-", .spatiotemporal, "", sep = "")
tag <- ""
dc <- here::here("data", "raw")
out_dir <- here::here("data-outputs", "bcgf-outputs")
fit_dir <- here::here("data-outputs", "bcgf-outputs", paste0("fits", tag))
dir.create(fit_dir, showWarnings = FALSE, recursive = TRUE)

# For fitting select species/regions/families
# ------------------------------------------------------------------------------
# regions <- c("SYN WCVI", "SYN QCS", "SYN HS", "SYN WCHG")
# families <- c("tweedie",
#   "delta-gamma", "delta-lognormal", "delta-gengamma")
#   #"delta-plink-lognormal", "delta-plink-gamma", "delta-plink-gengamma")
# tofit <- tidyr::expand_grid(.region = regions, .family = families,
#   .species = unique(spp_regions$.species)) |>
#   right_join(spp_regions) |>
#   mutate(.id = row_number())

# species <- spp_list[[4]]
# species <- "north pacific spiny dogfish"
# scaler <- 1e-2 # works
# region <- "SYN WCVI"

# species <- "walleye pollock"
# region <- "SYN WCVI"
# scaler <- 1e-3 # QCS needs 1e-1; nothing seems to fix SYN WCVI

# sp_file <- paste0(clean_name(species), ".rds")
# cutoff <- 8
# save_fits <- TRUE

# get_fit <- function(survey_dat, formula, region, family, species = NULL, cutoff = 8, time = "year",
#                    sp = "on", st = "iid", offset = "offset") {
#   survey_dat <- filter(survey_dat, survey_abbrev %in% region)

#   if (is.null(species)) {
#     species <- unique(survey_dat$species)
#   } else {
#     survey_dat <- filter(survey_dat, species == {{species}})
#   }
#   mesh <- make_mesh(survey_dat, xy_cols = c("X", "Y"), cutoff = cutoff)
#   fit <- tryCatch(
#     sdmTMB(
#       formula = formula,
#       data = survey_dat,
#       mesh = mesh,
#       time = "year",
#       spatial = sp,
#       spatiotemporal = st,
#       offset = "offset",
#       family = choose_family(family)
#     ),
#     error = function(e) paste(species, region, family, "\n\tError:", e, sep = " - ")
#   )

#   if (inherits(fit, "sdmTMB")) {
#     sanity_check <- all(unlist(sdmTMB::sanity(fit, gradient_thresh = 0.005)))
#   }
#   # Turn off spatial field if model does not fit and spatiotemporal == "off"
#   if ((!inherits(fit, "sdmTMB") | !sanity_check) & (sp == "on" & st == "off")) {
#     message("\tFitting: sp = ", sp, ", st = ", st, " for ", species, "-", region, "-", family, " failed")
#     message("\tUpdating with sp = off")
#     fit <- tryCatch(
#       update(fit, spatial = "off"),
#       error = function(e) paste(species, region, family, "\n\tError:", e, sep = " - ")
#     )
#   }
#   fit
# }

# regions <- c("SYN WCVI", "SYN QCS", "SYN HS")
# families <- c("tweedie", "delta-gamma", "delta-lognormal", "delta-gengamma")
# # "delta-plink-lognormal", "delta-plink-gamma", "delta-plink-gengamma")
# tofit <- tidyr::expand_grid(
#   .region = regions, .family = families,
#   .species = spp_list
# ) |>
#   mutate(.id = row_number())

# test_fit <-
#   tofit |>
#     filter(.species == species) |>
#     filter(.region == region) |>
#     filter(.family == "delta-gengamma") |>
#     purrr::pmap(\(.region, .family, .id, .species) {
#         get_fit(survey_dat = dat |> mutate(catch_weight = catch_weight * scaler),
#             formula = as.formula(catch_weight ~ 0 + as.factor(year)),
#             offset = "offset",
#             region = .region, family = .family, species = .species,
#             cutoff = cutoff, sp = "on", st = "off"
#             )
#     })
# beep()
# if (save_fits) {
#   message("Saving fits: ", sp_file)
#   saveRDS(fits, file.path(fit_dir, sp_file))
# }
# beep()
# ------------------------------------------------------------------------------

survey_dat <- dat
cutoff <- 8
regions <- c("SYN WCVI", "SYN QCS", "SYN HS")
families <- c(
  "tweedie",
  "delta-gamma", "delta-lognormal", "delta-gengamma",
  "delta-gamma-poisson-link", "delta-lognormal-poisson-link", "delta-gengamma-poisson-link"
)
# "delta-plink-lognormal", "delta-plink-gamma", "delta-plink-gengamma")
family_lu <- c(
  "tweedie" = "tweedie",
  "delta-gamma" = "Gamma",
  "delta-lognormal" = "lognormal",
  "delta-gengamma" = "gengamma",
  "delta-gamma-poisson-link" = "Gamma",
  "delta-lognormal-poisson-link" = "lognormal",
  "delta-gengamma-poisson-link" = "gengamma"
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

# Turn off spatial/spatiotemporal random fields if they collapse
# https://github.com/pbs-assess/sdmTMB/issues/263
# overwrite_cache = TRUE
overwrite_cache <- FALSE
# choose_multi_type(cores = 8)
future::plan(future::multisession, workers = 10)
progressr::with_progress({
  handler <- progressr::progressor(along = tofit$.id)
  fits <- tofit |>
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
