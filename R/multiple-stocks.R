library(dplyr)
#devtools::load_all('../sdmTMB')
#pak::pkg_install("pbs-assess/sdmTMB")
library(sdmTMB)
library(ggplot2)
theme_set(gfplot::theme_pbs())

source("R/00-utils.R")

# Load data
dat <- readRDS(file.path("data-outputs", "clean-survey-data.rds")) |>
  filter(survey_abbrev != "SYN WCHG")
spp_regions <- distinct(dat, species, survey_abbrev) |>
  rename(.species = 'species', .region = 'survey_abbrev')

# Simplified synoptic grid
syn_grid <- as_tibble(gfplot::synoptic_grid) |>
  select(survey, X, Y, depth)

#tag <- '_sp-on-st-iid'
#tag <- '_sp-on-st-off'
.spatial <- "off"
.spatiotemporal <- "off"
tag <- paste0("_sp-", .spatial, "-st-", .spatiotemporal, "", sep = "")
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

#species <- spp_list[[4]]
#species <- 'north pacific spiny dogfish'
#scaler <- 1e-2 # works
#region <- 'SYN WCVI'

# species <- 'walleye pollock'
# region <- 'SYN WCVI'
# scaler <- 1e-3 # QCS needs 1e-1; nothing seems to fix SYN WCVI

#sp_file <- paste0(clean_name(species), ".rds")
# cutoff <- 8
# save_fits <- TRUE

# test_fit <-
#   tofit |>
#     filter(.species == species) |>
#     filter(.region == region) |>
#     filter(.family == 'delta-gengamma') |>
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
  .species = unique(spp_regions$.species)) |>
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
#overwrite_cache = TRUE
overwrite_cache = FALSE
#choose_multi_type(cores = 8)
future::plan(future::multisession, workers = 10)
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
                cutoff = cutoff, sp = .spatial, st = .spatiotemporal,
                use_priors = FALSE)
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
f_name <- c(
  list.files(file.path(out_dir, 'fits_sp-off-st-off')),
  list.files(file.path(out_dir, 'fits_sp-on-st-off')),
  list.files(file.path(out_dir, 'fits_sp-on-st-on'))
)
f_name <- f_name[!grepl('WCHG', f_name)]
f <- file.path(fit_dir, f_name)

fits <-lapply(f, readRDS) |>
  setNames(stringr::str_extract(f_name, ("(.*)\\.rds"), group = 1))
beep()

#tag
tag <- "-random-effects-on-and-off"
sanity_filename <- paste0("sanity-df", tag, ".rds")
if (!file.exists(file.path(out_dir, sanity_filename))) {
sanity_df <- purrr::keep(fits, ~inherits(.x, "sdmTMB")) |>
  purrr::map_dfr(get_sanity_df, real_data = TRUE, silent = TRUE, .gradient_thresh = 0.005) |>
  left_join(lu_df, by = c("species" = ".species", "region" = ".region", "fit_family"))
saveRDS(sanity_df, file.path(out_dir, sanity_filename))
} else {message("\t Loading cached sanity-df")
  sanity_df <- readRDS(file.path(out_dir, sanity_filename))
}

ok_sanity <- sanity_df |> filter(all_ok == TRUE)

# Only predict and get index for fits that passed a sanity check
ok_fits <- fits[ok_sanity$fname]

predictions <- purrr::map(ok_fits, \(x) {
  region <- unique(x$data$survey_abbrev)
  years <- unique(x$data$year)
  sg <- syn_grid |> filter(survey %in% region)
  nd <- sdmTMB::replicate_df(
      dat = sg,
      time_name = "year",
      time_values = years)
  predict(x, newdata = nd, return_tmb_object = TRUE)
})

tictoc::tic()
inds3 <- purrr::map(predictions[301:472], get_index, bias_correct = TRUE) |>
  purrr::map2(.x = _, .y = names(predictions[301:472]), ~ mutate(.x, fname = .y)) |>
  bind_rows() #|>
tictoc::toc()
  left_join(lu_df) |>
  select(-fname)

index_filename <- paste0("index-df", tag, ".rds")
bind_rows(inds, inds2, inds3) |>
saveRDS(file.path(out_dir, index_filename))

fit_est_filename <- paste0("fitted-estimates", tag, ".rds")
fit_ests <- purrr::map(fits, get_fitted_estimates, real_data = TRUE) |>
  purrr::keep(is.data.frame) |>
  purrr::imap(~ mutate(.x, fname = .y)) |>
  bind_rows() |>
  left_join(lu_df) |>
  select(-fname)
beep()
saveRDS(fit_ests, file.path(out_dir, fit_est_filename))


test <- left_join(fit_ests, sanity_df)

test |> filter(all_ok == FALSE) |>
  View()

# gg_fits <- ok_fits[grepl('gengamma', names(ok_fits))] |>
#   purrr::keep(~ inherits(.x, 'sdmTMB'))
# gg_fit_ests <- gg_fits|>
#   purrr::map(get_fitted_estimates, real_data = TRUE) |>
#   purrr::map2(.x = _, .y = names(gg_fits), ~ mutate(.x, fname = .y)) |>
#   bind_rows() |>
#   left_join(lu_df) |>
#   select(-fname)
# beep()
# saveRDS(gg_fit_ests, file.path(out_dir, "gengamma-all-fitted-estimates-df.rds"))

gg_fit_ests

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

# Examine single species
sp_f <- f[grepl("north-pacific-spiny-dogfish", f)]
fa <- lapply(sp_f, readRDS) |>
  setNames(stringr::str_extract(sp_f, (".*fits/(.*)\\.rds"), group = 1))
sp_fits <- purrr::keep(fa, ~inherits(.x, "sdmTMB"))
test <- sp_fits |>
  purrr::map_dfr(get_sanity_df, real_data = TRUE, .gradient_thresh = 0.005)

g <- sp_fits |> purrr::map(~ .x$gradients)


sp_fits[["north-pacific-spiny-dogfish-SYN WCVI-delta-gengamma"]] |> names()

rqr_df <- purrr::map_dfr(1:length(sp_fits), ~ get_rqr(sp_fits[[.x]], id = names(sp_fits[.x])))

rqr_df |>
  mutate(region = stringr::str_extract_all(id, "(.*)-(SYN [A-Z]{2,4})-(.*)"))

rqr_df  |>
   mutate(species = stringr::str_extract(id, ".*(?=-SYN)"),
          region = stringr::str_extract(id, ".*-(SYN [A-Z]{2,4})-.*", group = 1),
          family = stringr::str_extract(id, ".*-SYN [A-Z]{2,4}-(.*)", group = 1)) |>
ggplot(aes(sample = r)) +
  geom_qq() +
  geom_abline(intercept = 0, slope = 1) +
  facet_grid(region ~ family)
  #facet_wrap(~ family)
ggsave(filename = file.path('figures', 'bc-gf-data', 'exploratory', 'dogfish-rqr.png'))


# Splitting the capture groups into separate columns




# RQR
# ------------------
rqr_df <- purrr::map_dfr(1:length(fits), ~ get_rqr(fits[[.x]], id = .x)) |>
  left_join(fit_df)
beep()