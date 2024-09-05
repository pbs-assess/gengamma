library(dplyr)
library(sdmTMB)

source("R/00-utils.R")

out_dir <- here::here("data-outputs", "multi-species")
fit_dir <- here::here("data-outputs", "multi-species", "fits")

lu_df <- readRDS(file.path(out_dir, "lu-df.rds")) |>
  janitor::clean_names()

# Simplified synoptic grid
syn_grid <- as_tibble(gfplot::synoptic_grid) |>
  select(survey, X, Y, depth)

# Load GoA grid
data(gulf_of_alaska_grid, package = "FishStatsUtils")
goa_grid <- as_tibble(gulf_of_alaska_grid) |>
  sdmTMB::add_utm_columns(c("Lon", "Lat"), utm_crs = 32609)
rm(list = "gulf_of_alaska_grid")

# ------------------------------------------------------------------------------
# Examine outputs
# ------------------------------------------------------------------------------
f <- list.files(file.path(out_dir, "fits"), full.names = TRUE)
# overwrite_outputs <- TRUE
overwrite_outputs <- FALSE
tag <- ""

fit_est_filename <- paste0("fitted-estimates", tag, ".rds")
pb <- cli::cli_progress_bar("Processing", total = length(f))
if (overwrite_outputs | !file.exists(file.path(out_dir, fit_est_filename))) {
  fit_ests <- f |>
    purrr::map(\(f) {
      cli::cli_progress_update(id = pb)
      fit <- readRDS(f)
      get_fitted_estimates(fit, real_data = TRUE)
    }
    ) |>
    purrr::keep(is.data.frame) |>
    bind_rows() |>
    mutate(type = gsub("poisson_link_delta", "poisson-link", type)) |>
    mutate(type = ifelse(fit_family == "tweedie", "standard", type))
  saveRDS(fit_ests, file.path(out_dir, fit_est_filename))
  beep()
} else {
  message("\t Loading cached fit_ests: ", fit_est_filename)
  fit_ests <- readRDS(file.path(out_dir, fit_est_filename))
}
cli::cli_progress_done(id = pb)

pb <- cli::cli_progress_bar("Processing", total = length(f))
sanity_filename <- paste0("sanity-df", tag, ".rds")
if (overwrite_outputs | !file.exists(file.path(out_dir, sanity_filename))) {
  sanity_df <- f |>
    purrr::map(\(f) {
      cli::cli_progress_update(id = pb)
      fit <- readRDS(f)
      if (inherits(fit, "sdmTMB")) {
        get_sanity_df(fit, real_data = TRUE, silent = TRUE, .gradient_thresh = 0.005)
      } else {
        fit
      }
    }) |>
    purrr::keep(~ is.data.frame(.x)) |>
    bind_rows() |>
    mutate(type = gsub("poisson_link_delta", "poisson-link", type)) |>
    mutate(type = ifelse(fit_family == "tweedie", "standard", type))
  saveRDS(sanity_df, file.path(out_dir, sanity_filename))
  beep()
} else {
  message("\t Loading cached sanity-df")
  sanity_df <- readRDS(file.path(out_dir, sanity_filename))
}
cli::cli_progress_done(id = pb)

ok_sanity <- sanity_df |>
  right_join(lu_df) |>
  filter(all_ok) |>
  mutate(fname = paste0(fname, '.rds'))

# Only predict and get index for fits that passed a sanity check

ff <- ok_sanity$fname
pb <- cli::cli_progress_bar("Processing", total = length(ff))
predictions <- purrr::map(ff, \(fit_file) {
  cli::cli_progress_update(id = pb)
  fit <- readRDS(file.path(fit_dir, fit_file))
  region <- unique(fit$data$region)
  years <- unique(fit$data$year)
  g <- switch(region,
  "GOA" = goa_grid,
  "HS-QCS" = syn_grid |> filter(survey %in% c("SYN HS", "SYN QCS")) |> mutate(region = "HS-QCS"),
  syn_grid |> filter(survey %in% region)
  )
  nd <- sdmTMB::replicate_df(
    dat = g,
    time_name = "year",
    time_values = years
  )
  predict(fit, newdata = nd, return_tmb_object = TRUE)
}) |>
  setNames(ff)
cli::cli_progress_done(id = pb)
beep()

# Slow
idx <- 1:length(predictions)
pb <- cli::cli_progress_bar("Processing", total = length(predictions[idx]))
ind <- purrr::map(predictions[idx], \(p) {
  cli::cli_progress_update(id = pb)
  area <- ifelse(unique(p$fit_obj$data$region) == "GOA", goa_grid$Area_in_survey_km2, 4)
  get_index(p, area = area, bias_correct = TRUE)
  }) |>
  purrr::map2(.x = _, .y = names(predictions[idx]), ~ mutate(.x, fname = gsub('\\.rds', '', .y)))
index_df <- bind_rows(ind)
cli::cli_progress_done(id = pb)
beep()

index_filename <- paste0("index-df", tag, ".rds")
saveRDS(index_df, file.path(out_dir, index_filename))


# ------------------------------------------------------------------------------
# Residuals
# -------------
ff <- ok_sanity$fname
pb <- cli::cli_progress_bar("Processing", total = length(ff))
rqr_catch <- purrr::map(ff, \(fit_file) {
  cli::cli_progress_update(id = pb)
  fit <- readRDS(file.path(fit_dir, fit_file))
  fname <- gsub("\\.rds", "", fit_file)
  get_rqr(fit, id = fname) |>
    rename(fname = "id")
})
cli::cli_progress_done(id = pb)
beep()
rqr_catch <- bind_rows(rqr_catch)
rqr_catch_filename <- paste0("rqr-catch-df", tag, ".rds")
saveRDS(rqr_catch, file.path(out_dir, rqr_catch_filename))

# ------------------------------------------------------------------------------
# Look at the cases where we tried using a prior
# ------------------------------------------------------------------------------
prior_fit_dir <- file.path(out_dir, "priors", "fits")

# Prior fits
f <- list.files(prior_fit_dir, full.names = TRUE)
p_fit_ests <- f |>
  purrr::map(\(f) {
    fit <- readRDS(f)
    ests <- get_fitted_estimates(readRDS((f)), real_data = TRUE)
    ests$prior_b_sd <- unique(fit$priors$b)[, 2]
    ests
    }) |>
  purrr::keep(is.data.frame) |>
  bind_rows() |>
  mutate(type = gsub("poisson_link_delta", "poisson-link", type)) |>
  mutate(type = ifelse(fit_family == "tweedie", "standard", type)) |>
  mutate(priors = TRUE)
saveRDS(p_fit_ests, file.path(out_dir, "priors", "fitted-estimates.rds"))

# Prior sanity
p_sanity_df <- f |> purrr::map(\(f) {
  fit <- readRDS(f)
  if (inherits(fit, "sdmTMB")) {
    get_sanity_df(fit, real_data = TRUE, silent = TRUE, .gradient_thresh = 0.005)
  }}) |>
  purrr::keep(~ is.data.frame(.x)) |>
  bind_rows() |>
  mutate(type = gsub("poisson_link_delta", "poisson-link", type),
         type = ifelse(fit_family == "tweedie", "standard", type),
         priors = TRUE)
saveRDS(p_sanity_df, file.path(out_dir, "priors", "sanity-df.rds"))

# Prior index (only for sanity checks that passed)
ff <- p_sanity_df |>
  right_join(lu_df) |>
  filter(all_ok) |>
  mutate(fname = paste0(fname, '.rds')) |>
  pull("fname")
predictions <- purrr::map(ff, \(fit_file) {
  message(fit_file)
  fit <- readRDS(file.path(prior_fit_dir, fit_file))
  region <- unique(fit$data$region)
  years <- unique(fit$data$year)
  g <- switch(region,
  "GOA" = goa_grid,
  "HS-QCS" = syn_grid |> filter(survey %in% c("SYN HS", "SYN QCS")) |> mutate(region = "HS-QCS"),
  syn_grid |> filter(survey %in% region)
  )
  nd <- sdmTMB::replicate_df(
    dat = g,
    time_name = "year",
    time_values = years
  )
  predict(fit, newdata = nd, return_tmb_object = TRUE)
}) |>
  setNames(ff)
beep()

# Slow (Should filter out petrale since petrale didn't fit...)
p_ind <- purrr::map(predictions, \(p) {
  area <- ifelse(unique(p$fit_obj$data$region) == "GOA", goa_grid$Area_in_survey_km2, 4)
  get_index(p, area = area, bias_correct = TRUE)
}) |>
  purrr::map2(.x = _, .y = names(predictions), ~ mutate(.x, fname = gsub('\\.rds', '', .y)))
p_index_df <- bind_rows(p_ind) |> mutate(priors = TRUE)

beep()
saveRDS(p_index_df, file.path(out_dir, "priors", "index-df.rds"))

p_rqr_catch <- purrr::map(ff, \(fit_file) {
  fit <- readRDS(file.path(prior_fit_dir, fit_file))
  fname <- gsub("\\.rds", "", fit_file)
  get_rqr(fit, id = fname) |> rename(fname = "id")
}) |>
  bind_rows() |>
  mutate(priors = TRUE)
beep()
saveRDS(p_rqr_catch, file.path(out_dir, "priors", "rqr-catch-df.rds"))

# ---

# Look at the cases where we tried scaling the response
# ------------------------------------------------------------------------------
scale_fit_dir <- file.path(out_dir, "scaled", "fits")

# Pollock scaled fits
f <- list.files(scale_fit_dir, full.names = TRUE)
s_fit_ests <- f |>
  purrr::map(\(f) {
    scale_factor <- stringr::str_extract(f, "(\\d*)\\.rds$", group = 1)
    get_fitted_estimates(readRDS((f)), real_data = TRUE) |>
      mutate(scale_factor = as.numeric(scale_factor))
    }) |>
  purrr::keep(is.data.frame) |>
  bind_rows() |>
  mutate(type = gsub("poisson_link_delta", "poisson-link", type)) |>
  mutate(type = ifelse(fit_family == "tweedie", "standard", type))
saveRDS(s_fit_ests, file.path(out_dir, "scaled", "fitted-estimates.rds"))

# Scaled sanity
s_sanity_df <- f |> purrr::map(\(f) {
  fit <- readRDS(f)
  scale_factor <- stringr::str_extract(f, "(\\d*)\\.rds$", group = 1)
  if (inherits(fit, "sdmTMB")) {
    get_sanity_df(fit, real_data = TRUE, silent = TRUE, .gradient_thresh = 0.005) |>
      mutate(scale_factor = as.numeric(scale_factor))
  }}) |>
  purrr::keep(~ is.data.frame(.x)) |>
  bind_rows() |>
  mutate(type = gsub("poisson_link_delta", "poisson-link", type),
         type = ifelse(fit_family == "tweedie", "standard", type))
saveRDS(s_sanity_df, file.path(out_dir, "scaled", "sanity-df.rds"))

# Sanity index (only for sanity checks that passed)
ff <- s_sanity_df |>
  right_join(lu_df) |>
  filter(all_ok) |>
  mutate(fname = paste0(fname, "-", scale_factor, '.rds')) |>
  pull("fname")
predictions <- purrr::map(ff, \(fit_file) {
  message(fit_file)
  fit <- readRDS(file.path(scale_fit_dir, fit_file))
  region <- unique(fit$data$region)
  years <- unique(fit$data$year)
  g <- switch(region,
  "GOA" = goa_grid,
  "HS-QCS" = syn_grid |> filter(survey %in% c("SYN HS", "SYN QCS")) |> mutate(region = "HS-QCS"),
  syn_grid |> filter(survey %in% region)
  )
  nd <- sdmTMB::replicate_df(
    dat = g,
    time_name = "year",
    time_values = years
  )
  predict(fit, newdata = nd, return_tmb_object = TRUE)
}) |>
  setNames(ff)
beep()

# Slow (Should filter out petrale since petrale didn't fit...)
s_ind <- purrr::map(predictions, \(p) {
  area <- ifelse(unique(p$fit_obj$data$region) == "GOA", goa_grid$Area_in_survey_km2, 4)
  get_index(p, area = area, bias_correct = TRUE)
}) |>
  purrr::map2(.x = _, .y = names(predictions), ~ mutate(.x, fname = gsub('\\.rds', '', .y)))
s_index_df <- bind_rows(s_ind) |>
  mutate(scale_factor = as.numeric(stringr::str_extract(fname, "\\d+$")),
    fname = stringr::str_remove(fname, "-\\d+$")
  )
beep()

saveRDS(s_index_df, file.path(out_dir, "scaled", "index-df.rds"))

s_rqr_catch <- purrr::map(ff, \(fit_file) {
  fit <- readRDS(file.path(scale_fit_dir, fit_file))
  fname <- gsub("\\.rds", "", fit_file)
  get_rqr(fit, id = fname) |> rename(fname = "id")
}) |>
  bind_rows() |>
  mutate(scale_factor = as.numeric(stringr::str_extract(fname, "\\d+$")),
    fname = stringr::str_remove(fname, "-\\d+$")
  )
beep()
saveRDS(s_rqr_catch, file.path(out_dir, "scaled", "rqr-catch-df.rds"))
# ----

# additional exclusion criterion is if index CV is > 1
i |>
  group_by(fname, id) |>
  summarise(mean_ind_cv = mean(sqrt(exp(se^2) - 1)))

# Frequency of estimated Q
.q <- fe |>
  filter(family == "delta-gengamma") |>
  select(species, region, est_q, est_qse) |>
  arrange(est_q) |>
  mutate(q_rank = row_number())

# ------
tr_aic_df <- fe |>
  group_by(species, region) |>
  mutate(
    min_aic = min(aic),
    daic = aic - min_aic,
    rel_lik = exp(-0.5 * daic)
  ) |>
  mutate(aic_w = rel_lik / sum(rel_lik)) |> # see: https://atsa-es.github.io/atsa-labs/sec-uss-comparing-models-with-aic-and-model-weights.html
  ungroup() |>
  left_join(.q) |>
  arrange(q_rank) |>
  mutate(priors = TRUE)
saveRDS(tr_aic_df, file.path(tr_df_dir, "aic-df.rds"))


scale_fit_dir <- file.path(out_dir, "scale", "fits")

# ---
# old scratch


tr_f <- left_join(s, lu_df) |> pull('fname')
tr_pred <- tr_f |>
  purrr::map(\(fit_file) {
    fit <- readRDS(file.path(tr_fit_dir, paste0(fit_file, '.rds')))
    region <- unique(fit$data$region)
    years <- unique(fit$data$year)
    sg <- syn_grid |> filter(survey %in% region)
    nd <- sdmTMB::replicate_df(
      dat = sg,
      time_name = "year",
      time_values = years
    )
    predict(fit, newdata = nd, return_tmb_object = TRUE)
}) |>
  setNames(tr_f)
beep()

i <- purrr::map(tr_pred, \(p) {
  get_index(p, area = 4, bias_correct = TRUE)
}) |>
  purrr::map2(.x = _, .y = names(tr_pred), ~ mutate(.x, fname = gsub('\\.rds', '', .y))) |>
  bind_rows()
i
saveRDS(i, file.path(out_dir,  'priors', "index-df.rds"))


tr_df_dir <- file.path(out_dir, 'priors')
# fit estimates of models that converged
fe <- readRDS(file.path(tr_df_dir, "fit-ests-df.rds")) |>
  janitor::clean_names() |>
  mutate(type = gsub("poisson_link_delta", "poisson-link", type)) |>
  filter(type == "standard") |>
  left_join(lu_df)

# sanity summary for all models fit (not necessarily converged)
s <- readRDS(file.path(tr_df_dir, "sanity-df.rds")) |>
  janitor::clean_names() |>
  filter(type == "standard") |>
  left_join(lu_df)

# indices of models that converged
i <- readRDS(file.path(tr_df_dir, "index-df.rds")) |>
  janitor::clean_names() |>
  filter(!stringr::str_detect(fname, 'poisson-link')) |>
  select(-type) |>
  left_join(lu_df) |>
  as_tibble()

# additional exclusion criterion is if index CV is > 1
i |>
  group_by(fname, id) |>
  summarise(mean_ind_cv = mean(sqrt(exp(se^2) - 1)))

# Frequency of estimated Q
.q <- fe |>
  filter(family == "delta-gengamma") |>
  select(species, region, est_q, est_qse) |>
  arrange(est_q) |>
  mutate(q_rank = row_number())

# ------
tr_aic_df <- fe |>
  group_by(species, region) |>
  mutate(
    min_aic = min(aic),
    daic = aic - min_aic,
    rel_lik = exp(-0.5 * daic)
  ) |>
  mutate(aic_w = rel_lik / sum(rel_lik)) |> # see: https://atsa-es.github.io/atsa-labs/sec-uss-comparing-models-with-aic-and-model-weights.html
  ungroup() |>
  left_join(.q) |>
  arrange(q_rank) |>
  mutate(priors = TRUE)
saveRDS(tr_aic_df, file.path(tr_df_dir, "aic-df.rds"))

# ------ functionify trouble fits
get_fit_outputs <- function(folder_tag) {
  tr_fit_dir <- file.path("data-outputs/multi-species/fits", folder_tag)
  dir.create(file.path("data-outputs/multi-species/", folder_tag), showWarnings = FALSE)
  tr_out_dir <- file.path("data-outputs/multi-species/", folder_tag)
  f <- list.files(tr_fit_dir, full.names = TRUE)

  fe <- purrr::map_df(f, \(file) {
      fit <- readRDS(file)
      get_fitted_estimates(fit, real_data = TRUE)
    }) |>
    mutate(type = gsub("poisson_link_delta", "poisson-link", type)) |>
    mutate(type = ifelse(fit_family == "tweedie", "standard", type))
  saveRDS(fe, file.path(tr_out_dir, folder_tag, 'fit-ests-df.rds'))

  s <- purrr::map_df(f, \(file) {
    fit <- readRDS(file)
    get_sanity_df(fit, real_data = TRUE)
    }) |>
    mutate(type = gsub("poisson_link_delta", "poisson-link", type)) |>
    mutate(type = ifelse(fit_family == "tweedie", "standard", type))
  saveRDS(s, file.path(tr_out_dir, folder_tag, 'sanity-df.rds'))

  lu_df <- readRDS("data-outputs/multi-species/lu-df.rds")
  tr_f <- left_join(s, lu_df) |> pull('fname')
  tr_pred <- tr_f |>
    purrr::map(\(fit_file) {
      fit <- readRDS(file.path(tr_fit_dir, paste0(fit_file, '.rds')))
      region <- unique(fit$data$region)
      years <- unique(fit$data$year)
      sg <- syn_grid |> filter(survey %in% region)
      nd <- sdmTMB::replicate_df(
        dat = sg,
        time_name = "year",
        time_values = years
      )
      predict(fit, newdata = nd, return_tmb_object = TRUE)
  }) |>
    setNames(tr_f)
  beep()

  i <- purrr::map(tr_pred, \(p) {
    get_index(p, area = 4, bias_correct = TRUE)
  }) |>
    purrr::map2(.x = _, .y = names(tr_pred), mutate(.x, fname = gsub('\\.rds', '', .y))) |>
    bind_rows()

  i |>
    group_by(fname, id) |>
    summarise(mean_ind_cv = mean(sqrt(exp(se^2) - 1)))
  saveRDS(i, file.path(tr_out_dir,  folder_tag, "index-df.rds"))


  tr_df_dir <- file.path(tr_out_dir, folder_tag)
  # fit estimates of models that converged
  fe <- readRDS(file.path(tr_df_dir, "fit-ests-df.rds")) |>
    janitor::clean_names() |>
    mutate(type = gsub("poisson_link_delta", "poisson-link", type)) |>
    filter(type == "standard") |>
    left_join(lu_df)

  # sanity summary for all models fit (not necessarily converged)
  s <- readRDS(file.path(tr_df_dir, "sanity-df.rds")) |>
    janitor::clean_names() |>
    filter(type == "standard") |>
    left_join(lu_df)

  # indices of models that converged
  i <- readRDS(file.path(tr_df_dir, "index-df.rds")) |>
    janitor::clean_names() |>
    filter(!stringr::str_detect(fname, 'poisson-link')) |>
    select(-type) |>
    left_join(lu_df) |>
    as_tibble()

  # additional exclusion criterion is if index CV is > 1
  i |>
    group_by(fname, id) |>
    summarise(mean_ind_cv = mean(sqrt(exp(se^2) - 1)))

  # Frequency of estimated Q
  .q <- fe |>
    filter(family == "delta-gengamma") |>
    select(species, region, est_q, est_qse) |>
    arrange(est_q) |>
    mutate(q_rank = row_number())

  # ------
  tr_aic_df <- fe |>
    group_by(species, region) |>
    mutate(
      min_aic = min(aic),
      daic = aic - min_aic,
      rel_lik = exp(-0.5 * daic)
    ) |>
    mutate(aic_w = rel_lik / sum(rel_lik)) |> # see: https://atsa-es.github.io/atsa-labs/sec-uss-comparing-models-with-aic-and-model-weights.html
    ungroup() |>
    left_join(.q) |>
    arrange(q_rank) |>
    mutate(priors = TRUE)
  saveRDS(tr_aic_df, file.path(tr_df_dir, "aic-df.rds"))
}