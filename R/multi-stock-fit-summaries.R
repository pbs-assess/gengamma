library(dplyr)
library(sdmTMB)

source("R/00-utils.R")

out_dir <- here::here("data-outputs", "bcgf-outputs")
fit_dir <- here::here("data-outputs", "bcgf-outputs", "fits")

lu_df <- readRDS(file.path(out_dir, "lu-df.rds")) |>
  janitor::clean_names()

syn_grid <- as_tibble(gfplot::synoptic_grid) |>
  select(survey, X, Y, depth)


# ------------------------------------------------------------------------------
# Examine outputs
# ------------------------------------------------------------------------------
f <- list.files(file.path(out_dir, "fits"), full.names = TRUE)
tag <- "-all"
overwrite_outputs <- TRUE

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
  left_join(lu_df) |>
  filter(all_ok) |>
  mutate(fname = paste0(fname, '.rds'))

# Only predict and get index for fits that passed a sanity check

ff <- ok_sanity$fname
pb <- cli::cli_progress_bar("Processing", total = length(ff))
predictions <- purrr::map(ff, \(fit_file) {
  cli::cli_progress_update(id = pb)
  fit <- readRDS(file.path(fit_dir, fit_file))
  region <- unique(fit$data$survey_abbrev)
  years <- unique(fit$data$year)
  sg <- syn_grid |> filter(survey %in% region)
  nd <- sdmTMB::replicate_df(
    dat = sg,
    time_name = "year",
    time_values = years
  )
  predict(fit, newdata = nd, return_tmb_object = TRUE)
}) |>
  setNames(ff)
cli::cli_progress_done(id = pb)
beep()

# Slow did it in chunks
idx <- 1:length(predictions)
#idx <- 201:length(predictions)
pb <- cli::cli_progress_bar("Processing", total = length(predictions[idx]))
ind <- purrr::map(predictions[idx], \(p) {
  cli::cli_progress_update(id = pb)
  get_index(p, area = 4, bias_correct = TRUE)
  }) |>
  purrr::map2(.x = _, .y = names(predictions[idx]), ~ mutate(.x, fname = gsub('\\.rds', '', .y)))
ind <- bind_rows(ind)
cli::cli_progress_done(id = pb)
beep()

index_filename <- paste0("index-df", tag, ".rds")
index_df <- bind_rows(ind)
saveRDS(index_df, file.path(out_dir, index_filename))


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

# DHARMA worth looking at? Transform to normal?
ff <- ok_sanity |>
  filter(fname %in% ff[grepl("lingcod-.*-.*", ff)]) |>
  pull("fname")
pb <- cli::cli_progress_bar("Processing", total = length(ff))
test <- ff |>
  purrr::map_dfr(\(fit_file) {
  cli::cli_progress_update(id = pb)
    fit <- readRDS(file.path(fit_dir, fit_file))
    fname <- gsub("\\.rds", "", fit_file)
    get_dr(fit, fit_file, seed = 100, nsim = 200)
  }) |>
  tibble::as_tibble() |>
  rename(fname = "id")
cli::cli_progress_done(id = pb)
beep()

test_f <- readRDS(file.path(fit_dir, ff[1]))
test_s <- simulate(test_f, nsim = 200, type = "mle-mvn")
dr <- dharma_residuals(test_s, test_f)

dr <- dharma_residuals(test_s, test_f, plot = T)
str(dr)


get_dr <- function(fit_obj, fit_id, nsim = 200, seed = sample.int(1e6, 1), type = "mle-mvn", ...) {
  set.seed(seed)
  simulate(fit_obj, nsim = nsim, type = type, ...) |>
    dharma_residuals(fit_obj, plot = FALSE, ...) |>
    mutate(id = fit_id, seed = seed)
}

DHARMa::testOutliers()
temp = testOutliers(simulationOutput, plot = F)
        legend("bottomright", c(paste("Outlier test: p=", round(temp$p.value,
            digits = 5)), paste("Deviation ", ifelse(temp$p.value <
            0.05, "significant", "n.s."))), text.col = ifelse(temp$p.value <
            0.05, "red", "black"), bty = "n")

test |> right_join(lu_df) |>
ggplot(data = _, aes(x = expected, y = observed)) +
  geom_point(shape = 21) +
  geom_abline(intercept = 0, slope = 1) +
  facet_grid(fit_family ~ region) +
  guides(colour = "none")

list.files(file.path(out_dir, "fits"))
test <- readRDS(file.path(out_dir, "fits", "arrowtooth-flounder-SYN HS-tweedie.rds"))
test2 <- get_fitted_estimates(test, real_data = TRUE)



# ------------------------------------------------------------------------------
# Look at the cases where we tried using a prior
# ------------------------------------------------------------------------------
tr_fit_dir <- file.path(fit_dir, 'priors')
dir.create(file.path(out_dir, 'priors'), showWarnings = FALSE)
tr_out_dir <- file.path(out_dir, 'priors')
f <- list.files(tr_fit_dir, full.names = TRUE)

fe <- purrr::map_df(f, \(file) {
  fit <- readRDS(file)
  get_fitted_estimates(fit, real_data = TRUE)
  }) |>
  mutate(type = gsub("poisson_link_delta", "poisson-link", type)) |>
  mutate(type = ifelse(fit_family == "tweedie", "standard", type))
saveRDS(fe, file.path(out_dir, 'priors', 'fit-ests-df.rds'))

s <- purrr::map_df(f, \(file) {
  fit <- readRDS(file)
  get_sanity_df(fit, real_data = TRUE)
  }) |>
  mutate(type = gsub("poisson_link_delta", "poisson-link", type)) |>
  mutate(type = ifelse(fit_family == "tweedie", "standard", type))
saveRDS(s, file.path(out_dir, 'priors', 'sanity-df.rds'))

#sable_gg <- readRDS(file.path(tr_fit_dir, "sablefish-SYN HS-delta-gengamma-poisson-link.rds"))

tr_f <- left_join(s, lu_df) |> pull('fname')
tr_pred <- tr_f |>
  purrr::map(\(fit_file) {
    fit <- readRDS(file.path(tr_fit_dir, paste0(fit_file, '.rds')))
    region <- unique(fit$data$survey_abbrev)
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
  tr_fit_dir <- file.path("data-outputs/bcgf-outputs/fits", folder_tag)
  dir.create(file.path("data-outputs/bcgf-outputs/", folder_tag), showWarnings = FALSE)
  tr_out_dir <- file.path("data-outputs/bcgf-outputs/", folder_tag)
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

  lu_df <- readRDS("data-outputs/bcgf-outputs/lu-df.rds")
  tr_f <- left_join(s, lu_df) |> pull('fname')
  tr_pred <- tr_f |>
    purrr::map(\(fit_file) {
      fit <- readRDS(file.path(tr_fit_dir, paste0(fit_file, '.rds')))
      region <- unique(fit$data$survey_abbrev)
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