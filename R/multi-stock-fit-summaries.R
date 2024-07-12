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

fit_est_filename <- paste0("fitted-estimates", tag, ".rds")
pb <- cli::cli_progress_bar("Processing", total = length(f))
if (!file.exists(file.path(out_dir, fit_est_filename))) {
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
} else {
  message("\t Loading cached fit_ests: ", fit_est_filename)
  fit_ests <- readRDS(file.path(out_dir, fit_est_filename))
}
beep()
cli::cli_progress_done(id = pb)

pb <- cli::cli_progress_bar("Processing", total = length(f))
sanity_filename <- paste0("sanity-df", tag, ".rds")
if (!file.exists(file.path(out_dir, sanity_filename))) {
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
} else {
  message("\t Loading cached sanity-df")
  sanity_df <- readRDS(file.path(out_dir, sanity_filename))
}
beep()
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
beepr::beep()
prediction_filename <- paste0("prediction-df", tag, ".rds")
saveRDS(predictions, file.path(out_dir, prediction_filename))

# Slow did it in chunks
idx <- 201:length(predictions)
pb <- cli::cli_progress_bar("Processing", total = length(predictions[idx]))
ind201 <- purrr::map(predictions[idx], \(p) {
  cli::cli_progress_update(id = pb)
  get_index(p, area = 4, bias_correct = TRUE)
  }) |>
  purrr::map2(.x = _, .y = names(predictions[idx]), ~ mutate(.x, fname = gsub('\\.rds', '', .y)))
ind201 <- bind_rows(ind201)
cli::cli_progress_done(id = pb)
beepr::beep()

index_filename <- paste0("index-df2", tag, ".rds")
index_df <- bind_rows(ind1_200, ind201)
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
beepr::beep()
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
beepr::beep()

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
