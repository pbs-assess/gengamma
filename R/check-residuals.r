library(sdmTMB)
library(dplyr)

#source(here::here('R', '00-utils.R'))

# fits <- readRDS(here::here('data-outputs', 'fits', '1-cv0.8-sigmao0.6.rds'))

# fit_df <- purrr::map_dfr(fits, get_fitted_estimates) |>
#   mutate(id = row_number())

# idx <- fit_df |>
#   filter(sim_family == "delta-gengamma", Q == 0.8, fit_family == "gengamma") |>
#   pull(id)
# test_fit <- fits[[idx]]

#saveRDS(test_fit, file.path(here::here('data-outputs', 'fits', 'test-fit.rds')))

test_fit <- readRDS(file.path(here::here('data-outputs', 'fits', 'test-fit.rds')))


seed <- sample.int(1e+06, 1)
print(seed)
set.seed(seed)
s <- simulate(test_fit, nsim = 200, type = 'mle-mvn')

stop()

source(here::here('R/check-residuals.r'))

r <- dharma_residuals(s, test_fit, plot = FALSE)
head(r)
plot(r$expected, r$observed)
abline(0, 1)

set.seed(1)
s <- simulate(test_fit, nsim = 200, type = 'mle-mvn')


set.seed(964601)
s <- simulate(test_fit, nsim = 200, type = 'mle-mvn')

set.seed(378612)
set.seed(970335)
s <- simulate(test_fit, nsim = 200, type = 'mle-mvn')

# Try Tweedie family:
fit <- sdmTMB(density ~ as.factor(year) + s(depth, k = 3),
  data = pcod_2011, mesh = pcod_mesh_2011,
  family = tweedie(link = "log"), spatial = "on")

# The `simulated_response` argument is first so the output from
# simulate() can be piped to dharma_residuals():

# not great:
simulate(fit, nsim = 200, type = "mle-mvn") |>
  dharma_residuals(fit)

# delta-lognormal looks better:
fit_dl <- update(fit, family = delta_lognormal())
simulate(fit_dl, nsim = 200, type = "mle-mvn") |>
  dharma_residuals(fit)

# or skip the pipe:
s <- simulate(fit_dl, nsim = 200, type = "mle-mvn")
# and manually plot it:
r <- dharma_residuals(s, fit_dl, plot = FALSE)
head(r)
plot(r$expected, r$observed)
abline(0, 1)

# return the DHARMa object and work with the DHARMa methods
ret <- simulate(fit_dl, nsim = 200, type = "mle-mvn") |>
  dharma_residuals(fit, return_DHARMa = TRUE)
plot(ret)

# Try gengamma
fit <- sdmTMB(density ~ as.factor(year) + s(depth, k = 3),
  data = pcod_2011 |> filter(density > 0),
  mesh = make_mesh(pcod_2011 |> filter(density > 0), xy_cols = c('X', 'Y'), cutoff = 8),
  family = gengamma(link = "log"), spatial = "on")

# not great:
simulate(fit, nsim = 200, type = "mle-mvn") |>
  dharma_residuals(fit)

