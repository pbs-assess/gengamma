library(dplyr)
# library(sdmTMB)
devtools::load_all("../sdmTMB") # gengamma RQR branch - not needed for dogfish test fitting
library(ggplot2)
theme_set(gfplot::theme_pbs())

source("R/00-utils.R")

# Load data
dat <- readRDS(file.path("data-outputs", "clean-survey-data.rds"))
spp_list <- unique(dat$species)
fig_dir <- here::here(file.path("figures", "bc-gf-data", "exploratory"))


# DOGFISH
# ------------------------------------------------------------------------------
test_dat <- dat |>
  filter(
    species == "north pacific spiny dogfish",
    survey_abbrev == "SYN WCVI"
  )

mesh <- make_mesh(test_dat, xy_cols = c("X", "Y"), cutoff = 8)

# catch_weight * 1e-2 works for spatial = 'on', spatiotemporal = 'iid'
# Also works without scaling and spatial = 'off', spatiotemporal = 'off'

fit1 <- sdmTMB(
  data = test_dat, #|> mutate(catch_weight = catch_weight * 1e-2),
  formula = as.formula(catch_weight ~ 0 + as.factor(year)),
  mesh = mesh,
  offset = "offset",
  spatial = "on",
  spatiotemporal = "off",
  time = "year",
  family = sdmTMB::delta_gengamma()
)

# make catch weight slightly smaller; converges:
fit2 <- update(fit1, data = test_dat |> mutate(catch_weight = catch_weight * 0.15))

# add light prior on year effects; converges:
fit3 <- update(fit1, priors = sdmTMBpriors(b = normal(rep(0, 10), rep(10, 10))))

# add light prior on year effects; converges:
fit3 <- update(
  fit1,
  priors = sdmTMBpriors(b = normal(rep(0, 10), rep(100, 10))),
  data = test_dat
)

fit4 <- update(
  fit1,
  control = sdmTMBcontrol(profile = TRUE),
  data = test_dat
)

# my local version only, this works too:
# move gengamma_Q only to inner optimization; converges:
PROFILE <- "gengamma_Q"
fit5 <- update(
  fit1,
  control = sdmTMBcontrol(profile = TRUE),
  data = test_dat
)

p <- get_pars(fit5)
p$gengamma_Q

# ----------


sp_f <- f[grepl("north-pacific-spiny-dogfish", f)]
fa <- lapply(sp_f, readRDS) |>
  setNames(stringr::str_extract(sp_f, (".*fits/(.*)\\.rds"), group = 1))
sp_fits <- purrr::keep(fa, ~ inherits(.x, "sdmTMB"))
test <- sp_fits |>
  purrr::map_dfr(get_sanity_df, real_data = TRUE, .gradient_thresh = 0.005)

g <- sp_fits |> purrr::map(~ .x$gradients)


sp_fits[["north-pacific-spiny-dogfish-SYN WCVI-delta-gengamma"]] |> names()

rqr_df <- purrr::map_dfr(1:length(sp_fits), ~ get_rqr(sp_fits[[.x]], id = names(sp_fits[.x])))

rqr_df |>
  mutate(region = stringr::str_extract_all(id, "(.*)-(SYN [A-Z]{2,4})-(.*)"))

rqr_df |>
  mutate(
    species = stringr::str_extract(id, ".*(?=-SYN)"),
    region = stringr::str_extract(id, ".*-(SYN [A-Z]{2,4})-.*", group = 1),
    family = stringr::str_extract(id, ".*-SYN [A-Z]{2,4}-(.*)", group = 1)
  ) |>
  ggplot(aes(sample = r)) +
  geom_qq() +
  geom_abline(intercept = 0, slope = 1) +
  facet_grid(region ~ family)
# facet_wrap(~ family)
ggsave(filename = file.path(fig_dir, "dogfish-rqr.png"))


filter(dat, species == "north pacific spiny dogfish", survey_abbrev == "SYN WCVI") |>
  ggplot() +
  geom_histogram(aes(x = catch_weight)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
ggsave(filename = file.path(fig_dir, "dogfish-catch-hist.png"))


regions <- c("SYN WCVI", "SYN QCS", "SYN HS", "SYN WCHG")
families <- c(
  "tweedie",
  "delta-gamma", "delta-lognormal", "delta-gengamma"
)
# "delta-plink-lognormal", "delta-plink-gamma", "delta-plink-gengamma")
tofit <- tidyr::expand_grid(
  .region = regions, .family = families,
  .species = spp_list
) |>
  mutate(.id = row_number())

species <- "north pacific spiny dogfish"
scaler <- 1e-2 # works
region <- "SYN WCVI"

sp_file <- paste0(clean_name(species), ".rds")
cutoff <- 8

test_fits <-
  tofit |>
  filter(.species == species) |>
  filter(.region == region) |>
  # filter(.family == 'delta-gengamma') |>
  purrr::pmap(\(.region, .family, .id, .species) {
    get_fit(
      survey_dat = dat |> mutate(catch_weight = catch_weight * scaler),
      formula = as.formula(catch_weight ~ 0 + as.factor(year)),
      offset = "offset",
      region = .region, family = .family, species = .species,
      cutoff = cutoff, sp = "on", st = "off"
    )
  })
beep()

fitted_ests <- purrr::map_dfr(test_fits, get_fitted_estimates, real_data = TRUE)

fitted_ests |>
  select(fit_family, est_Q, aic) |>
  flextable::flextable()



# species <- 'walleye pollock'
# region <- 'SYN WCVI'
# scaler <- 1e-3 # QCS needs 1e-1; nothing seems to fix SYN WCVI

filter(dat, species == "walleye pollock", survey_abbrev == "SYN WCVI") |>
  pull(catch_weight) |>
  hist()
