library(dplyr)
library(tidyr)
library(ggplot2)
source(here::here('R', '00-utils.R'))

sigma <- c(0.4, 1, 1.5)
Q <- c(-1, 0.01, 1)
x <- seq(0.01, 15, by = 0.01)
mu <- 1.8
fam <- c('lognormal', 'gamma', 'gengamma')

fig1_dat <- expand_grid(x, sigma, Q, mu) |>
  arrange(x, sigma, Q) |>
  mutate(gamma = dgamma(x, shape = 1 / sigma^2, rate = exp(-mu) / sigma^2),
         lognormal = dlnorm(x, meanlog = mu, sdlog = sigma),
         gengamma = flexsurv::dgengamma(x = x, mu = mu, sigma = sigma, Q = Q, log = FALSE)
  ) |>
  pivot_longer(gamma:gengamma, names_to = "family", values_to = "density") |>
  group_by(sigma, Q) |>
  mutate(panel_id = cur_group_id()) |>
  mutate(sigma_text = paste0("Sigma = ", sigma)) |>
  mutate(q_text = paste0("Q = ", Q)) |>
  mutate(family = factor(family, levels = fam))

family_colours_fig1 <- family_colours
names(family_colours_fig1) <- gsub("delta-", "", names(family_colours))
ggplot(data = fig1_dat, aes(x = x, y = density, colour = family)) +
  geom_line(aes(linetype = family), linewidth = 1.5) +
  facet_grid(sigma_text ~ q_text) +
  scale_y_continuous(limits = c(0, 0.2), breaks = c(0, 0.1, 0.2)) +
  scale_colour_manual(values = family_colours_fig1,
  labels = c('Lognormal', 'Gamma', 'Generalised Gamma')) +
  scale_linetype_manual(values = c(
    "lognormal" = "11", #' ocean green'
    "gamma" = "solid", # 'bamboo'
    "gengamma" = "2222"
  ), labels = c('Lognormal', 'Gamma', 'Generalised Gamma')) +
  gfplot::theme_pbs(base_size = 13) +
  theme(legend.position = "bottom",
        legend.margin = margin(-0.1, 0, 0, 0, 'cm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        strip.text = element_text(size = 11),
        axis.title = element_blank())

ggsave(width = 7.5, height = 5, filename = file.path("figures", "figure-1-probability-densities.png"))
