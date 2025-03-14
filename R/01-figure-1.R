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
  mutate(sigma_text = paste0("sigma==", sigma)) |>
  mutate(q_text = paste0("Q==", Q)) |>
  mutate(family = factor(family, levels = fam))

family_colours_fig1 <- family_colours
family_colours_fig1 <- c("lognormal" = "red", "gamma" = "dodgerblue", "gengamma" = "black")
names(family_colours_fig1) <- gsub("delta-", "", names(family_colours))
family_colours_fig1[['gengamma']] <- "black"
ggplot(data = fig1_dat, aes(x = x, y = density, colour = family)) +
  # geom_line(linewidth = 1.5) +
  geom_line(data = filter(fig1_dat, family != "gengamma"), aes(linetype = family), linewidth = 1.2, alpha = 0.6) +
  geom_line(data = filter(fig1_dat, family == "gengamma"), aes(linetype = family), linewidth = 1.2) +
  facet_grid(sigma_text ~ q_text, labeller = label_parsed) +
  scale_y_continuous(limits = c(0, 0.2), breaks = c(0, 0.1, 0.2)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_colour_manual(values = family_colours_fig1,
  labels = c('Lognormal', 'Gamma', 'Generalized Gamma')
  ) +
  scale_linetype_manual(values = c(
    "lognormal" = "solid", #' ocean green'
    "gamma" =  "solid",# 'bamboo'
    "gengamma" = "22"
  ), labels = c('Lognormal', 'Gamma', 'Generalized Gamma')) +
  gfplot::theme_pbs(base_size = 12) +
  theme(#axis.title = element_blank(),
        legend.position = "bottom",
        legend.margin = margin(0.05, 0, 0, 0, 'cm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        legend.key.width = unit(1, "cm"),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")
        # strip.text.y = ggtext::element_markdown(size = 11, angle = 0),
        ) +
  tagger::tag_facets(tag = "panel", position = list(x = 0.1, y = 0.9)) +
  theme(tagger.panel.tag.text = element_text(color = "black", size = 11)) +
  ylab("Density")

ggsave(width = 6.6, height = 4.4, filename = file.path("figures", "figure-1-probability-densities.png"))
ggsave(width = 6.6, height = 4.4, filename = file.path("figures", "figure-1-probability-densities.pdf"))
