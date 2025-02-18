library(dplyr)
library(forcats)
library(ggplot2)
library(patchwork)

theme_set(gfplot::theme_pbs(base_size = 12))
source(here::here("R", "00-utils.R"))

cv <- readRDS(here::here("data-outputs", "cv-summary-df.rds"))

family_levels <- c(
  "delta-lognormal", "delta-gamma", "Tweedie", "delta-gengamma"
)

aic_df <- readRDS(here::here("data-outputs", "multi-species", "aic-plot-df.rds")) |>
  filter(all_ok == TRUE) |>
  distinct(odd_species2, species, region, family, est_q, all_ok, aic_w) |>
  mutate(species = tolower(species)) |>
  mutate(family = factor(family, levels = family_levels)) |>
  bind_rows(
    tibble(species = c("walleye pollock", "walleye pollock"),
           region = c("GOA", "HS-QCS"),
           ))

e_df <- readRDS(here::here("data-outputs", "encounter-rates.rds"))

# qs <- aic_df |>
#   filter(family == "delta-gengamma") |>
#   mutate(.q = round(est_q, digits = 2)) |>
#   arrange(-.q) |>
#   distinct(fspecies2, .q)


ll_diff <- cv |>
  mutate(species = gsub("north ", "", species)) |>
  arrange(species, region, -loglik) |>
  group_by(species, region) |>
  mutate(max_ll = max(loglik, na.rm = TRUE), diff_ll = max_ll - loglik) |>
  ungroup() |>
  mutate(family = gsub("tweedie", "Tweedie", family)) |>
  mutate(family = factor(family, levels = family_levels)) |>
  left_join(aic_df) |>
  mutate(fregion = factor(gsub("SYN ", "", region), levels = c("GOA", "HS", "QCS", "HS-QCS", "WCVI"))) |>
  left_join(e_df)
  # filter(species != "walleye pollock")

bp <-
ggplot() +
  gfplot::theme_pbs(base_size = 11) +
  # geom_tile(
  #   data = aic_df,
  #   mapping = aes(
  #     x = 1, y = fspecies2,
  #     width = Inf, height = 1, fill = factor(odd_species2)
  #   )
  # ) +
  scale_fill_manual(values = c("grey95", "white")) +
  coord_cartesian(clip = "off") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(vjust = 0),
        axis.text = element_text(size = 9),
        axis.ticks.y = element_blank(),
        panel.spacing.x = unit(0.2, "lines"),
        panel.border = element_rect(fill = NA, linewidth = 0.5),
        plot.margin = margin(0, 0, 0, 1),
        legend.title = element_blank(),
        legend.position = "bottom"
        # legend.position = "top",
        # legend.margin = margin(3, 0, -0.3, 0, "cm")
        ) +
  guides(fill = "none")

plot_aic <- function(dat) {
  bp +
    geom_tile(
    data = dat,
    mapping = aes(
      x = 1, y = fspecies2,
      width = Inf, height = 1, fill = factor(odd_species2)
    )
  ) +
    geom_point(
      data = dat, #ll_diff |> filter(all_ok) |> filter(region == "GOA"), # models that converged
      aes(x = aic_w, y = fspecies2, colour = family, shape = family),
      stroke = 0.75, size = 1.75, position = ggstance::position_dodgev(height = 0.8)
    )  +
    scale_colour_manual(values = family_colours) +
    scale_shape_manual(values = family_shapes) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = scales::label_percent(suffix = ""),
      limits = c(0, 1)) +
      # limits = c(0, 1.27)) +
    labs(x = "AIC weight (%)", y = "Species",
      shape = "Family",
      colour = "Family") +
    guides(colour = "none", fill = "none", shape = "none") +
    facet_grid(rows = vars(fregion), switch = "y") +
    theme(strip.placement = "outside",
          strip.text.y = element_text(face = "bold"))
}

plot_cv <- function(dat) {
  bp +
    geom_tile(
      data = dat,
      mapping = aes(
        x = 1, y = fspecies2,
        width = Inf, height = 1, fill = factor(odd_species2)
    )) +
    geom_point(
      data = dat,
      aes(x = pmin(diff_ll, 500), y = fspecies2, colour = family, shape = family),
      stroke = 0.75, size = 1.75, position = ggstance::position_dodgev(height = 0.8)
    )  +
    scale_colour_manual(values = family_colours) +
    scale_shape_manual(values = family_shapes) +
    scale_x_continuous(breaks = c(400, 200, 0),
      labels = c(">400", "250", "0")) +#,
      # limits = c(0, 1)) +
    labs(x = "abs(sumLL diff)", y = "Species",
      shape = "Family",
      colour = "Family") +
    guides(colour = "none", fill = "none", shape = "none") +
    facet_grid(rows = vars(fregion), switch = "y") +
    theme(strip.text.y = element_blank())
}

plot_q <- function(dat) {
  bp +
    geom_tile(
      data = dat,
      mapping = aes(
        x = 1, y = fspecies2,
        width = Inf, height = 1, fill = factor(odd_species2)
    )) +
    geom_text(data = filter(dat, family == "delta-gengamma"),
    aes(x = 1, y = fspecies2, label = round(est_q, digits = 1)),
      hjust = 1, size = 3, colour = "grey50") +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()#,
          # axis.title.x = element_text()
          ) +
    scale_x_continuous(limits =c(0.3, 1.25)) +
    coord_cartesian(expand = FALSE) +
    guides(colour = "none", fill = "none", shape = "none") +
    facet_grid(rows = vars(fregion), switch = "y") +
    theme(strip.text.y = element_blank())
}

plot_n <- function(dat) {
  bp +
    geom_tile(
      data = dat,
      mapping = aes(
        x = 1, y = fspecies2,
        width = Inf, height = 1, fill = factor(odd_species2)
    )) +
    geom_text(data = filter(dat, family == "delta-gengamma"),
      # aes(x = 1, y = fspecies2, label = paste0("(", count, ")")),
      aes(x = 1, y = fspecies2, label = count),
      hjust = 1, size = 3, colour = "grey50") +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(colour = "white"), #element_blank(),
          axis.ticks.x = element_line(colour = "white"), # element_blank()) +
          axis.title.x = element_blank()
          ) +
    scale_x_continuous(limits =c(0.3, 1.25)) +
    coord_cartesian(expand = FALSE) +
    guides(colour = "none", fill = "none", shape = "none") +
    facet_grid(rows = vars(fregion), switch = "y") +
    theme(strip.text.y = element_blank())
}

no_x <- function() {
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()
  )
}

invis_x_axis <- function() {
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(colour = "white"),
    axis.ticks.x = element_line(colour = "white")
  )
}

get_q_rank <- function(dat, region) {
  browser()
  dat |>
    filter(family == "delta-gengamma") |>
    select(region, species, est_q) |>
    arrange(est_q) |>
    mutate(q_rank = row_number()) |>
  left_join(dat |> select(-est_q)) |>
  mutate(fspecies = forcats::fct_rev(species)) |>
  mutate(fspecies2 = forcats::fct_reorder(fspecies, as.numeric(q_rank))) |>
  mutate(odd_species2 = ifelse(as.numeric(fspecies2) %% 2 == 0, 1, 0))
}

get_legend <- function(plot) {
  gtable <- ggplotGrob(plot)
  legend <- gtable$grobs[[which(sapply(gtable$grobs, function(x) x$name) == "guide-box")]]
  return(legend)
}

legend_p <- ggplot(goa_df |> distinct(family), aes(x = 1, y = 1, color = family, shape = family)) +
  scale_colour_manual(values = family_colours) +
  scale_shape_manual(values = family_shapes) +
  geom_point(alpha = 0) +  # Invisible points
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  # theme_void() +  # Remove axes and background
  theme(plot.margin = margin(c(0, 0, 0, 0)),
        legend.margin = margin(c(0, 0, 0, 0)),
        legend.title = element_blank(),
        legend.position = "top")  # Adjust legend position

# Extract the legend
legend_grob <- get_legend(legend_p)
legend_panel <- wrap_elements(full = legend_grob)

goa_df <- ll_diff |>
  mutate(species = stringr::str_to_title(species)) |>
  filter(region == "GOA") |>
  get_q_rank(region = "GOA")
hs_df <- ll_diff |>
  mutate(species = stringr::str_to_title(species)) |>
  filter(region == "HS-QCS") |>
  get_q_rank(region = "HS-QCS")
wcvi_df <- ll_diff |>
  mutate(species = stringr::str_to_title(species)) |>
  filter(region == "SYN WCVI") |>
  get_q_rank(region == "SYN WCVI")


p1_goa <-
plot_aic(dat = goa_df) +
  no_x()
  # ggtitle("GOA") +
  # theme(plot.title.position = "plot")
p2_goa <-
plot_cv(dat = goa_df) +
  theme(axis.text.y = element_blank()) +
  no_x()
p3_goa <- plot_q(goa_df)
p4_goa <- plot_n(goa_df) +
 xlab("N")

# -- HS-QCS
p1_hs <-
plot_aic(hs_df) #+
  # ggtitle("HS-QCS") +
  # theme(plot.title.position = "plot")
p2_hs <-
plot_cv(hs_df) +
  theme(axis.text.y = element_blank())
p3_hs <- plot_q(hs_df) +
  xlab("Q")
p4_hs <- plot_n(hs_df) +
  xlab("N")

# -- WCVI
p1_wcvi <-
plot_aic(wcvi_df) +
  theme(axis.title.x = ggtext::element_markdown())
  # ggtitle("WCVI") +
  # theme(plot.title.position = "plot")
p2_wcvi <-
plot_cv(wcvi_df) +
  theme(axis.text.y = element_blank()) +
  xlab("&Delta; total CV log likelihood") +
  theme(axis.title.x = ggtext::element_markdown())
p3_wcvi <- plot_q(wcvi_df) +
  xlab("Q") +
  invis_x_axis()
  # theme(axis.text.x = element_text(colour = "black"),
  #       axis.title.x = element_text())
p4_wcvi <- plot_n(wcvi_df) +
  xlab("N") +
  theme(axis.text.x = element_text(colour = "black"),
        axis.ticks.x = element_line(colour = "black"),
        axis.title.x = element_text())

# Combination
pw1 <-
# p1_goa + p2_goa + p3_goa + p4_goa +
p1_goa + free(p2_goa, type = "label") + free(p3_goa, type = "label") +
  # free(p4_goa, type = "label") +
  # plot_layout(ncol = 4, widths = c(1, 1, 0.2, 0.2)) &
  plot_layout(ncol = 3, widths = c(1, 1, 0.2)) &
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank())
pw1 <- wrap_elements(full = pw1)

pw2 <-
# p1_hs + p2_hs + p3_hs + p4_hs +
p1_hs + free(p2_hs, type = "label") + free(p3_hs, type = "label") +
  # free(p4_hs, type = "label") +
  # plot_layout(ncol = 4, widths = c(1, 1, 0.2, 0.2)) &
  plot_layout(ncol = 3, widths = c(1, 1, 0.2)) &
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank())
pw2 <- wrap_elements(full = pw2)

pw3 <-
p1_wcvi + free(p2_wcvi, type = "label") + free(p3_wcvi, type = "label") +
  # free(p4_wcvi, type = "label") +
  # plot_layout(ncol = 4, widths = c(1, 1, 0.2, 0.2)) &
  plot_layout(ncol = 3, widths = c(1, 1, 0.2)) &
  theme(legend.position = "none")
pw3 <- wrap_elements(full = pw3)

legend_row <- plot_spacer() + legend_panel + plot_spacer() +
  plot_layout(ncol = 3, widths = c(0.4, 1, 0.1))
  # plot_layout(ncol = 3, widths = c(0.3, 1, 0.2))
legend_row <- wrap_elements(full = legend_row)

legend_row +
  pw1 +
  plot_spacer() +
  pw2 +
  plot_spacer() +
  pw3 +
  plot_spacer() +
  plot_layout(heights = c(0.001, 1, -0.116, 1, -0.117, 1.25, -0.18))

# ggsave(here::here("figures", "figure-4-option2-annual-encounter.png"), width = 7.86, height = 11.9)
ggsave(here::here("figures", "figure-4-option2.png"), width = 6.7, height = 9)

# ----