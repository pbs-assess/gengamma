library(dplyr)
library(forcats)
library(ggplot2)
library(patchwork)

theme_set(gfplot::theme_pbs(base_size = 12))
source(here::here("R", "00-utils.R"))
source(here::here("R", "00-cv-plotting-functions.R"))

# k <- 25
k <- 50

# sample_size <- TRUE
sample_size <- FALSE

filename <- paste0("figure-4-", k, "-fold", "-ss-", sample_size, ".png")

if (k == 25) {
  cv <- lapply(list.files(here::here("data-outputs", "cv"), full.names = TRUE), readRDS) |>
    bind_rows() |>
    filter(species != "longnose skate") |>
    rename(sumloglik = "loglik")
}

if (k == 50) {
  cv <- readRDS(here::here("data-outputs", "cv-summary-df-50.rds")) |>
    tidyr::separate_rows(foldloglik, sep = ";", convert = TRUE) |>
    filter(species != "longnose skate")
}

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


ll_diff <- cv |>
  distinct(family, region, species, sumloglik) |>
  mutate(species = gsub("north ", "", species)) |>
  arrange(species, region, -sumloglik) |>
  group_by(species, region) |>
  mutate(max_ll = max(sumloglik, na.rm = TRUE), diff_ll = max_ll - sumloglik) |>
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
        ) +
  guides(fill = "none")


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

centre_min <- FALSE
{
p1_goa <-
plot_aic(dat = goa_df) +
  no_x()
  # ggtitle("GOA") +
  # theme(plot.title.position = "plot")
p2_goa <-
plot_cv(dat = goa_df, centre_min) +
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
plot_cv(hs_df, centre_min) +
  theme(axis.text.y = element_blank())
p3_hs <- plot_q(hs_df) +
  xlab("Q")
p4_hs <- plot_n(hs_df) +
  xlab("N")

# -- WCVI
p1_wcvi <-
plot_aic(wcvi_df)
  # ggtitle("WCVI") +
  # theme(plot.title.position = "plot")
p2_wcvi <-
plot_cv(wcvi_df, centre_min) +
  theme(axis.text.y = element_blank()) +
  xlab(parse(text = "Delta~log(italic(L))~of~`out-of-sample`"))
p3_wcvi <- plot_q(wcvi_df) +
  xlab(parse(text = "Q")) +
  theme(axis.title.x = element_text(vjust = -1.5)) +
  invis_x_axis()
  # theme(axis.text.x = element_text(colour = "black"),
  #       axis.title.x = element_text())
p4_wcvi <- plot_n(wcvi_df) +
  xlab("N") +
  theme(axis.text.x = element_text(colour = "white"),
        axis.ticks.x = element_line(colour = "white"),
        axis.title.x = element_text())

# Combination
pw1 <-
p1_goa + free(p2_goa, type = "label") + free(p3_goa, type = "label") +
  plot_layout(ncol = 3, widths = c(1, 1, 0.2)) &
  # free(p4_goa, type = "label") +
  # plot_layout(ncol = 4, widths = c(1, 1, 0.2, 0.2)) &
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank())
pw1 <- wrap_elements(full = pw1)

pw2 <-
p1_hs + free(p2_hs, type = "label") + free(p3_hs, type = "label") +
  plot_layout(ncol = 3, widths = c(1, 1, 0.2)) &
  # free(p4_hs, type = "label") +
  # plot_layout(ncol = 4, widths = c(1, 1, 0.2, 0.2)) &
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank())
pw2 <- wrap_elements(full = pw2)

pw3 <-
p1_wcvi + free(p2_wcvi, type = "label") + free(p3_wcvi, type = "label") +
  plot_layout(ncol = 3, widths = c(1, 1, 0.2)) &
  # free(p4_wcvi, type = "label") +
  # plot_layout(ncol = 4, widths = c(1, 1, 0.2, 0.2)) &
  theme(legend.position = "none")
pw3 <- wrap_elements(full = pw3)

if (sample_size) {
pw1 <-
p1_goa + free(p2_goa, type = "label") + free(p3_goa, type = "label") +
  free(p4_goa, type = "label") +
  plot_layout(ncol = 4, widths = c(1, 1, 0.2, 0.2)) &
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank())
pw1 <- wrap_elements(full = pw1)

pw2 <-
p1_hs + free(p2_hs, type = "label") + free(p3_hs, type = "label") +
  free(p4_hs, type = "label") +
  plot_layout(ncol = 4, widths = c(1, 1, 0.2, 0.2)) &
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank())
pw2 <- wrap_elements(full = pw2)

pw3 <-
p1_wcvi + free(p2_wcvi, type = "label") + free(p3_wcvi, type = "label") +
  free(p4_wcvi, type = "label") +
  plot_layout(ncol = 4, widths = c(1, 1, 0.2, 0.2)) &
  theme(legend.position = "none")
pw3 <- wrap_elements(full = pw3)
}

legend_row <- plot_spacer() + legend_panel + plot_spacer() +
  plot_layout(ncol = 3, widths = c(0.4, 1, 0.1))
legend_row <- wrap_elements(full = legend_row)

legend_row +
  pw1 +
  plot_spacer() +
  pw2 +
  plot_spacer() +
  pw3 +
  plot_spacer() +
  plot_layout(heights = c(0.001, 1, -0.116, 1, -0.117, 1.25, -0.18))
}

# ggsave(here::here("figures", "figure-4-option2-annual-encounter.png"), width = 7.86, height = 11.9)
# ggsave(here::here("figures", filename), width = 6.7, height = 9)
ggsave(here::here("figures", "figure-4-50-fold.png"), width = 6.7, height = 9)
ggsave(here::here("figures", "figure-4-50-fold.pdf"), width = 6.7, height = 9)
# if (k == 50) ggsave(here::here("figures", "figure-4-50-fold.png"), width = 6.7, height = 9)


# ----

sds <- cv |>
  filter(!is.na(foldloglik)) |>
  group_by(family, region, species) |>
  summarise(mean_ll = round(mean(foldloglik), 1), sd_ll = round(sd(folIdloglik) / sqrt(n()), 1), k = n()) |>
  arrange(region, species, -mean_ll) |>
  filter(region != "GOA")
  #view_dt()

ll_diff |>
  filter(region != "GOA") |>
  arrange(region, species, diff_ll) |>
  left_join(sds) |>
  filter(family == "delta-gengamma") |>
  select(region, species, diff_ll, sd_ll) |>
  mutate(diff_ll = round(diff_ll, 1)) |>
  arrange(region, diff_ll) |>
  print(n = 27)

