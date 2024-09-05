# ------------------------------------------------------------------------------
# Plot map of survey regions

library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(sf)
library(patchwork)

dat <- readRDS(here::here("data-outputs", "data-used.rds")) |>
  filter(region %in% c("HS-QCS", "SYN WCVI", "GOA")) |>
  mutate(region = factor(region)) |>
  mutate(region = forcats::fct_recode(region, "WCVI" = "SYN WCVI"))

map_dat <- dat |>
  distinct(region, longitude, latitude) |>
  st_as_sf(coords = c("longitude", "latitude"), crs = "WGS84")

map_bbox <- st_bbox(c(xmin = -172, xmax = -120, ymin = 45, ymax = 63), crs = st_crs(4326))

coast_poly <- ne_countries(country = c("canada", "united states of america"), scale = 50) |>
  select(sovereignt) |>
  st_crop(map_bbox)

# region_colours <- rev(c("#fdcc8a", "#fc8d59", "#d7301f"))
# region_colours <- rev(c("#a1dab4", "#41b6c4", "#225ea8"))
region_colours <- rev(c("#fbb4b9", "#f768a1", "#ae017e"))

sampling_plot <- dat |>
  filter(species == "arrowtooth flounder") |>
  group_by(region, year) |>
  summarise(samples = n()) |>
  arrange(region, year) |>
ggplot(data = _, aes(x = year, y = samples, colour = region, shape = region)) +
  geom_point(size = 3) +
  scale_colour_manual(values = region_colours) +
  # viridis::scale_colour_viridis(discrete = TRUE) +
  guides(colour = "none", shape = "none") +
  labs(x = "Year", y = "Samples per year", title = "Annual data availability") +
  gfplot::theme_pbs(base_size = 12)

survey_map <- ggplot() +
  geom_sf(data = coast_poly, colour = "grey10", fill = "grey90") +
  geom_sf(data = map_dat, aes(colour = region, shape = region)) +
  scale_colour_manual(values = region_colours) +
  # viridis::scale_colour_viridis(discrete = TRUE) +
  coord_sf(ylim = c(46, 62), xlim = c(-172, -122), expand = FALSE) +
  scale_y_continuous(breaks = c(46, 50, 54, 58, 62)) +
  gfplot::theme_pbs(base_size = 12) +
  theme(legend.position = c(0.5, 0.1),
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.key.spacing.x = unit(1.5, "lines")) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  labs(x = "Longitude", y = "Latitude", title = "Data spatial coverage")
#survey_map

(wrap_elements(full = survey_map) / sampling_plot) +
  plot_layout(height = c(2, 1.5), width = c(1, 1)) +
  plot_annotation(tag_levels = "a", tag_suffix = ")") +
  theme(plot.tag.position = c(0, 0.99))

dir.create(here::here("figures", "supp"), showWarnings = FALSE, recursive = TRUE)
ggsave(width = 6, height = 6.8, filename = here::here("figures", "supp", "sampling-diff.png"))
