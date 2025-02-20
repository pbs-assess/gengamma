# What about generating a table for the supplement that shows the encounter rates
# and average encounter tows per year in the various surveys and then dropping
# petrale from the GOA panel and justifying that as being below 5%?
library(dplyr)
library(xtable)

dat0 <- readRDS(here::here("data-outputs", "data-used.rds")) |>
  filter(species != "longnose skate") |>
  mutate(species = gsub("north ", "", species))

dat <- dat0 |>
  group_by(region, species, year) |>
  summarise(e_rate = sum(present) / n(), e_count = sum(present)) |>
  group_by(species, region) |>
  summarise(rate = round(mean(e_rate), digits = 2),
    count = round(mean(e_count), digits = 0),
    .groups = "drop") |>
  arrange(region, -rate) |>
  relocate(region, species)
saveRDS(dat, here::here("data-outputs", "encounter-rates.rds"))

dat <- dat |>
  mutate(region = factor(region)) |>
  mutate(region = forcats::fct_recode(region, "WCVI" = "SYN WCVI")) |>
  mutate(species = stringr::str_to_title(species)) |>
  mutate(region = as.character(region)) |>
  mutate(region = ifelse(duplicated(region), "", region)) |>
  mutate(count = as.integer(count)) # to get rid of decimals in table

# Horizontal line breaks
line_pos <- c(0, 14, 28, 42)

# Convert to xtable
xt <- xtable(dat, align = c("l", "l", "l", "r", "r"))

# Output to use in overleaf
print(xt, include.rownames = FALSE, sanitize.text.function = identity,
  digits = c(0, 0, 2, 0, 0),
  add.to.row = list(
    pos = as.list(line_pos),  # Hard-coded positions for \hline
    command = rep("\\hline\n", length(line_pos))  # Insert horizontal lines
  )
)
