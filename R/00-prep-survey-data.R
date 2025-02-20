library(surveyjoin) # pak::pkg_install("DFO-NOAA-Pacific/surveyjoin")
library(dplyr)

# ------------------------------------------------------------------------------
# Query GFBio (used what is in gftrends); requires VPN access
# species_list <- c("Arrowtooth Flounder", "Big Skate", "Bocaccio", "Butter Sole",
#   "Canary Rockfish", "Copper Rockfish", "Darkblotched Rockfish",
#   "Dover Sole", "English Sole", "Flathead Sole", "Greenstriped Rockfish",
#   "Kelp Greenling", "Lingcod", "Longnose Skate", "Longspine Thornyhead",
#   "North Pacific Spiny Dogfish", "Pacific Cod", "Pacific Ocean Perch",
#   "Petrale Sole", "Quillback Rockfish", "Redstripe Rockfish", "Rex Sole",
#   "Rougheye/Blackspotted Rockfish Complex", "Sablefish", "Sandpaper Skate",
#   "Sharpchin Rockfish", "Shortraker Rockfish", "Shortspine Thornyhead",
#   "Silvergray Rockfish", "Southern Rock Sole", "Splitnose Rockfish",
#   "Spotted Ratfish", "Walleye Pollock", "Widow Rockfish", "Yelloweye Rockfish",
#   "Yellowmouth Rockfish", "Yellowtail Rockfish")
# species_list <- tolower(species_list)

# all_spp_trawl <- list()
# for (i in species_list) {
#   cat(i, "\n")
#   species <- i
#   all_spp_trawl[[i]] <- try({gfdata::get_survey_sets(species, ssid = c(1, 3, 4, 16))})
# }

# survey_dat <- do.call("rbind", all_spp_trawl)
# saveRDS(survey_dat, file = "data/survey-data.rds")

# survey_index <- list()
# for (i in species_list) {
#   cat(i, "\n")
#   species <- i
#   survey_index[[i]] <- try({gfdata::get_survey_index(species, ssid = c(1, 3, 4, 16))})
# }
# survey_index_df <- bind_rows(survey_index)
# saveRDS(survey_index_df, file = "data/survey-design-index.rds")

# ------------------------------------------------------------------------------
# Clean data used for analysis:
survey_dat <- readRDS("data/survey-data.rds")

clean_dat <- survey_dat |>
  filter(stringr::str_detect(survey_abbrev, 'SYN')) |>
  filter(!(survey_abbrev == 'SYN WCHG' & year == 2014)) |>
  sdmTMB::add_utm_columns(c("longitude", "latitude"), utm_crs = 32609) |>
  mutate(
    area_swept1 = doorspread_m * (speed_mpm * duration_min),
    area_swept2 = tow_length_m * doorspread_m,
    area_swept = ifelse(!is.na(area_swept2), area_swept2, area_swept1),
    present = ifelse(catch_weight > 0, 1, 0)) |>
  mutate(offset = log(area_swept / 1e5)) |>
  rename(species = "species_common_name") |>
  select(species, species_code, region = survey_abbrev, year, depth_m, X, Y, longitude, latitude,
    present, catch_weight, offset)

hs_qcs_dat <- clean_dat |>
  filter(region %in% c("SYN HS", "SYN QCS")) |>
  mutate(region = "HS-QCS")

clean_dat <- bind_rows(clean_dat, hs_qcs_dat)

mean_pos_sets <- clean_dat |>
  group_by(species, region, year) |>
  summarise(pos = sum(present), n_sets = n()) |>
  summarise(mean_pos = mean(pos), mean_sets = mean(n_sets), .groups = 'drop') |>
  mutate(mean_pos_sets = paste0(round(mean_pos), "/", round(mean_sets)),
         prop_pos = round(mean_pos / mean_sets, digits = 2))

dat_out <- left_join(clean_dat, mean_pos_sets)

saveRDS(dat_out, here::here("data", "clean-survey-data.rds"))

# ---------------------------
# Get the Gulf of Alaska data
# ---------------------------
spp <- c(
  "shortspine thornyhead",
  "rex sole",
  # "longnose skate", # never caught in GOA
  "arrowtooth flounder",
  "dover sole",
  "petrale sole",
  "english sole",
  "spotted ratfish",
  "flathead sole",
  "pacific ocean perch",
  "lingcod",
  "pacific cod",
  "sablefish",
  "walleye pollock",
  "north pacific spiny dogfish",
  "silvergray rockfish"
)

itis_lu <- gfsynopsis::get_spp_names() |>
  filter(species_common_name %in% spp) |>
  mutate(itis_tsn = replace(itis_tsn, itis_tsn == 160617, 160620)) |>
  select(species_common_name, species_code, itis_tsn)

# Open data base connection
cache_data()
load_sql_data()

d <- get_data(itis_id = itis_lu$itis_tsn, regions = "afsc") |>
  filter(survey_name == "Gulf of Alaska") |>
  tidyr::drop_na(common_name) # I think these are caused by the itis ids that are in the look up but not in the GoA data

# Exclude the years where there were no catch since we are fitting year as iid
no_catch_years <- d |>
  group_by(common_name, itis, year) |>
  summarise(annual_catch = sum(catch_weight)) |>
  filter(annual_catch == 0)

d2 <- d |>
  anti_join(no_catch_years |> select(-annual_catch)) |>
  sdmTMB::add_utm_columns(c("lon_start", "lat_start"), utm_crs = 32609) |>
  mutate(effort_km2 = effort * 0.01) |>
  mutate(offset = log(effort_km2),
         present = ifelse(catch_weight > 0, 1, 0)) |>
  select(common_name, itis, stratum, region, year, depth_m, X, Y, lat_start, lon_start,
         present, catch_weight, effort, offset) |>
  left_join(itis_lu, by = c("itis" = "itis_tsn")) |>
  arrange(species_common_name) |>
  select(-common_name) |>
  rename(species = "species_common_name") |>
  mutate(region = "GOA")
d2

goa_mean_pos_sets <- d2 |>
  group_by(species, region, year) |>
  summarise(pos = sum(present), n_sets = n()) |>
  summarise(mean_pos = mean(pos), mean_sets = mean(n_sets), .groups = 'drop') |>
  mutate(mean_pos_sets = paste0(round(mean_pos), "/", round(mean_sets)),
         prop_pos = round(mean_pos / mean_sets, digits = 2))

d2 <- left_join(d2, goa_mean_pos_sets)

saveRDS(d2, here::here("data", "clean-afsc-data.rds"))

# Get strata areas
shp_path <- system.file("extdata/goa_strata.shp", package = "surveyjoin")
s <- sf::read_sf(shp_path) |>
  janitor::clean_names()

s_area <- s |>
  sf::st_drop_geometry() |>
  group_by(stratum) |>
  summarise(area_km2 = sum(area_km2))

# goa_design <- d2 |>
#   left_join(s_area) |>
#   mutate(cpue_tow = catch_weight / (effort * 0.01)) |> # convert to kg / km2
#   group_by(species, year, stratum, area_km2) |>
#   summarise(c_ti = sum(cpue_tow) / n(), n_strata = n()) |>
#   group_by(species, year) |>
#   mutate(b_ti = sum(c_ti * area_km2)) |>
#   distinct(species, year, b_ti) |>
#   ungroup()


# Get bootstrap estimates:
# calculate design-based biomass estimate from output of get_survey_sets()
calc_bio <- function(dat, i = seq_len(nrow(dat))) {
  dat[i, ] |>
    group_by(year, area_km2, stratum) |>
    summarise(density = mean(catch_weight / (effort * 0.01)), .groups = "drop_last") |>
    group_by(year) |>
    summarise(biomass = sum(density * area_km2), .groups = "drop_last") |>
    pull(biomass)
}

boot_one_year <- function(x, reps) {
  b <- boot::boot(x, statistic = calc_bio, strata = x$stratum, R = reps)
  suppressWarnings(bci <- boot::boot.ci(b, type = "perc"))
  tibble::tibble(
    species = unique(b$data$species),
    index = mean(b$t),
    median_boot = median(b$t),
    lwr = bci$percent[[4]],
    upr = bci$percent[[5]],
    cv = sd(b$t) / mean(b$t),
    biomass = calc_bio(x)
  )
}

boot_all_years <- function(dat, reps) {
  out <- dat |>
    split(dat$year) |>
    purrr::map_dfr(boot_one_year, reps = reps, .id = "year")
  out$year <- as.numeric(out$year)
  out
}

boot_all_years_parallel <- function(dat, reps) {
  out <- dat |>
    split(dat$year) |>
    furrr::future_map_dfr(boot_one_year, reps = reps, .id = "year",
      .options = furrr::furrr_options(seed = TRUE))
  out$year <- as.numeric(out$year)
  out
}

# Slow, so only do the three focal species
d_species <- c("north pacific spiny dogfish", "pacific ocean perch", "arrowtooth flounder")

goa_design_boot <- d2 |>
  #filter(d2, species %in% d_species) |>
  left_join(s_area) |>
  mutate(density = catch_weight / (effort * 0.01)) |> # convert to kg / km2
  group_split(species) |>
  purrr::map_df(\(x) boot_all_years_parallel(x, reps = 1000))
beepr::beep()

saveRDS(goa_design_boot, here::here("data", "goa-design.rds"))


# ------------------------------------------------------------------------------
# What proportion of grids are surveyed in a given year?
#  - used to give me a sense of sampling amount to use in the simulation
# sp_dat <- filter(dat, species_common_name == "lingcod")
# n_trawls <-
#   sp_dat |>
#     filter(stringr::str_detect(survey_abbrev, 'SYN')) |>
#     group_by(survey_abbrev, year) |>
#     summarise(n_trawls = n(), .groups = "drop")
# n_cells <-
#   gfplot::synoptic_grid |>
#     group_by(survey) |>
#     summarise(n_cells = n(), .groups = "drop") |>
#     rename(survey_abbrev = 'survey')

# left_join(n_trawls, n_cells) |>
#   mutate(prop_surveyed = 100 * (n_trawls / n_cells)) |>
#   mutate(median_surveyed = median(prop_surveyed))

# ggplot(data = _, aes(x = year, y = prop_surveyed, colour = survey_abbrev)) +
#   geom_point()
# ~ 5% of survey grid is surveyed each year
# ------------------------------------------------------------------------------
