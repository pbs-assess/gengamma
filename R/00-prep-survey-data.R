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
  select(species, species_code, survey_abbrev, year, depth_m, X, Y, longitude, latitude,
    present, catch_weight, density_kgpm2, offset)

mean_pos_sets <- clean_dat |>
  group_by(species, survey_abbrev, year) |>
  summarise(pos = sum(present), n_sets = n()) |>
  summarise(mean_pos = mean(pos), mean_sets = mean(n_sets), .groups = 'drop') |>
  mutate(mean_pos_sets = paste0(round(mean_pos), "/", round(mean_sets)),
         prop_pos = round(mean_pos / mean_sets, digits = 2))

dat_out <- left_join(clean_dat, mean_pos_sets)

saveRDS(dat_out, file.path("data", "clean-survey-data.rds"))

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
