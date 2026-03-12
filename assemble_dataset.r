library(dplyr)
library(tidyr)
library(sf)
library(tigris)

screens <- read.csv("Screens # per 100K Table.csv", sep = "\t", fileEncoding = "UTF-16LE") %>%
  rename(zip = 2, ptsd_rate = X..per.100K) %>%
  mutate(zip = as.character(zip), ptsd_rate = ptsd_rate / 100000)

zip_pop <- read.csv("Chicago_Population_Counts.csv") %>%
  filter(Year == 2020) %>%
  transmute(zip = as.character(Geography), population = Population...Total)

covars_raw <- read.csv("Chicago Health Atlas Data Download.csv")
covars <- covars_raw[c(1, 5:62), c(2, 5, 7, 9, 11, 13, 15, 17)] %>%
  setNames(c("zip", "med_income", "uninsured_rate", "pm25", "hs_grad_rate", "unemployment_rate", "hardship_index", "prop_white")) %>%
  slice(-1) %>%
  mutate(across(-zip, as.numeric), zip = as.character(zip))

violent_types <- c("HOMICIDE", "ASSAULT", "BATTERY", "ROBBERY", "CRIMINAL SEXUAL ASSAULT")
zip_shapes <- zctas(cb = TRUE, year = 2020, progress_bar = FALSE) %>% st_transform(4326)

spatial_data <- read.csv("Crimes_-_2001_to_Present_20260307.csv") %>%
  filter(!is.na(Latitude), !is.na(Longitude), Primary.Type %in% violent_types) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_join(zip_shapes %>% select(zip = ZCTA5CE20)) %>%
  st_drop_geometry() %>%
  count(zip, Primary.Type) %>%
  pivot_wider(names_from = Primary.Type, values_from = n, values_fill = 0) %>%
  inner_join(zip_pop, by = "zip") %>%
  inner_join(screens, by = "zip") %>%
  inner_join(covars, by = "zip") %>%
  mutate(
    violent_rate = rowSums(across(all_of(violent_types))) / population,
    ptsd_cases = round(ptsd_rate * population),
    log_violent = log(violent_rate + 1)
  ) %>%
  inner_join(zip_shapes %>% select(zip = ZCTA5CE20), by = "zip") %>%
  st_as_sf()

coords <- st_coordinates(st_centroid(spatial_data))
spatial_data$lat <- coords[,2]
spatial_data$long <- coords[,1]
spatial_data <- as.data.frame(spatial_data)

spatial_data_scaled <- spatial_data %>%
  mutate(across(
    c(med_income, uninsured_rate, pm25, hs_grad_rate, 
      unemployment_rate, hardship_index, prop_white, 
      log_violent, lat, long), 
    ~ as.numeric(scale(.)), 
    .names = "{.col}_s"
  ))

summary(spatial_data_scaled %>% select(contains("_s")))
