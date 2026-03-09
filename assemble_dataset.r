library(dplyr)
library(tidyr)
library(sf)
library(tigris)

screens <- read.csv(
  "Screens # per 100K Table.csv",
  row.names = NULL,
  sep = "\t",
  fileEncoding = "UTF-16LE"
)
names(screens)[2] <- "zip"
zip_pop <- read.csv("Chicago_Population_Counts.csv")
crimes <- read.csv("Crimes_-_2001_to_Present_20260307.csv")
covars <- read.csv("Chicago Health Atlas Data Download.csv")

crime_data <- crimes %>%
  select(
    Primary.Type,
    Location.Description,
    Arrest,
    Domestic,
    Latitude,
    Longitude
  ) %>%
  filter(!is.na(Latitude), !is.na(Longitude)) %>%
  mutate(
    Primary.Type = as.factor(Primary.Type),
    Location.Description = as.factor(Location.Description),
    Arrest = as.factor(Arrest),
    Domestic = as.factor(Domestic)
  )

summary(crime_data)

zip_shapes <- zctas(cb = TRUE, year = 2020) %>%
  st_transform(4326) %>%
  select(zip = ZCTA5CE20)

crime_sf <- st_as_sf(
  crime_data,
  coords = c("Longitude", "Latitude"),
  crs = 4326,
  remove = FALSE
)

crime_zip <- st_join(crime_sf, zip_shapes)
crime_zip$zip <- as.character(crime_zip$zip)
sort(table(crime_zip$zip))

covars <- data.frame(covars[c(1,5:62), 2],
                     covars[c(1,5:62), 5],
                     covars[c(1,5:62), 7],
                     covars[c(1,5:62), 9],
                     covars[c(1,5:62), 11],
                     covars[c(1,5:62), 13],
                     covars[c(1,5:62), 15],
                     covars[c(1,5:62), 17])

names(covars) <- c("zip", "med_income", "uninsured_rate", "pm25", "hs_grad_rate", "unemployment_rate", "hardship_index", "prop_white")
covars <- covars[2:59, ]
covars <- as.data.frame(apply(covars, MARGIN = 2, FUN = as.numeric))
covars$zip <- as.factor(covars$zip)

violent_types <- c(
  "HOMICIDE",
  "ASSAULT",
  "BATTERY",
  "ROBBERY",
  "CRIMINAL SEXUAL ASSAULT"
)

zip_crime_summary <- crime_zip %>%
  filter(Primary.Type %in% violent_types) %>%
  st_drop_geometry() %>%
  count(zip, Primary.Type) %>%
  pivot_wider(
    names_from = Primary.Type,
    values_from = n,
    values_fill = 0
  )

zip_population <- zip_pop %>%
  filter(Year == 2020) %>%
  select(
    zip = Geography,
    population = Population...Total
  ) %>%
  slice(2:55) %>%
  mutate(zip = as.character(zip))

zip_crime_summary <- left_join(
  zip_crime_summary,
  zip_population,
  by = "zip"
)

zip_crime_summary <- zip_crime_summary %>%
  mutate(
    violent_rate = rowSums(across(all_of(violent_types))) / population
  )

ptsd_data <- screens %>%
  select(
    zip,
    ptsd_rate = "X..per.100K"
  ) %>%
  mutate(zip = as.character(zip)) %>%
  mutate(ptsd_rate = ptsd_rate/100000)

chicago_crime_ptsd <- zip_crime_summary %>%
  left_join(ptsd_data, by = "zip") %>%
  na.omit()

chicago_crime_ptsd <- chicago_crime_ptsd %>%
  left_join(covars, by = "zip") %>%
  na.omit()

names(chicago_crime_ptsd)
chicago_crime_ptsd
