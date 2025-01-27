
# libraries & working directory -------------------------------------------
## load libraries 
{
library(dplyr)
library(lubridate)
library(geosphere)  # For distance calculations (if needed)
library(ggplot2)
library(RNCEP)
library(maps)
library(leaflet)
library(geosphere)
library(RColorBrewer)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(foreach)
library(doParallel)
 }
# set working directory
setwd("/Users/tom/Desktop/Research_project_2/Track_data_montagu_Harrier/Data on interruptions migrating harrier/RAW data individual tracks")

# remove list 
rm(list=ls())

# Loading Data ------------------------------------------------------------
# load data
data <- read.csv("Data.All.Corry_autumn_2012.csv.NEW.csv") 

# Step 1: Filter for periods with at least 5 consecutive "fly"
fly_data <- data %>%
  mutate(fly_group = cumsum(speed.class != "fly")) %>%  # Group consecutive "fly" rows
  group_by(fly_group) %>%
  filter(speed.class == "fly" & n() >= 5) %>%          # Keep groups with at least 5 consecutive "fly"
  ungroup() %>%
  select(-fly_group)                                   # Remove helper column

# Step 2: Define custom 60-minute intervals starting from the first timestamp
filtered_data <- fly_data %>%
  mutate(
    date_time = as.POSIXct(date_time, format = "%Y-%m-%d %H:%M:%S"), # Ensure proper datetime format
    time_diff = as.numeric(difftime(date_time, min(date_time), units = "mins")), # Minutes since the first timestamp
    custom_hour_group = floor(time_diff / 60) # Create 60-minute groups
  )

# Step: Calculate orthodrome-based segment direction, length, and drift metrics
grouped_data <- grouped_data %>%
  arrange(date, custom_hour_group, date_time) %>% # Ensure proper order
  mutate(
    # Create a global row number across the entire dataset
    global_row_number = row_number(),
    
    # Skip the very first row of the entire dataset for segment length and speed calculation
    segment.length = if_else(
      global_row_number > 1 & 
        !is.na(latitude) & !is.na(longitude) & 
        !is.na(lag(latitude)) & !is.na(lag(longitude)),
      distVincentySphere(
        cbind(longitude, latitude), 
        cbind(lag(longitude), lag(latitude))
      ) / 1000, # Convert to kilometers
      NA_real_ # Assign NA for first row and invalid rows
    ),
    
    # Calculate segment direction (orthodrome-based bearing)
    segment.dir = if_else(
      global_row_number > 1 & 
        !is.na(latitude) & !is.na(longitude) & 
        !is.na(lag(latitude)) & !is.na(lag(longitude)),
      bearing(
        cbind(lag(longitude), lag(latitude)), 
        cbind(longitude, latitude)
      ),
      NA_real_ # Assign NA for the first row and invalid rows
    ),
    
    # Use existing dtime for speed calculation, skipping the very first row
    segment.speed = if_else(
      global_row_number > 1 & 
        !is.na(segment.length) & !is.na(dtime) & dtime > 0,
      segment.length / dtime,
      NA_real_ # Assign NA otherwise
    )
  ) %>%
  group_by(date, custom_hour_group) %>%
  mutate(
    # Reference coordinates and direction for each hour
    long.start.ref = first(longitude),
    lat.start.ref = first(latitude),
    long.end.ref = last(longitude),
    lat.end.ref = last(latitude),
    dir.ref = bearing(cbind(long.start.ref, lat.start.ref), cbind(long.end.ref, lat.end.ref)) # Hourly reference direction
  ) %>%
  ungroup() %>%
  mutate(
    # Movement components relative to the hourly reference direction
    dir.delta2 = dir.ref - segment.dir, # Directional difference
    forward.speed = cos(dir.delta2 / (180 / pi)) * segment.speed, # Forward component of speed
    perpen.speed = -sin(dir.delta2 / (180 / pi)) * segment.speed, # Perpendicular component of speed
    
    # Wind components relative to the hourly reference direction
    dir.delta1 = dir.ref - wind.dir, # Wind direction difference
    tailwind = cos(dir.delta1 / (180 / pi)) * wind.speed, # Tailwind component
    crosswind = -sin(dir.delta1 / (180 / pi)) * wind.speed # Crosswind component
  ) %>%
  select(-global_row_number) # Drop the global_row_number column after calculations







# Track plot, Whole track LEAFLET--------------------------------------------------------------
# Step 1: Prepare the data
map_data <- grouped_data %>%
  filter(!is.na(latitude) & !is.na(longitude)) %>%
  mutate(
    day_group = as.integer(difftime(date_time, min(date_time), units = "days")),
    custom_hour_group = hour(date_time)
  )

# Extract the beginpoint and endpoint of the bird track
beginpoint <- map_data %>%
  slice(1) %>%
  select(longitude, latitude) %>%
  unlist() %>%
  as.numeric()

endpoint <- map_data %>%
  slice(n()) %>%
  select(longitude, latitude) %>%
  unlist() %>%
  as.numeric()

# Step 2: Create a color palette for the hour groups (0 to 23)
color_palette <- colorFactor(
  palette = brewer.pal(n = 9, "Set1"),
  domain = map_data$custom_hour_group
)

# Step 3: Calculate the orthodrome (great-circle route)
orthodrome_route <- gcIntermediate(
  p1 = beginpoint, p2 = endpoint,
  n = 100, addStartEnd = TRUE, sp = TRUE
)

# Step 4: Create the Leaflet map
leaflet(map_data) %>%
  addTiles() %>%  # Add base map
  # Add the flight path from GPS data
  addPolylines(
    lng = ~longitude, lat = ~latitude,
    color = ~color_palette(custom_hour_group),
    weight = 2, opacity = 0.8,
    group = ~paste(day_group, custom_hour_group)
  ) %>%
  # Add markers for each GPS point
  addCircleMarkers(
    lng = ~longitude, lat = ~latitude,
    radius = 3, color = ~color_palette(custom_hour_group),
    fill = TRUE, fillOpacity = 0.6,
    popup = ~paste0(
      "<b>Time:</b> ", date_time, "<br>",
      "<b>Speed:</b> ", speed, " km/h<br>",
      "<b>Day Group:</b> ", day_group, "<br>",
      "<b>Hour:</b> ", custom_hour_group
    )
  ) %>%
  # Overlay the orthodrome
  addPolylines(
    data = orthodrome_route,
    color = "blue", weight = 2, opacity = 0.8,
    label = "Orthodrome (Great Circle)"
  ) %>%
  # Add a legend
  addLegend(
    position = "bottomright",
    pal = color_palette, values = map_data$custom_hour_group,
    title = "Hour of Day (0-23)",
    opacity = 0.8
  )


# Extract wind data -------------------------------------------------------
# Step 1: Ensure the Grouped_data has required columns
grouped_data <- grouped_data %>%
  mutate(
    date = as.Date(date_time),  # Extract the date from the datetime
    hour = hour(date_time)      # Extract the hour (0-23) from the datetime
  )

# Step 2: Combine 'date' and 'hour' into a single datetime column
grouped_data <- grouped_data %>%
  mutate(
    datetime = as.POSIXct(paste(date, sprintf("%02d:00:00", hour)),
                          format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    datetime = format(as.POSIXct(datetime), "%Y-%m-%d %H:%M:%S")  # Reformat to character
  )

# Step 3: Initialize empty vectors for wind data
uwind <- numeric(nrow(grouped_data))
vwind <- numeric(nrow(grouped_data))

# Step 4: Extract wind data for each hourly point (loop)
for (i in 1:nrow(grouped_data)) {
  tryCatch({
    uwind[i] <- NCEP.interp(
      variable = 'uwnd', level = 925, 
      lat = grouped_data$latitude[i], lon = grouped_data$longitude[i], 
      dt = grouped_data$datetime[i], reanalysis2 = TRUE, 
      keep.unpacking.info = TRUE
    )
    vwind[i] <- NCEP.interp(
      variable = 'vwnd', level = 925, 
      lat = grouped_data$latitude[i], lon = grouped_data$longitude[i], 
      dt = grouped_data$datetime[i], reanalysis2 = TRUE, 
      keep.unpacking.info = TRUE
    )
  }, error = function(e) {
    message("Error at index ", i, ": ", e$message)
  })
}

# Step 5: Add the extracted wind data to the dataframe
grouped_data <- grouped_data %>%
  mutate(
    uwind = uwind,
    vwind = vwind
  )

# Step 6: Calculate wind direction and speed
grouped_data <- grouped_data %>%
  mutate(
    wind.dir = ifelse((atan2(uwind, vwind) * 180 / pi) < 0,
                      (atan2(uwind, vwind) * 180 / pi) + 360,
                      (atan2(uwind, vwind) * 180 / pi)), 
    wind.speed = sqrt(uwind^2 + vwind^2)
  )
                             
problematic_rows <- grouped_data[uwind == 0 & vwind == 0, ]

# REFERENCE DIRECTION AND FW SPEED VS TAILWIND ----------------------------------------------
# Step 1: Calculate orthodrome metrics
data_ort <- grouped_data %>%
  arrange(custom_hour_group, date_time) %>% # Ensure proper ordering
  group_by(date, custom_hour_group) %>%
  mutate(
    # Reference coordinates and direction for each hour
    long.start.ref = first(longitude),
    lat.start.ref = first(latitude),
    long.end.ref = last(longitude),
    lat.end.ref = last(latitude),
    dir.ref = bearing(cbind(long.start.ref, lat.start.ref), cbind(long.end.ref, lat.end.ref)) # Hourly reference direction
  ) %>%
  ungroup() %>%
  mutate(
    # Calculate movement components relative to the hourly reference direction
    dir.delta2 = dir.ref - segment.dir, # Directional difference
    forward.speed = cos(dir.delta2 / (180 / pi)) * segment.speed, # Forward component of speed
    perpen.speed = -sin(dir.delta2 / (180 / pi)) * segment.speed, # Perpendicular component of speed
    
    # Wind components relative to the hourly reference direction
    dir.delta1 = dir.ref - wind.dir, # Wind direction difference
    tailwind = cos(dir.delta1 / (180 / pi)) * wind.speed, # Tailwind component
    crosswind = -sin(dir.delta1 / (180 / pi)) * wind.speed # Crosswind component
  )


# Plotting -----------------------------------------------------------------------
# Plot for data_ort: Tailwind vs Forward Speed
plot(data_ort$tailwind, data_ort$forward.speed, 
     pch = 19, col = 'firebrick', 
     xlab = "Tailwind (m/s)", ylab = "Forward Speed (km/h)",
     main = "Tailwind vs Forward Speed (Orthodrome)")
abline(lm(forward.speed ~ tailwind, data = data_ort), col = "blue", lwd = 2)

# Plot for data_ort: Crosswind vs Perpendicular Speed
plot(data_ort$crosswind, data_ort$perpen.speed, 
     pch = 19, col = 'firebrick', 
     xlab = "Crosswind (m/s)", ylab = "Perpendicular Speed (km/h)",
     main = "Crosswind vs Perpendicular Speed (Orthodrome)")
abline(lm(perpen.speed ~ crosswind, data = data_ort), col = "blue", lwd = 2)



