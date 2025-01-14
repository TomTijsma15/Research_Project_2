
# libraries & working directory -------------------------------------------
## load libraries 
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

# Step 3: View data grouped by custom intervals
grouped_data <- filtered_data %>%
  arrange(custom_hour_group, date_time) # Sort by group and time



# Track plot, Whole track LEAFLET--------------------------------------------------------------
# Step 1: Prepare the data
map_data <- grouped_data %>%
  filter(!is.na(latitude) & !is.na(longitude)) %>%
  mutate(
    # Create 'day_group' which counts days from the first observation
    day_group = as.integer(difftime(date_time, min(date_time), units = "days")),
    
    # Extract the hour of the day from the 'date_time' column (0 to 23)
    custom_hour_group = hour(date_time)  # Get hour from date_time (0 to 23)
  )

# Step 2: Create a color palette for the hour groups (0 to 23)
color_palette <- colorFactor(
  palette = brewer.pal(n = 9, "Set1"), # Set1 palette, adjust if you have more than 9 hours
  domain = map_data$custom_hour_group
)

# Step 3: Create the Leaflet map
leaflet(map_data) %>%
  addTiles() %>%  # Add base map
  addPolylines(
    lng = ~longitude, lat = ~latitude,  # Use longitude and latitude for the flight path
    color = ~color_palette(custom_hour_group),  # Color based on hour group (0-23)
    weight = 2, opacity = 0.8,
    group = ~paste(day_group, custom_hour_group)  # Group by day and hour
  ) %>%
  addCircleMarkers(
    lng = ~longitude, lat = ~latitude,  # Add markers for each GPS point
    radius = 3, color = ~color_palette(custom_hour_group), fill = TRUE, fillOpacity = 0.6,
    popup = ~paste0(
      "<b>Time:</b> ", date_time, "<br>",
      "<b>Speed:</b> ", speed, " km/h<br>",
      "<b>Day Group:</b> ", day_group, "<br>",
      "<b>Hour:</b> ", custom_hour_group
    )  # Add popups with additional information
  ) %>%
  addLegend(
    position = "bottomright",
    pal = color_palette, values = map_data$custom_hour_group,
    title = "Hour of Day (0-23)",
    opacity = 0.8
  )


# Track plot Per day GGPLOT ------------------------------------------------------
# Step 1: Prepare World Map
world_map <- ne_countries(scale = "medium", returnclass = "sf")

# Step 2: Prepare Bird Data
map_data <- map_data %>%
  filter(!is.na(latitude) & !is.na(longitude)) %>%
  mutate(
    longitude = ifelse(longitude > 180, longitude - 360, longitude)
  )

# Step 3: Calculate the Bounding Box (min/max lat/lon) for Each Day
map_data_limits <- map_data %>%
  group_by(date) %>%
  summarise(
    lon_min = min(longitude, na.rm = TRUE),
    lon_max = max(longitude, na.rm = TRUE),
    lat_min = min(latitude, na.rm = TRUE),
    lat_max = max(latitude, na.rm = TRUE)
  )

# Step 4: Create the Plot with Dynamic Zoom for Each Day
ggplot() +
  # Add the world map as the background
  geom_sf(data = world_map, fill = "lightgray", color = "white") +
  # Add the bird track with increased linewidth for visibility
  geom_path(
    data = map_data,
    aes(x = longitude, y = latitude, group = date),
    linewidth = 2, alpha = 0.8, color = "blue"
  ) +
  # Add points for the bird positions
  geom_point(
    data = map_data,
    aes(x = longitude, y = latitude),
    size = 3, alpha = 0.6, color = "red"
  ) +
  # Customize the appearance
  labs(
    title = "Bird Migration Tracks by Day",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  # Facet the plot by date
  #facet_wrap(~date) +
  # Apply coord_cartesian() to dynamically zoom for each day
  coord_cartesian(
    xlim = c(-180, 180), # Keep a global limit for longitude
    ylim = c(-90, 90),    # Keep a global limit for latitude
    expand = FALSE
  ) +
  theme(
    strip.text = element_text(size = 10), # Smaller facet labels
    axis.text = element_text(size = 8)    # Axis text size
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

# Step 2.1: Transform longitudes to the 0-360 range
grouped_data <- grouped_data %>%
  mutate(
    longitude = ifelse(longitude < 0, 360 + longitude, longitude)
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
                      (atan2(uwind, vwind) * 180 / pi))
  )
                             
problematic_rows <- grouped_data[uwind == 0 & vwind == 0, ]
print(problematic_rows)

# Plot with the wind vectors ----------------------------------------------

