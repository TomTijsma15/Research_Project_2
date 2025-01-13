
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



# Track plot --------------------------------------------------------------
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

