# LIBRARIES ---------------------------------------------------------------
{
  library(dplyr)  # For data manipulation
  library(lubridate)  # For date-time manipulation
  library(geosphere)  # For distance calculations (if needed)
  library(ggplot2)  # For plotting
  library(RNCEP)  # For wind data extraction
  library(maps)  # For map data
  library(leaflet)  # For interactive maps
  library(geosphere)  # For distance calculations
  library(RColorBrewer)  # For color palettes
  library(rnaturalearth)  # For map data
  library(rnaturalearthdata)  # For map data
  library(patchwork)  # For combining plots
  library(foreach)  # For parallel processing
  library(doParallel)  # For parallel processing
  library(progress)  # For progress bar
  library(progressr)  # For progress bar
  library(parallel)  # For parallel processing
  library(iterators)  # For parallel processing
  library(doFuture)  # For parallel processing
  library(future) # For parallel processing
  library(gridExtra)  # For arranging plots
  library(sp)
  library(viridis)
  library(sf)
  library(patchwork)
  library(maps)
  library(ggspatial)
  library(prettymapr)
  library(circular)
  library(tidyr)
  library(grid)
  library(lutz)
  library(ggspatial)
  library(readr)
}

# LOADING DATA ------------------------------------------------------------
# set working directory
setwd("/Users/tom/Documents/Masters/_Github/Research_Project_2/RAW data individual tracks")

# remove list
rm(list=ls())

# Select entire folder with data
files <- list.files(path = "/Users/tom/Desktop/Research_project_2/Track_data_montagu_Harrier/Data on interruptions migrating harrier/RAW data individual tracks", pattern = "\\.csv", full.names = TRUE)
# Read all files into a list
results <- list()

# Loop over each file
for (i in seq_along(files)) {
  file <- files[i]
  data <- read_csv(file)
  
  # Extract Bird_ID correctly
  bird_id <- sub("^Data\\.All\\.|\\.csv\\.NEW\\.csv$", "", basename(file))  
  bird_id <- sub("(\\w+_\\w+_\\d+).*", "\\1", bird_id)  # Keep only "Name_season_year"
  
  # Assign Bird_ID  to each row
  data$Bird_ID <- bird_id  
  
  # Step 1: Filter for periods with at least 5 consecutive "fly"
  fly_data1 <- data %>%
    mutate(fly_group = cumsum(speed.class != "fly")) %>%
    group_by(fly_group) %>%
    filter(speed.class == "fly" & n() >= 5) %>%
    ungroup() %>%
    select(-fly_group)
  
  # Step 2: Define custom 60-minute intervals starting from the first timestamp
  fly_data1 <- fly_data1 %>% 
    mutate(date_time = as.POSIXct(date_time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))
  
  start_time <- min(fly_data1$date_time, na.rm = TRUE)
  end_time   <- max(fly_data1$date_time, na.rm = TRUE)
  breaks <- seq(from = start_time, to = end_time, by = "60 min")
  
  data_60min <- fly_data1 %>%
    mutate(interval = cut(date_time, breaks = breaks, right = FALSE))
  
  hour_dat <- data_60min %>%
    group_by(interval) %>%
    summarise(
      date            = first(date),
      start_date_time = first(date_time),
      end_date_time   = last(date_time),
      start_longitude = first(longitude),
      start_latitude  = first(latitude),
      end_longitude   = last(longitude),
      end_latitude    = last(latitude),
      segment.dir     = bearing(cbind(first(longitude), first(latitude)),
                                cbind(last(longitude), last(latitude))),
      dir.ref         = bearingRhumb(cbind(first(longitude), first(latitude)),
                                     cbind(last(longitude), last(latitude)))  # Adjust as needed
    ) %>%
    ungroup() %>%
    arrange(start_date_time) %>%
    group_by(date) %>%
    mutate(
      dtime = as.numeric(difftime(start_date_time, lag(start_date_time), units = "hours")),
      segment.length = distVincentySphere(
        cbind(lag(start_longitude), lag(start_latitude)),
        cbind(start_longitude, start_latitude)
      ) / 1000,
      segment.speed = segment.length / dtime,
      dir.delta2    = dir.ref - segment.dir,
      forward.speed = cos(dir.delta2 * pi/180) * segment.speed,
      perpen.speed  = -sin(dir.delta2 * pi/180) * segment.speed
    ) %>%
    ungroup() %>%
    group_by(day = as.Date(start_date_time)) %>%
    mutate(
      segment.length = if_else(row_number() == 1, NA_real_, segment.length),
      segment.speed  = if_else(row_number() == 1, NA_real_, segment.speed),
      forward.speed  = if_else(row_number() == 1, NA_real_, forward.speed),
      perpen.speed   = if_else(row_number() == 1, NA_real_, perpen.speed)
    ) %>%
    ungroup() %>%
    select(-day) %>%
    mutate(forward.speed = abs(forward.speed)) %>%
    group_by(date) %>%
    filter(n() >= 4) %>%
    ungroup()
  
  # Add Bird_ID to the final dataset
  hour_dat$Bird_ID <- bird_id  
  
  # ----------------------------
  results[[i]] <- hour_dat  # Save processed data for this file
}

# Combine all processed data frames into one
final_data <- do.call(rbind, results)


# EXTRACT WIND DATA -------------------------------------------------------

# Step 1: Ensure final_data has required columns
final_data <- final_data %>%
  mutate(
    date = as.Date(start_date_time),  # Extract the date
    hour = hour(start_date_time)      # Extract the hour (0-23)
  )

# Step 2: Create standardized datetime column
final_data <- final_data %>%
  mutate(
    datetime = as.POSIXct(paste(date, sprintf("%02d:00:00", hour)),
                          format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    datetime = format(as.POSIXct(datetime), "%Y-%m-%d %H:%M:%S")  # Reformat to character
  )

# Step 3: Initialize empty vectors for wind data
uwind <- numeric(nrow(final_data))
vwind <- numeric(nrow(final_data))

# Step 4: Extract wind data (parallel processing)
{
  plan(multisession)  # Enable parallel execution
  registerDoFuture()
  total_rows <- nrow(final_data)
  start_time <- Sys.time()
}

progressr::with_progress({
  p <- progressr::progressor(total_rows)  # Initialize progress bar
  
  wind_data <- foreach(i = 1:total_rows, .combine = rbind, .packages = "RNCEP") %dopar% {
    result <- tryCatch({
      # Retrieve wind data at 925 hPa level
      uwind <- NCEP.interp(
        variable = 'uwnd', level = 925, 
        lat = final_data$start_latitude[i], lon = final_data$start_longitude[i], 
        dt = final_data$datetime[i], reanalysis2 = TRUE, 
        keep.unpacking.info = TRUE
      )
      vwind <- NCEP.interp(
        variable = 'vwnd', level = 925, 
        lat = final_data$start_latitude[i], lon = final_data$start_longitude[i], 
        dt = final_data$datetime[i], reanalysis2 = TRUE, 
        keep.unpacking.info = TRUE
      )
      c(uwind, vwind)
    }, error = function(e) {
      message("Error at index ", i, ": ", e$message)
      c(NA, NA)
    })
    
    p(sprintf("Processing row %d", i))  # Update progress bar
    result
  }
})

final_time <- Sys.time() - start_time
cat("Processing completed in", round(final_time, 2), "minutes\n")

# Step 5: Add extracted wind data
final_data$uwind <- wind_data[, 1]
final_data$vwind <- wind_data[, 2]

# Step 6: Calculate wind direction and speed
final_data <- final_data %>%
  mutate(
    wind.dir = ifelse((atan2(uwind, vwind) * 180 / pi) < 0,
                      (atan2(uwind, vwind) * 180 / pi) + 360,
                      (atan2(uwind, vwind) * 180 / pi)), 
    wind.speed = sqrt(uwind^2 + vwind^2) * 3.6  # Convert m/s to km/h
  )

# Step 7: Calculate wind components relative to movement direction
final_data <- final_data %>%
  mutate(
    dir.delta1 = dir.ref - wind.dir,
    tailwind   = cos(dir.delta1 * (pi / 180)) * wind.speed,
    crosswind  = -sin(dir.delta1 * (pi / 180)) * wind.speed
  ) %>% 
  filter(!is.na(segment.dir))

# Step 8: Compute heading and airspeed
deg_to_rad <- function(deg) deg * (pi / 180)
rad_to_deg <- function(rad) rad * (180 / pi)

final_data <- final_data %>%
  mutate(
    V_x_ground = segment.speed * cos(deg_to_rad(segment.dir)),  
    V_y_ground = segment.speed * sin(deg_to_rad(segment.dir)),  
    V_x_wind   = wind.speed * cos(deg_to_rad(wind.dir)),
    V_y_wind   = wind.speed * sin(deg_to_rad(wind.dir)),
    V_x_air    = V_x_ground - V_x_wind,
    V_y_air    = V_y_ground - V_y_wind,
    airspeed   = sqrt(V_x_air^2 + V_y_air^2),
    heading    = (rad_to_deg(atan2(V_y_air, V_x_air))) %% 360,
    ground_speed = sqrt(V_x_ground^2 + V_y_ground^2)
  )


