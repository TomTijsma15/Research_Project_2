
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

# LOAD DATA ---------------------------------------------------------------

# set working directory
setwd("/Users/tom/Documents/Masters/_Github/Research_Project_2/RAW data individual tracks")

# remove list
rm(list=ls())

# load data
data <- read_csv("Data.All.Ronny_autumn_2012.csv.NEW.csv") 

{
# Step 1: Filter for periods with at least 5 consecutive "fly"
fly_data1 <- data %>%
  mutate(fly_group = cumsum(speed.class != "fly")) %>%  # Group consecutive "fly" rows
  group_by(fly_group) %>%
  filter(speed.class == "fly" & n() >= 5) %>%          # Keep groups with at least 5 consecutive "fly"
  ungroup() %>%
  select(-fly_group)                                   # Remove helper column

# Step 2: Define custom 60-minute intervals starting from the first timestamp
filtered_data1 <- fly_data1 %>%
  mutate(
    date_time = as.POSIXct(date_time, format = "%Y-%m-%d %H:%M:%S"), # Ensure proper datetime format
    time_diff = as.numeric(difftime(date_time, min(date_time), units = "mins")), # Minutes since the first timestamp
    custom_hour_group = floor(time_diff / 60) # Create 60-minute groups
  )

#Step 3: Get global start/end reference points from the ENTIRE dataset (not grouped_data)
locs.start.ref <- cbind(data$longitude[1], data$latitude[1])
locs.end.ref <- cbind(data$longitude[nrow(data)], data$latitude[nrow(data)])

#Step 4: Calculate global reference direction (same for all rows)
global_dir.ref <- bearingRhumb(locs.start.ref, locs.end.ref)


# CREATE NEW DATAFRAME FOR HOURLY DATA ------------------------------------
# Ensure your date_time column is correctly converted to POSIXct.
fly_data1 <- fly_data1 %>% 
  mutate(date_time = as.POSIXct(date_time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))

# Define the start and end times for your data.
start_time <- min(fly_data1$date_time, na.rm = TRUE)
end_time   <- max(fly_data1$date_time, na.rm = TRUE)

# Create breaks for 60-minute intervals starting from the first timestamp.
breaks <- seq(from = start_time, to = end_time, by = "60 min")

# Create a new data frame with an "interval" column indicating each 60-minute period.
data_60min <- fly_data1 %>%
  mutate(interval = cut(date_time, breaks = breaks, right = FALSE))

# Summarise the data by 60-minute intervals.
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
    segment.dir     = bearing(
      cbind(first(longitude), first(latitude)),
      cbind(last(longitude), last(latitude))
    ),
    dir.ref         = global_dir.ref # Use the global reference direction
  ) %>%
  ungroup() %>%
  arrange(start_date_time) %>%
  # Group by date so we calculate differences within each day.
  group_by(date) %>%
  mutate(
    # Calculate time difference between consecutive hours using start_date_time.
    dtime = as.numeric(difftime(start_date_time, lag(start_date_time), units = "hours")),
    # Calculate segment length using the start coordinates of consecutive hours.
    segment.length = distVincentySphere(
      cbind(lag(start_longitude), lag(start_latitude)),
      cbind(start_longitude, start_latitude)
    ) / 1000,
    # Calculate segment speed.
    segment.speed = segment.length / dtime, # Speed in km/h
    dir.delta2    = dir.ref - segment.dir, # Direction difference
    forward.speed = cos(dir.delta2 * pi/180) * segment.speed, # Calculate forward speed
    perpen.speed  = -sin(dir.delta2 * pi/180) * segment.speed # Calculate perpendicular speed
  ) %>%
  ungroup() %>%
  # For the first row of each day, we don't have a previous hour so set the metrics to NA.
  group_by(day = as.Date(start_date_time)) %>%
  mutate(
    segment.length = if_else(row_number() == 1, NA_real_, segment.length),
    segment.speed  = if_else(row_number() == 1, NA_real_, segment.speed),
    forward.speed  = if_else(row_number() == 1, NA_real_, forward.speed),
    perpen.speed   = if_else(row_number() == 1, NA_real_, perpen.speed)
  ) %>%
  ungroup() %>%
  select(-day) %>%
  mutate(forward.speed = abs(forward.speed)) %>% # Ensure forward speed is always positive
  # Filter out days with less than 4 hours of flying
  group_by(date) %>%
  filter(n() >= 4) %>%
  ungroup()
}
####


# EXTRACT WIND DATA -------------------------------------------------------
# Step 1: Ensure the hour_dat has required columns
{
hour_dat <- hour_dat %>%
  mutate(
    date = as.Date(start_date_time),  # Extract the date from the datetime
    hour = hour(start_date_time)      # Extract the hour (0-23) from the datetime
  )

# Step 2: Combine 'date' and 'hour' into a single datetime column
hour_dat <- hour_dat %>%
  mutate(
    datetime = as.POSIXct(paste(date, sprintf("%02d:00:00", hour)),
                          format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    datetime = format(as.POSIXct(datetime), "%Y-%m-%d %H:%M:%S")  # Reformat to character
  )

# Step 3: Initialize empty vectors for wind data
uwind <- numeric(nrow(hour_dat))
vwind <- numeric(nrow(hour_dat))
}
# Step 4: Extract wind data for each hourly point (loop)
# Set up parallel backend with progress support
{
  plan(multisession)  # Uses available cores
  registerDoFuture()  # Register the future backend
  
  # Total number of rows
  total_rows <- nrow(hour_dat)
  
  # Record start time
  start_time <- Sys.time()
}
# Process with progress bar
progressr::with_progress({
  p <- progressr::progressor(total_rows)  # Initialize progressor
  
  wind_data <- foreach(i = 1:total_rows, .combine = rbind, .packages = "RNCEP") %dopar% {
    # Perform the task
    result <- tryCatch({
      # Your existing wind calculation code
      uwind <- NCEP.interp(
        variable = 'uwnd', level = 925, 
        lat = hour_dat$start_latitude[i], lon = hour_dat$start_longitude[i], 
        dt = hour_dat$datetime[i], reanalysis2 = TRUE, 
        keep.unpacking.info = TRUE
      )
      vwind <- NCEP.interp(
        variable = 'vwnd', level = 925, 
        lat = hour_dat$start_latitude[i], lon = hour_dat$start_longitude[i], 
        dt = hour_dat$datetime[i], reanalysis2 = TRUE, 
        keep.unpacking.info = TRUE
      )
      c(uwind, vwind)
    }, error = function(e) {
      message("Error at index ", i, ": ", e$message)
      c(NA, NA)
    })
    
    # Update progress bar (atomic operation)
    p(sprintf("Processing row %d", i))  # Updates both counter and message
    
    result
  }
})

# Print final time
{
  final_time <- Sys.time() - start_time
  cat("Processing completed in", round(final_time, 2), "minutes\n")
}

# Step 5: Add the extracted wind data to the dataframe
{ hour_dat$uwind <- wind_data[, 1]
hour_dat$vwind <- wind_data[, 2]
 } 
# Step 6: Calculate wind direction and speed
{
  hour_dat <- hour_dat %>%
  mutate(
    wind.dir = ifelse((atan2(uwind, vwind) * 180 / pi) < 0,
                      (atan2(uwind, vwind) * 180 / pi) + 360,
                      (atan2(uwind, vwind) * 180 / pi)), 
    wind.speed = sqrt(uwind^2 + vwind^2),
    wind.speed = wind.speed * 3.6  # Convert m/s to km/h
  )

# step 7: Convert wind components to speed/direction
hour_dat <- hour_dat %>%
  mutate(
    wind.dir = (atan2(uwind, vwind) * (180 / pi)) %% 360  # Wind direction (0-360°)
  ) %>%
  # Calculate tailwind/crosswind relative to reference direction (dir.ref)
  mutate(
    dir.delta1 = dir.ref - wind.dir,  # Direction difference
    tailwind = cos(dir.delta1 * (pi / 180)) * wind.speed,  # Tailwind component
    crosswind = -sin(dir.delta1 * (pi / 180)) * wind.speed  # Crosswind component
  )%>% 
  filter(!is.na(segment.dir))

# step 8: calculate heading and air speed 
# Convert degrees to radians and vice versa
deg_to_rad <- function(deg) deg * (pi / 180)
rad_to_deg <- function(rad) rad * (180 / pi)

hour_dat <- hour_dat %>%
  mutate(
    # Compute ground speed components
    V_x_ground = segment.speed * cos(deg_to_rad(segment.dir)),  
    V_y_ground = segment.speed * sin(deg_to_rad(segment.dir)),  
    
    # Compute wind speed components
    V_x_wind = wind.speed * cos(deg_to_rad(wind.dir)),
    V_y_wind = wind.speed * sin(deg_to_rad(wind.dir)),
    
    # Compute airspeed components
    V_x_air = V_x_ground - V_x_wind,
    V_y_air = V_y_ground - V_y_wind,
    
    # Compute airspeed magnitude
    airspeed = sqrt(V_x_air^2 + V_y_air^2),
    
    # Compute heading (direction of airspeed)
    heading = (rad_to_deg(atan2(V_y_air, V_x_air))) %% 360,
    
    # Ground speed (which you already have)
    ground_speed = sqrt(V_x_ground^2 + V_y_ground^2)
  )
 }
# CHECK THE VECOTRS -------------------------------------------------------
# For each observation, create a plot showing the airspeed, wind, and ground vectors
hour_dat %>%
  ggplot(aes(x = start_longitude, y = start_latitude)) +
  # Airspeed vector
  geom_segment(aes(xend = start_longitude + V_x_air, 
                   yend = start_latitude + V_y_air), 
               color = "goldenrod", arrow = arrow(length = unit(0.2, "cm")), linewidth = 1) +
  
  # Wind vector
  geom_segment(aes(xend = start_longitude + V_x_wind, 
                   yend = start_latitude + V_y_wind), 
               color = "deepskyblue", arrow = arrow(length = unit(0.2, "cm")), linewidth = 1) +
  
  # Ground vector
  geom_segment(aes(xend = start_longitude + V_x_ground, 
                   yend = start_latitude + V_y_ground), 
               color = "darkgreen", arrow = arrow(length = unit(0.2, "cm")), linewidth = 1) +
  
  # Parallelogram edges (optional)
  geom_segment(aes(x = start_longitude + V_x_air, y = start_latitude + V_y_air,
                   xend = start_longitude + V_x_wind + V_x_air, 
                   yend = start_latitude + V_y_wind + V_y_air),
               color = "gray", linetype = "dashed") +
  geom_segment(aes(x = start_longitude + V_x_wind, y = start_latitude + V_y_wind,
                   xend = start_longitude + V_x_ground, 
                   yend = start_latitude + V_y_ground),
               color = "gray", linetype = "dashed") +
  
  labs(title = "Airspeed, Wind, and Ground Vectors with Parallelogram",
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  coord_fixed(ratio = 1) + # Maintain aspect ratio
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),  # Remove ticks
    legend.title = element_blank()  # Remove legend title
  )



# SCATTER PLOTS TAILWIND VS FORWARD SPEED ---------------------------------

# Plot 1: Tailwind vs. Forward Speed
{
  p1 <- ggplot(hour_dat, aes(x = tailwind, y = forward.speed)) +
  geom_point(color = "firebrick", size = 2) +
  geom_smooth(method = "lm", color = "navy", se = FALSE, linewidth = 1) +
  labs(title = "Tailwind vs. Forward Speed",
       x = "Tailwind (km/h)",
       y = "Forward Speed (km/h)") +
  theme_minimal()

# Plot 2: Crosswind vs. Perpendicular Speed
p2 <- ggplot(hour_dat, aes(x = crosswind, y = perpen.speed)) +
  geom_point(color = "firebrick", size = 2) +
  geom_smooth(method = "lm", color = "navy", se = FALSE, linewidth = 1) +
  labs(title = "Crosswind vs. Perpendicular Speed",
       x = "Crosswind (km/h)",
       y = "Perpendicular Speed (km/h)") +
  theme_minimal()

# Combine the two plots side by side using patchwork
combined_plot <- p1 + p2

# Display the combined plot
print(combined_plot)
}


# TRACK PLOT PLUS WIND & AIRSPEED VECTORS ---------------------------------
# Prepare data for plotting
map_data_ggplot <- hour_dat %>%
  filter(!is.na(start_latitude) & !is.na(start_longitude)) %>%
  mutate(
    date_label = as.Date(datetime),
    hour_group = hour(datetime),
    
    # Calculate wind vector endpoints (0.01 scaling)
    wind_end_lon = start_longitude + 0.01 * wind.speed * sin(deg_to_rad(wind.dir)),
    wind_end_lat = start_latitude + 0.01 * wind.speed * cos(deg_to_rad(wind.dir)),
    
    # Calculate heading vector endpoints (using airspeed for length)
    heading_end_lon = start_longitude + 0.01 * airspeed * sin(deg_to_rad(heading)),
    heading_end_lat = start_latitude + 0.01 * airspeed * cos(deg_to_rad(heading))
  ) %>%
  mutate(day_facet = factor(as.character(date_label), 
                            levels = as.character(sort(unique(date_label)))))

# Create base plot
migration_plot <- ggplot(map_data_ggplot, aes(x = start_longitude, y = start_latitude)) +
  geom_path(aes(color = hour_group, group = date_label), 
            linewidth = 1.2, alpha = 1) +
  
  # Wind vector arrows (blue)
  geom_segment(aes(xend = wind_end_lon, yend = wind_end_lat,
                   alpha = wind.speed), 
               color = "deepskyblue", 
               arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +  
  
  # Bird airspeed vectors (goldenrod)
  geom_segment(aes(xend = heading_end_lon, yend = heading_end_lat,
                   linewidth = airspeed),  
               color = "goldenrod1", alpha = 0.8,
               arrow = arrow(length = unit(0.1, "cm"), type = "closed")) +  
  
  scale_color_viridis_c(name = "Hour of Day", option = "plasma") +
  scale_alpha_continuous(name = "Wind Speed (km/h)", range = c(0.3, 0.8)) +
  scale_linewidth_continuous(name = "Airspeed (km/h)", range = c(0.3, 0.8)) +  
  
  labs(title = "Migration Path with Wind Vectors and Bird Heading (Airspeed)",
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", colour = "gray40"),
    legend.position = "right",
    strip.background = element_rect(fill = "gray40")
  ) +
  # Add custom north arrow
  annotation_custom(
    grob = grid::gTree(children = grid::gList(
      # North Arrow
      grid::segmentsGrob(
        x0 = 0.05, y0 = 0.85, 
        x1 = 0.05, y1 = 0.95, 
        gp = grid::gpar(col = "black", lwd = 2),
        arrow = grid::arrow(
          type = "closed",
          length = unit(0.1, "inches")
        )
      ),
      # "N" label
      grid::textGrob(
        label = "N",
        x = 0.05, y = 0.82,               # Position below the arrow
        gp = grid::gpar(col = "black", fontface = "bold", fontsize = 10)
      )
    )),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) + 
  
  facet_wrap(~ day_facet, scales = 'free')

# Display the plot
print(migration_plot)

# WORKING UTM PROJECTION; CUSTOM PER TRACK --------------------------------------------------
#### CUSTOM UTM PROJECTION FROM TRANSFERSE MERCATOR PROJECTION +- MAX/MIN LATITUDE 
## THIS WORKS !!!!
{
  # --- 1. Prepare the Data ---
  map_data <- hour_dat %>%
    filter(!is.na(start_latitude) & !is.na(start_longitude)) %>%
    mutate(
      date_label = as.Date(datetime),
      hour_group = hour(datetime)
    )
  
  # --- 2. Convert to sf (WGS84) ---
  sf_data <- st_as_sf(map_data, coords = c("start_longitude", "start_latitude"), crs = 4326)
  
  # --- 3. Define a Custom Transverse Mercator Projection ---
  # Calculate the bounding box of your data and determine the custom central meridian
  bbox <- st_bbox(sf_data)
  custom_central_meridian <- mean(c(bbox["xmin"], bbox["xmax"]))
  
  # Build a custom CRS string (similar to UTM, but with our custom central meridian)
  custom_crs <- paste0(
    "+proj=tmerc ",
    "+lat_0=0 ",
    "+lon_0=", custom_central_meridian, " ",
    "+k=0.9996 ",          # scale factor (same as UTM)
    "+x_0=500000 ",        # false easting (same as UTM)
    "+y_0=0 ",
    "+datum=WGS84 ",
    "+units=m ",
    "+no_defs"
  )
  
  # --- 4. Transform Data to the Custom Projection ---
  sf_data_custom <- st_transform(sf_data, crs = custom_crs)
  
  # Extract UTM-like coordinates (in meters) and add them to map_data
  coords_custom <- st_coordinates(sf_data_custom)
  map_data <- map_data %>% 
    mutate(easting = coords_custom[, 1],
           northing = coords_custom[, 2])
  
  # --- 5. Compute Vector Endpoints in the Custom Projection ---
  # Helper function: convert degrees to radians
  deg_to_rad <- function(deg) { deg * pi / 180 }
  
  # Set a scaling factor for visualizing the vectors (speeds remain in km/h)
  scaling_factor_air <- 1500  # Adjust for visual effect
  scaling_factor_wind <- 1500  # Adjust for visual effect
  
  # Calculate wind and heading vector endpoints
  map_data <- map_data %>%
    mutate(
      wind_end_easting = ifelse(!is.na(wind.speed) & !is.na(wind.dir),
                                easting + scaling_factor_wind * wind.speed * sin(deg_to_rad(wind.dir)),
                                NA),
      wind_end_northing = ifelse(!is.na(wind.speed) & !is.na(wind.dir),
                                 northing + scaling_factor_wind * wind.speed * cos(deg_to_rad(wind.dir)),
                                 NA),
      heading_end_easting = ifelse(!is.na(airspeed) & !is.na(heading),
                                   easting + scaling_factor_air * airspeed * sin(deg_to_rad(heading)),
                                   NA),
      heading_end_northing = ifelse(!is.na(airspeed) & !is.na(heading),
                                    northing + scaling_factor_air * airspeed * cos(deg_to_rad(heading)),
                                    NA)
    )
  
  # --- 6. Prepare the World Map Background in the Custom Projection ---
  world_map <- ne_countries(scale = "medium", returnclass = "sf")
  world_map_custom <- st_transform(world_map, crs = custom_crs)
  
  # --- 7. Create a Plot for Each Day (with a 5° Margin) ---
  unique_days <- unique(map_data$date_label)
  plot_list <- list()
  
}
for(day in unique_days) {
  # Filter and sort the day's data by datetime (to preserve track order)
  day_data <- map_data %>% filter(date_label == day) %>% arrange(datetime)
  
  # -- Calculate Expanded Bounding Box in WGS84 --
  # Use original geographic coordinates to set a 5° margin
  lon_min <- min(day_data$start_longitude, na.rm = TRUE)
  lon_max <- max(day_data$start_longitude, na.rm = TRUE)
  lat_min <- min(day_data$start_latitude, na.rm = TRUE)
  lat_max <- max(day_data$start_latitude, na.rm = TRUE)
  
  # Expand by 1° on each side
  lon_min_exp <- lon_min - 1
  lon_max_exp <- lon_max + 1
  lat_min_exp <- lat_min - 1
  lat_max_exp <- lat_max + 1
  
  # Create a bounding box (polygon) in WGS84 and transform it to the custom projection
  bbox_polygon <- st_as_sfc(st_bbox(c(xmin = lon_min_exp, xmax = lon_max_exp,
                                      ymin = lat_min_exp, ymax = lat_max_exp),
                                    crs = 4326))
  bbox_projected <- st_transform(bbox_polygon, crs = custom_crs)
  bbox_proj <- st_bbox(bbox_projected)
  
  x_lim <- c(bbox_proj["xmin"], bbox_proj["xmax"])
  y_lim <- c(bbox_proj["ymin"], bbox_proj["ymax"])
  
  # -- Build the Plot --
  p <- ggplot() +
    # World map background in the custom projection
    geom_sf(data = world_map_custom, fill = "gray95", color = "gray50") +
    
    # Plot the day's track (colored by hour_group)
    geom_path(data = day_data, aes(x = easting, y = northing, color = hour_group),
              linewidth = 1.2, alpha = 1) +
    
    # Wind vector segments
    geom_segment(data = day_data %>% filter(!is.na(wind_end_easting)),
                 aes(x = easting, y = northing,
                     xend = wind_end_easting, yend = wind_end_northing,
                     alpha = wind.speed),
                 color = "deepskyblue",
                 arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
    
    # Airspeed (heading) vector segments
    geom_segment(data = day_data %>% filter(!is.na(heading_end_easting)),
                 aes(x = easting, y = northing,
                     xend = heading_end_easting, yend = heading_end_northing,
                     linewidth = airspeed),
                 color = "goldenrod1",
                 alpha = 0.8,
                 arrow = arrow(length = unit(0.1, "cm"), type = "closed")) +
    
    # Zoom to the expanded bounding box
    coord_sf(xlim = x_lim, ylim = y_lim) +
    
    # Scales and labels
    scale_color_viridis_c(name = "Hour of Day", option = "plasma") +
    scale_alpha_continuous(name = "Wind Speed (km/h)", range = c(0.3, 0.8)) +
    scale_linewidth_continuous(name = "Airspeed (km/h)", range = c(0.3, 0.8)) +
    labs(title = paste("Migration Track on", format(as.Date(day), "%d-%b-%Y")),
         x = "Easting (m)", y = "Northing (m)") +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) + 
    
    # Add custom north arrow
    annotation_custom(
      grob = grid::gTree(children = grid::gList(
        grid::segmentsGrob(
          x0 = 0.05, y0 = 0.85, 
          x1 = 0.05, y1 = 0.95, 
          gp = grid::gpar(col = "black", lwd = 2),
          arrow = grid::arrow(type = "closed", length = unit(0.1, "inches"))
        ),
        grid::textGrob(
          label = "N",
          x = 0.05, y = 0.82,
          gp = grid::gpar(col = "black", fontface = "bold", fontsize = 10)
        )
      )),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    )
  
  # Store the plot in the list
  plot_list[[as.character(day)]] <- p
}


# Print all plots
for (day in names(plot_list)) {
  print(plot_list[[day]])
}

# export the plots into a PDF 
# Specify the PDF file name and open the PDF device
pdf("Track_Marc_Spring_2013_3.pdf", height = 8, width = 8) # INSERT NAME YEAR AND SEASON OF THE BIRD 

# Loop through each plot in the plot_list and print it to the PDF
for (day in names(plot_list)) {
  print(plot_list[[day]])  # Print each plot to the PDF
}

# Close the PDF device
dev.off()



