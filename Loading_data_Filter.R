
# libraries & working directory -------------------------------------------
## load libraries 
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
}
# set working directory
setwd("/Users/tomtijsma/Documents/Masters/_Github/Research_Project_2/RAW data individual tracks")

# remove list 
rm(list=ls())

# Loading Data ------------------------------------------------------------
# load data
data <- read.csv("Data.All.Ronny_autumn_2014.csv.NEW.csv") 

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

#Step 3: Get global start/end reference points from the ENTIRE dataset (not grouped_data)
locs.start.ref <- cbind(data$longitude[1], data$latitude[1])
locs.end.ref <- cbind(data$longitude[nrow(data)], data$latitude[nrow(data)])

#Step 4: Calculate global reference direction (same for all rows)
global_dir.ref <- bearingRhumb(locs.start.ref, locs.end.ref)

#Step 5: Calculate the distance and direction between consecutive points
grouped_data <- filtered_data %>%
  arrange(date, date_time) %>% # Order by date and time within days
  group_by(date) %>% # Group by DAY instead of hour
  mutate(
    # Daily row counter (resets each day)
    daily_row_number = row_number(),
    
    # Segment length (daily-constrained)
    segment.length = if_else(
      daily_row_number > 1 & # Skip first row OF EACH DAY
        !is.na(latitude) & !is.na(longitude) &
        !is.na(lag(latitude)) & !is.na(lag(longitude)),
      distVincentySphere(
        cbind(longitude, latitude),
        cbind(lag(longitude), lag(latitude))
      ) / 1000, # km
      NA_real_
    ),
    
    # Segment direction (daily-constrained)
    segment.dir = if_else(
      daily_row_number > 1 & # Skip first row OF EACH DAY
        !is.na(latitude) & !is.na(longitude) &
        !is.na(lag(latitude)) & !is.na(lag(longitude)),
      bearing(
        cbind(lag(longitude), lag(latitude)),
        cbind(longitude, latitude)
      ),
      NA_real_
    ),
    
    # Segment speed (daily-constrained)
    segment.speed = if_else(
      daily_row_number > 1 &
        !is.na(segment.length) & !is.na(dtime) & dtime > 0,
      segment.length / dtime, # km/h
      NA_real_
    )
  ) %>%
  ungroup() %>%
  # Apply GLOBAL reference direction to all rows
  mutate(
    dir.ref = global_dir.ref, # Same value for every row
    # Movement components
    dir.delta2 = dir.ref - segment.dir,
    forward.speed = cos(dir.delta2 / (180 / pi)) * segment.speed,
    perpen.speed = -sin(dir.delta2 / (180 / pi)) * segment.speed
  ) %>%
  select(-daily_row_number)


# Track plot, Whole track LEAFLET--------------------------------------------------------------
# Step 1: Prepare the data for plotting
{
map_data <- grouped_data %>%
  filter(!is.na(latitude) & !is.na(longitude)) %>%
  mutate(
    day_group = as.integer(difftime(date_time, min(date_time), units = "days")),
    custom_hour_group = hour(date_time)
  )

# Step 2: Create orthodrome using global reference points
orthodrome_route <- gcIntermediate(
  p1 = as.numeric(locs.start.ref),  # Convert matrix to numeric vector
  p2 = as.numeric(locs.end.ref),    # Convert matrix to numeric vector
  n = 100, addStartEnd = TRUE, sp = TRUE
)

# Step 3: Create a color palette for the hour groups (0 to 23)
color_palette <- colorFactor(
  palette = brewer.pal(n = 9, "Set1"),
  domain = map_data$custom_hour_group
)

# plot the map
leaflet(map_data) %>%
  addTiles() %>%
  addPolylines(
    lng = ~longitude, lat = ~latitude,
    color = ~color_palette(custom_hour_group),
    weight = 2, opacity = 0.8,
    group = ~paste(day_group, custom_hour_group)
  ) %>%
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
  addPolylines(
    data = orthodrome_route,
    color = "blue", weight = 2, opacity = 0.8,
    label = "Orthodrome (Great Circle)"
  ) %>%
  addLegend(
    position = "bottomright",
    pal = color_palette, values = map_data$custom_hour_group,
    title = "Hour of Day (0-23)",
    opacity = 0.8
  )
}

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
# Set up parallel backend with progress support
{
  plan(multisession)  # Uses available cores
  registerDoFuture()  # Register the future backend

  # Total number of rows
  total_rows <- nrow(grouped_data)

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
        lat = grouped_data$latitude[i], lon = grouped_data$longitude[i], 
        dt = grouped_data$datetime[i], reanalysis2 = TRUE, 
        keep.unpacking.info = TRUE
      )
      vwind <- NCEP.interp(
        variable = 'vwnd', level = 925, 
        lat = grouped_data$latitude[i], lon = grouped_data$longitude[i], 
        dt = grouped_data$datetime[i], reanalysis2 = TRUE, 
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
grouped_data$uwind <- wind_data[, 1]
grouped_data$vwind <- wind_data[, 2]

# Step 6: Calculate wind direction and speed
grouped_data <- grouped_data %>%
  mutate(
    wind.dir = ifelse((atan2(uwind, vwind) * 180 / pi) < 0,
                      (atan2(uwind, vwind) * 180 / pi) + 360,
                      (atan2(uwind, vwind) * 180 / pi)), 
    wind.speed = sqrt(uwind^2 + vwind^2),
    wind.speed = wind.speed * 3.6  # Convert m/s to km/h
  )

# REFERENCE DIRECTION AND FW SPEED VS TAILWIND ----------------------------------------------
# step 1: Convert wind components to speed/direction
grouped_data <- grouped_data %>%
  mutate(
    wind.speed = sqrt(uwind^2 + vwind^2),  # Wind speed magnitude
    wind.dir = (atan2(uwind, vwind) * (180 / pi)) %% 360  # Wind direction (0-360°)
  ) %>%
  # Calculate tailwind/crosswind relative to reference direction (dir.ref)
  mutate(
    dir.delta1 = dir.ref - wind.dir,  # Direction difference
    tailwind = cos(dir.delta1 * (pi / 180)) * wind.speed,  # Tailwind component
    crosswind = -sin(dir.delta1 * (pi / 180)) * wind.speed  # Crosswind component
  )%>% 
  filter(!is.na(segment.dir))


# step 2: Plot the relationship between tailwind and forward speed
# Plot 1: Tailwind vs. Forward Speed (expect positive correlation)
plot(
  grouped_data$tailwind, grouped_data$forward.speed,
  main = "Tailwind vs. Forward Speed",
  xlab = "Tailwind (m/s)", ylab = "Forward Speed (km/h)",
  pch = 19, col = "firebrick"
)
abline(lm(forward.speed ~ tailwind, data = grouped_data), col = "navy", lwd = 2)

# Plot 2: Crosswind vs. Perpendicular Speed (drift check)
plot(
  grouped_data$crosswind, grouped_data$perpen.speed,
  main = "Crosswind vs. Perpendicular Speed",
  xlab = "Crosswind (m/s)", ylab = "Perpendicular Speed (km/h)",
  pch = 19, col = "firebrick"
)
abline(lm(perpen.speed ~ crosswind, data = grouped_data), col = "navy", lwd = 2)


# daily reference directions plots TW vs FS + CW vs PS----------------------------------------------
# Create new dataframe with daily reference directions
grouped_data_with_ref2 <- grouped_data %>%
  group_by(date) %>%  # Group by day
  mutate(
    # First coordinates of the day (start point)
    day_start_lon = first(longitude),
    day_start_lat = first(latitude),
    
    # Calculate daily reference direction (rhumb line to global endpoint)
    dir.ref2 = bearingRhumb(
      cbind(day_start_lon, day_start_lat),  # Start of the day
      matrix(locs.end.ref, ncol = 2)        # Global endpoint (replicated for each row)
    )
  ) %>%
  ungroup() %>%
  select(-day_start_lon, -day_start_lat)  # Remove temporary columns

# Recalculate delta using daily reference direction (dir.ref2)
grouped_data_with_ref2 <- grouped_data_with_ref2 %>%
  mutate(
    # Recalculate delta using daily reference direction (dir.ref2)
    dir.delta1 = dir.ref2 - wind.dir,
    tailwind = cos(dir.delta1 * (pi / 180)) * wind.speed,
    crosswind = -sin(dir.delta1 * (pi / 180)) * wind.speed
  )

# Split data by day (filtering NAs inline)
data_by_day <- grouped_data_with_ref2 %>%
  filter(!is.na(segment.dir)) %>%  # Remove NA segment.dir inline
  split(.$date)  # Split into daily dataframes

# Generate plots for each day
plot_list <- lapply(data_by_day, function(day_data) {
  # Plot 1: Tailwind vs. Forward Speed
  p1 <- ggplot(data_by_day, aes(tailwind, forward.speed)) +
    geom_point(color = "firebrick", alpha = 0.7) +
    geom_smooth(method = "lm", color = "navy", se = FALSE) +
    labs(
      title = paste("Day:", data_by_day$date[1]),
      x = "Tailwind (m/s)", 
      y = "Forward Speed (km/h)"
    ) +
    theme_minimal()
  
  # Plot 2: Crosswind vs. Perpendicular Speed
  p2 <- ggplot(data_by_day, aes(crosswind, perpen.speed)) +
    geom_point(color = "firebrick", alpha = 0.7) +
    geom_smooth(method = "lm", color = "navy", se = FALSE) +
    labs(
      title = paste("Day:", data_by_day$date[1]),
      x = "Crosswind (m/s)", 
      y = "Perpendicular Speed (km/h)"
    ) +
    theme_minimal()
  
  # Combine plots side-by-side
  grid.arrange(p1, p2, ncol = 2)
})

# Arrange all days in a vertical grid
final_pages <- marrangeGrob(
  grobs = plot_list, 
  nrow = 3,  # 3 days per page (vertical)
  ncol = 2,  # 2 plots per day (side-by-side)
  top = NULL
)
print(final_pages)


# overal daily plots ONE DIR.REF ------------------------------------------
# Calculate delta using the original dir.ref (no daily adjustment)
grouped_data_with_ref <- grouped_data %>%
  mutate(
    dir.delta1 = dir.ref - wind.dir,
    tailwind = cos(dir.delta1 * (pi / 180)) * wind.speed,
    crosswind = -sin(dir.delta1 * (pi / 180)) * wind.speed
  )

# Split data by day (filtering NAs inline)
data_by_day1 <- grouped_data_with_ref %>%
  filter(!is.na(segment.dir)) %>%  # Remove NA segment.dir
  split(.$date)  # Split into daily dataframes

# Generate plots for each day
plot_list <- lapply(data_by_day1, function(day_data) {
  # Plot 1: Tailwind vs. Forward Speed
  p1 <- ggplot(day_data, aes(tailwind, forward.speed)) +  # Fixed: use day_data
    geom_point(color = "firebrick", alpha = 0.7) +
    geom_smooth(method = "lm", color = "navy", se = FALSE) +
    labs(
      title = paste("Day:", day_data$date[1]),
      x = "Tailwind (m/s)", 
      y = "Forward Speed (km/h)"
    ) +
    theme_minimal()
  
  # Plot 2: Crosswind vs. Perpendicular Speed
  p2 <- ggplot(day_data, aes(crosswind, perpen.speed)) +  # Fixed: use day_data
    geom_point(color = "firebrick", alpha = 0.7) +
    geom_smooth(method = "lm", color = "navy", se = FALSE) +
    labs(
      title = paste("Day:", day_data$date[1]),
      x = "Crosswind (m/s)", 
      y = "Perpendicular Speed (km/h)"
    ) +
    theme_minimal()
  
  # Combine plots side-by-side
  grid.arrange(p1, p2, ncol = 2)
})

# Arrange all days in a vertical grid
final_pages <- marrangeGrob(
  grobs = plot_list, 
  nrow = 3,  # 3 days per page (vertical)
  ncol = 2,  # 2 plots per day (side-by-side)
  top = NULL
)
print(final_pages)



# track plot with wind vectors --------------------------------------------
map_data_ggplot <- grouped_data %>%
  filter(!is.na(latitude) & !is.na(longitude)) %>%
  mutate(
    date_label = as.Date(date_time),  # Extract actual date
    hour_group = hour(date_time),
    # Calculate wind vector endpoints (adjust 0.1 scaling as needed)
    wind_end_lon = longitude + 0.01 * wind.speed * sin(wind.dir * pi / 180),
    wind_end_lat = latitude + 0.01 * wind.speed * cos(wind.dir * pi / 180)
  ) %>%
  # Create ordered factor for proper facet labeling
  mutate(day_facet = factor(as.character(date_label), 
                            levels = as.character(sort(unique(date_label)))))

# Create base plot
migration_plot <- ggplot(map_data_ggplot, aes(x = longitude, y = latitude)) +
  # Bird migration path
  geom_path(aes(color = hour_group, group = date_label), 
            linewidth = 0.8, alpha = 0.8) +
  # Wind vectors
  geom_segment(aes(xend = wind_end_lon, yend = wind_end_lat,
                   alpha = wind.speed), 
               color = "darkred", 
               arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  # Styling
  scale_color_viridis_c(name = "Hour of Day", option = "plasma") +
  scale_alpha_continuous(name = "Wind Speed (km/h)", range = c(0.2, 0.8)) +
  labs(title = "Migration Path with Wind Vectors",
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "lightblue1", colour = "gray50"),
    panel.grid = element_line(color = "gray90"),
    legend.position = "right",
    strip.background = element_rect(fill = "gray80")
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
  
  facet_wrap(~ day_facet, scales = 'free')  # Facet by actual date

print(migration_plot)

# Heading calculation + Airspeed calculations -----------------------------------------------------
# Convert degrees to radians and vice versa
deg_to_rad <- function(deg) deg * (pi / 180)
rad_to_deg <- function(rad) rad * (180 / pi)

grouped_data <- grouped_data %>%
  mutate(
    # Convert wind speed to km/h and compute vector components
    V_x_ground = segment.speed * cos(deg_to_rad(segment.dir)),  
    V_y_ground = segment.speed * sin(deg_to_rad(segment.dir)),  
    
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


# plot heading on compass maps --------------------------------------------
# Reshape main data (segment, wind, heading) + change wind direction to "going to" direction
plot_data <- grouped_data %>%
  select(date_time, segment.dir, wind.dir, heading) %>%
  mutate(date = as.Date(date_time)) %>%  # Extract date from datetime
  select(date, segment.dir, wind.dir, heading) %>%  # Use extracted date
  pivot_longer(
    cols = -date,
    names_to = "type",
    values_to = "direction"
  ) %>%
  mutate(
    type = factor(
      type,
      levels = c("segment.dir", "wind.dir", "heading"),
      labels = c("Segment", "Wind", "Heading")
    ),
    # Convert wind direction to "going to" direction
    direction = ifelse(
      type == "Wind",
      (direction + 180) %% 360,  # Add 180° and wrap around 360
      direction  # Keep other directions unchanged
    )
  )

# Prepare reference direction data (one value per date)
ref_data <- grouped_data %>%
  mutate(date = as.Date(date_time)) %>%  # Extract date from datetime
  group_by(date) %>%
  summarise(direction = first(dir.ref)) %>%  # Assuming one dir.ref per date
  mutate(type = "Reference")  # Add type column


ggplot() +
  geom_segment(
    data = plot_data,
    aes(x = direction, xend = direction, y = 0, yend = 1, color = type),
    alpha = 0.3,
    linewidth = 0.5
  ) +
  geom_segment(
    data = ref_data,
    aes(x = direction, xend = direction, y = 0, yend = 1),
    color = "black",
    linewidth = 1.5,
    alpha = 0.8
  ) +
  coord_polar(start = -pi/2, direction = 1) +
  facet_wrap(~ date) +  
  scale_x_continuous(
    limits = c(0, 360),
    breaks = seq(0, 315, by = 45),
    labels = c("N", "45°", "E", "135°", "S", "225°", "W", "315°")
  ) +
  theme_minimal() +
  labs(
    title = "Heading, Segment, and Wind Directions (per Date)",
    subtitle = "Wind directions shown as flow direction (180° adjusted)",
    x = NULL,
    y = NULL,
    color = "Direction Type"
  ) +
  scale_color_manual(
    values = c(
      "Segment" = "blue",
      "Wind" = "red",
      "Heading" = "green"
    )
  )


# PLOT WIND PLUS HEADING  -------------------------------------------------
# Prepare data for plotting
map_data_ggplot1 <- grouped_data %>%
  filter(!is.na(latitude) & !is.na(longitude)) %>%
  mutate(
    date_label = as.Date(date_time),
    hour_group = hour(date_time),
    
    # Calculate wind vector endpoints (0.01 scaling)
    wind_end_lon = longitude + 0.01 * wind.speed * sin(deg_to_rad(wind.dir)),
    wind_end_lat = latitude + 0.01 * wind.speed * cos(deg_to_rad(wind.dir)),
    
    # Calculate heading vector endpoints (using airspeed for length)
    heading_end_lon = longitude + 0.01 * airspeed * sin(deg_to_rad(heading)),
    heading_end_lat = latitude + 0.01 * airspeed * cos(deg_to_rad(heading))
  ) %>%
  mutate(day_facet = factor(as.character(date_label), 
                            levels = as.character(sort(unique(date_label)))))

# Create base plot
# Create base plot
migration_plot1 <- ggplot(map_data_ggplot1, aes(x = longitude, y = latitude)) +
  geom_path(aes(color = hour_group, group = date_label), 
            linewidth = 1.2 , alpha = 1) +
  
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
  
  facet_wrap(~ day_facet, scales = 'free_y')

print(migration_plot1)

# Compute reference line start and end points
ref_lines <- map_data_ggplot1 %>%
  group_by(day_facet) %>%
  summarize(
    mean_lon = mean(longitude, na.rm = TRUE),
    mean_lat = mean(latitude, na.rm = TRUE),
    dir = mean(global_dir.ref, na.rm = TRUE), # Use the reference direction for the facet
    .groups = "drop"
  ) %>%
  mutate(
    dir_rad = (270 - dir) * pi / 180, # Convert to correct angle
    ref_end_lon = mean_lon + 1 * cos(dir_rad),  # Adjust length as needed
    ref_end_lat = mean_lat + 0.3 * sin(dir_rad)
  )

# Add reference direction to the plot
migration_plot2 <- migration_plot1 +
  geom_segment(data = ref_lines, aes(
    x = ref_end_lon, y = ref_end_lat,
    xend = mean_lon, yend = mean_lat
  ), 
  linetype = "dashed", 
  color = "black", 
  linewidth = 1, 
  arrow = arrow(length = unit(0.2, "cm"), type = "closed")  # Add arrowhead
  )
print(migration_plot2)

# optimal reference direction  --------------------------------------------
# Define the range of reference directions
directions <- seq(90, 300, by = 5)

# Create an empty dataframe to store R-squared values for each reference direction
r_squared_values <- data.frame(dir.ref = directions, R_squared = NA)

# Loop through each reference direction
for (direction in directions) {
  # Calculate forward.speed and tailwind for all rows for the current reference direction
  forward_speeds <- sapply(1:nrow(grouped_data), function(i) {
    row_data <- grouped_data[i, ]
    dir.delta2 <- direction - row_data$segment.dir
    forward.speed <- cos(dir.delta2 * (pi / 180)) * row_data$segment.speed
    return(forward.speed)
  })
  
  tailwinds <- sapply(1:nrow(grouped_data), function(i) {
    row_data <- grouped_data[i, ]
    dir.delta1 <- direction - row_data$wind.dir
    tailwind <- cos(dir.delta1 * (pi / 180)) * row_data$wind.speed
    return(tailwind)
  })
  
  # Fit a linear model and calculate the R-squared value for this reference direction
  model <- lm(forward_speeds ~ tailwinds)
  r_squared_values$R_squared[r_squared_values$dir.ref == direction] <- summary(model)$r.squared
}

# get the original reference direction
original_dir_ref <- grouped_data$dir.ref[1]

# Plot the average R-squared values
ggplot(r_squared_values, aes(x = dir.ref, y = R_squared)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = original_dir_ref, linetype = "dashed", color = "red") + # Add original reference direction line
  labs(x = "Reference Direction (°)", y = "Average R-squared of Tailwind vs Forward Speed",
       title = "Average R-squared of Tailwind vs Forward Speed for Various Reference Directions") +
  theme_minimal()
