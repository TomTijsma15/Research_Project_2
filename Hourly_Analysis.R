
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
}

# LOAD DATA ---------------------------------------------------------------

# set working directory
setwd("/Users/tom/Documents/Masters/_Github/Research_Project_2/RAW data individual tracks")

# remove list
rm(list=ls())

# load data
data <- read.csv("Data.All.Corry_autumn_2012.csv.NEW.csv") 

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
    # Use the first record’s date and time as the start, and the last as the end.
    date = first(date),
    start_date_time = first(date_time),
    end_date_time   = last(date_time),
    
    # Use the first record’s coordinates as the start and the last record’s as the end.
    start_longitude = first(longitude),
    start_latitude  = first(latitude),
    end_longitude   = last(longitude),
    end_latitude    = last(latitude),
    
    # Calculate segment length based on the start and end of the interval.
    segment.length = distVincentySphere(
      cbind(first(longitude), first(latitude)),
      cbind(last(longitude), last(latitude))
    ) / 1000,
    
    # Force dtime to 1 hour so that segment speed equals segment length (km/h).
    dtime = 1,
    segment.speed = segment.length,
    
    # Calculate the segment direction from the start to the end.
    segment.dir = bearing(
      cbind(first(longitude), first(latitude)),
      cbind(last(longitude), last(latitude))
    ),
    
    # Bring in the global reference direction.
    dir.ref = global_dir.ref
  ) %>%
  ungroup() %>%
  arrange(start_date_time) %>%
  # Adjust the start coordinates for intervals that have only one point.
  mutate(
    new_start_longitude = if_else(
      (start_longitude == end_longitude & start_latitude == end_latitude & row_number() > 1),
      lag(end_longitude),
      start_longitude
    ),
    new_start_latitude = if_else(
      (start_longitude == end_longitude & start_latitude == end_latitude & row_number() > 1),
      lag(end_latitude),
      start_latitude
    ),
    # Recalculate segment length using the new start coordinates (if adjusted) and current end coordinates.
    segment.length = distVincentySphere(
      cbind(new_start_longitude, new_start_latitude),
      cbind(end_longitude, end_latitude)
    ) / 1000,
    # Recalculate segment speed (still dtime = 1 hour).
    segment.speed = segment.length,
    
    # Calculate the difference from the reference direction and the movement components.
    dir.delta2    = dir.ref - segment.dir,
    forward.speed = cos(dir.delta2 * pi/180) * segment.speed,
    perpen.speed  = -sin(dir.delta2 * pi/180) * segment.speed
  ) %>%
  # For the first record of each day, remove segment length and speed since there is no prior point.
  group_by(day = as.Date(start_date_time)) %>%
  mutate(
    segment.length = if_else(row_number() == 1, NA_real_, segment.length),
    segment.speed  = if_else(row_number() == 1, NA_real_, segment.speed),
    forward.speed  = if_else(row_number() == 1, NA_real_, forward.speed),
    perpen.speed   = if_else(row_number() == 1, NA_real_, perpen.speed)
  ) %>%
  ungroup() %>%
  select(-day)

####

# EXTRACT WIND DATA -------------------------------------------------------
# Step 1: Ensure the hour_dat has required columns
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
hour_dat$uwind <- wind_data[, 1]
hour_dat$vwind <- wind_data[, 2]

# Step 6: Calculate wind direction and speed
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


# SCATTER PLOTS TAILWIND VS FORWARD SPEED ---------------------------------

# Plot 1: Tailwind vs. Forward Speed
{
  p1 <- ggplot(hour_dat, aes(x = tailwind, y = forward.speed)) +
  geom_point(color = "firebrick", size = 2) +
  geom_smooth(method = "lm", color = "navy", se = FALSE, linewidth = 1) +
  labs(title = "Tailwind vs. Forward Speed",
       x = "Tailwind (m/s)",
       y = "Forward Speed (km/h)") +
  theme_minimal()

# Plot 2: Crosswind vs. Perpendicular Speed
p2 <- ggplot(hour_dat, aes(x = crosswind, y = perpen.speed)) +
  geom_point(color = "firebrick", size = 2) +
  geom_smooth(method = "lm", color = "navy", se = FALSE, linewidth = 1) +
  labs(title = "Crosswind vs. Perpendicular Speed",
       x = "Crosswind (m/s)",
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



# WIND DATA PLOTS PDF EXPORT & SCALE CHANGES ------------------------------

# Prepare data for plotting
map_data_ggplot1 <- hour_dat %>%
  filter(!is.na(start_latitude) & !is.na(start_longitude)) %>%
  mutate(
    date_label = as.Date(datetime),  # Ensure it's a proper date format
    hour_group = lubridate::hour(datetime),  # Extract hour
    
    # Calculate wind vector endpoints (0.01 scaling)
    wind_end_lon = start_longitude + 0.01 * wind.speed * sin(deg_to_rad(wind.dir)),
    wind_end_lat = start_latitude + 0.01 * wind.speed * cos(deg_to_rad(wind.dir)),
    
    # Calculate heading vector endpoints (using airspeed for length)
    heading_end_lon = start_longitude + 0.01 * airspeed * sin(deg_to_rad(heading)),
    heading_end_lat = start_latitude + 0.01 * airspeed * cos(deg_to_rad(heading))
  )

# Get unique days for separate plots
unique_days <- unique(map_data_ggplot1$date_label)

# Loop over each unique day and create a separate plot
plot_list <- list()

for (day in unique_days) {
  day_data <- filter(map_data_ggplot1, date_label == day)
  
  # Ensure day is a proper Date before formatting
  day_formatted <- as.character(format(as.Date(day), "%d-%b-%y"))
  
  plot_list[[day_formatted]] <- ggplot(day_data, aes(x = start_longitude, y = start_latitude)) +
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
    
    labs(
      title = paste("Migration Path on", day_formatted),  # Corrected title formatting
      x = "Longitude", y = "Latitude"
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", colour = "gray40"),
      legend.position = "right",
      strip.background = element_rect(fill = "gray40"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")  # Center title and adjust size
    ) +
    
    # Fix aspect ratio for consistent grid size
    coord_fixed(1) +
    
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
    ) 
}

# Print all plots
for (day in names(plot_list)) {
  print(plot_list[[day]])
}

# export the plots into a PDF 
# Specify the PDF file name and open the PDF device
# INSERT NAME YEAR AND SEASON OF THE BIRD 
pdf("Ronny_autumn_2014.pdf", height = 8, width = 8)

# Loop through each plot in the plot_list and print it to the PDF
for (day in names(plot_list)) {
  print(plot_list[[day]])  # Print each plot to the PDF
}

# Close the PDF device
dev.off()


## PLOTS WIT VIXED AXIS FOR NO DISTORTION 
# Set fixed span (change this if you want a different range)
fixed_span <- 6      # total degrees (CHANGE AS NEEDED)
half_span <- fixed_span / 2  # half the span

for (day in unique_days) {
  day_data <- filter(map_data_ggplot1, date_label == day)
  
  # Calculate the center of the track
  center_lat <- mean(range(day_data$start_latitude))
  center_lon <- mean(range(day_data$start_longitude))
  
  # Format the day for title
  day_formatted <- as.character(format(as.Date(day), "%d-%b-%y"))
  
  plot_list[[day_formatted]] <- ggplot(day_data, aes(x = start_longitude, y = start_latitude)) +
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
    
    labs(
      title = paste("Migration Path on", day_formatted),
      x = "Longitude", y = "Latitude"
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", colour = "gray40"),
      legend.position = "right",
      strip.background = element_rect(fill = "gray40"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    
    # Fix aspect ratio and set fixed x and y limits based on the center
    coord_fixed(xlim = c(center_lon - half_span, center_lon + half_span),
                ylim = c(center_lat - half_span, center_lat + half_span)) +
    
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
}

# Print all plots
for (day in names(plot_list)) {
  print(plot_list[[day]])
}

# export the plots into a PDF 
# Specify the PDF file name and open the PDF device
pdf("Corry_autumn_2012_1.pdf", height = 8, width = 8) # INSERT NAME YEAR AND SEASON OF THE BIRD 

# Loop through each plot in the plot_list and print it to the PDF
for (day in names(plot_list)) {
  print(plot_list[[day]])  # Print each plot to the PDF
}

# Close the PDF device
dev.off()



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
