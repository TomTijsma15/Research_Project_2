

# relationship between crosswind and perpen.speed
# slope= Δcrosswind/ Δperpendicular_speed (n_hour - n_hour-1)

hour_dat <- hour_dat %>%
  arrange(date, hour) %>%  # Ensure data is sorted by date and hour
  group_by(date) %>%
  mutate(
    prev_perpen_speed = lag(perpen.speed),  # Get previous hour's perpen.speed
    prev_crosswind = lag(crosswind),        # Get previous hour's crosswind
    delta_perpen = perpen.speed - prev_perpen_speed,  # Change in perpen speed
    delta_crosswind = crosswind - prev_crosswind,     # Change in crosswind
    slope = ifelse(!is.na(delta_crosswind) & delta_crosswind != 0, 
                   delta_perpen / delta_crosswind, 
                   NA_real_),  # Avoid division by zero
    behavior = case_when(
      slope > 2  ~ "drifting",
      slope < -2 ~ "overcompensating",
      abs(slope) <= 2 ~ "compensating",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup()

## vector analysis
hour_dat <- hour_dat %>%
  group_by(date) %>% mutate(
    # Compute the angle (in degrees) for each vector relative to the x-axis
    angle_deg = atan2(perpen.speed, crosswind) * (180 / pi),
    # Adjust angle to be within [0, 180) degrees to handle all quadrants
    angle_mod = angle_deg %% 180,
    # Classify based on the adjusted angle
    behavior = case_when(
      # Drifting: angle ~45° (slope ~1)
      abs(angle_mod - 45) <= 5.7 ~ "drifting",
      # Compensating: angle near 0 or near 180
      (angle_mod <= 5.7 | angle_mod >= 180 - 5.7) ~ "compensating",
      # Overcompensating: angle ~135° (slope ~-1)
      abs(angle_mod - 135) <= 5.7 ~ "overcompensating",
      angle_mod > 5.7 & angle_mod < 39.3 ~ "partial_drift",  # Between compensation and drift
      angle_mod > 50.7 & angle_mod < 129.3 ~ "overdrift",     # Between drift and overcompensation
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(behavior)) %>% # Remove rows with NA behavior
  mutate(behavior = factor(behavior)) %>%
  ungroup()



# BAR PLOT THE BEHAVIORS
ggplot(hour_dat, aes(x = date, fill = behavior)) +
  geom_bar(position = "stack") +
  labs(title = "Daily Distribution of Bird Behaviors",
       x = "Date",
       y = "Count of Behaviors",
       fill = "Behavior") +
  scale_fill_manual(name = "Behavior",
                     values = c("drifting"         = "#4daf4a",   # clear green
                                "compensating"     = "#984ea3",   # strong purple
                                "overcompensating" = "#e41a1c",   # bold red
                                "partial_drift"    = "#FF69B4",   # hot pink
                                "overdrift"        = "#fdae61")) +   # soft orange
  theme_minimal()

# scatterplot perpen.speed Crosswind and tailwind, collored for behavior
ggplot(hour_dat, aes(x = crosswind, y = perpen.speed, color = behavior)) +
  geom_point(size = 2) +
  # Dotted lines for behavior boundaries
  geom_abline(slope = tan(5.7 * pi/180), intercept = 0, linetype = "dotted", color = "black") +
  geom_abline(slope = tan(-5.7 * pi/180), intercept = 0, linetype = "dotted", color = "black") +
  geom_abline(slope = tan(39.3 * pi/180), intercept = 0, linetype = "dotted", color = "black") +
  geom_abline(slope = tan(50.7 * pi/180), intercept = 0, linetype = "dotted", color = "black") +
  geom_abline(slope = tan(129.3 * pi/180), intercept = 0, linetype = "dotted", color = "black") +
  geom_abline(slope = tan(140.7 * pi/180), intercept = 0, linetype = "dotted", color = "black") +
  coord_equal() +
  labs(title = "Forward Speed vs. Crosswind",
       x = "Crosswind (km/h)",
       y = "Forward Speed (km/h)",
       color = "Behavior") +
  scale_color_manual(name = "Behavior",
                     values = c("drifting"         = "#4daf4a",   # clear green
                                "compensating"     = "#984ea3",   # strong purple
                                "overcompensating" = "#e41a1c",   # bold red
                                "partial_drift"    = "#FF69B4",   # hot pink
                                "overdrift"        = "#fdae61")) +   # soft orange
  theme_minimal()





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



# BEHAVIOR ON TRACK PLOTS -------------------------------------------------
{
  # --- 1. Prepare the Data ---
  map_data <- hour_dat %>%
    filter(!is.na(start_latitude) & !is.na(start_longitude)) %>%
    mutate(
      date_label = as.Date(datetime),
      hour_group = hour(datetime)  # still computed if needed elsewhere
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
  
  # --- 7. Create a Plot for Each Day (with a 1° Margin) ---
  unique_days <- unique(map_data$date_label)
  plot_list <- list()
  
}
for(day in unique_days) {
  # Filter and sort the day's data by datetime (to preserve track order)
  day_data <- map_data %>% filter(date_label == day) %>% arrange(datetime)
  
  # Create segments by computing the end points for each row
  day_data_seg <- day_data %>%
    mutate(easting_end = lead(easting),
           northing_end = lead(northing)) %>%
    filter(!is.na(easting_end) & !is.na(northing_end))
  
  # -- Calculate Expanded Bounding Box in WGS84 --
  # Use original geographic coordinates to set a 1° margin
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
    
    # Plot the track as segments, colored by behavior
    geom_segment(data = day_data_seg,
                 aes(x = easting, y = northing,
                     xend = easting_end, yend = northing_end,
                     color = behavior),
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
    scale_color_manual(name = "Behavior",
                       values = c("drifting"         = "#4daf4a",   # clear green
                                  "compensating"     = "#984ea3",   # strong purple
                                  "overcompensating" = "#e41a1c",   # bold red
                                  "partial_drift"    = "#FF69B4",   # hot pink
                                  "overdrift"        = "#fdae61")) +   # soft orange
    scale_alpha_continuous(name = "Wind Speed (km/h)", range = c(0.3, 0.8)) +
    scale_linewidth_continuous(name = "Airspeed (km/h)", range = c(0.3, 0.8)) +
    labs(title = paste("Migration Track on", format(as.Date(day), "%d-%b-%Y")),
         x = "Easting (m)", y = "Northing (m)") +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) + 
    
    # Add custom north arrow (using annotation_custom as before)
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
