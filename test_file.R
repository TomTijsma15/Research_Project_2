
# Create sample data
set.seed(123)
sample_data <- tibble(
  crosswind = runif(250, -30, 30),
  perpen.speed = runif(250, -30, 30)
) %>%
  mutate(
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
  )

# Create visualization
ggplot(sample_data, aes(x = crosswind, y = perpen.speed, color = behavior)) +
  geom_point() +
  # Dotted lines for behavior boundaries
  geom_abline(slope = tan(5.7 * pi/180), intercept = 0, linetype = "dotted", color = "black") +
  geom_abline(slope = tan(-5.7 * pi/180), intercept = 0, linetype = "dotted", color = "black") +
  geom_abline(slope = tan(39.3 * pi/180), intercept = 0, linetype = "dotted", color = "black") +
  geom_abline(slope = tan(50.7 * pi/180), intercept = 0, linetype = "dotted", color = "black") +
  geom_abline(slope = tan(129.3 * pi/180), intercept = 0, linetype = "dotted", color = "black") +
  geom_abline(slope = tan(140.7 * pi/180), intercept = 0, linetype = "dotted", color = "black") +
  coord_equal() +
  labs(title = "Behavior Classification Based on Wind Compensation",
       x = "Crosswind",
       y = "Perpendicular Speed") +
  scale_color_manual(name = "Behavior",
                     values = c("drifting"         = "#4daf4a",   # clear green
                                "compensating"     = "#984ea3",   # strong purple
                                "overcompensating" = "#e41a1c",   # bold red
                                "partial_drift"    = "#FF69B4",   # hot pink
                                "overdrift"        = "#fdae61")) +  # soft orange
  
  theme_minimal()


# filter for % of fly instead of >5 fly  ----------------------------------

fly_data2 <- data %>%
  # Ensure proper datetime format and compute time difference
  mutate(
    date_time = as.POSIXct(date_time, format = "%Y-%m-%d %H:%M:%S"),
    time_diff = as.numeric(difftime(date_time, min(date_time), units = "mins"))
  ) %>%
  # compute the hourly grouping based on time_diff
  mutate(custom_hour_group = floor(time_diff / 60)) %>%
  # Group by the custom hourly intervals
  group_by(custom_hour_group) %>%
  # Calculate the percentage of 'fly' entries in each hour
  mutate(fly_percentage = mean(speed.class == "fly") * 100) %>%
  # Keep only hours with >=75% 'fly' and filter to 'fly' entries
  filter(fly_percentage >= 75 & speed.class == "fly") %>%
  ungroup() %>%
  # Remove helper columns
  select(-time_diff, -custom_hour_group, -fly_percentage)

#Step 3: Get global start/end reference points from the ENTIRE dataset (not grouped_data)
locs.start.ref <- cbind(data$longitude[1], data$latitude[1])
locs.end.ref <- cbind(data$longitude[nrow(data)], data$latitude[nrow(data)])

#Step 4: Calculate global reference direction (same for all rows)
global_dir.ref <- bearingRhumb(locs.start.ref, locs.end.ref)

# Ensure your date_time column is correctly converted to POSIXct.
fly_data2 <- fly_data2 %>% 
  mutate(date_time = as.POSIXct(date_time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))

# Define the start and end times for your data.
start_time <- min(fly_data2$date_time, na.rm = TRUE)
end_time   <- max(fly_data2$date_time, na.rm = TRUE)

# Create breaks for 60-minute intervals starting from the first timestamp.
breaks <- seq(from = start_time, to = end_time, by = "60 min")

# Create a new data frame with an "interval" column indicating each 60-minute period.
data_60min1 <- fly_data2 %>%
  mutate(interval = cut(date_time, breaks = breaks, right = FALSE))

# Summarise the data by 60-minute intervals.
hour_dat2 <- data_60min1 %>%
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
