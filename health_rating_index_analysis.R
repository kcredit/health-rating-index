################################################################################
# Health Rating Index Analysis for Dublin, Ireland
# Authors: Credit, K., Damanpreet, K., and Eccles, E.
# Citation: Credit et al. 2025. "Exploring the transport-health-environment 
# nexus through a new 'Health Rating Index' for Dublin, Ireland." 
# Proceedings of the 33rd GISRUK Conference. 
# DOI: https://doi.org/10.5281/zenodo.15183740
################################################################################

# This script analyzes environmental health burdens and benefits to create
# a Health Rating Index (HRI) for small areas in Dublin, Ireland.
# Key components:
# 1. Multiple GP accessibility metrics (E2SFCA with different decay functions)
# 2. Blue space accessibility calculation
# 3. Green space accessibility with Gaussian decay
# 4. Years of Life Lost (YLL) calculations for air pollution and noise
# 5. Multiple index weighting approaches comparison
# 6. Spatial Cross-Validation Random Forest with comparison to spatial models
# 7. ALE and importance plots with confidence intervals
# 8. Interactive maps of results

################################################################################
# SETUP: Load Libraries
################################################################################

rm(list=ls())

# Required packages
packages.wanted <- c("stats", "sf", "vip", "spdep", "tidymodels", "vip", "pdp", 
                     "gridExtra", "magrittr", "randomForest", "conflicted", 
                     "finalfit", "leafgl", "mapview", "data.table", "s2", "sf", 
                     "r5r", "GGally", "MASS", "robustbase", "mvoutlier", 
                     "RColorBrewer", "rgdal", "ggplot2", "Hmisc", "purrr", 
                     "foreign", "stargazer", "tidyverse", "devtools", "otpr", 
                     "httr", "preogress", "sandwich", "lmtest", "googletraffic", 
                     "raster", "ranger", "mgcv", "xgboost", "e1071", "hydroGOF", 
                     "ALEPlot", "spatialreg", "blockCV", "leaflet", "leaflet.extras")

for (package in packages.wanted) require(package, character.only=TRUE)

# Resolve function conflicts
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("summarize", "dplyr")
conflicted::conflicts_prefer(pdp::partial)
conflicted::conflicts_prefer(vip::vi)

# Set working directory (update this to your local path)
setwd("~/Dropbox/Grants/New Foundations/Smart D8 Dashboard/Anaysis Data/health-rating-index") #Change for user directoy

################################################################################
# LOAD DATA
################################################################################

cat("Loading data files...\n")

sa <- read.csv("data/SA_Origin.csv")
gp <- read.csv("data/GP_Dest.csv")
pk <- read.csv("data/Park_Dest.csv")
bs <- read.csv("data/Blue_Dest.csv")

# Pre-calculated travel time matrices
dist <- read.csv("data/dist_sa_gp.csv")
distpk <- read.csv("data/dist_sa_pk1.csv")
disbs <- read.csv("data/dist_sa_bs.csv")

################################################################################
# OPTIONAL: Calculate Travel Time Matrices Using r5r
################################################################################
# NOTE: This section is commented out because it requires significant 
# computational time. Pre-calculated matrices are provided in the data folder.
# Uncomment and modify if you want to recalculate travel times.

# data_path <- "./data/dublinpoa"
# o_sa <- fread(file.path(data_path, "SA_Origin.csv"))
# d_gp <- fread(file.path(data_path, "GP_Dest.csv"))
# d_pk <- fread(file.path(data_path, "Park_Dest_Points.csv"))
# d_bs <- fread(file.path(data_path, "Blue_Dest.csv"))
# 
# options(java.parameters = "-Xmx10G")
# r5r_core <- setup_r5(data_path = data_path, verbose = FALSE)
# 
# walk_speed <- 5
# bike_speed <- 18
# max_trip_duration <- 45
# departure_datetime <- as.POSIXct("29-11-2022 09:00:00", 
#                                  format = "%d-%m-%Y %H:%M:%S")
# mode <- c("CAR")
# 
# dist <- travel_time_matrix(r5r_core = r5r_core,
#                            origins = o_sa, destinations = d_gp,
#                            mode = mode, departure_datetime = departure_datetime,
#                            time_window = 1, max_trip_duration = max_trip_duration,
#                            walk_speed = walk_speed, bike_speed = bike_speed,
#                            verbose = FALSE)
# 
# distbs <- travel_time_matrix(r5r_core = r5r_core,
#                              origins = o_sa, destinations = d_bs,
#                              mode = mode, departure_datetime = departure_datetime,
#                              time_window = 1, max_trip_duration = max_trip_duration,
#                              walk_speed = walk_speed, bike_speed = bike_speed,
#                              verbose = FALSE)
# 
# max_trip_duration <- 15
# distpk <- travel_time_matrix(r5r_core = r5r_core,
#                              origins = o_sa, destinations = d_pk,
#                              mode = mode, departure_datetime = departure_datetime,
#                              time_window = 1, max_trip_duration = max_trip_duration,
#                              walk_speed = walk_speed, bike_speed = bike_speed,
#                              verbose = FALSE)
# 
# stop_r5(r5r_core)
# rJava::.jgc(R.gc = TRUE)

################################################################################
# KEY COMPONENT 1: Enhanced Two-Step Floating Catchment Area (E2SFCA) 
# Multiple GP Accessibility Metrics
################################################################################

cat("\n=== CALCULATING GP ACCESSIBILITY METRICS ===\n")

# Distance decay function for accessibility calculations
calculate_decay_weight <- function(travel_time, method = "multi_zone") {
  
  if (method == "multi_zone") {
    # Three-zone decay (Luo & Qi 2009)
    weight <- case_when(
      travel_time <= 10 ~ 1.0,
      travel_time <= 20 ~ 0.68,
      travel_time <= 30 ~ 0.22,
      TRUE ~ 0
    )
  } else if (method == "gaussian") {
    # Gaussian decay
    d0 <- 20
    weight <- exp(-(travel_time^2) / (2*d0^2))
    weight <- ifelse(travel_time > 30, 0, weight)
  } else if (method == "kernel") {
    # Kernel density function
    d0 <- 30
    weight <- ifelse(travel_time < d0,
                     (3/4) * (1 - (travel_time/d0)^2),
                     0)
  } else {
    # Binary (traditional 2SFCA)
    weight <- ifelse(travel_time <= 30, 1, 0)
  }
  
  return(weight)
}

# Step 1: Calculate provider-to-population ratios
step1_provider_ratios <- function(travel_time_matrix, demand_data, supply_data, 
                                  max_travel_time = 30, decay_method = "multi_zone") {
  
  ttm_filtered <- travel_time_matrix %>%
    filter(travel_time <= max_travel_time)
  
  ttm_weighted <- ttm_filtered %>%
    mutate(weight = calculate_decay_weight(travel_time, method = decay_method))
  
  ttm_with_demand <- ttm_weighted %>%
    left_join(demand_data, by = c("origin_id" = "origin_id"))
  
  weighted_demand <- ttm_with_demand %>%
    group_by(destination_id) %>%
    summarise(
      total_weighted_demand = sum(population * weight, na.rm = TRUE),
      n_areas_in_catchment = n()
    ) %>%
    ungroup()
  
  provider_ratios <- supply_data %>%
    left_join(weighted_demand, by = "destination_id") %>%
    mutate(
      ratio = ifelse(total_weighted_demand > 0, 
                     capacity / total_weighted_demand, 
                     0)
    )
  
  return(provider_ratios)
}

# Step 2: Calculate accessibility scores
step2_accessibility <- function(travel_time_matrix, provider_ratios, 
                                max_travel_time = 30, decay_method = "multi_zone") {
  
  ttm_filtered <- travel_time_matrix %>%
    filter(travel_time <= max_travel_time)
  
  ttm_weighted <- ttm_filtered %>%
    mutate(weight = calculate_decay_weight(travel_time, method = decay_method))
  
  ttm_with_ratios <- ttm_weighted %>%
    left_join(
      provider_ratios %>% select(destination_id, ratio),
      by = "destination_id"
    )
  
  accessibility_scores <- ttm_with_ratios %>%
    group_by(origin_id) %>%
    summarise(
      accessibility = sum(ratio * weight, na.rm = TRUE),
      n_providers_in_catchment = n()
    ) %>%
    ungroup()
  
  return(accessibility_scores)
}

# Complete E2SFCA function
e2sfca_analysis <- function(travel_time_matrix, demand_data, supply_data,
                            max_travel_time = 30, decay_method = "multi_zone") {
  
  cat("  Calculating provider ratios...\n")
  provider_ratios <- step1_provider_ratios(
    travel_time_matrix, demand_data, supply_data,
    max_travel_time, decay_method
  )
  
  cat("  Calculating accessibility scores...\n")
  accessibility_scores <- step2_accessibility(
    travel_time_matrix, provider_ratios,
    max_travel_time, decay_method
  )
  
  final_results <- accessibility_scores %>%
    left_join(demand_data, by = "origin_id")
  
  return(list(
    accessibility = final_results,
    provider_ratios = provider_ratios
  ))
}

# Prepare data for E2SFCA
demand_data <- sa %>%
  select(id, T1_1AGETT) %>%
  rename(origin_id = id, population = T1_1AGETT)
demand_data$origin_id <- as.character(demand_data$origin_id)

supply_data <- gp %>%
  select(id, SUPPLY) %>%
  rename(destination_id = id, capacity = SUPPLY)
supply_data$destination_id <- as.character(supply_data$destination_id)

dist <- dist %>%
  rename(origin_id = from_id, destination_id = to_id, travel_time = travel_time_p50)
dist$origin_id <- as.character(dist$origin_id)
dist$destination_id <- as.character(dist$destination_id)

# Calculate accessibility using multiple decay methods
decay_methods <- c("multi_zone", "gaussian", "kernel", "binary")
gp_access_results <- list()

for (method in decay_methods) {
  cat(paste0("\nCalculating GP accessibility with ", method, " decay...\n"))
  results <- e2sfca_analysis(
    travel_time_matrix = dist,
    demand_data = demand_data,
    supply_data = supply_data,
    max_travel_time = 30,
    decay_method = method
  )
  
  suffix <- switch(method,
                   "multi_zone" = "mz",
                   "gaussian" = "gaus",
                   "kernel" = "kern",
                   "binary" = "bin")
  
  gp_access_results[[method]] <- results$accessibility %>%
    select(origin_id, accessibility) %>%
    rename(id = origin_id) %>%
    rename_with(~paste0(suffix, "_2s"), accessibility)
}

# Merge all GP accessibility metrics into sa dataframe
for (method in decay_methods) {
  sa <- merge(sa, gp_access_results[[method]], by = "id", all.x = TRUE)
}

cat("\nGP accessibility metrics calculated successfully.\n")

################################################################################
# KEY COMPONENT 2: Blue Space Accessibility
################################################################################

cat("\n=== CALCULATING BLUE SPACE ACCESSIBILITY ===\n")

# Find nearest coastline point for each small area
nearestbs <- disbs %>%
  group_by(from_id) %>%
  slice_min(travel_time_p50, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(id = from_id)

# Apply Gaussian decay
sigma <- 20
nearestbs$bs_access <- exp(-(nearestbs$travel_time_p50^2) / (2*sigma^2))

nearestbs <- nearestbs %>%
  select(id, bs_access)

sa <- merge(sa, nearestbs, by = "id", all.x = TRUE)

cat("Blue space accessibility calculated successfully.\n")

################################################################################
# KEY COMPONENT 3: Green Space Accessibility
################################################################################

cat("\n=== CALCULATING GREEN SPACE ACCESSIBILITY ===\n")

# Load park polygons and origin points
o_sa <- read.csv("data/SA_Origin.csv")
parks_poly <- st_read("data/DCC Parks.shp", quiet = TRUE)

df_sf <- st_as_sf(o_sa, coords = c("lon", "lat"), crs = 4326)
parks_poly <- st_transform(parks_poly, 4326)
parks_poly <- st_make_valid(parks_poly)

# Assign maximum accessibility (5) to points within parks
within_park <- lengths(st_intersects(df_sf, parks_poly)) > 0
df_sf$gs_access <- NA
df_sf$gs_access[within_park] <- 5

# Calculate accessibility for points outside parks
outside_ids <- df_sf$id[!within_park]

if (length(outside_ids) > 0) {
  distpk_outside <- distpk %>%
    filter(from_id %in% outside_ids)
  
  pk <- pk %>%
    select(id, Name)
  pk$id <- as.character(pk$id)
  
  distpk_outside <- merge(distpk_outside, pk, by.x = "to_id", by.y = "id")
  distpk_outside$Name[distpk_outside$Name == ""] <- "NA"
  
  # Find minimum time to each park
  min_times <- distpk_outside %>%
    group_by(from_id, Name) %>%
    summarise(min_time = min(travel_time_p50), .groups = 'drop')
  
  # Find the 5 nearest parks per origin
  nearest_5 <- min_times %>%
    group_by(from_id) %>%
    arrange(min_time) %>%
    slice_head(n = 5) %>%
    ungroup()
  
  # Apply Gaussian decay
  sigma <- 5
  nearest_5 <- nearest_5 %>%
    mutate(decayed = exp(-(min_time^2) / (2 * sigma^2)))
  
  # Sum decayed values per origin
  gs_access_calculated <- nearest_5 %>%
    group_by(from_id) %>%
    summarise(gs_access = sum(decayed), .groups = 'drop')
  
  # Merge calculated accessibility
  df_sf <- df_sf %>%
    left_join(gs_access_calculated, by = c("id" = "from_id"))
  
  # Combine: keep 5 for parks, use calculated for others
  df_sf <- df_sf %>%
    mutate(gs_access = coalesce(gs_access.x, gs_access.y)) %>%
    select(-gs_access.x, -gs_access.y)
}

o_sa <- st_drop_geometry(df_sf)
sa <- merge(sa, o_sa %>% select(id, gs_access), by = "id", all.x = TRUE)

cat("Green space accessibility calculated successfully.\n")

#Save accessibility results
sa_access <- sa %>%
  select(id, SA_PUB202, mz_2s, gaus_2s, kern_2s, bin_2s, bs_access, gs_access)

################################################################################
# KEY COMPONENT 4: Years of Life Lost (YLL) Calculations
################################################################################

cat("\n=== CALCULATING YEARS OF LIFE LOST ===\n")

# Load noise and air quality data
noise <- read.csv("data/NoiseContours_2011_Dissolved_SA2022_Intersection_Area.csv")
air <- read.csv("data/SA2022_Dublin_Postcode_greenR_AirView_Noise_GS_distroad_GP_NO_TUNNEL_airviewmedians.csv")

## AIR POLLUTION YLL CALCULATION ##
cat("\nCalculating air pollution-related years of life lost...\n")

# Calculate population aged 30+
air$POP30P <- air$T1_1AGE3_1 + air$T1_1AGE3_2 + air$T1_1AGE4_1 + 
  air$T1_1AGE4_2 + air$T1_1AGE5_1 + air$T1_1AGE5_2 + 
  air$T1_1AGE6_1 + air$T1_1AGE6_2 + air$T1_1AGE7_1 + 
  air$T1_1AGE7_2 + air$T1_1AGE8_1 + air$T1_1AGEG_1

# Relative risks using WHO guidelines (WHO 2022)
air$RRPM25 <- ifelse(!is.na(air$PM25_media) & air$PM25_media >= 5,
                     (((air$PM25_media - 5) / 10) * 0.08) + 1, 1)

air$RRNO2 <- ifelse(!is.na(air$NO2_median) & air$NO2_median >= 10,
                    (((air$NO2_median - 10) / 10) * 0.02) + 1, 1)

air$RRO3 <- ifelse(!is.na(air$O3_median) & air$O3_median >= 71.05,
                   (((air$O3_median - 71.05) / 10) * 0.0043) + 1, 1)

# Population attributable fractions
air$PAFPM25 <- (air$RRPM25-1) / air$RRPM25
air$PAFNO2 <- (air$RRNO2-1) / air$RRNO2
air$PAFO3 <- (air$RRO3-1) / air$RRO3

# All-cause mortality rates (from Irish vital statistics)
CDRAll17 <- 0.005917
CDRAll17_f <- 0.0058954
CDRAll17_m <- 0.0059396
CDR3017 <- 0.009650213
CDR3017_f <- 0.009460165
CDR3017_m <- 0.009858216

## Calculate life expectancy for exposed populations ##

# Age midpoints for calculation
age_midpoints_all <- c(.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 
                       10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 
                       19.5, 22, 27, 32, 37, 42, 47, 52, 57, 62, 67, 72, 77, 
                       82, 87)
age_midpoints_30P <- c(32, 37, 42, 47, 52, 57, 62, 67, 72, 77, 82, 87)

# Define age columns
start_col_all_f <- "T1_1AGE0F"
start_col_30P_f <- "T1_1AGE301"
end_col_f <- "T1_1AGEGE1"

start_col_all_m <- "T1_1AGE0M"
start_col_30P_m <- "T1_1AGE30_"
end_col_m <- "T1_1AGEGE_"

age_cols_all_f <- names(air)[which(names(air) == start_col_all_f):
                                which(names(air) == end_col_f)]
age_cols_30P_f <- names(air)[which(names(air) == start_col_30P_f):
                                which(names(air) == end_col_f)]
age_cols_all_m <- names(air)[which(names(air) == start_col_all_m):
                                which(names(air) == end_col_m)]
age_cols_30P_m <- names(air)[which(names(air) == start_col_30P_m):
                                which(names(air) == end_col_m)]

# Weight ages by population
age_weighted_all_f <- air[, age_cols_all_f] * 
                      rep(age_midpoints_all, each = nrow(air))
age_weighted_30P_f <- air[, age_cols_30P_f] * 
                      rep(age_midpoints_30P, each = nrow(air))
age_weighted_all_m <- air[, age_cols_all_m] * 
                      rep(age_midpoints_all, each = nrow(air))
age_weighted_30P_m <- air[, age_cols_30P_m] * 
                      rep(age_midpoints_30P, each = nrow(air))

weighted_sum_all_f <- rowSums(age_weighted_all_f, na.rm = TRUE)
weighted_sum_30P_f <- rowSums(age_weighted_30P_f, na.rm = TRUE)
weighted_sum_all_m <- rowSums(age_weighted_all_m, na.rm = TRUE)
weighted_sum_30P_m <- rowSums(age_weighted_30P_m, na.rm = TRUE)

air$total_pop_all_f <- rowSums(air[, age_cols_all_f], na.rm = TRUE)
air$total_pop_30P_f <- rowSums(air[, age_cols_30P_f], na.rm = TRUE)
air$total_pop_all_m <- rowSums(air[, age_cols_all_m], na.rm = TRUE)
air$total_pop_30P_m <- rowSums(air[, age_cols_30P_m], na.rm = TRUE)

# Calculate average age by sex and age group
air$avg_age_all_f <- weighted_sum_all_f / air$total_pop_all_f
air$avg_age_30P_f <- weighted_sum_30P_f / air$total_pop_30P_f
air$avg_age_all_m <- weighted_sum_all_m / air$total_pop_all_m
air$avg_age_30P_m <- weighted_sum_30P_m / air$total_pop_30P_m

# Load Irish life table
life_table <- read.csv("data/irish_life_table.csv")
life_table <- life_table[, c("Age", "Sex", "LifeExpectancy")]

# Function to interpolate life expectancy
get_life_expectancy <- function(avg_age, sex, life_table) {
  if (is.na(avg_age)) return(NA)
  
  lt_sex <- life_table[life_table$Sex == sex, ]
  lt_sex$Age <- as.numeric(lt_sex$Age)
  lt_sex$LifeExpectancy <- as.numeric(lt_sex$LifeExpectancy)
  
  if (avg_age <= min(lt_sex$Age, na.rm = TRUE)) {
    return(lt_sex$LifeExpectancy[which.min(lt_sex$Age)])
  }
  if (avg_age >= max(lt_sex$Age, na.rm = TRUE)) {
    return(lt_sex$LifeExpectancy[which.max(lt_sex$Age)])
  }
  
  lower_idx <- max(which(lt_sex$Age <= avg_age))
  upper_idx <- min(which(lt_sex$Age >= avg_age))
  
  age_lower <- lt_sex$Age[lower_idx]
  age_upper <- lt_sex$Age[upper_idx]
  le_lower <- lt_sex$LifeExpectancy[lower_idx]
  le_upper <- lt_sex$LifeExpectancy[upper_idx]
  
  le_interp <- le_lower + (le_upper - le_lower) * 
               (avg_age - age_lower) / (age_upper - age_lower)
  
  return(le_interp)
}

# Calculate life expectancy for each small area
air$LE_all_m <- sapply(air$avg_age_all_m, get_life_expectancy, 
                       sex="Male", life_table=life_table)
air$LE_all_f <- sapply(air$avg_age_all_f, get_life_expectancy, 
                       sex="Female", life_table=life_table)
air$LE_30P_m <- sapply(air$avg_age_30P_m, get_life_expectancy, 
                       sex="Male", life_table=life_table)
air$LE_30P_f <- sapply(air$avg_age_30P_f, get_life_expectancy, 
                       sex="Female", life_table=life_table)

# Calculate Years of Life Lost (YLL)
air$YLLPM25 <- (air$PAFPM25 * air$total_pop_30P_f * CDR3017_f * air$LE_30P_f) +
               (air$PAFPM25 * air$total_pop_30P_m * CDR3017_m * air$LE_30P_m)

air$YLLNO2 <- (air$PAFNO2 * air$total_pop_30P_f * CDR3017_f * air$LE_30P_f) +
              (air$PAFNO2 * air$total_pop_30P_m * CDR3017_m * air$LE_30P_m)

air$YLLO3 <- (air$PAFO3 * air$total_pop_all_f * CDRAll17_f * air$LE_all_f) +
             (air$PAFO3 * air$total_pop_all_m * CDRAll17_m * air$LE_all_m)

air$YLLAQ <- air$YLLPM25 + air$YLLNO2 + air$YLLO3

# Simplify air data
air <- air %>%
  select(SA_PUB2022, POP30P, SHAPE_Area, T1_1AGETT, total_pop_all_f, total_pop_all_m, LE_all_m, 
         LE_all_f, LE_30P_m, LE_30P_f,YLLPM25, YLLNO2, YLLO3, YLLAQ, 
         T10_4_OD_2, T10_4_HD_2, T10_4_PDT,T10_4_DT, T10_4_TT, T11_1_FT, 
         T11_1_BIT,T11_1_BUT, T11_1_TDLT,T11_1_TT, T12_3_BT, T12_3_VBT, T12_3_TT)

cat("  Total YLL from air pollution:", round(sum(air$YLLAQ, na.rm=TRUE), 2), "\n")

## NOISE YLL CALCULATION ##
cat("\nCalculating noise-related years of life lost...\n")

noise_d <- noise %>%
  spread(DBs, area) %>%
  select(SA_PUB2022, `55`, `60`, `65`, `70`, `75`)

noise_w <- noise_d %>%
  group_by(SA_PUB2022) %>% 
  summarise(across(everything(), ~ reduce(., coalesce)), .groups = 'drop')

air <- merge(air, noise_w, by = "SA_PUB2022", all.x = TRUE)

# Calculate proportion of area exposed to each noise level
air$D55_P <- air$`55`/air$SHAPE_Area
air$D60_P <- air$`60`/air$SHAPE_Area
air$D65_P <- air$`65`/air$SHAPE_Area
air$D70_P <- air$`70`/air$SHAPE_Area
air$D75_P <- air$`75`/air$SHAPE_Area

# Replace NA with 0
air$D55_P[is.na(air$D55_P)] <- 0
air$D60_P[is.na(air$D60_P)] <- 0
air$D65_P[is.na(air$D65_P)] <- 0
air$D70_P[is.na(air$D70_P)] <- 0
air$D75_P[is.na(air$D75_P)] <- 0

# Odds ratios for noise levels (WHO 2022)
air$OR55 <- ifelse(air$D55_P==0, 1, 1.000195)
air$OR60 <- ifelse(air$D60_P==0, 1, 1.01296)
air$OR65 <- ifelse(air$D65_P==0, 1, 1.061315)
air$OR70 <- ifelse(air$D70_P==0, 1, 1.15078)
air$OR75 <- ifelse(air$D75_P==0, 1, 1.286875)

# Background IHD rates (from Irish vital statistics)
Ik <- 0.000654
Ik_f <- 0.000449
Ik_m <- 0.000868

# Noise YLL
air$YLL55 <- ((air$OR55-1) * air$total_pop_all_f * Ik_f * air$D55_P * air$LE_all_f) +
             ((air$OR55-1) * air$total_pop_all_m * Ik_m * air$D55_P * air$LE_all_m)
air$YLL60 <- ((air$OR60-1) * air$total_pop_all_f * Ik_f * air$D60_P * air$LE_all_f) +
             ((air$OR60-1) * air$total_pop_all_m * Ik_m * air$D60_P * air$LE_all_m)
air$YLL65 <- ((air$OR65-1) * air$total_pop_all_f * Ik_f * air$D65_P * air$LE_all_f) +
             ((air$OR65-1) * air$total_pop_all_m * Ik_m * air$D65_P * air$LE_all_m)
air$YLL70 <- ((air$OR70-1) * air$total_pop_all_f * Ik_f * air$D70_P * air$LE_all_f) +
             ((air$OR70-1) * air$total_pop_all_m * Ik_m * air$D70_P * air$LE_all_m)
air$YLL75 <- ((air$OR75-1) * air$total_pop_all_f * Ik_f * air$D75_P * air$LE_all_f) +
             ((air$OR75-1) * air$total_pop_all_m * Ik_m * air$D75_P * air$LE_all_m)

air$YLLN <- air$YLL55 + air$YLL60 + air$YLL65 + air$YLL70 + air$YLL75

cat("  Total YLL from noise pollution:", round(sum(air$YLLN, na.rm=TRUE), 2), "\n")

# Replace NA with 0
air$YLLAQ[is.na(air$YLLAQ)] <- 0
air$YLLN[is.na(air$YLLN)] <- 0

# OUTPUT: Save total YLL
cat("\n=== OUTPUT 1: TOTAL YEARS OF LIFE LOST ===\n")
total_yll <- data.frame(
  Source = c("Air Pollution (PM2.5, NO2, O3)", "Road Noise (55-75+ dB)"),
  Total_YLL = c(sum(air$YLLAQ, na.rm=TRUE), sum(air$YLLN, na.rm=TRUE))
)
print(total_yll)
write.csv(total_yll, "output/total_years_life_lost.csv", row.names = FALSE)

################################################################################
# KEY COMPONENT 5: Multiple Index Weighting Approaches
################################################################################

cat("\n=== TESTING MULTIPLE INDEX WEIGHTING APPROACHES ===\n")

# Load spatial data
SA <- st_read("data/SA2022_Dublin_AllData3.shp", quiet = TRUE) %>%
  select(SA_PUB2022, RoutingKey, In22_ED, SOUTH, log_dist)

SA <- SA %>%
  left_join(sa_access, by = c("SA_PUB2022" = "SA_PUB202")) %>%
  left_join(air, by = "SA_PUB2022")

summary(SA[, c("YLLAQ", "YLLN")])

# Create rate variables
SA$YLLAQPC <- (SA$YLLAQ / SA$T1_1AGETT) * 100000
SA$YLLNPC <- (SA$YLLN / SA$T1_1AGETT) * 100000

# Function to calculate naive weights
calc_naive <- function(df, gp_var) {
  df_data <- st_drop_geometry(df)
  df_scaled <- df_data %>%
    select(all_of(c(gp_var, "gs_access", "bs_access", "YLLAQPC", "YLLNPC"))) %>%
    mutate(across(everything(), ~scale(.)[,1]))
  
  hri <- (1/4) * df_scaled[[gp_var]] + 
         (1/8) * df_scaled$gs_access + 
         (1/8) * df_scaled$bs_access - 
         (1/4) * df_scaled$YLLAQPC - 
         (1/4) * df_scaled$YLLNPC
  
  (hri - min(hri)) / (max(hri) - min(hri))
}

# Function to calculate entropy weights
calc_entropy <- function(df, gp_var) {
  df_data <- st_drop_geometry(df)
  df_norm <- df_data %>%
    select(all_of(c(gp_var, "gs_access", "bs_access", "YLLAQPC", "YLLNPC"))) %>%
    mutate(across(everything(), ~ (. - min(.)) / (max(.) - min(.)))) %>%
    mutate(YLLAQPC = 1 - YLLAQPC, YLLNPC = 1 - YLLNPC)
  
  p <- df_norm / rowSums(df_norm)
  n <- nrow(df_norm)
  k <- 1 / log(n)
  p[p == 0] <- 1e-10
  
  e <- -k * colSums(p * log(p))
  d <- 1 - e
  w <- d / sum(d)
  
  index <- as.numeric(as.matrix(df_norm) %*% w)
  (index - min(index)) / (max(index) - min(index))
}

# Function to calculate PCA weights
calc_pca <- function(df, gp_var) {
  df_data <- st_drop_geometry(df)
  df_pca <- df_data %>%
    select(all_of(c(gp_var, "gs_access", "bs_access", "YLLAQPC", "YLLNPC"))) %>%
    mutate(YLLAQPC = max(YLLAQPC, na.rm = TRUE) - YLLAQPC,
           YLLNPC = max(YLLNPC, na.rm = TRUE) - YLLNPC) %>%
    na.omit()
  
  row_ids <- as.numeric(rownames(df_pca))
  df_pca_scaled <- as.data.frame(scale(df_pca))
  
  pca_res <- prcomp(df_pca_scaled, center = FALSE, scale. = FALSE)
  eig_vals <- (pca_res$sdev)^2
  prop_variance <- eig_vals / sum(eig_vals)
  retain <- which(eig_vals > 1)
  
  variable_weights <- matrix(0, nrow = ncol(df_pca_scaled), ncol = 1)
  for (i in seq_along(retain)) {
    pc_idx <- retain[i]
    variable_weights <- variable_weights + 
      (pca_res$rotation[, pc_idx]^2) * prop_variance[pc_idx]
  }
  variable_weights <- variable_weights / sum(variable_weights)
  
  hri <- as.numeric(as.matrix(df_pca_scaled) %*% variable_weights)
  hri_scaled <- (hri - min(hri)) / (max(hri) - min(hri))
  
  result <- data.frame(row_id = row_ids, hri = hri_scaled)
  df_data <- st_drop_geometry(df)
  df_data %>%
    mutate(row_id = row_number()) %>%
    left_join(result, by = "row_id") %>%
    pull(hri)
}

# Generate all 12 indices (4 GP methods × 3 weighting approaches)
gp_vars <- c("mz_2s", "gaus_2s", "kern_2s", "bin_2s")
methods <- c("naive", "entropy", "pca")

SA_results <- SA
for (gp_var in gp_vars) {
  for (method in methods) {
    col_name <- paste0("HRI_", substr(gp_var, 1, nchar(gp_var)-3), "_", 
                      substr(method, 1, 1))
    
    SA_results[[col_name]] <- switch(method,
                                     "naive" = calc_naive(SA, gp_var),
                                     "entropy" = calc_entropy(SA, gp_var),
                                     "pca" = calc_pca(SA, gp_var)
    )
  }
}

# OUTPUT: Calculate correlation matrix
cat("\n=== OUTPUT 2: CORRELATION MATRIX OF HRI VARIANTS ===\n")
index_cols <- grep("^HRI_", names(SA_results), value = TRUE)
cor_matrix <- cor(st_drop_geometry(SA_results)[, index_cols], 
                  method = "spearman", 
                  use = "pairwise.complete.obs")

col_means <- colMeans(cor_matrix, na.rm = TRUE)

# Print mean values with names (to test stability)
print(col_means)

write.csv(round(cor_matrix, 3), "output/hri_correlation_matrix.csv")

################################################################################
# KEY COMPONENT 6: Spatial Cross-Validation Random Forest
################################################################################

cat("\n=== SPATIAL CROSS-VALIDATION RANDOM FOREST ===\n")

# Prepare model data
SA_results$Bach_p <- (SA_results$T10_4_OD_2 + SA_results$T10_4_HD_2 + 
                      SA_results$T10_4_PDT + SA_results$T10_4_DT) / 
                      SA_results$T10_4_TT
SA_results$NoAuto_p <- (SA_results$T11_1_FT + SA_results$T11_1_BIT + 
                        SA_results$T11_1_BUT + SA_results$T11_1_TDLT) / 
                        SA_results$T11_1_TT
SA_results$BVBHth_p <- (SA_results$T12_3_BT + SA_results$T12_3_VBT) / 
                        SA_results$T12_3_TT
SA_results$log_POP <- log(SA_results$T1_1AGETT)
SA_results$POPD <- SA_results$T1_1AGETT / SA_results$SHAPE_Area
SA_results$Bach_p[is.na(SA_results$Bach_p)] <- 0

# Create spatial neighbors (20 nearest neighbors)
coords <- st_centroid(st_geometry(SA_results), of_largest_polygon=TRUE)
sar.nb20 <- knearneigh(coords, k=20)
sar.nb20 <- knn2nb(sar.nb20)
sar.wt <- nb2listw(sar.nb20, style="W")

# Calculate spatial lags
SA_results$w_HRI_gaus_p <- lag.listw(sar.wt, as.numeric(SA_results$HRI_gaus_p))
SA_results$w_Bach_p <- lag.listw(sar.wt, as.numeric(SA_results$Bach_p))
SA_results$w_NoAuto_p <- lag.listw(sar.wt, SA_results$NoAuto_p)
SA_results$w_BVBHth_p <- lag.listw(sar.wt, SA_results$BVBHth_p)
SA_results$w_POPD <- lag.listw(sar.wt, SA_results$POPD)

# Create postcode dummy variables
SA <- as.data.frame(SA_results)
SA$D01 <- ifelse(SA$RoutingKey=="DUBLIN 1", 1, 0)
SA$D02 <- ifelse(SA$RoutingKey=="DUBLIN 2", 1, 0)
SA$D03 <- ifelse(SA$RoutingKey=="DUBLIN 3", 1, 0)
SA$D04 <- ifelse(SA$RoutingKey=="DUBLIN 4", 1, 0)
SA$D05 <- ifelse(SA$RoutingKey=="DUBLIN 5", 1, 0)
SA$D06 <- ifelse(SA$RoutingKey=="DUBLIN 6", 1, 0)
SA$D6W <- ifelse(SA$RoutingKey=="DUBLIN 6W", 1, 0)
SA$D07 <- ifelse(SA$RoutingKey=="DUBLIN 7", 1, 0)
SA$D08 <- ifelse(SA$RoutingKey=="DUBLIN 8", 1, 0)
SA$D09 <- ifelse(SA$RoutingKey=="DUBLIN 9", 1, 0)
SA$D10 <- ifelse(SA$RoutingKey=="DUBLIN 10", 1, 0)
SA$D11 <- ifelse(SA$RoutingKey=="DUBLIN 11", 1, 0)
SA$D12 <- ifelse(SA$RoutingKey=="DUBLIN 12", 1, 0)
SA$D13 <- ifelse(SA$RoutingKey=="DUBLIN 13", 1, 0)
SA$D17 <- ifelse(SA$RoutingKey=="DUBLIN 17", 1, 0)
SA$D20 <- ifelse(SA$RoutingKey=="DUBLIN 20", 1, 0)

## Spatial Cross-Validation Setup ##

model_data <- SA_results %>%
  filter(!is.na(HRI_gaus_p)) %>%
  select(HRI_gaus_p, In22_ED, NoAuto_p, BVBHth_p, POPD, log_dist)

model_data_spatial <- st_transform(model_data, 29902)

model_data <- model_data %>%
  mutate(across(-geometry, ~scale(.)[,1]))

# Create spatial blocks for cross-validation
set.seed(1111)
sb <- spatialBlock(
  speciesData = model_data,
  theRange = 1000,
  k = 5,
  selection = "random",
  iteration = 100
)

# Function to compute spatial lag
compute_spatial_lag <- function(data, neighbors) {
  y <- data$HRI_gaus_p
  lag_y <- lag.listw(nb2listw(neighbors, style = "W"), y)
  return(lag_y)
}

# Spatial CV with spatial lags
spatial_cv_results <- list()
predictors <- c("In22_ED", "NoAuto_p", "BVBHth_p", "POPD", "log_dist")

cv_models <- list()
train_data_list <- list()

cat("Running spatial cross-validation...\n")

for (fold in 1:5) {
  train_idx <- sb$folds[[fold]][[1]]
  test_idx <- sb$folds[[fold]][[2]]
  
  train_data <- model_data[train_idx, ]
  test_data <- model_data[test_idx, ]
  
  train_coords <- st_coordinates(st_centroid(st_geometry(train_data)))
  train_nb <- knn2nb(knearneigh(train_coords, k = 20))
  
  train_data_df <- st_drop_geometry(train_data)
  train_lag <- compute_spatial_lag(train_data_df, train_nb)
  train_data_df$spatial_lag <- train_lag
  
  rf_model <- ranger(
    HRI_gaus_p ~ .,
    data = train_data_df[, c("HRI_gaus_p", predictors, "spatial_lag")],
    importance = "permutation",
    num.trees = 500,
    mtry = 3,
    seed = 1111
  )
  
  test_data_df <- st_drop_geometry(test_data)[, c("HRI_gaus_p", predictors)]
  test_data_df$spatial_lag <- 0
  
  predictions <- predict(rf_model, test_data_df)$predictions
  
  spatial_cv_results[[fold]] <- data.frame(
    fold = fold,
    observed = test_data_df$HRI_gaus_p,
    predicted = predictions,
    residual = test_data_df$HRI_gaus_p - predictions
  )
  
  cv_models[[fold]] <- rf_model
  train_data_list[[fold]] <- train_data_df
}

cv_combined <- do.call(rbind, spatial_cv_results)

spatial_cv_rmse <- sqrt(mean(cv_combined$residual^2, na.rm = TRUE))
spatial_cv_r2 <- cor(cv_combined$observed, cv_combined$predicted, 
                     use = "complete.obs")^2

cat("\nSpatial CV Results:\n")
cat("  RMSE:", round(spatial_cv_rmse, 4), "\n")
cat("  R-squared:", round(spatial_cv_r2, 4), "\n")

## Moran's I on Residuals ##

model_data_spatial <- model_data %>%
  na.omit() %>%
  mutate(residuals = cv_combined$residual)

coords <- st_coordinates(st_centroid(st_geometry(model_data_spatial)))
nb <- knn2nb(knearneigh(coords, k = 20))
lw <- nb2listw(nb, style = "W")

moran_rf <- moran.test(model_data_spatial$residuals, lw)

cat("\nMoran's I on RF Residuals:\n")
print(moran_rf)

## Comparison to Spatial Econometric Models ##

cat("\n=== COMPARING TO SPATIAL ECONOMETRIC MODELS ===\n")

sar_formula <- as.formula("HRI_gaus_p ~ In22_ED + NoAuto_p + BVBHth_p + POPD + log_dist")
sar_model <- lagsarlm(sar_formula, data = model_data, listw = lw)
sem_model <- errorsarlm(sar_formula, data = model_data, listw = lw)
sac_model <- sacsarlm(sar_formula, data = model_data, listw = lw, type = "sac")
ols_model <- lm(sar_formula, data = model_data)

# Calculate metrics
calc_metrics <- function(observed, predicted) {
  residuals <- observed - predicted
  rmse <- sqrt(mean(residuals^2))
  mae <- mean(abs(residuals))
  r2 <- cor(observed, predicted)^2
  return(c(RMSE = rmse, MAE = mae, R2 = r2))
}

pred_sar <- predict(sar_model)
pred_sem <- predict(sem_model)
pred_sac <- sac_model$fitted.values
pred_ols <- predict(ols_model)

moran_ols <- moran.test(residuals(ols_model), lw)
moran_sar <- moran.test(residuals(sar_model), lw)
moran_sem <- moran.test(residuals(sem_model), lw)
moran_sac <- moran.test(residuals(sac_model), lw)

# OUTPUT: Model comparison table
cat("\n=== OUTPUT 3: MODEL COMPARISON ===\n")
comparison <- data.frame(
  Model = c("OLS", "SAR", "SEM", "SAC", "Random Forest Spatial CV"),
  RMSE = c(
    calc_metrics(model_data$HRI_gaus_p, pred_ols)[1],
    calc_metrics(model_data$HRI_gaus_p, pred_sar)[1],
    calc_metrics(model_data$HRI_gaus_p, pred_sem)[1],
    calc_metrics(model_data$HRI_gaus_p, pred_sac)[1],
    spatial_cv_rmse
  ),
  R2 = c(
    calc_metrics(model_data$HRI_gaus_p, pred_ols)[3],
    calc_metrics(model_data$HRI_gaus_p, pred_sar)[3],
    calc_metrics(model_data$HRI_gaus_p, pred_sem)[3],
    calc_metrics(model_data$HRI_gaus_p, pred_sac)[3],
    spatial_cv_r2
  ),
  Morans_I = c(
    moran_ols$estimate[1],
    moran_sar$estimate[1],
    moran_sem$estimate[1],
    moran_sac$estimate[1],
    moran_rf$estimate[1]
  ),
  Morans_pval = c(
    moran_ols$p.value,
    moran_sar$p.value,
    moran_sem$p.value,
    moran_sac$p.value,
    moran_rf$p.value
  )
)

print(comparison)
write.csv(comparison, "output/model_comparison.csv", row.names = FALSE)

# OUTPUT: SAC model summary
cat("\n=== OUTPUT 4: SPATIAL LAG + ERROR MODEL (SAC) SUMMARY ===\n")
sac_summary <- summary(sac_model)
print(sac_summary)

# Save SAC model summary to file
sink("output/sac_model_summary.txt")
print(sac_summary)
sink()

################################################################################
# KEY COMPONENT 7: Bootstrap Spatial CV for ALE and Importance with CIs
################################################################################

cat("\n=== BOOTSTRAP SPATIAL CV FOR UNCERTAINTY QUANTIFICATION ===\n")

# Function for spatial CV bootstrap
spatial_cv_bootstrap_rf <- function(data, folds, predictors, target, 
                                    n_iterations = 20, 
                                    compute_spatial_lag = TRUE, 
                                    k_neighbors = 20) {
  
  n <- nrow(data)
  cv_results <- list()
  result_counter <- 1
  
  has_spatial_lag <- "HRI_lag" %in% predictors
  
  for (iter in 1:n_iterations) {
    for (fold in 1:length(folds$folds)) {
      
      train_idx <- folds$folds[[fold]][[1]]
      test_idx <- folds$folds[[fold]][[2]]
      
      train_data_spatial <- data[train_idx, ]
      test_data_spatial <- data[test_idx, ]
      
      if (has_spatial_lag && compute_spatial_lag) {
        train_coords <- st_coordinates(st_centroid(st_geometry(train_data_spatial)))
        train_nb <- knn2nb(knearneigh(train_coords, k = k_neighbors))
        train_lw <- nb2listw(train_nb, style = "W", zero.policy = TRUE)
        
        train_data_df <- st_drop_geometry(train_data_spatial)
        train_lag <- lag.listw(train_lw, train_data_df[[target]], 
                               zero.policy = TRUE)
        train_data_df$HRI_lag <- train_lag
        
        test_coords <- st_coordinates(st_centroid(st_geometry(test_data_spatial)))
        test_lags <- numeric(nrow(test_data_spatial))
        for (i in 1:nrow(test_data_spatial)) {
          dists <- sqrt(rowSums((t(train_coords) - test_coords[i, ])^2))
          nearest_idx <- order(dists)[1:min(k_neighbors, length(dists))]
          test_lags[i] <- mean(train_data_df[[target]][nearest_idx], na.rm = TRUE)
        }
        
        test_data_df <- st_drop_geometry(test_data_spatial)
        test_data_df$HRI_lag <- test_lags
        
      } else {
        train_data_df <- st_drop_geometry(train_data_spatial)
        test_data_df <- st_drop_geometry(test_data_spatial)
      }
      
      rf_model <- ranger(
        formula = as.formula(paste(target, "~", paste(predictors, collapse = " + "))),
        data = train_data_df[, c(target, predictors)],
        num.trees = 500,
        importance = "permutation",
        seed = 1111 + iter * 10 + fold
      )
      
      if (has_spatial_lag && compute_spatial_lag) {
        full_coords <- st_coordinates(st_centroid(st_geometry(data)))
        full_nb <- knn2nb(knearneigh(full_coords, k = k_neighbors))
        full_lw <- nb2listw(full_nb, style = "W", zero.policy = TRUE)
        full_data_df <- st_drop_geometry(data)
        full_lag <- lag.listw(full_lw, full_data_df[[target]], 
                              zero.policy = TRUE)
        full_data_df$HRI_lag <- full_lag
      } else {
        full_data_df <- st_drop_geometry(data)
      }
      
      predictions_full <- predict(rf_model, full_data_df[, predictors])$predictions
      
      cv_results[[result_counter]] <- list(
        predictions = data.frame(
          cv_iter = result_counter,
          fold = fold,
          iteration = iter,
          row_id = 1:n,
          prediction = predictions_full,
          in_training = 1:n %in% train_idx
        ),
        importance = rf_model$variable.importance,
        fold_id = fold,
        iteration_id = iter,
        spatial_lag_computed = has_spatial_lag
      )
      
      result_counter <- result_counter + 1
    }
  }
  
  return(cv_results)
}

# Add spatial lag to model data
model_data_spatial$HRI_lag <- lag.listw(lw, model_data_spatial$HRI_gaus_p)

predictors <- c("In22_ED", "NoAuto_p", "BVBHth_p", "POPD", "log_dist", "HRI_lag")
target <- "HRI_gaus_p"

cat(paste0("Running spatial CV bootstrap: 5 folds × 20 iterations = 100 models\n"))
cat("This may take several minutes...\n")

cv_boot_results <- spatial_cv_bootstrap_rf(
  data = model_data_spatial,
  folds = sb,
  predictors = predictors,
  target = target,
  n_iterations = 20,
  compute_spatial_lag = TRUE,
  k_neighbors = 20
)

## Extract Variable Importance with CIs ##

importance_list <- lapply(cv_boot_results, function(x) x$importance)

importance_boot <- bind_rows(
  lapply(seq_along(importance_list), function(i) {
    df <- as.data.frame(t(importance_list[[i]]))
    df$cv_iter <- i
    df$fold <- cv_boot_results[[i]]$fold_id
    df$iteration <- cv_boot_results[[i]]$iteration_id
    return(df)
  })
) %>%
  pivot_longer(-c(cv_iter, fold, iteration), names_to = "variable", 
               values_to = "importance")

importance_summary <- importance_boot %>%
  group_by(variable) %>%
  summarise(
    mean = mean(importance),
    sd = sd(importance),
    lower = quantile(importance, 0.025),
    upper = quantile(importance, 0.975),
    cv = sd(importance) / mean(importance)
  ) %>%
  ungroup()

# OUTPUT: Variable importance with CIs
cat("\n=== OUTPUT 5: VARIABLE IMPORTANCE WITH 95% CONFIDENCE INTERVALS ===\n")
print(importance_summary %>% arrange(desc(mean)))
write.csv(importance_summary, "output/variable_importance_with_ci.csv", 
          row.names = FALSE)

# Plot importance with error bars
p_importance <- ggplot(importance_summary, 
                       aes(x = reorder(variable, mean), y = mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  coord_flip() +
  labs(
    title = "Variable Importance with 95% Confidence Intervals",
    subtitle = "Based on Spatial Cross-Validation Bootstrap",
    x = "Variable",
    y = "Importance (Permutation)"
  ) +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave("output/variable_importance_plot.png", p_importance, 
       width = 8, height = 6, dpi = 300)

## Calculate ALE Plots with Confidence Ribbons ##

# Function to calculate ALE for each CV fold
ale_spatial_cv <- function(data, cv_results, folds, predictor, predictors, 
                           target, compute_spatial_lag = TRUE, k_neighbors = 8) {
  
  ale_cv_list <- list()
  has_spatial_lag <- "HRI_lag" %in% predictors
  
  for (i in seq_along(cv_results)) {
    
    fold_id <- cv_results[[i]]$fold_id
    iter_id <- cv_results[[i]]$iteration_id
    
    train_idx <- folds$folds[[fold_id]][[1]]
    train_data_spatial <- data[train_idx, ]
    
    if (has_spatial_lag && compute_spatial_lag) {
      train_coords <- st_coordinates(st_centroid(st_geometry(train_data_spatial)))
      train_nb <- knn2nb(knearneigh(train_coords, k = k_neighbors))
      train_lw <- nb2listw(train_nb, style = "W", zero.policy = TRUE)
      
      train_data_df <- st_drop_geometry(train_data_spatial)
      train_lag <- lag.listw(train_lw, train_data_df[[target]], 
                             zero.policy = TRUE)
      train_data_df$HRI_lag <- train_lag
    } else {
      train_data_df <- st_drop_geometry(train_data_spatial)
    }
    
    rf_model <- ranger(
      formula = as.formula(paste(target, "~", paste(predictors, collapse = " + "))),
      data = train_data_df[, c(target, predictors)],
      num.trees = 500,
      seed = 1111 + iter_id * 10 + fold_id
    )
    
    yhat <- function(X.model, newdata) {
      predict(X.model, newdata)$predictions
    }
    
    ale_result <- tryCatch({
      ALEPlot(
        X = train_data_df[, predictors],
        X.model = rf_model,
        pred.fun = yhat,
        J = which(predictors == predictor),
        K = 40
      )
    }, error = function(e) NULL)
    
    if (!is.null(ale_result)) {
      ale_cv_list[[i]] <- data.frame(
        x = ale_result$x.values,
        ale = ale_result$f.values,
        cv_iter = i,
        fold = fold_id,
        iteration = iter_id
      )
    }
  }
  
  ale_combined <- bind_rows(ale_cv_list)
  
  ale_summary <- ale_combined %>%
    group_by(x) %>%
    summarise(
      ale_mean = mean(ale),
      ale_sd = sd(ale),
      ale_lower = quantile(ale, 0.025),
      ale_upper = quantile(ale, 0.975),
      n_folds = n()
    ) %>%
    filter(n_folds >= 10)
  
  return(list(
    summary = ale_summary,
    raw = ale_combined
  ))
}

# Fit main model for variable ranking
rf_main <- ranger(
  formula = as.formula(paste(target, "~", paste(predictors, collapse = " + "))),
  data = st_drop_geometry(model_data_spatial)[, c(target, predictors)],
  num.trees = 500,
  importance = "permutation"
)

top_vars <- names(sort(rf_main$variable.importance, decreasing = TRUE))[1:6]

cat("\nCalculating ALE plots with spatial CV uncertainty for top variables...\n")

ale_results <- list()
for (var in top_vars) {
  cat(paste0("  Processing: ", var, "\n"))
  ale_results[[var]] <- ale_spatial_cv(
    data = model_data_spatial,
    cv_results = cv_boot_results,
    folds = sb,
    predictor = var,
    predictors = predictors,
    target = target,
    compute_spatial_lag = TRUE,
    k_neighbors = 8
  )
}

# OUTPUT: ALE plots with confidence ribbons
cat("\n=== OUTPUT 6: ALE PLOTS WITH 95% CONFIDENCE RIBBONS ===\n")

ale_plots <- list()

for (var in names(ale_results)) {
  
  ale_data <- ale_results[[var]]$summary
  
  if (nrow(ale_data) < 5) {
    cat(paste0("Skipping ", var, " - insufficient data points\n"))
    next
  }
  
  # LOESS smoothing
  loess_mean <- loess(ale_mean ~ x, data = ale_data, span = 0.3)
  loess_lower <- loess(ale_lower ~ x, data = ale_data, span = 0.3)
  loess_upper <- loess(ale_upper ~ x, data = ale_data, span = 0.3)
  
  ale_data$ale_mean_smooth <- predict(loess_mean)
  ale_data$ale_lower_smooth <- predict(loess_lower)
  ale_data$ale_upper_smooth <- predict(loess_upper)
  
  p <- ggplot(ale_data, aes(x = x)) +
    geom_ribbon(aes(ymin = ale_lower_smooth, ymax = ale_upper_smooth), 
                alpha = 0.3, fill = "steelblue") +
    geom_line(aes(y = ale_mean_smooth), size = 1.2, color = "darkblue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      title = paste0("ALE: ", var),
      x = var,
      y = "Effect on HRI"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      axis.title = element_text(size = 9)
    )
  
  ale_plots[[var]] <- p
}

# Combine all ALE plots into one figure
library(gridExtra)
ale_combined_plot <- do.call(grid.arrange, c(ale_plots, ncol = 2))
ggsave("output/ale_plots_combined.png", ale_combined_plot, 
       width = 12, height = 12, dpi = 300)

################################################################################
# KEY COMPONENT 8: Interactive Maps
################################################################################

cat("\n=== CREATING INTERACTIVE MAPS ===\n")

SA_results <- st_transform(SA_results, 4326)

# Map 1: Health Rating Index
palhri <- colorQuantile("RdBu", SA_results$HRI_gaus_p, n = 5)
map_hri <- leaflet(SA_results) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolygons(
    fillColor = ~palhri(HRI_gaus_p),
    fillOpacity = 0.5,
    color = "white",
    weight = 0.1,
    label = ~paste0("HRI: ", round(HRI_gaus_p, 2))
  ) %>%
  addLegend(
    pal = palhri,
    values = ~HRI_gaus_p,
    title = "Health Rating Index",
    position = "bottomright"
  )

# Map 2: Green Space Access
palgs <- colorQuantile("YlGn", SA_results$gs_access, n = 5)
map_gs <- leaflet(SA_results) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolygons(
    fillColor = ~palgs(gs_access),
    fillOpacity = 0.5,
    color = "white",
    weight = 0.1,
    label = ~paste0("Park access: ", round(gs_access, 2))
  ) %>%
  addLegend(
    pal = palgs,
    values = ~gs_access,
    title = "Park Access",
    position = "bottomright"
  )

# Map 3: Blue Space Access
palbs <- colorQuantile("Blues", SA_results$bs_access, n = 5)
map_bs <- leaflet(SA_results) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolygons(
    fillColor = ~palbs(bs_access),
    fillOpacity = 0.5,
    color = "white",
    weight = 0.1,
    label = ~paste0("Blue space access: ", round(bs_access, 2))
  ) %>%
  addLegend(
    pal = palbs,
    values = ~bs_access,
    title = "Blue Space Access",
    position = "bottomright"
  )

# Map 4: GP Access
palgp <- colorQuantile("magma", SA_results$gaus_2s, n = 5)
map_gp <- leaflet(SA_results) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolygons(
    fillColor = ~palgp(gaus_2s),
    fillOpacity = 0.5,
    color = "white",
    weight = 0.1,
    label = ~paste0("GP access: ", round(gaus_2s, 4))
  ) %>%
  addLegend(
    pal = palgp,
    values = ~gaus_2s,
    title = "GP Access",
    position = "bottomright"
  )

# Map 5: Air Pollution YLL
palaq <- colorQuantile("YlOrBr", SA_results$YLLAQPC, n = 5)
map_aq <- leaflet(SA_results) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolygons(
    fillColor = ~palaq(YLLAQPC),
    fillOpacity = 0.5,
    color = "white",
    weight = 0.1,
    label = ~paste0("YLL (per 100K) air pollution: ", round(YLLAQPC, 2))
  ) %>%
  addLegend(
    pal = palaq,
    values = ~YLLAQPC,
    title = "YLL Air Pollution",
    position = "bottomright"
  )

# Map 6: Noise YLL
SA_results_filtered <- SA_results %>%
  filter(!is.na(YLLNPC) & YLLNPC > 0)

paln <- colorQuantile("Reds", SA_results_filtered$YLLNPC, n = 5)
map_noise <- leaflet(SA_results_filtered) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolygons(
    fillColor = ~paln(YLLNPC),
    fillOpacity = 0.5,
    color = "white",
    weight = 0.1,
    label = ~paste0("YLL (per 100K) noise: ", round(YLLNPC, 2))
  ) %>%
  addLegend(
    pal = paln,
    values = ~YLLNPC,
    title = "YLL Road Noise",
    position = "bottomright"
  )

# Save maps as HTML widgets
library(htmlwidgets)

saveWidget(map_hri, "output/map_hri.html")
saveWidget(map_gs, "output/map_green_space.html")
saveWidget(map_bs, "output/map_blue_space.html")
saveWidget(map_gp, "output/map_gp_access.html")
saveWidget(map_aq, "output/map_air_pollution_yll.html")
saveWidget(map_noise, "output/map_noise_yll.html")

cat("\n=== OUTPUT 7: MAPS SAVED ===\n")
cat("Maps have been saved as HTML files in the output/ directory\n")

################################################################################
# FINAL SUMMARY
################################################################################

cat("\n" , rep("=", 80), "\n", sep = "")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("Key outputs saved:\n")
cat("  1. Total years of life lost: output/total_years_life_lost.csv\n")
cat("  2. HRI correlation matrix: output/hri_correlation_matrix.csv\n")
cat("  3. Model comparison: output/model_comparison.csv\n")
cat("  4. SAC model summary: output/sac_model_summary.txt\n")
cat("  5. Variable importance with CIs: output/variable_importance_with_ci.csv\n")
cat("  6. ALE plots: output/ale_plots_combined.png\n")
cat("  7. Interactive maps: output/map_*.html (6 files)\n\n")

cat("For questions or issues, please refer to:\n")
cat("Credit, K., Damanpreet, K., and Eccles, E. 2025.\n")
cat("DOI: https://doi.org/10.5281/zenodo.15183740\n\n")

cat(rep("=", 80), "\n", sep = "")
