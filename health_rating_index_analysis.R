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

# Install SArf locally
# Set working directory to your package
setwd("~/Dropbox/Packages/SArf/SArf/")
devtools::document()
devtools::install()
# Restart R session (important!)
# In RStudio: Session > Restart R


# Install SArf from GitHub
devtools::install_github("kcredit/SArf", auth_token = NULL, force = TRUE)

# Load ONLY what you need for SArf
library(ggplot2)  # Load first
library(sf)
library(dplyr)
library(tidyr)
library(purrr)
library(spdep)
library(leaflet)
library(leaflet.extras)
library(blockCV)  # Load after ggplot2
library(SArf)     # Load last


# Resolve conflicts
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")


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
  select(SA_PUB2022, POP30P, SHAPE_Area, T1_1AGETT, T6_2_TH, total_pop_all_f, total_pop_all_m, LE_all_m, 
         LE_all_f, LE_30P_m, LE_30P_f,YLLPM25, YLLNO2, YLLO3, YLLAQ, 
         T10_4_OD_2, T10_4_HD_2, T10_4_PDT,T10_4_DT, T10_4_TT, T11_1_FT, 
         T11_1_BIT,T11_1_BUT, T11_1_TDLT,T11_1_TT, T12_3_BT, T12_3_VBT, 
         T12_3_TT,T1_1AGE6_1,T1_1AGE6_2,T1_1AGE7_1, T1_1AGE7_2,
         T1_1AGE8_1,T1_1AGEG_1,T2_1IEC,T2_1TC)

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
  select(SA_PUB2022, RoutingKey, In22_ED, SOUTH, log_dist, total)

SA <- SA %>%
  left_join(sa_access, by = c("SA_PUB2022" = "SA_PUB202")) %>%
  left_join(air, by = "SA_PUB2022")

summary(SA[, c("YLLAQ", "YLLN")])

# Create rate variables
SA$YLLAQPC <- (SA$YLLAQ / SA$T1_1AGETT) * 100000
SA$YLLNPC <- (SA$YLLN / SA$T1_1AGETT) * 100000
SA$ph_rate <- SA$total/SA$T6_2_TH
SA$ph_rate[is.infinite(SA$ph_rate) | is.nan(SA$ph_rate)] <- 0
SA$ph_rate <- pmin(SA$ph_rate, 1)

# Function to calculate naive weights
calc_naive <- function(df, gp_var) {
  df_data <- st_drop_geometry(df)
  
  # Normalize to 0-1 first
  df_norm <- df_data %>%
    select(all_of(c(gp_var, "gs_access", "bs_access", "YLLAQPC", "YLLNPC", "ph_rate"))) %>%
    mutate(across(everything(), ~ (. - min(.)) / (max(.) - min(.))))
  
  # Flip the bad indicators
  df_norm <- df_norm %>%
    mutate(
      YLLAQPC = 1 - YLLAQPC,
      YLLNPC = 1 - YLLNPC,
      ph_rate = 1 - ph_rate
    )
  
  # Simple average
  hri <- (1/6) * (df_norm[[gp_var]] + df_norm$gs_access + df_norm$bs_access + 
                    df_norm$YLLAQPC + df_norm$YLLNPC + df_norm$ph_rate)
  
  # Scale to 0-1
  (hri - min(hri)) / (max(hri) - min(hri))
}

# Function to calculate entropy weights
calc_entropy <- function(df, gp_var) {
  df_data <- st_drop_geometry(df)
  df_norm <- df_data %>%
    select(all_of(c(gp_var, "gs_access", "bs_access", "YLLAQPC", "YLLNPC", "ph_rate"))) %>%  # Added ph_rate
    mutate(across(everything(), ~ (. - min(.)) / (max(.) - min(.)))) %>%
    mutate(YLLAQPC = 1 - YLLAQPC, 
           YLLNPC = 1 - YLLNPC,
           ph_rate = 1 - ph_rate)  # Reverse ph_rate (higher = worse)
  
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
    select(all_of(c(gp_var, "gs_access", "bs_access", "YLLAQPC", "YLLNPC", "ph_rate"))) %>%  # Added ph_rate
    mutate(YLLAQPC = max(YLLAQPC, na.rm = TRUE) - YLLAQPC,
           YLLNPC = max(YLLNPC, na.rm = TRUE) - YLLNPC,
           ph_rate = max(ph_rate, na.rm = TRUE) - ph_rate) %>%  # Reverse ph_rate
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
# KEY COMPONENTS 6-7: Spatial Cross-Validation Random Forest, and 
# Bootstrap Spatial CV for ALE and Importance with CIs, using the 
# SArf package (https://github.com/kcredit/SArf)
################################################################################

cat("\n=== SPATIAL CROSS-VALIDATION RANDOM FOREST ===\n")

# Prepare model data
SA_results$NoAuto_p <- (SA_results$T11_1_FT + SA_results$T11_1_BIT + 
                        SA_results$T11_1_BUT + SA_results$T11_1_TDLT) / 
                        SA_results$T11_1_TT
SA_results$BVBHth_p <- (SA_results$T12_3_BT + SA_results$T12_3_VBT) / 
                        SA_results$T12_3_TT
SA_results$ov60 <- (SA_results$T1_1AGE6_1+SA_results$T1_1AGE6_2+SA_results$T1_1AGE7_1+
                      SA_results$T1_1AGE7_2+SA_results$T1_1AGE8_1+SA_results$T1_1AGEG_1) / SA_results$T1_1AGETT
SA_results$nonIrish <- 1-(SA_results$T2_1IEC / SA_results$T2_1TC)
SA_results$POPD <- SA_results$T1_1AGETT / SA_results$SHAPE_Area

# Create spatial neighbors (20 nearest neighbors)
coords <- st_centroid(st_geometry(SA_results), of_largest_polygon=TRUE)
sar.nb20 <- spdep::knearneigh(coords, k=20)
sar.nb20 <- spdep::knn2nb(sar.nb20)
sar.wt <- spdep::nb2listw(sar.nb20, style="W")

## Spatial Cross-Validation Setup ##

model_data <- SA_results %>%
  filter(!is.na(HRI_gaus_e)) %>%
  select(HRI_gaus_e, In22_ED, NoAuto_p, POPD, log_dist, ov60, nonIrish)

model_data <- model_data %>%
  mutate(across(-geometry, ~scale(.)[,1]))

# Export
# write_sf(model_data, "model_data.shp")

#Smaller sample for testing
# sampled_data <- model_data %>% slice_sample(n = 500)

# Run SArf
results <- SArf(
  HRI_gaus_e ~ In22_ED + NoAuto_p + POPD + log_dist + ov60 + nonIrish,
  data = model_data,
  k_neighbors = 20,            # Number of neighbors for spatial weights
  n_folds = 5,                 # Spatial CV folds
  n_bootstrap = 20,            # Bootstrap iterations
  num_trees = 500,             # Trees in random forest
  include_naive_rf = TRUE,     # For comparison with "naive" random forest
  naive_test_fraction = 0.2,   
  create_map = TRUE,           # Spatial folds map
  block_range = 5000,          # Size of the spatial blocks in units of metres
  verbose = TRUE               # Print progress
)

# print(results) # Shows all results

# View specific results outputs
results$moran_test              # Moran's I test for outcome variable
results$moran_plot              # Moran's I scatter plot
results$model_comparison        # RF vs OLS/SAR/SEM/SAC table
show_models(results)            # Full spatial econometric model results
results$variable_importance     # Importance with CIs
results$importance_plot         # Importance bar chart
results$ale_results             # ALE data
results$ale_plots               # ALE plots with CIs

# =================================================================
# VERIFICATION CODE TO ENSURE SPATIAL LAGS ARE CALCULATED CORRECTLY
# =================================================================

# Extract spatial_weights from results
spatial_weights <- results$spatial_weights

# ========================================
# 1. Get fold structure for verification
# ========================================
fold_1 <- results$spatial_cv_results$predictions %>%
  filter(fold == 1, iteration == 1)

test_idx <- fold_1 %>% filter(in_training == FALSE) %>% pull(row_id)
train_idx <- fold_1 %>% filter(in_training == TRUE) %>% pull(row_id)

cat("Fold 1, Iteration 1:\n")
cat("Training observations:", length(train_idx), "\n")
cat("Test observations:", length(test_idx), "\n\n")

# ========================================
# 2. Get coordinates and extract k
# ========================================
all_coords <- st_coordinates(st_centroid(st_geometry(model_data)))
train_coords <- all_coords[train_idx, , drop = FALSE]

# Extract k from spatial_weights
k_neighbors <- round(mean(sapply(spatial_weights$neighbours, length)))
cat("k (number of neighbors):", k_neighbors, "\n\n")

# ========================================
# 3. For EACH test observation, verify neighbors are ALL from training set
# ========================================
cat("=== Verifying Spatial Lag Calculation ===\n\n")

# Get response variable name from results
response_var <- all.vars(results$formula)[1]

verification_results <- data.frame(
  test_obs = test_idx,
  stored_lag = numeric(length(test_idx)),
  calculated_lag = numeric(length(test_idx)),
  match = logical(length(test_idx)),
  n_neighbors = numeric(length(test_idx)),
  any_from_test = logical(length(test_idx))
)

for (i in 1:length(test_idx)) {
  obs <- test_idx[i]
  
  # Get stored spatial lag from predictions
  stored_lag <- fold_1 %>% filter(row_id == obs) %>% pull(spatial_lag)
  
  # Calculate what spatial lag SHOULD be (using training neighbors only)
  obs_coord <- all_coords[obs, , drop = FALSE]
  
  # Find distances to ALL training observations
  dists <- sqrt(rowSums(sweep(train_coords, 2, obs_coord, "-")^2))
  
  # Find k nearest TRAINING neighbors
  k_to_use <- min(k_neighbors, length(dists))
  nearest_indices <- order(dists)[1:k_to_use]
  neighbor_rows <- train_idx[nearest_indices]
  
  # Calculate spatial lag with W-style weights (simple mean)
  calculated_lag <- mean(model_data[[response_var]][neighbor_rows])
  
  # Check if any neighbors are from test set (should be NONE!)
  any_from_test <- any(neighbor_rows %in% test_idx)
  
  # Store results
  verification_results$stored_lag[i] <- stored_lag
  verification_results$calculated_lag[i] <- calculated_lag
  verification_results$match[i] <- abs(stored_lag - calculated_lag) < 0.001
  verification_results$n_neighbors[i] <- k_to_use
  verification_results$any_from_test[i] <- any_from_test
}

# ========================================
# 4. Summary Statistics
# ========================================
cat("=== VERIFICATION RESULTS ===\n\n")

cat("Total test observations checked:", nrow(verification_results), "\n")
cat("Stored values match calculated:", sum(verification_results$match), "/", nrow(verification_results), "\n")
cat("Observations with neighbors from test set:", sum(verification_results$any_from_test), "\n\n")

if (sum(verification_results$any_from_test) > 0) {
  cat("❌ PROBLEM: Some test observations have neighbors from test set!\n")
  cat("This indicates data leakage.\n\n")
  
  # Show which observations have test neighbors
  problem_obs <- verification_results %>% filter(any_from_test == TRUE)
  cat("Problem observations:\n")
  print(problem_obs %>% select(test_obs, any_from_test))
  
} else {
  cat("✅ SUCCESS: NO test observations have neighbors from test set!\n")
  cat("All spatial lags use only training neighbors.\n\n")
}

if (all(verification_results$match)) {
  cat("✅ SUCCESS: All stored spatial_lag values match manual calculations!\n")
  cat("The order() bug is fixed.\n\n")
} else {
  cat("❌ PROBLEM: Some stored values don't match calculations!\n")
  cat("Number of mismatches:", sum(!verification_results$match), "\n\n")
  
  # Show mismatches
  mismatches <- verification_results %>% filter(match == FALSE)
  cat("Mismatched observations (first 10):\n")
  print(head(mismatches %>% select(test_obs, stored_lag, calculated_lag), 10))
}

# ========================================
# 5. Detailed Check for Random Sample
# ========================================
cat("\n=== DETAILED CHECK: 5 Random Test Observations ===\n\n")

sample_obs <- sample(test_idx, min(5, length(test_idx)))

for (obs in sample_obs) {
  cat("Observation:", obs, "\n")
  
  # Get stored value
  stored <- fold_1 %>% filter(row_id == obs) %>% pull(spatial_lag)
  
  # Calculate manually
  obs_coord <- all_coords[obs, , drop = FALSE]
  dists <- sqrt(rowSums(sweep(train_coords, 2, obs_coord, "-")^2))
  k_to_use <- min(k_neighbors, length(dists))
  nearest_indices <- order(dists)[1:k_to_use]
  neighbor_rows <- train_idx[nearest_indices]
  calculated <- mean(model_data[[response_var]][neighbor_rows])
  
  cat("  Stored spatial_lag:    ", round(stored, 6), "\n")
  cat("  Calculated spatial_lag:", round(calculated, 6), "\n")
  cat("  Match:                 ", abs(stored - calculated) < 0.001, "\n")
  cat("  Number of neighbors:   ", k_to_use, "\n")
  cat("  All from training set: ", all(neighbor_rows %in% train_idx), "\n")
  cat("  Any from test set:     ", any(neighbor_rows %in% test_idx), "\n")
  
  # Show neighbor IDs
  cat("  Neighbor row IDs:      ", paste(head(neighbor_rows, 10), collapse = ", "))
  if (k_to_use > 10) cat(", ...")
  cat("\n\n")
}

# ========================================
# 6. Final Summary
# ========================================
cat("=== FINAL ASSESSMENT ===\n\n")

if (sum(verification_results$any_from_test) == 0 && all(verification_results$match)) {
  cat("✅✅✅ PERFECT! ✅✅✅\n")
  cat("1. NO data leakage (all neighbors from training set)\n")
  cat("2. Storage is correct (stored values match calculations)\n")
  cat("3. Your spatial CV implementation is working correctly!\n")
} else {
  if (sum(verification_results$any_from_test) > 0) {
    cat("❌ Data leakage detected\n")
  }
  if (!all(verification_results$match)) {
    cat("❌ Storage bug detected\n")
  }
  cat("\nSee details above for specific issues.\n")
}

################################################################################
# KEY COMPONENT 8: Interactive Maps
################################################################################

cat("\n=== CREATING INTERACTIVE MAPS ===\n")

SA_results <- st_transform(SA_results, 4326)

# Map 1: Health Rating Index
palhri <- colorQuantile("RdYlBu", SA_results$HRI_gaus_e, n = 5)
map_hri <- leaflet(SA_results) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolygons(
    fillColor = ~palhri(HRI_gaus_e),
    fillOpacity = 0.5,
    color = "white",
    weight = 0.1,
    label = ~paste0("HRI: ", round(HRI_gaus_e, 2))
  ) %>%
  addLegend(
    pal = palhri,
    values = ~HRI_gaus_e,
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

# Map 7: Poor quality housing
SA_results$ph_rate[is.na(SA_results$ph_rate) | is.infinite(SA_results$ph_rate)] <- 0

palhq <- colorQuantile("viridis", SA_results$ph_rate, n = 5)
map_hq <- leaflet(SA_results) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolygons(
    fillColor = ~palhq(ph_rate),
    fillOpacity = 0.5,
    color = "white",
    weight = 0.1,
    label = ~paste0("Rate of housing with BER of E or worse: ", round(ph_rate, 2))
  ) %>%
  addLegend(
    pal = palhq,
    values = ~ph_rate,
    title = "Poor quality housing rate",
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
saveWidget(map_hq, "output/map_housing_quality.html")

cat("\n=== OUTPUT 7: MAPS SAVED ===\n")
cat("Maps have been saved as HTML files in the output/ directory\n")
