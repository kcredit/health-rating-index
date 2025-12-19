################################################################################
# Health Rating Index (HRI) Analysis for Dublin, Ireland
# 
# Description:
#   Analyzes environmental health burdens and benefits to create a Health 
#   Rating Index for small areas in Dublin, Ireland. Components include:
#   - GP accessibility (E2SFCA with multiple decay functions)
#   - Green and blue space accessibility
#   - Years of Life Lost from air pollution (PM2.5, NO2, O3) and road noise
#   - Poor quality housing indicators
#   
#   Methods:
#   - Multiple index weighting approaches (naive, entropy, PCA)
#   - Spatial cross-validation random forest with bootstrap confidence intervals
#   - Comparison with spatial econometric models (OLS, SAR, SEM, SAC)
#   - Interactive web maps of results
################################################################################

# Clear workspace
rm(list=ls())

################################################################################
# 1. LOAD LIBRARIES
################################################################################

# Core packages
library(ggplot2)
library(sf)
library(dplyr)
library(tidyr)
library(purrr)

# Spatial analysis
library(spdep)
library(spatialreg)
library(blockCV)
library(SArf)  # Install: devtools::install_github("kcredit/SArf")

# Visualization
library(leaflet)
library(leaflet.extras)
library(htmlwidgets)

# Resolve namespace conflicts
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

################################################################################
# 2. SETUP
################################################################################

# Set working directory (update to your local path)
setwd("~/Dropbox/Grants/New Foundations/Smart D8 Dashboard/Anaysis Data/health-rating-index")

# Create output directory
if (!dir.exists("output")) dir.create("output")

################################################################################
# 3. LOAD DATA
################################################################################

cat("Loading data...\n")

# Origin and destination data
sa <- read.csv("data/SA_Origin.csv")
gp <- read.csv("data/GP_Dest.csv")
pk <- read.csv("data/Park_Dest.csv")
bs <- read.csv("data/Blue_Dest.csv")

# Pre-calculated travel time matrices
dist <- read.csv("data/dist_sa_gp.csv")
distpk <- read.csv("data/dist_sa_pk1.csv")
disbs <- read.csv("data/dist_sa_bs.csv")

# Environmental health data
noise <- read.csv("data/NoiseContours_2011_Dissolved_SA2022_Intersection_Area.csv")
air <- read.csv("data/SA2022_Dublin_Postcode_greenR_AirView_Noise_GS_distroad_GP_NO_TUNNEL_airviewmedians.csv")

# Irish life table for YLL calculations
life_table <- read.csv("data/irish_life_table.csv")
life_table <- life_table[, c("Age", "Sex", "LifeExpectancy")]

################################################################################
# 4. GP ACCESSIBILITY (Enhanced 2SFCA)
################################################################################

cat("\n=== CALCULATING GP ACCESSIBILITY ===\n")

# Adjust GP supply by county (per-capita ratios)
gp <- gp %>%
  mutate(SUPPLY = case_when(
    COUNTY_ENG == "SOUTH DUBLIN" ~ 2.106383676,
    COUNTY_ENG == "FINGAL" ~ 2.259737083,
    COUNTY_ENG == "DUBLIN CITY" ~ 1.330671623,
    COUNTY_ENG == "DUN LAOGHAIRE/RATHDOWN" ~ 1.730009988,
    TRUE ~ SUPPLY
  ))

# Distance decay functions
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
    # Gaussian decay (Bauer & Groneberg 2016)
    d0 <- 20
    weight <- exp(-(travel_time^2) / (2*d0^2))
    weight <- ifelse(travel_time > 30, 0, weight)
  } else if (method == "kernel") {
    # Epanechnikov kernel (Dai & Wang 2011)
    d0 <- 30
    weight <- ifelse(travel_time < d0, (3/4) * (1 - (travel_time/d0)^2), 0)
  } else {
    # Binary (traditional 2SFCA - Luo & Wang 2003)
    weight <- ifelse(travel_time <= 30, 1, 0)
  }
  return(weight)
}

# E2SFCA Step 1: Provider-to-population ratios
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
    mutate(ratio = ifelse(total_weighted_demand > 0, 
                         capacity / total_weighted_demand, 0))
  
  return(provider_ratios)
}

# E2SFCA Step 2: Accessibility scores
step2_accessibility <- function(travel_time_matrix, provider_ratios, 
                                max_travel_time = 30, decay_method = "multi_zone") {
  ttm_filtered <- travel_time_matrix %>%
    filter(travel_time <= max_travel_time)
  
  ttm_weighted <- ttm_filtered %>%
    mutate(weight = calculate_decay_weight(travel_time, method = decay_method))
  
  ttm_with_ratios <- ttm_weighted %>%
    left_join(provider_ratios %>% select(destination_id, ratio), by = "destination_id")
  
  accessibility_scores <- ttm_with_ratios %>%
    group_by(origin_id) %>%
    summarise(
      accessibility = sum(ratio * weight, na.rm = TRUE),
      n_providers_in_catchment = n()
    ) %>%
    ungroup()
  
  return(accessibility_scores)
}

# Complete E2SFCA workflow
e2sfca_analysis <- function(travel_time_matrix, demand_data, supply_data,
                            max_travel_time = 30, decay_method = "multi_zone") {
  cat("  Calculating provider ratios...\n")
  provider_ratios <- step1_provider_ratios(
    travel_time_matrix, demand_data, supply_data, max_travel_time, decay_method
  )
  
  cat("  Calculating accessibility scores...\n")
  accessibility_scores <- step2_accessibility(
    travel_time_matrix, provider_ratios, max_travel_time, decay_method
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
  rename(origin_id = id, population = T1_1AGETT) %>%
  mutate(origin_id = as.character(origin_id))

supply_data <- gp %>%
  select(id, SUPPLY) %>%
  rename(destination_id = id, capacity = SUPPLY) %>%
  mutate(destination_id = as.character(destination_id))

dist <- dist %>%
  rename(origin_id = from_id, destination_id = to_id, travel_time = travel_time_p50) %>%
  mutate(origin_id = as.character(origin_id), 
         destination_id = as.character(destination_id))

# Calculate accessibility using multiple decay methods
decay_methods <- c("multi_zone", "gaussian", "kernel", "binary")
gp_access_results <- list()

for (method in decay_methods) {
  cat(paste0("\nCalculating GP accessibility: ", method, " decay\n"))
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

# Merge all GP accessibility metrics
for (method in decay_methods) {
  sa <- merge(sa, gp_access_results[[method]], by = "id", all.x = TRUE)
}

cat("GP accessibility metrics calculated.\n")

################################################################################
# 5. BLUE SPACE ACCESSIBILITY
################################################################################

cat("\n=== CALCULATING BLUE SPACE ACCESSIBILITY ===\n")

# Find nearest coastline point
nearestbs <- disbs %>%
  group_by(from_id) %>%
  slice_min(travel_time_p50, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(id = from_id)

# Apply Gaussian decay (sigma = 20 minutes)
sigma <- 20
nearestbs$bs_access <- exp(-(nearestbs$travel_time_p50^2) / (2*sigma^2))

nearestbs <- nearestbs %>% select(id, bs_access)
sa <- merge(sa, nearestbs, by = "id", all.x = TRUE)

cat("Blue space accessibility calculated.\n")

################################################################################
# 6. GREEN SPACE ACCESSIBILITY
################################################################################

cat("\n=== CALCULATING GREEN SPACE ACCESSIBILITY ===\n")

# Load park polygons and origin points
o_sa <- read.csv("data/SA_Origin.csv")
parks_poly <- st_read("data/DCC Parks.shp", quiet = TRUE)

df_sf <- st_as_sf(o_sa, coords = c("lon", "lat"), crs = 4326)
parks_poly <- st_transform(parks_poly, 4326)
parks_poly <- st_make_valid(parks_poly)

# Maximum accessibility for points within parks
within_park <- lengths(st_intersects(df_sf, parks_poly)) > 0
df_sf$gs_access <- NA
df_sf$gs_access[within_park] <- 5

# Calculate accessibility for points outside parks
outside_ids <- df_sf$id[!within_park]

if (length(outside_ids) > 0) {
  distpk_outside <- distpk %>% filter(from_id %in% outside_ids)
  
  pk <- pk %>%
    select(id, Name) %>%
    mutate(id = as.character(id))
  
  distpk_outside <- merge(distpk_outside, pk, by.x = "to_id", by.y = "id")
  distpk_outside$Name[distpk_outside$Name == ""] <- "NA"
  
  # Minimum time to each unique park
  min_times <- distpk_outside %>%
    group_by(from_id, Name) %>%
    summarise(min_time = min(travel_time_p50), .groups = 'drop')
  
  # Find 5 nearest parks per origin
  nearest_5 <- min_times %>%
    group_by(from_id) %>%
    arrange(min_time) %>%
    slice_head(n = 5) %>%
    ungroup()
  
  # Apply Gaussian decay (sigma = 5 minutes)
  sigma <- 5
  nearest_5 <- nearest_5 %>%
    mutate(decayed = exp(-(min_time^2) / (2 * sigma^2)))
  
  # Sum decayed values
  gs_access_calculated <- nearest_5 %>%
    group_by(from_id) %>%
    summarise(gs_access = sum(decayed), .groups = 'drop')
  
  # Merge calculated accessibility
  df_sf <- df_sf %>%
    left_join(gs_access_calculated, by = c("id" = "from_id"))
  
  # Combine: 5 for parks, calculated for others
  df_sf <- df_sf %>%
    mutate(gs_access = coalesce(gs_access.x, gs_access.y)) %>%
    select(-gs_access.x, -gs_access.y)
}

o_sa <- st_drop_geometry(df_sf)
sa <- merge(sa, o_sa %>% select(id, gs_access), by = "id", all.x = TRUE)

cat("Green space accessibility calculated.\n")

# Save accessibility results
sa_access <- sa %>%
  select(id, SA_PUB202, mz_2s, gaus_2s, kern_2s, bin_2s, bs_access, gs_access)

################################################################################
# 7. YEARS OF LIFE LOST (YLL) CALCULATIONS
################################################################################

cat("\n=== CALCULATING YEARS OF LIFE LOST ===\n")

## AIR POLLUTION YLL ##
cat("Air pollution YLL...\n")

# Relative risks (WHO 2021 Air Quality Guidelines)
air$RRPM25 <- ifelse(!is.na(air$PM25_media) & air$PM25_media >= 5,
                     (((air$PM25_media - 5) / 10) * 0.08) + 1, 1)

air$RRNO2 <- ifelse(!is.na(air$NO2_median) & air$NO2_median >= 10,
                    (((air$NO2_median - 10) / 10) * 0.02) + 1, 1)

air$RRO3 <- ifelse(!is.na(air$O3_median) & air$O3_median >= 71.05,
                   (((air$O3_median - 71.05) / 10) * 0.0043) + 1, 1)

# Population attributable fractions
air$PAFPM25 <- (air$RRPM25 - 1) / air$RRPM25
air$PAFNO2 <- (air$RRNO2 - 1) / air$RRNO2
air$PAFO3 <- (air$RRO3 - 1) / air$RRO3

# Life expectancy interpolation
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

# Age-sex specific mortality rates (Ireland 2017, CSO DHA13)
mortality_rates <- data.frame(
  AgeGroup = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", 
               "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", 
               "75-79", "80-84", "85+"),
  Female = c(0.0006601, 0.0000886, 0.0000521, 0.0000519, 0.0000925, 0.0002121,
             0.0002898, 0.00054, 0.0007498, 0.0011977, 0.0019599, 0.0032209,
             0.0064141, 0.0080121, 0.015255, 0.0281099, 0.0499962, 0.13867),
  Male = c(0.0005456, 0, 0.0000996, 0.0001769, 0.0004129, 0.0004458,
           0.0004142, 0.0009302, 0.0011782, 0.0022411, 0.0037875, 0.0054263,
           0.008684, 0.0146503, 0.0226823, 0.0425332, 0.0705163, 0.1520202)
)

# Map census columns to mortality bands
age_band_mapping <- list(
  "0-4" = list(
    female_cols = c("T1_1AGE0F", "T1_1AGE1F", "T1_1AGE2F", "T1_1AGE3F", "T1_1AGE4F"),
    male_cols = c("T1_1AGE0M", "T1_1AGE1M", "T1_1AGE2M", "T1_1AGE3M", "T1_1AGE4M"),
    age_midpoint = 2.5
  ),
  "5-9" = list(
    female_cols = c("T1_1AGE5F", "T1_1AGE6F", "T1_1AGE7F", "T1_1AGE8F", "T1_1AGE9F"),
    male_cols = c("T1_1AGE5M", "T1_1AGE6M", "T1_1AGE7M", "T1_1AGE8M", "T1_1AGE9M"),
    age_midpoint = 7.5
  ),
  "10-14" = list(
    female_cols = c("T1_1AGE10F", "T1_1AGE11F", "T1_1AGE12F", "T1_1AGE13F", "T1_1AGE14F"),
    male_cols = c("T1_1AGE10M", "T1_1AGE11M", "T1_1AGE12M", "T1_1AGE13M", "T1_1AGE14M"),
    age_midpoint = 12.5
  ),
  "15-19" = list(
    female_cols = c("T1_1AGE15F", "T1_1AGE16F", "T1_1AGE17F", "T1_1AGE18F", "T1_1AGE19F"),
    male_cols = c("T1_1AGE15M", "T1_1AGE16M", "T1_1AGE17M", "T1_1AGE18M", "T1_1AGE19M"),
    age_midpoint = 17.5
  ),
  "20-24" = list(
    female_cols = "T1_1AGE201",
    male_cols = "T1_1AGE20_",
    age_midpoint = 22.5
  ),
  "25-29" = list(
    female_cols = "T1_1AGE251",
    male_cols = "T1_1AGE25_",
    age_midpoint = 27.5
  ),
  "30-34" = list(
    female_cols = "T1_1AGE301",
    male_cols = "T1_1AGE30_",
    age_midpoint = 32.5
  ),
  "35-39" = list(
    female_cols = "T1_1AGE351",
    male_cols = "T1_1AGE35_",
    age_midpoint = 37.5
  ),
  "40-44" = list(
    female_cols = "T1_1AGE401",
    male_cols = "T1_1AGE40_",
    age_midpoint = 42.5
  ),
  "45-49" = list(
    female_cols = "T1_1AGE451",
    male_cols = "T1_1AGE45_",
    age_midpoint = 47.5
  ),
  "50-54" = list(
    female_cols = "T1_1AGE501",
    male_cols = "T1_1AGE50_",
    age_midpoint = 52.5
  ),
  "55-59" = list(
    female_cols = "T1_1AGE551",
    male_cols = "T1_1AGE55_",
    age_midpoint = 57.5
  ),
  "60-64" = list(
    female_cols = "T1_1AGE601",
    male_cols = "T1_1AGE60_",
    age_midpoint = 62.5
  ),
  "65-69" = list(
    female_cols = "T1_1AGE651",
    male_cols = "T1_1AGE65_",
    age_midpoint = 67.5
  ),
  "70-74" = list(
    female_cols = "T1_1AGE701",
    male_cols = "T1_1AGE70_",
    age_midpoint = 72.5
  ),
  "75-79" = list(
    female_cols = "T1_1AGE751",
    male_cols = "T1_1AGE75_",
    age_midpoint = 77.5
  ),
  "80-84" = list(
    female_cols = "T1_1AGE801",
    male_cols = "T1_1AGE80_",
    age_midpoint = 82.5
  ),
  "85+" = list(
    female_cols = "T1_1AGEGE1",
    male_cols = "T1_1AGEGE_",
    age_midpoint = 87
  )
)

# Calculate stratified YLL
calculate_stratified_yll <- function(data, paf_col, age_range = "all") {
  total_yll <- numeric(nrow(data))
  
  # Determine age bands
  if (age_range == "30+") {
    bands_to_use <- mortality_rates$AgeGroup[mortality_rates$AgeGroup >= "30-34"]
  } else {
    bands_to_use <- mortality_rates$AgeGroup
  }
  
  # Loop through each age-sex band
  for (i in 1:nrow(mortality_rates)) {
    age_group <- mortality_rates$AgeGroup[i]
    
    if (!(age_group %in% bands_to_use)) next
    
    # Get mortality rates
    mort_f <- mortality_rates$Female[i]
    mort_m <- mortality_rates$Male[i]
    
    # Get life expectancy
    age_mid <- age_band_mapping[[age_group]]$age_midpoint
    le_f <- get_life_expectancy(age_mid, "Female", life_table)
    le_m <- get_life_expectancy(age_mid, "Male", life_table)
    
    # Get population
    female_cols <- age_band_mapping[[age_group]]$female_cols
    male_cols <- age_band_mapping[[age_group]]$male_cols
    
    pop_f <- rowSums(data[, female_cols, drop = FALSE], na.rm = TRUE)
    pop_m <- rowSums(data[, male_cols, drop = FALSE], na.rm = TRUE)
    
    # Calculate YLL for this stratum
    yll_stratum <- (pop_f * mort_f * le_f) + (pop_m * mort_m * le_m)
    total_yll <- total_yll + yll_stratum
  }
  
  # Apply PAF
  total_yll <- total_yll * data[[paf_col]]
  
  return(total_yll)
}

# Calculate air pollution YLL (age 30+ for PM2.5/NO2, all ages for O3)
air$YLLPM25 <- calculate_stratified_yll(air, "PAFPM25", age_range = "30+")
air$YLLNO2 <- calculate_stratified_yll(air, "PAFNO2", age_range = "30+")
air$YLLO3 <- calculate_stratified_yll(air, "PAFO3", age_range = "all")

air$YLLAQ <- air$YLLPM25 + air$YLLNO2 + air$YLLO3

cat("  Total YLL from air pollution:", round(sum(air$YLLAQ, na.rm=TRUE), 2), "\n")

## NOISE YLL ##
cat("Noise YLL...\n")

# Process noise data
noise_d <- noise %>%
  spread(DBs, area) %>%
  select(SA_PUB2022, `55`, `60`, `65`, `70`, `75`)

noise_w <- noise_d %>%
  group_by(SA_PUB2022) %>% 
  summarise(across(everything(), ~ reduce(., coalesce)), .groups = 'drop')

air <- merge(air, noise_w, by = "SA_PUB2022", all.x = TRUE)

# Proportion of area exposed to each noise level
air$D55_P <- air$`55`/air$SHAPE_Area
air$D60_P <- air$`60`/air$SHAPE_Area
air$D65_P <- air$`65`/air$SHAPE_Area
air$D70_P <- air$`70`/air$SHAPE_Area
air$D75_P <- air$`75`/air$SHAPE_Area

air[, c("D55_P", "D60_P", "D65_P", "D70_P", "D75_P")][
  is.na(air[, c("D55_P", "D60_P", "D65_P", "D70_P", "D75_P")])] <- 0

# Odds ratios for noise levels (WHO 2018)
air$OR55 <- ifelse(air$D55_P==0, 1, 1.000195)
air$OR60 <- ifelse(air$D60_P==0, 1, 1.01296)
air$OR65 <- ifelse(air$D65_P==0, 1, 1.061315)
air$OR70 <- ifelse(air$D70_P==0, 1, 1.15078)
air$OR75 <- ifelse(air$D75_P==0, 1, 1.286875)

# IHD mortality rates (Ireland 2017, CSO DHA11)
ihd_mortality_rates <- data.frame(
  AgeGroup = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", 
               "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", 
               "75-79", "80-84", "85+"),
  Female = c(0, 0, 0, 0, 0, 0, 0.0000055, 0.000005, 0.0000273, 0.0000539,
             0.0001371, 0.0001954, 0.0004435, 0.0006484, 0.0012399, 0.0028503,
             0.0061533, 0.0210593),
  Male = c(0, 0, 0, 0, 0, 0, 0.000018, 0.0000727, 0.0001722, 0.0003309,
           0.0005262, 0.000979, 0.0014019, 0.002289, 0.0038114, 0.0069162,
           0.0112589, 0.0286301)
)

# Calculate stratified noise YLL
calculate_noise_yll <- function(data, or_col, prop_exposed_col) {
  total_yll <- numeric(nrow(data))
  
  for (i in 1:nrow(ihd_mortality_rates)) {
    age_group <- ihd_mortality_rates$AgeGroup[i]
    
    # Get IHD mortality rates
    mort_ihd_f <- ihd_mortality_rates$Female[i]
    mort_ihd_m <- ihd_mortality_rates$Male[i]
    
    if (mort_ihd_f == 0 && mort_ihd_m == 0) next
    
    # Get life expectancy
    age_mid <- age_band_mapping[[age_group]]$age_midpoint
    le_f <- get_life_expectancy(age_mid, "Female", life_table)
    le_m <- get_life_expectancy(age_mid, "Male", life_table)
    
    # Get population
    female_cols <- age_band_mapping[[age_group]]$female_cols
    male_cols <- age_band_mapping[[age_group]]$male_cols
    
    pop_f <- rowSums(data[, female_cols, drop = FALSE], na.rm = TRUE)
    pop_m <- rowSums(data[, male_cols, drop = FALSE], na.rm = TRUE)
    
    # Calculate YLL: (OR - 1) × IHD_mortality × Population × Proportion_exposed × LE
    yll_stratum <- ((data[[or_col]] - 1) * mort_ihd_f * pop_f * data[[prop_exposed_col]] * le_f) +
      ((data[[or_col]] - 1) * mort_ihd_m * pop_m * data[[prop_exposed_col]] * le_m)
    
    total_yll <- total_yll + yll_stratum
  }
  
  return(total_yll)
}

# Calculate noise YLL for each exposure level
air$YLL55 <- calculate_noise_yll(air, "OR55", "D55_P")
air$YLL60 <- calculate_noise_yll(air, "OR60", "D60_P")
air$YLL65 <- calculate_noise_yll(air, "OR65", "D65_P")
air$YLL70 <- calculate_noise_yll(air, "OR70", "D70_P")
air$YLL75 <- calculate_noise_yll(air, "OR75", "D75_P")

air$YLLN <- air$YLL55 + air$YLL60 + air$YLL65 + air$YLL70 + air$YLL75

cat("  Total YLL from noise pollution:", round(sum(air$YLLN, na.rm=TRUE), 2), "\n")

# Replace NA with 0
air$YLLAQ[is.na(air$YLLAQ)] <- 0
air$YLLN[is.na(air$YLLN)] <- 0

# Save total YLL summary
total_yll <- data.frame(
  Source = c("Air Pollution (PM2.5, NO2, O3)", "Road Noise (55-75+ dB)"),
  Total_YLL = c(sum(air$YLLAQ, na.rm=TRUE), sum(air$YLLN, na.rm=TRUE))
)
print(total_yll)
write.csv(total_yll, "output/total_years_life_lost.csv", row.names = FALSE)

# Calculate population-weighted average age
age_cols_all <- c(paste0("T1_1AGE", 0:19, "M"), 
                  paste0("T1_1AGE", c(20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80), "_"),
                  "T1_1AGEGE_",
                  paste0("T1_1AGE", 0:19, "F"),
                  paste0("T1_1AGE", c(20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80), "1"),
                  "T1_1AGEGE1")

age_midpoints_all <- c(0:19, 22, 27, 32, 37, 42, 47, 52, 57, 62, 67, 72, 77, 82, 87,
                       0:19, 22, 27, 32, 37, 42, 47, 52, 57, 62, 67, 72, 77, 82, 87)

age_weighted <- air[, age_cols_all] * rep(age_midpoints_all, each = nrow(air))
weighted_sum <- rowSums(age_weighted, na.rm = TRUE)
total_pop <- rowSums(air[, age_cols_all], na.rm = TRUE)

air$avg_age <- weighted_sum / total_pop

# Simplify air data for merging
air <- air %>%
  select(SA_PUB2022, SHAPE_Area, T1_1AGETT, T6_2_TH,
         YLLPM25, YLLNO2, YLLO3, YLLAQ, YLL55, YLL60, YLL65, YLL70, YLL75, YLLN,
         T10_4_OD_2, T10_4_HD_2, T10_4_PDT, T10_4_DT, T10_4_TT, T11_1_FT,
         T11_1_BIT, T11_1_BUT, T11_1_TDLT, T11_1_TT, T12_3_BT, T12_3_VBT,
         T12_3_TT, T1_1AGE6_1, T1_1AGE6_2, T1_1AGE7_1, T1_1AGE7_2,
         T1_1AGE8_1, T1_1AGEG_1, T2_1IEC, T2_1TC, avg_age)

################################################################################
# 8. INDEX WEIGHTING APPROACHES
################################################################################

cat("\n=== CALCULATING HEALTH RATING INDEX ===\n")

# Load spatial data
SA <- st_read("data/SA2022_Dublin_AllData3.shp", quiet = TRUE) %>%
  select(SA_PUB2022, RoutingKey, In22_ED, SOUTH, log_dist, total)

SA <- SA %>%
  left_join(sa_access, by = c("SA_PUB2022" = "SA_PUB202")) %>%
  left_join(air, by = "SA_PUB2022")

# Create rate variables
SA$YLLAQPC <- (SA$YLLAQ / SA$T1_1AGETT) * 100000
SA$YLLNPC <- (SA$YLLN / SA$T1_1AGETT) * 100000
SA$YLLPC <- ((SA$YLLN + SA$YLLAQ) / SA$T1_1AGETT) * 100000
SA$ph_rate <- SA$total / SA$T6_2_TH
SA$ph_rate[is.infinite(SA$ph_rate) | is.nan(SA$ph_rate)] <- 0
SA$ph_rate <- pmin(SA$ph_rate, 1)

# Calculate excess YLL (age-adjusted burden)
SA$avg_age_z <- scale(SA$avg_age)
SA$YLLPC_z <- scale(SA$YLLPC)
SA$excess_YLLPC <- SA$YLLPC_z - SA$avg_age_z

# Function: Naive weights (equal)
calc_naive <- function(df, gp_var) {
  df_data <- st_drop_geometry(df)
  
  df_norm <- df_data %>%
    select(all_of(c(gp_var, "gs_access", "bs_access", "YLLAQPC", "YLLNPC", "ph_rate"))) %>%
    mutate(across(everything(), ~ (. - min(.)) / (max(.) - min(.))))
  
  df_norm <- df_norm %>%
    mutate(YLLAQPC = 1 - YLLAQPC, YLLNPC = 1 - YLLNPC, ph_rate = 1 - ph_rate)
  
  hri <- (1/6) * (df_norm[[gp_var]] + df_norm$gs_access + df_norm$bs_access + 
                    df_norm$YLLAQPC + df_norm$YLLNPC + df_norm$ph_rate)
  
  (hri - min(hri)) / (max(hri) - min(hri))
}

# Function: Entropy weights (Karagiannis & Karagiannis 2020)
calc_entropy <- function(df, gp_var) {
  df_data <- st_drop_geometry(df)
  df_raw <- df_data %>%
    select(all_of(c(gp_var, "gs_access", "bs_access", "YLLAQPC", "YLLNPC", "ph_rate")))
  
  # Robust normalization (5th-95th percentile)
  df_norm <- df_raw %>%
    mutate(across(everything(), ~ {
      p05 <- quantile(., 0.05, na.rm = TRUE)
      p95 <- quantile(., 0.95, na.rm = TRUE)
      x_trimmed <- pmax(pmin(., p95), p05)
      (x_trimmed - p05) / (p95 - p05)
    }))
  
  df_norm <- df_norm %>%
    mutate(YLLAQPC = 1 - YLLAQPC, YLLNPC = 1 - YLLNPC, ph_rate = 1 - ph_rate)
  
  # Shannon entropy weights
  p <- sweep(df_norm, 2, colSums(df_norm), "/")
  
  n <- nrow(df_norm)
  k <- 1 / log(n)
  p[p == 0] <- 1e-10
  
  e <- -k * colSums(p * log(p))
  d <- 1 - e
  w <- d / sum(d)
  
  index <- as.numeric(as.matrix(df_norm) %*% w)
  (index - min(index)) / (max(index) - min(index))
}

# Function: PCA variance weights (OECD 2008)
calc_pca <- function(df, gp_var) {
  df_data <- st_drop_geometry(df)
  df_pca <- df_data %>%
    select(all_of(c(gp_var, "gs_access", "bs_access", "YLLAQPC", "YLLNPC", "ph_rate"))) %>%
    mutate(YLLAQPC = max(YLLAQPC, na.rm = TRUE) - YLLAQPC,
           YLLNPC = max(YLLNPC, na.rm = TRUE) - YLLNPC,
           ph_rate = max(ph_rate, na.rm = TRUE) - ph_rate) %>%
    na.omit()
  
  row_ids <- as.numeric(rownames(df_pca))
  df_pca_scaled <- as.data.frame(scale(df_pca))
  
  pca_res <- prcomp(df_pca_scaled, center = FALSE, scale. = FALSE)
  eig_vals <- (pca_res$sdev)^2
  prop_variance <- eig_vals / sum(eig_vals)
  retain <- which(eig_vals > 1)  # Kaiser criterion
  
  # Calculate variable weights
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
  df_data %>%
    mutate(row_id = row_number()) %>%
    left_join(result, by = "row_id") %>%
    pull(hri)
}

# Generate all 12 HRI variants (4 GP methods × 3 weighting approaches)
gp_vars <- c("mz_2s", "gaus_2s", "kern_2s", "bin_2s")
methods <- c("naive", "entropy", "pca")

SA_results <- SA
for (gp_var in gp_vars) {
  for (method in methods) {
    col_name <- paste0("HRI_", substr(gp_var, 1, nchar(gp_var)-3), "_", 
                       substr(method, 1, 1))
    
    cat(paste0("  ", col_name, "\n"))
    
    SA_results[[col_name]] <- switch(method,
                                     "naive" = calc_naive(SA, gp_var),
                                     "entropy" = calc_entropy(SA, gp_var),
                                     "pca" = calc_pca(SA, gp_var))
  }
}

# Calculate correlation matrix
cat("\nCalculating HRI variant correlations...\n")
index_cols <- grep("^HRI_", names(SA_results), value = TRUE)
cor_matrix <- cor(st_drop_geometry(SA_results)[, index_cols], 
                  method = "spearman", 
                  use = "pairwise.complete.obs")

write.csv(round(cor_matrix, 3), "output/hri_correlation_matrix.csv")

cat("Mean Spearman correlations:\n")
print(round(colMeans(cor_matrix, na.rm = TRUE), 3))

################################################################################
# 9. SPATIAL CROSS-VALIDATION RANDOM FOREST
################################################################################

cat("\n=== SPATIAL RANDOM FOREST MODEL ===\n")

# Create derived variables
SA_results$NoAuto_p <- (SA_results$T11_1_FT + SA_results$T11_1_BIT + 
                        SA_results$T11_1_BUT + SA_results$T11_1_TDLT) / 
                        SA_results$T11_1_TT
SA_results$BVBHth_p <- (SA_results$T12_3_BT + SA_results$T12_3_VBT) / 
                        SA_results$T12_3_TT
SA_results$ov60 <- (SA_results$T1_1AGE6_1 + SA_results$T1_1AGE6_2 + 
                    SA_results$T1_1AGE7_1 + SA_results$T1_1AGE7_2 +
                    SA_results$T1_1AGE8_1 + SA_results$T1_1AGEG_1) / SA_results$T1_1AGETT
SA_results$nonIrish <- 1 - (SA_results$T2_1IEC / SA_results$T2_1TC)
SA_results$POPD <- SA_results$T1_1AGETT / SA_results$SHAPE_Area

# Prepare model data
model_data <- SA_results %>%
  filter(!is.na(HRI_gaus_p)) %>%
  select(HRI_gaus_p, In22_ED, NoAuto_p, POPD, log_dist, ov60, nonIrish)

# Standardize variables
model_data <- model_data %>%
  mutate(across(-geometry, ~scale(.)[,1]))

# Run SArf (Spatial Random Forest with within-fold spatial CV)
cat("Running spatial cross-validation random forest...\n")
cat("This may take several minutes...\n\n")

results <- SArf(
  HRI_gaus_p ~ In22_ED + NoAuto_p + POPD + log_dist + ov60 + nonIrish,
  data = model_data,
  k_neighbors = 20,            # Neighbors for spatial lag
  n_folds = 5,                 # Spatial CV folds
  n_bootstrap = 20,            # Bootstrap iterations for CIs
  num_trees = 500,             # Trees in random forest
  include_naive_rf = TRUE,     # Compare with naive RF
  naive_test_fraction = 0.2,   
  create_map = TRUE,           # Map of spatial folds
  block_range = 5000,          # Block size (meters)
  verbose = TRUE
)

# Display results
cat("\n=== MODEL RESULTS ===\n")
print(results$model_comparison)

cat("\n=== VARIABLE IMPORTANCE (with 95% CIs) ===\n")
# print(results$variable_importance)

# View plots
# print(results$moran_plot)
print(results$importance_plot)
print(results$ale_plots)

# Save detailed spatial model results
show_models(results)

################################################################################
# 10. INTERACTIVE MAPS
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
    label = ~paste0("YLL (per 100K): ", round(YLLAQPC, 2))
  ) %>%
  addLegend(
    pal = palaq,
    values = ~YLLAQPC,
    title = "YLL Air Pollution",
    position = "bottomright"
  )

# Map 6: Noise YLL
SA_results_filtered <- SA_results %>% filter(!is.na(YLLNPC) & YLLNPC > 0)

paln <- colorQuantile("Reds", SA_results_filtered$YLLNPC, n = 5)
map_noise <- leaflet(SA_results_filtered) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolygons(
    fillColor = ~paln(YLLNPC),
    fillOpacity = 0.5,
    color = "white",
    weight = 0.1,
    label = ~paste0("YLL (per 100K): ", round(YLLNPC, 2))
  ) %>%
  addLegend(
    pal = paln,
    values = ~YLLNPC,
    title = "YLL Road Noise",
    position = "bottomright"
  )

# Map 7: Poor Quality Housing
SA_results$ph_rate[is.na(SA_results$ph_rate) | is.infinite(SA_results$ph_rate)] <- 0

palhq <- colorQuantile("viridis", SA_results$ph_rate, n = 5)
map_hq <- leaflet(SA_results) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolygons(
    fillColor = ~palhq(ph_rate),
    fillOpacity = 0.5,
    color = "white",
    weight = 0.1,
    label = ~paste0("BER E or worse: ", round(ph_rate, 2))
  ) %>%
  addLegend(
    pal = palhq,
    values = ~ph_rate,
    title = "Poor Quality Housing",
    position = "bottomright"
  )

# Map 8: Excess YLL Risk (age-adjusted)
sd_bins <- c(-Inf, -2, -1, 0, 1, 2, Inf)
palen <- colorBin("RdYlBu", SA_results$excess_YLLPC, bins = sd_bins, reverse = TRUE)
map_excess <- leaflet(SA_results) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolygons(
    fillColor = ~palen(excess_YLLPC),
    fillOpacity = 0.5,
    color = "white",
    weight = 0.1,
    label = ~paste0("Excess YLL: ", round(excess_YLLPC, 2))
  ) %>%
  addLegend(
    pal = palen,
    values = ~excess_YLLPC,
    title = "Excess YLL",
    position = "bottomright"
  )

# Save maps
saveWidget(map_hri, "output/map_hri.html")
saveWidget(map_gs, "output/map_green_space.html")
saveWidget(map_bs, "output/map_blue_space.html")
saveWidget(map_gp, "output/map_gp_access.html")
saveWidget(map_aq, "output/map_air_pollution_yll.html")
saveWidget(map_noise, "output/map_noise_yll.html")
saveWidget(map_hq, "output/map_housing_quality.html")
saveWidget(map_excess, "output/map_excess_yll.html")

################################################################################
# ANALYSIS COMPLETE
################################################################################

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Results saved to output/ directory:\n")
cat("  - total_years_life_lost.csv\n")
cat("  - hri_correlation_matrix.csv\n")
cat("  - Interactive HTML maps\n")
cat("  - Model results (from SArf package)\n")
