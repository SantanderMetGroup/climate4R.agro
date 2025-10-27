################################################################################
# CREATE SYNTHETIC PRECIPITATION AND TEMPERATURE DATA
# 6 complete years (2015-2020) for testing agroindexGrid
################################################################################

library(transformeR)

message("[", Sys.time(), "] Creating synthetic climate data...")

# Set working directory
setwd("C:/Users/pablo/Desktop/TrabajoSara/example_data")

# Load original grid to get spatial structure
pr_original <- readRDS("pr_mildiu_obs_vid.rds")

################################################################################
# CREATE COMPLETE DAILY TIME SERIES (2015-2020)
################################################################################

# Create complete daily time series
start_date <- as.Date("2015-01-01")
end_date <- as.Date("2020-12-31")
dates_complete <- seq(start_date, end_date, by = "day")
n_days <- length(dates_complete)

message("  - Total days: ", n_days, " (", as.character(start_date), " to ", as.character(end_date), ")")

# Get spatial dimensions from original grid
n_lat <- length(pr_original$xyCoords$y)
n_lon <- length(pr_original$xyCoords$x)
message("  - Grid dimensions: ", n_lat, " lat x ", n_lon, " lon")

################################################################################
# GENERATE SYNTHETIC PRECIPITATION
################################################################################

message("\n[", Sys.time(), "] Generating synthetic precipitation...")

set.seed(123)  # For reproducibility
pr_data <- array(0, dim = c(n_days, n_lat, n_lon))

for (i in 1:n_lat) {
  for (j in 1:n_lon) {
    # Seasonal rain probability (higher in winter)
    doy <- as.numeric(format(dates_complete, "%j"))
    seasonal_factor <- 0.3 + 0.2 * cos((doy - 1) * 2 * pi / 365)  # 0.1 to 0.5
    
    # Random rain occurrence
    rain_days <- runif(n_days) < seasonal_factor
    
    # Rain amounts: exponential distribution with mean ~10mm
    rain_amounts <- rexp(n_days, rate = 1/10)
    
    # Combine: rain only on rain days
    pr_data[, i, j] <- ifelse(rain_days, rain_amounts, 0)
  }
}

# Add spatial variability
for (i in 1:n_lat) {
  lat_factor <- 0.8 + 0.4 * sin(i * pi / n_lat)  # Varies from 0.8 to 1.2
  pr_data[, i, ] <- pr_data[, i, ] * lat_factor
}

message("  - Mean precipitation: ", round(mean(pr_data), 2), " mm/day")
message("  - Max precipitation: ", round(max(pr_data), 2), " mm/day")
message("  - Days with pr >= 1mm: ", round(100 * sum(pr_data >= 1) / length(pr_data), 1), "%")
message("  - Days with pr >= 10mm: ", round(100 * sum(pr_data >= 10) / length(pr_data), 1), "%")

# Create precipitation grid
pr <- pr_original
pr$Data <- pr_data
pr$Dates$start <- dates_complete
pr$Dates$end <- dates_complete
attr(pr$Data, "dimensions") <- c("time", "lat", "lon")

################################################################################
# GENERATE SYNTHETIC TEMPERATURE (tasmin and tasmax)
################################################################################

message("\n[", Sys.time(), "] Generating synthetic temperature...")

tasmin_data <- array(0, dim = c(n_days, n_lat, n_lon))
tasmax_data <- array(0, dim = c(n_days, n_lat, n_lon))

for (i in 1:n_lat) {
  for (j in 1:n_lon) {
    # Day of year for seasonal cycle
    doy <- as.numeric(format(dates_complete, "%j"))
    
    # Seasonal temperature cycle (winter ~0°C, summer ~22°C for tasmin)
    # This ensures we have cold periods below 5°C for GSL calculation
    seasonal_mean <- 11 + 11 * sin((doy - 80) * 2 * pi / 365)
    
    # Add random daily variation
    daily_noise <- rnorm(n_days, mean = 0, sd = 2.5)
    
    # Add spatial variation (cooler at higher latitudes)
    lat_effect <- -(i - n_lat/2) * 0.3
    
    # Tasmin
    tasmin_data[, i, j] <- seasonal_mean + daily_noise + lat_effect
    
    # Tasmax (typically 8-12°C warmer than tasmin)
    diurnal_range <- 8 + rnorm(n_days, mean = 2, sd = 1.5)
    tasmax_data[, i, j] <- tasmin_data[, i, j] + diurnal_range
  }
}

message("  - Mean tasmin: ", round(mean(tasmin_data), 2), " °C")
message("  - Mean tasmax: ", round(mean(tasmax_data), 2), " °C")
message("  - Mean diurnal range: ", round(mean(tasmax_data - tasmin_data), 2), " °C")
message("  - Tasmin range: ", round(min(tasmin_data), 2), " to ", round(max(tasmin_data), 2), " °C")
message("  - Tasmax range: ", round(min(tasmax_data), 2), " to ", round(max(tasmax_data), 2), " °C")

# Create temperature grids
tasmin <- pr_original  # Use same structure
tasmin$Data <- tasmin_data
tasmin$Dates$start <- dates_complete
tasmin$Dates$end <- dates_complete
tasmin$Variable$varName <- "tasmin"
attr(tasmin$Data, "dimensions") <- c("time", "lat", "lon")

tasmax <- pr_original
tasmax$Data <- tasmax_data
tasmax$Dates$start <- dates_complete
tasmax$Dates$end <- dates_complete
tasmax$Variable$varName <- "tasmax"
attr(tasmax$Data, "dimensions") <- c("time", "lat", "lon")

################################################################################
# VERIFY DATA HAS COMPLETE YEARS
################################################################################

message("\n[", Sys.time(), "] Verifying complete years...")

years <- unique(format(dates_complete, "%Y"))
for(y in years) {
  jan1 <- paste0(y, "-01-01")
  dec31 <- paste0(y, "-12-31")
  has_jan1 <- jan1 %in% as.character(dates_complete)
  has_dec31 <- dec31 %in% as.character(dates_complete)
  n_days_year <- sum(format(dates_complete, "%Y") == y)
  message("  Year ", y, ": ", n_days_year, " days (Jan 1: ", has_jan1, ", Dec 31: ", has_dec31, ")")
}

################################################################################
# SAVE SYNTHETIC DATA
################################################################################

message("\n[", Sys.time(), "] Saving synthetic data...")

saveRDS(pr, "pr_synthetic_6years.rds")
saveRDS(tasmin, "tasmin_synthetic_6years.rds")
saveRDS(tasmax, "tasmax_synthetic_6years.rds")

message("[", Sys.time(), "] Done! Synthetic data saved:")
message("  - pr_synthetic_6years.rds")
message("  - tasmin_synthetic_6years.rds")
message("  - tasmax_synthetic_6years.rds")
message("\nUse these files in test.R by loading:")
message('  pr <- readRDS("pr_synthetic_6years.rds")')
message('  tasmin <- readRDS("tasmin_synthetic_6years.rds")')
message('  tasmax <- readRDS("tasmax_synthetic_6years.rds")')

