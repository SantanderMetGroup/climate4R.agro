################################################################################
# TEST SCRIPT FOR agroindexGrid.R
# Comprehensive validation of agroclimatic indices calculation
# Author: P. Lavin
# Date: 2025
################################################################################

library(transformeR)
library(magrittr)
library(dplyr)
library(abind)

source("C:/Users/pablo/Desktop/climate4R.agro/R/create_synthetic_data.R")
source("C:/Users/pablo/Desktop/TrabajoSara/example_data/create_realistic_temp_data.R")

# Source the climate4R.agro functions (development version)
message("Loading climate4R.agro functions...")
source("C:/Users/pablo/Desktop/climate4R.agro/R/agroindexGrid.R")
source("C:/Users/pablo/Desktop/climate4R.agro/R/indicesFAO_tier1.R")
source("C:/Users/pablo/Desktop/climate4R.agro/R/indicesFAO.R")
source("C:/Users/pablo/Desktop/climate4R.agro/R/cdi_cei_functions.R")

# Set working directory to data folder
setwd("C:/Users/pablo/Desktop/TrabajoSara/example_data")

################################################################################
# 1. LOAD DATA
################################################################################

message("[", Sys.time(), "] Loading climate data...")

# Load precipitation and temperature data (SYNTHETIC with complete years)
pr <- readRDS("pr_synthetic_6years.rds")
# No need for unit conversion - synthetic data already in mm

# Check actual dates
dates_vec <- as.Date(pr$Dates$start)
cat("First date:", as.character(head(dates_vec, 1)), "\n")
cat("Last date:", as.character(tail(dates_vec, 1)), "\n")
cat("First few dates:", paste(head(dates_vec, 10), collapse=", "), "\n")
cat("Unique years:", paste(unique(format(dates_vec, "%Y")), collapse=", "), "\n")

# Check if Jan 1 exists for each year
years <- unique(format(dates_vec, "%Y"))
for(y in years) {
  jan1 <- paste0(y, "-01-01")
  dec31 <- paste0(y, "-12-31")
  has_jan1 <- jan1 %in% as.character(dates_vec)
  has_dec31 <- dec31 %in% as.character(dates_vec)
  cat("Year", y, "- Has Jan 1:", has_jan1, ", Has Dec 31:", has_dec31, "\n")
}

tasmin <- readRDS("tasmin_synthetic_6years.rds")
message("    Temperature range: ", paste(round(range(tasmin$Data, na.rm=TRUE), 1), collapse=" to "), " °C")

# Fix dimensions if needed (remove extra singleton dimensions)
if (length(dim(tasmin$Data)) > 3) {
  message("    Fixing tasmin dimensions from ", length(dim(tasmin$Data)), " to 3...")
  tasmin_data_fixed <- array(tasmin$Data, dim = c(length(tasmin$Dates$start), 
                                                    length(tasmin$xyCoords$y), 
                                                    length(tasmin$xyCoords$x)))
  tasmin$Data <- tasmin_data_fixed
  attr(tasmin$Data, "dimensions") <- c("time", "lat", "lon")
  message("    Fixed dimensions: ", paste(dim(tasmin$Data), collapse=" x "))
}

# Inspect data structure
message("\n=== DATA STRUCTURE ===")
message("Precipitation structure:")
str(pr, max.level = 2)
message("\nMinimum temperature structure:")
str(tasmin, max.level = 2)

# Check dimensions
message("\nPrecipitation dimensions: ", paste(dim(pr$Data), collapse = " x "))
message("Tasmin dimensions: ", paste(dim(tasmin$Data), collapse = " x "))

# Check time period
message("\nTime period:")
message("Start: ", pr$Dates$start[1])
message("End: ", tail(pr$Dates$end, 1))
message("Number of days: ", length(pr$Dates$start))

# Check spatial extent
message("\nSpatial extent:")
message("Latitudes: ", paste(range(pr$xyCoords$y), collapse = " to "))
message("Longitudes: ", paste(range(pr$xyCoords$x), collapse = " to "))
message("Grid points: ", length(pr$xyCoords$y), " x ", length(pr$xyCoords$x))

################################################################################
# 2. TEST FAO TIER1 INDICES
################################################################################

message("\n\n=== TESTING FAO TIER1 INDICES ===\n")

# -----------------------------------------------------------------------------
# Test 2.1: Number of Rainy Days (NRD)
# -----------------------------------------------------------------------------
message("[", Sys.time(), "] Testing index: NRD (Number of Rainy Days)")

# Using agroindexGrid
nrd_grid <- agroindexGrid(
  index.code = "nrd",
  pr = pr,
  index.arg.list = list(
    wet.threshold = 1  # days with >= 1mm
  )
)

message("NRD grid structure:")
str(nrd_grid, max.level = 2)
message("NRD dimensions: ", paste(dim(nrd_grid$Data), collapse = " x "))
message("NRD range: ", paste(round(range(nrd_grid$Data, na.rm = TRUE), 2), collapse = " to "))

# Manual verification for one grid point with actual data (indices [8, 14])
message("\n--- Manual verification for NRD at test grid point [8, 14] ---")

# Use direct indices
lats <- pr$xyCoords$y
lons <- pr$xyCoords$x
lat_idx <- 8
lon_idx <- 14

message("Test point: lat_idx=", lat_idx, " (lat=", lats[lat_idx], "), lon_idx=", lon_idx, " (lon=", lons[lon_idx], ")")
summary(pr$Data[, lat_idx, lon_idx])
# Extract data for test grid point
pr_point_test <- pr$Data[, lat_idx, lon_idx]  # time x lat x lon
dates_mat <- cbind(
  as.numeric(format(as.Date(pr$Dates$start), "%Y")),
  as.numeric(format(as.Date(pr$Dates$start), "%m")),
  as.numeric(format(as.Date(pr$Dates$start), "%d"))
)

message("First 10 pr values at test point: ", paste(head(pr_point_test, 10), collapse = ", "))
message("Number of non-NA values: ", sum(!is.na(pr_point_test)))

# Calculate NRD manually
nrd_manual <- nrd(pr = pr_point_test, dates = dates_mat, wet.threshold = 1)
nrd_grid_test <- nrd_grid$Data[, lat_idx, lon_idx]

message("Manual calculation result: ", paste(nrd_manual, collapse = ", "))
message("agroindexGrid result: ", paste(nrd_grid_test, collapse = ", "))
message("Match: ", all.equal(nrd_manual, nrd_grid_test, tolerance = 1e-6))

# -----------------------------------------------------------------------------
# Test 2.2: Total Precipitation in Wet Days (PRCPTOT)
# -----------------------------------------------------------------------------
message("\n[", Sys.time(), "] Testing index: PRCPTOT (Total Precipitation)")

prcptot_grid <- agroindexGrid(
  index.code = "prcptot",
  pr = pr,
  index.arg.list = list(
    wet.threshold = 1
  )
)

message("PRCPTOT range: ", paste(round(range(prcptot_grid$Data, na.rm = TRUE), 2), collapse = " to "), " mm")

# Manual verification at test point
prcptot_manual <- prcptot(pr = pr_point_test, dates = dates_mat, wet.threshold = 1)
prcptot_grid_test <- prcptot_grid$Data[, lat_idx, lon_idx]

message("Manual calculation: ", paste(round(prcptot_manual, 2), collapse = ", "))
message("agroindexGrid result: ", paste(round(prcptot_grid_test, 2), collapse = ", "))
message("Match: ", all.equal(prcptot_manual, prcptot_grid_test, tolerance = 1e-6))

# -----------------------------------------------------------------------------
# Test 2.3: Simple Daily Intensity Index (SDII)
# -----------------------------------------------------------------------------
message("\n[", Sys.time(), "] Testing index: SDII (Precipitation Intensity)")

sdii_grid <- agroindexGrid(
  index.code = "sdii",
  pr = pr,
  index.arg.list = list(
    wet.threshold = 1
  )
)

message("SDII range: ", paste(round(range(sdii_grid$Data, na.rm = TRUE), 2), collapse = " to "), " mm/day")

# Manual verification at test point
sdii_manual <- sdii(pr = pr_point_test, dates = dates_mat, wet.threshold = 1)
sdii_grid_test <- sdii_grid$Data[, lat_idx, lon_idx]

message("Manual calculation: ", paste(round(sdii_manual, 2), collapse = ", "))
message("agroindexGrid result: ", paste(round(sdii_grid_test, 2), collapse = ", "))
message("Match: ", all.equal(sdii_manual, sdii_grid_test, tolerance = 1e-6))

################################################################################
# 3. TEST COMPOSITE INDICES (requiring multiple variables)
################################################################################

message("\n\n=== TESTING COMPOSITE INDICES ===\n")

# Load tasmax from synthetic data
tasmax <- readRDS("tasmax_synthetic_6years.rds")

# Fix dimensions if needed (remove extra singleton dimensions)
if (length(dim(tasmax$Data)) > 3) {
  message("    Fixing tasmax dimensions from ", length(dim(tasmax$Data)), " to 3...")
  tasmax_data_fixed <- array(tasmax$Data, dim = c(length(tasmax$Dates$start), 
                                                    length(tasmax$xyCoords$y), 
                                                    length(tasmax$xyCoords$x)))
  tasmax$Data <- tasmax_data_fixed
  attr(tasmax$Data, "dimensions") <- c("time", "lat", "lon")
  message("    Fixed dimensions: ", paste(dim(tasmax$Data), collapse=" x "))
}

# -----------------------------------------------------------------------------
# Test 3.1: Diurnal Temperature Range (DR)
# -----------------------------------------------------------------------------
message("[", Sys.time(), "] Testing index: DR (Diurnal Range)")

dr_grid <- agroindexGrid(
  index.code = "dr",
  tn = tasmin,
  tx = tasmax
)

message("DR range: ", paste(round(range(dr_grid$Data, na.rm = TRUE), 2), collapse = " to "), " °C")

# Manual verification at test point
tn_point_test <- tasmin$Data[, lat_idx, lon_idx]
tx_point_test <- tasmax$Data[, lat_idx, lon_idx]

dr_manual <- dr(tx = tx_point_test, tn = tn_point_test, dates = dates_mat)
dr_grid_test <- dr_grid$Data[, lat_idx, lon_idx]

message("Manual calculation: ", paste(round(dr_manual, 2), collapse = ", "))
message("agroindexGrid result: ", paste(round(dr_grid_test, 2), collapse = ", "))
message("Match: ", all.equal(dr_manual, dr_grid_test, tolerance = 1e-6))

################################################################################
# 4. TEST FAO AGRONOMIC SEASON INDICES
################################################################################

message("\n\n=== TESTING FAO AGRONOMIC SEASON INDICES ===\n")

# -----------------------------------------------------------------------------
# Test 4.1: Growing Season Length (GSL)
# -----------------------------------------------------------------------------
message("[", Sys.time(), "] Testing index: GSL (Growing Season Length)")

t2m <- readRDS("C:/Users/pablo/Desktop/TrabajoSara/example_data/t2m_synthetic_6years.rds")

# DIAGNOSTIC: Check t2m structure
message("=== DIAGNOSING T2M STRUCTURE ===")
message("Data array dimensions: ", paste(dim(t2m$Data), collapse=" x "))
message("Data class: ", class(t2m$Data))
if (!is.null(attr(t2m$Data, "dimensions"))) {
  message("Dimension attribute: ", paste(attr(t2m$Data, "dimensions"), collapse=", "))
}
message("Grid type: ", typeofGrid(t2m))
str(t2m, max.level = 1)

# Fix: Remove extra dimensions by recreating the data array
message("\nFixing data array structure...")
t2m_data_fixed <- array(t2m$Data, dim = c(length(t2m$Dates$start), 
                                            length(t2m$xyCoords$y), 
                                            length(t2m$xyCoords$x)))
t2m$Data <- t2m_data_fixed
attr(t2m$Data, "dimensions") <- c("time", "lat", "lon")

message("Fixed dimensions: ", paste(dim(t2m$Data), collapse=" x "))
message("================================\n")

gsl_grid <- agroindexGrid(
  index.code = "gsl",
  t2m = t2m,
  index.arg.list = list(
    pnan = 25  # maximum % of missing data allowed
  )
)

message("GSL structure:")
str(gsl_grid, max.level = 2)
message("GSL range: ", paste(round(range(gsl_grid$Data, na.rm = TRUE), 2), collapse = " to "), " days")

# Manual verification at test point [8, 14]
tm_point_test <- t2m$Data[, lat_idx, lon_idx]
lat_test <- t2m$xyCoords$y[lat_idx]

# DIAGNOSTIC: Check temperature patterns for GSL
message("\n--- DIAGNOSTIC: Temperature patterns for GSL ---")
message("Temperature range: ", paste(round(range(tm_point_test, na.rm = TRUE), 2), collapse = " to "), " °C")
message("Mean temperature: ", round(mean(tm_point_test, na.rm = TRUE), 2), " °C")
message("Days with temp > 5°C: ", sum(tm_point_test > 5, na.rm = TRUE))
message("Days with temp < 5°C: ", sum(tm_point_test < 5, na.rm = TRUE))

# Check for 6-day spells
warm_spell <- binSpell(tm_point_test > 5)
cold_spell <- binSpell(tm_point_test < 5)
message("Longest warm spell (>5°C): ", max(warm_spell$len[warm_spell$val], na.rm = TRUE), " days")
message("Longest cold spell (<5°C): ", max(cold_spell$len[cold_spell$val], na.rm = TRUE), " days")
message("Number of 6+ day warm spells: ", sum(warm_spell$len[warm_spell$val] >= 6, na.rm = TRUE))
message("Number of 6+ day cold spells: ", sum(cold_spell$len[cold_spell$val] >= 6, na.rm = TRUE))

# Check if the pattern exists in the correct order (for Northern Hemisphere)
if (lat_test >= 0) {
  message("\nNorthern Hemisphere - checking for pattern:")
  message("  Need: 6-day warm spell before July 1, then 6-day cold spell after July 1")
  
  # Find July 1 index for first year
  year_2015 <- which(dates_mat[,1] == 2015)
  if (length(year_2015) > 0) {
    july1_idx <- which(dates_mat[,1] == 2015 & dates_mat[,2] == 7 & dates_mat[,3] == 1)
    if (length(july1_idx) > 0) {
      july1_pos <- july1_idx - min(year_2015) + 1
      message("  July 1 is at position: ", july1_pos, " of ", length(year_2015), " days in 2015")
      
      first_half <- tm_point_test[year_2015[1]:year_2015[july1_pos-1]]
      second_half <- tm_point_test[year_2015[july1_pos]:year_2015[length(year_2015)]]
      
      message("  First half (Jan-Jun) temp range: ", paste(round(range(first_half, na.rm=TRUE), 1), collapse=" to "))
      message("  Second half (Jul-Dec) temp range: ", paste(round(range(second_half, na.rm=TRUE), 1), collapse=" to "))
    }
  }
}
message("-----------------------------------------------\n")

gsl_manual <- gsl(tm = tm_point_test, dates = dates_mat, lat = lat_test, pnan = 25)
message("Manual calculation (GSL): ", paste(gsl_manual$GSL, collapse = ", "))
message("agroindexGrid result (GSL): ", paste(gsl_grid$Data[, lat_idx, lon_idx], collapse = ", "))
message("Match: ", all.equal(gsl_manual$GSL, gsl_grid$Data[, lat_idx, lon_idx], tolerance = 1e-6))

# -----------------------------------------------------------------------------
# Test 4.2: Agronomic Season Length (dl_agsn)
# -----------------------------------------------------------------------------
message("\n[", Sys.time(), "] Testing index: DL_AGSN (Agronomic Season Length)")

dl_agsn_grid <- agroindexGrid(
  index.code = "dl_agsn",
  tn = tasmin,
  tx = tasmax,
  pr = pr,
  index.arg.list = list(
    pnan = 25,
    shc = 100,   # soil holding capacity
    rndy = 2.5,  # rainy day threshold
    rnlg = 50    # large rain threshold
  )
)

message("DL_AGSN range: ", paste(round(range(dl_agsn_grid$Data, na.rm = TRUE), 2), collapse = " to "), " days")

# Manual verification at test point
dl_agsn_manual <- agroindexFAO(
  lat = lat_test,
  dates = dates_mat,
  index.code = "dl_agsn",
  pr = pr_point_test,
  tx = tx_point_test,
  tn = tn_point_test,
  pnan = 25,
  shc = 100,
  rndy = 2.5,
  rnlg = 50
)

message("Manual calculation: ", paste(dl_agsn_manual, collapse = ", "))
message("agroindexGrid result: ", paste(dl_agsn_grid$Data[, lat_idx, lon_idx], collapse = ", "))
message("Match: ", all.equal(dl_agsn_manual, dl_agsn_grid$Data[, lat_idx, lon_idx], tolerance = 1e-6))

################################################################################
# 5. SPATIAL VISUALIZATION
################################################################################

message("\n\n=== SPATIAL ANALYSIS ===\n")

if (requireNamespace("visualizeR", quietly = TRUE)) {
  library(visualizeR)
  
  # Plot NRD spatial distribution
  message("Plotting NRD spatial distribution...")
  spatialPlot(
    climatology(nrd_grid),
    main = "Mean Annual Number of Rainy Days",
    backdrop.theme = "countries"
  )
  
  # Plot PRCPTOT spatial distribution
  message("Plotting PRCPTOT spatial distribution...")
  spatialPlot(
    climatology(prcptot_grid),
    main = "Mean Annual Precipitation Total (mm)",
    backdrop.theme = "countries"
  )
  
  # Plot GSL spatial distribution
  message("Plotting GSL spatial distribution...")
  spatialPlot(
    climatology(gsl_grid),
    main = "Mean Growing Season Length (days)",
    backdrop.theme = "countries"
  )
} else {
  message("visualizeR not available for plotting")
}

################################################################################
# 6. VALIDATION SUMMARY
################################################################################

message("\n\n=== VALIDATION SUMMARY ===\n")

# Create summary table
validation_results <- data.frame(
  Index = c("NRD", "PRCPTOT", "SDII", "DR", "GSL", "DL_AGSN"),
  Description = c(
    "Number of Rainy Days",
    "Total Precipitation",
    "Precipitation Intensity",
    "Diurnal Range",
    "Growing Season Length",
    "Agronomic Season Length"
  ),
  Status = c(
    ifelse(exists("nrd_manual") && all.equal(nrd_manual, nrd_grid_point1, tolerance = 1e-6) == TRUE, "PASS", "CHECK"),
    ifelse(exists("prcptot_manual") && all.equal(prcptot_manual, prcptot_grid_point1, tolerance = 1e-6) == TRUE, "PASS", "CHECK"),
    ifelse(exists("sdii_manual") && all.equal(sdii_manual, sdii_grid_point1, tolerance = 1e-6) == TRUE, "PASS", "CHECK"),
    ifelse(exists("dr_manual") && all.equal(dr_manual, dr_grid_point1, tolerance = 1e-6) == TRUE, "PASS", "CHECK"),
    "MANUAL_CHECK",
    ifelse(exists("dl_agsn_manual") && all.equal(dl_agsn_manual, dl_agsn_grid$Data[, 1, 1], tolerance = 1e-6) == TRUE, "PASS", "CHECK")
  ),
  stringsAsFactors = FALSE
)

print(validation_results)

################################################################################
# 7. ADDITIONAL DIAGNOSTICS
################################################################################

message("\n\n=== DIAGNOSTICS ===\n")

# Check for NA values
message("NA values in results:")
message("NRD: ", sum(is.na(nrd_grid$Data)), " / ", length(nrd_grid$Data))
message("PRCPTOT: ", sum(is.na(prcptot_grid$Data)), " / ", length(prcptot_grid$Data))
message("SDII: ", sum(is.na(sdii_grid$Data)), " / ", length(sdii_grid$Data))
message("DR: ", sum(is.na(dr_grid$Data)), " / ", length(dr_grid$Data))
message("GSL: ", sum(is.na(gsl_grid$Data)), " / ", length(gsl_grid$Data))
message("DL_AGSN: ", sum(is.na(dl_agsn_grid$Data)), " / ", length(dl_agsn_grid$Data))

# Check output grid structure compliance
message("\nGrid structure validation:")
message("NRD has Variable metadata: ", !is.null(nrd_grid$Variable))
message("NRD has Data array: ", !is.null(nrd_grid$Data))
message("NRD has xyCoords: ", !is.null(nrd_grid$xyCoords))
message("NRD has Dates: ", !is.null(nrd_grid$Dates))
message("NRD varName: ", nrd_grid$Variable$varName)
message("NRD units: ", attr(nrd_grid$Variable, "units"))
message("NRD description: ", attr(nrd_grid$Variable, "description"))

# Save results for later inspection
message("\n\nSaving results...")
saveRDS(nrd_grid, "test_results_nrd.rds")
saveRDS(prcptot_grid, "test_results_prcptot.rds")
saveRDS(sdii_grid, "test_results_sdii.rds")
saveRDS(gsl_grid, "test_results_gsl.rds")
saveRDS(dl_agsn_grid, "test_results_dl_agsn.rds")

message("\n[", Sys.time(), "] TEST COMPLETED SUCCESSFULLY!")
message("\nResults saved to: ", getwd())

