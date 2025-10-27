# Test Suite for agroindexGrid.R
# Comprehensive testing of all FAO indices, CDI, CEI and error handling

# Load required libraries
library(transformeR)
library(visualizeR)
library(dplyr)
library(magrittr)
library(abind)

# Source the functions to test
source("R/agroindexGrid.R")
source("R/indicesFAO_tier1.R")
source("R/indicesFAO.R")
source("R/cdi_cei_functions.R")

# =============================================================================
# 1. TEST DATA SETUP
# =============================================================================

cat("=== Setting up test data ===\n")

# Create synthetic climate data for testing
create_test_data <- function(n_days = 365, n_lat = 3, n_lon = 3, start_year = 2020) {
  # Create dates
  dates <- seq(as.Date(paste0(start_year, "-01-01")), 
               by = "day", length.out = n_days)
  
  # Create coordinates
  lats <- seq(40, 42, length.out = n_lat)
  lons <- seq(-3, -1, length.out = n_lon)
  
  # Create synthetic data with realistic patterns
  set.seed(123)
  
  # Temperature data (seasonal pattern)
  base_temp <- 15 + 10 * sin(2 * pi * (1:n_days) / 365)
  tx_data <- array(NA, dim = c(n_days, n_lat, n_lon))
  tn_data <- array(NA, dim = c(n_days, n_lat, n_lon))
  t2m_data <- array(NA, dim = c(n_days, n_lat, n_lon))
  
  for (i in 1:n_lat) {
    for (j in 1:n_lon) {
      # Add spatial variation
      spatial_var <- (i-1) * 0.5 + (j-1) * 0.3
      tx_data[, i, j] <- base_temp + spatial_var + rnorm(n_days, 0, 2)
      tn_data[, i, j] <- base_temp + spatial_var - 5 + rnorm(n_days, 0, 1.5)
      t2m_data[, i, j] <- (tx_data[, i, j] + tn_data[, i, j]) / 2
    }
  }
  
  # Precipitation data (more variable, some dry periods)
  pr_data <- array(NA, dim = c(n_days, n_lat, n_lon))
  for (i in 1:n_lat) {
    for (j in 1:n_lon) {
      # Seasonal precipitation pattern
      seasonal_pr <- 2 + 3 * sin(2 * pi * (1:n_days) / 365 + pi/2)
      pr_data[, i, j] <- pmax(0, seasonal_pr + rnorm(n_days, 0, 2))
      # Add some dry periods
      dry_periods <- sample(1:n_days, size = n_days * 0.3)
      pr_data[dry_periods, i, j] <- 0
    }
  }
  
  # Relative humidity (inverse relationship with temperature)
  hurs_data <- array(NA, dim = c(n_days, n_lat, n_lon))
  for (i in 1:n_lat) {
    for (j in 1:n_lon) {
      hurs_data[, i, j] <- pmax(20, pmin(100, 80 - (t2m_data[, i, j] - 15) * 2 + rnorm(n_days, 0, 5)))
    }
  }
  
  # Create climate4R grid objects
  tx_grid <- makeGrid(tx_data, lats, lons, dates, "tx")
  tn_grid <- makeGrid(tn_data, lats, lons, dates, "tn")
  pr_grid <- makeGrid(pr_data, lats, lons, dates, "pr")
  t2m_grid <- makeGrid(t2m_data, lats, lons, dates, "t2m")
  hurs_grid <- makeGrid(hurs_data, lats, lons, dates, "hurs")
  
  return(list(
    tx = tx_grid,
    tn = tn_grid,
    pr = pr_grid,
    t2m = t2m_grid,
    hurs = hurs_grid,
    dates = dates,
    lats = lats,
    lons = lons
  ))
}

# Helper function to create proper climate4R grid objects
makeGrid <- function(data, lats, lons, dates, varName = "test") {
  # Create proper climate4R grid structure
  grid <- list()
  
  # Data array with proper dimensions
  grid$Data <- data
  attr(grid$Data, "dimensions") <- c("time", "lat", "lon")
  
  # Spatial coordinates
  grid$xyCoords <- list(x = lons, y = lats)
  attr(grid$xyCoords, "resX") <- abs(lons[2] - lons[1])
  attr(grid$xyCoords, "projection") <- "LatLonProjection"
  attr(grid$xyCoords, "interpolation") <- "bilinear"
  attr(grid$xyCoords, "subset") <- "subsetSpatial"
  
  # Temporal coordinates
  grid$Dates <- list(start = dates, end = dates)
  
  # Variable metadata
  grid$Variable <- list(varName = varName, level = NULL)
  attr(grid$Variable, "use_dictionary") <- FALSE
  attr(grid$Variable, "description") <- paste("Daily", varName)
  attr(grid$Variable, "units") <- ifelse(varName %in% c("tx", "tn", "t2m"), "degC", 
                                        ifelse(varName == "pr", "mm", 
                                              ifelse(varName == "hurs", "%", "units")))
  attr(grid$Variable, "longname") <- paste("Daily", varName)
  attr(grid$Variable, "daily_agg_cellfun") <- "none"
  attr(grid$Variable, "monthly_agg_cellfun") <- "none"
  attr(grid$Variable, "verification_time") <- "none"
  
  # Additional metadata
  grid$InitializationDates <- NULL
  grid$Members <- NULL
  
  # Dataset attributes
  attr(grid, "dataset") <- "synthetic_test_data"
  attr(grid, "R_package_desc") <- "climate4R.agro-test"
  attr(grid, "R_package_URL") <- "https://github.com/SantanderMetGroup/climate4R.agro"
  attr(grid, "R_package_ref") <- "test_data"
  
  # Set class
  class(grid) <- c("grid", "regular")
  
  return(grid)
}

# Create test data
test_data <- create_test_data(n_days = 730, n_lat = 3, n_lon = 3)  # 2 years of data

cat("Test data created: 730 days, 3x3 grid\n")
cat("Date range:", as.character(min(test_data$dates)), "to", as.character(max(test_data$dates)), "\n")
cat("Latitude range:", min(test_data$lats), "to", max(test_data$lats), "\n")
cat("Longitude range:", min(test_data$lons), "to", max(test_data$lons), "\n\n")

# =============================================================================
# 2. TEST FAO TIER1 INDICES
# =============================================================================

cat("=== Testing FAO Tier1 Indices ===\n")

test_fao_tier1_indices <- function() {
  results <- list()
  
  # Test GSL (Growing Season Length) - requires t2m and lat
  cat("Testing GSL (Growing Season Length)...\n")
  tryCatch({
    gsl_result <- agroindexGrid(
      index.code = "gsl",
      t2m = test_data$t2m,
      index.arg.list = list(lat = 40),
      time.resolution = "year"
    )
    results$gsl <- list(success = TRUE, result = gsl_result)
    cat("  ✓ GSL test passed\n")
  }, error = function(e) {
    results$gsl <- list(success = FALSE, error = e$message)
    cat("  ✗ GSL test failed:", e$message, "\n")
  })
  
  # Test AVG (Average Temperature) - requires t2m
  cat("Testing AVG (Average Temperature)...\n")
  tryCatch({
    avg_result <- agroindexGrid(
      index.code = "avg",
      t2m = test_data$t2m,
      time.resolution = "year"
    )
    results$avg <- list(success = TRUE, result = avg_result)
    cat("  ✓ AVG test passed\n")
  }, error = function(e) {
    results$avg <- list(success = FALSE, error = e$message)
    cat("  ✗ AVG test failed:", e$message, "\n")
  })
  
  # Test ND_THRE (Number of Days Exceeding Threshold) - requires t2m
  cat("Testing ND_THRE (Number of Days Exceeding Threshold)...\n")
  tryCatch({
    nd_thre_result <- agroindexGrid(
      index.code = "nd_thre",
      t2m = test_data$t2m,
      index.arg.list = list(threshold = 25, direction = "geq"),
      time.resolution = "year"
    )
    results$nd_thre <- list(success = TRUE, result = nd_thre_result)
    cat("  ✓ ND_THRE test passed\n")
  }, error = function(e) {
    results$nd_thre <- list(success = FALSE, error = e$message)
    cat("  ✗ ND_THRE test failed:", e$message, "\n")
  })
  
  # Test NHW (Number of Heatwaves) - requires tx
  cat("Testing NHW (Number of Heatwaves)...\n")
  tryCatch({
    nhw_result <- agroindexGrid(
      index.code = "nhw",
      tx = test_data$tx,
      index.arg.list = list(threshold = 30, duration = 3),
      time.resolution = "year"
    )
    results$nhw <- list(success = TRUE, result = nhw_result)
    cat("  ✓ NHW test passed\n")
  }, error = function(e) {
    results$nhw <- list(success = FALSE, error = e$message)
    cat("  ✗ NHW test failed:", e$message, "\n")
  })
  
  # Test DR (Diurnal Temperature Range) - requires tx and tn
  cat("Testing DR (Diurnal Temperature Range)...\n")
  tryCatch({
    dr_result <- agroindexGrid(
      index.code = "dr",
      tx = test_data$tx,
      tn = test_data$tn,
      time.resolution = "year"
    )
    results$dr <- list(success = TRUE, result = dr_result)
    cat("  ✓ DR test passed\n")
  }, error = function(e) {
    results$dr <- list(success = FALSE, error = e$message)
    cat("  ✗ DR test failed:", e$message, "\n")
  })
  
  # Test PRCPTOT (Total Precipitation Wet Days) - requires pr
  cat("Testing PRCPTOT (Total Precipitation Wet Days)...\n")
  tryCatch({
    prcptot_result <- agroindexGrid(
      index.code = "prcptot",
      pr = test_data$pr,
      time.resolution = "year"
    )
    results$prcptot <- list(success = TRUE, result = prcptot_result)
    cat("  ✓ PRCPTOT test passed\n")
  }, error = function(e) {
    results$prcptot <- list(success = FALSE, error = e$message)
    cat("  ✗ PRCPTOT test failed:", e$message, "\n")
  })
  
  # Test NRD (Number of Rainy Days) - requires pr
  cat("Testing NRD (Number of Rainy Days)...\n")
  tryCatch({
    nrd_result <- agroindexGrid(
      index.code = "nrd",
      pr = test_data$pr,
      time.resolution = "year"
    )
    results$nrd <- list(success = TRUE, result = nrd_result)
    cat("  ✓ NRD test passed\n")
  }, error = function(e) {
    results$nrd <- list(success = FALSE, error = e$message)
    cat("  ✗ NRD test failed:", e$message, "\n")
  })
  
  # Test LDS (Length of Dry Spell) - requires pr
  cat("Testing LDS (Length of Dry Spell)...\n")
  tryCatch({
    lds_result <- agroindexGrid(
      index.code = "lds",
      pr = test_data$pr,
      index.arg.list = list(wet.threshold = 1, length.spell = "mean"),
      time.resolution = "year"
    )
    results$lds <- list(success = TRUE, result = lds_result)
    cat("  ✓ LDS test passed\n")
  }, error = function(e) {
    results$lds <- list(success = FALSE, error = e$message)
    cat("  ✗ LDS test failed:", e$message, "\n")
  })
  
  # Test SDII (Simple Daily Intensity Index) - requires pr
  cat("Testing SDII (Simple Daily Intensity Index)...\n")
  tryCatch({
    sdii_result <- agroindexGrid(
      index.code = "sdii",
      pr = test_data$pr,
      time.resolution = "year"
    )
    results$sdii <- list(success = TRUE, result = sdii_result)
    cat("  ✓ SDII test passed\n")
  }, error = function(e) {
    results$sdii <- list(success = FALSE, error = e$message)
    cat("  ✗ SDII test failed:", e$message, "\n")
  })
  
  # Test PRCPTOT_THRE (Total Precipitation Above Threshold) - requires pr
  cat("Testing PRCPTOT_THRE (Total Precipitation Above Threshold)...\n")
  tryCatch({
    prcptot_thre_result <- agroindexGrid(
      index.code = "prcptot_thre",
      pr = test_data$pr,
      index.arg.list = list(threshold = 10),
      time.resolution = "year"
    )
    results$prcptot_thre <- list(success = TRUE, result = prcptot_thre_result)
    cat("  ✓ PRCPTOT_THRE test passed\n")
  }, error = function(e) {
    results$prcptot_thre <- list(success = FALSE, error = e$message)
    cat("  ✗ PRCPTOT_THRE test failed:", e$message, "\n")
  })
  
  # Test NS (Number of Spells) - requires pr
  cat("Testing NS (Number of Spells)...\n")
  tryCatch({
    ns_result <- agroindexGrid(
      index.code = "ns",
      pr = test_data$pr,
      index.arg.list = list(wet.threshold = 1, duration = 3, type.spell = "dry"),
      time.resolution = "year"
    )
    results$ns <- list(success = TRUE, result = ns_result)
    cat("  ✓ NS test passed\n")
  }, error = function(e) {
    results$ns <- list(success = FALSE, error = e$message)
    cat("  ✗ NS test failed:", e$message, "\n")
  })
  
  return(results)
}

# Run FAO tier1 tests
fao_tier1_results <- test_fao_tier1_indices()

# =============================================================================
# 3. TEST FAO AGRONOMIC SEASON INDICES
# =============================================================================

cat("\n=== Testing FAO Agronomic Season Indices ===\n")

test_fao_agronomic_indices <- function() {
  results <- list()
  
  # All agronomic indices require tn, tx, pr and lat
  agronomic_indices <- c(
    "dt_st_rnagsn", "nm_flst_rnagsn", "dt_fnst_rnagsn", "dt_ed_rnagsn",
    "dl_agsn", "dc_agsn", "rn_agsn", "avrn_agsn", "dc_rnlg_agsn",
    "tm_agsn", "dc_txh_agsn", "dc_tnh_agsn"
  )
  
  for (idx in agronomic_indices) {
    cat("Testing", idx, "...\n")
    tryCatch({
      result <- agroindexGrid(
        index.code = idx,
        tn = test_data$tn,
        tx = test_data$tx,
        pr = test_data$pr,
        index.arg.list = list(lat = 40),
        time.resolution = "year"
      )
      results[[idx]] <- list(success = TRUE, result = result)
      cat("  ✓", idx, "test passed\n")
    }, error = function(e) {
      results[[idx]] <- list(success = FALSE, error = e$message)
      cat("  ✗", idx, "test failed:", e$message, "\n")
    })
  }
  
  return(results)
}

# Run FAO agronomic tests
fao_agronomic_results <- test_fao_agronomic_indices()

# =============================================================================
# 4. TEST CDI AND CEI STRESS INDICES
# =============================================================================

cat("\n=== Testing CDI and CEI Stress Indices ===\n")

test_cdi_cei_indices <- function() {
  results <- list()
  
  # Test CDI (Condition Duration Index) - requires tx, t2m, hurs
  cat("Testing CDI (Condition Duration Index)...\n")
  tryCatch({
    cdi_result <- agroindexGrid(
      index.code = "CDI",
      tx = test_data$tx,
      t2m = test_data$t2m,
      hurs = test_data$hurs,
      index.arg.list = list(
        bounds = data.frame(
          var = c("tx", "hurs"),
          lower = c(30, 0),
          upper = c(Inf, 40)
        ),
        combiner = "all",
        min_duration = 3
      ),
      time.resolution = "year"
    )
    results$CDI <- list(success = TRUE, result = cdi_result)
    cat("  ✓ CDI test passed\n")
  }, error = function(e) {
    results$CDI <- list(success = FALSE, error = e$message)
    cat("  ✗ CDI test failed:", e$message, "\n")
  })
  
  # Test CEI (Condition Excess Index) - requires tx, t2m, hurs
  cat("Testing CEI (Condition Excess Index)...\n")
  tryCatch({
    cei_result <- agroindexGrid(
      index.code = "CEI",
      tx = test_data$tx,
      t2m = test_data$t2m,
      hurs = test_data$hurs,
      index.arg.list = list(
        x_col = "tx",
        lower = 25,
        upper = Inf,
        min_duration = 2
      ),
      time.resolution = "year"
    )
    results$CEI <- list(success = TRUE, result = cei_result)
    cat("  ✓ CEI test passed\n")
  }, error = function(e) {
    results$CEI <- list(success = FALSE, error = e$message)
    cat("  ✗ CEI test failed:", e$message, "\n")
  })
  
  # Test CDI with different combiner options
  cat("Testing CDI with 'any' combiner...\n")
  tryCatch({
    cdi_any_result <- agroindexGrid(
      index.code = "CDI",
      tx = test_data$tx,
      t2m = test_data$t2m,
      hurs = test_data$hurs,
      index.arg.list = list(
        bounds = data.frame(
          var = c("tx", "hurs"),
          lower = c(25, 30),
          upper = c(35, 80)
        ),
        combiner = "any",
        min_duration = 2
      ),
      time.resolution = "year"
    )
    results$CDI_any <- list(success = TRUE, result = cdi_any_result)
    cat("  ✓ CDI 'any' combiner test passed\n")
  }, error = function(e) {
    results$CDI_any <- list(success = FALSE, error = e$message)
    cat("  ✗ CDI 'any' combiner test failed:", e$message, "\n")
  })
  
  # Test CDI with k_of_n combiner
  cat("Testing CDI with 'k_of_n' combiner...\n")
  tryCatch({
    cdi_kofn_result <- agroindexGrid(
      index.code = "CDI",
      tx = test_data$tx,
      t2m = test_data$t2m,
      hurs = test_data$hurs,
      index.arg.list = list(
        bounds = data.frame(
          var = c("tx", "t2m", "hurs"),
          lower = c(20, 15, 40),
          upper = c(40, 30, 80)
        ),
        combiner = "k_of_n",
        k = 2,
        min_duration = 3
      ),
      time.resolution = "year"
    )
    results$CDI_kofn <- list(success = TRUE, result = cdi_kofn_result)
    cat("  ✓ CDI 'k_of_n' combiner test passed\n")
  }, error = function(e) {
    results$CDI_kofn <- list(success = FALSE, error = e$message)
    cat("  ✗ CDI 'k_of_n' combiner test failed:", e$message, "\n")
  })
  
  return(results)
}

# Run CDI/CEI tests
cdi_cei_results <- test_cdi_cei_indices()

# =============================================================================
# 5. TEST ERROR HANDLING AND EDGE CASES
# =============================================================================

cat("\n=== Testing Error Handling and Edge Cases ===\n")

test_error_handling <- function() {
  results <- list()
  
  # Test 1: Missing required variables
  cat("Testing missing required variables...\n")
  tryCatch({
    agroindexGrid(index.code = "gsl", time.resolution = "year")
    results$missing_vars <- list(success = FALSE, error = "Should have failed but didn't")
    cat("  ✗ Missing variables test failed - should have thrown error\n")
  }, error = function(e) {
    results$missing_vars <- list(success = TRUE, error = e$message)
    cat("  ✓ Missing variables test passed - correctly threw error\n")
  })
  
  # Test 2: Invalid index code
  cat("Testing invalid index code...\n")
  tryCatch({
    agroindexGrid(index.code = "invalid_index", t2m = test_data$t2m, time.resolution = "year")
    results$invalid_index <- list(success = FALSE, error = "Should have failed but didn't")
    cat("  ✗ Invalid index test failed - should have thrown error\n")
  }, error = function(e) {
    results$invalid_index <- list(success = TRUE, error = e$message)
    cat("  ✓ Invalid index test passed - correctly threw error\n")
  })
  
  # Test 3: Invalid time resolution
  cat("Testing invalid time resolution...\n")
  tryCatch({
    agroindexGrid(index.code = "avg", t2m = test_data$t2m, time.resolution = "invalid")
    results$invalid_time_res <- list(success = FALSE, error = "Should have failed but didn't")
    cat("  ✗ Invalid time resolution test failed - should have thrown error\n")
  }, error = function(e) {
    results$invalid_time_res <- list(success = TRUE, error = e$message)
    cat("  ✓ Invalid time resolution test passed - correctly threw error\n")
  })
  
  # Test 4: Non-daily data (create monthly data)
  cat("Testing non-daily data...\n")
  monthly_data <- test_data$t2m
  # Modify dates to be monthly
  monthly_dates <- seq(as.Date("2020-01-01"), by = "month", length.out = 24)
  monthly_data$Dates$start <- monthly_dates
  monthly_data$Dates$end <- monthly_dates
  
  tryCatch({
    agroindexGrid(index.code = "avg", t2m = monthly_data, time.resolution = "year")
    results$non_daily <- list(success = FALSE, error = "Should have failed but didn't")
    cat("  ✗ Non-daily data test failed - should have thrown error\n")
  }, error = function(e) {
    results$non_daily <- list(success = TRUE, error = e$message)
    cat("  ✓ Non-daily data test passed - correctly threw error\n")
  })
  
  # Test 5: CDI with missing bounds
  cat("Testing CDI with missing bounds...\n")
  tryCatch({
    agroindexGrid(
      index.code = "CDI",
      tx = test_data$tx,
      t2m = test_data$t2m,
      hurs = test_data$hurs,
      time.resolution = "year"
    )
    results$cdi_missing_bounds <- list(success = FALSE, error = "Should have failed but didn't")
    cat("  ✗ CDI missing bounds test failed - should have thrown error\n")
  }, error = function(e) {
    results$cdi_missing_bounds <- list(success = TRUE, error = e$message)
    cat("  ✓ CDI missing bounds test passed - correctly threw error\n")
  })
  
  # Test 6: CEI with missing x_col
  cat("Testing CEI with missing x_col...\n")
  tryCatch({
    agroindexGrid(
      index.code = "CEI",
      tx = test_data$tx,
      t2m = test_data$t2m,
      hurs = test_data$hurs,
      index.arg.list = list(lower = 25),
      time.resolution = "year"
    )
    results$cei_missing_xcol <- list(success = FALSE, error = "Should have failed but didn't")
    cat("  ✗ CEI missing x_col test failed - should have thrown error\n")
  }, error = function(e) {
    results$cei_missing_xcol <- list(success = TRUE, error = e$message)
    cat("  ✓ CEI missing x_col test passed - correctly threw error\n")
  })
  
  # Test 7: Grid consistency (different spatial grids)
  cat("Testing grid consistency...\n")
  # Create grids with different spatial dimensions
  small_grid <- makeGrid(
    array(rnorm(730 * 2 * 2), dim = c(730, 2, 2)),
    c(40, 41), c(-3, -2),
    test_data$dates
  )
  
  tryCatch({
    agroindexGrid(
      index.code = "dr",
      tx = test_data$tx,  # 3x3 grid
      tn = small_grid,    # 2x2 grid
      time.resolution = "year"
    )
    results$grid_consistency <- list(success = FALSE, error = "Should have failed but didn't")
    cat("  ✗ Grid consistency test failed - should have thrown error\n")
  }, error = function(e) {
    results$grid_consistency <- list(success = TRUE, error = e$message)
    cat("  ✓ Grid consistency test passed - correctly threw error\n")
  })
  
  return(results)
}

# Run error handling tests
error_results <- test_error_handling()

# =============================================================================
# 6. TEST INPUT VALIDATION AND CONSISTENCY CHECKS
# =============================================================================

cat("\n=== Testing Input Validation and Consistency Checks ===\n")

test_validation <- function() {
  results <- list()
  
  # Test 1: Parallel processing
  cat("Testing parallel processing...\n")
  tryCatch({
    parallel_result <- agroindexGrid(
      index.code = "avg",
      t2m = test_data$t2m,
      time.resolution = "year",
      parallel = TRUE,
      max.ncores = 2
    )
    results$parallel <- list(success = TRUE, result = parallel_result)
    cat("  ✓ Parallel processing test passed\n")
  }, error = function(e) {
    results$parallel <- list(success = FALSE, error = e$message)
    cat("  ✗ Parallel processing test failed:", e$message, "\n")
  })
  
  # Test 2: Different time resolutions
  cat("Testing different time resolutions...\n")
  time_resolutions <- c("month", "year", "climatology")
  
  for (tr in time_resolutions) {
    cat("  Testing time resolution:", tr, "\n")
    tryCatch({
      tr_result <- agroindexGrid(
        index.code = "avg",
        t2m = test_data$t2m,
        time.resolution = tr
      )
      results[[paste0("time_res_", tr)]] <- list(success = TRUE, result = tr_result)
      cat("    ✓", tr, "test passed\n")
    }, error = function(e) {
      results[[paste0("time_res_", tr)]] <- list(success = FALSE, error = e$message)
      cat("    ✗", tr, "test failed:", e$message, "\n")
    })
  }
  
  # Test 3: Station data (simulate by creating single point)
  cat("Testing station data...\n")
  station_data <- makeGrid(
    array(rnorm(730), dim = c(730, 1, 1)),
    c(40), c(-3),
    test_data$dates
  )
  class(station_data) <- c("station", "grid")
  
  tryCatch({
    station_result <- agroindexGrid(
      index.code = "avg",
      t2m = station_data,
      time.resolution = "year"
    )
    results$station <- list(success = TRUE, result = station_result)
    cat("  ✓ Station data test passed\n")
  }, error = function(e) {
    results$station <- list(success = FALSE, error = e$message)
    cat("  ✗ Station data test failed:", e$message, "\n")
  })
  
  # Test 4: Data with NAs
  cat("Testing data with NAs...\n")
  na_data <- test_data$t2m
  # Introduce some NAs
  na_data$Data[1:10, 1, 1] <- NA
  
  tryCatch({
    na_result <- agroindexGrid(
      index.code = "avg",
      t2m = na_data,
      time.resolution = "year"
    )
    results$na_data <- list(success = TRUE, result = na_result)
    cat("  ✓ NA data test passed\n")
  }, error = function(e) {
    results$na_data <- list(success = FALSE, error = e$message)
    cat("  ✗ NA data test failed:", e$message, "\n")
  })
  
  return(results)
}

# Run validation tests
validation_results <- test_validation()

# =============================================================================
# 7. SUMMARY AND RESULTS
# =============================================================================

cat("\n=== TEST SUMMARY ===\n")

# Count successes and failures
count_results <- function(results) {
  successes <- sum(sapply(results, function(x) if(is.list(x) && "success" %in% names(x)) x$success else FALSE))
  failures <- sum(sapply(results, function(x) if(is.list(x) && "success" %in% names(x)) !x$success else TRUE))
  return(list(successes = successes, failures = failures))
}

# FAO Tier1 summary
fao_tier1_summary <- count_results(fao_tier1_results)
cat("FAO Tier1 Indices:", fao_tier1_summary$successes, "passed,", fao_tier1_summary$failures, "failed\n")

# FAO Agronomic summary
fao_agronomic_summary <- count_results(fao_agronomic_results)
cat("FAO Agronomic Indices:", fao_agronomic_summary$successes, "passed,", fao_agronomic_summary$failures, "failed\n")

# CDI/CEI summary
cdi_cei_summary <- count_results(cdi_cei_results)
cat("CDI/CEI Indices:", cdi_cei_summary$successes, "passed,", cdi_cei_summary$failures, "failed\n")

# Error handling summary
error_summary <- count_results(error_results)
cat("Error Handling Tests:", error_summary$successes, "passed,", error_summary$failures, "failed\n")

# Validation summary
validation_summary <- count_results(validation_results)
cat("Validation Tests:", validation_summary$successes, "passed,", validation_summary$failures, "failed\n")

# Overall summary
total_successes <- fao_tier1_summary$successes + fao_agronomic_summary$successes + 
                  cdi_cei_summary$successes + error_summary$successes + validation_summary$successes
total_failures <- fao_tier1_summary$failures + fao_agronomic_summary$failures + 
                  cdi_cei_summary$failures + error_summary$failures + validation_summary$failures

cat("\nOVERALL RESULTS:", total_successes, "passed,", total_failures, "failed\n")

if (total_failures == 0) {
  cat("🎉 ALL TESTS PASSED! 🎉\n")
} else {
  cat("⚠️  Some tests failed. Check the detailed output above.\n")
}

# =============================================================================
# 8. DETAILED ERROR REPORTING
# =============================================================================

if (total_failures > 0) {
  cat("\n=== DETAILED ERROR REPORT ===\n")
  
  # Report FAO Tier1 failures
  if (fao_tier1_summary$failures > 0) {
    cat("\nFAO Tier1 Failures:\n")
    for (name in names(fao_tier1_results)) {
      if (!fao_tier1_results[[name]]$success) {
        cat("  -", name, ":", fao_tier1_results[[name]]$error, "\n")
      }
    }
  }
  
  # Report FAO Agronomic failures
  if (fao_agronomic_summary$failures > 0) {
    cat("\nFAO Agronomic Failures:\n")
    for (name in names(fao_agronomic_results)) {
      if (!fao_agronomic_results[[name]]$success) {
        cat("  -", name, ":", fao_agronomic_results[[name]]$error, "\n")
      }
    }
  }
  
  # Report CDI/CEI failures
  if (cdi_cei_summary$failures > 0) {
    cat("\nCDI/CEI Failures:\n")
    for (name in names(cdi_cei_results)) {
      if (!cdi_cei_results[[name]]$success) {
        cat("  -", name, ":", cdi_cei_results[[name]]$error, "\n")
      }
    }
  }
  
  # Report error handling failures
  if (error_summary$failures > 0) {
    cat("\nError Handling Failures:\n")
    for (name in names(error_results)) {
      if (!error_results[[name]]$success) {
        cat("  -", name, ":", error_results[[name]]$error, "\n")
      }
    }
  }
  
  # Report validation failures
  if (validation_summary$failures > 0) {
    cat("\nValidation Failures:\n")
    for (name in names(validation_results)) {
      if (!validation_results[[name]]$success) {
        cat("  -", name, ":", validation_results[[name]]$error, "\n")
      }
    }
  }
}

cat("\n=== Test completed ===\n")
