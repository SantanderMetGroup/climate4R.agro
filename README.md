# climate4R.agro

A comprehensive R package for calculating agroclimatic indices from climate data, built on the `climate4R` framework.

## Overview

`climate4R.agro` provides a unified interface for calculating various agroclimatic indices from climate data grids. The package supports FAO Tier1 indices, FAO Agronomic Season indices, and stress indices (CDI/CEI) for agricultural climate analysis.

## Features

### FAO Tier1 Indices
- **GSL** - Growing Season Length
- **AVG** - Average Temperature
- **ND_THRE** - Number of Days Exceeding Threshold
- **NHW** - Number of Heatwaves
- **DR** - Diurnal Temperature Range
- **PRCPTOT** - Total Precipitation Wet Days
- **NRD** - Number of Rainy Days
- **LDS** - Length of Dry Spell
- **SDII** - Simple Daily Intensity Index
- **PRCPTOT_THRE** - Total Precipitation Above Threshold
- **NS** - Number of Spells

### FAO Agronomic Season Indices
- **dt_st_rnagsn** - Start date of rainy agronomic season
- **nm_flst_rnagsn** - Number of false starts of rainy agronomic season
- **dt_fnst_rnagsn** - First start date of rainy agronomic season
- **dt_ed_rnagsn** - End date of rainy agronomic season
- **dl_agsn** - Duration of agronomic season
- **dc_agsn** - Duration of crop season
- **rn_agsn** - Rainy season
- **avrn_agsn** - Average rainy season
- **dc_rnlg_agsn** - Duration of rainy season length
- **tm_agsn** - Temperature of agronomic season
- **dc_txh_agsn** - Duration of temperature extremes high
- **dc_tnh_agsn** - Duration of temperature extremes low

### Stress Indices
- **CDI** - Condition Duration Index (multi-variable stress analysis)
- **CEI** - Condition Excess Index (multi-variable stress analysis)

## Installation

```r
# Install from GitHub
devtools::install_github("SantanderMetGroup/climate4R.agro")
```

## Usage

### Basic Usage

```r
library(climate4R.agro)

# Calculate a single index
result <- agroindexGrid(
  index.code = "gsl",
  tx = tx_grid,  # Maximum temperature grid
  tn = tn_grid, # Minimum temperature grid
  time.resolution = "year"
)

# Calculate multiple indices
indices <- c("gsl", "avg", "dr", "prcptot")
results <- lapply(indices, function(idx) {
  agroindexGrid(
    index.code = idx,
    tx = tx_grid,
    tn = tn_grid,
    pr = pr_grid,
    time.resolution = "year"
  )
})
```

### Advanced Usage

```r
# CDI with custom bounds
cdi_result <- agroindexGrid(
  index.code = "CDI",
  tx = tx_grid,
  tn = tn_grid,
  pr = pr_grid,
  bounds = list(tx = c(25, 35), tn = c(10, 20), pr = c(0, 5)),
  combiner = "any"
)

# FAO Agronomic Season indices
agsn_result <- agroindexGrid(
  index.code = "dt_st_rnagsn",
  tx = tx_grid,
  tn = tn_grid,
  pr = pr_grid,
  time.resolution = "year"
)
```

## Input Data Requirements

The package expects `climate4R` grid objects with the following structure:

```r
# Example grid structure
str(grid_object)
# List of 4
#  $ Variable: List of 2
#    ..$ varName: chr "tasmax"
#    ..$ level : NULL
#  $ Data    : num [1:2192, 1:10, 1:16] 15.1 15.3 11.7 14.4 13.1 ...
#  $ xyCoords: List of 2
#    ..$ x: num [1:16] -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 ...
#    ..$ y: num [1:10] 35 36 37 38 39 40 41 42 43 44 ...
#  $ Dates   : List of 2
#    ..$ start: Date[1:2192], format: "2015-01-01" "2015-01-02" ...
#    ..$ end  : Date[1:2192], format: "2015-01-01" "2015-01-02" ...
```

### Required Variables by Index Type

- **Temperature indices**: `tx` (max temp), `tn` (min temp), `t2m` (mean temp)
- **Precipitation indices**: `pr` (precipitation)
- **Multi-variable indices**: Combination of temperature and precipitation variables

## Dependencies

- `transformeR` - Grid manipulation and processing
- `visualizeR` - Visualization utilities
- `dplyr` - Data manipulation
- `magrittr` - Pipe operator
- `abind` - Array binding

## Recent Updates

### Version 0.1
- **Comprehensive test suite**: Added complete testing coverage for all indices
- **Bug fixes**: Fixed critical issues in FAO Tier1 index calculations
- **Improved robustness**: Enhanced error handling and edge case management
- **Grid compatibility**: Fixed climate4R grid object structure compatibility
- **Performance improvements**: Optimized parallel processing and grid operations

## Citation

If you use this package in your research, please cite:

```r
citation("climate4R.agro")
```

## License

This package is part of the climate4R framework developed by the Santander Meteorology Group.

## Contact

- **Repository**: https://github.com/SantanderMetGroup/climate4R.agro
- **Documentation**: See package vignettes and help files
- **Issues**: Report bugs and feature requests on GitHub
