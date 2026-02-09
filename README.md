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

The recommended procedure for installing the package is using the devtools package.

```r
# Install from GitHub (also available )
devtools::install_github("SantanderMetGroup/climate4R.agro")
```

If you are experiencing issues when using `devtools` or `remotes`, we recommend manual installation using the source file:

1. Download the `climate4R.agro_1.0.0.tar.gz` file from (https://github.com/SantanderMetGroup/climate4R.agro/releases).

2. Install it from the R console:

```R
install.packages("climate4R.agro_1.0.0.tar.gz", repos = NULL, type = "source")
```

A list of all available indices and the atomic functions calculating them is printed on screen with:

```r
library(climate4R.agro)
agroindexShow()
?agroindexGrid   # see the examples 
```

## Usage

```r
library(climate4R.agro)

# Calculate GSL (Growing Season Length)
gsl_result <- agroindexGrid(
  index.code = "gsl",
  tm = tmean_grid,  # Mean temperature
)

# Calculate Average Temperature with custom period
avg_result <- agroindexGrid(
  index.code = "avg",
  tm = tmean_grid,
)

# Calculate Number of Heatwaves
nhw_result <- agroindexGrid(
  index.code = "nhw",
  tx = tx_grid,  # Maximum temperature
  index.arg.list = list(
    threshold = 35,  # Temperature threshold in degrees C
    duration = 3     # Minimum duration in days
  )
)

# Calculate multiple indices
indices <- c("gsl", "avg", "dr", "prcptot")
results <- lapply(indices, function(idx) {
  agroindexGrid(
    index.code = idx,
    tx = tx_grid,
    tn = tn_grid,
    pr = pr_grid,
    tm = tmean_grid
    )
  }
)

# CDI with custom bounds
cdi_result <- agroindexGrid(
  index.code = "CDI",
  tx = tx_grid,
  pr = pr_grid,
  index.arg.list = list(
    season_start  = "01-01",
    season_end    = "12-31",
    bounds = data.frame(
      var   = c("tx", "pr"),
      lower = c(25, 0),
      upper = c(35, 10)
      ),
      combiner = "any",
      min_duration = 3
      ),
  parallel = TRUE
)

# Number of Days Exceeding Threshold
nd_thre_result <- agroindexGrid(
  index.code = "nd_thre",
  tm = tmean_grid,
  index.arg.list = list(
    threshold = 30,
    direction = "geq"  # "geq" for >=, "leq" for <=
  )
)

# Length of Dry Spell
lds_result <- agroindexGrid(
  index.code = "lds",
  pr = pr_grid,
  index.arg.list = list(
    wet.threshold = 1,      # Minimum precipitation for wet day (mm)
    spell.length = "max"    # "max", "mean", or numeric value
  )
)

# FAO Agronomic Season indices
agsn_result <- agroindexGrid(
  index.code = "dt_st_rnagsn",
  tx = tx_grid,
  tn = tn_grid,
  pr = pr_grid
)
```

## Input Data Requirements

The package expects `climate4R` grid objects with daily temporal resolution. Grids should have the following structure:

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

- **Temperature indices**: `tx` (max temp), `tn` (min temp), `tm` (mean temp)
- **Precipitation indices**: `pr` (precipitation)
- **Multi-variable indices**: Combination of temperature and precipitation variables

### Index-Specific Arguments (`index.arg.list`)

Many indices require additional parameters passed via `index.arg.list`:

- **lat**: Latitude in degrees (optional)
- **threshold**: Threshold value for indices like `nd_thre`, `nhw`, `prcptot_thre`
- **duration**: Duration in days for indices like `nhw`, `ns`
- **direction**: "geq" or "leq" for `nd_thre` (greater/less than threshold)
- **wet.threshold**: Minimum precipitation for wet day (default: 1 mm)
- **spell.length**: For `lds` and `ns` (e.g., "max", "mean", or numeric)
- **spell.type**: "wet" or "dry" for `ns`
- **year.start**, **year.end**: Date strings for custom periods (e.g., "2000-06-01")
- **pnan**: Maximum percentage of missing data allowed (default: 25%)

For FAO agronomic indices (`agroindexFAO`), additional arguments include: `shc`, `rndy`, `rnlg`, `txh`, `tnh`.

For CDI/CEI indices, arguments include: `bounds` (data.frame), `combiner`, `min_duration`.

See the help files for individual index functions (e.g., `?gsl`, `?avg`) for complete parameter lists.

## Important Notes

The `time.resolution` parameter is ignored for FAO indices, and the output will always be yearly data regardless of the `time.resolution` setting. This also applies to the stress indices (CDI/CEI).

## Dependencies

- `transformeR` - Grid manipulation and processing
- `dplyr` - Data manipulation
- `magrittr` - Pipe operator
- `abind` - Array binding
- `parallel` - Parallel computing utilities
- `tibble` - Modern data frames
- `rlang` - Tools for working with R expressions and environments

## License

This package is part of the climate4R framework developed by the Santander Meteorology Group.

## Contact

- **Repository**: https://github.com/SantanderMetGroup/climate4R.agro
- **Documentation**: See package vignettes and help files
- **Issues**: Report bugs and feature requests on GitHub
