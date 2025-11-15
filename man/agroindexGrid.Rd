% Generated manually; do not regenerate with roxygen2
\name{agroindexGrid}
\alias{agroindexGrid}
\title{Agroclimatic Indices in Climate4R}
\description{
Compute agroclimatic indices (FAO tier1, FAO agronomic season, and CDI/CEI stress indices)
from daily climate4R grid objects, returning a grid aggregated to the requested temporal
resolution.
}
\usage{
agroindexGrid(
  index.code,
  tn = NULL,
  tx = NULL,
  pr = NULL,
  tm = NULL,
  hurs = NULL,
  sfcwind = NULL,
  ssrd = NULL,
  index.arg.list = list(),
  cal = "365_day",
  time.resolution = "year",
  parallel = FALSE,
  max.ncores = 16,
  ncores = NULL,
  verbose = FALSE,
  ._retry = FALSE
)
}
\arguments{
\item{index.code}{Character string identifying the index to compute. Supported values include
the FAO tier1 family (\code{gsl}, \code{avg}, \code{nd_thre}, \code{nhw}, \code{dr},
\code{prcptot}, \code{nrd}, \code{lds}, \code{sdii}, \code{prcptot_thre}, \code{ns}),
the FAO agronomic-season set (\code{dt_st_rnagsn}, \code{nm_flst_rnagsn}, \code{dt_fnst_rnagsn}, \code{dt_ed_rnagsn},
\code{dl_agsn}, \code{dc_agsn}, \code{rn_agsn}, \code{avrn_agsn}, \code{dc_rnlg_agsn}, \code{tm_agsn},
\code{dc_txh_agsn}, \code{dc_tnh_agsn}), and the stress indices \code{CDI} and \code{CEI}.}

\item{tn, tx, pr, tm, hurs, sfcwind, ssrd}{Optional daily climate4R grids providing minimum and maximum
temperature, precipitation, mean temperature, relative humidity, wind speed, and surface solar radiation,
respectively. Supply the variables required by the selected index.}

\item{index.arg.list}{Named list of additional arguments passed to the specific index function.
Common entries include \code{lat}, \code{threshold}, \code{duration}, \code{direction}, \code{wet.threshold},
\code{spell.length}, \code{spell.type}, \code{year.start}, \code{year.end}, and \code{pnan}.
See the help page of the chosen index for the full list.}

\item{cal}{Character string describing the calendar type to assume when working with FAO agronomic
indices. Defaults to \code{"365_day"}.}

\item{time.resolution}{Temporal aggregation for the output grid; one of \code{"month"},
\code{"year"} (default), or \code{"climatology"}. FAO indices are always computed yearly.}

\item{parallel}{Logical flag enabling member-wise parallel execution when multiple ensemble members are present.}

\item{max.ncores}{Upper bound on the number of cores used when \code{parallel = TRUE}.}

\item{ncores}{Exact number of cores to use for parallel loops (overrides the automatic detection
when supplied).}

\item{verbose}{Logical flag that prints per-member progress and debugging output, including detailed
status messages for CDI/CEI parameter lists. Defaults to \code{FALSE}.}

\item{._retry}{Internal logical guard that avoids infinite recursion when the function retries
sequential execution to mitigate parallel connection errors. Users should not set this argument
directly.}
}
\details{
The helper \code{\link{agroindexShow}} lists all available indices and their requirements.
When multiple input grids are provided, the function checks temporal and spatial compatibility,
interpolates where needed, and aggregates to the requested resolution. Stress indices (CDI/CEI)
operate on generic bounds supplied through \code{index.arg.list}; FAO agronomic indices fall back
to the \code{agroindexFAO} implementation.
}
\section{Parallel Processing}{
Member-level parallelisation uses the \pkg{parallel} package. The default strategy relies on
socket clusters; on Windows, only single-member execution is available. Use \code{max.ncores}
and \code{ncores} to control resource usage.
}
\examples{
\dontrun{
require(transformeR)
require(visualizeR)

# Load temperature and precipitation data
data("tasmax.grid")
data("tasmin.grid")
data("pr.grid")

## Example 1: Growing Season Length (GSL)
gsl.grid <- agroindexGrid(tm = tasmax.grid, 
                          index.code = "gsl",
                          time.resolution = "year",
                          index.arg.list = list(lat = 40))

## Example 2: Average temperature with custom season
avg.grid <- agroindexGrid(tm = tasmax.grid,
                          index.code = "avg",
                          time.resolution = "year",
                          index.arg.list = list(
                            year.start = "2000-06-01",
                            year.end = "2000-08-31"))

## Example 3: Number of heat waves
nhw.grid <- agroindexGrid(tx = tasmax.grid,
                          index.code = "nhw",
                          time.resolution = "year",
                          index.arg.list = list(
                            threshold = 35,
                            duration = 3))

## Example 4: CDI - Condition Duration Index (multi-variable stress)
# Requires bounds specification for each variable
cdi.grid <- agroindexGrid(tx = tasmax.grid,
                          hurs = hurs.grid,
                          index.code = "CDI",
                          time.resolution = "year",
                          index.arg.list = list(
                            bounds = data.frame(
                              var = c("tx", "hurs"),
                              lower = c(30, 0),
                              upper = c(Inf, 40)),
                            combiner = "all",
                            min_duration = 3))

}
}
\author{
Pablo Lavin Pellon
}

