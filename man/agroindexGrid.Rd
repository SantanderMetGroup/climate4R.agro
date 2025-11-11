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
  ncores = NULL
)
}
\arguments{
\item{index.code}{Character string identifying the index to compute. Supported values include
the FAO tier1 family (gsl, avg, nd\_thre, nhw, dr, prcptot, nrd, lds, sdii, prcptot\_thre, ns),
the FAO agronomic-season set (dt\_st\_rnagsn, nm\_flst\_rnagsn, dt\_fnst\_rnagsn, dt\_ed\_rnagsn,
dl\_agsn, dc\_agsn, rn\_agsn, avrn\_agsn, dc\_rnlg\_agsn, tm\_agsn, dc\_txh\_agsn, dc\_tnh\_agsn),
and the stress indices CDI and CEI.}

\item{tn, tx, pr, tm, hurs, sfcwind, ssrd}{Optional daily climate4R grids providing minimum and maximum
temperature, precipitation, mean temperature, relative humidity, wind speed, and surface solar radiation,
respectively. Supply the variables required by the selected index.}

\item{index.arg.list}{Named list of additional arguments passed to the specific index function.
Common entries include \code{lat}, \code{threshold}, \code{duration}, \code{direction}, \code{wet.threshold},
\code{spell.length}, \code{spell.type}, \code{year.start}, \code{year.end}, and \code{pnan}.
See the help page of the chosen index for the full list.}

\item{cal}{Calendar identifier (defaults to \code{"365_day"}).}

\item{time.resolution}{Temporal aggregation for the output grid; one of \code{"month"},
\code{"year"} (default), or \code{"climatology"}. FAO indices are always computed yearly.}

\item{parallel}{Logical flag enabling member-wise parallel execution when multiple ensemble members are present.}

\item{max.ncores}{Upper bound on the number of cores used when \code{parallel = TRUE}.}

\item{ncores}{Exact number of cores to use for parallel loops (overrides the automatic detection
when supplied).}
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
data("tasmax.grid", package = "climate4R.base")
data("tasmin.grid", package = "climate4R.base")

gsl.grid <- agroindexGrid(
  tm = tasmax.grid,
  index.code = "gsl",
  time.resolution = "year",
  index.arg.list = list(lat = 40)
)
}
}
\author{
Pablo Lavin Pellon
}

