#     agroindexGrid.R Agroclimatic Index in Climate4R
#
#     Copyright (C) 2018 (original) and 2025 (modifications) 
#     Santander Meteorology Group (http://www.meteo.unican.es)
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program. If not, see <http://www.gnu.org/licenses/>.

#' @title Agroclimatic Indices in Climate4R
#' @description Calculation of agroclimatic indices including FAO tier1 indices and CDI/CEI stress indices.
#' The function processes climate4R grid objects for seamless integration.
#' @param tn A climate4R dataset of daily minimum temperature (degrees C)
#' @param tx A climate4R dataset of daily maximum temperature (degrees C)
#' @param pr A climate4R dataset of daily precipitation (mm)
#' @param tm A climate4R dataset of daily mean air temperature (degrees C).
#' This parameter matches FAO naming conventions (tm = temperature mean).
#' For indices that don't require tm specifically, it can be derived from tx and tn if provided.
#' @param hurs A climate4R dataset of daily relative humidity (\%)
#' @param sfcwind A climate4R dataset of daily surface wind speed (m/s)
#' @param sp A climate4R dataset of daily surface pressure (Pa)
#' @param ssrd A climate4R dataset of daily surface solar radiation downwards (J/m^2)
#' @param pvpot A climate4R dataset of daily vapor pressure potential (Pa)
#' @param cal A calendar definition. Default to 365-day calendar (not used by FAO/CDI/CEI indices, kept for compatibility).
#' @param time.resolution Output time resolution. Choices are "month", "year" (default) and "climatology".
#' Note: FAO indices are calculated year by year by definition and will ignore this parameter.
#' @param index.code Character string, indicating the specific code of the agroclimatic index. 
#' Options include FAO tier1 indices (gsl, avg, nd_thre, nhw, dr, prcptot, nrd, lds, sdii, prcptot_thre, ns),
#' FAO agronomic season indices (dt_st_rnagsn, nm_flst_rnagsn, dt_fnst_rnagsn, dt_ed_rnagsn, dl_agsn, 
#' dc_agsn, rn_agsn, avrn_agsn, dc_rnlg_agsn, tm_agsn, dc_txh_agsn, dc_tnh_agsn), 
#' and stress indices (CDI, CEI). See Details.
#' @param index.arg.list Optional list of index-specific arguments. Depending on the index, this may be required.
#' Common arguments include:
#' \itemize{
#'   \item \strong{lat}: Latitude in degrees (required for GSL, optional for others)
#'   \item \strong{threshold}: Threshold value for indices like nd_thre, nhw, prcptot_thre
#'   \item \strong{duration}: Duration in days for indices like nhw, ns
#'   \item \strong{direction}: "geq" or "leq" for nd_thre (greater/less than threshold)
#'   \item \strong{wet.threshold}: Minimum precipitation for wet day (default: 1 mm)
#'   \item \strong{spell.length}: For lds and ns (e.g., "max", "mean", or numeric)
#'   \item \strong{spell.type}: "wet" or "dry" for ns
#'   \item \strong{year.start}, \strong{year.end}: Date strings for custom periods (e.g., "2000-06-01")
#'   \item \strong{pnan}: Maximum percentage of missing data allowed (default: 25\%)
#' }
#' For FAO agronomic indices (agroindexFAO), additional arguments include: shc, rndy, rnlg, txh, tnh.
#' For CDI/CEI indices, arguments include: bounds (data.frame), combiner, min_duration.
#' See the help files for individual index functions (e.g., \code{?gsl}, \code{?avg}) for complete parameter lists.
#' @template templateParallelParams
#' @import transformeR
#' @importFrom parallel stopCluster
#' @importFrom magrittr %>% %<>% extract2
#' @importFrom utils head
#' @importFrom abind abind
#' @details \code{\link{agroindexShow}} will display on screen a full list of available agroclimatic indices and their codes.
#' 
#' \strong{Index Groups}
#' \itemize{
#'   \item \strong{FAO Tier1 indices}: Basic agroclimatic indices (gsl, avg, nd_thre, nhw, dr, prcptot, nrd, lds, sdii, prcptot_thre, ns)
#'   \item \strong{FAO Agronomic Season indices}: Indices based on water balance and agronomic season definition
#'   \item \strong{Stress indices}: CDI (Condition Duration Index) and CEI (Condition Excess Index) for multi-variable stress analysis
#' }
#' 
#' For index-specific arguments, refer to the help files of individual index functions.
#'
#' @template templateParallel
#'
#' @examples \dontrun{
#' require(transformeR)
#' require(visualizeR)
#' # Load temperature and precipitation data
#' data("tasmax.grid")
#' data("tasmin.grid")
#' data("pr.grid")
#' 
#' ## Example 1: Growing Season Length (GSL)
#' gsl.grid <- agroindexGrid(tm = tasmax.grid, 
#'                           index.code = "gsl",
#'                           time.resolution = "year",
#'                           index.arg.list = list(lat = 40))
#' 
#' ## Example 2: Average temperature with custom season
#' avg.grid <- agroindexGrid(tm = tasmax.grid,
#'                           index.code = "avg",
#'                           time.resolution = "year",
#'                           index.arg.list = list(
#'                             year.start = "2000-06-01",
#'                             year.end = "2000-08-31"))
#' 
#' ## Example 3: Number of heat waves
#' nhw.grid <- agroindexGrid(tx = tasmax.grid,
#'                           index.code = "nhw",
#'                           time.resolution = "year",
#'                           index.arg.list = list(
#'                             threshold = 35,
#'                             duration = 3))
#' 
#' ## Example 4: CDI - Condition Duration Index (multi-variable stress)
#' # Requires bounds specification for each variable
#' cdi.grid <- agroindexGrid(tx = tasmax.grid,
#'                           hurs = hurs.grid,
#'                           index.code = "CDI",
#'                           time.resolution = "year",
#'                           index.arg.list = list(
#'                             bounds = data.frame(
#'                               var = c("tx", "hurs"),
#'                               lower = c(30, 0),
#'                               upper = c(Inf, 40)),
#'                             combiner = "all",
#'                             min_duration = 3))
#' }
#'
#' @author J. Bedia is the original author of climdexGrid.R 
#' @author Pablo Lavin is the author of the modifications in main function agroindexGrid
#' @export

agroindexGrid <- function(index.code,
                        tn = NULL,
                        tx = NULL,
                        pr = NULL,
                        tm = NULL,
                        hurs = NULL,
                        sfcwind = NULL,
                        sp = NULL,
                        ssrd = NULL,
                        pvpot = NULL,
                        index.arg.list = list(),
                        cal = "365_day",
                        time.resolution = "year",
                        parallel = FALSE,
                        max.ncores = 16,
                        ncores = NULL) {
    if (!is.null(tn)) {
        if (getTimeResolution(tn) != "DD") stop("Daily data is required as input", call. = FALSE)
    }
    if (!is.null(tx)) {
        if (getTimeResolution(tx) != "DD") stop("Daily data is required as input", call. = FALSE)
    }
    if (!is.null(pr)) {
        if (getTimeResolution(pr) != "DD") stop("Daily data is required as input", call. = FALSE)
    }
    if (!is.null(tm)) {
        if (getTimeResolution(tm) != "DD") stop("Daily data is required as input", call. = FALSE)
    }
    if (!is.null(hurs)) {
        if (getTimeResolution(hurs) != "DD") stop("Daily data is required as input", call. = FALSE)
    }
    if (!is.null(sfcwind)) {
        if (getTimeResolution(sfcwind) != "DD") stop("Daily data is required as input", call. = FALSE)
    }
    if (!is.null(sp)) {
        if (getTimeResolution(sp) != "DD") stop("Daily data is required as input", call. = FALSE)
    }
    if (!is.null(ssrd)) {
        if (getTimeResolution(ssrd) != "DD") stop("Daily data is required as input", call. = FALSE)
    }
    if (!is.null(pvpot)) {
        if (getTimeResolution(pvpot) != "DD") stop("Daily data is required as input", call. = FALSE)
    }
    
    # Validate time resolution
    time.resolution <- match.arg(time.resolution, choices = c("month", "year", "climatology"))
    
    index.code <- match.arg(index.code,
                            choices = c("gsl", "avg", "nd_thre", "nhw", "dr", 
                                        "prcptot", "nrd", "lds", "sdii", "prcptot_thre", "ns",
                                        "dt_st_rnagsn", "nm_flst_rnagsn", "dt_fnst_rnagsn", 
                                        "dt_ed_rnagsn", "dl_agsn", "dc_agsn", "rn_agsn", 
                                        "avrn_agsn", "dc_rnlg_agsn", "tm_agsn", 
                                        "dc_txh_agsn", "dc_tnh_agsn", "CDI", "CEI"))
    

    # Enhanced dimension conversion function
    convert_dimension_order <- function(grid, var_name) {
        if (is.null(grid)) return(NULL)
        
        data_dims <- dim(grid$Data)
        dim_names <- attr(grid$Data, "dimensions")
        
        # Only process if we have 4 dimensions
        if (length(data_dims) != 4) return(grid)
        
        # Check if we have the required dimensions
        required_dims <- c("member", "time", "lat", "lon")
        if (!all(required_dims %in% dim_names)) return(grid)
        
        # Get current dimension order
        current_order <- match(required_dims, dim_names)
        
        # Check if member is first (user's format: member x time x lat x lon)
        if (current_order[1] == 1) {
            message("[", Sys.time(), "] Converting ", var_name, " from (member x time x lat x lon) to (time x lat x lon x member)")
            
            # Transpose to climate4R standard: (time x lat x lon x member)
            target_order <- c(2, 3, 4, 1)  # time, lat, lon, member
            grid$Data <- aperm(grid$Data, target_order)
            attr(grid$Data, "dimensions") <- c("time", "lat", "lon", "member")
            
            # Ensure member information is properly set
            if (is.null(grid$Members)) {
                n_members <- data_dims[1]
                grid$Members <- paste0("member_", 1:n_members)
            }
            
            return(grid)
        }
        
        # Check if already in climate4R format (time x lat x lon x member)
        if (all(current_order == c(2, 3, 4, 1))) {
            # Already in correct format, no conversion needed
            return(grid)
        }
        
        # Handle other dimension orders if needed
        if (any(current_order != c(2, 3, 4, 1))) {
            warning("Unusual dimension order detected for ", var_name, 
                    ": ", paste(dim_names, collapse = " x "), 
                    ". Attempting to convert to (time x lat x lon x member)")
            
            # Try to convert to standard order
            target_order <- order(current_order)
            grid$Data <- aperm(grid$Data, target_order)
            attr(grid$Data, "dimensions") <- c("time", "lat", "lon", "member")
        }
        
        return(grid)
    }
    
    # Apply dimension conversion to all input grids
    if (!is.null(tn)) tn <- convert_dimension_order(tn, "tn")
    if (!is.null(tx)) tx <- convert_dimension_order(tx, "tx")
    if (!is.null(pr)) pr <- convert_dimension_order(pr, "pr")
    if (!is.null(tm)) tm <- convert_dimension_order(tm, "tm")
    if (!is.null(hurs)) hurs <- convert_dimension_order(hurs, "hurs")
    if (!is.null(sfcwind)) sfcwind <- convert_dimension_order(sfcwind, "sfcwind")
    if (!is.null(sp)) sp <- convert_dimension_order(sp, "sp")
    if (!is.null(ssrd)) ssrd <- convert_dimension_order(ssrd, "ssrd")
    if (!is.null(pvpot)) pvpot <- convert_dimension_order(pvpot, "pvpot")
    aux <- read.master()
    metadata <- aux[grep(paste0("^", index.code, "$"), aux$code, fixed = FALSE), ]
    a <- c(!is.null(tn), !is.null(tx), !is.null(pr), !is.null(tm), !is.null(hurs), 
           !is.null(sfcwind), !is.null(sp), !is.null(ssrd), !is.null(pvpot)) %>% as.numeric()
    b <- metadata[ , 4:12] %>% as.numeric()
    if (any(b - a > 0)) {
        stop("The required input variable(s) for ", index.code,
             " index calculation are missing\nType \'?",
             metadata$indexfun, "\' for help", call. = FALSE)
    }
    # Remove any possible uneeded input grid
    if (any(a - b > 0)) {
        ind <- which((a - b) > 0)
        rem <- c("tn", "tx", "pr", "tm", "hurs", "sfcwind", "sp", "ssrd", "pvpot")[ind]
        sapply(rem, function(x) assign(x, NULL)) %>% invisible()
        message("NOTE: some input grids provided for ", index.code,
                " index calculation are not required and were removed")
    }
    
    # Grid consistency checks for multi-variable inputs
    grid.list <- list("tn" = tn, "tx" = tx, "pr" = pr, "tm" = tm, 
                     "hurs" = hurs, "sfcwind" = sfcwind, "sp" = sp, 
                     "ssrd" = ssrd, "pvpot" = pvpot)
    grid.list <- grid.list[!sapply(grid.list, is.null)]
    
    if (length(grid.list) > 1) {
        # Check grid consistency
        locs <- lapply(grid.list, isRegular)
        if (!any(sum(unlist(locs)) != 0, sum(unlist(locs)) != length(grid.list))) {
            stop("Regular and Irregular grids can not be combined. See function interpGrid")
        }
        
        # Intersect grids temporally
        original_names <- names(grid.list)
        grid.list <- intersectGrid(grid.list, type = "temporal", which.return = 1:length(grid.list))
        namesgridlist <- original_names
        
        # Check spatial consistency and interpolate if needed
        refgrid <- getGrid(grid.list[[1]])
        indinterp <- which(isFALSE(unlist(lapply(2:length(grid.list), function(i) 
            identical(refgrid, getGrid(grid.list[[i]]))))))
        
        if (length(indinterp) > 0) {
            grid.list.aux <- suppressMessages(lapply(grid.list[indinterp], function(i) 
                interpGrid(i, getGrid(grid.list[[1]]))))
            grid.list[indinterp] <- grid.list.aux
            grid.list.aux <- NULL
        }
        
        # Update variables with processed grids
        for (i in 1:length(grid.list)) {
            assign(namesgridlist[i], grid.list[[i]])
        }
    }
    # Ensure member is present in data structures / handle station <--> grid
    if (!is.null(tx)) {
        station <- ifelse(typeofGrid(tx) == "station", TRUE, FALSE)
        tx %<>% redim(member = TRUE, var = FALSE)
    }
    if (!is.null(tn)) {
        station <- ifelse(typeofGrid(tn) == "station", TRUE, FALSE)
        tn %<>% redim(member = TRUE, var = FALSE)
    }
    if (!is.null(pr)) {
        station <- ifelse(typeofGrid(pr) == "station", TRUE, FALSE)
        pr %<>% redim(member = TRUE, var = FALSE)
    }
    if (!is.null(tm)) {
        station <- ifelse(typeofGrid(tm) == "station", TRUE, FALSE)
        tm %<>% redim(member = TRUE, var = FALSE)
    }
    if (!is.null(hurs)) {
        station <- ifelse(typeofGrid(hurs) == "station", TRUE, FALSE)
        hurs %<>% redim(member = TRUE, var = FALSE)
    }
    if (!is.null(sfcwind)) {
        station <- ifelse(typeofGrid(sfcwind) == "station", TRUE, FALSE)
        sfcwind %<>% redim(member = TRUE, var = FALSE)
    }
    if (!is.null(sp)) {
        station <- ifelse(typeofGrid(sp) == "station", TRUE, FALSE)
        sp %<>% redim(member = TRUE, var = FALSE)
    }
    if (!is.null(ssrd)) {
        station <- ifelse(typeofGrid(ssrd) == "station", TRUE, FALSE)
        ssrd %<>% redim(member = TRUE, var = FALSE)
    }
    if (!is.null(pvpot)) {
        station <- ifelse(typeofGrid(pvpot) == "station", TRUE, FALSE)
        pvpot %<>% redim(member = TRUE, var = FALSE)
    }
    # Check structural consistency of data arrays when multiple
    if (index.code == "DTR" | index.code == "GSL") {
        sapply(list(tn, tx), "checkDim", dimensions = c("member", "time", "lat", "lon")) %>% invisible()
    }
    # name of the reference grid
    refGridName <- c("tn","tx","pr","tm","hurs","sfcwind","sp","ssrd","pvpot")[which(c(!is.null(tn),
                                                                                            !is.null(tx),
                                                                                            !is.null(pr),
                                                                                            !is.null(tm),
                                             !is.null(hurs),
                                             !is.null(sfcwind),
                                             !is.null(sp),
                                             !is.null(ssrd),
                                             !is.null(pvpot)) %>% as.numeric() != 0)] %>% head(1)
    assign("refGrid", get(refGridName))
    # Member dimension consistency check (similar to indexGrid.R)
    grid.list.check <- list("tn" = tn, "tx" = tx, "pr" = pr, "tm" = tm, 
                           "hurs" = hurs, "sfcwind" = sfcwind, "sp" = sp, 
                           "ssrd" = ssrd, "pvpot" = pvpot)
    grid.list.check <- grid.list.check[!sapply(grid.list.check, is.null)]
    if (length(grid.list.check) > 1) {
        ns.mem <- lapply(grid.list.check, function(r) getShape(r)[["member"]])
        # Handle NA (no member dimension) as 1
        ns.mem <- lapply(ns.mem, function(x) if (is.na(x)) 1 else x)
        if (sum(unlist(ns.mem) - rep(ns.mem[[1]], length(ns.mem))) != 0) {
            stop("Number of members is different across input grids")
        }
        n.mem <- unique(unlist(ns.mem))
    } else {
        n.mem <- getShape(refGrid, "member")
        if (is.na(n.mem)) n.mem <- 1  # Handle grids without member dimension
    }
    # Get reference dates (as Date, FAO/CDI/CEI)
    refDates <- getRefDates(refGrid)
    message("[", Sys.time(), "] Calculating ", index.code, " ...")
    if (n.mem > 1) {
        parallel.pars <- parallelCheck(parallel, max.ncores, ncores)
        apply_fun <- selectPar.pplyFun(parallel.pars, .pplyFUN = "lapply")
        if (parallel.pars$hasparallel) on.exit(parallel::stopCluster(parallel.pars$cl))
    } else {
        if (isTRUE(parallel)) message("NOTE: Parallel processing was skipped (unable to parallelize one single member)")
        apply_fun <- lapply
    }
    out.list <- apply_fun(1:n.mem, function(x) {
        tryCatch({
        # if (n.mem > 1) message("[", Sys.time(), "] Calculating ",
                               # index.code, " for member ", x, " ...")
        # Extract data as 3D arrays (time x lat x lon) following agroclimGrid pattern
        aux.tx <- aux.tn <- aux.pr <- aux.tm <- aux.hurs <- NULL
        aux.sfcwind <- aux.sp <- aux.ssrd <- aux.pvpot <- NULL
        if (!is.null(tx)) {
            if (n.mem > 1) {
                tmp <- subsetGrid(tx, members = x, drop = TRUE) %>% redim(loc = FALSE, member = FALSE)
            } else {
                tmp <- tx %>% redim(loc = FALSE, member = FALSE)
            }
            aux.tx <- tmp[["Data"]]
        }
        if (!is.null(tn)) {
            if (n.mem > 1) {
                tmp <- subsetGrid(tn, members = x, drop = TRUE) %>% redim(loc = FALSE, member = FALSE)
            } else {
                tmp <- tn %>% redim(loc = FALSE, member = FALSE)
            }
            aux.tn <- tmp[["Data"]]
        }
        if (!is.null(pr)) {
            if (n.mem > 1) {
                tmp <- subsetGrid(pr, members = x, drop = TRUE) %>% redim(loc = FALSE, member = FALSE)
            } else {
                tmp <- pr %>% redim(loc = FALSE, member = FALSE)
            }
            aux.pr <- tmp[["Data"]]
        }
        if (!is.null(tm)) {
            if (n.mem > 1) {
                tmp <- subsetGrid(tm, members = x, drop = TRUE) %>% redim(loc = FALSE, member = FALSE)
            } else {
                tmp <- tm %>% redim(loc = FALSE, member = FALSE)
            }
            aux.tm <- tmp[["Data"]]
        }
        if (!is.null(hurs)) {
            if (n.mem > 1) {
                tmp <- subsetGrid(hurs, members = x, drop = TRUE) %>% redim(loc = FALSE, member = FALSE)
            } else {
                tmp <- hurs %>% redim(loc = FALSE, member = FALSE)
            }
            aux.hurs <- tmp[["Data"]]
        }
        if (!is.null(sfcwind)) {
            if (n.mem > 1) {
                tmp <- subsetGrid(sfcwind, members = x, drop = TRUE) %>% redim(loc = FALSE, member = FALSE)
            } else {
                tmp <- sfcwind %>% redim(loc = FALSE, member = FALSE)
            }
            aux.sfcwind <- tmp[["Data"]]
        }
        if (!is.null(sp)) {
            if (n.mem > 1) {
                tmp <- subsetGrid(sp, members = x, drop = TRUE) %>% redim(loc = FALSE, member = FALSE)
            } else {
                tmp <- sp %>% redim(loc = FALSE, member = FALSE)
            }
            aux.sp <- tmp[["Data"]]
        }
        if (!is.null(ssrd)) {
            if (n.mem > 1) {
                tmp <- subsetGrid(ssrd, members = x, drop = TRUE) %>% redim(loc = FALSE, member = FALSE)
            } else {
                tmp <- ssrd %>% redim(loc = FALSE, member = FALSE)
            }
            aux.ssrd <- tmp[["Data"]]
        }
        if (!is.null(pvpot)) {
            if (n.mem > 1) {
                tmp <- subsetGrid(pvpot, members = x, drop = TRUE) %>% redim(loc = FALSE, member = FALSE)
            } else {
                tmp <- pvpot %>% redim(loc = FALSE, member = FALSE)
            }
            aux.pvpot <- tmp[["Data"]]
        }
        # EXCEPTION for FAO INDICES (require lat, dates, and NO temporal subsetting)
        # Check for FAO agronomic indices or tier1 indices (by function name or index code)
        if (metadata$indexfun %in% c("agroindexFAO", "agroindexFAO_tier1") || 
            index.code %in% c("gsl", "avg", "nd_thre", "nhw", "dr", "prcptot", "nrd", "lds", "sdii", "prcptot_thre", "ns")) {
            if (time.resolution != "year") {
                message(index.code, " is calculated year by year by definition. argument time.resolution ignored.")
            }
            
            # Build grid.list.aux with processed grids for this member (consistent with indexGrid.R)
            # Always use subsetGrid to extract member x - it handles grids with/without member dimension
            grid.list.aux <- list()
            if (!is.null(tx)) {
                grid.list.aux[["tx"]] <- subsetGrid(tx, members = x) %>% redim(loc = FALSE, member = FALSE)
            }
            if (!is.null(tn)) {
                grid.list.aux[["tn"]] <- subsetGrid(tn, members = x) %>% redim(loc = FALSE, member = FALSE)
            }
            if (!is.null(pr)) {
                grid.list.aux[["pr"]] <- subsetGrid(pr, members = x) %>% redim(loc = FALSE, member = FALSE)
            }
            if (!is.null(tm)) {
                grid.list.aux[["tm"]] <- subsetGrid(tm, members = x) %>% redim(loc = FALSE, member = FALSE)
            }
            
            # Create output grid structure using aggregateGrid on sugrid (sets yearly structure)
            out.aux <- suppressMessages(aggregateGrid(grid.list.aux[[1]], aggr.y = list(FUN = "mean", na.rm = TRUE)))
            
            # Extract data arrays (3D: time x lat x lon) from processed grids
            input.arg.list <- lapply(grid.list.aux, function(d) d[["Data"]])
            
            # Get dates from the processed grid (convert to date matrix format)
            datess <- as.Date(grid.list.aux[[1]][["Dates"]][["start"]])
            datess <- cbind(
                as.numeric(format(datess, "%Y")),
                as.numeric(format(datess, "%m")),
                as.numeric(format(datess, "%d"))
            )
            
            # Get coordinates
            lats <- grid.list.aux[[1]][["xyCoords"]][["y"]]
            n_lons <- getShape(grid.list.aux[[1]])["lon"]
            
            # Prepare arguments
            index.arg.list[["dates"]] <- datess
            # Only add index.code for agroindexFAO, not for tier1 functions
            if (metadata$indexfun == "agroindexFAO") {
                index.arg.list[["index.code"]] <- index.code
            }
            
            # Process lat-lon loops
            latloop <- lapply(1:length(lats), function(l) {
                lonloop <- lapply(1:n_lons, function(lo) {
                    index.arg.list[["lat"]] <- lats[l]
                    # Extract time series for this point from each data array
                    point_data <- lapply(input.arg.list, function(z) z[, l, lo])
                    
                    # Call FAO function with point data and index arguments
                    tryCatch({
                        result <- do.call(metadata$indexfun, c(point_data, index.arg.list))
                        # GSL returns a list, extract the GSL component
                        if (index.code == "gsl" && is.list(result)) {
                            result$GSL
                        } else {
                            result
                        }
                    }, error = function(e) {
                        if (l <= 3 && lo <= 3) {
                            warning("Error at lat=", lats[l], ", lon=", lo, ": ", e$message)
                        }
                        NA
                    })
                })
                do.call(abind::abind, list(lonloop, along = 0))
            })
            
            # Reconstruct output array
            # latloop_bound has dimensions (lat, lon, time) from abind
            # We need (time, lat, lon), so use aperm(..., c(3, 1, 2))
            latloop_bound <- do.call(abind::abind, list(latloop, along = 0))
            out.aux[["Data"]] <- unname(aperm(latloop_bound, c(3, 1, 2)))
            attr(out.aux[["Data"]], "dimensions") <- c("time", "lat", "lon")
            
            # Dates are already correctly set from aggregateGrid (yearly structure)
            tx <- tn <- pr <- tm <- hurs <- sfcwind <- sp <- ssrd <- pvpot <- NULL
            return(out.aux)
        }
        
        # Set dates in index.arg.list if not already set
        if (is.null(index.arg.list[["dates"]])) {
            index.arg.list[["dates"]] <- as.Date(refDates)
        }
        
        # Iterate over latitudes
        latloop <- lapply(1:length(lats), function(l) {
            # Iterate over longitudes
            lonloop <- lapply(1:length(lons), function(lo) {
                # Extract time series for this lat-lon point
                result <- tryCatch({
                    if (index.code %in% c("CDI", "CEI")) {
                        # CDI/CEI require dataframe format
                        df <- data.frame(
                            id = 1,
                            date = as.Date(refDates)
                        )
                        if (!is.null(aux.tx)) df$tx <- aux.tx[, l, lo]
                        if (!is.null(aux.tn)) df$tn <- aux.tn[, l, lo]
                        if (!is.null(aux.pr)) df$pr <- aux.pr[, l, lo]
                        if (!is.null(aux.tm)) df$tm <- aux.tm[, l, lo]
                        if (!is.null(aux.hurs)) df$hurs <- aux.hurs[, l, lo]
                        if (!is.null(aux.sfcwind)) df$sfcwind <- aux.sfcwind[, l, lo]
                        if (!is.null(aux.sp)) df$sp <- aux.sp[, l, lo]
                        if (!is.null(aux.ssrd)) df$ssrd <- aux.ssrd[, l, lo]
                        if (!is.null(aux.pvpot)) df$pvpot <- aux.pvpot[, l, lo]
                        
                        args_list <- c(list(df = df, id = "id", start_date = min(df$date)), index.arg.list)
                        cdi_cei_result <- do.call(metadata$indexfun, args_list)
                        # Extract final cumulative value per season
                        cdi_cei_result %>% group_by(season_id) %>% 
                            summarize(value = max(if (index.code == "CDI") cum_days else cum_excess, na.rm = TRUE)) %>% 
                            pull(value)
                    } else {
                        # Other indices (non-FAO, non-CDI/CEI)
                        dates_mat <- cbind(
                            as.numeric(format(as.Date(refDates), "%Y")),
                            as.numeric(format(as.Date(refDates), "%m")),
                            as.numeric(format(as.Date(refDates), "%d"))
                        )
                        args_list <- index.arg.list[!names(index.arg.list) %in% c("dates")]
                        args_list[["dates"]] <- dates_mat
                        
                        if (index.code == "gsl") {
                            args_list[["tm"]] <- if (!is.null(aux.tm)) aux.tm[, l, lo] else (aux.tx[, l, lo] + aux.tn[, l, lo]) / 2
                            args_list[["lat"]] <- lats[l]
                        } else if (index.code == "avg") {
                            args_list[["tm"]] <- if (!is.null(aux.tm)) aux.tm[, l, lo] else if (!is.null(aux.tx)) aux.tx[, l, lo] else NULL
                        } else if (index.code == "nd_thre") {
                            args_list[["any"]] <- if (!is.null(aux.tm)) aux.tm[, l, lo] else if (!is.null(aux.tx)) aux.tx[, l, lo] else NULL
                        } else if (index.code == "nhw") {
                            args_list[["tx"]] <- aux.tx[, l, lo]
                        } else if (index.code == "dr") {
                            args_list[["tx"]] <- aux.tx[, l, lo]
                            args_list[["tn"]] <- aux.tn[, l, lo]
                        } else if (index.code %in% c("prcptot", "nrd", "lds", "sdii", "prcptot_thre", "ns")) {
                            args_list[["pr"]] <- aux.pr[, l, lo]
                        }
                        
                        # Call the index function
                        # For tier1 indices, call the function directly by name
                        if (index.code %in% c("gsl", "avg", "nd_thre", "nhw", "dr", "prcptot", "nrd", "lds", "sdii", "prcptot_thre", "ns")) {
                            index_result <- do.call(index.code, args_list)
                        } else {
                            index_result <- do.call(metadata$indexfun, args_list)
                        }
                        
                        # GSL returns a list, extract the GSL component
                        if (index.code == "gsl" && is.list(index_result)) {
                            index_result$GSL
                        } else {
                            index_result
                        }
                    }
                }, error = function(e) {
                    # Return NA on error - show first 10 errors for debugging
                    if (l <= 3 && lo <= 3) {
                        warning("Error at lat=", lats[l], ", lon=", lons[lo], ": ", e$message)
                    }
                    NA
                })
                
                return(result)
            })
            do.call(abind::abind, list(lonloop, along = 0))
        })
        
        # Reconstruct data array with proper dimensions
        latloop_bound <- do.call(abind::abind, list(latloop, along = 0))
        
        out_array <- unname(aperm(latloop_bound, c(3, 1, 2)))
        
        # Update dates to match the output time dimension (years for FAO indices)
        n_output_times <- dim(out_array)[1]  # First dimension is time
        if (n_output_times != length(refDates)) {
            # FAO tier1 and agronomic indices return yearly values
            # Extract unique years and create year boundaries
            years <- unique(format(as.Date(refDates), "%Y"))
            if (n_output_times == length(years)) {
                # Create start (Jan 1) and end (Dec 31) for each year
                start_dates <- as.Date(paste0(years, "-01-01"))
                end_dates <- as.Date(paste0(years, "-12-31"))
                message("[", Sys.time(), "] Converting daily dates to yearly dates (", n_output_times, " years)")
            } else {
                # For other cases, take evenly spaced dates
                indices <- round(seq(1, length(refDates), length.out = n_output_times))
                start_dates <- as.Date(refDates)[indices]
                end_dates <- start_dates
                message("[", Sys.time(), "] Adjusting dates to match output time steps (", n_output_times, ")")
            }
        } else {
            # Daily data preserved (for CDI/CEI or daily indices)
            start_dates <- as.Date(refDates)
            end_dates <- start_dates
        }
        
        # Transform to climate4R grid structure
        refGrid <- if (!is.null(tx)) {
            if (n.mem > 1) {
                tmp <- subsetGrid(tx, members = x, drop = TRUE) %>% redim(loc = FALSE, member = FALSE)
            } else {
                tmp <- tx %>% redim(loc = FALSE, member = FALSE)
            }
        } else if (!is.null(tn)) {
            if (n.mem > 1) {
                tmp <- subsetGrid(tn, members = x, drop = TRUE) %>% redim(loc = FALSE, member = FALSE)
            } else {
                tmp <- tn %>% redim(loc = FALSE, member = FALSE)
            }
        } else if (!is.null(pr)) {
            if (n.mem > 1) {
                tmp <- subsetGrid(pr, members = x, drop = TRUE) %>% redim(loc = FALSE, member = FALSE)
            } else {
                tmp <- pr %>% redim(loc = FALSE, member = FALSE)
            }
        } else if (!is.null(tm)) {
            if (n.mem > 1) {
                tmp <- subsetGrid(tm, members = x, drop = TRUE) %>% redim(loc = FALSE, member = FALSE)
            } else {
                tmp <- tm %>% redim(loc = FALSE, member = FALSE)
            }
        } else if (!is.null(hurs)) {
            if (n.mem > 1) {
                tmp <- subsetGrid(hurs, members = x, drop = TRUE) %>% redim(loc = FALSE, member = FALSE)
            } else {
                tmp <- hurs %>% redim(loc = FALSE, member = FALSE)
            }
        } else if (!is.null(sfcwind)) {
            if (n.mem > 1) {
                tmp <- subsetGrid(sfcwind, members = x, drop = TRUE) %>% redim(loc = FALSE, member = FALSE)
            } else {
                tmp <- sfcwind %>% redim(loc = FALSE, member = FALSE)
            }
        } else if (!is.null(sp)) {
            if (n.mem > 1) {
                tmp <- subsetGrid(sp, members = x, drop = TRUE) %>% redim(loc = FALSE, member = FALSE)
            } else {
                tmp <- sp %>% redim(loc = FALSE, member = FALSE)
            }
        } else if (!is.null(ssrd)) {
            if (n.mem > 1) {
                tmp <- subsetGrid(ssrd, members = x, drop = TRUE) %>% redim(loc = FALSE, member = FALSE)
            } else {
                tmp <- ssrd %>% redim(loc = FALSE, member = FALSE)
            }
        } else {
            if (n.mem > 1) {
                tmp <- subsetGrid(pvpot, members = x, drop = TRUE) %>% redim(loc = FALSE, member = FALSE)
            } else {
                tmp <- pvpot %>% redim(loc = FALSE, member = FALSE)
            }
        }
        
        # Update refGrid with calculated data and dates
        refGrid[["Data"]] <- out_array
        attr(refGrid[["Data"]], "dimensions") <- c("time", "lat", "lon")
        refGrid[["Dates"]][["start"]] <- start_dates
        refGrid[["Dates"]][["end"]] <- end_dates
        
        # Memory cleanup
        tx <- tn <- pr <- tm <- hurs <- sfcwind <- sp <- ssrd <- pvpot <- NULL
        gc()  # Force garbage collection
        return(refGrid)
        }, error = function(e) {
            warning("Error processing member ", x, ": ", e$message)
            # Return a grid with NA values but same structure
            refGrid[["Data"]][] <- NA
            return(refGrid)
        })
    })
    message("[", Sys.time(), "] Done")
    out <- if (length(out.list) == 1) {
        out.list %>% extract2(1) %>% redim(drop = TRUE)
    } else {
        suppressMessages(suppressWarnings(do.call("bindGrid", c(out.list, dimension = "member"))))
    }
    # Update Variable metadata after binding (preserve level attribute)
    out[["Variable"]] <- list("varName" = index.code, 
                              "level" = out[["Variable"]][["level"]])
    attr(out[["Variable"]], "description") <- metadata$description
    attr(out[["Variable"]], "units") <- metadata$units
    attr(out[["Variable"]], "longname") <- metadata$longname
    if (station) out %<>% redim(drop = FALSE, loc = TRUE, member = FALSE)
    invisible(out)
}





#' @title List all available Agroclimatic Indices
#' @description Print a table with a summary of the available agroclimatic indices including FAO tier1 indices and stress indices
#' @return Print a table on the screen with the following columns:
#' \itemize{
#' \item \strong{code}: Code of the index. This is the character string used as input value
#' for the argument \code{index.code} in \code{\link{agroindexGrid}}
#' \item \strong{longname}: Long description of the index
#' \item \strong{index.fun}: The name of the internal function used to calculate it
#' \item \strong{tn,tx,tm,pr,hurs,sfcwind,sp,ssrd,pvpot}: Logical values (0/1) indicating the input variables required for index calculation
#' \item \strong{units}: The units of the index (when different from those of the input variable)
#' }
#' @references FAO agroclimatic indices documentation
#' @author J. Bedia (original)
#' @export
#' @importFrom magrittr %>%

agroindexShow <- function() {
    read.master()
}



#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom utils read.table

read.master <- function() {
    # Try installed package first, then development location
    master_file <- system.file("master", package = "climate4R.agro")
    
    # If package not installed, try development paths
    if (master_file == "" || !file.exists(master_file)) {
        # Try common development locations
        dev_paths <- c(
            "C:/Users/pablo/Desktop/climate4R.agro/inst/master",  # absolute path for development
            file.path(getwd(), "inst", "master"),  # relative to working directory
            file.path(getwd(), "..", "inst", "master"),  # one level up
            "inst/master",  # relative from package root
            "../inst/master"  # relative from R/
        )
        
        for (path in dev_paths) {
            if (file.exists(path)) {
                master_file <- path
                message("Using master file from: ", path)
                break
            }
        }
    }
    
    if (master_file == "" || !file.exists(master_file)) {
        stop("Cannot find master file. Tried paths:\n",
             "  - system.file('master', package = 'climate4R.agro')\n",
             "  - C:/Users/pablo/Desktop/climate4R.agro/inst/master\n",
             "  - inst/master (relative paths)\n",
             "Please ensure climate4R.agro package is installed or master file exists in inst/ folder.")
    }
    
    read.table(master_file, 
               header = TRUE,
                                                                        sep = ";",
                                                                        stringsAsFactors = FALSE,
                                                                        na.strings = "")
}



