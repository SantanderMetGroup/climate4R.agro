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
#' @importFrom dplyr group_by summarize pull
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
    

    # Constants
    MAX_ERROR_DISPLAY <- 10  # Maximum number of errors to display per member
    
    # Helper function to extract 1D time series from 3D arrays (time x lat x lon)
    # After redim(member = FALSE), arrays are always 3D in climate4R
    extract_ts <- function(arr, l, lo) {
        if (is.null(arr)) return(NULL)
        # Arrays are already 3D (time x lat x lon) after redim(member = FALSE)
        return(arr[, l, lo])
    }
    
    # Helper function to extract member grid and remove member dimension
    # Returns grid with 3D Data array: (time x lat x lon)
    extract_member_grid <- function(grid, member_idx, n_mem) {
        if (is.null(grid)) return(NULL)
        if (n_mem > 1) {
            result <- subsetGrid(grid, members = member_idx, drop = TRUE) %>% redim(loc = FALSE, member = FALSE)
        } else {
            result <- grid %>% redim(loc = FALSE, member = FALSE)
        }
        
        # Ensure member dimension is actually removed from Data array
        # Sometimes redim doesn't fully remove the dimension if it's size 1
        if (!is.null(result) && !is.null(result[["Data"]])) {
            data_arr <- result[["Data"]]
            dim_names <- attr(data_arr, "dimensions")
            
            # If we still have 4 dimensions and one is member (size 1), remove it manually
            if (length(dim(data_arr)) == 4) {
                if (!is.null(dim_names) && "member" %in% dim_names) {
                    member_dim_idx <- which(dim_names == "member")
                    if (dim(data_arr)[member_dim_idx] == 1) {
                        # Select the first (and only) slice along the member dimension
                        if (member_dim_idx == 1) {
                            result[["Data"]] <- data_arr[1, , , , drop = TRUE]
                        } else if (member_dim_idx == 2) {
                            result[["Data"]] <- data_arr[, 1, , , drop = TRUE]
                        } else if (member_dim_idx == 3) {
                            result[["Data"]] <- data_arr[, , 1, , drop = TRUE]
                        } else if (member_dim_idx == 4) {
                            result[["Data"]] <- data_arr[, , , 1, drop = TRUE]
                        }
                        # Update dimensions attribute
                        new_dim_names <- dim_names[-member_dim_idx]
                        attr(result[["Data"]], "dimensions") <- new_dim_names
                    }
                } else if (dim(data_arr)[1] == 1) {
                    # Assume first dimension is member if size 1
                    result[["Data"]] <- array(data_arr, dim = dim(data_arr)[-1])
                }
            }
        }
        
        return(result)
    }
    
    # Helper function to extract data array from grid (returns 3D array)
    extract_data_array <- function(grid, member_idx, n_mem) {
        tmp <- extract_member_grid(grid, member_idx, n_mem)
        if (is.null(tmp)) return(NULL)
        tmp[["Data"]]
    }
    aux <- read.master()
    metadata <- aux[grep(paste0("^", index.code, "$"), aux$code, fixed = FALSE), ]
    a <- c(!is.null(tn), !is.null(tx), !is.null(pr), !is.null(tm), !is.null(hurs), 
           !is.null(sfcwind), !is.null(sp), !is.null(ssrd), !is.null(pvpot)) %>% as.numeric()
    b <- metadata[ , 4:12] %>% as.numeric()
    
    # Skip variable requirement check for CDI/CEI (they accept any variables via df)
    if (!(index.code %in% c("CDI", "CEI"))) {
        if (any(b - a > 0)) {
            stop("The required input variable(s) for ", index.code,
                 " index calculation are missing\nType \'?",
                 metadata$indexfun, "\' for help", call. = FALSE)
        }
    }
    # Remove any possible uneeded input grid
    # Skip this for CDI/CEI since they accept any variables via bounds/x_col
    if (!(index.code %in% c("CDI", "CEI"))) {
        if (any(a - b > 0)) {
            ind <- which((a - b) > 0)
            rem <- c("tn", "tx", "pr", "tm", "hurs", "sfcwind", "sp", "ssrd", "pvpot")[ind]
            sapply(rem, function(x) assign(x, NULL)) %>% invisible()
            message("NOTE: some input grids provided for ", index.code,
                    " index calculation are not required and were removed")
        }
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
        # Use helper function to simplify repetitive code
        aux.tx <- extract_data_array(tx, x, n.mem)
        aux.tn <- extract_data_array(tn, x, n.mem)
        aux.pr <- extract_data_array(pr, x, n.mem)
        aux.tm <- extract_data_array(tm, x, n.mem)
        aux.hurs <- extract_data_array(hurs, x, n.mem)
        aux.sfcwind <- extract_data_array(sfcwind, x, n.mem)
        aux.sp <- extract_data_array(sp, x, n.mem)
        aux.ssrd <- extract_data_array(ssrd, x, n.mem)
        aux.pvpot <- extract_data_array(pvpot, x, n.mem)
        # EXCEPTION for FAO AGRONOMIC INDICES (require lat, dates, and NO temporal subsetting)
        # Note: Tier1 indices (gsl, avg, etc.) are handled in the non-FAO path below
        # because they use direct function calls, not the agroindexFAO wrapper
        if (metadata$indexfun == "agroindexFAO") {
            if (time.resolution != "year") {
                message(index.code, " is calculated year by year by definition. argument time.resolution ignored.")
            }
            
            # Build grid.list.aux with processed grids for this member (consistent with indexGrid.R)
            # Extract member x and remove member dimension to get (time, lat, lon)
            # Use helper function to simplify repetitive code
            grid.list.aux <- list()
            if (!is.null(tx)) grid.list.aux[["tx"]] <- extract_member_grid(tx, x, n.mem)
            if (!is.null(tn)) grid.list.aux[["tn"]] <- extract_member_grid(tn, x, n.mem)
            if (!is.null(pr)) grid.list.aux[["pr"]] <- extract_member_grid(pr, x, n.mem)
            if (!is.null(tm)) grid.list.aux[["tm"]] <- extract_member_grid(tm, x, n.mem)
            
            # Ensure temporal consistency after member extraction (critical for multi-member grids)
            # After subsetGrid, grids might have slightly different temporal structures
            if (length(grid.list.aux) > 1) {
                original_names_aux <- names(grid.list.aux)
                # Check if grids have the same dates before intersecting
                dates_list <- lapply(grid.list.aux, function(g) {
                    if (!is.null(g) && !is.null(g$Dates) && !is.null(g$Dates$start)) {
                        return(as.Date(g$Dates$start))
                    }
                    return(NULL)
                })
                dates_list <- dates_list[!sapply(dates_list, is.null)]
                
                # Only intersect if dates differ
                dates_differ <- FALSE
                if (length(dates_list) > 1) {
                    first_dates <- dates_list[[1]]
                    for (i in 2:length(dates_list)) {
                        if (!identical(first_dates, dates_list[[i]])) {
                            dates_differ <- TRUE
                            break
                        }
                    }
                }
                
                if (dates_differ) {
                    tryCatch({
                        grid.list.aux <- intersectGrid(grid.list.aux, type = "temporal", which.return = 1:length(grid.list.aux))
                        names(grid.list.aux) <- original_names_aux
                    }, error = function(e) {
                        # If intersection fails, try to align dates manually
                        message("[", Sys.time(), "] Warning: Temporal intersection failed, attempting manual alignment: ", e$message)
                        # Find the grid with the shortest time period
                        time_lengths <- sapply(grid.list.aux, function(g) {
                            if (!is.null(g) && !is.null(g$Dates) && !is.null(g$Dates$start)) {
                                return(length(g$Dates$start))
                            }
                            return(Inf)
                        })
                        min_length <- min(time_lengths)
                        ref_grid_idx <- which(time_lengths == min_length)[1]
                        
                        # Extract reference dates
                        ref_dates <- grid.list.aux[[ref_grid_idx]]$Dates$start
                        ref_dates_end <- grid.list.aux[[ref_grid_idx]]$Dates$end
                        
                        # Align all grids to reference dates
                        for (i in 1:length(grid.list.aux)) {
                            if (i != ref_grid_idx && !is.null(grid.list.aux[[i]])) {
                                current_dates <- as.Date(grid.list.aux[[i]]$Dates$start)
                                # Find overlapping dates
                                date_match <- match(ref_dates, current_dates)
                                valid_dates <- !is.na(date_match)
                                
                                if (any(valid_dates)) {
                                    # Subset grid to matching dates
                                    date_indices <- date_match[valid_dates]
                                    # Get dimensions to handle correctly
                                    data_dims <- dim(grid.list.aux[[i]]$Data)
                                    if (length(data_dims) == 3) {
                                        # 3D: time x lat x lon
                                        grid.list.aux[[i]]$Data <- grid.list.aux[[i]]$Data[date_indices, , , drop = FALSE]
                                    } else if (length(data_dims) == 2) {
                                        # 2D: time x space
                                        grid.list.aux[[i]]$Data <- grid.list.aux[[i]]$Data[date_indices, , drop = FALSE]
                                    } else {
                                        stop("Unexpected data dimensions: ", paste(data_dims, collapse=" x "))
                                    }
                                    grid.list.aux[[i]]$Dates$start <- ref_dates[valid_dates]
                                    grid.list.aux[[i]]$Dates$end <- ref_dates_end[valid_dates]
                                } else {
                                    stop("No overlapping dates found between grids after member extraction")
                                }
                            }
                        }
                    })
                }
            }
            
            # Create output grid structure using aggregateGrid on sugrid (sets yearly structure)
            out.aux <- suppressMessages(aggregateGrid(grid.list.aux[[1]], aggr.y = list(FUN = "mean", na.rm = TRUE)))
            
            # Extract data arrays (3D: time x lat x lon) from processed grids
            input.arg.list <- lapply(grid.list.aux, function(d) {
                if (is.null(d)) return(NULL)
                data_arr <- d[["Data"]]
                data_dims <- dim(data_arr)
                n_dims <- length(data_dims)
                
                # Handle different dimension cases
                if (n_dims == 4) {
                    # 4D array: likely member x time x lat x lon or time x lat x lon x something
                    # Check dimension names if available
                    dim_names <- attr(data_arr, "dimensions")
                    if (!is.null(dim_names)) {
                        # If member dimension is present and size 1, drop it
                        if ("member" %in% dim_names && data_dims[which(dim_names == "member")] == 1) {
                            # Remove member dimension by selecting first (and only) member
                            member_idx <- which(dim_names == "member")
                            if (member_idx == 1) {
                                data_arr <- data_arr[1, , , , drop = TRUE]
                            } else if (member_idx == 2) {
                                data_arr <- data_arr[, 1, , , drop = TRUE]
                            } else if (member_idx == 3) {
                                data_arr <- data_arr[, , 1, , drop = TRUE]
                            } else if (member_idx == 4) {
                                data_arr <- data_arr[, , , 1, drop = TRUE]
                            }
                        }
                    } else {
                        # No dimension names - assume first dimension is member if it's size 1
                        if (data_dims[1] == 1) {
                            data_arr <- array(data_arr, dim = data_dims[-1])
                        } else {
                            stop("Cannot handle 4D array without dimension names. Dimensions: ", 
                                 paste(data_dims, collapse=" x "))
                        }
                    }
                } else if (n_dims == 3) {
                    # Already 3D, keep as is
                    # No action needed
                } else if (n_dims == 2) {
                    # 2D array might be time x space (for station data)
                    # This is acceptable
                } else {
                    stop("Unexpected array dimensions: ", n_dims, " dimensions (", 
                         paste(data_dims, collapse=" x "), 
                         "). Expected 3D (time x lat x lon) or 4D with removable member dimension.")
                }
                
                # Final check - ensure we have at least 2 dimensions
                if (length(dim(data_arr)) < 2) {
                    stop("Array collapsed to fewer than 2 dimensions. Original dimensions: ", 
                         paste(data_dims, collapse=" x "))
                }
                
                return(data_arr)
            })
            
            # Verify all arrays have the same dimensions
            dims_list <- lapply(input.arg.list, function(x) if (!is.null(x)) dim(x) else NULL)
            dims_list <- dims_list[!sapply(dims_list, is.null)]
            if (length(dims_list) > 1) {
                first_dims <- dims_list[[1]]
                for (i in 2:length(dims_list)) {
                    if (!identical(dims_list[[i]], first_dims)) {
                        stop("Data arrays have inconsistent dimensions. Expected: ", 
                             paste(first_dims, collapse=" x "), 
                             " but got: ", paste(dims_list[[i]], collapse=" x "))
                    }
                }
            }
            
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
            
            # Process lat-lon loops with error tracking
            error_count_fao <- 0
            error_display_count_fao <- 0
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
                        error_count_fao <<- error_count_fao + 1
                        if (error_display_count_fao < MAX_ERROR_DISPLAY) {
                            warning("Error at lat=", lats[l], ", lon=", lo, ": ", e$message)
                            error_display_count_fao <<- error_display_count_fao + 1
                        }
                        NA
                    })
                })
                do.call(abind::abind, list(lonloop, along = 0))
            })
            
            # Report error summary for FAO indices
            if (error_count_fao > 0) {
                message("[", Sys.time(), "] FAO indices: ", error_count_fao, " grid points encountered errors")
            }
            
            # Reconstruct output array
            # latloop_bound has dimensions (lat, lon, time) from abind
            # We need (time, lat, lon), so use aperm(..., c(3, 1, 2))
            latloop_bound <- do.call(abind::abind, list(latloop, along = 0))
            out.aux[["Data"]] <- unname(aperm(latloop_bound, c(3, 1, 2)))
            attr(out.aux[["Data"]], "dimensions") <- c("time", "lat", "lon")
            
            # Update dates to yearly (FAO agronomic indices return one value per year)
            # The output should have one value per year
            n_output_times <- dim(out.aux[["Data"]])[1]
            
            # Extract unique years from the date matrix (datess has columns: year, month, day)
            unique_years <- unique(datess[, 1])
            
            if (n_output_times == length(unique_years)) {
                # Create yearly dates (Jan 1 to Dec 31 for each year)
                years_str <- sprintf("%04d", unique_years)  # Ensure 4-digit years
                start_dates <- as.Date(paste0(years_str, "-01-01"))
                end_dates <- as.Date(paste0(years_str, "-12-31"))
                out.aux[["Dates"]][["start"]] <- start_dates
                out.aux[["Dates"]][["end"]] <- end_dates
                message("[", Sys.time(), "] FAO index ", index.code, 
                       ": Set yearly dates (", n_output_times, " years from ", 
                       min(unique_years), " to ", max(unique_years), ")")
            } else {
                # Dimension mismatch - this indicates a problem
                # Check if maybe all values are NA or there's an error in calculation
                warning("FAO index ", index.code, 
                       ": Mismatch between output time dimension (", n_output_times, 
                       ") and number of unique years (", length(unique_years), 
                       "). Output may have incorrect time structure.")
                # Try to use first n_output_times years
                if (n_output_times <= length(unique_years)) {
                    years_str <- sprintf("%04d", sort(unique_years)[1:n_output_times])
                    start_dates <- as.Date(paste0(years_str, "-01-01"))
                    end_dates <- as.Date(paste0(years_str, "-12-31"))
                    out.aux[["Dates"]][["start"]] <- start_dates
                    out.aux[["Dates"]][["end"]] <- end_dates
                }
            }
            
            # Memory cleanup
            rm(list = c("latloop", "latloop_bound", "input.arg.list", "grid.list.aux"), 
               envir = environment(), inherits = FALSE)
            tx <- tn <- pr <- tm <- hurs <- sfcwind <- sp <- ssrd <- pvpot <- NULL
            return(out.aux)
        }
        
        # Set dates in index.arg.list if not already set
        if (is.null(index.arg.list[["dates"]])) {
            index.arg.list[["dates"]] <- as.Date(refDates)
        }
        
        # Get coordinates (needed for non-FAO indices)
        lats <- getGrid(refGrid)$y
        lons <- getGrid(refGrid)$x
        
        # Create dates matrix (ndates x 3): year, month, day
        dates_mat <- cbind(
            as.numeric(format(as.Date(refDates), "%Y")),
            as.numeric(format(as.Date(refDates), "%m")),
            as.numeric(format(as.Date(refDates), "%d"))
        )
        dates_mat <- as.matrix(dates_mat)
        colnames(dates_mat) <- NULL
        
        # Prepare base args_list (without dates)
        base_args_list <- index.arg.list[!names(index.arg.list) %in% c("dates")]
        
        # Error tracking for non-FAO indices
        error_count <- 0
        error_display_count <- 0
        
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
                        
                        # Extract time series for all variables
                        if (!is.null(aux.tx)) df$tx <- extract_ts(aux.tx, l, lo)
                        if (!is.null(aux.tn)) df$tn <- extract_ts(aux.tn, l, lo)
                        if (!is.null(aux.pr)) df$pr <- extract_ts(aux.pr, l, lo)
                        if (!is.null(aux.tm)) df$tm <- extract_ts(aux.tm, l, lo)
                        if (!is.null(aux.hurs)) df$hurs <- extract_ts(aux.hurs, l, lo)
                        if (!is.null(aux.sfcwind)) df$sfcwind <- extract_ts(aux.sfcwind, l, lo)
                        if (!is.null(aux.sp)) df$sp <- extract_ts(aux.sp, l, lo)
                        if (!is.null(aux.ssrd)) df$ssrd <- extract_ts(aux.ssrd, l, lo)
                        if (!is.null(aux.pvpot)) df$pvpot <- extract_ts(aux.pvpot, l, lo)
                        
                        # Remove 'dates' from args as CDI/CEI don't accept it
                        args_clean <- index.arg.list[!names(index.arg.list) %in% c("dates")]
                        args_list <- c(list(df = df, id = "id", start_date = min(df$date)), args_clean)
                        cdi_cei_result <- do.call(metadata$indexfun, args_list)
                        # Extract final cumulative value per season
                        cdi_cei_result %>% group_by(season_id) %>% 
                            summarize(value = max(if (index.code == "CDI") cum_days else cum_excess, na.rm = TRUE)) %>% 
                            pull(value)
                    } else {
                        # Other indices (non-FAO, non-CDI/CEI)
                        args_list <- base_args_list
                        
                        # Build arguments for tier1 indices with proper NULL checks
                        if (index.code == "gsl") {
                            if (!is.null(aux.tm)) {
                                args_list[["tm"]] <- aux.tm[, l, lo]
                            } else if (!is.null(aux.tx) && !is.null(aux.tn)) {
                                args_list[["tm"]] <- (aux.tx[, l, lo] + aux.tn[, l, lo]) / 2
                            } else {
                                stop("GSL requires tm or both tx and tn")
                            }
                            args_list[["lat"]] <- lats[l]
                        } else if (index.code == "avg") {
                            if (!is.null(aux.tm)) {
                                args_list[["tm"]] <- aux.tm[, l, lo]
                            } else if (!is.null(aux.tx)) {
                                args_list[["tm"]] <- aux.tx[, l, lo]
                            } else {
                                stop("avg requires tm or tx")
                            }
                        } else if (index.code == "nd_thre") {
                            if (!is.null(aux.tm)) {
                                args_list[["any"]] <- aux.tm[, l, lo]
                            } else if (!is.null(aux.tx)) {
                                args_list[["any"]] <- aux.tx[, l, lo]
                            } else {
                                stop("nd_thre requires tm or tx")
                            }
                        } else if (index.code == "nhw") {
                            if (is.null(aux.tx)) stop("nhw requires tx")
                            args_list[["tx"]] <- aux.tx[, l, lo]
                        } else if (index.code == "dr") {
                            if (is.null(aux.tx)) stop("dr requires tx")
                            if (is.null(aux.tn)) stop("dr requires tn")
                            args_list[["tx"]] <- aux.tx[, l, lo]
                            args_list[["tn"]] <- aux.tn[, l, lo]
                        } else if (index.code %in% c("prcptot", "nrd", "lds", "sdii", "prcptot_thre", "ns")) {
                            if (is.null(aux.pr)) stop(index.code, " requires pr")
                            args_list[["pr"]] <- aux.pr[, l, lo]
                        }
                        
                        # Add dates matrix to args_list
                        args_list[["dates"]] <- dates_mat
                        
                        # Call the index function
                        # For tier1 indices, use agroindexFAO_tier1 wrapper which handles function lookup
                        if (index.code %in% c("gsl", "avg", "nd_thre", "nhw", "dr", "prcptot", "nrd", "lds", "sdii", "prcptot_thre", "ns")) {
                            # Use the wrapper function which internally calls the correct function
                            # This is more robust than direct function lookup
                            tryCatch({
                                tier1_wrapper <- get("agroindexFAO_tier1", mode = "function", inherits = TRUE)
                                # The wrapper expects index.code as first arg, then ...
                                args_with_code <- c(list(index.code = index.code), args_list)
                                index_result <- do.call(tier1_wrapper, args_with_code)
                            }, error = function(e) {
                                # Fallback: try direct function lookup
                                tier1_fun <- NULL
                                tryCatch({
                                    tier1_fun <- get(index.code, mode = "function", inherits = TRUE)
                                }, error = function(e1) {
                                    tryCatch({
                                        tier1_fun <<- get(index.code, envir = asNamespace("climate4R.agro"), inherits = FALSE)
                                    }, error = function(e2) {
                                        stop("Could not find function '", index.code, 
                                             "' or wrapper 'agroindexFAO_tier1'. ",
                                             "Make sure the climate4R.agro package is loaded or source 'R/indicesFAO_tier1.R' first.")
                                    })
                                })
                                if (is.null(tier1_fun) || !is.function(tier1_fun)) {
                                    stop("Could not find function '", index.code, "'")
                                }
                                index_result <<- do.call(tier1_fun, args_list)
                            })
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
                    error_count <<- error_count + 1
                    if (error_display_count < MAX_ERROR_DISPLAY) {
                        warning("Error at lat=", lats[l], ", lon=", lons[lo], ": ", e$message)
                        error_display_count <<- error_display_count + 1
                    }
                    years <- unique(format(as.Date(refDates), "%Y"))
                    rep(NA, length(years))
                })
                
                return(result)
            })
            do.call(abind::abind, list(lonloop, along = 0))
        })
        
        # Report error summary for non-FAO indices
        if (error_count > 0) {
            total_points <- length(lats) * length(lons)
            message("[", Sys.time(), "] Non-FAO indices: ", error_count, " out of ", 
                    total_points, " grid points encountered errors")
        }
        
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
        
        # Transform to climate4R grid structure using helper function
        # Find first available grid to use as template
        refGrid <- NULL
        grid_vars <- list(tx = tx, tn = tn, pr = pr, tm = tm, hurs = hurs, 
                         sfcwind = sfcwind, sp = sp, ssrd = ssrd, pvpot = pvpot)
        for (var_name in names(grid_vars)) {
            if (!is.null(grid_vars[[var_name]])) {
                refGrid <- extract_member_grid(grid_vars[[var_name]], x, n.mem)
                break
            }
        }
        
        if (is.null(refGrid)) {
            stop("No input grid available to use as template")
        }
        
        # Update refGrid with calculated data and dates
        refGrid[["Data"]] <- out_array
        attr(refGrid[["Data"]], "dimensions") <- c("time", "lat", "lon")
        refGrid[["Dates"]][["start"]] <- start_dates
        refGrid[["Dates"]][["end"]] <- end_dates
        
        # Memory cleanup
        rm(list = c("aux.tx", "aux.tn", "aux.pr", "aux.tm", "aux.hurs", 
                   "aux.sfcwind", "aux.sp", "aux.ssrd", "aux.pvpot", 
                   "latloop", "latloop_bound", "out_array"), 
           envir = environment(), inherits = FALSE)
        tx <- tn <- pr <- tm <- hurs <- sfcwind <- sp <- ssrd <- pvpot <- NULL
        gc()  # Force garbage collection
        return(refGrid)
        }, error = function(e) {
            warning("Error processing member ", x, ": ", e$message)
            # Create a grid with NA values but same structure as other members
            # Use helper function to simplify code
            grid_vars <- list(tx = tx, tn = tn, pr = pr, tm = tm, hurs = hurs, 
                             sfcwind = sfcwind, sp = sp, ssrd = ssrd, pvpot = pvpot)
            refGrid_member <- NULL
            for (var_name in names(grid_vars)) {
                if (!is.null(grid_vars[[var_name]])) {
                    refGrid_member <- extract_member_grid(grid_vars[[var_name]], x, n.mem)
                    break
                }
            }
            
            if (is.null(refGrid_member)) {
                # Fallback to refGrid from outer scope
                refGrid_member <- if (n.mem > 1) {
                    subsetGrid(refGrid, members = x, drop = TRUE) %>% redim(loc = FALSE, member = FALSE)
                } else {
                    refGrid %>% redim(loc = FALSE, member = FALSE)
                }
            }
            
            # Fill with NA and set appropriate dates
            refGrid_member[["Data"]][] <- NA
            refGrid_member[["Dates"]][["start"]] <- as.Date(refDates)
            refGrid_member[["Dates"]][["end"]] <- as.Date(refDates)
            return(refGrid_member)
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
            '/lustre/gmeteo/WORK/lavinp/test/inst/master',
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



