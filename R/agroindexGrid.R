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
#' @param ssrd A climate4R dataset of daily surface solar radiation downwards (J/m^2)
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
#' @importFrom utils getExportedValue
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
#' @author Pablo Lavin Pellon is the author of agroindexGrid.R
#' @export

agroindexGrid <- function(index.code,
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
                        ncores = NULL) {
                            
    # Validate time resolution for all input grids
    var_names <- c("tn", "tx", "pr", "tm", "hurs", "sfcwind", "ssrd")
    for (var_name in var_names) {
        var_obj <- get(var_name)
        if (!is.null(var_obj) && getTimeResolution(var_obj) != "DD") {
            stop("Daily data is required as input", call. = FALSE)
        }
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
    
    # Helper function to create grid list (used multiple times)
    create_grid_list <- function() {
        grid_list <- list()
        for (var_name in var_names) {
            var_obj <- get(var_name)
            if (!is.null(var_obj)) {
                grid_list[[var_name]] <- var_obj
            }
        }
        return(grid_list)
    }
    
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
    
    # Helper function to get grid variables as a list (used in error handlers)
    get_grid_vars_list <- function() {
        list(tx = tx, tn = tn, pr = pr, tm = tm, hurs = hurs,
             sfcwind = sfcwind, ssrd = ssrd)
    }

    first_non_null <- function(items) {
        for (obj in items) {
            if (!is.null(obj)) {
                return(obj)
            }
        }
        NULL
    }

    get_member_template_grid <- function(member_idx, n_mem) {
        grid_vars <- get_grid_vars_list()
        candidate <- first_non_null(grid_vars)
        if (is.null(candidate)) return(NULL)
        extract_member_grid(candidate, member_idx, n_mem)
    }
    
    # Helper function to load function from various paths
    load_function_from_paths <- function(fun_name, possible_paths) {
        is_absolute_path <- function(path) {
            grepl("^([A-Za-z]:|/|\\\\)", path)
        }
        extra_bases <- getOption("climate4R.agro.extra_paths")
        expanded_paths <- unique(unlist(lapply(possible_paths, function(path) {
            if (is.null(path) || path == "") {
                return(character(0))
            }
            candidates <- path
            if (!is.null(extra_bases) && length(extra_bases) > 0 && !is_absolute_path(path)) {
                candidates <- c(candidates, file.path(extra_bases, path))
            }
            candidates
        })))
        for (path in expanded_paths) {
            if (!is.null(path) && path != "") {
                actual_path <- tryCatch(normalizePath(path, winslash = "\\/", mustWork = FALSE),
                                        error = function(...) path)
                if (file.exists(actual_path)) {
                    source(actual_path, local = FALSE)
                    if (exists(fun_name, mode = "function", inherits = TRUE)) {
                        fun_obj <- get(fun_name, mode = "function", inherits = TRUE)
                        attr(fun_obj, "source_path") <- actual_path
                        return(fun_obj)
                    }
                }
            }
        }
        return(NULL)
    }
    cluster_assign_function <- function(cl, fun_name, fun_obj) {
        if (is.null(cl) || !inherits(cl, "cluster")) return(invisible(NULL))
        if (is.null(fun_obj) || !is.function(fun_obj)) return(invisible(NULL))
        tryCatch({
            parallel::clusterCall(cl, function(fn_name, fn_obj) {
                assign(fn_name, fn_obj, envir = .GlobalEnv)
                invisible(NULL)
            }, fun_name, fun_obj)
        }, error = function(e) {
            # Silently ignore export errors (e.g., serialization issues)
        })
        invisible(NULL)
    }
    cluster_load_packages <- function(cl, pkgs) {
        if (is.null(cl) || !inherits(cl, "cluster")) return(invisible(NULL))
        pkgs <- unique(pkgs)
        pkgs <- pkgs[pkgs != ""]
        if (length(pkgs) == 0) return(invisible(NULL))
        missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
        if (length(missing) > 0) {
            stop("Required package(s) not installed: ", paste(missing, collapse = ", "))
        }
        var_name <- ".climate4R_agro_pkgs"
        pkg_env <- new.env(parent = emptyenv())
        assign(var_name, pkgs, envir = pkg_env)
        parallel::clusterExport(cl, varlist = var_name, envir = pkg_env)
        results <- parallel::clusterEvalQ(cl, {
            if (!exists(".climate4R_agro_pkgs", envir = .GlobalEnv, inherits = FALSE)) {
                stop("Package list not found on worker")
            }
            pkgs_vec <- get(".climate4R_agro_pkgs", envir = .GlobalEnv, inherits = FALSE)
            ok <- sapply(pkgs_vec, function(pkg) {
                suppressPackageStartupMessages(library(pkg, character.only = TRUE))
                TRUE
            })
            rm(list = ".climate4R_agro_pkgs", envir = .GlobalEnv)
            ok
        })
        failed <- vapply(results, function(res) any(!res), logical(1))
        if (any(failed)) {
            stop("Failed to load packages on worker(s): ", paste(pkgs, collapse = ", "))
        }
        invisible(NULL)
    }
    combine_member_results <- function(member_list, station_flag, metadata_row, idx_code) {
        if (length(member_list) == 0) {
            stop("No member results to combine")
        }
        valid_members <- member_list[!sapply(member_list, is.null)]
        if (length(valid_members) == 0) {
            stop("All member results are NULL")
        }
        combined <- if (length(valid_members) == 1) {
            tryCatch({
                valid_members %>% extract2(1) %>% redim(drop = TRUE)
            }, error = function(e) {
                valid_members[[1]]
            })
        } else {
            member_arrays <- lapply(valid_members, function(g) {
                if (is.null(g) || is.null(g$Data)) return(NULL)
                g$Data
            })
            valid_indices <- !sapply(member_arrays, is.null)
            if (sum(valid_indices) == 0) {
                stop("Member grids lack Data arrays")
            }
            member_arrays <- member_arrays[valid_indices]
            valid_members <- valid_members[valid_indices]
            first_dims <- dim(member_arrays[[1]])
            for (i in seq_along(member_arrays)) {
                if (!identical(dim(member_arrays[[i]]), first_dims)) {
                    stop("Member grids have inconsistent dimensions. Expected ",
                         paste(first_dims, collapse = " x "), " but got ",
                         paste(dim(member_arrays[[i]]), collapse = " x "))
                }
            }
            combined_array <- do.call(abind::abind, c(member_arrays, list(along = 0)))
            template <- valid_members[[1]]
            template[["Data"]] <- unname(combined_array)
            data_dims <- length(dim(template[["Data"]]))
            if (data_dims == 4) {
                attr(template[["Data"]], "dimensions") <- c("member", "time", "lat", "lon")
            } else if (data_dims == 3) {
                attr(template[["Data"]], "dimensions") <- c("member", "time", "loc")
            }
            tryCatch({
                template %<>% redim(member = TRUE, drop = FALSE)
            }, error = function(e) {
                # Ignore redim issues; template already has member dimension in Data
            })
            template
        }
        level_value <- NULL
        if (!is.null(combined[["Variable"]]) && !is.null(combined[["Variable"]][["level"]])) {
            level_value <- combined[["Variable"]][["level"]]
        }
        combined[["Variable"]] <- list("varName" = idx_code,
                                        "level" = level_value)
        attr(combined[["Variable"]], "description") <- metadata_row$description
        attr(combined[["Variable"]], "units") <- metadata_row$units
        attr(combined[["Variable"]], "longname") <- metadata_row$longname
        if (station_flag) combined %<>% redim(drop = FALSE, loc = TRUE, member = FALSE)
        invisible(combined)
    }
    ensure_function_available <- function(fun_name, possible_paths) {
        if (exists(fun_name, mode = "function", inherits = TRUE)) {
            return(get(fun_name, mode = "function", inherits = TRUE))
        }
        fun_obj <- NULL
        tryCatch({
            fun_obj <<- get(fun_name, envir = asNamespace("climate4R.agro"), inherits = FALSE)
        }, error = function(...) {})
        if (is.null(fun_obj)) {
            tryCatch({
                if (requireNamespace("climate4R.agro", quietly = TRUE)) {
                    fun_obj <<- getExportedValue("climate4R.agro", fun_name)
                }
            }, error = function(...) {})
        }
        if (is.null(fun_obj)) {
            fun_obj <- load_function_from_paths(fun_name, possible_paths)
        }
        if (!is.null(fun_obj) && is.function(fun_obj) && !exists(fun_name, mode = "function", inherits = TRUE)) {
            assign(fun_name, fun_obj, envir = .GlobalEnv)
        }
        fun_obj
    }
    aux <- read.master()
    metadata <- aux[grep(paste0("^", index.code, "$"), aux$code, fixed = FALSE), ]
    # Check which variables are provided
    a <- sapply(var_names, function(v) as.numeric(!is.null(get(v))))
    missing_meta_cols <- setdiff(var_names, names(metadata))
    if (length(missing_meta_cols) > 0) {
        stop("Master metadata is missing required columns: ",
             paste(missing_meta_cols, collapse = ", "))
    }
    b <- metadata[, var_names, drop = FALSE] %>% unlist(use.names = FALSE) %>% as.numeric()
    
    # Skip variable requirement check for CDI/CEI (they accept any variables via df)
    if (!(index.code %in% c("CDI", "CEI"))) {
        if (any(b - a > 0)) {
            stop("The required input variable(s) for ", index.code,
                 " index calculation are missing\nType \'?",
                 metadata$indexfun, "\' for help", call. = FALSE)
        }
    }
    # Remove any possible unneeded input grid
    # Skip this for CDI/CEI since they accept any variables via bounds/x_col
    if (!(index.code %in% c("CDI", "CEI")) && any(a - b > 0)) {
        ind <- which((a - b) > 0)
        rem <- var_names[ind]
        for (var_name in rem) assign(var_name, NULL)
        message("NOTE: some input grids provided for ", index.code,
                " index calculation are not required and were removed")
    }
    
    # Grid consistency checks for multi-variable inputs
    grid.list <- create_grid_list()
    
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
    station <- FALSE
    for (var_name in var_names) {
        var_obj <- get(var_name)
        if (!is.null(var_obj)) {
            if (!station) station <- typeofGrid(var_obj) == "station"
            assign(var_name, redim(var_obj, member = TRUE, var = FALSE))
        }
    }
    # Get reference grid name (first non-null grid)
    refGridName <- var_names[sapply(var_names, function(v) !is.null(get(v)))][1]
    assign("refGrid", get(refGridName))
    
    # Member dimension consistency check (similar to indexGrid.R)
    grid.list.check <- create_grid_list()
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
    ref_dates_date <- as.Date(refDates)
    ref_years <- unique(format(ref_dates_date, "%Y"))
    message("[", Sys.time(), "] Calculating ", index.code, " ...")
    
    # Pre-load functions before parallel processing (for multi-member grids)
    # This ensures functions are available in parallel worker environments
    cdi_cei_fun_loaded <- NULL
    build_seasons_map_loaded <- NULL
    fao_fun_loaded <- NULL
    computeET0_loaded <- NULL
    binSpell_loaded <- NULL
    
    # Pre-load agroindexFAO function if needed (for FAO agronomic indices)
    if (metadata$indexfun == "agroindexFAO") {
        fun_name <- "agroindexFAO"
        if (exists(fun_name, mode = "function", inherits = TRUE)) {
            fao_fun_loaded <- get(fun_name, mode = "function", inherits = TRUE)
        } else {
            tryCatch({
                fao_fun_loaded <- get(fun_name, envir = asNamespace("climate4R.agro"), inherits = FALSE)
            }, error = function(e1) {
                tryCatch({
                    if (requireNamespace("climate4R.agro", quietly = TRUE)) {
                        fao_fun_loaded <<- getExportedValue("climate4R.agro", fun_name)
                    }
                }, error = function(e2) {
                    # Try to source the file directly (for development)
                    possible_paths <- list(
                        "R/indicesFAO.R",
                        system.file("R", "indicesFAO.R", package = "climate4R.agro"),
                        file.path(system.file(package = "climate4R.agro"), "..", "R", "indicesFAO.R"),
                        file.path(getwd(), "R", "indicesFAO.R")
                    )
                    fao_fun_loaded <<- load_function_from_paths(fun_name, possible_paths)
                })
            })
        }
        
        if (is.null(fao_fun_loaded) || !is.function(fao_fun_loaded)) {
            stop("Could not find function 'agroindexFAO'. ",
                 "Make sure the climate4R.agro package is loaded or source 'R/indicesFAO.R' first.")
        }
        computeET0_loaded <- ensure_function_available(
            "computeET0",
            list(
                "R/indicesFAO.R",
                system.file("R", "indicesFAO.R", package = "climate4R.agro"),
                file.path(system.file(package = "climate4R.agro"), "..", "R", "indicesFAO.R"),
                file.path(getwd(), "R", "indicesFAO.R")
            )
        )
        if (is.null(computeET0_loaded) || !is.function(computeET0_loaded)) {
            stop("Could not find function 'computeET0'. Make sure the climate4R.agro package is loaded or source 'R/indicesFAO.R'.")
        }
        binSpell_loaded <- ensure_function_available(
            "binSpell",
            list(
                "R/indicesFAO_tier1.R",
                system.file("R", "indicesFAO_tier1.R", package = "climate4R.agro"),
                file.path(system.file(package = "climate4R.agro"), "..", "R", "indicesFAO_tier1.R"),
                file.path(getwd(), "R", "indicesFAO_tier1.R")
            )
        )
        if (is.null(binSpell_loaded) || !is.function(binSpell_loaded)) {
            stop("Could not find function 'binSpell'. Make sure the climate4R.agro package is loaded or source 'R/indicesFAO_tier1.R'.")
        }
    }
    
    if (index.code %in% c("CDI", "CEI")) {
        fun_name <- index.code
        # Try to find the function using multiple methods
        if (exists(fun_name, mode = "function", inherits = TRUE)) {
            cdi_cei_fun_loaded <- get(fun_name, mode = "function", inherits = TRUE)
        } else {
            tryCatch({
                cdi_cei_fun_loaded <- get(fun_name, envir = asNamespace("climate4R.agro"), inherits = FALSE)
            }, error = function(e1) {
                tryCatch({
                    if (requireNamespace("climate4R.agro", quietly = TRUE)) {
                        cdi_cei_fun_loaded <<- getExportedValue("climate4R.agro", fun_name)
                    }
                }, error = function(e2) {
                    # Try to source the file directly (for development)
                    possible_paths <- list(
                        if (fun_name == "CDI") "R/indexCDI.R" else if (fun_name == "CEI") "R/indexCEI.R" else paste0("R/index", fun_name, ".R"),
                        system.file("R", paste0("index", fun_name, ".R"), package = "climate4R.agro"),
                        file.path(system.file(package = "climate4R.agro"), "..", "R", paste0("index", fun_name, ".R")),
                        file.path(getwd(), "R", paste0("index", fun_name, ".R"))
                    )
                    cdi_cei_fun_loaded <<- load_function_from_paths(fun_name, possible_paths)
                })
            })
        }
        
        # Also load build_seasons_map
        if (!exists("build_seasons_map", mode = "function", inherits = TRUE)) {
            build_seasons_paths <- list(
                "R/build_seasons_map.R",
                system.file("R", "build_seasons_map.R", package = "climate4R.agro"),
                file.path(system.file(package = "climate4R.agro"), "..", "R", "build_seasons_map.R"),
                file.path(getwd(), "R", "build_seasons_map.R")
            )
            build_seasons_map_loaded <- load_function_from_paths("build_seasons_map", build_seasons_paths)
        } else {
            build_seasons_map_loaded <- get("build_seasons_map", mode = "function", inherits = TRUE)
        }
        
        # Verify functions are loaded
        if (is.null(cdi_cei_fun_loaded) || !is.function(cdi_cei_fun_loaded)) {
            stop("Could not find function '", fun_name, 
                 "'. Make sure the climate4R.agro package is loaded or source 'R/index", 
                 fun_name, ".R' first. Also ensure build_seasons_map is available.")
        }
    }
    allow_member_parallel <- TRUE
    if (n.mem > 1) {
        effective_member_parallel <- isTRUE(parallel) && allow_member_parallel
        if (!allow_member_parallel && isTRUE(parallel)) {
            tryCatch({
                message("[", Sys.time(), "] Member-level parallelization disabled for FAO agronomic indices; using per-member serial execution")
            }, error = function(e) {
                # Ignore connection write errors
            })
        }
        parallel.pars <- parallelCheck(effective_member_parallel, max.ncores, ncores)
        apply_fun <- selectPar.pplyFun(parallel.pars, .pplyFUN = "lapply")
        if (parallel.pars$hasparallel) {
            # Clean up parallel cluster with proper error handling to avoid connection errors
            on.exit({
                tryCatch({
                    # Suppress warnings during cleanup to avoid connection write errors
                    suppressWarnings({
                        if (exists("parallel.pars") && !is.null(parallel.pars$cl)) {
                            parallel::stopCluster(parallel.pars$cl)
                        }
                    })
                }, error = function(e) {
                    # Silently ignore connection errors during cleanup
                    # (these occur when connection is already closed)
                })
            }, add = TRUE)

            if (!is.null(parallel.pars$cl)) {
                cluster_load_packages(parallel.pars$cl, c("transformeR", "magrittr", "abind", "dplyr"))
                export_fun_to_cluster <- function(name) {
                    if (exists(name, mode = "function", inherits = TRUE)) {
                        cluster_assign_function(parallel.pars$cl, name, get(name, mode = "function", inherits = TRUE))
                    }
                }
                if (metadata$indexfun == "agroindexFAO") {
                    if (!is.null(fao_fun_loaded)) {
                        cluster_assign_function(parallel.pars$cl, "agroindexFAO", fao_fun_loaded)
                    } else {
                        export_fun_to_cluster("agroindexFAO")
                    }
                    if (!is.null(computeET0_loaded)) {
                        cluster_assign_function(parallel.pars$cl, "computeET0", computeET0_loaded)
                    } else {
                        export_fun_to_cluster("computeET0")
                    }
                    if (!is.null(binSpell_loaded)) {
                        cluster_assign_function(parallel.pars$cl, "binSpell", binSpell_loaded)
                    } else {
                        export_fun_to_cluster("binSpell")
                    }
                }
                if (index.code %in% c("CDI", "CEI")) {
                    export_fun_to_cluster(index.code)
                    export_fun_to_cluster("build_seasons_map")
                }
                if (!is.null(cdi_cei_fun_loaded) && is.function(cdi_cei_fun_loaded)) {
                    cluster_assign_function(parallel.pars$cl, index.code, cdi_cei_fun_loaded)
                }
                if (!is.null(build_seasons_map_loaded) && is.function(build_seasons_map_loaded)) {
                    cluster_assign_function(parallel.pars$cl, "build_seasons_map", build_seasons_map_loaded)
                }
            }
        }
    } else {
        if (isTRUE(parallel)) message("NOTE: Parallel processing was skipped (unable to parallelize one single member)")
        apply_fun <- lapply
    }
    # Suppress output in parallel workers to avoid SIGPIPE and connection errors
    # Save current options and set to suppress warnings/messages in parallel mode
    old_options <- NULL
    if (n.mem > 1 && exists("parallel.pars") && parallel.pars$hasparallel) {
        old_options <- options(warn = -1)  # Suppress warnings
        # Note: We can't suppress messages with options, but we handle them via tryCatch
    }
    on.exit({
        if (!is.null(old_options)) {
            options(old_options)
        }
    }, add = TRUE)
    
    # Wrap the entire parallel computation to catch SIGPIPE and connection errors
    # This is necessary because SIGPIPE can occur at the system level
    # Use withCallingHandlers to catch warnings (like "ignoring SIGPIPE signal")
    # and tryCatch to catch errors
    out.list <- suppressWarnings({
        withCallingHandlers({
            tryCatch({
                apply_fun(1:n.mem, function(x) {
            # Comprehensive output suppression for parallel workers to avoid connection errors.
            # Only engage full suppression if member-level parallelism is actually active.
            if (n.mem > 1 && isTRUE(parallel.pars$hasparallel)) {
            # Save current state
            old_warn <- options(warn = -1)$warn
            old_message <- getOption("message", FALSE)
            
            # Suppress all output types aggressively
            options(warn = -1)  # Suppress warnings
            options(message = FALSE)  # Suppress messages
            
            # Set up sink for all output types - do this unconditionally for parallel workers
            # Use tryCatch to handle any connection errors during sink setup
            sinks_opened <- list(output = FALSE, message = FALSE)
            
            tryCatch({
                # Always sink output in parallel workers, even if sink is already active
                # This ensures we catch all output
                sink_count_output <- sink.number(type = "output")
                sink_count_message <- sink.number(type = "message")
                
                # Open sinks to nullfile to suppress all output
                if (sink_count_output == 0) {
                    sink(nullfile(), type = "output")
                    sinks_opened$output <- TRUE
                } else {
                    # If sink is already open, open another one (nested sink)
                    sink(nullfile(), type = "output", append = TRUE)
                    sinks_opened$output <- TRUE
                }
                
                if (sink_count_message == 0) {
                    sink(nullfile(), type = "message")
                    sinks_opened$message <- TRUE
                } else {
                    # If sink is already open, open another one (nested sink)
                    sink(nullfile(), type = "message", append = TRUE)
                    sinks_opened$message <- TRUE
                }
            }, error = function(e) {
                # If sink setup fails, continue anyway - we've at least suppressed via options
                sinks_opened$output <- FALSE
                sinks_opened$message <- FALSE
            })
            
            # Create safe wrapper functions that suppress all output
            safe_message <- function(...) {
                # Do nothing in parallel workers
            }
            safe_warning <- function(...) {
                # Do nothing in parallel workers
            }
            safe_cat <- function(...) {
                # Do nothing in parallel workers
            }
            safe_print <- function(...) {
                # Do nothing in parallel workers
            }
            
            # Override message/warning functions locally
            assign("message", safe_message, envir = environment())
            assign("warning", safe_warning, envir = environment())
            assign("cat", safe_cat, envir = environment())
            assign("print", safe_print, envir = environment())
            
            on.exit({
                # Restore sinks
                tryCatch({
                    if (sinks_opened$message) {
                        sink(type = "message")
                    }
                }, error = function(e) {})
                
                tryCatch({
                    if (sinks_opened$output) {
                        sink(type = "output")
                    }
                }, error = function(e) {})
                
                # Restore options
                options(warn = old_warn)
                if (!is.null(old_message)) {
                    options(message = old_message)
                }
            }, add = TRUE)
            } else {
            old_warn <- options(warn = -1)$warn
            on.exit(options(warn = old_warn), add = TRUE)
        }
        
        # Wrap entire computation in connection error handler for parallel workers
        # This catches SIGPIPE and connection read/write errors
        result <- tryCatch({
            # Inner tryCatch for the actual computation
            tryCatch({
            # Extract data as 3D arrays (time x lat x lon) following agroclimGrid pattern
            # Use helper function to simplify repetitive code
            aux.tx <- extract_data_array(tx, x, n.mem)
        aux.tn <- extract_data_array(tn, x, n.mem)
        aux.pr <- extract_data_array(pr, x, n.mem)
        aux.tm <- extract_data_array(tm, x, n.mem)
        aux.hurs <- extract_data_array(hurs, x, n.mem)
        aux.sfcwind <- extract_data_array(sfcwind, x, n.mem)
        aux.ssrd <- extract_data_array(ssrd, x, n.mem)
        # EXCEPTION for FAO AGRONOMIC INDICES (require lat, dates, and NO temporal subsetting)
        # Note: Tier1 indices (gsl, avg, etc.) are handled in the non-FAO path below
        # because they use direct function calls, not the agroindexFAO wrapper
        if (metadata$indexfun == "agroindexFAO") {
            if (time.resolution != "year") {
                # Suppress message in parallel workers to avoid connection errors
                if (n.mem == 1) {
                    message(index.code, " is calculated year by year by definition. argument time.resolution ignored.")
                }
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
                        # Suppress message in parallel workers
                        if (n.mem == 1) {
                            tryCatch({
                                message("[", Sys.time(), "] Warning: Temporal intersection failed, attempting manual alignment: ", e$message)
                            }, error = function(e2) {
                                # Ignore connection write errors
                            })
                        }
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
            lons_vec <- tryCatch(grid.list.aux[[1]][["xyCoords"]][["x"]], error = function(...) NULL)
            
            # Prepare arguments
            index.arg.list[["dates"]] <- datess
            # Only add index.code for agroindexFAO, not for tier1 functions
            if (metadata$indexfun == "agroindexFAO") {
                index.arg.list[["index.code"]] <- index.code
            }
            
            # Process lat-lon loops with error tracking and parallelization
            # Create all lat/lon combinations for parallel processing
            grid_points <- expand.grid(l = 1:length(lats), lo = 1:n_lons)
            total_points <- nrow(grid_points)
            
            # Set up parallel processing for FAO agronomic indices (lat/lon loops)
            # Parallelize at grid point level if:
            # - parallel is requested AND
            # - we have many grid points (>10) AND  
            # - (single member OR not already parallelizing at member level)
            # Note: We avoid nested parallelization (members already parallel) to prevent overhead
            use_fao_parallel <- (parallel && total_points > 10 && n.mem == 1)
            if (use_fao_parallel) {
                fao_parallel.pars <- parallelCheck(parallel, max.ncores, ncores)
                fao_apply_fun <- selectPar.pplyFun(fao_parallel.pars, .pplyFUN = "lapply")
                if (fao_parallel.pars$hasparallel) {
                    # Only print message if not in parallel member processing
                    if (n.mem == 1) {
                        core_count <- fao_parallel.pars$ncores
                        if (is.null(core_count) || is.na(core_count)) {
                            if (!is.null(fao_parallel.pars$cl)) {
                                core_count <- length(fao_parallel.pars$cl)
                            } else {
                                detected <- tryCatch(parallel::detectCores(), error = function(...) NA_integer_)
                                fallback <- if (!is.null(max.ncores) && !is.na(max.ncores)) max.ncores else detected
                                if (is.na(fallback)) fallback <- 1L
                                core_count <- min(fallback, total_points)
                            }
                        }
                        tryCatch({
                            message("[", Sys.time(), "] Parallelizing FAO agronomic index calculation over ", 
                                   total_points, " grid points using ", core_count, " cores")
                        }, error = function(e) {
                            # Ignore connection write errors
                        })
                    }
                    if (!is.null(fao_parallel.pars$cl)) {
                        cluster_load_packages(fao_parallel.pars$cl, c("transformeR", "magrittr", "abind", "dplyr"))
                        export_fun_to_cluster <- function(name) {
                            if (exists(name, mode = "function", inherits = TRUE)) {
                                cluster_assign_function(fao_parallel.pars$cl, name,
                                                        get(name, mode = "function", inherits = TRUE))
                            }
                        }
                        if (!is.null(fao_fun_loaded) && is.function(fao_fun_loaded)) {
                            cluster_assign_function(fao_parallel.pars$cl, "agroindexFAO", fao_fun_loaded)
                        } else {
                            export_fun_to_cluster("agroindexFAO")
                        }
                        if (!is.null(computeET0_loaded) && is.function(computeET0_loaded)) {
                            cluster_assign_function(fao_parallel.pars$cl, "computeET0", computeET0_loaded)
                        } else {
                            export_fun_to_cluster("computeET0")
                        }
                        if (!is.null(binSpell_loaded) && is.function(binSpell_loaded)) {
                            cluster_assign_function(fao_parallel.pars$cl, "binSpell", binSpell_loaded)
                        } else {
                            export_fun_to_cluster("binSpell")
                        }
                    }
                } else {
                    if (n.mem == 1) {
                        tryCatch({
                            message("[", Sys.time(), "] Parallel processing requested but not available (Windows or single core). Using sequential processing for ", 
                                   total_points, " grid points.")
                        }, error = function(e) {
                            # Ignore connection write errors
                        })
                    }
                    fao_apply_fun <- lapply
                }
            } else {
                fao_apply_fun <- lapply
                fao_parallel.pars <- list(hasparallel = FALSE, cl = NULL)
            }
            
            # Process all grid points (parallelized if enabled)
            error_info <- list(count = 0, display_count = 0, first_error = NULL, first_all_na = NULL)
            
            # For multi-member grids, process in smaller batches to reduce memory usage
            # Use smaller batch sizes for multi-member to prevent memory allocation failures
            if (n.mem > 1 && total_points > 50) {
                # Smaller batches for multi-member grids (50 instead of 100)
                batch_size <- min(50, total_points)
                n_batches <- ceiling(total_points / batch_size)
                grid_results <- vector("list", total_points)  # Pre-allocate
                
                for (batch_idx in 1:n_batches) {
                    start_idx <- (batch_idx - 1) * batch_size + 1
                    end_idx <- min(batch_idx * batch_size, total_points)
                    batch_indices <- start_idx:end_idx
                    
                    batch_results <- lapply(batch_indices, function(i) {
                        l <- grid_points$l[i]
                        lo <- grid_points$lo[i]
                        index.arg.list[["lat"]] <- lats[l]
                        # Extract time series for this point from each data array
                        point_data <- lapply(input.arg.list, function(z) {
                            if (is.null(z)) return(NULL)
                            z[, l, lo]
                        })
                        # Remove NULLs
                        point_data <- point_data[!sapply(point_data, is.null)]
                        
                        # Call FAO function with point data and index arguments
                        tryCatch({
                            # Use pre-loaded function if available, otherwise try lookup
                            fun_obj <- fao_fun_loaded
                            if (is.null(fun_obj) || !is.function(fun_obj)) {
                                # Fallback: try to get the function
                                fun_obj <- tryCatch({
                                    get(metadata$indexfun, mode = "function", inherits = TRUE)
                                }, error = function(e1) {
                                    tryCatch({
                                        get(metadata$indexfun, envir = asNamespace("climate4R.agro"), inherits = FALSE)
                                    }, error = function(e2) {
                                        # Last resort: try to get from global environment
                                        get(metadata$indexfun, envir = .GlobalEnv, inherits = TRUE)
                                    })
                                })
                            }
                            
                            result <- do.call(fun_obj, c(point_data, index.arg.list))
                            # GSL returns a list, extract the GSL component
                            if (index.code == "gsl" && is.list(result)) {
                                result <- result$GSL
                            }
                            if (all(is.na(result))) {
                                if (is.null(error_info$first_all_na)) {
                                    year_counts <- tryCatch(table(datess[, 1]), error = function(...) NULL)
                                    var_na_counts <- tryCatch({
                                        lapply(point_data, function(vec) sum(is.na(vec)))
                                    }, error = function(...) NULL)
                                    error_info$first_all_na <<- list(
                                        lat = lats[l],
                                        lon = if (!is.null(lons_vec)) lons_vec[lo] else NA_real_,
                                        year_counts = year_counts,
                                        var_na = var_na_counts
                                    )
                                }
                            }
                            result
                        }, error = function(e) {
                            error_info$count <<- error_info$count + 1
                            # Store first error for debugging (especially useful for multi-member)
                            if (is.null(error_info$first_error)) {
                                error_info$first_error <<- paste0("Error at lat=", lats[l], ", lon=", lo, ": ", e$message)
                            }
                            if (error_info$display_count < MAX_ERROR_DISPLAY) {
                                # Only show warnings for single-member processing
                                if (n.mem == 1) {
                                    tryCatch({
                                        warning("Error at lat=", lats[l], ", lon=", lo, ": ", e$message)
                                    }, error = function(e2) {
                                        # Ignore connection write errors during warning
                                    })
                                }
                                error_info$display_count <<- error_info$display_count + 1
                            }
                            # Return a vector of NAs with correct length (number of years)
                            n_years_expected <- length(unique(datess[, 1]))
                            rep(NA_real_, n_years_expected)
                        })
                    })
                    
                    # Store batch results efficiently
                    for (i in seq_along(batch_indices)) {
                        grid_results[[batch_indices[i]]] <- batch_results[[i]]
                    }
                    
                    # Clean up batch results and force garbage collection more aggressively
                    rm(batch_results)
                    if (batch_idx %% 3 == 0) {  # More frequent GC (every 3 batches instead of 5)
                        gc(verbose = FALSE)
                    }
                }
            } else {
                # Original processing for single-member or small grids
                grid_results <- fao_apply_fun(1:total_points, function(i) {
                    l <- grid_points$l[i]
                    lo <- grid_points$lo[i]
                    index.arg.list[["lat"]] <- lats[l]
                    # Extract time series for this point from each data array
                    point_data <- lapply(input.arg.list, function(z) z[, l, lo])
                    
                    # Call FAO function with point data and index arguments
                    tryCatch({
                        result <- do.call(metadata$indexfun, c(point_data, index.arg.list))
                        # GSL returns a list, extract the GSL component
                        if (index.code == "gsl" && is.list(result)) {
                            result <- result$GSL
                        }
                        if (all(is.na(result))) {
                            if (is.null(error_info$first_all_na)) {
                                year_counts <- tryCatch(table(datess[, 1]), error = function(...) NULL)
                                var_na_counts <- tryCatch({
                                    lapply(point_data, function(vec) sum(is.na(vec)))
                                }, error = function(...) NULL)
                                error_info$first_all_na <<- list(
                                    lat = lats[l],
                                    lon = if (!is.null(lons_vec)) lons_vec[lo] else NA_real_,
                                    year_counts = year_counts,
                                    var_na = var_na_counts
                                )
                            }
                        }
                        result
                    }, error = function(e) {
                        error_info$count <<- error_info$count + 1
                        # Store first error for debugging (especially useful for multi-member)
                        if (is.null(error_info$first_error)) {
                            error_info$first_error <<- paste0("Error at lat=", lats[l], ", lon=", lo, ": ", e$message)
                        }
                        if (error_info$display_count < MAX_ERROR_DISPLAY) {
                            # Only show warnings for single-member processing
                            if (n.mem == 1) {
                                tryCatch({
                                    warning("Error at lat=", lats[l], ", lon=", lo, ": ", e$message)
                                }, error = function(e2) {
                                    # Ignore connection write errors during warning
                                })
                            }
                            error_info$display_count <<- error_info$display_count + 1
                        }
                        # Return a vector of NAs with correct length (number of years)
                        n_years_expected <- length(unique(datess[, 1]))
                        rep(NA_real_, n_years_expected)
                    })
                })
            }
            
            # Clean up parallel cluster if it was created specifically for FAO grid points
            if (use_fao_parallel && fao_parallel.pars$hasparallel && !is.null(fao_parallel.pars$cl)) {
                tryCatch({
                    suppressWarnings(parallel::stopCluster(fao_parallel.pars$cl))
                }, error = function(e) {
                    # Ignore connection errors during cleanup
                })
            }
            
            # Report error summary for FAO indices
            if (error_info$count > 0) {
                tryCatch({
                    if (n.mem == 1) {
                        message("[", Sys.time(), "] FAO indices: ", error_info$count, " out of ", 
                               total_points, " grid points encountered errors")
                    } else {
                        # For multi-member, report first error to help debug
                        # Only report if this is NOT a connection error (to avoid nested connection issues)
                        if (!is.null(error_info$first_error)) {
                            # Check if first error is memory-related
                            if (grepl("cannot allocate|no se puede ubicar", error_info$first_error, ignore.case = TRUE)) {
                                message("[", Sys.time(), "] FAO indices (member ", x, "): MEMORY ERROR - ", 
                                       error_info$count, " out of ", total_points, " grid points failed")
                                message("  Error: ", error_info$first_error)
                            } else {
                                message("[", Sys.time(), "] FAO indices (member ", x, "): ", error_info$count, 
                                       " out of ", total_points, " grid points encountered errors")
                                message("  First error: ", error_info$first_error)
                            }
                        }
                    }
                }, error = function(e) {
                    # Ignore connection write errors
                })
            }
            
            # Check if ALL results are NA (indicates a systematic problem)
            all_na_count <- sum(sapply(grid_results, function(r) all(is.na(r))))
            if (all_na_count == total_points) {
                tryCatch({
                    if (n.mem == 1) {
                        warning("All grid points returned NA for index ", index.code, 
                               ". This may indicate a data or parameter issue.")
                    } else {
                        message("[", Sys.time(), "] FAO indices (member ", x, "): all grid points returned NA.")
                    }
                    if (!is.null(error_info$first_all_na)) {
                        diag <- error_info$first_all_na
                        message("  Example grid point lat=", diag$lat, ", lon=", diag$lon)
                        if (!is.null(diag$year_counts)) {
                            message("  Daily records per year: ",
                                    paste(names(diag$year_counts), diag$year_counts, sep = "=", collapse = "; "))
                        }
                        if (!is.null(diag$var_na)) {
                            na_str <- paste(
                                names(diag$var_na),
                                unlist(diag$var_na),
                                sep = "=",
                                collapse = "; "
                            )
                            message("  NA counts (per input series): ", na_str)
                        }
                    }
                }, error = function(e) {})
            }
            
            # Reconstruct output array: convert list of results to 3D array (time x lat x lon)
            # Vectorized approach: build array directly from results using index mapping
            # Each element in grid_results is a vector of length = number of years
            
            # Find first valid result to determine dimensions
            first_valid_idx <- which(sapply(grid_results, function(r) length(r) > 1))[1]
            if (is.na(first_valid_idx)) {
                # All results are scalar NA - something went very wrong
                # Use expected number of years from dates
                n_years <- length(unique(datess[, 1]))
                if (n.mem == 1) {
                    tryCatch({
                        warning("All grid points returned scalar NA for FAO index ", index.code, 
                               ". Using expected number of years (", n_years, ") for output structure.")
                    }, error = function(e) {})
                }
            } else {
                n_years <- length(grid_results[[first_valid_idx]])
            }
            
            # Pre-allocate output array (time x lat x lon)
            result_array <- array(NA_real_, dim = c(n_years, length(lats), n_lons))
            
            # Fill array using vectorized indexing
            # grid_points is in expand.grid order: all lons for each lat sequentially
            # Map each grid point result to correct array position
            for (i in 1:total_points) {
                l <- grid_points$l[i]
                lo <- grid_points$lo[i]
                current_result <- grid_results[[i]]
                # Handle both vector and scalar NA results
                if (length(current_result) == n_years) {
                    result_array[, l, lo] <- current_result
                } else if (length(current_result) == 1 && is.na(current_result)) {
                    # Scalar NA - fill with NA vector
                    result_array[, l, lo] <- rep(NA_real_, n_years)
                } else {
                    # Length mismatch - log warning and fill with NAs
                    if (n.mem == 1 && i <= 5) {  # Only report first 5 mismatches
                        tryCatch({
                            warning("Result length mismatch at grid point ", i, 
                                   " (lat=", lats[l], ", lon=", lo, "): expected ", n_years, 
                                   ", got ", length(current_result))
                        }, error = function(e) {})
                    }
                    result_array[, l, lo] <- rep(NA_real_, n_years)
                }
            }
            
            # Clean up large intermediate objects after reconstruction
            rm(grid_points)
            gc(verbose = FALSE)
            
            out.aux[["Data"]] <- unname(result_array)
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
                if (n.mem == 1) {
                    tryCatch({
                        message("[", Sys.time(), "] FAO index ", index.code, 
                               ": Set yearly dates (", n_output_times, " years from ", 
                               min(unique_years), " to ", max(unique_years), ")")
                    }, error = function(e) {
                        # Ignore connection write errors
                    })
                }
            } else {
                # Dimension mismatch - this indicates a problem
                # Check if maybe all values are NA or there's an error in calculation
                # Suppress warning output to avoid connection errors in batch jobs
                if (n.mem == 1) {
                    tryCatch({
                        warning("FAO index ", index.code, 
                               ": Mismatch between output time dimension (", n_output_times, 
                               ") and number of unique years (", length(unique_years), 
                               "). Output may have incorrect time structure.")
                    }, error = function(e2) {
                        # Ignore connection write errors during warning
                    })
                }
                # Try to use first n_output_times years
                if (n_output_times <= length(unique_years)) {
                    years_str <- sprintf("%04d", sort(unique_years)[1:n_output_times])
                    start_dates <- as.Date(paste0(years_str, "-01-01"))
                    end_dates <- as.Date(paste0(years_str, "-12-31"))
                    out.aux[["Dates"]][["start"]] <- start_dates
                    out.aux[["Dates"]][["end"]] <- end_dates
                }
            }
            
            # Safety check: ensure out.aux is valid before cleanup
            if (is.null(out.aux) || is.null(out.aux[["Data"]])) {
                # If out.aux is invalid, this is a serious error - let the error handler deal with it
                stop("Failed to create valid output grid for FAO agronomic index. ",
                     "This may indicate an issue with input data or grid structure.")
            }
            
            # Memory cleanup
            rm(list = c("grid_results", "result_array", "input.arg.list", "grid.list.aux"), 
               envir = environment(), inherits = FALSE)
            # Explicit garbage collection to free memory early (important for multi-member)
            gc(verbose = FALSE)
            
            return(out.aux)
        }
        
        # Set dates in index.arg.list if not already set
        if (is.null(index.arg.list[["dates"]])) {
            index.arg.list[["dates"]] <- ref_dates_date
        }
        
        # Get coordinates from refGrid (set before member processing)
        # This ensures we get the full grid coordinates before any transformations
        # refGrid is assigned from one of the original input grids before member extraction
        lats <- refGrid$xyCoords$y
        lons <- refGrid$xyCoords$x
        
        # Verify dimensions match the extracted data arrays
        # Find first available data array to verify dimensions
        ref_data_array <- NULL
        for (var_name in c("tx", "tn", "pr", "tm", "hurs", "sfcwind", "ssrd")) {
            var_arr <- get(paste0("aux.", var_name), envir = environment())
            if (!is.null(var_arr) && length(dim(var_arr)) == 3) {
                ref_data_array <- var_arr
                break
            }
        }
        
        if (is.null(ref_data_array)) {
            stop("No valid data array found to determine grid dimensions")
        }
        
        # Final verification - coordinates should match data array dimensions
        if (length(lats) != dim(ref_data_array)[2] || length(lons) != dim(ref_data_array)[3]) {
            stop("Could not determine correct coordinates. ",
                 "Data array dimensions: ", paste(dim(ref_data_array), collapse=" x "),
                 " (time x lat x lon). Expected lat=", dim(ref_data_array)[2], 
                 ", lon=", dim(ref_data_array)[3],
                 " but got lat=", length(lats), 
                 ", lon=", length(lons),
                 ". This may indicate an issue with grid extraction.")
        }
        
        # Create dates matrix (ndates x 3): year, month, day
        dates_mat <- cbind(
            as.numeric(format(ref_dates_date, "%Y")),
            as.numeric(format(ref_dates_date, "%m")),
            as.numeric(format(ref_dates_date, "%d"))
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
                            date = ref_dates_date
                        )
                        
                        # Extract time series for all variables
                        if (!is.null(aux.tx)) df$tx <- extract_ts(aux.tx, l, lo)
                        if (!is.null(aux.tn)) df$tn <- extract_ts(aux.tn, l, lo)
                        if (!is.null(aux.pr)) df$pr <- extract_ts(aux.pr, l, lo)
                        if (!is.null(aux.tm)) df$tm <- extract_ts(aux.tm, l, lo)
                        if (!is.null(aux.hurs)) df$hurs <- extract_ts(aux.hurs, l, lo)
                        if (!is.null(aux.sfcwind)) df$sfcwind <- extract_ts(aux.sfcwind, l, lo)
                        if (!is.null(aux.ssrd)) df$ssrd <- extract_ts(aux.ssrd, l, lo)
                        
                        # Remove 'dates' from args as CDI/CEI don't accept it
                        args_clean <- index.arg.list[!names(index.arg.list) %in% c("dates")]
                        args_list <- c(list(df = df, id = "id", start_date = min(df$date)), args_clean)
                        
                        # Use pre-loaded function if available, otherwise try lookup
                        cdi_cei_fun <- cdi_cei_fun_loaded
                        if (is.null(cdi_cei_fun) || !is.function(cdi_cei_fun)) {
                            # Fallback: try lookup (original code)
                            fun_name <- metadata$indexfun
                            if (exists(fun_name, mode = "function", inherits = TRUE)) {
                                cdi_cei_fun <- get(fun_name, mode = "function", inherits = TRUE)
                            } else {
                                tryCatch({
                                    cdi_cei_fun <- get(fun_name, envir = asNamespace("climate4R.agro"), inherits = FALSE)
                                }, error = function(e1) {
                                    tryCatch({
                                        if (requireNamespace("climate4R.agro", quietly = TRUE)) {
                                            cdi_cei_fun <<- getExportedValue("climate4R.agro", fun_name)
                                        }
                                    }, error = function(e2) {
                                        # Try sourcing in worker (less ideal but fallback)
                                        possible_paths <- list(
                                            if (fun_name == "CDI") "R/indexCDI.R" else if (fun_name == "CEI") "R/indexCEI.R" else paste0("R/index", fun_name, ".R"),
                                            system.file("R", paste0("index", fun_name, ".R"), package = "climate4R.agro"),
                                            file.path(system.file(package = "climate4R.agro"), "..", "R", paste0("index", fun_name, ".R")),
                                            file.path(getwd(), "R", paste0("index", fun_name, ".R"))
                                        )
                                        source_file <- NULL
                                        for (path in possible_paths) {
                                            if (!is.null(path) && path != "" && file.exists(path)) {
                                                source_file <- path
                                                break
                                            }
                                        }
                                        if (!is.null(source_file)) {
                                            source(source_file, local = FALSE)
                                            if (exists(fun_name, mode = "function", inherits = TRUE)) {
                                                cdi_cei_fun <<- get(fun_name, mode = "function", inherits = TRUE)
                                            }
                                        }
                                    })
                                })
                            }
                        }
                        
                        if (is.null(cdi_cei_fun) || !is.function(cdi_cei_fun)) {
                            # Suppress error message in parallel workers to avoid connection errors
                            error_msg <- paste0("Could not find function '", fun_name, 
                                              "'. Make sure the climate4R.agro package is loaded or source 'R/index", 
                                              fun_name, ".R' first. Also ensure build_seasons_map is available.")
                            if (n.mem == 1) {
                                stop(error_msg)
                            } else {
                                # In parallel workers, just return NA and log error silently
                                return(rep(NA, length(ref_years)))
                            }
                        }
                        
                        # Call CDI/CEI function - only suppress output for multi-member to avoid connection errors
                        if (n.mem > 1) {
                            cdi_cei_result <- suppressMessages(suppressWarnings({
                                do.call(cdi_cei_fun, args_list)
                            }))
                        } else {
                            cdi_cei_result <- do.call(cdi_cei_fun, args_list)
                        }
                        
                        # Validate CDI/CEI result structure
                        if (is.null(cdi_cei_result) || !is.data.frame(cdi_cei_result)) {
                            # Return NAs if result is invalid
                            expected_length <- length(ref_years)
                            result <- rep(NA_real_, expected_length)
                        } else {
                            # Check that required columns exist
                            has_season_id <- "season_id" %in% names(cdi_cei_result)
                            if (index.code == "CDI") {
                                has_value_col <- "cum_days" %in% names(cdi_cei_result)
                            } else {
                                has_value_col <- "cum_excess" %in% names(cdi_cei_result)
                            }
                            
                            if (!has_season_id || !has_value_col) {
                                # Missing required columns
                                if (n.mem == 1) {
                                    tryCatch({
                                        warning("CDI/CEI result missing required columns. ",
                                               "Has season_id: ", has_season_id, 
                                               ", Has value column: ", has_value_col)
                                    }, error = function(e) {})
                                }
                                expected_length <- length(ref_years)
                                result <- rep(NA_real_, expected_length)
                            } else {
                                # Extract final cumulative value per season
                                # For CDI: max of cum_days per season
                                # For CEI: max of cum_excess per season
                                if (index.code == "CDI") {
                                    result <- cdi_cei_result %>% 
                                        group_by(season_id) %>% 
                                        summarize(value = max(cum_days, na.rm = TRUE), .groups = "drop") %>% 
                                        pull(value)
                                } else {
                                    # CEI
                                    result <- cdi_cei_result %>% 
                                        group_by(season_id) %>% 
                                        summarize(value = max(cum_excess, na.rm = TRUE), .groups = "drop") %>% 
                                        pull(value)
                                }
                                
                                # Check if result length matches expected (should be number of seasons/years)
                                expected_length <- length(ref_years)
                                if (length(result) != expected_length) {
                                    # Length mismatch - this can happen if season definitions differ
                                    # Only warn for single-member and first few points to avoid spam
                                    if (n.mem == 1 && l <= 3 && lo <= 3) {
                                        tryCatch({
                                            message("Warning at lat=", lats[l], ", lon=", lons[lo], 
                                                   ": CDI/CEI result length (", length(result), 
                                                   ") doesn't match expected years (", expected_length, "). ",
                                                   "Unique season_ids: ", length(unique(cdi_cei_result$season_id)))
                                        }, error = function(e) {})
                                    }
                                }
                            }
                        }
                        
                        result
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
                        # Suppress warning output to avoid connection errors in batch jobs
                        # Only show warnings for single-member processing
                        if (n.mem == 1) {
                            tryCatch({
                                warning("Error at lat=", lats[l], ", lon=", lons[lo], ": ", e$message)
                            }, error = function(e2) {
                                # Ignore connection write errors during warning
                            })
                        }
                        error_display_count <<- error_display_count + 1
                    }
                    years <- ref_years
                    rep(NA, length(years))
                })
                
                return(result)
            })
            do.call(abind::abind, list(lonloop, along = 0))
        })
        
        # Report error summary for non-FAO indices
        if (error_count > 0 && n.mem == 1) {
            tryCatch({
                message("[", Sys.time(), "] Non-FAO indices: ", error_count, " out of ", 
                        length(lats) * length(lons), " grid points encountered errors")
            }, error = function(e) {
                # Ignore connection write errors
            })
        }
        
        # Reconstruct data array with proper dimensions
        latloop_bound <- do.call(abind::abind, list(latloop, along = 0))
        out_array <- unname(aperm(latloop_bound, c(3, 1, 2)))
        
        # Update dates to match the output time dimension (years for FAO indices)
        n_output_times <- dim(out_array)[1]  # First dimension is time
        if (n_output_times != length(refDates)) {
            # FAO tier1 and agronomic indices return yearly values
            # Extract unique years and create year boundaries
            years <- ref_years
            if (n_output_times == length(years)) {
                # Create start (Jan 1) and end (Dec 31) for each year
                start_dates <- as.Date(paste0(years, "-01-01"))
                end_dates <- as.Date(paste0(years, "-12-31"))
                if (n.mem == 1) {
                    tryCatch({
                        message("[", Sys.time(), "] Converting daily dates to yearly dates (", n_output_times, " years)")
                    }, error = function(e) {
                        # Ignore connection write errors
                    })
                }
            } else {
                # For other cases, take evenly spaced dates
                indices <- round(seq(1, length(refDates), length.out = n_output_times))
                start_dates <- ref_dates_date[indices]
                end_dates <- start_dates
                if (n.mem == 1) {
                    tryCatch({
                        message("[", Sys.time(), "] Adjusting dates to match output time steps (", n_output_times, ")")
                    }, error = function(e) {
                        # Ignore connection write errors
                    })
                }
            }
        } else {
            # Daily data preserved (for CDI/CEI or daily indices)
            start_dates <- ref_dates_date
            end_dates <- start_dates
        }
        
        # Transform to climate4R grid structure using helper function
        # Extract member grid to use as template (needed for proper grid structure)
        refGrid_member <- NULL
        refGrid_member <- get_member_template_grid(x, n.mem)
        
        if (is.null(refGrid_member)) {
            stop("No input grid available to use as template")
        }
        
        # Use refGrid_member as template for output grid structure
        refGrid_template <- refGrid_member
        
        # Update template grid with calculated data and dates
        # Ensure the template grid has the correct spatial dimensions
        refGrid_template[["Data"]] <- out_array
        attr(refGrid_template[["Data"]], "dimensions") <- c("time", "lat", "lon")
        refGrid_template[["Dates"]][["start"]] <- start_dates
        refGrid_template[["Dates"]][["end"]] <- end_dates
        
        # Update xyCoords to match the actual data dimensions if needed
        if (length(refGrid_template$xyCoords$y) != length(lats) || 
            length(refGrid_template$xyCoords$x) != length(lons)) {
            refGrid_template$xyCoords$y <- lats
            refGrid_template$xyCoords$x <- lons
        }
        
        # Memory cleanup
        rm(list = c("aux.tx", "aux.tn", "aux.pr", "aux.tm", "aux.hurs", 
                   "aux.sfcwind", "aux.ssrd", 
                   "latloop", "latloop_bound", "out_array"), 
           envir = environment(), inherits = FALSE)
        # Explicit garbage collection to free memory early (important for multi-member)
        gc(verbose = FALSE)
        refGrid_template  # Return the result
        }, error = function(e) {
            # Suppress warning output to avoid connection errors in batch jobs
            # The error is already captured in the error handler
            if (n.mem == 1) {
                tryCatch({
                    warning("Error processing member ", x, ": ", e$message)
                }, error = function(e2) {
                    # Ignore connection write errors during warning
                })
            }
            # Create a grid with NA values but same structure as other members
            # Use helper function to simplify code
            refGrid_member <- get_member_template_grid(x, n.mem)

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
            refGrid_member[["Dates"]][["start"]] <- ref_dates_date
            refGrid_member[["Dates"]][["end"]] <- ref_dates_date
            refGrid_member  # Return error grid
            })
        }, error = function(e) {
            # Outer error handler: catch all errors and return a valid grid structure
            # This ensures we always return something that can be combined, even if it's all NAs
            error_msg <- if (is.character(e$message)) e$message else as.character(e)
            
            # Check if this is a connection/SIGPIPE error
            is_connection_error <- (
                grepl("conexi[oó]n", error_msg, ignore.case = TRUE) ||
                grepl("connection", error_msg, ignore.case = TRUE) ||
                grepl("SIGPIPE", error_msg, ignore.case = TRUE) ||
                grepl("broken pipe", error_msg, ignore.case = TRUE) ||
                grepl("ignoring SIGPIPE", error_msg, ignore.case = TRUE)
            )
            
            # For multi-member grids, always try to return a valid grid structure (even if all NAs)
            # This allows the computation to continue even if one worker has issues
            # SIGPIPE errors are common in parallel processing and should be handled gracefully
            if (n.mem > 1 || is_connection_error) {
                tryCatch({
                    refGrid_member <- get_member_template_grid(x, n.mem)
                    
                    if (is.null(refGrid_member)) {
                        # Fallback: try to use refGrid from outer scope
                        tryCatch({
                            refGrid_member <- if (n.mem > 1) {
                                subsetGrid(refGrid, members = x, drop = TRUE) %>% redim(loc = FALSE, member = FALSE)
                            } else {
                                refGrid %>% redim(loc = FALSE, member = FALSE)
                            }
                        }, error = function(e3) {
                            # If that fails, try to get any available grid
                            if (!is.null(tx)) {
                                refGrid_member <<- extract_member_grid(tx, x, n.mem)
                            } else if (!is.null(tn)) {
                                refGrid_member <<- extract_member_grid(tn, x, n.mem)
                            } else if (!is.null(pr)) {
                                refGrid_member <<- extract_member_grid(pr, x, n.mem)
                            } else if (!is.null(tm)) {
                                refGrid_member <<- extract_member_grid(tm, x, n.mem)
                            }
                        })
                    }
                    
                    if (!is.null(refGrid_member) && !is.null(refGrid_member[["Data"]])) {
                        # Create a grid with NA values but same structure
                        # For FAO agronomic indices, we need yearly structure
                        # Check if metadata is available and if this is an FAO agronomic index
                        is_fao_agro <- tryCatch({
                            !is.null(metadata) && !is.null(metadata$indexfun) && metadata$indexfun == "agroindexFAO"
                        }, error = function(e) {
                            # If metadata check fails, assume it's not FAO agronomic
                            FALSE
                        })
                        
                        if (is_fao_agro) {
                            # Get unique years from refDates
                            years <- ref_years
                            n_years <- length(years)
                            # Get spatial dimensions
                            if (length(dim(refGrid_member[["Data"]])) == 3) {
                                # 3D: time x lat x lon
                                n_lat <- dim(refGrid_member[["Data"]])[2]
                                n_lon <- dim(refGrid_member[["Data"]])[3]
                                refGrid_member[["Data"]] <- array(NA, dim = c(n_years, n_lat, n_lon))
                                attr(refGrid_member[["Data"]], "dimensions") <- c("time", "lat", "lon")
                            } else if (length(dim(refGrid_member[["Data"]])) == 2) {
                                # 2D: time x space (station data)
                                n_space <- dim(refGrid_member[["Data"]])[2]
                                refGrid_member[["Data"]] <- array(NA, dim = c(n_years, n_space))
                                attr(refGrid_member[["Data"]], "dimensions") <- c("time", "loc")
                            }
                            # Set yearly dates
                            start_dates <- as.Date(paste0(years, "-01-01"))
                            end_dates <- as.Date(paste0(years, "-12-31"))
                            refGrid_member[["Dates"]][["start"]] <- start_dates
                            refGrid_member[["Dates"]][["end"]] <- end_dates
                        } else {
                            # For other indices, use original dates
                            refGrid_member[["Data"]][] <- NA
                            refGrid_member[["Dates"]][["start"]] <- ref_dates_date
                            refGrid_member[["Dates"]][["end"]] <- ref_dates_date
                        }
                        return(refGrid_member)
                    } else {
                        # Last resort: return NULL (will be handled by combination code)
                        return(NULL)
                    }
                }, error = function(e2) {
                    # If even error recovery fails, return NULL and let the outer code handle it
                    # But log the error for debugging (only for single-member to avoid connection issues)
                    if (n.mem == 1) {
                        tryCatch({
                            warning("Error recovery failed for member ", x, ": ", e2$message, 
                                   " (original error: ", error_msg, ")")
                        }, error = function(e3) {})
                    }
                    return(NULL)
                })
            } else {
                # For single-member, re-throw the error so user can see what went wrong
                stop(e)
            }
        })
        
        # Return the result from the inner computation
        return(result)
        })
            }, error = function(e) {
        # Top-level error handler for parallel computation
        # This should rarely be reached since inner error handlers should catch most errors
        # SIGPIPE and connection errors should be caught by the outer error handler
        error_msg <- if (is.character(e$message)) e$message else as.character(e)
        
        # Check if this is a connection/SIGPIPE error
        is_connection_error <- (
            grepl("conexi[oó]n", error_msg, ignore.case = TRUE) ||
            grepl("connection", error_msg, ignore.case = TRUE) ||
            grepl("SIGPIPE", error_msg, ignore.case = TRUE) ||
            grepl("broken pipe", error_msg, ignore.case = TRUE) ||
            grepl("ignoring SIGPIPE", error_msg, ignore.case = TRUE)
        )
        
        # If this is a SIGPIPE/connection error and we have multi-member grids,
        # the outer error handler should have caught it. If we reach here, it means
        # the error occurred before the outer handler could catch it.
        # For SIGPIPE errors, try to create valid grids for all members
        if (is_connection_error && n.mem > 1) {
            # Create a list of valid (but empty) grids for all members
            # This allows the computation to continue
            tryCatch({
                ref_grid_template <- get_member_template_grid(1, n.mem)
                
                if (is.null(ref_grid_template)) {
                    ref_grid_template <- if (n.mem > 1) {
                        subsetGrid(refGrid, members = 1, drop = TRUE) %>% redim(loc = FALSE, member = FALSE)
                    } else {
                        refGrid %>% redim(loc = FALSE, member = FALSE)
                    }
                }
                
                # Create grids for all members with NA values
                out_list_recovery <- lapply(1:n.mem, function(mem_idx) {
                    mem_grid <- if (mem_idx == 1) {
                        ref_grid_template
                    } else {
                        # For other members, try to extract their grid
                        tryCatch({
                            member_grid_candidate <- get_member_template_grid(mem_idx, n.mem)
                            if (!is.null(member_grid_candidate)) member_grid_candidate else ref_grid_template
                        }, error = function(e) {
                            ref_grid_template
                        })
                    }
                    
                    # Fill with NAs and set appropriate structure
                    if (!is.null(mem_grid) && !is.null(mem_grid[["Data"]])) {
                        # Check if this is FAO agronomic index
                        is_fao_agro <- tryCatch({
                            !is.null(metadata) && !is.null(metadata$indexfun) && metadata$indexfun == "agroindexFAO"
                        }, error = function(e) FALSE)
                        
                        if (is_fao_agro) {
                            # Yearly structure for FAO agronomic indices
                            years <- ref_years
                            n_years <- length(years)
                            if (length(dim(mem_grid[["Data"]])) == 3) {
                                n_lat <- dim(mem_grid[["Data"]])[2]
                                n_lon <- dim(mem_grid[["Data"]])[3]
                                mem_grid[["Data"]] <- array(NA, dim = c(n_years, n_lat, n_lon))
                                attr(mem_grid[["Data"]], "dimensions") <- c("time", "lat", "lon")
                            }
                            start_dates <- as.Date(paste0(years, "-01-01"))
                            end_dates <- as.Date(paste0(years, "-12-31"))
                            mem_grid[["Dates"]][["start"]] <- start_dates
                            mem_grid[["Dates"]][["end"]] <- end_dates
                        } else {
                            mem_grid[["Data"]][] <- NA
                            mem_grid[["Dates"]][["start"]] <- ref_dates_date
                            mem_grid[["Dates"]][["end"]] <- ref_dates_date
                        }
                    }
                    return(mem_grid)
                })
                
                return(out_list_recovery)
            }, error = function(e2) {
                # If recovery fails, re-throw the original error
                stop("Error during parallel computation: ", error_msg, 
                     "\nRecovery attempt also failed: ", e2$message)
            })
        } else {
            # For non-connection errors or single-member, re-throw
            stop("Error during parallel computation: ", error_msg, 
                 "\nThis error occurred at the top level of member processing. ",
                 "Check that all required functions and data are available in parallel workers.")
        }
            })
        }, warning = function(w) {
            # Catch SIGPIPE warnings specifically in withCallingHandlers
            warn_msg <- if (is.character(w$message)) w$message else as.character(w)
            if (grepl("SIGPIPE", warn_msg, ignore.case = TRUE) && n.mem > 1) {
                # Suppress SIGPIPE warnings in parallel mode
                invokeRestart("muffleWarning")
            }
        }, message = function(m) {
            # Catch SIGPIPE messages (sometimes printed as messages, not warnings)
            msg_text <- if (is.character(m$message)) m$message else as.character(m)
            if (grepl("SIGPIPE", msg_text, ignore.case = TRUE) && n.mem > 1) {
                # Suppress SIGPIPE messages in parallel mode
                invokeRestart("muffleMessage")
            }
        })
    })
    
    # Suppress message output to avoid connection errors in batch jobs
    tryCatch({
        message("[", Sys.time(), "] Done")
    }, error = function(e) {
        # Ignore connection write errors during message
    })
    
    out <- combine_member_results(out.list, station, metadata, index.code)
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
#' \item \strong{tn,tx,tm,pr,hurs,sfcwind,ssrd}: Logical values (0/1) indicating the input variables required for index calculation
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
    
    master_dir <- dirname(master_file)
    extra_paths <- unique(na.omit(c(getOption("climate4R.agro.extra_paths"),
                                    master_dir,
                                    normalizePath(file.path(master_dir, ".."), mustWork = FALSE),
                                    normalizePath(file.path(master_dir, "..", ".."), mustWork = FALSE))))
    options(climate4R.agro.extra_paths = extra_paths)
    read.table(master_file, 
               header = TRUE,
               sep = ";",
               stringsAsFactors = FALSE,
               na.strings = "")
}



