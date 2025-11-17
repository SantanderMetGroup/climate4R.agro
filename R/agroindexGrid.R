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

#' @noRd
#' @template templateParallelParams
#' @import transformeR
#' @importFrom parallel stopCluster
#' @importFrom utils head capture.output str globalVariables
#' @importFrom abind abind
#' @details \code{\link{agroindexShow}} will display on screen a full list of
#' available agroclimatic indices and their codes. Index groups include:
#' (i) FAO Tier1 indices such as gsl, avg, nd_thre, nhw, dr, prcptot, nrd,
#' lds, sdii, prcptot_thre, and ns; (ii) FAO agronomic season indices derived
#' from the water balance (dt_st_rnagsn, nm_flst_rnagsn, dt_fnst_rnagsn,
#' dt_ed_rnagsn, dl_agsn, dc_agsn, rn_agsn, avrn_agsn, dc_rnlg_agsn, tm_agsn,
#' dc_txh_agsn, dc_tnh_agsn); and (iii) stress indices CDI (Condition Duration
#' Index) and CEI (Condition Excess Index). For index-specific arguments, refer
#' to the help files of individual index functions.
#'
#' @template templateParallel
#'
#' @examples
#' \dontrun{
#'   library(transformeR)
#'   library(visualizeR)
#'   # Load temperature and precipitation data
#'   data("tasmax.grid")
#'   data("tasmin.grid")
#'   data("pr.grid")
#'
#'   ## Example 1: Growing Season Length (GSL)
#'   gsl.grid <- agroindexGrid(
#'     tm = tasmax.grid,
#'     index.code = "gsl",
#'     time.resolution = "year",
#'     index.arg.list = list(lat = 40)
#'   )
#
#'   ## Example 2: Average temperature with custom season
#'   avg.grid <- agroindexGrid(
#'     tm = tasmax.grid,
#'     index.code = "avg",
#'     time.resolution = "year",
#'     index.arg.list = list(
#'       year.start = "2000-06-01",
#'       year.end = "2000-08-31"
#'     )
#'   )
#'
#'   ## Example 3: Number of heat waves
#'   nhw.grid <- agroindexGrid(
#'     tx = tasmax.grid,
#'     index.code = "nhw",
#'     time.resolution = "year",
#'     index.arg.list = list(
#'       threshold = 35,
#'       duration = 3
#'     )
#'   )
#'
#'   ## Example 4: CDI - Condition Duration Index (multi-variable stress)
#'   ## Requires bounds specification for each variable
#'   cdi.grid <- agroindexGrid(
#'     tx = tasmax.grid,
#'     hurs = hurs.grid,
#'     index.code = "CDI",
#'     time.resolution = "year",
#'     index.arg.list = list(
#'       bounds = data.frame(
#'         var = c("tx", "hurs"),
#'         lower = c(30, 0),
#'         upper = c(Inf, 40)
#'       ),
#'       combiner = "all",
#'       min_duration = 3
#'     )
#'   )
#' }
#'
#' @author Pablo Lavin Pellon is the author of agroindexGrid.R

#' @name get_time_resolution
#' @title Get time resolution of a grid
#' @description Wrapper function for \code{\link[transformeR]{getTimeResolution}} 
#' to extract the time resolution from a climate grid object.
#' @param ... Arguments passed to \code{\link[transformeR]{getTimeResolution}}
#' @return A character string indicating the time resolution (e.g., "DD" for daily, 
#' "MM" for monthly, "YY" for yearly)
#' @seealso \code{\link[transformeR]{getTimeResolution}}
#' @export
get_time_resolution <- function(...) {
  transformeR::getTimeResolution(...)
}

#' @title Get shape/dimensions of a grid
#' @description Wrapper function for \code{\link[transformeR]{getShape}} 
#' to extract the dimensions (shape) of a climate grid object.
#' @param ... Arguments passed to \code{\link[transformeR]{getShape}}
#' @return A list or vector containing the dimensions of the grid (e.g., time, 
#' lat, lon, member)
#' @seealso \code{\link[transformeR]{getShape}}
#' @export
get_shape <- function(...) {
  transformeR::getShape(...)
}

#' @title Redimension a grid
#' @description Wrapper function for \code{\link[transformeR]{redim}} 
#' to reshape or redimension a climate grid object by adding or removing dimensions.
#' @param ... Arguments passed to \code{\link[transformeR]{redim}}
#' @return A redimensioned grid object
#' @seealso \code{\link[transformeR]{redim}}
#' @export
redim_tr <- function(...) {
  transformeR::redim(...)
}

subset_grid <- transformeR::subsetGrid 
subset_dimension <- transformeR::subsetDimension 
get_grid <- transformeR::getGrid 
interp_grid <- transformeR::interpGrid 
typeof_grid <- transformeR::typeofGrid 
get_ref_dates <- transformeR::getRefDates 
is_regular <- transformeR::isRegular 
intersect_grid <- transformeR::intersectGrid 
parallel_check <- transformeR::parallelCheck 
select_par_pply_fun <- transformeR::selectPar.pplyFun 

# Initialize symbols used via NSE to avoid R CMD check notes
tn <- tx <- pr <- tm <- hurs <- sfcwind <- ssrd <- NULL

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "season_id", "cum_days", "cum_excess", "value", "tx"
  ))
}

#' @title Agroclimatic Index Calculation for Grid Data
#' @description Calculate agroclimatic indices from climate grid data
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
                          ncores = NULL,
                          verbose = FALSE,
                          ._retry = FALSE) {
    original_call <- match.call()

    # Validate time resolution for all input grids
    var_names <- c("tn", "tx", "pr", "tm", "hurs", "sfcwind", "ssrd")
    for (var_name in var_names) {
        var_obj <- get(var_name)
        if (!is.null(var_obj) && get_time_resolution(var_obj) != "DD") {
            stop("Daily data is required as input", call. = FALSE)
        }
    }
    
    # Validate time resolution
    time.resolution <- match.arg(time.resolution, choices = c("month", "year", "climatology"))
    original_call$time.resolution <- time.resolution
    
    index.code <- match.arg(index.code,
                            choices = c("gsl", "avg", "nd_thre", "nhw", "dr", 
                                        "prcptot", "nrd", "lds", "sdii", "prcptot_thre", "ns",
                                        "dt_st_rnagsn", "nm_flst_rnagsn", "dt_fnst_rnagsn", 
                                        "dt_ed_rnagsn", "dl_agsn", "dc_agsn", "rn_agsn", 
                                        "avrn_agsn", "dc_rnlg_agsn", "tm_agsn", 
                                        "dc_txh_agsn", "dc_tnh_agsn", "CDI", "CEI"))
    original_call$index.code <- index.code
    original_call$verbose <- verbose
    

    # Constants
    MAX_ERROR_DISPLAY <- 10  # Maximum number of errors to display per member
    
    # Helper to normalise date metadata (drop names, ensure Date class)
    normalize_grid_dates <- function(grid_obj) {
        if (is.null(grid_obj) || is.null(grid_obj$Dates)) return(grid_obj)
        for (slot in c("start", "end")) {
            dates_vec <- grid_obj$Dates[[slot]]
            if (!is.null(dates_vec)) {
                grid_obj$Dates[[slot]] <- as.Date(dates_vec)
                if (!is.null(names(grid_obj$Dates[[slot]]))) {
                    names(grid_obj$Dates[[slot]]) <- NULL
                }
            }
        }
        grid_obj
    }
    
    # Helper function to create grid list (used multiple times)
    create_grid_list <- function() {
        grid_list <- list()
        for (var_name in var_names) {
            var_obj <- get(var_name)
            if (!is.null(var_obj)) {
                grid_list[[var_name]] <- normalize_grid_dates(var_obj)
            }
        }
        return(grid_list)
    }
    # Helper to ensure member dimension with memory-safe fallback
    ensure_member_dimension <- function(grid_obj) {
        if (is.null(grid_obj)) return(NULL)
        member_size <- tryCatch({
            get_shape(grid_obj, "member")
        }, error = function(...) NA_real_)
        if (!is.na(member_size) && member_size > 0) {
            return(normalize_grid_dates(grid_obj))
        }
        tryCatch({
            redim_tr(grid_obj, member = TRUE, var = FALSE)
        }, error = function(e) {
            msg <- conditionMessage(e)
            mem_error <- grepl("cannot allocate", msg, ignore.case = TRUE) ||
                         grepl("no se puede ubicar", msg, ignore.case = TRUE)
            if (!mem_error) stop(e)
            data_arr <- grid_obj[["Data"]]
            if (is.null(data_arr)) stop(e)
            original_dims <- dim(data_arr)
            new_data <- array(data_arr, dim = c(1, original_dims))
            dim_names <- attr(data_arr, "dimensions")
            if (is.null(dim_names)) {
                default_names <- switch(length(original_dims),
                    `3` = c("time", "lat", "lon"),
                    `2` = c("time", "loc"),
                    `1` = c("time"),
                    NULL)
                dim_names <- default_names
            }
            attr(new_data, "dimensions") <- c("member", dim_names)
            grid_obj[["Data"]] <- new_data
            if (!is.null(grid_obj[["Members"]])) {
                if (length(grid_obj[["Members"]]) == 0) {
                    grid_obj[["Members"]] <- list("member" = "member_1")
                }
            } else {
                grid_obj[["Members"]] <- list("member" = "member_1")
            }
            normalize_grid_dates(grid_obj)
        })
    }
    
    # Helper function to extract 1D time series from 3D arrays (time x lat x lon)
    # After redim(member = FALSE), arrays are always 3D in climate4R
    extract_ts <- function(arr, l, lo) {
        if (is.null(arr)) return(NULL)
        # Arrays are already 3D (time x lat x lon) after redim(member = FALSE)
        return(arr[, l, lo])
    }
    
    # Helper function to safely drop dimensions from a grid
    # Only drops dimensions that actually exist to avoid warnings
    safe_redim_drop <- function(grid, drop_member = TRUE, drop_loc = TRUE) {
        if (is.null(grid)) return(NULL)
        
        # Check which dimensions exist before trying to drop them
        grid_dims <- tryCatch({
            get_shape(grid)
        }, error = function(e) NULL)
        
        # Build redim arguments based on existing dimensions
        redim_args <- list()
        if (!is.null(grid_dims)) {
            # Only drop dimensions that exist
            if (drop_member && "member" %in% names(grid_dims)) {
                redim_args[["member"]] <- FALSE
            }
            if (drop_loc && "loc" %in% names(grid_dims)) {
                redim_args[["loc"]] <- FALSE
            }
        } else {
            # Fallback: try to drop requested dimensions (may generate warning)
            if (drop_member) redim_args[["member"]] <- FALSE
            if (drop_loc) redim_args[["loc"]] <- FALSE
        }
        
        if (length(redim_args) > 0) {
            do.call(redim_tr, c(list(grid), redim_args))
        } else {
            grid
        }
    }
    
    # Helper function to extract member grid and remove member dimension
    # Returns grid with 3D Data array: (time x lat x lon)
    extract_member_grid <- function(grid, member_idx, n_mem) {
        if (is.null(grid)) return(NULL)
        
        if (n_mem > 1) {
            grid_subset <- subset_grid(grid, members = member_idx, drop = TRUE)
            result <- safe_redim_drop(grid_subset, drop_member = TRUE, drop_loc = TRUE)
        } else {
            result <- safe_redim_drop(grid, drop_member = TRUE, drop_loc = TRUE)
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
        result <- normalize_grid_dates(result)
        if (isTRUE(verbose)) {
            res_class <- if (is.null(result)) "NULL" else paste(class(result), collapse = "/")
            message("[", Sys.time(), "] Member ", member_idx, ": finished with result class ", res_class)
            if (is.null(result)) {
                message("[", Sys.time(), "] Member ", member_idx, ": result is NULL (check member_errors).")
                err_attr <- attr(result, ".agroindex_error", exact = TRUE)
                if (is.null(err_attr)) {
                    err_attr <- attr(refGrid, ".agroindex_error", exact = TRUE)
                }
                if (!is.null(err_attr)) {
                    message("[", Sys.time(), "] Member ", member_idx, ": error detail -> ", err_attr$message)
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
    get_grid_vars_list <- function(tx, tn, pr, tm, hurs, sfcwind, ssrd) {
        list(tx = tx, tn = tn, pr = pr, tm = tm, hurs = hurs,
             sfcwind = sfcwind, ssrd = ssrd)
    }
    # Helper function to get first non-null item
    first_non_null <- function(items) {
        for (obj in items) {
            if (!is.null(obj)) {
                return(obj)
            }
        }
        NULL
    }
    # Helper function to get member template grid
    get_member_template_grid <- function(member_idx, n_mem,
                                         tx, tn, pr, tm, hurs, sfcwind, ssrd) {
        grid_vars <- get_grid_vars_list(tx, tn, pr, tm, hurs, sfcwind, ssrd)
        candidate <- first_non_null(grid_vars)
        if (is.null(candidate)) return(NULL)
        extract_member_grid(candidate, member_idx, n_mem)
    }
    
    # Helper function to assign function to cluster
    cluster_assign_function <- function(cl, fun_name, fun_obj) {
        if (is.null(cl) || !inherits(cl, "cluster")) return(invisible(NULL))
        if (is.null(fun_obj) || !is.function(fun_obj)) return(invisible(NULL))
        tryCatch({
            tmp_env <- new.env(parent = emptyenv())
            assign(fun_name, fun_obj, envir = tmp_env)
            parallel::clusterExport(cl, varlist = fun_name, envir = tmp_env)
        }, error = function(e) {
            # Silently ignore export errors (e.g., serialization issues)
        })
        invisible(NULL)
    }
    # Helper function to load packages on cluster
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
    # Helper function to verify that required functions exist on workers
    cluster_verify_functions <- function(cl, fun_names) {
        if (is.null(cl) || !inherits(cl, "cluster")) return(invisible(NULL))
        fun_names <- unique(fun_names)
        fun_names <- fun_names[fun_names != ""]
        if (length(fun_names) == 0) return(invisible(NULL))
        var_name <- ".climate4R_required_functions"
        fun_env <- new.env(parent = emptyenv())
        assign(var_name, fun_names, envir = fun_env)
        parallel::clusterExport(cl, varlist = var_name, envir = fun_env)
        results <- parallel::clusterEvalQ(cl, {
            if (!exists(".climate4R_required_functions", envir = .GlobalEnv, inherits = FALSE)) {
                stop("Function list not found on worker")
            }
            fn_list <- get(".climate4R_required_functions", envir = .GlobalEnv, inherits = FALSE)
            ok <- sapply(fn_list, function(fn) exists(fn, mode = "function", inherits = TRUE))
            rm(list = ".climate4R_required_functions", envir = .GlobalEnv)
            ok
        })
        missing <- character(0)
        if (is.list(results)) {
            missing <- unique(unlist(lapply(results, function(vals) {
                if (is.null(vals)) return(fun_names)
                fun_names[!vals]
            })))
        } else if (is.atomic(results)) {
            missing <- fun_names[!as.logical(results)]
        }
        if (length(missing) > 0) {
            stop("Cluster workers are missing required function(s): ",
                 paste(missing, collapse = ", "),
                 ". Ensure the necessary packages are loaded on each worker.")
        }
        invisible(NULL)
    }
    # Helper function to combine member results
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
                redim_tr(valid_members[[1]], drop = TRUE)
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
                template <- redim_tr(template, member = TRUE, drop = FALSE)
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
        if (station_flag) combined <- redim_tr(combined, drop = FALSE, loc = TRUE, member = FALSE)
        invisible(combined)
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
    b <- as.numeric(unlist(metadata[, var_names, drop = FALSE], use.names = FALSE))
    
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
        locs <- lapply(grid.list, is_regular)
        if (!any(sum(unlist(locs)) != 0, sum(unlist(locs)) != length(grid.list))) {
            stop("Regular and Irregular grids can not be combined. See function interpGrid")
        }
        
        # Intersect grids temporally
        original_names <- names(grid.list)
        safe_temporal_intersection <- function(gl) {
            tryCatch({
                intersect_grid(gl, type = "temporal", which.return = seq_along(gl))
            }, error = function(e) {
                msg <- conditionMessage(e)
                mem_error <- grepl("cannot allocate", msg, ignore.case = TRUE) ||
                             grepl("no se puede ubicar", msg, ignore.case = TRUE)
                if (!mem_error) stop(e)
                common_dates <- Reduce(
                    intersect,
                    lapply(gl, function(g) {
                        if (is.null(g)) return(NULL)
                        as.Date(g$Dates$start)
                    })
                )
                if (length(common_dates) == 0) {
                    stop("Temporal intersection resulted in an empty set of dates after memory-safe fallback.")
                }
                fall_back <- lapply(gl, function(g) {
                    if (is.null(g)) return(NULL)
                    grid_dates <- as.Date(g$Dates$start)
                    idx <- match(common_dates, grid_dates)
                    idx <- idx[!is.na(idx)]
                    if (length(idx) == 0) {
                        stop("Fallback intersection found no matching dates for at least one grid.")
                    }
                    normalize_grid_dates(
                        subset_dimension(g, dimension = "time", indices = idx)
                    )
                })
                fall_back
            })
        }
        grid.list <- safe_temporal_intersection(grid.list)
        grid.list <- lapply(grid.list, normalize_grid_dates)
        namesgridlist <- original_names
        
        # Check spatial consistency and interpolate if needed
        refgrid <- get_grid(grid.list[[1]])
        indinterp <- which(isFALSE(unlist(lapply(seq_along(grid.list)[-1], function(i) 
            identical(refgrid, get_grid(grid.list[[i]]))))))
        
        if (length(indinterp) > 0) {
            grid.list.aux <- suppressMessages(lapply(grid.list[indinterp], function(i) 
                interp_grid(i, get_grid(grid.list[[1]]))))
            grid.list[indinterp] <- lapply(grid.list.aux, normalize_grid_dates)
            grid.list.aux <- NULL
        }
        
        # Update variables with processed grids
        for (i in seq_along(grid.list)) {
            assign(namesgridlist[i], grid.list[[i]])
        }
    }
    # Ensure member is present in data structures / handle station <--> grid
    station <- FALSE
    for (var_name in var_names) {
        var_obj <- get(var_name)
        if (!is.null(var_obj)) {
            if (!station) station <- typeof_grid(var_obj) == "station"
            assign(var_name, ensure_member_dimension(var_obj))
        }
    }
    # Get reference grid name (first non-null grid)
    refGridName <- var_names[sapply(var_names, function(v) !is.null(get(v)))][1]
    assign("refGrid", get(refGridName))
    
    # Member dimension consistency check (similar to indexGrid.R)
    grid.list.check <- create_grid_list()
    if (length(grid.list.check) > 1) {
        ns.mem <- lapply(grid.list.check, function(r) get_shape(r)[["member"]])
        # Handle NA (no member dimension) as 1
        ns.mem <- lapply(ns.mem, function(x) if (is.na(x)) 1 else x)
        if (sum(unlist(ns.mem) - rep(ns.mem[[1]], length(ns.mem))) != 0) {
            stop("Number of members is different across input grids")
        }
        n.mem <- unique(unlist(ns.mem))
    } else {
        n.mem <- get_shape(refGrid, "member")
        if (is.na(n.mem)) n.mem <- 1  # Handle grids without member dimension
    }
    # Get reference dates (as Date, FAO/CDI/CEI)
    refDates <- get_ref_dates(refGrid)
    ref_dates_date <- as.Date(refDates)
    ref_years <- unique(format(ref_dates_date, "%Y"))
    message("[", Sys.time(), "] Calculating ", index.code, " ...")
    member_error_summary <- list()
    
    # Pre-load functions before parallel processing (for multi-member grids)
    # This ensures functions are available in parallel worker environments
    cdi_cei_fun_loaded <- NULL
    build_seasons_map_loaded <- NULL
    fao_fun_loaded <- NULL
    computeET0_loaded <- NULL
    binSpell_loaded <- NULL
    tier1_fun_loaded <- NULL
    tier1_index_fun_loaded <- NULL
    yearStartEnd_loaded <- NULL
    
    # Tier1 indices that need special function exports (defined once)
    tier1_indices <- c("gsl", "avg", "nd_thre", "nhw", "dr", "prcptot", "nrd", "lds", "sdii", "prcptot_thre", "ns")
    is_tier1 <- index.code %in% tier1_indices
    is_cdi_cei <- index.code %in% c("CDI", "CEI")
    is_fao_agro <- metadata$indexfun == "agroindexFAO"
    
    # Helper function to load functions efficiently (avoids repetitive tryCatch)
    load_function <- function(fun_name, required = TRUE) {
        fun <- tryCatch({
            get(fun_name, mode = "function", inherits = TRUE)
        }, error = function(...) {
            ns_loaded <- "climate4R.agro" %in% loadedNamespaces()
            if (ns_loaded) {
                utils::getFromNamespace(fun_name, ns = "climate4R.agro")
            } else if (required) {
                stop("Function '", fun_name, "' not found. Please ensure 'climate4R.agro' is installed and loaded.")
            } else {
                NULL
            }
        })
        fun
    }
    
    # Pre-load functions only if needed (avoid redundant loading)
    if (is_tier1) {
        tier1_fun_loaded <- load_function("agroindexFAO_tier1")
        tier1_index_fun_loaded <- load_function(index.code)
        binSpell_loaded <- load_function("binSpell")
        yearStartEnd_loaded <- load_function("yearStartEnd")
    }
    
    if (is_fao_agro) {
        fao_fun_loaded <- agroindexFAO
        computeET0_loaded <- load_function("computeET0")
        # binSpell may already be loaded for Tier1, avoid redundant loading
        if (is.null(binSpell_loaded)) {
            binSpell_loaded <- load_function("binSpell")
        }
    }
    
    if (is_cdi_cei) {
        cdi_cei_fun_loaded <- if (index.code == "CDI") CDI else CEI
        build_seasons_map_loaded <- build_seasons_map
    }
    allow_member_parallel <- TRUE
    member_parallel_active <- FALSE
    if (n.mem > 1) {
        effective_member_parallel <- isTRUE(parallel) && allow_member_parallel
        parallel.pars <- parallel_check(effective_member_parallel, max.ncores, ncores)
        apply_fun <- select_par_pply_fun(parallel.pars, .pplyFUN = "lapply")
        member_parallel_active <- isTRUE(parallel.pars$hasparallel)
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
                cluster_load_packages(parallel.pars$cl, c("transformeR", "abind", "climate4R.agro"))
                
                # Export functions to parallel workers (only what's needed)
                verify_funs_list <- character(0)
                
                if (is_tier1) {
                    # Export Tier1 functions to parallel workers
                    if (is.null(tier1_fun_loaded) || is.null(tier1_index_fun_loaded) || 
                        is.null(binSpell_loaded) || is.null(yearStartEnd_loaded)) {
                        stop("Tier1 functions failed to load. Please ensure 'climate4R.agro' is properly installed.")
                    }
                    # For GSL (and other Tier1 indices that use binSpell), ensure binSpell is available
                    # Export binSpell first so it's in .GlobalEnv before gsl is exported
                    cluster_assign_function(parallel.pars$cl, "binSpell", binSpell_loaded)
                    cluster_assign_function(parallel.pars$cl, "yearStartEnd", yearStartEnd_loaded)
                    cluster_assign_function(parallel.pars$cl, "agroindexFAO_tier1", tier1_fun_loaded)
                    
                    # For gsl specifically: create a wrapper that ensures binSpell is accessible
                    # When gsl calls binSpell internally, it needs to find it in its environment
                    # Since exported functions may lose their namespace context, we wrap gsl
                    if (index.code == "gsl") {
                        # Create a wrapper function with binSpell in its closure
                        # This ensures that when gsl calls binSpell, it can find it
                        gsl_wrapper <- (function(original_gsl, binSpell_fn) {
                            function(tm, dates, lat, pnan = 25) {
                                # Ensure binSpell is in .GlobalEnv as fallback
                                if (!exists("binSpell", mode = "function", envir = .GlobalEnv, inherits = FALSE)) {
                                    assign("binSpell", binSpell_fn, envir = .GlobalEnv)
                                }
                                # Try to set binSpell in gsl's environment before calling
                                gsl_env <- environment(original_gsl)
                                if (!is.null(gsl_env)) {
                                    tryCatch({
                                        assign("binSpell", binSpell_fn, envir = gsl_env, inherits = FALSE)
                                    }, error = function(e) {
                                        # If assignment fails, binSpell is still in .GlobalEnv
                                    })
                                }
                                # Call original gsl function
                                original_gsl(tm = tm, dates = dates, lat = lat, pnan = pnan)
                            }
                        })(tier1_index_fun_loaded, binSpell_loaded)
                        # Export wrapper instead of original
                        cluster_assign_function(parallel.pars$cl, index.code, gsl_wrapper)
                    } else {
                        cluster_assign_function(parallel.pars$cl, index.code, tier1_index_fun_loaded)
                    }
                    
                    # Ensure binSpell is accessible in .GlobalEnv on all workers (for gsl and other functions)
                    parallel::clusterEvalQ(parallel.pars$cl, {
                        if (!exists("binSpell", mode = "function", envir = .GlobalEnv, inherits = FALSE)) {
                            if ("climate4R.agro" %in% loadedNamespaces()) {
                                assign("binSpell", 
                                       utils::getFromNamespace("binSpell", ns = "climate4R.agro"),
                                       envir = .GlobalEnv)
                            }
                        }
                    })
                    verify_funs_list <- c(verify_funs_list, "agroindexFAO_tier1", index.code, "binSpell", "yearStartEnd")
                }
                
                if (is_fao_agro) {
                    cluster_assign_function(parallel.pars$cl, "agroindexFAO", fao_fun_loaded)
                    cluster_assign_function(parallel.pars$cl, "computeET0", computeET0_loaded)
                    verify_funs_list <- c(verify_funs_list, "agroindexFAO", "computeET0")
                    # binSpell only exported once if needed by both Tier1 and FAO
                    if (!is.null(binSpell_loaded) && !("binSpell" %in% verify_funs_list)) {
                        cluster_assign_function(parallel.pars$cl, "binSpell", binSpell_loaded)
                        verify_funs_list <- c(verify_funs_list, "binSpell")
                    }
                }
                
                if (is_cdi_cei) {
                    cluster_assign_function(parallel.pars$cl, index.code, cdi_cei_fun_loaded)
                    cluster_assign_function(parallel.pars$cl, "build_seasons_map", build_seasons_map_loaded)
                    verify_funs_list <- c(verify_funs_list, index.code, "build_seasons_map")
                }
                
                # Verify all functions at once (more efficient)
                if (length(verify_funs_list) > 0) {
                    cluster_verify_functions(parallel.pars$cl, verify_funs_list)
                }
            }
        }
    } else {
        if (isTRUE(parallel)) message("NOTE: Parallel processing was skipped (unable to parallelize one single member)")
        apply_fun <- lapply
        member_parallel_active <- FALSE
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
            apply_fun(1:n.mem, function(x) {
            if (isTRUE(verbose)) {
                message("[", Sys.time(), "] Member ", x, ": starting ", index.code)
                if (x == 1 && index.code %in% c("CDI", "CEI")) {
                    debug_lines <- capture.output(str(index.arg.list))
                    message("[", Sys.time(), "] index.arg.list:\n", paste(debug_lines, collapse = "\n"))
                }
            }
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
            # Cache Tier1 index function lookup once per member (avoids repeated get() calls in loops)
            tier1_index_fun_cached <- NULL
            if (is_tier1 && n.mem > 1 && member_parallel_active) {
                # Try multiple methods to find the function in parallel workers
                tier1_index_fun_cached <- tryCatch({
                    # First try: get from global environment (where cluster_assign_function exports it)
                    get(index.code, envir = .GlobalEnv, mode = "function", inherits = FALSE)
                }, error = function(e1) {
                    tryCatch({
                        # Second try: get with inheritance (might be in parent env)
                        get(index.code, mode = "function", inherits = TRUE)
                    }, error = function(e2) {
                        tryCatch({
                            # Third try: get from namespace
                            if ("climate4R.agro" %in% loadedNamespaces()) {
                                utils::getFromNamespace(index.code, ns = "climate4R.agro")
                            } else {
                                NULL
                            }
                        }, error = function(e3) {
                            NULL
                        })
                    })
                })
                # If still NULL, this is a problem - but we'll handle it in the call site
                if (is.null(tier1_index_fun_cached)) {
                    # Try one more time with explicit environment search
                    tier1_index_fun_cached <- tryCatch({
                        eval(parse(text = index.code), envir = .GlobalEnv)
                    }, error = function(e) NULL)
                }
            }
            
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
                    for (i in seq_along(dates_list)[-1]) {
                        if (!identical(first_dates, dates_list[[i]])) {
                            dates_differ <- TRUE
                            break
                        }
                    }
                }
                
                if (dates_differ) {
                    tryCatch({
                        grid.list.aux <- intersect_grid(grid.list.aux, type = "temporal", which.return = seq_along(grid.list.aux))
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
                        for (i in seq_along(grid.list.aux)) {
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
            
            fao_out <- local({
                base_grid <- grid.list.aux[[1]]
                base_dates <- as.Date(base_grid[["Dates"]][["start"]])
                dates_mat <- cbind(
                    as.numeric(format(base_dates, "%Y")),
                    as.numeric(format(base_dates, "%m")),
                    as.numeric(format(base_dates, "%d"))
                )
                years_unique <- sort(unique(as.integer(format(base_dates, "%Y"))))
                n_years <- length(years_unique)
                if (n_years == 0) {
                    stop("No yearly data available for FAO agronomic index.")
                }

            input.arg.list <- lapply(grid.list.aux, function(d) {
                if (is.null(d)) return(NULL)
                data_arr <- d[["Data"]]
                    dims <- dim(data_arr)
                    dim_names <- attr(data_arr, "dimensions")
                    if (length(dims) == 4 && !is.null(dim_names) && "member" %in% dim_names && dims[which(dim_names == "member")] == 1) {
                        member_dim_idx <- which(dim_names == "member")
                        if (member_dim_idx == 1) {
                                data_arr <- data_arr[1, , , , drop = TRUE]
                        } else if (member_dim_idx == 2) {
                                data_arr <- data_arr[, 1, , , drop = TRUE]
                        } else if (member_dim_idx == 3) {
                                data_arr <- data_arr[, , 1, , drop = TRUE]
                        } else {
                                data_arr <- data_arr[, , , 1, drop = TRUE]
                            }
                    } else if (length(dims) == 4 && dims[1] == 1) {
                        data_arr <- array(data_arr, dim = dims[-1])
                    } else if (length(dims) < 2) {
                        stop("Unexpected data dimensions for FAO agronomic processing.")
                    }
                    data_arr
                })
                dims_list <- lapply(input.arg.list, function(arr) if (!is.null(arr)) dim(arr) else NULL)
            dims_list <- dims_list[!sapply(dims_list, is.null)]
            if (length(dims_list) > 1) {
                    ref_dims <- dims_list[[1]]
                for (i in seq_along(dims_list)[-1]) {
                        if (!identical(dims_list[[i]], ref_dims)) {
                            stop("Data arrays have inconsistent dimensions. Expected ",
                                 paste(ref_dims, collapse = " x "), " but got ",
                                 paste(dims_list[[i]], collapse = " x "))
                        }
                    }
                }

                lats <- base_grid[["xyCoords"]][["y"]]
                lons_vec <- tryCatch(base_grid[["xyCoords"]][["x"]], error = function(...) NULL)
                n_lat <- length(lats)
                n_lon <- if (!is.null(lons_vec)) length(lons_vec) else get_shape(base_grid)[["lon"]]
                if (is.null(n_lon) || is.na(n_lon)) n_lon <- dim(base_grid[["Data"]])[3]
                if (is.null(n_lat) || is.null(n_lon) || n_lat == 0 || n_lon == 0) {
                    stop("Unable to determine spatial dimensions for FAO agronomic index.")
                }

                fao_args_static <- index.arg.list
                fao_args_static[c("dates", "lat", "index.code")] <- NULL
                
                result_array <- array(NA_real_, dim = c(n_years, n_lat, n_lon))
                error_counter <- 0L
                first_error_msg <- NULL
                first_all_na <- NULL
                
                point_base_args <- c(
                    list(dates = dates_mat, index.code = index.code),
                    fao_args_static
                )
                
                fao_runner <- if (n.mem > 1 && member_parallel_active) {
                    # In parallel mode, use the exported function name directly to avoid serialization issues
                    function(args_point) {
                        withCallingHandlers(
                            suppressMessages(suppressWarnings(
                                do.call("agroindexFAO", args_point)
                            )),
                            warning = function(w) invokeRestart("muffleWarning"),
                            message = function(m) invokeRestart("muffleMessage")
                        )
                    }
                } else {
                    # In serial mode, use the local variable
                    function(args_point) {
                        do.call(fao_fun_loaded, args_point)
                    }
                }
                
                # Use lapply instead of nested for loops for better performance
                # Iterate over latitudes
                latloop_results <- lapply(seq_len(n_lat), function(lat_idx) {
                    lat_val <- lats[lat_idx]
                    # Iterate over longitudes
                    lonloop_results <- lapply(seq_len(n_lon), function(lon_idx) {
                        args_point <- point_base_args
                        args_point$lat <- lat_val
                        if (!is.null(input.arg.list$pr)) args_point$pr <- input.arg.list$pr[, lat_idx, lon_idx]
                        if (!is.null(input.arg.list$tx)) args_point$tx <- input.arg.list$tx[, lat_idx, lon_idx]
                        if (!is.null(input.arg.list$tn)) args_point$tn <- input.arg.list$tn[, lat_idx, lon_idx]
                        if (!is.null(input.arg.list$tm)) args_point$tm <- input.arg.list$tm[, lat_idx, lon_idx]

                        point_result <- tryCatch({
                            fao_runner(args_point)
                        }, error = function(e) {
                            # Return error info along with NA result
                            lon_val <- if (!is.null(lons_vec) && length(lons_vec) >= lon_idx) lons_vec[lon_idx] else lon_idx
                            return(list(result = rep(NA_real_, n_years), 
                                       error = paste0("lat=", lat_val, ", lon=", lon_val, ": ", e$message),
                                       is_error = TRUE))
                        })

                        # Handle error return structure
                        is_error <- is.list(point_result) && !is.null(point_result$is_error) && point_result$is_error
                        if (is_error) {
                            error_info <- point_result
                            point_result <- error_info$result
                        } else {
                            error_info <- NULL
                        }

                        if (is.list(point_result) && index.code == "gsl" && "GSL" %in% names(point_result)) {
                            point_result <- point_result$GSL
                        }

                        if (length(point_result) != n_years) {
                            if (length(point_result) == 1 && is.na(point_result)) {
                                point_result <- rep(NA_real_, n_years)
                            } else {
                                lon_val <- if (!is.null(lons_vec) && length(lons_vec) >= lon_idx) lons_vec[lon_idx] else lon_idx
                                error_msg <- paste0("lat=", lat_val, ", lon=", lon_val,
                                                   ": unexpected result length ", length(point_result),
                                                   " (expected ", n_years, ")")
                                point_result <- rep(NA_real_, n_years)
                                return(list(result = point_result, error = error_msg, is_error = TRUE))
                            }
                        }

                        # Check for all NA and collect info
                        all_na_info <- NULL
                        if (all(is.na(point_result))) {
                            lon_val <- if (!is.null(lons_vec) && length(lons_vec) >= lon_idx) lons_vec[lon_idx] else lon_idx
                            na_counts <- lapply(args_point[c("pr", "tx", "tn", "tm")], function(vec) {
                                if (is.null(vec)) return(NA_integer_)
                                sum(is.na(vec))
                            })
                            all_na_info <- list(lat = lat_val, lon = lon_val, na = na_counts)
                        }

                        # Return result with metadata
                        list(result = point_result, 
                             error = if (!is.null(error_info)) error_info$error else NULL,
                             all_na = all_na_info,
                             is_error = is_error)
                    })
                    # Return list of results for this latitude
                    lonloop_results
                })
                
                # Process results and populate result_array, track errors
                for (lat_idx in seq_len(n_lat)) {
                    for (lon_idx in seq_len(n_lon)) {
                        res <- latloop_results[[lat_idx]][[lon_idx]]
                        result_array[, lat_idx, lon_idx] <- res$result
                        
                        if (!is.null(res$error)) {
                            error_counter <- error_counter + 1L
                            if (is.null(first_error_msg)) {
                                first_error_msg <- res$error
                            }
                        }
                        
                        if (!is.null(res$all_na) && is.null(first_all_na)) {
                            first_all_na <- res$all_na
                        }
                    }
                }

                if (error_counter > 0 && n.mem == 1) {
                    tryCatch({
                        message("[", Sys.time(), "] FAO agronomic index (member ", x, "): ",
                                error_counter, " grid point issue(s)",
                                if (!is.null(first_error_msg)) paste0(". First: ", first_error_msg) else "")
                        if (!is.null(first_all_na)) {
                            message("  Example NA location lat=", first_all_na$lat, ", lon=", first_all_na$lon)
                        }
                    }, error = function(e) {})
                } else if (error_counter > 0 && n.mem > 1) {
                    tryCatch({
                        message("[", Sys.time(), "] FAO agronomic index (member ", x, "): first error -> ",
                                if (is.null(first_error_msg)) "(missing message)" else first_error_msg)
                    }, error = function(e) {})
                }
            
                out_grid <- base_grid
                out_grid[["Data"]] <- unname(result_array)
                attr(out_grid[["Data"]], "dimensions") <- c("time", "lat", "lon")
                years_str <- sprintf("%04d", years_unique)
                out_grid[["Dates"]][["start"]] <- as.Date(paste0(years_str, "-01-01"))
                out_grid[["Dates"]][["end"]] <- as.Date(paste0(years_str, "-12-31"))
                if (length(out_grid$xyCoords$y) != n_lat) out_grid$xyCoords$y <- lats
                if (!is.null(lons_vec) && length(out_grid$xyCoords$x) != length(lons_vec)) {
                    out_grid$xyCoords$x <- lons_vec
                }

                out_grid
            })
            return(fao_out)
        }
        
        # NON-FAO INDICES (CDI, CEI, and Tier1 indices)
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
        first_error_msg <- NULL
    first_error_location <- NULL
    first_error_call <- NULL
    first_error_trace <- NULL
        
        # Iterate over latitudes
        latloop <- lapply(seq_along(lats), function(l) {
            # Iterate over longitudes
            lonloop <- lapply(seq_along(lons), function(lo) {
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
                        
                        # In parallel mode, use the exported function name directly to avoid serialization issues
                        # In serial mode, use the local variable
                        cdi_cei_fun <- if (n.mem > 1 && member_parallel_active) {
                            index.code  # Use the exported function name (CDI or CEI)
                        } else {
                            cdi_cei_fun_loaded  # Use local variable in serial mode
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
                                compute_season_max <- function(column_name) {
                                    seasons <- cdi_cei_result$season_id
                                    values_vec <- cdi_cei_result[[column_name]]
                                    if (is.null(seasons) || length(seasons) == 0) {
                                        return(numeric(0))
                                    }
                                    unique_seasons <- unique(seasons)
                                    vapply(
                                        unique_seasons,
                                        function(season_val) {
                                            season_values <- values_vec[seasons == season_val]
                                            if (length(season_values) == 0 || all(is.na(season_values))) {
                                                return(NA_real_)
                                            }
                                            max(season_values, na.rm = TRUE)
                                        },
                                        numeric(1)
                                    )
                                }
                                if (index.code == "CDI") {
                                    result <- compute_season_max("cum_days")
                                } else {
                                    # CEI
                                    result <- compute_season_max("cum_excess")
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
                        if (is_tier1) {
                            # In parallel mode, use cached function to avoid repeated lookups
                            if (n.mem > 1 && member_parallel_active) {
                                if (!is.null(tier1_index_fun_cached) && is.function(tier1_index_fun_cached)) {
                                    index_result <- do.call(tier1_index_fun_cached, args_list)
                                } else {
                                    # Fallback: try to use agroindexFAO_tier1 if direct function not found
                                    # This should work since we exported it to workers
                                    args_with_code <- c(list(index.code = index.code), args_list)
                                    index_result <- tryCatch({
                                        do.call("agroindexFAO_tier1", args_with_code)
                                    }, error = function(e) {
                                        # Last resort: try to get function again
                                        fun <- tryCatch({
                                            get(index.code, mode = "function", inherits = TRUE)
                                        }, error = function(e2) {
                                            if ("climate4R.agro" %in% loadedNamespaces()) {
                                                utils::getFromNamespace(index.code, ns = "climate4R.agro")
                                            } else {
                                                stop("Function ", index.code, " not found in parallel worker: ", e$message)
                                            }
                                        })
                                        do.call(fun, args_list)
                                    })
                                }
                            } else {
                                # In serial mode, use agroindexFAO_tier1 wrapper (works fine)
                                args_with_code <- c(list(index.code = index.code), args_list)
                                index_result <- do.call("agroindexFAO_tier1", args_with_code)
                            }
                        } else {
                            index_result <- do.call(metadata$indexfun, args_list)
                        }
                        
                        # GSL returns a list, extract the GSL component
                        if (index.code == "gsl") {
                            # Ensure index_result is a list with GSL component
                            if (!is.list(index_result)) {
                                # If not a list, it's an error - return NAs for all years
                                gsl_result <- rep(NA_real_, length(ref_years))
                            } else if (!"GSL" %in% names(index_result)) {
                                # If list but no GSL component, return NAs for all years
                                gsl_result <- rep(NA_real_, length(ref_years))
                            } else {
                                gsl_result <- index_result$GSL
                                # Ensure result is a vector with correct length (one value per year)
                                if (is.null(gsl_result)) {
                                    # If GSL is NULL, return NAs for all years
                                    gsl_result <- rep(NA_real_, length(ref_years))
                                } else if (!is.vector(gsl_result) || !is.numeric(gsl_result)) {
                                    # If GSL is not a numeric vector, return NAs
                                    gsl_result <- rep(NA_real_, length(ref_years))
                                } else if (length(gsl_result) != length(ref_years)) {
                                    # If length doesn't match, try to align with ref_years
                                    # gsl function returns values for years in the data, which might differ
                                    if ("year" %in% names(index_result) && !is.null(index_result$year)) {
                                        # Map GSL values to ref_years
                                        gsl_years <- index_result$year
                                        # Convert gsl_years to character for matching (ref_years is character)
                                        gsl_years_char <- as.character(gsl_years)
                                        aligned_result <- rep(NA_real_, length(ref_years))
                                        # Match ref_years (character) to gsl_years_char
                                        year_match <- match(ref_years, gsl_years_char)
                                        valid_match <- !is.na(year_match) & year_match <= length(gsl_result)
                                        if (any(valid_match)) {
                                            aligned_result[valid_match] <- gsl_result[year_match[valid_match]]
                                        }
                                        gsl_result <- aligned_result
                                    } else {
                                        # Can't align, return NAs
                                        gsl_result <- rep(NA_real_, length(ref_years))
                                    }
                                }
                            }
                            gsl_result
                        } else {
                            # For other Tier1 indices, ensure result is a vector
                            if (is.vector(index_result) && length(index_result) != length(ref_years)) {
                                # Try to align if possible, otherwise return NAs
                                if (length(index_result) > 0 && !all(is.na(index_result))) {
                                    # If result is shorter, pad with NAs; if longer, truncate
                                    if (length(index_result) < length(ref_years)) {
                                        result_padded <- rep(NA_real_, length(ref_years))
                                        result_padded[seq_along(index_result)] <- index_result
                                        index_result <- result_padded
                                    } else {
                                        index_result <- index_result[seq_len(length(ref_years))]
                                    }
                                } else {
                                    index_result <- rep(NA_real_, length(ref_years))
                                }
                            }
                            index_result
                        }
                    }
                }, error = function(e) {
                    error_count <<- error_count + 1
                    error_msg <- conditionMessage(e)
                    error_location <- paste0("lat=", lats[l], ", lon=", lons[lo])
                    error_call <- conditionCall(e)
                    if (is.null(error_call)) {
                        error_call_str <- "NULL"
                    } else {
                        error_call_str <- paste(deparse(error_call), collapse = "")
                    }
                    if (is.null(first_error_trace)) {
                        first_error_trace <<- capture.output({
                            traceback_calls <- sys.calls()
                            print(traceback_calls)
                        })
                    }
                    if (is.null(first_error_call)) {
                        first_error_call <<- error_call_str
                    }
                    
                    # Store first error for reporting
                    if (is.null(first_error_msg)) {
                        first_error_msg <<- error_msg
                        first_error_location <<- error_location
                    }
                    
                    if (error_display_count < MAX_ERROR_DISPLAY) {
                        # Show error messages - try to display even in parallel mode (with safe handling)
                        if (n.mem == 1) {
                            tryCatch({
                                warning("Error at ", error_location, ": ", error_msg)
                            }, error = function(e2) {
                                # Ignore connection write errors during warning
                            })
                        } else {
                            # In parallel mode, try to log first few errors (may fail due to connection issues)
                            if (error_display_count == 0) {
                                tryCatch({
                                    message("[", Sys.time(), "] Member ", x, " - First error at ", error_location, ": ", error_msg)
                                }, error = function(e2) {
                                    # Ignore connection write errors in parallel workers
                                })
                            }
                        }
                        error_display_count <<- error_display_count + 1
                    }
                    # Return NAs with correct length for yearly data
                    rep(NA_real_, length(ref_years))
                })
                
                return(result)
            })
            do.call(abind::abind, list(lonloop, along = 0))
        })
        
        # Report error summary for non-FAO indices
        if (error_count > 0) {
            total_points <- length(lats) * length(lons)
            err_msg <- paste0("Non-FAO indices: ", error_count, " out of ", total_points, " grid point(s) failed.")
            if (!is.null(first_error_location) && !is.null(first_error_msg)) {
                err_msg <- paste0(err_msg, " First error at ", first_error_location, ": ", first_error_msg)
            } else if (!is.null(first_error_msg)) {
                err_msg <- paste0(err_msg, " First error: ", first_error_msg)
            } else if (error_display_count > 0) {
                err_msg <- paste0(err_msg, " Check earlier warnings for details.")
            }
            if (!is.null(first_error_call)) {
                err_msg <- paste0(err_msg, "\nCall: ", first_error_call)
            }
            if (!is.null(first_error_trace)) {
                err_msg <- paste0(
                    err_msg,
                    "\nTraceback snapshot:\n",
                    paste(first_error_trace, collapse = "\n")
                )
            }
            if (n.mem > 1) {
                stop(err_msg)
            } else {
                tryCatch({
                    message("[", Sys.time(), "] ", err_msg)
                }, error = function(e) {})
            }
        }
        
        # Reconstruct data array with proper dimensions
        latloop_bound <- do.call(abind::abind, list(latloop, along = 0))
        out_array <- unname(aperm(latloop_bound, c(3, 1, 2)))
        
        # Update dates to match the output time dimension (years for FAO indices)
        n_output_times <- dim(out_array)[1]  # First dimension is time
        # Tier1 indices (including gsl) always return yearly values, even if n_output_times == length(refDates)
        # This can happen if there's a bug or if the data happens to have the same number of days as years
        if (is_tier1) {
            # FAO tier1 indices always return yearly values
            # Extract unique years and create year boundaries
            years <- ref_years
            if (n_output_times == length(years)) {
                # Create start (Jan 1) and end (Dec 31) for each year
                start_dates <- as.Date(paste0(years, "-01-01"))
                end_dates <- as.Date(paste0(years, "-12-31"))
                if (n.mem == 1 && n_output_times != length(refDates)) {
                    tryCatch({
                        message("[", Sys.time(), "] Converting daily dates to yearly dates (", n_output_times, " years)")
                    }, error = function(e) {
                        # Ignore connection write errors
                    })
                }
            } else {
                # Output doesn't match expected number of years - might be an error
                # Still try to create yearly dates if possible
                if (n_output_times <= length(years)) {
                    # Take first n_output_times years
                    years_subset <- years[seq_len(n_output_times)]
                    start_dates <- as.Date(paste0(years_subset, "-01-01"))
                    end_dates <- as.Date(paste0(years_subset, "-12-31"))
                } else {
                    # More output times than years - take evenly spaced dates as fallback
                    indices <- round(seq(1, length(refDates), length.out = n_output_times))
                    start_dates <- ref_dates_date[indices]
                    end_dates <- start_dates
                }
                if (n.mem == 1) {
                    tryCatch({
                        message("[", Sys.time(), "] Warning: Output time steps (", n_output_times, 
                               ") don't match expected years (", length(years), ")")
                    }, error = function(e) {
                        # Ignore connection write errors
                    })
                }
            }
        } else if (n_output_times != length(refDates)) {
            # For FAO agronomic indices (non-tier1), check if output matches years
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
        refGrid_member <- get_member_template_grid(x, n.mem,
                                                   tx, tn, pr, tm, hurs, sfcwind, ssrd)
        
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
            error_msg <- conditionMessage(e)
            if (n.mem == 1) {
                tryCatch({
                    warning("Error processing member ", x, ": ", error_msg)
                }, error = function(e2) {
                    # Ignore connection write errors during warning
                })
            }
            # Create a grid with NA values but same structure as other members
            # Use helper function to simplify code
            refGrid_member <- get_member_template_grid(x, n.mem,
                                                       tx, tn, pr, tm, hurs, sfcwind, ssrd)

            if (is.null(refGrid_member)) {
                # Fallback to refGrid from outer scope
                refGrid_member <- if (n.mem > 1) {
                    safe_redim_drop(subset_grid(refGrid, members = x, drop = TRUE), drop_member = TRUE, drop_loc = TRUE)
                } else {
                    safe_redim_drop(refGrid, drop_member = TRUE, drop_loc = TRUE)
                }
            }
            
            # Fill with NA and set appropriate dates
            if (!is.null(refGrid_member) && !is.null(refGrid_member[["Data"]])) {
                refGrid_member[["Data"]][] <- NA
                refGrid_member[["Dates"]][["start"]] <- ref_dates_date
                refGrid_member[["Dates"]][["end"]] <- ref_dates_date
                attr(refGrid_member, ".agroindex_error") <- list(
                    member = x,
                    stage = "member_computation",
                    message = error_msg
                )
            }
            refGrid_member  # Return error grid
        })
        }, error = function(e) {
            # Outer error handler: catch all errors and return a valid grid structure
            # This ensures we always return something that can be combined, even if it's all NAs
            error_msg <- if (is.character(e$message)) e$message else as.character(e)
            
            # Check if this is a connection/SIGPIPE error
            is_connection_error <- (
                grepl("connection", error_msg, ignore.case = TRUE) ||
                grepl("conexion", error_msg, ignore.case = TRUE) ||
                grepl("conex", error_msg, ignore.case = TRUE) ||
                grepl("SIGPIPE", error_msg, ignore.case = TRUE) ||
                grepl("broken pipe", error_msg, ignore.case = TRUE) ||
                grepl("ignoring SIGPIPE", error_msg, ignore.case = TRUE)
            )
            
            # For multi-member grids, always try to return a valid grid structure (even if all NAs)
            # This allows the computation to continue even if one worker has issues
            # SIGPIPE errors are common in parallel processing and should be handled gracefully
            if (member_parallel_active) {
                tryCatch({
                    refGrid_member <- get_member_template_grid(x, n.mem,
                                                               tx, tn, pr, tm, hurs, sfcwind, ssrd)
                    
                    if (is.null(refGrid_member)) {
                        # Fallback: try to use refGrid from outer scope
                        tryCatch({
                            refGrid_member <- if (n.mem > 1) {
                                safe_redim_drop(subset_grid(refGrid, members = x, drop = TRUE), drop_member = TRUE, drop_loc = TRUE)
                            } else {
                                safe_redim_drop(refGrid, drop_member = TRUE, drop_loc = TRUE)
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
                        attr(refGrid_member, ".agroindex_error") <- list(
                            member = x,
                            stage = if (is_connection_error) "connection_recovery" else "member_computation",
                            message = error_msg
                        )
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
                # Re-throw to expose underlying error when not parallelizing
                stop(e)
            }
        })
        
        # Return the result from the inner computation
        return(result)
            })
    }, warning = function(w) {
        # Catch SIGPIPE warnings specifically in withCallingHandlers
        warn_msg <- if (is.character(w$message)) {
            w$message
        } else {
            as.character(w)
        }
        if (grepl("SIGPIPE", warn_msg, ignore.case = TRUE) && n.mem > 1) {
            # Suppress SIGPIPE warnings in parallel mode
            invokeRestart("muffleWarning")
        }
    }, message = function(m) {
        # Catch SIGPIPE messages (sometimes printed as messages, not warnings)
        msg_text <- if (is.character(m$message)) {
            m$message
        } else {
            as.character(m)
        }
        if (grepl("SIGPIPE", msg_text, ignore.case = TRUE) && n.mem > 1) {
            # Suppress SIGPIPE messages in parallel mode
            invokeRestart("muffleMessage")
        }
    })
})
    
    # Collect error information returned by workers (if any)
    if (is.list(out.list)) {
        for (idx in seq_along(out.list)) {
            err_attr <- attr(out.list[[idx]], ".agroindex_error", exact = TRUE)
            if (!is.null(err_attr)) {
                if (is.null(err_attr$member)) {
                    err_attr$member <- idx
                }
                err_attr$return_index <- idx
                member_error_summary[[length(member_error_summary) + 1]] <- err_attr
                attr(out.list[[idx]], ".agroindex_error") <- NULL
            }
        }
    } else {
        err_attr <- attr(out.list, ".agroindex_error", exact = TRUE)
        if (!is.null(err_attr)) {
            if (is.null(err_attr$member)) {
                err_attr$member <- 1
            }
            err_attr$return_index <- 1
            member_error_summary[[length(member_error_summary) + 1]] <- err_attr
            attr(out.list, ".agroindex_error") <- NULL
        }
    }
    
    if (length(member_error_summary) > 0) {
        member_error_names <- vapply(member_error_summary, function(err) {
            member_val <- err$member
            if (is.null(member_val) || length(member_val) == 0 || is.na(member_val)) {
                return("member_unknown")
            }
            paste0("member_", member_val)
        }, character(1))
        names(member_error_summary) <- member_error_names
    }
    
    connection_patterns <- c("SIGPIPE", "broken pipe", "connection", "conexion", "conex", "error al escribir en una conexion")
    connection_only <- length(member_error_summary) > 0 &&
        all(vapply(member_error_summary, function(err) {
            msg <- err$message
            if (is.null(msg)) return(FALSE)
            any(vapply(connection_patterns, function(pat) grepl(pat, msg, ignore.case = TRUE), logical(1)))
        }, logical(1)))
    
    if (connection_only && !isTRUE(._retry)) {
        tryCatch({
            message("[", Sys.time(), "] Worker connection errors detected; retrying sequentially to expose underlying error.")
        }, error = function(e) {
            return(NULL)
        })
        original_call$parallel <- FALSE
        original_call$ncores <- 1
        original_call[["._retry"]] <- TRUE
        return(eval(original_call, envir = parent.frame()))
    }
    
    # Suppress message output to avoid connection errors in batch jobs
    tryCatch({
        message("[", Sys.time(), "] Done")
    }, error = function(e) {
        # Ignore connection write errors during message
    })
    
    out <- combine_member_results(out.list, station, metadata, index.code)
    if (length(member_error_summary) > 0) {
        first_err <- member_error_summary[[1]]
        member_label <- if (!is.null(first_err$member)) first_err$member else first_err$return_index
        err_msg <- if (!is.null(first_err$message)) first_err$message else "(missing error message)"
        issues_count <- length(member_error_summary)
        issues_suffix <- if (issues_count == 1) " issue" else " issues"
        message_text <- paste0(
            "[", Sys.time(), "] Worker issue detected (member ", member_label, "): ", err_msg,
            "\nInspect attr(result, \"member_errors\") for details on ", issues_count, issues_suffix, "."
        )
        tryCatch({
            message(message_text)
        }, error = function(e) {
            # Ignore connection write errors during message
        })
        attr(out, "member_errors") <- member_error_summary
    }
    invisible(out)
}                          


#' @title List all available Agroclimatic Indices
#' @description Print a table with a summary of the available agroclimatic indices including FAO tier1 indices and stress indices
#' @return Print a table on the screen with the following columns:
#' \itemize{
#' \item \strong{code}: Code of the index. This is the character string used as input value for the argument \code{index.code} in \code{\link{agroindexGrid}}
#' \item \strong{longname}: Long description of the index
#' \item \strong{index.fun}: The name of the internal function used to calculate it
#' \item \strong{tn,tx,tm,pr,hurs,sfcwind,ssrd}: Logical values (0/1) indicating the input variables required for index calculation
#' \item \strong{units}: The units of the index (when different from those of the input variable)
#' }
#' @references FAO agroclimatic indices documentation
#' @author J. Bedia (original)
#' @export
#'
agroindexShow <- function() {
    read.master()
}



#' @keywords internal
#' @importFrom utils read.table

read.master <- function() {
  # Try installed package first, then development location
  master_file <- system.file("master", package = "climate4R.agro")

  # If package not installed, try development paths
  if (master_file == "" || !file.exists(master_file)) {
    # Common development locations
    dev_paths <- c(
      "C:/Users/pablo/Desktop/climate4R.agro/inst/master",
      "/lustre/gmeteo/WORK/lavinp/test/inst/master",
      file.path(getwd(), "inst", "master"),
      file.path(getwd(), "..", "inst", "master"),
      "inst/master",
      "../inst/master"
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
    stop(
      "Cannot find master file. Tried paths:\n",
      "  - system.file(\"master\", package = \"climate4R.agro\")\n",
      "  - C:/Users/pablo/Desktop/climate4R.agro/inst/master\n",
      "  - inst/master (relative paths)\n",
      "Please ensure climate4R.agro package is installed or master file exists ",
      "in inst/ folder."
    )
  }

  master_dir <- dirname(master_file)
  extra_paths <- unique(na.omit(c(
    getOption("climate4R.agro.extra_paths"),
    master_dir,
    normalizePath(file.path(master_dir, ".."), mustWork = FALSE),
    normalizePath(file.path(master_dir, "..", ".."), mustWork = FALSE)
  )))
  options(climate4R.agro.extra_paths = extra_paths)
  read.table(
    master_file,
    header = TRUE,
    sep = ";",
    stringsAsFactors = FALSE,
    na.strings = ""
  )
}



