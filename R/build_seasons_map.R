#' Build season mapping for a daily-regular series
#'
#' @description
#' Given the number of daily observations and a start date, builds:
#' (i) a table of non-overlapping seasons defined by `season_start` and `season_end`
#' (which may cross calendar years), and (ii) a vector `season_id` mapping each
#' day index (1..n_days) to its season (NA if outside all seasons).
#'
#' @param n_days Integer, length of the daily series.
#' @param start_date Date, date of the first observation.
#' @param season_start Character "mm-dd", season start (default "07-01").
#' @param season_end Character "mm-dd", season end (default "06-30").
#'
#' @return A list with:
#' \itemize{
#'   \item `seasons`: tibble with columns `start`, `end`, `season_id`, `season_label`.
#'   \item `season_id`: integer vector length `n_days`, mapping day index -> season id (or NA).
#'   \item `dates`: Date vector length `n_days`.
#' }
#'
#' @author Sara Herrera
#'
#' @examples
#' SM <- build_seasons_map(730, as.Date("2019-01-01"))
#'
#' @importFrom dplyr filter mutate row_number
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#' @export
build_seasons_map <- function(n_days, start_date,
                              season_start = "07-01", season_end = "06-30") {
  if (!inherits(start_date, "Date")) start_date <- as.Date(start_date)
  dates <- seq.Date(start_date, by = "day", length.out = n_days)

  # Validate and parse season_start and season_end format (must be "mm-dd")
  # Check format: should be 5 characters with "-" in the middle
  validate_season_format <- function(season_str, param_name) {
    if (nchar(season_str) != 5 || substr(season_str, 3, 3) != "-") {
      stop("Invalid format for ", param_name, ": '", season_str, 
           "'. Expected format: 'mm-dd' (e.g., '11-01', '04-30')")
    }
    month <- as.integer(substr(season_str, 1, 2))
    day <- as.integer(substr(season_str, 4, 5))
    if (is.na(month) || is.na(day) || month < 1 || month > 12 || day < 1 || day > 31) {
      stop("Invalid date values in ", param_name, ": '", season_str, 
           "'. Month must be 01-12, day must be 01-31")
    }
    # Check if it looks like dd-mm format (day > 12 suggests wrong format)
    if (month > 12) {
      stop("Invalid format for ", param_name, ": '", season_str, 
           "'. It appears to be in 'dd-mm' format. Expected 'mm-dd' format (e.g., '04-30' not '30-04')")
    }
    return(list(month = month, day = day))
  }
  
  season_start_parsed <- validate_season_format(season_start, "season_start")
  season_end_parsed <- validate_season_format(season_end, "season_end")
  
  # Does the season cross calendar year? (e.g., Jul->Jun)
  ms <- season_start_parsed$month
  ds <- season_start_parsed$day
  me <- season_end_parsed$month
  de <- season_end_parsed$day
  crosses <- (me < ms) || (me == ms && de < ds)

  # Build many candidate seasons around the data range
  y_min <- as.integer(format(min(dates), "%Y")) - 1L
  y_max <- as.integer(format(max(dates), "%Y")) + 1L

  S <- tibble::tibble(
    start = as.Date(paste0(seq.int(y_min, y_max), "-", season_start)),
    end   = as.Date(paste0(seq.int(y_min, y_max) + crosses, "-", season_end))
  ) %>%
    dplyr::filter(start <= end) %>%
    dplyr::mutate(
      season_id    = dplyr::row_number(),
      season_label = paste0(format(start, "%Y"), "-", format(end, "%Y"))
    )

  # Map each date to its season (seasons are ordered and non-overlapping)
  idx <- findInterval(dates, S$start)
  # Ensure idx is within valid range before accessing S$end
  # findInterval returns 0 for dates before first season, and can return length(S$start) 
  # for dates after last season. We need to check that idx is valid (1 to length(S$end))
  # and that the date is within the season end date.
  n_seasons <- length(S$end)
  # idx must be > 0 and <= n_seasons to be valid
  # Also, if idx == n_seasons, we need to check that the date is actually within that season
  valid_idx <- idx > 0 & idx <= n_seasons
  # For valid indices, check if date falls within the season (date <= season end)
  # Use vectorized comparison only for valid indices to avoid out-of-bounds access
  in_season <- rep(FALSE, length(dates))
  if (any(valid_idx)) {
    # Safely access S$end only for valid indices
    season_ends <- S$end[idx[valid_idx]]
    in_season[valid_idx] <- dates[valid_idx] <= season_ends
  }
  season_id <- ifelse(valid_idx & in_season, idx, NA_integer_)

  list(seasons = S, season_id = season_id, dates = dates)
}

