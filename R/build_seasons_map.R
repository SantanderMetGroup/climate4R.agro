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

  # Does the season cross calendar year? (e.g., Jul->Jun)
  ms <- as.integer(substr(season_start, 1, 2))
  ds <- as.integer(substr(season_start, 4, 5))
  me <- as.integer(substr(season_end, 1, 2))
  de <- as.integer(substr(season_end, 4, 5))
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
  season_id <- ifelse(idx > 0 & dates <= S$end[pmax(idx, 1)], idx, NA_integer_)

  list(seasons = S, season_id = season_id, dates = dates)
}

