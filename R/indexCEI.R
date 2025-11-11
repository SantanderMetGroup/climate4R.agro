#' Condition Excess Index
#'
#' @description
#' Daily accumulator of *excess over lower* for a single variable.
#' For each `id` and season, computes a *daily* cumulative series `cum_excess`
#' from a **single climatic variable** `x`. A day contributes the value
#' `excess = (x - lower)` **only if** `x` lies within the interval
#' `[lower, upper]` (with configurable inclusivity). The contribution is further
#' gated by a **minimum run-length** requirement `min_duration`: a day only
#' counts if it belongs to a TRUE-run whose length has reached at least
#' `min_duration`.
#'
#' @param df Tibble/data.frame with one row per day and per `id`. All `id`s must
#'   have the same number of rows (daily-regular series).
#' @param id Character, name of the `id` column (default `"id"`).
#' @param x_col Character, name of the **single** numeric variable to analyse.
#' @param start_date Date of the first row for each `id` (used to build the
#'   calendar and season mapping).
#' @param season_start,season_end Character `"mm-dd"` giving the season window.
#'   The window may **cross years** (e.g., `"07-01"` -> `"06-30"`).
#'
#' @param lower,upper Numeric bounds defining the valid interval for `x`.
#'   \describe{
#'     \item{`lower`}{**Required** and must be finite. The excess is defined as
#'       `x - lower` (capped at 0 below `lower`).}
#'     \item{`upper`}{Optional (default `Inf`). Use to restrict the valid range
#'       above; values beyond `upper` contribute 0 and break runs depending on
#'       `na_action`.}
#'   }
#' @param inc_lower,inc_upper Logical flags controlling interval inclusivity
#'   (both default `TRUE`):
#'   \describe{
#'     \item{`inc_lower = TRUE`}{`x >= lower` qualifies; `FALSE` means
#'     `x > lower`.}
#'     \item{`inc_upper = TRUE`}{`x <= upper` qualifies; `FALSE` means
#'     `x < upper`.}
#'   }
#'
#' @param min_duration Integer >= 1. Minimum number of **consecutive qualifying
#'   days** required before the run contributes to the cumulative:
#'   \itemize{
#'     \item If `min_duration = 1`, every qualifying day contributes its excess
#'           immediately.
#'     \item If `min_duration = m > 1`, the first day that brings the run length
#'           to `m` contributes its own excess **plus** a one-off **bonus** equal
#'           to the sum of the previous `m - 1` days' excesses in that run (bias
#'           correction). Subsequent days in the run contribute their own excess
#'           normally.
#'   }
#'
#' @param na_action Character controlling how missing values in `x` are handled:
#'   \describe{
#'     \item{`"false"`}{Treat `NA` as out-of-range (i.e., `cond = FALSE`). This
#'       **breaks runs** and contributes 0 excess on those days.}
#'     \item{`"skip_days"`}{Drop days with `NA` in `x` **before** evaluating
#'       conditions. Consecutivity is then defined over the remaining observed
#'       days.}
#'   }
#'
#' @return A tibble with daily rows and columns:
#' \itemize{
#'   \item `id` - identifier (station/cell).
#'   \item `day_idx` - 1..n index within the full time series per `id`.
#'   \item `season_id`, `season_label` - season identifiers.
#'   \item `cond` - logical; `TRUE` when `x` lies within the specified
#'          interval.
#'   \item `valid` - logical; `TRUE` once the current TRUE-run reaches
#'          `min_duration`.
#'   \item `excess` - numeric; `pmax(x - lower, 0)` on qualifying days,
#'          else 0.
#'   \item `cum_excess` - season-wise cumulative of excess with per-run
#'          correction.
#' }
#' 
#' @author Sara Herrera
#'
#' @details
#' **Units.** `excess` and `cum_excess` are in the same units as `x`.
#' **Capping.** By default, `excess` is not capped at `(upper - lower)`; if you
#' want to cap it, post-process with `excess <- pmin(excess, upper - lower)`
#' before accumulation.
#' **Semantics.** The online bias-corrected cumulative equals the retrospective
#' (full-season) total for runs that have already completed, while remaining
#' conservative for runs that are still in progress at the season boundary.
#'
#' @examples
#' set.seed(42)
#' df <- tibble::tibble(
#'   id = rep(c("A", "B"), each = 10),
#'   date = as.Date("2000-01-01") + rep(0:9, times = 2),
#'   temp = runif(20, 20, 40)
#' )
#'
#' # Accumulate excess temperature above 25 degC whenever 25 <= temp (no upper cap),
#' # with runs of at least 3 consecutive qualifying days.
#' CEI(
#'   df,
#'   id = "id", x_col = "temp",
#'   start_date = min(df$date),
#'   lower = 25, upper = Inf,
#'   min_duration = 3
#' )
#'
#' # Same but only count when 25 <= temp <= 35, exclusive upper bound
#' CEI(
#'   df,
#'   id = "id", x_col = "temp",
#'   start_date = min(df$date),
#'   lower = 25, upper = 35,
#'   inc_upper = FALSE,
#'   min_duration = 2
#' )
#'
#' @importFrom dplyr group_by mutate ungroup filter select arrange left_join if_else lag row_number n
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
CEI <- function(
    df,
    id = "id",
    x_col,
    start_date,
    season_start = "07-01",
    season_end = "06-30",
    lower,
    upper = Inf,
    inc_lower = TRUE,
    inc_upper = TRUE,
    min_duration = 1L,
    na_action = c("false", "skip_days")[1]) {
  
  # ---- 0) Validation & setup
  
  stopifnot(is.finite(lower), min_duration >= 1L)
  na_action <- match.arg(na_action, choices = c("false", "skip_days"))
  stopifnot(id %in% names(df), x_col %in% names(df))
  
  df <- tibble::as_tibble(df) %>%
    dplyr::group_by(.data[[as.character(id)]]) %>%
    dplyr::mutate(.n = dplyr::n()) %>%
    dplyr::ungroup()
  nd <- df$.n[1]
  if (any(df$.n != nd)) stop("All groups (id) must share the same number of rows.")
  df <- dplyr::select(df, -".n")
  
  SM <- build_seasons_map(nd, start_date, season_start, season_end)
  S <- SM$seasons
  sid <- SM$season_id
  if (all(is.na(sid))) stop("No days fall within any season window.")
  
  df <- df %>%
    dplyr::group_by(.data[[as.character(id)]]) %>%
    dplyr::mutate(
      day_idx   = dplyr::row_number(),
      season_id = sid[day_idx]
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(season_id))
  
  
  # ---- 1) In-range condition and per-day excess
  
  if (na_action == "skip_days") df <- dplyr::filter(df, !is.na(.data[[as.character(x_col)]]))
  
  x <- df[[x_col]]
  lo_ok <- if (inc_lower) x >= lower else x > lower
  hi_ok <- if (inc_upper) x <= upper else x < upper
  cond <- lo_ok & hi_ok
  if (na_action == "false") cond[is.na(cond)] <- FALSE
  
  df <- df %>%
    dplyr::mutate(
      cond   = cond,
      excess = dplyr::if_else(cond, pmax(.data[[as.character(x_col)]] - lower, 0), 0)
    )
  
  # ---- 2) Run-length logic (online with per-run correction)
  
  k <- max(min_duration - 1L, 0L)
  
  df <- df %>%
    dplyr::arrange(.data[[as.character(id)]], season_id, day_idx) %>%
    dplyr::group_by(.data[[as.character(id)]], season_id) %>%
    dplyr::mutate(
      # TRUE-run counter via grouping trick (no self-reference)
      grp = cumsum(!cond),
      run_len = dplyr::if_else(cond,
                               as.integer(ave(seq_along(cond), grp, FUN = seq_along)),
                               0L
      ),
      valid = cond & (run_len >= min_duration),
      # cumulative excess within group
      cum_exc = cumsum(excess),
      prev_cum = dplyr::lag(cum_exc, 1, default = 0),
      prev_k = if (k > 0L) dplyr::lag(prev_cum, k, default = 0) else 0,
      # bonus = sum of previous (min_duration - 1) excesses when run hits threshold
      bonus = dplyr::if_else(run_len == min_duration, prev_cum - prev_k, 0),
      cum_excess = cumsum((as.numeric(valid) * excess) + bonus)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(S %>% dplyr::select(season_id, season_label), by = "season_id") %>%
    dplyr::arrange(.data[[as.character(id)]], season_id, day_idx) %>%
    dplyr::select(!!id, day_idx, season_id, season_label, cond, valid, excess, cum_excess)
  
  return(df)
}

