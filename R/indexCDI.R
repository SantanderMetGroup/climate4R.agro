#' Condition Duration Index
#'
#' @description
#' Daily accumulator for multi-variable day counts.
#' For each `id` and season, compute a *daily* cumulative counter `cum_days`
#' of days that satisfy a combined condition across multiple variables, with
#' a minimum-run constraint `min_duration`. The function supports flexible
#' logic to combine conditions across variables and to control how missing
#' values are treated. It uses an **online** algorithm with **per-run bias
#' correction**, ensuring that cumulative values match retrospective counts
#' for completed runs, while remaining conservative for ongoing (incomplete)
#' runs.
#'
#' @param df Tibble/data.frame with one row per day and per `id`. Must have
#'   the same number of rows for each `id`. Variables referenced in `bounds`
#'   must be numeric columns in `df`.
#' @param id Character, name of the id column (default `"id"`).
#' @param start_date Date of the first row for each `id` (daily-regular series).
#' @param season_start,season_end Character `"mm-dd"`, defining the start and
#'   end of the season window. The window may cross calendar years (e.g.
#'   `"07-01"` to `"06-30"`).
#'
#' @param bounds Tibble or data frame defining **per-variable thresholds**.
#'   It must contain one row per variable with at least the column:
#'   \describe{
#'     \item{`var`}{Character: name of the variable in `df` to evaluate.}
#'   }
#'   Optional columns allow specifying the interval and inclusivity criteria
#'   for each variable:
#'   \describe{
#'     \item{`lower`}{Numeric lower bound (default `-Inf`).}
#'     \item{`upper`}{Numeric upper bound (default `Inf`).}
#'     \item{`inc_lower`}{Logical; whether the lower bound is inclusive (default `TRUE`).}
#'     \item{`inc_upper`}{Logical; whether the upper bound is inclusive (default `TRUE`).}
#'   }
#'   For example:
#'   ```
#'   bounds <- tibble::tibble(
#'     var = c("temp","humidity"),
#'     lower = c(10, 40),
#'     upper = c(30, 80),
#'     inc_lower = TRUE,
#'     inc_upper = TRUE
#'   )
#'   ```
#'   A day is considered valid for a given variable if its value lies within
#'   the specified range (taking inclusivity into account).
#'
#' @param combiner Character specifying how to **combine** per-variable conditions
#'   into a single daily condition (`cond`):
#'   \describe{
#'     \item{`"all"`}{All variables must meet their condition for the day to be
#'       considered valid. Equivalent to a logical AND.}
#'     \item{`"any"`}{At least one variable meets its condition. Equivalent to
#'       a logical OR.}
#'     \item{`"k_of_n"`}{At least `k` of the `n` variables meet their condition.
#'       Useful for indices that require simultaneous satisfaction of several,
#'       but not necessarily all, criteria (e.g., "at least 2 of 3 stressors
#'       below threshold"). Requires argument `k`.}
#'   }
#'
#' @param k Integer >= 1, required only when `combiner = "k_of_n"`. It defines
#'   the minimum number of variables that must pass their individual conditions
#'   for a day to be marked as valid (`cond = TRUE`). For example, with
#'   three variables and `k = 2`, a day is valid if at least two of them are
#'   within their specified ranges.
#'
#' @param min_duration Integer >= 1. Minimum number of **consecutive days**
#'   that must satisfy the combined daily condition (`cond = TRUE`) before
#'   the run starts contributing to the cumulative count.
#'   - If `min_duration = 1`, every valid day contributes immediately.
#'   - If `min_duration = 5`, only when a run of 5 consecutive valid days is
#'     reached do those days begin to count, and the function automatically
#'     adds the previous 4 days at that point (bias correction).
#'
#' @param na_action Character, controlling how **missing values** (`NA`) are
#'   handled. Options:
#'   \describe{
#'     \item{`"false"`}{Treats `NA` as *not satisfying* the condition.
#'       Missing values break runs and are equivalent to days that fail the
#'       thresholds.}
#'     \item{`"skip_vars"`}{When combining across variables, ignore `NA`
#'       for that variable only. For `"all"` or `"k_of_n"`, this means a day
#'       can still be valid if enough other variables are available and meet
#'       their conditions.}
#'     \item{`"skip_days"`}{Drops entire days that contain at least one `NA`
#'       among the variables referenced in `bounds`, before evaluating
#'       conditions. Consecutivity is preserved only across remaining (non-NA)
#'       days.}
#'   }
#'
#' @return Tibble with daily rows and columns:
#' \itemize{
#'   \item `id` - station or spatial cell identifier.
#'   \item `day_idx` - daily index within the full time series.
#'   \item `season_id`, `season_label` - season grouping identifiers.
#'   \item `cond` - logical vector of days meeting the combined condition.
#'   \item `valid` - logical vector indicating if the day is part of a run
#'         that reached the required `min_duration`.
#'   \item `cum_days` - cumulative number of valid days since season start,
#'         corrected for incomplete runs.
#' }
#' 
#' @author Sara Herrera
#'
#' @details
#' The algorithm proceeds day by day ("online"), tracking consecutive valid
#' runs and applying a **per-run correction** once a run first reaches
#' `min_duration`. This ensures that the cumulative counts per season match
#' what would be obtained with a complete retrospective analysis, while
#' remaining unbiased for incomplete ongoing runs (e.g., when the current
#' season is still in progress).
#'
#' @examples
#' set.seed(123)
#' df <- tibble::tibble(
#'   id = rep(c("A", "B"), each = 10),
#'   date = as.Date("2000-01-01") + rep(0:9, times = 2),
#'   temp = runif(20, 5, 35),
#'   humidity = runif(20, 20, 90)
#' )
#'
#' bounds <- tibble::tibble(
#'   var = c("temp", "humidity"),
#'   lower = c(10, 40),
#'   upper = c(30, 80)
#' )
#'
#' # Count days where both variables meet their conditions
#' CDI(df,
#'   id = "id", start_date = min(df$date),
#'   bounds = bounds, combiner = "all", min_duration = 3
#' )
#'
#' # Require at least 2 of 3 variables to meet their conditions
#' CDI(df,
#'   id = "id", start_date = min(df$date),
#'   bounds = bounds, combiner = "k_of_n", k = 2
#' )
#'
#' @importFrom dplyr group_by mutate ungroup filter select arrange left_join across all_of if_else if_all n row_number
#' @importFrom tibble as_tibble tibble
#' @importFrom magrittr %>%
#' @importFrom rlang .data :=
#' @importFrom stats ave
#' @export
CDI <- function(
    df,
    id = "id",
    start_date,
    season_start = "07-01",
    season_end = "06-30",
    bounds,
    combiner = c("all", "any", "k_of_n")[1],
    k = NULL,
    min_duration = 1,
    na_action = c("false", "skip_vars", "skip_days")[1]) {
  
  # ---- 0) Validation & setup
  
  stopifnot(min_duration >= 1L)
  combiner <- match.arg(combiner, choices = c("all", "any", "k_of_n"))
  na_action <- match.arg(na_action, choices = c("false", "skip_vars", "skip_days"))

  df <- tibble::as_tibble(df) %>%
    dplyr::group_by(.data[[as.character(id)]]) %>%
    dplyr::mutate(.n = dplyr::n()) %>%
    dplyr::ungroup()
  nd <- df$.n[1]
  if (any(df$.n != nd)) stop("All groups (id) must share the same number of rows.")
  df <- dplyr::select(df, -".n")

  # Seasons (built once using the common length)
  SM <- build_seasons_map(nd, start_date, season_start, season_end)
  S <- SM$seasons
  sid <- SM$season_id
  if (all(is.na(sid))) stop("No days fall within any season window.")

  # Index days per id and attach season_id
  df <- df %>%
    dplyr::group_by(.data[[as.character(id)]]) %>%
    dplyr::mutate(
      day_idx   = dplyr::row_number(),
      season_id = sid[day_idx]
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(season_id))

  
  # ---- 1) Per-variable conditions
  
  B <- tibble::as_tibble(bounds)
  if (!"lower" %in% names(B)) B <- dplyr::mutate(B, lower = -Inf)
  if (!"upper" %in% names(B)) B <- dplyr::mutate(B, upper = Inf)
  if (!"inc_lower" %in% names(B)) B <- dplyr::mutate(B, inc_lower = TRUE)
  if (!"inc_upper" %in% names(B)) B <- dplyr::mutate(B, inc_upper = TRUE)
  if (!all(B$var %in% names(df))) stop("Some bounds$var not found in df.")

  in_interval <- function(v, lo, hi, il, iu) {
    lo_ok <- if (il) v >= lo else v > lo
    hi_ok <- if (iu) v <= hi else v < hi
    lo_ok & hi_ok
  }

  # Optionally drop days with any NA on referenced variables
  if (na_action == "skip_days") {
    df <- df %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(B$var), ~ !is.na(.x)))
  }

  # Build per-variable logical columns cond__var
  for (i in seq_len(nrow(B))) {
    v <- B$var[i]
    df <- df %>%
      dplyr::mutate(
        !!paste0("cond__", v) :=
          in_interval(
            .data[[as.character(v)]], B$lower[i], B$upper[i],
            isTRUE(B$inc_lower[i]), isTRUE(B$inc_upper[i])
          )
      )
  }
  cond_cols <- paste0("cond__", B$var)

  # NA policy per-variable
  if (na_action == "false") {
    df <- df %>%
      dplyr::mutate(dplyr::across(
        dplyr::all_of(cond_cols),
        ~ dplyr::if_else(is.na(.x), FALSE, .x)
      ))
  }

  # ---- 2) Combine variables into a single daily condition

  # Count TRUEs (optionally ignoring NA with na.rm)
  df <- df %>%
    dplyr::mutate(
      n_true = rowSums(
        dplyr::across(dplyr::all_of(cond_cols),
          ~ as.integer(.x %in% TRUE),
          .names = NULL
        ),
        na.rm = (na_action == "skip_vars")
      ),
      n_vars = if (na_action == "skip_vars") {
        rowSums(dplyr::across(dplyr::all_of(cond_cols), ~ !is.na(.x), .names = NULL))
      } else {
        length(cond_cols)
      }
    )

  # Build 'cond' by chosen combiner
  if (combiner == "all") {
    df <- df %>% dplyr::mutate(cond = n_true == n_vars)
  } else if (combiner == "any") {
    df <- df %>% dplyr::mutate(cond = n_true >= 1L)
  } else {
    if (is.null(k)) stop("Provide k for combiner='k_of_n'")
    kk <- as.integer(k)
    df <- df %>% dplyr::mutate(cond = n_true >= !!kk)
  }

  # ---- 3) Run-length logic (online with per-run correction)
  
  # We avoid self-reference by computing a TRUE-run counter via grouping trick.
  df <- df %>%
    dplyr::arrange(.data[[id]], season_id, day_idx) %>%
    dplyr::group_by(.data[[id]], season_id) %>%
    dplyr::mutate(
      # Group index that increases on FALSE -> gives 1,2,3,... within TRUE runs
      grp = cumsum(!cond),
      run_len = dplyr::if_else(cond,
        as.integer(ave(seq_along(cond), grp, FUN = seq_along)),
        0L
      ),
      valid = cond & (run_len >= min_duration),
      # Bonus: add min_duration-1 days exactly when a run first hits the threshold
      bonus = as.integer(run_len == min_duration) * pmax(min_duration - 1L, 0L),
      cum_days = cumsum(as.integer(valid) + bonus)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(S %>% dplyr::select(season_id, season_label), by = "season_id") %>%
    dplyr::arrange(.data[[as.character(id)]], season_id, day_idx) %>%
    dplyr::select(!!id, day_idx, season_id, season_label, cond, valid, cum_days)

  return(df)
}

