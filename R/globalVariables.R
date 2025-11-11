#' Internal global variable declarations for R CMD check
#'
#' Declares symbols used in tidy evaluation pipelines to silence
#' "no visible binding for global variable" notes during R CMD check.
#'
#' @noRd
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ".data", "bonus", "cond", "cum_days", "cum_exc", "cum_excess",
    "day_idx", "end", "excess", "grp", "n_true", "n_vars", "prev_cum",
    "prev_k", "run_len", "season_id", "season_label", "start", "valid",
    "value"
  ))
}


