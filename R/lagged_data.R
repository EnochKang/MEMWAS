#' @title Create within-subject lagged panel variables
#' @description Sort panel data by subject and time, then create one-period lagged response and predictor columns without crossing subject boundaries.
#' @param data data.frame. Long-format panel data.
#' @param subject_id character. Name of the subject identifier column.
#' @param time_var character. Name of the measurement-time column used for ordering rows within subject.
#' @param response_var character. Name of the response variable to lag.
#' @param predictor_vars character. Names of predictor variables to lag.
#' @param remove_NA logical. Whether to remove rows whose newly created lag columns are missing.
#' @returns data.frame. The input data sorted by subject and time with added `lag_*` columns.
#' @examples
#' \dontrun{
#' dat <- simulate_panel_data(n_id = 10, n_time = 4)
#' lagged <- lagged_data(dat, "id", "time", "var_1", c("var_2", "var_3"), TRUE)
#' }
#'
#' @noRd
.lagged_data <- function(
    data,
    subject_id,
    time_var,
    response_var,
    predictor_vars,
    remove_NA = FALSE) {

  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  # 1. Input validation -----
  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }

  if (!is.character(subject_id) || length(subject_id) != 1 || is.na(subject_id)) {
    stop("`subject_id` must be a single character string.", call. = FALSE)
  }

  if (!is.character(time_var) || length(time_var) != 1 || is.na(time_var)) {
    stop("`time_var` must be a single character string.", call. = FALSE)
  }

  if (!is.character(response_var) || length(response_var) != 1 || is.na(response_var)) {
    stop("`response_var` must be a single character string.", call. = FALSE)
  }

  if (!is.character(predictor_vars) || length(predictor_vars) == 0 || anyNA(predictor_vars)) {
    stop("`predictor_vars` must be a non-empty character vector.", call. = FALSE)
  }

  if (!is.logical(remove_NA) || length(remove_NA) != 1 || is.na(remove_NA)) {
    stop("`remove_NA` must be TRUE or FALSE.", call. = FALSE)
  }

  input_vars <- c(subject_id, time_var, response_var, predictor_vars)
  missing_vars <- setdiff(input_vars, names(data))

  if (length(missing_vars) > 0) {
    stop(
      "The following variables are not found in `data`: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }

  if (anyDuplicated(predictor_vars)) {
    stop("`predictor_vars` contains duplicated variable names.", call. = FALSE)
  }

  if (anyNA(data[[subject_id]])) {
    stop("`subject_id` contains missing values.", call. = FALSE)
  }

  if (anyNA(data[[time_var]])) {
    stop("`time_var` contains missing values.", call. = FALSE)
  }

  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  # 2. Sort by subject ID and time -----
  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  data <- data[order(data[[subject_id]], data[[time_var]]), , drop = FALSE]
  rownames(data) <- NULL

  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  # 3. Create lag names -----
  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  lag_dep <- paste0("lag_", response_var)
  lag_inds <- paste0("lag_", predictor_vars)
  lag_vars <- c(lag_dep, lag_inds)

  existing_lags <- intersect(lag_vars, names(data))

  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  # 4. Initialize lag columns with NA -----
  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  data[[lag_dep]] <- data[[response_var]][rep(NA_integer_, nrow(data))]

  for (v in seq_along(predictor_vars)) {
    data[[lag_inds[v]]] <- data[[predictor_vars[v]]][rep(NA_integer_, nrow(data))]
  }

  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  # 5. Process lags per subject ID to avoid cross-subject contamination -----
  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  ids <- unique(data[[subject_id]])

  for (id in ids) {
    idx <- which(data[[subject_id]] == id)

    if (length(idx) > 1) {
      data[idx[-1], lag_dep] <- data[idx[-length(idx)], response_var]

      for (v in seq_along(predictor_vars)) {
        data[idx[-1], lag_inds[v]] <- data[idx[-length(idx)], predictor_vars[v]]
      }
    }
  }

  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  # 6. Remove rows where lags are NA if requested -----
  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  if (remove_NA) {
    data <- data[complete.cases(data[, lag_vars, drop = FALSE]), , drop = FALSE]
    rownames(data) <- NULL
  }

  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  # 7. Return lagged data -----
  # --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
  return(data)
}


# --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
# FUNCTION: Main function of MEMWAS -----
# --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
