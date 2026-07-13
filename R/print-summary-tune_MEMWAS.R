#' @title Print a MEMWAS tuning summary
#' @description Print the summary produced by `summary.tune_MEMWAS()`, including best final-model formulas, coefficients, penalties, and autocorrelation information when available.
#' @param x summary.tune_MEMWAS. Summary object to print.
#' @param ... list. Additional arguments; currently unused.
#' @returns invisible. Invisibly returns `x`.
#' @examples
#' \dontrun{
#' print(summary(tuned))
#' }
#' @export
print.summary.tune_MEMWAS <- function(x, ...) {
  cat("\n================ MEMWAS Tuning Summary ================\n")

  cat("\nCall:\n")
  print(x$call)

  cat("\nApproximation:\n")
  cat("  Final: ", x$approximation_label %||% x$approximation %||% "NA", "\n", sep = "")
  if (!is.na(x$init_approximation %||% NA_character_)) cat("  Initialization: ", .memwas_approximation_label(x$init_approximation), "\n", sep = "")
  if (!is.na(x$se_method %||% NA_character_)) cat("  SE method: ", paste(x$se_method, collapse = ", "), "\n", sep = "")

  cat("\nBest hyperparameters and CV performance:\n")
  b <- x$best[1L, , drop = FALSE]
  keep <- intersect(c("lambda", "log_lambda", "alpha", "L1_penalty", "L2_penalty",
                      "median_metric", "mean_metric", "stability", "selection_score",
                      "n_successful_folds", "n_folds", "failed_folds",
                      "final_fit_available", "n_fixed_effects", "n_selected_fixed_effects"),
                    names(b))
  print(b[, keep, drop = FALSE], row.names = FALSE)

  cat("\nBest final model formulas:\n")
  if (is.null(x$formulas)) {
    cat("  Structured best-model formulas are not available.\n")
  } else {
    .memwas_cat_wrapped("original formula: ", x$formulas$original)
    .memwas_cat_wrapped("formal formula: ", x$formulas$formal)
    .memwas_cat_wrapped("fixed effects: ", x$formulas$fixed_effects_with_coefficients)
    .memwas_cat_wrapped("random effects: ", x$formulas$random_effects_with_coefficients)
    .memwas_cat_wrapped("fixed penalty: ", x$formulas$fixed_effect_penalty)
    .memwas_cat_wrapped("autocorrelation: ", x$formulas$autocorrelation)
    .memwas_cat_wrapped("autocor penalty: ", x$formulas$autocorrelation_penalty)
  }

  cat("\nFixed-effect coefficients:\n")
  if (is.data.frame(x$fixed_effects) && nrow(x$fixed_effects) > 0L) {
    print(x$fixed_effects, row.names = FALSE)
  } else {
    cat("  Not available. Refit the final model with `refit_final = TRUE` to store coefficients.\n")
  }

  cat("\nRandom effects and covariance:\n")
  if (!is.null(x$random_covariance)) {
    cat("  Random-effect covariance matrix:\n")
    print(x$random_covariance)
  } else {
    cat("  Random-effect covariance matrix not available.\n")
  }
  if (is.data.frame(x$random_effects_preview) && nrow(x$random_effects_preview) > 0L) {
    cat("  Random-effect BLUP preview (", nrow(x$random_effects_preview), " of ",
        x$random_effects_n_rows, " rows):\n", sep = "")
    print(x$random_effects_preview, row.names = FALSE)
  } else {
    cat("  No subject-specific random-effect BLUP rows stored.\n")
  }

  cat("\nPenalties:\n")
  if (!is.null(x$penalties)) {
    cat("  fixed penalty value:", .memwas_format_number(x$penalties$fixed_effect_value), "\n")
    cat("  autocorrelation penalty value:", .memwas_format_number(x$penalties$autocorrelation_value), "\n")
  } else {
    cat("  Penalty details are not available.\n")
  }

  cat("\nAutocorrelation estimates:\n")
  if (!is.null(x$autocorrelation) && is.data.frame(x$autocorrelation$estimates) && nrow(x$autocorrelation$estimates) > 0L) {
    print(x$autocorrelation$estimates, row.names = FALSE)
  } else {
    cat("  No autocorrelation parameter estimates stored, or structure is independence.\n")
  }

  cat("\nTop observed tuning results:\n")
  if (is.data.frame(x$top_results) && nrow(x$top_results) > 0L) {
    print(x$top_results, row.names = FALSE)
  } else {
    cat("  No tuning results available.\n")
  }

  cat("\nConvergence:\n")
  print(x$convergence)

  if (!is.null(x$final_fit_error) && length(x$final_fit_error) > 0L &&
      !is.na(x$final_fit_error[1L]) && nzchar(x$final_fit_error[1L])) {
    cat("\nFinal fit note/error:\n", x$final_fit_error[1L], "\n", sep = "")
  }

  cat("\n======================================================\n")
  invisible(x)
}
