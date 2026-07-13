#' @title Print a MEMWAS model summary
#' @description Print the summary produced by `summary.MEMWAS_fit()`.
#' @param x summary.MEMWAS_fit. Summary object to print.
#' @param ... list. Additional arguments; currently unused.
#' @returns invisible. Invisibly returns `x`.
#' @examples
#' \dontrun{
#' print(summary(fit))
#' }
#' @export
print.summary.MEMWAS_fit <- function(x, ...) {

  cat("\n================ MEMWAS Model Summary ================\n")

  cat("\nCall:\n")
  print(x$call)

  cat("\nModel:\n")
  cat("  Family :", x$family, "\n")
  cat("  Final approximation:", x$approximation_label %||% x$approximation %||% "NA", "\n")
  if (!is.na(x$init_approximation %||% NA_character_)) cat("  Initialization:", x$init_approximation_label %||% x$init_approximation, "\n")
  if (!is.na(x$se_source %||% NA_character_)) cat("  SE source:", x$se_source, "\n")
  if (!is.na(x$family_link %||% NA_character_)) cat("  Link:", x$family_link, "\n")
  cat("  Formula:", deparse(x$formula), "\n")

  cat("\n--- Fixed Effects ---\n")
  print(x$coefficients, row.names = FALSE)

  cat("\n--- Random Effects (BLUPs) ---\n")
  if (is.data.frame(x$random_effects_preview) && nrow(x$random_effects_preview) > 0L) {
    cat("Preview of ", nrow(x$random_effects_preview), " of ", x$random_effects_n_rows,
        " BLUP row(s). Intervals use level ", x$random_effects_ci_level, ".\n", sep = "")
    print(x$random_effects_preview, row.names = FALSE)
  } else {
    cat("No subject-specific random-effect BLUP rows stored.\n")
  }
  if (!is.null(x$random_covariance)) {
    cat("Random-effect covariance matrix:\n")
    print(x$random_covariance)
  }

  cat("\n--- Model Fit ---\n")
  print(x$metrics)

  cat("\n--- GLMM Assumption Screening Methods ---\n")
  if (!is.null(x$assumption_check_methods) && is.data.frame(x$assumption_check_methods) &&
      nrow(x$assumption_check_methods) > 0L) {
    print(x$assumption_check_methods, row.names = FALSE)
  } else {
    cat("No GLMM assumption-screening methods stored.\n")
  }

  cat("\n--- GLMM Assumption Screening Results ---\n")
  if (!is.null(x$assumption_check_model) && !is.na(x$assumption_check_model)) {
    cat("Screening model:", x$assumption_check_model, "\n")
  }
  .memwas_print_assumption_results(x$assumption_check_table)
  if (!is.null(x$assumption_check_note) && !is.na(x$assumption_check_note) && nzchar(x$assumption_check_note)) {
    cat("Note:", x$assumption_check_note, "\n")
  }

  if (is.data.frame(x$profile_ci) && nrow(x$profile_ci) > 0L) {
    cat("\n--- Profile likelihood CIs for variance/dispersion parameters ---\n")
    print(x$profile_ci, row.names = FALSE)
  }

  cat("\n--- Non-linear Effects ---\n")

  if (length(x$spline_variables) == 0) {
    cat("None\n")
  } else {
    for (v in x$spline_variables) {

      cat("\nVariable:", v, "\n")

      tp <- x$turning_points[[v]]
      cat("Turning points:", paste(signif(tp, 5), collapse = ", "), "\n")

      cat("Segment-specific slopes:\n")
      print(x$segment_effects[[v]], row.names = FALSE)
    }
  }

  cat("\n=====================================================\n")
}


# --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
# S3 METHOD: predict for MEMWAS_fit objects -----
# --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
