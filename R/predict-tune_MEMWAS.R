#' @title Predict from the best tuned MEMWAS model
#' @description Generate link-scale or response-scale predictions from the final model refit stored inside a `tune_MEMWAS` object.
#' @param object tune_MEMWAS. Tuning object returned by `tune_MEMWAS()`.
#' @param newdata data.frame. Optional new data; when omitted, predictions are computed for the final model's original data.
#' @param type character. Prediction scale: `"link"` or `"response"`.
#' @param include_random logical. Whether to include random effects. This argument is passed to `predict.MEMWAS_fit()`; fixed-effect predictions are used for new data.
#' @param ... list. Additional arguments passed to `predict.MEMWAS_fit()`.
#' @returns numeric. Predicted values from the best final MEMWAS model on the requested scale.
#' @examples
#' \dontrun{
#' pred <- predict(tuned, type = "response")
#' }
#' @export
predict.tune_MEMWAS <- function(object, newdata = NULL,
                                type = c("link", "response"),
                                include_random = FALSE, ...) {
  if (!inherits(object, "tune_MEMWAS")) {
    stop("`object` must be a tune_MEMWAS object.", call. = FALSE)
  }
  type <- match.arg(type)
  if (is.null(object$best_fit)) {
    msg <- "No final MEMWAS fit is stored in this tune_MEMWAS object."
    if (!is.null(object$final_fit_error) && length(object$final_fit_error) > 0L &&
        !is.na(object$final_fit_error[1L]) && nzchar(object$final_fit_error[1L])) {
      msg <- paste(msg, object$final_fit_error[1L])
    } else {
      msg <- paste(msg, "Run tune_MEMWAS() with `refit_final = TRUE` to enable prediction.")
    }
    stop(msg, call. = FALSE)
  }
  stats::predict(object$best_fit, newdata = newdata, type = type,
                 include_random = include_random, ...)
}
