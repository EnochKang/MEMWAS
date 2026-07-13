#' @title Predict from a fitted MEMWAS model
#' @description Generate link-scale or response-scale fixed-effect predictions from a `MEMWAS_fit` object, rebuilding stored spline terms for new data when needed.
#' @param object MEMWAS_fit. A fitted model object returned by `fit_MEMWAS()`.
#' @param newdata data.frame. Optional new data; when omitted, predictions are computed for the original model matrix.
#' @param type character. Prediction scale: `"link"` or `"response"`.
#' @param include_random logical. Whether to include random effects; currently only fixed-effect prediction is implemented for new data.
#' @param ... list. Additional arguments; currently unused.
#' @returns numeric. Predicted values on the requested scale.
#' @examples
#' \dontrun{
#' pred <- predict(fit, type = "response")
#' }
#' @export
predict.MEMWAS_fit <- function(object, newdata = NULL,
                               type = c("link", "response"),
                               include_random = FALSE, ...) {

  rcs_basis <- function(x, knots) {
    x <- as.numeric(x)
    knots <- sort(unique(as.numeric(knots)))
    K <- length(knots)
    if (K < 3L) stop("At least three unique knots are required for restricted cubic splines.", call. = FALSE)
    tp <- function(z, k) pmax(z - k, 0)^3
    denom <- max(knots[K] - knots[K - 1L], 1e-12)
    B <- matrix(NA_real_, nrow = length(x), ncol = K - 2L)
    for (j in seq_len(K - 2L)) {
      B[, j] <- tp(x, knots[j]) - tp(x, knots[K - 1L]) * ((knots[K] - knots[j]) / denom) +
        tp(x, knots[K]) * ((knots[K - 1L] - knots[j]) / denom)
    }
    scale <- max(diff(range(knots)), 1e-12)^3
    B <- B / scale
    colnames(B) <- paste0("spline", seq_len(ncol(B)))
    B
  }

  type <- match.arg(type)

  if (is.null(newdata)) {
    X <- object$X
  } else {

    data_new <- newdata

    # rebuild spline terms
    for (v in object$spline_variables) {

      info <- object$settings$spline_info[[v]]
      knots <- info$knots
      basis_names <- info$basis_names

      B <- rcs_basis(data_new[[v]], knots)

      for (j in seq_len(ncol(B))) {
        data_new[[basis_names[j]]] <- B[, j]
      }
    }

    X <- stats::model.matrix(stats::delete.response(stats::terms(object$formal_formula)), data_new)
  }

  beta <- object$coefficients

  eta <- as.vector(X %*% beta)

  if (include_random && !is.null(object$random_effects)) {
    warning("Random effects prediction only implemented for original data.")
  }

  if (type == "link") return(eta)

  # response transformation
  if (object$family == "gaussian") return(eta)

  if (object$family == "binomial") {
    return(1 / (1 + exp(-eta)))
  }

  if (object$family %in% c("poisson", "negative_binomial", "gamma", "exponential")) {
    return(exp(pmin(pmax(eta, -30), 30)))
  }

  stop("Unsupported family stored in MEMWAS fit: ", object$family, call. = FALSE)
}

# --- # --- # --- # --- # --- # --- # --- # --- # --- # ---
# FUNCTION: Hyperparameter tuning for MEMWAS by grouped K-fold CV -----
# --- # --- # --- # --- # --- # --- # --- # --- # --- # ---

# tune_MEMWAS() tunes the MEMWAS fixed-effect elastic-net penalties.
# It uses lambda as the overall regularization strength and alpha as the
# L1/L2 mixing parameter:
#   L1_penalty = alpha * lambda
#   L2_penalty = (1 - alpha) * lambda
#
# The function expects fit_MEMWAS() to be available in the environment. It can be
# called with a MEMWAS settings object produced by set_MEMWAS(), or with a
# formula plus data/id/time arguments, matching the revised fit_MEMWAS() interface.
