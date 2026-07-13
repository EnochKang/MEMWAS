# Internal helpers used by tune_MEMWAS(). These functions are intentionally not exported.

.memwas_augment_control <- function(control, autocor_regularization = NULL,
                                    autocor_regularization_type = "l2",
                                    autocor_regularization_alpha = 0) {
  control <- control %||% list()
  if (!is.list(control)) stop("`control` must be a list.", call. = FALSE)
  if (!is.null(autocor_regularization)) {
    if (is.list(autocor_regularization)) {
      reg <- autocor_regularization
      if (is.null(reg$enabled)) reg$enabled <- TRUE
      if (is.null(reg$type)) reg$type <- autocor_regularization_type
      if (is.null(reg$alpha)) reg$alpha <- autocor_regularization_alpha
    } else {
      if (!is.numeric(autocor_regularization) || length(autocor_regularization) != 1L ||
          is.na(autocor_regularization) || autocor_regularization < 0) {
        stop("`autocor_regularization` must be NULL, a non-negative numeric scalar, or a list.", call. = FALSE)
      }
      reg <- list(enabled = TRUE,
                  lambda = as.numeric(autocor_regularization),
                  type = autocor_regularization_type,
                  alpha = autocor_regularization_alpha)
    }
    control$autocor_regularization <- reg
    control$autocor_penalty <- reg$lambda %||% control$autocor_penalty
    control$autocor_regularization_type <- reg$type %||% autocor_regularization_type
    control$autocor_regularization_alpha <- reg$alpha %||% autocor_regularization_alpha
  }
  control
}


.memwas_make_id_fold_split <- function(id, K, engine = "R") {
  backend <- .memwas_normalize_engine(engine)

  row_folds <- .memwas_make_fold_assignment(id, K = K, engine = backend)
  id_chr <- as.character(id)
  unique_ids <- unique(id_chr)
  id_fold <- vapply(unique_ids, function(sid) {
    u <- unique(row_folds[id_chr == sid])
    if (length(u) != 1L) {
      stop("Internal error: subject ID `", sid, "` was assigned to multiple folds.", call. = FALSE)
    }
    as.integer(u[1L])
  }, integer(1L))
  id_n <- vapply(unique_ids, function(sid) sum(id_chr == sid), integer(1L))
  id_table <- data.frame(id = unique_ids, fold = as.integer(id_fold), n_rows = as.integer(id_n),
                         stringsAsFactors = FALSE, check.names = FALSE)
  list(row_folds = as.integer(row_folds), id_folds = id_table)
}

.memwas_response_from_data <- function(formula, data, family, mape_epsilon = .Machine$double.eps) {
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass,
                           drop.unused.levels = TRUE)
  y <- stats::model.response(mf)
  if (is.factor(y)) y <- as.numeric(y) - 1
  if (is.logical(y)) y <- as.numeric(y)
  if (is.matrix(y) && family == "binomial" && ncol(y) == 2L) {
    denom <- rowSums(y)
    y <- y[, 1L] / pmax(denom, mape_epsilon)
  }
  as.numeric(y)
}

.memwas_complete_validation_data <- function(data_i, settings) {
  random_formula <- .memwas_parse_random_formula(settings$random)
  serial_value <- settings$serial %||% NULL
  if (is.null(serial_value) && exists(".memwas_is_multi_serial_autocor", mode = "function") &&
      .memwas_is_multi_serial_autocor(settings$autocor %||% NULL)) {
    serial_value <- settings$autocor
  }
  serial_vars <- character(0L)
  if (!is.null(serial_value) && exists(".memwas_normalize_serial_components", mode = "function")) {
    defs <- try(.memwas_normalize_serial_components(serial_value, data_i, parent.frame()), silent = TRUE)
    if (!inherits(defs, "try-error")) {
      serial_vars <- unique(unlist(lapply(defs, function(z) all.vars(z$design_formula)), use.names = FALSE))
    }
  }
  needed <- unique(c(all.vars(settings$formula), all.vars(random_formula), serial_vars,
                     settings$id, settings$time, settings$spline_variables %||% character(0L)))
  needed <- needed[needed %in% names(data_i)]
  cc <- stats::complete.cases(data_i[, needed, drop = FALSE])
  data_i[cc, , drop = FALSE]
}

.memwas_lambda_from_log <- function(log_lambda, lambda_base = 10) {
  lambda_base^log_lambda
}

.memwas_log_from_lambda <- function(lambda, lambda_base = 10) {
  log(lambda) / log(lambda_base)
}

.memwas_selection_score <- function(med, disp, use_stability_penalty = TRUE,
                                    stability_weight = 1) {
  med + if (isTRUE(use_stability_penalty)) stability_weight * disp else 0
}

.memwas_lhs_design <- function(n, d) {
  n <- as.integer(n)
  d <- as.integer(d)
  if (n <= 0L) return(matrix(numeric(0L), nrow = 0L, ncol = d))
  X <- matrix(NA_real_, nrow = n, ncol = d)
  for (j in seq_len(d)) {
    X[, j] <- (sample.int(n, n, replace = FALSE) - stats::runif(n)) / n
  }
  X
}

.memwas_scale_inputs <- function(X, lambda_log_range, alpha_range) {
  X <- as.matrix(X)
  X[, 1L] <- (X[, 1L] - lambda_log_range[1L]) / diff(lambda_log_range)
  X[, 2L] <- (X[, 2L] - alpha_range[1L]) / diff(alpha_range)
  X
}

.memwas_dedupe_points <- function(points, lambda_log_range, alpha_range, digits = 12L) {
  points <- as.data.frame(points)
  points$log_lambda <- as.numeric(points$log_lambda)
  points$alpha <- as.numeric(points$alpha)
  points <- points[is.finite(points$log_lambda) & is.finite(points$alpha), , drop = FALSE]
  if (nrow(points) == 0L) return(points[, c("log_lambda", "alpha"), drop = FALSE])
  points$log_lambda <- .memwas_clamp(points$log_lambda, lambda_log_range[1L], lambda_log_range[2L])
  points$alpha <- .memwas_clamp(points$alpha, alpha_range[1L], alpha_range[2L])
  key <- paste(round(points$log_lambda, digits), round(points$alpha, digits), sep = "_")
  points <- points[!duplicated(key), c("log_lambda", "alpha"), drop = FALSE]
  row.names(points) <- NULL
  points
}

.memwas_parse_initial_points <- function(points, lambda_log_range, alpha_range,
                                         lambda_base = 10) {
  if (is.null(points)) return(data.frame(log_lambda = numeric(0L), alpha = numeric(0L)))
  points <- as.data.frame(points)
  if (!"log_lambda" %in% names(points)) {
    if ("lambda" %in% names(points)) {
      points$log_lambda <- .memwas_log_from_lambda(as.numeric(points$lambda), lambda_base = lambda_base)
    } else {
      stop("`initial_points` must contain either `log_lambda` or `lambda`.", call. = FALSE)
    }
  }
  if (!"alpha" %in% names(points)) stop("`initial_points` must contain `alpha`.", call. = FALSE)
  out <- data.frame(log_lambda = as.numeric(points$log_lambda),
                    alpha = as.numeric(points$alpha))
  if (any(!is.finite(out$log_lambda)) || any(!is.finite(out$alpha))) {
    stop("`initial_points` contains non-finite values.", call. = FALSE)
  }
  if (any(out$log_lambda < lambda_log_range[1L] | out$log_lambda > lambda_log_range[2L]) ||
      any(out$alpha < alpha_range[1L] | out$alpha > alpha_range[2L])) {
    stop("`initial_points` must lie inside `lambda_log_range` and `alpha_range`.", call. = FALSE)
  }
  out
}

.memwas_make_initial_points <- function(lambda_log_range, alpha_range,
                                        initial_design = "lhs", n_initial = 10L,
                                        initial_points = NULL, lambda_base = 10) {
  lo <- lambda_log_range[1L]
  hi <- lambda_log_range[2L]
  alo <- alpha_range[1L]
  ahi <- alpha_range[2L]
  mid <- mean(lambda_log_range)
  amid <- mean(alpha_range)
  boundary <- data.frame(log_lambda = c(lo, lo, hi, hi, mid),
                         alpha = c(alo, ahi, alo, ahi, amid))
  if (initial_design == "lhs") {
    n_lhs <- max(0L, n_initial - nrow(boundary))
    U <- .memwas_lhs_design(n_lhs, 2L)
    lhs <- if (nrow(U) > 0L) {
      data.frame(log_lambda = lo + U[, 1L] * (hi - lo),
                 alpha = alo + U[, 2L] * (ahi - alo))
    } else {
      data.frame(log_lambda = numeric(0L), alpha = numeric(0L))
    }
    pts <- rbind(boundary, lhs,
                 .memwas_parse_initial_points(initial_points, lambda_log_range, alpha_range, lambda_base = lambda_base))
  } else {
    pts <- rbind(boundary,
                 .memwas_parse_initial_points(initial_points, lambda_log_range, alpha_range, lambda_base = lambda_base))
  }
  .memwas_dedupe_points(pts, lambda_log_range, alpha_range)
}

.memwas_make_candidate_points <- function(n, lambda_log_range, alpha_range) {
  lo <- lambda_log_range[1L]
  hi <- lambda_log_range[2L]
  alo <- alpha_range[1L]
  ahi <- alpha_range[2L]
  U <- .memwas_lhs_design(n, 2L)
  pts <- data.frame(log_lambda = lo + U[, 1L] * (hi - lo),
                    alpha = alo + U[, 2L] * (ahi - alo))
  fixed <- data.frame(log_lambda = c(lo, lo, hi, hi, mean(lambda_log_range)),
                      alpha = c(alo, ahi, alo, ahi, mean(alpha_range)))
  .memwas_dedupe_points(rbind(fixed, pts), lambda_log_range, alpha_range)
}

.memwas_not_observed <- function(candidates, observed) {
  if (nrow(observed) == 0L || nrow(candidates) == 0L) return(candidates)
  key_c <- paste(round(candidates$log_lambda, 10L), round(candidates$alpha, 10L), sep = "_")
  key_o <- paste(round(observed$log_lambda, 10L), round(observed$alpha, 10L), sep = "_")
  candidates[!key_c %in% key_o, , drop = FALSE]
}

.memwas_distance_uncertainty <- function(Xnew, Xobs, yobs, lambda_log_range, alpha_range) {
  Xn <- .memwas_scale_inputs(Xnew, lambda_log_range, alpha_range)
  Xo <- .memwas_scale_inputs(Xobs, lambda_log_range, alpha_range)
  yscale <- stats::sd(yobs[is.finite(yobs)], na.rm = TRUE)
  if (!is.finite(yscale) || yscale <= 0) yscale <- max(abs(yobs[is.finite(yobs)]), na.rm = TRUE)
  if (!is.finite(yscale) || yscale <= 0) yscale <- 1
  apply(Xn, 1L, function(x) {
    d <- sqrt(min(rowSums((t(t(Xo) - x))^2)))
    yscale * max(d, 0.05)
  })
}

.memwas_fit_qrs <- function(X, y) {
  X <- as.matrix(X)
  df <- data.frame(y = as.numeric(y), log_lambda = X[, 1L], alpha = X[, 2L])
  fit <- try(stats::lm(y ~ log_lambda + alpha + I(log_lambda^2) + I(alpha^2) + I(log_lambda * alpha),
                       data = df), silent = TRUE)
  if (inherits(fit, "try-error")) fit <- NULL
  list(type = "QRS", fit = fit, X = X, y = y)
}

.memwas_predict_qrs <- function(model, Xnew, lambda_log_range, alpha_range) {
  Xnew <- as.matrix(Xnew)
  mu0 <- mean(model$y, na.rm = TRUE)
  if (!is.finite(mu0)) mu0 <- 0
  se0 <- .memwas_distance_uncertainty(Xnew, model$X, model$y, lambda_log_range, alpha_range)
  if (is.null(model$fit)) return(list(mean = rep(mu0, nrow(Xnew)), se = se0))
  newdata <- data.frame(log_lambda = Xnew[, 1L], alpha = Xnew[, 2L])
  pr <- try(stats::predict(model$fit, newdata = newdata, se.fit = TRUE), silent = TRUE)
  if (inherits(pr, "try-error")) return(list(mean = rep(mu0, nrow(Xnew)), se = se0))
  mu <- as.numeric(pr$fit)
  se <- as.numeric(pr$se.fit)
  mu[!is.finite(mu)] <- mu0
  se[!is.finite(se) | se <= 0] <- se0[!is.finite(se) | se <= 0]
  se <- pmax(se, 0.25 * se0)
  list(mean = mu, se = se)
}

.memwas_fit_gp <- function(X, y, lambda_log_range, alpha_range, engine = "R") {
  backend <- .memwas_normalize_engine(engine)

  X <- as.matrix(X)
  Xn <- .memwas_scale_inputs(X, lambda_log_range, alpha_range)
  y <- as.numeric(y)
  y_mean <- mean(y, na.rm = TRUE)
  y_sd <- stats::sd(y, na.rm = TRUE)
  if (!is.finite(y_sd) || y_sd <= 0) y_sd <- 1
  ys <- (y - y_mean) / y_sd
  n <- length(ys)

  nll <- function(par) {
    ell <- exp(par[seq_len(ncol(Xn))])
    noise <- exp(par[ncol(Xn) + 1L])
    Kmat <- .memwas_kernel_matrix(Xn, Xn, ell, engine = backend) + diag(noise^2 + 1e-8, n)
    R <- try(chol(Kmat), silent = TRUE)
    if (inherits(R, "try-error")) return(1e100)
    alpha <- backsolve(R, forwardsolve(t(R), ys))
    val <- 0.5 * sum(ys * alpha) + sum(log(diag(R))) + 0.5 * n * log(2 * pi)
    if (!is.finite(val)) 1e100 else val
  }

  par0 <- log(c(rep(0.35, ncol(Xn)), 0.05))
  opt <- try(stats::optim(par0, nll, method = "L-BFGS-B",
                          lower = log(c(rep(0.02, ncol(Xn)), 1e-6)),
                          upper = log(c(rep(5, ncol(Xn)), 1))), silent = TRUE)
  if (inherits(opt, "try-error") || !is.finite(opt$value)) opt <- list(par = par0, convergence = 1L)
  ell <- exp(opt$par[seq_len(ncol(Xn))])
  noise <- exp(opt$par[ncol(Xn) + 1L])
  Kmat <- .memwas_kernel_matrix(Xn, Xn, ell, engine = backend) + diag(noise^2 + 1e-8, n)
  R <- try(chol(Kmat), silent = TRUE)
  if (inherits(R, "try-error")) return(.memwas_fit_qrs(X, y))
  alpha <- backsolve(R, forwardsolve(t(R), ys))
  list(type = "GP", X = X, Xn = Xn, y = y, y_mean = y_mean, y_sd = y_sd,
       ell = ell, noise = noise, R = R, alpha = alpha, convergence = opt$convergence)
}

.memwas_predict_gp <- function(model, Xnew, lambda_log_range, alpha_range, engine = "R") {
  backend <- .memwas_normalize_engine(engine)

  if (!identical(model$type, "GP")) return(.memwas_predict_qrs(model, Xnew, lambda_log_range, alpha_range))
  Xnew <- as.matrix(Xnew)
  Xnewn <- .memwas_scale_inputs(Xnew, lambda_log_range, alpha_range)
  Ks <- .memwas_kernel_matrix(Xnewn, model$Xn, model$ell, engine = backend)
  mu_s <- as.vector(Ks %*% model$alpha)
  v <- forwardsolve(t(model$R), t(Ks))
  var_s <- pmax(1 - colSums(v^2), 1e-10)
  mu <- model$y_mean + model$y_sd * mu_s
  se <- model$y_sd * sqrt(var_s)
  bad <- !is.finite(se) | se <= 0
  if (any(bad)) {
    se[bad] <- .memwas_distance_uncertainty(Xnew, model$X, model$y, lambda_log_range, alpha_range)[bad]
  }
  list(mean = as.numeric(mu), se = as.numeric(se))
}

.memwas_fit_surrogate <- function(X, y, surrogate = "QRS", lambda_log_range,
                                  alpha_range, engine = "R") {
  backend <- .memwas_normalize_engine(engine)

  X <- as.matrix(X)
  y <- as.numeric(y)
  ok <- is.finite(y) & is.finite(X[, 1L]) & is.finite(X[, 2L])
  X <- X[ok, , drop = FALSE]
  y <- y[ok]
  if (length(y) == 0L) stop("No finite cross-validation results are available for surrogate fitting.", call. = FALSE)
  if (surrogate == "GP" && length(y) >= 2L) {
    out <- try(.memwas_fit_gp(X, y, lambda_log_range, alpha_range, engine = backend), silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
  }
  .memwas_fit_qrs(X, y)
}

.memwas_predict_surrogate <- function(model, Xnew, lambda_log_range, alpha_range,
                                      engine = "R") {
  backend <- .memwas_normalize_engine(engine)

  if (identical(model$type, "GP")) {
    .memwas_predict_gp(model, Xnew, lambda_log_range, alpha_range, engine = backend)
  } else {
    .memwas_predict_qrs(model, Xnew, lambda_log_range, alpha_range)
  }
}

.memwas_make_settings_for_pair <- function(settings, data_i, L1_penalty, L2_penalty,
                                           engine = "R", autocor_regularization = NULL,
                                           autocor_regularization_type = "l2",
                                           autocor_regularization_alpha = 0) {
  s <- settings
  s$data <- data_i
  s$L1_penalty <- L1_penalty
  s$L2_penalty <- L2_penalty
  s$control <- .memwas_augment_control(s$control,
                                       autocor_regularization = autocor_regularization,
                                       autocor_regularization_type = autocor_regularization_type,
                                       autocor_regularization_alpha = autocor_regularization_alpha)
  s$engine <- engine
  s$verbose <- FALSE
  if (is.null(s$spline_variables)) s$spline_variables <- character(0L)
  if (is.null(s$spline_info)) s$spline_info <- list()
  class(s) <- "MEMWAS"
  s
}

.memwas_evaluate_penalty_pair <- function(log_lambda, alpha, iteration, stage,
                                          settings, row_folds, K,
                                          metric, stability_metric,
                                          lambda_base, alpha_range,
                                          use_stability_penalty, stability_weight,
                                          mape_epsilon, fail_value,
                                          keep_fold_fits, engine,
                                          autocor_regularization = NULL,
                                          autocor_regularization_type = "l2",
                                          autocor_regularization_alpha = 0,
                                          ...) {
  backend <- .memwas_normalize_engine(engine)

  lambda <- .memwas_lambda_from_log(log_lambda, lambda_base = lambda_base)
  alpha <- .memwas_clamp(alpha, alpha_range[1L], alpha_range[2L])
  L1_penalty <- alpha * lambda
  L2_penalty <- (1 - alpha) * lambda
  fold_rows <- vector("list", K)
  fold_fits <- if (isTRUE(keep_fold_fits)) vector("list", K) else NULL

  for (fold in seq_len(K)) {
    train_data <- settings$data[row_folds != fold, , drop = FALSE]
    valid_data <- settings$data[row_folds == fold, , drop = FALSE]
    valid_data <- .memwas_complete_validation_data(valid_data, settings)

    fold_metric <- fail_value
    fold_n <- nrow(valid_data)
    fold_error <- ""
    fit_convergence <- NA_integer_
    pred <- y <- numeric(0L)

    if (fold_n == 0L) {
      fold_error <- "No complete validation rows in this fold."
    } else {
      train_settings <- .memwas_make_settings_for_pair(settings, train_data, L1_penalty, L2_penalty,
                                                       engine = backend,
                                                       autocor_regularization = autocor_regularization,
                                                       autocor_regularization_type = autocor_regularization_type,
                                                       autocor_regularization_alpha = autocor_regularization_alpha)
      fit <- try(fit_MEMWAS(train_settings, verbose = FALSE, ...), silent = TRUE)
      if (inherits(fit, "try-error")) {
        fold_error <- paste0("MEMWAS fit failed: ", as.character(fit))
      } else {
        if (isTRUE(keep_fold_fits)) fold_fits[[fold]] <- fit
        fit_convergence <- as.integer(fit$convergence %||% NA_integer_)
        y <- try(.memwas_response_from_data(settings$formula, valid_data, settings$family,
                                            mape_epsilon = mape_epsilon), silent = TRUE)
        pred <- try(stats::predict(fit, newdata = valid_data, type = "response"), silent = TRUE)
        if (inherits(y, "try-error")) {
          fold_error <- paste0("Response extraction failed: ", as.character(y))
        } else if (inherits(pred, "try-error")) {
          fold_error <- paste0("Prediction failed: ", as.character(pred))
        } else {
          y <- as.numeric(y)
          pred <- as.numeric(pred)
          if (length(pred) != length(y)) {
            fold_error <- paste0("Prediction length mismatch: got ", length(pred), " predictions for ", length(y), " validation rows.")
          } else {
            fold_metric <- .memwas_metric_value(y, pred, metric = metric, epsilon = mape_epsilon,
                                                fail_value = fail_value, engine = backend)
          }
        }
      }
    }

    fold_rows[[fold]] <- data.frame(
      iteration = iteration,
      stage = stage,
      fold = fold,
      log_lambda = log_lambda,
      lambda = lambda,
      alpha = alpha,
      L1_penalty = L1_penalty,
      L2_penalty = L2_penalty,
      metric = fold_metric,
      n_validation = fold_n,
      convergence = fit_convergence,
      error = fold_error,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  folds <- do.call(rbind, fold_rows)
  finite_metrics <- folds$metric[is.finite(folds$metric)]
  median_metric <- if (length(finite_metrics)) stats::median(finite_metrics) else fail_value
  mean_metric <- if (length(finite_metrics)) mean(finite_metrics) else fail_value
  disp <- if (length(finite_metrics)) {
    .memwas_stability_value(finite_metrics, metric = stability_metric, engine = backend)
  } else {
    fail_value
  }
  score <- .memwas_selection_score(median_metric, disp,
                                   use_stability_penalty = use_stability_penalty,
                                   stability_weight = stability_weight)

  row <- data.frame(
    iteration = iteration,
    stage = stage,
    log_lambda = log_lambda,
    lambda = lambda,
    alpha = alpha,
    L1_penalty = L1_penalty,
    L2_penalty = L2_penalty,
    median_metric = median_metric,
    mean_metric = mean_metric,
    stability = disp,
    selection_score = score,
    n_successful_folds = length(finite_metrics),
    n_folds = K,
    min_fold_metric = if (length(finite_metrics)) min(finite_metrics) else fail_value,
    max_fold_metric = if (length(finite_metrics)) max(finite_metrics) else fail_value,
    failed_folds = sum(!is.finite(folds$metric)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  list(row = row, folds = folds, fold_fits = fold_fits)
}
