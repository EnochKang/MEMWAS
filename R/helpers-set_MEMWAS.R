# Internal helpers used by set_MEMWAS(). These functions are intentionally not exported.

.memwas_choose_spline_probs <- function(K, spline_probs = NULL) {
  if (!is.null(spline_probs)) {
    if (!is.numeric(spline_probs) || length(spline_probs) != K || any(spline_probs <= 0 | spline_probs >= 1)) {
      stop("`spline_probs` must contain exactly `spline_knots` probabilities strictly between 0 and 1.",
           call. = FALSE)
    }
    return(sort(spline_probs))
  }
  if (K == 3L) return(c(0.10, 0.50, 0.90))
  if (K == 4L) return(c(0.05, 0.35, 0.65, 0.95))
  if (K == 5L) return(c(0.05, 0.275, 0.50, 0.725, 0.95))
  seq(0.05, 0.95, length.out = K)
}

.memwas_get_spline_knots <- function(x, K, spline_probs = NULL) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(unique(x)) < K) return(sort(unique(x)))
  probs <- .memwas_choose_spline_probs(K, spline_probs = spline_probs)
  sort(unique(as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 7))))
}

.memwas_candidate_predictors <- function(formula, data, id, time, nonlinear_predictors = NULL) {
  labels <- attr(stats::terms(formula), "term.labels")
  simple <- labels[labels %in% names(data)]
  simple <- setdiff(simple, c(id, time))
  if (!is.null(nonlinear_predictors)) simple <- intersect(simple, nonlinear_predictors)
  is_num <- vapply(simple, function(v) is.numeric(data[[v]]) || is.integer(data[[v]]), logical(1L))
  simple[is_num]
}

.memwas_make_basis_names <- function(var, n_basis, existing_names) {
  root <- paste0(".MEMWAS_spl_", make.names(var), "_")
  nm <- paste0(root, seq_len(n_basis))
  if (!any(nm %in% existing_names)) return(nm)
  suffix <- 1L
  repeat {
    nm <- paste0(root, seq_len(n_basis), "_", suffix)
    if (!any(nm %in% existing_names)) return(nm)
    suffix <- suffix + 1L
  }
}

.memwas_find_turning_points <- function(grid, effect) {
  ok <- is.finite(grid) & is.finite(effect)
  grid <- grid[ok]
  effect <- effect[ok]
  if (length(grid) < 5L) return(numeric(0L))
  d <- diff(effect) / pmax(diff(grid), 1e-12)
  s <- sign(d)
  s[s == 0] <- NA
  if (all(is.na(s))) return(numeric(0L))
  for (i in seq_along(s)) {
    if (is.na(s[i])) s[i] <- if (i == 1L) 0 else s[i - 1L]
  }
  s[s == 0] <- 1
  idx <- which(diff(s) != 0L)
  if (length(idx) == 0L) return(numeric(0L))
  unique(round((grid[idx] + grid[idx + 1L]) / 2, 8))
}

.memwas_spline_turning_points <- function(var, fit, knots, basis_names, data, grid_n,
                                          engine = "R") {
  backend <- .memwas_normalize_engine(engine)

  co <- fit$coefficients
  x <- as.numeric(data[[var]])
  grid <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = grid_n)
  B <- .memwas_rcs_basis(grid, knots, engine = backend)
  beta_linear <- if (var %in% names(co)) co[var] else 0
  beta_basis <- rep(0, length(basis_names))
  names(beta_basis) <- basis_names
  common <- intersect(basis_names, names(co))
  beta_basis[common] <- co[common]
  effect <- beta_linear * grid + as.vector(B %*% beta_basis)
  .memwas_find_turning_points(grid, effect)
}

.memwas_fit_screening_model <- function(formula_i, data_i, family, id, time, random,
                                        autocor, serial = NULL,
                                        L1_penalty, L2_penalty, control,
                                        screen_method = "ML", random_cov = "diagonal",
                                        approximation = "laplace",
                                        init_approximation = "variational_inference",
                                        se_method = "hessian",
                                        engine = "R") {
  backend <- .memwas_normalize_engine(engine)

  tmp <- list(formula = formula_i, family = family, data = data_i, id = id, time = time,
              random = random, autocor = autocor, serial = serial,
              L1_penalty = L1_penalty,
              L2_penalty = L2_penalty, control = control, method = screen_method,
              random_cov = random_cov, approximation = approximation, init_approximation = init_approximation, se_method = se_method, engine = backend,
              spline_variables = character(0L), spline_info = list(),
              verbose = FALSE)
  class(tmp) <- "MEMWAS"
  fit_MEMWAS(tmp, verbose = FALSE)
}

.memwas_spline_failure_row <- function(pred, message, n) {
  data.frame(predictor = pred, n = n, knots = NA_character_,
             linear_logLik = NA_real_, spline_logLik = NA_real_,
             linear_AIC = NA_real_, spline_AIC = NA_real_,
             linear_BIC = NA_real_, spline_BIC = NA_real_,
             delta_logLik = NA_real_, delta_AIC = NA_real_, delta_BIC = NA_real_,
             LRT = NA_real_, df = NA_integer_, p_value = NA_real_, nonlinear = FALSE,
             turning_points = NA_character_, convergence_linear = NA_integer_,
             convergence_spline = NA_integer_, message = message,
             stringsAsFactors = FALSE, check.names = FALSE)
}
