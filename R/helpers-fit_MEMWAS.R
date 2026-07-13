# Internal helpers used by fit_MEMWAS(). These functions are intentionally not exported.

.memwas_make_spd_corr <- function(R, eps = 1e-7) {
  R <- as.matrix(R)
  R <- (R + t(R)) / 2
  diag(R) <- 1
  if (nrow(R) <= 1L) return(matrix(1, 1L, 1L))
  eg <- eigen(R, symmetric = TRUE)
  if (min(eg$values) < eps || any(!is.finite(eg$values))) {
    vals <- pmax(eg$values, eps)
    R <- eg$vectors %*% diag(vals, nrow = length(vals)) %*% t(eg$vectors)
    R <- (R + t(R)) / 2
    d <- sqrt(pmax(diag(R), eps))
    R <- R / outer(d, d)
    diag(R) <- 1
  }
  R
}

.memwas_prepare_model_data <- function(formula, data, id, time, random, family) {
  if (!inherits(formula, "formula")) stop("`formula` must be an R formula.", call. = FALSE)
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  if (!is.character(id) || length(id) != 1L || !id %in% names(data)) {
    stop("`id` must be the name of an existing subject-id column in `data`.", call. = FALSE)
  }
  if (!is.character(time) || length(time) != 1L || !time %in% names(data)) {
    stop("`time` must be the name of an existing measurement-time column in `data`.", call. = FALSE)
  }
  random_formula <- .memwas_parse_random_formula(random)
  needed <- unique(c(all.vars(formula), all.vars(random_formula), id, time))
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0L) {
    stop("The following required variable(s) are absent from `data`: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  cc <- stats::complete.cases(data[, needed, drop = FALSE])
  if (!any(cc)) stop("No complete rows remain after removing missing values.", call. = FALSE)
  d <- droplevels(data[cc, , drop = FALSE])
  ord <- order(as.factor(d[[id]]), d[[time]])
  d <- d[ord, , drop = FALSE]
  mf <- stats::model.frame(formula, data = d, na.action = stats::na.pass,
                           drop.unused.levels = TRUE)
  y <- stats::model.response(mf)
  if (is.factor(y)) y <- as.numeric(y) - 1
  if (is.logical(y)) y <- as.numeric(y)
  if (is.matrix(y) && family == "binomial" && ncol(y) == 2L) {
    denom <- rowSums(y)
    if (any(denom <= 0)) stop("Binomial two-column response has non-positive trial totals.", call. = FALSE)
    y <- y[, 1L] / denom
    warning("Two-column binomial response was converted to proportions; trial weights are not currently modeled.",
            call. = FALSE)
  }
  y <- as.numeric(y)
  if (any(!is.finite(y))) stop("Response contains non-finite values.", call. = FALSE)
  if (family == "binomial" && (any(y < 0 | y > 1))) {
    stop("For `family = 'binomial'`, the response must be coded as 0/1 or proportions in [0, 1].",
         call. = FALSE)
  }
  if (family %in% c("poisson", "negative_binomial") &&
      (any(y < 0) || any(abs(y - round(y)) > 1e-8))) {
    stop("For `family = 'poisson'` or `family = 'negative_binomial'`, the response must contain non-negative counts.",
         call. = FALSE)
  }
  if (family %in% c("gamma", "exponential") && any(y <= 0)) {
    stop("For `family = 'gamma'` or `family = 'exponential'`, the response must contain strictly positive values.",
         call. = FALSE)
  }
  terms_obj <- attr(mf, "terms")
  X <- stats::model.matrix(terms_obj, data = mf)
  x_assign <- attr(X, "assign")
  term_labels <- attr(terms_obj, "term.labels")
  model_terms <- rep("(Intercept)", ncol(X))
  if (!is.null(x_assign) && length(x_assign) == ncol(X) && length(term_labels) > 0L) {
    non_intercept <- x_assign > 0L
    model_terms[non_intercept] <- term_labels[x_assign[non_intercept]]
  } else {
    model_terms <- colnames(X)
  }
  coefficient_term_map <- data.frame(
    term = colnames(X),
    model_term = model_terms,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  Z <- stats::model.matrix(random_formula, data = d)
  ids <- as.factor(d[[id]])
  groups <- split(seq_along(y), ids, drop = TRUE)
  list(data = d, y = y, X = X, Z = Z, id = ids, time = d[[time]], groups = groups,
       random_formula = random_formula, dropped_rows = sum(!cc), formula = formula,
       family = family, coefficient_term_map = coefficient_term_map)
}

.memwas_normalize_autocor <- function(autocor) {
  if (is.null(autocor)) return(list(type = "NONE", p = 0L, label = "NONE"))
  a <- toupper(gsub("\\s+", "", as.character(autocor)[1L]))
  if (a %in% c("", "NONE", "NULL", "I", "INDEPENDENCE", "INDEPENDENT")) {
    return(list(type = "NONE", p = 0L, label = "NONE"))
  }
  if (a %in% c("AR1", "AR(1)")) return(list(type = "AR1", p = 1L, label = "AR(1)"))
  if (grepl("^AR\\([0-9]+\\)$", a)) {
    p <- as.integer(sub("^AR\\(([0-9]+)\\)$", "\\1", a))
    if (p <= 1L) return(list(type = "AR1", p = 1L, label = "AR(1)"))
    return(list(type = "ARP", p = p, label = paste0("AR(", p, ")")))
  }
  if (a %in% c("ARMA", "ARMA(1,1)", "ARMA11")) return(list(type = "ARMA11", p = 2L, label = "ARMA(1,1)"))
  if (a %in% c("CS", "EXCHANGEABLE", "COMPOUNDSYMMETRY", "COMPOUND_SYMMETRY")) {
    return(list(type = "CS", p = 1L, label = "CS"))
  }
  if (a %in% c("TOEP", "TOEPLITZ")) return(list(type = "TOEP", p = NA_integer_, label = "TOEP"))
  if (a %in% c("UN", "UNSTRUCTURED")) return(list(type = "UN", p = NA_integer_, label = "UN"))
  stop("Unsupported autocorrelation structure: ", autocor, call. = FALSE)
}

.memwas_build_autocor_info <- function(autocor, groups, control) {
  ac <- .memwas_normalize_autocor(autocor)
  max_group <- max(vapply(groups, length, integer(1L)))
  max_lag <- max(0L, max_group - 1L)
  continuous_time <- isTRUE(control$continuous_time)
  max_un_dim <- control$max_unstructured_dim %||% 20L
  if (ac$type == "NONE") npar <- 0L
  if (ac$type == "AR1") npar <- 1L
  if (ac$type == "ARP") npar <- ac$p
  if (ac$type == "ARMA11") npar <- 2L
  if (ac$type == "CS") npar <- 1L
  if (ac$type == "TOEP") {
    npar <- as.integer(control$toep_lags %||% max_lag)
    npar <- max(0L, min(npar, max_lag))
  }
  if (ac$type == "UN") {
    if (max_group > max_un_dim) {
      stop("`autocor = 'UN'` requires ", max_group, " within-subject positions; this exceeds ",
           "`control$max_unstructured_dim = ", max_un_dim, "`. Increase this control value or choose a parsimonious structure.",
           call. = FALSE)
    }
    npar <- max_group * (max_group - 1L) / 2L
  }
  list(type = ac$type, label = ac$label, p = ac$p, npar = npar, max_lag = max_lag,
       max_group = max_group, continuous_time = continuous_time)
}

.memwas_make_R <- function(info, raw, ti, engine = "R") {
  backend <- .memwas_normalize_engine(engine)

  n <- length(ti)
  if (n <= 1L || info$type == "NONE" || info$npar == 0L) return(diag(1, n))
  lag <- .memwas_lag_matrix(ti, continuous = info$continuous_time, engine = backend)
  if (info$type == "AR1") {
    if (info$continuous_time) {
      rho <- .memwas_bounded_logistic(raw[1L], low = 0, high = 0.99, engine = backend)
    } else {
      rho <- tanh(raw[1L]) * 0.99
    }
    R <- rho^lag
    return(.memwas_make_spd_corr(R))
  }
  if (info$type == "ARP") {
    pacf <- tanh(raw[seq_len(info$p)]) * 0.95
    phi <- .memwas_pacf_to_ar(pacf, engine = backend)
    acf <- .memwas_ar_acf(phi, max_lag = info$max_lag, engine = backend)
    lag_int <- pmin(round(lag), info$max_lag)
    R <- matrix(acf[lag_int + 1L], n, n)
    return(.memwas_make_spd_corr(R))
  }
  if (info$type == "ARMA11") {
    phi <- tanh(raw[1L]) * 0.95
    theta <- tanh(raw[2L]) * 0.95
    gamma0 <- (1 + theta^2 + 2 * phi * theta) / max(1 - phi^2, 1e-8)
    rho1 <- (phi * gamma0 + theta) / max(gamma0, 1e-8)
    acf <- numeric(info$max_lag + 1L)
    acf[1L] <- 1
    if (info$max_lag >= 1L) acf[2L] <- rho1
    if (info$max_lag >= 2L) {
      for (k in 2L:info$max_lag) acf[k + 1L] <- (phi^(k - 1L)) * rho1
    }
    lag_int <- pmin(round(lag), info$max_lag)
    R <- matrix(pmax(pmin(acf[lag_int + 1L], 0.999), -0.999), n, n)
    return(.memwas_make_spd_corr(R))
  }
  if (info$type == "CS") {
    lower <- if (info$max_group > 1L) -1 / (info$max_group - 1L) + 1e-5 else -0.99
    rho <- .memwas_bounded_logistic(raw[1L], low = lower, high = 0.99, engine = backend)
    R <- matrix(rho, n, n)
    diag(R) <- 1
    return(.memwas_make_spd_corr(R))
  }
  if (info$type == "TOEP") {
    rho <- tanh(raw[seq_len(info$npar)]) * 0.95
    lag_int <- round(lag)
    vals <- rep(0, length(lag_int))
    ok <- lag_int > 0L & lag_int <= length(rho)
    vals[lag_int == 0L] <- 1
    vals[ok] <- rho[lag_int[ok]]
    R <- matrix(vals, n, n)
    diag(R) <- 1
    return(.memwas_make_spd_corr(R))
  }
  if (info$type == "UN") {
    m <- info$max_group
    L <- diag(1, m)
    low <- which(lower.tri(L), arr.ind = TRUE)
    if (nrow(low) > 0L) L[low] <- raw[seq_len(nrow(low))]
    C <- L %*% t(L)
    d <- sqrt(pmax(diag(C), 1e-10))
    C <- C / outer(d, d)
    diag(C) <- 1
    R <- C[seq_len(n), seq_len(n), drop = FALSE]
    return(.memwas_make_spd_corr(R))
  }
  stop("Internal error: unknown autocorrelation type.", call. = FALSE)
}

.memwas_describe_autocor <- function(info, raw, engine = "R") {
  backend <- .memwas_normalize_engine(engine)

  if (info$type == "NONE" || info$npar == 0L) return(list(type = "NONE"))
  if (info$type == "AR1") {
    return(list(type = "AR(1)", rho = if (info$continuous_time) .memwas_bounded_logistic(raw[1L], low = 0, high = 0.99, engine = backend) else tanh(raw[1L]) * 0.99))
  }
  if (info$type == "ARP") {
    pacf <- tanh(raw[seq_len(info$p)]) * 0.95
    return(list(type = paste0("AR(", info$p, ")"), partial_autocorrelations = pacf,
                ar_coefficients = .memwas_pacf_to_ar(pacf, engine = backend)))
  }
  if (info$type == "ARMA11") return(list(type = "ARMA(1,1)", phi = tanh(raw[1L]) * 0.95,
                                         theta = tanh(raw[2L]) * 0.95))
  if (info$type == "CS") {
    lower <- if (info$max_group > 1L) -1 / (info$max_group - 1L) + 1e-5 else -0.99
    return(list(type = "CS", rho = .memwas_bounded_logistic(raw[1L], low = lower, high = 0.99, engine = backend)))
  }
  if (info$type == "TOEP") return(list(type = "TOEP", lag_correlations = tanh(raw[seq_len(info$npar)]) * 0.95))
  if (info$type == "UN") {
    m <- info$max_group
    L <- diag(1, m)
    low <- which(lower.tri(L), arr.ind = TRUE)
    if (nrow(low) > 0L) L[low] <- raw[seq_len(nrow(low))]
    C <- L %*% t(L)
    d <- sqrt(pmax(diag(C), 1e-10))
    C <- C / outer(d, d)
    diag(C) <- 1
    return(list(type = "UN", correlation = C))
  }
  list(type = info$type, raw = raw)
}

.memwas_random_npar <- function(q, random_cov) {
  if (q == 0L) return(0L)
  if (random_cov == "diagonal") return(q)
  q * (q + 1L) / 2L
}

.memwas_build_random_cov <- function(par, q, random_cov, engine = "R") {
  backend <- .memwas_normalize_engine(engine)

  if (q == 0L) return(matrix(0, 0L, 0L))
  if (random_cov == "diagonal") {
    return(.memwas_diag_vec(pmax(.memwas_safe_exp(par[seq_len(q)], engine = backend), 1e-8), engine = backend))
  }
  L <- matrix(0, q, q)
  low <- which(lower.tri(L, diag = TRUE), arr.ind = TRUE)
  for (k in seq_len(nrow(low))) {
    i <- low[k, 1L]
    j <- low[k, 2L]
    L[i, j] <- if (i == j) pmax(.memwas_safe_exp(par[k], engine = backend), 1e-8) else par[k]
  }
  D <- L %*% t(L)
  (D + t(D)) / 2
}

.memwas_extract_optim_control <- function(control) {
  if (!is.list(control)) return(list(maxit = 200L))
  if (!is.null(control$optim_control)) return(control$optim_control)
  custom <- c("continuous_time", "max_unstructured_dim", "toep_lags", "optim_method",
              "optim_control", "nongaussian_optim_control", "approximation_optim_control",
              "nongaussian_optim_method", "approximation_optim_method",
              "pql_maxit", "pql_tol", "cd_maxit", "cd_tol",
              "vi_maxit", "vi_tol", "laplace_maxit", "laplace_tol",
              "saddlepoint_maxit", "saddlepoint_tol",
              "skew_laplace_maxit", "skew_laplace_tol",
              "agq_maxit", "agq_tol", "aghq_maxit", "aghq_tol",
              "agq_nodes", "aghq_nodes", "quadrature_nodes",
              "agq_max_dim", "aghq_max_dim", "agq_max_nodes", "aghq_max_nodes",
              "mode_maxit", "conditional_mode_maxit", "mode_tol", "conditional_mode_tol",
              "hessian_eps", "se_method", "bootstrap_B", "parametric_bootstrap_B",
              "profile_ci", "profile_likelihood_ci", "profile_alpha", "profile_grid",
              "estimate_dispersion", "estimate_nb_theta", "estimate_gamma_shape",
              "negative_binomial_theta", "nb_theta", "theta", "size",
              "gamma_shape", "shape", "gamma_k", "gamma_alpha",
              "random_effect_ci_level", "blup_ci_level",
              "max_toep_lag", "spline_grid", "min_unique_nonlinear")
  out <- control[setdiff(names(control), custom)]
  if (is.null(out$maxit)) out$maxit <- 200L
  out
}

.memwas_estimate_beta <- function(Xw, yw, lambda1, lambda2, xnames, control,
                                  engine = "R") {
  backend <- .memwas_normalize_engine(engine)

  pf <- .memwas_penalty_factor(xnames, engine = backend)
  if (lambda1 > 0) {
    return(.memwas_coordinate_descent_enet(
      X = Xw,
      y = yw,
      lambda1 = lambda1,
      lambda2 = lambda2,
      pf = pf,
      maxit = control$cd_maxit %||% 1000L,
      tol = control$cd_tol %||% 1e-7,
      engine = backend
    ))
  }
  H <- crossprod(Xw) + .memwas_diag_vec(lambda2 * pf + 1e-10, engine = backend)
  as.vector(.memwas_safe_solve(H, crossprod(Xw, yw), engine = backend))
}

.memwas_autocor_penalty_value <- function(auto_raw, control) {
  if (length(auto_raw) == 0L || is.null(control) || !is.list(control)) return(0)
  reg <- control$autocor_regularization %||% NULL
  lambda <- control$autocor_penalty %||% 0
  type <- control$autocor_regularization_type %||% "L2"
  alpha <- control$autocor_regularization_alpha %||% 0
  enabled <- TRUE

  if (is.list(reg)) {
    enabled <- isTRUE(reg$enabled %||% TRUE)
    lambda <- reg$lambda %||% lambda
    type <- reg$type %||% type
    alpha <- reg$alpha %||% alpha
  } else if (is.numeric(reg) && length(reg) >= 1L) {
    lambda <- reg[1L]
  }

  if (!isTRUE(enabled)) return(0)
  lambda <- suppressWarnings(as.numeric(lambda)[1L])
  alpha <- suppressWarnings(as.numeric(alpha)[1L])
  type <- tolower(as.character(type)[1L])
  if (!is.finite(lambda) || lambda <= 0 || type == "none") return(0)
  if (!is.finite(alpha)) alpha <- 0
  alpha <- pmin(pmax(alpha, 0), 1)

  if (type == "l1") return(lambda * sum(abs(auto_raw)))
  if (type == "elasticnet") {
    return(lambda * (alpha * sum(abs(auto_raw)) + 0.5 * (1 - alpha) * sum(auto_raw^2)))
  }
  0.5 * lambda * sum(auto_raw^2)
}

.memwas_fit_gaussian_lmm <- function(prep, autocor, random_cov, L1_penalty, L2_penalty,
                                     control, method = "ML", obs_var = NULL,
                                     theta_start = NULL, working_label = "gaussian",
                                     engine = "R") {
  backend <- .memwas_normalize_engine(engine)

  y <- prep$y
  X <- prep$X
  Z <- prep$Z
  groups <- prep$groups
  time <- prep$time
  n <- length(y)
  p <- ncol(X)
  q <- ncol(Z)
  if (p == 0L) stop("The fixed-effects design matrix has zero columns.", call. = FALSE)
  if (is.null(obs_var)) obs_var <- rep(0, n)
  obs_var <- pmax(as.numeric(obs_var), 0)
  random_cov <- match.arg(random_cov, c("diagonal", "unstructured"))
  if ((L1_penalty > 0 || L2_penalty > 0) && toupper(method) == "REML") method <- "ML"

  ac_info <- .memwas_build_autocor_info(autocor, groups, control)
  rn <- .memwas_random_npar(q, random_cov)
  s2_init <- stats::var(y)
  lm0 <- try(stats::lm.fit(x = X, y = y), silent = TRUE)
  if (!inherits(lm0, "try-error") && !is.null(lm0$residuals)) {
    s2_lm <- stats::var(as.numeric(lm0$residuals))
    if (is.finite(s2_lm) && s2_lm > 0) s2_init <- s2_lm
  }
  if (!is.finite(s2_init) || s2_init <= 0) s2_init <- 1
  rand0 <- if (rn == 0L) numeric(0L) else {
    if (random_cov == "diagonal") rep(log(max(0.10 * s2_init, 1e-4)), q) else {
      vals <- numeric(rn)
      low <- which(lower.tri(matrix(0, q, q), diag = TRUE), arr.ind = TRUE)
      for (k in seq_len(nrow(low))) vals[k] <- if (low[k, 1L] == low[k, 2L]) log(sqrt(max(0.10 * s2_init, 1e-4))) else 0
      vals
    }
  }
  theta0 <- c(rand0, log(sqrt(max(0.90 * s2_init, 1e-4))), rep(0, ac_info$npar))
  if (!is.null(theta_start) && length(theta_start) == length(theta0) && all(is.finite(theta_start))) {
    theta0 <- theta_start
  }
  idx_random <- if (rn > 0L) seq_len(rn) else integer(0L)
  idx_sigma <- rn + 1L
  idx_auto <- if (ac_info$npar > 0L) seq(from = rn + 2L, length.out = ac_info$npar) else integer(0L)
  opt_control <- .memwas_extract_optim_control(control)
  opt_method <- control$optim_method %||% "BFGS"

  evaluate_theta <- function(theta, return_all = FALSE) {
    D <- .memwas_build_random_cov(theta[idx_random], q, random_cov, engine = backend)
    sigma <- pmax(.memwas_safe_exp(theta[idx_sigma], engine = backend), 1e-8)
    sigma2 <- sigma^2
    auto_raw <- if (length(idx_auto)) theta[idx_auto] else numeric(0L)
    Xw_list <- vector("list", length(groups))
    yw_list <- vector("list", length(groups))
    chol_list <- vector("list", length(groups))
    R_list <- vector("list", length(groups))
    V_list <- vector("list", length(groups))
    logdetV <- 0
    gnames <- names(groups)
    for (g in seq_along(groups)) {
      ii <- groups[[g]]
      Zi <- Z[ii, , drop = FALSE]
      Xi <- X[ii, , drop = FALSE]
      Ri <- .memwas_make_R(ac_info, auto_raw, time[ii], engine = backend)
      Vi <- sigma2 * Ri
      if (q > 0L) Vi <- Vi + Zi %*% D %*% t(Zi)
      if (any(obs_var[ii] > 0)) Vi <- Vi + .memwas_diag_vec(obs_var[ii], engine = backend)
      Ui <- .memwas_safe_chol(Vi, engine = backend)
      logdetV <- logdetV + 2 * sum(log(diag(Ui)))
      Xw_list[[g]] <- forwardsolve(t(Ui), Xi)
      yw_list[[g]] <- as.vector(forwardsolve(t(Ui), matrix(y[ii], ncol = 1L)))
      chol_list[[g]] <- Ui
      R_list[[g]] <- Ri
      V_list[[g]] <- Vi
    }
    names(chol_list) <- gnames
    names(R_list) <- gnames
    names(V_list) <- gnames
    Xw <- do.call(rbind, Xw_list)
    yw <- unlist(yw_list, use.names = FALSE)
    beta <- .memwas_estimate_beta(Xw, yw, L1_penalty, L2_penalty, colnames(X), control, engine = backend)
    resw <- as.vector(yw - Xw %*% beta)
    rss <- sum(resw^2)
    pf <- .memwas_penalty_factor(colnames(X), engine = backend)
    pen <- .memwas_fixed_effect_penalty(beta, pf, L1_penalty, L2_penalty, engine = backend)
    autocor_pen <- .memwas_autocor_penalty_value(auto_raw, control)
    if (toupper(method) == "REML" && L1_penalty == 0 && L2_penalty == 0) {
      XtVX <- crossprod(Xw)
      rankX <- qr(XtVX)$rank
      Ux <- .memwas_safe_chol(XtVX + diag(1e-10, ncol(XtVX)), engine = backend)
      logdetX <- 2 * sum(log(diag(Ux)))
      nll_unpen <- 0.5 * ((n - rankX) * log(2 * pi) + logdetV + logdetX + rss)
    } else {
      rankX <- qr(X)$rank
      logdetX <- NA_real_
      nll_unpen <- 0.5 * (n * log(2 * pi) + logdetV + rss)
    }
    nll <- nll_unpen + pen + autocor_pen
    if (!is.finite(nll)) nll <- 1e100
    if (!return_all) return(nll)
    list(nll = nll, nll_unpen = nll_unpen, penalty = pen, autocor_penalty = autocor_pen, beta = beta, D = D,
         sigma = sigma, sigma2 = sigma2, auto_raw = auto_raw, ac_info = ac_info,
         logdetV = logdetV, logdetX = logdetX, rss = rss, rankX = rankX,
         Xw = Xw, yw = yw, chol = chol_list, R = R_list, V = V_list)
  }

  obj <- function(theta) {
    val <- try(evaluate_theta(theta, return_all = FALSE), silent = TRUE)
    if (inherits(val, "try-error") || length(val) != 1L || !is.finite(val)) return(1e100)
    val
  }

  opt <- try(stats::optim(theta0, obj, method = opt_method, control = opt_control), silent = TRUE)
  if (inherits(opt, "try-error") || !is.finite(opt$value)) {
    opt <- try(stats::optim(theta0, obj, method = "Nelder-Mead", control = opt_control), silent = TRUE)
  }
  if (inherits(opt, "try-error")) {
    stop("Optimization failed in the mixed-effects covariance model.", call. = FALSE)
  }

  final <- evaluate_theta(opt$par, return_all = TRUE)
  beta <- final$beta
  names(beta) <- colnames(X)
  pf <- .memwas_penalty_factor(colnames(X), engine = backend)
  Hbeta <- crossprod(final$Xw) + .memwas_diag_vec(L2_penalty * pf + 1e-10, engine = backend)
  vcov_beta <- .memwas_safe_solve(Hbeta, engine = backend)
  dimnames(vcov_beta) <- list(colnames(X), colnames(X))
  se <- sqrt(pmax(diag(vcov_beta), 0))
  zval <- beta / se
  pval <- 2 * stats::pnorm(abs(zval), lower.tail = FALSE)
  coef_table <- data.frame(term = names(beta), estimate = as.numeric(beta), std_error = se,
                           statistic = as.numeric(zval), p_value = as.numeric(pval),
                           row.names = NULL, check.names = FALSE)

  fitted_link <- as.vector(X %*% beta)
  residuals <- y - fitted_link
  ci_level <- control$random_effect_ci_level %||% control$blup_ci_level %||% 0.95
  ci_level <- suppressWarnings(as.numeric(ci_level)[1L])
  if (!is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) {
    stop("`control$random_effect_ci_level` must be a numeric value strictly between 0 and 1.", call. = FALSE)
  }
  ci_mult <- stats::qnorm(0.5 + ci_level / 2)
  blup <- vector("list", length(groups))
  blup_se <- vector("list", length(groups))
  blup_ci <- vector("list", length(groups))
  blup_cov <- vector("list", length(groups))
  names(blup) <- names(groups)
  names(blup_se) <- names(groups)
  names(blup_ci) <- names(groups)
  names(blup_cov) <- names(groups)
  blup_rows <- vector("list", length(groups))
  if (q > 0L) {
    for (g in seq_along(groups)) {
      ii <- groups[[g]]
      Zi <- Z[ii, , drop = FALSE]
      Xi <- X[ii, , drop = FALSE]
      ri <- matrix(y[ii] - Xi %*% beta, ncol = 1L)
      Ui <- final$chol[[g]]
      vinv_r <- backsolve(Ui, forwardsolve(t(Ui), ri))
      bi <- final$D %*% t(Zi) %*% vinv_r
      bi <- as.numeric(bi)
      names(bi) <- colnames(Z)

      vinv_ZD <- backsolve(Ui, forwardsolve(t(Ui), Zi %*% final$D))
      cond_cov <- final$D - final$D %*% t(Zi) %*% vinv_ZD
      cond_cov <- (cond_cov + t(cond_cov)) / 2
      dimnames(cond_cov) <- list(colnames(Z), colnames(Z))
      se_i <- sqrt(pmax(diag(cond_cov), 0))
      names(se_i) <- colnames(Z)
      ci_i <- cbind(conf_low = bi - ci_mult * se_i,
                    conf_high = bi + ci_mult * se_i)
      rownames(ci_i) <- colnames(Z)

      blup[[g]] <- bi
      blup_se[[g]] <- se_i
      blup_ci[[g]] <- ci_i
      blup_cov[[g]] <- cond_cov
      blup_rows[[g]] <- data.frame(
        id = rep(names(groups)[g], length(bi)),
        effect = colnames(Z),
        estimate = as.numeric(bi),
        std_error = as.numeric(se_i),
        conf_low = as.numeric(ci_i[, "conf_low"]),
        conf_high = as.numeric(ci_i[, "conf_high"]),
        ci_level = rep(ci_level, length(bi)),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }
  blup_rows <- blup_rows[!vapply(blup_rows, is.null, logical(1L))]
  blup_table <- if (length(blup_rows)) {
    out_re <- do.call(rbind, blup_rows)
    row.names(out_re) <- NULL
    out_re
  } else {
    data.frame(id = character(0L), effect = character(0L), estimate = numeric(0L),
               std_error = numeric(0L), conf_low = numeric(0L), conf_high = numeric(0L),
               ci_level = numeric(0L), stringsAsFactors = FALSE, check.names = FALSE)
  }

  df_fixed <- if (L1_penalty > 0) sum(abs(beta) > 1e-8) else qr(X)$rank
  df_total <- df_fixed + length(opt$par)
  logLik <- -final$nll_unpen
  metrics <- list(logLik = as.numeric(logLik), penalized_logLik = as.numeric(-final$nll),
                  AIC = as.numeric(-2 * logLik + 2 * df_total),
                  BIC = as.numeric(-2 * logLik + log(n) * df_total),
                  deviance = as.numeric(-2 * logLik), df = as.numeric(df_total),
                  nobs = n, method = toupper(method), penalty = final$penalty,
                  autocorrelation_penalty = final$autocor_penalty %||% 0,
                  working_label = working_label)

  list(call = match.call(), formula = prep$formula, data = prep$data, family = prep$family,
       y = y, X = X, Z = Z, groups = groups, coefficients = beta,
       coefficient_table = coef_table, coefficient_term_map = prep$coefficient_term_map,
       vcov = vcov_beta, fitted = fitted_link,
       residuals = residuals, random_effects = blup,
       random_effects_se = blup_se, random_effects_ci = blup_ci,
       random_effects_covariance = blup_cov, random_effects_table = blup_table,
       random_effects_ci_level = ci_level,
       random_effects_note = "BLUP standard errors and confidence intervals use the conditional random-effect covariance D - D Z' V^{-1} Z D.",
       random_covariance = final$D,
       residual_sigma = final$sigma, autocorrelation = .memwas_describe_autocor(final$ac_info, final$auto_raw, engine = backend),
       raw_theta = opt$par, metrics = metrics, convergence = opt$convergence,
       optim_message = opt$message %||% "", approximate = FALSE,
       dropped_rows = prep$dropped_rows)
}

.memwas_control_positive <- function(control, names, default, lower = .Machine$double.eps) {
  control <- control %||% list()
  for (nm in names) {
    if (!is.null(control[[nm]])) {
      val <- suppressWarnings(as.numeric(control[[nm]][1L]))
      if (is.finite(val) && val > lower) return(val)
    }
  }
  default
}

.memwas_estimate_nb_theta <- function(y, control) {
  val <- .memwas_control_positive(control,
                                  c("negative_binomial_theta", "nb_theta", "theta", "size"),
                                  default = NA_real_, lower = 1e-8)
  if (is.finite(val)) return(val)
  y <- as.numeric(y)
  m <- mean(y, na.rm = TRUE)
  v <- stats::var(y, na.rm = TRUE)
  if (is.finite(m) && is.finite(v) && m > 0 && v > m) {
    return(pmax(m^2 / (v - m), 1e-4))
  }
  1000
}

.memwas_estimate_gamma_shape <- function(y, control, default = NA_real_) {
  val <- .memwas_control_positive(control,
                                  c("gamma_shape", "shape", "gamma_k", "gamma_alpha"),
                                  default = default, lower = 1e-8)
  if (is.finite(val)) return(val)
  y <- as.numeric(y)
  m <- mean(y, na.rm = TRUE)
  v <- stats::var(y, na.rm = TRUE)
  if (is.finite(m) && is.finite(v) && m > 0 && v > 0) {
    return(pmax(m^2 / v, 1e-4))
  }
  1
}

.memwas_family_parts <- function(family, y = NULL, control = list()) {
  family <- .memwas_normalize_family(family)
  clamp_eta <- function(eta) pmin(pmax(eta, -30), 30)
  safe_mu <- function(eta) pmax(exp(clamp_eta(eta)), 1e-8)
  loglik_sum <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    sum(x)
  }

  if (family == "binomial") {
    return(list(
      family = family,
      link = "logit",
      parameters = list(),
      linkinv = function(eta) pmin(pmax(stats::plogis(clamp_eta(eta)), 1e-6), 1 - 1e-6),
      mu_eta = function(eta) {
        mu <- pmin(pmax(stats::plogis(clamp_eta(eta)), 1e-6), 1 - 1e-6)
        pmax(mu * (1 - mu), 1e-8)
      },
      variance = function(mu) pmax(mu * (1 - mu), 1e-8),
      weight = function(mu) pmax(mu * (1 - mu), 1e-8),
      skewness = function(mu) (1 - 2 * mu) / sqrt(pmax(mu * (1 - mu), 1e-8)),
      excess_kurtosis = function(mu) (1 - 6 * mu * (1 - mu)) / pmax(mu * (1 - mu), 1e-8),
      loglik = function(y, mu) loglik_sum(stats::dbinom(round(y), size = 1, prob = mu, log = TRUE))
    ))
  }

  if (family == "poisson") {
    return(list(
      family = family,
      link = "log",
      parameters = list(),
      linkinv = safe_mu,
      mu_eta = safe_mu,
      variance = function(mu) pmax(mu, 1e-8),
      weight = function(mu) pmax(mu, 1e-8),
      skewness = function(mu) 1 / sqrt(pmax(mu, 1e-8)),
      excess_kurtosis = function(mu) 1 / pmax(mu, 1e-8),
      loglik = function(y, mu) loglik_sum(stats::dpois(round(y), lambda = pmax(mu, 1e-8), log = TRUE))
    ))
  }

  if (family == "negative_binomial") {
    theta <- .memwas_estimate_nb_theta(y, control)
    return(list(
      family = family,
      link = "log",
      parameters = list(theta = theta, size = theta),
      linkinv = safe_mu,
      mu_eta = safe_mu,
      variance = function(mu) pmax(mu + mu^2 / theta, 1e-8),
      weight = function(mu) pmax(mu^2 / (mu + mu^2 / theta), 1e-8),
      skewness = function(mu) (1 + 2 * mu / theta) / sqrt(pmax(mu + mu^2 / theta, 1e-8)),
      excess_kurtosis = function(mu) 6 / pmax(theta, 1e-8) +
        (1 / pmax(mu, 1e-8)) * (1 + 6 * mu / pmax(theta, 1e-8) + 6 * mu^2 / pmax(theta, 1e-8)^2) /
        pmax((1 + mu / pmax(theta, 1e-8))^2, 1e-8),
      loglik = function(y, mu) loglik_sum(stats::dnbinom(round(y), size = theta, mu = pmax(mu, 1e-8), log = TRUE))
    ))
  }

  if (family == "gamma") {
    shape <- .memwas_estimate_gamma_shape(y, control)
    return(list(
      family = family,
      link = "log",
      parameters = list(shape = shape),
      linkinv = safe_mu,
      mu_eta = safe_mu,
      variance = function(mu) pmax(mu^2 / shape, 1e-8),
      weight = function(mu) rep(shape, length(mu)),
      skewness = function(mu) rep(2 / sqrt(pmax(shape, 1e-8)), length(mu)),
      excess_kurtosis = function(mu) rep(6 / pmax(shape, 1e-8), length(mu)),
      loglik = function(y, mu) loglik_sum(stats::dgamma(y, shape = shape, rate = shape / pmax(mu, 1e-8), log = TRUE))
    ))
  }

  if (family == "exponential") {
    return(list(
      family = family,
      link = "log",
      parameters = list(shape = 1),
      linkinv = safe_mu,
      mu_eta = safe_mu,
      variance = function(mu) pmax(mu^2, 1e-8),
      weight = function(mu) rep(1, length(mu)),
      skewness = function(mu) rep(2, length(mu)),
      excess_kurtosis = function(mu) rep(6, length(mu)),
      loglik = function(y, mu) loglik_sum(stats::dexp(y, rate = 1 / pmax(mu, 1e-8), log = TRUE))
    ))
  }

  stop("Internal error: unsupported non-Gaussian family.", call. = FALSE)
}

.memwas_initial_glm_beta <- function(X, y, family, control = list()) {
  family <- .memwas_normalize_family(family)
  fam <- switch(family,
                binomial = stats::binomial(),
                poisson = stats::poisson(),
                negative_binomial = stats::poisson(),
                gamma = stats::Gamma(link = "log"),
                exponential = stats::Gamma(link = "log"),
                stats::gaussian())
  out <- try(stats::glm.fit(x = X, y = y, family = fam, control = list(maxit = 50)), silent = TRUE)
  if (!inherits(out, "try-error") && length(out$coefficients) == ncol(X) && all(is.finite(out$coefficients))) {
    beta <- out$coefficients
    beta[!is.finite(beta)] <- 0
    return(as.numeric(beta))
  }
  beta <- rep(0, ncol(X))
  if ("(Intercept)" %in% colnames(X)) {
    idx <- match("(Intercept)", colnames(X))
    if (family == "binomial") beta[idx] <- stats::qlogis(pmin(pmax(mean(y), 1e-4), 1 - 1e-4))
    if (family %in% c("poisson", "negative_binomial", "gamma", "exponential")) {
      beta[idx] <- log(pmax(mean(y), 1e-4))
    }
  }
  beta
}

.memwas_logdet_spd <- function(M, engine = "R") {
  backend <- .memwas_normalize_engine(engine)

  M <- as.matrix(M)
  if (nrow(M) == 0L) return(0)
  U <- try(.memwas_safe_chol((M + t(M)) / 2, engine = backend), silent = TRUE)
  if (inherits(U, "try-error")) {
    M <- (M + t(M)) / 2 + diag(1e-6, nrow(M))
    U <- chol(M)
  }
  2 * sum(log(pmax(diag(U), 1e-300)))
}

.memwas_logsumexp <- function(x) {
  x <- as.numeric(x)
  if (!length(x)) return(-Inf)
  m <- max(x, na.rm = TRUE)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

.memwas_gauss_hermite_nodes <- function(n = 7L) {
  n <- as.integer(n[1L])
  if (!is.finite(n) || n < 3L) n <- 7L
  if (n > 25L) n <- 25L
  J <- matrix(0, n, n)
  if (n > 1L) {
    for (i in seq_len(n - 1L)) {
      val <- sqrt(i / 2)
      J[i, i + 1L] <- val
      J[i + 1L, i] <- val
    }
  }
  eg <- eigen(J, symmetric = TRUE)
  ord <- order(eg$values)
  nodes <- as.numeric(eg$values[ord])
  weights <- as.numeric(sqrt(pi) * eg$vectors[1L, ord]^2)
  list(nodes = nodes, weights = weights)
}

.memwas_glmm_group_matrices <- function(prep, fit, group_index, beta, parts, engine = "R") {
  ii <- prep$groups[[group_index]]
  Xi <- prep$X[ii, , drop = FALSE]
  Zi <- prep$Z[ii, , drop = FALSE]
  yi <- prep$y[ii]
  q <- ncol(Zi)
  b <- if (q > 0L && !is.null(fit$random_effects[[group_index]])) {
    as.numeric(fit$random_effects[[group_index]])
  } else rep(0, q)
  eta <- as.vector(Xi %*% beta + if (q > 0L) Zi %*% b else 0)
  mu <- parts$linkinv(eta)
  w <- pmax(as.numeric(parts$weight(mu)), 1e-8)
  list(ii = ii, X = Xi, Z = Zi, y = yi, b = b, eta = eta, mu = mu, w = w)
}

.memwas_glmm_logpost <- function(y, X, Z, beta, b, D, Dinv, logdetD, parts) {
  q <- length(b)
  eta <- as.vector(X %*% beta + if (q > 0L) Z %*% b else 0)
  mu <- parts$linkinv(eta)
  cond <- parts$loglik(y, mu)
  if (q == 0L) return(cond)
  prior <- -0.5 * (q * log(2 * pi) + logdetD + as.numeric(crossprod(b, Dinv %*% b)))
  cond + prior
}

.memwas_agq_group_loglik <- function(gm, beta, D, Dinv, logdetD, parts,
                                     n_nodes = 7L, max_nodes = 50000L,
                                     engine = "R") {
  backend <- .memwas_normalize_engine(engine)

  q <- ncol(gm$Z)
  if (q == 0L) return(list(logLik = parts$loglik(gm$y, gm$mu), used = TRUE))
  H <- Dinv + crossprod(gm$Z, gm$Z * gm$w)
  S <- try(.memwas_safe_solve(H + diag(1e-8, q), engine = backend), silent = TRUE)
  if (inherits(S, "try-error") || any(!is.finite(S))) {
    return(list(logLik = NA_real_, used = FALSE))
  }
  U <- try(.memwas_safe_chol((S + t(S)) / 2, engine = backend), silent = TRUE)
  if (inherits(U, "try-error")) return(list(logLik = NA_real_, used = FALSE))
  gh <- .memwas_gauss_hermite_nodes(n_nodes)
  total_nodes <- length(gh$nodes)^q
  if (total_nodes > max_nodes) return(list(logLik = NA_real_, used = FALSE))
  idx_grid <- do.call(expand.grid, c(rep(list(seq_along(gh$nodes)), q), KEEP.OUT.ATTRS = FALSE))
  log_terms <- numeric(nrow(idx_grid))
  for (r in seq_len(nrow(idx_grid))) {
    idx <- as.integer(idx_grid[r, ])
    z <- gh$nodes[idx]
    b <- gm$b + as.vector(sqrt(2) * t(U) %*% z)
    logw <- sum(log(pmax(gh$weights[idx], 1e-300)))
    log_terms[r] <- logw + sum(z^2) +
      .memwas_glmm_logpost(gm$y, gm$X, gm$Z, beta, b, D, Dinv, logdetD, parts)
  }
  logdetL <- sum(log(pmax(diag(U), 1e-300)))
  list(logLik = 0.5 * q * log(2) + logdetL + .memwas_logsumexp(log_terms), used = TRUE)
}

.memwas_glmm_approx_loglik <- function(fit, prep, family, parts, approximation,
                                       control = list(), engine = "R") {
  backend <- .memwas_normalize_engine(engine)

  approximation <- .memwas_normalize_approximation(approximation)
  beta <- as.numeric(fit$coefficients)
  names(beta) <- colnames(prep$X)
  q <- ncol(prep$Z)
  if (q == 0L) {
    eta <- as.vector(prep$X %*% beta)
    mu <- parts$linkinv(eta)
    ll <- parts$loglik(prep$y, mu)
    return(list(logLik = as.numeric(ll), approximation = approximation,
                label = .memwas_approximation_label(approximation),
                n_groups = length(prep$groups), n_quadrature_nodes = 0L,
                fallback_groups = 0L,
                note = "No random effects were present; family log-likelihood was evaluated at the fixed-effect linear predictor.",
                details = data.frame()))
  }
  D <- fit$random_covariance
  if (is.null(D) || any(!is.finite(D))) {
    return(list(logLik = NA_real_, approximation = approximation,
                label = .memwas_approximation_label(approximation),
                n_groups = length(prep$groups), n_quadrature_nodes = NA_integer_,
                fallback_groups = length(prep$groups),
                note = "Approximate family-scale log-likelihood was not computed because the random-effect covariance was unavailable.",
                details = data.frame()))
  }
  D <- as.matrix(D)
  D <- (D + t(D)) / 2 + diag(1e-8, nrow(D))
  Dinv <- .memwas_safe_solve(D, engine = backend)
  logdetD <- .memwas_logdet_spd(D, engine = backend)
  n_nodes <- as.integer(control$aghq_nodes %||% control$agq_nodes %||% control$quadrature_nodes %||% 7L)
  if (!is.finite(n_nodes) || n_nodes < 3L) n_nodes <- 7L
  max_dim <- as.integer(control$aghq_max_dim %||% control$agq_max_dim %||% 2L)
  if (!is.finite(max_dim) || max_dim < 1L) max_dim <- 1L
  max_nodes <- as.integer(control$aghq_max_nodes %||% control$agq_max_nodes %||% 50000L)
  if (!is.finite(max_nodes) || max_nodes < 100L) max_nodes <- 50000L

  rows <- vector("list", length(prep$groups))
  fallback <- 0L
  for (g in seq_along(prep$groups)) {
    gm <- .memwas_glmm_group_matrices(prep, fit, g, beta, parts, engine = backend)
    H <- Dinv + crossprod(gm$Z, gm$Z * gm$w)
    logdetH <- .memwas_logdet_spd(H + diag(1e-8, q), engine = backend)
    S <- .memwas_safe_solve(H + diag(1e-8, q), engine = backend)
    logdetS <- .memwas_logdet_spd(S + diag(1e-10, q), engine = backend)
    cond_ll <- parts$loglik(gm$y, gm$mu)
    prior_ll <- -0.5 * (q * log(2 * pi) + logdetD + as.numeric(crossprod(gm$b, Dinv %*% gm$b)))
    laplace <- cond_ll + prior_ll + 0.5 * q * log(2 * pi) - 0.5 * logdetH

    diag_zsz <- rowSums((gm$Z %*% S) * gm$Z)
    vi_loglik_corr <- -0.5 * sum(gm$w * pmax(diag_zsz, 0))
    prior_expect <- -0.5 * (q * log(2 * pi) + logdetD +
                              as.numeric(crossprod(gm$b, Dinv %*% gm$b)) +
                              sum(Dinv * t(S)))
    entropy <- 0.5 * (q * (1 + log(2 * pi)) + logdetS)
    vi <- cond_ll + vi_loglik_corr + prior_expect + entropy

    skew <- parts$skewness(gm$mu)
    kurt <- parts$excess_kurtosis(gm$mu)
    leverage <- pmin(pmax(diag_zsz, 0), 1)
    saddle_corr <- -0.125 * sum((skew^2) * leverage^2, na.rm = TRUE)
    saddle_corr <- .memwas_clamp(saddle_corr, -10, 10)
    skew_corr <- sum(skew * leverage^(3 / 2), na.rm = TRUE) / 6
    skew_corr <- .memwas_clamp(skew_corr, -10, 10)

    used_agq <- FALSE
    agq <- NA_real_
    if (approximation %in% c("adaptive_gaussian_quadrature", "adaptive_gauss_hermite_quadrature") &&
        q <= max_dim) {
      agq_res <- .memwas_agq_group_loglik(gm, beta, D, Dinv, logdetD, parts,
                                          n_nodes = n_nodes, max_nodes = max_nodes,
                                          engine = backend)
      agq <- agq_res$logLik
      used_agq <- isTRUE(agq_res$used) && is.finite(agq)
    }
    if (approximation %in% c("adaptive_gaussian_quadrature", "adaptive_gauss_hermite_quadrature") && !used_agq) {
      fallback <- fallback + 1L
    }

    value <- switch(approximation,
                    variational_inference = vi,
                    laplace = laplace,
                    saddlepoint = laplace + saddle_corr,
                    skew_corrected_laplace = laplace + skew_corr,
                    adaptive_gaussian_quadrature = if (used_agq) agq else laplace,
                    adaptive_gauss_hermite_quadrature = if (used_agq) agq else laplace,
                    pql = cond_ll,
                    vi)
    rows[[g]] <- data.frame(group = names(prep$groups)[g] %||% as.character(g),
                            logLik = as.numeric(value),
                            conditional_logLik = as.numeric(cond_ll),
                            laplace_logLik = as.numeric(laplace),
                            variational_elbo = as.numeric(vi),
                            saddlepoint_correction = as.numeric(saddle_corr),
                            skew_laplace_correction = as.numeric(skew_corr),
                            quadrature_used = used_agq,
                            stringsAsFactors = FALSE,
                            check.names = FALSE)
  }
  detail <- do.call(rbind, rows)
  finite_ll <- detail$logLik[is.finite(detail$logLik)]
  total <- if (length(finite_ll)) sum(finite_ll) else NA_real_
  note <- switch(approximation,
                 variational_inference = "Mean-field variational lower-bound diagnostic using the working quadratic model around subject-specific modes.",
                 laplace = "Laplace diagnostic using a second-order expansion around subject-specific random-effect modes.",
                 saddlepoint = "Saddlepoint diagnostic applies a bounded cumulant correction to the Laplace diagnostic.",
                 skew_corrected_laplace = "Skew-corrected Laplace diagnostic applies a bounded residual-skewness correction to the Laplace diagnostic.",
                 adaptive_gaussian_quadrature = "Adaptive Gaussian quadrature diagnostic; groups beyond the configured dimension/node limits fall back to Laplace.",
                 adaptive_gauss_hermite_quadrature = "Adaptive Gauss-Hermite quadrature diagnostic; groups beyond the configured dimension/node limits fall back to Laplace.",
                 pql = "Legacy penalized quasi-likelihood diagnostic evaluated at fitted conditional means.",
                 "Approximate family-scale log-likelihood diagnostic.")
  list(logLik = as.numeric(total), approximation = approximation,
       label = .memwas_approximation_label(approximation),
       n_groups = length(prep$groups), n_quadrature_nodes = if (approximation %in% c("adaptive_gaussian_quadrature", "adaptive_gauss_hermite_quadrature")) n_nodes else 0L,
       fallback_groups = fallback, note = note, details = detail)
}

.memwas_fit_pql_glmm <- function(prep, family, autocor, random_cov, L1_penalty, L2_penalty,
                                 control, method = "ML", approximation = "variational_inference",
                                 engine = "R") {
  backend <- .memwas_normalize_engine(engine)

  approximation <- .memwas_normalize_approximation(approximation)
  parts <- .memwas_family_parts(family, y = prep$y, control = control)
  X <- prep$X
  y <- prep$y
  beta <- .memwas_initial_glm_beta(X, y, family, control = control)
  theta_start <- NULL
  maxit_name <- switch(approximation,
                       variational_inference = "vi_maxit",
                       laplace = "laplace_maxit",
                       saddlepoint = "saddlepoint_maxit",
                       skew_corrected_laplace = "skew_laplace_maxit",
                       adaptive_gaussian_quadrature = "agq_maxit",
                       adaptive_gauss_hermite_quadrature = "aghq_maxit",
                       pql = "pql_maxit",
                       "pql_maxit")
  tol_name <- switch(approximation,
                     variational_inference = "vi_tol",
                     laplace = "laplace_tol",
                     saddlepoint = "saddlepoint_tol",
                     skew_corrected_laplace = "skew_laplace_tol",
                     adaptive_gaussian_quadrature = "agq_tol",
                     adaptive_gauss_hermite_quadrature = "aghq_tol",
                     pql = "pql_tol",
                     "pql_tol")
  maxit <- max(1L, as.integer(control[[maxit_name]] %||% control$pql_maxit %||% 20L))
  tol <- control[[tol_name]] %||% control$pql_tol %||% 1e-5
  fit <- NULL
  converged <- FALSE
  for (iter in seq_len(maxit)) {
    eta <- as.vector(X %*% beta)
    mu <- parts$linkinv(eta)
    dmu <- parts$mu_eta(eta)
    varmu <- parts$variance(mu)
    w <- pmax((dmu^2) / varmu, 1e-6)
    z <- eta + (y - mu) / pmax(dmu, 1e-8)
    prep_z <- prep
    prep_z$y <- as.numeric(z)
    fit <- .memwas_fit_gaussian_lmm(prep_z, autocor = autocor, random_cov = random_cov,
                                    L1_penalty = L1_penalty, L2_penalty = L2_penalty,
                                    control = control, method = "ML", obs_var = 1 / w,
                                    theta_start = theta_start,
                                    working_label = paste0(family, "_", approximation),
                                    engine = backend)
    beta_new <- as.numeric(fit$coefficients)
    theta_start <- fit$raw_theta
    if (max(abs(beta_new - beta), na.rm = TRUE) < tol) {
      beta <- beta_new
      converged <- TRUE
      break
    }
    beta <- beta_new
  }
  names(beta) <- colnames(X)
  eta <- as.vector(X %*% beta)
  mu <- parts$linkinv(eta)
  fit$family <- family
  fit$family_link <- parts$link
  fit$family_parameters <- parts$parameters
  fit$coefficients <- beta
  fit$working_fitted <- fit$fitted
  fit$working_residuals <- fit$residuals
  fit$fitted_link <- eta
  fit$fitted_response <- mu
  fit$fitted <- mu
  fit$residuals <- y - mu
  fit$approximation <- approximation
  fit$approximation_label <- .memwas_approximation_label(approximation)
  fit$pql <- list(iterations = iter, converged = converged,
                  note = paste0("Working Gaussian iterations used for ",
                                .memwas_approximation_label(approximation),
                                "; exact high-dimensional GLMM likelihood is not directly optimized."))
  fit$approximation_details <- .memwas_glmm_approx_loglik(fit, prep, family, parts,
                                                          approximation = approximation,
                                                          control = control,
                                                          engine = backend)
  fit$approximate <- TRUE
  fit$metrics$family_logLik_fixed_only <- as.numeric(parts$loglik(y, mu))
  fit$metrics$approximate_marginal_logLik <- as.numeric(fit$approximation_details$logLik)
  df_total <- fit$metrics$df %||% (length(beta) + length(fit$raw_theta %||% numeric(0L)))
  n <- length(y)
  fit$metrics$AIC_approximation <- if (is.finite(fit$metrics$approximate_marginal_logLik)) {
    -2 * fit$metrics$approximate_marginal_logLik + 2 * df_total
  } else NA_real_
  fit$metrics$BIC_approximation <- if (is.finite(fit$metrics$approximate_marginal_logLik)) {
    -2 * fit$metrics$approximate_marginal_logLik + log(n) * df_total
  } else NA_real_
  fit$metrics$approximation <- approximation
  fit$metrics$approximation_label <- fit$approximation_label
  fit$metrics$approximation_note <- fit$approximation_details$note
  fit$metrics$approximation_fallback_groups <- fit$approximation_details$fallback_groups
  fit$metrics$family_link <- parts$link
  fit$metrics$family_parameters <- parts$parameters
  fit$convergence <- if (isTRUE(converged) && fit$convergence == 0L) 0L else 1L
  fit
}

.memwas_fit_model <- function(formula, data, id, time, random, family, autocor,
                              L1_penalty, L2_penalty, control, method, random_cov,
                              approximation = "laplace",
                              init_approximation = "variational_inference",
                              se_method = "hessian",
                              engine = "R") {
  backend <- .memwas_normalize_engine(engine)

  family <- .memwas_normalize_family(family)
  approximation <- .memwas_normalize_approximation(approximation)
  init_approximation <- .memwas_normalize_init_approximation(init_approximation)
  se_method <- .memwas_normalize_se_method(se_method)
  prep <- .memwas_prepare_model_data(formula, data, id, time, random, family)
  if (family == "gaussian") {
    fit <- .memwas_fit_gaussian_lmm(prep, autocor = autocor, random_cov = random_cov,
                                    L1_penalty = L1_penalty, L2_penalty = L2_penalty,
                                    control = control, method = method, engine = backend)
    fit$approximation <- "exact_gaussian"
    fit$approximation_label <- "Exact Gaussian marginal likelihood"
    fit$metrics$approximation <- "exact_gaussian"
    fit$metrics$approximation_label <- fit$approximation_label
    return(fit)
  }
  .memwas_fit_nongaussian_marginal(prep, family = family, autocor = autocor, random_cov = random_cov,
                                      L1_penalty = L1_penalty, L2_penalty = L2_penalty,
                                      control = control, method = method,
                                      approximation = approximation,
                                      init_approximation = init_approximation,
                                      se_method = se_method,
                                      id = id, time = time, random = random,
                                      engine = backend)
}

.memwas_add_stored_splines <- function(data, formula, spline_variables, spline_info,
                                       engine = "R") {
  backend <- .memwas_normalize_engine(engine)

  added <- character(0L)
  if (length(spline_variables) == 0L) return(list(data = data, formula = formula, added_terms = added))
  for (v in spline_variables) {
    if (!v %in% names(data)) stop("Spline-selected predictor `", v, "` is not present in `data`.", call. = FALSE)
    info <- spline_info[[v]]
    if (is.null(info$knots) || is.null(info$basis_names)) {
      stop("Missing stored spline information for predictor `", v, "`.", call. = FALSE)
    }
    B <- .memwas_rcs_basis(data[[v]], info$knots, engine = backend)
    if (ncol(B) != length(info$basis_names)) {
      stop("Stored spline basis for `", v, "` has inconsistent dimension.", call. = FALSE)
    }
    for (j in seq_len(ncol(B))) data[[info$basis_names[j]]] <- B[, j]
    added <- c(added, info$basis_names)
  }
  list(data = data, formula = .memwas_add_terms_to_formula(formula, added), added_terms = added)
}
