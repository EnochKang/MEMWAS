# Internal non-Gaussian approximation helpers for MEMWAS.
# Subject-level integration, conditional-mode solving, and approximation kernels
# are evaluated in native C++ through .Call(). R orchestrates formula handling,
# optimizer calls, S3 object assembly, and optional resampling/profiling.

.memwas_default_approximation <- function() "laplace"
.memwas_default_init_approximation <- function() "variational_inference"

.memwas_normalize_init_approximation <- function(init_approximation = "variational_inference") {
  if (is.null(init_approximation) || length(init_approximation) == 0L ||
      is.na(init_approximation[1L]) || !nzchar(as.character(init_approximation[1L]))) {
    init_approximation <- "variational_inference"
  }
  .memwas_normalize_approximation(init_approximation)
}

.memwas_supported_se_methods <- function() {
  c("hessian", "cluster_sandwich", "parametric_bootstrap", "profile")
}

.memwas_normalize_se_method <- function(se_method = "hessian") {
  if (is.null(se_method) || length(se_method) == 0L ||
      all(is.na(se_method)) || !any(nzchar(as.character(se_method)))) {
    se_method <- "hessian"
  }
  key <- tolower(trimws(as.character(se_method)))
  key <- gsub("[[:space:]-]+", "_", key)
  key <- gsub("\\.+", "_", key)
  aliases <- c(
    observed_hessian = "hessian", hessian = "hessian", ml_hessian = "hessian",
    cluster = "cluster_sandwich", sandwich = "cluster_sandwich",
    cluster_sandwich = "cluster_sandwich", robust = "cluster_sandwich",
    cluster_robust = "cluster_sandwich",
    bootstrap = "parametric_bootstrap", parametric_bootstrap = "parametric_bootstrap",
    param_bootstrap = "parametric_bootstrap",
    profile = "profile", profile_likelihood = "profile", profile_ci = "profile"
  )
  if (any(!key %in% names(aliases))) {
    bad <- unique(key[!key %in% names(aliases)])
    stop("Unsupported `se_method`. Supported values are 'hessian', 'cluster_sandwich', ",
         "'parametric_bootstrap', and 'profile'. Invalid value(s): ",
         paste(bad, collapse = ", "), call. = FALSE)
  }
  unique(unname(aliases[key]))
}

.memwas_select_default_approximation <- function(family, q, approximation = NULL, control = list()) {
  family <- .memwas_normalize_family(family)
  if (family == "gaussian") return("exact_gaussian")
  if (!is.null(approximation) && length(approximation) > 0L &&
      !is.na(approximation[1L]) && nzchar(as.character(approximation[1L]))) {
    key <- tolower(trimws(as.character(approximation[1L])))
    if (!key %in% c("auto", "automatic", "default")) return(.memwas_normalize_approximation(approximation))
  }
  q <- as.integer(q %||% 0L)
  aghq_max_dim <- as.integer(control$aghq_max_dim %||% control$agq_max_dim %||% 5L)
  if (!is.finite(aghq_max_dim) || aghq_max_dim < 1L) aghq_max_dim <- 5L
  prefer_aghq <- isTRUE(control$auto_aghq %||% FALSE)
  if (isTRUE(prefer_aghq) && q > 0L && q <= aghq_max_dim) return("adaptive_gauss_hermite_quadrature")
  "laplace"
}

.memwas_family_npar <- function(family, control = list()) {
  family <- .memwas_normalize_family(family)
  if (family == "negative_binomial") {
    if (isFALSE(control$estimate_dispersion %||% TRUE) || isFALSE(control$estimate_nb_theta %||% TRUE)) return(0L)
    return(1L)
  }
  if (family == "gamma") {
    if (isFALSE(control$estimate_dispersion %||% TRUE) || isFALSE(control$estimate_gamma_shape %||% TRUE)) return(0L)
    return(1L)
  }
  0L
}

.memwas_random_par_from_cov <- function(D, random_cov = "diagonal") {
  D <- as.matrix(D)
  q <- nrow(D)
  if (q == 0L) return(numeric(0L))
  random_cov <- match.arg(random_cov, c("diagonal", "unstructured"))
  if (random_cov == "diagonal") return(log(pmax(diag(D), 1e-8)))
  U <- try(chol((D + t(D)) / 2 + diag(1e-8, q)), silent = TRUE)
  if (inherits(U, "try-error")) {
    U <- diag(sqrt(pmax(diag(D), 1e-8)), q)
  }
  L <- t(U)
  idx <- which(lower.tri(L, diag = TRUE), arr.ind = TRUE)
  out <- numeric(nrow(idx))
  for (k in seq_len(nrow(idx))) {
    i <- idx[k, 1L]
    j <- idx[k, 2L]
    out[k] <- if (i == j) log(pmax(L[i, j], 1e-8)) else L[i, j]
  }
  out
}

.memwas_random_cov_from_par <- function(par, q, random_cov = "diagonal", engine = "R") {
  backend <- .memwas_normalize_engine(engine)

  .memwas_build_random_cov(par, q, match.arg(random_cov, c("diagonal", "unstructured")), engine = backend)
}

.memwas_family_par_start <- function(family, y, control = list()) {
  family <- .memwas_normalize_family(family)
  npar <- .memwas_family_npar(family, control)
  if (npar == 0L) return(numeric(0L))
  if (family == "negative_binomial") return(log(.memwas_estimate_nb_theta(y, control)))
  if (family == "gamma") return(log(.memwas_estimate_gamma_shape(y, control)))
  numeric(0L)
}

.memwas_family_params_from_par <- function(family, par_tail, y, control = list()) {
  family <- .memwas_normalize_family(family)
  npar <- .memwas_family_npar(family, control)
  if (family == "negative_binomial") {
    theta <- if (npar >= 1L && length(par_tail) >= 1L) exp(par_tail[1L]) else .memwas_estimate_nb_theta(y, control)
    theta <- pmax(theta, 1e-6)
    return(list(theta = theta, size = theta, dispersion_estimated = npar >= 1L))
  }
  if (family == "gamma") {
    shape <- if (npar >= 1L && length(par_tail) >= 1L) exp(par_tail[1L]) else .memwas_estimate_gamma_shape(y, control)
    shape <- pmax(shape, 1e-6)
    return(list(shape = shape, dispersion_estimated = npar >= 1L))
  }
  if (family == "exponential") return(list(shape = 1, dispersion_estimated = FALSE))
  list()
}

.memwas_gh_nodes_weights <- function(n = 7L) {
  gh <- .memwas_gauss_hermite_nodes(n)
  list(nodes = as.numeric(gh$nodes), weights = as.numeric(gh$weights))
}

.memwas_cpp_glmm_approximation <- function(par, prep, family, approximation,
                                           random_cov = "diagonal", control = list(),
                                           return_details = FALSE) {
  family <- .memwas_normalize_family(family)
  approximation <- .memwas_normalize_approximation(approximation)
  random_cov <- match.arg(random_cov, c("diagonal", "unstructured"))
  p <- ncol(prep$X)
  q <- ncol(prep$Z)
  rn <- .memwas_random_npar(q, random_cov)
  fn <- .memwas_family_npar(family, control)
  base_parts <- .memwas_family_parts(family, y = prep$y, control = control)
  fixed_theta <- base_parts$parameters$theta %||% base_parts$parameters$size %||% .memwas_estimate_nb_theta(prep$y, control)
  fixed_shape <- base_parts$parameters$shape %||% .memwas_estimate_gamma_shape(prep$y, control)
  n_nodes <- as.integer(control$aghq_nodes %||% control$agq_nodes %||% control$quadrature_nodes %||% 7L)
  if (!is.finite(n_nodes) || n_nodes < 3L) n_nodes <- 7L
  gh <- .memwas_gh_nodes_weights(n_nodes)
  max_dim <- as.integer(control$aghq_max_dim %||% control$agq_max_dim %||% 5L)
  if (!is.finite(max_dim) || max_dim < 1L) max_dim <- 5L
  max_nodes <- as.integer(control$aghq_max_nodes %||% control$agq_max_nodes %||% 50000L)
  if (!is.finite(max_nodes) || max_nodes < 100L) max_nodes <- 50000L
  mode_maxit <- as.integer(control$mode_maxit %||% control$conditional_mode_maxit %||% 50L)
  if (!is.finite(mode_maxit) || mode_maxit < 1L) mode_maxit <- 50L
  mode_tol <- as.numeric(control$mode_tol %||% control$conditional_mode_tol %||% 1e-7)
  if (!is.finite(mode_tol) || mode_tol <= 0) mode_tol <- 1e-7

  storage.mode(prep$X) <- "double"
  storage.mode(prep$Z) <- "double"
  .Call("memwas_glmm_approximation",
        as.numeric(par), prep$X, prep$Z, as.numeric(prep$y), prep$groups,
        as.character(family), as.character(approximation), as.character(random_cov),
        as.integer(p), as.integer(q), as.integer(rn), as.integer(fn),
        as.numeric(fixed_theta), as.numeric(fixed_shape),
        as.numeric(gh$nodes), as.numeric(gh$weights),
        as.integer(max_dim), as.integer(max_nodes), as.integer(mode_maxit),
        as.numeric(mode_tol), as.logical(return_details), PACKAGE = "MEMWAS")
}

.memwas_numeric_hessian <- function(par, fn, eps = NULL) {
  par <- as.numeric(par)
  k <- length(par)
  if (k == 0L) return(matrix(0, 0L, 0L))
  if (is.null(eps)) eps <- pmax(1e-4, abs(par) * 1e-4)
  eps <- rep(as.numeric(eps), length.out = k)
  H <- matrix(NA_real_, k, k)
  f0 <- fn(par)
  for (i in seq_len(k)) {
    ei <- rep(0, k); ei[i] <- eps[i]
    fpi <- fn(par + ei)
    fmi <- fn(par - ei)
    H[i, i] <- (fpi - 2 * f0 + fmi) / (eps[i]^2)
    if (i < k) {
      for (j in (i + 1L):k) {
        ej <- rep(0, k); ej[j] <- eps[j]
        fpp <- fn(par + ei + ej)
        fpm <- fn(par + ei - ej)
        fmp <- fn(par - ei + ej)
        fmm <- fn(par - ei - ej)
        H[i, j] <- H[j, i] <- (fpp - fpm - fmp + fmm) / (4 * eps[i] * eps[j])
      }
    }
  }
  H[!is.finite(H)] <- NA_real_
  H
}

.memwas_hessian_vcov <- function(par, objective, p, control = list(), engine = "R") {
  backend <- .memwas_normalize_engine(engine)

  h_eps <- control$hessian_eps %||% NULL
  H <- try(.memwas_numeric_hessian(par, objective, eps = h_eps), silent = TRUE)
  if (inherits(H, "try-error") || any(!is.finite(H))) return(list(vcov = NULL, hessian = NULL, positive_definite = FALSE))
  H <- (H + t(H)) / 2
  eg <- try(eigen(H, symmetric = TRUE, only.values = TRUE), silent = TRUE)
  pd <- !inherits(eg, "try-error") && all(is.finite(eg$values)) && min(eg$values) > 1e-8
  Vinv <- try(.memwas_safe_solve(H + diag(1e-8, nrow(H)), engine = backend), silent = TRUE)
  if (inherits(Vinv, "try-error") || any(!is.finite(Vinv))) return(list(vcov = NULL, hessian = H, positive_definite = FALSE))
  Vb <- Vinv[seq_len(p), seq_len(p), drop = FALSE]
  list(vcov = (Vb + t(Vb)) / 2, hessian = H, positive_definite = pd)
}

.memwas_cluster_sandwich_vcov <- function(par, prep, family, approximation, random_cov,
                                          control, bread, engine = "R") {
  p <- ncol(prep$X)
  if (is.null(bread) || any(!is.finite(bread))) return(NULL)
  eps <- pmax(1e-4, abs(par[seq_len(p)]) * 1e-4)
  base <- .memwas_cpp_glmm_approximation(par, prep, family, approximation,
                                         random_cov = random_cov, control = control,
                                         return_details = TRUE)
  if (is.null(base$details) || !is.data.frame(base$details)) return(NULL)
  G <- nrow(base$details)
  scores <- matrix(0, G, p)
  for (j in seq_len(p)) {
    pp <- par; pm <- par
    pp[j] <- pp[j] + eps[j]
    pm[j] <- pm[j] - eps[j]
    lp <- .memwas_cpp_glmm_approximation(pp, prep, family, approximation, random_cov, control, TRUE)$details$logLik
    lm <- .memwas_cpp_glmm_approximation(pm, prep, family, approximation, random_cov, control, TRUE)$details$logLik
    scores[, j] <- (lp - lm) / (2 * eps[j])
  }
  meat <- crossprod(scores)
  V <- bread %*% meat %*% bread
  dimnames(V) <- list(colnames(prep$X), colnames(prep$X))
  (V + t(V)) / 2
}

.memwas_response_simulate <- function(family, mu, params) {
  family <- .memwas_normalize_family(family)
  mu <- as.numeric(mu)
  if (family == "binomial") return(stats::rbinom(length(mu), size = 1L, prob = pmin(pmax(mu, 1e-6), 1 - 1e-6)))
  if (family == "poisson") return(stats::rpois(length(mu), lambda = pmax(mu, 1e-8)))
  if (family == "negative_binomial") {
    theta <- params$theta %||% params$size %||% 1
    return(stats::rnbinom(length(mu), size = pmax(theta, 1e-6), mu = pmax(mu, 1e-8)))
  }
  if (family == "gamma") {
    shape <- params$shape %||% 1
    return(stats::rgamma(length(mu), shape = pmax(shape, 1e-6), rate = pmax(shape, 1e-6) / pmax(mu, 1e-8)))
  }
  if (family == "exponential") return(stats::rexp(length(mu), rate = 1 / pmax(mu, 1e-8)))
  stats::rnorm(length(mu), mean = mu, sd = 1)
}

.memwas_parametric_bootstrap_vcov <- function(fit, prep, formula, id, time, random, autocor,
                                              family, random_cov, L1_penalty, L2_penalty,
                                              control, method, approximation,
                                              init_approximation, engine = "R") {
  backend <- .memwas_normalize_engine(engine)

  B <- as.integer(control$bootstrap_B %||% control$parametric_bootstrap_B %||% 0L)
  if (!is.finite(B) || B <= 1L) return(NULL)
  if (is.null(id) || is.null(time) || is.null(random)) return(NULL)
  p <- length(fit$coefficients)
  boot <- matrix(NA_real_, B, p)
  colnames(boot) <- names(fit$coefficients)
  response_name <- all.vars(formula)[1L]
  for (b in seq_len(B)) {
    dat_b <- prep$data
    dat_b[[response_name]] <- .memwas_response_simulate(family, fit$fitted_response %||% fit$fitted,
                                                        fit$family_parameters %||% list())
    ctrl_b <- control
    ctrl_b$bootstrap_B <- 0L
    ctrl_b$parametric_bootstrap_B <- 0L
    fit_b <- try(.memwas_fit_model(formula = formula, data = dat_b, id = id, time = time,
                                   random = random, family = family, autocor = autocor,
                                   L1_penalty = L1_penalty, L2_penalty = L2_penalty,
                                   control = ctrl_b, method = method, random_cov = random_cov,
                                   approximation = approximation,
                                   init_approximation = init_approximation,
                                   se_method = "hessian", engine = backend), silent = TRUE)
    if (!inherits(fit_b, "try-error") && length(fit_b$coefficients) == p) {
      boot[b, ] <- as.numeric(fit_b$coefficients)
    }
  }
  ok <- stats::complete.cases(boot)
  if (sum(ok) < 2L) return(NULL)
  stats::cov(boot[ok, , drop = FALSE])
}

.memwas_profile_ci <- function(par, objective, p, rn, fn, random_cov, family,
                               control = list(), engine = "R") {
  do_profile <- isTRUE(control$profile_ci %||% FALSE) || isTRUE(control$profile_likelihood_ci %||% FALSE)
  if (!do_profile) return(data.frame())
  alpha <- control$profile_alpha %||% 0.05
  alpha <- suppressWarnings(as.numeric(alpha)[1L])
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) alpha <- 0.05
  cutoff <- stats::qchisq(1 - alpha, df = 1)
  n_grid <- as.integer(control$profile_grid %||% 21L)
  if (!is.finite(n_grid) || n_grid < 9L) n_grid <- 21L
  idx <- seq.int(from = p + 1L, length.out = rn + fn)
  if (!length(idx)) return(data.frame())
  nll0 <- objective(par)
  rows <- vector("list", length(idx))
  for (m in seq_along(idx)) {
    j <- idx[m]
    grid <- seq(par[j] - 2, par[j] + 2, length.out = n_grid)
    prof <- rep(NA_real_, n_grid)
    for (g in seq_along(grid)) {
      pg <- par; pg[j] <- grid[g]
      prof[g] <- objective(pg)
    }
    keep <- is.finite(prof) & (2 * (prof - nll0) <= cutoff)
    if (any(keep)) {
      low_raw <- min(grid[keep]); high_raw <- max(grid[keep])
    } else {
      low_raw <- high_raw <- NA_real_
    }
    nm <- if (m <= rn) paste0("random_covariance_parameter_", m) else {
      if (family == "negative_binomial") "negative_binomial_theta" else if (family == "gamma") "gamma_shape" else paste0("family_parameter_", m - rn)
    }
    rows[[m]] <- data.frame(parameter = nm,
                            estimate_raw = par[j],
                            conf_low_raw = low_raw,
                            conf_high_raw = high_raw,
                            estimate = exp(par[j]),
                            conf_low = if (is.finite(low_raw)) exp(low_raw) else NA_real_,
                            conf_high = if (is.finite(high_raw)) exp(high_raw) else NA_real_,
                            alpha = alpha,
                            method = "profile_likelihood_grid",
                            stringsAsFactors = FALSE,
                            check.names = FALSE)
  }
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}


.memwas_make_nongaussian_fit_shell <- function(prep, family, autocor, method = "ML") {
  ac <- .memwas_normalize_autocor(autocor)
  list(call = match.call(), formula = prep$formula, data = prep$data, family = family,
       family_link = NA_character_, family_parameters = list(), y = prep$y,
       X = prep$X, Z = prep$Z, id = prep$id, time = prep$time,
       groups = prep$groups, coefficients = numeric(0L), coefficient_table = data.frame(),
       coefficient_term_map = prep$coefficient_term_map, vcov = matrix(NA_real_, 0L, 0L),
       fitted = numeric(0L), fitted_link = numeric(0L), fitted_response = numeric(0L),
       residuals = numeric(0L), random_effects = setNames(vector("list", length(prep$groups)), names(prep$groups)),
       random_effects_se = setNames(vector("list", length(prep$groups)), names(prep$groups)),
       random_effects_ci = setNames(vector("list", length(prep$groups)), names(prep$groups)),
       random_effects_covariance = setNames(vector("list", length(prep$groups)), names(prep$groups)),
       random_effects_table = data.frame(id = character(0L), effect = character(0L), estimate = numeric(0L),
                                         std_error = numeric(0L), conf_low = numeric(0L), conf_high = numeric(0L),
                                         ci_level = numeric(0L), stringsAsFactors = FALSE, check.names = FALSE),
       random_effects_ci_level = NA_real_, random_covariance = matrix(0, ncol(prep$Z), ncol(prep$Z)),
       residual_sigma = NA_real_,
       autocorrelation = list(type = ac$label,
                              note = "For non-Gaussian marginal fits, MEMWAS carries the requested autocorrelation setting in the model object; the native likelihood approximation conditions on the subject-specific random-effect integration."),
       raw_theta = numeric(0L), metrics = list(method = toupper(method), nobs = length(prep$y)),
       convergence = NA_integer_, optim_message = "", approximate = TRUE,
       dropped_rows = prep$dropped_rows)
}

.memwas_cpp_initialize_nongaussian <- function(prep, family, random_cov = "diagonal", control = list(),
                                               init_approximation = "variational_inference",
                                               engine = "cpp") {
  family <- .memwas_normalize_family(family)
  init_approximation <- .memwas_normalize_init_approximation(init_approximation)
  random_cov <- match.arg(random_cov, c("diagonal", "unstructured"))
  p <- ncol(prep$X)
  q <- ncol(prep$Z)
  rn <- .memwas_random_npar(q, random_cov)
  fn <- .memwas_family_npar(family, control)
  beta0 <- .memwas_initial_glm_beta(prep$X, prep$y, family, control)
  if (length(beta0) != p || any(!is.finite(beta0))) beta0 <- rep(0, p)
  names(beta0) <- colnames(prep$X)
  rpar0 <- if (rn > 0L) rep(log(0.1), rn) else numeric(0L)
  fpar0 <- .memwas_family_par_start(family, prep$y, control)
  if (length(fpar0) != fn || any(!is.finite(fpar0))) fpar0 <- rep(0, fn)
  par0 <- c(beta0, rpar0, fpar0)
  if (!length(par0)) {
    return(list(par = par0, convergence = 0L, message = "No parameters required initialization.",
                method = init_approximation, optimizer_method = NA_character_, objective = NA_real_,
                backend = "native_cpp_subject_level_kernels", used = FALSE))
  }

  init_objective <- function(par) {
    out <- try(.memwas_cpp_glmm_approximation(par, prep, family, init_approximation,
                                              random_cov = random_cov, control = control,
                                              return_details = FALSE), silent = TRUE)
    if (inherits(out, "try-error") || is.null(out$negLogLik) || !is.finite(out$negLogLik)) return(1e100)
    as.numeric(out$negLogLik)
  }
  opt_control <- control$init_optim_control %||% control$vi_optim_control %||%
    control$initialization_optim_control %||% list(maxit = as.integer(control$vi_maxit %||% 100L))
  opt_method <- control$init_optim_method %||% control$vi_optim_method %||%
    control$initialization_optim_method %||% "BFGS"
  opt <- try(stats::optim(par0, init_objective, method = opt_method, control = opt_control), silent = TRUE)
  used_method <- opt_method
  if (inherits(opt, "try-error") || !is.finite(opt$value)) {
    opt <- try(stats::optim(par0, init_objective, method = "Nelder-Mead", control = opt_control), silent = TRUE)
    used_method <- "Nelder-Mead"
  }
  if (inherits(opt, "try-error") || !is.finite(opt$value)) {
    return(list(par = par0, convergence = 1L,
                message = paste("Native-C++", .memwas_approximation_label(init_approximation),
                                "initialization optimizer failed; GLM/dispersion starts were used."),
                method = init_approximation, optimizer_method = used_method, objective = init_objective(par0),
                backend = "native_cpp_subject_level_kernels", used = FALSE))
  }
  list(par = as.numeric(opt$par), convergence = opt$convergence,
       message = opt$message %||% "", method = init_approximation,
       optimizer_method = used_method, objective = as.numeric(opt$value),
       backend = "native_cpp_subject_level_kernels", used = TRUE)
}

.memwas_build_random_effect_outputs_from_cpp <- function(cpp, groups, z_names, ci_level = 0.95) {
  q <- length(z_names)
  gnames <- names(groups)
  if (q == 0L || is.null(cpp$random_effects)) {
    empty <- data.frame(id = character(0L), effect = character(0L), estimate = numeric(0L),
                        std_error = numeric(0L), conf_low = numeric(0L), conf_high = numeric(0L),
                        ci_level = numeric(0L), stringsAsFactors = FALSE, check.names = FALSE)
    return(list(random_effects = setNames(vector("list", length(groups)), gnames),
                random_effects_se = setNames(vector("list", length(groups)), gnames),
                random_effects_ci = setNames(vector("list", length(groups)), gnames),
                random_effects_covariance = setNames(vector("list", length(groups)), gnames),
                random_effects_table = empty))
  }
  B <- as.matrix(cpp$random_effects)
  S <- as.matrix(cpp$random_effects_se)
  cov_list <- cpp$random_effects_covariance %||% vector("list", nrow(B))
  ci_mult <- stats::qnorm(0.5 + ci_level / 2)
  re <- se <- ci <- cov <- vector("list", nrow(B))
  names(re) <- names(se) <- names(ci) <- names(cov) <- gnames
  rows <- vector("list", nrow(B))
  for (g in seq_len(nrow(B))) {
    bi <- as.numeric(B[g, seq_len(q), drop = TRUE]); names(bi) <- z_names
    si <- as.numeric(S[g, seq_len(q), drop = TRUE]); names(si) <- z_names
    cii <- cbind(conf_low = bi - ci_mult * si, conf_high = bi + ci_mult * si)
    rownames(cii) <- z_names
    cvi <- if (length(cov_list) >= g && is.matrix(cov_list[[g]])) cov_list[[g]] else diag(si^2, q)
    dimnames(cvi) <- list(z_names, z_names)
    re[[g]] <- bi; se[[g]] <- si; ci[[g]] <- cii; cov[[g]] <- cvi
    rows[[g]] <- data.frame(id = rep(gnames[g], q), effect = z_names,
                            estimate = bi, std_error = si,
                            conf_low = cii[, "conf_low"], conf_high = cii[, "conf_high"],
                            ci_level = rep(ci_level, q), stringsAsFactors = FALSE,
                            check.names = FALSE)
  }
  table <- do.call(rbind, rows); row.names(table) <- NULL
  list(random_effects = re, random_effects_se = se, random_effects_ci = ci,
       random_effects_covariance = cov, random_effects_table = table)
}


.memwas_fit_nongaussian_marginal <- function(prep, family, autocor, random_cov, L1_penalty, L2_penalty,
                                             control, method = "ML", approximation = "laplace",
                                             init_approximation = "variational_inference",
                                             se_method = "hessian", id = NULL, time = NULL,
                                             random = NULL, engine = "cpp") {
  backend <- .memwas_normalize_engine(engine)

  family <- .memwas_normalize_family(family)
  approximation <- .memwas_normalize_approximation(approximation %||% "laplace")
  init_approximation <- .memwas_normalize_init_approximation(init_approximation %||% "variational_inference")
  se_method <- .memwas_normalize_se_method(se_method %||% "hessian")
  random_cov <- match.arg(random_cov, c("diagonal", "unstructured"))
  p <- ncol(prep$X)
  q <- ncol(prep$Z)
  rn <- .memwas_random_npar(q, random_cov)
  fn <- .memwas_family_npar(family, control)

  init_start <- .memwas_cpp_initialize_nongaussian(prep, family = family,
                                                   random_cov = random_cov,
                                                   control = control,
                                                   init_approximation = init_approximation,
                                                   engine = backend)
  par0 <- init_start$par

  objective <- function(par) {
    out <- try(.memwas_cpp_glmm_approximation(par, prep, family, approximation,
                                              random_cov = random_cov, control = control,
                                              return_details = FALSE), silent = TRUE)
    if (inherits(out, "try-error") || is.null(out$negLogLik) || !is.finite(out$negLogLik)) return(1e100)
    as.numeric(out$negLogLik)
  }
  opt_control <- control$nongaussian_optim_control %||% control$approximation_optim_control %||% .memwas_extract_optim_control(control)
  opt_method <- control$nongaussian_optim_method %||% control$approximation_optim_method %||% "BFGS"
  opt <- try(stats::optim(par0, objective, method = opt_method, control = opt_control), silent = TRUE)
  used_opt_method <- opt_method
  if (inherits(opt, "try-error") || !is.finite(opt$value)) {
    opt <- try(stats::optim(par0, objective, method = "Nelder-Mead", control = opt_control), silent = TRUE)
    used_opt_method <- "Nelder-Mead"
  }
  if (inherits(opt, "try-error") || !is.finite(opt$value)) {
    warning("Native non-Gaussian final approximation optimizer failed; returning a legacy PQL fallback fit. Error: ",
            if (inherits(opt, "try-error")) as.character(opt) else "non-finite objective value", call. = FALSE)
    fallback_fit <- .memwas_fit_pql_glmm(prep, family = family, autocor = autocor, random_cov = random_cov,
                                         L1_penalty = L1_penalty, L2_penalty = L2_penalty,
                                         control = control, method = method,
                                         approximation = "pql", engine = backend)
    fallback_fit$approximation <- "pql"
    fallback_fit$approximation_label <- .memwas_approximation_label("pql")
    fallback_fit$init_approximation <- init_approximation
    fallback_fit$init_approximation_label <- .memwas_approximation_label(init_approximation)
    fallback_fit$initialization <- init_start
    fallback_fit$final_approximation_failed <- TRUE
    fallback_fit$se_method <- se_method
    fallback_fit$approximation_strategy <- list(initial = init_approximation, final = "pql",
                                                backend = "legacy_R_PQL_fallback_after_native_failure",
                                                random_effect_dimension = q,
                                                random_cov = random_cov)
    return(fallback_fit)
  }

  final_cpp <- .memwas_cpp_glmm_approximation(opt$par, prep, family, approximation,
                                              random_cov = random_cov, control = control,
                                              return_details = TRUE)
  beta <- as.numeric(opt$par[seq_len(p)]); names(beta) <- colnames(prep$X)
  rand_par <- if (rn > 0L) opt$par[p + seq_len(rn)] else numeric(0L)
  fam_tail <- if (fn > 0L) opt$par[p + rn + seq_len(fn)] else numeric(0L)
  D <- .memwas_random_cov_from_par(rand_par, q, random_cov, engine = backend)
  dimnames(D) <- list(colnames(prep$Z), colnames(prep$Z))
  fam_params <- .memwas_family_params_from_par(family, fam_tail, prep$y, control)
  control_final <- modifyList(control, fam_params)
  parts <- .memwas_family_parts(family, y = prep$y, control = control_final)
  eta <- as.vector(prep$X %*% beta)
  mu <- parts$linkinv(eta)

  hess <- .memwas_hessian_vcov(opt$par, objective, p = p, control = control, engine = backend)
  vcov_beta <- hess$vcov
  se_source <- "observed_hessian"
  if ("cluster_sandwich" %in% se_method && !is.null(hess$vcov)) {
    Vrob <- .memwas_cluster_sandwich_vcov(opt$par, prep, family, approximation,
                                          random_cov, control, bread = hess$vcov,
                                          engine = backend)
    if (!is.null(Vrob)) {
      vcov_beta <- Vrob
      se_source <- "cluster_sandwich_by_id"
    }
  }
  if ("parametric_bootstrap" %in% se_method) {
    sim_fit <- list(coefficients = beta, fitted_response = mu, family_parameters = fam_params)
    Vboot <- .memwas_parametric_bootstrap_vcov(sim_fit, prep, prep$formula, id,
                                               time, random, autocor, family, random_cov,
                                               L1_penalty, L2_penalty, control, method,
                                               approximation, init_approximation, engine = backend)
    if (!is.null(Vboot)) {
      vcov_beta <- Vboot
      se_source <- "parametric_bootstrap"
    }
  }
  if (is.null(vcov_beta)) {
    vcov_beta <- matrix(NA_real_, p, p)
    se_source <- "unavailable"
  }
  dimnames(vcov_beta) <- list(colnames(prep$X), colnames(prep$X))
  se <- sqrt(pmax(diag(vcov_beta), 0))
  zval <- beta / pmax(se, 1e-12)
  pval <- 2 * stats::pnorm(abs(zval), lower.tail = FALSE)
  coef_table <- data.frame(term = names(beta), estimate = as.numeric(beta),
                           std_error = as.numeric(se), statistic = as.numeric(zval),
                           p_value = as.numeric(pval), row.names = NULL,
                           check.names = FALSE)

  ci_level <- control$random_effect_ci_level %||% control$blup_ci_level %||% 0.95
  ci_level <- suppressWarnings(as.numeric(ci_level)[1L])
  if (!is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) ci_level <- 0.95
  re_out <- .memwas_build_random_effect_outputs_from_cpp(final_cpp, prep$groups, colnames(prep$Z), ci_level)
  profile_control <- control
  if ("profile" %in% se_method) profile_control$profile_ci <- TRUE
  profile_ci <- .memwas_profile_ci(opt$par, objective, p = p, rn = rn, fn = fn,
                                   random_cov = random_cov, family = family,
                                   control = profile_control, engine = backend)

  n <- length(prep$y)
  df_total <- p + rn + fn
  logLik <- as.numeric(final_cpp$logLik)
  metrics <- list()
  metrics$logLik <- logLik
  metrics$penalized_logLik <- logLik
  metrics$AIC <- -2 * logLik + 2 * df_total
  metrics$BIC <- -2 * logLik + log(n) * df_total
  metrics$deviance <- -2 * logLik
  metrics$df <- df_total
  metrics$nobs <- n
  metrics$method <- "ML"
  metrics$approximate_marginal_logLik <- logLik
  metrics$AIC_approximation <- metrics$AIC
  metrics$BIC_approximation <- metrics$BIC
  metrics$approximation <- approximation
  metrics$approximation_label <- .memwas_approximation_label(approximation)
  metrics$init_approximation <- init_approximation
  metrics$init_approximation_label <- .memwas_approximation_label(init_approximation)
  metrics$family_link <- parts$link
  metrics$family_parameters <- fam_params
  metrics$se_source <- se_source
  metrics$hessian_positive_definite <- isTRUE(hess$positive_definite)
  metrics$approximation_fallback_groups <- final_cpp$fallback_groups %||% 0L

  fit <- .memwas_make_nongaussian_fit_shell(prep, family = family, autocor = autocor, method = method)
  fit$family <- family
  fit$family_link <- parts$link
  fit$family_parameters <- fam_params
  fit$coefficients <- beta
  fit$coefficient_table <- coef_table
  fit$vcov <- vcov_beta
  fit$fitted_link <- eta
  fit$fitted_response <- mu
  fit$fitted <- mu
  fit$residuals <- prep$y - mu
  fit$random_covariance <- D
  fit$random_effects <- re_out$random_effects
  fit$random_effects_se <- re_out$random_effects_se
  fit$random_effects_ci <- re_out$random_effects_ci
  fit$random_effects_covariance <- re_out$random_effects_covariance
  fit$random_effects_table <- re_out$random_effects_table
  fit$random_effects_ci_level <- ci_level
  fit$random_effects_note <- "Non-Gaussian BLUP uncertainty uses the inverse native-C++ negative Hessian of the conditional log posterior at each subject-specific mode."
  fit$raw_theta <- opt$par
  fit$raw_theta_layout <- list(p = p, random_npar = rn, family_npar = fn,
                               beta = seq_len(p),
                               random = if (rn) p + seq_len(rn) else integer(0L),
                               family = if (fn) p + rn + seq_len(fn) else integer(0L))
  fit$metrics <- metrics
  fit$convergence <- opt$convergence
  fit$optim_message <- opt$message %||% ""
  fit$approximation <- approximation
  fit$approximation_label <- .memwas_approximation_label(approximation)
  fit$init_approximation <- init_approximation
  fit$init_approximation_label <- .memwas_approximation_label(init_approximation)
  fit$initialization <- init_start
  fit$approximation_strategy <- list(initial = init_approximation, final = approximation,
                                     backend = "native_cpp_subject_level_kernels",
                                     initialization_optimizer = init_start$optimizer_method,
                                     final_optimizer = used_opt_method,
                                     random_effect_dimension = q,
                                     random_cov = random_cov)
  fit$approximation_details <- list(logLik = logLik,
                                    approximation = approximation,
                                    label = .memwas_approximation_label(approximation),
                                    init_approximation = init_approximation,
                                    n_groups = length(prep$groups),
                                    n_quadrature_nodes = final_cpp$n_quadrature_nodes %||% 0L,
                                    fallback_groups = final_cpp$fallback_groups %||% 0L,
                                    note = final_cpp$note %||% "Native C++ non-Gaussian marginal approximation.",
                                    details = final_cpp$details)
  fit$approximation_diagnostics <- list(
    family = family,
    link = parts$link,
    approximation = approximation,
    approximation_label = .memwas_approximation_label(approximation),
    init_approximation = init_approximation,
    init_approximation_label = .memwas_approximation_label(init_approximation),
    initialization_convergence = init_start$convergence,
    initialization_message = init_start$message %||% "",
    random_effect_dimension = q,
    n_groups = length(prep$groups),
    nobs = n,
    converged = identical(opt$convergence, 0L),
    optimizer_convergence = opt$convergence,
    optimizer_message = opt$message %||% "",
    hessian_positive_definite = isTRUE(hess$positive_definite),
    se_source = se_source,
    fallback_groups = final_cpp$fallback_groups %||% 0L,
    quadrature_nodes = final_cpp$n_quadrature_nodes %||% 0L,
    profile_ci_computed = nrow(profile_ci) > 0L,
    details = final_cpp$details
  )
  fit$se_method <- se_method
  fit$se_source <- se_source
  fit$observed_hessian <- hess$hessian
  fit$profile_ci <- profile_ci
  fit$variance_parameter_ci <- profile_ci
  fit$approximate <- TRUE
  fit$pql <- list(iterations = NA_integer_, converged = NA,
                  note = paste("Legacy PQL was not used for the default path. Initialization used",
                               .memwas_approximation_label(init_approximation),
                               "through native-C++ kernels before the native-C++",
                               .memwas_approximation_label(approximation), "final fit."))
  fit
}

