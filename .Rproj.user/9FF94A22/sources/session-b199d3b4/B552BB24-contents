# Multi-family likelihood for MEMWAS
#
# It not only support Gaussian identity-link implementation as the default
# Gaussian distribution, but also support non-Gaussian distribution with a coherent
# Laplace/low-dimensional quadrature likelihood. The same predictor, offsets,
# grouped-binomial totals, penalties, random effects and serial effects are used
# by fitting, prediction and residual calculations.
#
# These functions are intentionally not exported and do not add package dependencies.

.memwas_legacy_fit_MEMWAS <- if (exists("fit_MEMWAS", inherits = FALSE)) fit_MEMWAS else NULL
.memwas_legacy_fit_formals <- if (!is.null(.memwas_legacy_fit_MEMWAS)) names(formals(.memwas_legacy_fit_MEMWAS)) else character()

.memwas_normalize_fit_call <- function(raw, envir) {
  if (!length(raw)) return(raw)
  nm <- names(raw); if (is.null(nm)) nm <- rep("", length(raw))
  unnamed <- which(nm == "")
  # Prefer semantic positional recognition, allowing both formula-data and data-formula APIs.
  vals <- lapply(raw[unnamed[seq_len(min(length(unnamed), 3L))]], function(z) try(eval(z, envir), silent = TRUE))
  isform <- vapply(vals, inherits, logical(1), "formula")
  isdata <- vapply(vals, function(z) is.data.frame(z) || is.matrix(z), logical(1))
  if (length(unnamed) >= 2L && isform[1L] && isdata[2L]) {
    nm[unnamed[1L]] <- "formula"; nm[unnamed[2L]] <- "data"
  } else if (length(unnamed) >= 2L && isdata[1L] && isform[2L]) {
    nm[unnamed[1L]] <- "data"; nm[unnamed[2L]] <- "formula"
  } else if (length(unnamed) >= 3L && isform[1L] && isform[2L] && isdata[3L]) {
    nm[unnamed[1L]] <- "formula"; nm[unnamed[2L]] <- "random"; nm[unnamed[3L]] <- "data"
  } else if (length(unnamed) >= 3L && isdata[1L] && isform[2L] && isform[3L]) {
    nm[unnamed[1L]] <- "data"; nm[unnamed[2L]] <- "formula"; nm[unnamed[3L]] <- "random"
  }
  # Map remaining unnamed arguments according to the legacy public formals.
  remaining <- which(nm == "")
  formal_names <- setdiff(.memwas_legacy_fit_formals, "...")
  used <- unique(nm[nzchar(nm)])
  available <- formal_names[!formal_names %in% used]
  if (length(remaining) && length(available)) {
    take <- seq_len(min(length(remaining), length(available)))
    nm[remaining[take]] <- available[take]
  }
  names(raw) <- nm
  raw
}

.memwas_null_coalesce <- function(x, y) if (is.null(x)) y else x
.memwas_is_scalar_character <- function(x) is.character(x) && length(x) == 1L && !is.na(x)
.memwas_safe_name <- function(x) gsub("[^A-Za-z0-9_.]+", "_", x)

.memwas_normalize_random_covariance <- function(x) {
  if (!.memwas_is_scalar_character(x) || !nzchar(trimws(x))) {
    stop("random_cov must be one non-empty character value.", call. = FALSE)
  }
  key <- tolower(gsub("[ _-]+", "", x))
  if (key %in% c("diag", "diagonal", "independent")) return("diagonal")
  if (key %in% c("un", "unstructured", "full")) return("unstructured")
  stop("random_cov must be 'diagonal' or 'unstructured'.", call. = FALSE)
}

.memwas_normalize_likelihood_method <- function(x) {
  if (!.memwas_is_scalar_character(x) || !nzchar(trimws(x))) {
    stop("method must be 'ML' or 'REML'.", call. = FALSE)
  }
  method <- toupper(trimws(x))
  if (!method %in% c("ML", "REML")) {
    stop("method must be 'ML' or 'REML'.", call. = FALSE)
  }
  method
}

.memwas_validate_quadrature_points <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < 3 || x != as.integer(x)) {
    stop("quadrature_points must be one integer of at least 3.", call. = FALSE)
  }
  as.integer(x)
}

.memwas_validate_raw_start_vector <- function(x, expected_length, label) {
  if (!is.numeric(x) || length(x) != expected_length || any(!is.finite(x))) {
    stop(label, " must be a finite numeric vector of length ", expected_length, ".", call. = FALSE)
  }
  as.numeric(x)
}

.memwas_family_logsumexp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}
.memwas_inverse_logit <- function(x) {
  out <- numeric(length(x))
  pos <- x >= 0
  out[pos] <- 1 / (1 + exp(-x[pos]))
  ex <- exp(x[!pos])
  out[!pos] <- ex / (1 + ex)
  out
}
.memwas_chol_pd <- function(x, jitter = 1e-9, max_tries = 9L) {
  x <- (x + t(x)) / 2
  if (!length(x)) return(matrix(numeric(), 0L, 0L))
  scale <- max(1, max(abs(diag(x)), na.rm = TRUE))
  for (k in 0:max_tries) {
    ans <- try(chol(x + diag(jitter * (10^k) * scale, nrow(x))), silent = TRUE)
    if (!inherits(ans, "try-error")) return(ans)
  }
  stop("A covariance or Hessian matrix is not positive definite.", call. = FALSE)
}
.memwas_solve_pd <- function(x, b = NULL) {
  r <- .memwas_chol_pd(x)
  if (is.null(b)) return(chol2inv(r))
  backsolve(r, forwardsolve(t(r), b))
}
.memwas_logdet_pd <- function(x) 2 * sum(log(diag(.memwas_chol_pd(x))))
.memwas_block_diag <- function(a, b) {
  if (!length(a)) return(b)
  if (!length(b)) return(a)
  out <- matrix(0, nrow(a) + nrow(b), ncol(a) + ncol(b))
  out[seq_len(nrow(a)), seq_len(ncol(a))] <- a
  ii <- nrow(a) + seq_len(nrow(b)); jj <- ncol(a) + seq_len(ncol(b))
  out[ii, jj] <- b
  out
}
.memwas_strip_parentheses <- function(x) {
  while (is.call(x) && identical(x[[1L]], as.name("("))) x <- x[[2L]]
  x
}
.memwas_find_named_argument <- function(raw, aliases) {
  nm <- names(raw); if (is.null(nm)) nm <- rep("", length(raw))
  idx <- which(nm %in% aliases)
  if (length(idx)) raw[[idx[[1L]]]] else NULL
}
.memwas_unnamed_arguments <- function(raw) {
  nm <- names(raw); if (is.null(nm)) nm <- rep("", length(raw))
  raw[nm == ""]
}
.memwas_evaluate_argument <- function(expr, envir, data = NULL, column_nse = FALSE) {
  if (is.null(expr)) return(NULL)
  if (column_nse && is.symbol(expr) && !is.null(data) && as.character(expr) %in% names(data)) {
    return(as.character(expr))
  }
  if (column_nse && is.character(expr) && length(expr) == 1L) return(expr)
  ans <- try(eval(expr, envir = envir), silent = TRUE)
  if (!inherits(ans, "try-error")) return(ans)
  if (!is.null(data)) {
    ans <- try(eval(expr, envir = data, enclos = envir), silent = TRUE)
    if (!inherits(ans, "try-error")) return(ans)
  }
  stop(sprintf("Could not evaluate argument: %s", paste(deparse(expr), collapse = " ")), call. = FALSE)
}
.memwas_resolve_data_value <- function(x, data, envir, what) {
  if (is.null(x)) return(NULL)
  if (.memwas_is_scalar_character(x)) {
    if (!x %in% names(data)) stop(sprintf("%s column '%s' was not found in data.", what, x), call. = FALSE)
    return(data[[x]])
  }
  if (is.symbol(x)) {
    nm <- as.character(x)
    if (nm %in% names(data)) return(data[[nm]])
  }
  if (is.language(x)) return(eval(x, data, envir))
  if (length(x) == 1L && is.atomic(x)) {
    nm <- as.character(x)
    if (nm %in% names(data) && length(data[[nm]]) == nrow(data)) return(data[[nm]])
  }
  x
}
.memwas_as_family_spec <- function(family = stats::gaussian(), link = NULL, theta = NULL, shape = NULL) {
  original <- family
  if (!is.null(theta) && (!is.numeric(theta) || length(theta) != 1L || !is.finite(theta) || theta <= 0)) stop("theta/size must be one positive number.", call. = FALSE)
  if (!is.null(shape) && (!is.numeric(shape) || length(shape) != 1L || !is.finite(shape) || shape <= 0)) stop("shape must be one positive number.", call. = FALSE)
  if (is.function(family) && !inherits(family, "family")) family <- family()
  if (inherits(family, "family")) {
    nm0 <- tolower(family$family)
    if (grepl("^quasi", nm0)) stop("Quasi families do not define a likelihood and cannot be used with Laplace or quadrature integration.", call. = FALSE)
    link_name <- family$link
    if (!is.null(link) && !identical(link, link_name)) {
      # Rebuild standard families with the requested link rather than discarding it.
      if (grepl("binomial", nm0) && !grepl("negative", nm0)) family <- stats::binomial(link)
      else if (grepl("poisson", nm0)) family <- stats::poisson(link)
      else if (grepl("gamma", nm0)) family <- stats::Gamma(link)
      else if (grepl("gaussian", nm0)) family <- stats::gaussian(link)
      else stop("The supplied custom family does not support replacing its link.", call. = FALSE)
      nm0 <- tolower(family$family); link_name <- family$link
    }
    if (grepl("negative[ _-]*binomial", nm0)) {
      parsed <- suppressWarnings(as.numeric(sub(".*\\(([^)]+)\\).*", "\\1", family$family)))
      fixed_theta <- if (is.finite(parsed)) parsed else theta
      if (!is.null(fixed_theta) && (!is.finite(fixed_theta) || fixed_theta <= 0)) stop("Negative-binomial size/theta must be positive.", call. = FALSE)
      lk <- stats::make.link(link_name)
      return(list(name = "negative_binomial", label = family$family, link = link_name,
                  linkfun = lk$linkfun, linkinv = lk$linkinv, mu.eta = lk$mu.eta,
                  valideta = lk$valideta, family_object = family,
                  fixed_theta = fixed_theta, fixed_shape = NULL))
    }
    if (grepl("exponential", nm0)) {
      lk <- stats::make.link(link_name)
      return(list(name = "exponential", label = family$family, link = link_name,
                  linkfun = lk$linkfun, linkinv = lk$linkinv, mu.eta = lk$mu.eta,
                  valideta = lk$valideta, family_object = family,
                  fixed_theta = NULL, fixed_shape = 1))
    }
    name <- if (grepl("binomial", nm0)) "binomial" else
      if (grepl("poisson", nm0)) "poisson" else
      if (grepl("gamma", nm0)) "gamma" else
      if (grepl("gaussian", nm0)) "gaussian" else NA_character_
    if (is.na(name)) stop(sprintf("Unsupported family '%s'.", family$family), call. = FALSE)
    return(list(name = name, label = family$family, link = family$link,
                linkfun = family$linkfun, linkinv = family$linkinv,
                mu.eta = family$mu.eta, valideta = family$valideta,
                family_object = family, fixed_theta = theta,
                fixed_shape = if (name == "gamma") shape else NULL))
  }
  if (!.memwas_is_scalar_character(family)) stop("family must be a family object, family function, or character string.", call. = FALSE)
  key <- tolower(gsub("[ .-]+", "_", family))
  key <- sub("^negativebinomial$", "negative_binomial", key)
  if (key %in% c("bernoulli", "binary")) key <- "binomial"
  if (key %in% c("nb", "nb2", "negbin", "negative_binomial_2")) key <- "negative_binomial"
  if (key %in% c("exp")) key <- "exponential"
  default_link <- switch(key, gaussian = "identity", binomial = "logit", poisson = "log",
                         negative_binomial = "log", gamma = "log", exponential = "log", NULL)
  if (is.null(default_link)) stop(sprintf("Unsupported family '%s'.", family), call. = FALSE)
  link <- .memwas_null_coalesce(link, default_link)
  if (key == "gaussian") return(.memwas_as_family_spec(stats::gaussian(link), theta = theta, shape = shape))
  if (key == "binomial") return(.memwas_as_family_spec(stats::binomial(link), theta = theta, shape = shape))
  if (key == "poisson") return(.memwas_as_family_spec(stats::poisson(link), theta = theta, shape = shape))
  if (key == "gamma") return(.memwas_as_family_spec(stats::Gamma(link), theta = theta, shape = shape))
  lk <- stats::make.link(link)
  list(name = key, label = family, link = link, linkfun = lk$linkfun,
       linkinv = lk$linkinv, mu.eta = lk$mu.eta, valideta = lk$valideta,
       family_object = structure(list(family = family, link = link,
                                      linkfun = lk$linkfun, linkinv = lk$linkinv,
                                      mu.eta = lk$mu.eta, valideta = lk$valideta),
                                 class = "family"),
       fixed_theta = theta,
       fixed_shape = if (key == "exponential") 1 else shape)
}
.memwas_family_spec_from_call <- function(raw, envir) {
  fexpr <- .memwas_find_named_argument(raw, c("family", "distribution"))
  lexpr <- .memwas_find_named_argument(raw, c("link"))
  thetaexpr <- .memwas_find_named_argument(raw, c("theta", "size", "nb_size"))
  shapeexpr <- .memwas_find_named_argument(raw, c("shape", "gamma_shape"))
  f <- if (is.null(fexpr)) stats::gaussian() else .memwas_evaluate_argument(fexpr, envir)
  link <- if (is.null(lexpr)) NULL else .memwas_evaluate_argument(lexpr, envir)
  theta <- if (is.null(thetaexpr)) NULL else .memwas_evaluate_argument(thetaexpr, envir)
  shape <- if (is.null(shapeexpr)) NULL else .memwas_evaluate_argument(shapeexpr, envir)
  .memwas_as_family_spec(f, link, theta, shape)
}
.memwas_parse_fit_call <- function(raw, envir, call) {
  unnamed <- .memwas_unnamed_arguments(raw)
  fexpr <- .memwas_find_named_argument(raw, c("formula", "fixed", "fixed_formula", "fixed.effects", "fixed_effects", "fixed_effects_formula", "fixed.formula", "model_formula"))
  dexpr <- .memwas_find_named_argument(raw, c("data", "data_frame", "dataset", "data_long", "long_data"))
  if (is.null(fexpr) && length(unnamed) >= 1L) fexpr <- unnamed[[1L]]
  if (is.null(dexpr) && length(unnamed) >= 2L) dexpr <- unnamed[[2L]]
  if (is.null(fexpr)) stop("A fixed-effects formula is required.", call. = FALSE)
  if (is.null(dexpr)) stop("A data frame is required.", call. = FALSE)
  data <- .memwas_evaluate_argument(dexpr, envir)
  if (!is.data.frame(data)) data <- as.data.frame(data)
  formula <- .memwas_evaluate_argument(fexpr, envir, data)
  if (!inherits(formula, "formula")) stop("formula/fixed must be a formula.", call. = FALSE)
  getv <- function(aliases, default = NULL, column_nse = FALSE) {
    z <- .memwas_find_named_argument(raw, aliases)
    if (is.null(z)) return(default)
    .memwas_evaluate_argument(z, envir, data, column_nse = column_nse)
  }
  control <- getv(c("control"), list())
  if (is.null(control)) control <- list()
  if (!is.list(control)) stop("control must be a list.", call. = FALSE)
  list(
    formula = formula, data = data,
    family = .memwas_family_spec_from_call(raw, envir),
    random = getv(c("random", "random_formula", "random.effects", "random_effects", "random_effects_formula", "random.formula", "re_formula")),
    id = getv(c("id", "subject", "subject_id", "cluster", "group", "subject_col", "subject_var", "id_var", "group_var", "cluster_id"), NULL, TRUE),
    time = getv(c("time", "time_var", "time_variable", "visit", "time_col", "visit_time", "occasion"), NULL, TRUE),
    autocor = getv(c("autocor", "autocorrelation", "correlation", "cor_struct", "correlation_structure", "correlation_type", "correlation_struct"), "NONE"),
    random_cov = getv(c("random_cov", "random_covariance", "covariance"), "unstructured"),
    approximation = getv(c("approximation"), "laplace"),
    method = getv(c("method"), "ML"),
    quadrature_points = getv(c("quadrature_points", "nAGQ", "nodes"), 9L),
    L1_penalty = getv(c("L1_penalty", "lambda1", "l1"), 0),
    L2_penalty = getv(c("L2_penalty", "lambda2", "l2"), 0),
    offset = getv(c("offset"), NULL, TRUE), weights = getv(c("weights", "prior_weights"), NULL, TRUE),
    subset = getv(c("subset"), NULL), na.action = getv(c("na.action"), stats::na.omit),
    theta = getv(c("theta", "size", "nb_size"), NULL),
    shape = getv(c("shape", "gamma_shape"), NULL),
    start = getv(c("start", "initial"), NULL), control = control,
    call = call, environment = envir
  )
}
.memwas_random_effects_spec <- function(random, id, data, envir) {
  if (is.list(random) && !inherits(random, "formula")) {
    if (!is.null(random$formula)) {
      if (is.null(id)) id <- .memwas_null_coalesce(random$group, random$id)
      random <- random$formula
    } else if (length(random) == 1L && inherits(random[[1L]], "formula")) {
      if (is.null(id) && nzchar(names(random)[1L])) id <- names(random)[1L]
      random <- random[[1L]]
    }
  }
  if (is.character(random) && length(random) > 1L) random <- stats::reformulate(random)
  if (.memwas_is_scalar_character(random) && grepl("~", random, fixed = TRUE)) random <- stats::as.formula(random, env = envir)
  if (!is.null(random) && !inherits(random, "formula")) stop("random must be a one-sided formula, a bar formula, or NULL.", call. = FALSE)
  design <- random; group_expr <- NULL
  if (!is.null(random)) {
    rhs <- .memwas_strip_parentheses(random[[length(random)]])
    if (is.call(rhs) && identical(rhs[[1L]], as.name("|"))) {
      group_expr <- rhs[[3L]]
      design <- stats::as.formula(call("~", rhs[[2L]]), env = environment(random))
    }
  }
  if (is.null(id) && !is.null(group_expr)) {
    id <- if (is.symbol(group_expr) && as.character(group_expr) %in% names(data)) as.character(group_expr) else eval(group_expr, data, envir)
  }
  list(formula = design, id = id)
}
.memwas_serial_covariance_spec <- function(x, id, time, control) {
  if (is.list(x)) {
    type <- .memwas_null_coalesce(x$type, .memwas_null_coalesce(x$structure, .memwas_null_coalesce(x$name, "NONE")))
    order <- .memwas_null_coalesce(x$order, x$p)
  } else { type <- x; order <- NULL }
  if (is.null(type)) type <- "NONE"
  type0 <- toupper(gsub("[ _-]+", "", as.character(type)[1L]))
  p <- NULL
  if (grepl("^AR\\([0-9]+\\)$", type0)) {
    p <- as.integer(sub("^AR\\(([0-9]+)\\)$", "\\1", type0)); type0 <- if (p == 1L) "AR1" else "ARP"
  }
  if (type0 %in% c("AR", "AUTOREGRESSIVE")) { p <- as.integer(.memwas_null_coalesce(order, 1L)); type0 <- if (p == 1L) "AR1" else "ARP" }
  if (type0 %in% c("AR1", "AR(1)")) type0 <- "AR1"
  if (type0 %in% c("ARMA11", "ARMA(1,1)")) type0 <- "ARMA11"
  if (type0 %in% c("COMPOUNDSYMMETRY", "EXCHANGEABLE")) type0 <- "CS"
  if (type0 %in% c("TOEPLITZ")) type0 <- "TOEP"
  if (type0 %in% c("UNSTRUCTURED")) type0 <- "UN"
  if (type0 %in% c("INDEPENDENT", "IID", "NULL", "NO", "NONE")) type0 <- "NONE"
  if (type0 %in% c("EXPONENTIAL", "OU")) type0 <- "EXP"
  allowed <- c("NONE", "AR1", "ARP", "ARMA11", "CS", "TOEP", "UN", "EXP")
  if (!type0 %in% allowed) stop(sprintf("Unsupported autocorrelation structure '%s'.", type), call. = FALSE)
  sizes <- table(id); max_m <- if (length(sizes)) max(sizes) else 1L
  if (is.null(time)) time <- ave(seq_along(id), id, FUN = seq_along)
  time <- as.numeric(time)
  gaps <- unlist(lapply(split(time, id), function(z) diff(sort(unique(z)))))
  positive <- gaps[is.finite(gaps) & gaps > 0]
  scale <- if (length(positive)) min(positive) else 1
  equal_grid <- all(vapply(split(time, id), function(z) {
    zz <- (z - min(z, na.rm = TRUE)) / scale
    all(abs(zz - round(zz)) < 1e-7)
  }, logical(1)))
  if (type0 %in% c("ARP", "ARMA11", "TOEP") && !equal_grid)
    stop(sprintf("%s requires equally spaced observation times; use AR(1)/EXP for irregular time.", type0), call. = FALSE)
  if (type0 == "ARP") p <- as.integer(.memwas_null_coalesce(p, .memwas_null_coalesce(order, 2L)))
  if (type0 == "TOEP") p <- as.integer(.memwas_null_coalesce(order, min(5L, max(1L, max_m - 1L))))
  if (!is.null(p) && (p < 1L || p > 20L)) stop("AR/Toeplitz order must be between 1 and 20.", call. = FALSE)
  levels_time <- sort(unique(time))
  max_un <- as.integer(.memwas_null_coalesce(control$max_unstructured_visits, 10L))
  if (type0 == "UN" && length(levels_time) > max_un)
    stop(sprintf("Unstructured serial covariance has %d visit times; increase control$max_unstructured_visits (currently %d) or use a parsimonious structure.", length(levels_time), max_un), call. = FALSE)
  list(type = type0, order = p, max_m = max_m, scale = scale,
       equal_grid = equal_grid, allow_negative_ar1 = equal_grid,
       time_levels = levels_time)
}
.memwas_family_pacf_to_ar <- function(kappa) {
  p <- length(kappa); if (!p) return(numeric())
  phi <- kappa[1L]
  if (p == 1L) return(phi)
  for (m in 2:p) {
    old <- phi
    phi <- numeric(m)
    phi[m] <- kappa[m]
    phi[seq_len(m - 1L)] <- old - kappa[m] * rev(old)
  }
  phi
}
.memwas_prepare_single_serial <- function(spec) {
  data <- spec$data
  if (!is.null(spec$subset)) {
    ss <- spec$subset
    if (is.language(ss)) ss <- eval(ss, data, spec$environment)
    data <- data[ss, , drop = FALSE]
  }
  mf <- stats::model.frame(spec$formula, data = data, na.action = stats::na.pass, drop.unused.levels = TRUE)
  tt <- stats::terms(mf)
  X <- stats::model.matrix(tt, mf)
  yraw <- stats::model.response(mf)
  n <- nrow(mf)
  off <- stats::model.offset(mf); if (is.null(off)) off <- rep(0, n)
  if (!is.null(spec$offset)) {
    oo <- .memwas_resolve_data_value(spec$offset, data, spec$environment, "offset")
    if (length(oo) == 1L) oo <- rep(oo, n)
    if (length(oo) != n) stop("offset must have one value per data row.", call. = FALSE)
    off <- off + as.numeric(oo)
  }
  w <- if (is.null(spec$weights)) rep(1, n) else .memwas_resolve_data_value(spec$weights, data, spec$environment, "weights")
  if (length(w) == 1L) w <- rep(w, n)
  if (length(w) != n || any(w < 0, na.rm = TRUE)) stop("weights must be nonnegative with one value per row.", call. = FALSE)
  rs <- .memwas_random_effects_spec(spec$random, spec$id, data, spec$environment)
  id <- rs$id
  if (!is.null(id)) id <- .memwas_resolve_data_value(id, data, spec$environment, "id")
  Z <- matrix(numeric(), n, 0L)
  random_terms <- NULL; random_xlevels <- NULL; random_contrasts <- NULL
  if (!is.null(rs$formula)) {
    random_mf <- stats::model.frame(rs$formula, data = data, na.action = stats::na.pass, drop.unused.levels = TRUE)
    random_terms <- stats::terms(random_mf)
    Z <- stats::model.matrix(random_terms, random_mf)
    random_xlevels <- stats::.getXlevels(random_terms, random_mf)
    random_contrasts <- attr(Z, "contrasts")
  }
  if (ncol(Z) && is.null(id)) stop("A subject/group id is required when random effects are specified.", call. = FALSE)
  if (is.null(id)) id <- rep(".all", n)
  if (length(id) != n) stop("id must have one value per data row.", call. = FALSE)
  time <- if (is.null(spec$time)) NULL else .memwas_resolve_data_value(spec$time, data, spec$environment, "time")
  if (!is.null(time) && length(time) != n) stop("time must have one value per data row.", call. = FALSE)
  # Family-aware response representation.
  fam <- spec$family
  freq <- as.numeric(w); trials <- rep(1, n); y <- yraw
  if (fam$name == "binomial") {
    if (is.character(yraw)) yraw <- factor(yraw)
    if (is.matrix(yraw) && ncol(yraw) == 2L) {
      y <- as.numeric(yraw[, 1L]); trials <- rowSums(yraw); freq <- as.numeric(w)
    } else if (is.factor(yraw)) {
      if (nlevels(yraw) != 2L) stop("Binomial factor responses must have exactly two levels.", call. = FALSE)
      y <- as.numeric(yraw) - 1L
    } else {
      y <- as.numeric(yraw)
      fractional <- is.finite(y) & y > 0 & y < 1
      if (any(fractional)) {
        if (any(abs(w - round(w)) > 1e-7 | w < 1, na.rm = TRUE))
          stop("Binomial proportions require integer trial totals supplied through weights or a cbind(success, failure) response.", call. = FALSE)
        trials <- as.numeric(w); y <- y * trials; freq <- rep(1, n)
      }
    }
  } else y <- as.numeric(yraw)
  if (is.null(time)) {
    time <- unsplit(lapply(split(seq_len(n), id), seq_along), id)
    time <- as.numeric(time)
  } else time <- as.numeric(time)
  ok <- stats::complete.cases(X, y, trials, freq, off, id, time) & freq > 0
  if (ncol(Z)) ok <- ok & stats::complete.cases(Z)
  if (!all(ok)) {
    X <- X[ok, , drop = FALSE]; Z <- Z[ok, , drop = FALSE]
    y <- y[ok]; trials <- trials[ok]; freq <- freq[ok]; off <- off[ok]
    id <- id[ok]; time <- time[ok]
  }
  if (!length(y)) stop("No complete observations remain after missing-value handling.", call. = FALSE)
  if (any(freq < 0) || any(trials <= 0)) stop("Frequency weights must be nonnegative and binomial totals positive.", call. = FALSE)
  if (fam$name == "binomial") {
    if (any(y < -1e-7 | y > trials + 1e-7) || any(abs(y - round(y)) > 1e-7) || any(abs(trials - round(trials)) > 1e-7))
      stop("Binomial successes and trial totals must be integer-valued with 0 <= successes <= trials.", call. = FALSE)
    y <- round(y); trials <- round(trials)
  }
  if (fam$name %in% c("poisson", "negative_binomial") && (any(y < 0) || any(abs(y - round(y)) > 1e-7)))
    stop(sprintf("%s responses must be nonnegative counts.", fam$label), call. = FALSE)
  if (fam$name %in% c("gamma", "exponential") && any(y <= 0))
    stop(sprintf("%s responses must be strictly positive.", fam$label), call. = FALSE)
  id <- factor(id)
  serial <- .memwas_serial_covariance_spec(spec$autocor, id, time, spec$control)
  list(y = y, yraw = yraw, trials = trials, freq = freq, X = X, Z = Z,
       offset = off, id = id, time = time, keep = which(ok), data = data,
       fixed_terms = tt, fixed_xlevels = stats::.getXlevels(tt, mf),
       fixed_contrasts = attr(X, "contrasts"), random_terms = random_terms,
       random_xlevels = random_xlevels, random_contrasts = random_contrasts,
       family = fam, serial = serial, id_spec = rs$id, time_spec = spec$time,
       response_name = as.character(attr(tt, "variables")[[2L]]))
}
.memwas_single_serial_layout <- function(prep, spec) {
  p <- ncol(prep$X); q <- ncol(prep$Z)
  random_cov <- .memwas_normalize_random_covariance(spec$random_cov)
  rn <- if (!q) 0L else if (random_cov == "diagonal") q else q * (q + 1L) / 2L
  st <- prep$serial$type
  sn <- if (identical(st, "NONE")) {
    0L
  } else if (st %in% c("AR1", "EXP", "CS")) {
    2L
  } else if (identical(st, "ARMA11")) {
    3L
  } else if (st %in% c("ARP", "TOEP")) {
    1L + prep$serial$order
  } else if (identical(st, "UN")) {
    m <- length(prep$serial$time_levels)
    as.integer(m * (m + 1L) / 2L)
  } else {
    stop("Unknown serial structure in layout: ", st, call. = FALSE)
  }
  fn <- if (prep$family$name == "gaussian") 1L else
    if (prep$family$name == "gamma" && is.null(prep$family$fixed_shape)) 1L else
    if (prep$family$name == "negative_binomial" && is.null(prep$family$fixed_theta)) 1L else 0L
  idx <- list(beta = seq_len(p))
  at <- p
  idx$random <- if (rn) at + seq_len(rn) else integer(); at <- at + rn
  idx$serial <- if (sn) at + seq_len(sn) else integer(); at <- at + sn
  idx$family <- if (fn) at + seq_len(fn) else integer(); at <- at + fn
  nms <- colnames(prep$X)
  if (rn) nms <- c(nms, paste0("cov_random_", seq_len(rn)))
  if (sn) nms <- c(nms, paste0("cov_serial_", seq_len(sn)))
  if (fn) nms <- c(nms, switch(prep$family$name, gaussian = "log_sigma", gamma = "log_shape", negative_binomial = "log_size"))
  list(p = p, q = q, random_cov = random_cov, random_n = rn, serial_n = sn,
       family_n = fn, idx = idx, npar = at, names = nms)
}
.memwas_decode_cholesky_parameters <- function(x, q, diagonal = FALSE) {
  if (!q) return(matrix(numeric(), 0L, 0L))
  if (diagonal) return(diag(exp(2 * x), q))
  L <- matrix(0, q, q); k <- 0L
  for (i in seq_len(q)) for (j in seq_len(i)) {
    k <- k + 1L; L[i, j] <- if (i == j) exp(x[k]) else x[k]
  }
  L %*% t(L)
}
.memwas_decode_single_serial_parameters <- function(par, prep, spec, layout) {
  beta <- par[layout$idx$beta]
  D <- .memwas_decode_cholesky_parameters(par[layout$idx$random], layout$q, layout$random_cov == "diagonal")
  sp <- par[layout$idx$serial]; ss <- prep$serial
  serial <- list(type = ss$type, time_levels = ss$time_levels, scale = ss$scale, order = ss$order)
  if (ss$type %in% c("AR1", "EXP")) {
    serial$sd <- exp(sp[1L])
    if (ss$type == "EXP") {
      serial$range <- exp(sp[2L])
      serial$rho <- exp(-1 / serial$range)
    } else {
      serial$rho <- if (ss$allow_negative_ar1) tanh(sp[2L]) else .memwas_inverse_logit(sp[2L])
    }
  } else if (ss$type == "CS") {
    lower <- if (ss$max_m <= 1L) 0 else -1 / (ss$max_m - 1L) + 1e-6
    serial$sd <- exp(sp[1L]); serial$rho <- lower + (0.999 - lower) * .memwas_inverse_logit(sp[2L])
  } else if (ss$type == "ARMA11") {
    serial$sd <- exp(sp[1L]); serial$phi <- tanh(sp[2L]); serial$theta <- tanh(sp[3L])
  } else if (ss$type %in% c("ARP", "TOEP")) {
    serial$sd <- exp(sp[1L]); serial$pacf <- tanh(sp[-1L]); serial$ar <- .memwas_family_pacf_to_ar(serial$pacf)
  } else if (ss$type == "UN") {
    m <- length(ss$time_levels); serial$Sigma_full <- .memwas_decode_cholesky_parameters(sp, m, FALSE)
  }
  fampar <- list()
  if (prep$family$name == "gaussian") fampar$sigma <- exp(par[layout$idx$family])
  if (prep$family$name == "gamma") fampar$shape <- if (is.null(prep$family$fixed_shape)) exp(par[layout$idx$family]) else prep$family$fixed_shape
  if (prep$family$name == "negative_binomial") fampar$theta <- if (is.null(prep$family$fixed_theta)) exp(par[layout$idx$family]) else prep$family$fixed_theta
  if (prep$family$name == "exponential") fampar$shape <- 1
  list(beta = beta, D = D, serial = serial, family = fampar)
}
.memwas_serial_covariance_matrix <- function(s, time1, time2 = time1) {
  n1 <- length(time1); n2 <- length(time2)
  if (s$type == "NONE") return(matrix(numeric(), n1, n2))
  if (s$type == "UN") {
    i <- match(time1, s$time_levels); j <- match(time2, s$time_levels)
    if (anyNA(i) || anyNA(j)) stop("Prediction at unseen times is unavailable for an unstructured serial covariance.", call. = FALSE)
    return(s$Sigma_full[i, j, drop = FALSE])
  }
  d <- abs(outer(time1, time2, "-") / s$scale)
  if (s$type == "EXP") return(s$sd^2 * exp(-d / max(s$range, 1e-8)))
  if (s$type == "AR1") {
    rho <- s$rho
    if (rho < 0) {
      di <- round(d)
      R <- (abs(rho)^d) * ifelse(di %% 2L, -1, 1)
    } else R <- rho^d
    return(s$sd^2 * R)
  }
  if (s$type == "CS") return(s$sd^2 * ifelse(d < 1e-12, 1, s$rho))
  lag <- as.integer(round(d)); maxlag <- max(lag)
  if (s$type == "ARMA11") ac <- stats::ARMAacf(ar = s$phi, ma = s$theta, lag.max = maxlag)
  else ac <- stats::ARMAacf(ar = s$ar, lag.max = maxlag)
  s$sd^2 * matrix(ac[lag + 1L], n1, n2)
}
.memwas_family_mean <- function(eta, fam) {
  mu <- suppressWarnings(fam$linkinv(eta))
  valid_eta <- try(fam$valideta(eta), silent = TRUE)
  if (!inherits(valid_eta, "try-error") && length(valid_eta) == length(mu)) mu[!valid_eta] <- NA_real_
  if (fam$name == "binomial") {
    mu[!is.finite(mu) | mu < 0 | mu > 1] <- NA_real_
    good <- is.finite(mu)
    mu[good] <- pmin(pmax(mu[good], 1e-12), 1 - 1e-12)
  }
  if (fam$name %in% c("poisson", "negative_binomial", "gamma", "exponential"))
    mu[!is.finite(mu) | mu <= 0] <- NA_real_
  mu
}
.memwas_observation_loglikelihood <- function(eta, y, trials, freq, fam, fp) {
  mu <- .memwas_family_mean(eta, fam)
  bad <- !is.finite(mu)
  ans <- switch(fam$name,
    gaussian = stats::dnorm(y, mean = mu, sd = fp$sigma, log = TRUE),
    binomial = stats::dbinom(y, size = trials, prob = mu, log = TRUE),
    poisson = stats::dpois(y, lambda = mu, log = TRUE),
    negative_binomial = stats::dnbinom(y, size = fp$theta, mu = mu, log = TRUE),
    gamma = stats::dgamma(y, shape = fp$shape, scale = mu / fp$shape, log = TRUE),
    exponential = stats::dexp(y, rate = 1 / mu, log = TRUE))
  ans <- freq * ans
  ans[bad | !is.finite(ans)] <- -1e100
  ans
}
.memwas_score_weight <- function(eta, y, trials, freq, fam, fp) {
  mu <- .memwas_family_mean(eta, fam)
  dm <- suppressWarnings(fam$mu.eta(eta)); dm[!is.finite(dm)] <- 0
  tiny <- 1e-12
  score <- weight <- rep(NA_real_, length(eta))
  if (fam$name == "gaussian") {
    score <- freq * (y - mu) * dm / fp$sigma^2
    weight <- freq * dm^2 / fp$sigma^2
  } else if (fam$name == "binomial") {
    v <- pmax(mu * (1 - mu), tiny)
    score <- freq * (y - trials * mu) * dm / v
    weight <- freq * trials * dm^2 / v
  } else if (fam$name == "poisson") {
    score <- freq * (y - mu) * dm / pmax(mu, tiny)
    weight <- freq * dm^2 / pmax(mu, tiny)
  } else if (fam$name == "negative_binomial") {
    den <- pmax(mu * (fp$theta + mu), tiny)
    score <- freq * fp$theta * (y - mu) * dm / den
    weight <- freq * fp$theta * dm^2 / den
  } else if (fam$name == "gamma") {
    score <- freq * fp$shape * (y - mu) * dm / pmax(mu^2, tiny)
    weight <- freq * fp$shape * dm^2 / pmax(mu^2, tiny)
  } else {
    score <- freq * (y - mu) * dm / pmax(mu^2, tiny)
    weight <- freq * dm^2 / pmax(mu^2, tiny)
  }
  score[!is.finite(score)] <- 0
  weight[!is.finite(weight) | weight < 1e-10] <- 1e-10
  list(score = score, weight = pmin(weight, 1e12), mu = mu)
}
.memwas_observed_curvature <- function(eta, y, trials, freq, fam, fp, fallback) {
  h <- 1e-4 * (1 + abs(eta))
  f0 <- .memwas_observation_loglikelihood(eta, y, trials, freq, fam, fp)
  fpv <- .memwas_observation_loglikelihood(eta + h, y, trials, freq, fam, fp)
  fmv <- .memwas_observation_loglikelihood(eta - h, y, trials, freq, fam, fp)
  curv <- -(fpv - 2 * f0 + fmv) / h^2
  bad <- !is.finite(curv) | abs(curv) > 1e12
  curv[bad] <- fallback[bad]
  curv
}
.memwas_optimize_latent_mode <- function(base_eta, A, C, y, trials, freq, fam, fp, start = NULL, control = list()) {
  k <- ncol(A)
  if (!k) return(list(mode = numeric(), joint = -.memwas_null_coalesce(sum(.memwas_observation_loglikelihood(base_eta, y, trials, freq, fam, fp)), -Inf), H = matrix(numeric(),0,0), eta = base_eta, converged = TRUE))
  Q <- .memwas_solve_pd(C); logdetC <- .memwas_logdet_pd(C)
  a <- if (is.null(start) || length(start) != k || any(!is.finite(start))) rep(0, k) else start
  joint <- function(a) {
    eta <- base_eta + as.vector(A %*% a)
    ll <- sum(.memwas_observation_loglikelihood(eta, y, trials, freq, fam, fp))
    if (!is.finite(ll)) return(1e100)
    -ll + 0.5 * sum(a * (Q %*% a)) + 0.5 * logdetC + 0.5 * k * log(2 * pi)
  }
  maxit <- as.integer(.memwas_null_coalesce(control$inner_maxit, 60L)); tol <- .memwas_null_coalesce(control$inner_tol, 1e-8)
  old <- joint(a); converged <- FALSE
  for (it in seq_len(maxit)) {
    eta <- base_eta + as.vector(A %*% a)
    sw <- .memwas_score_weight(eta, y, trials, freq, fam, fp)
    g <- -as.vector(crossprod(A, sw$score)) + as.vector(Q %*% a)
    Hf <- crossprod(A, A * sw$weight) + Q
    step <- try(.memwas_solve_pd(Hf, g), silent = TRUE)
    if (inherits(step, "try-error") || any(!is.finite(step))) break
    if (max(abs(g)) < tol || max(abs(step)) < tol * (1 + max(abs(a)))) { converged <- TRUE; break }
    alpha <- 1
    improved <- FALSE
    for (ls in 0:20) {
      cand <- a - alpha * step; val <- joint(cand)
      if (is.finite(val) && val <= old - 1e-4 * alpha * sum(g * step)) { a <- cand; old <- val; improved <- TRUE; break }
      alpha <- alpha / 2
    }
    if (!improved) break
  }
  if (!converged) {
    oo <- try(stats::optim(a, joint, method = "BFGS", control = list(maxit = maxit, reltol = tol)), silent = TRUE)
    if (!inherits(oo, "try-error") && is.finite(oo$value)) { a <- oo$par; old <- oo$value; converged <- oo$convergence == 0L }
  }
  eta <- base_eta + as.vector(A %*% a)
  sw <- .memwas_score_weight(eta, y, trials, freq, fam, fp)
  curv <- .memwas_observed_curvature(eta, y, trials, freq, fam, fp, sw$weight)
  H <- crossprod(A, A * curv) + Q
  exact_pd <- !inherits(try(chol((H + t(H))/2), silent = TRUE), "try-error")
  if (!exact_pd) H <- crossprod(A, A * sw$weight) + Q
  .memwas_chol_pd(H)
  list(mode = a, joint = joint(a), H = H, eta = eta, converged = converged,
       hessian = if (exact_pd) "observed" else "fisher_fallback")
}
.memwas_gauss_hermite_rule <- function(n) {
  n <- as.integer(n); if (n < 3L) n <- 3L
  J <- matrix(0, n, n)
  if (n > 1L) for (i in seq_len(n - 1L)) J[i, i + 1L] <- J[i + 1L, i] <- sqrt(i / 2)
  ee <- eigen(J, symmetric = TRUE)
  ord <- order(ee$values)
  list(nodes = ee$values[ord], weights = sqrt(pi) * ee$vectors[1L, ord]^2)
}
.memwas_quadrature_negative_loglikelihood <- function(base_eta, A, C, y, trials, freq, fam, fp, nodes,
                                           mode_info = NULL, control = list()) {
  k <- ncol(A)
  if (!k) return(-sum(.memwas_observation_loglikelihood(base_eta, y, trials, freq, fam, fp)))
  if (k > 2L) stop("Adaptive Gauss-Hermite quadrature is limited to at most two latent dimensions per subject; use approximation='laplace' for serial latent processes or larger random-effects vectors.", call. = FALSE)
  if (is.null(mode_info))
    mode_info <- .memwas_optimize_latent_mode(base_eta, A, C, y, trials, freq, fam, fp, control = control)
  gh <- .memwas_gauss_hermite_rule(nodes)
  Q <- .memwas_solve_pd(C); logdetC <- .memwas_logdet_pd(C)
  RH <- .memwas_chol_pd(mode_info$H)
  Hinvhalf <- backsolve(RH, diag(k))
  logdetH <- 2 * sum(log(diag(RH)))
  log_joint <- function(a) {
    eta <- base_eta + as.vector(A %*% a)
    sum(.memwas_observation_loglikelihood(eta, y, trials, freq, fam, fp)) -
      0.5 * sum(a * (Q %*% a)) - 0.5 * logdetC - 0.5 * k * log(2*pi)
  }
  if (k == 1L) {
    vals <- vapply(seq_along(gh$nodes), function(i) {
      x <- gh$nodes[i]
      a <- mode_info$mode + sqrt(2) * as.vector(Hinvhalf %*% x)
      log(gh$weights[i]) + log_joint(a) + x^2
    }, numeric(1))
  } else {
    grid <- expand.grid(i = seq_along(gh$nodes), j = seq_along(gh$nodes))
    vals <- vapply(seq_len(nrow(grid)), function(r) {
      ii <- grid$i[r]; jj <- grid$j[r]; x <- c(gh$nodes[ii], gh$nodes[jj])
      a <- mode_info$mode + sqrt(2) * as.vector(Hinvhalf %*% x)
      log(gh$weights[ii]) + log(gh$weights[jj]) + log_joint(a) + sum(x^2)
    }, numeric(1))
  }
  -(0.5 * k * log(2) - 0.5 * logdetH + .memwas_family_logsumexp(vals))
}

.memwas_single_serial_subject_components <- function(prep, dec, idx) {
  Zi <- prep$Z[idx, , drop = FALSE]; ti <- prep$time[idx]
  mats <- list(); covs <- list()
  if (ncol(Zi)) { mats[[length(mats)+1L]] <- Zi; covs[[length(covs)+1L]] <- dec$D }
  if (dec$serial$type != "NONE") {
    mats[[length(mats)+1L]] <- diag(length(idx))
    covs[[length(covs)+1L]] <- .memwas_serial_covariance_matrix(dec$serial, ti)
  }
  A <- if (length(mats)) do.call(cbind, mats) else matrix(numeric(), length(idx), 0L)
  C <- matrix(numeric(), 0L, 0L)
  for (cc in covs) C <- .memwas_block_diag(C, cc)
  list(A = A, C = C, q = ncol(Zi), r = if (dec$serial$type == "NONE") 0L else length(idx))
}
.memwas_single_serial_initial_values <- function(prep, spec, layout) {
  fam <- prep$family
  beta <- rep(0, layout$p); names(beta) <- colnames(prep$X)
  fit0 <- try({
    if (fam$name == "binomial") {
      yy <- cbind(prep$y, prep$trials - prep$y)
      stats::glm.fit(prep$X, yy, family = stats::binomial(fam$link), offset = prep$offset, weights = prep$freq)
    } else {
      gf <- if (fam$name == "gaussian") stats::gaussian(fam$link) else if (fam$name == "poisson") stats::poisson(fam$link) else if (fam$name == "gamma" || fam$name == "exponential") stats::Gamma(fam$link) else stats::poisson("log")
      stats::glm.fit(prep$X, prep$y, family = gf, offset = prep$offset, weights = prep$freq)
    }
  }, silent = TRUE)
  if (!inherits(fit0, "try-error") && length(fit0$coefficients) == layout$p) {
    beta <- fit0$coefficients; beta[!is.finite(beta)] <- 0
  }
  par <- numeric(layout$npar); par[layout$idx$beta] <- beta
  if (layout$random_n) {
    if (layout$random_cov == "diagonal") par[layout$idx$random] <- log(0.5)
    else {
      k <- 0L; vals <- numeric(layout$random_n)
      for (i in seq_len(layout$q)) for (j in seq_len(i)) { k <- k + 1L; vals[k] <- if (i == j) log(0.5) else 0 }
      par[layout$idx$random] <- vals
    }
  }
  if (layout$serial_n) {
    st <- prep$serial$type
    if (st == "UN") {
      m <- length(prep$serial$time_levels); vals <- numeric(layout$serial_n); k <- 0L
      for (i in seq_len(m)) for (j in seq_len(i)) { k <- k + 1L; vals[k] <- if (i == j) log(0.35) else 0 }
      par[layout$idx$serial] <- vals
    } else { par[layout$idx$serial] <- 0; par[layout$idx$serial[1L]] <- log(0.35) }
  }
  if (layout$family_n) {
    if (fam$name == "gaussian") par[layout$idx$family] <- log(max(stats::sd(prep$y), 0.1))
    if (fam$name == "gamma") par[layout$idx$family] <- log(2)
    if (fam$name == "negative_binomial") par[layout$idx$family] <- log(5)
  }
  if (!is.null(spec$start)) {
    st <- spec$start
    if (is.list(st)) {
      if (!is.null(st$beta)) {
        par[layout$idx$beta] <- .memwas_validate_raw_start_vector(
          st$beta, layout$p, "start$beta"
        )
      }
    } else {
      par <- .memwas_validate_raw_start_vector(st, length(par), "start")
    }
  }
  names(par) <- layout$names
  par
}
.memwas_fit_single_serial_engine <- function(spec) {
  if (is.null(spec$control)) spec$control <- list()
  if (!is.list(spec$control)) stop("control must be a list.", call. = FALSE)
  .memwas_check_scalar_nonnegative(spec$L1_penalty, "L1_penalty")
  .memwas_check_scalar_nonnegative(spec$L2_penalty, "L2_penalty")
  spec$quadrature_points <- .memwas_validate_quadrature_points(spec$quadrature_points)
  method <- .memwas_normalize_likelihood_method(.memwas_null_coalesce(spec$method, "ML"))
  if (!identical(method, "ML")) {
    stop(
      "The joint single-serial likelihood engine supports method = 'ML' only; ",
      "Gaussian REML remains available through the legacy Gaussian path.",
      call. = FALSE
    )
  }
  prep <- .memwas_prepare_single_serial(spec)
  layout <- .memwas_single_serial_layout(prep, spec)
  if (!.memwas_is_scalar_character(spec$approximation)) {
    stop("approximation must be 'laplace' or 'quadrature'.", call. = FALSE)
  }
  approx <- tolower(gsub("[ _-]+", "", spec$approximation))
  if (approx %in% c("laplace", "lap")) approx <- "laplace"
  else if (approx %in% c("quadrature", "aghq", "adaptivegausshermite",
                         "adaptivegausshermitequadrature",
                         "adaptivegaussianquadrature",
                         "gausshermite", "ghq")) approx <- "quadrature"
  else stop("approximation must be 'laplace' or 'quadrature'.", call. = FALSE)
  groups <- split(seq_along(prep$y), prep$id)
  max_latent <- max(vapply(groups, function(ii) layout$q + if (prep$serial$type == "NONE") 0L else length(ii), integer(1)))
  if (approx == "quadrature" && max_latent > 2L)
    stop("Quadrature would require more than two latent dimensions per subject. Use approximation='laplace'; the request is rejected rather than silently fitting a different model.", call. = FALSE)
  init <- .memwas_single_serial_initial_values(prep, spec, layout)
  cache <- new.env(parent = emptyenv()); cache$modes <- vector("list", length(groups)); names(cache$modes) <- names(groups)
  penalty_idx <- setdiff(seq_len(layout$p), which(colnames(prep$X) == "(Intercept)"))
  if (isTRUE(spec$control$penalize_intercept)) penalty_idx <- seq_len(layout$p)
  evaluate <- function(par, details = FALSE, penalized = TRUE) {
    dec <- try(.memwas_decode_single_serial_parameters(par, prep, spec, layout), silent = TRUE)
    if (inherits(dec, "try-error")) return(if (details) NULL else 1e50)
    total <- 0; out <- vector("list", length(groups)); names(out) <- names(groups)
    for (g in seq_along(groups)) {
      ii <- groups[[g]]; sc <- try(.memwas_single_serial_subject_components(prep, dec, ii), silent = TRUE)
      if (inherits(sc, "try-error")) return(if (details) NULL else 1e50)
      base_eta <- as.vector(prep$X[ii, , drop = FALSE] %*% dec$beta) + prep$offset[ii]
      if (approx == "laplace") {
        im <- try(.memwas_optimize_latent_mode(base_eta, sc$A, sc$C, prep$y[ii], prep$trials[ii], prep$freq[ii], prep$family, dec$family, cache$modes[[g]], spec$control), silent = TRUE)
        if (inherits(im, "try-error")) return(if (details) NULL else 1e50)
        val <- im$joint + if (ncol(sc$A)) 0.5 * .memwas_logdet_pd(im$H) - 0.5 * ncol(sc$A) * log(2*pi) else 0
        cache$modes[[g]] <- im$mode
      } else {
        im <- try(.memwas_optimize_latent_mode(base_eta, sc$A, sc$C, prep$y[ii], prep$trials[ii], prep$freq[ii], prep$family, dec$family, cache$modes[[g]], spec$control), silent = TRUE)
        if (inherits(im, "try-error")) return(if (details) NULL else 1e50)
        val <- try(.memwas_quadrature_negative_loglikelihood(base_eta, sc$A, sc$C, prep$y[ii], prep$trials[ii], prep$freq[ii], prep$family, dec$family, spec$quadrature_points, mode_info = im, control = spec$control), silent = TRUE)
        if (inherits(val, "try-error") || !is.finite(val)) return(if (details) NULL else 1e50)
        cache$modes[[g]] <- im$mode
      }
      if (!is.finite(val)) return(if (details) NULL else 1e50)
      total <- total + val; out[[g]] <- c(im, list(idx = ii, components = sc))
    }
    penalty <- spec$L1_penalty * sum(abs(dec$beta[penalty_idx])) + 0.5 * spec$L2_penalty * sum(dec$beta[penalty_idx]^2)
    if (details) list(nll = total, penalty = penalty, objective = total + penalty, decoded = dec, subjects = out) else total + if (penalized) penalty else 0
  }
  maxit <- as.integer(.memwas_null_coalesce(spec$control$outer_maxit, 150L))
  evalmax <- as.integer(.memwas_null_coalesce(spec$control$outer_eval_max, max(500L, 5L * maxit)))
  reltol <- .memwas_null_coalesce(spec$control$outer_tol, 1e-7)
  opt <- try(stats::nlminb(init, objective = evaluate,
                           control = list(iter.max = maxit, eval.max = evalmax, rel.tol = reltol, x.tol = reltol)), silent = TRUE)
  if (inherits(opt, "try-error")) {
    oo <- stats::optim(init, evaluate, method = "BFGS", control = list(maxit = maxit, reltol = reltol))
    opt <- list(par = oo$par, objective = oo$value, convergence = oo$convergence, message = oo$message, iterations = oo$counts)
  }
  details <- evaluate(opt$par, details = TRUE)
  if (is.null(details)) stop("The final likelihood evaluation failed.", call. = FALSE)
  dec <- details$decoded
  eta_fixed <- as.vector(prep$X %*% dec$beta) + prep$offset
  eta_cond <- eta_fixed; random_modes <- matrix(NA_real_, nlevels(prep$id), layout$q, dimnames = list(levels(prep$id), colnames(prep$Z)))
  serial_modes <- rep(0, length(prep$y))
  for (g in seq_along(details$subjects)) {
    ss <- details$subjects[[g]]; ii <- ss$idx; a <- ss$mode
    if (layout$q) random_modes[g, ] <- a[seq_len(layout$q)]
    if (prep$serial$type != "NONE") serial_modes[ii] <- a[layout$q + seq_along(ii)]
    if (length(a)) eta_cond[ii] <- eta_fixed[ii] + as.vector(ss$components$A %*% a)
  }
  mu_fixed <- .memwas_family_mean(eta_fixed, prep$family); mu_cond <- .memwas_family_mean(eta_cond, prep$family)
  V <- matrix(NA_real_, layout$npar, layout$npar, dimnames = list(layout$names, layout$names))
  compute_vcov <- isTRUE(.memwas_null_coalesce(spec$control$compute_vcov, layout$npar <= 30L && spec$L1_penalty == 0))
  if (compute_vcov) {
    hh <- try(stats::optimHess(opt$par, evaluate), silent = TRUE)
    if (!inherits(hh, "try-error") && all(is.finite(hh))) {
      vv <- try(.memwas_solve_pd(hh), silent = TRUE)
      if (!inherits(vv, "try-error")) V <- vv
    }
  }
  loglik <- -details$nll; df <- layout$npar
  structure(list(
    call = spec$call, formula = spec$formula, random = spec$random,
    family = prep$family$family_object, family_name = prep$family$name, link = prep$family$link,
    method = method, approximation = approx, autocor = prep$serial$type,
    serial_structure = prep$serial$type,
    coefficients = stats::setNames(dec$beta, colnames(prep$X)), fixed_effects = stats::setNames(dec$beta, colnames(prep$X)), beta = dec$beta,
    vcov = V[layout$idx$beta, layout$idx$beta, drop = FALSE], vcov_full = V,
    random_effects = random_modes, ranef = random_modes, serial_effects = serial_modes,
    random_covariance = dec$D, serial_parameters = dec$serial, family_parameters = dec$family,
    linear_predictor = eta_cond, linear_predictor_fixed = eta_fixed,
    fitted.values = mu_cond, fitted_values = mu_cond, fitted_fixed = mu_fixed,
    residuals = if (prep$family$name == "binomial") prep$y / prep$trials - mu_cond else prep$y - mu_cond,
    response = prep$y, trials = prep$trials, frequency_weights = prep$freq, offset = prep$offset,
    id = prep$id, time = prep$time, X = prep$X, Z = prep$Z,
    logLik = loglik, negLogLik = details$nll, objective = details$objective, penalty = details$penalty,
    AIC = -2 * loglik + 2 * df, BIC = -2 * loglik + log(length(prep$y)) * df,
    df = df, nobs = length(prep$y), convergence = opt$convergence,
    converged = isTRUE(opt$convergence == 0L), message = .memwas_null_coalesce(opt$message, ""), iterations = opt$iterations,
    L1_penalty = spec$L1_penalty, L2_penalty = spec$L2_penalty,
    terms = prep$fixed_terms, xlevels = prep$fixed_xlevels, contrasts = prep$fixed_contrasts,
    random_terms = prep$random_terms, random_xlevels = prep$random_xlevels, random_contrasts = prep$random_contrasts,
    training_data = prep$data, kept_rows = prep$keep, engine = "single_serial_joint_likelihood",
    .prep = prep, .layout = layout, .decoded = dec, .spec = spec
  ), class = c("MEMWAS_family_fit", "MEMWAS_fit"))
}

# Gaussian identity-link requests continue through the original covariance
# engine unless control$force_repaired = TRUE.
fit_MEMWAS <- function(...) {
  envir <- parent.frame()
  raw <- .memwas_normalize_fit_call(as.list(match.call(expand.dots = FALSE)$...), envir)
  fam <- .memwas_family_spec_from_call(raw, envir)
  cexpr <- .memwas_find_named_argument(raw, c("control")); ctl <- if (is.null(cexpr)) list() else try(.memwas_evaluate_argument(cexpr, envir), silent = TRUE)
  force <- is.list(ctl) && isTRUE(ctl$force_repaired)
  legacy_ok <- !force && fam$name == "gaussian" && identical(fam$link, "identity") && !is.null(.memwas_legacy_fit_MEMWAS)
  if (legacy_ok) return(.memwas_legacy_fit_MEMWAS(...))
  spec <- .memwas_parse_fit_call(raw, envir, match.call())
  .memwas_fit_single_serial_engine(spec)
}

#' @export
coef.MEMWAS_family_fit <- function(object, ...) object$coefficients

#' @export
vcov.MEMWAS_family_fit <- function(object, full = FALSE, ...) if (isTRUE(full)) object$vcov_full else object$vcov

#' @export
logLik.MEMWAS_family_fit <- function(object, ...) structure(object$logLik, df = object$df, nobs = object$nobs, class = "logLik")

.nobs.MEMWAS_family_fit <- function(object, ...) object$nobs

#' @export
fitted.MEMWAS_family_fit <- function(object, fixed_only = FALSE, ...) if (isTRUE(fixed_only)) object$fitted_fixed else object$fitted.values

.ranef.MEMWAS_family_fit <- function(object, ...) object$random_effects

#' @export
print.MEMWAS_family_fit <- function(x, ...) {
  cat("MEMWAS corrected multi-family fit\n")
  cat(" Family:", x$family_name, " (", x$link, ")\n", sep = "")
  cat(" Approximation:", x$approximation, "  Serial:", x$autocor, "\n")
  cat(" Observations:", x$nobs, "  logLik:", format(x$logLik, digits = 6), "\n")
  cat(" Converged:", x$converged, if (nzchar(x$message)) paste0(" (", x$message, ")") else "", "\n")
  print(x$coefficients)
  invisible(x)
}

#' @export
summary.MEMWAS_family_fit <- function(object, ...) {
  se <- sqrt(pmax(diag(object$vcov), 0)); est <- object$coefficients
  if (length(se) != length(est)) se <- rep(NA_real_, length(est))
  z <- est / se; tab <- cbind(Estimate = est, `Std. Error` = se, `z value` = z, `Pr(>|z|)` = 2 * stats::pnorm(abs(z), lower.tail = FALSE))
  out <- list(call = object$call, family = object$family_name, link = object$link,
              approximation = object$approximation, autocor = object$autocor,
              coefficients = tab, logLik = object$logLik, AIC = object$AIC, BIC = object$BIC,
              convergence = object$convergence, message = object$message,
              random_covariance = object$random_covariance,
              serial_parameters = object$serial_parameters,
              family_parameters = object$family_parameters)
  class(out) <- "summary.MEMWAS_family_fit"; out
}

#' @export
print.summary.MEMWAS_family_fit <- function(x, ...) {
  cat("MEMWAS corrected multi-family model\n")
  cat("Family:", x$family, "  Link:", x$link, "  Approximation:", x$approximation, "\n")
  cat("Serial structure:", x$autocor, "\n\nFixed effects:\n")
  stats::printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
  cat("\nlogLik:", format(x$logLik, digits = 6), " AIC:", format(x$AIC, digits = 6), " BIC:", format(x$BIC, digits = 6), "\n")
  invisible(x)
}

#' @export
residuals.MEMWAS_family_fit <- function(object, type = c("response", "pearson", "deviance"), fixed_only = FALSE, ...) {
  type <- match.arg(type); mu <- if (fixed_only) object$fitted_fixed else object$fitted.values
  y <- object$response; fam <- object$family_name
  if (fam == "binomial") { obs <- y / object$trials; raw <- obs - mu; variance <- mu * (1 - mu) / object$trials }
  else { raw <- y - mu; variance <- switch(fam, gaussian = object$family_parameters$sigma^2,
      poisson = mu, negative_binomial = mu + mu^2 / object$family_parameters$theta,
      gamma = mu^2 / object$family_parameters$shape, exponential = mu^2) }
  if (type == "response") return(raw)
  if (type == "pearson") return(raw / sqrt(pmax(variance, 1e-12)))
  # Signed square-root likelihood deviance, computed from saturated likelihood where available.
  if (fam == "gaussian") return(raw / object$family_parameters$sigma)
  if (fam == "poisson") {
    d <- 2 * ifelse(y == 0, mu, y * log(y / mu) - (y - mu))
  } else if (fam == "binomial") {
    n <- object$trials; p <- pmin(pmax(mu, 1e-12), 1-1e-12)
    d <- 2 * (ifelse(y == 0, 0, y * log(y/(n*p))) + ifelse(y == n, 0, (n-y)*log((n-y)/(n*(1-p)))))
  } else return(raw / sqrt(pmax(variance, 1e-12)))
  sign(raw) * sqrt(pmax(d, 0))
}
.memwas_new_model_design <- function(object, newdata) {
  mf <- stats::model.frame(stats::delete.response(object$terms), newdata,
                           na.action = stats::na.pass, xlev = object$xlevels)
  X <- stats::model.matrix(stats::delete.response(object$terms), mf, contrasts.arg = object$contrasts)
  off <- stats::model.offset(mf); if (is.null(off)) off <- rep(0, nrow(X))
  Z <- matrix(numeric(), nrow(X), 0L)
  if (!is.null(object$random_terms)) {
    rmf <- stats::model.frame(object$random_terms, newdata, na.action = stats::na.pass, xlev = object$random_xlevels)
    Z <- stats::model.matrix(object$random_terms, rmf, contrasts.arg = object$random_contrasts)
  }
  list(X = X, Z = Z, offset = off)
}

#' @export
predict.MEMWAS_family_fit <- function(object, newdata = NULL, type = c("response", "link"),
                                      include_random = FALSE, include_serial = include_random,
                                      id = NULL, time = NULL, offset = NULL, allow_new_levels = TRUE, ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    eta <- object$linear_predictor_fixed
    if (include_random && ncol(object$Z)) {
      b_by_row <- object$random_effects[as.character(object$id), , drop = FALSE]
      eta <- eta + rowSums(object$Z * b_by_row)
    }
    if (include_serial && object$autocor != "NONE") eta <- eta + object$serial_effects
    return(if (type == "link") eta else .memwas_family_mean(eta, object$.prep$family))
  }
  nd <- as.data.frame(newdata); des <- .memwas_new_model_design(object, nd)
  eta <- as.vector(des$X %*% object$coefficients) + des$offset
  offset_source <- .memwas_null_coalesce(offset, object$.spec$offset)
  if (!is.null(offset_source)) {
    oo <- try(.memwas_resolve_data_value(offset_source, nd, parent.frame(), "offset"), silent = TRUE)
    if (inherits(oo, "try-error")) stop("Could not evaluate the newdata offset; supply offset= explicitly.", call. = FALSE)
    if (length(oo) == 1L) oo <- rep(oo, nrow(nd))
    if (length(oo) != nrow(nd)) stop("offset must have one value per newdata row.", call. = FALSE)
    eta <- eta + as.numeric(oo)
  }
  idv <- id
  id_source <- .memwas_null_coalesce(object$.spec$id, object$.prep$id_spec)
  if (is.null(idv) && !is.null(id_source)) idv <- try(.memwas_resolve_data_value(id_source, nd, parent.frame(), "id"), silent = TRUE)
  if (inherits(idv, "try-error")) idv <- NULL
  if ((include_random || include_serial) && is.null(idv)) stop("Conditional prediction requires id values for newdata.", call. = FALSE)
  if (!is.null(idv) && length(idv) != nrow(nd)) stop("id must have one value per newdata row.", call. = FALSE)
  if (include_random && ncol(des$Z)) {
    known <- match(as.character(idv), rownames(object$random_effects)); good <- !is.na(known)
    if (any(!good) && !allow_new_levels) stop("newdata contains unseen subject levels.", call. = FALSE)
    if (any(good)) eta[good] <- eta[good] + rowSums(des$Z[good,,drop=FALSE] * object$random_effects[known[good],,drop=FALSE])
  }
  if (include_serial && object$autocor != "NONE") {
    tv <- time
    time_source <- .memwas_null_coalesce(object$.spec$time, object$.prep$time_spec)
    if (is.null(tv) && !is.null(time_source)) tv <- try(.memwas_resolve_data_value(time_source, nd, parent.frame(), "time"), silent = TRUE)
    if (inherits(tv, "try-error") || is.null(tv)) stop("Serial conditional prediction requires time values for newdata.", call. = FALSE)
    for (lev in unique(as.character(idv))) {
      ni <- which(as.character(idv) == lev); oi <- which(as.character(object$id) == lev)
      if (!length(oi)) next
      Koo <- .memwas_serial_covariance_matrix(object$serial_parameters, object$time[oi])
      Kno <- .memwas_serial_covariance_matrix(object$serial_parameters, as.numeric(tv[ni]), object$time[oi])
      eta[ni] <- eta[ni] + as.vector(Kno %*% .memwas_solve_pd(Koo, object$serial_effects[oi]))
    }
  }
  if (type == "link") eta else .memwas_family_mean(eta, object$.prep$family)
}

#' @export
AIC.MEMWAS_family_fit <- function(object, ..., k = 2) -2 * object$logLik + k * object$df

.BIC.MEMWAS_family_fit <- function(object, ...) -2 * object$logLik + log(object$nobs) * object$df

.simulate.MEMWAS_family_fit <- function(object, nsim = 1L, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  nsim <- as.integer(nsim)
  if (nsim < 1L) stop("nsim must be positive.", call. = FALSE)
  prep <- object$.prep; dec <- object$.decoded
  groups <- split(seq_len(object$nobs), object$id)
  ans <- matrix(NA_real_, object$nobs, nsim)
  for (ss in seq_len(nsim)) {
    eta <- object$linear_predictor_fixed
    for (ii in groups) {
      if (ncol(object$Z)) {
        Lb <- t(.memwas_chol_pd(dec$D)); b <- as.vector(Lb %*% stats::rnorm(ncol(object$Z)))
        eta[ii] <- eta[ii] + as.vector(object$Z[ii,,drop=FALSE] %*% b)
      }
      if (object$autocor != "NONE") {
        Ku <- .memwas_serial_covariance_matrix(dec$serial, object$time[ii])
        u <- as.vector(t(.memwas_chol_pd(Ku)) %*% stats::rnorm(length(ii)))
        eta[ii] <- eta[ii] + u
      }
    }
    mu <- .memwas_family_mean(eta, prep$family)
    ans[,ss] <- switch(object$family_name,
      gaussian = stats::rnorm(object$nobs, mu, dec$family$sigma),
      binomial = stats::rbinom(object$nobs, object$trials, mu),
      poisson = stats::rpois(object$nobs, mu),
      negative_binomial = stats::rnbinom(object$nobs, size = dec$family$theta, mu = mu),
      gamma = stats::rgamma(object$nobs, shape = dec$family$shape, scale = mu / dec$family$shape),
      exponential = stats::rexp(object$nobs, rate = 1 / mu))
  }
  ans <- as.data.frame(ans)
  names(ans) <- paste0("sim_", seq_len(nsim))
  ans
}
