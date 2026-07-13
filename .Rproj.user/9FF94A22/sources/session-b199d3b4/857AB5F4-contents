# Multiple predictor-specific serial processes
#
# ---------------
# For subject i, the linear predictor is
#
#   eta_i = X_i beta + Z_i b_i + sum_k diag(h_ik) u_ik + offset_i,
#
# where each u_ik has its own covariance structure and is independent of the
# other serial components and of b_i.  h_ik is a single design column, usually
# 1 for a residual process or the observed value of one predictor for a
# predictor-specific time-varying coefficient process.
# ---------------
#
# This file loads after helpers-family.R.  It leaves single-autocorrelation intact
# and extend engine only when `serial` or a multi-component `autocor` specification
# is supplied.
#
# These functions are intentionally not exported and do not add package dependencies.

#' Define one MEMWAS serial-process component
#'
#' @param structure Character scalar naming the covariance structure. Supported
#'   values are `"NONE"`, `"AR(1)"`, `"EXP"`/`"OU"`, `"AR(p)"`,
#'   `"ARMA(1,1)"`, `"CS"`, `"TOEP"`, and `"UN"`.
#' @param predictor Optional character scalar naming one numeric predictor. It
#'   is shorthand for `design = ~ 0 + predictor`.
#' @param design Optional one-sided formula that must create exactly one design
#'   column. Use `~ 1` for an outcome-level residual process and `~ 0 + x` for
#'   a predictor-specific process.
#' @param order Optional positive integer for AR(p) or TOEP structures.
#' @param name Optional descriptive component name. A name supplied by the
#'   enclosing `serial = list(...)` takes precedence.
#' @param control Optional component-specific controls. These override matching
#'   entries in the model-level `control` list for this component.
#'
#' @return An object of class `MEMWAS_serial_component` for use in
#'   `fit_MEMWAS(serial = ...)` or `set_MEMWAS(serial = ...)`.
#' @export
serial_component <- function(structure = "AR(1)", predictor = NULL,
                             design = NULL, order = NULL, name = NULL,
                             control = list()) {
  if (!is.character(structure) || length(structure) != 1L || is.na(structure) || !nzchar(structure)) {
    stop("`structure` must be one non-empty character value.", call. = FALSE)
  }
  if (!is.null(predictor) && !is.null(design)) {
    stop("Supply either `predictor` or `design`, not both.", call. = FALSE)
  }
  if (!is.null(predictor)) {
    if (!is.character(predictor) || length(predictor) != 1L || is.na(predictor) || !nzchar(predictor)) {
      stop("`predictor` must be one non-empty column name.", call. = FALSE)
    }
  }
  if (!is.null(design)) {
    if (is.character(design) && length(design) == 1L) {
      design <- if (grepl("~", design, fixed = TRUE)) {
        stats::as.formula(design, env = parent.frame())
      } else {
        stats::reformulate(design, intercept = FALSE, env = parent.frame())
      }
    }
    if (!inherits(design, "formula") || length(design) != 2L) {
      stop("`design` must be a one-sided formula such as `~ 1` or `~ 0 + x`.", call. = FALSE)
    }
  }
  if (!is.null(order)) {
    if (!is.numeric(order) || length(order) != 1L || !is.finite(order) || order < 1 || order != as.integer(order)) {
      stop("`order` must be one positive integer.", call. = FALSE)
    }
    order <- as.integer(order)
  }
  if (!is.null(name) && (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name))) {
    stop("`name` must be NULL or one non-empty character value.", call. = FALSE)
  }
  if (!is.list(control)) stop("`control` must be a list.", call. = FALSE)
  base::structure(
    list(structure = structure, predictor = predictor, design = design,
         order = order, name = name, control = control),
    class = "MEMWAS_serial_component"
  )
}

#' @noRd
#' @export
print.MEMWAS_serial_component <- function(x, ...) {
  target <- if (!is.null(x$predictor)) {
    paste0("predictor ", sQuote(x$predictor))
  } else if (!is.null(x$design)) {
    paste(deparse(x$design), collapse = " ")
  } else {
    "~ 1"
  }
  cat("MEMWAS serial component:", x$structure, "on", target, "\n")
  invisible(x)
}

.memwas_serial_component_fields <- c(
  "structure", "type", "autocor", "correlation", "design", "formula",
  "predictor", "term", "order", "p", "name", "label", "control"
)

.memwas_residual_component_names <- c(
  "residual", "residuals", "intercept", "(intercept)",
  "outcome", "response", "error", "errors", "epsilon"
)

.memwas_is_serial_component_spec <- function(x) {
  inherits(x, "MEMWAS_serial_component") ||
    (is.list(x) && length(intersect(names(x), .memwas_serial_component_fields)) > 0L)
}

.memwas_is_none_serial <- function(x) {
  if (is.null(x)) return(TRUE)
  if (is.character(x) && length(x) == 1L) {
    key <- toupper(gsub("[ _-]+", "", x))
    return(key %in% c("NONE", "NO", "NULL", "IID", "INDEPENDENT"))
  }
  FALSE
}

.memwas_is_multi_serial_autocor <- function(x) {
  if (is.null(x)) return(FALSE)
  if (inherits(x, "MEMWAS_serial_component")) return(TRUE)
  if (is.character(x)) {
    nm <- names(x)
    return(length(x) > 1L || (!is.null(nm) && any(nzchar(nm))))
  }
  if (!is.list(x)) return(FALSE)
  if (.memwas_is_serial_component_spec(x)) {
    return(!is.null(x$design) || !is.null(x$formula) || !is.null(x$predictor) || !is.null(x$term))
  }
  length(x) > 0L
}

.memwas_validate_serial_syntax <- function(serial) {
  if (is.null(serial)) return(invisible(TRUE))
  ok <- inherits(serial, "MEMWAS_serial_component") || is.character(serial) || is.list(serial)
  if (!ok) {
    stop("`serial` must be a serial_component(), a character specification, or a list of components.", call. = FALSE)
  }
  invisible(TRUE)
}

.memwas_collect_serial_components <- function(serial) {
  .memwas_validate_serial_syntax(serial)
  if (is.null(serial)) return(list())
  if (inherits(serial, "MEMWAS_serial_component") || .memwas_is_serial_component_spec(serial)) {
    return(list(serial))
  }
  if (is.character(serial)) {
    nm <- names(serial)
    if (length(serial) == 1L && (is.null(nm) || !nzchar(nm[1L]))) return(list(serial))
    if (is.null(nm) || any(!nzchar(nm))) {
      stop("A character vector with multiple serial structures must be fully named by component or predictor.", call. = FALSE)
    }
    return(as.list(serial))
  }
  if (!is.list(serial)) stop("Invalid `serial` specification.", call. = FALSE)
  if (length(serial) > 1L && (is.null(names(serial)) || any(!nzchar(names(serial))))) {
    stop("Every component in a multi-component `serial` list must have a unique name.", call. = FALSE)
  }
  serial
}

.memwas_normalize_serial_components <- function(serial, data, envir) {
  entries <- .memwas_collect_serial_components(serial)
  if (!length(entries)) return(list())
  entry_names <- names(entries)
  if (is.null(entry_names)) entry_names <- rep("", length(entries))
  defs <- vector("list", length(entries))
  used <- character()
  for (k in seq_along(entries)) {
    entry <- entries[[k]]
    outer_name <- entry_names[k]
    if (is.character(entry)) {
      if (length(entry) != 1L || is.na(entry)) stop("Each character serial component must name exactly one structure.", call. = FALSE)
      entry <- list(structure = entry)
    }
    if (!is.list(entry)) stop("Each serial component must be a list, character scalar, or serial_component().", call. = FALSE)
    structure_name <- .memwas_null_coalesce(entry$structure,
      .memwas_null_coalesce(entry$type,
        .memwas_null_coalesce(entry$autocor, entry$correlation)))
    if (is.null(structure_name)) stop("Each serial component requires `structure` (or `type`).", call. = FALSE)
    if (!is.character(structure_name) || length(structure_name) != 1L ||
        is.na(structure_name) || !nzchar(structure_name)) {
      stop("Each serial component structure must be one non-empty character value.", call. = FALSE)
    }
    predictor <- .memwas_null_coalesce(entry$predictor, entry$term)
    design <- .memwas_null_coalesce(entry$design, entry$formula)
    if (!is.null(predictor) && !is.null(design)) {
      stop("Serial component ", sQuote(outer_name), " supplies both predictor and design.", call. = FALSE)
    }
    inner_name <- .memwas_null_coalesce(entry$name, entry$label)
    outer_ok <- length(outer_name) == 1L && !is.na(outer_name) && nzchar(outer_name)
    component_name <- if (outer_ok) outer_name else if (!is.null(inner_name)) as.character(inner_name)[1L] else ""
    if (length(component_name) != 1L || is.na(component_name)) {
      stop("Each serial component name must be one non-missing character value.", call. = FALSE)
    }
    component_key <- tolower(trimws(component_name))
    reserved_residual_name <- nzchar(component_key) &&
      component_key %in% .memwas_residual_component_names
    if (is.null(predictor) && is.null(design) && nzchar(component_name) &&
        !reserved_residual_name && component_name %in% names(data)) {
      predictor <- component_name
    }
    if (is.null(predictor) && is.null(design)) {
      design <- stats::as.formula("~ 1", env = envir)
    }
    if (!is.null(predictor)) {
      if (!is.character(predictor) || length(predictor) != 1L || is.na(predictor) || !nzchar(predictor)) {
        stop("A serial component predictor must be one non-empty column name.", call. = FALSE)
      }
      if (!predictor %in% names(data)) {
        stop("Serial component predictor ", sQuote(predictor), " was not found in data.", call. = FALSE)
      }
      design <- stats::reformulate(predictor, intercept = FALSE, env = envir)
    }
    if (is.character(design) && length(design) == 1L) {
      design <- if (grepl("~", design, fixed = TRUE)) {
        stats::as.formula(design, env = envir)
      } else {
        stats::reformulate(design, intercept = FALSE, env = envir)
      }
    }
    if (!inherits(design, "formula") || length(design) != 2L) {
      stop("Each serial design must be one-sided, for example `~ 1` or `~ 0 + x`.", call. = FALSE)
    }
    if (!nzchar(component_name)) {
      component_name <- if (!is.null(predictor)) predictor else if (length(entries) == 1L) "residual" else paste0("serial", k)
    }
    component_name <- as.character(component_name)[1L]
    if (component_name %in% used) stop("Serial component names must be unique: ", component_name, call. = FALSE)
    used <- c(used, component_name)
    component_control <- entry$control
    if (is.null(component_control)) component_control <- list()
    if (!is.list(component_control)) stop("Component-specific control must be a list.", call. = FALSE)
    defs[[k]] <- list(
      name = component_name,
      structure = structure_name,
      order = .memwas_null_coalesce(entry$order, entry$p),
      predictor = predictor,
      design_formula = design,
      control = component_control
    )
  }
  names(defs) <- vapply(defs, `[[`, character(1), "name")
  defs
}

.memwas_serial_design_component <- function(def, data) {
  mf <- stats::model.frame(def$design_formula, data = data,
                           na.action = stats::na.pass,
                           drop.unused.levels = TRUE)
  tt <- stats::terms(mf)
  mm <- stats::model.matrix(tt, mf)
  if (nrow(mm) != nrow(data)) {
    stop("Serial design ", sQuote(def$name), " did not produce one row per data row.", call. = FALSE)
  }
  if (ncol(mm) != 1L) {
    stop(
      "Serial component ", sQuote(def$name), " produced ", ncol(mm),
      " design columns (", paste(colnames(mm), collapse = ", "), "). ",
      "Each component must produce exactly one column. Use `~ 1` for a residual process ",
      "or `~ 0 + predictor` for a predictor-specific process.",
      call. = FALSE
    )
  }
  list(
    values = as.numeric(mm[, 1L]),
    terms = tt,
    xlevels = stats::.getXlevels(tt, mf),
    contrasts = attr(mm, "contrasts"),
    design_column = colnames(mm)[1L]
  )
}

.memwas_subset_prepared_data <- function(prep, keep2) {
  row_vector <- c("y", "trials", "freq", "offset", "id", "time")
  for (nm in row_vector) prep[[nm]] <- prep[[nm]][keep2]
  prep$X <- prep$X[keep2, , drop = FALSE]
  prep$Z <- prep$Z[keep2, , drop = FALSE]
  prep$keep <- prep$keep[keep2]
  prep$id <- droplevels(factor(prep$id))
  prep
}

.memwas_serial_component_signature <- function(comp) {
  paste(comp$spec$type, .memwas_null_coalesce(comp$spec$order, ""), sep = ":")
}

.memwas_check_serial_identifiability <- function(H, components, control, Z = NULL,
                                                family_name = NULL, time = NULL) {
  if (!ncol(H)) return(invisible(TRUE))
  tol <- as.numeric(.memwas_null_coalesce(control$serial_design_tolerance, 1e-9))
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("control$serial_design_tolerance must be one positive finite number.", call. = FALSE)
  }
  for (j in seq_len(ncol(H))) {
    hj <- H[, j]
    if (!all(is.finite(hj))) {
      stop("Serial component ", sQuote(colnames(H)[j]), " contains non-finite design values after missing-value handling.", call. = FALSE)
    }
    if (sqrt(sum(hj^2)) <= tol) {
      stop("Serial component ", sQuote(colnames(H)[j]), " has an all-zero design and is not identifiable.", call. = FALSE)
    }
  }

  # A CS process on a design already represented by a random coefficient has
  # an exactly confounded rank-one covariance component.  An UN process can
  # also absorb that random-effect covariance.  Reject these decompositions
  # instead of relying on an optimizer to choose an arbitrary split.
  if (!is.null(Z) && ncol(Z)) {
    qrz <- qr(Z, tol = tol)
    for (j in seq_len(ncol(H))) {
      hj <- H[, j]
      resid <- qr.resid(qrz, hj)
      relative <- sqrt(sum(resid^2)) / max(sqrt(sum(hj^2)), tol)
      st <- components[[j]]$spec$type
      if (is.finite(relative) && relative <= tol && st %in% c("CS", "UN")) {
        stop(
          "Serial component ", sQuote(colnames(H)[j]), " uses ", st,
          " and its design lies in the random-effects design space. ",
          "The corresponding covariance component is not separately identifiable; ",
          "remove the overlapping random coefficient or choose a structured serial covariance.",
          call. = FALSE
        )
      }
      if (is.finite(relative) && relative <= tol &&
          isTRUE(.memwas_null_coalesce(control$warn_serial_random_overlap, FALSE))) {
        warning(
          "Serial component ", sQuote(colnames(H)[j]),
          " has a design in the random-effects design space; inspect covariance profiles for weak identification.",
          call. = FALSE
        )
      }
    }
  }

  if (ncol(H) >= 2L) {
    for (a in seq_len(ncol(H) - 1L)) for (b in (a + 1L):ncol(H)) {
      x <- H[, a]; y <- H[, b]
      den <- sum(x * x)
      if (!is.finite(den) || den <= tol) next
      scale <- sum(x * y) / den
      relative <- sqrt(sum((y - scale * x)^2)) / max(sqrt(sum(y^2)), tol)
      if (is.finite(relative) && relative <= tol) {
        same <- identical(.memwas_serial_component_signature(components[[a]]),
                          .memwas_serial_component_signature(components[[b]]))
        has_un <- any(vapply(components[c(a, b)], function(z) identical(z$spec$type, "UN"), logical(1)))
        msg <- paste0(
          "Serial designs ", sQuote(colnames(H)[a]), " and ", sQuote(colnames(H)[b]),
          " are proportional"
        )
        if (same) {
          stop(msg, " and use the same covariance family; their variance components are not separately identifiable.", call. = FALSE)
        }
        if (has_un) {
          stop(msg, " and one component is unstructured; the UN component subsumes the other covariance on this design.", call. = FALSE)
        }
        warning(msg, "; different covariance shapes may be weakly identifiable. Inspect profile likelihoods and convergence diagnostics.", call. = FALSE)
      }
    }
  }

  # In Gaussian models the observation nugget is a free multiple of I.  The
  # diagonal parameters of one or more UN components can absorb that nugget
  # whenever, at every visit time, the all-ones vector lies in the column space
  # of the squared UN design columns.  This includes a residual UN process and
  # predictor designs whose magnitude is constant within visit time.  Detect
  # the full multi-component span rather than checking components one by one.
  if (identical(family_name, "gaussian")) {
    un_idx <- which(vapply(components, function(z) identical(z$spec$type, "UN"), logical(1)))
    if (length(un_idx)) {
      if (is.null(time) || length(time) != nrow(H) || any(!is.finite(time))) {
        stop("Finite observation times are required to check Gaussian UN/nugget identifiability.", call. = FALSE)
      }
      strata <- split(seq_len(nrow(H)), time)
      nugget_in_span <- vapply(strata, function(ii) {
        G <- H[ii, un_idx, drop = FALSE]^2
        keep <- vapply(seq_len(ncol(G)), function(j) sqrt(sum(G[, j]^2)) > tol, logical(1))
        G <- G[, keep, drop = FALSE]
        if (!ncol(G)) return(FALSE)
        target <- rep(1, length(ii))
        residual <- try(qr.resid(qr(G, tol = tol), target), silent = TRUE)
        if (inherits(residual, "try-error")) return(FALSE)
        sqrt(sum(residual^2)) / max(sqrt(sum(target^2)), tol) <= tol
      }, logical(1))
      if (length(nugget_in_span) && all(nugget_in_span)) {
        component_labels <- paste(sQuote(names(components)[un_idx]), collapse = ", ")
        msg <- paste0(
          "The Gaussian observation variance is not separately identifiable from the diagonal covariance of UN serial component(s) ",
          component_labels, ". At every observed visit time, their squared design columns span a constant nugget. ",
          "Use a structured covariance, add within-visit design variation, remove the Gaussian nugget through a dedicated model reformulation, or set ",
          "control$allow_unstructured_nugget_confounding = TRUE only for sensitivity analysis."
        )
        if (!isTRUE(control$allow_unstructured_nugget_confounding)) stop(msg, call. = FALSE)
        warning(msg, call. = FALSE)
      }
    }
  }
  invisible(TRUE)
}

.memwas_prepare_multi_serial <- function(spec) {
  spec0 <- spec
  spec0$autocor <- "NONE"
  prep <- .memwas_prepare_single_serial(spec0)
  finite_blocks <- c(prep$y, prep$trials, prep$freq, prep$offset, prep$time,
                     as.numeric(prep$X), as.numeric(prep$Z))
  if (any(!is.finite(finite_blocks))) {
    stop("The multiple-serial engine requires finite response, design, offset, weight, and time values.", call. = FALSE)
  }
  defs <- .memwas_normalize_serial_components(spec$serial, prep$data, spec$environment)
  max_components_raw <- .memwas_null_coalesce(spec$control$max_serial_components, 8L)
  if (!is.numeric(max_components_raw) || length(max_components_raw) != 1L ||
      !is.finite(max_components_raw) || max_components_raw < 1 ||
      max_components_raw != as.integer(max_components_raw)) {
    stop("control$max_serial_components must be one positive integer.", call. = FALSE)
  }
  max_components <- as.integer(max_components_raw)
  if (length(defs) > max_components) {
    stop("The model requests ", length(defs), " serial components; increase control$max_serial_components (currently ", max_components, ") only after checking identifiability.", call. = FALSE)
  }
  if (!length(defs)) {
    stop("`serial` must contain at least one component; use serial = NULL for a model without serial processes.", call. = FALSE)
  }

  # Normalize structures once so NONE components can be discarded before their
  # design is allowed to affect complete-case filtering.
  first_specs <- lapply(defs, function(def) {
    cc <- utils::modifyList(spec$control, def$control)
    .memwas_serial_covariance_spec(list(type = def$structure, order = def$order), prep$id, prep$time, cc)
  })
  active <- !vapply(first_specs, function(x) identical(x$type, "NONE"), logical(1))
  if (any(!active)) {
    warning("Ignoring serial component(s) with structure NONE: ", paste(names(defs)[!active], collapse = ", "), call. = FALSE)
    defs <- defs[active]
  }
  if (!length(defs)) {
    stop("All supplied serial components use structure NONE; use serial = NULL instead.", call. = FALSE)
  }

  design_meta <- lapply(defs, .memwas_serial_design_component, data = prep$data)
  H <- do.call(cbind, lapply(design_meta, function(z) z$values[prep$keep]))
  if (is.null(dim(H))) H <- matrix(H, ncol = 1L)
  colnames(H) <- names(defs)
  okH <- stats::complete.cases(H)
  if (!all(okH)) {
    prep <- .memwas_subset_prepared_data(prep, okH)
    H <- H[okH, , drop = FALSE]
  }
  if (!nrow(H)) stop("No complete observations remain after evaluating serial-component designs.", call. = FALSE)

  components <- vector("list", length(defs)); names(components) <- names(defs)
  for (k in seq_along(defs)) {
    def <- defs[[k]]; dm <- design_meta[[k]]
    cc <- utils::modifyList(spec$control, def$control)
    ss <- .memwas_serial_covariance_spec(list(type = def$structure, order = def$order), prep$id, prep$time, cc)
    components[[k]] <- c(def, dm[c("terms", "xlevels", "contrasts", "design_column")],
                            list(spec = ss, design = H[, k], control = cc))
  }
  .memwas_check_serial_identifiability(
    H, components, spec$control, Z = prep$Z,
    family_name = prep$family$name, time = prep$time
  )
  prep$serial_components <- components
  prep$serial_H <- H
  prep$serial <- if (length(components) == 1L) components[[1L]]$spec else list(type = "MULTIPLE")
  prep
}

.memwas_serial_parameter_count <- function(ss) {
  type <- ss$type
  if (identical(type, "NONE")) return(0L)
  if (type %in% c("AR1", "EXP", "CS")) return(2L)
  if (identical(type, "ARMA11")) return(3L)
  if (type %in% c("ARP", "TOEP")) return(1L + ss$order)
  if (identical(type, "UN")) {
    m <- length(ss$time_levels)
    return(as.integer(m * (m + 1L) / 2L))
  }
  stop("Unknown serial structure in layout: ", type, call. = FALSE)
}

.memwas_serial_parameter_names <- function(component) {
  ss <- component$spec
  type <- ss$type
  suffix <- if (identical(type, "NONE")) {
    character()
  } else if (identical(type, "AR1") || identical(type, "CS")) {
    c("log_sd", "correlation_raw")
  } else if (identical(type, "EXP")) {
    c("log_sd", "log_range")
  } else if (identical(type, "ARMA11")) {
    c("log_sd", "phi_raw", "theta_raw")
  } else if (type %in% c("ARP", "TOEP")) {
    c("log_sd", paste0("pacf_raw_", seq_len(ss$order)))
  } else if (identical(type, "UN")) {
    m <- length(ss$time_levels)
    out <- character(m * (m + 1L) / 2L)
    z <- 0L
    for (i in seq_len(m)) {
      for (j in seq_len(i)) {
        z <- z + 1L
        out[z] <- paste0("chol_", i, "_", j)
      }
    }
    out
  } else {
    stop("Unknown serial structure when naming parameters: ", type, call. = FALSE)
  }
  paste0("serial_", .memwas_safe_name(component$name), "_", suffix)
}

.memwas_multi_serial_layout <- function(prep, spec) {
  p <- ncol(prep$X); q <- ncol(prep$Z)
  random_cov <- .memwas_normalize_random_covariance(spec$random_cov)
  rn <- if (!q) 0L else if (random_cov == "diagonal") q else q * (q + 1L) / 2L
  serial_sizes <- vapply(prep$serial_components, function(z) .memwas_serial_parameter_count(z$spec), integer(1))
  total_serial <- sum(serial_sizes)
  max_serial_parameters_raw <- .memwas_null_coalesce(spec$control$max_total_serial_parameters, 100L)
  if (!is.numeric(max_serial_parameters_raw) || length(max_serial_parameters_raw) != 1L ||
      !is.finite(max_serial_parameters_raw) || max_serial_parameters_raw < 1 ||
      max_serial_parameters_raw != as.integer(max_serial_parameters_raw)) {
    stop("control$max_total_serial_parameters must be one positive integer.", call. = FALSE)
  }
  max_serial_parameters <- as.integer(max_serial_parameters_raw)
  if (total_serial > max_serial_parameters) {
    stop("The requested serial components require ", total_serial,
         " covariance parameters; increase control$max_total_serial_parameters (currently ",
         max_serial_parameters, ") only with adequate replication.", call. = FALSE)
  }
  fn <- if (prep$family$name == "gaussian") 1L else
    if (prep$family$name == "gamma" && is.null(prep$family$fixed_shape)) 1L else
      if (prep$family$name == "negative_binomial" && is.null(prep$family$fixed_theta)) 1L else 0L
  idx <- list(beta = seq_len(p))
  at <- p
  idx$random <- if (rn) at + seq_len(rn) else integer(); at <- at + rn
  idx$serial <- vector("list", length(serial_sizes)); names(idx$serial) <- names(serial_sizes)
  serial_names <- character()
  if (length(serial_sizes)) {
    for (k in seq_along(serial_sizes)) {
      nk <- serial_sizes[k]
      idx$serial[[k]] <- if (nk) at + seq_len(nk) else integer()
      at <- at + nk
      serial_names <- c(serial_names, .memwas_serial_parameter_names(prep$serial_components[[k]]))
    }
  }
  idx$serial_flat <- unlist(idx$serial, use.names = FALSE)
  idx$family <- if (fn) at + seq_len(fn) else integer(); at <- at + fn
  nms <- colnames(prep$X)
  if (rn) nms <- c(nms, paste0("cov_random_", seq_len(rn)))
  nms <- c(nms, serial_names)
  if (fn) nms <- c(nms, switch(prep$family$name,
    gaussian = "log_sigma", gamma = "log_shape", negative_binomial = "log_size"))
  list(
    p = p, q = q, random_cov = random_cov, random_n = rn,
    serial_n = total_serial, serial_sizes = serial_sizes,
    serial_component_n = length(serial_sizes), serial_names = names(serial_sizes),
    family_n = fn, idx = idx, npar = at, names = nms
  )
}

.memwas_decode_serial_component <- function(sp, component) {
  ss <- component$spec
  serial <- list(
    name = component$name,
    type = ss$type,
    time_levels = ss$time_levels,
    scale = ss$scale,
    order = ss$order,
    design_column = component$design_column,
    predictor = component$predictor
  )
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
    serial$sd <- exp(sp[1L])
    serial$rho <- lower + (0.999 - lower) * .memwas_inverse_logit(sp[2L])
  } else if (ss$type == "ARMA11") {
    serial$sd <- exp(sp[1L]); serial$phi <- tanh(sp[2L]); serial$theta <- tanh(sp[3L])
  } else if (ss$type %in% c("ARP", "TOEP")) {
    serial$sd <- exp(sp[1L]); serial$pacf <- tanh(sp[-1L]); serial$ar <- .memwas_family_pacf_to_ar(serial$pacf)
  } else if (ss$type == "UN") {
    m <- length(ss$time_levels)
    serial$Sigma_full <- .memwas_decode_cholesky_parameters(sp, m, FALSE)
  }
  serial
}

.memwas_decode_multi_serial_parameters <- function(par, prep, spec, layout) {
  beta <- par[layout$idx$beta]
  D <- .memwas_decode_cholesky_parameters(par[layout$idx$random], layout$q, layout$random_cov == "diagonal")
  serial_components <- vector("list", length(prep$serial_components))
  names(serial_components) <- names(prep$serial_components)
  for (k in seq_along(serial_components)) {
    serial_components[[k]] <- .memwas_decode_serial_component(par[layout$idx$serial[[k]]], prep$serial_components[[k]])
  }
  fampar <- list()
  if (prep$family$name == "gaussian") fampar$sigma <- exp(par[layout$idx$family])
  if (prep$family$name == "gamma") fampar$shape <- if (is.null(prep$family$fixed_shape)) exp(par[layout$idx$family]) else prep$family$fixed_shape
  if (prep$family$name == "negative_binomial") fampar$theta <- if (is.null(prep$family$fixed_theta)) exp(par[layout$idx$family]) else prep$family$fixed_theta
  if (prep$family$name == "exponential") fampar$shape <- 1
  list(
    beta = beta, D = D,
    serial_components = serial_components,
    serial = if (length(serial_components) == 1L) serial_components[[1L]] else serial_components,
    family = fampar
  )
}

.memwas_multi_serial_subject_components <- function(prep, dec, idx) {
  Zi <- prep$Z[idx, , drop = FALSE]
  ti <- prep$time[idx]
  mats <- list(); covs <- list(); slices <- list()
  cursor <- 0L
  if (ncol(Zi)) {
    mats[[length(mats) + 1L]] <- Zi
    covs[[length(covs) + 1L]] <- dec$D
    cursor <- ncol(Zi)
  }
  for (nm in names(dec$serial_components)) {
    s <- dec$serial_components[[nm]]
    h <- prep$serial_H[idx, nm]
    m <- length(idx)
    Ak <- diag(as.numeric(h), nrow = m, ncol = m)
    Kk <- .memwas_serial_covariance_matrix(s, ti)
    mats[[length(mats) + 1L]] <- Ak
    covs[[length(covs) + 1L]] <- Kk
    slices[[nm]] <- cursor + seq_len(m)
    cursor <- cursor + m
  }
  A <- if (length(mats)) do.call(cbind, mats) else matrix(numeric(), length(idx), 0L)
  C <- matrix(numeric(), 0L, 0L)
  for (cc in covs) C <- .memwas_block_diag(C, cc)
  list(A = A, C = C, q = ncol(Zi), serial_slices = slices,
       r = length(idx) * length(dec$serial_components))
}

.memwas_multi_serial_initial_values <- function(prep, spec, layout) {
  fam <- prep$family
  beta <- rep(0, layout$p); names(beta) <- colnames(prep$X)
  fit0 <- try({
    if (fam$name == "binomial") {
      yy <- cbind(prep$y, prep$trials - prep$y)
      stats::glm.fit(prep$X, yy, family = stats::binomial(fam$link),
                     offset = prep$offset, weights = prep$freq)
    } else {
      gf <- if (fam$name == "gaussian") stats::gaussian(fam$link) else
        if (fam$name == "poisson") stats::poisson(fam$link) else
          if (fam$name %in% c("gamma", "exponential")) stats::Gamma(fam$link) else stats::poisson("log")
      stats::glm.fit(prep$X, prep$y, family = gf,
                     offset = prep$offset, weights = prep$freq)
    }
  }, silent = TRUE)
  if (!inherits(fit0, "try-error") && length(fit0$coefficients) == layout$p) {
    beta <- fit0$coefficients; beta[!is.finite(beta)] <- 0
  }
  par <- numeric(layout$npar); par[layout$idx$beta] <- beta
  if (layout$random_n) {
    if (layout$random_cov == "diagonal") {
      par[layout$idx$random] <- log(0.5)
    } else {
      z <- 0L; vals <- numeric(layout$random_n)
      for (i in seq_len(layout$q)) for (j in seq_len(i)) {
        z <- z + 1L; vals[z] <- if (i == j) log(0.5) else 0
      }
      par[layout$idx$random] <- vals
    }
  }
  if (layout$serial_n) {
    for (k in seq_along(prep$serial_components)) {
      ss <- prep$serial_components[[k]]$spec
      ii <- layout$idx$serial[[k]]
      vals <- numeric(length(ii))
      if (ss$type == "UN") {
        m <- length(ss$time_levels); z <- 0L
        for (i in seq_len(m)) for (j in seq_len(i)) {
          z <- z + 1L; vals[z] <- if (i == j) log(0.35) else 0
        }
      } else if (length(vals)) {
        vals[1L] <- log(0.35)
      }
      par[ii] <- vals
    }
  }
  if (layout$family_n) {
    if (fam$name == "gaussian") {
      sigma0 <- stats::sd(prep$y)
      if (!is.finite(sigma0) || sigma0 <= 0) sigma0 <- max(abs(prep$y - mean(prep$y)), 0.1)
      if (!is.finite(sigma0) || sigma0 <= 0) sigma0 <- 0.1
      par[layout$idx$family] <- log(sigma0)
    }
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
      sstart <- .memwas_null_coalesce(st$serial, st$serial_parameters)
      if (!is.null(sstart)) {
        if (is.numeric(sstart) && length(prep$serial_components) == 1L) sstart <- list(sstart)
        if (!is.list(sstart)) stop("start$serial must be a named list of raw parameter vectors.", call. = FALSE)
        if (is.null(names(sstart)) && length(sstart) == length(prep$serial_components)) names(sstart) <- names(prep$serial_components)
        if (is.null(names(sstart)) || any(!nzchar(names(sstart)))) {
          stop("start$serial must be named by serial component.", call. = FALSE)
        }
        unknown_start <- setdiff(names(sstart), names(layout$idx$serial))
        if (length(unknown_start)) {
          stop("Unknown start$serial component(s): ", paste(unknown_start, collapse = ", "), call. = FALSE)
        }
        for (nm in names(sstart)) {
          par[layout$idx$serial[[nm]]] <- .memwas_validate_raw_start_vector(
            sstart[[nm]],
            length(layout$idx$serial[[nm]]),
            paste0("start$serial[[", sQuote(nm), "]]")
          )
        }
      }
    } else {
      par <- .memwas_validate_raw_start_vector(st, length(par), "start")
    }
  }
  names(par) <- layout$names
  par
}

.memwas_normalize_multi_serial_approximation <- function(x) {
  if (!.memwas_is_scalar_character(x) || !nzchar(trimws(x))) {
    stop("For multiple serial components, approximation must be 'laplace' or low-dimensional 'quadrature'.", call. = FALSE)
  }
  key <- tolower(gsub("[ _-]+", "", x))
  if (key %in% c("laplace", "lap", "skewcorrectedlaplace", "saddlepoint", "variationalinference", "pql")) {
    if (!key %in% c("laplace", "lap")) {
      warning("Multiple serial components currently use the full joint Laplace objective; approximation '", x, "' was mapped to Laplace.", call. = FALSE)
    }
    return("laplace")
  }
  if (key %in% c("quadrature", "aghq", "adaptivegausshermite", "adaptivegausshermitequadrature", "adaptivegaussianquadrature", "gausshermite", "ghq")) return("quadrature")
  stop("For multiple serial components, approximation must be 'laplace' or low-dimensional 'quadrature'.", call. = FALSE)
}

.memwas_fit_multi_serial_engine <- function(spec) {
  if (is.null(spec$control)) spec$control <- list()
  if (!is.list(spec$control)) stop("control must be a list.", call. = FALSE)
  .memwas_check_scalar_nonnegative(spec$L1_penalty, "L1_penalty")
  .memwas_check_scalar_nonnegative(spec$L2_penalty, "L2_penalty")
  spec$quadrature_points <- .memwas_validate_quadrature_points(spec$quadrature_points)
  cross_value <- .memwas_null_coalesce(
    spec$control$serial_cross_covariance, "independent"
  )
  if (!.memwas_is_scalar_character(cross_value) || !nzchar(trimws(cross_value))) {
    stop("control$serial_cross_covariance must be one non-empty character value.", call. = FALSE)
  }
  cross <- tolower(gsub("[ _-]+", "", cross_value))
  if (!cross %in% c("independent", "none", "zero")) {
    stop("Only independent serial components are currently identified. Set control$serial_cross_covariance = 'independent'.", call. = FALSE)
  }
  prep <- .memwas_prepare_multi_serial(spec)
  layout <- .memwas_multi_serial_layout(prep, spec)
  approx <- .memwas_normalize_multi_serial_approximation(spec$approximation)
  method <- .memwas_normalize_likelihood_method(.memwas_null_coalesce(spec$method, "ML"))
  if (method == "REML") {
    stop("REML is not available for predictor-specific multiple serial components. Use method = 'ML'; single-structure Gaussian models may still use the legacy REML engine.", call. = FALSE)
  }
  groups <- split(seq_along(prep$y), prep$id)
  kserial <- length(prep$serial_components)
  max_latent <- max(vapply(groups, function(ii) layout$q + kserial * length(ii), integer(1)))
  if (approx == "quadrature" && max_latent > 2L) {
    stop("Quadrature would require more than two latent dimensions per subject. Use approximation='laplace'; the requested multi-component model was not simplified silently.", call. = FALSE)
  }
  warn_dim_raw <- .memwas_null_coalesce(spec$control$warn_laplace_latent_dim, 80L)
  if (!is.numeric(warn_dim_raw) || length(warn_dim_raw) != 1L ||
      !is.finite(warn_dim_raw) || warn_dim_raw < 1 ||
      warn_dim_raw != as.integer(warn_dim_raw)) {
    stop("control$warn_laplace_latent_dim must be one positive integer.", call. = FALSE)
  }
  warn_dim <- as.integer(warn_dim_raw)
  if (approx == "laplace" && max_latent > warn_dim) {
    warning("The largest subject-specific latent vector has dimension ", max_latent,
            ". Multiple serial components can be computationally expensive; check convergence and sensitivity.", call. = FALSE)
  }
  init <- .memwas_multi_serial_initial_values(prep, spec, layout)
  cache <- new.env(parent = emptyenv())
  cache$modes <- vector("list", length(groups)); names(cache$modes) <- names(groups)
  penalty_idx <- setdiff(seq_len(layout$p), which(colnames(prep$X) == "(Intercept)"))
  if (isTRUE(spec$control$penalize_intercept)) penalty_idx <- seq_len(layout$p)

  evaluate <- function(par, details = FALSE, penalized = TRUE) {
    dec <- try(.memwas_decode_multi_serial_parameters(par, prep, spec, layout), silent = TRUE)
    if (inherits(dec, "try-error")) return(if (details) NULL else 1e50)
    total <- 0; out <- vector("list", length(groups)); names(out) <- names(groups)
    for (g in seq_along(groups)) {
      ii <- groups[[g]]
      sc <- try(.memwas_multi_serial_subject_components(prep, dec, ii), silent = TRUE)
      if (inherits(sc, "try-error")) return(if (details) NULL else 1e50)
      base_eta <- as.vector(prep$X[ii, , drop = FALSE] %*% dec$beta) + prep$offset[ii]
      im <- try(.memwas_optimize_latent_mode(
        base_eta, sc$A, sc$C, prep$y[ii], prep$trials[ii], prep$freq[ii],
        prep$family, dec$family, cache$modes[[g]], spec$control
      ), silent = TRUE)
      if (inherits(im, "try-error")) return(if (details) NULL else 1e50)
      if (approx == "laplace") {
        val <- im$joint + if (ncol(sc$A)) {
          0.5 * .memwas_logdet_pd(im$H) - 0.5 * ncol(sc$A) * log(2 * pi)
        } else 0
      } else {
        val <- try(.memwas_quadrature_negative_loglikelihood(
          base_eta, sc$A, sc$C, prep$y[ii], prep$trials[ii], prep$freq[ii],
          prep$family, dec$family, spec$quadrature_points,
          mode_info = im, control = spec$control
        ), silent = TRUE)
        if (inherits(val, "try-error") || !is.finite(val)) return(if (details) NULL else 1e50)
      }
      cache$modes[[g]] <- im$mode
      if (!is.finite(val)) return(if (details) NULL else 1e50)
      total <- total + val
      out[[g]] <- c(im, list(idx = ii, components = sc))
    }
    penalty <- spec$L1_penalty * sum(abs(dec$beta[penalty_idx])) +
      0.5 * spec$L2_penalty * sum(dec$beta[penalty_idx]^2)
    if (details) {
      list(nll = total, penalty = penalty, objective = total + penalty,
           decoded = dec, subjects = out)
    } else {
      total + if (penalized) penalty else 0
    }
  }

  maxit_raw <- .memwas_null_coalesce(spec$control$outer_maxit, 150L)
  if (!is.numeric(maxit_raw) || length(maxit_raw) != 1L ||
      !is.finite(maxit_raw) || maxit_raw < 1 || maxit_raw != as.integer(maxit_raw)) {
    stop("control$outer_maxit must be one positive integer.", call. = FALSE)
  }
  maxit <- as.integer(maxit_raw)
  evalmax_raw <- .memwas_null_coalesce(spec$control$outer_eval_max, max(500L, 5L * maxit))
  if (!is.numeric(evalmax_raw) || length(evalmax_raw) != 1L ||
      !is.finite(evalmax_raw) || evalmax_raw < 1 ||
      evalmax_raw != as.integer(evalmax_raw)) {
    stop("control$outer_eval_max must be one positive integer.", call. = FALSE)
  }
  evalmax <- as.integer(evalmax_raw)
  if (evalmax < maxit) {
    stop("control$outer_eval_max must be at least control$outer_maxit.", call. = FALSE)
  }
  reltol <- .memwas_null_coalesce(spec$control$outer_tol, 1e-7)
  if (!is.numeric(reltol) || length(reltol) != 1L ||
      !is.finite(reltol) || reltol <= 0) {
    stop("control$outer_tol must be one positive finite number.", call. = FALSE)
  }
  opt <- try(stats::nlminb(
    init, objective = evaluate,
    control = list(iter.max = maxit, eval.max = evalmax,
                   rel.tol = reltol, x.tol = reltol)
  ), silent = TRUE)
  if (inherits(opt, "try-error")) {
    oo <- stats::optim(init, evaluate, method = "BFGS",
                       control = list(maxit = maxit, reltol = reltol))
    opt <- list(par = oo$par, objective = oo$value,
                convergence = oo$convergence,
                message = .memwas_null_coalesce(oo$message, ""),
                iterations = oo$counts)
  }
  details <- evaluate(opt$par, details = TRUE)
  if (is.null(details)) stop("The final multiple-serial likelihood evaluation failed.", call. = FALSE)
  dec <- details$decoded
  eta_fixed <- as.vector(prep$X %*% dec$beta) + prep$offset
  eta_cond <- eta_fixed
  random_modes <- matrix(NA_real_, nlevels(prep$id), layout$q,
                         dimnames = list(levels(prep$id), colnames(prep$Z)))
  serial_modes <- matrix(0, length(prep$y), kserial,
                         dimnames = list(NULL, names(prep$serial_components)))
  serial_contributions <- serial_modes
  for (g in seq_along(details$subjects)) {
    ss <- details$subjects[[g]]; ii <- ss$idx; a <- ss$mode
    if (layout$q) random_modes[g, ] <- a[seq_len(layout$q)]
    if (kserial) {
      for (nm in names(ss$components$serial_slices)) {
        sl <- ss$components$serial_slices[[nm]]
        serial_modes[ii, nm] <- a[sl]
        serial_contributions[ii, nm] <- prep$serial_H[ii, nm] * a[sl]
      }
    }
    if (length(a)) eta_cond[ii] <- eta_fixed[ii] + as.vector(ss$components$A %*% a)
  }
  serial_effects <- if (kserial) rowSums(serial_contributions) else rep(0, length(prep$y))
  mu_fixed <- .memwas_family_mean(eta_fixed, prep$family)
  mu_cond <- .memwas_family_mean(eta_cond, prep$family)

  V <- matrix(NA_real_, layout$npar, layout$npar,
              dimnames = list(layout$names, layout$names))
  compute_vcov <- isTRUE(.memwas_null_coalesce(
    spec$control$compute_vcov,
    layout$npar <= 30L && spec$L1_penalty == 0
  ))
  if (compute_vcov) {
    hh <- try(stats::optimHess(opt$par, evaluate), silent = TRUE)
    if (!inherits(hh, "try-error") && all(is.finite(hh))) {
      vv <- try(.memwas_solve_pd(hh), silent = TRUE)
      if (!inherits(vv, "try-error")) V <- vv
    }
  }
  loglik <- -details$nll; df <- layout$npar
  aic <- -2 * loglik + 2 * df
  bic <- -2 * loglik + log(length(prep$y)) * df
  structures <- vapply(dec$serial_components, `[[`, character(1), "type")
  autocor_label <- if (!length(structures)) "NONE" else if (length(structures) == 1L) structures[[1L]] else "MULTIPLE"
  serial_parameters <- if (length(dec$serial_components) == 1L) dec$serial_components[[1L]] else dec$serial_components
  beta_named <- stats::setNames(dec$beta, colnames(prep$X))
  beta_vcov <- V[layout$idx$beta, layout$idx$beta, drop = FALSE]
  se <- sqrt(pmax(diag(beta_vcov), 0)); z <- beta_named / se
  coefficient_table <- data.frame(
    term = names(beta_named), estimate = as.numeric(beta_named),
    std_error = as.numeric(se), statistic = as.numeric(z),
    p_value = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  if (all(is.na(beta_vcov))) coefficient_table[, c("std_error", "statistic", "p_value")] <- NA_real_

  fit <- structure(list(
    call = spec$call, formula = spec$formula, random = spec$random,
    serial = spec$serial,
    family = prep$family$family_object,
    family_name = prep$family$name, link = prep$family$link,
    method = method, approximation = approx,
    autocor = autocor_label,
    serial_structure = structures,
    serial_component_names = names(structures),
    coefficients = beta_named, fixed_effects = beta_named, beta = dec$beta,
    coefficient_table = coefficient_table,
    vcov = beta_vcov, vcov_full = V,
    random_effects = random_modes, ranef = random_modes,
    serial_modes = serial_modes,
    serial_latent_modes = serial_modes,
    serial_contributions = serial_contributions,
    serial_effects = serial_effects,
    random_covariance = dec$D,
    serial_parameters = serial_parameters,
    serial_parameters_by_component = dec$serial_components,
    family_parameters = dec$family,
    linear_predictor = eta_cond,
    linear_predictor_fixed = eta_fixed,
    fitted.values = mu_cond, fitted_values = mu_cond,
    fitted_fixed = mu_fixed,
    residuals = if (prep$family$name == "binomial") prep$y / prep$trials - mu_cond else prep$y - mu_cond,
    response = prep$y, trials = prep$trials,
    frequency_weights = prep$freq, offset = prep$offset,
    id = prep$id, time = prep$time,
    X = prep$X, Z = prep$Z,
    serial_design = prep$serial_H,
    serial_component_info = prep$serial_components,
    logLik = loglik, negLogLik = details$nll,
    objective = details$objective, penalty = details$penalty,
    AIC = aic, BIC = bic,
    metrics = list(logLik = loglik, AIC = aic, BIC = bic),
    df = df, nobs = length(prep$y),
    convergence = opt$convergence,
    converged = isTRUE(opt$convergence == 0L),
    message = .memwas_null_coalesce(opt$message, ""),
    iterations = opt$iterations,
    L1_penalty = spec$L1_penalty,
    L2_penalty = spec$L2_penalty,
    terms = prep$fixed_terms,
    xlevels = prep$fixed_xlevels,
    contrasts = prep$fixed_contrasts,
    random_terms = prep$random_terms,
    random_xlevels = prep$random_xlevels,
    random_contrasts = prep$random_contrasts,
    training_data = prep$data,
    kept_rows = prep$keep,
    approximate = !(prep$family$name == "gaussian" && identical(prep$family$link, "identity")),
    approximation_label = if (prep$family$name == "gaussian" && identical(prep$family$link, "identity")) {
      "exact Gaussian latent integration evaluated through the Laplace identity"
    } else {
      "joint Laplace approximation"
    },
    engine = "multiple_serial_joint_likelihood",
    autocorrelation = list(type = autocor_label, components = dec$serial_components),
    settings = .memwas_null_coalesce(spec$.settings, NULL),
    .prep = prep, .layout = layout,
    .decoded = dec, .spec = spec
  ), class = c("MEMWAS_multi_serial_fit", "MEMWAS_family_fit", "MEMWAS_fit"))
  fit
}

.memwas_select_serial_components <- function(include_serial, available) {
  if (is.null(include_serial) || (is.logical(include_serial) && length(include_serial) == 1L && !isTRUE(include_serial))) return(character())
  if (is.logical(include_serial) && length(include_serial) == 1L && isTRUE(include_serial)) return(available)
  if (is.character(include_serial)) {
    unknown <- setdiff(include_serial, available)
    if (length(unknown)) stop("Unknown serial component(s): ", paste(unknown, collapse = ", "), call. = FALSE)
    return(unique(include_serial))
  }
  stop("`include_serial` must be FALSE, TRUE, or a character vector of component names.", call. = FALSE)
}

.memwas_validate_serial_prediction_times <- function(serial, time_new, time_old,
                                                tolerance = 1e-7) {
  discrete <- serial$type %in% c("ARP", "ARMA11", "TOEP") ||
    (identical(serial$type, "AR1") && is.finite(serial$rho) && serial$rho < 0)
  if (!discrete) return(invisible(TRUE))
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance <= 0) {
    stop("The serial prediction-time tolerance must be one positive finite number.", call. = FALSE)
  }
  scale <- serial$scale
  if (!is.numeric(scale) || length(scale) != 1L || !is.finite(scale) || scale <= 0) {
    stop("The fitted serial component has an invalid time scale.", call. = FALSE)
  }
  lag <- outer(as.numeric(time_new), as.numeric(time_old), "-") / scale
  if (any(abs(lag - round(lag)) > tolerance)) {
    stop(
      "This discrete-time serial component cannot predict at off-grid times. ",
      "Use times aligned with the fitted grid, or use EXP/OU for continuous-time prediction.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.memwas_new_serial_design <- function(object, newdata) {
  out <- matrix(numeric(), nrow(newdata), 0L)
  if (!length(object$serial_component_info)) return(out)
  out <- matrix(NA_real_, nrow(newdata), length(object$serial_component_info),
                dimnames = list(NULL, names(object$serial_component_info)))
  for (nm in names(object$serial_component_info)) {
    comp <- object$serial_component_info[[nm]]
    mf <- stats::model.frame(comp$terms, newdata,
                             na.action = stats::na.pass,
                             xlev = comp$xlevels)
    mm <- stats::model.matrix(comp$terms, mf,
                              contrasts.arg = comp$contrasts)
    if (ncol(mm) != 1L || nrow(mm) != nrow(newdata)) {
      stop("New data did not reproduce the one-column design for serial component ", sQuote(nm), ".", call. = FALSE)
    }
    out[, nm] <- as.numeric(mm[, 1L])
  }
  out
}

#' @export
predict.MEMWAS_multi_serial_fit <- function(
    object, newdata = NULL, type = c("response", "link"),
    include_random = FALSE, include_serial = include_random,
    id = NULL, time = NULL, offset = NULL,
    allow_new_levels = TRUE, ...) {
  type <- match.arg(type)
  if (!is.logical(include_random) || length(include_random) != 1L || is.na(include_random)) {
    stop("`include_random` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(allow_new_levels) || length(allow_new_levels) != 1L || is.na(allow_new_levels)) {
    stop("`allow_new_levels` must be TRUE or FALSE.", call. = FALSE)
  }
  selected <- .memwas_select_serial_components(include_serial, names(object$serial_parameters_by_component))
  if (is.null(newdata)) {
    eta <- object$linear_predictor_fixed
    if (include_random && ncol(object$Z)) {
      b_by_row <- object$random_effects[as.character(object$id), , drop = FALSE]
      eta <- eta + rowSums(object$Z * b_by_row)
    }
    if (length(selected)) {
      eta <- eta + rowSums(object$serial_contributions[, selected, drop = FALSE])
    }
    return(if (type == "link") eta else .memwas_family_mean(eta, object$.prep$family))
  }

  nd <- as.data.frame(newdata)
  des <- .memwas_new_model_design(object, nd)
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
  if (is.character(idv) && length(idv) == 1L && idv %in% names(nd)) idv <- nd[[idv]]
  id_source <- .memwas_null_coalesce(object$.spec$id, object$.prep$id_spec)
  if (is.null(idv) && !is.null(id_source)) idv <- try(.memwas_resolve_data_value(id_source, nd, parent.frame(), "id"), silent = TRUE)
  if (inherits(idv, "try-error")) idv <- NULL
  if ((include_random || length(selected)) && is.null(idv)) stop("Conditional prediction requires id values for newdata.", call. = FALSE)
  if (!is.null(idv) && length(idv) != nrow(nd)) stop("id must have one value per newdata row.", call. = FALSE)
  if ((include_random || length(selected)) && anyNA(idv)) stop("Conditional prediction does not allow missing id values.", call. = FALSE)
  idc <- if (is.null(idv)) rep(NA_character_, nrow(nd)) else as.character(idv)

  if (include_random && ncol(des$Z)) {
    known <- match(idc, rownames(object$random_effects)); good <- !is.na(known)
    if (any(!good) && !allow_new_levels) stop("newdata contains unseen subject levels.", call. = FALSE)
    if (any(good)) eta[good] <- eta[good] + rowSums(des$Z[good, , drop = FALSE] * object$random_effects[known[good], , drop = FALSE])
  }
  if (length(selected)) {
    tv <- time
    if (is.character(tv) && length(tv) == 1L && tv %in% names(nd)) tv <- nd[[tv]]
    time_source <- .memwas_null_coalesce(object$.spec$time, object$.prep$time_spec)
    if (is.null(tv) && !is.null(time_source)) tv <- try(.memwas_resolve_data_value(time_source, nd, parent.frame(), "time"), silent = TRUE)
    if (inherits(tv, "try-error") || is.null(tv)) stop("Serial conditional prediction requires time values for newdata.", call. = FALSE)
    if (length(tv) != nrow(nd)) stop("time must have one value per newdata row.", call. = FALSE)
    if (anyNA(tv)) stop("Serial conditional prediction does not allow missing time values.", call. = FALSE)
    tv <- as.numeric(tv)
    if (any(!is.finite(tv))) stop("Serial prediction times must be finite numeric values.", call. = FALSE)
    Hnew <- .memwas_new_serial_design(object, nd)
    if (any(!is.finite(Hnew[, selected, drop = FALSE]))) {
      stop("Selected serial-component designs contain missing or non-finite values in newdata.", call. = FALSE)
    }
    for (lev in unique(idc)) {
      ni <- which(idc == lev)
      oi <- which(as.character(object$id) == lev)
      if (!length(oi)) {
        if (!allow_new_levels) stop("newdata contains unseen subject level ", sQuote(lev), ".", call. = FALSE)
        next
      }
      for (nm in selected) {
        s <- object$serial_parameters_by_component[[nm]]
        pred_u <- try({
          prediction_tolerance <- .memwas_null_coalesce(
            object$.spec$control$prediction_time_tolerance, 1e-7
          )
          .memwas_validate_serial_prediction_times(
            s, tv[ni], object$time[oi], prediction_tolerance
          )
          Koo <- .memwas_serial_covariance_matrix(s, object$time[oi])
          Kno <- .memwas_serial_covariance_matrix(s, tv[ni], object$time[oi])
          as.vector(Kno %*% .memwas_solve_pd(Koo, object$serial_modes[oi, nm]))
        }, silent = TRUE)
        if (inherits(pred_u, "try-error")) {
          stop("Serial prediction failed for component ", sQuote(nm), ": ", as.character(pred_u), call. = FALSE)
        }
        eta[ni] <- eta[ni] + Hnew[ni, nm] * pred_u
      }
    }
  }
  if (type == "link") eta else .memwas_family_mean(eta, object$.prep$family)
}

#' Simulate from a fitted multiple-serial MEMWAS model
#'
#' Draw responses from the fitted fixed effects, conventional random effects,
#' every named serial component, and the selected response family.
#'
#' @param object A fitted `MEMWAS_multi_serial_fit` object.
#' @param nsim Positive integer number of response vectors to simulate.
#' @param seed Optional random-number seed.
#' @param ... Reserved for compatibility with the `simulate()` generic.
#'
#' @return A data frame with one simulated response vector per column.
#' @noRd
.simulate.MEMWAS_multi_serial_fit <- function(object, nsim = 1L, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  if (!is.numeric(nsim) || length(nsim) != 1L || !is.finite(nsim) ||
      nsim < 1 || nsim != as.integer(nsim)) {
    stop("`nsim` must be one positive integer.", call. = FALSE)
  }
  nsim <- as.integer(nsim)
  prep <- object$.prep; dec <- object$.decoded
  groups <- split(seq_len(object$nobs), object$id)
  ans <- matrix(NA_real_, object$nobs, nsim)
  for (ss in seq_len(nsim)) {
    eta <- object$linear_predictor_fixed
    for (ii in groups) {
      if (ncol(object$Z)) {
        Lb <- t(.memwas_chol_pd(dec$D))
        b <- as.vector(Lb %*% stats::rnorm(ncol(object$Z)))
        eta[ii] <- eta[ii] + as.vector(object$Z[ii, , drop = FALSE] %*% b)
      }
      for (nm in names(dec$serial_components)) {
        Ku <- .memwas_serial_covariance_matrix(dec$serial_components[[nm]], object$time[ii])
        u <- as.vector(t(.memwas_chol_pd(Ku)) %*% stats::rnorm(length(ii)))
        eta[ii] <- eta[ii] + object$serial_design[ii, nm] * u
      }
    }
    mu <- .memwas_family_mean(eta, prep$family)
    ans[, ss] <- switch(object$family_name,
      gaussian = stats::rnorm(object$nobs, mu, dec$family$sigma),
      binomial = stats::rbinom(object$nobs, object$trials, mu),
      poisson = stats::rpois(object$nobs, mu),
      negative_binomial = stats::rnbinom(object$nobs, size = dec$family$theta, mu = mu),
      gamma = stats::rgamma(object$nobs, shape = dec$family$shape, scale = mu / dec$family$shape),
      exponential = stats::rexp(object$nobs, rate = 1 / mu)
    )
  }
  ans <- as.data.frame(ans)
  names(ans) <- paste0("sim_", seq_len(nsim))
  ans
}

#' @export
print.MEMWAS_multi_serial_fit <- function(x, ...) {
  cat("MEMWAS mixed-effects fit with multiple serial components\n")
  cat(" Family:", x$family_name, "(", x$link, ")\n")
  cat(" Approximation:", x$approximation, " Method:", x$method, "\n")
  cat(" Observations:", x$nobs, " logLik:", format(x$logLik, digits = 6), "\n")
  cat(" Serial components:\n")
  for (nm in names(x$serial_structure)) {
    info <- x$serial_component_info[[nm]]
    cat("  -", nm, ":", x$serial_structure[[nm]], "on", paste(deparse(info$design_formula), collapse = " "), "\n")
  }
  cat(" Converged:", x$converged, if (nzchar(x$message)) paste0(" (", x$message, ")") else "", "\n")
  print(x$coefficients)
  invisible(x)
}

#' @export
summary.MEMWAS_multi_serial_fit <- function(object, ...) {
  base <- summary.MEMWAS_family_fit(object, ...)
  component_table <- data.frame(
    component = names(object$serial_structure),
    structure = unname(object$serial_structure),
    design = vapply(object$serial_component_info, function(z) paste(deparse(z$design_formula), collapse = " "), character(1)),
    stringsAsFactors = FALSE
  )
  base$serial_components <- component_table
  class(base) <- c("summary.MEMWAS_multi_serial_fit", class(base))
  base
}

#' @export
print.summary.MEMWAS_multi_serial_fit <- function(x, ...) {
  cat("MEMWAS mixed-effects model with predictor-specific serial processes\n")
  cat("Family:", x$family, " Link:", x$link, " Approximation:", x$approximation, "\n\n")
  cat("Serial components:\n")
  print(x$serial_components, row.names = FALSE)
  cat("\nFixed effects:\n")
  stats::printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
  cat("\nlogLik:", format(x$logLik, digits = 6),
      " AIC:", format(x$AIC, digits = 6),
      " BIC:", format(x$BIC, digits = 6), "\n")
  invisible(x)
}

.memwas_get_raw_argument <- function(raw, aliases, envir, data = NULL,
                                   default = NULL, column_nse = FALSE) {
  expr <- .memwas_find_named_argument(raw, aliases)
  if (is.null(expr)) return(default)
  .memwas_evaluate_argument(expr, envir, data, column_nse = column_nse)
}

.memwas_multi_serial_spec_from_settings <- function(settings, raw, envir, call, serial_value) {
  formula_override <- !is.null(.memwas_find_named_argument(raw, c("formula", "fixed", "fixed_formula")))
  data_override <- !is.null(.memwas_find_named_argument(raw, c("data", "dataset", "data_frame")))
  data <- .memwas_get_raw_argument(raw, c("data", "dataset", "data_frame"), envir,
                                 default = settings$data)
  if (!is.data.frame(data)) data <- as.data.frame(data)
  formula <- .memwas_get_raw_argument(raw, c("formula", "fixed", "fixed_formula"), envir, data,
                                    default = settings$formula)
  if (!formula_override && !data_override && length(settings$spline_variables) &&
      exists(".memwas_add_stored_splines", mode = "function")) {
    aug <- try(.memwas_add_stored_splines(
      settings$data, settings$formula,
      settings$spline_variables, settings$spline_info,
      engine = .memwas_null_coalesce(settings$engine, "R")
    ), silent = TRUE)
    if (!inherits(aug, "try-error")) {
      data <- aug$data; formula <- aug$formula
    }
  }
  if (!inherits(formula, "formula")) stop("A settings object used with multiple serial components must contain a formula.", call. = FALSE)
  family_value <- .memwas_get_raw_argument(raw, c("family", "distribution"), envir, data,
                                         default = settings$family)
  link_value <- .memwas_get_raw_argument(raw, c("link"), envir, data, default = NULL)
  theta <- .memwas_get_raw_argument(raw, c("theta", "size", "nb_size"), envir, data,
                                  default = .memwas_null_coalesce(settings$theta, NULL))
  shape <- .memwas_get_raw_argument(raw, c("shape", "gamma_shape"), envir, data,
                                  default = .memwas_null_coalesce(settings$shape, NULL))
  control <- .memwas_get_raw_argument(raw, c("control"), envir, data,
                                    default = .memwas_null_coalesce(settings$control, list()))
  if (!is.list(control)) stop("control must be a list.", call. = FALSE)
  spec <- list(
    formula = formula, data = data,
    family = .memwas_as_family_spec(family_value, link_value, theta, shape),
    random = .memwas_get_raw_argument(raw, c("random", "random_formula", "re_formula"), envir, data,
                                    default = settings$random),
    id = .memwas_get_raw_argument(raw, c("id", "subject", "subject_id", "cluster", "group"), envir, data,
                                default = settings$id, column_nse = TRUE),
    time = .memwas_get_raw_argument(raw, c("time", "time_var", "visit", "occasion"), envir, data,
                                  default = settings$time, column_nse = TRUE),
    autocor = "NONE", serial = serial_value,
    random_cov = .memwas_get_raw_argument(raw, c("random_cov", "random_covariance", "covariance"), envir, data,
                                        default = .memwas_null_coalesce(settings$random_cov, "unstructured")),
    approximation = .memwas_get_raw_argument(raw, c("approximation"), envir, data,
                                           default = .memwas_null_coalesce(settings$approximation, "laplace")),
    method = .memwas_get_raw_argument(raw, c("method"), envir, data,
                                    default = .memwas_null_coalesce(settings$method, "ML")),
    quadrature_points = .memwas_get_raw_argument(raw, c("quadrature_points", "nAGQ", "nodes"), envir, data,
                                                        default = 9L),
    L1_penalty = .memwas_get_raw_argument(raw, c("L1_penalty", "lambda1", "l1"), envir, data,
                                                 default = .memwas_null_coalesce(settings$L1_penalty, 0)),
    L2_penalty = .memwas_get_raw_argument(raw, c("L2_penalty", "lambda2", "l2"), envir, data,
                                                 default = .memwas_null_coalesce(settings$L2_penalty, 0)),
    offset = .memwas_get_raw_argument(raw, c("offset"), envir, data,
                                    default = .memwas_null_coalesce(settings$offset, NULL), column_nse = TRUE),
    weights = .memwas_get_raw_argument(raw, c("weights", "prior_weights"), envir, data,
                                     default = .memwas_null_coalesce(settings$weights, NULL), column_nse = TRUE),
    subset = .memwas_get_raw_argument(raw, c("subset"), envir, data,
                                    default = .memwas_null_coalesce(settings$subset, NULL)),
    na.action = .memwas_get_raw_argument(raw, c("na.action"), envir, data,
                                      default = .memwas_null_coalesce(settings$na.action, stats::na.omit)),
    theta = theta, shape = shape,
    start = .memwas_get_raw_argument(raw, c("start", "initial"), envir, data,
                                   default = .memwas_null_coalesce(settings$start, NULL)),
    control = control, call = call, environment = envir,
    .settings = settings
  )
  spec
}

.memwas_previous_fit_MEMWAS <- fit_MEMWAS

# Public dispatcher multiple-serial processes. The complete formal interface is
# retained so positional matching, help pages, and inherited non-standard evaluation
# continue to behave as the core MEMWAS did with single-autocorrelation.
fit_MEMWAS <- function(
    object = NULL,
    formula = NULL,
    family = NULL,
    data = NULL, id = NULL, time = NULL,
    random = NULL,
    autocor = NULL,
    serial = NULL,
    L1_penalty = NULL, L2_penalty = NULL,
    control = NULL,
    method = NULL,
    random_cov = NULL,
    approximation = NULL,
    init_approximation = NULL,
    se_method = NULL,
    dot_predictors = NULL, dot_alternative = NULL,
    dot_threshold = 0, dot_alpha = 0.05,
    spline_variables = NULL, spline_info = NULL,
    turning_points = NULL, nonlinear_summary = NULL,
    all_screened_spline_info = NULL,
    baseline_screen_metrics = NULL,
    engine = NULL,
    verbose = NULL, ...) {
  envir <- parent.frame()
  original_call <- sys.call()
  public_call <- match.call(expand.dots = TRUE)
  raw <- .memwas_normalize_fit_call(as.list(original_call)[-1L], envir)

  serial_expr <- .memwas_find_named_argument(raw, c(
    "serial", "serial_components", "predictor_serial",
    "autocor_by_predictor", "predictor_autocor"
  ))
  autocor_expr <- .memwas_find_named_argument(raw, c(
    "autocor", "autocorrelation", "correlation", "cor_struct",
    "correlation_structure", "correlation_type", "correlation_struct"
  ))
  object_expr <- .memwas_find_named_argument(raw, "object")
  formula_expr <- .memwas_find_named_argument(raw, c(
    "formula", "fixed", "fixed_formula", "fixed.effects", "fixed_effects"
  ))
  data_expr <- .memwas_find_named_argument(raw, c(
    "data", "data_frame", "dataset", "data_long", "long_data"
  ))

  object_value <- NULL
  if (!is.null(object_expr)) {
    object_value <- try(.memwas_evaluate_argument(object_expr, envir), silent = TRUE)
    if (inherits(object_value, "try-error")) object_value <- NULL
  }

  # A formula in the core argument of set_MEMWAS and fit_MEMWAS is the long-
  # standing direct-formula linking.
  if (inherits(object_value, "formula")) {
    formula_value <- if (is.null(formula_expr)) NULL else
      try(.memwas_evaluate_argument(formula_expr, envir), silent = TRUE)
    if (!is.null(formula_expr) &&
        (inherits(formula_value, "data.frame") || is.matrix(formula_value)) &&
        is.null(data_expr)) {
      names(raw)[which(names(raw) == "formula")[1L]] <- "data"
      names(raw)[which(names(raw) == "object")[1L]] <- "formula"
      formula_expr <- object_expr
      data_expr <- .memwas_find_named_argument(raw, "data")
      object_value <- NULL
    } else if (is.null(formula_expr)) {
      names(raw)[which(names(raw) == "object")[1L]] <- "formula"
      formula_expr <- object_expr
      object_value <- NULL
    } else {
      stop("A model formula was supplied both as `object` and `formula`.", call. = FALSE)
    }
  }

  serial_present <- !is.null(serial_expr)
  serial_value <- if (serial_present) .memwas_evaluate_argument(serial_expr, envir) else NULL
  autocor_present <- !is.null(autocor_expr)
  autocor_value <- if (autocor_present) .memwas_evaluate_argument(autocor_expr, envir) else NULL

  # Explicit serial = NULL deliberately disables a serial definition stored in
  # a settings object.  Otherwise inherit either the new `serial` field or a
  # multi-valued backward-compatible autocor field.
  if (!serial_present && inherits(object_value, "MEMWAS")) {
    if (!is.null(object_value$serial)) {
      serial_value <- object_value$serial
    } else if (!is.null(object_value$extra$serial)) {
      # Compatibility with settings objects where `serial` was captured
      # through `...`.
      serial_value <- object_value$extra$serial
    } else if (.memwas_is_multi_serial_autocor(object_value$autocor)) {
      serial_value <- object_value$autocor
    }
  }
  if (is.null(serial_value) && .memwas_is_multi_serial_autocor(autocor_value)) {
    serial_value <- autocor_value
  }
  if (!is.null(serial_value) && autocor_present &&
      !.memwas_is_multi_serial_autocor(autocor_value) &&
      !.memwas_is_none_serial(autocor_value)) {
    stop(
      "Supply all serial processes through `serial`; do not also supply a non-NONE `autocor`.",
      call. = FALSE
    )
  }

  multi_requested <- !is.null(serial_value)
  if (!multi_requested) {
    # Reconstruct the delegated call from the normalized, still-unevaluated
    # expressions.  In particular, this preserves the historical
    # fit_MEMWAS(y ~ x, data = dat, ...) positional formula API for the
    # single-structure family wrapper, whose own `...` call otherwise labels
    # that first expression as `object` before its raw parser sees it.
    legacy_call <- as.call(c(
      list(as.name(".memwas_previous_fit_MEMWAS")), raw
    ))
    eval_env <- new.env(parent = envir)
    eval_env$.memwas_previous_fit_MEMWAS <- .memwas_previous_fit_MEMWAS
    return(eval(legacy_call, envir = eval_env))
  }

  .memwas_validate_serial_syntax(serial_value)
  if (inherits(object_value, "MEMWAS")) {
    spec <- .memwas_multi_serial_spec_from_settings(
      object_value, raw, envir, public_call, serial_value
    )
  } else {
    spec <- .memwas_parse_fit_call(raw, envir, public_call)
    spec$serial <- serial_value
    spec$autocor <- "NONE"
    approx_expr <- .memwas_find_named_argument(raw, "approximation")
    spec$approximation <- if (is.null(approx_expr)) "laplace" else
      .memwas_evaluate_argument(approx_expr, envir)
    method_expr <- .memwas_find_named_argument(raw, "method")
    spec$method <- if (is.null(method_expr)) "ML" else
      .memwas_evaluate_argument(method_expr, envir)
  }
  fit <- .memwas_fit_multi_serial_engine(spec)
  if (inherits(object_value, "MEMWAS")) {
    fit$settings <- object_value
    fit$formal_formula <- spec$formula
  }
  fit
}
