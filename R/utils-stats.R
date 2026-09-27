# Shared statistical helpers -------------------------------------------------
#
# One implementation of each diagnostic, so that (for example) the sample-size
# guard on Shapiro-Wilk behaves identically in every function that calls it.

#' Shapiro-Wilk test of residual normality with sample-size guards
#'
#' \code{stats::shapiro.test()} errors below 3 and above 5000 observations.
#' Rather than let that error escape from a diagnostic nobody asked for, this
#' returns \code{NULL} and a note.
#'
#' @param x numeric vector, or a matrix (which is flattened explicitly rather
#'   than left to \code{shapiro.test}'s recycling behaviour).
#' @return list with \code{test} (an \code{htest} or \code{NULL}), \code{n} and
#'   \code{note} (character, possibly empty).
#' @noRd
.check_normality <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]                  # drops NA, NaN and +/-Inf alike
  n <- length(x)
  if (n < 3L) {
    return(list(test = NULL, n = n, note = sprintf(
      "Shapiro-Wilk skipped: %d usable residual(s), at least 3 are required.", n)))
  }
  if (n > 5000L) {
    return(list(test = NULL, n = n, note = sprintf(
      "Shapiro-Wilk skipped: %d residuals exceeds the 5000 the test supports. Inspect the Q-Q plot instead.", n)))
  }
  sdx <- stats::sd(x)
  if (!is.finite(sdx) || sdx == 0) {
    return(list(test = NULL, n = n,
                note = "Shapiro-Wilk skipped: the values have no usable variation."))
  }
  list(test = stats::shapiro.test(x), n = n, note = character(0))
}

#' Normality of the response within each group
#'
#' The right diagnostic for a method that permits unequal variances: pooling
#' residuals across groups with different spreads produces a mixture that fails
#' normality tests even when every group is normal.
#'
#' @param values numeric response.
#' @param cell factor of group memberships.
#' @noRd
.normality_by_group <- function(values, cell) {
  cell <- droplevels(as.factor(cell))
  rows <- lapply(levels(cell), function(lv) {
    v <- values[cell == lv]
    nm <- .check_normality(v)
    data.frame(
      group     = lv,
      n         = nm$n,        # the number the test actually used, not length(v)
      statistic = if (is.null(nm$test)) NA_real_ else unname(nm$test$statistic),
      p_value   = if (is.null(nm$test)) NA_real_ else unname(nm$test$p.value),
      # Why a group has no test (too few or too many values, or no variation);
      # an NA with no reason reads as a failure.
      note      = if (length(nm$note) == 0L) NA_character_ else nm$note,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

#' A p-value for a sentence: "p = 0.012", or "p < 2.2e-16" at machine precision
#' @noRd
.fmt_p <- function(p, digits = 2L) {
  eps <- .Machine$double.eps
  ifelse(!is.na(p) & p < eps,
         paste("p <", format(eps, digits = digits)),
         paste("p =", format(signif(p, digits))))
}

#' Is a Pearson dispersion above 1 by more than chance?
#'
#' Under the model the Pearson statistic is roughly chi-square on the residual
#' degrees of freedom, so a dispersion a little above 1 is expected half the
#' time. A note about overdispersion is worth giving only when the statistic
#' is beyond what that distribution produces.
#' @return the upper-tail p-value, or \code{NA} when it cannot be computed.
#' @noRd
.overdispersion_p <- function(dispersion, df) {
  if (!is.finite(dispersion) || is.null(df) || !is.finite(df) || df < 1) {
    return(NA_real_)
  }
  stats::pchisq(dispersion * df, df, lower.tail = FALSE)
}

#' Pearson dispersion statistic for a fitted GLM
#'
#' Pearson rather than the deviance ratio, which is biased for small counts, and
#' the same statistic in every function that reports dispersion.
#' @noRd
.dispersion <- function(fit) {
  pr <- stats::residuals(fit, type = "pearson")
  dfr <- stats::df.residual(fit)
  if (is.null(dfr) || dfr <= 0) return(NA_real_)
  sum(pr^2, na.rm = TRUE) / dfr
}

#' Partial eta squared and partial omega squared from an ANOVA table
#'
#' Drops both the intercept row and the residual row. Leaving the residual row
#' in produces an entry of exactly 0.5 for every model, which is meaningless.
#'
#' Both measures are \emph{partial}. The classical (non-partial) omega squared,
#' \code{(SS - df * MSE) / (SS_total + MSE)}, is only a proportion of variance
#' when the effect sums of squares partition the total, which Type II and Type
#' III sums of squares do not: on a near-aliased design it can sum well above 1.
#' Partial omega squared,
#' \code{df (MS - MSE) / (df MS + (N - df) MSE)}, is defined per effect and
#' needs no such partition, and it is the companion measure to partial eta
#' squared rather than a different kind of quantity sitting in the next column.
#'
#' @param aov_table a \code{car::Anova} table (or any data.frame with a
#'   \code{"Sum Sq"} column and a \code{Residuals} row).
#' @param ss_error residual sum of squares; taken from the table when absent.
#' @param df_error residual degrees of freedom, needed for omega squared.
#' @param n_obs number of observations, needed for omega squared.
#' @noRd
.partial_eta_squared <- function(aov_table, ss_error = NULL, df_error = NULL,
                                 n_obs = NULL) {
  tab <- as.data.frame(aov_table)
  if (!"Sum Sq" %in% names(tab)) {
    return(data.frame(term = character(0), df = numeric(0),
                      sum_sq = numeric(0), partial_eta_sq = numeric(0),
                      stringsAsFactors = FALSE))
  }
  rn <- rownames(tab)
  if (is.null(ss_error) && "Residuals" %in% rn) {
    ss_error <- tab[["Sum Sq"]][rn == "Residuals"]
  }
  if (is.null(df_error) && "Residuals" %in% rn && "Df" %in% names(tab)) {
    df_error <- tab[["Df"]][rn == "Residuals"]
  }
  drop_rows <- rn %in% c("(Intercept)", "Residuals")
  eff <- tab[!drop_rows, , drop = FALSE]
  if (nrow(eff) == 0L || is.null(ss_error)) {
    return(data.frame(term = character(0), df = numeric(0),
                      sum_sq = numeric(0), partial_eta_sq = numeric(0),
                      stringsAsFactors = FALSE))
  }
  ss <- eff[["Sum Sq"]]
  df <- if ("Df" %in% names(eff)) eff[["Df"]] else rep(NA_real_, length(ss))
  out <- data.frame(
    term           = gsub("`", "", rownames(eff), fixed = TRUE),
    df             = as.numeric(df),
    sum_sq         = as.numeric(ss),
    partial_eta_sq = as.numeric(ss / (ss + ss_error)),
    stringsAsFactors = FALSE,
    row.names      = NULL
  )
  if (!is.null(df_error) && !is.null(n_obs) && is.finite(df_error) &&
      df_error > 0 && is.finite(n_obs) && n_obs > 0) {
    mse <- ss_error / df_error
    ms  <- ss / out$df
    om  <- (out$df * (ms - mse)) / (out$df * ms + (n_obs - out$df) * mse)
    out$partial_omega_sq <- pmax(0, pmin(1, om))
    if (any(om < 0, na.rm = TRUE)) {
      attr(out, "omega_floored") <- out$term[which(om < 0)]
    }
  }
  out
}

#' Group-wise mean, sd, n, se and confidence interval
#' @noRd
.summary_stats <- function(values, cell, conf_level = 0.95) {
  cell <- droplevels(as.factor(cell))
  lv <- levels(cell)
  rows <- lapply(lv, function(g) {
    v <- values[cell == g]
    v <- v[!is.na(v)]
    n <- length(v)
    m <- if (n > 0L) mean(v) else NA_real_
    s <- if (n > 1L) stats::sd(v) else NA_real_
    se <- if (n > 1L) s / sqrt(n) else NA_real_
    tcrit <- if (n > 1L) stats::qt(1 - (1 - conf_level) / 2, df = n - 1L) else NA_real_
    data.frame(group = g, n = n, mean = m, sd = s, se = se,
               conf_low  = if (is.na(se)) NA_real_ else m - tcrit * se,
               conf_high = if (is.na(se)) NA_real_ else m + tcrit * se,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$group <- factor(out$group, levels = lv)
  row.names(out) <- NULL
  out
}

#' Mardia's tests of multivariate skewness and kurtosis
#'
#' Implemented directly rather than taken from a package, so that multivariate
#' normality can be assessed without adding a dependency.
#'
#' Both statistics are functions of the Mahalanobis forms
#' \eqn{d_{ij} = y_i' S^{-1} y_j} (ML covariance, divisor \eqn{n}): skewness
#' \eqn{b_1 = \sum_{ij} d_{ij}^3 / n^2} and kurtosis
#' \eqn{b_2 = \sum_i d_{ii}^2 / n}. The \eqn{n \times n} matrix of forms is
#' never built: with whitened scores \eqn{Z} (so that \eqn{d_{ij} = z_i'z_j}),
#' \eqn{b_1 = \sum_{abc} (\sum_i z_{ia} z_{ib} z_{ic})^2 / n^2}, which takes
#' \eqn{O(n p^3)} time and \eqn{O(n p)} memory. The forms do not depend on the
#' scale of the columns, so the columns are standardised first, which keeps
#' responses on very different scales from making \eqn{S} look singular.
#'
#' Applied to model residuals, the statistics carry no information when the
#' residual degrees of freedom equal the number of columns (the forms are then
#' a fixed function of the design), and close to it they are dominated by the
#' design; the test is skipped, with a note, unless the residual degrees of
#' freedom exceed \code{p} by at least \code{margin}.
#'
#' @param Y numeric matrix, rows are observations (typically residuals).
#' @param df_resid residual degrees of freedom of the model that produced
#'   \code{Y}, or \code{NULL} for raw data (\code{n - 1}).
#' @param margin how far \code{df_resid} must exceed the number of columns.
#' @return list with \code{table} (data.frame with one row per test, or
#'   \code{NULL}) and \code{note}.
#' @references Mardia, K. V. (1970). Measures of multivariate skewness and
#'   kurtosis with applications. \emph{Biometrika}, 57(3), 519-530.
#' @noRd
.mardia_test <- function(Y, df_resid = NULL, margin = 10L) {
  Y <- as.matrix(Y)
  Y <- Y[stats::complete.cases(Y), , drop = FALSE]
  n <- nrow(Y); p <- ncol(Y)
  dfr <- if (is.null(df_resid)) n - 1L else df_resid
  if (n < 4L || !is.finite(dfr) || dfr < p + margin) {
    return(list(table = NULL, note = sprintf(
      "Mardia's tests were not computed: the residuals have %s degree(s) of freedom for %d responses, and at least %d are needed. At %d they would not depend on the data at all, and close to it they are dominated by the design.",
      format(dfr), p, p + margin, p)))
  }
  singular <- list(table = NULL, note = "Mardia's tests were not computed: the residual covariance matrix is singular (the responses are collinear).")
  Yc <- sweep(Y, 2L, colMeans(Y), "-")
  s <- sqrt(colSums(Yc^2) / n)
  if (any(!is.finite(s) | s <= 0)) return(singular)
  Yc <- sweep(Yc, 2L, s, "/")                  # the forms are scale-invariant
  S  <- crossprod(Yc) / n                      # ML covariance, Mardia's divisor
  ev <- eigen((S + t(S)) / 2, symmetric = TRUE)
  if (min(ev$values) <= max(ev$values) * 1e-12) return(singular)
  # Whitened scores: Z Z' = Yc S^-1 Yc', the matrix of Mahalanobis forms.
  Z <- Yc %*% ev$vectors %*% diag(1 / sqrt(ev$values), p)

  b1 <- 0
  for (a in seq_len(p)) b1 <- b1 + sum(crossprod(Z, Z * Z[, a])^2)
  b1 <- b1 / n^2
  b2 <- sum(rowSums(Z^2)^2) / n

  skew_stat <- n * b1 / 6
  skew_df   <- p * (p + 1) * (p + 2) / 6
  skew_p    <- stats::pchisq(skew_stat, df = skew_df, lower.tail = FALSE)

  kurt_stat <- (b2 - p * (p + 2)) / sqrt(8 * p * (p + 2) / n)
  kurt_p    <- 2 * stats::pnorm(-abs(kurt_stat))

  list(table = data.frame(
    test      = c("Mardia skewness", "Mardia kurtosis"),
    statistic = c(skew_stat, kurt_stat),
    df        = c(skew_df, NA_real_),
    p_value   = c(skew_p, kurt_p),
    stringsAsFactors = FALSE
  ), note = character(0))
}

#' Box's M test for homogeneity of covariance matrices
#'
#' @param Y numeric matrix of responses.
#' @param group factor of group memberships.
#' @return list with \code{statistic}, \code{df}, \code{p_value} and
#'   \code{note}, or \code{NULL} when it cannot be computed.
#' @references Box, G. E. P. (1949). A general distribution theory for a class
#'   of likelihood criteria. \emph{Biometrika}, 36(3/4), 317-346.
#' @noRd
.box_m <- function(Y, group) {
  Y <- as.matrix(Y)
  group <- droplevels(as.factor(group))
  ok <- stats::complete.cases(Y) & !is.na(group)
  Y <- Y[ok, , drop = FALSE]; group <- droplevels(group[ok])
  p <- ncol(Y); k <- nlevels(group); N <- nrow(Y)
  if (k < 2L) return(NULL)
  ns <- as.integer(table(group))
  if (any(ns <= p)) {
    return(list(statistic = NA_real_, df = NA_real_, p_value = NA_real_,
                note = sprintf(
                  "Box's M not computed: every group needs more than %d observations (smallest group has %d).",
                  p, min(ns))))
  }
  covs <- lapply(levels(group), function(g) stats::cov(Y[group == g, , drop = FALSE]))
  pooled <- Reduce(`+`, Map(function(S, n) (n - 1) * S, covs, ns)) / (N - k)
  dets <- vapply(covs, function(S) determinant(S, logarithm = TRUE)$modulus, numeric(1))
  det_pooled <- determinant(pooled, logarithm = TRUE)$modulus
  if (!all(is.finite(c(dets, det_pooled)))) {
    return(list(statistic = NA_real_, df = NA_real_, p_value = NA_real_,
                note = "Box's M not computed: at least one covariance matrix is singular."))
  }
  M  <- (N - k) * as.numeric(det_pooled) - sum((ns - 1) * as.numeric(dets))
  c1 <- (sum(1 / (ns - 1)) - 1 / (N - k)) *
        (2 * p^2 + 3 * p - 1) / (6 * (p + 1) * (k - 1))
  stat <- M * (1 - c1)
  df <- (k - 1) * p * (p + 1) / 2
  list(statistic = stat, df = df,
       p_value = stats::pchisq(stat, df = df, lower.tail = FALSE),
       note = character(0))
}

#' Canonical discriminant analysis from hypothesis and error SSCP matrices
#'
#' Solves the generalised eigenproblem \eqn{H v = \lambda E v} and projects
#' the centred responses onto the leading discriminant axes. The problem is
#' solved on responses scaled to unit error variance
#' (\eqn{D^{-1} E D^{-1}} and \eqn{D^{-1} H D^{-1}}, with \eqn{D} the square
#' roots of \code{diag(E)}), which leaves the eigenvalues unchanged and keeps
#' responses on very different scales from making \code{E} look singular; the
#' eigenvectors are transformed back.
#'
#' @param Y numeric matrix of responses. It must live in the same space as
#'   \code{H} and \code{E}: with covariates, or with other terms in the model,
#'   their fitted effects must already have been removed, so that the
#'   deviations of \code{Y} from its means within the levels of the described
#'   term are the model residuals. Projecting the raw responses onto
#'   eigenvectors derived from an adjusted model produces scores that
#'   discriminate far worse than the eigenvalue printed beside them.
#' @param H hypothesis SSCP matrix.
#' @param E error SSCP matrix.
#' @param df_h degrees of freedom of the hypothesis term.
#' @param df_e residual degrees of freedom, used to scale the axes to unit
#'   pooled within-group variance (\eqn{v'Ev / df_e = 1}), which is the
#'   convention that makes the spread of the group centroids on the plot mean
#'   what it appears to mean.
#' @param group optional factor of the described term's levels, returned with
#'   the scores for plotting.
#' @return list with \code{scores}, \code{group}, \code{canonical} (axis,
#'   eigenvalue, canonical correlation, proportion of the term's between-group
#'   variation) and \code{structure} (a \code{response} column and one column
#'   per axis: the pooled within-group correlations between each response and
#'   each axis, \eqn{diag(E)^{-1/2} E v / \sqrt{v'Ev}}, SAS's "pooled within
#'   canonical structure"); or, on failure, a list whose \code{canonical} is
#'   \code{NULL} and whose \code{reason} says why.
#' @noRd
.canonical_discriminant <- function(Y, H, E, df_h, df_e = NULL, group = NULL) {
  Y <- as.matrix(Y)
  fail <- function(reason) list(canonical = NULL, reason = reason)
  if (is.null(H) || is.null(E)) {
    return(fail("its hypothesis or error matrix is unavailable"))
  }
  if (!is.finite(df_h) || df_h < 1) {
    return(fail("the term has no estimable degrees of freedom"))
  }
  s <- sqrt(diag(E))
  if (any(!is.finite(s) | s <= 0)) return(fail("the error matrix is singular"))
  Es <- E / outer(s, s)
  Hs <- H / outer(s, s)
  ee <- eigen((Es + t(Es)) / 2, symmetric = TRUE)
  if (min(ee$values) <= max(ee$values) * 1e-12) {
    return(fail("the error matrix is singular"))
  }
  # Symmetric square root of the scaled error matrix, and its inverse.
  Es_half  <- ee$vectors %*% (t(ee$vectors) * sqrt(ee$values))
  Es_ihalf <- ee$vectors %*% (t(ee$vectors) / sqrt(ee$values))
  M <- Es_ihalf %*% Hs %*% Es_ihalf
  ev <- eigen((M + t(M)) / 2, symmetric = TRUE)
  vals <- ev$values
  keep <- min(df_h, ncol(Y), sum(vals > .Machine$double.eps^0.5))
  if (keep < 1L) return(fail("the term does not separate the groups on any axis"))
  vals <- vals[seq_len(keep)]
  W <- ev$vectors[, seq_len(keep), drop = FALSE]
  # Structure coefficients diag(E)^-1/2 E v / sqrt(v'Ev), which for these
  # eigenvectors is Es^1/2 w. The sign of an eigenvector is arbitrary: make
  # each axis's largest structure coefficient positive, so the orientation of
  # the plot is reproducible.
  struct <- Es_half %*% W
  big <- max.col(t(abs(struct)), ties.method = "first")
  flip <- sign(struct[cbind(big, seq_len(keep))])
  flip[flip == 0] <- 1
  W <- sweep(W, 2L, flip, "*")
  struct <- sweep(struct, 2L, flip, "*")
  # Discriminant coefficients on the original scale, normalised to v'Ev = 1.
  V <- (Es_ihalf %*% W) / s

  dfe <- if (is.null(df_e) || !is.finite(df_e) || df_e <= 0) {
    max(nrow(Y) - 1L, 1L)
  } else df_e
  Yc <- sweep(Y, 2L, colMeans(Y), "-")
  # v'Ev / df_e is the pooled within-group variance of an axis, so with
  # v'Ev = 1 multiplying by sqrt(df_e) gives unit pooled within-group variance.
  scores <- (Yc %*% V) * sqrt(dfe)
  axes <- paste0("Can", seq_len(keep))
  colnames(scores) <- axes
  colnames(struct) <- axes
  rn <- colnames(Y) %||% rownames(E) %||% paste0("y", seq_len(ncol(Y)))

  list(
    scores = scores,
    group = if (is.null(group)) NULL else droplevels(as.factor(group)),
    canonical = data.frame(
      axis = axes,
      eigenvalue = vals,
      canonical_r = sqrt(vals / (1 + vals)),
      prop_variance = vals / sum(vals),
      stringsAsFactors = FALSE
    ),
    structure = data.frame(response = rn, struct, stringsAsFactors = FALSE,
                           row.names = NULL, check.names = FALSE),
    reason = character(0)
  )
}

#' Detect separation in a binary GLM, or an all-zero cell in a count GLM
#'
#' Wald statistics collapse on the boundary of the parameter space: the odds
#' ratio or rate ratio explodes, the interval covers everything, and the
#' p-value approaches 1. \code{glm()} reports \code{converged = TRUE}
#' regardless, so this has to be checked explicitly rather than inferred from
#' the fit.
#'
#' The check is made on the fitted linear predictor of each distinct cell of
#' the design, not on the coefficients. A coefficient depends on the coding: under
#' sum-to-zero contrasts the divergence of one homogeneous level is spread
#' across the intercept and every deviation coefficient, none of which need
#' look extreme on its own. A cell's linear predictor and its standard error do
#' not depend on the coding. A cell is flagged when its linear predictor has
#' run away (beyond +/-10 on the logit scale, or below -10 on the log scale)
#' together with its standard error (above 10). A homogeneous cell of an
#' additive model whose predictor is estimable from the other cells has a
#' moderate standard error and is not flagged, and neither is a rare but
#' well-determined event rate in a large sample.
#'
#' @param fit a fitted binomial, Poisson or negative binomial GLM.
#' @return a character note, or \code{character(0)}.
#' @noRd
.check_separation <- function(fit) {
  fam <- tryCatch(stats::family(fit)$family, error = function(e) "")
  is_count <- inherits(fit, "negbin") ||
    fam %in% c("poisson", "quasipoisson") || grepl("^Negative Binomial", fam)
  co <- stats::coef(fit)
  ok <- !is.na(co)
  X <- tryCatch(stats::model.matrix(fit), error = function(e) NULL)
  V <- tryCatch(suppressWarnings(stats::vcov(fit)), error = function(e) NULL)
  if (is.null(X) || is.null(V) || !any(ok)) return(character(0))
  keep <- intersect(names(co)[ok], intersect(colnames(X), colnames(V)))
  X <- X[, keep, drop = FALSE]
  V <- V[keep, keep, drop = FALSE]
  key <- do.call(paste, c(as.data.frame(X), list(sep = "\r")))
  first <- !duplicated(key)
  Xu <- X[first, , drop = FALSE]
  eta <- drop(Xu %*% co[keep])
  se_eta <- sqrt(pmax(rowSums((Xu %*% V) * Xu), 0))
  flag <- is.finite(eta) & is.finite(se_eta) & se_eta > 10 &
    (if (is_count) eta < -10 else abs(eta) > 10)
  if (!any(flag)) return(character(0))

  mf <- tryCatch(stats::model.frame(fit), error = function(e) NULL)
  cells <- if (!is.null(mf)) {
    facs <- names(mf)[-1L][vapply(mf[-1L], is.factor, logical(1))]
    if (length(facs) > 0L) {
      sub <- mf[first, facs, drop = FALSE][flag, , drop = FALSE]
      apply(sub, 1L, function(r) paste(sprintf("%s = %s", facs, r), collapse = ", "))
    } else character(0)
  } else character(0)
  where <- if (length(cells) > 0L) .abbrev(unname(cells), 5L) else
    sprintf("%d cell(s) of the design", sum(flag))

  co_ok <- co[keep]
  se_co <- sqrt(pmax(diag(V), 0))
  # Coefficients are worth naming only under treatment coding, where each one
  # is a level against its reference; a sum-to-zero or polynomial coefficient
  # names nothing the reader can find in the data. The cells are named anyway.
  ctr <- fit$contrasts
  treatment <- is.null(ctr) || all(vapply(ctr, function(x)
    identical(x, "contr.treatment"), logical(1)))
  runaway <- if (treatment) {
    setdiff(names(co_ok)[abs(co_ok) > 10 & se_co > 10], "(Intercept)")
  } else character(0)
  coefs <- if (length(runaway) > 0L) {
    sprintf(" The affected coefficient(s): %s.",
            paste(gsub("`", "", runaway, fixed = TRUE), collapse = ", "))
  } else ""

  if (is_count) {
    return(paste0(
      "No events were observed in ", where, ", so the fitted rate there is ",
      "numerically zero.", coefs, " A rate ratio involving such a cell is not ",
      "identified: it will be near zero or enormous, its interval unbounded on one ",
      "side, and its Wald p-value near 1 no matter how large the difference is; ",
      "the same goes for that cell's marginal mean and the comparisons involving it. ",
      "The likelihood-ratio omnibus test remains usable. Consider collapsing the ",
      "level, or an exact or penalised method for the affected comparisons."))
  }
  paste0(
    "Complete or quasi-complete separation detected: the fitted probability is ",
    "numerically 0 or 1 in ", where, ".", coefs,
    " An affected odds ratio is not identified: it will be enormous, its ",
    "interval will be unbounded on one side, and its Wald p-value will be near 1 ",
    "no matter how strong the association is; the same goes for that cell's ",
    "marginal probability and the comparisons involving it. With events this sparse the ",
    "likelihood-ratio omnibus test is also liberal. Consider a penalised fit such ",
    "as logistf::logistf(), or collapsing the offending level.")
}
