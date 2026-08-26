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
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
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
    term           = rownames(eff),
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
      attr(out, "omega_floored") <- rownames(eff)[which(om < 0)]
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
#' @param Y numeric matrix, rows are observations.
#' @return data.frame with one row per test, or \code{NULL} if \code{Y} is too
#'   small or singular.
#' @references Mardia, K. V. (1970). Measures of multivariate skewness and
#'   kurtosis with applications. \emph{Biometrika}, 57(3), 519-530.
#' @noRd
.mardia_test <- function(Y) {
  Y <- as.matrix(Y)
  Y <- Y[stats::complete.cases(Y), , drop = FALSE]
  n <- nrow(Y); p <- ncol(Y)
  if (n < p + 1L || n < 4L) return(NULL)
  Yc <- scale(Y, center = TRUE, scale = FALSE)
  S  <- crossprod(Yc) / n                      # ML covariance, Mardia's divisor
  Si <- tryCatch(solve(S), error = function(e) NULL)
  if (is.null(Si)) return(NULL)

  D  <- Yc %*% Si %*% t(Yc)                    # n x n matrix of Mahalanobis forms
  b1 <- sum(D^3) / (n^2)
  b2 <- sum(diag(D)^2) / n

  skew_stat <- n * b1 / 6
  skew_df   <- p * (p + 1) * (p + 2) / 6
  skew_p    <- stats::pchisq(skew_stat, df = skew_df, lower.tail = FALSE)

  kurt_stat <- (b2 - p * (p + 2)) / sqrt(8 * p * (p + 2) / n)
  kurt_p    <- 2 * stats::pnorm(-abs(kurt_stat))

  data.frame(
    test      = c("Mardia skewness", "Mardia kurtosis"),
    statistic = c(skew_stat, kurt_stat),
    df        = c(skew_df, NA_real_),
    p_value   = c(skew_p, kurt_p),
    stringsAsFactors = FALSE
  )
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
#' Solves the generalised eigenproblem for \code{solve(E) \%*\% H} and projects
#' the centred responses onto the leading discriminant axes.
#'
#' @param Y numeric matrix of responses. For a MANCOVA this must be the
#'   \emph{residuals} of the responses on the covariates, so that the scores
#'   live in the same adjusted space as \code{H} and \code{E}. Projecting the
#'   raw responses onto eigenvectors derived from an adjusted model produces
#'   scores that discriminate far worse than the eigenvalue printed beside
#'   them, and structure coefficients that can have the wrong sign.
#' @param H hypothesis SSCP matrix.
#' @param E error SSCP matrix.
#' @param df_h degrees of freedom of the hypothesis term.
#' @param df_e residual degrees of freedom, used to scale the axes to unit
#'   pooled within-group variance, which is the convention that makes the
#'   spread of the group centroids on the plot mean what it appears to mean.
#' @param group optional factor, used for the pooled within-group structure
#'   coefficients that CDA conventionally reports.
#' @return list with \code{scores}, \code{canonical} and \code{structure}, or
#'   \code{NULL} on failure.
#' @noRd
.canonical_discriminant <- function(Y, H, E, df_h, df_e = NULL, group = NULL) {
  Y <- as.matrix(Y)
  Einv <- tryCatch(solve(E), error = function(e) NULL)
  if (is.null(Einv)) return(NULL)
  ev <- tryCatch(eigen(Einv %*% H, symmetric = FALSE), error = function(e) NULL)
  if (is.null(ev)) return(NULL)
  vals <- Re(ev$values)
  vecs <- Re(ev$vectors)
  ord <- order(vals, decreasing = TRUE)
  vals <- vals[ord]; vecs <- vecs[, ord, drop = FALSE]
  keep <- min(df_h, ncol(Y), sum(vals > .Machine$double.eps^0.5))
  if (keep < 1L) return(NULL)
  vals <- vals[seq_len(keep)]
  vecs <- vecs[, seq_len(keep), drop = FALSE]

  Yc <- scale(Y, center = TRUE, scale = FALSE)
  scores <- Yc %*% vecs

  # Unit pooled within-group variance: v' E v / df_e is the within-group
  # variance of the axis, so dividing by its square root gives unit scale.
  wsd <- sqrt(pmax(diag(t(vecs) %*% E %*% vecs), 0) /
                (if (is.null(df_e) || !is.finite(df_e) || df_e <= 0)
                  max(nrow(Y) - 1L, 1L) else df_e))
  wsd[!is.finite(wsd) | wsd == 0] <- 1
  scores <- sweep(scores, 2L, wsd, "/")
  colnames(scores) <- paste0("Can", seq_len(keep))

  # Pooled within-group correlations between responses and axes: the quantity
  # conventionally called the structure coefficient.
  struct <- if (!is.null(group)) {
    g <- droplevels(as.factor(group))
    resid_within <- function(M) {
      M <- as.matrix(M)
      for (j in seq_len(ncol(M))) M[, j] <- M[, j] - stats::ave(M[, j], g)
      M
    }
    suppressWarnings(stats::cor(resid_within(Y), resid_within(scores)))
  } else {
    suppressWarnings(stats::cor(Y, scores))
  }

  list(
    scores = scores,
    group = if (is.null(group)) NULL else droplevels(as.factor(group)),
    canonical = data.frame(
      axis = colnames(scores),
      eigenvalue = vals,
      canonical_r = sqrt(vals / (1 + vals)),
      prop_variance = vals / sum(vals),
      stringsAsFactors = FALSE
    ),
    structure = as.data.frame(struct)
  )
}

#' Detect complete or quasi-complete separation in a binomial GLM
#'
#' Wald statistics collapse under separation: the odds ratio explodes, the
#' interval covers everything, and the p-value approaches 1. \code{glm()}
#' reports \code{converged = TRUE} regardless, so this has to be checked
#' explicitly rather than inferred from the fit.
#'
#' The test is on the fitted model, not on the design. A homogeneous cell of a
#' crossed grouping is \emph{not} evidence of separation when the model is
#' additive: that cell's linear predictor is estimable from the others, and a
#' cell-based rule flags close to half of all clean additive fits. Separation is
#' therefore identified from the coefficients themselves — one that has run away
#' on the log-odds scale together with a standard error that has run away with
#' it — and from fitted probabilities that have reached the numerical boundary.
#' The offending coefficients are named, since it is the coefficients, not the
#' cells, that are unidentified.
#'
#' @param fit a fitted binomial GLM.
#' @return a character note, or \code{character(0)}.
#' @noRd
.check_separation <- function(fit) {
  co <- stats::coef(fit)
  se <- suppressWarnings(sqrt(diag(stats::vcov(fit))))
  keep <- intersect(names(co), names(se))
  co <- co[keep]; se <- se[keep]

  runaway <- names(co)[is.finite(co) & is.finite(se) &
                         abs(co) > 10 & se > 10]
  mu <- stats::fitted(fit)
  boundary <- sum(mu < 1e-10 | mu > 1 - 1e-10, na.rm = TRUE)

  if (length(runaway) == 0L && boundary == 0L) return(character(0))

  detail <- if (length(runaway) > 0L) {
    sprintf(" The affected coefficient(s): %s.",
            paste(setdiff(runaway, "(Intercept)"), collapse = ", "))
  } else {
    sprintf(" %d fitted probabilit%s numerically 0 or 1.",
            boundary, if (boundary == 1L) "y is" else "ies are")
  }
  paste0(
    "Complete or quasi-complete separation detected.", detail,
    " An affected odds ratio is not identified: it will be enormous, its ",
    "interval will be unbounded on one side, and its Wald p-value will be near 1 ",
    "no matter how strong the association is. Consider a penalised fit such as ",
    "logistf::logistf(), or collapsing the offending level.")
}
