# The shared return class ----------------------------------------------------

#' The object returned by every anovakit function
#'
#' All eight analysis functions return an object of class
#' \code{anovakit_fit}. It is a plain list, so \code{$} extraction works as
#' usual, with \code{print}, \code{summary} and \code{plot} methods for
#' convenience.
#'
#' @section Components:
#' \describe{
#'   \item{method}{Character. The analysis that was run, as printed.}
#'   \item{call}{The matched call, so a result can always say how it was made.}
#'   \item{model}{The fitted model object: an \code{lm} from
#'     \code{\link{anova_ancova}}, a \code{glm} from \code{\link{anova_bin}},
#'     \code{\link{anova_count}} and \code{\link{anova_glm}}, a
#'     \code{manova} from \code{\link{anova_manova}}, an \code{afex_aov} from
#'     \code{\link{anova_rm}}, and the \code{htest} returned by
#'     \code{\link[stats]{oneway.test}} or \code{\link[stats]{kruskal.test}}
#'     from \code{\link{anova_welch}} and \code{\link{anova_kw}}, which fit no
#'     model.}
#'   \item{anova}{Data frame. The omnibus test table.}
#'   \item{effect_sizes}{Data frame, or \code{NULL}. What it holds depends on
#'     the method: partial eta squared with partial omega squared
#'     (\code{anova_ancova}, \code{anova_manova}, and \code{anova_glm} on a
#'     Gaussian family), Hedges' g with intervals (\code{anova_welch}), epsilon
#'     squared and eta squared for the rank statistic (\code{anova_kw}), partial
#'     and generalised eta squared (\code{anova_rm}), odds ratios
#'     (\code{anova_bin}), incidence rate ratios (\code{anova_count}), and the
#'     deviance explained with McFadden's pseudo R squared (\code{anova_glm} on
#'     any other family). Both variance measures are \emph{partial}: the
#'     classical omega squared is only a proportion of variance when the effect
#'     sums of squares partition the total, which Type II and Type III sums of
#'     squares do not.}
#'   \item{emmeans}{Data frame, or \code{NULL}. Estimated marginal means with
#'     intervals at \code{conf_level}.}
#'   \item{emmeans_object}{The \code{emmGrid}, or \code{NULL}. Pass this to
#'     \pkg{emmeans} for contrasts the wrapper does not cover. It carries
#'     \pkg{emmeans}' own default confidence level rather than
#'     \code{conf_level}, so give \code{level =} when you summarise it. It is
#'     \code{NULL} for \code{\link{anova_welch}} and \code{\link{anova_kw}},
#'     which fit no model \pkg{emmeans} can use, and for
#'     \code{\link{anova_manova}}, where there is one grid per response under
#'     \code{$univariate[[response]]$emmeans_object}.}
#'   \item{posthoc}{Data frame, or \code{NULL}. Pairwise comparisons.}
#'   \item{assumptions}{Named list of assumption checks. Contents vary by
#'     method; each element is a test object, a data frame or \code{NULL}. It is
#'     empty for \code{\link{anova_kw}} unless \code{diagnostics = TRUE}, since
#'     the test assumes no distribution.}
#'   \item{plots}{Named list of \pkg{ggplot2} objects. Empty when
#'     \code{plots = FALSE}. Nothing is ever drawn as a side effect.}
#'   \item{data_used}{Data frame. The rows and columns the model was fitted on.}
#'   \item{n_removed}{Integer. Input rows that are not in \code{data_used}:
#'     dropped for missing or infinite values or, in \code{\link{anova_rm}},
#'     for an incomplete within-subject design or because \pkg{afex}
#'     aggregated duplicated subject-by-cell rows into their means.
#'     \code{nrow(data_used) + n_removed} is always the number of rows given.}
#'   \item{conf_level}{Numeric. The level used for every interval returned.}
#'   \item{notes}{Character vector. Everything the function decided on your
#'     behalf, could not compute, or thinks you should know -- including
#'     anything \pkg{car}, \code{\link[stats]{glm}} or \pkg{afex} said while
#'     the model was being fitted. Always read this.}
#' }
#'
#' @section Components individual functions add:
#' Each function returns everything above plus whatever its own method
#' produces. \code{names(fit)} lists them all. \code{\link{anova_ancova}} adds
#' \code{$slopes_test}, \code{$simple_slopes} and \code{$covariate_means};
#' \code{\link{anova_rm}} adds \code{$sphericity}, \code{$subjects_dropped} and
#' the breakdown of \code{$n_removed}; \code{\link{anova_manova}} adds
#' \code{$multivariate}, \code{$univariate}, \code{$canonical} and
#' \code{$canonical_term}; \code{\link{anova_count}} adds \code{$model_type},
#' \code{$dispersion} and \code{$model_dispersion};
#' \code{\link{anova_bin}} and \code{\link{anova_glm}} add
#' \code{$model_stats}. Each is documented on the function that produces it.
#'
#' @name anovakit_fit
#' @seealso \code{\link{print.anovakit_fit}},
#'   \code{\link{summary.anovakit_fit}}, \code{\link{plot.anovakit_fit}}
NULL

#' Build an anovakit_fit object
#' @noRd
.new_fit <- function(method, call, model = NULL, anova = NULL,
                     effect_sizes = NULL, emmeans = NULL, emmeans_object = NULL,
                     posthoc = NULL, assumptions = list(), plots = list(),
                     data_used = NULL, n_removed = 0L, conf_level = 0.95,
                     notes = character(0), extra = list()) {
  out <- c(
    list(
      method         = method,
      call           = call,
      model          = model,
      anova          = anova,
      effect_sizes   = effect_sizes,
      emmeans        = emmeans,
      emmeans_object = emmeans_object,
      posthoc        = posthoc,
      assumptions    = assumptions,
      plots          = plots,
      data_used      = data_used,
      n_removed      = as.integer(n_removed),
      conf_level     = conf_level,
      notes          = unique(notes)
    ),
    extra
  )
  class(out) <- "anovakit_fit"
  out
}

#' Print an anovakit result
#'
#' Shows the method, the sample size, the omnibus table and any notes. Use
#' \code{\link{summary.anovakit_fit}} for effect sizes, assumption checks
#' and post-hoc comparisons.
#'
#' @param x An \code{\link{anovakit_fit}} object.
#' @param digits Number of significant digits for the printed tables. Default
#'   \code{4}. Unlike most printing in R this is not capped by
#'   \code{getOption("digits")}.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @examples
#' set.seed(1)
#' d <- data.frame(g = rep(c("a", "b", "c"), each = 20), y = rnorm(60))
#' fit <- anova_welch(d, "y", "g", plots = FALSE)
#' print(fit)
#' @export
print.anovakit_fit <- function(x, digits = 4L, ...) {
  digits <- .check_digits(digits)
  cat(x$method, "\n")
  cat(strrep("-", nchar(x$method)), "\n", sep = "")
  if (!is.null(x$call)) {
    cat("Call: ", paste(deparse(x$call), collapse = " "), "\n", sep = "")
  }
  n <- if (!is.null(x$data_used)) nrow(x$data_used) else NA_integer_
  cat(sprintf("Observations used: %s%s\n",
              if (is.na(n)) "unknown" else format(n),
              if (x$n_removed > 0L)
                sprintf("  (%d input row(s) not in $data_used; see $notes)",
                        x$n_removed)
              else ""))
  if (!is.null(x$anova) && NROW(x$anova) > 0L) {
    cat("\n", .omnibus_header(x$anova), "\n", sep = "")
    .print_df(x$anova, digits)
  }
  if (length(x$notes) > 0L) {
    cat("\nNotes\n")
    cat(paste0("  - ", x$notes, collapse = "\n"), "\n", sep = "")
  }
  if (length(x$plots) > 0L) {
    cat("\nPlots available: ", paste(names(x$plots), collapse = ", "),
        "\n  (use plot(x, which = \"", names(x$plots)[1L], "\"))\n", sep = "")
  }
  invisible(x)
}

#' The heading of the omnibus table, naming what it contains
#'
#' The statistic (F, likelihood-ratio or Wald chi-square), the type of sums of
#' squares and the sphericity correction all change what the numbers mean, and
#' none of them can be read off the table itself.
#' @noRd
.omnibus_header <- function(tab) {
  bits <- c(
    if (!is.null(attr(tab, "ss_type"))) sprintf("Type %s", attr(tab, "ss_type")),
    attr(tab, "statistic"),
    if (!is.null(attr(tab, "correction"))) {
      corr <- attr(tab, "correction")
      if (identical(corr, "none")) "no sphericity correction" else
        sprintf("%s-corrected degrees of freedom", corr)
    })
  if (length(bits) == 0L) "Omnibus test" else
    sprintf("Omnibus test (%s)", paste(bits, collapse = ", "))
}

#' Summarise an anovakit result
#'
#' Everything \code{print} shows, plus assumption checks, effect sizes,
#' marginal means, post-hoc comparisons and whatever further tables the method
#' produces (simple slopes, sphericity, the multivariate tests and univariate
#' follow-ups, the canonical axes).
#'
#' @param object An \code{\link{anovakit_fit}} object.
#' @param digits Number of significant digits for the printed tables. Default
#'   \code{4}.
#' @param ... Ignored.
#' @return An object of class \code{summary.anovakit_fit}, which prints the
#'   summary. Its \code{$fit} element is \code{object}.
#' @examples
#' set.seed(1)
#' d <- data.frame(g = rep(c("a", "b", "c"), each = 20), y = rnorm(60))
#' summary(anova_welch(d, "y", "g", plots = FALSE))
#' @export
summary.anovakit_fit <- function(object, digits = 4L, ...) {
  structure(list(fit = object, digits = .check_digits(digits)),
            class = "summary.anovakit_fit")
}

#' @rdname summary.anovakit_fit
#' @param x A \code{summary.anovakit_fit} object.
#' @export
print.summary.anovakit_fit <- function(x, digits = x$digits, ...) {
  object <- x$fit
  digits <- .check_digits(digits)
  print(object, digits = digits)

  if (length(object$assumptions) > 0L) {
    cat("\nAssumption checks\n")
    for (nm in names(object$assumptions)) {
      el <- object$assumptions[[nm]]
      if (is.null(el)) next
      cat("\n  ", nm, "\n", sep = "")
      if (inherits(el, "htest")) {
        cat(sprintf("    %s: statistic = %s, p = %s\n", el$method,
                    format(signif(unname(el$statistic), digits)),
                    format.pval(unname(el$p.value), digits = digits,
                                eps = .Machine$double.eps)))
      } else if (is.data.frame(el)) {
        .print_df(el, digits, row.names = FALSE)
      } else if (is.numeric(el) && length(el) == 1L) {
        cat("    ", format(signif(el, digits)), "\n", sep = "")
      } else {
        print(el, digits = digits)
      }
    }
  }
  if (!is.null(object$effect_sizes) && NROW(object$effect_sizes) > 0L) {
    cat("\nEffect sizes\n")
    .print_df(object$effect_sizes, digits, row.names = FALSE)
  }
  if (!is.null(object$emmeans) && NROW(object$emmeans) > 0L) {
    header <- if (is.null(object$emmeans_object)) "Group summaries" else
      "Estimated marginal means"
    cat(sprintf("\n%s (%s%% intervals)\n", header, .pct(object$conf_level)))
    .print_df(object$emmeans, digits, row.names = FALSE)
  }
  if (!is.null(object$posthoc) && NROW(object$posthoc) > 0L) {
    cat("\nPairwise comparisons\n")
    .print_df(object$posthoc, digits, row.names = FALSE)
  }
  extra <- c(sphericity = "Sphericity", slopes_test = "Homogeneity of slopes",
             simple_slopes = "Covariate slopes by group",
             multivariate = "Multivariate tests", canonical = "Canonical axes")
  for (nm in names(extra)) {
    el <- object[[nm]]
    if (is.data.frame(el) && NROW(el) > 0L) {
      cat("\n", extra[[nm]], "\n", sep = "")
      .print_df(el, digits, row.names = FALSE)
    } else if (is.list(el) && is.data.frame(el$table) && NROW(el$table) > 0L) {
      cat("\n", extra[[nm]], "\n", sep = "")
      .print_df(el$table, digits, row.names = FALSE)
    }
  }
  if (is.list(object$univariate) && length(object$univariate) > 0L) {
    for (r in names(object$univariate)) {
      tab <- object$univariate[[r]]$anova
      if (is.data.frame(tab) && NROW(tab) > 0L) {
        cat("\nUnivariate follow-up: ", r, "\n", sep = "")
        .print_df(tab, digits, row.names = FALSE)
      }
    }
  }
  invisible(x)
}

#' Plot an anovakit result
#'
#' Returns one of the \pkg{ggplot2} objects the analysis built. The available
#' names are listed by \code{print()} and stored in \code{x$plots}.
#'
#' @param x An \code{\link{anovakit_fit}} object.
#' @param which Name or index of the plot. Defaults to the first one.
#' @param ... Ignored.
#' @return A \pkg{ggplot2} object, invisibly returning \code{NULL} when the fit
#'   holds no plots.
#' @examples
#' set.seed(1)
#' d <- data.frame(g = rep(c("a", "b", "c"), each = 20), y = rnorm(60))
#' fit <- anova_welch(d, "y", "g")
#' p <- plot(fit, which = "means")
#' @export
plot.anovakit_fit <- function(x, which = 1L, ...) {
  if (length(x$plots) == 0L) {
    message("This fit holds no plots (it was called with plots = FALSE).")
    return(invisible(NULL))
  }
  if (is.null(which) || length(which) != 1L || is.na(which)) {
    .stopf("`which` must be a single plot name or index. Available: %s.",
           paste(names(x$plots), collapse = ", "))
  }
  if (is.character(which)) {
    if (!which %in% names(x$plots)) {
      .stopf("No plot called \"%s\". Available: %s.",
             which, paste(names(x$plots), collapse = ", "))
    }
    return(x$plots[[which]])
  }
  if (!is.numeric(which) || which != round(which) ||
      which < 1L || which > length(x$plots)) {
    .stopf("`which` must be a whole number between 1 and %d, or one of: %s.",
           length(x$plots), paste(names(x$plots), collapse = ", "))
  }
  x$plots[[as.integer(which)]]
}

#' Validate a digits argument
#' @noRd
.check_digits <- function(digits) {
  if (!is.numeric(digits) || length(digits) != 1L || !is.finite(digits)) {
    .stopf("`digits` must be a single whole number between 1 and 22.")
  }
  as.integer(min(max(round(digits), 1L), 22L))
}

#' A confidence level as a percentage, without spurious rounding
#' @noRd
.pct <- function(conf_level) {
  format(conf_level * 100, digits = 8L, drop0trailing = TRUE, trim = TRUE)
}

#' Round the numeric columns of a data frame for printing
#'
#' Counts and degrees of freedom are whole numbers; rounding them to
#' significant digits would print n = 12345 as 12340, so whole-number columns
#' are left alone. p-value columns are formatted with \code{format.pval()},
#' which prints a p-value below machine precision as "< 2.2e-16" rather than as
#' zero.
#' @noRd
.round_df <- function(df, digits = 4L) {
  df <- as.data.frame(df)
  for (nm in names(df)) {
    v <- df[[nm]]
    if (!is.numeric(v)) next
    if (.is_p_column(nm)) {
      df[[nm]] <- ifelse(is.na(v), NA_character_,
                         format.pval(v, digits = digits, eps = .Machine$double.eps))
      next
    }
    if (is.integer(v) || all(is.na(v) | !is.finite(v) | v == round(v))) next
    df[[nm]] <- signif(v, digits)
  }
  df
}

#' Is this column a p-value?
#' @noRd
.is_p_column <- function(nm) {
  grepl("^p_value$|^p_adjusted$|^p_(gg|hf|value_.*)$|^p$", nm)
}

#' Print a table at the requested number of digits
#'
#' \code{print.data.frame} re-formats with \code{getOption("digits")}, so
#' rounding alone lets a small global setting override the argument the caller
#' asked for. Setting the option for the duration of the call fixes that.
#' @noRd
.print_df <- function(df, digits = 4L, ...) {
  digits <- .check_digits(digits)
  old <- options(digits = digits)
  on.exit(options(old), add = TRUE)
  print(.round_df(df, digits), ...)
  invisible(NULL)
}
