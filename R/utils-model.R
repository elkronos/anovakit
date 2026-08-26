# Model fitting helpers ------------------------------------------------------

#' Fit a model under the contrasts the requested ANOVA type needs
#'
#' Type III sums of squares are only meaningful when the factors are coded with
#' contrasts that sum to zero. Setting \code{options(contrasts = ...)} after a
#' model has been fitted does nothing: \code{car::Anova()} reads the contrasts
#' stored on the object, so the model has to be fitted under them.
#'
#' @param fit_fun a function of no arguments that fits and returns the model.
#' @param type \code{"I"}, \code{"II"} or \code{"III"}.
#' @return the fitted model.
#' @noRd
.fit_with_contrasts <- function(fit_fun, type = "II") {
  if (identical(type, "III")) {
    old <- options(contrasts = c("contr.sum", "contr.poly"))
    on.exit(options(old), add = TRUE)
  }
  fit_fun()
}

#' Run car::Anova with a validated type and test statistic
#'
#' @param model a fitted model, already fitted under the right contrasts.
#' @param type \code{"II"} or \code{"III"}.
#' @param test_statistic passed to \code{car::Anova}.
#' @return list with \code{table} (a tidied data.frame), \code{raw} (the
#'   \code{car} object) and \code{note}.
#' @noRd
.car_anova <- function(model, type = "II", test_statistic = NULL, ...) {
  args <- list(mod = model, type = if (type == "III") 3L else 2L)
  if (!is.null(test_statistic)) args$test.statistic <- test_statistic
  args <- c(args, list(...))
  # car emits "Note: model has aliased coefficients" through message(), which
  # suppressWarnings() does not catch. The package promises not to write to the
  # console, so anything it says is collected into $notes instead.
  got <- .collect_conditions(do.call(car::Anova, args))
  raw <- got$value
  if (inherits(raw, "error")) {
    return(list(table = NULL, raw = NULL, note = sprintf(
      "The Type %s ANOVA table could not be computed: %s",
      type, conditionMessage(raw))))
  }
  list(table = .tidy_anova(raw, type), raw = raw,
       note = .said_note(got$said, "car::Anova()"))
}

#' Evaluate an expression, capturing anything it says instead of printing it
#'
#' Third-party warnings and messages are informative but they are console
#' output, and the package guarantees it produces none. They are captured here
#' and surfaced in \code{$notes}, where the rest of what the function decided
#' already lives.
#'
#' @return \code{list(value, said)}; \code{value} is the result, or the
#'   condition object if the expression errored.
#' @noRd
.collect_conditions <- function(expr) {
  said <- character(0)
  value <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = function(w) {
      said <<- c(said, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      said <<- c(said, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  # Condition messages may span several lines. $notes is printed one bullet per
  # element, so flatten them to a single line.
  said <- trimws(gsub("[[:space:]]+", " ", said))
  list(value = value, said = unique(said[nzchar(said)]))
}

#' Turn captured conditions into a single note
#'
#' Routine chatter is dropped. A note is for something the caller might act on;
#' a package announcing a default it applied -- which the wrapper set
#' deliberately, and documents -- is not that, and printing it on every call
#' teaches people to stop reading \code{$notes}.
#' @noRd
.said_note <- function(said, who) {
  boring <- paste(c(
    "^Contrasts set to contr\\.sum",   # afex, on every between-subjects fit
    "^Setting contrasts",              # afex, older versions
    "^NOTE: Results may be misleading due to involvement in interactions",
    "^Note: D\\.f\\. calculations"
  ), collapse = "|")
  said <- said[!grepl(boring, said)]
  if (length(said) == 0L) return(character(0))
  sprintf("%s reported: %s", who, paste(said, collapse = " / "))
}

#' Turn a car::Anova or stats::anova table into a consistent data frame
#' @noRd
.tidy_anova <- function(tab, type = NA_character_) {
  df <- as.data.frame(tab)
  out <- data.frame(term = rownames(df), stringsAsFactors = FALSE)
  ren <- c("Sum Sq" = "sum_sq", "Df" = "df", "F value" = "statistic",
           "F values" = "statistic", "LR Chisq" = "statistic",
           "Chisq" = "statistic", "Wald Chisq" = "statistic",
           "Deviance" = "deviance", "Resid. Df" = "resid_df",
           "Resid. Dev" = "resid_dev", "Mean Sq" = "mean_sq")
  for (nm in names(df)) {
    target <- if (nm %in% names(ren)) unname(ren[nm]) else
      if (grepl("^Pr\\(", nm)) "p_value" else
        gsub("[^A-Za-z0-9]+", "_", tolower(nm))
    if (!target %in% names(out)) out[[target]] <- df[[nm]]
  }
  out <- out[!out$term %in% "(Intercept)", , drop = FALSE]
  row.names(out) <- NULL
  if (!is.na(type)) attr(out, "ss_type") <- type
  out
}

#' Build the right-hand side of a model formula from grouping variables
#'
#' @param groups character vector of grouping columns.
#' @param interaction \code{FALSE} for an additive model, \code{TRUE} for the
#'   full factorial, or an integer giving the highest order of interaction to
#'   include.
#' @return character vector of term labels, already backtick-quoted.
#' @noRd
.group_terms <- function(groups, interaction = FALSE) {
  q <- .bq(groups)
  if (length(q) == 1L) return(q)
  if (isFALSE(interaction)) return(q)
  order_max <- if (isTRUE(interaction)) length(q) else as.integer(interaction)
  if (!is.finite(order_max) || order_max < 1L) {
    .stopf("`interaction` must be TRUE, FALSE, or a whole number of at least 1.")
  }
  order_max <- min(order_max, length(q))
  if (order_max == 1L) return(q)
  terms <- unlist(lapply(seq_len(order_max), function(k) {
    utils::combn(q, k, FUN = function(z) paste(z, collapse = ":"))
  }), use.names = FALSE)
  terms
}

#' Validate the interaction argument
#' @noRd
.check_interaction <- function(interaction, n_groups) {
  if (is.logical(interaction) && length(interaction) == 1L && !is.na(interaction)) {
    return(invisible(interaction))
  }
  if (is.numeric(interaction) && length(interaction) == 1L &&
      is.finite(interaction) && interaction >= 1 &&
      interaction == round(interaction)) {
    if (interaction > n_groups) {
      .stopf("`interaction` = %d exceeds the number of grouping variables (%d).",
             as.integer(interaction), n_groups)
    }
    return(invisible(interaction))
  }
  .stopf("`interaction` must be TRUE, FALSE, or a whole number between 1 and %d.",
         n_groups)
}

#' Warn when a factorial model is too rich for the data
#' @noRd
.check_model_size <- function(fit) {
  dfr <- stats::df.residual(fit)
  np  <- length(stats::coef(fit))
  notes <- character(0)
  if (!is.null(dfr) && is.finite(dfr) && dfr < 1) {
    notes <- c(notes, sprintf(
      "The model has %d parameters and no residual degrees of freedom: it is saturated, so there is nothing left to test against. Every p-value, interval and effect size below is undefined or degenerate (partial eta squared is exactly 1 by construction). Fewer groups, a simpler interaction structure, or more data are needed. Zero prior weights count against the sample size here.",
      np))
  } else if (!is.null(dfr) && is.finite(dfr) && dfr < np) {
    notes <- c(notes, sprintf(
      "The model has %d parameters and only %d residual degrees of freedom, so it is close to saturated. Consider fewer groups, a simpler interaction structure, or more data.",
      np, dfr))
  }
  if (any(is.na(stats::coef(fit)))) {
    bad <- names(stats::coef(fit))[is.na(stats::coef(fit))]
    notes <- c(notes, sprintf(
      "%d coefficient(s) could not be estimated because the corresponding cells are empty: %s.",
      length(bad), paste(utils::head(bad, 5L), collapse = ", ")))
  }
  notes
}

#' Resolve an optional weights column into a vector
#'
#' Weights are named as a column rather than passed through \code{...}. An
#' expression evaluated in the caller is evaluated against the \emph{original}
#' frame, so as soon as one row is dropped for missing values the weights no
#' longer line up with the model frame and \code{glm()} errors on the length
#' mismatch -- or, worse, does not.
#'
#' @param data the prepared (complete-case) frame.
#' @param weights a single column name, or \code{NULL}.
#' @noRd
.resolve_weights <- function(data, weights) {
  if (is.null(weights)) return(NULL)
  .check_weights_column(data, weights)
  w <- data[[weights]]
  # Non-finite weights are dropped with their rows by .prepare_frame(), so by
  # the time we get here every weight is finite and non-negative. What can
  # still be wrong is that nothing is left to fit on: glm.fit() with no
  # positive weight fails with "object 'fit' not found", which names nothing
  # the caller can act on.
  pos <- sum(w > 0)
  if (pos == 0L) {
    .stopf("Every weight in `%s` is zero, so no observation would contribute to the fit.",
           weights)
  }
  if (pos < 2L) {
    .stopf("Only %d observation has a non-zero weight in `%s`; a model cannot be fitted on it.",
           pos, weights)
  }
  w
}

#' Validate a weights column on the data as supplied
#'
#' Run before rows are dropped, so that an unusable weight is reported rather
#' than silently removed along with a row that was missing something else.
#' Non-finite weights are deliberately \emph{not} an error here: they are
#' treated as missing, dropped by \code{.prepare_frame()} and counted in
#' \code{$n_removed}, which is what happens to a non-finite value in any other
#' analysed column.
#' @noRd
.check_weights_column <- function(data, weights, response = NULL,
                                  groups = NULL, offset = NULL) {
  if (is.null(weights)) return(invisible(NULL))
  .check_name(weights, "weights")
  .check_columns(data, weights, "Weights column")
  # A column cannot be both a predictor and the prior weights. Catching it here
  # by name gives a message the caller can act on; left to the fitting call it
  # surfaces as a type error about a factor conversion the function itself
  # performed, or -- for the response -- not at all.
  role <- c(response = "the response", groups = "a grouping variable",
            offset = "the exposure offset")
  for (nm in names(role)) {
    if (weights %in% get(nm)) {
      .stopf("Weights column `%s` is already %s; a column cannot be both.",
             weights, role[[nm]])
    }
  }
  w <- data[[weights]]
  if (!is.numeric(w)) {
    .stopf("Weights column `%s` must be numeric; it is %s.",
           weights, paste(class(w), collapse = "/"))
  }
  neg <- which(is.finite(w) & w < 0)
  if (length(neg) > 0L) {
    .stopf("Weights column `%s` contains %d negative value(s) (row%s %s); weights must be non-negative.",
           weights, length(neg), if (length(neg) == 1L) "" else "s",
           .abbrev(neg))
  }
  invisible(TRUE)
}

#' Does this model estimate its dispersion parameter?
#'
#' \code{summary.glm()} uses a t reference distribution whenever the dispersion
#' is estimated (Gaussian, Gamma, inverse Gaussian, every quasi family) and a
#' normal one only when it is fixed at 1 (binomial, Poisson). Using the normal
#' everywhere makes the p-values of a Gamma fit wrong by orders of magnitude.
#' @noRd
.estimates_dispersion <- function(model) {
  if (!inherits(model, "glm")) return(TRUE)              # lm, and lm-alikes
  fam <- tryCatch(stats::family(model)$family, error = function(e) NA_character_)
  if (is.na(fam)) return(TRUE)
  # MASS::glm.nb() estimates theta but then fixes the dispersion at 1, and
  # summary.negbin() reports a z value. Match on the CLASS, not the family
  # string: glm(family = MASS::negative.binomial(2)) produces the same
  # "Negative Binomial(2)" string but is an ordinary glm whose dispersion
  # summary.glm() estimates, and which therefore takes a t reference.
  if (inherits(model, "negbin")) return(FALSE)
  !(fam %in% c("binomial", "poisson"))
}

#' The reference distribution's critical value for a coefficient interval
#' @noRd
.coef_crit <- function(model, conf_level) {
  p <- 1 - (1 - conf_level) / 2
  if (!.estimates_dispersion(model)) return(stats::qnorm(p))
  dfr <- stats::df.residual(model)
  # A saturated model has no residual degrees of freedom, and qt(df = 0) warns
  # "NaNs produced" straight past every suppressor in the package. There is no
  # inference to be had there; return NA and let the table say so.
  if (is.null(dfr) || !is.finite(dfr) || dfr < 1) return(NA_real_)
  stats::qt(p, df = dfr)
}

#' Coefficient table with model-based or robust standard errors
#'
#' Computed directly from the coefficients and the covariance matrix, so no
#' external package and no defensive column matching are involved.
#'
#' @noRd
.coef_table <- function(model, vcov_matrix = NULL) {
  est <- stats::coef(model)
  V <- if (is.null(vcov_matrix)) stats::vcov(model) else vcov_matrix
  se <- sqrt(diag(V))
  keep <- intersect(names(est), colnames(V))
  est <- est[keep]; se <- se[keep]
  stat <- est / se
  use_t <- .estimates_dispersion(model)
  p <- if (use_t) {
    2 * stats::pt(-abs(stat), df = stats::df.residual(model))
  } else {
    2 * stats::pnorm(-abs(stat))
  }
  data.frame(term = names(est), estimate = unname(est), se = unname(se),
             statistic = unname(stat), p_value = unname(p),
             distribution = if (use_t) "t" else "z",
             stringsAsFactors = FALSE, row.names = NULL)
}

#' Fit a GLM with prior weights that cannot be shadowed by a data column
#'
#' \code{glm()} looks \code{weights} up in \code{data} first, so a column
#' called \code{wts} would capture a local variable of that name. Adding the
#' weights to the frame under a generated name and referring to that name
#' removes the ambiguity entirely.
#'
#' @param fml model formula.
#' @param data the prepared frame.
#' @param family a family object, or \code{NULL} for \code{MASS::glm.nb}.
#' @param wts numeric weights, or \code{NULL}.
#' @noRd
.fit_glm <- function(fml, data, family = NULL, wts = NULL) {
  if (is.null(wts)) {
    if (is.null(family)) return(MASS::glm.nb(fml, data = data))
    return(stats::glm(fml, data = data, family = family))
  }
  wt_name <- .safe_name(data, ".prior_weights")
  data[[wt_name]] <- wts
  call <- if (is.null(family)) {
    bquote(MASS::glm.nb(.(fml), data = data, weights = .(as.name(wt_name))))
  } else {
    bquote(stats::glm(.(fml), data = data, family = family,
                      weights = .(as.name(wt_name))))
  }
  eval(call)
}
