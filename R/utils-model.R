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
.car_anova <- function(model, type = "II", test_statistic = NULL,
                       vcov_matrix = NULL, ...) {
  # car::Anova() refuses a linear model whose residual sum of squares is below
  # an ABSOLUTE tolerance (sqrt(.Machine$double.eps)), so a response measured in
  # small units -- metres instead of millimetres -- loses its whole table to a
  # false "residual sum of squares is 0". The test statistics do not depend on
  # the units, so the model is rescaled for car and the sums of squares are
  # rescaled back. A genuinely perfect fit (RSS tiny relative to the total) is
  # left alone and still refused.
  k <- .lm_rescale_factor(model)
  mod <- if (k != 1) .rescale_lm(model, k) else model
  args <- list(mod = mod, type = if (type == "III") 3L else 2L)
  if (!is.null(test_statistic)) args$test.statistic <- test_statistic
  # A covariance for the coefficients rescales with them (by k^2), so Wald
  # tests built from it are unchanged by the rescaling too.
  if (!is.null(vcov_matrix)) args$vcov. <- vcov_matrix * k^2
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
  if (k != 1 && "Sum Sq" %in% names(raw)) raw[["Sum Sq"]] <- raw[["Sum Sq"]] / k^2
  said <- got$said
  if (!is.null(vcov_matrix)) {
    # With a user-supplied covariance car's lm method names its F column "F"
    # and announces the covariance by deparsing it into a message.
    if ("F" %in% names(raw)) names(raw)[names(raw) == "F"] <- "F value"
    said <- said[!grepl("^Coefficient covariances computed by", said)]
  }
  list(table = .tidy_anova(raw, type), raw = raw,
       note = .said_note(said, "car::Anova()"))
}

#' Say whether the omnibus test uses a robust covariance
#'
#' Likelihood-ratio and F tests compare deviances or residual sums of squares
#' and cannot use a sandwich covariance; a Wald test can. With robust standard
#' errors requested the reader needs to know which kind the table holds.
#' @noRd
.omnibus_vcov_note <- function(test_statistic, vcov_type, robust_used) {
  if (robust_used) {
    return(sprintf("The omnibus Wald test uses the robust (%s) covariance, as the coefficient table, the marginal means and the comparisons do.",
                   vcov_type))
  }
  sprintf("The omnibus %s test compares deviances, so it rests on the model-based variance; the robust (%s) covariance is used only by the coefficient table, the marginal means and the comparisons. For an omnibus test that uses it, set test_statistic = \"Wald\".",
          if (identical(test_statistic, "F")) "F" else "likelihood-ratio",
          vcov_type)
}

#' The factor that brings a small-unit linear model above car's tolerance
#' @noRd
.lm_rescale_factor <- function(model) {
  if (!inherits(model, "lm") || inherits(model, "glm") || inherits(model, "mlm")) {
    return(1)
  }
  r <- stats::residuals(model)
  rss <- sum(r^2)
  y <- stats::fitted(model) + r
  tss <- sum((y - mean(y))^2)
  if (!is.finite(rss) || !is.finite(tss) || rss <= 0 || tss <= 0) return(1)
  if (rss >= 1e-4 || rss / tss < 1e-12) return(1)
  1 / stats::sd(y)
}

#' Multiply the response of a fitted linear model by a constant
#'
#' Every quantity car::Anova() reads from an lm -- coefficients, residuals,
#' fitted values, effects and the response in the model frame -- scales by k,
#' so its F statistics and p-values are unchanged and its sums of squares
#' scale by k^2.
#' @noRd
.rescale_lm <- function(model, k) {
  model$coefficients <- model$coefficients * k
  model$residuals <- model$residuals * k
  model$fitted.values <- model$fitted.values * k
  if (!is.null(model$effects)) model$effects <- model$effects * k
  if (!is.null(model$model)) model$model[[1L]] <- model$model[[1L]] * k
  model
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
  # Non-syntactic column names reach the model formula backtick-quoted, and the
  # term labels keep the quotes. The package forbids backticks inside column
  # names (see .bq()), so stripping them is lossless.
  out <- data.frame(term = gsub("`", "", rownames(df), fixed = TRUE),
                    stringsAsFactors = FALSE)
  stat_name <- intersect(c("F value", "F values", "LR Chisq", "Wald Chisq", "Chisq"),
                         names(df))[1L]
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
  if (!is.na(stat_name)) {
    attr(out, "statistic") <- switch(stat_name,
      "F value" = , "F values" = "F", "LR Chisq" = "likelihood-ratio chi-square",
      "Wald Chisq" = , "Chisq" = "Wald chi-square")
  }
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

#' Number of model coefficients the grouping terms need
#'
#' Counted from the levels, without building the model matrix: every term
#' contributes the product of its factors' (levels - 1), whatever the
#' contrasts. Three 22-level factors in a full factorial need 10647 on top of
#' the intercept, and \code{lm()} spends minutes pivoting a model matrix that
#' wide however few rows there are.
#' @noRd
.n_group_coefs <- function(d, groups, interaction = FALSE) {
  k <- vapply(groups, function(g) nlevels(droplevels(as.factor(d[[g]]))) - 1,
              numeric(1))
  order_max <- if (length(groups) == 1L || isFALSE(interaction)) 1L else
    if (isTRUE(interaction)) length(groups) else
      min(as.integer(interaction), length(groups))
  sum(vapply(seq_len(order_max), function(m) {
    # combn() on the indices: given a single number it would count up to it.
    sum(utils::combn(seq_along(k), m, FUN = function(i) prod(k[i])))
  }, numeric(1)))
}

#' Does an interaction argument give a model with no interaction terms?
#' @noRd
.is_additive <- function(interaction) {
  isFALSE(interaction) || (is.numeric(interaction) && interaction == 1)
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
#'
#' What a saturated model loses depends on the family. When the dispersion is
#' estimated (Gaussian, Gamma, quasi families) there is no residual variance to
#' test against, so nothing is testable. When it is fixed (binomial, Poisson,
#' negative binomial) the likelihood-ratio and Wald tests remain valid -- a
#' saturated log-linear model is how a G-test of independence is computed --
#' and only goodness of fit and overdispersion become uncheckable.
#'
#' Aliased coefficients are reported as aliased, without guessing why: empty
#' cells of a crossed design and collinear covariates both produce them.
#' @noRd
.check_model_size <- function(fit) {
  dfr <- stats::df.residual(fit)
  np  <- sum(!is.na(stats::coef(fit)))
  fixed <- !.estimates_dispersion(fit)
  gaussian_lm <- inherits(fit, "lm") &&
    (!inherits(fit, "glm") || identical(stats::family(fit)$family, "gaussian"))
  notes <- character(0)
  if (!is.null(dfr) && is.finite(dfr) && dfr < 1) {
    notes <- c(notes, if (fixed) {
      sprintf("The model has %d parameters and no residual degrees of freedom: it is saturated. Its likelihood-ratio and Wald tests are still valid, but goodness of fit and overdispersion cannot be checked.",
              np)
    } else {
      sprintf("The model has %d parameters and no residual degrees of freedom: it is saturated, so there is nothing left to test against. Every p-value, interval and effect size below is undefined or degenerate%s. Fewer groups, a simpler interaction structure, or more data are needed.",
              np, if (gaussian_lm) " (partial eta squared is exactly 1 by construction)" else "")
    })
  } else if (!is.null(dfr) && is.finite(dfr) && dfr < np && !fixed) {
    notes <- c(notes, sprintf(
      "The model has %d parameters and only %d residual degrees of freedom, so it is close to saturated. Consider fewer groups, a simpler interaction structure, or more data.",
      np, dfr))
  }
  if (any(is.na(stats::coef(fit)))) {
    bad <- gsub("`", "", names(stats::coef(fit))[is.na(stats::coef(fit))], fixed = TRUE)
    notes <- c(notes, sprintf(
      "%d coefficient(s) are aliased and could not be estimated: %s. This happens when level combinations are empty or when predictors are collinear; the affected terms are tested with fewer degrees of freedom than their levels suggest.",
      length(bad), paste(utils::head(bad, 5L), collapse = ", ")))
  }
  notes
}

#' Drop rows whose prior weight is zero
#'
#' A zero-weight row contributes nothing to the likelihood, but it still
#' counts as an observation for anything that counts rows: the sandwich
#' covariance divides by it (so every robust standard error shrinks), the
#' information criteria and the omega squared sample size include it, and a
#' group whose rows all have weight zero passes the check that every group has
#' data. Removing those rows up front makes every downstream count agree with
#' the fit. They are counted in \code{$n_removed}, like any other row that is
#' not in \code{$data_used}.
#'
#' @param data the prepared (complete-case) frame.
#' @param weights a single column name, or \code{NULL}.
#' @return list with \code{data}, \code{n_removed} and \code{notes}.
#' @noRd
.drop_zero_weights <- function(data, weights) {
  if (is.null(weights)) return(list(data = data, n_removed = 0L, notes = character(0)))
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
  zero <- w == 0
  if (!any(zero)) return(list(data = data, n_removed = 0L, notes = character(0)))
  out <- data[!zero, , drop = FALSE]
  row.names(out) <- NULL
  list(data = out, n_removed = sum(zero), notes = sprintf(
    "Dropped %d row(s) whose weight in `%s` is zero: they contribute nothing to the fit, and keeping them would make row counts, robust standard errors and information criteria disagree with it.",
    sum(zero), weights))
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
  dfr <- stats::df.residual(model)
  p <- if (use_t) {
    # A saturated model has no residual degrees of freedom; pt(df = 0) warns
    # "NaNs produced" to the console. There is no inference to report.
    if (is.null(dfr) || !is.finite(dfr) || dfr < 1) rep(NA_real_, length(stat))
    else 2 * stats::pt(-abs(stat), df = dfr)
  } else {
    2 * stats::pnorm(-abs(stat))
  }
  data.frame(term = names(est), estimate = unname(est), se = unname(se),
             statistic = unname(stat), p_value = unname(p),
             distribution = if (use_t) "t" else "z",
             stringsAsFactors = FALSE, row.names = NULL)
}

#' An environment that holds only the analysed frame
#'
#' A formula remembers the environment it was created in. Created inside an
#' analysis function, that is the function's frame, which holds the caller's
#' raw \code{data}: the fitted model then drags the whole input along when it
#' is saved, and anything that re-evaluates the model's call there --
#' \pkg{emmeans} rebuilding the data for an offset model, for one -- finds the
#' raw input instead of the rows that were analysed. Fitting in a small
#' environment that binds \code{data_used} to the analysed frame avoids both.
#' Its parent is the stats namespace, so \code{offset()} and \code{log()} in a
#' formula still resolve.
#' @noRd
.model_env <- function(data) {
  env <- new.env(parent = asNamespace("stats"))
  env$data_used <- data
  env
}

#' A reference to an object of a model environment that evaluates anywhere
#'
#' \code{update()}, \code{lmtest::lrtest()}, \code{step()} and
#' \code{MASS::stepAIC()} re-evaluate a model's call in the caller's frame,
#' where a bare \code{data_used} does not exist. A call to \code{get()} that
#' carries the environment itself finds it from any frame, so those functions
#' work on \code{fit$model} as they would on a model the user fitted.
#' @noRd
.env_ref <- function(env, name = "data_used") {
  call("get", name, envir = env)
}

#' A family as a call that reads well and re-evaluates anywhere
#' @noRd
.family_call <- function(family) {
  std <- c("binomial", "quasibinomial", "poisson", "quasipoisson", "gaussian",
           "Gamma", "inverse.gaussian")
  if (family$family %in% std) {
    return(call(family$family, link = family$link))
  }
  NULL
}

#' Fit a linear model whose call and formula refer to the analysed frame
#' @noRd
.fit_lm <- function(fml, data) {
  env <- .model_env(data)
  environment(fml) <- env
  eval(bquote(stats::lm(formula = .(fml), data = .(.env_ref(env)))), env)
}

#' Fit a GLM (or a negative binomial GLM) on the analysed frame
#'
#' The prior weights are named by their column, so they are subsetted with the
#' data and cannot be shadowed by a local variable. The stored call reaches
#' the frame the model was fitted on through \code{.env_ref()}, so
#' \code{update(fit$model, . ~ . - g)} refits it from anywhere.
#'
#' @param fml model formula.
#' @param data the prepared frame.
#' @param family a family object, or \code{NULL} for \code{MASS::glm.nb}.
#' @param weights name of the prior-weights column in \code{data}, or
#'   \code{NULL}.
#' @noRd
.fit_glm <- function(fml, data, family = NULL, weights = NULL) {
  env <- .model_env(data)
  environment(fml) <- env
  wsym <- if (is.null(weights)) NULL else as.name(weights)
  dref <- .env_ref(env)
  if (is.null(family)) {
    call <- if (is.null(wsym)) {
      bquote(MASS::glm.nb(formula = .(fml), data = .(dref)))
    } else {
      bquote(MASS::glm.nb(formula = .(fml), data = .(dref), weights = .(wsym)))
    }
  } else {
    fam <- .family_call(family)
    if (is.null(fam)) {
      env$family_used <- family
      fam <- .env_ref(env, "family_used")
    }
    call <- if (is.null(wsym)) {
      bquote(stats::glm(formula = .(fml), family = .(fam), data = .(dref)))
    } else {
      bquote(stats::glm(formula = .(fml), family = .(fam), data = .(dref),
                        weights = .(wsym)))
    }
  }
  eval(call, env)
}

#' Coefficient names that a model produces more than once
#' @noRd
.dup_coef_names <- function(fit) {
  co <- stats::coef(fit)
  nm <- if (is.matrix(co)) rownames(co) else names(co)
  unique(nm[duplicated(nm)])
}

#' Give every grouping factor contrasts whose column labels cannot collide
#'
#' A covariate named like a contrast coefficient -- \code{dose1} beside a
#' grouping factor \code{dose} -- makes the model produce two coefficients
#' of the same name, and emmeans then reads the wrong column. Used by
#' \code{\link{anova_ancova}} and \code{\link{anova_manova}} to refit when that
#' happens.
#'
#' The contrasts are the ones the fit would have used anyway -- sum-to-zero
#' under Type III, the session default otherwise, polynomial for an ordered
#' factor -- with their column labels wrapped in brackets, so a coefficient
#' reads \code{dose[1]} rather than \code{dose1}. No estimate changes.
#' @noRd
.bracket_group_contrasts <- function(d, groups, type) {
  opts <- if (identical(type, "III")) {
    c("contr.sum", "contr.poly")
  } else {
    getOption("contrasts", c("contr.treatment", "contr.poly"))
  }
  for (g in groups) {
    f <- d[[g]]
    fun <- if (is.ordered(f)) opts[[2L]] else opts[[1L]]
    M <- match.fun(fun)(levels(f))
    cn <- colnames(M)
    if (is.null(cn)) cn <- seq_len(ncol(M))
    colnames(M) <- paste0("[", cn, "]")
    stats::contrasts(f, ncol(M)) <- M
    d[[g]] <- f
  }
  d
}

#' Put an analysis-of-deviance table on the observations behind frequency weights
#'
#' car::Anova() refits reduced models, whose residual degrees of freedom glm()
#' counts in rows, so it is given the fit as glm() made it. Every statistic that
#' car scales by an estimated dispersion (all of them for a quasi family, the F
#' test for any family) is then too small by the factor rows / observations:
#' the dispersion is the Pearson statistic divided by the residual degrees of
#' freedom. The statistics are multiplied back, and F tests are referred to the
#' observation-based residual degrees of freedom. The result equals car's table
#' on the data expanded to one row per observation.
#'
#' @param av the result of \code{.car_anova()}.
#' @param df_row,df_obs residual degrees of freedom in rows and in observations.
#' @param dispersion_based does the table's statistic use an estimated
#'   dispersion?
#' @noRd
.freq_weight_anova <- function(av, df_row, df_obs, dispersion_based) {
  tab <- av$table
  if (is.null(tab) || !dispersion_based || !is.finite(df_row) || df_row < 1 ||
      !is.finite(df_obs) || df_obs < 1) {
    return(av)
  }
  k <- df_obs / df_row
  eff <- tab$term != "Residuals"
  tab$statistic[eff] <- tab$statistic[eff] * k
  if (identical(attr(tab, "statistic"), "F")) {
    tab$p_value[eff] <- stats::pf(tab$statistic[eff], tab$df[eff], df_obs,
                                  lower.tail = FALSE)
    tab$df[!eff] <- df_obs
  } else {
    tab$p_value[eff] <- stats::pchisq(tab$statistic[eff], tab$df[eff],
                                      lower.tail = FALSE)
  }
  av$table <- tab
  av
}
