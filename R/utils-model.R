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
  # A Wald chi-square test of a fixed-dispersion GLM uses only the
  # coefficients and their covariance, so it stays valid on a saturated fit
  # (see .check_model_size()); car's linearHypothesis() refuses residual df 0
  # all the same, and is given infinite residual df, which a chi-square test
  # does not read.
  if (identical(test_statistic, "Wald") && inherits(mod, "glm") &&
      !.estimates_dispersion(mod) && isTRUE(stats::df.residual(mod) < 1)) {
    mod$df.residual <- Inf
  }
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
    # car's error messages can span lines; $notes is one bullet per element.
    return(list(table = NULL, raw = NULL, note = sprintf(
      "The Type %s ANOVA table could not be computed: %s",
      type, trimws(gsub("[[:space:]]+", " ", conditionMessage(raw))))))
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

#' Fall back to Type II tests when Type III ones are undefined
#'
#' A model with aliased coefficients (an empty cell of the design, or
#' collinear predictors) has no Type III tests: car refuses them, and the
#' function would otherwise return no omnibus table at all -- by default, in
#' anova_ancova(). Type II tests are defined, so they are computed instead,
#' as anova_manova() does, and the method label and table heading say so.
#' @return list with \code{type} (the type to use) and \code{note}.
#' @noRd
.type_for_model <- function(model, type) {
  if (identical(type, "III") && anyNA(stats::coef(model))) {
    return(list(type = "II", note = "The model has aliased coefficients (an empty cell of the design, or collinear predictors), and Type III tests are not defined for it: car refuses them. Type II tests were computed instead."))
  }
  list(type = type, note = character(0))
}

#' Say whether the omnibus test uses a robust covariance
#'
#' Likelihood-ratio and F tests compare deviances or residual sums of squares
#' and cannot use a sandwich covariance; a Wald test can. With robust standard
#' errors requested the reader needs to know which kind the table holds. An F
#' test divides by a dispersion estimated from the Pearson residuals, so it
#' does not fix the dispersion at 1; what it does not use is the sandwich.
#' On a Poisson, binomial or negative binomial fit the note on the F test says
#' that instead, and this one is not added.
#' @param table what the function calls its table of coefficients: "the
#'   coefficient table" in \code{anova_glm()}, "$effect_sizes" in
#'   \code{anova_bin()} and \code{anova_count()}, which have no other.
#' @noRd
.omnibus_vcov_note <- function(test_statistic, vcov_type, robust_used,
                               table = "the coefficient table") {
  if (robust_used) {
    return(sprintf("The omnibus Wald test uses the robust (%s) covariance, as %s, the marginal means and the comparisons do.",
                   vcov_type, table))
  }
  if (identical(test_statistic, "F")) {
    return(sprintf("The omnibus F test compares deviances scaled by a dispersion estimated from the Pearson residuals, so it uses the model's variance function, not the robust (%s) covariance, which only %s, the marginal means and the comparisons use. For an omnibus test that uses it, set test_statistic = \"Wald\".",
                   vcov_type, table))
  }
  sprintf("The omnibus likelihood-ratio test compares deviances, so it rests on the model-based variance; the robust (%s) covariance is used only by %s, the marginal means and the comparisons. For an omnibus test that uses it, set test_statistic = \"Wald\".",
          vcov_type, table)
}

#' What overdispersion does to a fixed-dispersion fit with robust errors
#'
#' With a robust (sandwich) covariance the coefficient table, the marginal
#' means and the comparisons do not assume a dispersion of 1, so they are not
#' too narrow under overdispersion. A robust Wald omnibus test does not assume
#' it either, and an F test estimates the dispersion; only a likelihood-ratio
#' omnibus test still rests on it. Under an F test the note on that test says
#' which covariance each result uses, so this text does not repeat it.
#' @param size_text what the dispersion does to the size of a nominal 5\%
#'   test, as a clause.
#' @param remedy the argument that gives quasi-likelihood inference.
#' @param table as in \code{.omnibus_vcov_note()}.
#' @noRd
.robust_overdispersion_text <- function(test_statistic, vcov_type, size_text,
                                        remedy, table = "the coefficient table") {
  if (identical(test_statistic, "Wald")) {
    return(sprintf(
      "None of the results reported here rests on that assumption: %s, the marginal means and the comparisons use the robust (%s) covariance, and the omnibus Wald test uses it too.",
      table, vcov_type))
  }
  if (identical(test_statistic, "F")) {
    return(sprintf(
      "None of the results reported here rests on that assumption: the robust (%s) covariance does not, and the omnibus F test estimates the dispersion instead of fixing it at 1.",
      vcov_type))
  }
  sprintf(
    "%s, the marginal means and the comparisons use the robust (%s) covariance, which does not rest on that assumption, so their standard errors and intervals are not too small. Only the omnibus likelihood-ratio test does: its statistic is inflated by about that factor, and %s. test_statistic = \"Wald\" gives an omnibus test with the robust covariance, and %s gives quasi-likelihood inference throughout.",
    if (startsWith(table, "the ")) sub("^the", "The", table) else table,
    vcov_type, size_text, remedy)
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
  # car orders its columns by statistic (Df first for Wald, last for LR); one
  # order whatever the test keeps positional access and printing stable.
  lead <- intersect(c("term", "sum_sq", "mean_sq", "df", "statistic"), names(out))
  tail <- intersect("p_value", names(out))
  out <- out[c(lead, setdiff(names(out), c(lead, tail)), tail)]
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
#' Aliased coefficients are reported with their cause: the empty level
#' combinations of a crossed design, named by their levels, or collinearity,
#' naming the covariates involved (see \code{.aliased_notes()}).
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
  cf <- stats::coef(fit)
  na <- if (is.matrix(cf)) apply(is.na(cf), 1L, any) else is.na(cf)
  if (any(na)) notes <- c(notes, .aliased_notes(fit, na))
  notes
}

#' NA inference for a saturated GLM that estimates its dispersion
#'
#' With no residual degrees of freedom the dispersion is 0/0. car::Anova()
#' then reports an F statistic of 0 for some families and NaN for others,
#' emmeans gives NaN standard errors, and base R warns "NaNs produced" -- all
#' of it saying, less clearly, what the saturated-model note of
#' \code{.check_model_size()} says. The tests (\code{tests = TRUE}) become
#' NA, every NaN in the tables becomes NA, and the relayed "NaNs produced"
#' warnings are dropped.
#' @param out a named list of the result's tables and its notes.
#' @noRd
.saturated_na <- function(out, tests = TRUE) {
  nan_to_na <- function(tab) {
    if (!is.data.frame(tab)) return(tab)
    for (j in seq_along(tab)) {
      if (is.double(tab[[j]])) tab[[j]][is.nan(tab[[j]])] <- NA_real_
    }
    tab
  }
  for (nm in setdiff(names(out), "notes")) out[[nm]] <- nan_to_na(out[[nm]])
  tab <- out$anova
  if (tests && is.data.frame(tab) && all(c("statistic", "p_value") %in% names(tab))) {
    eff <- tab$term != "Residuals"
    tab$statistic[eff] <- NA_real_
    tab$p_value[eff] <- NA_real_
    out$anova <- tab
  }
  # A note from .said_note() is "<who> reported: <a> / <b>".
  notes <- out$notes
  said <- grepl(" reported: ", notes, fixed = TRUE)
  for (i in which(said)) {
    who <- sub(" reported: .*$", "", notes[i])
    parts <- strsplit(sub("^.*? reported: ", "", notes[i], perl = TRUE), " / ", fixed = TRUE)[[1L]]
    if (!"NaNs produced" %in% parts) next
    parts <- parts[parts != "NaNs produced"]
    notes[i] <- if (length(parts)) .said_note(parts, who) else NA_character_
  }
  out$notes <- notes[!is.na(notes)]
  out
}

#' Say why coefficients are aliased
#'
#' A coefficient of a term made only of grouping factors, some of whose level
#' combinations hold no data, is aliased because of those empty cells: the
#' note names the cells by their levels, not the coefficient, whose name under
#' sum-to-zero contrasts points at a different cell. The same holds for a
#' covariate-by-group term, whose slope in an empty combination of its grouping
#' factors cannot be estimated. Any other aliased
#' coefficient is collinear with the rest of the model, and the note names it
#' and any covariate in its term.
#' @param fit a fitted model.
#' @param na named logical, one per coefficient: is it aliased?
#' @noRd
.aliased_notes <- function(fit, na) {
  bad <- gsub("`", "", names(na)[na], fixed = TRUE)
  tail <- "the affected terms are tested with fewer degrees of freedom than their levels suggest."
  generic <- sprintf(
    "%d coefficient(s) are aliased and could not be estimated: %s. This happens when level combinations are empty or when predictors are collinear; %s",
    length(bad), .abbrev(bad, 5L), tail)
  info <- tryCatch(list(
    asg = attr(stats::model.matrix(fit), "assign"),
    fac = attr(stats::terms(fit), "factors"),
    labs = attr(stats::terms(fit), "term.labels"),
    mf = stats::model.frame(fit)), error = function(e) NULL)
  if (is.null(info) || length(info$asg) != length(na) || !is.matrix(info$fac)) {
    return(generic)
  }
  term_vars <- function(j) {
    if (j < 1L) return(character(0))
    gsub("`", "", rownames(info$fac)[info$fac[, j] > 0L], fixed = TRUE)
  }
  is_group <- function(v) {
    x <- info$mf[[v]]
    !is.null(x) && (is.factor(x) || is.character(x) || is.logical(x))
  }
  out <- character(0)
  other <- logical(length(na))
  for (j in unique(info$asg[na])) {
    idx <- na & info$asg == j
    vars <- term_vars(j)
    grp <- vars[vapply(vars, is_group, logical(1))]
    covs <- setdiff(vars, grp)
    if (length(grp) >= 2L) {
      ec <- .empty_cells(info$mf, grp)
      # A covariate-by-group term has one slope per combination of its
      # grouping factors, so an empty combination aliases at most one
      # coefficient per covariate product; more aliased than that is
      # collinearity as well, and is reported as such.
      if (ec$n_empty > 0 &&
          (length(covs) == 0L || sum(idx) <= ec$n_empty)) {
        term <- gsub("`", "", info$labs[j], fixed = TRUE)
        out <- c(out, if (length(covs) == 0L) {
          sprintf(
            "%d coefficient(s) of the %s term are aliased and could not be estimated because %s of its %s level combination(s) contain no data: %s; %s",
            sum(idx), term, .count_text(ec$n_empty), .count_text(ec$n_possible),
            .empty_cells_text(ec), tail)
        } else {
          sprintf(
            "%d coefficient(s) of the %s term are aliased and could not be estimated because %s of the %s level combination(s) of %s contain no data, so the slope of %s cannot be estimated in them: %s; %s",
            sum(idx), term, .count_text(ec$n_empty), .count_text(ec$n_possible),
            paste(grp, collapse = " and "),
            paste(sprintf("`%s`", covs), collapse = ", "),
            .empty_cells_text(ec), tail)
        })
        next
      }
    }
    other <- other | idx
  }
  if (any(other)) {
    vars <- unique(unlist(lapply(unique(info$asg[other]), term_vars)))
    covs <- vars[!vapply(vars, is_group, logical(1))]
    obad <- gsub("`", "", names(na)[other], fixed = TRUE)
    out <- c(out, sprintf(
      "%d coefficient(s) are aliased and could not be estimated because they are collinear with other terms of the model, not because of empty cells: %s.%s The affected terms are tested with fewer degrees of freedom than their levels suggest.",
      length(obad), .abbrev(obad, 5L),
      if (length(covs) > 0L) sprintf(
        " They involve the covariate(s) %s: a covariate that is a linear combination of other predictors, or constant within a group it interacts with, does this.",
        paste(sprintf("`%s`", covs), collapse = ", ")) else
        " Their grouping variables are confounded with other grouping variables."))
  }
  out
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
#' A group whose rows all have weight zero is lost with them; like a group
#' lost to missing values in \code{.prepare_frame()}, it is named in a note
#' and its level is dropped, so \code{$data_used} matches the fit.
#'
#' @param data the prepared (complete-case) frame.
#' @param weights a single column name, or \code{NULL}.
#' @param groups the grouping columns (factors) of \code{data}.
#' @return list with \code{data}, \code{n_removed} and \code{notes}.
#' @noRd
.drop_zero_weights <- function(data, weights, groups = character(0)) {
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
  for (g in groups) out[[g]] <- droplevels(out[[g]])
  list(data = out, n_removed = sum(zero), notes = c(sprintf(
    "Dropped %d row(s) whose weight in `%s` is zero: they contribute nothing to the fit, and keeping them would make row counts, robust standard errors and information criteria disagree with it.",
    sum(zero), weights),
    .lost_level_notes(data, !zero, groups,
                      what = sprintf("rows with a non-zero weight in `%s`", weights))))
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
  # A family whose initialisation was changed (see .null_start_family()) is
  # not the one the bare call would build.
  if (family$family %in% std && !isTRUE(attr(family, "null_start"))) {
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
#' @param control a \code{glm.control()} list, or \code{NULL} for the
#'   default. It is kept on the model, so \pkg{car}'s and \code{profile()}'s
#'   refits use it too.
#' @noRd
.fit_glm <- function(fml, data, family = NULL, weights = NULL, control = NULL) {
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
    if (!is.null(control)) {
      env$control_used <- control
      call$control <- .env_ref(env, "control_used")
    }
  }
  eval(call, env)
}

#' A binomial family that starts its iterations from the overall proportion
#'
#' \code{binomial()$initialize} starts each row at \code{(w y + 0.5) / (w + 1)}.
#' For a 0/1 response with frequency weights that puts every row next to 0 or
#' 1, and under a cloglog link the iterations then run off to coefficients of
#' 1e15 and stop at a deviance many times the true minimum, sometimes
#' reporting convergence. Starting every row at the weighted overall
#' proportion -- the fitted value of the intercept-only model -- avoids that.
#' The start is part of the family, so \pkg{car}'s and \code{drop1()}'s refits
#' of reduced models, which call \code{glm.fit()} with the model's family and
#' no starting values, start there too. Starting values given to the fit
#' (\code{mustart}, \code{etastart}, \code{start}) still take precedence.
#' @noRd
.null_start_family <- function(family) {
  init <- family$initialize
  family$initialize <- bquote({
    .(init)
    .w_pos <- weights > 0
    .m_all <- sum((weights * y)[.w_pos]) / sum(weights[.w_pos])
    if (is.finite(.m_all) && .m_all > 0 && .m_all < 1) {
      mustart <- rep(.m_all, length(y))
    }
    rm(.w_pos, .m_all)
  })
  attr(family, "null_start") <- TRUE
  family
}

#' Is a GLM's response constant, or fitted exactly?
#'
#' A family that estimates its dispersion divides every test statistic by it;
#' when the fit reproduces the response to rounding error so is the
#' dispersion, and the F statistics come out as 1e12, NaN or negative. A fit
#' is exact when the response takes one value in every row that carries
#' weight and the offset is constant too (a constant count over varying
#' exposures is a varying rate, which the model need not fit), or when every
#' fitted value equals its response to within 1e-13 of its size: the
#' rounding error of the fit, below real variation of even 1e-12 of the
#' response's level (counts of 1e12 differing by a few), which lm() and car
#' resolve.
#'
#' A response constant within every group but 0 (or 1, for a proportion) in
#' some group is fitted only in the limit: glm() stops with fitted values of
#' about 1e-11 there and F statistics of 1e12 to 1e17. When some response is
#' at such a boundary the fit is therefore also exact when every boundary
#' value is fitted to within 1e-6 of the response's range and every other
#' value to within 1e-11 of its own size (only for families bounded there:
#' Poisson, negative binomial and binomial, with their quasi forms). A
#' response constant within every design cell of a model saturated in its
#' cells (\code{.constant_within_cells()}) and fitted to within 1e-10 of its
#' size is exact too.
#'
#' The deviance is not used here: on a large level it can be mostly rounding
#' error (or negative) while the fit leaves real residual variation, and a
#' fixed fraction of the null deviance refuses real variation in a group of
#' small values beside a group of large ones. Whether the deviance resolves
#' the variation the fit leaves is judged test by test in
#' \code{.deviance_precision()}.
#' @return \code{NULL} when the fit leaves residual variation; otherwise the
#'   error message.
#' @noRd
.glm_exact_fit_message <- function(model, response) {
  y <- as.numeric(model$y)
  mu <- as.numeric(model$fitted.values)
  off <- model$offset
  keep <- rep(TRUE, length(y))
  w <- model$prior.weights
  if (!is.null(w) && length(w) == length(y)) keep <- w > 0
  y <- y[keep]
  mu <- if (length(mu) == length(keep)) mu[keep] else NULL
  off <- if (length(off) == length(keep)) as.numeric(off[keep]) else NULL
  if (length(y) == 0L || anyNA(y)) return(NULL)
  const_off <- is.null(off) || anyNA(off) ||
    all(abs(off - off[1L]) <= 1e-12 * max(1, abs(off)))
  # The same value to a few units in the last place: whole-number counts of
  # 1e12 that differ by 1 are different values.
  if (all(abs(y - y[1L]) <= 8 * .Machine$double.eps * max(abs(y))) &&
      (const_off || all(y == 0))) {
    return(sprintf(
      "Response `%s` takes the same value (%s) in every analysed row, so the model fits it exactly and leaves no residual variation: the dispersion cannot be estimated, and every test, interval and comparison would be a ratio of rounding errors. There is no difference between the groups to test.",
      response, format(signif(y[1L], 6))))
  }
  what <- NULL
  if (length(mu) == length(y) && all(is.finite(mu))) {
    # A response at a boundary of its family (0, or 1 for a proportion) is
    # fitted by glm() only in the limit, to about 1e-11: those rows are
    # measured against the response's range; every other row against its own
    # size, so that a constant group of large values cannot hide real
    # variation in a group of small ones.
    famobj <- tryCatch(stats::family(model), error = function(e) NULL)
    fam <- if (is.null(famobj)) "" else famobj$family
    binom <- fam %in% c("binomial", "quasibinomial") ||
      (identical(fam, "quasi") && identical(famobj$varfun, "mu(1-mu)"))
    # Only a family whose mean is bounded at 0 (and 1, for a proportion) is
    # fitted there in the limit; a Gaussian 0 is an ordinary value.
    bounded <- binom || fam %in% c("poisson", "quasipoisson") ||
      startsWith(fam, "Negative Binomial") ||
      (identical(fam, "quasi") && identical(famobj$varfun, "mu"))
    bnd <- bounded & (y == 0 | (binom & y == 1))
    spread <- diff(range(y))
    # Each residual is measured against its own row's size, so that a
    # constant group of large values cannot hide real variation in a group of
    # small ones; but a residual within a few dozen units in the last place
    # of the response's overall level is rounding error whatever its row
    # (a Gaussian group at 0 beside groups at 1e9 is fitted to about 1e-7).
    # (Only where the fit is linear in the response -- identity link,
    # constant variance -- is its rounding set by the overall level; a log
    # link fits each group on its own scale.)
    gauss_like <- identical(famobj$link, "identity") &&
      (identical(fam, "gaussian") ||
         (identical(fam, "quasi") && identical(famobj$varfun, "constant")))
    floor_ <- if (gauss_like) 32 * .Machine$double.eps * max(abs(y)) else 0
    rel_of <- function(mu) {
      res <- abs(y - mu)
      scl <- pmax(abs(y), abs(mu))
      r <- ifelse(res <= floor_, 0, ifelse(scl > 0, res / scl, Inf))
      if (any(!bnd)) max(r[!bnd]) else NA_real_
    }
    rel_i <- rel_of(mu)
    # An iteratively fitted family reaches an exact response only to its
    # convergence (about 3e-12 for a quasi-Poisson fit near 1), not to
    # rounding: when the residuals are that small, the fit is pushed to full
    # convergence before it is judged. Real variation stays where it is.
    # The trigger is wide (1e-4): glm() stops by a change in deviance
    # relative to the deviance plus 0.1, so on a small level an exact
    # response can be left at 1e-5. The refit's tolerance is made relative
    # to the deviance for the same reason.
    if (!gauss_like && isTRUE(rel_i > 1e-13 && rel_i <= 1e-4)) {
      X <- tryCatch(stats::model.matrix(model), error = function(e) NULL)
      cf <- stats::coef(model)
      if (!is.null(X) && nrow(X) == length(keep)) {
        X <- X[, !is.na(cf), drop = FALSE]
        o <- model$offset
        fr <- tryCatch(suppressWarnings(stats::glm.fit(
          X, as.numeric(model$y), weights = model$prior.weights,
          start = cf[!is.na(cf)],
          offset = if (length(o) == length(keep)) o else NULL,
          family = famobj,
          control = list(epsilon = 1e-15 * min(1, model$deviance /
                                                  (model$deviance + 0.1)),
                         maxit = 100, trace = FALSE))),
          error = function(e) NULL)
        if (!is.null(fr) && all(is.finite(fr$fitted.values))) {
          mu <- as.numeric(fr$fitted.values)[keep]
          rel_i <- rel_of(mu)
        }
      }
    }
    res <- abs(y - mu)
    if (any(bnd) && spread > 0) {
      rel_b <- max(res[bnd]) / spread
      if (rel_b <= 1e-6 && (!any(!bnd) || isTRUE(rel_i <= 1e-11))) {
        what <- if (any(!bnd)) {
          sprintf("every fitted value equals the response to within %s of the response's range where the response is %s, and to within %s of its own size elsewhere, which is rounding error",
                  format(signif(rel_b, 3)), if (binom) "0 or 1" else "0",
                  format(signif(rel_i, 3)))
        } else {
          sprintf("every response is %s, and every fitted value equals it to within %s of the response's range, which is rounding error",
                  if (binom) "0 or 1" else "0", format(signif(rel_b, 3)))
        }
      }
    } else if (!any(bnd) && isTRUE(rel_i <= 1e-10) &&
               .constant_within_cells(model, keep, off)) {
      what <- "the response is constant within every design cell and every fitted value equals it to rounding error"
    } else if (!any(bnd) && isTRUE(if (gauss_like) rel_i == 0 else rel_i <= 1e-13)) {
      # The fit reproduces every response to the rounding error of the fit
      # itself (an exactly additive two-factor response, say) -- for a fit
      # linear in the response, to the floor above; otherwise to 1e-13 of
      # each row's size: the Pearson
      # residuals, and so the estimated dispersion, are rounding error. Real
      # variation of even 1e-12 of the response's level (counts of 1e12
      # differing by a few) is above this; whether the deviance can resolve
      # it is a separate question, which .deviance_precision() answers for
      # each test.
      what <- "every fitted value equals the response to the rounding error of the fit"
    }
  }
  if (is.null(what)) return(NULL)
  sprintf(
    "The model fits response `%s` exactly: %s, as when the response is constant within every group. There is no residual variation to estimate the dispersion from, so the F tests, intervals and comparisons are undefined (each F statistic would be a ratio of rounding errors).",
    response, what)
}

#' Is a GLM's response constant within every design cell?
#'
#' A design cell is a set of analysed rows (\code{keep}) that share their row
#' of the model matrix and their offset, so the model cannot tell them apart.
#' The response is constant within a cell when its values agree to within a
#' few units in the last place, which a constant entered or computed in
#' double precision meets and real variation, even 1 on a level of 1e14,
#' does not. FALSE unless the model has a
#' coefficient for every cell and some cell has two rows or more.
#' @noRd
.constant_within_cells <- function(model, keep, off = NULL) {
  X <- tryCatch(stats::model.matrix(model), error = function(e) NULL)
  y <- as.numeric(model$y)
  if (is.null(X) || nrow(X) != length(keep) || length(y) != length(keep)) {
    return(FALSE)
  }
  X <- X[keep, , drop = FALSE]
  y <- y[keep]
  key <- do.call(paste, c(lapply(seq_len(ncol(X)), function(j)
    format(X[, j], digits = 15)), sep = "\r"))
  if (!is.null(off)) key <- paste(key, format(off, digits = 15), sep = "\r")
  cell <- match(key, unique(key))
  # Only a model with a coefficient for every cell reproduces a response
  # constant within its cells; an additive model of the same cells leaves
  # real lack of fit, and cells of one row are trivially constant.
  if (model$rank < max(cell) || !any(tabulate(cell) >= 2L)) return(FALSE)
  lo <- stats::ave(y, cell, FUN = min)
  hi <- stats::ave(y, cell, FUN = max)
  all(hi - lo <= 8 * .Machine$double.eps * pmax(abs(lo), abs(hi)))
}

#' How precisely a GLM's deviance is computed
#'
#' A unit deviance is evaluated from \code{y} and the fitted mean \code{mu}.
#' Where the family's deviance is written as a difference of nearly equal
#' terms -- \code{y * log(y / mu) - (y - mu)} for the Poisson and negative
#' binomial, \code{log(y / mu) - (y - mu) / mu} for the Gamma -- each term
#' carries a rounding error of about \code{eps * y} (Poisson) or \code{eps}
#' (Gamma), whatever the size of the true unit deviance, which is about
#' \code{(y - mu)^2 / V(mu)}. On a large level with little variation (counts
#' of 1e9 varying by a few units, epoch seconds varying by a minute) the
#' deviance is then mostly rounding error, and so is every deviance
#' difference that a likelihood-ratio or F test is built on: car returns F
#' statistics that are arbitrary or negative. The Gaussian and inverse
#' Gaussian deviances are computed from \code{y - mu} directly and do not
#' cancel in this way.
#'
#' The fitted means are themselves rounded, to about \code{16 * eps * |mu|}
#' (a few units in the last place, times the conditioning of the fit), which
#' moves each unit deviance by about \code{2 w |y - mu| / V(mu)} times that.
#' Both are summed in quadrature: the errors of different rows are
#' independent, and the deviance is stationary at the fit, so a shared error
#' in the coefficients moves it only to second order. Each is about twice the
#' error a careful analysis gives, so it is a bound to compare with, not an
#' estimate. The attribute \code{cancels} is FALSE for the Gaussian and
#' inverse Gaussian deviances, which are computed from \code{y - mu}
#' directly; their rounding is that of the residuals themselves, which the
#' exact-fit guard covers.
#' @return the bound on the absolute error of \code{deviance(model)}.
#' @noRd
.deviance_rounding <- function(model) {
  y <- as.numeric(model$y)
  mu <- as.numeric(model$fitted.values)
  w <- model$prior.weights
  if (is.null(w) || length(w) != length(y)) w <- rep(1, length(y))
  ok <- w > 0 & is.finite(y) & is.finite(mu)
  y <- y[ok]; mu <- mu[ok]; w <- w[ok]
  if (length(y) == 0L) return(0)
  eps <- .Machine$double.eps
  fam <- tryCatch(stats::family(model), error = function(e) NULL)
  fname <- if (is.null(fam)) "" else fam$family
  vtype <- if (identical(fname, "quasi")) fam$varfun else
    if (fname == "gaussian") "constant" else
    if (fname %in% c("poisson", "quasipoisson")) "mu" else
    if (fname == "Gamma") "mu^2" else
    if (fname == "inverse.gaussian") "mu^3" else
    if (fname %in% c("binomial", "quasibinomial")) "mu(1-mu)" else
    if (startsWith(fname, "Negative Binomial")) "nb" else "mu"
  big <- pmax(abs(y), abs(mu))
  theta <- if (vtype == "nb") {
    th <- model$theta
    if (is.null(th)) th <- suppressWarnings(as.numeric(
      sub("^Negative Binomial\\((.*)\\)$", "\\1", fname)))
    if (length(th) == 1L && is.finite(th)) th else 0
  } else 0
  cancel <- switch(vtype,
    "constant" = , "mu^3" = rep(0, length(y)),
    "mu^2" = , "mu(1-mu)" = 2 * eps * w,
    "nb" = 2 * eps * w * (2 * big + theta),
    2 * eps * w * big)
  V <- if (is.null(fam)) rep(NA_real_, length(mu)) else
    tryCatch(fam$variance(mu), error = function(e) rep(NA_real_, length(mu)))
  fit_err <- ifelse(is.finite(V) & V > 0,
                    2 * w * abs(y - mu) / V * 16 * eps * abs(mu), 0)
  out <- sqrt(sum(cancel^2) + sum(fit_err^2))
  attr(out, "cancels") <- !vtype %in% c("constant", "mu^3")
  out
}

#' Make NA the deviance tests that rounding error swamps
#'
#' A likelihood-ratio or F test of a GLM compares two deviances; their
#' difference is uncertain by about twice \code{.deviance_rounding()} (the
#' reduced model's deviance carries the same rounding as the full one's). In
#' the units of the reported statistic that is the difference itself for a
#' fixed-dispersion likelihood-ratio test, and the difference divided by the
#' Pearson dispersion (and, for F, the term's degrees of freedom) when car
#' scales by one. The bound is about three times the error seen in practice.
#' A test whose statistic is uncertain by more than 20\% of itself or of its
#' size under the null hypothesis, whichever is larger (its degrees of
#' freedom for a chi-square, 1 for an F), cannot be read, and is made NA
#' with a note; between 2\% and 20\% the test is kept and a note gives the
#' uncertainty. A negative statistic within its rounding error is shown as 0,
#' with p = 1. Families whose deviance does not cancel (Gaussian, inverse
#' Gaussian) are left alone. A Wald test does not use the
#' deviance and can still be computed. Applied to car's table
#' before any frequency-weight rescaling, which scales the statistic and its
#' error alike.
#' @param av the result of \code{.car_anova()}.
#' @param scaled does the table divide the deviance by an estimated
#'   dispersion?
#' @return \code{av}, with the affected rows NA and a note added to
#'   \code{av$note}.
#' @noRd
.deviance_precision <- function(av, model, scaled, k = 1) {
  tab <- av$table
  stat <- attr(tab, "statistic")
  if (is.null(tab) || nrow(tab) == 0L || is.null(stat) ||
      !(identical(stat, "F") || startsWith(stat, "likelihood-ratio") ||
        startsWith(stat, "scaled-deviance"))) {
    return(av)
  }
  rnd <- .deviance_rounding(model)
  if (!isTRUE(attr(rnd, "cancels"))) return(av)
  err_d <- 2 * as.numeric(rnd)
  phi <- if (scaled || identical(stat, "F")) .dispersion(model) else 1
  if (!is.finite(err_d) || err_d <= 0 || !isTRUE(is.finite(phi) && phi > 0)) {
    return(av)
  }
  eff <- tab$term != "Residuals" & is.finite(tab$statistic)
  dfs <- pmax(tab$df, 1)
  # On the scale of the statistic as reported: frequency weights multiply a
  # dispersion-based statistic, and so its error, by k (.freq_weight_anova()).
  err <- k * err_d / phi / if (identical(stat, "F")) dfs else 1
  orig <- k * tab$statistic
  # The statistic's size under the null hypothesis: its degrees of freedom
  # for a chi-square, 1 for an F (already divided by them).
  scale_ <- pmax(abs(orig), if (identical(stat, "F")) 1 else dfs)
  bad <- eff & err > 0.2 * scale_
  shaky <- eff & !bad & err > 0.02 * scale_
  # A negative statistic within its rounding error is rounding error: it is
  # shown as 0 (a deviance cannot rise when terms are added), with p = 1,
  # and noted only when the error itself matters on the statistic's scale.
  neg <- eff & !bad & tab$statistic < 0 & abs(orig) < err
  if (any(neg)) {
    tab$statistic[neg] <- 0
    if ("sum_sq" %in% names(tab)) tab$sum_sq[neg] <- 0
    tab$p_value[neg] <- 1
  }
  if (any(shaky)) {
    av$note <- c(av$note, sprintf(
      "The response varies little relative to its level, so this family's deviance is computed with limited precision: from rounding alone, the statistic(s) of %s are uncertain by about %s (on %s). The tests are still usable; test_statistic = \"Wald\", which uses the estimates and their covariance rather than the deviance, avoids this.",
      paste(tab$term[shaky], collapse = ", "),
      paste(format(signif(err[shaky], 2)), collapse = ", "),
      paste(format(signif(orig[shaky], 3)), collapse = ", ")))
  }
  if (!any(bad)) {
    av$table <- tab
    return(av)
  }
  tab$statistic[bad] <- NA_real_
  tab$p_value[bad] <- NA_real_
  if ("sum_sq" %in% names(tab)) tab$sum_sq[bad] <- NA_real_
  av$table <- tab
  av$note <- c(av$note, sprintf(
    "The test(s) of %s are NA: the response varies too little relative to its level for this family's deviance to resolve it. Each deviance is uncertain by about %s from rounding alone, so the change in deviance behind each test is uncertain by about %s on the scale of the statistic, which swamps it (the statistic came out as %s). This is a limit of computing the deviance in double precision, not a property of the data. test_statistic = \"Wald\", which uses the estimates and their covariance rather than the deviance, can still test them.",
    paste(tab$term[bad], collapse = ", "), format(signif(k * err_d / 2, 2)),
    paste(format(signif(err[bad], 2)), collapse = ", "),
    paste(format(signif(orig[bad], 3)), collapse = ", ")))
  av
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

#' Refit a GLM whose grouping factors produce one coefficient name twice
#'
#' Two grouping columns can collide on their own: a factor \code{a} with a
#' level \code{b1} and a factor \code{ab} with a level \code{1} both give a
#' coefficient \code{ab1}. The fit is still right, but every result looked up
#' by coefficient name -- the odds ratios or rate ratios, the marginal means --
#' then reads the wrong column. Used by \code{\link{anova_bin}} and
#' \code{\link{anova_count}}: when that happens the grouping factors get
#' bracketed contrast labels (see \code{.bracket_group_contrasts()}) and the
#' model is refitted, which changes no estimate.
#'
#' @param fitted the result of \code{.collect_conditions()} on the fit.
#' @param refit_fun function of the data that refits the model the same way,
#'   returning a \code{.collect_conditions()} result.
#' @param table,what where the caller reports the coefficients, and what they
#'   are, for the note.
#' @return list with \code{fitted} (the original or the refit), \code{data}
#'   (with bracketed contrasts when it was refitted), \code{bracketed} and
#'   \code{note}.
#' @noRd
.unclash_group_coefs <- function(fitted, refit_fun, d, groups, type,
                                 table = "$effect_sizes",
                                 what = "the effect sizes") {
  out <- list(fitted = fitted, data = d, bracketed = FALSE, note = character(0))
  if (inherits(fitted$value, "error")) return(out)
  dup <- .dup_coef_names(fitted$value)
  if (length(dup) == 0L) return(out)
  dup_lab <- gsub("`", "", dup, fixed = TRUE)
  involved <- groups[vapply(groups, function(g) any(startsWith(dup_lab, g)),
                            logical(1))]
  d2 <- .bracket_group_contrasts(d, groups, type)
  refit <- refit_fun(d2)
  if (inherits(refit$value, "error") ||
      length(.dup_coef_names(refit$value)) > 0L) {
    .stopf("Column names collide with the names of the model's coefficients (%s is produced twice), so %s and marginal means would be computed from the wrong columns. Rename the column(s) involved: %s.",
           paste(dup_lab, collapse = ", "), what,
           paste(sprintf("`%s`", involved), collapse = ", "))
  }
  relabelled <- grep("[", names(stats::coef(refit$value)), fixed = TRUE,
                     value = TRUE)
  out$note <- sprintf(
    "The names and levels of grouping columns %s produce the same coefficient name (%s would appear twice), which would make %s and the marginal means read the wrong model column. The grouping factors' contrast labels were put in brackets in $model and %s (%s, ...) to keep the names distinct; no estimate or test is affected.",
    paste(sprintf("`%s`", involved), collapse = ", "),
    paste(dup_lab, collapse = ", "), what, table,
    paste(gsub("`", "", utils::head(relabelled, 2L), fixed = TRUE),
          collapse = ", "))
  out$fitted <- refit
  out$data <- d2
  out$bracketed <- TRUE
  out
}

#' A note when weights that are not frequencies do not average 1
#'
#' glm() prior weights that are not whole numbers are not frequencies, but
#' the likelihood still multiplies by them, so every model-based test and
#' interval scales with their arbitrary scale: weights normalised to sum to 1
#' make any data look like a handful of observations. Only the robust
#' covariance treats each row as one unit. Used by \code{\link{anova_bin}}
#' and \code{\link{anova_count}}.
#' @noRd
.weight_scale_note <- function(d, weights, tol = 0.1) {
  if (is.null(weights)) return(character(0))
  w <- as.numeric(d[[weights]])
  if (all(abs(w - round(w)) < sqrt(.Machine$double.eps))) return(character(0))
  m <- mean(w)
  if (!is.finite(m) || m <= 0 || abs(log(m)) <= log(1 + tol)) {
    return(character(0))
  }
  sprintf(
    "The weights in `%s` are not all whole numbers, so they are not read as frequencies, and they average %s rather than 1. The fit, the tests and the intervals use them as given, as glm() uses prior weights, so their scale matters: multiplying every weight by a constant multiplies each likelihood-ratio and Wald statistic by about that constant. For sampling weights, rescale them to average 1, or use a design-based method such as survey::svyglm().",
    weights, format(signif(m, 3)))
}

#' The note for effect sizes whose treatment-coded labels were bracketed
#'
#' A main fit under sum-to-zero contrasts can have distinct coefficient names
#' while the treatment coding of \code{$effect_sizes} does not; only the
#' latter is then relabelled.
#' @noRd
.bracketed_terms_note <- function(terms) {
  sprintf(
    "The names and levels of the grouping columns paste to the same coefficient name under treatment coding, so the contrast labels in $effect_sizes were put in brackets (%s, ...) to keep the names distinct; no estimate is affected.",
    paste(gsub("`", "", utils::head(grep("[", terms, fixed = TRUE, value = TRUE), 2L),
               fixed = TRUE), collapse = ", "))
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

#' A call fit to keep and print
#'
#' \code{match.call()} records the arguments as they were written. Through
#' \code{do.call(anova_glm, list(data = d, ...))} they were written as values:
#' the function itself instead of its name, and the whole data frame instead
#' of the name of one. Kept as they are, the fit would carry the raw input
#' (every column, analysed or not) and \code{print()} would deparse all of it
#' into its \code{Call:} line. The function is put back as its name, and a
#' data frame, a family object or a function given as a value is replaced by
#' a placeholder naming what it was.
#' @param cl the result of \code{match.call()}.
#' @param name the name of the function called.
#' @noRd
.compact_call <- function(cl, name) {
  if (!is.call(cl)) return(cl)
  if (!is.name(cl[[1L]]) && !is.call(cl[[1L]])) cl[[1L]] <- as.name(name)
  for (i in seq_along(cl)[-1L]) {
    a <- cl[[i]]
    if (is.data.frame(a)) {
      cl[[i]] <- as.name(sprintf("<data frame: %d rows>", nrow(a)))
    } else if (inherits(a, "family")) {
      cl[[i]] <- as.name(sprintf("<family: %s(link = %s)>", a$family, a$link))
    } else if (is.function(a)) {
      cl[[i]] <- as.name("<function>")
    }
  }
  cl
}

#' How far weights used as given inflate a GLM's model-based tests
#'
#' Prior weights that are not frequencies (sampling or importance weights,
#' say) enter a fixed-dispersion GLM's likelihood as given, and its tests
#' assume Var(y) = V(mu) / w. When the response's variance does not depend on
#' the weights -- Var(y) = phi V(mu) -- the estimates have the sandwich
#' covariance A^-1 B A^-1, with A = X' diag(W) X the model's information
#' (W the working weights) and B = phi X' diag(w W) X, against the A^-1 the
#' tests use. For each term the Wald statistic is then a sum of
#' lambda_j chi-square_1, where the lambda_j are the eigenvalues of the
#' term's block of A^-1 B A^-1 relative to that of A^-1; Satterthwaite's
#' approximation, c chi-square_nu with c = sum(lambda^2) / sum(lambda) and
#' nu = sum(lambda)^2 / sum(lambda^2), gives the size of a nominal 5\% test.
#' With equal weights w every lambda is phi w; with unequal ones the factor
#' depends on where the heavy rows fall, which no single summary of the
#' weights captures (weights of 9.1 on a tenth of the rows inflate a test by
#' 8 if they fall in every cell, and by less if they fall in a few).
#' Likelihood-ratio statistics behave the same way to first order. Weights
#' that are inverse variances would not inflate the tests, so callers say
#' which reading they assume.
#' @param phi the dispersion of the response's variance relative to the
#'   family's (estimated with the weights rescaled to average 1); divided by
#'   the dispersion an F test estimates, for the size of that test.
#' @param type the type of the tests: under Type II a term is tested in the
#'   model without the terms that contain it.
#' @return a data frame with one row per term: \code{term}, \code{df},
#'   \code{factor} (mean lambda), \code{size} (of a nominal 5\% test), or
#'   \code{NULL} when it cannot be computed.
#' @noRd
.weight_inflation <- function(model, phi, type = "II") {
  X <- tryCatch(stats::model.matrix(model), error = function(e) NULL)
  W <- model$weights
  w <- model$prior.weights
  if (is.null(X) || is.null(W) || is.null(w) || !is.finite(phi) || phi <= 0 ||
      nrow(X) != length(W) || length(w) != length(W)) return(NULL)
  keep <- !is.na(stats::coef(model))
  asg <- attr(X, "assign")[keep]
  X <- X[, keep, drop = FALSE]
  ok <- is.finite(W) & W > 0 & is.finite(w) & w > 0
  y <- as.numeric(model$y)
  off <- model$offset
  # lambda for the columns J of a design Xm with working weights Wm
  lam_for <- function(Xm, Wm, J) {
    A <- crossprod(Xm * sqrt(Wm))
    B <- phi * crossprod(Xm * sqrt(w[ok] * Wm))
    Ai <- tryCatch(solve(A), error = function(e) NULL)
    if (is.null(Ai)) return(NA_real_)
    Vt <- Ai %*% B %*% Ai
    tryCatch(Re(eigen(solve(Ai[J, J, drop = FALSE], Vt[J, J, drop = FALSE]),
                      only.values = TRUE)$values), error = function(e) NA_real_)
  }
  tt <- stats::terms(model)
  labs <- attr(tt, "term.labels")
  fac <- attr(tt, "factors")
  # Under Type II a term is tested in the model without the terms that
  # contain it (car's marginality), so its lambda comes from that fit; under
  # Type III (sum-to-zero contrasts) from the full model's coefficients.
  containing <- function(k) {
    if (is.null(fac) || ncol(fac) < length(labs)) return(integer(0))
    vars <- fac[, k] > 0
    which(vapply(seq_along(labs), function(j)
      j != k && all(fac[vars, j] > 0), logical(1)))
  }
  out <- lapply(sort(unique(asg[asg > 0])), function(k) {
    drop <- if (identical(type, "II")) containing(k) else integer(0)
    if (length(drop)) {
      cols <- !asg %in% drop
      Xr <- X[ok, cols, drop = FALSE]
      fr <- tryCatch(suppressWarnings(stats::glm.fit(
        Xr, y[ok], weights = w[ok],
        offset = if (length(off) == length(y)) off[ok] else NULL,
        family = stats::family(model))), error = function(e) NULL)
      if (is.null(fr) || any(!is.finite(fr$weights))) return(NULL)
      lam <- lam_for(Xr, fr$weights, which(asg[cols] == k))
    } else {
      lam <- lam_for(X[ok, , drop = FALSE], W[ok], which(asg == k))
    }
    if (anyNA(lam) || any(lam <= 0)) return(NULL)
    cc <- sum(lam^2) / sum(lam)
    nu <- sum(lam)^2 / sum(lam^2)
    df <- length(lam)
    data.frame(term = gsub("`", "", labs[k], fixed = TRUE), df = df,
               factor = mean(lam),
               size = stats::pchisq(stats::qchisq(0.95, df) / cc, nu,
                                    lower.tail = FALSE),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  if (is.null(out) || nrow(out) == 0L) NULL else out
}

#' Clauses describing what weights used as given do to the tests
#' @param wi the result of \code{.weight_inflation()}.
#' @return list with \code{by} (\"by about ... (why)\", to follow
#'   \"inflated\") and \code{size} (the size of the worst-affected test).
#' @noRd
.weight_inflation_text <- function(wi) {
  rng <- range(wi$factor)
  fac <- if (diff(rng) < 0.05 * rng[2L]) sprintf("about %.2f", mean(wi$factor)) else
    sprintf("between about %.2f and %.2f, depending on the term", rng[1L], rng[2L])
  worst <- wi[which.max(wi$size), ]
  list(
    by = sprintf("by %s (the weights are used as given and are not frequencies; computed from the sandwich covariance, assuming the weights are not inverse variances, which would not inflate them)", fac),
    size = sprintf("a nominal 5%% test of %s rejects roughly %.0f%% of true null hypotheses",
                   worst$term, 100 * worst$size))
}

#' One note on what unequal weights that are not frequencies do to the tests
#'
#' Shared by \code{\link{anova_glm}} and \code{\link{anova_bin}}, so that the
#' note always describes the test that \code{$anova} reports: none when that
#' is a robust Wald test (it does not rest on the weights' reading); the F
#' test's own inflation when an estimated dispersion divides it (that
#' removes the average inflation but not the unevenness); otherwise the
#' likelihood-ratio and Wald tests'. The response's own dispersion phi0 is 1
#' for a 0/1 response and a negative binomial, and the Pearson dispersion of
#' the weighted fit over the mean weight otherwise. No note is given when no
#' test is inflated to a size of 7.5\% or more.
#' @param wnum the weights of the analysed rows.
#' @param binary is the response 0/1?
#' @param robust the HC type in use, or \code{NULL}.
#' @return list with \code{note} and \code{applies} (the weights are unequal
#'   sampling weights, so the caller's overdispersion note, which assumes
#'   otherwise, should give way).
#' @noRd
.weight_note <- function(model, wnum, weights, binary, test_statistic, robust,
                         type) {
  out <- list(note = character(0), applies = FALSE)
  wp <- wnum[wnum > 0]
  if (length(wp) < 2L || stats::var(wp) <= 1e-4 * mean(wp)^2) return(out)
  out$applies <- TRUE
  if (!is.null(robust) && identical(test_statistic, "Wald")) return(out)
  fam <- tryCatch(stats::family(model)$family, error = function(e) "")
  disp_w <- .dispersion(model)
  phi0 <- if (binary || startsWith(fam, "Negative Binomial")) 1 else
    if (isTRUE(is.finite(disp_w) && disp_w > 0)) disp_w / mean(wp) else NA_real_
  f_kind <- identical(test_statistic, "F") || .estimates_dispersion(model)
  phi_t <- if (f_kind) {
    if (isTRUE(is.finite(disp_w) && disp_w > 0)) phi0 / disp_w else NA_real_
  } else phi0
  if (!is.finite(phi_t)) return(out)
  wi <- .weight_inflation(model, phi_t, type)
  if (is.null(wi) || !isTRUE(max(wi$size) >= 0.075)) return(out)
  wt <- .weight_inflation_text(wi)
  out$note <- sprintf(
    "The weights in `%s` are unequal and not frequencies, so the %s in $anova %s inflated %s: %s.%s",
    weights,
    if (f_kind) "F tests" else "likelihood-ratio and Wald tests",
    "are", wt$by, wt$size,
    if (is.null(robust)) {
      if (f_kind && !identical(test_statistic, "F")) {
        " The estimated dispersion removes the average inflation, not that of unequal weights, and the standard errors elsewhere are affected alike; vcov_type = \"HC3\" gives ones that are not."
      } else if (f_kind) {
        " The F test's estimated dispersion removes the average inflation, not that of unequal weights."
      } else {
        " The standard errors in the other tables are too small for the same reason; vcov_type = \"HC3\" gives ones that are not."
      }
    } else {
      sprintf(" The robust (%s) standard errors in the other tables do not depend on this, and test_statistic = \"Wald\" gives an omnibus test that does not either.", robust)
    })
  out
}
