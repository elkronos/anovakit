#' Analysis of Deviance for a Binary Response
#'
#' Fits a logistic regression for a binary response across one or more grouping
#' variables and reports an analysis of deviance table, odds ratios with
#' confidence intervals, estimated marginal probabilities, and pairwise
#' comparisons on the odds-ratio scale.
#'
#' @details
#' \strong{Interactions.} By default the model is additive. Set
#' \code{interaction = TRUE} for the full factorial, or to a whole number for
#' the highest order of interaction to include. An additive model cannot detect
#' an interaction between grouping variables, so if you have more than one
#' grouping variable this is a decision worth making deliberately.
#'
#' \strong{Type III sums of squares.} When \code{type = "III"} the model is
#' \emph{fitted} under sum-to-zero contrasts, which is what makes Type III
#' tests meaningful. Setting the global contrast option after fitting has no
#' effect on an existing model.
#'
#' \strong{Odds ratios.} Under sum-to-zero or polynomial contrasts an
#' exponentiated coefficient is not an odds ratio against anything a reader
#' would recognise. \code{$effect_sizes} is therefore computed from a
#' reference-coded (treatment-coded) fit of the same model, with every grouping
#' factor treated as unordered and \code{reference} honoured, whatever
#' \code{type}, the global \code{contrasts} option or the class of the factor.
#' \code{$anova} still comes from the Type II or Type III fit, which describes
#' the same fitted probabilities.
#'
#' \strong{Separation.} Complete and quasi-complete separation are detected and
#' reported in \code{$notes}. Under separation, Wald odds ratios and their
#' intervals are meaningless even though \code{glm()} reports convergence, so
#' \code{ci_method = "profile"} is the default. The profile-likelihood interval
#' of an odds ratio that separation drives to infinity (or to zero) is open on
#' that side, and is reported with \code{conf_high = Inf} (or
#' \code{conf_low = 0}); its finite end is found by profiling the likelihood
#' directly, since \code{stats::confint()} cannot step from a diverged estimate,
#' and is \code{NA} when the likelihood rules out no value on that side either.
#'
#' \strong{Sparse data.} When fewer than about 5 events or non-events are
#' expected in a cell under the null hypothesis, the likelihood-ratio omnibus
#' test is liberal (its false-positive rate can be well above the nominal
#' level) and \code{$notes} says so. \code{test_statistic = "Wald"} is not a
#' remedy; an exact test, or the score test, holds its level better.
#'
#' \strong{Weights.} The prior weights enter the fit, the model statistics and
#' the observed proportions alike. A row with weight \code{w} counts as
#' \code{w} identical observations, so data aggregated to one row per cell and
#' outcome, with the count as the weight, give the same results as the
#' individual rows (apart from BIC, which R bases on the number of rows).
#'
#' @param data A data frame, or anything inheriting from one, such as a
#'   \code{data.table} or a tibble.
#' @param response Character. Name of the binary response column. May be a
#'   factor with two observed levels, a logical, a numeric 0/1 vector, or a
#'   character vector with two distinct values (whose levels are sorted in the
#'   C locale, so the default success level does not depend on the session).
#' @param groups Character vector. One or more grouping columns.
#' @param success Optional. A single value: the response level to model as the
#'   "success". By default the second level of the factor, which is what
#'   \code{glm()} uses. Works for ordered factors too.
#' @param reference Optional named list of reference levels for the grouping
#'   factors, one level per factor, for example \code{list(site = "north")}.
#'   An ordered factor has no reference level (its contrasts are polynomial
#'   trends), so naming one here makes the function treat that factor as
#'   unordered, with a note; its levels are never reordered as a scale.
#' @param interaction \code{FALSE} (additive, the default), \code{TRUE} (full
#'   factorial), or a whole number giving the highest interaction order.
#' @param type Character. \code{"II"} (default) or \code{"III"} sums of squares.
#' @param test_statistic Character. \code{"LR"} (default) or \code{"Wald"},
#'   passed to \code{\link[car]{Anova}}.
#' @param conf_level Numeric in (0, 1). Level for every interval returned.
#'   Default \code{0.95}.
#' @param ci_method Character. \code{"profile"} (default) for
#'   profile-likelihood intervals, or \code{"wald"}.
#' @param vcov_type Character. \code{"model"} (default) or an HC type
#'   (\code{"HC0"} to \code{"HC4"}) for robust standard errors, which requires
#'   the \pkg{sandwich} package. They reach the coefficient table, the
#'   marginal means and the comparisons; the omnibus table uses them only with
#'   \code{test_statistic = "Wald"}, since a likelihood-ratio test compares
#'   deviances. \code{$notes} says which applies.
#' @param adjust Character. Multiplicity adjustment for the pairwise
#'   comparisons, passed to \code{\link[emmeans]{contrast}}. Default
#'   \code{"tukey"}.
#' @param posthoc Logical. Compute pairwise comparisons. There are
#'   \code{choose(k, 2)} of them, so this is worth turning off when the number
#'   of cells is large; \code{$notes} records that they were skipped. Default
#'   \code{TRUE}.
#' @param plots Logical. Build \pkg{ggplot2} objects. They are returned in
#'   \code{$plots}, never drawn. Default \code{TRUE}.
#' @param verbose Logical. Emit progress through \code{\link[base]{message}}.
#'   Default \code{FALSE}.
#' @param weights Optional character. Name of a numeric column of prior
#'   weights. Given as a column name rather than a vector so that it is
#'   subsetted with the data: a vector supplied by the caller is evaluated
#'   against the original frame and silently misaligns as soon as one row is
#'   dropped for missing values. Whole-number weights (some above 1) are read
#'   as frequency weights, each row standing for that many identical
#'   observations, and the robust covariance (\code{vcov_type}) is then that
#'   of the data expanded to one row per observation; other weights are
#'   treated as sampling weights, with each row as one unit.
#'
#' @return An \code{\link{anovakit_fit}} object. Besides the standard
#'   components:
#'   \describe{
#'     \item{\code{$model_stats}}{A list: \code{AIC}, \code{BIC} and
#'       \code{logLik}; \code{mcfadden_r2}, McFadden's pseudo R squared,
#'       \code{1 - deviance / null deviance}; \code{success_level}, the level
#'       modelled as the success; and \code{lr_vs_null}, the model-versus-null
#'       likelihood-ratio test (\code{statistic} = null deviance minus
#'       deviance, \code{df}, \code{p_value}). All of them use the prior
#'       weights. When any weight is not a whole number, \code{AIC},
#'       \code{BIC} and \code{logLik} are \code{NA} with a note, because R's
#'       binomial log-likelihood rounds the weights; the deviance-based
#'       statistics are unaffected. \code{BIC} uses \code{nobs()}, the number
#'       of rows with a non-zero weight.}
#'     \item{\code{$effect_sizes}}{One row per coefficient of the
#'       reference-coded model (see Details): \code{term} (the coefficient),
#'       \code{factor} (the grouping variable, or \code{"g1:g2"} for an
#'       interaction), \code{comparison} (for example \code{"b vs a"}, a level
#'       against the reference level; in a model with interactions a lower-order
#'       odds ratio holds at the reference level of the factors it interacts
#'       with, \code{"b vs a at site = north"}, and an interaction row such as
#'       \code{"(b vs a) x (south vs north)"} is a ratio of odds ratios),
#'       \code{odds_ratio}, \code{conf_low}, \code{conf_high},
#'       \code{se_log_or}, \code{statistic}, \code{p_value} and
#'       \code{ci_method}.}
#'     \item{\code{$assumptions}}{\code{dispersion}, which is \code{NA}: the
#'       Pearson dispersion of a 0/1 response carries no information about
#'       overdispersion, so it is not reported (a note says why); and
#'       \code{proportions}, the observed proportion of each response level in
#'       each cell, with columns for the cell, the response level, \code{n}
#'       (the weighted count), \code{group_total} (the weighted cell total) and
#'       \code{proportion}; when \code{weights} is given, \code{rows} and
#'       \code{group_rows} give the number of data rows behind them. A column
#'       whose name the response already uses gets a suffix, e.g.
#'       \code{n.1}.}
#'   }
#'   \code{$emmeans} holds the estimated marginal probabilities and
#'   \code{$posthoc} the pairwise odds ratios.
#'
#' @seealso \code{\link{anova_count}} for counts, \code{\link{anova_glm}} for
#'   other families.
#'
#' @examples
#' set.seed(1)
#' n <- 300
#' d <- data.frame(g = factor(sample(c("a", "b", "c"), n, replace = TRUE)))
#' d$y <- rbinom(n, 1, c(a = 0.2, b = 0.5, c = 0.8)[as.character(d$g)])
#' fit <- anova_bin(d, "y", "g")
#' fit
#' fit$effect_sizes
#'
#' # Two factors, with their interaction
#' d$site <- factor(sample(c("north", "south"), n, replace = TRUE))
#' anova_bin(d, "y", c("g", "site"), interaction = TRUE, plots = FALSE)$anova
#'
#' # Aggregated data: one row per group and outcome, the count as the weight
#' agg <- data.frame(g = rep(c("a", "b", "c"), each = 2), y = rep(0:1, 3),
#'                   n = c(80, 20, 50, 50, 20, 80))
#' anova_bin(agg, "y", "g", weights = "n", plots = FALSE)$model_stats
#'
#' @export
anova_bin <- function(data, response, groups,
                      success = NULL,
                      reference = NULL,
                      interaction = FALSE,
                      type = c("II", "III"),
                      test_statistic = c("LR", "Wald"),
                      conf_level = 0.95,
                      ci_method = c("profile", "wald"),
                      vcov_type = "model",
                      adjust = "tukey",
                      posthoc = TRUE,
                      plots = TRUE,
                      weights = NULL,
                      verbose = FALSE) {
  cl <- match.call()
  type <- match.arg(type)
  test_statistic <- match.arg(test_statistic)
  ci_method <- match.arg(ci_method)
  .check_adjust(vcov_type, .vcov_choices(), "vcov_type")
  .check_adjust(adjust, .adjust_choices(), "adjust")
  data <- .as_df(data)
  .check_name(response, "response")
  .check_names(groups, "groups")
  .check_columns(data, response, "Response column")
  .check_columns(data, groups, "Grouping column(s)")
  .check_conf_level(conf_level)
  .check_interaction(interaction, length(groups))
  .check_roles(list(`the response` = response, `a grouping variable` = groups,
                    `the weights` = weights))
  .check_flag(posthoc, "posthoc")
  .check_flag(plots, "plots")
  .check_flag(verbose, "verbose")
  .check_success(success)
  .check_reference(reference, groups)

  .say(verbose, "Preparing data.")
  .check_weights_column(data, weights, response, groups)
  prep <- .prepare_frame(data, c(response, groups, weights), factors = groups,
                         keep_all = FALSE)
  d <- prep$data
  notes <- prep$notes
  zw <- .drop_zero_weights(d, weights)
  d <- zw$data
  n_removed <- prep$n_removed + zw$n_removed
  notes <- c(notes, zw$notes)
  .check_groups(d, groups, min_levels = 2L, min_n = 1L)

  resp <- .prep_binary(d[[response]], response, success)
  d[[response]] <- resp$values
  success <- resp$success
  notes <- c(notes, resp$notes, sprintf(
    "Modelling P(%s = %s); the other level is the baseline.", response, success))

  if (!is.null(reference)) {
    refd <- .set_references(d, groups, reference)
    d <- refd$data
    notes <- c(notes, refd$notes)
  }

  ## Fit ---------------------------------------------------------------------
  terms_rhs <- .group_terms(groups, interaction)
  fml <- stats::reformulate(terms_rhs, response = .bq(response))
  .say(verbose, "Fitting %s", paste(deparse(fml), collapse = " "))
  fitted <- .collect_conditions(.fit_with_contrasts(
    function() .fit_glm(fml, d, stats::binomial(), weights),
    type = type))
  model <- fitted$value
  if (inherits(model, "error")) {
    .stopf("The logistic regression could not be fitted: %s",
           conditionMessage(model))
  }
  # Non-integer prior weights make glm() warn about "non-integer #successes".
  # That is worth knowing and it is not an error, so it goes to $notes rather
  # than to the console.
  notes <- c(notes, .said_note(fitted$said, "stats::glm()"))
  if (!isTRUE(model$converged)) {
    notes <- c(notes, "The logistic regression did not converge; every result below is unreliable.")
  }
  notes <- c(notes, .check_model_size(model))

  # The odds ratios come from a reference-coded fit of the same model, so that
  # each one is a level against its reference level whatever `type`, the
  # contrasts option or the class of the factor. The separation check is
  # coding-invariant for the cells it flags; run on this fit, the coefficients
  # it names are the ones in $effect_sizes.
  tf <- .treatment_fit(model, fml, d, groups, weights)
  trt <- tf$model
  sep_note <- .check_separation(trt)
  notes <- c(notes, sep_note)

  ## Analysis of deviance ----------------------------------------------------
  # Whole-number weights on a 0/1 response count identical observations; the
  # robust covariance is then that of the expanded data.
  freq_rows <- .frequency_weights(d, weights)
  rv <- .robust_vcov(model, vcov_type, freq_weights = freq_rows)
  notes <- c(notes, rv$note)
  robust_wald <- !is.null(rv$matrix) && identical(test_statistic, "Wald")
  av <- .car_anova(model, type = type, test_statistic = test_statistic,
                   vcov_matrix = if (robust_wald) rv$matrix else NULL)
  notes <- c(notes, av$note)
  if (!is.null(rv$matrix) && !is.null(av$table)) {
    notes <- c(notes, .omnibus_vcov_note(test_statistic, vcov_type, robust_wald))
  }

  ## Model versus null, from the model's own (weighted) deviances -------------
  ms <- .bin_model_stats(model, success)
  notes <- c(notes, ms$note)

  ## Odds ratios -------------------------------------------------------------
  if (!is.null(rv$matrix) && ci_method == "profile") {
    ci_method <- "wald"
    notes <- c(notes, "Robust standard errors were requested, so the odds-ratio intervals are Wald intervals built from them; profile-likelihood intervals cannot use a sandwich covariance.")
  }
  rv_trt <- if (is.null(rv$matrix) || !tf$refit) rv$matrix else
    .robust_vcov(trt, vcov_type, freq_weights = freq_rows)$matrix
  eff <- .odds_ratios(trt, conf_level = conf_level, ci_method = ci_method,
                      vcov_matrix = rv_trt, separated = length(sep_note) > 0L)
  notes <- c(notes, eff$note)

  ## Marginal probabilities and pairwise odds ratios -------------------------
  emm <- .emm_block(model, groups, additive = .is_additive(interaction),
                    type = "response", vcov_matrix = rv$matrix,
                    conf_level = conf_level, adjust = adjust, posthoc = posthoc,
                    link = "logit", data = d)
  notes <- c(notes, emm$notes)
  if (!posthoc) {
    notes <- c(notes, "$effect_sizes still reports the odds ratios from the model.")
  }

  ## Proportions, sparse cells and plot ---------------------------------------
  added <- .add_cell(d, groups)
  prop <- .proportion_table(added$data, response, added$cell, weights)
  notes <- c(notes, .sparse_bin_note(d, response, groups, weights,
                                     additive = .is_additive(interaction)))
  # The Pearson X^2 / df of ungrouped 0/1 data is not a diagnostic: in a model
  # saturated in the cells it equals N / (N - k) whatever the data, and it
  # changes when the same trials are aggregated into weighted rows.
  notes <- c(notes, "$assumptions$dispersion is NA: the Pearson dispersion of a 0/1 response carries no information about overdispersion (in a model saturated in the cells it equals N / (N - k) whatever the data).")
  plot_list <- list()
  if (plots) {
    .say(verbose, "Building plots.")
    plot_list$proportions <- .bin_proportion_plot(
      prop, added$cell, response, success,
      xlab = paste(groups, collapse = " : "))
    if (!is.null(emm$table)) {
      plot_list$emmeans <- .plot_emmeans(
        emm$table, groups, estimate = "estimate",
        lower = "conf_low", upper = "conf_high", conf_level = conf_level,
        title = sprintf("Estimated probability of %s = %s", response, success),
        ylab = "Probability")
    }
    if (is.null(plot_list$emmeans)) notes <- c(notes, .no_emm_plot_note())
  }

  .new_fit(
    method       = "Analysis of deviance for a binary response (logistic regression)",
    call         = cl,
    model        = model,
    anova        = av$table,
    effect_sizes = eff$table,
    emmeans      = emm$table,
    emmeans_object = emm$grid,
    posthoc      = emm$posthoc,
    assumptions  = list(dispersion = NA_real_,
                        proportions = prop$table),
    plots        = plot_list,
    data_used    = d,
    n_removed    = n_removed,
    conf_level   = conf_level,
    notes        = notes,
    extra        = list(model_stats = ms$stats)
  )
}

#' Validate `success` before anything is fitted
#' @noRd
.check_success <- function(success) {
  if (is.null(success)) return(invisible(NULL))
  if (!is.atomic(success) || length(success) != 1L || is.na(success)) {
    .stopf("`success` must be a single, non-missing value: the one response level to model as the success.")
  }
  invisible(TRUE)
}

#' Validate `reference` before anything is fitted
#' @noRd
.check_reference <- function(reference, groups) {
  if (is.null(reference)) return(invisible(NULL))
  nms <- names(reference)
  if (!is.list(reference) || is.null(nms) || any(is.na(nms) | !nzchar(nms))) {
    .stopf("`reference` must be a named list, for example list(%s = \"...\").",
           groups[1L])
  }
  if (anyDuplicated(nms)) {
    .stopf("`reference` names each grouping variable at most once; %s %s repeated.",
           paste(sprintf("`%s`", unique(nms[duplicated(nms)])), collapse = ", "),
           if (sum(duplicated(nms)) == 1L) "is" else "are")
  }
  unknown <- setdiff(nms, groups)
  if (length(unknown) > 0L) {
    .stopf("`reference` names must be grouping variables; %s %s not.",
           paste(sprintf("`%s`", unknown), collapse = ", "),
           if (length(unknown) == 1L) "is" else "are")
  }
  for (g in nms) {
    r <- reference[[g]]
    if (!is.atomic(r) || length(r) != 1L || is.na(r)) {
      .stopf("`reference$%s` must be a single, non-missing level of `%s`.", g, g)
    }
  }
  invisible(TRUE)
}

#' Coerce and validate a binary response
#'
#' A numeric vector containing only zeros still produces a two-level factor when
#' the levels are supplied, which is how an all-one-outcome column can slip
#' through unnoticed. Levels are therefore counted after dropping empty ones.
#' A character response is sorted in the C locale, like the grouping columns,
#' so which level is the default success does not depend on the session.
#' @noRd
.prep_binary <- function(x, name, success = NULL) {
  notes <- character(0)
  if (is.logical(x)) {
    x <- factor(x, levels = c(FALSE, TRUE), labels = c("FALSE", "TRUE"))
  } else if (is.numeric(x)) {
    vals <- sort(unique(x[!is.na(x)]))
    if (!all(vals %in% c(0, 1))) {
      .stopf("Response `%s` is numeric but contains values other than 0 and 1: %s.",
             name, paste(utils::head(vals, 5L), collapse = ", "))
    }
    x <- factor(x, levels = c(0, 1))
  } else if (is.character(x)) {
    x <- .as_group_factor(x, name)
  } else if (!is.factor(x)) {
    x <- factor(x)
  }
  observed <- levels(droplevels(x))
  if (length(observed) != 2L) {
    .stopf("Response `%s` must take exactly two distinct values; %d observed (%s).",
           name, length(observed),
           if (length(observed) == 0L) "none" else paste(observed, collapse = ", "))
  }
  x <- droplevels(x)
  if (!is.null(success)) {
    success <- as.character(success)
    if (!success %in% observed) {
      .stopf("`success` = \"%s\" is not one of the observed levels of `%s`: %s.",
             success, name, paste(observed, collapse = ", "))
    }
    # Rebuilt rather than relevel()ed: relevel() refuses ordered factors, and
    # glm() only needs the level order (the second level is the success).
    x <- factor(as.character(x), levels = c(setdiff(observed, success), success))
  }
  list(values = x, success = levels(x)[2L], notes = notes)
}

#' Apply user-supplied reference levels to grouping factors
#'
#' An ordered factor has no reference level: it is coded by polynomial trends,
#' and moving a level to the front would only permute the scale those trends
#' run over. Such a factor is made unordered first, and a note says so.
#' @return list with \code{data} and \code{notes}.
#' @noRd
.set_references <- function(d, groups, reference) {
  .check_reference(reference, groups)
  notes <- character(0)
  for (g in names(reference)) {
    ref <- as.character(reference[[g]])
    if (!ref %in% levels(d[[g]])) {
      .stopf("Reference level \"%s\" not found in `%s`. Levels: %s.",
             ref, g, paste(levels(d[[g]]), collapse = ", "))
    }
    if (is.ordered(d[[g]])) {
      d[[g]] <- factor(d[[g]], levels = levels(d[[g]]), ordered = FALSE)
      notes <- c(notes, sprintf(
        "`%s` is an ordered factor, which has no reference level (its contrasts are polynomial trends). Because reference = list(%s = \"%s\") was given, it was treated as unordered, with \"%s\" as the reference level; its trend terms are not estimated.",
        g, g, ref, ref))
    }
    d[[g]] <- stats::relevel(d[[g]], ref = ref)
  }
  list(data = d, notes = notes)
}

#' The same model, fitted under reference (treatment) coding
#'
#' Every grouping factor is made unordered and fitted under
#' \code{contr.treatment}, so each coefficient is a level against the first
#' (reference) level. When the model already is coded that way it is returned
#' as it is.
#' @return list with \code{model} and \code{refit} (whether it was refitted).
#' @noRd
.treatment_fit <- function(model, fml, d, groups, weights) {
  ctr <- model$contrasts
  already <- !any(vapply(groups, function(g) is.ordered(d[[g]]), logical(1))) &&
    all(vapply(ctr, function(x) identical(x, "contr.treatment"), logical(1)))
  if (already) return(list(model = model, refit = FALSE))
  dt <- d
  for (g in groups) dt[[g]] <- factor(dt[[g]], levels = levels(dt[[g]]),
                                      ordered = FALSE)
  old <- options(contrasts = c("contr.treatment", "contr.poly"))
  on.exit(options(old), add = TRUE)
  # Any warning repeats one the main fit has already reported.
  refit <- .collect_conditions(.fit_glm(fml, dt, stats::binomial(), weights))$value
  if (inherits(refit, "error")) list(model = model, refit = FALSE) else
    list(model = refit, refit = TRUE)
}

#' Model-versus-null statistics from the fitted model's own deviances
#'
#' For a 0/1 response the saturated log-likelihood is zero, so the deviance is
#' -2 log L and McFadden's 1 - log L / log L0 equals 1 - D / D0. Using the
#' deviances keeps the prior weights -- aggregated counts, sampling weights --
#' exactly as the fit used them. R's binomial log-likelihood, and hence AIC and
#' BIC, rounds non-integer weights, so those are NA then.
#' @noRd
.bin_model_stats <- function(model, success) {
  w <- model$prior.weights
  whole <- is.null(w) || all(abs(w - round(w)) < 1e-8)
  dev <- model$deviance
  dev0 <- model$null.deviance
  lr <- max(dev0 - dev, 0)
  lr_df <- as.integer(model$df.null - model$df.residual)
  out <- list(
    AIC = if (whole) stats::AIC(model) else NA_real_,
    BIC = if (whole) stats::BIC(model) else NA_real_,
    logLik = if (whole) as.numeric(stats::logLik(model)) else NA_real_,
    mcfadden_r2 = if (is.finite(dev0) && dev0 > 0) 1 - dev / dev0 else NA_real_,
    success_level = success,
    lr_vs_null = data.frame(
      statistic = lr, df = lr_df,
      p_value = if (lr_df > 0L) stats::pchisq(lr, df = lr_df, lower.tail = FALSE)
                else NA_real_,
      stringsAsFactors = FALSE)
  )
  note <- if (whole) character(0) else
    "AIC, BIC and the log-likelihood in $model_stats are NA: some weights are not whole numbers, and R's binomial log-likelihood rounds them, so it would describe different data from the fit. McFadden's R squared and the model-versus-null test are computed from the deviances, which use the weights as given."
  list(stats = out, note = note)
}

#' Readable labels for the coefficients of a reference-coded model
#'
#' @return data.frame with \code{term}, \code{factor} and \code{comparison},
#'   one row per column of the model matrix.
#' @noRd
.or_labels <- function(model) {
  X <- stats::model.matrix(model)
  asg <- attr(X, "assign")
  tt <- stats::terms(model)
  labs <- attr(tt, "term.labels")
  fac <- attr(tt, "factors")
  xl <- model$xlevels
  unq <- function(x) gsub("`", "", x, fixed = TRUE)
  names(xl) <- unq(names(xl))
  out <- data.frame(term = colnames(X), factor = NA_character_,
                    comparison = NA_character_, stringsAsFactors = FALSE)
  term_vars <- lapply(labs, function(tl) rownames(fac)[fac[, tl] > 0])
  for (a in setdiff(unique(asg), 0L)) {
    vars <- term_vars[[a]]
    codes <- fac[vars, labs[a]]
    uv <- unq(vars)
    if (!all(uv %in% names(xl))) next
    lv <- lapply(seq_along(uv), function(i) {
      l <- xl[[uv[i]]]
      if (codes[i] == 1) l[-1L] else l
    })
    refs <- vapply(uv, function(v) xl[[v]][1L], character(1))
    combos <- expand.grid(lv, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    idx <- which(asg == a)
    if (nrow(combos) != length(idx)) next
    # model.matrix() names a column by pasting variable and level; check that
    # the combinations line up with the columns before labelling them.
    expect <- do.call(paste, c(lapply(seq_along(vars), function(i)
      paste0(unq(vars[i]), combos[[i]])), list(sep = ":")))
    if (!identical(expect, unq(colnames(X)[idx]))) next
    # Factors this term interacts with elsewhere in the model are held at
    # their reference level.
    higher <- unique(unlist(term_vars[vapply(term_vars, function(v)
      all(vars %in% v) && length(v) > length(vars), logical(1))]))
    others <- setdiff(higher, vars)
    at <- if (length(others) > 0L) {
      uo <- unq(others)
      sprintf(" at %s", paste(sprintf("%s = %s", uo,
                                      vapply(uo, function(v) xl[[v]][1L], character(1))),
                              collapse = ", "))
    } else ""
    comp <- vapply(seq_len(nrow(combos)), function(r) {
      lvl <- unlist(combos[r, ], use.names = FALSE)
      parts <- ifelse(codes == 1, sprintf("%s vs %s", lvl, refs), lvl)
      if (length(parts) == 1L) parts else
        paste(sprintf("(%s)", parts), collapse = " x ")
    }, character(1))
    out$factor[idx] <- unq(labs[a])
    out$comparison[idx] <- paste0(comp, at)
  }
  out
}

#' Odds ratios with profile-likelihood or Wald intervals
#'
#' @param model a reference-coded binomial GLM.
#' @param separated whether separation was detected. A coefficient it drives
#'   to infinity then gets its profile interval computed directly.
#' @noRd
.odds_ratios <- function(model, conf_level = 0.95, ci_method = "profile",
                         vcov_matrix = NULL, separated = FALSE) {
  ct <- .coef_table(model, vcov_matrix)
  notes <- character(0)
  if (ci_method == "profile") {
    ci <- tryCatch(
      suppressWarnings(suppressMessages(
        stats::confint(model, level = conf_level))),
      error = function(e) e)
    if (inherits(ci, "error")) {
      notes <- c(notes, sprintf(
        "Profile-likelihood intervals failed, Wald intervals were used instead: %s",
        conditionMessage(ci)))
      ci <- NULL
    } else {
      ci <- as.matrix(ci)
      if (ncol(ci) != 2L) ci <- matrix(ci, ncol = 2L,
                                       dimnames = list(ct$term, NULL))
    }
  } else {
    ci <- NULL
  }
  if (is.null(ci)) {
    z <- .coef_crit(model, conf_level)
    ci <- cbind(ct$estimate - z * ct$se, ct$estimate + z * ct$se)
    rownames(ci) <- ct$term
    method_used <- "wald"
  } else {
    method_used <- "profile"
  }
  idx <- match(ct$term, rownames(ci))
  lo <- ci[idx, 1L]
  hi <- ci[idx, 2L]

  # Under separation stats::confint() profiles from the diverged estimate with
  # steps scaled by its enormous standard error, and returns bounds that are
  # meaningless on both sides. A standard error above 10 on the log-odds scale
  # does not occur for a coefficient the data identify.
  if (method_used == "profile" && separated) {
    flag <- which(ct$term != "(Intercept)" & is.finite(ct$estimate) &
                    is.finite(ct$se) & ct$se > 10)
    for (k in flag) {
      b <- .separated_profile(model, ct$term[k], conf_level)
      lo[k] <- b[1L]
      hi[k] <- b[2L]
    }
    if (length(flag) > 0L) {
      notes <- c(notes, sprintf(
        "Separation drives the odds ratio(s) for %s to infinity or to zero, so the profile-likelihood interval is open on that side and is reported as conf_high = Inf (or conf_low = 0). Its finite end was found by profiling the likelihood directly, as stats::confint() cannot step from a diverged estimate; it is NA where the likelihood rules out no value on that side either.",
        paste(gsub("`", "", ct$term[flag], fixed = TRUE), collapse = ", ")))
    }
  }

  lab <- .or_labels(model)
  li <- match(ct$term, lab$term)
  out <- data.frame(
    term       = gsub("`", "", ct$term, fixed = TRUE),
    factor     = lab$factor[li],
    comparison = lab$comparison[li],
    odds_ratio = exp(ct$estimate),
    conf_low   = exp(lo),
    conf_high  = exp(hi),
    se_log_or  = ct$se,
    statistic  = ct$statistic,
    p_value    = ct$p_value,
    ci_method  = method_used,
    stringsAsFactors = FALSE, row.names = NULL
  )
  out <- out[out$term != "(Intercept)", , drop = FALSE]
  row.names(out) <- NULL
  list(table = out, note = notes)
}

#' Profile-likelihood interval of a coefficient that separation has diverged
#'
#' The profile deviance at a fixed value b of the coefficient is the deviance
#' of the model refitted with b times its column as an offset. The interval
#' holds every b whose profile deviance, divided by the dispersion, is within
#' \code{cutoff} of the minimum. In the direction the estimate diverged the
#' profile never rises that far, so that end is infinite; the other end is
#' found by stepping away from the estimate until the profile crosses the
#' cut-off and then solving for the crossing.
#'
#' A logistic model is refitted by a monotone Newton descent (see
#' \code{.logit_profile_deviance()}); any other GLM by \code{glm.fit()}, whose
#' IRLS starts from the data rather than from the diverged fit.
#'
#' @param cutoff the squared critical value: \code{qchisq(conf_level, 1)}
#'   for a fixed dispersion, the squared t quantile when it is estimated.
#' @param dispersion what the deviance is divided by: 1 for a fixed
#'   dispersion, the Pearson estimate otherwise, as \code{profile.glm()} does.
#' @return the two ends on the linear-predictor scale: +/-Inf for an open end,
#'   NA for a finite end that cannot be found.
#' @noRd
.separated_profile <- function(model, term, conf_level,
                               cutoff = stats::qchisq(conf_level, 1),
                               dispersion = 1) {
  X <- stats::model.matrix(model)
  co <- stats::coef(model)
  keep <- names(co)[!is.na(co)]
  X <- X[, keep, drop = FALSE]
  j <- match(term, colnames(X))
  est <- unname(co[term])
  y <- model$y
  w <- model$prior.weights
  off0 <- if (is.null(model$offset)) 0 else model$offset
  Xo <- X[, -j, drop = FALSE]
  start <- unname(co[keep][-j])
  fam <- stats::family(model)
  logit <- fam$family %in% c("binomial", "quasibinomial") && fam$link == "logit"
  if (logit) {
    dmin <- .logit_deviance(drop(X %*% co[keep]) + off0, y, w)
    dev_at <- function(off) .logit_profile_deviance(Xo, y, w, off, start)
  } else {
    dmin <- stats::deviance(model)
    dev_at <- function(off) {
      f <- tryCatch(suppressWarnings(stats::glm.fit(
        Xo, y, weights = w, offset = off, family = fam,
        control = stats::glm.control(maxit = 100L))), error = function(e) NULL)
      if (is.null(f) || !isTRUE(f$converged)) NA_real_ else f$deviance
    }
  }
  gap <- function(b) (dev_at(off0 + b * X[, j]) - dmin) / dispersion - cutoff
  find_end <- function(dir) {
    prev <- est
    step <- 0.5
    while (step <= 128) {
      b <- est + dir * step
      g <- gap(b)
      if (is.finite(g) && g > 0) {
        r <- tryCatch(stats::uniroot(gap, sort(c(prev, b)), tol = 1e-8)$root,
                      error = function(e) NA_real_)
        return(r)
      }
      prev <- b
      step <- step * 2
    }
    NA_real_
  }
  up <- if (est >= 0) 1 else -1           # the direction it diverged in
  open_end <- find_end(up)
  far_end <- find_end(-up)
  if (is.na(open_end)) open_end <- up * Inf
  if (up > 0) c(far_end, open_end) else c(open_end, far_end)
}

#' Deviance of a logistic model for a 0/1 response, from its linear predictor
#'
#' The saturated log-likelihood of a 0/1 response is zero, so the deviance is
#' -2 log L. Computed on the log scale so that a fitted probability at the
#' boundary does not underflow.
#' @noRd
.logit_deviance <- function(eta, y, w) {
  -2 * sum(w * (y * stats::plogis(eta, log.p = TRUE) +
                  (1 - y) * stats::plogis(-eta, log.p = TRUE)))
}

#' Minimum deviance of a logistic model with a fixed offset
#'
#' Newton-Raphson with step halving, so the deviance never increases. Plain
#' IRLS (glm.fit()) has no line search and, started far from the solution --
#' as it is when the offset moves a diverged coefficient a long way -- can
#' overshoot and settle on a deviance thousands of units too high, which would
#' put a spurious end on a profile interval. The log-likelihood is concave, so
#' a monotone descent from any start reaches the minimum (or, when some
#' coefficient diverges, its infimum).
#' @noRd
.logit_profile_deviance <- function(X, y, w, offset, start, maxit = 500L) {
  beta <- start
  eta <- offset + drop(X %*% beta)
  dev <- .logit_deviance(eta, y, w)
  p <- ncol(X)
  for (it in seq_len(maxit)) {
    mu <- stats::plogis(eta)
    score <- crossprod(X, w * (y - mu))
    info <- crossprod(X, X * (w * mu * (1 - mu)))
    ridge <- 1e-10 * max(1, diag(info))
    delta <- tryCatch(drop(solve(info + diag(ridge, p), score)),
                      error = function(e) NULL)
    if (is.null(delta) || !all(is.finite(delta))) break
    # Far from the solution the information is nearly zero and the Newton
    # step enormous; halving it as often as it takes is what makes the descent
    # work from any start.
    step <- 1
    repeat {
      nb <- beta + step * delta
      ne <- offset + drop(X %*% nb)
      nd <- .logit_deviance(ne, y, w)
      if (is.finite(nd) && nd <= dev) break
      step <- step / 2
      if (step < 2^-80) break
    }
    if (!is.finite(nd) || nd > dev) break
    done <- (dev - nd) < 1e-10 * (abs(nd) + 0.1)
    beta <- nb
    eta <- ne
    dev <- nd
    if (done) break
  }
  dev
}

#' Observed proportions of each response level within each cell
#'
#' Counts and totals are sums of the prior weights, so aggregated data (one
#' row per cell and outcome, weighted by the count) give the same table as the
#' individual rows, and the table agrees with the weighted fit beside it.
#' Cells are matched by position, never by name: an empty-string level is a
#' valid name for a group but not for indexing.
#'
#' @return list with \code{table} (the returned table) and \code{plot_data}
#'   (cell, level and proportion under the names the plot needs).
#' @noRd
.proportion_table <- function(d, response, cell, weights = NULL) {
  cf <- d[[cell]]
  rf <- d[[response]]
  w <- if (is.null(weights)) rep(1, nrow(d)) else as.numeric(d[[weights]])
  cl <- levels(cf)
  rl <- levels(rf)
  ci <- as.integer(cf)
  ri <- as.integer(rf)
  i <- rep(seq_along(cl), times = length(rl))   # the cell varies fastest
  k <- rep(seq_along(rl), each = length(cl))
  rows <- tabulate(ci + (ri - 1L) * length(cl), nbins = length(cl) * length(rl))
  wsum <- vapply(seq_along(i), function(m) sum(w[ci == i[m] & ri == k[m]]),
                 numeric(1))
  wtot <- vapply(seq_along(cl), function(m) sum(w[ci == m]), numeric(1))
  rtot <- tabulate(ci, nbins = length(cl))
  if (is.null(weights)) {
    n <- as.integer(rows)
    total <- as.integer(rtot)[i]
  } else {
    n <- wsum
    total <- wtot[i]
  }

  taken <- c(cell, response)
  pick <- function(candidate) {
    nm <- .safe_name(stats::setNames(vector("list", length(taken)), taken),
                     candidate)
    taken <<- c(taken, nm)
    nm
  }
  out <- data.frame(factor(cl[i], levels = cl), rl[k],
                    stringsAsFactors = FALSE)
  names(out) <- c(cell, response)
  nm_n <- pick("n")
  nm_total <- pick("group_total")
  nm_prop <- pick("proportion")
  out[[nm_n]] <- n
  out[[nm_total]] <- total
  prop <- ifelse(total > 0, n / total, NA_real_)
  out[[nm_prop]] <- prop
  if (!is.null(weights)) {
    out[[pick("rows")]] <- as.integer(rows)
    out[[pick("group_rows")]] <- as.integer(rtot)[i]
  }

  # The plot helper reads a column called "proportion" and one named after the
  # response; give it a frame in which those cannot be the same column.
  lvl <- if (response %in% c("proportion", cell)) ".level" else response
  plot_data <- data.frame(out[[cell]], out[[response]], prop,
                          stringsAsFactors = FALSE)
  names(plot_data) <- c(cell, lvl, "proportion")
  list(table = out, plot_data = plot_data, level = lvl)
}

#' The stacked proportion plot, labelled with the response name
#' @noRd
.bin_proportion_plot <- function(prop, cell, response, success, xlab) {
  p <- .plot_proportions(prop$plot_data, cell, prop$level, success, xlab = xlab)
  if (!identical(prop$level, response)) {
    p <- p + ggplot2::labs(
      title = sprintf("Proportion of %s by group", response),
      subtitle = sprintf("Modelled outcome: %s = %s", response, success),
      fill = response)
  }
  p
}

#' A note when too few events or non-events are expected in a cell
#'
#' The likelihood-ratio chi-square is liberal with sparse binary data: with 2
#' to 4 expected events per group its Type I error is 0.07 to 0.09 at a nominal
#' 0.05, mostly from outcomes in which a group has no events. Expected counts
#' are taken under the null hypothesis (the pooled proportion), as in Cochran's
#' rule, over the cells of the model: each factor's levels for an additive
#' model, the combined cells otherwise.
#' @noRd
.sparse_bin_note <- function(d, response, groups, weights, additive,
                             threshold = 5) {
  w <- if (is.null(weights)) rep(1, nrow(d)) else as.numeric(d[[weights]])
  y <- as.integer(d[[response]]) == 2L
  pbar <- sum(w * y) / sum(w)
  sets <- if (additive && length(groups) > 1L) as.list(groups) else list(groups)
  where <- character(0)
  smallest <- Inf
  for (gs in sets) {
    cell <- .cell_vector(d, gs)
    tot <- vapply(seq_len(nlevels(cell)),
                  function(m) sum(w[as.integer(cell) == m]), numeric(1))
    expct <- pmin(tot * pbar, tot * (1 - pbar))
    low <- expct < threshold
    if (any(low)) {
      lab <- if (length(gs) == 1L) sprintf("%s = %s", gs, levels(cell)[low]) else
        sprintf("%s = %s", paste(gs, collapse = " : "), levels(cell)[low])
      where <- c(where, lab)
      smallest <- min(smallest, expct[low])
    }
  }
  if (length(where) == 0L) return(character(0))
  sprintf(
    "Sparse data: fewer than %g events or non-events are expected under the null hypothesis in %s (smallest %s). With data this sparse the likelihood-ratio omnibus test is liberal, rejecting a true null more often than its nominal level. test_statistic = \"Wald\" is not a remedy: Wald tests are less reliable still here. An exact test (for one grouping variable, fisher.test() on the group-by-outcome table) or the score test, anova(fit$model, test = \"Rao\"), holds its level better.",
    threshold, .abbrev(unique(where), 5L), format(signif(smallest, 3)))
}
