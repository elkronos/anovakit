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
#' directly, since \code{stats::confint()} cannot step from a diverged estimate.
#' When the likelihood rules out no value on that side either (for example an
#' interaction whose cells are all separated) the interval is \code{(0, Inf)}:
#' the data do not identify that odds ratio at all, its point estimate is
#' arbitrary, and \code{$notes} says so. An end the search cannot find is
#' \code{NA}.
#'
#' \strong{Interval cost.} A profile-likelihood interval refits the model
#' several times for each coefficient. The model has only grouping factors,
#' so the refits run over the data collapsed to one weighted row per cell and
#' outcome, which has the same likelihood: the number of rows does not add to
#' the cost, and individual rows and the same data aggregated with counts as
#' weights get the same intervals. The number of coefficients does: a factor
#' with 200 levels can take minutes where Wald intervals take about a second.
#' Left at its default, \code{ci_method} therefore gives Wald intervals, with
#' a note, when the model has more than 50 coefficients, as
#' \code{\link{anova_glm}} does. An odds ratio that separation drives to
#' infinity or to zero still gets its profile-likelihood interval then, since
#' a Wald interval means nothing there. Wald and profile intervals are close
#' when every cell has many events and non-events, and can differ noticeably
#' when events are rare. Give \code{ci_method = "profile"} explicitly to
#' profile every coefficient.
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
#' individual rows, BIC included: with frequency weights its penalty uses the
#' number of observations they represent, not the number of rows.
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
#'   A model with aliased coefficients (an empty cell of the design, or
#'   collinear predictors) has no Type III tests; Type II tests are then
#'   computed instead, and \code{$notes}, the method and the table heading
#'   say so.
#' @param test_statistic Character. \code{"LR"} (default) or \code{"Wald"},
#'   passed to \code{\link[car]{Anova}}.
#' @param conf_level Numeric in (0, 1). Level for every interval returned.
#'   Default \code{0.95}.
#' @param ci_method Character. \code{"profile"} (default) for
#'   profile-likelihood intervals, or \code{"wald"}. Profile intervals are
#'   slow with many coefficients; left at the default, a model with more
#'   than 50 coefficients gets Wald intervals instead, with a note, except
#'   for odds ratios that separation drives to infinity or zero (see
#'   Details). The p-values in \code{$effect_sizes} are Wald p-values, so in
#'   small or sparse samples a profile interval can exclude an odds ratio of 1
#'   while the p-value beside it is above \code{1 - conf_level}, or the
#'   reverse; \code{$notes} names such rows, and \code{"wald"} intervals
#'   always agree with the p-values.
#' @param vcov_type Character. \code{"model"} (default) or an HC type
#'   (\code{"HC0"} to \code{"HC4"}) for robust standard errors, which requires
#'   the \pkg{sandwich} package. They reach \code{$effect_sizes}, the
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
#'   of the data expanded to one row per observation. Other weights are used
#'   as given, as \code{glm()} uses prior weights: the fit, the tests, the
#'   intervals and the marginal means treat them as multipliers of each row,
#'   like frequencies, so their scale matters (weights averaging 0.01 make
#'   300 rows count as 3 observations), and \code{$notes} says so when they
#'   do not average 1. Only the robust covariance, and the sparse-data check
#'   (which rescales them to average 1), treat each row as one unit. Sampling
#'   weights should therefore
#'   be rescaled to average 1, or analysed with a design-based method such as
#'   \code{survey::svyglm()}.
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
#'       statistics are unaffected. \code{BIC} uses \code{log(nobs())}, the
#'       number of rows with a non-zero weight, or, with frequency weights,
#'       the number of observations they represent (their sum), so that it
#'       equals the BIC of the expanded data.}
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
#'       \code{ci_method} (\code{"profile"} or \code{"wald"}, row by row:
#'       see Details).}
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
  cl <- .compact_call(match.call(), "anova_bin")
  type <- match.arg(type)
  test_statistic <- match.arg(test_statistic)
  ci_default <- missing(ci_method)
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
  data <- .as_plain_columns(data, list(Response = response,
                                       `Grouping column` = groups,
                                       `Weights column` = weights))
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
  zw <- .drop_zero_weights(d, weights, groups)
  d <- zw$data
  n_removed <- prep$n_removed + zw$n_removed
  notes <- c(notes, zw$notes, .weight_scale_note(d, weights))
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
  fit_on <- function(dd) .collect_conditions(.fit_with_contrasts(
    function() .fit_glm(fml, dd, stats::binomial(), weights),
    type = type))
  fitted <- fit_on(d)
  # Grouping columns whose names and levels paste to the same coefficient
  # name would make the odds ratios and marginal means read the wrong column.
  uc <- .unclash_group_coefs(fitted, fit_on, d, groups, type)
  fitted <- uc$fitted
  d <- uc$data
  notes <- c(notes, uc$note)
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
  tt <- .type_for_model(model, type)
  type <- tt$type
  notes <- c(notes, tt$note)

  # The odds ratios come from a reference-coded fit of the same model, so that
  # each one is a level against its reference level whatever `type`, the
  # contrasts option or the class of the factor. The separation check is
  # coding-invariant for the cells it flags; run on this fit, the coefficients
  # it names are the ones in $effect_sizes.
  tf <- .treatment_fit(model, fml, d, groups, weights)
  trt <- tf$model
  if (isTRUE(tf$bracketed) && !uc$bracketed) {
    notes <- c(notes, .bracketed_terms_note(names(stats::coef(trt))))
  }
  sep_note <- .check_separation(trt, test_statistic)
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
    notes <- c(notes, .omnibus_vcov_note(test_statistic, vcov_type, robust_wald,
                                         table = "$effect_sizes"))
  }
  # Unequal weights that are not whole numbers are sampling weights, used as
  # given: the model-based tests are inflated by an amount that depends on
  # where the heavy rows fall (.weight_note()). Rescaling them to average 1
  # removes only their scale, not this.
  if (!is.null(weights) && !freq_rows &&
      any(abs(d[[weights]] - round(d[[weights]])) > sqrt(.Machine$double.eps))) {
    notes <- c(notes, .weight_note(model, d[[weights]], weights, TRUE,
                                   test_statistic,
                                   if (!is.null(rv$matrix)) vcov_type, type)$note)
  }

  ## Model versus null, from the model's own (weighted) deviances -------------
  ms <- .bin_model_stats(model, success,
                         freq_n = if (freq_rows) sum(d[[weights]]) else NULL)
  notes <- c(notes, ms$note)

  ## Odds ratios -------------------------------------------------------------
  if (!is.null(rv$matrix) && ci_method == "profile") {
    ci_method <- "wald"
    notes <- c(notes, "Robust standard errors were requested, so the odds-ratio intervals are Wald intervals built from them; profile-likelihood intervals cannot use a sandwich covariance.")
  }
  # Profiling refits the model several times per coefficient. The refits run
  # over the data collapsed to one row per cell and outcome (see
  # .collapsed_bin_fit()), so the number of rows does not add to the cost,
  # but the number of coefficients does. Left at the default, a model with
  # many coefficients gets Wald intervals instead, as in anova_glm(), except
  # for the coefficients separation diverges, whose Wald intervals mean
  # nothing; asked for by name, profiling is honoured. The note comes before
  # the cost would have been incurred.
  n_or <- sum(!is.na(stats::coef(trt)))
  if (ci_default && ci_method == "profile" &&
      n_or > .max_profiled_coefficients) {
    ci_method <- if (length(sep_note) == 0L) "wald" else "separated"
    notes <- c(notes, sprintf(
      "The model has %d coefficients, and profiling the likelihood of each refits the model many times, so the odds-ratio intervals are Wald intervals%s; give ci_method = \"profile\" explicitly to profile them all.",
      n_or, if (ci_method == "separated") ", except for any odds ratio that separation drives to infinity or zero, whose profile-likelihood interval is computed directly (a Wald interval means nothing there)" else ""))
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
  if (posthoc && NROW(emm$table) == 2L && !isTRUE(emm$per_factor) &&
      !is.null(emm$posthoc)) {
    notes <- c(notes, .two_cell_note(model, av$table, test_statistic,
                                     robust = !is.null(rv$matrix),
                                     vcov_type = vcov_type))
  }

  ## Proportions, sparse cells and plot ---------------------------------------
  added <- .add_cell(d, groups)
  prop <- .proportion_table(added$data, response, added$cell, weights)
  notes <- c(notes, .sparse_bin_note(d, response, groups, weights,
                                     additive = .is_additive(interaction),
                                     frequencies = freq_rows))
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
        ylab = "Probability",
        ordered = is.ordered(data[[groups[1L]]]))
    }
    if (is.null(plot_list$emmeans)) notes <- c(notes, .no_emm_plot_note())
  }

  .new_fit(
    method       = sprintf(
      "Analysis of deviance for a binary response (logistic regression, Type %s)",
      type),
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
#' as it is. When grouping columns' coefficient names collide under treatment
#' coding (see \code{.unclash_group_coefs()}) the refit's contrasts get
#' bracketed labels, so its coefficient names stay distinct.
#' @return list with \code{model}, \code{refit} (whether it was refitted) and
#'   \code{bracketed} (whether its labels were bracketed).
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
  refit_on <- function(dd)
    .collect_conditions(.fit_glm(fml, dd, stats::binomial(), weights))$value
  refit <- refit_on(dt)
  # Grouping columns whose names and levels paste to the same coefficient
  # name under treatment coding (which a sum-to-zero main fit can avoid)
  # get bracketed labels, as in .unclash_group_coefs().
  bracketed <- !inherits(refit, "error") && length(.dup_coef_names(refit)) > 0L
  if (bracketed) refit <- refit_on(.bracket_group_contrasts(dt, groups, "II"))
  if (inherits(refit, "error")) list(model = model, refit = FALSE) else
    list(model = refit, refit = TRUE, bracketed = bracketed)
}

#' Model-versus-null statistics from the fitted model's own deviances
#'
#' For a 0/1 response the saturated log-likelihood is zero, so the deviance is
#' -2 log L and McFadden's 1 - log L / log L0 equals 1 - D / D0. Using the
#' deviances keeps the prior weights -- aggregated counts, sampling weights --
#' exactly as the fit used them. R's binomial log-likelihood, and hence AIC and
#' BIC, rounds non-integer weights, so those are NA then.
#'
#' \code{stats::BIC()} takes \code{log(nobs())}, and \code{nobs()} of a
#' weighted glm counts rows. With frequency weights the penalty is that of the
#' observations they represent (\code{freq_n}), as on the expanded data and as
#' in \code{anova_glm()}; the log-likelihood, and so AIC, already is.
#' @param freq_n the number of observations frequency weights represent, or
#'   \code{NULL} when the weights are not frequencies.
#' @noRd
.bin_model_stats <- function(model, success, freq_n = NULL) {
  w <- model$prior.weights
  whole <- is.null(w) || all(abs(w - round(w)) < 1e-8)
  dev <- model$deviance
  dev0 <- model$null.deviance
  lr <- max(dev0 - dev, 0)
  lr_df <- as.integer(model$df.null - model$df.residual)
  out <- list(
    AIC = if (whole) stats::AIC(model) else NA_real_,
    BIC = if (!whole) NA_real_ else if (is.null(freq_n)) stats::BIC(model) else {
      ll <- stats::logLik(model)
      -2 * as.numeric(ll) + log(freq_n) * attr(ll, "df")
    },
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
    # model.matrix() names a column by pasting variable and level (inside
    # brackets when the labels were bracketed to keep them distinct); check
    # that the combinations line up with the columns before labelling them.
    expect <- function(open, close) do.call(paste, c(lapply(seq_along(vars),
      function(i) paste0(unq(vars[i]), open, combos[[i]], close)), list(sep = ":")))
    got <- unq(colnames(X)[idx])
    if (!identical(expect("", ""), got) && !identical(expect("[", "]"), got)) next
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
#' @param ci_method \code{"profile"}, \code{"wald"}, or \code{"separated"}:
#'   Wald intervals except for the coefficients separation diverges, which
#'   get their profile intervals computed directly.
#' @param separated whether separation was detected. A coefficient it drives
#'   to infinity then gets its profile interval computed directly.
#' @noRd
.odds_ratios <- function(model, conf_level = 0.95, ci_method = "profile",
                         vcov_matrix = NULL, separated = FALSE) {
  ct <- .coef_table(model, vcov_matrix)
  notes <- character(0)
  # Profiles are computed on the same model fitted to the data collapsed to
  # one row per cell and outcome: the likelihood, and so every profile, is
  # the same, and the refits no longer run over every row.
  pm <- if (ci_method %in% c("profile", "separated")) .collapsed_bin_fit(model)
  if (ci_method == "profile") {
    ci <- tryCatch(
      suppressWarnings(suppressMessages(
        stats::confint(pm, level = conf_level))),
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
    # Each profile refits the model several times, so the cost grows with
    # the number of coefficients times the cost of a fit; past about 30 that
    # is seconds to minutes, which the caller should know is avoidable.
    n_prof <- sum(ct$term != "(Intercept)")
    if (n_prof > 30L) {
      notes <- c(notes, sprintf(
        "Profile-likelihood intervals were computed for %d odds ratios. Each refits the model several times, which is slow with this many coefficients; ci_method = \"wald\" is much faster, and adequate when there is no separation.",
        n_prof))
    }
  }
  idx <- match(ct$term, rownames(ci))
  lo <- ci[idx, 1L]
  hi <- ci[idx, 2L]

  # Under separation stats::confint() profiles from the diverged estimate with
  # steps scaled by its enormous standard error, and returns bounds that are
  # meaningless on both sides. A standard error above 10 on the log-odds scale
  # does not occur for a coefficient the data identify.
  method <- rep(method_used, nrow(ct))
  if ((method_used == "profile" || ci_method == "separated") && separated) {
    flag <- which(ct$term != "(Intercept)" & is.finite(ct$estimate) &
                    is.finite(ct$se) & ct$se > 10)
    method[flag] <- "profile"
    for (k in flag) {
      b <- .separated_profile(pm, ct$term[k], conf_level)
      lo[k] <- b[1L]
      hi[k] <- b[2L]
    }
    if (length(flag) > 0L) {
      notes <- c(notes, sprintf(
        "Separation drives the odds ratio(s) for %s to infinity or to zero, so the profile-likelihood interval is open on that side and is reported as conf_high = Inf (or conf_low = 0). Its finite end was found by profiling the likelihood directly, as stats::confint() cannot step from a diverged estimate%s%s.",
        paste(gsub("`", "", ct$term[flag], fixed = TRUE), collapse = ", "),
        .unidentified_clause(ct$term[flag], lo[flag], hi[flag], "odds ratio"),
        if (anyNA(c(lo[flag], hi[flag]))) "; it is NA where that search failed" else ""))
    }
  }

  if (method_used == "profile") {
    # A diverged odds ratio has its own separation note.
    keep <- !(separated & is.finite(ct$se) & ct$se > 10)
    notes <- c(notes, .profile_wald_note(ct$term[keep], lo[keep], hi[keep],
                                         ct$p_value[keep], conf_level,
                                         "odds ratio of 1"))
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
    ci_method  = method,
    stringsAsFactors = FALSE, row.names = NULL
  )
  out <- out[out$term != "(Intercept)", , drop = FALSE]
  row.names(out) <- NULL
  list(table = out, note = notes)
}

#' The same logistic model fitted to its data collapsed by cell and outcome
#'
#' The model has only grouping factors, so rows with the same cell and the
#' same outcome differ in nothing but their weight, and the binomial
#' log-likelihood sum(w * log p) is unchanged when they are merged into one
#' row carrying the sum of their weights (whole or not). The collapsed fit
#' therefore has the same coefficients, standard errors, deviance and profile
#' likelihoods, at a cost set by the number of cells rather than of rows. The
#' model itself is returned when collapsing would not save much, when it has
#' an offset, or when the refit does not reproduce it.
#' @noRd
.collapsed_bin_fit <- function(model) {
  mf <- model$model
  if (is.null(mf) || !is.null(stats::model.offset(mf))) return(model)
  w <- stats::model.weights(mf)
  if (is.null(w)) w <- rep(1, nrow(mf))
  keep <- setdiff(names(mf), "(weights)")
  key <- do.call(paste, c(lapply(mf[keep], function(x) as.integer(factor(x))),
                          sep = "\r"))
  first <- !duplicated(key)
  if (sum(first) > 0.5 * nrow(mf)) return(model)
  agg <- mf[first, keep, drop = FALSE]
  wsum <- rowsum(as.numeric(w), key, reorder = FALSE)
  fit <- tryCatch(suppressWarnings(do.call(stats::glm, list(
    formula = stats::formula(model), family = stats::family(model),
    data = agg, weights = as.numeric(wsum[key[first], 1L]),
    contrasts = model$contrasts, control = model$control))),
    error = function(e) NULL)
  if (is.null(fit)) return(model)
  b0 <- stats::coef(model)
  b1 <- stats::coef(fit)
  se0 <- sqrt(pmax(diag(stats::vcov(model)), 0))
  ok <- identical(names(b0), names(b1)) && identical(is.na(b0), is.na(b1)) &&
    abs(fit$deviance - model$deviance) <= 1e-6 * max(1, model$deviance)
  if (ok) {
    # Coefficients that separation diverges stop wherever the iterations did;
    # the others are the same estimates.
    id <- !is.na(b0) & is.finite(se0[names(b0)]) & se0[names(b0)] <= 10
    ok <- isTRUE(all.equal(unname(b0[id]), unname(b1[id]), tolerance = 1e-6))
  }
  if (ok) fit else model
}

#' The part of a separation note about intervals open on both sides
#'
#' A coefficient whose profile likelihood rules out no value in either
#' direction is not identified at all: its interval is (0, Inf) and its point
#' estimate is wherever the fit happened to stop.
#' @param lo,hi the interval ends on the linear-predictor scale.
#' @noRd
.unidentified_clause <- function(terms, lo, hi, what) {
  both <- !is.na(lo) & !is.na(hi) & lo == -Inf & hi == Inf
  if (!any(both)) return("")
  sprintf("; the interval for %s is open on both sides (0 to Inf), because the likelihood rules out no value either way: the data do not identify that %s at all, and its point estimate is arbitrary",
          paste(gsub("`", "", terms[both], fixed = TRUE), collapse = ", "), what)
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
#' @param design optional list with \code{X} and \code{coef}: the same fit on
#'   another coding (see the code).
#' @return the two ends on the linear-predictor scale: +/-Inf for an open end,
#'   NA for a finite end that cannot be found. The far end is open too when
#'   the profile stays below the cut-off on that side as well (a coefficient
#'   the data do not identify at all). The \code{"flat"} attribute (lower end,
#'   upper end) marks an end that is open because the profile stayed flat out
#'   to the search limit, rather than because the estimate diverged that way.
#' @noRd
.separated_profile <- function(model, term, conf_level,
                               cutoff = stats::qchisq(conf_level, 1),
                               dispersion = 1, design = NULL) {
  # `design` reparameterises the same fit (a treatment-coded design matrix
  # `X` spanning the model's column space, with the coefficients `coef` it
  # gives): the likelihood, and so the profile, are unchanged, but `term`
  # can then be a coefficient the model itself was not fitted with.
  X <- if (is.null(design)) stats::model.matrix(model) else design$X
  co <- if (is.null(design)) stats::coef(model) else design$coef
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
    refit <- function(off, s) {
      f <- tryCatch(suppressWarnings(stats::glm.fit(
        Xo, y, weights = w, offset = off, family = fam, start = s,
        control = stats::glm.control(maxit = 100L))), error = function(e) NULL)
      if (is.null(f) || !isTRUE(f$converged) || !is.finite(f$deviance)) return(NULL)
      # glm.fit()'s convergence flag says only that the deviance stopped
      # changing. Far from the estimate, started from the data, it can stop at
      # a deviance hundreds of units above the minimum, and whether it does
      # depends on the platform. The likelihood is concave, so a refit is at
      # the minimum exactly when its Newton decrement (in deviance units) is
      # zero; one that is not cannot mark where the profile crosses.
      mu <- f$fitted.values
      me <- fam$mu.eta(f$linear.predictors)
      va <- fam$variance(mu)
      if (!all(is.finite(c(me, va))) || any(va <= 0)) return(NULL)
      score <- crossprod(Xo, w * (y - mu) * me / va)
      info <- crossprod(Xo, Xo * (w * me^2 / va))
      dec <- tryCatch(drop(crossprod(score, solve(
        info + diag(1e-10 * max(1, diag(info)), ncol(Xo)), score))),
        error = function(e) Inf)
      if (!is.finite(dec) || dec > 1e-3) NULL else f
    }
  }
  # The search moves out from the estimate in steps. Each refit is tried from
  # the solution at the previous point (continuation, as profile.glm() does),
  # from the fitted coefficients and from the data, and the lowest deviance
  # of the refits that reached their minimum is the profile there.
  fit_at <- function(b, s) {
    off <- off0 + b * X[, j]
    if (logit) {
      return(list(gap = (dev_at(off) - dmin) / dispersion - cutoff, coef = s))
    }
    fits <- list(refit(off, s), if (!identical(s, start)) refit(off, start),
                 refit(off, NULL))
    fits <- Filter(Negate(is.null), fits)
    if (length(fits) == 0L) return(list(gap = NA_real_, coef = s))
    f <- fits[[which.min(vapply(fits, `[[`, numeric(1), "deviance"))]]
    cf <- unname(f$coefficients)
    list(gap = (f$deviance - dmin) / dispersion - cutoff,
         coef = if (anyNA(cf)) s else cf)
  }
  # A profile that stays below the cut-off out to 32 or more units from the
  # estimate on the linear-predictor scale (a ratio of e^32, about 1e14) rules
  # out no value on that side: the end is open, not missing, whichever way
  # the arbitrary estimate of an unidentified coefficient happens to point.
  # NA is kept for a search that failed (no finite profile that far out).
  find_end <- function(dir) {
    prev <- est
    s_prev <- start
    step <- 0.5
    reached <- 0
    while (step <= 128) {
      b <- est + dir * step
      at <- fit_at(b, s_prev)
      g <- at$gap
      if (is.finite(g) && g > 0) {
        s0 <- s_prev
        r <- tryCatch(suppressWarnings(stats::uniroot(
          function(x) fit_at(x, s0)$gap, sort(c(prev, b)), tol = 1e-8))$root,
          error = function(e) NA_real_)
        return(structure(r, flat = FALSE))
      }
      if (is.finite(g)) {
        reached <- step
        s_prev <- at$coef
      }
      prev <- b
      step <- step * 2
    }
    if (reached >= 32) structure(dir * Inf, flat = TRUE) else
      structure(NA_real_, flat = FALSE)
  }
  up <- if (est >= 0) 1 else -1           # the direction it diverged in
  open_end <- find_end(up)
  far_end <- find_end(-up)
  far_flat <- isTRUE(attr(far_end, "flat"))
  open_end <- as.numeric(open_end)
  far_end <- as.numeric(far_end)
  if (is.na(open_end)) open_end <- up * Inf
  if (up > 0) {
    structure(c(far_end, open_end), flat = c(far_flat, FALSE))
  } else {
    structure(c(open_end, far_end), flat = c(FALSE, far_flat))
  }
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
#'
#' Frequency weights count observations and are used as they are. Other
#' weights have an arbitrary scale -- weights averaging 0.01 would make 100
#' rows look like one expected event -- so they are rescaled to average 1,
#' which keeps each row one unit, as the dispersion of anova_count() does.
#' @param frequencies are the weights frequencies (see .frequency_weights())?
#' @noRd
.sparse_bin_note <- function(d, response, groups, weights, additive,
                             threshold = 5, frequencies = FALSE) {
  w <- if (is.null(weights)) rep(1, nrow(d)) else as.numeric(d[[weights]])
  if (!frequencies) w <- w / mean(w)
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
