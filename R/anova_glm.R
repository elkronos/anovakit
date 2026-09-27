#' Analysis of Deviance for a Generalised Linear Model
#'
#' Fits a generalised linear model of a response on one or more grouping
#' variables in any \code{\link[stats]{family}}, and reports an analysis of
#' deviance table, effect sizes, estimated marginal means and pairwise
#' comparisons.
#'
#' @details
#' \strong{Type III sums of squares.} When \code{type = "III"} the model is
#' \emph{fitted} under sum-to-zero contrasts. This is the whole point: setting
#' \code{options(contrasts = )} after a model has been fitted does not change
#' the contrasts stored on it, so \code{car::Anova(type = 3)} would silently
#' report simple effects at the reference level under the label "Type III".
#'
#' \strong{Post-hoc.} Marginal means and comparisons come from \pkg{emmeans}
#' and are on the response scale. When the grouping variables enter the model
#' additively (the default with more than one), they are reported for each
#' factor on its own, averaged over the others; with an interaction they are
#' made across the cells of the grid. With a log or logit link a comparison is
#' a ratio (of means, or of odds); with any other non-identity link (inverse,
#' probit, square root, ...) it is a difference between the marginal means on
#' the response scale. Families that estimate their dispersion use t with the
#' residual degrees of freedom, as the coefficient table does.
#'
#' \strong{Separation.} For binomial, quasi-binomial, Poisson, quasi-Poisson
#' and negative binomial families the fit is checked for groups whose fitted
#' probability is numerically 0 or 1 (or whose fitted rate is numerically 0),
#' where coefficients, Wald tests and comparisons are not identified;
#' \code{$notes} says so when it happens.
#'
#' For binary responses see \code{\link{anova_bin}}, and for counts see
#' \code{\link{anova_count}}: both add family-specific diagnostics and effect
#' sizes that this general wrapper does not.
#'
#' @param data A data frame, or anything inheriting from one, such as a
#'   \code{data.table} or a tibble.
#' @param response Character. Name of the response column. For a binomial or
#'   quasi-binomial family it may be 0/1, logical, a proportion (with
#'   \code{weights} giving the numbers of trials), or a factor or character
#'   column with exactly two observed values; the model is then for the
#'   second level (in C-locale order for a character column), and
#'   \code{$notes} names it.
#' @param groups Character vector. One or more grouping columns.
#' @param family A \code{\link[stats]{family}} object, a family function, or the
#'   name of one as a character string. Default \code{stats::gaussian()}.
#' @param interaction \code{FALSE} (additive, the default), \code{TRUE} (full
#'   factorial), or a whole number giving the highest interaction order.
#' @param type Character. \code{"II"} (default) or \code{"III"} sums of squares.
#' @param test_statistic Character or \code{NULL}. One of \code{"LR"},
#'   \code{"Wald"} or \code{"F"}, passed to \code{\link[car]{Anova}}. When
#'   \code{NULL} (the default) it is chosen after fitting: an F test when the
#'   fitted model estimates its dispersion -- the Gaussian, Gamma, inverse
#'   Gaussian, quasi, quasi-Poisson and quasi-binomial families, and
#'   \code{MASS::negative.binomial(theta)}, whose coefficient table also uses
#'   t -- and a likelihood ratio test when the dispersion is fixed (binomial,
#'   Poisson).
#' @param conf_level Numeric in (0, 1). Level for every interval returned.
#'   Default \code{0.95}.
#' @param ci_method Character. \code{"profile"} (default) or \code{"wald"} for
#'   the coefficient intervals. Both use the reference distribution of the
#'   coefficient table's p-values. For a family with a fixed dispersion
#'   (binomial, Poisson) the profile interval is \code{stats::confint()}'s,
#'   which cuts the signed-root deviance profile at normal quantiles. For a
#'   family that estimates its dispersion the same profile is cut at
#'   \code{qt(., df.residual)} instead, so that the interval agrees with the
#'   t test beside it; for the Gaussian identity model it equals
#'   \code{confint(lm())}. If profiling fails, Wald intervals are used, with a
#'   note. A bound the profile cannot reach (as under separation) is
#'   \code{NA}, also with a note.
#' @param vcov_type Character. \code{"model"} (default) or an HC type
#'   (\code{"HC0"} to \code{"HC4"}), which requires the \pkg{sandwich} package.
#'   The robust covariance is used by the coefficient table (with Wald
#'   intervals), the marginal means and the comparisons. The omnibus
#'   likelihood ratio and F tests are computed from the likelihood and cannot
#'   use it; with \code{test_statistic = "Wald"} the omnibus test uses it too.
#'   \code{$notes} says which applies.
#' @param adjust Character. Multiplicity adjustment for the pairwise
#'   comparisons, applied by \pkg{emmeans}. One of \code{"tukey"},
#'   \code{"sidak"}, \code{"scheffe"}, \code{"dunnettx"}, \code{"bonferroni"},
#'   \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"BH"},
#'   \code{"BY"}, \code{"fdr"} or \code{"none"}. Default \code{"tukey"}.
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
#'   dropped for missing values. Rows with weight zero contribute nothing to
#'   the fit; they are dropped before fitting and counted in
#'   \code{$n_removed}. Whole-number weights (some above 1) on a count or 0/1
#'   response are read as frequency weights, each row standing for that many
#'   identical observations: the robust covariance (\code{vcov_type}) and, for
#'   the Poisson families, the residual degrees of freedom are then those of
#'   the data expanded to one row per observation. Other weights (precision
#'   weights, or the trials behind a proportion) keep each row as one unit.
#'
#' @return An \code{\link{anovakit_fit}} object, with \code{$family} (the
#'   family object actually used) and \code{$model_stats} (AIC, BIC, the null
#'   and residual deviances, and the residual degrees of freedom).
#'   \code{$assumptions} holds \code{coefficients} -- the coefficient table,
#'   with the intervals \code{ci_method} selects -- and \code{dispersion}.
#'   \code{$effect_sizes} depends on the family:
#'   \itemize{
#'     \item Gaussian: partial eta squared and partial omega squared for each
#'       term, from the F table's sums of squares whatever
#'       \code{test_statistic} is. Partial omega squared is floored at 0, with
#'       a note. When the ANOVA table cannot be computed (a Type III table
#'       with an empty cell, say) the table is empty and a note says why.
#'     \item Any other family: \code{deviance_explained},
#'       \code{1 - deviance / null deviance}, and \code{mcfadden_r2},
#'       McFadden's pseudo R squared \code{1 - logLik(model) /
#'       logLik(intercept-only model)}. McFadden's measure is reported only for
#'       the binomial, Poisson and negative binomial families, whose
#'       log-likelihoods are probabilities; for a binomial response it is
#'       computed from the Bernoulli log-likelihood, so it is the same whether
#'       the data are given one row per trial or aggregated as proportions
#'       with trial weights (and for a 0/1 response it equals the deviance
#'       explained). For the Gamma, inverse Gaussian and quasi families
#'       it is \code{NA} with a note: a density's log-likelihood changes with
#'       the units of the response, and a quasi family has none. The deviance
#'       explained is measured against the saturated model, so for binomial
#'       data it does depend on the aggregation.
#'   }
#'
#' @seealso \code{\link{anova_bin}} and \code{\link{anova_count}}, which add
#'   family-specific diagnostics and effect sizes.
#'
#' @examples
#' set.seed(42)
#' d <- data.frame(
#'   A = factor(sample(c("a1", "a2"), 200, replace = TRUE)),
#'   B = factor(sample(c("b1", "b2"), 200, replace = TRUE))
#' )
#' d$y <- 3 + 2 * (d$A == "a2") + (d$B == "b2") +
#'   4 * (d$A == "a2") * (d$B == "b2") + rnorm(200)
#'
#' # Type III tests, with the model fitted under sum-to-zero contrasts
#' fit <- anova_glm(d, "y", c("A", "B"), interaction = TRUE, type = "III")
#' fit$anova
#'
#' # The coefficient table, with profile-likelihood intervals by default
#' fit$assumptions$coefficients
#'
#' # A family may be given as a string, as stats::glm() allows. For a binary
#' # response anova_bin() is usually the better entry point.
#' d$hit <- rbinom(200, 1, plogis(-0.5 + 1.2 * (d$A == "a2")))
#' bin <- anova_glm(d, "hit", "A", family = "binomial", plots = FALSE)
#' bin$anova
#' bin$effect_sizes
#'
#' # Heteroskedasticity-consistent standard errors, which flow through to the
#' # marginal means and the comparisons
#' anova_glm(d, "y", "A", vcov_type = "HC3", plots = FALSE)$emmeans
#'
#' @export
anova_glm <- function(data, response, groups,
                      family = stats::gaussian(),
                      interaction = FALSE,
                      type = c("II", "III"),
                      test_statistic = NULL,
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
  ci_method <- match.arg(ci_method)
  # Validated here rather than left to car::Anova(), whose error would be
  # caught and turned into a note on a fit with no omnibus table.
  if (!is.null(test_statistic)) {
    .check_adjust(test_statistic, c("LR", "Wald", "F"), "test_statistic")
  }
  .check_adjust(vcov_type, .vcov_choices(), "vcov_type")
  .check_adjust(adjust, .adjust_choices(), "adjust")
  family <- .as_family(family)
  data <- .as_df(data)
  .check_name(response, "response")
  .check_names(groups, "groups")
  .check_columns(data, response, "Response column")
  .check_columns(data, groups, "Grouping column(s)")
  .check_conf_level(conf_level)
  .check_interaction(interaction, length(groups))
  .check_roles(list(`the response` = response, `a grouping variable` = groups,
                    `the weights` = weights))
  # The column's type is checked on the input; the levels and values only on
  # the rows that are analysed, after missing values and zero weights are gone.
  .check_response_for_family(data, response, family, values = FALSE)
  .check_flag(posthoc, "posthoc")
  .check_flag(plots, "plots")
  .check_flag(verbose, "verbose")

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
  .check_response_for_family(d, response, family)
  if (family$family %in% c("binomial", "quasibinomial") &&
      (is.factor(d[[response]]) || is.character(d[[response]]))) {
    # glm() models the probability of the second level of a factor response.
    # A character response is not accepted by glm() at all; it becomes a
    # factor with its levels in C-locale order, as grouping columns do.
    was_char <- is.character(d[[response]])
    if (was_char) d[[response]] <- .as_group_factor(d[[response]], response)
    d[[response]] <- droplevels(d[[response]])
    notes <- c(notes, sprintf(
      "Response `%s` is %s: the model is for the probability that it equals \"%s\" (its second level), against \"%s\".",
      response,
      if (was_char) "a character column, taken as a factor with its values in sorted order" else "a factor",
      levels(d[[response]])[2L], levels(d[[response]])[1L]))
  }

  terms_rhs <- .group_terms(groups, interaction)
  fml <- stats::reformulate(terms_rhs, response = .bq(response))
  .say(verbose, "Fitting %s with family %s.",
       paste(deparse(fml), collapse = " "), family$family)

  fitted <- .collect_conditions(.fit_with_contrasts(
    function() .fit_glm(fml, d, family, weights),
    type = type))
  model <- fitted$value
  if (inherits(model, "error")) {
    .stopf("The model could not be fitted: %s", conditionMessage(model))
  }
  notes <- c(notes, .said_note(fitted$said, "stats::glm()"))
  if (!isTRUE(model$converged)) {
    notes <- c(notes, "The model did not converge; every result below is unreliable.")
  }
  fam_name <- family$family
  # A Poisson frequency table (whole-number weights, some above 1) holds one
  # row per distinct count. glm() counts residual degrees of freedom in rows;
  # the data the table stands for have them in observations. Corrected on the
  # fit, so the dispersion, F tests, t intervals and marginal means all equal
  # those of the expanded data, as in anova_count().
  freq_n <- if (fam_name %in% c("poisson", "quasipoisson") &&
                .frequency_weights(d, weights)) sum(d[[weights]]) else NULL
  row_model <- model
  if (!is.null(freq_n) && freq_n - model$rank >= 1) {
    model$df.residual <- freq_n - model$rank
  }
  notes <- c(notes, .check_model_size(model))
  sep_note <- if (fam_name %in% c("binomial", "quasibinomial", "poisson", "quasipoisson") ||
                  grepl("^Negative Binomial", fam_name)) {
    .check_separation(model)
  } else character(0)
  notes <- c(notes, sep_note)

  # The default follows what the fitted model does with its dispersion, which
  # is also what decides t or z in the coefficient table: a family string
  # cannot tell glm(family = MASS::negative.binomial(theta)), whose dispersion
  # summary.glm() estimates, from a fixed-dispersion fit.
  if (is.null(test_statistic)) {
    test_statistic <- if (.estimates_dispersion(model)) "F" else "LR"
  }

  ## Robust covariance --------------------------------------------------------
  # Whole-number weights on a count or 0/1 response count identical
  # observations; the robust covariance is then that of the expanded data.
  # Other weights (precision weights, trials behind a proportion) keep the
  # row as the unit.
  freq_rows <- .frequency_weights(d, weights) && (
    !is.null(freq_n) || grepl("^Negative Binomial", fam_name) ||
      (fam_name %in% c("binomial", "quasibinomial") && all(model$y %in% c(0, 1))))
  rv <- .robust_vcov(model, vcov_type, freq_weights = freq_rows)
  notes <- c(notes, rv$note)

  ## Analysis of deviance -----------------------------------------------------
  av <- .glm_anova(row_model, type = type, test_statistic = test_statistic,
                   vcov_matrix = rv$matrix)
  if (!is.null(freq_n) && !(identical(test_statistic, "Wald") && !is.null(rv$matrix))) {
    av <- .freq_weight_anova(av, stats::df.residual(row_model),
                             stats::df.residual(model),
                             .estimates_dispersion(model) || identical(test_statistic, "F"))
  }
  notes <- c(notes, av$note)
  if (!is.null(rv$matrix) && !is.null(av$table)) {
    notes <- c(notes, .omnibus_vcov_note(test_statistic, vcov_type,
                                         identical(test_statistic, "Wald")))
  }

  ## Effect sizes -------------------------------------------------------------
  eff <- .glm_effect_sizes(model, av$raw, type = type,
                           test_statistic = test_statistic)
  notes <- c(notes, eff$notes)

  ## Dispersion ---------------------------------------------------------------
  dispersion <- .dispersion(model)
  if (!is.null(freq_n)) {
    notes <- c(notes, sprintf(
      "The weights in `%s` are whole numbers, so they are treated as frequency weights: the residual degrees of freedom are counted in the %s observations they represent, not the %d rows, so the dispersion, F tests and intervals equal those of the data expanded to one row per observation.",
      weights, format(freq_n), nrow(d)))
  }
  # A binary response cannot be overdispersed, and with frequency weights its
  # Pearson statistic grows with the weights: the check only means something
  # for counts and for proportions out of several trials.
  binary <- fam_name %in% c("binomial", "quasibinomial") && all(model$y %in% c(0, 1))
  if (binary) {
    # For 0/1 data the Pearson statistic says nothing about overdispersion (in a
    # one-way model it is N / (N - k) whatever the data), and with frequency
    # weights it grows with them.
    dispersion <- NA_real_
    notes <- c(notes, "The dispersion is reported as NA: for a 0/1 response the Pearson statistic carries no information about overdispersion, which binary data cannot show.")
  }
  fixed_family <- fam_name %in% c("poisson", "binomial") && !binary
  p_over <- .overdispersion_p(dispersion, stats::df.residual(model))
  if (fixed_family && is.finite(dispersion) && dispersion > 1 &&
      is.finite(p_over) && p_over < 0.05) {
    # Under a true dispersion phi a test statistic of this family is about phi
    # times its nominal chi-square; say what that does rather than only that
    # it happened.
    size <- stats::pchisq(stats::qchisq(0.95, 1) / dispersion, 1,
                          lower.tail = FALSE)
    notes <- c(notes, sprintf(
      "Pearson dispersion is %.2f (Pearson test of overdispersion: %s), where a %s model assumes 1: its test statistics are inflated by about that factor, its standard errors are about %.2f times too small, and a nominal 5%% test on 1 degree of freedom rejects roughly %.0f%% of true null hypotheses. %s",
      dispersion, .fmt_p(p_over), fam_name, sqrt(dispersion), 100 * size,
      if (identical(test_statistic, "F")) {
        sprintf("The F table already scales by the dispersion, but the coefficient table, the marginal means and the comparisons do not: use family = \"quasi%s\" to scale them too.",
                fam_name)
      } else {
        sprintf("Consider family = \"quasi%s\", or %s.", fam_name,
                if (fam_name == "poisson") "anova_count(model = \"negbin\")" else
                  "a beta-binomial model")
      }))
  }
  if (identical(test_statistic, "F") && fam_name %in% c("poisson", "binomial")) {
    notes <- c(notes, sprintf(
      "test_statistic = \"F\" on a %s fit is a quasi-likelihood test: car::Anova() estimates the dispersion from the Pearson statistic (%.3g) instead of fixing it at 1. $effect_sizes, the coefficient table, $emmeans and $posthoc still assume a dispersion of 1; family = \"quasi%s\" makes all of them consistent with the F table.",
      fam_name, .dispersion(model), fam_name))
  }

  ## Coefficient intervals ----------------------------------------------------
  if (!is.null(rv$matrix) && ci_method == "profile") {
    notes <- c(notes, sprintf(
      "Robust standard errors (%s) were requested, so the coefficient intervals are Wald intervals built from them; profile-likelihood intervals cannot use a sandwich covariance.",
      vcov_type))
  }
  ci <- .coef_intervals(model, conf_level = conf_level,
                        ci_method = if (is.null(rv$matrix)) ci_method else "wald",
                        vcov_matrix = rv$matrix, known = fitted$said,
                        separated = length(sep_note) > 0L)
  notes <- c(notes, ci$note)

  ## Marginal means and comparisons ------------------------------------------
  emm <- .emm_block(model, groups, additive = .is_additive(interaction),
                    type = "response", vcov_matrix = rv$matrix,
                    conf_level = conf_level, adjust = adjust, posthoc = posthoc,
                    link = family$link, data = d)
  notes <- c(notes, emm$notes)
  if (!is.null(emm$posthoc) &&
      !family$link %in% c("identity", "log", "logit")) {
    notes <- c(notes, sprintf(
      "With the %s link the pairwise comparisons are differences between the marginal means on the response scale; only log and logit links give ratios.",
      family$link))
  }
  if (posthoc && !is.null(emm$table) && !emm$per_factor && !is.null(emm$grid)) {
    if (NROW(emm$table) < 2L) {
      notes <- c(notes, "No pairwise comparisons: fewer than two estimable cells.")
    } else if (NROW(emm$table) == 2L && !is.null(emm$posthoc)) {
      notes <- c(notes, .two_cell_note(model, av$table, test_statistic,
                                       robust = !is.null(rv$matrix),
                                       vcov_type = vcov_type))
    }
  }

  ## Plots --------------------------------------------------------------------
  plot_list <- list()
  if (plots) {
    .say(verbose, "Building plots.")
    added <- .add_cell(d, groups)
    if (is.numeric(d[[response]])) {
      plot_list$box <- .plot_box(added$data, response, added$cell,
                                 title = "Observed response by group",
                                 xlab = paste(groups, collapse = " : "))
    }
    plot_list$residuals <- .plot_resid_fitted(
      stats::fitted(model), stats::residuals(model, type = "deviance"),
      title = "Deviance residuals vs fitted", ylab = "Deviance residuals")
    plot_list$qq <- .plot_qq(stats::residuals(model, type = "deviance"),
                             title = "Normal Q-Q plot of deviance residuals")
    if (!is.null(emm$table)) {
      plot_list$emmeans <- .plot_emmeans(
        emm$table, groups, estimate = "estimate",
        lower = "conf_low", upper = "conf_high", conf_level = conf_level,
        title = "Estimated marginal means",
        ylab = sprintf("Estimated %s", response))
    }
    if (is.null(plot_list$emmeans)) notes <- c(notes, .no_emm_plot_note())
  }

  .new_fit(
    method       = sprintf("Analysis of deviance (%s family, %s link, Type %s)",
                           family$family, family$link, type),
    call         = cl,
    model        = model,
    anova        = av$table,
    effect_sizes = eff$table,
    emmeans      = emm$table,
    emmeans_object = emm$grid,
    posthoc      = emm$posthoc,
    assumptions  = list(dispersion = dispersion,
                        coefficients = ci$table),
    plots        = plot_list,
    data_used    = d,
    n_removed    = n_removed,
    conf_level   = conf_level,
    notes        = notes,
    extra        = list(
      family = family,
      model_stats = list(AIC = .quiet_number(stats::AIC(model)),
                         BIC = .quiet_number(stats::BIC(model)),
                         null_deviance = model$null.deviance,
                         residual_deviance = model$deviance,
                         df_residual = model$df.residual))
  )
}

#' A number from an expression that may warn or fail, without console output
#' @noRd
.quiet_number <- function(expr) {
  got <- .collect_conditions(expr)
  v <- got$value
  if (inherits(v, "error") || !is.numeric(v) || length(v) != 1L) NA_real_ else
    as.numeric(v)
}

#' The analysis of deviance table, with a robust covariance for a Wald test
#'
#' A likelihood-ratio or F test is computed from deviances and has no place
#' for a sandwich covariance; a Wald test does, and car accepts one through
#' \code{vcov.}. car announces that with a message that deparses its
#' argument -- the whole matrix, when it is passed as a value -- so that
#' message is dropped here and the caller says it in words instead.
#' @noRd
.glm_anova <- function(model, type = "II", test_statistic = "LR",
                       vcov_matrix = NULL) {
  robust <- !is.null(vcov_matrix) && identical(test_statistic, "Wald")
  .car_anova(model, type = type, test_statistic = test_statistic,
             vcov_matrix = if (robust) vcov_matrix else NULL)
}

#' The note for a design with exactly two cells
#'
#' The single comparison is the omnibus test only for a Gaussian identity
#' model with model-based standard errors and an F test on one term with one
#' degree of freedom, where t squared equals F. Anything else -- a likelihood
#' ratio omnibus test against a Wald comparison, or a robust covariance on one
#' side only -- can disagree, sometimes badly (under separation, say).
#' @noRd
.two_cell_note <- function(model, anova_table, test_statistic, robust = FALSE,
                           vcov_type = "model") {
  fam <- stats::family(model)
  terms <- if (is.null(anova_table)) NULL else
    anova_table[anova_table$term != "Residuals", , drop = FALSE]
  same <- identical(fam$family, "gaussian") && identical(fam$link, "identity") &&
    !robust && identical(test_statistic, "F") && NROW(terms) == 1L &&
    isTRUE(all.equal(as.numeric(terms$df), 1))
  if (same) {
    return("With two cells the single pairwise comparison is the omnibus test; no multiplicity adjustment was applied.")
  }
  stat <- switch(test_statistic, LR = "likelihood-ratio", F = "F", Wald = "Wald")
  sprintf("With two cells there is a single pairwise comparison, so no multiplicity adjustment was applied. It is a Wald test%s and can differ from the omnibus %s test.",
          if (robust) sprintf(" with the robust (%s) covariance", vcov_type) else "",
          stat)
}

#' Effect sizes appropriate to the family
#'
#' Partial eta and omega squared for a Gaussian model, where sums of squares
#' exist. For any other family, the deviance explained and, where the
#' log-likelihood is a probability (binomial, Poisson, negative binomial),
#' McFadden's pseudo R squared.
#'
#' @return list with \code{table} and \code{notes}.
#' @noRd
.glm_effect_sizes <- function(model, anova_raw = NULL, type = "II",
                              test_statistic = "F") {
  fam <- stats::family(model)$family
  if (identical(fam, "gaussian")) {
    has_ss <- function(r) !is.null(r) && "Sum Sq" %in% names(as.data.frame(r))
    raw <- anova_raw
    # The sums of squares do not depend on the test statistic, but only car's
    # F table reports them: a likelihood-ratio or Wald table has none.
    if (!has_ss(raw) && !identical(test_statistic, "F")) {
      raw <- .car_anova(model, type = type, test_statistic = "F")$raw
    }
    if (!has_ss(raw)) {
      empty <- data.frame(term = character(0), df = numeric(0),
                          sum_sq = numeric(0), partial_eta_sq = numeric(0),
                          partial_omega_sq = numeric(0),
                          stringsAsFactors = FALSE)
      return(list(table = empty, notes = "Partial eta and omega squared could not be computed because the ANOVA table could not be computed, so $effect_sizes is empty. (McFadden's pseudo R squared is no substitute for a Gaussian model: its value depends on the units of the response.)"))
    }
    eff <- .partial_eta_squared(raw, df_error = stats::df.residual(model),
                                n_obs = stats::nobs(model))
    notes <- character(0)
    if (!is.null(attr(eff, "omega_floored"))) {
      notes <- sprintf(
        "Partial omega squared was negative for %s and has been floored at 0, which happens when an effect explains less than its degrees of freedom would by chance.",
        paste(attr(eff, "omega_floored"), collapse = ", "))
    }
    return(list(table = eff, notes = notes))
  }
  dev_explained <- if (isTRUE(is.finite(model$null.deviance)) &&
                       model$null.deviance != 0) {
    1 - model$deviance / model$null.deviance
  } else NA_real_
  mc <- .mcfadden_r2(model)
  list(table = data.frame(
    measure  = c("deviance_explained", "mcfadden_r2"),
    estimate = c(dev_explained, mc$value),
    stringsAsFactors = FALSE),
    notes = mc$note)
}

#' McFadden's pseudo R squared for a GLM without an offset
#'
#' \code{1 - logLik(model) / logLik(null)}, where the null model has only an
#' intercept. Without an offset its fitted mean is the weighted mean of the
#' response under any link, so it is computed directly rather than refitted
#' (a refit repeats every warning the fit gave, straight to the console).
#'
#' The ratio is only meaningful when the log-likelihood has no additive
#' constant unrelated to the fit. That holds for the binomial family when it
#' is written as a Bernoulli likelihood -- \code{logLik()} adds
#' \code{log(choose(n, k))} for aggregated data, which makes the measure
#' depend on how identical data are grouped -- and for the Poisson and
#' negative binomial, whose log-likelihoods are log probabilities. For a
#' continuous family the log-likelihood is a log density that shifts with the
#' units of the response, and a quasi family has none: both are \code{NA}.
#'
#' @return list with \code{value} and \code{note}.
#' @noRd
.mcfadden_r2 <- function(model) {
  fam <- stats::family(model)
  fname <- fam$family
  nb <- grepl("^Negative Binomial", fname)
  if (!fname %in% c("binomial", "poisson") && !nb) {
    why <- if (grepl("^quasi", fname)) {
      "a quasi family has no likelihood. The deviance explained is reported instead"
    } else {
      sprintf("the log-likelihood of a %s family is a log density, which shifts with the units of the response, so the ratio has no unit-free meaning. The deviance explained does not depend on the units and is reported instead", fname)
    }
    return(list(value = NA_real_, note = sprintf(
      "McFadden's pseudo R squared is not reported (NA): %s.", why)))
  }
  y <- as.numeric(model$y)
  w <- as.numeric(model$prior.weights)
  mu <- as.numeric(stats::fitted(model))
  if (length(y) != length(mu) || length(w) != length(mu) || sum(w) <= 0) {
    return(list(value = NA_real_, note = character(0)))
  }
  mu0 <- sum(w * y) / sum(w)
  note <- character(0)
  ll <- if (fname == "binomial") {
    if (any(y > 0 & y < 1)) {
      note <- "The deviance explained is measured against the saturated model, so for binomial data given as proportions with trial weights it is larger than for the same data given one row per trial. McFadden's pseudo R squared is computed from the Bernoulli log-likelihood and is the same either way."
    }
    function(m) {
      m <- rep_len(m, length(y))
      sum(w * (ifelse(y > 0, y * log(m), 0) + ifelse(y < 1, (1 - y) * log1p(-m), 0)))
    }
  } else if (fname == "poisson") {
    if (any(abs(y - round(y)) > 1e-8)) {
      return(list(value = NA_real_, note = "McFadden's pseudo R squared is not reported (NA): the response is not integer-valued, and the Poisson log-likelihood is only defined for counts."))
    }
    function(m) sum(w * stats::dpois(round(y), rep_len(m, length(y)), log = TRUE))
  } else {
    function(m) -fam$aic(y, rep(1, length(y)), rep_len(m, length(y)), w, 0) / 2
  }
  got <- .collect_conditions(1 - ll(mu) / ll(mu0))
  v <- got$value
  if (inherits(v, "error") || !is.numeric(v) || !is.finite(v)) v <- NA_real_
  list(value = v, note = note)
}

#' Coefficient estimates with intervals
#'
#' Profile intervals cut the signed-root deviance profile of each coefficient
#' at the quantiles of the distribution its test uses. \code{stats::confint()}
#' always uses normal quantiles, which is right when the dispersion is fixed
#' (binomial, Poisson) but too narrow when it is estimated: the profile is
#' then scaled by the estimated dispersion and its reference distribution is
#' t on the residual degrees of freedom -- exactly so for the Gaussian
#' identity model, where the result is \code{confint(lm())}.
#' @noRd
.coef_intervals <- function(model, conf_level = 0.95, ci_method = "profile",
                            vcov_matrix = NULL, known = character(0),
                            separated = FALSE) {
  ct <- .coef_table(model, vcov_matrix)
  notes <- character(0)
  ci <- NULL
  open <- character(0)
  crit <- .coef_crit(model, conf_level)
  # A diverged coefficient cannot be profiled by stats::profile() -- in a
  # Poisson model the attempt fails for every coefficient -- and gets its own
  # direct profile below, so it is left out here.
  diverged <- if (separated) {
    ct$term[is.finite(ct$se) & ct$se > 10]
  } else character(0)
  if (ci_method == "profile" && is.finite(crit)) {
    got <- .collect_conditions(.profile_ci(model, conf_level, crit,
                                           skip = diverged))
    tried <- got$value
    if (inherits(tried, "error")) {
      notes <- c(notes, sprintf(
        "Profile-likelihood intervals failed, Wald intervals were used instead: %s",
        conditionMessage(tried)))
    } else {
      ci <- tried
      # Refitting repeats whatever the fit itself said (non-integer counts,
      # say); only what is new is worth a note.
      notes <- c(notes, .said_note(setdiff(got$said, known), "Profiling the likelihood"))
      est <- stats::coef(model)
      open <- rownames(ci)[!is.na(est[rownames(ci)]) &
                             (is.na(ci[, 1L]) | is.na(ci[, 2L]))]
    }
  }
  method_used <- "profile"
  if (is.null(ci)) {
    ci <- cbind(ct$estimate - crit * ct$se, ct$estimate + crit * ct$se)
    rownames(ci) <- ct$term
    method_used <- "wald"
  }
  # Under separation (or a zero-count cell) a coefficient has run off towards
  # infinity, and confint()'s spline extrapolates the profile into junk on
  # both sides. The side it diverged in is open; the other end is found by
  # profiling the likelihood directly.
  if (separated) {
    diverged <- intersect(diverged, rownames(ci))
    for (tm in diverged) {
      # A Wald end is as meaningless as the Wald standard error behind it, and
      # confint()'s spline through a diverged profile is junk: both ends come
      # from profiling the likelihood directly, with the cutoff and dispersion
      # the other coefficients' profile intervals use.
      disp <- if (.estimates_dispersion(model)) summary(model)$dispersion else 1
      ci[tm, ] <- .separated_profile(model, tm, conf_level, cutoff = crit^2,
                                     dispersion = disp)
    }
    if (length(diverged) > 0L) {
      notes <- c(notes, sprintf(
        "The estimate of %s has diverged (separation or a zero-count cell), so its interval is open on that side (Inf or -Inf). The finite end was found by profiling the likelihood directly, as stats::confint() cannot step from a diverged estimate%s.",
        paste(gsub("`", "", diverged, fixed = TRUE), collapse = ", "),
        if (anyNA(ci[diverged, ])) "; it is NA where that search failed" else ""))
    }
    open <- setdiff(open, diverged)
  }
  if (length(open) > 0L) {
    notes <- c(notes, sprintf(
      "The profile-likelihood interval of %s is unbounded on at least one side (NA): the profile never reaches the cutoff there, which happens when an estimate runs off to the boundary of the parameter space (see any separation note), and the Wald test of such a coefficient is unreliable too.",
      paste(gsub("`", "", open, fixed = TRUE), collapse = ", ")))
  }
  idx <- match(ct$term, rownames(ci))
  out <- data.frame(ct, conf_low = ci[idx, 1L], conf_high = ci[idx, 2L],
                    ci_method = method_used, stringsAsFactors = FALSE,
                    row.names = NULL)
  list(table = out, note = notes)
}

#' Profile-likelihood intervals cut at a given critical value
#'
#' \code{profile()} on a GLM steps each coefficient away from its estimate,
#' refits the others, and records the signed root of the deviance increase
#' divided by the dispersion (Pearson, as \code{summary.glm()} uses). The
#' interval is where a spline through that profile crosses \code{-crit} and
#' \code{crit}, which is what \code{MASS}'s \code{confint.profile.glm()} does
#' with \code{crit = qnorm()}. The profile is computed with
#' \code{alpha = (1 - conf_level) / 4}, as \code{confint()} does, which for
#' an estimated dispersion carries it out to \code{qt(1 - alpha / 2)} and so
#' past the t cutoff.
#' @param skip coefficients not to profile (left NA).
#' @return a two-column matrix with one row per coefficient.
#' @noRd
.profile_ci <- function(model, conf_level, crit, skip = character(0)) {
  if (getRversion() < "4.4.0" && !requireNamespace("MASS", quietly = TRUE)) {
    .stopf("profile-likelihood intervals for a GLM need the MASS package.")
  }
  pn <- names(stats::coef(model))
  ci <- matrix(NA_real_, nrow = length(pn), ncol = 2L,
               dimnames = list(pn, c("lower", "upper")))
  which <- which(!is.na(stats::coef(model)) & !pn %in% skip)
  if (length(which) == 0L) return(ci)
  prof <- stats::profile(model, which = which, alpha = (1 - conf_level) / 4)
  for (nm in intersect(names(prof), pn)) {
    pro <- prof[[nm]]
    if (is.null(pro) || NROW(pro) < 2L) next
    pv <- pro[["par.vals"]]
    x <- if (is.matrix(pv)) pv[, nm] else pv
    z <- pro[[1L]]
    sp <- stats::spline(x = x, y = z)
    ci[nm, ] <- stats::approx(sp$y, sp$x, xout = c(-crit, crit))$y
  }
  ci
}
