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
#' \strong{Post-hoc.} Comparisons come from \pkg{emmeans} and are available for
#' every design, including factorial ones. When more than one grouping variable
#' is present, comparisons are made across the cells of the full grid.
#'
#' For binary responses see \code{\link{anova_bin}}, and for counts see
#' \code{\link{anova_count}}: both add family-specific diagnostics and effect
#' sizes that this general wrapper does not.
#'
#' @param data A data frame, or anything inheriting from one, such as a
#'   \code{data.table} or a tibble.
#' @param response Character. Name of the response column.
#' @param groups Character vector. One or more grouping columns.
#' @param family A \code{\link[stats]{family}} object, a family function, or the
#'   name of one as a character string. Default \code{stats::gaussian()}.
#' @param interaction \code{FALSE} (additive, the default), \code{TRUE} (full
#'   factorial), or a whole number giving the highest interaction order.
#' @param type Character. \code{"II"} (default) or \code{"III"} sums of squares.
#' @param test_statistic Character or \code{NULL}. Passed to
#'   \code{\link[car]{Anova}}. When \code{NULL} (the default) an F test is used
#'   for the Gaussian, Gamma, inverse Gaussian, quasi, quasi-Poisson and
#'   quasi-binomial families -- every family whose dispersion is estimated --
#'   and a likelihood ratio test otherwise.
#' @param conf_level Numeric in (0, 1). Level for every interval returned.
#'   Default \code{0.95}.
#' @param ci_method Character. \code{"profile"} (default) or \code{"wald"} for
#'   the coefficient intervals.
#' @param vcov_type Character. \code{"model"} (default) or an HC type
#'   (\code{"HC0"} to \code{"HC4"}), which requires the \pkg{sandwich} package.
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
#'   dropped for missing values.
#'
#' @return An \code{\link{anovakit_fit}} object, with \code{$family} (the
#'   family object actually used) and \code{$model_stats} (AIC, BIC, the null
#'   and residual deviances, and the residual degrees of freedom).
#'   \code{$assumptions} holds \code{coefficients} -- the coefficient table,
#'   with the intervals \code{ci_method} selects -- and \code{dispersion}.
#'   \code{$effect_sizes} is partial eta squared and partial omega squared for a
#'   Gaussian family, and the deviance explained with McFadden's pseudo R
#'   squared otherwise.
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
  .check_response_for_family(data, response, family)
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
  if (family$family %in% c("binomial", "quasibinomial") &&
      (is.factor(d[[response]]) || is.character(d[[response]]))) {
    # glm() models the probability of the second level of a factor response.
    # A character response is not accepted by glm() at all; it becomes a
    # factor with its levels in C-locale order, as grouping columns do.
    if (is.character(d[[response]])) d[[response]] <- .as_group_factor(d[[response]], response)
    d[[response]] <- droplevels(d[[response]])
    notes <- c(notes, sprintf(
      "Response `%s` is a factor: the model is for the probability that it equals \"%s\" (its second level), against \"%s\".",
      response, levels(d[[response]])[2L], levels(d[[response]])[1L]))
  }

  if (is.null(test_statistic)) {
    test_statistic <- if (family$family %in% c("gaussian", "quasipoisson",
                                               "quasibinomial", "quasi",
                                               "Gamma", "inverse.gaussian")) {
      "F"
    } else {
      "LR"
    }
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
  notes <- c(notes, .check_model_size(model))

  av <- .car_anova(model, type = type, test_statistic = test_statistic)
  notes <- c(notes, av$note)

  ## Effect sizes -------------------------------------------------------------
  eff <- .glm_effect_sizes(model, family, av$raw)

  ## Dispersion ---------------------------------------------------------------
  dispersion <- .dispersion(model)
  if (family$family %in% c("poisson", "binomial") &&
      is.finite(dispersion) && dispersion > 1.5) {
    notes <- c(notes, sprintf(
      "Pearson dispersion is %.2f. For a %s family this suggests overdispersion; consider quasi%s, or %s.",
      dispersion, family$family, family$family,
      if (family$family == "poisson") "anova_count(model = \"negbin\")" else
        "a beta-binomial model"))
  }

  ## Coefficient intervals ----------------------------------------------------
  rv <- .robust_vcov(model, vcov_type)
  notes <- c(notes, rv$note)
  if (!is.null(rv$matrix) && ci_method == "profile") {
    notes <- c(notes, sprintf(
      "Robust standard errors (%s) were requested, so the coefficient intervals are Wald intervals built from them; profile-likelihood intervals cannot use a sandwich covariance.",
      vcov_type))
  }
  ci <- .coef_intervals(model, conf_level = conf_level,
                        ci_method = if (is.null(rv$matrix)) ci_method else "wald",
                        vcov_matrix = rv$matrix)
  notes <- c(notes, ci$note)

  ## Marginal means and comparisons ------------------------------------------
  emm <- .emm_block(model, groups, additive = .is_additive(interaction),
                    type = "response", vcov_matrix = rv$matrix,
                    conf_level = conf_level, adjust = adjust, posthoc = posthoc,
                    link = family$link, data = d)
  notes <- c(notes, emm$notes)
  if (posthoc && !emm$per_factor && !is.null(emm$table)) {
    if (NROW(emm$table) < 2L) {
      notes <- c(notes, "No pairwise comparisons: fewer than two estimable cells.")
    } else if (NROW(emm$table) == 2L && !is.null(emm$posthoc)) {
      notes <- c(notes, "With two cells the single pairwise comparison is the omnibus test; no multiplicity adjustment was applied.")
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
  }

  .new_fit(
    method       = sprintf("Analysis of deviance (%s family, %s link, Type %s)",
                           family$family, family$link, type),
    call         = cl,
    model        = model,
    anova        = av$table,
    effect_sizes = eff,
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
      model_stats = list(AIC = tryCatch(stats::AIC(model),
                                        error = function(e) NA_real_),
                         BIC = tryCatch(stats::BIC(model),
                                        error = function(e) NA_real_),
                         null_deviance = model$null.deviance,
                         residual_deviance = model$deviance,
                         df_residual = model$df.residual))
  )
}

#' Effect sizes appropriate to the family
#'
#' Partial eta squared for Gaussian models (where sums of squares exist);
#' McFadden's pseudo R squared and the deviance explained otherwise.
#' @noRd
.glm_effect_sizes <- function(model, family, anova_raw) {
  if (family$family == "gaussian" && !is.null(anova_raw) &&
      "Sum Sq" %in% names(as.data.frame(anova_raw))) {
    y <- stats::model.response(stats::model.frame(model))
    n_obs <- length(y)
    return(.partial_eta_squared(anova_raw, df_error = stats::df.residual(model),
                                n_obs = n_obs))
  }
  dev_explained <- if (isTRUE(is.finite(model$null.deviance)) &&
                       model$null.deviance != 0) {
    1 - model$deviance / model$null.deviance
  } else NA_real_
  # The intercept-only model is fitted directly rather than through update().
  # The fitting call is built inside .fit_glm(), so its `family` and `data`
  # are local to that frame and update() cannot re-evaluate them: it fails,
  # and McFadden's R squared came back NA for every non-Gaussian family.
  mcfadden <- tryCatch({
    ll <- as.numeric(stats::logLik(model))
    mf <- stats::model.frame(model)
    y0 <- stats::model.response(mf)
    # A two-column matrix response carries its own totals. Passing the prior
    # weights as well makes binomial()$initialize multiply them in a second
    # time, which inflates the null log-likelihood by a factor of n.
    w0 <- if (is.matrix(y0)) NULL else stats::weights(model, "prior")
    null_fit <- stats::glm(y0 ~ 1, family = family, weights = w0,
                           offset = stats::model.offset(mf))
    1 - ll / as.numeric(stats::logLik(null_fit))
  }, error = function(e) NA_real_)
  data.frame(
    measure  = c("deviance_explained", "mcfadden_r2"),
    estimate = c(dev_explained, mcfadden),
    stringsAsFactors = FALSE
  )
}

#' Coefficient estimates with intervals
#' @noRd
.coef_intervals <- function(model, conf_level = 0.95, ci_method = "profile",
                            vcov_matrix = NULL) {
  ct <- .coef_table(model, vcov_matrix)
  notes <- character(0)
  ci <- NULL
  if (ci_method == "profile") {
    tried <- tryCatch(
      suppressWarnings(suppressMessages(stats::confint(model, level = conf_level))),
      error = function(e) e)
    if (inherits(tried, "error")) {
      notes <- c(notes, sprintf(
        "Profile-likelihood intervals failed, Wald intervals were used instead: %s",
        conditionMessage(tried)))
    } else {
      ci <- matrix(as.numeric(tried), ncol = 2L,
                   dimnames = list(if (is.matrix(tried)) rownames(tried) else ct$term,
                                   NULL))
    }
  }
  method_used <- "profile"
  if (is.null(ci)) {
    crit <- .coef_crit(model, conf_level)
    ci <- cbind(ct$estimate - crit * ct$se, ct$estimate + crit * ct$se)
    rownames(ci) <- ct$term
    method_used <- "wald"
  }
  idx <- match(ct$term, rownames(ci))
  out <- data.frame(ct, conf_low = ci[idx, 1L], conf_high = ci[idx, 2L],
                    ci_method = method_used, stringsAsFactors = FALSE,
                    row.names = NULL)
  list(table = out, note = notes)
}
