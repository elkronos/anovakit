#' Analysis of Deviance for a Count Response
#'
#' Fits a Poisson, quasi-Poisson or negative binomial regression for a count
#' response across one or more grouping variables, checks for overdispersion,
#' and reports an analysis of deviance table, incidence rate ratios, estimated
#' marginal rates and pairwise comparisons.
#'
#' @details
#' \strong{Overdispersion.} The Pearson dispersion statistic is computed from
#' the Poisson fit. With \code{model = "auto"} (the default) a dispersion above
#' \code{overdispersion_threshold} switches the model to negative binomial when
#' \pkg{MASS} is installed, and to quasi-Poisson otherwise (or when the
#' negative binomial fit fails). The threshold is a fixed cut-off, not a test,
#' and a Poisson model kept below it is still affected by a dispersion
#' \eqn{\phi} above 1: its likelihood-ratio and Wald statistics are inflated by
#' about \eqn{\phi} and its intervals are too narrow, so at \eqn{\phi = 1.5} a
#' nominal 5\% test on 1 degree of freedom rejects about 11\% of true null
#' hypotheses, however large the sample. \code{$notes} then gives the
#' dispersion and the size of the effect; use \code{model = "quasipoisson"} or
#' \code{"negbin"} when that matters. When the Poisson model has no residual
#' degrees of freedom the dispersion cannot be computed, and \code{$notes} says
#' that overdispersion was not checked. Which model was used, and why, is
#' always recorded in \code{$notes} and in \code{$model_type}. Set \code{model}
#' explicitly to take the decision out of the function's hands: an explicit
#' \code{model = "poisson"} is kept on overdispersed data (with a note), and an
#' explicit \code{model = "negbin"} that cannot be fitted is an error rather
#' than a substitution.
#'
#' \strong{Negative binomial inference.} The analysis of deviance treats the
#' dispersion parameter theta as known, and its tests and the Wald intervals
#' are asymptotic. With small groups they are anti-conservative, and not only
#' because theta is treated as known: with fewer than about 10 observations
#' per group a nominal 5\% test rejects roughly 10-15\% of true null hypotheses
#' and 95\% intervals cover about 90\%; the excess fades by about 30 per group.
#' The note on theta states the size of the problem for the smallest group in
#' the data, whether theta is well determined, and whether the fit converged.
#' Convergence is read from the fitted object, so the note is the same in
#' every session language.
#'
#' \strong{Exposure.} An \code{offset} column enters the model as
#' \code{offset(log(exposure))} in the formula, not through \code{glm()}'s
#' \code{offset} argument. That matters: an offset supplied through the argument
#' is invisible to \code{terms()}, which means \pkg{emmeans} cannot see it and
#' \code{predict(newdata = )} silently recycles it. Estimated marginal means are
#' reported at an exposure of 1, so they are rates per unit of exposure.
#'
#' \strong{Weights.} Prior weights multiply each row's contribution to the
#' likelihood, which is exactly what a frequency weight does: a row with
#' weight 3 counts as three identical observations, and a table with one row
#' per distinct count and a column of frequencies gives the same estimates and
#' likelihood-ratio tests as the expanded data. When every weight is a whole
#' number and some exceed 1, the weights are taken to be frequencies and the
#' Pearson dispersion (the one the model choice is made on, reported in
#' \code{$dispersion}) divides by the number of observations they represent
#' minus the number of parameters rather than by the number of rows, which
#' would inflate it by up to the ratio of the two. A quasi-Poisson fit and an
#' F test estimate their dispersion from the rows all the same, so for those
#' expand a frequency table to one row per observation; \code{$notes} says so.
#'
#' \strong{Incidence rate ratios.} \code{$effect_sizes} is computed from a
#' treatment-coded parameterisation of the fitted model, whatever \code{type}
#' and the global \code{contrasts} option are, and with ordered factors treated
#' as unordered. Each ratio compares one level of a factor with its first
#' (reference) level, as the \code{contrast} column says; in a model with
#' interactions a main-effect ratio is taken at the reference levels of the
#' factors it interacts with, and an interaction term is a ratio of rate
#' ratios. The coefficients and their covariance (model-based or robust) are
#' mapped exactly onto that coding, which gives what a refit with
#' \code{contr.treatment} would give.
#'
#' \strong{Interactions.} The model is additive by default. Set
#' \code{interaction = TRUE} for the full factorial, or to a whole number for
#' the highest interaction order. With three or more grouping variables the full
#' factorial is often unestimable, so this is deliberately not the default.
#'
#' @param data A data frame, or anything inheriting from one, such as a
#'   \code{data.table} or a tibble.
#' @param response Character. Name of the count column. Must be non-negative
#'   whole numbers, not all zero.
#' @param groups Character vector. One or more grouping columns.
#' @param offset Optional character. Name of a strictly positive exposure
#'   column, entered as a log offset. It cannot also be the response, a
#'   grouping column or the weights.
#' @param interaction \code{FALSE} (additive, the default), \code{TRUE} (full
#'   factorial), or a whole number giving the highest interaction order.
#' @param model Character. \code{"auto"} (default), \code{"poisson"},
#'   \code{"negbin"} or \code{"quasipoisson"}.
#' @param overdispersion_threshold Numeric. Pearson dispersion above which
#'   \code{model = "auto"} moves away from Poisson. Default \code{1.5}.
#' @param type Character. \code{"II"} (default) or \code{"III"} sums of squares.
#'   With \code{"III"} the model is fitted under sum-to-zero contrasts; this
#'   changes \code{$anova} only.
#' @param test_statistic Character. \code{"LR"} (default), \code{"Wald"} or
#'   \code{"F"}, passed to \code{\link[car]{Anova}}. \code{"F"} is the right
#'   choice for a quasi-Poisson model, which has no likelihood and is given it
#'   in place of \code{"LR"}. For a Poisson or negative binomial model,
#'   \code{"F"} makes \code{car::Anova()} estimate a dispersion from the Pearson
#'   residuals, so the table becomes a quasi-likelihood F test while
#'   \code{$effect_sizes} and \code{$posthoc} still assume a dispersion of 1;
#'   \code{$notes} gives the dispersion the F test used. Use
#'   \code{model = "quasipoisson"} for quasi-likelihood inference throughout.
#' @param conf_level Numeric in (0, 1). Level for every interval returned.
#'   Default \code{0.95}.
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
#'   weights, which multiply each row's log-likelihood as frequency weights
#'   would (see Details). Given as a column name rather than a vector so that
#'   it is subsetted with the data: a vector supplied by the caller is
#'   evaluated against the original frame and silently misaligns as soon as
#'   one row is dropped for missing values.
#'
#' @return An \code{\link{anovakit_fit}} object, with \code{$model_type}
#'   naming the model actually fitted, \code{$dispersion} the Pearson dispersion
#'   of the \emph{Poisson} fit (the statistic the model choice was made on) and
#'   \code{$model_dispersion} that of the model actually returned (both with
#'   the frequency-weight denominator when the weights are frequencies).
#'   \code{$effect_sizes} holds incidence rate ratios against the reference
#'   level, with columns \code{term} (the treatment-coded coefficient),
#'   \code{contrast} (what it compares), \code{IRR}, \code{conf_low},
#'   \code{conf_high} and \code{p_value}. \code{$emmeans} holds the estimated
#'   marginal rates.
#'
#' @seealso \code{\link{anova_bin}} for binary responses.
#'
#' @examples
#' set.seed(1)
#' n <- 400
#' d <- data.frame(
#'   g1 = factor(rep(c("A", "B"), each = n / 2)),
#'   g2 = factor(rep(rep(c("X", "Y"), each = n / 4), 2))
#' )
#' lambda <- with(d, ifelse(g1 == "A" & g2 == "X", 5,
#'                   ifelse(g1 == "A" & g2 == "Y", 10,
#'                   ifelse(g1 == "B" & g2 == "X", 15, 20))))
#' d$count <- rpois(n, lambda)
#' fit <- anova_count(d, "count", c("g1", "g2"), interaction = TRUE)
#' fit
#' fit$emmeans
#' fit$posthoc
#'
#' # With an exposure offset, the marginal means are rates per unit exposure
#' d$hours <- runif(n, 0.5, 4)
#' anova_count(d, "count", "g1", offset = "hours", plots = FALSE)$emmeans
#'
#' @export
anova_count <- function(data, response, groups,
                        offset = NULL,
                        interaction = FALSE,
                        model = c("auto", "poisson", "negbin", "quasipoisson"),
                        overdispersion_threshold = 1.5,
                        type = c("II", "III"),
                        test_statistic = c("LR", "Wald", "F"),
                        conf_level = 0.95,
                        vcov_type = "model",
                        adjust = "tukey",
                        posthoc = TRUE,
                        plots = TRUE,
                        weights = NULL,
                        verbose = FALSE) {
  cl <- match.call()
  model <- match.arg(model)
  type <- match.arg(type)
  test_statistic <- match.arg(test_statistic)
  .check_adjust(vcov_type, .vcov_choices(), "vcov_type")
  .check_adjust(adjust, .adjust_choices(), "adjust")
  data <- .as_df(data)
  .check_name(response, "response")
  .check_names(groups, "groups")
  .check_columns(data, response, "Response column")
  .check_columns(data, groups, "Grouping column(s)")
  if (!is.null(offset)) {
    .check_name(offset, "offset")
    .check_columns(data, offset, "Offset column")
  }
  .check_conf_level(conf_level)
  .check_interaction(interaction, length(groups))
  .check_roles(list(`the response` = response, `a grouping variable` = groups,
                    `the exposure offset` = offset, `the weights` = weights))
  .check_flag(posthoc, "posthoc")
  .check_flag(plots, "plots")
  .check_flag(verbose, "verbose")
  if (!is.numeric(overdispersion_threshold) ||
      length(overdispersion_threshold) != 1L ||
      !is.finite(overdispersion_threshold) || overdispersion_threshold <= 0) {
    .stopf("`overdispersion_threshold` must be a single positive number.")
  }

  .check_counts(data[[response]], response)

  .say(verbose, "Preparing data.")
  .check_weights_column(data, weights, response, groups, offset)
  cols <- c(response, groups, offset, weights)
  prep <- .prepare_frame(data, cols, factors = groups, keep_all = FALSE)
  d <- prep$data
  notes <- prep$notes
  zw <- .drop_zero_weights(d, weights)
  d <- zw$data
  n_removed <- prep$n_removed + zw$n_removed
  notes <- c(notes, zw$notes)
  .check_groups(d, groups, min_levels = 2L, min_n = 1L)
  # The rows that are analysed can be all zero even when the input was not.
  .check_counts(d[[response]], response)

  if (!is.null(offset)) {
    if (!is.numeric(d[[offset]]) || any(d[[offset]] <= 0)) {
      .stopf("Offset column `%s` must be numeric and strictly positive; a log offset is undefined otherwise.",
             offset)
    }
  }

  freq_w <- .frequency_weights(d, weights)
  notes <- c(notes, .sparse_cell_note(d, groups, response, interaction,
                                      if (freq_w) weights else NULL))

  ## Build the formula, with the offset inside it -----------------------------
  terms_rhs <- .group_terms(groups, interaction)
  if (!is.null(offset)) {
    terms_rhs <- c(terms_rhs, sprintf("offset(log(%s))", .bq(offset)))
  }
  fml <- stats::reformulate(terms_rhs, response = .bq(response))
  .say(verbose, "Fitting %s", paste(deparse(fml), collapse = " "))

  ## Poisson first, to measure dispersion -------------------------------------
  # Wrapped like every other fitting call: glm() warns about numerically zero
  # fitted rates, and that belongs in $notes rather than on the console.
  got <- .collect_conditions(.fit_with_contrasts(
    function() .fit_glm(fml, d, stats::poisson(link = "log"), weights),
    type = type))
  pois <- got$value
  if (inherits(pois, "error")) {
    .stopf("The Poisson model could not be fitted: %s", conditionMessage(pois))
  }
  notes <- c(notes, .said_note(got$said, "stats::glm()"))
  # With frequency weights the residual degrees of freedom are counted in
  # observations, not rows (see .count_dispersion()).
  n_obs <- if (freq_w) sum(d[[weights]]) else NULL
  dispersion <- .count_dispersion(pois, n_obs)

  chosen <- .choose_count_model(model, dispersion, overdispersion_threshold)
  notes <- c(notes, chosen$notes)
  if (freq_w) {
    notes <- c(notes, sprintf(
      "The weights in `%s` are whole numbers, so they are treated as frequency weights: each row stands for that many observations. The Pearson dispersion divides by the %s observations they represent minus the %d parameters, not by the %d rows, which would inflate it.",
      weights, format(n_obs, big.mark = ","), pois$rank, nrow(d)))
  }

  fit <- pois
  model_type <- "poisson"
  if (chosen$model == "negbin") {
    # glm.nb() captures `link` with substitute(), so it must not be passed
    # through do.call() as a value. The default is already log.
    nbg <- .collect_conditions(.fit_with_contrasts(
      function() .fit_glm(fml, d, NULL, weights),
      type = type))
    nb <- nbg$value
    if (inherits(nb, "error")) {
      if (model == "negbin") {
        # Asked for by name: substituting another model would contradict the
        # request, so say what failed instead.
        .stopf("The negative binomial model could not be fitted: %s%s. Use model = \"poisson\" or model = \"quasipoisson\" instead.",
               conditionMessage(nb),
               if (!isTRUE(stats::df.residual(pois) >= 1))
                 " (theta cannot be estimated when the model has no residual degrees of freedom)"
               else "")
      }
      notes <- c(notes, sprintf(
        "The negative binomial fit failed (%s); a quasi-Poisson model was used instead.",
        conditionMessage(nb)), .said_note(nbg$said, "MASS::glm.nb()"))
      qp <- .collect_conditions(.fit_with_contrasts(
        function() .fit_glm(fml, d, stats::quasipoisson(), weights),
        type = type))
      fit <- qp$value
      if (inherits(fit, "error")) {
        .stopf("Neither a negative binomial nor a quasi-Poisson model could be fitted: %s",
               conditionMessage(fit))
      }
      notes <- c(notes, .said_note(qp$said, "stats::glm()"))
      model_type <- "quasipoisson"
    } else {
      fit <- nb
      model_type <- "negbin"
      attr(fit, "nb_warnings") <- nbg$said
      notes <- c(notes, .theta_note(
        fit, .smallest_group(d, groups, interaction,
                             if (freq_w) weights else NULL)))
      # A theta that hit its limit is reported by .theta_note(); MASS warns
      # with the same (translated) text, which would only repeat it.
      th_warn <- trimws(gsub("[[:space:]]+", " ", as.character(fit$th.warn)))
      notes <- c(notes, .said_note(setdiff(nbg$said, th_warn), "MASS::glm.nb()"))
    }
  } else if (chosen$model == "quasipoisson") {
    qp <- .collect_conditions(.fit_with_contrasts(
      function() .fit_glm(fml, d, stats::quasipoisson(), weights),
      type = type))
    fit <- qp$value
    if (inherits(fit, "error")) {
      .stopf("The quasi-Poisson model could not be fitted: %s",
             conditionMessage(fit))
    }
    notes <- c(notes, .said_note(qp$said, "stats::glm()"))
    model_type <- "quasipoisson"
  }

  if (model_type == "quasipoisson" && test_statistic == "LR") {
    test_statistic <- "F"
    notes <- c(notes, "A quasi-Poisson model has no likelihood, so the analysis of deviance uses an F test rather than a likelihood ratio test.")
  }
  notes <- c(notes, .check_model_size(fit))
  notes <- c(notes, .check_separation(fit))

  ## Analysis of deviance -----------------------------------------------------
  av <- .car_anova(fit, type = type, test_statistic = test_statistic)
  notes <- c(notes, av$note)
  if (test_statistic == "F" && model_type != "quasipoisson") {
    notes <- c(notes, .f_test_note(model_type, .dispersion(fit)))
  }
  if (freq_w && (model_type == "quasipoisson" || test_statistic == "F")) {
    quasi <- model_type == "quasipoisson"
    notes <- c(notes, sprintf(
      "%s estimates its dispersion from the %d rows, not from the %s observations the frequency weights represent, so %s can be far too conservative here. Expand the table to one row per observation for %s.",
      if (quasi) "A quasi-Poisson model" else "The F test",
      nrow(d), format(n_obs, big.mark = ","),
      if (quasi) "its tests and intervals" else "it",
      if (quasi) "quasi-Poisson inference" else "an F test"))
  }

  ## Marginal rates and pairwise incidence rate ratios ------------------------
  rv <- .robust_vcov(fit, vcov_type)
  notes <- c(notes, rv$note)
  emm_extra <- if (!is.null(offset)) list(offset = 0) else list()
  emm <- do.call(.emm_block, c(list(
    model = fit, groups = groups, additive = .is_additive(interaction),
    type = "response", vcov_matrix = rv$matrix, conf_level = conf_level,
    adjust = adjust, posthoc = posthoc, link = "log", data = d), emm_extra))
  notes <- c(notes, emm$notes)
  if (!is.null(offset) && !is.null(emm$table)) {
    notes <- c(notes, sprintf(
      "Estimated marginal means are rates at %s = 1 (one unit of exposure).", offset))
  }
  if (!posthoc) {
    notes <- c(notes, "$effect_sizes still reports the incidence rate ratios from the model.")
  }
  if (!is.null(emm$posthoc) && "ratio" %in% names(emm$posthoc)) {
    names(emm$posthoc)[names(emm$posthoc) == "ratio"] <- "IRR"
  }

  ## Incidence rate ratios against the reference level -----------------------
  irr <- .irr_table(fit, groups, rv$matrix, conf_level)
  notes <- c(notes, irr$note)

  ## Plots --------------------------------------------------------------------
  plot_list <- list()
  if (plots && !is.null(emm$table)) {
    .say(verbose, "Building plots.")
    plot_list$emmeans <- .plot_emmeans(
      emm$table, groups, estimate = "estimate",
      lower = "conf_low", upper = "conf_high", conf_level = conf_level,
      title = sprintf("Estimated %s rate by group", response),
      ylab = sprintf("Estimated %s", response))
    added <- .add_cell(d, groups)
    plot_list$observed <- .plot_box(added$data, response, added$cell,
                                    title = "Observed counts by group",
                                    xlab = paste(groups, collapse = " : "))
  } else if (plots) {
    notes <- c(notes, "Plots were skipped because the estimated marginal means could not be computed.")
  }

  model_dispersion <- .count_dispersion(fit, n_obs)
  .new_fit(
    method       = sprintf("Analysis of deviance for a count response (%s regression)",
                           model_type),
    call         = cl,
    model        = fit,
    anova        = av$table,
    effect_sizes = irr$table,
    emmeans      = emm$table,
    emmeans_object = emm$grid,
    posthoc      = emm$posthoc,
    assumptions  = list(poisson_dispersion = dispersion,
                        model_dispersion = model_dispersion,
                        cell_counts = .cell_counts(d, groups)),
    plots        = plot_list,
    data_used    = d,
    n_removed    = n_removed,
    conf_level   = conf_level,
    notes        = notes,
    extra        = list(model_type = model_type,
                        dispersion = dispersion,
                        model_dispersion = model_dispersion)
  )
}

#' Validate a count response
#' @noRd
.check_counts <- function(x, name) {
  if (!is.numeric(x)) {
    .stopf("Response `%s` must be numeric (non-negative whole numbers); it is %s.",
           name, paste(class(x), collapse = "/"))
  }
  # Non-finite values are dropped by .prepare_frame() and counted in $n_removed,
  # so they must not be judged here: abs(Inf - round(Inf)) is NaN, and the
  # whole-number test below would fail with "missing value where TRUE/FALSE
  # needed" before the row ever reached the code that removes it.
  v <- x[is.finite(x)]
  if (length(v) == 0L) {
    .stopf("Response `%s` has no finite values; every row is missing or infinite.",
           name)
  }
  if (any(v < 0)) {
    .stopf("Response `%s` contains negative values; counts cannot be negative.", name)
  }
  if (any(abs(v - round(v)) > .Machine$double.eps^0.5)) {
    .stopf("Response `%s` contains non-integer values; counts must be whole numbers.",
           name)
  }
  # With no events at all every log rate is minus infinity: glm() stops at an
  # arbitrary point, the likelihood-ratio statistic is numerically zero, and
  # nothing reported would mean anything.
  if (all(v == 0)) {
    .stopf("Response `%s` is zero in every row: with no events there is no rate to estimate or compare.",
           name)
  }
  invisible(TRUE)
}

#' Are the prior weights frequencies?
#'
#' Whole-number weights of which some exceed 1 are read as a frequency table:
#' each row stands for that many identical observations. For a Poisson or
#' negative binomial likelihood that reading changes nothing in the estimates,
#' only the number of observations behind the residual degrees of freedom.
#' @noRd
.frequency_weights <- function(d, weights) {
  if (is.null(weights)) return(FALSE)
  w <- d[[weights]]
  all(abs(w - round(w)) < sqrt(.Machine$double.eps)) && any(w > 1)
}

#' Pearson dispersion, with the frequency-weight denominator when needed
#'
#' With frequency weights the Pearson sum of squares already equals that of
#' the expanded data, but \code{df.residual()} counts rows. Dividing by rows
#' inflates the statistic by (observations - p) / (rows - p), which can switch
#' Poisson data to a negative binomial model.
#' @param fit a fitted count GLM.
#' @param n_obs the number of observations the weights represent, or
#'   \code{NULL} when the weights are not frequencies.
#' @noRd
.count_dispersion <- function(fit, n_obs = NULL) {
  if (is.null(n_obs)) return(.dispersion(fit))
  dfr <- n_obs - fit$rank
  if (!is.finite(dfr) || dfr <= 0) return(NA_real_)
  sum(stats::residuals(fit, type = "pearson")^2, na.rm = TRUE) / dfr
}

#' What a Poisson fit at dispersion phi does to its tests and intervals
#'
#' Under a true dispersion phi a Poisson likelihood-ratio or Wald statistic is
#' about phi times its nominal chi-square, so the size of a nominal 5\% test
#' on 1 degree of freedom is P(chisq_1 > qchisq(0.95, 1) / phi).
#' @noRd
.inflation_note <- function(phi) {
  size <- stats::pchisq(stats::qchisq(0.95, 1) / phi, 1, lower.tail = FALSE)
  sprintf("A Poisson model assumes a dispersion of 1; at %.2f its likelihood-ratio and Wald statistics are inflated by about that factor and its standard errors are about %.2f times too small, so its intervals are too narrow and a nominal 5%% test on 1 degree of freedom rejects roughly %.0f%% of true null hypotheses. Use model = \"quasipoisson\" or \"negbin\" when that matters.",
          phi, sqrt(phi), 100 * size)
}

#' Decide which count model to fit
#' @noRd
.choose_count_model <- function(model, dispersion, threshold) {
  if (model == "negbin" && !requireNamespace("MASS", quietly = TRUE)) {
    .stopf("model = \"negbin\" requires the {MASS} package. Install it, or use model = \"quasipoisson\".")
  }
  if (model %in% c("negbin", "quasipoisson")) {
    return(list(model = model, notes = character(0)))
  }
  finite <- is.finite(dispersion)
  if (model == "poisson") {
    # Requested by name, so kept; but the same diagnosis anova_glm() gives a
    # Poisson family belongs here too.
    notes <- if (finite && dispersion > threshold) {
      sprintf("Pearson dispersion is %.2f, above the threshold of %.2f: the counts are overdispersed relative to Poisson. model = \"poisson\" was requested, so the Poisson model was kept. %s",
              dispersion, threshold, .inflation_note(dispersion))
    } else if (finite && dispersion > 1) {
      sprintf("Pearson dispersion is %.2f. %s", dispersion, .inflation_note(dispersion))
    } else character(0)
    return(list(model = "poisson", notes = notes))
  }
  if (!finite) {
    return(list(model = "poisson", notes =
      "The Pearson dispersion cannot be computed because the Poisson model has no residual degrees of freedom, so overdispersion could not be checked and a Poisson model was kept."))
  }
  if (dispersion <= threshold) {
    note <- sprintf(
      "Pearson dispersion is %.2f, at or below the threshold of %.2f, so a Poisson model was kept.",
      dispersion, threshold)
    if (dispersion > 1) note <- paste(note, .inflation_note(dispersion))
    return(list(model = "poisson", notes = note))
  }
  if (requireNamespace("MASS", quietly = TRUE)) {
    return(list(model = "negbin", notes = sprintf(
      "Pearson dispersion is %.2f, above the threshold of %.2f, so a negative binomial model was fitted instead of Poisson. Set model = \"poisson\" to override.",
      dispersion, threshold)))
  }
  list(model = "quasipoisson", notes = sprintf(
    "Pearson dispersion is %.2f, above the threshold of %.2f. The {MASS} package is not installed, so a quasi-Poisson model was fitted instead of a negative binomial one.",
    dispersion, threshold))
}

#' The dispersion behind an F test on a fixed-dispersion count model
#'
#' car::Anova() computes its F test for any glm with the Pearson dispersion
#' estimated from the residuals, which turns a Poisson or negative binomial
#' analysis of deviance into a quasi-likelihood test, while the coefficient
#' table and the emmeans comparisons keep the dispersion at 1.
#' @noRd
.f_test_note <- function(model_type, dispersion) {
  what <- if (model_type == "negbin") "negative binomial" else "Poisson"
  disp <- if (is.finite(dispersion)) sprintf("%.2f", dispersion) else
    "which cannot be estimated here, as there are no residual degrees of freedom"
  sprintf("test_statistic = \"F\" was requested for a %s model, so car::Anova() estimated a dispersion from the Pearson residuals (%s) and divided the deviance by it: the analysis of deviance is a quasi-likelihood F test. $effect_sizes and $posthoc still assume a dispersion of 1, so they can disagree with it. Use model = \"quasipoisson\" for quasi-likelihood inference throughout.",
          what, disp)
}

#' Note sparse or empty cells
#'
#' A count cell's information is in its events, not its rows: one row of
#' aggregated data holding 80 events pins its rate down well, while a hundred
#' rows holding two events do not. Sparse cells are therefore judged on the
#' events they hold (weighted by frequency weights when there are some).
#'
#' Whether an empty combination of levels is estimable depends on the model:
#' the full factorial cannot estimate it, an additive model extrapolates it
#' from the main effects.
#' @noRd
.sparse_cell_note <- function(d, groups, response, interaction = FALSE,
                              weights = NULL) {
  levs <- lapply(groups, function(g) levels(droplevels(as.factor(d[[g]]))))
  n_possible <- prod(lengths(levs))
  cell <- .cell_vector(d, groups)
  n_observed <- nlevels(cell)
  notes <- character(0)
  empty <- n_possible - n_observed
  if (empty > 0L) {
    full <- isTRUE(interaction) ||
      (is.numeric(interaction) && interaction >= length(groups))
    notes <- c(notes, sprintf(
      "%d of the %d possible group combination(s) contain no observations; %s",
      empty, n_possible,
      if (.is_additive(interaction)) {
        "the additive model still estimates them from the main effects, so marginal means that average over them rely on there being no interaction."
      } else if (full) {
        "contrasts involving them are not estimable."
      } else {
        "they are estimable only where the lower-order terms determine them (see the notes on the marginal means)."
      }))
  }
  w <- if (is.null(weights)) 1 else d[[weights]]
  events <- tapply(w * d[[response]], cell, sum)
  small <- sum(events < 5)
  if (small > 0L) {
    notes <- c(notes, sprintf(
      "%d group combination(s) hold fewer than 5 events in total; their rates, and any rate ratio involving them, rest on very little information, so their Wald intervals and p-values are unreliable.",
      small))
  }
  notes
}

#' The number of observations behind the smallest group
#'
#' For an additive model each factor's effects are estimated from its levels,
#' so the smallest level of any factor is what matters; with interactions each
#' cell has parameters of its own, so the smallest populated cell is.
#' @noRd
.smallest_group <- function(d, groups, interaction = FALSE, weights = NULL) {
  w <- if (is.null(weights)) rep(1, nrow(d)) else d[[weights]]
  if (.is_additive(interaction)) {
    return(min(vapply(groups, function(g) {
      min(tapply(w, droplevels(as.factor(d[[g]])), sum))
    }, numeric(1))))
  }
  min(tapply(w, .cell_vector(d, groups), sum))
}

#' Describe the negative binomial dispersion parameter honestly
#'
#' Non-convergence is read from the fit (\code{th.warn}, which
#' \code{MASS::glm.nb()} sets whenever theta hits its iteration limit, the
#' alternation limit, or zero; and \code{converged}), never from warning text,
#' which MASS translates into the session language.
#'
#' A theta in the hundreds or thousands is not an estimate of overdispersion, it
#' is the fit telling you there is none: the negative binomial has collapsed
#' towards the Poisson it generalises. A small theta is strong overdispersion
#' (variance mu + mu^2 / theta) and is imprecise in small samples for that very
#' reason, so an imprecise theta is not evidence of weak overdispersion.
#'
#' The calibration statements come from simulation (3 groups, theta 1 and 2,
#' 3000 replicates per cell): nominal 5\% tests rejected 14-15\% of true nulls
#' at 5 per group, 10-12\% at 8, 9-10\% at 10, 7-8\% at 15-20, 6-7\% at 30 and
#' 5.6-6.1\% at 50; 95\% intervals for rate ratios covered 91-92\% at 8, 92.5-93\%
#' at 10 and 94\% at 15, and Tukey intervals jointly 89-90\% at 8.
#' @param nb a \code{negbin} fit.
#' @param n_min observations in the smallest group (see .smallest_group()).
#' @noRd
.theta_note <- function(nb, n_min = NA_real_) {
  th <- nb$theta; se <- nb$SE.theta
  if (!is.null(nb$th.warn) || !isTRUE(nb$converged)) {
    why <- if (!is.null(nb$th.warn)) {
      sprintf("MASS::glm.nb() stopped with \"%s\"", nb$th.warn[1L])
    } else {
      "the iteratively reweighted least squares fit did not converge"
    }
    return(sprintf(
      "The negative binomial fit did not converge (%s; theta reached %.3g). Treat the model as unreliable and compare it with model = \"poisson\" or model = \"quasipoisson\".",
      why, th))
  }
  if (!is.finite(th) || th > 1000) {
    return(sprintf(
      "The negative binomial dispersion parameter is very large (theta = %.3g), which means the data are not meaningfully overdispersed relative to Poisson and the negative binomial has collapsed towards it. model = \"poisson\" would give essentially the same answer more simply.",
      th))
  }
  out <- sprintf(
    "Negative binomial dispersion parameter theta = %.3f (SE %.3f), so the variance is about mu + mu^2/%.3g.",
    th, se, th)
  if (is.finite(se) && se > 0 && th / se < 2) {
    out <- paste(out, if (th > 10) {
      "Theta is poorly determined, but it is large: the data are only weakly overdispersed, and model = \"poisson\" would give a similar answer."
    } else {
      "Theta is poorly determined by data this sparse, although a value this small means substantial overdispersion, not weak."
    })
  }
  size <- if (is.finite(n_min)) sprintf("%s observation%s",
                                        format(n_min, big.mark = ","),
                                        if (n_min == 1) "" else "s") else NA
  calib <- if (is.finite(n_min) && n_min < 10) {
    sprintf("With %s in the smallest group, the likelihood-ratio tests and Wald intervals are markedly anti-conservative: at 5-10 observations per group, nominal 5%% tests reject roughly 10-15%% of true null hypotheses (2-3 times too often) and 95%% intervals, including the pairwise comparisons, cover only about 89-93%%. Treating theta as known explains only part of this; the rest is small-sample asymptotics. model = \"quasipoisson\" (F tests) is better calibrated at this size.",
            size)
  } else if (is.finite(n_min) && n_min < 30) {
    sprintf("With %s in the smallest group, the likelihood-ratio tests and Wald intervals are moderately anti-conservative: at 10-30 observations per group, nominal 5%% tests reject roughly 6-10%% of true null hypotheses and 95%% intervals cover about 91-95%%.",
            size)
  } else {
    "The analysis of deviance treats theta as known, so its p-values are mildly anti-conservative (roughly 5-7% rejection of true null hypotheses at a nominal 5% level with 30 or more observations per group)."
  }
  paste(out, calib)
}

#' Incidence rate ratios against the reference level
#'
#' \code{exp(coef(fit))} is a rate ratio against a reference level only under
#' treatment contrasts. Under the sum-to-zero coding that \code{type = "III"}
#' needs, or a global \code{contrasts} option, or polynomial contrasts for an
#' ordered factor, it is something else. The fitted coefficients and their
#' covariance are therefore mapped onto a treatment-coded design of the same
#' model: both designs span the same column space, so \eqn{X_o = X_t A} for a
#' square matrix \eqn{A}, and the treatment-coded coefficients are
#' \eqn{A \beta} with covariance \eqn{A V A'}. This is exact and is what a
#' refit with \code{contr.treatment} gives (for a quasi-likelihood or negative
#' binomial fit too, whose estimates are equivariant under reparameterisation),
#' and it works with robust covariances and with aliased coefficients.
#'
#' @param fit the fitted count model.
#' @param groups the grouping factors.
#' @param vcov_matrix a robust covariance, or \code{NULL} for model-based.
#' @return list with \code{table} and \code{note}.
#' @noRd
.irr_table <- function(fit, groups, vcov_matrix = NULL, conf_level = 0.95) {
  got <- tryCatch(.treatment_coefficients(fit, groups, vcov_matrix),
                  error = function(e) e)
  if (inherits(got, "error")) {
    return(list(table = NULL, note = sprintf(
      "Incidence rate ratios could not be computed: %s", conditionMessage(got))))
  }
  crit <- .coef_crit(fit, conf_level)
  use_t <- .estimates_dispersion(fit)
  dfr <- stats::df.residual(fit)
  stat <- got$estimate / got$se
  p <- if (use_t) {
    if (is.null(dfr) || !is.finite(dfr) || dfr < 1) rep(NA_real_, length(stat))
    else 2 * stats::pt(-abs(stat), df = dfr)
  } else {
    2 * stats::pnorm(-abs(stat))
  }
  tab <- data.frame(
    term = got$term, contrast = got$contrast,
    IRR = exp(got$estimate),
    conf_low = exp(got$estimate - crit * got$se),
    conf_high = exp(got$estimate + crit * got$se),
    p_value = p, stringsAsFactors = FALSE)
  tab <- tab[tab$term != "(Intercept)", , drop = FALSE]
  row.names(tab) <- NULL
  list(table = tab, note = character(0))
}

#' Map a fitted model's coefficients onto treatment coding
#' @noRd
.treatment_coefficients <- function(fit, groups, vcov_matrix = NULL) {
  co <- stats::coef(fit)
  V <- if (is.null(vcov_matrix)) stats::vcov(fit) else vcov_matrix
  keep <- names(co)[!is.na(co)]
  keep <- keep[keep %in% colnames(V)]
  Xo <- stats::model.matrix(fit)[, keep, drop = FALSE]
  bo <- co[keep]
  Vo <- V[keep, keep, drop = FALSE]

  mf <- stats::model.frame(fit)
  tt <- stats::terms(fit)
  facs <- groups[groups %in% names(mf)]
  mf[facs] <- lapply(mf[facs], function(x) factor(x, levels = levels(x), ordered = FALSE))
  ca <- stats::setNames(rep(list("contr.treatment"), length(facs)), facs)
  Xt <- stats::model.matrix(tt, mf, contrasts.arg = ca)

  # The treatment-coded columns that are estimable, in the order glm() would
  # keep them when some are aliased (empty cells of a factorial).
  qx <- qr(Xt)
  est <- sort(qx$pivot[seq_len(qx$rank)])
  Xe <- Xt[, est, drop = FALSE]
  A <- qr.coef(qr(Xe), Xo)
  if (length(est) != length(keep) || anyNA(A) ||
      max(abs(Xe %*% A - Xo)) > 1e-7 * max(1, abs(Xo))) {
    stop("the treatment-coded design does not reproduce the fitted one", call. = FALSE)
  }
  b <- drop(A %*% bo)
  Vt <- A %*% Vo %*% t(A)

  estimate <- se <- rep(NA_real_, ncol(Xt))
  estimate[est] <- b
  se[est] <- sqrt(pmax(diag(Vt), 0))
  list(term = gsub("`", "", colnames(Xt), fixed = TRUE),
       contrast = .treatment_labels(Xt, tt, mf, facs),
       estimate = estimate, se = se)
}

#' Say what each treatment-coded coefficient compares
#'
#' Under treatment coding a model matrix column for a term of factors
#' f1, ..., fm corresponds to one non-reference level of each, with f1
#' varying fastest (the order of \code{expand.grid()}). A main effect is a
#' ratio against the reference level, taken at the reference levels of any
#' factor it interacts with; a two-way term is a ratio of two such ratios.
#' @noRd
.treatment_labels <- function(Xt, tt, mf, facs) {
  labels <- gsub("`", "", colnames(Xt), fixed = TRUE)
  asg <- attr(Xt, "assign")
  fmat <- attr(tt, "factors")
  if (is.null(fmat) || length(fmat) == 0L) return(labels)
  vars <- gsub("`", "", rownames(fmat), fixed = TRUE)
  term_vars <- lapply(seq_len(ncol(fmat)), function(j) vars[fmat[, j] > 0])
  ref <- vapply(facs, function(f) levels(mf[[f]])[1L], character(1))
  for (j in seq_along(term_vars)) {
    vs <- term_vars[[j]]
    cols <- which(asg == j)
    if (!all(vs %in% facs) || length(cols) == 0L) next
    combos <- expand.grid(lapply(vs, function(f) levels(mf[[f]])[-1L]),
                          stringsAsFactors = FALSE)
    if (nrow(combos) != length(cols)) next
    # Factors this term is conditioned on: those it interacts with in a
    # higher-order term, held at their reference level.
    higher <- vapply(term_vars, function(tv) all(vs %in% tv) &&
                       length(tv) > length(vs), logical(1))
    cond <- setdiff(unique(unlist(term_vars[higher])), vs)
    at <- if (length(cond) > 0L) {
      paste(sprintf("%s = %s", cond, ref[cond]), collapse = ", ")
    } else ""
    name <- paste(vs, collapse = ":")
    for (r in seq_len(nrow(combos))) {
      lv <- unlist(combos[r, ], use.names = FALSE)
      labels[cols[r]] <- if (length(vs) == 1L) {
        sprintf("%s: %s / %s%s", name, lv, ref[vs],
                if (nzchar(at)) sprintf(" at %s", at) else "")
      } else if (length(vs) == 2L) {
        sprintf("%s: (%s / %s at %s = %s) / (%s / %s at %s = %s)%s",
                name, lv[1L], ref[vs[1L]], vs[2L], lv[2L],
                lv[1L], ref[vs[1L]], vs[2L], ref[vs[2L]],
                if (nzchar(at)) sprintf(", at %s", at) else "")
      } else {
        sprintf("%s: interaction ratio of rate ratios for %s%s", name,
                paste(sprintf("%s = %s vs %s", vs, lv, ref[vs]), collapse = ", "),
                if (nzchar(at)) sprintf(", at %s", at) else "")
      }
    }
  }
  labels
}
