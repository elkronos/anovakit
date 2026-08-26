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
#' \pkg{MASS} is installed, and to quasi-Poisson otherwise. Which model was
#' used, and why, is always recorded in \code{$notes} and in
#' \code{$model_type}. Set \code{model} explicitly to take the decision out of
#' the function's hands.
#'
#' \strong{Exposure.} An \code{offset} column enters the model as
#' \code{offset(log(exposure))} in the formula, not through \code{glm()}'s
#' \code{offset} argument. That matters: an offset supplied through the argument
#' is invisible to \code{terms()}, which means \pkg{emmeans} cannot see it and
#' \code{predict(newdata = )} silently recycles it. Estimated marginal means are
#' reported at an exposure of 1, so they are rates per unit of exposure.
#'
#' \strong{Interactions.} The model is additive by default. Set
#' \code{interaction = TRUE} for the full factorial, or to a whole number for
#' the highest interaction order. With three or more grouping variables the full
#' factorial is often unestimable, so this is deliberately not the default.
#'
#' @param data A data frame, or anything inheriting from one, such as a
#'   \code{data.table} or a tibble.
#' @param response Character. Name of the count column. Must be non-negative
#'   whole numbers.
#' @param groups Character vector. One or more grouping columns.
#' @param offset Optional character. Name of a strictly positive exposure
#'   column, entered as a log offset.
#' @param interaction \code{FALSE} (additive, the default), \code{TRUE} (full
#'   factorial), or a whole number giving the highest interaction order.
#' @param model Character. \code{"auto"} (default), \code{"poisson"},
#'   \code{"negbin"} or \code{"quasipoisson"}.
#' @param overdispersion_threshold Numeric. Pearson dispersion above which
#'   \code{model = "auto"} moves away from Poisson. Default \code{1.5}.
#' @param type Character. \code{"II"} (default) or \code{"III"} sums of squares.
#'   With \code{"III"} the model is fitted under sum-to-zero contrasts.
#' @param test_statistic Character. \code{"LR"} (default), \code{"Wald"} or
#'   \code{"F"}, passed to \code{\link[car]{Anova}}. \code{"F"} is the right
#'   choice for a quasi-Poisson model.
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
#'   weights. Given as a column name rather than a vector so that it is
#'   subsetted with the data: a vector supplied by the caller is evaluated
#'   against the original frame and silently misaligns as soon as one row is
#'   dropped for missing values.
#'
#' @return An \code{\link{anovakit_fit}} object, with \code{$model_type}
#'   naming the model actually fitted, \code{$dispersion} the Pearson dispersion
#'   of the \emph{Poisson} fit (the statistic the model choice was made on) and
#'   \code{$model_dispersion} that of the model actually returned.
#'   \code{$effect_sizes} holds incidence rate ratios, \code{$emmeans} the
#'   estimated marginal rates.
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
  .check_groups(d, groups, min_levels = 2L, min_n = 1L)

  if (!is.null(offset)) {
    if (!is.numeric(d[[offset]]) || any(d[[offset]] <= 0)) {
      .stopf("Offset column `%s` must be numeric and strictly positive; a log offset is undefined otherwise.",
             offset)
    }
  }

  wts <- .resolve_weights(d, weights)
  notes <- c(notes, .sparse_cell_note(d, groups))

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
    function() .fit_glm(fml, d, stats::poisson(link = "log"), wts),
    type = type))
  pois <- got$value
  if (inherits(pois, "error")) {
    .stopf("The Poisson model could not be fitted: %s", conditionMessage(pois))
  }
  notes <- c(notes, .said_note(got$said, "stats::glm()"))
  dispersion <- .dispersion(pois)

  chosen <- .choose_count_model(model, dispersion, overdispersion_threshold)
  notes <- c(notes, chosen$notes)

  fit <- pois
  model_type <- "poisson"
  if (chosen$model == "negbin") {
    # glm.nb() captures `link` with substitute(), so it must not be passed
    # through do.call() as a value. The default is already log.
    nb_warnings <- character(0)
    nb <- withCallingHandlers(
      tryCatch(
        .fit_with_contrasts(
          function() .fit_glm(fml, d, NULL, wts),
          type = type),
        error = function(e) e),
      warning = function(w) {
        nb_warnings <<- c(nb_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    if (inherits(nb, "error")) {
      notes <- c(notes, sprintf(
        "The negative binomial fit failed (%s); a quasi-Poisson model was used instead.",
        conditionMessage(nb)))
      qp <- .collect_conditions(.fit_with_contrasts(
        function() .fit_glm(fml, d, stats::quasipoisson(), wts),
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
      attr(fit, "nb_warnings") <- nb_warnings
      # Read from `fit`, which carries the attribute: `nb` is the pre-copy
      # object and never has it, so the non-convergence branch was dead.
      notes <- c(notes, .theta_note(fit))
    }
  } else if (chosen$model == "quasipoisson") {
    qp <- .collect_conditions(.fit_with_contrasts(
      function() .fit_glm(fml, d, stats::quasipoisson(), wts),
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

  ## Analysis of deviance -----------------------------------------------------
  av <- .car_anova(fit, type = type, test_statistic = test_statistic)
  notes <- c(notes, av$note)

  ## Marginal rates and pairwise incidence rate ratios ------------------------
  rv <- .robust_vcov(fit, vcov_type)
  notes <- c(notes, rv$note)
  emm_args <- list(model = fit, specs = groups, type = "response",
                   vcov_matrix = rv$matrix)
  if (!is.null(offset)) emm_args$offset <- 0
  emm <- do.call(.emmeans_grid, emm_args)
  notes <- c(notes, emm$note)
  if (!is.null(offset) && !is.null(emm$grid)) {
    notes <- c(notes, sprintf(
      "Estimated marginal means are rates at %s = 1 (one unit of exposure).", offset))
  }
  emm_tab <- .emmeans_table(emm$grid, conf_level, protect = groups)
  notes <- c(notes, emm_tab$note)
  ph <- if (posthoc) {
    .emmeans_pairs(emm$grid, adjust = adjust, conf_level = conf_level,
                   ratios = TRUE)
  } else list(table = NULL,
              note = .no_posthoc_note(NROW(emm_tab$table),
                                      extra = "$effect_sizes still reports the incidence rate ratios from the model."))
  notes <- c(notes, ph$note)
  if (!is.null(ph$table) && "ratio" %in% names(ph$table)) {
    names(ph$table)[names(ph$table) == "ratio"] <- "IRR"
  }

  ## Coefficient-scale incidence rate ratios ---------------------------------
  ct <- .coef_table(fit, rv$matrix)
  z <- .coef_crit(fit, conf_level)
  effect_sizes <- data.frame(
    term = ct$term, IRR = exp(ct$estimate),
    conf_low = exp(ct$estimate - z * ct$se),
    conf_high = exp(ct$estimate + z * ct$se),
    p_value = ct$p_value, stringsAsFactors = FALSE)
  effect_sizes <- effect_sizes[effect_sizes$term != "(Intercept)", , drop = FALSE]
  row.names(effect_sizes) <- NULL

  ## Plots --------------------------------------------------------------------
  plot_list <- list()
  if (plots && !is.null(emm_tab$table)) {
    .say(verbose, "Building plots.")
    plot_list$emmeans <- .plot_emmeans(
      emm_tab$table, groups, estimate = "estimate",
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

  .new_fit(
    method       = sprintf("Analysis of deviance for a count response (%s regression)",
                           model_type),
    call         = cl,
    model        = fit,
    anova        = av$table,
    effect_sizes = effect_sizes,
    emmeans      = emm_tab$table,
    emmeans_object = emm$grid,
    posthoc      = ph$table,
    assumptions  = list(poisson_dispersion = dispersion,
                        model_dispersion = .dispersion(fit),
                        cell_counts = .cell_counts(d, groups)),
    plots        = plot_list,
    data_used    = d,
    n_removed    = prep$n_removed,
    conf_level   = conf_level,
    notes        = notes,
    extra        = list(model_type = model_type,
                        dispersion = dispersion,
                        model_dispersion = .dispersion(fit))
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
  invisible(TRUE)
}

#' Decide which count model to fit
#' @noRd
.choose_count_model <- function(model, dispersion, threshold) {
  notes <- character(0)
  if (model != "auto") {
    if (model == "negbin" && !requireNamespace("MASS", quietly = TRUE)) {
      .stopf("model = \"negbin\" requires the {MASS} package. Install it, or use model = \"quasipoisson\".")
    }
    return(list(model = model, notes = notes))
  }
  overdispersed <- is.finite(dispersion) && dispersion > threshold
  if (!overdispersed) {
    if (is.finite(dispersion)) {
      notes <- c(notes, sprintf(
        "Pearson dispersion is %.2f, at or below the threshold of %.2f, so a Poisson model was kept.",
        dispersion, threshold))
    }
    return(list(model = "poisson", notes = notes))
  }
  if (requireNamespace("MASS", quietly = TRUE)) {
    notes <- c(notes, sprintf(
      "Pearson dispersion is %.2f, above the threshold of %.2f, so a negative binomial model was fitted instead of Poisson. Set model = \"poisson\" to override.",
      dispersion, threshold))
    return(list(model = "negbin", notes = notes))
  }
  notes <- c(notes, sprintf(
    "Pearson dispersion is %.2f, above the threshold of %.2f. The {MASS} package is not installed, so a quasi-Poisson model was fitted instead of a negative binomial one.",
    dispersion, threshold))
  list(model = "quasipoisson", notes = notes)
}

#' Note sparse or empty cells
#' @noRd
.sparse_cell_note <- function(d, groups) {
  levs <- lapply(groups, function(g) levels(droplevels(as.factor(d[[g]]))))
  names(levs) <- groups
  full <- expand.grid(levs, stringsAsFactors = FALSE)
  observed <- .cell_counts(d, groups)
  n_possible <- nrow(full)
  n_observed <- nrow(observed)
  notes <- character(0)
  empty <- n_possible - n_observed
  small <- sum(observed$n < 5L)
  if (empty > 0L) {
    notes <- c(notes, sprintf(
      "%d of the %d possible group combination(s) contain no observations; contrasts involving them are not estimable.",
      empty, n_possible))
  }
  if (small > 0L) {
    notes <- c(notes, sprintf(
      "%d group combination(s) have fewer than 5 observations; inference for those cells is unstable.",
      small))
  }
  notes
}

#' Describe the negative binomial dispersion parameter honestly
#'
#' A theta in the hundreds or thousands is not an estimate of overdispersion, it
#' is the fit telling you there is none: the negative binomial has collapsed
#' towards the Poisson it generalises. Reporting it to three decimals as though
#' it were a parameter estimate invites the wrong reading.
#' @noRd
.theta_note <- function(nb) {
  th <- nb$theta; se <- nb$SE.theta
  warns <- attr(nb, "nb_warnings")
  if (!is.null(warns) && any(grepl("iteration limit|alternation limit", warns))) {
    return(sprintf(
      "The negative binomial dispersion parameter did not converge (theta reached %.3g before the iteration limit). Treat the model as unreliable and compare it with model = \"poisson\".",
      th))
  }
  if (!is.finite(th) || th > 1000) {
    return(sprintf(
      "The negative binomial dispersion parameter is very large (theta = %.3g), which means the data are not meaningfully overdispersed relative to Poisson and the negative binomial has collapsed towards it. model = \"poisson\" would give essentially the same answer more simply.",
      th))
  }
  if (is.finite(se) && se > 0 && th / se < 2) {
    return(sprintf(
      "Negative binomial dispersion parameter theta = %.3f (SE %.3f), which is not well determined: the data are only weakly overdispersed. The analysis of deviance treats theta as known, so its p-values are mildly anti-conservative.",
      th, se))
  }
  sprintf(
    "Negative binomial dispersion parameter theta = %.3f (SE %.3f). The analysis of deviance treats theta as known, so its p-values are mildly anti-conservative.",
    th, se)
}
