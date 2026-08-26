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
#' \strong{Separation.} Complete and quasi-complete separation are detected and
#' reported in \code{$notes}. Under separation, Wald odds ratios and their
#' intervals are meaningless even though \code{glm()} reports convergence, so
#' \code{ci_method = "profile"} is the default.
#'
#' @param data A data frame, or anything inheriting from one, such as a
#'   \code{data.table} or a tibble.
#' @param response Character. Name of the binary response column. May be a
#'   factor with two observed levels, a logical, or a numeric 0/1 vector.
#' @param groups Character vector. One or more grouping columns.
#' @param success Optional. The response level to model as the "success". By
#'   default the second level of the factor, which is what \code{glm()} uses.
#' @param reference Optional named list of reference levels for the grouping
#'   factors, for example \code{list(site = "north")}.
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
#'   the \pkg{sandwich} package.
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
#'   dropped for missing values.
#'
#' @return An \code{\link{anovatoolbox_fit}} object with \code{$model_stats}
#'   (AIC, BIC, the log-likelihood, McFadden's pseudo R squared, the level
#'   treated as a success, and \code{lr_vs_null}, the model-versus-null
#'   likelihood ratio test), and \code{$assumptions} holding
#'   \code{dispersion} and \code{proportions} (the observed proportion in each
#'   cell). \code{$effect_sizes} holds the odds ratios, \code{$emmeans} the
#'   estimated marginal probabilities and \code{$posthoc} the pairwise odds
#'   ratios.
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
  .check_flag(posthoc, "posthoc")
  .check_flag(plots, "plots")
  .check_flag(verbose, "verbose")

  .say(verbose, "Preparing data.")
  .check_weights_column(data, weights, response, groups)
  prep <- .prepare_frame(data, c(response, groups, weights), factors = groups,
                         keep_all = FALSE)
  d <- prep$data
  notes <- prep$notes
  .check_groups(d, groups, min_levels = 2L, min_n = 1L)

  wts <- .resolve_weights(d, weights)
  resp <- .prep_binary(d[[response]], response, success)
  d[[response]] <- resp$values
  success <- resp$success
  notes <- c(notes, resp$notes, sprintf(
    "Modelling P(%s = %s); the other level is the baseline.", response, success))

  if (!is.null(reference)) {
    d <- .set_references(d, groups, reference)
  }

  ## Fit ---------------------------------------------------------------------
  terms_rhs <- .group_terms(groups, interaction)
  fml <- stats::reformulate(terms_rhs, response = .bq(response))
  .say(verbose, "Fitting %s", paste(deparse(fml), collapse = " "))
  fitted <- .collect_conditions(.fit_with_contrasts(
    function() .fit_glm(fml, d, stats::binomial(), wts),
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
  notes <- c(notes, .check_separation(model))

  ## Analysis of deviance ----------------------------------------------------
  av <- .car_anova(model, type = type, test_statistic = test_statistic)
  notes <- c(notes, av$note)

  ## Null model, fitted once, used for both the LR test and pseudo R squared --
  null_model <- stats::glm(stats::reformulate("1", response = .bq(response)),
                           data = d, family = stats::binomial())
  ll_full <- stats::logLik(model); ll_null <- stats::logLik(null_model)
  lr_stat <- as.numeric(2 * (ll_full - ll_null))
  lr_df   <- attr(ll_full, "df") - attr(ll_null, "df")
  model_stats <- list(
    AIC = stats::AIC(model), BIC = stats::BIC(model),
    logLik = as.numeric(ll_full),
    mcfadden_r2 = as.numeric(1 - ll_full / ll_null),
    success_level = success,
    lr_vs_null = data.frame(
      statistic = lr_stat, df = lr_df,
      p_value = stats::pchisq(lr_stat, df = lr_df, lower.tail = FALSE),
      stringsAsFactors = FALSE)
  )

  ## Odds ratios -------------------------------------------------------------
  rv <- .robust_vcov(model, vcov_type)
  notes <- c(notes, rv$note)
  if (!is.null(rv$matrix) && ci_method == "profile") {
    ci_method <- "wald"
    notes <- c(notes, "Robust standard errors were requested, so the odds-ratio intervals are Wald intervals built from them; profile-likelihood intervals cannot use a sandwich covariance.")
  }
  eff <- .odds_ratios(model, conf_level = conf_level, ci_method = ci_method,
                      vcov_matrix = rv$matrix)
  notes <- c(notes, eff$note)

  ## Marginal probabilities and pairwise odds ratios -------------------------
  emm <- .emmeans_grid(model, groups, type = "response",
                       vcov_matrix = rv$matrix)
  notes <- c(notes, emm$note)
  emm_tab <- .emmeans_table(emm$grid, conf_level, protect = groups)
  notes <- c(notes, emm_tab$note)
  ph <- if (posthoc) {
    .emmeans_pairs(emm$grid, adjust = adjust, conf_level = conf_level,
                   ratios = TRUE)
  } else list(table = NULL,
              note = .no_posthoc_note(NROW(emm_tab$table),
                                      extra = "$effect_sizes still reports the odds ratios from the model."))
  notes <- c(notes, ph$note)

  ## Proportions and plot ----------------------------------------------------
  added <- .add_cell(d, groups)
  prop_tab <- .proportion_table(added$data, response, added$cell)
  plot_list <- list()
  if (plots) {
    .say(verbose, "Building plots.")
    plot_list$proportions <- .plot_proportions(
      prop_tab, added$cell, response, success,
      xlab = paste(groups, collapse = " : "))
    if (!is.null(emm_tab$table)) {
      plot_list$emmeans <- .plot_emmeans(
        emm_tab$table, groups, estimate = "estimate",
        lower = "conf_low", upper = "conf_high", conf_level = conf_level,
        title = sprintf("Estimated probability of %s = %s", response, success),
        ylab = "Probability")
    }
  }

  .new_fit(
    method       = "Analysis of deviance for a binary response (logistic regression)",
    call         = cl,
    model        = model,
    anova        = av$table,
    effect_sizes = eff$table,
    emmeans      = emm_tab$table,
    emmeans_object = emm$grid,
    posthoc      = ph$table,
    assumptions  = list(dispersion = .dispersion(model),
                        proportions = prop_tab),
    plots        = plot_list,
    data_used    = d,
    n_removed    = prep$n_removed,
    conf_level   = conf_level,
    notes        = notes,
    extra        = list(model_stats = model_stats)
  )
}

#' Coerce and validate a binary response
#'
#' A numeric vector containing only zeros still produces a two-level factor when
#' the levels are supplied, which is how an all-one-outcome column can slip
#' through unnoticed. Levels are therefore counted after dropping empty ones.
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
    x <- stats::relevel(x, ref = setdiff(observed, success))
  }
  list(values = x, success = levels(x)[2L], notes = notes)
}

#' Apply user-supplied reference levels to grouping factors
#' @noRd
.set_references <- function(d, groups, reference) {
  if (!is.list(reference) || is.null(names(reference))) {
    .stopf("`reference` must be a named list, for example list(%s = \"...\").",
           groups[1L])
  }
  unknown <- setdiff(names(reference), groups)
  if (length(unknown) > 0L) {
    .stopf("`reference` names must be grouping variables; %s %s not.",
           paste(sprintf("`%s`", unknown), collapse = ", "),
           if (length(unknown) == 1L) "is" else "are")
  }
  for (g in names(reference)) {
    ref <- as.character(reference[[g]])
    if (!ref %in% levels(d[[g]])) {
      .stopf("Reference level \"%s\" not found in `%s`. Levels: %s.",
             ref, g, paste(levels(d[[g]]), collapse = ", "))
    }
    if (is.ordered(d[[g]])) {
      lv <- levels(d[[g]])
      d[[g]] <- factor(as.character(d[[g]]), levels = c(ref, setdiff(lv, ref)),
                       ordered = TRUE)
    } else {
      d[[g]] <- stats::relevel(d[[g]], ref = ref)
    }
  }
  d
}

#' Odds ratios with profile-likelihood or Wald intervals
#' @noRd
.odds_ratios <- function(model, conf_level = 0.95, ci_method = "profile",
                         vcov_matrix = NULL) {
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
  out <- data.frame(
    term       = ct$term,
    odds_ratio = exp(ct$estimate),
    conf_low   = exp(ci[idx, 1L]),
    conf_high  = exp(ci[idx, 2L]),
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

#' Observed proportions of each response level within each cell
#' @noRd
.proportion_table <- function(d, response, cell) {
  tab <- table(d[[cell]], d[[response]])
  totals <- rowSums(tab)
  out <- as.data.frame(tab, stringsAsFactors = FALSE)
  names(out) <- c(cell, response, "n")
  out$group_total <- as.integer(totals[out[[cell]]])
  out$proportion <- ifelse(out$group_total > 0L, out$n / out$group_total, NA_real_)
  out[[cell]] <- factor(out[[cell]], levels = levels(d[[cell]]))
  out
}
