#' Welch's Analysis of Variance
#'
#' Compares the mean of a numeric response across groups without assuming equal
#' variances, using \code{\link[stats]{oneway.test}}. When more than two groups
#' are present, pairwise Welch t-tests with a multiplicity adjustment and
#' standardised mean differences are also returned.
#'
#' @details
#' Welch's test is the appropriate default for comparing means: it is barely
#' less powerful than the classical F test when variances are equal, and it
#' keeps its nominal error rate when they are not.
#'
#' Because the method explicitly permits unequal variances, normality is
#' assessed \emph{within each group} rather than on pooled residuals. Pooling
#' residuals from groups with different spreads produces a mixture distribution
#' that fails normality tests even when every group is perfectly normal, which
#' would argue against the very method being used.
#'
#' When several grouping variables are supplied they are combined into a single
#' cell factor and unused level combinations are dropped, so the number of
#' groups reported is the number that actually contain data.
#'
#' @param data A data frame, or anything inheriting from one, such as a
#'   \code{data.table} or a tibble.
#' @param response Character. Name of the numeric response column.
#' @param groups Character vector. One or more grouping columns.
#' @param conf_level Numeric in (0, 1). Level for every interval returned:
#'   group means, pairwise differences and standardised effect sizes. Default
#'   \code{0.95}.
#' @param adjust Character. Multiplicity adjustment for the pairwise
#'   comparisons, passed to \code{\link[stats]{p.adjust}}. One of
#'   \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"},
#'   \code{"BH"}, \code{"BY"}, \code{"fdr"} or \code{"none"}. Default
#'   \code{"holm"}. Note that \code{"tukey"} is not available here, as it is on
#'   the functions whose comparisons come from \pkg{emmeans}: these are Welch
#'   t-tests, not linear contrasts on a common error term.
#' @param hedges_correction Logical. Apply Hedges' small-sample correction to
#'   the standardised mean differences. Default \code{TRUE}, in which case the
#'   effect size column is named \code{hedges_g}; otherwise \code{cohens_d}.
#' @param posthoc Logical. Compute the pairwise comparisons. There are
#'   \code{choose(k, 2)} of them, so this is worth turning off when the number
#'   of groups is large. Default \code{TRUE}. Hedges' g is computed pair by
#'   pair alongside them, so \code{$effect_sizes} is empty when this is
#'   \code{FALSE}; \code{$notes} says so.
#' @param plots Logical. Build \pkg{ggplot2} objects. They are returned in
#'   \code{$plots}, never drawn. Default \code{TRUE}.
#' @param verbose Logical. Emit progress through \code{\link[base]{message}}.
#'   Default \code{FALSE}.
#'
#' @return An \code{\link{anovakit_fit}} object. \code{$anova} holds the
#'   Welch test, \code{$emmeans} the group means with intervals,
#'   \code{$posthoc} the pairwise comparisons with adjusted p-values and
#'   intervals, \code{$effect_sizes} the standardised mean differences (both
#'   \code{NULL} when \code{posthoc = FALSE}), and
#'   \code{$assumptions} a per-group Shapiro-Wilk table in \code{normality} and
#'   the ratio of the largest to the smallest group variance in
#'   \code{variance_ratio}. \code{$model} is the \code{htest} returned by
#'   \code{\link[stats]{oneway.test}}; there is no fitted model object and
#'   \code{$emmeans_object} is \code{NULL}.
#'
#' @seealso \code{\link{anova_kw}} for a rank-based alternative,
#'   \code{\link{anova_ancova}} to adjust for a covariate.
#'
#' @examples
#' set.seed(123)
#' d <- data.frame(
#'   group = rep(c("A", "B", "C"), each = 30),
#'   value = c(rnorm(30, 10, 2), rnorm(30, 12, 2.5), rnorm(30, 9, 1.5))
#' )
#' fit <- anova_welch(d, "value", "group")
#' fit
#' fit$posthoc
#'
#' # Two grouping variables are combined into one cell factor
#' d$site <- rep(c("north", "south"), 45)
#' anova_welch(d, "value", c("group", "site"), plots = FALSE)$anova
#'
#' @export
anova_welch <- function(data, response, groups,
                        conf_level = 0.95,
                        adjust = c("holm", "hochberg", "hommel", "bonferroni",
                                   "BH", "BY", "fdr", "none"),
                        hedges_correction = TRUE,
                        posthoc = TRUE,
                        plots = TRUE,
                        verbose = FALSE) {
  cl <- match.call()
  adjust <- match.arg(adjust)
  data <- .as_df(data)
  .check_name(response, "response")
  .check_names(groups, "groups")
  .check_columns(data, response, "Response column")
  .check_columns(data, groups, "Grouping column(s)")
  .check_numeric_col(data, response)
  .check_conf_level(conf_level)
  .check_flag(hedges_correction, "hedges_correction")
  .check_flag(posthoc, "posthoc")
  .check_flag(plots, "plots")
  .check_flag(verbose, "verbose")

  .say(verbose, "Preparing data.")
  prep <- .prepare_frame(data, c(response, groups), factors = groups,
                         keep_all = FALSE)
  d <- prep$data
  notes <- prep$notes
  .check_groups(d, groups, min_levels = 1L, min_n = 2L)

  added <- .add_cell(d, groups)
  d <- added$data
  cell <- added$cell
  n_groups <- nlevels(d[[cell]])
  if (n_groups < 2L) {
    .stopf("At least two populated groups are required; found %d.", n_groups)
  }

  group_var <- tapply(d[[response]], d[[cell]], stats::var)
  degenerate <- names(group_var)[!is.finite(group_var) | group_var == 0]
  if (length(degenerate) > 0L) {
    .stopf("Welch's test needs a positive variance in every group, and %s none: %s. Every observation there takes the same value, so consider dropping the group or using anova_kw().",
           if (length(degenerate) == 1L) "this group has" else "these groups have",
           paste(degenerate, collapse = ", "))
  }

  ## Omnibus test ------------------------------------------------------------
  .say(verbose, "Running Welch's ANOVA on %d groups.", n_groups)
  fml <- .formula(response, cell)
  welch <- stats::oneway.test(fml, data = d, var.equal = FALSE)
  anova_tab <- data.frame(
    term      = paste(groups, collapse = " : "),
    statistic = unname(welch$statistic),
    num_df    = unname(welch$parameter[1L]),
    den_df    = unname(welch$parameter[2L]),
    p_value   = unname(welch$p.value),
    stringsAsFactors = FALSE
  )

  ## Group summaries ---------------------------------------------------------
  summary_stats <- .summary_stats(d[[response]], d[[cell]], conf_level)

  ## Assumptions -------------------------------------------------------------
  normality <- .normality_by_group(d[[response]], d[[cell]])
  if (all(is.na(normality$p_value))) {
    notes <- c(notes, "Per-group Shapiro-Wilk could not be computed for any group; see $assumptions$normality for the group sizes.")
  }
  sds <- summary_stats$sd
  variance_ratio <- if (length(sds) >= 2L && all(is.finite(sds))) {
    if (min(sds) > 0) (max(sds) / min(sds))^2 else Inf
  } else NA_real_
  if (is.finite(variance_ratio) && variance_ratio > 4) {
    notes <- c(notes, sprintf(
      "Largest group variance is %.1f times the smallest. Welch's test handles this; the classical F test would not.",
      variance_ratio))
  }

  ## Pairwise comparisons ----------------------------------------------------
  posthoc_tab <- NULL
  effect_sizes <- NULL
  if (!posthoc) {
    # Hedges' g is computed pair by pair alongside the comparisons, so
    # suppressing them suppresses it too. Say so: an $effect_sizes of NULL is
    # otherwise indistinguishable from one that could not be computed.
    notes <- c(notes, .no_posthoc_note(
      n_groups,
      extra = "$effect_sizes is empty for the same reason: Hedges' g is a per-pair quantity, computed with the comparisons."))
  } else if (n_groups >= 2L) {
    pw <- .welch_pairwise(d[[response]], d[[cell]], conf_level = conf_level,
                          adjust = adjust,
                          hedges_correction = hedges_correction)
    posthoc_tab <- pw$posthoc
    effect_sizes <- pw$effect_sizes
    if (n_groups == 2L) {
      notes <- c(notes, "With two groups the pairwise comparison is the omnibus test; no multiplicity adjustment was needed.")
    }
  }

  ## Plots -------------------------------------------------------------------
  plot_list <- list()
  if (plots) {
    .say(verbose, "Building plots.")
    plot_list$means <- .plot_means(summary_stats, response, conf_level,
                                   xlab = paste(groups, collapse = " : "))
    plot_list$box <- .plot_box(d, response, cell,
                               xlab = paste(groups, collapse = " : "))
    resid_within <- unlist(lapply(levels(d[[cell]]), function(g) {
      v <- d[[response]][d[[cell]] == g]
      v - mean(v, na.rm = TRUE)
    }), use.names = FALSE)
    plot_list$qq <- .plot_qq(
      resid_within,
      title = "Normal Q-Q plot of within-group deviations",
      subtitle = "Deviations from each group's own mean, as Welch's test assumes")
  }

  .new_fit(
    method       = "Welch's analysis of variance",
    call         = cl,
    model        = welch,
    anova        = anova_tab,
    effect_sizes = effect_sizes,
    emmeans      = summary_stats,
    posthoc      = posthoc_tab,
    assumptions  = list(normality = normality,
                        variance_ratio = variance_ratio),
    plots        = plot_list,
    data_used    = d,
    n_removed    = prep$n_removed,
    conf_level   = conf_level,
    notes        = notes
  )
}

#' Pairwise Welch t-tests with intervals and standardised differences
#' @noRd
.welch_pairwise <- function(values, cell, conf_level = 0.95, adjust = "holm",
                            hedges_correction = TRUE) {
  lv <- levels(droplevels(as.factor(cell)))
  combos <- utils::combn(lv, 2L, simplify = FALSE)
  rows <- lapply(combos, function(pair) {
    x <- values[cell == pair[1L]]; x <- x[!is.na(x)]
    y <- values[cell == pair[2L]]; y <- y[!is.na(y)]
    n1 <- length(x); n2 <- length(y)
    if (n1 < 2L || n2 < 2L) {
      return(data.frame(group1 = pair[1L], group2 = pair[2L],
                        difference = mean(x) - mean(y),
                        conf_low = NA_real_, conf_high = NA_real_,
                        statistic = NA_real_, df = NA_real_,
                        p_value = NA_real_,
                        effect = NA_real_, effect_low = NA_real_,
                        effect_high = NA_real_,
                        stringsAsFactors = FALSE))
    }
    tt <- stats::t.test(x, y, var.equal = FALSE, conf.level = conf_level)
    es <- .std_mean_diff(x, y, conf_level = conf_level,
                         hedges_correction = hedges_correction)
    data.frame(
      group1 = pair[1L], group2 = pair[2L],
      difference = unname(diff(rev(tt$estimate))),
      conf_low = tt$conf.int[1L], conf_high = tt$conf.int[2L],
      statistic = unname(tt$statistic), df = unname(tt$parameter),
      p_value = unname(tt$p.value),
      effect = es$estimate, effect_low = es$conf_low,
      effect_high = es$conf_high,
      stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  tab$p_adjusted <- stats::p.adjust(tab$p_value, method = adjust)
  tab$adjustment <- adjust

  es_name <- if (hedges_correction) "hedges_g" else "cohens_d"
  effect_sizes <- data.frame(
    group1 = tab$group1, group2 = tab$group2,
    estimate = tab$effect, conf_low = tab$effect_low,
    conf_high = tab$effect_high,
    magnitude = .effect_magnitude(tab$effect),
    stringsAsFactors = FALSE)
  names(effect_sizes)[names(effect_sizes) == "estimate"] <- es_name

  posthoc <- tab[, c("group1", "group2", "difference", "conf_low", "conf_high",
                     "statistic", "df", "p_value", "p_adjusted", "adjustment")]
  row.names(posthoc) <- NULL
  row.names(effect_sizes) <- NULL
  list(posthoc = posthoc, effect_sizes = effect_sizes)
}

#' Standardised mean difference with an approximate confidence interval
#'
#' Cohen's d, optionally with Hedges' small-sample correction. The interval uses
#' the standard large-sample variance of d.
#' @noRd
.std_mean_diff <- function(x, y, conf_level = 0.95, hedges_correction = TRUE) {
  n1 <- length(x); n2 <- length(y)
  s_pooled <- sqrt(((n1 - 1) * stats::var(x) + (n2 - 1) * stats::var(y)) /
                     (n1 + n2 - 2))
  if (!is.finite(s_pooled) || s_pooled == 0) {
    return(list(estimate = NA_real_, conf_low = NA_real_, conf_high = NA_real_))
  }
  d <- (mean(x) - mean(y)) / s_pooled
  J <- if (hedges_correction) 1 - 3 / (4 * (n1 + n2) - 9) else 1
  est <- d * J
  # Var(g) = J^2 * Var(d), and Var(d) is a function of the UNCORRECTED d.
  se <- J * sqrt((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2)))
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  list(estimate = est, conf_low = est - z * se, conf_high = est + z * se)
}

#' Conventional magnitude labels for a standardised mean difference
#' @noRd
.effect_magnitude <- function(d) {
  out <- rep(NA_character_, length(d))
  ok <- is.finite(d)
  a <- abs(d[ok])
  out[ok] <- ifelse(a < 0.2, "negligible",
             ifelse(a < 0.5, "small",
             ifelse(a < 0.8, "medium", "large")))
  out
}
