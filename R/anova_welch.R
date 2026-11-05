#' Welch's Analysis of Variance
#'
#' Compares the mean of a numeric response across groups without assuming equal
#' variances, using \code{\link[stats]{oneway.test}}. Pairwise Welch t-tests
#' with a multiplicity adjustment, and standardised mean differences, are also
#' returned for any number of groups from two upwards; with two groups the
#' single comparison is the omnibus test.
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
#' would argue against the very method being used. For the same reason the Q-Q
#' plot shows each observation's deviation from its group mean divided by its
#' group's standard deviation, so every group contributes on its own scale.
#' A group whose Shapiro-Wilk test could not be run (fewer than 3 or more than
#' 5000 observations) is named in \code{$notes}, with the reason.
#'
#' \strong{Standardised mean differences.} The pairwise effect sizes do not
#' assume equal variances either. Each difference in means is divided by the
#' average-variance standard deviation of its two groups,
#' \eqn{s^* = \sqrt{(s_1^2 + s_2^2)/2}}{s* = sqrt((s1^2 + s2^2) / 2)}, rather
#' than by the pooled standard deviation, whose meaning and sampling variance
#' both depend on equal variances. The interval is for the population value
#' \eqn{\delta^* = (\mu_1 - \mu_2) / \sqrt{(\sigma_1^2 + \sigma_2^2)/2}}{delta*
#' = (mu1 - mu2) / sqrt((sigma1^2 + sigma2^2) / 2)}. It inverts the noncentral
#' t distribution of Welch's statistic on Welch's degrees of freedom and
#' rescales the bounds on the noncentrality to the average-variance standard
#' deviation, which is the interval \code{effectsize::cohens_d(pooled_sd =
#' FALSE)} reports. In simulation (variance ratios from 1/4 to 16 in either
#' group, group sizes from 2 to 50, true differences of 0 and 0.8) a 95\%
#' interval covered 93-97\% whenever the smaller group had at least five
#' observations, 93-98\% for equal groups of two or three, and 86-95\% when
#' a group of two or three sat beside one ten or more times its size, where
#' no interval does well: with two observations a group's standard deviation
#' is barely estimated. Bonett's (2008) normal-theory interval for the same
#' quantity was also examined; it was wider than needed for equal small groups
#' (99-100\% at two per group) and covered less in the unbalanced cases
#' (81-94\%).
#'
#' With \code{hedges_correction = TRUE} the estimate is multiplied by the
#' small-sample bias correction \eqn{J(\nu) = \Gamma(\nu/2) / (\sqrt{\nu/2}\,
#' \Gamma((\nu-1)/2))}{J(nu) = gamma(nu/2) / (sqrt(nu/2) gamma((nu-1)/2))},
#' evaluated at the Satterthwaite degrees of freedom of the standardiser,
#' \eqn{\nu = (s_1^2 + s_2^2)^2 / (s_1^4/(n_1-1) + s_2^4/(n_2-1))}{nu = (s1^2 +
#' s2^2)^2 / (s1^4/(n1-1) + s2^4/(n2-1))} (Delacre et al., 2021). When the two
#' sample variances are equal and the groups are the same size this is the
#' familiar \eqn{n_1 + n_2 - 2}{n1 + n2 - 2}; when the variances differ it
#' removes the bias that the pooled-variance correction leaves. The correction is a property of the point estimate, so
#' the interval, which is for \eqn{\delta^*}{delta*} itself, is the same with or
#' without it. The corrected estimate can fall outside that interval at any
#' \code{conf_level} when \eqn{\nu}{nu} is near 1, which happens when a group
#' of two or three holds most of the variance: \eqn{J(\nu)}{J(nu)} then
#' approaches 0 and pulls the estimate towards zero. With few observations and
#' a low \code{conf_level} the estimate, corrected or not, can also fall
#' outside, because the interval inverts a skewed distribution. \code{$notes}
#' says so when either happens.
#'
#' The interval is found with an exact noncentral t distribution function:
#' beyond a noncentrality of 30 it is computed by numerical integration,
#' because \code{stats::pt()} switches to a normal approximation for a large
#' noncentrality that is badly wrong at small degrees of freedom.
#'
#' A pair whose standard error is below the floating-point resolution of its
#' means, which \code{stats::t.test()} refuses as "essentially constant", is
#' reported as \code{NA} with a note.
#'
#' The standardiser is recorded in the \code{standardiser} column of
#' \code{$effect_sizes}. The magnitude labels use Cohen's conventional
#' thresholds of 0.2, 0.5 and 0.8 on the absolute value of the estimate; when
#' the estimate lies outside its own interval the label is taken from the
#' nearer interval bound instead, so that it does not describe a value the
#' interval rules out.
#'
#' \strong{Several grouping variables.} They are combined into a single cell
#' factor and unused level combinations are dropped, so the number of groups
#' reported is the number that actually contain data. The result is one
#' omnibus test across all populated cells, labelled \code{"A x B cells"} in
#' \code{$anova}: it cannot separate a main effect of one factor from a main
#' effect of another, and it cannot test an interaction. A note to that effect
#' is added to the returned object. For a factorial analysis use
#' \code{\link{anova_glm}} with \code{interaction = TRUE}.
#'
#' @param data A data frame, or anything inheriting from one, such as a
#'   \code{data.table} or a tibble.
#' @param response Character. Name of the numeric response column.
#' @param groups Character vector. One or more grouping columns.
#' @param conf_level Numeric in (0, 1). Level for every interval returned:
#'   group means, pairwise differences and standardised effect sizes. Default
#'   \code{0.95}.
#' @param adjust Character. Multiplicity adjustment for the pairwise
#'   comparisons' p-values, passed to \code{\link[stats]{p.adjust}}. One of
#'   \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"},
#'   \code{"BH"}, \code{"BY"}, \code{"fdr"} or \code{"none"}, spelled out in
#'   full. Default \code{"holm"}. Note that \code{"tukey"} is not available
#'   here, as it is on the functions whose comparisons come from
#'   \pkg{emmeans}: these are Welch t-tests, not linear contrasts on a common
#'   error term.
#' @param hedges_correction Logical. Apply the small-sample bias correction to
#'   the standardised mean differences (see Details). Default \code{TRUE}, in
#'   which case the effect size column is named \code{hedges_g}; otherwise
#'   \code{cohens_d}. Both use the average-variance standardiser.
#' @param posthoc Logical. Compute the pairwise comparisons. There are
#'   \code{choose(k, 2)} of them, so this is worth turning off when the number
#'   of groups is large. Default \code{TRUE}. The standardised mean
#'   differences are computed pair by pair alongside them, so
#'   \code{$effect_sizes} is empty when this is \code{FALSE}; \code{$notes}
#'   says so.
#' @param plots Logical. Build \pkg{ggplot2} objects. They are returned in
#'   \code{$plots}, never drawn. Default \code{TRUE}.
#' @param verbose Logical. Emit progress through \code{\link[base]{message}}.
#'   Default \code{FALSE}.
#'
#' @return An \code{\link{anovakit_fit}} object.
#'   \itemize{
#'     \item \code{$anova}: the Welch test, one row. Its \code{term} is the
#'       grouping column, or \code{"A x B cells"} when several were combined.
#'     \item \code{$emmeans}: the group means with t intervals (raw group
#'       summaries; no model is fitted).
#'     \item \code{$posthoc}: one row per pair of groups, with \code{group1},
#'       \code{group2}, \code{difference} (the mean of \code{group1} minus the
#'       mean of \code{group2}), \code{conf_low} and \code{conf_high} (the
#'       Welch interval for that difference at \code{conf_level}; an
#'       unadjusted, per-comparison interval, not a simultaneous one),
#'       \code{statistic} (Welch's t, same sign as \code{difference}),
#'       \code{df}, \code{p_value} (unadjusted), \code{p_adjusted} (adjusted by
#'       \code{adjust}) and \code{adjustment}. This is the same column layout
#'       as the functions built on \pkg{emmeans}, whose intervals are however
#'       adjusted.
#'     \item \code{$effect_sizes}: one row per pair, with \code{group1},
#'       \code{group2}, \code{hedges_g} (or \code{cohens_d}), the same
#'       direction as \code{difference}, its interval \code{conf_low} and
#'       \code{conf_high}, \code{magnitude}, and \code{standardiser}, which
#'       names the standard deviation used.
#'     \item \code{$assumptions}: a per-group Shapiro-Wilk table in
#'       \code{normality}, with a \code{note} column giving the reason for any
#'       group that was not tested, and the ratio of the largest to the
#'       smallest group variance in \code{variance_ratio}.
#'     \item \code{$model}: the \code{htest} returned by
#'       \code{\link[stats]{oneway.test}}; there is no fitted model object and
#'       \code{$emmeans_object} is \code{NULL}.
#'   }
#'   \code{$posthoc} and \code{$effect_sizes} are \code{NULL} when
#'   \code{posthoc = FALSE}.
#'
#' @references
#' Bonett, D. G. (2008). Confidence intervals for standardized linear
#' contrasts of means. \emph{Psychological Methods}, 13(2), 99-109.
#'
#' Delacre, M., Lakens, D., Ley, C., Liu, L., & Leys, C. (2021). Why Hedges'
#' g*s based on the non-pooled standard deviation should be reported with
#' Welch's t-test. \emph{PsyArXiv}.
#'
#' Welch, B. L. (1951). On the comparison of several mean values: an
#' alternative approach. \emph{Biometrika}, 38(3/4), 330-336.
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
#' # Standardised mean differences on the average-variance standard deviation,
#' # with an interval that does not assume equal variances
#' fit$effect_sizes
#'
#' # Two grouping variables are combined into one cell factor
#' d$site <- rep(c("north", "south"), 45)
#' anova_welch(d, "value", c("group", "site"), plots = FALSE)$anova
#'
#' @export
anova_welch <- function(data, response, groups,
                        conf_level = 0.95,
                        adjust = "holm",
                        hedges_correction = TRUE,
                        posthoc = TRUE,
                        plots = TRUE,
                        verbose = FALSE) {
  cl <- match.call()
  .check_adjust(adjust, stats::p.adjust.methods, "adjust")
  data <- .as_df(data)
  .check_name(response, "response")
  .check_names(groups, "groups")
  .check_columns(data, response, "Response column")
  .check_columns(data, groups, "Grouping column(s)")
  .check_numeric_col(data, response)
  .check_roles(list(`the response` = response, `a grouping variable` = groups))
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
  # A zero or non-finite variance means constant data only when the values
  # really are all equal; distinct values far from 1 in magnitude (around
  # 1e-200 or 1e200) can also give one, because squaring them under- or
  # overflows. That needs a rescale, not dropping the group.
  n_distinct <- tapply(d[[response]], d[[cell]], function(v) length(unique(v)))
  bad_var <- !is.finite(group_var) | group_var == 0
  degenerate <- levels(d[[cell]])[bad_var & n_distinct == 1L]
  if (length(degenerate) > 0L) {
    .stopf("Welch's test needs a positive variance in every group, and %s none: %s. Every observation there takes the same value, so consider dropping the group or using anova_kw().",
           if (length(degenerate) == 1L) "this group has" else "these groups have",
           paste(sprintf("\"%s\"", degenerate), collapse = ", "))
  }
  unrepresentable <- levels(d[[cell]])[bad_var & n_distinct > 1L]
  if (length(unrepresentable) > 0L) {
    .stopf("The variance of %s %s cannot be represented in double precision: the values differ, but at their scale the variance underflows to zero or overflows to infinity. Rescale the response (multiply or divide it by a power of ten) and rerun.",
           if (length(unrepresentable) == 1L) "group" else "groups",
           paste(sprintf("\"%s\"", unrepresentable), collapse = ", "))
  }

  term_label <- .cells_label(groups)
  if (length(groups) > 1L) {
    notes <- c(notes, sprintf(
      "The %d grouping variables were combined into %d cells and tested with a single omnibus Welch test. This cannot separate main effects or test an interaction; for that, use anova_glm() with interaction = TRUE.",
      length(groups), n_groups))
  }

  ## Omnibus test ------------------------------------------------------------
  .say(verbose, "Running Welch's ANOVA on %d groups.", n_groups)
  fml <- .formula(response, cell)
  welch <- stats::oneway.test(fml, data = d, var.equal = FALSE)
  # oneway.test names the internal cell column; name what the user supplied.
  welch$data.name <- sprintf("%s and %s", response, term_label)
  anova_tab <- data.frame(
    term      = term_label,
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
  notes <- c(notes, .normality_skip_notes(normality))
  group_var <- as.numeric(group_var)
  variance_ratio <- if (length(group_var) >= 2L && all(is.finite(group_var))) {
    if (min(group_var) > 0) max(group_var) / min(group_var) else Inf
  } else NA_real_
  if (is.finite(variance_ratio) && variance_ratio > 4) {
    notes <- c(notes, sprintf(
      "Largest group variance is %.1f times the smallest. Welch's test handles this; the classical F test would not.",
      variance_ratio))
  }

  ## Pairwise comparisons ----------------------------------------------------
  posthoc_tab <- NULL
  effect_sizes <- NULL
  es_label <- if (hedges_correction) "Hedges' g" else "Cohen's d"
  if (!posthoc) {
    # The standardised differences are computed pair by pair alongside the
    # comparisons, so suppressing them suppresses those too. Say so: an
    # $effect_sizes of NULL is otherwise indistinguishable from one that could
    # not be computed.
    notes <- c(notes, .no_posthoc_note(
      n_groups,
      extra = sprintf("$effect_sizes is empty for the same reason: %s is a per-pair quantity, computed with the comparisons.",
                      es_label)))
  } else {
    pw <- .welch_pairwise(d[[response]], d[[cell]], conf_level = conf_level,
                          adjust = adjust,
                          hedges_correction = hedges_correction)
    posthoc_tab <- pw$posthoc
    effect_sizes <- pw$effect_sizes
    notes <- c(notes, pw$notes)
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
    plot_list$qq <- .plot_qq(
      .standardise_within(d[[response]], d[[cell]]),
      title = "Normal Q-Q plot of within-group standardised deviations",
      subtitle = "(value - group mean) / group SD: each group on its own scale, as Welch's test allows")
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

#' The label of an omnibus test across the combined cells of several factors
#'
#' \code{"A : B"} is the conventional notation for an interaction term, which
#' a one-way test across the cells is not.
#' @noRd
.cells_label <- function(groups) {
  if (length(groups) == 1L) return(groups)
  sprintf("%s cells", paste(groups, collapse = " x "))
}

#' Deviations from each group's mean in units of that group's SD
#'
#' Non-finite values and groups with fewer than two finite values or no
#' spread contribute nothing.
#' @noRd
.standardise_within <- function(values, cell) {
  cell <- droplevels(as.factor(cell))
  unlist(lapply(levels(cell), function(g) {
    v <- values[cell == g]
    v <- v[is.finite(v)]
    if (length(v) < 2L) return(numeric(0))
    s <- stats::sd(v)
    if (!is.finite(s) || s == 0) return(numeric(0))
    (v - mean(v)) / s
  }), use.names = FALSE)
}

#' Notes naming the groups whose Shapiro-Wilk test was skipped, and why
#' @noRd
.normality_skip_notes <- function(normality) {
  if (is.null(normality) || !"note" %in% names(normality)) return(character(0))
  skipped <- !is.na(normality$note)
  if (!any(skipped)) return(character(0))
  reasons <- sub("^Shapiro-Wilk skipped: ", "", normality$note[skipped])
  grp <- as.character(normality$group[skipped])
  vapply(unique(reasons), function(r) {
    g <- grp[reasons == r]
    sprintf("Shapiro-Wilk was not run for group%s %s: %s",
            if (length(g) == 1L) "" else "s",
            .abbrev(sprintf("\"%s\"", g)), r)
  }, character(1), USE.NAMES = FALSE)
}

#' Pairwise Welch t-tests with intervals and standardised differences
#'
#' Welch's t, degrees of freedom, p-value and interval are computed for every
#' pair at once, in the closed form \code{stats::t.test()} uses, from each
#' group's mean, variance and size; calling \code{t.test()} and building a
#' data frame pair by pair made many groups very slow. Only the effect-size
#' interval, which inverts a noncentral t, is found pair by pair.
#'
#' \code{t.test()} refuses a pair whose standard error is below
#' \code{10 * .Machine$double.eps} times the larger absolute mean ("data are
#' essentially constant"): the spread is then at the floating-point
#' resolution of the means and any t computed from it is rounding noise. Such
#' a pair is left as \code{NA} with a note rather than given that noise.
#' @noRd
.welch_pairwise <- function(values, cell, conf_level = 0.95, adjust = "holm",
                            hedges_correction = TRUE) {
  cell <- droplevels(as.factor(cell))
  lv <- levels(cell)
  # Split by position, not by name: a level can be the empty string, which
  # cannot be used to index a named list.
  by_group <- lapply(split(values, cell), function(v) v[!is.na(v)])
  n <- vapply(by_group, length, numeric(1), USE.NAMES = FALSE)
  m <- vapply(by_group, mean, numeric(1), USE.NAMES = FALSE)
  v <- vapply(by_group, function(x) if (length(x) < 2L) NA_real_ else stats::var(x),
              numeric(1), USE.NAMES = FALSE)
  idx <- utils::combn(length(lv), 2L)
  i <- idx[1L, ]; j <- idx[2L, ]

  vi <- v[i] / n[i]; vj <- v[j] / n[j]
  se <- sqrt(vi + vj)
  difference <- m[i] - m[j]
  statistic <- difference / se
  df <- se^4 / (vi^2 / (n[i] - 1) + vj^2 / (n[j] - 1))
  p_value <- 2 * stats::pt(-abs(statistic), df)
  q <- stats::qt(1 - (1 - conf_level) / 2, df)
  conf_low <- (statistic - q) * se
  conf_high <- (statistic + q) * se

  sized <- n[i] >= 2L & n[j] >= 2L
  unresolved <- sized & se < 10 * .Machine$double.eps * pmax(abs(m[i]), abs(m[j]))
  ok <- sized & !unresolved
  conf_low[!ok] <- NA_real_; conf_high[!ok] <- NA_real_
  statistic[!ok] <- NA_real_; df[!ok] <- NA_real_; p_value[!ok] <- NA_real_

  effect <- effect_low <- effect_high <- rep(NA_real_, length(i))
  for (k in which(ok)) {
    es <- .std_mean_diff(by_group[[i[k]]], by_group[[j[k]]],
                         conf_level = conf_level,
                         hedges_correction = hedges_correction)
    effect[k] <- es$estimate
    effect_low[k] <- es$conf_low
    effect_high[k] <- es$conf_high
  }

  tab <- data.frame(group1 = lv[i], group2 = lv[j],
                    difference = difference,
                    conf_low = conf_low, conf_high = conf_high,
                    statistic = statistic, df = df, p_value = p_value,
                    stringsAsFactors = FALSE)
  tab$p_adjusted <- stats::p.adjust(tab$p_value, method = adjust)
  tab$adjustment <- adjust

  # An estimate outside its own interval would get a magnitude label that the
  # interval rules out ("negligible" beside an interval excluding everything
  # smaller than 1.8 in absolute value), so such a label is taken from the
  # interval bound nearest the estimate.
  outside <- which(effect < effect_low | effect > effect_high)
  labelled <- effect
  labelled[outside] <- pmin(pmax(effect[outside], effect_low[outside]),
                            effect_high[outside])

  es_name <- if (hedges_correction) "hedges_g" else "cohens_d"
  effect_sizes <- data.frame(
    group1 = tab$group1, group2 = tab$group2,
    estimate = effect, conf_low = effect_low,
    conf_high = effect_high,
    magnitude = .effect_magnitude(labelled),
    standardiser = "sqrt((s1^2 + s2^2) / 2)",
    stringsAsFactors = FALSE)
  names(effect_sizes)[names(effect_sizes) == "estimate"] <- es_name

  notes <- character(0)
  if (any(unresolved)) {
    notes <- c(notes, sprintf(
      "For %s, the spread within the groups is at the floating-point resolution of their means, so a t statistic computed from it would be rounding noise; that comparison and its effect size are left as NA. Subtracting a constant close to the values from the response before the analysis avoids this.",
      .abbrev(sprintf("%s vs %s", tab$group1[unresolved],
                      tab$group2[unresolved]), 5L)))
  }
  if (length(outside) > 0L) {
    pairs <- .abbrev(sprintf("%s vs %s", tab$group1[outside],
                             tab$group2[outside]), 5L)
    notes <- c(notes, if (hedges_correction) {
      sprintf(
        "For %s, Hedges' g lies outside its %s%% interval, so the magnitude label is taken from the nearer interval bound. The bias correction shrinks the estimate towards zero, and shrinks it strongly when the standardiser's degrees of freedom are near 1 (a group of two or three holding most of the variance), while the interval is for the population value and is not shrunk; with few observations and a low conf_level the interval, which inverts a skewed noncentral t distribution, can also miss the estimate. Cohen's d (hedges_correction = FALSE) is the uncorrected estimate.",
        pairs, .pct(conf_level))
    } else {
      sprintf(
        "For %s, Cohen's d lies outside its %s%% interval, so the magnitude label is taken from the nearer interval bound. The interval inverts the skewed noncentral t distribution of Welch's statistic, so with few degrees of freedom and a low conf_level it need not contain the point estimate.",
        pairs, .pct(conf_level))
    })
  }

  posthoc <- tab[, c("group1", "group2", "difference", "conf_low", "conf_high",
                     "statistic", "df", "p_value", "p_adjusted", "adjustment")]
  row.names(posthoc) <- NULL
  row.names(effect_sizes) <- NULL
  list(posthoc = posthoc, effect_sizes = effect_sizes, notes = notes)
}

#' Standardised mean difference that does not assume equal variances
#'
#' The difference in means over the average-variance standard deviation
#' \code{sqrt((s1^2 + s2^2) / 2)}, with the noncentral-t interval of Welch's
#' statistic rescaled to it, for the population value. With \code{hedges_correction} the
#' estimate (not the interval) is multiplied by the bias correction J at the
#' Satterthwaite degrees of freedom of the standardiser (Delacre et al.,
#' 2021).
#' @return list with \code{estimate}, \code{conf_low}, \code{conf_high},
#'   \code{se} (the scale of Welch's statistic on the standardised scale,
#'   \code{se_w / s*}), \code{standardiser} and
#'   \code{df} (of the standardiser).
#' @noRd
.std_mean_diff <- function(x, y, conf_level = 0.95, hedges_correction = TRUE) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  n1 <- length(x); n2 <- length(y)
  na <- list(estimate = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
             se = NA_real_, standardiser = NA_real_, df = NA_real_)
  if (n1 < 2L || n2 < 2L) return(na)
  v1 <- stats::var(x); v2 <- stats::var(y)
  df1 <- n1 - 1; df2 <- n2 - 1
  s <- sqrt((v1 + v2) / 2)
  if (!is.finite(s) || s == 0) return(na)
  d <- (mean(x) - mean(y)) / s
  # The interval inverts the noncentral t distribution of Welch's statistic
  # on Welch's degrees of freedom, and rescales the bounds on its noncentrality
  # to the average-variance standardiser: d* = t * se_w / s*.
  se_w <- sqrt(v1 / n1 + v2 / n2)
  t_w <- (mean(x) - mean(y)) / se_w
  df_w <- se_w^4 / ((v1 / n1)^2 / df1 + (v2 / n2)^2 / df2)
  ci <- .ncp_t_interval(t_w, df_w, conf_level) * se_w / s
  nu <- (v1 + v2)^2 / (v1^2 / df1 + v2^2 / df2)
  J <- if (hedges_correction) .hedges_j(nu) else 1
  list(estimate = J * d, conf_low = ci[1L], conf_high = ci[2L],
       se = se_w / s, standardiser = s, df = nu)
}

#' Confidence interval for the noncentrality of a t statistic
#'
#' The bounds are the noncentrality parameters at which the observed t is the
#' upper and the lower (1 - level) / 2 quantile. The distribution function is
#' \code{.pnt()}, not \code{stats::pt()}, whose noncentral algorithm is not
#' accurate for a large noncentrality. \code{pt()} also warns that full
#' precision may not have been achieved when the search strays far into a
#' tail, where the probability is far from the target and the shortfall does
#' not move the root, so those warnings are not passed on. They are muffled
#' once per root rather than per evaluation, which is much cheaper.
#' @noRd
.ncp_t_interval <- function(t, df, conf_level) {
  a <- (1 - conf_level) / 2
  f <- function(ncp, p) .pnt(t, df, ncp) - p
  # Start from the normal approximation to each root, t - z_p * s; uniroot()
  # widens the bracket when that misses, as it does at small df. Starting close
  # saves evaluations, which matters with many pairs.
  s <- sqrt(1 + t^2 / (2 * df))
  root <- function(p) suppressWarnings(tryCatch(
    stats::uniroot(f, t - stats::qnorm(p) * s + c(-1, 1) * s, p = p,
                   extendInt = "yes", tol = 1e-10)$root,
    error = function(e) NA_real_))
  c(root(1 - a), root(a))
}

#' Distribution function of the noncentral t, accurate for any noncentrality
#'
#' For |ncp| > 37.62 \code{stats::pt()} switches to a normal approximation
#' that is badly wrong at small df (and discontinuous there), and with tens of
#' thousands of df it is already wrong a little below that. Up to |ncp| = 30
#' it agrees with the integral below to about 1e-11 at every df, and it is
#' much faster, so it is used there.
#'
#' Beyond that, P(T <= t) for T = (Z + ncp) / sqrt(V / df) is integrated over
#' the standard normal Z, conditioning on which leaves a chi-square
#' probability for V: for t < 0 the event needs Z + ncp < 0 and
#' V <= df ((Z + ncp) / t)^2; for t > 0 it holds whenever Z + ncp <= 0, and
#' otherwise needs V >= df ((Z + ncp) / t)^2. The integral is split where the
#' integrand changes fastest, so \code{integrate()} sees every feature.
#' @noRd
.pnt <- function(t, df, ncp) {
  if (abs(ncp) <= 30) return(stats::pt(t, df, ncp))
  if (t == 0) return(stats::pnorm(-ncp))
  # dnorm() is zero to double precision beyond +-38.5.
  if (t < 0) {
    lo <- -38.5; hi <- -ncp; base <- 0; lower <- TRUE
  } else {
    lo <- -ncp; hi <- 38.5; base <- stats::pnorm(-ncp); lower <- FALSE
  }
  if (lo >= hi) return(base)
  g <- function(z) {
    stats::pchisq(df * ((z + ncp) / t)^2, df, lower.tail = lower) *
      stats::dnorm(z)
  }
  # The chi-square probability changes fastest around z = t - ncp, over a
  # width of about |t| sqrt(2 / df); the normal density over a width of 1.
  h <- abs(t) * sqrt(2 / df)
  br <- sort(unique(c(lo, hi, (t - ncp) + c(-6, -2, 0, 2, 6) * h,
                      -8, -4, -2, 0, 2, 4, 8)))
  br <- br[br >= lo & br <= hi]
  pieces <- vapply(seq_len(length(br) - 1L), function(k) {
    stats::integrate(g, br[k], br[k + 1L], rel.tol = 1e-10, abs.tol = 1e-15,
                     subdivisions = 1000L)$value
  }, numeric(1))
  min(1, base + sum(pieces))
}

#' Hedges' small-sample bias correction for a standardiser on \code{df}
#' degrees of freedom, in its exact gamma-function form
#' @noRd
.hedges_j <- function(df) {
  exp(lgamma(df / 2) - log(sqrt(df / 2)) - lgamma((df - 1) / 2))
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
