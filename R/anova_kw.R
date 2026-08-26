#' Kruskal-Wallis Test with Dunn Post-hoc Comparisons
#'
#' Compares the distribution of a numeric response across groups without
#' assuming normality, using \code{\link[stats]{kruskal.test}}, followed by
#' Dunn's pairwise rank-sum comparisons with a multiplicity adjustment. An
#' ordered factor must be converted first, with \code{as.integer()}.
#'
#' @details
#' Dunn's test is computed directly from the rank sums, with the standard tie
#' correction, so the package does not depend on an external implementation and
#' the full range of \code{\link[stats]{p.adjust}} methods is available.
#'
#' Normality is not required and is not tested. When \code{diagnostics = TRUE} a
#' Q-Q plot of within-group deviations is produced for context only; it plays no
#' part in the inference.
#'
#' When several grouping variables are supplied they are combined into a single
#' cell factor. The result is then one omnibus test across all populated cells:
#' it cannot separate a main effect of one factor from a main effect of another,
#' and it cannot test an interaction. A note to that effect is added to the
#' returned object. For a factorial rank-based analysis, consider the
#' Scheirer-Ray-Hare extension or an aligned rank transform.
#'
#' @param data A data frame, or anything inheriting from one, such as a
#'   \code{data.table} or a tibble.
#' @param response Character. Name of the numeric response column.
#' @param groups Character vector. One or more grouping columns.
#' @param conf_level Numeric in (0, 1). Used for the group summary intervals.
#'   Default \code{0.95}.
#' @param adjust Character. Multiplicity adjustment for Dunn's comparisons,
#'   passed to \code{\link[stats]{p.adjust}}. One of \code{"BH"},
#'   \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"},
#'   \code{"BY"}, \code{"fdr"} or \code{"none"}. Default \code{"BH"}. Note
#'   that \code{"tukey"} is not available here, as it is on the functions whose
#'   comparisons come from \pkg{emmeans}: Dunn's test is a rank comparison, not
#'   a linear contrast.
#' @param posthoc Logical. Compute Dunn's pairwise comparisons. There are
#'   \code{choose(k, 2)} of them, so this is worth turning off when the number
#'   of groups is large. Default \code{TRUE}.
#' @param diagnostics Logical. Compute the optional normality diagnostics into
#'   \code{$assumptions$normality}, with a Q-Q plot. They are contextual only
#'   and are skipped by default, since Kruskal-Wallis makes no distributional
#'   assumption. Default \code{FALSE}.
#' @param plots Logical. Build \pkg{ggplot2} objects. They are returned in
#'   \code{$plots}, never drawn. Default \code{TRUE}.
#' @param verbose Logical. Emit progress through \code{\link[base]{message}}.
#'   Default \code{FALSE}.
#'
#' @return An \code{\link{anovatoolbox_fit}} object. \code{$anova} holds the
#'   Kruskal-Wallis test, \code{$posthoc} Dunn's comparisons with raw and
#'   adjusted p-values, \code{$effect_sizes} epsilon squared and eta squared for
#'   the rank statistic, and \code{$emmeans} one row per group with its size,
#'   median, mean rank (the quantity Dunn's comparisons are built from), mean,
#'   standard deviation and interval. \code{$assumptions$normality} is present
#'   only when \code{diagnostics = TRUE}.
#'
#' @seealso \code{\link{anova_welch}} for the parametric equivalent.
#'
#' @references
#' Dunn, O. J. (1964). Multiple comparisons using rank sums.
#' \emph{Technometrics}, 6(3), 241-252.
#'
#' @examples
#' set.seed(123)
#' d <- data.frame(
#'   group = rep(c("A", "B", "C"), each = 25),
#'   value = c(rnorm(25, 5), rnorm(25, 7, 1.5), rnorm(25, 4, 0.8))
#' )
#' fit <- anova_kw(d, "value", "group")
#' fit
#' fit$posthoc
#'
#' # Epsilon squared and eta squared for the rank statistic, and the group
#' # summary the comparisons are built from
#' fit$effect_sizes
#' fit$emmeans
#'
#' # Any p.adjust() method may be used; the default is "BH"
#' anova_kw(d, "value", "group", adjust = "bonferroni",
#'          plots = FALSE)$posthoc$p_adjusted
#'
#' # Optional normality diagnostics, which play no part in the test itself
#' anova_kw(d, "value", "group", diagnostics = TRUE,
#'          plots = FALSE)$assumptions$normality
#'
#' # Small samples work: there is no minimum group size
#' small <- data.frame(g = rep(c("a", "b"), each = 3), y = c(1, 2, 3, 4, 5, 6))
#' anova_kw(small, "y", "g", plots = FALSE)$anova
#'
#' @export
anova_kw <- function(data, response, groups,
                     conf_level = 0.95,
                     adjust = c("BH", "holm", "hochberg", "hommel",
                                "bonferroni", "BY", "fdr", "none"),
                     diagnostics = FALSE,
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
  .check_flag(diagnostics, "diagnostics")
  .check_flag(posthoc, "posthoc")
  .check_flag(plots, "plots")
  .check_flag(verbose, "verbose")

  .say(verbose, "Preparing data.")
  prep <- .prepare_frame(data, c(response, groups), factors = groups,
                         keep_all = FALSE)
  d <- prep$data
  notes <- prep$notes
  .check_groups(d, groups, min_levels = 1L, min_n = 1L)

  added <- .add_cell(d, groups)
  d <- added$data
  cell <- added$cell
  k <- nlevels(d[[cell]])
  if (k < 2L) .stopf("At least two populated groups are required; found %d.", k)

  if (length(groups) > 1L) {
    notes <- c(notes, sprintf(
      "The %d grouping variables were combined into %d cells and tested with a single omnibus Kruskal-Wallis test. This cannot separate main effects or test an interaction; for that, use a Scheirer-Ray-Hare or aligned-rank-transform analysis.",
      length(groups), k))
  }

  ## Omnibus test ------------------------------------------------------------
  .say(verbose, "Running Kruskal-Wallis across %d groups.", k)
  kw <- stats::kruskal.test(.formula(response, cell), data = d)
  anova_tab <- data.frame(
    term      = paste(groups, collapse = " : "),
    statistic = unname(kw$statistic),
    df        = unname(kw$parameter),
    p_value   = unname(kw$p.value),
    stringsAsFactors = FALSE
  )

  ## Effect sizes ------------------------------------------------------------
  n <- nrow(d)
  H <- unname(kw$statistic)
  effect_sizes <- data.frame(
    measure  = c("epsilon_squared", "eta_squared_H"),
    estimate = c(H / ((n^2 - 1) / (n + 1)), (H - k + 1) / (n - k)),
    stringsAsFactors = FALSE
  )
  effect_sizes$estimate <- pmax(0, pmin(1, effect_sizes$estimate))

  ## Post-hoc ----------------------------------------------------------------
  posthoc_tab <- if (posthoc) {
    .dunn_test(d[[response]], d[[cell]], adjust = adjust)
  } else {
    notes <- c(notes, .no_posthoc_note(
      k,
      extra = "$effect_sizes still reports epsilon squared and eta squared, which come from the omnibus test."))
    NULL
  }

  ## Group summaries ---------------------------------------------------------
  summary_stats <- .summary_stats(d[[response]], d[[cell]], conf_level)
  summary_stats$median <- as.numeric(tapply(d[[response]], d[[cell]],
                                            stats::median, na.rm = TRUE)[
                                              as.character(summary_stats$group)])
  summary_stats$mean_rank <- as.numeric(
    tapply(rank(d[[response]]), d[[cell]], mean)[as.character(summary_stats$group)])
  summary_stats <- summary_stats[, c("group", "n", "median", "mean_rank",
                                     "mean", "sd", "se", "conf_low", "conf_high")]

  ## Optional diagnostics ----------------------------------------------------
  assumptions <- list()
  if (diagnostics) {
    .say(verbose, "Computing optional normality diagnostics.")
    assumptions$normality <- .normality_by_group(d[[response]], d[[cell]])
    notes <- c(notes, "The normality diagnostics are contextual only. Kruskal-Wallis does not assume normality and the test does not depend on them.")
  }

  ## Plots -------------------------------------------------------------------
  plot_list <- list()
  if (plots) {
    .say(verbose, "Building plots.")
    plot_list$box <- .plot_box(d, response, cell,
                               xlab = paste(groups, collapse = " : "))
    if (diagnostics) {
      resid_within <- unlist(lapply(levels(d[[cell]]), function(g) {
        v <- d[[response]][d[[cell]] == g]
        v - stats::median(v, na.rm = TRUE)
      }), use.names = FALSE)
      plot_list$qq <- .plot_qq(
        resid_within,
        title = "Normal Q-Q plot of within-group deviations",
        subtitle = "Diagnostic context only; Kruskal-Wallis assumes no distribution")
    }
  }

  .new_fit(
    method       = "Kruskal-Wallis rank sum test",
    call         = cl,
    model        = kw,
    anova        = anova_tab,
    effect_sizes = effect_sizes,
    emmeans      = summary_stats,
    posthoc      = posthoc_tab,
    assumptions  = assumptions,
    plots        = plot_list,
    data_used    = d,
    n_removed    = prep$n_removed,
    conf_level   = conf_level,
    notes        = notes
  )
}

#' Dunn's test of pairwise rank sums, with tie correction
#'
#' Computed directly rather than via an external package, so that the console
#' output of that package cannot leak into the caller's session and so that any
#' \code{stats::p.adjust} method may be used.
#'
#' @param values numeric response.
#' @param cell factor of group memberships.
#' @param adjust a \code{stats::p.adjust} method.
#' @noRd
.dunn_test <- function(values, cell, adjust = "BH") {
  cell <- droplevels(as.factor(cell))
  ok <- !is.na(values) & !is.na(cell)
  values <- values[ok]; cell <- droplevels(cell[ok])
  N <- length(values)
  r <- rank(values)
  ns <- as.integer(table(cell))
  names(ns) <- levels(cell)
  rbar <- tapply(r, cell, mean)

  # Tie correction: sum over tied groups of (t^3 - t), divided by 12(N - 1)
  tie_sizes <- as.integer(table(values))
  tie_term <- sum(tie_sizes^3 - tie_sizes) / (12 * (N - 1))
  sigma_base <- (N * (N + 1) / 12) - tie_term

  combos <- utils::combn(levels(cell), 2L, simplify = FALSE)
  rows <- lapply(combos, function(pair) {
    i <- pair[1L]; j <- pair[2L]
    se <- sqrt(sigma_base * (1 / ns[[i]] + 1 / ns[[j]]))
    z <- if (is.finite(se) && se > 0) (rbar[[i]] - rbar[[j]]) / se else NA_real_
    data.frame(group1 = i, group2 = j,
               mean_rank_diff = rbar[[i]] - rbar[[j]],
               z = z,
               p_value = if (is.na(z)) NA_real_ else 2 * stats::pnorm(-abs(z)),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$p_adjusted <- stats::p.adjust(out$p_value, method = adjust)
  out$adjustment <- adjust
  row.names(out) <- NULL
  out
}
