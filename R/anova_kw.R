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
#' \strong{Ties.} Values are ranked exactly as stored, and ties are counted on
#' those ranks, both for the omnibus test and for Dunn's comparisons, so the
#' tie correction always matches the ranks used. Two values that differ only
#' through floating-point error (\code{0.1 + 0.2} and \code{0.3}, say) are
#' therefore distinct, not tied; \code{$notes} says so when the data contain
#' such values. Round the response first if they are meant to be equal.
#' (\code{\link[stats]{kruskal.test}} itself counts ties on the printed values,
#' so in that one situation its statistic differs slightly from the one
#' reported here, which equals \code{kruskal.test()} applied to the ranks.)
#'
#' \strong{Infinite values} of the response are kept, not dropped: they rank
#' as the most extreme observations, which is how a rank analysis treats, for
#' example, non-completers coded as \code{Inf}. A group containing one has an
#' infinite (or, with both signs, undefined) mean and no standard deviation, so
#' those columns of \code{$emmeans} are reported as \code{Inf}, \code{-Inf} or
#' \code{NA}. Its mean rank is exact, as is its median unless that falls
#' between \code{-Inf} and \code{Inf} (then \code{NA}). The box plot cannot
#' draw infinite values. Missing values are dropped and counted as usual.
#'
#' \strong{Effect sizes.} With \eqn{H} the tie-corrected statistic, \eqn{n}
#' observations and \eqn{k} groups, \code{epsilon_squared} is
#' \eqn{H / (n - 1)} and \code{eta_squared_H} is
#' \eqn{(H - k + 1) / (n - k)}, the names used by Tomczak and Tomczak (2014)
#' and \pkg{effectsize}. The first is exactly the proportion of rank variance
#' between groups (the \eqn{R^2} of a one-way analysis of variance on the
#' ranks), and the second is its bias-adjusted counterpart, so
#' \code{epsilon_squared} is never smaller than \code{eta_squared_H}: the
#' reverse of what the same names mean for a parametric analysis. A negative
#' \code{eta_squared_H} (less between-group variation than chance produces) is
#' reported as 0, with a note giving the unfloored value; with one observation
#' per group (\eqn{n = k}) it is undefined and reported as \code{NA}, also
#' with a note.
#'
#' Normality is not required and is not tested. When \code{diagnostics = TRUE}
#' a Q-Q plot of within-group standardised deviations (each value minus its
#' group mean, divided by its group standard deviation) is produced for
#' context only; it plays no part in the inference.
#'
#' When several grouping variables are supplied they are combined into a single
#' cell factor. The result is then one omnibus test across all populated cells,
#' labelled \code{"A x B cells"} in \code{$anova}: it cannot separate a main
#' effect of one factor from a main effect of another, and it cannot test an
#' interaction. A note to that effect is added to the returned object. For a
#' factorial rank-based analysis, consider the Scheirer-Ray-Hare extension or
#' an aligned rank transform.
#'
#' A response in which every value is the same is refused: there is nothing to
#' rank.
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
#'   \code{"BY"}, \code{"fdr"} or \code{"none"}, spelled out in full. Default
#'   \code{"BH"}. Note that \code{"tukey"} is not available here, as it is on
#'   the functions whose comparisons come from \pkg{emmeans}: Dunn's test is a
#'   rank comparison, not a linear contrast.
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
#' @return An \code{\link{anovakit_fit}} object.
#'   \itemize{
#'     \item \code{$anova}: the Kruskal-Wallis test, one row. Its \code{term}
#'       is the grouping column, or \code{"A x B cells"} when several were
#'       combined.
#'     \item \code{$posthoc}: Dunn's comparisons, one row per pair of groups,
#'       with \code{group1}, \code{group2}, \code{mean_rank_diff} (the mean
#'       rank of \code{group1} minus that of \code{group2}), \code{z} (the
#'       same difference over its standard error, so with the same sign),
#'       \code{p_value} (two-sided, unadjusted), \code{p_adjusted} (adjusted
#'       by \code{adjust}) and \code{adjustment}: the same p-value columns as
#'       every other function in the package. \code{rstatix::dunn_test()}
#'       reports \code{group2} minus \code{group1}, so its statistic has the
#'       opposite sign.
#'     \item \code{$effect_sizes}: \code{epsilon_squared} and
#'       \code{eta_squared_H}, defined in Details.
#'     \item \code{$emmeans}: one row per group with its size, median, mean
#'       rank (the quantity Dunn's comparisons are built from), mean, standard
#'       deviation and t interval for the mean. These are raw group summaries;
#'       no model is fitted.
#'     \item \code{$assumptions$normality}: present only when
#'       \code{diagnostics = TRUE}.
#'     \item \code{$model}: the \code{htest} from
#'       \code{\link[stats]{kruskal.test}}, computed on the ranks; its
#'       \code{data.name} names the response and the grouping column(s).
#'   }
#'
#' @seealso \code{\link{anova_welch}} for the parametric equivalent.
#'
#' @references
#' Cohen, B. H. (2008). \emph{Explaining Psychological Statistics} (3rd ed.).
#' Wiley.
#'
#' Dunn, O. J. (1964). Multiple comparisons using rank sums.
#' \emph{Technometrics}, 6(3), 241-252.
#'
#' Tomczak, M., & Tomczak, E. (2014). The need to report effect size
#' estimates revisited. An overview of some recommended measures of effect
#' size. \emph{Trends in Sport Sciences}, 21(1), 19-25.
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
#' # Infinite values are kept and ranked as the most extreme outcomes
#' worst <- data.frame(g = rep(c("a", "b"), each = 5),
#'                     y = c(1:5, 6, 7, Inf, Inf, Inf))
#' anova_kw(worst, "y", "g", plots = FALSE)$anova
#'
#' @export
anova_kw <- function(data, response, groups,
                     conf_level = 0.95,
                     adjust = "BH",
                     diagnostics = FALSE,
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
  .check_flag(diagnostics, "diagnostics")
  .check_flag(posthoc, "posthoc")
  .check_flag(plots, "plots")
  .check_flag(verbose, "verbose")

  .say(verbose, "Preparing data.")
  # Ranks can use infinite values, so they are kept in the response.
  prep <- .prepare_frame(data, c(response, groups), factors = groups,
                         keep_all = FALSE, keep_infinite = response)
  d <- prep$data
  notes <- prep$notes
  .check_groups(d, groups, min_levels = 1L, min_n = 1L)

  added <- .add_cell(d, groups)
  d <- added$data
  cell <- added$cell
  k <- nlevels(d[[cell]])
  if (k < 2L) .stopf("At least two populated groups are required; found %d.", k)

  y <- as.numeric(d[[response]])
  if (length(unique(y)) == 1L) {
    .stopf("The Kruskal-Wallis test needs variation in the response, but every one of the %d observations of `%s` takes the same value (%s). There is nothing to rank.",
           length(y), response, format(y[1L]))
  }
  n_infinite <- sum(is.infinite(y))
  if (n_infinite > 0L) {
    notes <- c(notes, sprintf(
      "%d infinite value(s) of `%s` were kept and ranked as the most extreme observations. Group means and standard deviations that involve them are reported as Inf, -Inf or NA; the mean ranks are exact. The box plot does not draw them.",
      n_infinite, response))
  }
  if (.has_near_ties(y)) {
    notes <- c(notes, sprintf(
      "Some values of `%s` differ only in their last few significant digits (as 0.1 + 0.2 and 0.3 do) and were ranked as distinct values, not as ties. If they are meant to be equal, round the response before the analysis.",
      response))
  }

  term_label <- .cells_label(groups)
  if (length(groups) > 1L) {
    notes <- c(notes, sprintf(
      "The %d grouping variables were combined into %d cells and tested with a single omnibus Kruskal-Wallis test. This cannot separate main effects or test an interaction; for that, use a Scheirer-Ray-Hare or aligned-rank-transform analysis.",
      length(groups), k))
  }

  ## Omnibus test ------------------------------------------------------------
  .say(verbose, "Running Kruskal-Wallis across %d groups.", k)
  # kruskal.test() counts ties on the printed values (table(x)) but ranks the
  # exact ones. Given the ranks, it counts ties on those, which is always
  # consistent with the ranking, and is what Dunn's comparisons use too.
  r <- rank(y)
  kw <- stats::kruskal.test(r, d[[cell]])
  kw$data.name <- sprintf("%s by %s", response, term_label)
  anova_tab <- data.frame(
    term      = term_label,
    statistic = unname(kw$statistic),
    df        = unname(kw$parameter),
    p_value   = unname(kw$p.value),
    stringsAsFactors = FALSE
  )

  ## Effect sizes ------------------------------------------------------------
  n <- length(y)
  H <- unname(kw$statistic)
  epsilon_sq <- min(1, max(0, H / (n - 1)))
  eta_sq_h <- if (n > k) (H - k + 1) / (n - k) else NA_real_
  if (is.na(eta_sq_h)) {
    notes <- c(notes, "eta_squared_H is undefined with one observation per group (its denominator, n - k, is zero) and is reported as NA.")
  } else if (eta_sq_h < 0) {
    notes <- c(notes, sprintf(
      "eta_squared_H was negative (%s) and has been floored at 0, which happens when the groups differ less than their degrees of freedom would by chance.",
      format(signif(eta_sq_h, 3))))
    eta_sq_h <- 0
  } else if (eta_sq_h > 1) {
    eta_sq_h <- 1
  }
  effect_sizes <- data.frame(
    measure  = c("epsilon_squared", "eta_squared_H"),
    estimate = c(epsilon_sq, eta_sq_h),
    stringsAsFactors = FALSE
  )

  ## Post-hoc ----------------------------------------------------------------
  posthoc_tab <- if (posthoc) {
    .dunn_test(y, d[[cell]], adjust = adjust)
  } else {
    notes <- c(notes, .no_posthoc_note(
      k,
      extra = "$effect_sizes still reports epsilon squared and eta squared, which come from the omnibus test."))
    NULL
  }

  ## Group summaries ---------------------------------------------------------
  # tapply() over the cell factor returns one value per level, in level order,
  # which is the order of .summary_stats(); indexing by name would lose a
  # level called "".
  summary_stats <- .summary_stats(y, d[[cell]], conf_level)
  summary_stats$median <- as.numeric(tapply(y, d[[cell]], stats::median))
  summary_stats$mean_rank <- as.numeric(tapply(r, d[[cell]], mean))
  summary_stats <- summary_stats[, c("group", "n", "median", "mean_rank",
                                     "mean", "sd", "se", "conf_low", "conf_high")]
  if (n_infinite > 0L) {
    # Inf - Inf is NaN: an undefined quantity, reported as NA.
    for (col in c("median", "mean", "sd", "se", "conf_low", "conf_high")) {
      v <- summary_stats[[col]]
      v[is.nan(v)] <- NA_real_
      summary_stats[[col]] <- v
    }
    has_inf <- as.logical(tapply(is.infinite(y), d[[cell]], any))
    summary_stats[has_inf, c("sd", "se", "conf_low", "conf_high")] <- NA_real_
  }

  ## Optional diagnostics ----------------------------------------------------
  assumptions <- list()
  if (diagnostics) {
    .say(verbose, "Computing optional normality diagnostics.")
    assumptions$normality <- .normality_by_group(y, d[[cell]])
    notes <- c(notes, "The normality diagnostics are contextual only. Kruskal-Wallis does not assume normality and the test does not depend on them.",
               .normality_skip_notes(assumptions$normality))
  }

  ## Plots -------------------------------------------------------------------
  plot_list <- list()
  if (plots) {
    .say(verbose, "Building plots.")
    plot_list$box <- .plot_box(d, response, cell,
                               xlab = paste(groups, collapse = " : "))
    if (diagnostics) {
      plot_list$qq <- .plot_qq(
        .standardise_within(y, d[[cell]]),
        title = "Normal Q-Q plot of within-group standardised deviations",
        subtitle = "(value - group mean) / group SD. Context only; Kruskal-Wallis assumes no distribution")
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

#' Are there distinct finite values that differ only by rounding error?
#'
#' Adjacent distinct values within 64 units in the last place of each other
#' are almost always one value computed two ways (0.1 + 0.2 and 0.3 are one
#' unit apart), not two measurements: continuous data would need millions of
#' observations to produce such a gap by chance.
#' @noRd
.has_near_ties <- function(x) {
  u <- sort(unique(x[is.finite(x)]))
  if (length(u) < 2L) return(FALSE)
  gap <- diff(u)
  scale <- pmax(abs(u[-1L]), abs(u[-length(u)]))
  any(gap <= 64 * .Machine$double.eps * scale)
}

#' Dunn's test of pairwise rank sums, with tie correction
#'
#' Computed directly rather than via an external package, so that the console
#' output of that package cannot leak into the caller's session and so that any
#' \code{stats::p.adjust} method may be used. Groups are handled by position,
#' so a level may be any string, including the empty one. Ties are counted on
#' the ranks, so the correction always matches the ranking.
#'
#' @param values numeric response.
#' @param cell factor of group memberships.
#' @param adjust a \code{stats::p.adjust} method.
#' @return data.frame with \code{group1}, \code{group2}, \code{mean_rank_diff}
#'   and \code{z} (both \code{group1} minus \code{group2}), \code{p_value},
#'   \code{p_adjusted} and \code{adjustment}.
#' @noRd
.dunn_test <- function(values, cell, adjust = "BH") {
  cell <- droplevels(as.factor(cell))
  ok <- !is.na(values) & !is.na(cell)
  values <- values[ok]; cell <- droplevels(cell[ok])
  lv <- levels(cell)
  N <- length(values)
  r <- rank(values)
  ns <- as.vector(table(cell))
  rbar <- as.vector(tapply(r, cell, mean))

  # Tie correction: sum over tied groups of (t^3 - t), divided by 12(N - 1).
  # Tied values share a rank, and distinct values never do.
  tie_sizes <- tabulate(match(r, unique(r)))
  tie_term <- sum(tie_sizes^3 - tie_sizes) / (12 * (N - 1))
  sigma_base <- (N * (N + 1) / 12) - tie_term

  idx <- utils::combn(length(lv), 2L)
  i <- idx[1L, ]; j <- idx[2L, ]
  diff_rank <- rbar[i] - rbar[j]
  se <- sqrt(sigma_base * (1 / ns[i] + 1 / ns[j]))
  z <- ifelse(is.finite(se) & se > 0, diff_rank / se, NA_real_)
  out <- data.frame(group1 = lv[i], group2 = lv[j],
                    mean_rank_diff = diff_rank,
                    z = z,
                    p_value = 2 * stats::pnorm(-abs(z)),
                    stringsAsFactors = FALSE)
  out$p_adjusted <- stats::p.adjust(out$p_value, method = adjust)
  out$adjustment <- adjust
  row.names(out) <- NULL
  out
}
