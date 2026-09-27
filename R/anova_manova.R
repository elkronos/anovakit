#' Multivariate Analysis of Variance and Covariance
#'
#' Compares two or more numeric responses jointly across groups, optionally
#' adjusting for covariates. Reports the multivariate test, follow-up univariate
#' analyses with effect sizes, covariate-adjusted estimated marginal means,
#' assumption checks (Mardia's multivariate normality tests, Box's M and, with
#' covariates, a multivariate test of homogeneous regression slopes), and a
#' canonical discriminant analysis with its plot.
#'
#' @details
#' \strong{The multivariate tests} are Type II or Type III tests, as
#' \code{type} says, computed by \code{car::Manova()} on the multivariate
#' linear model \code{lm(cbind(<responses>) ~ <covariates> + <groups>)}. Type
#' II tests each term adjusted for every other term that does not contain it,
#' so the result does not depend on the order of \code{groups} or
#' \code{covariates}; Type III fits under sum-to-zero contrasts and tests each
#' term adjusted for all others. (\code{summary(stats::manova())} would give
#' sequential, order-dependent tests instead.) When the model has aliased
#' coefficients -- an empty cell of a crossed design, or collinear covariates
#' -- \code{car::Manova()} refuses the model, so the Type II tests are
#' computed by the equivalent model comparisons, and a note says so. Type III
#' tests are not defined when coefficients are aliased; Type II tests are then
#' reported throughout, with a note.
#'
#' With \code{test = "Roy"} the F statistic for a term with more than one
#' degree of freedom (and more than one response) is an upper bound, so its
#' p-value is a lower bound: it is anti-conservative, and a note says so.
#'
#' \strong{Covariates} are mean-centred before fitting (the means are in
#' \code{$covariate_means} and in a note). Centring changes no test; it keeps a
#' covariate with a large offset and a small spread (a time stamp, say) from
#' being mistaken for a constant by the fitting routine, and it means the
#' adjusted marginal means are read at the covariate means. A MANCOVA assumes
#' that each covariate has the same slope in every group. That assumption is
#' tested by comparing the model with one in which each covariate interacts
#' with every grouping term, using the chosen multivariate statistic; the
#' result is in \code{$slopes_test}, and a note says when it is rejected at
#' the 0.05 level, since the common-slope adjusted means and comparisons are
#' then not interpretable as constant group differences.
#'
#' \strong{Assumption checks.} Mardia's tests of multivariate skewness and
#' kurtosis (on the residuals of the multivariate model) and Box's M test of
#' equality of covariance matrices (across the cells of the grouping
#' variables, with the covariate effects removed) are computed from first
#' principles, so no additional package is required. Mardia's tests are
#' skipped, with a note, unless the residual degrees of freedom exceed the
#' number of responses by at least 10: when they equal it the statistics do not
#' depend on the data at all, and close to it they are dominated by the
#' design. Both are asymptotic tests, and Mardia's kurtosis test in particular
#' rejects too often when the sample is small relative to the square of the
#' number of responses. Box's M is notoriously sensitive to non-normality:
#' treat a small p-value as a prompt to look at the group covariances rather
#' than as a verdict.
#'
#' \strong{Univariate follow-ups.} One Type II or III ANOVA is fitted per
#' response and kept under \code{$univariate}. Their p-values are \emph{not}
#' adjusted across responses and are reported whether or not the multivariate
#' test is significant; \code{adjust} applies only to the pairwise comparisons
#' within each response. A note says when a term's multivariate test is not
#' significant at the 0.05 level, since its univariate follow-ups are then not
#' protected by it.
#'
#' \strong{Canonical discriminant analysis.} The hypothesis and error sum of
#' squares and cross-products matrices of one term are taken from the same
#' \code{car::Manova()} tests as \code{$multivariate}, and the generalised
#' eigenproblem is solved directly (on responses scaled to unit error
#' variance, so responses on very different scales do not make it singular).
#' \code{$canonical} reports the eigenvalues, canonical correlations and the
#' proportion of the term's between-group variation on each axis;
#' \code{$plots$canonical} plots the first two axes with group centroids, and
#' \code{$assumptions$structure_coefficients} gives the pooled within-group
#' correlations between each response and each axis, computed from the error
#' matrix. The analysis describes \emph{one} term: the full interaction of the
#' grouping variables when it is in the model (\code{interaction = TRUE}, or an
#' order equal to the number of grouping variables); otherwise the main effect
#' of the \emph{first} grouping variable, \code{groups[1]} -- list another
#' variable first to describe it instead. \code{$canonical_term} names the
#' term, and a note says which term was used and why. The scores are computed
#' after removing the fitted effects of every other term in the model (the
#' covariates, and the other grouping terms, under sum-to-zero coding), so
#' that the plot shows the described term alone. Axes are scaled to unit
#' pooled within-group variance (the error matrix divided by its degrees of
#' freedom), so the spread of the centroids on the plot means what it appears
#' to mean; in a balanced design the between-group to within-group ratio of the
#' plotted scores equals the eigenvalue exactly, and in an unbalanced one
#' approximately.
#'
#' \strong{Estimated marginal means} come from \pkg{emmeans} applied to each
#' univariate model, so with covariates present they are adjusted means, not
#' raw group means. When several grouping variables enter additively they are
#' reported for each factor separately, averaged over the others, with a
#' \code{term} column.
#'
#' Responses and covariates are named as columns rather than composed into a
#' formula. If you want a transformed response, create the column first.
#'
#' @param data A data frame, or anything inheriting from one, such as a
#'   \code{data.table} or a tibble.
#' @param responses Character vector of at least two numeric response columns.
#'   A single response is accepted and falls back to a univariate analysis,
#'   with a note.
#' @param groups Character vector. One or more grouping columns.
#' @param covariates Character vector or \code{NULL}. Numeric covariates, which
#'   turn the analysis into a MANCOVA. They are mean-centred before fitting.
#' @param test Character. Multivariate statistic to report: \code{"Pillai"}
#'   (default), \code{"Wilks"}, \code{"Hotelling-Lawley"} or \code{"Roy"}.
#'   Pillai's trace is the most robust to departures from the assumptions.
#'   Roy's largest root has only an upper bound for its F statistic when the
#'   hypothesis has more than one degree of freedom, so its p-value is then a
#'   lower bound (anti-conservative).
#' @param interaction \code{FALSE} (additive, the default), \code{TRUE} (full
#'   factorial across the grouping variables), or a whole number giving the
#'   highest interaction order.
#' @param type Character. \code{"II"} (default) or \code{"III"}: the type of
#'   both the multivariate tests and the univariate follow-up tests. With
#'   \code{"III"} the multivariate model and the univariate ones are fitted
#'   under sum-to-zero contrasts.
#' @param conf_level Numeric in (0, 1). Level for every interval returned.
#'   Default \code{0.95}.
#' @param adjust Character. Multiplicity adjustment for the pairwise
#'   comparisons within each response. It does not adjust across responses.
#'   Default \code{"tukey"}.
#' @param assumptions Logical. Compute Mardia's tests and Box's M. Default
#'   \code{TRUE}. The homogeneity-of-slopes test is computed whenever there are
#'   covariates.
#' @param posthoc Logical. Compute pairwise comparisons within each response.
#'   Default \code{TRUE}.
#' @param plots Logical. Build \pkg{ggplot2} objects. They are returned in
#'   \code{$plots}, never drawn. Default \code{TRUE}.
#' @param verbose Logical. Emit progress through \code{\link[base]{message}}.
#'   Default \code{FALSE}.
#'
#' @return An \code{\link{anovakit_fit}} object. Its \code{$model} is the
#'   multivariate linear model (class \code{mlm}; pass it to
#'   \code{car::Manova()} for further tests), and its \code{$anova} is the
#'   multivariate table. Extra components: \code{$univariate} (a named list of
#'   per-response fits, each holding that response's own \code{$model},
#'   \code{$anova}, \code{$effect_sizes}, \code{$emmeans},
#'   \code{$emmeans_object} and \code{$posthoc}; the univariate p-values are
#'   not adjusted across responses), \code{$multivariate} (the Type II or III
#'   multivariate table, one row per term, with columns \code{term}, \code{df},
#'   \code{statistic}, \code{approx_f}, \code{num_df}, \code{den_df} and
#'   \code{p_value}), \code{$canonical} (the discriminant axes),
#'   \code{$canonical_term} (the term they describe, spelled as in
#'   \code{$multivariate$term}), \code{$slopes_test} (with covariates, the
#'   multivariate test of homogeneous regression slopes; otherwise
#'   \code{NULL}), \code{$covariate_means} (the means the covariates were
#'   centred at) and \code{$test} (the multivariate statistic used).
#'   \code{$assumptions$structure_coefficients} has a \code{response} column
#'   and one column per axis. The top-level \code{$emmeans} and
#'   \code{$posthoc} stack every response's table with a \code{response}
#'   column (\code{response.1} if a grouping column is already called
#'   \code{response}); \code{$emmeans_object} is \code{NULL}, because there is
#'   one grid per response -- take them from
#'   \code{$univariate[[r]]$emmeans_object} (a named list of grids, one per
#'   factor, when the factors enter additively).
#'
#' @seealso \code{\link{anova_ancova}} for one response with covariates,
#'   \code{\link{anova_glm}} for one response in any family.
#'
#' @references
#' Mardia, K. V. (1970). Measures of multivariate skewness and kurtosis with
#' applications. \emph{Biometrika}, 57(3), 519-530.
#'
#' Box, G. E. P. (1949). A general distribution theory for a class of likelihood
#' criteria. \emph{Biometrika}, 36(3/4), 317-346.
#'
#' @examples
#' set.seed(1)
#' n <- 120
#' d <- data.frame(g = factor(rep(c("a", "b", "c"), each = n / 3)),
#'                 age = rnorm(n, 40, 8))
#' d$score1 <- 5 + 2 * (d$g == "b") + 4 * (d$g == "c") + 0.1 * d$age + rnorm(n)
#' d$score2 <- 3 + 1 * (d$g == "b") + 2 * (d$g == "c") + 0.05 * d$age + rnorm(n)
#'
#' fit <- anova_manova(d, c("score1", "score2"), "g")
#' fit
#'
#' # Assumption checks, computed directly rather than through a dependency
#' fit$assumptions$mardia
#' fit$assumptions$box_m
#'
#' # The discriminant axes, the term they describe, and how each response
#' # loads on them
#' fit$canonical
#' fit$canonical_term
#' fit$assumptions$structure_coefficients
#'
#' # Follow-up tests, one per response, are kept separately as well as stacked
#' fit$univariate$score1$anova
#'
#' # As a MANCOVA: the marginal means are adjusted for age, and the
#' # common-slope assumption is tested
#' mc <- anova_manova(d, c("score1", "score2"), "g",
#'                    covariates = "age", plots = FALSE)
#' mc$emmeans
#' mc$slopes_test
#'
#' @export
anova_manova <- function(data, responses, groups, covariates = NULL,
                         test = c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"),
                         interaction = FALSE,
                         type = c("II", "III"),
                         conf_level = 0.95,
                         adjust = "tukey",
                         assumptions = TRUE,
                         posthoc = TRUE,
                         plots = TRUE,
                         verbose = FALSE) {
  cl <- match.call()
  test <- match.arg(test)
  type <- match.arg(type)
  .check_adjust(adjust, .adjust_choices(), "adjust")
  data <- .as_df(data)
  .check_names(responses, "responses")
  .check_names(groups, "groups")
  if (!is.null(covariates)) .check_names(covariates, "covariates")
  .check_columns(data, responses, "Response column(s)")
  .check_columns(data, groups, "Grouping column(s)")
  if (!is.null(covariates)) .check_columns(data, covariates, "Covariate column(s)")
  .check_roles(list(`a response` = responses, `a grouping variable` = groups,
                    `a covariate` = covariates))
  for (r in responses) .check_numeric_col(data, r)
  if (!is.null(covariates)) {
    for (cv in covariates) .check_numeric_col(data, cv, role = "Covariate")
  }
  .check_conf_level(conf_level)
  .check_interaction(interaction, length(groups))
  .check_flag(assumptions, "assumptions")
  .check_flag(posthoc, "posthoc")
  .check_flag(plots, "plots")
  .check_flag(verbose, "verbose")

  .say(verbose, "Preparing data.")
  cols <- c(responses, groups, covariates)
  prep <- .prepare_frame(data, cols, factors = groups, keep_all = FALSE)
  d <- prep$data
  notes <- prep$notes
  .check_groups(d, groups, min_levels = 2L, min_n = 1L)

  multivariate <- length(responses) >= 2L
  if (!multivariate) {
    notes <- c(notes, "Only one response was supplied, so this is a univariate analysis; there is no multivariate test and no canonical discriminant analysis.")
  }

  ## Centre the covariates ----------------------------------------------------
  # Centring leaves every test and every adjusted mean unchanged (no covariate
  # enters an interaction), but lm()'s pivoting QR declares a covariate whose
  # coefficient of variation is below about 1e-7 -- a time stamp spanning a few
  # minutes -- aliased with the intercept, and silently drops it.
  covariate_means <- NULL
  if (!is.null(covariates)) {
    covariate_means <- vapply(covariates, function(cv) mean(d[[cv]]), numeric(1))
    for (cv in covariates) d[[cv]] <- d[[cv]] - covariate_means[[cv]]
    notes <- c(notes, sprintf(
      "Covariate(s) mean-centred before fitting (%s), so $emmeans are adjusted means at the covariate mean(s). Centring changes no test.",
      paste(sprintf("%s: %.6g", covariates, covariate_means), collapse = "; ")))
  }

  cellv <- .cell_vector(d, groups)
  g_terms <- .group_terms(groups, interaction)
  rhs <- c(.bq(covariates), g_terms)
  lhs <- if (multivariate) {
    sprintf("cbind(%s)", paste(.bq(responses), collapse = ", "))
  } else .bq(responses)

  ## The model the tests are computed from -------------------------------------
  .say(verbose, "Fitting the %s model.",
       if (multivariate) "multivariate" else "univariate")
  fml_main <- stats::reformulate(rhs, response = lhs)
  fit <- .fit_with_contrasts(function() .fit_lm(fml_main, d), type)
  # A covariate named like a contrast coefficient (`dose1` next to a factor
  # `dose` under sum-to-zero contrasts) gives the model two coefficients of
  # the same name, and emmeans, which matches by name, then silently reads the
  # wrong column for every response.
  dup <- .dup_coef_names(fit)
  if (length(dup) > 0L) {
    d <- .bracket_group_contrasts(d, groups, type)
    fit <- .fit_with_contrasts(function() .fit_lm(fml_main, d), type)
    dup_lab <- gsub("`", "", dup, fixed = TRUE)
    involved <- unique(c(
      covariates[covariates %in% dup_lab],
      groups[vapply(groups, function(g) any(startsWith(dup_lab, g)), logical(1))]))
    if (length(.dup_coef_names(fit)) > 0L) {
      .stopf("Column names collide with the names of the model's coefficients (%s is produced twice), so the marginal means would be computed from the wrong columns. Rename the column(s) involved: %s.",
             paste(dup_lab, collapse = ", "),
             paste(sprintf("`%s`", involved), collapse = ", "))
    }
    relabelled <- grep("[", rownames(as.matrix(stats::coef(fit))), fixed = TRUE,
                       value = TRUE)
    notes <- c(notes, sprintf(
      "The names of columns %s collide with the coefficient names of a grouping factor's contrasts (%s would appear twice), which would make emmeans read the wrong model column. The grouping factors' contrast labels were put in brackets in the models (%s, ...) to keep the names distinct; no estimate or test is affected.",
      paste(sprintf("`%s`", involved), collapse = ", "),
      paste(dup_lab, collapse = ", "),
      paste(gsub("`", "", utils::head(relabelled, 2L), fixed = TRUE),
            collapse = ", ")))
  }
  aliased <- anyNA(stats::coef(fit))
  type_used <- type
  if (aliased && identical(type, "III")) {
    type_used <- "II"
    notes <- c(notes, "The model has aliased coefficients (an empty cell of the design, or collinear covariates), and Type III tests are not defined for it: car refuses them. Type II tests were computed instead, for the multivariate table and the univariate follow-ups alike.")
  }

  ## Multivariate tests ---------------------------------------------------------
  mv_table <- NULL
  ssp <- NULL
  if (multivariate) {
    ssp <- .manova_ssp(fit, type_used)
    notes <- c(notes, ssp$note)
    Y <- as.matrix(d[, responses, drop = FALSE])
    tss <- colSums(sweep(Y, 2L, colMeans(Y))^2)
    mv_table <- tryCatch(.manova_table(ssp, test, tss), error = function(e) e)
    if (inherits(mv_table, "error")) {
      .stopf("The multivariate test could not be computed: %s. This usually means the responses are collinear, or there are too few residual degrees of freedom for %d responses.",
             conditionMessage(mv_table), length(responses))
    }
    p <- length(responses)
    roy_bound <- mv_table$term[is.finite(mv_table$df) & pmin(mv_table$df, p) > 1]
    if (identical(test, "Roy") && length(roy_bound) > 0L) {
      notes <- c(notes, sprintf(
        "For %s, Roy's largest root has only an upper bound for its F statistic (the hypothesis has more than one degree of freedom and there are %d responses), so the reported p-value is a lower bound: it is anti-conservative. Pillai's trace or Wilks' lambda give calibrated p-values.",
        paste(sprintf("`%s`", roy_bound), collapse = ", "), p))
    }
  }

  ## Homogeneity of regression slopes ------------------------------------------
  slopes_test <- NULL
  if (!is.null(covariates)) {
    sl <- .slopes_test(d, lhs, rhs, covariates, g_terms, test, multivariate)
    slopes_test <- sl$table
    notes <- c(notes, sl$note)
  }

  ## Univariate follow-ups ---------------------------------------------------
  .say(verbose, "Fitting %d univariate model(s).", length(responses))
  univariate <- stats::setNames(vector("list", length(responses)), responses)
  eff_rows <- list()
  emm_rows <- list()
  ph_rows <- list()
  resp_col <- "response"
  term_col <- NULL
  for (r in responses) {
    m <- if (!multivariate) fit else {
      fml_r <- stats::reformulate(rhs, response = .bq(r))
      .fit_with_contrasts(function() .fit_lm(fml_r, d), type)
    }
    av <- .car_anova(m, type = type_used)
    # Aliasing is a property of the design, so it is the same for every
    # response; report it once rather than per response.
    if (r == responses[[1L]]) notes <- c(notes, .check_model_size(m))
    notes <- c(notes, av$note)
    es <- .partial_eta_squared(av$raw, df_error = stats::df.residual(m),
                               n_obs = nrow(d))
    if (!is.null(attr(es, "omega_floored"))) {
      notes <- c(notes, sprintf(
        "Partial omega squared was negative for %s in the model for `%s` and has been floored at 0.",
        paste(attr(es, "omega_floored"), collapse = ", "), r))
    }
    if (nrow(es) > 0L) eff_rows[[r]] <- .cbind_label(es, "response", r)

    emm <- .emm_block(m, groups, additive = .is_additive(interaction),
                      type = "link", conf_level = conf_level, adjust = adjust,
                      posthoc = posthoc, link = "identity", data = d,
                      n_each = length(responses))
    notes <- c(notes, emm$notes)
    if (!is.null(emm$table)) {
      term_col <- attr(emm$table, "term_col")
      # The label column is called `response` unless a grouping column already
      # is; the plot must facet on the name actually used.
      resp_col <- .safe_name(emm$table, "response")
      emm_rows[[r]] <- .cbind_label(emm$table, "response", r)
    }
    if (!is.null(emm$posthoc)) ph_rows[[r]] <- .cbind_label(emm$posthoc, "response", r)

    univariate[[r]] <- list(model = m, anova = av$table, effect_sizes = es,
                            emmeans = emm$table, emmeans_object = emm$grid,
                            posthoc = emm$posthoc)
  }
  stack <- function(rows) {
    if (length(rows) == 0L) return(NULL)
    out <- do.call(rbind, unname(rows))
    row.names(out) <- NULL
    out
  }
  effect_sizes <- stack(eff_rows)
  emm_all <- stack(emm_rows)
  ph_all  <- stack(ph_rows)

  if (multivariate) {
    ns <- mv_table$term[is.finite(mv_table$p_value) & mv_table$p_value >= 0.05]
    if (length(ns) > 0L) {
      notes <- c(notes, sprintf(
        "The multivariate test is not significant at the 0.05 level for %s, yet $univariate reports a follow-up test of %s for each of the %d responses. Those follow-ups are not protected by the multivariate test and not adjusted across responses, so a small univariate p-value there is not a finding on its own.",
        paste(sprintf("`%s`", ns), collapse = ", "),
        if (length(ns) == 1L) "it" else "them", length(responses)))
    }
  }

  ## Responses adjusted for the other terms -------------------------------------
  # The assumption checks and the canonical scores live in the space the
  # multivariate test operates in, so the fitted contributions of the
  # covariates (and, for the canonical scores, of the other grouping terms)
  # are removed. They are taken from the full model under sum-to-zero coding,
  # where a term's contribution is its effect averaged over the others.
  Y <- as.matrix(d[, responses, drop = FALSE])
  canonical_term <- NA_character_
  canon_label <- NULL
  canon_why <- NULL
  Ycov <- Y
  Ycda <- Y
  term_labels <- attr(stats::terms(fit), "term.labels")
  if (multivariate) {
    ct <- .choose_canonical_term(fit, groups)
    canon_label <- ct$label
    canon_why <- ct$why
    canonical_term <- gsub("`", "", canon_label, fixed = TRUE)
    part <- .partial_out(fit, Y, covariates, canon_label)
    if (is.null(part)) {
      notes <- c(notes, "The fitted effects of the covariates and of the other terms could not be removed from the responses, so Box's M and the canonical scores below describe the unadjusted responses.")
    } else {
      Ycov <- part$covariates
      Ycda <- part$canonical
      if (!is.null(covariates) && assumptions) {
        notes <- c(notes, "Box's M compares the cell covariance matrices of the responses with the covariate effects removed (as estimated in the full model), and Mardia's tests use the model residuals: the conditional distribution the MANCOVA assumes.")
      }
    }
  }

  ## Assumptions -------------------------------------------------------------
  assum <- list()
  if (assumptions && multivariate) {
    .say(verbose, "Checking multivariate assumptions.")
    mt <- .mardia_test(stats::residuals(fit), df_resid = ssp$error_df)
    assum$mardia <- mt$table
    notes <- c(notes, mt$note)
    if (!is.null(mt$table) && any(mt$table$p_value < 0.05, na.rm = TRUE)) {
      notes <- c(notes, paste(
        "Mardia's tests reject multivariate normality of the residuals.",
        if (identical(test, "Pillai")) {
          "The multivariate table uses Pillai's trace, the most robust of the four statistics to this."
        } else {
          "Pillai's trace is the most robust of the four multivariate statistics to this; consider test = \"Pillai\"."
        }))
    }
    bm <- .box_m(Ycov, cellv)
    if (is.null(bm)) {
      notes <- c(notes, "Box's M could not be computed: fewer than two cells have complete data.")
    } else {
      assum$box_m <- data.frame(statistic = bm$statistic, df = bm$df,
                                p_value = bm$p_value, stringsAsFactors = FALSE)
      notes <- c(notes, bm$note)
      if (!is.na(bm$p_value) && bm$p_value < 0.001) {
        n_cell <- table(cellv)
        balanced <- length(unique(n_cell[n_cell > 0L])) == 1L
        pillai <- identical(test, "Pillai")
        notes <- c(notes, paste(
          "Box's M rejects equality of the group covariance matrices. The test is very sensitive to non-normality, so inspect the group covariances before acting on it.",
          if (balanced) {
            paste("With equal group sizes the multivariate tests are fairly robust to unequal covariances, Pillai's trace most of all",
                  if (pillai) "(it is the statistic used here)." else "(consider test = \"Pillai\").")
          } else {
            paste("With unequal group sizes as well, the multivariate tests can be too liberal or too conservative;",
                  if (pillai) "Pillai's trace, used here, is the least affected." else "prefer test = \"Pillai\".")
          }))
      }
    }
  }
  assum$cell_counts <- .cell_counts(d, groups)

  ## Canonical discriminant analysis -----------------------------------------
  canonical <- NULL
  cds <- NULL
  if (multivariate) {
    term_vars <- ct$vars
    cds <- .canonical_discriminant(
      Ycda, ssp$SSP[[canon_label]], ssp$SSPE,
      df_h = ssp$df[[canon_label]], df_e = ssp$error_df,
      group = .cell_vector(d, term_vars))
    if (is.null(cds$canonical)) {
      notes <- c(notes, sprintf(
        "The canonical discriminant analysis of `%s` could not be computed: %s.",
        canonical_term, cds$reason))
      cds <- NULL
    } else {
      canonical <- cds$canonical
      attr(canonical, "term") <- canonical_term
      assum$structure_coefficients <- cds$structure
      if (length(term_labels) > 1L) {
        notes <- c(notes, sprintf(
          "The canonical discriminant analysis describes the `%s` term only: %s. The model has %d terms (%s); the canonical scores have the fitted effects of the other terms removed, so they show `%s` alone.",
          canonical_term, canon_why, length(term_labels),
          paste(gsub("`", "", term_labels, fixed = TRUE), collapse = ", "),
          canonical_term))
      }
    }
  }

  ## Plots -------------------------------------------------------------------
  plot_list <- list()
  if (plots) {
    .say(verbose, "Building plots.")
    for (r in responses) {
      m <- univariate[[r]]$model
      plot_list[[paste0("residuals_", r)]] <- .plot_resid_fitted(
        stats::fitted(m), stats::residuals(m),
        title = sprintf("Residuals vs fitted: %s", r))
      plot_list[[paste0("qq_", r)]] <- .plot_qq(
        stats::residuals(m), title = sprintf("Normal Q-Q plot: %s", r))
    }
    if (!is.null(emm_all)) {
      pe <- .plot_manova_emmeans(emm_all, groups, conf_level,
                                 response_col = resp_col, term_col = term_col)
      if (!is.null(pe)) plot_list$emmeans <- pe
    }
    if (is.null(plot_list$emmeans)) notes <- c(notes, .no_emm_plot_note())
    if (!is.null(cds)) {
      plot_list$canonical <- .plot_canonical(
        cds$scores, cds$group, paste(ct$vars, collapse = " : "),
        context = if (length(term_labels) > 1L) {
          sprintf("Axes for the `%s` term; other terms' effects removed",
                  canonical_term)
        } else NULL)
    }
  }

  .new_fit(
    method       = if (multivariate) {
      sprintf("Multivariate analysis of %s (%s test, Type %s)",
              if (is.null(covariates)) "variance" else "covariance", test,
              type_used)
    } else {
      sprintf("Univariate analysis of %s (Type %s)",
              if (is.null(covariates)) "variance" else "covariance", type_used)
    },
    call         = cl,
    model        = fit,
    anova        = mv_table %||% univariate[[1L]]$anova,
    effect_sizes = effect_sizes,
    emmeans      = emm_all,
    posthoc      = ph_all,
    assumptions  = assum,
    plots        = plot_list,
    data_used    = d,
    n_removed    = prep$n_removed,
    conf_level   = conf_level,
    notes        = notes,
    extra        = list(univariate = univariate,
                        multivariate = mv_table,
                        canonical = canonical,
                        canonical_term = canonical_term,
                        slopes_test = slopes_test,
                        covariate_means = covariate_means,
                        test = test)
  )
}

#' Hypothesis and error SSCP matrices of the Type II or III multivariate tests
#'
#' \code{car::Manova()} is used whenever it accepts the model. It refuses a
#' model with aliased coefficients; the Type II matrices are then computed by
#' the model comparisons that define them -- for each term, the error SSCP of
#' the model without the term and the terms that contain it, minus that of the
#' model without only the terms that contain it -- which is what
#' \code{car::Anova()} does for a univariate model in the same situation.
#' Degrees of freedom are rank differences, so a fully aliased term has none.
#'
#' @param fit the multivariate \code{mlm}, fitted under the contrasts
#'   \code{type} needs.
#' @param type \code{"II"} or \code{"III"}.
#' @return list with \code{SSP} (named list of hypothesis matrices, keyed by
#'   the model's term labels), \code{SSPE}, \code{df} (named), \code{error_df},
#'   \code{terms} and \code{note}.
#' @noRd
.manova_ssp <- function(fit, type = "II") {
  note <- character(0)
  if (!anyNA(stats::coef(fit))) {
    got <- .collect_conditions(
      car::Manova(fit, type = if (identical(type, "III")) 3L else 2L))
    M <- got$value
    if (!inherits(M, "error")) {
      keep <- M$terms != "(Intercept)"
      tl <- M$terms[keep]
      return(list(SSP = M$SSP[tl], SSPE = M$SSPE,
                  df = stats::setNames(as.numeric(M$df[keep]), tl),
                  error_df = M$error.df, terms = tl, type = type,
                  note = .said_note(got$said, "car::Manova()")))
    }
    note <- sprintf("car::Manova() failed (%s), so the multivariate tests were computed by model comparison.",
                    conditionMessage(M))
  } else {
    note <- "The model has aliased coefficients, which car::Manova() does not accept, so the Type II multivariate tests were computed by the equivalent model comparisons; a term that is entirely aliased has no degrees of freedom and no test."
  }

  X <- stats::model.matrix(fit)
  asg <- attr(X, "assign")
  Y <- stats::fitted(fit) + stats::residuals(fit)
  tt <- stats::terms(fit)
  tl <- attr(tt, "term.labels")
  fac <- attr(tt, "factors") > 0
  sspe <- function(cols) {
    q <- qr(X[, cols, drop = FALSE])
    r <- qr.resid(q, Y)
    list(E = crossprod(r), rank = q$rank)
  }
  full <- sspe(seq_len(ncol(X)))
  SSP <- stats::setNames(vector("list", length(tl)), tl)
  df <- stats::setNames(numeric(length(tl)), tl)
  for (j in seq_along(tl)) {
    # Terms that contain term j (its higher-order relatives).
    rel <- which(vapply(seq_along(tl), function(k) {
      k != j && all(fac[, k] | !fac[, j])
    }, logical(1)))
    if (identical(type, "III")) {
      a <- sspe(which(asg != j))
      b <- full
    } else {
      a <- sspe(which(!asg %in% c(j, rel)))
      b <- sspe(which(!asg %in% rel))
    }
    H <- a$E - b$E
    SSP[[j]] <- (H + t(H)) / 2
    df[[j]] <- b$rank - a$rank
  }
  list(SSP = SSP, SSPE = full$E, df = df,
       error_df = nrow(X) - full$rank, terms = tl, type = type, note = note)
}

#' The multivariate test table from hypothesis and error SSCP matrices
#'
#' The eigenvalues of \code{solve(E, H)} are computed on responses scaled to
#' unit error variance, which leaves them unchanged but keeps responses on
#' very different scales from making \code{E} look singular. The F
#' approximations are the standard ones (as in \code{summary.manova} and
#' \pkg{car}).
#' @noRd
.manova_table <- function(ssp, test, tss = NULL) {
  E <- ssp$SSPE
  p <- ncol(E)
  dfe <- ssp$error_df
  if (!is.finite(dfe) || dfe < p) {
    stop(sprintf("the residuals have %s degree(s) of freedom, fewer than the %d responses",
                 format(dfe), p), call. = FALSE)
  }
  # A response the model fits exactly leaves an error matrix of rounding
  # errors, which the scaling below would otherwise dress up as information.
  if (!is.null(tss) && any(diag(E) <= 1e-12 * tss)) {
    exact <- colnames(E)[diag(E) <= 1e-12 * tss]
    stop(sprintf("the model fits %s exactly, leaving no residual variation",
                 if (length(exact) > 0L) paste(sprintf("`%s`", exact), collapse = ", ")
                 else "a response"), call. = FALSE)
  }
  w <- .whiten_sscp(E)
  if (is.null(w)) stop("the error sum of squares and cross-products matrix is singular", call. = FALSE)
  rows <- lapply(ssp$terms, function(tm) {
    q <- ssp$df[[tm]]
    H <- ssp$SSP[[tm]]
    if (!is.finite(q) || q < 1 || is.null(H)) return(rep(NA_real_, 5L))
    Hs <- H / outer(w$scale, w$scale)
    S <- w$isqrt %*% Hs %*% w$isqrt
    eig <- eigen((S + t(S)) / 2, symmetric = TRUE, only.values = TRUE)$values
    st <- .mv_statistic(eig, q, dfe, test)
    ok <- all(is.finite(st)) && st[2L] >= 0 && st[3L] > 0 && st[4L] > 0
    c(st, if (ok) stats::pf(st[2L], st[3L], st[4L], lower.tail = FALSE) else NA_real_)
  })
  m <- do.call(rbind, rows)
  out <- data.frame(term = gsub("`", "", ssp$terms, fixed = TRUE),
                    df = unname(as.numeric(ssp$df[ssp$terms])),
                    statistic = m[, 1L], approx_f = m[, 2L],
                    num_df = m[, 3L], den_df = m[, 4L], p_value = m[, 5L],
                    stringsAsFactors = FALSE, row.names = NULL)
  attr(out, "test") <- test
  attr(out, "statistic") <- sprintf("%s, approximate F", switch(test,
    Pillai = "Pillai's trace", Wilks = "Wilks' lambda",
    `Hotelling-Lawley` = "Hotelling-Lawley trace", Roy = "Roy's largest root"))
  attr(out, "ss_type") <- ssp$type
  out
}

#' Scale an SSCP matrix to unit diagonal and take its inverse square root
#' @return list with \code{scale} (the square roots of the diagonal) and
#'   \code{isqrt}, or \code{NULL} when the matrix is singular.
#' @noRd
.whiten_sscp <- function(E) {
  s <- sqrt(diag(E))
  if (any(!is.finite(s) | s <= 0)) return(NULL)
  Es <- E / outer(s, s)
  ev <- eigen((Es + t(Es)) / 2, symmetric = TRUE)
  if (min(ev$values) <= max(ev$values) * 1e-12) return(NULL)
  isqrt <- ev$vectors %*% (t(ev$vectors) / sqrt(ev$values))
  list(scale = s, isqrt = isqrt)
}

#' A multivariate test statistic and its F approximation
#'
#' @param eig eigenvalues of \code{solve(E, H)}.
#' @param q hypothesis degrees of freedom.
#' @param df_res error degrees of freedom.
#' @return numeric: statistic, approximate F, numerator df, denominator df.
#' @noRd
.mv_statistic <- function(eig, q, df_res, test) {
  p <- length(eig)
  switch(test,
    Pillai = {
      stat <- sum(eig / (1 + eig))
      s <- min(p, q)
      n <- 0.5 * (df_res - p - 1)
      m <- 0.5 * (abs(p - q) - 1)
      t1 <- 2 * m + s + 1
      t2 <- 2 * n + s + 1
      c(stat, (t2 / t1 * stat) / (s - stat), s * t1, s * t2)
    },
    Wilks = {
      stat <- prod(1 / (1 + eig))
      t1 <- df_res - 0.5 * (p - q + 1)
      t2 <- (p * q - 2) / 4
      t3 <- p^2 + q^2 - 5
      t3 <- if (t3 > 0) sqrt(((p * q)^2 - 4) / t3) else 1
      c(stat, ((stat^(-1 / t3) - 1) * (t1 * t3 - 2 * t2)) / p / q,
        p * q, t1 * t3 - 2 * t2)
    },
    `Hotelling-Lawley` = {
      stat <- sum(eig)
      m <- 0.5 * (abs(p - q) - 1)
      n <- 0.5 * (df_res - p - 1)
      s <- min(p, q)
      t1 <- 2 * m + s + 1
      t2 <- 2 * (s * n + 1)
      c(stat, (t2 * stat) / s / s / t1, s * t1, t2)
    },
    Roy = {
      stat <- max(eig)
      t1 <- max(p, q)
      t2 <- df_res - t1 + q
      c(stat, (t2 * stat) / t1, t1, t2)
    })
}

#' Which term the canonical discriminant analysis describes
#'
#' The full interaction of the grouping variables when the model contains it,
#' otherwise the main effect of the first grouping variable. Terms are matched
#' through the model's own \code{factors} attribute, so the backtick quoting of
#' non-syntactic names cannot make a match fail.
#' @return list with \code{label} (the model's term label), \code{vars} (the
#'   grouping variables in it) and \code{why}.
#' @noRd
.choose_canonical_term <- function(fit, groups) {
  tt <- stats::terms(fit)
  tl <- attr(tt, "term.labels")
  fac <- attr(tt, "factors")
  vars_of <- lapply(tl, function(t) {
    gsub("`", "", rownames(fac)[fac[, t] > 0], fixed = TRUE)
  })
  full <- which(vapply(vars_of, function(v) setequal(v, groups), logical(1)))
  if (length(full) > 0L) {
    why <- if (length(groups) > 1L) {
      "it is the full interaction of the grouping variables, the term whose levels are the cells of the design"
    } else "it is the grouping term"
    return(list(label = tl[full[1L]], vars = groups, why = why))
  }
  first <- which(vapply(vars_of, function(v) identical(v, groups[1L]), logical(1)))
  list(label = tl[first[1L]], vars = groups[1L],
       why = sprintf("the model does not contain the full interaction of the grouping variables, so the main effect of the first grouping variable is used (list another variable first in `groups` to describe it instead)"))
}

#' Remove the fitted contributions of model terms from the responses
#'
#' The coefficients are those of the full model under sum-to-zero coding, so
#' removing a grouping term takes out its effect averaged over the others, and
#' removing a (centred) covariate reads the responses at the covariate mean.
#' Aliased coefficients count as zero, which is what the fit itself does.
#'
#' @return list with \code{covariates} (responses with the covariate effects
#'   removed) and \code{canonical} (with every term but \code{keep_term}
#'   removed), or \code{NULL} on failure.
#' @noRd
.partial_out <- function(fit, Y, covariates, keep_term) {
  tryCatch({
    tt <- stats::terms(fit)
    mf <- stats::model.frame(fit)
    # Contrasts stored on a factor would override the sum-to-zero coding.
    for (j in seq_along(mf)) {
      if (is.factor(mf[[j]])) attr(mf[[j]], "contrasts") <- NULL
    }
    X <- .fit_with_contrasts(function() stats::model.matrix(tt, mf), "III")
    asg <- attr(X, "assign")
    tl <- attr(tt, "term.labels")
    fac <- attr(tt, "factors")
    B <- qr.coef(qr(X), Y)
    B <- as.matrix(B)
    B[is.na(B)] <- 0
    cov_terms <- which(vapply(tl, function(t) {
      v <- gsub("`", "", rownames(fac)[fac[, t] > 0], fixed = TRUE)
      length(v) > 0L && all(v %in% covariates)
    }, logical(1)))
    drop_cols <- function(cols) {
      if (!any(cols)) return(Y)
      Y - X[, cols, drop = FALSE] %*% B[cols, , drop = FALSE]
    }
    Ycov <- drop_cols(asg %in% cov_terms)
    Ycda <- drop_cols(asg != 0L & asg != match(keep_term, tl))
    if (!all(is.finite(Ycov)) || !all(is.finite(Ycda))) stop("non-finite")
    dimnames(Ycov) <- dimnames(Ycda) <- list(NULL, colnames(Y))
    list(covariates = Ycov, canonical = Ycda)
  }, error = function(e) NULL)
}

#' Multivariate (or univariate) test of homogeneous regression slopes
#'
#' Compares the model with a common slope for each covariate against one in
#' which each covariate interacts with every grouping term of the model.
#' @return list with \code{table} and \code{note}.
#' @noRd
.slopes_test <- function(d, lhs, rhs, covariates, g_terms, test, multivariate) {
  int_terms <- as.vector(outer(.bq(covariates), g_terms, paste, sep = ":"))
  fml_add <- stats::reformulate(rhs, response = lhs)
  fml_int <- stats::reformulate(c(rhs, int_terms), response = lhs)
  got <- .collect_conditions({
    m_add <- .fit_lm(fml_add, d)
    m_int <- .fit_lm(fml_int, d)
    # With no residual variation left there is nothing to compare: the test
    # statistic would be a ratio of rounding errors.
    res <- as.matrix(stats::residuals(m_add))
    y <- as.matrix(stats::fitted(m_add)) + res
    tss <- colSums(sweep(y, 2L, colMeans(y))^2)
    if (any(colSums(res^2) <= 1e-12 * tss)) {
      stop("the common-slope model fits a response exactly", call. = FALSE)
    }
    if (multivariate) stats::anova(m_add, m_int, test = test) else
      stats::anova(m_add, m_int)
  })
  cmp <- got$value
  comparison <- "common slopes vs covariate-by-group slopes"
  empty <- function() {
    data.frame(comparison = comparison, df = NA_real_, statistic = NA_real_,
               approx_f = NA_real_, num_df = NA_real_, den_df = NA_real_,
               p_value = NA_real_, homogeneous = NA, stringsAsFactors = FALSE)
  }
  if (inherits(cmp, "error")) {
    return(list(table = empty(), note = sprintf(
      "The homogeneity-of-slopes test could not be computed (%s). Check that each covariate varies within every group.",
      conditionMessage(cmp))))
  }
  cmp <- as.data.frame(cmp)
  row <- cmp[2L, , drop = FALSE]
  if (multivariate) {
    tab <- data.frame(comparison = comparison, df = abs(row[["Df"]]),
                      statistic = row[[test]], approx_f = row[["approx F"]],
                      num_df = row[["num Df"]], den_df = row[["den Df"]],
                      p_value = row[["Pr(>F)"]], stringsAsFactors = FALSE)
  } else {
    tab <- data.frame(comparison = comparison, df = abs(row[["Df"]]),
                      statistic = row[["F"]], approx_f = row[["F"]],
                      num_df = abs(row[["Df"]]), den_df = row[["Res.Df"]],
                      p_value = row[["Pr(>F)"]], stringsAsFactors = FALSE)
  }
  tab$homogeneous <- if (is.na(tab$p_value)) NA else tab$p_value >= 0.05
  row.names(tab) <- NULL
  attr(tab, "test") <- if (multivariate) test else "F"
  note <- .said_note(got$said, "the homogeneity-of-slopes test")
  if (isFALSE(tab$homogeneous)) {
    note <- c(note, sprintf(
      "Homogeneity of regression slopes is rejected (%s = %.4g, approximate F(%g, %g) = %.4g, p = %.3g): the covariate slopes differ across groups. The common-slope %s, its adjusted marginal means and their comparisons then describe groups at the covariate mean only, not constant adjusted differences, and rejections by Mardia's tests or Box's M may be artefacts of the misspecified slopes. Model the covariate-by-group interaction (e.g. with car::Manova on your own model), or analyse each response with anova_ancova().",
      if (multivariate) test else "F", tab$statistic, tab$num_df, tab$den_df,
      tab$approx_f, tab$p_value, if (multivariate) "MANCOVA" else "ANCOVA"))
  }
  list(table = tab, note = note)
}

#' Faceted marginal-means plot across responses
#'
#' One panel per response. A single grouping variable is drawn on the x axis;
#' with several, the remaining ones are mapped to colour and dodged (a cell
#' grid), or, when the factors enter additively and the means are per factor,
#' one column of panels per factor.
#'
#' @param emm_all the stacked table.
#' @param response_col the column naming the response.
#' @param term_col the column naming the factor of a per-factor table, or
#'   \code{NULL}.
#' @noRd
.plot_manova_emmeans <- function(emm_all, groups, conf_level,
                                 response_col = "response", term_col = NULL) {
  df <- as.data.frame(emm_all)
  protect <- c(groups, response_col, term_col)
  est <- .emm_col(df, "estimate", protect)
  lo  <- .emm_col(df, "conf_low",  protect)
  hi  <- .emm_col(df, "conf_high", protect)
  if (is.null(est) || !response_col %in% names(df) ||
      !all(groups %in% names(df))) return(NULL)
  has_ci <- !is.null(lo) && !is.null(hi)
  resp <- .safe_name(df, ".response")
  df[[resp]] <- factor(df[[response_col]], levels = unique(df[[response_col]]))
  subtitle <- if (has_ci) {
    sprintf("Error bars are %s%% confidence intervals", .pct(conf_level))
  } else NULL
  theme <- ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.position = "bottom")

  if (!is.null(term_col) && term_col %in% names(df)) {
    # Per-factor means of an additive model: each row belongs to one factor.
    lvl <- rep(NA_character_, nrow(df))
    for (g in groups) {
      v <- as.character(df[[g]])
      lvl[!is.na(v)] <- v[!is.na(v)]
    }
    lvl_col <- .safe_name(df, ".level")
    df[[lvl_col]] <- factor(lvl, levels = unique(lvl))
    trm <- .safe_name(df, ".term")
    df[[trm]] <- factor(df[[term_col]], levels = unique(df[[term_col]]))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[lvl_col]], y = .data[[est]])) +
      ggplot2::geom_point(size = 2.5)
    if (has_ci) {
      p <- p + ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data[[lo]], ymax = .data[[hi]]), width = 0.15)
    }
    return(p +
      ggplot2::facet_grid(rows = ggplot2::vars(.data[[resp]]),
                          cols = ggplot2::vars(.data[[trm]]), scales = "free") +
      ggplot2::labs(title = "Estimated marginal means by response",
                    subtitle = subtitle, x = NULL,
                    y = "Estimated marginal mean") +
      theme)
  }

  xvar <- groups[1L]
  df[[xvar]] <- as.factor(df[[xvar]])
  series <- .safe_name(df, ".series")
  if (length(groups) > 1L) {
    df[[series]] <- droplevels(interaction(df[groups[-1L]], drop = TRUE, sep = " : "))
    dodge <- ggplot2::position_dodge(width = 0.3)
    p <- ggplot2::ggplot(df, ggplot2::aes(
      x = .data[[xvar]], y = .data[[est]],
      colour = .data[[series]], group = .data[[series]])) +
      ggplot2::geom_line(linewidth = 0.7, position = dodge) +
      ggplot2::geom_point(size = 2.5, position = dodge)
    if (has_ci) {
      p <- p + ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data[[lo]], ymax = .data[[hi]]), width = 0.15,
        position = dodge)
    }
    lab_colour <- paste(groups[-1L], collapse = " : ")
  } else {
    # One factor: a line only when its levels are ordered. Between unordered
    # categories it would suggest a trend, as in the other functions' plots.
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[xvar]], y = .data[[est]],
                                          group = 1L))
    if (is.ordered(df[[xvar]])) {
      p <- p + ggplot2::geom_line(linewidth = 0.7, colour = "grey40")
    }
    p <- p + ggplot2::geom_point(size = 2.5)
    if (has_ci) {
      p <- p + ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data[[lo]], ymax = .data[[hi]]), width = 0.15)
    }
    lab_colour <- NULL
  }
  p +
    ggplot2::facet_wrap(ggplot2::vars(.data[[resp]]), scales = "free_y") +
    ggplot2::labs(
      title = "Estimated marginal means by response",
      subtitle = subtitle, x = xvar, y = "Estimated marginal mean",
      colour = lab_colour) +
    theme
}
