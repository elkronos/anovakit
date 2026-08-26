#' Multivariate Analysis of Variance and Covariance
#'
#' Compares two or more numeric responses jointly across groups, optionally
#' adjusting for covariates. Reports the multivariate test, follow-up univariate
#' analyses with effect sizes, covariate-adjusted estimated marginal means,
#' assumption checks (Mardia's multivariate normality tests and Box's M), and a
#' canonical discriminant analysis with its plot.
#'
#' @details
#' \strong{Assumption checks.} Mardia's tests of multivariate skewness and
#' kurtosis and Box's M test of equality of covariance matrices are computed
#' from first principles, so no additional package is required. Box's M is
#' notoriously sensitive to non-normality: treat a small p-value as a prompt to
#' look at the group covariances rather than as a verdict.
#'
#' \strong{Canonical discriminant analysis.} The hypothesis and error sum of
#' squares and cross-products matrices are taken from the fitted model and the
#' generalised eigenproblem is solved directly. \code{$canonical} reports the
#' eigenvalues, canonical correlations and the proportion of between-group
#' variance on each axis; \code{$plots$canonical} plots the first two axes with
#' group centroids, and \code{$assumptions$structure_coefficients} gives the
#' pooled within-group correlations between each response and each axis. The
#' analysis describes \emph{one} term: the full grouping interaction when it is
#' in the model, and otherwise the grouping term itself. \code{$canonical_term}
#' names it, and a note repeats it. Axes are scaled to unit pooled within-group
#' variance, so the spread of the centroids on the plot means what it appears
#' to mean.
#'
#' \strong{Estimated marginal means} come from \pkg{emmeans} applied to each
#' univariate model, so with covariates present they are adjusted means, not
#' raw group means.
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
#'   turn the analysis into a MANCOVA.
#' @param test Character. Multivariate statistic to report: \code{"Pillai"}
#'   (default), \code{"Wilks"}, \code{"Hotelling-Lawley"} or \code{"Roy"}.
#'   Pillai's trace is the most robust to departures from the assumptions.
#' @param interaction \code{FALSE} (additive, the default), \code{TRUE} (full
#'   factorial across the grouping variables), or a whole number giving the
#'   highest interaction order.
#' @param type Character. \code{"II"} (default) or \code{"III"} sums of squares
#'   for the follow-up univariate tests. With \code{"III"} both the multivariate
#'   model and the univariate ones are fitted under sum-to-zero contrasts.
#' @param conf_level Numeric in (0, 1). Level for every interval returned.
#'   Default \code{0.95}.
#' @param adjust Character. Multiplicity adjustment for the pairwise
#'   comparisons within each response. Default \code{"tukey"}.
#' @param assumptions Logical. Compute Mardia's tests and Box's M. Default
#'   \code{TRUE}.
#' @param posthoc Logical. Compute pairwise comparisons within each response.
#'   Default \code{TRUE}.
#' @param plots Logical. Build \pkg{ggplot2} objects. They are returned in
#'   \code{$plots}, never drawn. Default \code{TRUE}.
#' @param verbose Logical. Emit progress through \code{\link[base]{message}}.
#'   Default \code{FALSE}.
#'
#' @return An \code{\link{anovakit_fit}} object with extra components
#'   \code{$univariate} (a named list of per-response fits, each holding that
#'   response's own \code{$model}, \code{$anova}, \code{$effect_sizes},
#'   \code{$emmeans}, \code{$emmeans_object} and \code{$posthoc}),
#'   \code{$multivariate} (the full multivariate table), \code{$canonical} (the
#'   discriminant axes), \code{$canonical_term} (the term they describe) and
#'   \code{$test} (the multivariate statistic used). The top-level
#'   \code{$emmeans} and \code{$posthoc} stack every response's table with a
#'   \code{response} column; \code{$emmeans_object} is \code{NULL}, because
#'   there is one grid per response -- take them from
#'   \code{$univariate[[r]]$emmeans_object}.
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
#' # As a MANCOVA: the marginal means are adjusted for age
#' anova_manova(d, c("score1", "score2"), "g",
#'              covariates = "age", plots = FALSE)$emmeans
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
  clash <- intersect(responses, c(groups, covariates))
  if (length(clash) > 0L) {
    .stopf("A column cannot be both a response and a predictor: %s.",
           paste(clash, collapse = ", "))
  }
  overlap <- intersect(groups, covariates)
  if (length(overlap) > 0L) {
    .stopf("A column cannot be both a group and a covariate: %s.",
           paste(overlap, collapse = ", "))
  }
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
  prep <- .prepare_frame(data, cols, factors = groups, keep_all = TRUE)
  d <- prep$data
  notes <- prep$notes
  .check_groups(d, groups, min_levels = 2L, min_n = 1L)

  multivariate <- length(responses) >= 2L
  if (!multivariate) {
    notes <- c(notes, "Only one response was supplied, so this is a univariate analysis; there is no multivariate test and no canonical discriminant analysis.")
  }

  added <- .add_cell(d, groups)
  d <- added$data
  cell <- added$cell

  rhs <- c(.bq(covariates), .group_terms(groups, interaction))

  ## Multivariate fit --------------------------------------------------------
  fit <- NULL
  mv_table <- NULL
  if (multivariate) {
    .say(verbose, "Fitting the multivariate model.")
    Y <- as.matrix(d[, responses, drop = FALSE])
    mf <- d
    mf[[".Y"]] <- Y
    fml <- stats::reformulate(rhs, response = ".Y")
    fit <- .fit_with_contrasts(function() stats::manova(fml, data = mf), type)
    smry <- tryCatch(summary(fit, test = test), error = function(e) e)
    if (inherits(smry, "error")) {
      .stopf("The multivariate test could not be computed: %s. This usually means the responses are collinear, or there are too few residual degrees of freedom for %d responses.",
             conditionMessage(smry), length(responses))
    }
    st <- as.data.frame(smry$stats)
    mv_table <- data.frame(
      term = rownames(st), df = st[["Df"]],
      statistic = st[[test]],
      approx_f = st[["approx F"]],
      num_df = st[["num Df"]], den_df = st[["den Df"]],
      p_value = st[["Pr(>F)"]],
      stringsAsFactors = FALSE, row.names = NULL)
    mv_table <- mv_table[mv_table$term != "Residuals", , drop = FALSE]
    attr(mv_table, "test") <- test
  }

  ## Univariate follow-ups ---------------------------------------------------
  .say(verbose, "Fitting %d univariate model(s).", length(responses))
  univariate <- stats::setNames(vector("list", length(responses)), responses)
  eff_rows <- list()
  emm_rows <- list()
  ph_rows <- list()
  for (r in responses) {
    fml_r <- stats::reformulate(rhs, response = .bq(r))
    m <- .fit_with_contrasts(function() stats::lm(fml_r, data = d), type)
    av <- .car_anova(m, type = type)
    # Aliasing is a property of the design, so it is the same for every
    # response; report it once rather than per response.
    if (r == responses[[1L]]) notes <- c(notes, .check_model_size(m))
    notes <- c(notes, av$note)
    y <- d[[r]]
    es <- .partial_eta_squared(av$raw, df_error = stats::df.residual(m),
                               n_obs = length(y))
    if (!is.null(attr(es, "omega_floored"))) {
      notes <- c(notes, sprintf(
        "Partial omega squared was negative for %s in the model for `%s` and has been floored at 0.",
        paste(attr(es, "omega_floored"), collapse = ", "), r))
    }
    if (nrow(es) > 0L) eff_rows[[r]] <- .cbind_label(es, "response", r)

    emm <- .emmeans_grid(m, groups, type = "link")
    notes <- c(notes, emm$note)
    tab <- .emmeans_table(emm$grid, conf_level, protect = groups)
    notes <- c(notes, tab$note)
    if (!is.null(tab$table)) emm_rows[[r]] <- .cbind_label(tab$table, "response", r)
    ph <- if (posthoc) {
      .emmeans_pairs(emm$grid, adjust = adjust, conf_level = conf_level)
    } else list(table = NULL,
                # One note for the whole fit, but $posthoc would have been the
                # comparisons for EVERY response stacked together, so the count
                # is per-response times the number of responses.
                note = if (r == responses[[1L]])
                  .no_posthoc_note(NROW(tab$table), each = length(responses))
                else character(0))
    notes <- c(notes, ph$note)
    if (!is.null(ph$table)) ph_rows[[r]] <- .cbind_label(ph$table, "response", r)

    univariate[[r]] <- list(model = m, anova = av$table, effect_sizes = es,
                            emmeans = tab$table, emmeans_object = emm$grid,
                            posthoc = ph$table)
  }
  effect_sizes <- if (length(eff_rows)) do.call(rbind, eff_rows) else NULL
  emm_all <- if (length(emm_rows)) do.call(rbind, emm_rows) else NULL
  ph_all  <- if (length(ph_rows))  do.call(rbind, ph_rows)  else NULL
  for (x in list(effect_sizes, emm_all, ph_all)) if (!is.null(x)) row.names(x) <- NULL

  ## Assumptions -------------------------------------------------------------
  # With covariates present the assumptions are about the CONDITIONAL
  # distribution, so both Mardia and Box's M are computed on the responses with
  # the covariates partialled out. Using the raw responses instead rejects
  # homogeneity whenever a covariate's spread differs by group, on data where
  # the assumption holds exactly.
  Yadj <- as.matrix(d[, responses, drop = FALSE])
  if (!is.null(covariates) && multivariate && !is.null(fit)) {
    # Subtract the covariate contribution as estimated in the FULL model, with
    # the groups in it. Regressing the responses on the covariates alone would
    # also strip the between-group signal the covariates happen to carry,
    # leaving scores that discriminate far less than their own eigenvalue.
    adj <- tryCatch({
      B <- stats::coef(fit)
      cov_rows <- intersect(rownames(B), .bq(covariates))
      if (length(cov_rows) == 0L) cov_rows <- intersect(rownames(B), covariates)
      Xc <- as.matrix(d[, covariates, drop = FALSE])
      Xc <- sweep(Xc, 2L, colMeans(Xc), "-")
      Yadj - Xc %*% B[cov_rows, , drop = FALSE]
    }, error = function(e) NULL)
    if (is.null(adj) || !all(dim(adj) == dim(Yadj))) {
      notes <- c(notes, "The covariates could not be partialled out of the responses, so the assumption checks below describe the unadjusted responses.")
    } else {
      Yadj <- adj
      colnames(Yadj) <- responses
      notes <- c(notes, "Assumption checks and the canonical discriminant analysis use the responses with the covariates partialled out, taking the covariate effect as estimated in the full model, which is the space the multivariate test operates in.")
    }
  }

  assum <- list()
  if (assumptions && multivariate) {
    .say(verbose, "Checking multivariate assumptions.")
    Ymat <- Yadj
    resid_mat <- if (!is.null(fit)) stats::residuals(fit) else Ymat
    assum$mardia <- .mardia_test(resid_mat)
    if (is.null(assum$mardia)) {
      notes <- c(notes, "Mardia's tests could not be computed (too few observations, or a singular covariance matrix).")
    } else if (any(assum$mardia$p_value < 0.05, na.rm = TRUE)) {
      notes <- c(notes, "Mardia's tests reject multivariate normality of the residuals. Pillai's trace is the most robust of the four multivariate statistics; consider it if you have not already.")
    }
    bm <- .box_m(Ymat, d[[cell]])
    assum$box_m <- if (is.null(bm)) NULL else
      data.frame(statistic = bm$statistic, df = bm$df, p_value = bm$p_value,
                 stringsAsFactors = FALSE)
    if (!is.null(bm)) notes <- c(notes, bm$note)
    if (!is.null(bm) && !is.na(bm$p_value) && bm$p_value < 0.001) {
      notes <- c(notes, "Box's M rejects equality of the group covariance matrices. The test is very sensitive to non-normality, so inspect the group covariances before acting on it; with unequal group sizes as well, prefer Pillai's trace.")
    }
  }
  assum$cell_counts <- .cell_counts(d, groups)

  ## Canonical discriminant analysis -----------------------------------------
  canonical <- NULL
  cds <- NULL
  canonical_term <- NA_character_
  if (multivariate && !is.null(fit)) {
    smry_ss <- smry
    # The SSCP list and the stats table are keyed by the model's own term
    # labels, which are backtick-quoted whenever a column name is not
    # syntactic. Match on those rather than on the bare column name.
    avail <- setdiff(names(smry_ss$SS), "Residuals")
    wanted <- c(paste(.bq(groups), collapse = ":"),
                paste(groups, collapse = ":"), .bq(groups)[1L], groups[1L])
    canonical_term <- utils::head(intersect(wanted, avail), 1L)
    if (length(canonical_term) == 0L) canonical_term <- utils::head(avail, 1L)
    H <- smry_ss$SS[[canonical_term]]
    E <- smry_ss$SS[["Residuals"]]
    df_h <- smry_ss$stats[canonical_term, "Df"]
    df_e <- smry_ss$stats["Residuals", "Df"]
    if (!is.null(H) && !is.null(E)) {
      # Colour by the factor(s) the hypothesis matrix actually describes, not
      # by the full cell grid, which would label the plot with a grouping the
      # axes were not derived from.
      term_vars <- intersect(
        unlist(strsplit(gsub("`", "", canonical_term, fixed = TRUE), ":",
                        fixed = TRUE)), groups)
      if (length(term_vars) == 0L) term_vars <- groups
      cds <- .canonical_discriminant(Yadj, H, E, df_h, df_e = df_e,
                                     group = .cell_vector(d, term_vars))
      if (is.null(cds)) {
        notes <- c(notes, "The canonical discriminant analysis could not be computed (singular error matrix).")
      } else {
        canonical <- cds$canonical
        attr(canonical, "term") <- canonical_term
        assum$structure_coefficients <- cds$structure
        if (length(avail) > 1L) {
          notes <- c(notes, sprintf(
            "The canonical discriminant analysis describes the `%s` term only; the model has %d terms (%s).",
            canonical_term, length(avail), paste(avail, collapse = ", ")))
        }
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
      plot_list$emmeans <- .plot_manova_emmeans(emm_all, groups, conf_level)
    }
    if (!is.null(cds)) {
      plot_list$canonical <- .plot_canonical(
        cds$scores, cds$group,
        gsub("`", "", canonical_term, fixed = TRUE))
    }
  }

  .new_fit(
    method       = if (multivariate) {
      sprintf("Multivariate analysis of %s (%s test, Type %s follow-ups)",
              if (is.null(covariates)) "variance" else "covariance", test, type)
    } else {
      sprintf("Univariate analysis of %s (Type %s)",
              if (is.null(covariates)) "variance" else "covariance", type)
    },
    call         = cl,
    model        = fit %||% univariate[[1L]]$model,
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
                        test = test)
  )
}

#' Faceted marginal-means plot across responses
#' @noRd
.plot_manova_emmeans <- function(emm_all, groups, conf_level) {
  df <- as.data.frame(emm_all)
  xvar <- groups[1L]
  est <- .emm_col(df, "estimate", c(groups, "response"))
  lo  <- .emm_col(df, "conf_low",  c(groups, "response"))
  hi  <- .emm_col(df, "conf_high", c(groups, "response"))
  if (!xvar %in% names(df) || is.null(est)) return(NULL)
  df[[xvar]] <- as.factor(df[[xvar]])
  has_ci <- !is.null(lo) && !is.null(hi)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[xvar]], y = .data[[est]],
                                        group = 1L)) +
    ggplot2::geom_line(linewidth = 0.7, colour = "grey40") +
    ggplot2::geom_point(size = 2.5)
  if (has_ci) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data[[lo]], ymax = .data[[hi]]), width = 0.15)
  }
  p +
    ggplot2::facet_wrap(~ .data[["response"]], scales = "free_y") +
    ggplot2::labs(
      title = "Estimated marginal means by response",
      subtitle = if (has_ci) sprintf("Error bars are %g%% confidence intervals",
                                     round(conf_level * 100, 1)) else NULL,
      x = xvar, y = "Estimated marginal mean") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
