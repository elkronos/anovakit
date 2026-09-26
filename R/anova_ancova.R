#' Analysis of Covariance
#'
#' Compares the mean of a numeric response across groups while adjusting for one
#' or more numeric covariates. Tests the homogeneity-of-slopes assumption,
#' checks residual normality and equality of variances, reports Type II or Type
#' III sums of squares with partial eta squared, and returns covariate-adjusted
#' estimated marginal means.
#'
#' @details
#' \strong{Grouping structure.} With two or more grouping variables the default,
#' \code{interaction = TRUE}, fits their full factorial -- \code{y ~ x + A * B}
#' -- which is what a factorial ANCOVA is. \code{interaction = FALSE} enters
#' them additively and a whole number keeps interactions up to that order, as
#' in \code{\link{anova_glm}}. The structure chosen is the group part of both
#' candidate models below. Under the full factorial \code{$emmeans} and
#' \code{$posthoc} compare the cells, and a cell with no data is left out with a
#' note; with the grouping variables entered additively they are reported for
#' each factor separately, averaged over the others.
#'
#' \strong{Covariates are centred.} Every covariate is mean-centred before the
#' model is fitted. That moves the intercept and, when the covariate-by-group
#' terms are in the model and \code{type = "III"}, the group row of the ANOVA
#' table: that row tests the group difference \emph{at covariate = 0}, which for
#' an uncentred covariate may be a point nowhere near the data, and the p-value
#' can be anything. On the example below -- a baseline distributed around 100
#' -- the same data gives p = 0.35 uncentred against p = 5.5e-21 centred. In the
#' additive model, and under Type II, the group row is the same either way.
#' \code{$emmeans}, \code{$posthoc} and \code{$simple_slopes} do not depend on
#' centring at all: they are always evaluated with every covariate held at its
#' mean, including a covariate with only two values (a 0/1 indicator, say),
#' which is held at its mean rather than averaged over its two values. Set
#' \code{center_covariates = FALSE} only if the covariate's own zero is the
#' point you want the intercept (and the Type III group row) to refer to. When
#' the covariate mean lies outside the range a group was observed over, that
#' group's adjusted mean is an extrapolation along the fitted slope, and
#' \code{$notes} names the group.
#'
#' A covariate that is constant, that is determined by the grouping variables
#' (constant within every group, such as a group-level attribute), or that is
#' an exact linear combination of the grouping variables and the other
#' covariates is refused with an error naming it: its effect cannot be
#' separated from the group effect, so there is nothing to adjust for.
#'
#' \strong{Homogeneity of slopes.} The two candidate models are the grouping
#' structure with the covariates entered additively, and the same model plus
#' the product of every covariate with every grouping term in it, so that under
#' the full factorial each cell has its own slope. They are compared by one
#' nested F test, joint over all covariates: heterogeneity in any one covariate
#' retains the slope terms for all of them. \code{homogeneity_alpha} is its
#' level. The test is reported as not computable (\code{NA}) when the
#' covariate-by-group terms add nothing estimable, or when the model without
#' them already fits the data exactly. If the slopes differ, the interaction
#' model is used, and a note explains that the group effect is then a
#' comparison at one point on the covariate rather than a constant difference,
#' and that \code{$simple_slopes} is the thing to read.
#'
#' Choosing the model with a test on the same data has a cost, in both
#' directions. When the slopes really differ but the test misses the
#' difference, the adjusted comparisons of the additive model are biased by the
#' slope difference times the difference in covariate means. That vanishes on
#' average when the groups are randomised but not otherwise, and the test often
#' lacks the power to
#' find a slope difference large enough to matter. In simulations with three
#' groups of 20, slopes 0.5, 1 and 1.5, and covariate means at -1, 0 and 1
#' within-group SDs, the default pipeline rejected a true null of equal
#' adjusted means for the group row about 11--13\% of the time at the 5\% level
#' (about 34\% with 10 per group and means at -2, 0 and 2), against about 5\%
#' with \code{force_interaction = TRUE}; always fitting the additive model
#' (\code{force_interaction = FALSE}) was worse still. Randomised groups were
#' not affected. So when the additive model is used and the covariate means
#' differ between groups by more than half a pooled within-group SD,
#' \code{$notes} says so. Conversely, when the test retains the interaction, the p-values do
#' not allow for that choice and can also be too small. \code{force_interaction
#' = TRUE} avoids the selection step: with centred covariates and \code{type =
#' "III"}, the group row is then a valid comparison of the groups at the
#' covariate mean whether or not the slopes differ, which makes it the safer
#' choice when the groups were not randomised. \code{$notes} always reports
#' which model was fitted and why.
#'
#' \strong{Robust standard errors.} With \code{vcov_type} other than
#' \code{"model"}, the F tests in \code{$anova} are Wald F tests
#' (\code{car::Anova(vcov. = )}) built from the heteroscedasticity-consistent
#' covariance, on the model's residual degrees of freedom, so that table has no
#' sums of squares. \code{$emmeans}, \code{$posthoc} and \code{$simple_slopes}
#' use the same covariance. \code{$effect_sizes} are still computed from the
#' model-based sums of squares, which the covariance does not change, and the
#' homogeneity-of-slopes test stays the model-based F test.
#'
#' @param data A data frame, or anything inheriting from one, such as a
#'   \code{data.table} or a tibble.
#' @param response Character. Name of the numeric response column.
#' @param groups Character vector. One or more grouping columns.
#' @param covariates Character vector. One or more numeric covariate columns.
#' @param center_covariates Logical. Mean-centre the covariates before fitting.
#'   Default \code{TRUE}. See Details.
#' @param force_interaction Logical or \code{NULL}. \code{NULL} (default) lets
#'   the slopes test decide; \code{TRUE} always fits the covariate-by-group
#'   interaction; \code{FALSE} always fits the additive model. See Details for
#'   what the test-based choice costs.
#' @param homogeneity_alpha Numeric in (0, 1). Significance level for the slopes
#'   test. Default \code{0.05}.
#' @param interaction How the grouping variables combine: \code{TRUE} (the
#'   default, their full factorial), \code{FALSE} (additive), or a whole number
#'   giving the highest order of interaction among them. It has no effect with
#'   a single grouping variable. The covariate slopes are allowed to differ
#'   across every grouping term this includes. See Details.
#' @param type Character. \code{"III"} (default) or \code{"II"} sums of squares.
#'   The model is fitted under sum-to-zero contrasts when \code{"III"}.
#' @param conf_level Numeric in (0, 1). Level for every interval returned,
#'   including the bands of the covariate plot. Default \code{0.95}.
#' @param vcov_type Character. \code{"model"} (default) or an HC type
#'   (\code{"HC0"} to \code{"HC4"}), which requires the \pkg{sandwich} package.
#'   An HC type reaches \code{$anova} (as Wald F tests), \code{$emmeans},
#'   \code{$posthoc} and \code{$simple_slopes}; see Details.
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
#'
#' @return An \code{\link{anovakit_fit}} object, with five extra components:
#'   \code{$slopes_test} (the homogeneity-of-slopes comparison, also in
#'   \code{$assumptions$slopes}); \code{$simple_slopes} (the slope of each
#'   covariate in each cell of the grouping structure, from
#'   \code{\link[emmeans]{emtrends}}, with intervals at \code{conf_level} and
#'   the covariance chosen by \code{vcov_type}; \code{NULL} when the additive
#'   model is used); \code{$covariate_means} (the mean of each covariate over
#'   the analysed rows, returned whether or not the covariates were centred:
#'   they are the values subtracted when \code{center_covariates = TRUE}, and
#'   the point at which \code{$emmeans} are evaluated either way); and
#'   \code{$model_additive} and \code{$model_interaction}, both candidate fits,
#'   so the one that was not chosen is still available. \code{$assumptions}
#'   also holds \code{normality} (Shapiro-Wilk on the residuals) and
#'   \code{levene} (Levene's test on the residuals across the cells).
#'   \code{$data_used} holds the analysed columns and the cell factor, with the
#'   covariates centred when \code{center_covariates = TRUE}. The covariate
#'   plot shows each covariate on its original scale.
#'
#' @seealso \code{\link{anova_welch}} when there is no covariate,
#'   \code{\link{anova_manova}} for several responses at once.
#'
#' @examples
#' set.seed(7)
#' n <- 150
#' d <- data.frame(
#'   grp = factor(rep(c("control", "treated"), each = n / 2)),
#'   baseline = rnorm(n, mean = 100, sd = 5)
#' )
#' d$score <- 2 * d$baseline + ifelse(d$grp == "treated", 6, 0) + rnorm(n, 0, 3)
#'
#' fit <- anova_ancova(d, "score", "grp", "baseline")
#' fit
#'
#' # Group means adjusted to the covariate mean, not the raw group means
#' fit$emmeans
#' fit$effect_sizes
#'
#' # Why the additive model was chosen, and what centring did
#' fit$slopes_test
#' fit$notes
#'
#' # Fix the model in advance instead of letting the slopes test choose. With
#' # the interaction in the model, centring is what keeps the group row
#' # interpretable: uncentred, it tests the groups at baseline = 0.
#' int <- anova_ancova(d, "score", "grp", "baseline",
#'                     force_interaction = TRUE, plots = FALSE)
#' int$simple_slopes
#' int$anova$p_value[int$anova$term == "grp"]
#' unc <- anova_ancova(d, "score", "grp", "baseline", force_interaction = TRUE,
#'                     center_covariates = FALSE, plots = FALSE)
#' unc$anova$p_value[unc$anova$term == "grp"]
#'
#' # Two grouping variables: their interaction is part of the model
#' d$site <- factor(rep(c("north", "south"), times = n / 2))
#' two <- anova_ancova(d, "score", c("grp", "site"), "baseline", plots = FALSE)
#' two$anova
#'
#' @export
anova_ancova <- function(data, response, groups, covariates,
                         center_covariates = TRUE,
                         force_interaction = NULL,
                         homogeneity_alpha = 0.05,
                         interaction = TRUE,
                         type = c("III", "II"),
                         conf_level = 0.95,
                         vcov_type = "model",
                         adjust = "tukey",
                         posthoc = TRUE,
                         plots = TRUE,
                         verbose = FALSE) {
  cl <- match.call()
  type <- match.arg(type)
  .check_adjust(vcov_type, .vcov_choices(), "vcov_type")
  .check_adjust(adjust, .adjust_choices(), "adjust")
  data <- .as_df(data)
  .check_name(response, "response")
  .check_names(groups, "groups")
  .check_names(covariates, "covariates")
  .check_columns(data, response, "Response column")
  .check_columns(data, groups, "Grouping column(s)")
  .check_columns(data, covariates, "Covariate column(s)")
  .check_roles(list(`the response` = response, `a grouping variable` = groups,
                    `a covariate` = covariates))
  .check_numeric_col(data, response)
  for (cv in covariates) .check_numeric_col(data, cv, role = "Covariate")
  .check_conf_level(conf_level)
  .check_interaction(interaction, length(groups))
  .check_flag(center_covariates, "center_covariates")
  .check_flag(posthoc, "posthoc")
  .check_flag(plots, "plots")
  .check_flag(verbose, "verbose")
  if (!is.null(force_interaction)) .check_flag(force_interaction, "force_interaction")
  if (!is.numeric(homogeneity_alpha) || length(homogeneity_alpha) != 1L ||
      is.na(homogeneity_alpha) || homogeneity_alpha <= 0 ||
      homogeneity_alpha >= 1) {
    .stopf("`homogeneity_alpha` must be a single number strictly between 0 and 1.")
  }

  .say(verbose, "Preparing data.")
  cols <- c(response, groups, covariates)
  prep <- .prepare_frame(data, cols, factors = groups, keep_all = FALSE)
  d <- prep$data
  notes <- prep$notes
  .check_groups(d, groups, min_levels = 2L, min_n = 1L)
  added <- .add_cell(d, groups)
  d <- added$data
  cell <- added$cell

  ## Covariates ----------------------------------------------------------------
  for (cv in covariates) {
    if (length(unique(d[[cv]])) < 2L) {
      .stopf("Covariate `%s` has the same value (%s) in every analysed row, so there is nothing to adjust for. Drop it from `covariates`.",
             cv, format(d[[cv]][1L]))
    }
  }
  covariate_means <- vapply(covariates, function(cv) mean(d[[cv]]), numeric(1))
  covariate_ranges <- vapply(covariates, function(cv) range(d[[cv]]), numeric(2))
  # The plots show the covariates in the units they were measured in; the
  # centred copy is only what the model is fitted on.
  d_raw <- d
  # The note is deferred: whether centring moves the group row at all depends
  # on the model that is chosen further down. It shifts the row only when the
  # covariate-by-group term is in the model AND type = "III"; in the additive
  # model, and under Type II, the group row is bit-identical either way, and
  # telling the user to distrust it would send them away from a correct result.
  if (center_covariates) {
    for (cv in covariates) d[[cv]] <- d[[cv]] - covariate_means[[cv]]
  }

  g_terms <- .group_terms(groups, interaction)
  cv_terms <- .bq(covariates)
  .check_covariates_identified(d, covariates, g_terms)

  ## Homogeneity of slopes ---------------------------------------------------
  # The slopes may differ across every grouping term in the model: under the
  # full factorial that is one slope per cell, not slopes forced to be additive
  # across factors.
  additive_terms <- c(cv_terms, g_terms)
  inter_terms <- c(additive_terms,
                   as.vector(outer(cv_terms, g_terms, paste, sep = ":")))
  fml_add <- stats::reformulate(additive_terms, response = .bq(response))
  fml_int <- stats::reformulate(inter_terms, response = .bq(response))

  fits <- .fit_ancova_pair(d, fml_add, fml_int, groups, covariates, type)
  d <- fits$data
  m_add <- fits$additive
  m_int <- fits$interaction
  notes <- c(notes, fits$notes)

  slopes <- .slopes_comparison(m_add, m_int)
  notes <- c(notes, slopes$note)
  slopes_p <- slopes$p_value
  slopes_test <- data.frame(
    comparison = "additive vs covariate-by-group interaction",
    df = slopes$df,
    statistic = slopes$statistic,
    p_value = slopes_p,
    # NA, not FALSE, when the test could not be computed: asserting
    # heterogeneity is a different claim from having no evidence either way.
    homogeneous = if (is.na(slopes_p)) NA else slopes_p >= homogeneity_alpha,
    stringsAsFactors = FALSE
  )

  use_interaction <- if (!is.null(force_interaction)) {
    force_interaction
  } else if (is.na(slopes_p)) {
    FALSE
  } else {
    !slopes_test$homogeneous
  }

  model <- if (use_interaction) m_int else m_add
  notes <- c(notes, .check_model_size(model))

  # Now that the model is known, say what centring did to this fit. $emmeans
  # are evaluated at the covariate mean whether or not it was centred; only
  # the intercept (and, with the interaction under Type III, the group row)
  # refer to the covariate's zero.
  shifts_group_row <- use_interaction && identical(type, "III")
  notes <- c(notes, if (center_covariates) {
    sprintf("Covariate(s) mean-centred before fitting (%s), so %s to the covariate mean rather than to covariate = 0. $emmeans are evaluated at the covariate mean either way.",
            paste(sprintf("%s: %.4g", covariates, covariate_means),
                  collapse = "; "),
            if (shifts_group_row) "the intercept and the Type III group row refer"
            else "the intercept refers")
  } else {
    outside <- covariates[covariate_ranges[1L, ] > 0 | covariate_ranges[2L, ] < 0]
    sprintf("Covariate(s) were NOT centred, so the intercept refers to covariate = 0; $emmeans and $posthoc are still evaluated at the covariate mean (%s).%s%s Observed ranges: %s.",
            paste(sprintf("%s: %.4g", covariates, covariate_means),
                  collapse = "; "),
            if (shifts_group_row)
              " The Type III group row of the ANOVA table also tests the groups at covariate = 0; read $emmeans instead."
            else
              " The group row of the ANOVA table is unaffected: centring moves the intercept, not that test.",
            if (length(outside) > 0L) sprintf(
              " Zero lies outside the observed range of %s, so the intercept%s describes a point the data never reach.",
              paste(outside, collapse = ", "),
              if (shifts_group_row) " (and that group row)" else "") else "",
            paste(sprintf("%s [%.4g, %.4g]", covariates,
                          covariate_ranges[1L, ], covariate_ranges[2L, ]),
                  collapse = "; "))
  })
  .say(verbose, "Using the %s model.",
       if (use_interaction) "covariate-by-group interaction" else "additive")

  # Exactly one note explains the model choice, on every path, and it is true
  # of the fit it describes: with force_interaction the choice was the
  # caller's, so reporting the slopes test as the reason would state something
  # the test did not find.
  slopes_txt <- if (!is.na(slopes_p)) sprintf("gave p = %.4g", slopes_p) else
    if (slopes$exact) "could not be computed (the fit is exact)" else
      "could not be computed"
  where <- if (!identical(type, "III")) "that ignores the interaction (Type II)" else
    if (center_covariates) "at the covariate mean" else "at covariate = 0"
  read_instead <- sprintf(
    "The group row of the ANOVA table is then a comparison %s, not a constant adjusted difference: read $simple_slopes and $emmeans instead.",
    where)
  notes <- c(notes, if (use_interaction) {
    if (is.null(force_interaction)) {
      sprintf("Slopes differ across groups (p = %.4g), so the covariate-by-group interaction was retained. %s",
              slopes_p, read_instead)
    } else {
      sprintf("The covariate-by-group interaction was fitted because force_interaction = TRUE; the homogeneity-of-slopes test itself %s. %s",
              slopes_txt, read_instead)
    }
  } else if (is.null(force_interaction)) {
    if (is.na(slopes_p) && slopes$exact) {
      "The homogeneity-of-slopes test could not be computed because the model without covariate-by-group terms already fits the data exactly (its residual sum of squares is rounding error), so the additive model was used."
    } else if (is.na(slopes_p)) {
      "The homogeneity-of-slopes test could not be computed, so the additive model was used. Check that the covariate varies within every group."
    } else {
      sprintf("Slopes are consistent with being equal across groups (p = %.4g), so the additive model was used. A non-significant test is not proof of equal slopes, only an absence of evidence against it.",
              slopes_p)
    }
  } else {
    sprintf("The additive model was fitted because force_interaction = FALSE; the homogeneity-of-slopes test %s.",
            slopes_txt)
  })
  if (use_interaction && is.null(force_interaction)) {
    notes <- c(notes, paste0(
      "The model was chosen by a test on these same data, and the p-values below do not allow for that choice: given that the test retained the slope terms, they can be noticeably too small, especially in small samples. ",
      if (identical(type, "III") && center_covariates)
        "Setting force_interaction = TRUE in advance avoids the selection step; the group row is then a valid comparison at the covariate mean whether or not the slopes differ."
      else
        "Set force_interaction in advance to avoid the selection step."))
  }
  if (!use_interaction) {
    notes <- c(notes, .imbalance_note(d_raw, covariates, cell))
  }
  notes <- c(notes, .extrapolation_note(d_raw, covariates, cell,
                                        covariate_means))

  ## Covariance ----------------------------------------------------------------
  rv <- .robust_vcov(model, vcov_type)
  notes <- c(notes, rv$note)
  robust <- !is.null(rv$matrix)

  ## ANOVA table and effect sizes --------------------------------------------
  av <- .car_anova(model, type = type)
  notes <- c(notes, av$note)
  anova_tab <- av$table
  if (robust) {
    rob <- .robust_anova(model, type, rv$matrix, vcov_type,
                         fallback = !is.null(av$table))
    notes <- c(notes, rob$note)
    if (!is.null(rob$table)) anova_tab <- rob$table
  }
  y <- stats::model.response(stats::model.frame(model))
  eff <- .partial_eta_squared(av$raw, df_error = stats::df.residual(model),
                              n_obs = length(y))
  if (!is.null(attr(eff, "omega_floored"))) {
    notes <- c(notes, sprintf(
      "Partial omega squared was negative for %s and has been floored at 0, which happens when an effect explains less than its degrees of freedom would by chance.",
      paste(attr(eff, "omega_floored"), collapse = ", ")))
  }

  ## Assumptions -------------------------------------------------------------
  res <- stats::residuals(model)
  norm <- .check_normality(res)
  notes <- c(notes, norm$note)
  lev <- .levene_on_model(model, d, cell,
                          vcov_type = if (robust) vcov_type else "model")
  notes <- c(notes, lev$note)

  ## Marginal means, comparisons and slopes ----------------------------------
  emm <- .emm_block(model, groups, additive = .is_additive(interaction),
                    type = "link", vcov_matrix = rv$matrix,
                    conf_level = conf_level, adjust = adjust, posthoc = posthoc,
                    link = "identity", data = d)
  notes <- c(notes, emm$notes)

  simple_slopes <- NULL
  if (use_interaction) {
    parts <- lapply(covariates, function(cv) {
      ss <- .emtrends_table(model, groups, cv, conf_level,
                            vcov_matrix = rv$matrix, protect = groups)
      notes <<- c(notes, ss$note)
      if (is.null(ss$table)) NULL else .cbind_label(ss$table, "covariate", cv)
    })
    parts <- parts[!vapply(parts, is.null, logical(1))]
    if (length(parts) > 0L) {
      simple_slopes <- do.call(rbind, parts)
      row.names(simple_slopes) <- NULL
    }
  }

  ## Plots -------------------------------------------------------------------
  plot_list <- list()
  if (plots) {
    .say(verbose, "Building plots.")
    plot_list$residuals <- .plot_resid_fitted(stats::fitted(model), res)
    plot_list$qq <- .plot_qq(res)
    for (cv in covariates) {
      nm <- if (length(covariates) == 1L) "covariate" else paste0("covariate_", cv)
      plot_list[[nm]] <- .plot_covariate(d_raw, response, cv, cell,
                                         conf_level = conf_level)
    }
    if (!is.null(emm$table)) {
      plot_list$emmeans <- .plot_emmeans(
        emm$table, groups, estimate = "estimate",
        lower = "conf_low", upper = "conf_high", conf_level = conf_level,
        title = "Covariate-adjusted marginal means",
        ylab = sprintf("Adjusted %s", response))
    }
  }

  .new_fit(
    method       = sprintf("Analysis of covariance (Type %s)", type),
    call         = cl,
    model        = model,
    anova        = anova_tab,
    effect_sizes = eff,
    emmeans      = emm$table,
    emmeans_object = emm$grid,
    posthoc      = emm$posthoc,
    assumptions  = list(normality = norm$test,
                        levene = lev$table,
                        slopes = slopes_test),
    plots        = plot_list,
    data_used    = d,
    n_removed    = prep$n_removed,
    conf_level   = conf_level,
    notes        = notes,
    extra        = list(slopes_test = slopes_test,
                        simple_slopes = simple_slopes,
                        covariate_means = covariate_means,
                        model_additive = m_add,
                        model_interaction = m_int)
  )
}

#' Refuse covariates whose effect cannot be separated from the groups
#'
#' A covariate that is a function of the grouping variables -- a group-level
#' attribute, or anything constant within every group -- or an exact linear
#' combination of them and the other covariates is aliased in the model. The
#' coefficient \code{lm()} drops is then whichever column comes last, which
#' may be a group contrast rather than the covariate, and the Type III table,
#' the adjusted means and the comparisons all degrade without saying why.
#' @noRd
.check_covariates_identified <- function(d, covariates, g_terms) {
  G <- stats::model.matrix(stats::reformulate(g_terms), data = d)
  C <- as.matrix(d[covariates])
  # The same tolerance lm() uses to declare a column aliased.
  rk <- function(M) qr(M, tol = 1e-7)$rank
  r_g <- rk(G)
  r_all <- rk(cbind(G, C))
  if (r_all == r_g + ncol(C)) return(invisible(TRUE))
  by_groups <- character(0)
  by_others <- character(0)
  for (j in seq_along(covariates)) {
    if (rk(cbind(G, C[, j])) == r_g) {
      by_groups <- c(by_groups, covariates[j])
    } else if (rk(cbind(G, C[, -j, drop = FALSE])) == r_all) {
      by_others <- c(by_others, covariates[j])
    }
  }
  msg <- c(
    if (length(by_groups) > 0L) sprintf(
      "Covariate(s) %s are determined by the grouping variable(s) (each is constant within every group), so their effect cannot be separated from the group effect and there is nothing to adjust for.",
      paste(sprintf("`%s`", by_groups), collapse = ", ")),
    if (length(by_others) > 0L) sprintf(
      "Covariate(s) %s are an exact linear combination of the grouping variable(s) and the other covariates, so their effects cannot be estimated separately.",
      paste(sprintf("`%s`", by_others), collapse = ", ")))
  if (length(msg) == 0L) {
    msg <- "The covariates are collinear with the grouping variables or with each other, so their effects cannot be estimated separately."
  }
  .stopf("%s Drop the offending covariate(s) from `covariates`.",
         paste(msg, collapse = " "))
}

#' Fit both candidate models, keeping their coefficient names distinct
#'
#' A column whose name equals a grouping factor's name plus one of its contrast
#' labels (a covariate \code{dose1} next to a factor \code{dose} fitted under
#' sum-to-zero contrasts, say) gives the model two coefficients of the same
#' name. \pkg{car} indexes terms by position and is unaffected, but
#' \pkg{emmeans} matches the model matrix to the coefficients by name and
#' silently reads the wrong column. When that happens the models are refitted
#' with the grouping factors' contrasts relabelled (\code{[1]}, \code{[2]},
#' ...), which changes no estimate; if even that collides, the call is refused.
#' @return list with \code{additive}, \code{interaction}, \code{data} (the
#'   frame the models were fitted on) and \code{notes}.
#' @noRd
.fit_ancova_pair <- function(d, fml_add, fml_int, groups, covariates, type) {
  fit_one <- function(fml, dd) {
    got <- .collect_conditions(.fit_with_contrasts(
      function() .fit_lm(fml, dd), type = type))
    if (inherits(got$value, "error")) {
      .stopf("The model could not be fitted: %s", conditionMessage(got$value))
    }
    got
  }
  fit_both <- function(dd) {
    a <- fit_one(fml_add, dd)
    i <- fit_one(fml_int, dd)
    list(additive = a$value, interaction = i$value, said = c(a$said, i$said))
  }
  clashes <- function(f) {
    per_model <- list(names(stats::coef(f$additive)),
                      names(stats::coef(f$interaction)))
    unique(unlist(lapply(per_model, function(n) n[duplicated(n)])))
  }
  fits <- fit_both(d)
  dup <- clashes(fits)
  notes <- character(0)
  if (length(dup) > 0L) {
    dup_lab <- gsub("`", "", dup, fixed = TRUE)
    involved <- unique(c(
      covariates[covariates %in% dup_lab],
      groups[vapply(groups, function(g) any(startsWith(dup_lab, g)), logical(1))]))
    d2 <- .bracket_group_contrasts(d, groups, type)
    fits2 <- fit_both(d2)
    if (length(clashes(fits2)) > 0L) {
      .stopf("Column names collide with the names of the model's coefficients (%s is produced twice), so the marginal means would be computed from the wrong columns. Rename the column(s) involved: %s.",
             paste(dup_lab, collapse = ", "),
             paste(sprintf("`%s`", involved), collapse = ", "))
    }
    relabelled <- grep("[", names(stats::coef(fits2$additive)), fixed = TRUE,
                       value = TRUE)
    notes <- sprintf(
      "The names of columns %s collide with the coefficient names of a grouping factor's contrasts (%s would appear twice), which would make emmeans read the wrong model column. The grouping factors' contrast labels were put in brackets in $model (%s, ...) to keep the names distinct; no estimate or test is affected.",
      paste(sprintf("`%s`", involved), collapse = ", "),
      paste(dup_lab, collapse = ", "),
      paste(gsub("`", "", utils::head(relabelled, 2L), fixed = TRUE),
            collapse = ", "))
    d <- d2
    fits <- fits2
  }
  list(additive = fits$additive, interaction = fits$interaction, data = d,
       notes = c(.said_note(unique(fits$said), "stats::lm()"), notes))
}

#' The homogeneity-of-slopes F test, guarded against an exact fit
#'
#' \code{stats::anova()} compares two residual sums of squares with no check
#' that there is any residual variation to compare: on data that lie exactly on
#' parallel lines both are rounding error, their ratio is noise, and a
#' reordering of the rows can flip the conclusion. The test is reported as not
#' computable instead.
#' @noRd
.slopes_comparison <- function(m_add, m_int) {
  rss_add <- stats::deviance(m_add)
  y <- stats::model.response(stats::model.frame(m_add))
  tss <- sum((y - mean(y))^2)
  exact <- is.finite(rss_add) && is.finite(tss) && rss_add <= 1e-12 * tss
  got <- .collect_conditions(stats::anova(m_add, m_int))
  cmp <- got$value
  if (inherits(cmp, "error")) {
    return(list(df = NA_real_, statistic = NA_real_, p_value = NA_real_,
                exact = exact, note = sprintf(
                  "The homogeneity-of-slopes comparison failed: %s",
                  conditionMessage(cmp))))
  }
  out <- list(df = cmp[["Df"]][2L], statistic = cmp[["F"]][2L],
              p_value = cmp[["Pr(>F)"]][2L], exact = exact,
              note = character(0))
  if (exact) {
    out$statistic <- NA_real_
    out$p_value <- NA_real_
  }
  out
}

#' Robust Wald F tests for the ANOVA table
#'
#' \code{car::Anova()} given a covariance through \code{vcov.} tests each term
#' with a Wald F statistic on the model's residual degrees of freedom. The
#' table has no sums of squares, since those are properties of the fit rather
#' than of the covariance.
#' @noRd
.robust_anova <- function(model, type, V, vcov_type, fallback = TRUE) {
  # Through .car_anova(), so a response in small units is rescaled for car's
  # absolute tolerance here as well, with the covariance rescaled to match.
  got <- .car_anova(model, type = type, vcov_matrix = V)
  if (is.null(got$table)) {
    return(list(table = NULL, note = sprintf(
      "Robust (%s) Wald F tests could not be computed%s: %s",
      vcov_type,
      if (fallback) ", so the F tests in $anova are model-based" else "",
      sub("^The Type [I]+ ANOVA table could not be computed: ", "", got$note))))
  }
  tab <- got$table
  attr(tab, "statistic") <- sprintf("Wald F with %s covariance", vcov_type)
  list(table = tab, note = c(
    sprintf("vcov_type = \"%s\": the F tests in $anova are Wald F tests built from that covariance (car::Anova(vcov. = )), on the model's residual degrees of freedom, so the table has no sums of squares. $effect_sizes are still computed from the model-based sums of squares, which the covariance does not change, and the homogeneity-of-slopes test in $slopes_test is the model-based F test.",
            vcov_type),
    got$note))
}

#' A note when the covariate means differ materially between groups
#'
#' With the slopes assumed equal, the adjusted difference between two groups
#' is biased by (slope difference) x (covariate-mean difference). That product
#' vanishes when the covariate is balanced across groups, and can be large
#' when it is not, which is when a slopes test with little power matters.
#' @noRd
.imbalance_note <- function(d, covariates, cell, threshold = 0.5) {
  g <- droplevels(as.factor(d[[cell]]))
  gap <- vapply(covariates, function(cv) {
    x <- d[[cv]]
    dfw <- length(x) - nlevels(g)
    if (dfw < 1L) return(NA_real_)
    sp <- sqrt(sum((x - stats::ave(x, g))^2) / dfw)
    if (!is.finite(sp) || sp <= 0) return(NA_real_)
    m <- tapply(x, g, mean)
    (max(m) - min(m)) / sp
  }, numeric(1))
  big <- !is.na(gap) & gap > threshold
  if (!any(big)) return(character(0))
  sprintf("The covariate means differ between groups (largest difference: %s pooled within-group SDs). With the slopes assumed equal, any real difference in slopes biases the adjusted comparisons by the slope difference times the covariate-mean difference, and the homogeneity-of-slopes test often lacks the power to detect a difference large enough to matter, so the group row and $posthoc can reject far more often than their nominal level. force_interaction = TRUE does not assume equal slopes.",
          paste(sprintf("%s %.2g", covariates[big], gap[big]), collapse = ", "))
}

#' A note when the covariate mean lies outside a group's observed range
#' @noRd
.extrapolation_note <- function(d, covariates, cell, ref) {
  g <- droplevels(as.factor(d[[cell]]))
  bits <- character(0)
  for (cv in covariates) {
    lo <- tapply(d[[cv]], g, min)
    hi <- tapply(d[[cv]], g, max)
    out <- which(ref[[cv]] < lo | ref[[cv]] > hi)
    if (length(out) > 0L) {
      bits <- c(bits, sprintf(
        "%s = %.4g lies outside the range observed in %s", cv, ref[[cv]],
        .abbrev(sprintf("%s [%.4g, %.4g]", names(lo)[out], lo[out], hi[out]))))
    }
  }
  if (length(bits) == 0L) return(character(0))
  sprintf("$emmeans are evaluated at the covariate mean, but %s. The adjusted means there are extrapolations along the fitted slope to covariate values none of those observations have, so they and the comparisons that involve them rest on the model's linearity rather than on data.",
          paste(bits, collapse = "; "))
}

#' Levene's test on the residuals of a fitted model
#'
#' The residuals are taken from the model's own frame and the grouping factor is
#' rebuilt from the same rows, so the two can never fall out of alignment when
#' the model has dropped incomplete cases.
#'
#' Levene's statistic is an ANOVA of the absolute deviations from the group
#' medians. When those do not vary within any group -- every group has at most
#' two observations, whose deviations from their median are equal -- that
#' ANOVA has no residual variation and its F is a ratio of rounding errors
#' (about 1e31), so the test is skipped.
#' @param vcov_type the covariance used for the other results, which decides
#'   what a rejection note recommends.
#' @noRd
.levene_on_model <- function(model, d, cell, vcov_type = "model") {
  res <- stats::residuals(model)
  rows <- as.integer(names(res))
  if (anyNA(rows) || length(rows) != length(res)) {
    rows <- seq_len(nrow(d))[stats::complete.cases(stats::model.frame(model))]
  }
  g <- droplevels(as.factor(d[[cell]][rows]))
  if (length(g) != length(res)) {
    return(list(table = NULL, note =
      "Levene's test was skipped: the residuals and the grouping factor could not be aligned."))
  }
  if (nlevels(g) < 2L) {
    return(list(table = NULL, note =
      "Levene's test was skipped: fewer than two populated groups."))
  }
  r <- as.numeric(res)
  z <- abs(r - stats::ave(r, g, FUN = stats::median))
  ss_within <- sum((z - stats::ave(z, g))^2)
  ss_z <- sum(z^2)
  if (!is.finite(ss_within) || !is.finite(ss_z) || ss_z == 0 ||
      ss_within <= 1e-10 * ss_z) {
    return(list(table = NULL, note =
      "Levene's test was skipped: the absolute deviations from the group medians do not vary within any group (as happens when no group has more than two observations), so its F statistic is undefined."))
  }
  tmp <- data.frame(res = r, grp = g)
  got <- .collect_conditions(car::leveneTest(res ~ grp, data = tmp))
  out <- got$value
  if (inherits(out, "error")) {
    return(list(table = NULL, note = sprintf(
      "Levene's test failed: %s", conditionMessage(out))))
  }
  tab <- as.data.frame(out)
  res_tab <- data.frame(
    df1 = tab[["Df"]][1L], df2 = tab[["Df"]][2L],
    statistic = tab[["F value"]][1L], p_value = tab[["Pr(>F)"]][1L],
    stringsAsFactors = FALSE)
  note <- .said_note(got$said, "car::leveneTest()")
  if (!is.na(res_tab$p_value) && res_tab$p_value < 0.05) {
    note <- c(note, if (identical(vcov_type, "model")) {
      sprintf("Levene's test rejects equality of residual variances across groups (p = %.3g). The F tests in $anova and the intervals and tests in $emmeans, $posthoc and $simple_slopes all use the pooled residual variance; vcov_type = \"HC3\" bases them on a heteroscedasticity-consistent covariance instead (the homogeneity-of-slopes test stays model-based).",
              res_tab$p_value)
    } else {
      sprintf("Levene's test rejects equality of residual variances across groups (p = %.3g); the %s covariance requested through vcov_type already allows for this in $anova, $emmeans, $posthoc and $simple_slopes (the homogeneity-of-slopes test is model-based).",
              res_tab$p_value, vcov_type)
    })
  }
  list(table = res_tab, note = note)
}

#' Covariate slopes in each cell of the grouping structure
#'
#' Evaluated like the marginal means: global \code{emm_options()} are ignored,
#' every other covariate is held at its mean (\code{cov.keep = character(0)}),
#' and the robust covariance is used when one was requested.
#'
#' Only the trend column emmeans produces is renamed to \code{slope}, located
#' by position after the grid factors, so a grouping column called
#' \code{slope} (or \code{<covariate>.trend}) keeps its labels and the slopes
#' take the name \code{emm_slope} instead, as other result tables do.
#' @noRd
.emtrends_table <- function(model, groups, covariate, conf_level,
                            vcov_matrix = NULL, protect = character(0)) {
  old <- .emm_options_guard()
  on.exit(options(emmeans = old), add = TRUE)
  args <- list(object = model, specs = .formula(NULL, groups), var = covariate,
               cov.keep = character(0))
  if (!is.null(vcov_matrix)) args$vcov. <- vcov_matrix
  # Captured, not suppressed: emtrends() reaches summary.lm(), which warns on a
  # near-perfect fit, and that would otherwise reach the console.
  got <- .collect_conditions(do.call(emmeans::emtrends, args))
  out <- got$value
  if (inherits(out, "error")) {
    return(list(table = NULL, note = sprintf(
      "Per-group covariate slopes could not be computed: %s",
      conditionMessage(out))))
  }
  said <- .said_note(got$said, "emmeans::emtrends()")
  s <- tryCatch(
    suppressWarnings(summary(out, infer = c(TRUE, TRUE), level = conf_level)),
    error = function(e) e)
  if (inherits(s, "error")) {
    return(list(table = NULL, note = sprintf(
      "Per-group covariate slopes could not be summarised: %s",
      conditionMessage(s))))
  }
  pri <- attr(s, "pri.vars")
  est <- attr(s, "estName")
  df <- as.data.frame(s)
  nm <- names(df)
  n_pri <- length(pri)
  pos <- n_pri + match(est, nm[-seq_len(n_pri)])
  if (is.na(pos)) {
    return(list(table = NULL, note =
      "Per-group covariate slopes could not be summarised: emmeans returned no trend column."))
  }
  # A placeholder no other column can have, so nothing is indexed by a name
  # that may be duplicated.
  placeholder <- .safe_name(df, ".slope")
  names(df)[pos] <- placeholder
  df <- .normalise_emm_names(df, protect = protect)
  renamed <- attr(df, "renamed")
  target <- "slope"
  if (target %in% setdiff(names(df), placeholder)) {
    k <- 1L
    target <- "emm_slope"
    while (target %in% setdiff(names(df), placeholder)) {
      k <- k + 1L
      target <- paste0("emm_slope", k)
    }
    renamed <- c(renamed, sprintf("slope -> %s", target))
  }
  names(df)[names(df) == placeholder] <- target
  attr(df, "renamed") <- renamed
  # A cell with no data has no estimable slope.
  miss <- is.na(df[[target]])
  note <- c(said, .renamed_note(df))
  if (any(miss)) {
    lab <- do.call(paste, c(lapply(df[miss, intersect(groups, names(df)), drop = FALSE],
                                   as.character), sep = " : "))
    note <- c(note, sprintf(
      "The slope of %s is not estimable in %d cell(s), which are left out of $simple_slopes: %s.",
      covariate, sum(miss), .abbrev(lab)))
    df <- df[!miss, , drop = FALSE]
    row.names(df) <- NULL
    attr(df, "renamed") <- renamed
  }
  list(table = df, note = note)
}
