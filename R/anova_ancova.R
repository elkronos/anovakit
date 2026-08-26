#' Analysis of Covariance
#'
#' Compares the mean of a numeric response across groups while adjusting for one
#' or more numeric covariates. Tests the homogeneity-of-slopes assumption,
#' checks residual normality and equality of variances, reports Type II or Type
#' III sums of squares with partial eta squared, and returns covariate-adjusted
#' estimated marginal means.
#'
#' @details
#' \strong{Covariates are centred.} Every covariate is mean-centred before the
#' model is fitted. This matters most when the covariate-by-group interaction is
#' in the model and \code{type = "III"}: that Type III row tests the group
#' difference \emph{at covariate = 0}, which for an uncentred covariate may be a
#' point nowhere near the data, and the p-value can be anything. On the example
#' below -- a baseline distributed around 100 -- the same data gives p = 0.35
#' uncentred against p = 5.5e-21 centred. In the additive model the group row is
#' the same either
#' way; what centring buys there is that the intercept and \code{$emmeans} are
#' read at the covariate mean rather than at zero. Set
#' \code{center_covariates = FALSE} only if the covariate's own zero is the
#' point you want.
#'
#' \strong{Homogeneity of slopes.} The interaction between each covariate and
#' the grouping factor is tested against the additive model. If the slopes
#' differ, the interaction model is used and a note explains that the group
#' effect is now a comparison at one point on the covariate rather than a
#' constant difference, and that \code{$simple_slopes} is the thing to read.
#' Because the model is selected from a test on the same data, the subsequent
#' p-values are mildly optimistic; use \code{force_interaction} to fix the model
#' in advance and avoid this. \code{$notes} always reports which model was
#' fitted and why.
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
#'   interaction; \code{FALSE} always fits the additive model.
#' @param homogeneity_alpha Numeric in (0, 1). Significance level for the slopes
#'   test. Default \code{0.05}.
#' @param type Character. \code{"III"} (default) or \code{"II"} sums of squares.
#'   The model is fitted under sum-to-zero contrasts when \code{"III"}.
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
#'
#' @return An \code{\link{anovakit_fit}} object, with five extra components:
#'   \code{$slopes_test} (the homogeneity-of-slopes comparison, also in
#'   \code{$assumptions$slopes}), \code{$simple_slopes} (per-group covariate
#'   slopes from \code{\link[emmeans]{emtrends}}, \code{NULL} when the additive
#'   model is used), \code{$covariate_means} (the means subtracted when
#'   centring, so a slope can be read back at the original scale), and
#'   \code{$model_additive} and \code{$model_interaction}, both candidate fits,
#'   so the one that was not chosen is still available. \code{$assumptions} also
#'   holds \code{normality} (Shapiro-Wilk on the residuals) and \code{levene}.
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
#' @export
anova_ancova <- function(data, response, groups, covariates,
                         center_covariates = TRUE,
                         force_interaction = NULL,
                         homogeneity_alpha = 0.05,
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
  .check_numeric_col(data, response)
  for (cv in covariates) .check_numeric_col(data, cv, role = "Covariate")
  .check_conf_level(conf_level)
  .check_flag(center_covariates, "center_covariates")
  .check_flag(posthoc, "posthoc")
  .check_flag(plots, "plots")
  .check_flag(verbose, "verbose")
  if (!is.null(force_interaction)) .check_flag(force_interaction, "force_interaction")
  if (!is.numeric(homogeneity_alpha) || length(homogeneity_alpha) != 1L ||
      homogeneity_alpha <= 0 || homogeneity_alpha >= 1) {
    .stopf("`homogeneity_alpha` must be a single number strictly between 0 and 1.")
  }
  overlap <- intersect(groups, covariates)
  if (length(overlap) > 0L) {
    .stopf("A column cannot be both a group and a covariate: %s.",
           paste(overlap, collapse = ", "))
  }
  if (response %in% c(groups, covariates)) {
    .stopf("The response `%s` cannot also be a %s.", response,
           if (response %in% groups) "grouping variable" else "covariate")
  }

  .say(verbose, "Preparing data.")
  cols <- c(response, groups, covariates)
  prep <- .prepare_frame(data, cols, factors = groups, keep_all = TRUE)
  d <- prep$data
  notes <- prep$notes
  .check_groups(d, groups, min_levels = 2L, min_n = 1L)

  ## Centre the covariates ---------------------------------------------------
  # The note is deferred: whether centring moves the group row at all depends
  # on the model that is chosen further down. It shifts the row only when the
  # covariate-by-group term is in the model AND type = "III"; in the additive
  # model, and under Type II, the group row is bit-identical either way, and
  # telling the user to distrust it would send them away from a correct result.
  covariate_means <- vapply(covariates, function(cv) mean(d[[cv]]), numeric(1))
  covariate_ranges <- vapply(covariates, function(cv) range(d[[cv]]), numeric(2))
  if (center_covariates) {
    for (cv in covariates) d[[cv]] <- d[[cv]] - covariate_means[[cv]]
  }

  added <- .add_cell(d, groups)
  d <- added$data
  cell <- added$cell

  ## Homogeneity of slopes ---------------------------------------------------
  g_terms <- .bq(groups)
  cv_terms <- .bq(covariates)
  additive_terms <- c(cv_terms, g_terms)
  inter_terms <- c(additive_terms,
                   as.vector(outer(cv_terms, g_terms, paste, sep = ":")))
  fml_add <- stats::reformulate(additive_terms, response = .bq(response))
  fml_int <- stats::reformulate(inter_terms, response = .bq(response))

  m_add <- .fit_with_contrasts(function() stats::lm(fml_add, data = d), type)
  m_int <- .fit_with_contrasts(function() stats::lm(fml_int, data = d), type)
  cmp <- stats::anova(m_add, m_int)
  slopes_p <- cmp[["Pr(>F)"]][2L]
  slopes_test <- data.frame(
    comparison = "additive vs covariate-by-group interaction",
    df = cmp[["Df"]][2L],
    statistic = cmp[["F"]][2L],
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

  # Now that the model is known, say what centring did to this fit.
  shifts_group_row <- use_interaction && identical(type, "III")
  notes <- c(notes, if (center_covariates) {
    sprintf("Covariate(s) mean-centred before fitting (%s), so $emmeans and the intercept are read at the covariate mean%s.",
            paste(sprintf("%s: %.4g", covariates, covariate_means),
                  collapse = "; "),
            if (shifts_group_row)
              ", and the Type III group row is a comparison there rather than at covariate = 0"
            else "")
  } else {
    outside <- covariates[covariate_ranges[1L, ] > 0 | covariate_ranges[2L, ] < 0]
    sprintf("Covariate(s) were NOT centred, so $emmeans and the intercept are read at covariate = 0.%s%s Observed ranges: %s.",
            if (shifts_group_row)
              " The group row of the ANOVA table tests the groups at that point too; read $emmeans instead."
            else
              " The group row of the ANOVA table is unaffected: centring moves the intercept, not that test.",
            if (length(outside) > 0L) sprintf(
              " Zero lies outside the observed range of %s, so anything read there describes a point the data never reach.",
              paste(outside, collapse = ", ")) else "",
            paste(sprintf("%s [%.4g, %.4g]", covariates,
                          covariate_ranges[1L, ], covariate_ranges[2L, ]),
                  collapse = "; "))
  })
  .say(verbose, "Using the %s model.",
       if (use_interaction) "covariate-by-group interaction" else "additive")

  # Say why the model was chosen, and say it truthfully. With
  # force_interaction the choice was the caller's, so reporting the slopes
  # test as the reason states something the test did not find.
  # Exactly one note explains the model choice, on every path, and it is true
  # of the fit it describes: with force_interaction the choice was the
  # caller's, so reporting the slopes test as the reason would state something
  # the test did not find.
  slopes_txt <- if (is.na(slopes_p)) "could not be computed" else
    sprintf("gave p = %.4g", slopes_p)
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
    if (is.na(slopes_p)) {
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
    notes <- c(notes, "The model was chosen by a test on these same data, so the p-values below are mildly optimistic. Set force_interaction to fix the model in advance.")
  }

  ## ANOVA table and effect sizes --------------------------------------------
  av <- .car_anova(model, type = type)
  notes <- c(notes, av$note)
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
  lev <- .levene_on_model(model, d, cell)
  notes <- c(notes, lev$note)

  ## Marginal means, comparisons and slopes ----------------------------------
  rv <- .robust_vcov(model, vcov_type)
  notes <- c(notes, rv$note)
  emm <- .emmeans_grid(model, groups, type = "link", vcov_matrix = rv$matrix)
  notes <- c(notes, emm$note)
  emm_tab <- .emmeans_table(emm$grid, conf_level, protect = groups)
  notes <- c(notes, emm_tab$note)
  ph <- if (posthoc) {
    .emmeans_pairs(emm$grid, adjust = adjust, conf_level = conf_level)
  } else list(table = NULL,
              note = .no_posthoc_note(NROW(emm_tab$table)))
  notes <- c(notes, ph$note)

  simple_slopes <- NULL
  if (use_interaction) {
    parts <- lapply(covariates, function(cv) {
      ss <- .emtrends_table(model, groups, cv, conf_level, protect = groups)
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
      plot_list[[nm]] <- .plot_covariate(d, response, cv, cell)
    }
    if (!is.null(emm_tab$table)) {
      plot_list$emmeans <- .plot_emmeans(
        emm_tab$table, groups, estimate = "estimate",
        lower = "conf_low", upper = "conf_high", conf_level = conf_level,
        title = "Covariate-adjusted marginal means",
        ylab = sprintf("Adjusted %s", response))
    }
  }

  .new_fit(
    method       = sprintf("Analysis of covariance (Type %s)", type),
    call         = cl,
    model        = model,
    anova        = av$table,
    effect_sizes = eff,
    emmeans      = emm_tab$table,
    emmeans_object = emm$grid,
    posthoc      = ph$table,
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

#' Levene's test on the residuals of a fitted model
#'
#' The residuals are taken from the model's own frame and the grouping factor is
#' rebuilt from the same rows, so the two can never fall out of alignment when
#' the model has dropped incomplete cases.
#' @noRd
.levene_on_model <- function(model, d, cell) {
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
  tmp <- data.frame(res = as.numeric(res), grp = g)
  out <- tryCatch(car::leveneTest(res ~ grp, data = tmp), error = function(e) e)
  if (inherits(out, "error")) {
    return(list(table = NULL, note = sprintf(
      "Levene's test failed: %s", conditionMessage(out))))
  }
  tab <- as.data.frame(out)
  res_tab <- data.frame(
    df1 = tab[["Df"]][1L], df2 = tab[["Df"]][2L],
    statistic = tab[["F value"]][1L], p_value = tab[["Pr(>F)"]][1L],
    stringsAsFactors = FALSE)
  note <- if (!is.na(res_tab$p_value) && res_tab$p_value < 0.05) {
    "Levene's test rejects equality of residual variances across groups. The ANCOVA F tests assume equal variances; consider vcov_type = \"HC3\"."
  } else character(0)
  list(table = res_tab, note = note)
}

#' Per-group covariate slopes
#' @noRd
.emtrends_table <- function(model, groups, covariate, conf_level,
                            protect = character(0)) {
  # Captured, not suppressed: emtrends() reaches summary.lm(), which warns on a
  # near-perfect fit, and that would otherwise reach the console.
  got <- .collect_conditions(emmeans::emtrends(
    model, specs = .formula(NULL, groups), var = covariate))
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
  df <- as.data.frame(s)
  names(df)[grepl("\\.trend$", names(df))] <- "slope"
  tab <- .normalise_emm_names(df, protect = protect)
  list(table = tab, note = c(said, .renamed_note(tab)))
}
