#' Repeated Measures Analysis of Variance
#'
#' Fits a repeated measures ANOVA to data in long format using
#' \code{\link[afex]{aov_ez}}, reports Mauchly's test of sphericity with the
#' Greenhouse-Geisser and Huynh-Feldt corrections, generalised and partial eta
#' squared, residual diagnostics, and estimated marginal means.
#'
#' @details
#' \strong{What is fixed rather than inherited.} The sums of squares are Type
#' III unless \code{type = 2} is given through \code{...}; the between-subject
#' factors are coded with sum-to-zero contrasts, which Type III needs; and the
#' marginal means and comparisons use \pkg{afex}'s multivariate model, whose
#' standard errors and degrees of freedom come from each within-subject cell's
#' own variance rather than a pooled error term. None of these is read from
#' \code{afex::afex_options()}, so a global setting made for some other
#' analysis cannot change the result. The type is recorded in
#' \code{attr(fit$anova, "ss_type")} (and printed), the emmeans model in
#' \code{attr(fit$emmeans, "emmeans_model")} and
#' \code{attr(fit$posthoc, "emmeans_model")}.
#'
#' \strong{Sphericity.} Mauchly's test and both epsilon estimates are read
#' from the multivariate model \pkg{afex} fitted (through
#' \code{car::summary.Anova.mlm()}). Each term's error matrix is first
#' rescaled to unit average variance, so a response measured in small units
#' does not fall below \pkg{car}'s \emph{absolute} singularity tolerance and
#' lose its corrections. When the error matrix is genuinely singular -- fewer
#' subjects than the term has contrasts -- \pkg{car} cannot give the
#' corrections and Mauchly's test is undefined. The Greenhouse-Geisser
#' epsilon, \eqn{tr(S)^2 / ((k-1) tr(S^2))}, is still defined, where \eqn{S}
#' is the covariance matrix of the orthonormal within-subject contrasts; it
#' and the Huynh-Feldt epsilon (with Lecoutre's correction) are then computed
#' directly, the requested correction is applied with them, and
#' \code{$notes} says so. The Huynh-Feldt estimate can exceed 1; it is used,
#' and reported in \code{$sphericity$hf_epsilon}, capped at 1, with the
#' uncapped value in \code{$sphericity$hf_epsilon_raw}. The corrected degrees
#' of freedom are the uncorrected ones multiplied by epsilon, as in
#' \pkg{afex}.
#'
#' \strong{Unbalanced subjects.} A repeated measures ANOVA needs every subject
#' to appear in every within-subject cell. Subjects who do not are removed
#' before fitting, and both the count and the identifiers (in the order of the
#' subject column's levels) are reported in \code{$subjects_dropped} and in
#' \code{$notes}.
#'
#' \strong{Repeated rows.} A subject with more than one row in the same
#' within-subject cell has those rows combined into one value before fitting,
#' with \code{fun_aggregate} (the mean unless another function is given). The
#' note says how many subject-by-cell combinations were affected and how many
#' extra rows were combined away; those rows are counted in \code{$n_removed}.
#'
#' \strong{Column names and labels.} \pkg{afex} builds a model formula from
#' the column names and wide-format column names from the within-subject
#' levels, so a column name that is not a syntactic R name (\code{"subject
#' id"}), or level labels such as \code{0, 4, 12}, would be mangled or parsed
#' as code. \pkg{afex} is therefore given internal names (\code{Y}, \code{ID},
#' \code{W1}, \code{W2}, ..., \code{B1}, ...) and level labels (\code{L01},
#' \code{L02}, ... in level order), and every table, \code{$data_used},
#' \code{$emmeans_object} and the plots are mapped back to your names and
#' labels. \code{$model} is the \pkg{afex} fit itself and keeps the internal
#' names; \code{$internal_names} gives the correspondence.
#'
#' \strong{Normality.} The F tests for within-subject effects depend on the
#' within-subject errors: what is left of each value after its subject's mean
#' and its cell's mean (within each between-subject group) are removed,
#' \eqn{y_{ij} - \bar{y}_{i\cdot} - \bar{y}_{g(i)\cdot j} +
#' \bar{y}_{g(i)\cdot\cdot}}{y[ij] - mean(subject i) - mean(cell j in group
#' g(i)) + mean(group g(i))}. These are \code{$residuals}, and Shapiro-Wilk
#' is run on them separately in each within-subject cell
#' (\code{$assumptions$normality}), as \code{\link{anova_welch}} tests each
#' group: the cells' variances usually differ, and pooling values of
#' different spread gives a mixture that fails a normality test even when
#' every cell is normal. With a single two-level within-subject factor the
#' residuals are plus and minus half of each subject's centred difference, so
#' both rows test the normality of the differences, which is exactly the
#' assumption of the equivalent paired t-test. The Q-Q plot shows the same
#' residuals standardised within each cell. Normality of the subject means,
#' which the between-subject tests rely on, is not tested.
#'
#' \strong{A factor response.} With \code{factorize = TRUE}, a response stored
#' as a factor is converted with \code{as.numeric(as.character(x))}, which
#' recovers the numbers. A bare \code{as.numeric()} on a factor returns the
#' level \emph{indices}, so a response of 11, 17, 21 would silently become
#' 1, 2, 3.
#'
#' Requires the \pkg{afex} package.
#'
#' @param data A data frame in long format: one row per subject per cell.
#' @param response Character. Name of the numeric response column.
#' @param subject Character. Name of the subject identifier column.
#' @param within Character vector. One or more within-subject factors.
#' @param between Character vector or \code{NULL}. Between-subject factors,
#'   each constant within a subject.
#' @param factorize Logical. Convert the response to numeric through its labels.
#'   Default \code{TRUE}; with \code{FALSE} a non-numeric response is an error.
#'   The subject, within- and between-subject columns are always coerced to
#'   factors, whatever this is set to.
#' @param emm_specs Character vector or \code{NULL}. Factors to compute
#'   estimated marginal means over. Defaults to all within-subject factors;
#'   when it leaves out a factor of the model, the means are averaged over
#'   that factor's levels and \code{$notes} says so.
#' @param adjust Character. Multiplicity adjustment for the pairwise
#'   comparisons, applied by \pkg{emmeans}. One of \code{"tukey"},
#'   \code{"sidak"}, \code{"scheffe"}, \code{"dunnettx"}, \code{"bonferroni"},
#'   \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"BH"},
#'   \code{"BY"}, \code{"fdr"} or \code{"none"}. Default \code{"tukey"}.
#' @param conf_level Numeric in (0, 1). Level for every interval returned.
#'   Default \code{0.95}.
#' @param correction Character. Sphericity correction applied to the ANOVA
#'   table: \code{"GG"} (default), \code{"HF"} or \code{"none"}. It is
#'   recorded in \code{attr(fit$anova, "correction")} and printed.
#' @param posthoc Logical. Compute pairwise comparisons. There are
#'   \code{choose(k, 2)} of them, so this is worth turning off when the number
#'   of cells is large; \code{$notes} records that they were skipped. Default
#'   \code{TRUE}.
#' @param plots Logical. Build \pkg{ggplot2} objects. They are returned in
#'   \code{$plots}, never drawn. Default \code{TRUE}.
#' @param verbose Logical. Emit progress through \code{\link[base]{message}}.
#'   Default \code{FALSE}.
#' @param ... Only these arguments of \code{\link[afex]{aov_ez}} are
#'   accepted; anything else is an error rather than silently ignored.
#'   \describe{
#'     \item{\code{fun_aggregate}}{A function combining the rows of a subject
#'       that share a within-subject cell into one number. Default
#'       \code{mean}.}
#'     \item{\code{observed}}{Character. Factors that were observed (measured)
#'       rather than manipulated. Generalised eta squared depends on this
#'       (Olejnik and Algina, 2003).}
#'     \item{\code{type}}{Sums of squares: \code{3} (or \code{"III"}, the
#'       default) or \code{2} (\code{"II"}).}
#'     \item{\code{anova_table}}{A list holding only \code{p_adjust_method},
#'       one of \code{\link[stats]{p.adjust.methods}}, which adjusts the
#'       p-values of \code{$anova} across its terms. The \code{p_gg} and
#'       \code{p_hf} columns of \code{$sphericity} stay unadjusted.}
#'   }
#'   \code{covariate} is refused: \code{anova_rm()} does not fit covariates.
#'   \code{id}, \code{dv}, \code{return}, \code{include_aov} and
#'   \code{print.formula} are set by the function itself, and
#'   \code{transformation} is refused (transform the response column
#'   instead).
#'
#' @return An \code{\link{anovakit_fit}} object. \code{$anova} has the
#'   corrected degrees of freedom, the MSE, F, partial eta squared and the
#'   p-value for each term; \code{$effect_sizes} has partial and generalised
#'   eta squared, the latter honouring \code{observed}. The function adds
#'   \describe{
#'     \item{\code{$sphericity}}{Mauchly's test (\code{mauchly_w},
#'       \code{p_value}, \code{NA} where it is undefined), the
#'       Greenhouse-Geisser and Huynh-Feldt epsilons (\code{gg_epsilon},
#'       \code{hf_epsilon}, the latter capped at 1), the p-values under each
#'       correction (\code{p_gg}, \code{p_hf}) and the uncapped Huynh-Feldt
#'       estimate (\code{hf_epsilon_raw}), one row per within-subject term
#'       with more than one degree of freedom.}
#'     \item{\code{$residuals}}{The within-subject residuals (see Details),
#'       one per row of \code{$data_used} and in the same order, so
#'       \code{$residuals[i]} belongs to \code{$data_used[i, ]}.}
#'     \item{\code{$subjects_dropped}}{The subjects removed for incomplete
#'       within-subject cells.}
#'     \item{\code{$internal_names}}{A data frame with the \code{role},
#'       your \code{name} and the \code{internal} name each column has in
#'       \code{$model}.}
#'     \item{\code{$n_removed_missing}, \code{$n_removed_unbalanced},
#'       \code{$n_removed_aggregated}}{\code{$n_removed} broken into its three
#'       causes.}
#'   }
#'   \code{$data_used} has one row per subject and within-subject cell, in the
#'   order the rows came in (a combined row takes the place of its first
#'   row). \code{$assumptions$normality} has one Shapiro-Wilk test per
#'   within-subject cell, with a \code{note} giving the reason when a cell
#'   could not be tested.
#'
#' @references Olejnik, S., & Algina, J. (2003). Generalized eta and omega
#'   squared statistics: measures of effect size for some common research
#'   designs. \emph{Psychological Methods}, 8(4), 434-447.
#'
#'   Lecoutre, B. (1991). A correction for the epsilon-tilde approximate test
#'   in repeated measures designs with two or more independent groups.
#'   \emph{Journal of Educational Statistics}, 16(4), 371-372.
#'
#' @seealso \code{\link{anova_welch}} for independent groups,
#'   \code{\link[afex]{aov_ez}} for the fit itself.
#'
#' @examples
#' if (requireNamespace("afex", quietly = TRUE)) {
#'   set.seed(1)
#'   d <- expand.grid(id = factor(1:24), time = factor(c("t1", "t2", "t3")))
#'   d$arm <- factor(rep(rep(c("ctrl", "trt"), each = 12), 3))
#'   d$score <- 10 + 2 * as.numeric(d$time) +
#'     1.5 * (d$arm == "trt") + rnorm(nrow(d), 0, 2)
#'
#'   fit <- anova_rm(d, "score", subject = "id",
#'                   within = "time", between = "arm")
#'   print(fit)
#'
#'   # Mauchly's test with both epsilon corrections, read off the fitted model
#'   print(fit$sphericity)
#'
#'   # Generalised eta squared alongside partial eta squared, and the marginal
#'   # means over the within-subject factor (averaged over arm, as $notes says)
#'   print(fit$effect_sizes)
#'   print(fit$emmeans)
#'
#'   # A factor that was measured rather than assigned (an age group, say) is
#'   # declared with `observed`, which changes generalised eta squared
#'   d$age_grp <- factor(rep(rep(c("younger", "older"), 12), 3))
#'   anova_rm(d, "score", subject = "id", within = "time",
#'            between = c("arm", "age_grp"), observed = "age_grp",
#'            plots = FALSE)$effect_sizes
#' }
#'
#' @export
anova_rm <- function(data, response, subject, within, between = NULL,
                     factorize = TRUE,
                     emm_specs = NULL,
                     adjust = "tukey",
                     conf_level = 0.95,
                     correction = c("GG", "HF", "none"),
                     posthoc = TRUE,
                     plots = TRUE,
                     verbose = FALSE,
                     ...) {
  cl <- match.call()
  correction <- match.arg(correction)
  .check_adjust(adjust, .adjust_choices(), "adjust")
  if (!requireNamespace("afex", quietly = TRUE)) {
    .stopf("anova_rm() requires the {afex} package. Install it with install.packages(\"afex\").")
  }
  data <- .as_df(data)
  .check_name(response, "response")
  .check_name(subject, "subject")
  .check_names(within, "within")
  if (!is.null(between)) .check_names(between, "between")
  .check_roles(list(`the response` = response,
                    `the subject identifier` = subject,
                    `a within-subject factor` = within,
                    `a between-subject factor` = between))
  .check_columns(data, response, "Response column")
  .check_columns(data, subject, "Subject column")
  .check_columns(data, within, "Within-subject factor(s)")
  if (!is.null(between)) .check_columns(data, between, "Between-subject factor(s)")
  .check_conf_level(conf_level)
  .check_flag(factorize, "factorize")
  .check_flag(posthoc, "posthoc")
  .check_flag(plots, "plots")
  .check_flag(verbose, "verbose")
  opts <- .rm_dots(list(...), as.list(substitute(list(...)))[-1L],
                   within, between)
  specs <- emm_specs %||% within
  .check_names(specs, "emm_specs")
  unknown <- setdiff(specs, c(within, between))
  if (length(unknown) > 0L) {
    .stopf("`emm_specs` must name within- or between-subject factors; %s did not.",
           paste(sprintf("`%s`", unknown), collapse = ", "))
  }

  .say(verbose, "Preparing data.")
  notes <- character(0)
  if (factorize) {
    conv <- .coerce_response_numeric(data[[response]], response)
    data[[response]] <- conv$values
    notes <- c(notes, conv$notes)
  } else {
    .check_numeric_col(data, response)
  }

  cols <- c(response, subject, within, between)
  # A subject identifier stored as numbers is a label, not a covariate, so the
  # "many distinct values" warning meant for grouping columns does not apply.
  prep <- .prepare_frame(data, cols,
                         factors = c(subject, within, between),
                         keep_all = FALSE, quiet_factors = subject)
  d <- prep$data
  notes <- c(notes, prep$notes)

  for (w in within) {
    if (nlevels(d[[w]]) < 2L) {
      .stopf("Within-subject factor `%s` has %d level(s); at least 2 are required.",
             w, nlevels(d[[w]]))
    }
  }
  .rm_check_design(d, subject, within, between)

  ## Balance the design ------------------------------------------------------
  bal <- .balance_subjects(d, subject, within)
  d <- bal$data
  if (length(bal$dropped) > 0L) {
    notes <- c(notes, sprintf(
      "%d subject(s) were removed because they are missing at least one within-subject cell: %s%s. A repeated measures ANOVA requires a complete design.",
      length(bal$dropped),
      paste(utils::head(bal$dropped, 10L), collapse = ", "),
      if (length(bal$dropped) > 10L) ", ..." else ""))
  }
  if (nlevels(d[[subject]]) < 2L) {
    .stopf("Fewer than two subjects remain with a complete set of within-subject cells.")
  }
  for (w in within) d[[w]] <- droplevels(d[[w]])
  for (b in between) {
    d[[b]] <- droplevels(d[[b]])
    if (nlevels(d[[b]]) < 2L) {
      .stopf("Between-subject factor `%s` has %d level(s) among the subjects with a complete set of within-subject cells; at least 2 are required.",
             b, nlevels(d[[b]]))
    }
  }
  n_after_balance <- nrow(d)

  ## Combine repeated subject-by-cell rows -----------------------------------
  agg <- .rm_aggregate(d, response, subject, within, opts$fun_aggregate)
  d <- agg$data
  if (agg$n_extra > 0L) {
    notes <- c(notes, sprintf(
      "%d subject-by-cell combination(s) have more than one row (%d extra row(s)); each was aggregated into one value with %s before fitting. $data_used holds the aggregated rows, and the %d row(s) aggregated away are included in $n_removed.",
      agg$n_combinations, agg$n_extra,
      opts$fun_label %||% "the mean (the default; pass fun_aggregate to choose another function)",
      agg$n_extra))
  }
  data_used <- d

  ## Fit ---------------------------------------------------------------------
  .say(verbose, "Fitting the repeated measures model with afex::aov_ez().")
  map <- .rm_internal_frame(d, response, subject, within, between)
  obs_int <- if (length(opts$observed) > 0L) unname(map$to_internal[opts$observed]) else NULL
  afex_args <- list(id = "ID", dv = "Y", data = map$data,
                    within = unname(map$to_internal[within]),
                    between = if (length(between) > 0L) unname(map$to_internal[between]) else NULL,
                    observed = obs_int,
                    type = if (opts$type == "II") 2L else 3L,
                    factorize = FALSE, return = "afex_aov", include_aov = FALSE)
  if (!is.null(opts$p_adjust_method)) {
    afex_args$anova_table <- list(p_adjust_method = opts$p_adjust_method)
  }
  got <- .collect_conditions(do.call(afex::aov_ez, afex_args))
  fit <- got$value
  if (inherits(fit, "error")) {
    .stopf("afex::aov_ez() failed: %s",
           .rm_unmap_text(conditionMessage(fit), map$to_user))
  }
  notes <- c(notes, .said_note(.rm_unmap_text(got$said, map$to_user),
                               "afex::aov_ez()"))

  ## ANOVA table and sphericity ----------------------------------------------
  base <- .rm_base_table(fit, obs_int, map$to_user)
  notes <- c(notes, attr(base, "note"))
  sph <- .rm_sphericity(fit)
  notes <- c(notes, .rm_unmap_text(sph$said, map$to_user))
  sph_tab <- NULL
  if (!is.null(sph$table)) {
    sph_tab <- sph$table
    idx <- match(sph_tab$term, base$term_int)
    f_stat <- base$statistic[idx]
    df1 <- base$num_df[idx]
    df2 <- base$den_df[idx]
    sph_tab$hf_epsilon <- pmin(1, sph_tab$hf_epsilon_raw)
    sph_tab$p_gg <- .rm_pf(f_stat, df1 * sph_tab$gg_epsilon,
                           df2 * sph_tab$gg_epsilon)
    sph_tab$p_hf <- .rm_pf(f_stat, df1 * sph_tab$hf_epsilon,
                           df2 * sph_tab$hf_epsilon)
    sph_tab$term_user <- base$term[idx]
  }
  notes <- c(notes, .rm_sphericity_notes(sph, sph_tab, correction))

  anova_tab <- .rm_corrected_table(base, sph_tab, correction)
  attr(anova_tab, "ss_type") <- opts$type
  if (!is.null(opts$p_adjust_method) && opts$p_adjust_method != "none") {
    anova_tab$p_value <- stats::p.adjust(anova_tab$p_value,
                                         method = opts$p_adjust_method)
    attr(anova_tab, "p_adjust_method") <- opts$p_adjust_method
    notes <- c(notes, sprintf(
      "The p-values in $anova are adjusted across its %d terms with the %s method (anova_table = list(p_adjust_method = \"%s\")); the p_gg and p_hf columns of $sphericity are not adjusted.",
      nrow(anova_tab), opts$p_adjust_method, opts$p_adjust_method))
  }
  notes <- c(notes, .rm_degenerate_notes(base))

  eff <- data.frame(term = base$term,
                    partial_eta_sq = base$partial_eta_sq,
                    generalised_eta_sq = base$generalised_eta_sq,
                    stringsAsFactors = FALSE, row.names = NULL)
  if (length(opts$observed) > 0L) attr(eff, "observed") <- opts$observed

  if (!is.null(sph_tab)) {
    bad <- sph_tab$term_user[!is.na(sph_tab$p_value) & sph_tab$p_value < 0.05]
    if (length(bad) > 0L) {
      notes <- c(notes, if (correction == "none") {
        sprintf("Mauchly's test rejects sphericity for %s; the reported table is uncorrected (correction = \"none\"), so its p-values for %s may be too small.",
                paste(bad, collapse = ", "),
                if (length(bad) == 1L) "that term" else "those terms")
      } else {
        sprintf("Mauchly's test rejects sphericity for %s; the reported table uses the %s correction.",
                paste(bad, collapse = ", "),
                if (correction == "GG") "Greenhouse-Geisser" else "Huynh-Feldt")
      })
    }
    sph_tab <- data.frame(term = sph_tab$term_user,
                          mauchly_w = sph_tab$mauchly_w,
                          p_value = sph_tab$p_value,
                          gg_epsilon = sph_tab$gg_epsilon,
                          hf_epsilon = sph_tab$hf_epsilon,
                          p_gg = sph_tab$p_gg,
                          p_hf = sph_tab$p_hf,
                          hf_epsilon_raw = sph_tab$hf_epsilon_raw,
                          stringsAsFactors = FALSE, row.names = NULL)
  }

  ## Residual diagnostics ----------------------------------------------------
  res <- .rm_within_residuals(data_used, response, subject, within, between)
  norm <- .normality_by_group(res$residuals, res$cell)
  names(norm)[names(norm) == "group"] <- "cell"
  skipped <- is.na(norm$p_value)
  if (any(skipped)) {
    notes <- c(notes, sprintf(
      "Shapiro-Wilk skipped for %d of %d within-subject cell(s); see $assumptions$normality$note (%s).",
      sum(skipped), nrow(norm),
      paste(unique(sub("\\.$", "", sub("^Shapiro-Wilk skipped: ", "",
                                        norm$note[skipped]))),
            collapse = " / ")))
  }

  ## Marginal means ----------------------------------------------------------
  # The emmeans model is passed explicitly: left to afex_options("emmeans_model")
  # a session setting could switch the comparisons to the pooled-error
  # univariate model, with different standard errors and degrees of freedom.
  emm <- .emmeans_grid(fit, unname(map$to_internal[specs]), type = "link",
                       model = "multivariate")
  notes <- c(notes, .rm_unmap_text(emm$note, map$to_user))
  grid <- emm$grid
  if (!is.null(grid)) {
    rl <- .rm_relabel_grid(grid, map)
    grid <- rl$grid
    notes <- c(notes, rl$note)
  }
  emm_tab <- .emmeans_table(grid, conf_level, protect = specs)
  notes <- c(notes, emm_tab$note)
  notes <- c(notes, .rm_averaged_note(emm_tab$mesg, specs, within, between,
                                      anova_tab, base$term_int, map))
  n_cells <- NROW(emm_tab$table)
  n_comp <- if (n_cells >= 2L) choose(n_cells, 2L) else 0
  ph <- if (!posthoc) {
    list(table = NULL, note = .no_posthoc_note(n_cells))
  } else if (n_comp > 5000) {
    list(table = NULL, note = .too_many_note(n_comp, 5000L))
  } else if (n_cells >= 2L) {
    .emmeans_pairs(grid, adjust = adjust, conf_level = conf_level)
  } else list(table = NULL, note = character(0))
  notes <- c(notes, ph$note)
  emm_out <- emm_tab$table
  if (!is.null(emm_out)) attr(emm_out, "emmeans_model") <- "multivariate"
  ph_out <- ph$table
  if (!is.null(ph_out)) attr(ph_out, "emmeans_model") <- "multivariate"

  ## Plots -------------------------------------------------------------------
  plot_list <- list()
  if (plots) {
    .say(verbose, "Building plots.")
    plot_list$residuals <- .plot_resid_fitted(
      res$fitted, res$residuals,
      title = "Residuals vs fitted (within-subject)")
    plot_list$qq <- .plot_qq(
      res$standardised,
      title = "Normal Q-Q plot of within-subject residuals",
      subtitle = "Standardised within each within-subject cell, as tested in $assumptions$normality")
    plot_list$index <- .plot_resid_index(
      res$residuals,
      title = "Within-subject residuals in the row order of $data_used")
    if (!is.null(emm_out)) {
      plot_list$emmeans <- .plot_emmeans(
        emm_out, specs, estimate = "estimate",
        lower = "conf_low", upper = "conf_high", conf_level = conf_level,
        title = "Estimated marginal means",
        ylab = sprintf("Estimated %s", response))
    }
  }

  anova_tab$term_int <- NULL
  n_input <- nrow(data)
  n_unbalanced <- n_input - prep$n_removed - n_after_balance
  .new_fit(
    method       = "Repeated measures analysis of variance",
    call         = cl,
    model        = fit,
    anova        = anova_tab,
    effect_sizes = eff,
    emmeans      = emm_out,
    emmeans_object = grid,
    posthoc      = ph_out,
    assumptions  = list(sphericity = sph_tab, normality = norm),
    plots        = plot_list,
    data_used    = data_used,
    n_removed    = n_input - nrow(data_used),
    conf_level   = conf_level,
    notes        = notes,
    extra        = list(subjects_dropped = bal$dropped,
                        n_removed_missing = prep$n_removed,
                        n_removed_unbalanced = as.integer(n_unbalanced),
                        n_removed_aggregated = agg$n_extra,
                        sphericity = sph_tab,
                        residuals = res$residuals,
                        internal_names = map$table)
  )
}

#' Validate the arguments passed through `...` to afex::aov_ez()
#'
#' afex's own \code{...} only reaches \code{fun_aggregate}, so anything it does
#' not name is accepted and ignored. Only the arguments the wrapper can honour
#' pass; everything else is refused by name.
#' @param dots \code{list(...)}.
#' @param exprs the unevaluated expressions of \code{...}, used to name the
#'   aggregation function in the note.
#' @return list with \code{fun_aggregate}, \code{fun_label}, \code{observed},
#'   \code{type} (\code{"II"} or \code{"III"}) and \code{p_adjust_method}.
#' @noRd
.rm_dots <- function(dots, exprs, within, between) {
  out <- list(fun_aggregate = NULL, fun_label = NULL, observed = NULL,
              type = "III", p_adjust_method = NULL)
  if (length(dots) == 0L) return(out)
  nms <- names(dots) %||% rep("", length(dots))
  nms[is.na(nms)] <- ""
  if (any(!nzchar(nms))) {
    .stopf("Every argument passed through `...` must be named; anova_rm() received %d unnamed argument(s).",
           sum(!nzchar(nms)))
  }
  if (anyDuplicated(nms)) {
    .stopf("Argument(s) given more than once: %s.",
           paste(sprintf("`%s`", unique(nms[duplicated(nms)])), collapse = ", "))
  }
  allowed <- c("fun_aggregate", "observed", "type", "anova_table")
  set_here <- c(id = "subject", dv = "response", return = NA, include_aov = NA,
                print.formula = NA)
  unsupported <- c("covariate", "transformation")
  unknown <- setdiff(nms, c(allowed, names(set_here), unsupported))
  if (length(unknown) > 0L) {
    .stopf("anova_rm() has no argument %s, so it would be ignored. The arguments passed on to afex::aov_ez() through `...` may only be: %s.",
           paste(sprintf("`%s`", unknown), collapse = ", "),
           paste(allowed, collapse = ", "))
  }
  clash <- intersect(nms, names(set_here))
  if (length(clash) > 0L) {
    how <- vapply(clash, function(a) {
      if (is.na(set_here[[a]])) sprintf("`%s` is set by anova_rm() itself", a) else
        sprintf("`%s` is set from `%s`", a, set_here[[a]])
    }, character(1))
    .stopf("%s; do not pass %s through `...`.", paste(how, collapse = "; "),
           if (length(clash) == 1L) "it" else "them")
  }
  if ("covariate" %in% nms) {
    .stopf("anova_rm() does not fit covariates (`covariate`). A between-subject covariate changes what the within-subject tests mean unless it is centred, and it needs checks this function does not make; fit afex::aov_ez() directly, or analyse the subject-level summaries with anova_ancova().")
  }
  if ("transformation" %in% nms) {
    .stopf("`transformation` is not supported; transform the response column before calling anova_rm(), so that $data_used and the marginal means are on the scale you analysed.")
  }

  if ("fun_aggregate" %in% nms) {
    f <- dots$fun_aggregate
    if (!is.function(f)) {
      .stopf("`fun_aggregate` must be a function, such as mean or median.")
    }
    out$fun_aggregate <- f
    e <- exprs[["fun_aggregate"]]
    out$fun_label <- if (is.symbol(e) ||
                         (is.call(e) && as.character(e[[1L]])[1L] %in% c("::", ":::"))) {
      sprintf("`%s`", paste(deparse(e), collapse = ""))
    } else {
      "the function given as fun_aggregate"
    }
  }
  if ("observed" %in% nms) {
    o <- dots$observed
    if (!is.null(o)) {
      if (!is.character(o) || anyNA(o) || length(o) == 0L) {
        .stopf("`observed` must be a character vector naming within- or between-subject factors.")
      }
      bad <- setdiff(o, c(within, between))
      if (length(bad) > 0L) {
        .stopf("`observed` must name within- or between-subject factors; %s did not.",
               paste(sprintf("`%s`", bad), collapse = ", "))
      }
      out$observed <- unique(o)
    }
  }
  if ("type" %in% nms) {
    t <- dots$type
    ok <- length(t) == 1L && !is.na(t) && as.character(t) %in% c("2", "3", "II", "III")
    if (!ok) {
      .stopf("`type` must be 3 (or \"III\", the default) or 2 (or \"II\").")
    }
    out$type <- if (as.character(t) %in% c("2", "II")) "II" else "III"
  }
  if ("anova_table" %in% nms) {
    at <- dots$anova_table
    if (!is.list(at)) {
      .stopf("`anova_table` must be a list, such as list(p_adjust_method = \"holm\").")
    }
    at_names <- names(at) %||% rep("", length(at))
    at_names[is.na(at_names)] <- ""
    extra <- at_names[at_names != "p_adjust_method"]
    if (length(extra) > 0L) {
      .stopf("`anova_table` may hold only `p_adjust_method`; %s cannot be honoured. Use the `correction` argument for the sphericity correction and `observed` for observed factors; partial and generalised eta squared and the MSE are always reported.",
             paste(unique(ifelse(nzchar(extra), sprintf("`%s`", extra),
                                 "an unnamed element")), collapse = ", "))
    }
    m <- at$p_adjust_method
    if (!is.null(m)) {
      if (!is.character(m) || length(m) != 1L || is.na(m) ||
          !m %in% stats::p.adjust.methods) {
        .stopf("`anova_table$p_adjust_method` must be one of: %s.",
               paste(stats::p.adjust.methods, collapse = ", "))
      }
      out$p_adjust_method <- m
    }
  }
  out
}

#' Coerce a response to numeric without turning a factor into level indices
#' @noRd
.coerce_response_numeric <- function(x, name) {
  notes <- character(0)
  if (is.numeric(x)) return(list(values = x, notes = notes))
  if (is.factor(x) || is.character(x)) {
    chr <- as.character(x)
    num <- suppressWarnings(as.numeric(chr))
    n_bad <- sum(is.na(num) & !is.na(chr))
    if (n_bad > 0L) {
      bad <- unique(chr[is.na(num) & !is.na(chr)])
      .stopf("Response `%s` is %s and %d value(s) cannot be read as numbers (for example: %s).",
             name, if (is.factor(x)) "a factor" else "character", n_bad,
             paste(utils::head(bad, 3L), collapse = ", "))
    }
    notes <- c(notes, sprintf(
      "Response `%s` was stored as %s and was converted to numeric via its labels, not its level codes.",
      name, if (is.factor(x)) "a factor" else "character"))
    return(list(values = num, notes = notes))
  }
  if (is.logical(x)) return(list(values = as.numeric(x), notes = notes))
  .stopf("Response `%s` must be numeric; it is %s.",
         name, paste(class(x), collapse = "/"))
}

#' Check that the design can be fitted as a repeated measures design
#'
#' A between-subject factor must be constant within each subject, and the
#' within-subject factors must be fully crossed. afex refuses both, but in its
#' own terms and with the internal column names.
#' @noRd
.rm_check_design <- function(d, subject, within, between) {
  for (b in between) {
    n_val <- tapply(as.integer(d[[b]]), d[[subject]],
                    function(z) length(unique(z)))
    bad <- names(n_val)[!is.na(n_val) & n_val > 1L]
    if (length(bad) > 0L) {
      .stopf("Between-subject factor `%s` must take one value per subject, but it varies within %d subject(s): %s.",
             b, length(bad), .abbrev(bad))
    }
  }
  if (length(within) > 1L) {
    codes <- do.call(paste, c(lapply(within, function(w) as.integer(d[[w]])),
                              list(sep = "\r")))
    n_obs <- length(unique(codes))
    n_all <- prod(vapply(within, function(w) nlevels(d[[w]]), integer(1)))
    if (n_obs < n_all) {
      .stopf("The within-subject factors %s are not fully crossed: %d of their %d level combinations never occur. A repeated measures ANOVA needs every combination.",
             paste(sprintf("`%s`", within), collapse = ", "), n_all - n_obs, n_all)
    }
  }
  invisible(TRUE)
}

#' Drop subjects that do not appear in every within-subject cell
#'
#' The dropped subjects are returned in the order of the subject column's
#' levels, so which ones a note names does not depend on the row order.
#' @noRd
.balance_subjects <- function(d, subject, within) {
  # A separator that cannot appear in a level. Pasting with sep = "" would
  # make subject "a" in cell "bc" indistinguishable from subject "ab" in
  # cell "c", miscounting the duplicates the note reports. Written as an
  # escape rather than a literal control byte, which is invisible in source.
  cell <- .cell_vector(d, within, sep = "\001")
  n_cells <- nlevels(cell)
  subj <- as.factor(d[[subject]])
  complete_by_subject <- tapply(as.character(cell), subj,
                                function(z) length(unique(z)))
  complete_by_subject[is.na(complete_by_subject)] <- 0L
  keep_ids <- names(complete_by_subject)[complete_by_subject == n_cells]
  present <- levels(droplevels(subj))
  dropped <- present[!present %in% keep_ids]
  out <- d[as.character(subj) %in% keep_ids, , drop = FALSE]
  out[[subject]] <- droplevels(as.factor(out[[subject]]))
  row.names(out) <- NULL
  key <- paste(out[[subject]], .cell_vector(out, within, sep = "\001"),
               sep = "\001")
  list(data = out, dropped = dropped, duplicates = sum(duplicated(key)))
}

#' Combine a subject's repeated rows within a within-subject cell
#'
#' @param fun the aggregation function, or \code{NULL} for the mean.
#' @return list with \code{data} (one row per subject and cell, in the order
#'   of each combination's first row), \code{n_combinations} (combinations
#'   that had more than one row) and \code{n_extra} (rows combined away).
#' @noRd
.rm_aggregate <- function(d, response, subject, within, fun = NULL) {
  key <- do.call(paste, c(lapply(c(subject, within), function(v) as.integer(d[[v]])),
                          list(sep = "\r")))
  first <- !duplicated(key)
  n_extra <- sum(!first)
  if (n_extra == 0L) {
    return(list(data = d, n_combinations = 0L, n_extra = 0L))
  }
  fun <- fun %||% mean
  counts <- table(key)
  groups <- split(d[[response]], factor(key, levels = unique(key)))
  vals <- lapply(groups, function(v) {
    r <- tryCatch(fun(v), error = function(e) e)
    if (inherits(r, "error")) {
      .stopf("`fun_aggregate` failed on a subject's repeated rows: %s",
             conditionMessage(r))
    }
    if (!is.numeric(r) && !is.logical(r) || length(r) != 1L || !is.finite(r)) {
      .stopf("`fun_aggregate` must return one finite number for each subject-by-cell combination; it returned %s.",
             if (length(r) != 1L) sprintf("%d values", length(r)) else
               paste(format(r), collapse = ""))
    }
    as.numeric(r)
  })
  out <- d[first, , drop = FALSE]
  out[[response]] <- unlist(vals, use.names = FALSE)[match(key[first], names(groups))]
  row.names(out) <- NULL
  list(data = out, n_combinations = sum(counts > 1L), n_extra = n_extra)
}

#' The frame afex is given: internal column names and level labels
#'
#' afex pastes column names into a formula it parses, and runs within-subject
#' levels through make.names() to name its wide columns, so a non-syntactic
#' name fails (or is evaluated as code) and numeric levels come back as X0,
#' X4, X12. Every column gets a fixed syntactic name and every level a label
#' L01, L02, ... in level order.
#' @return list with \code{data}, \code{to_internal} and \code{to_user} (named
#'   character maps), \code{levels} (for each internal factor name, the user's
#'   labels in internal order) and \code{table} (for \code{$internal_names}).
#' @noRd
.rm_internal_frame <- function(d, response, subject, within, between) {
  w_int <- sprintf("W%d", seq_along(within))
  b_int <- if (length(between) > 0L) sprintf("B%d", seq_along(between)) else character(0)
  user <- c(response, subject, within, between)
  int <- c("Y", "ID", w_int, b_int)
  recode <- function(x, prefix) {
    x <- droplevels(as.factor(x))
    n <- nlevels(x)
    lab <- sprintf("%s%0*d", prefix, nchar(as.character(n)), seq_len(n))
    f <- factor(as.integer(x), levels = seq_len(n), labels = lab)
    list(f = f, user = levels(x))
  }
  di <- data.frame(Y = as.numeric(d[[response]]))
  levs <- list()
  s <- recode(d[[subject]], "S")
  di$ID <- s$f
  for (i in seq_along(within)) {
    r <- recode(d[[within[i]]], "L")
    di[[w_int[i]]] <- r$f
    levs[[w_int[i]]] <- r$user
  }
  for (i in seq_along(between)) {
    r <- recode(d[[between[i]]], "L")
    di[[b_int[i]]] <- r$f
    levs[[b_int[i]]] <- r$user
  }
  # Sum-to-zero contrasts, which Type III tests need. afex sets them too, but
  # only while afex_options("check_contrasts") is TRUE.
  for (v in c(w_int, b_int)) stats::contrasts(di[[v]]) <- "contr.sum"
  role <- c("response", "subject", rep("within", length(within)),
            rep("between", length(between)))
  list(data = di,
       to_internal = stats::setNames(int, user),
       to_user = stats::setNames(user, int),
       levels = levs,
       table = data.frame(role = role, name = user, internal = int,
                          stringsAsFactors = FALSE))
}

#' Map afex's term labels (or any text) back to the user's names
#' @noRd
.rm_unmap_terms <- function(terms, to_user) {
  vapply(strsplit(as.character(terms), ":", fixed = TRUE), function(p) {
    hit <- p %in% names(to_user)
    p[hit] <- unname(to_user[p[hit]])
    paste(p, collapse = ":")
  }, character(1))
}

#' In one pass, so that a user's column called, say, "B1" is not translated
#' a second time, and without treating the user's names as regular
#' expressions.
#' @noRd
.rm_unmap_text <- function(x, to_user) {
  if (length(x) == 0L) return(x)
  pat <- sprintf("\\b(%s)\\b", paste(names(to_user), collapse = "|"))
  m <- gregexpr(pat, x, perl = TRUE)
  regmatches(x, m) <- lapply(regmatches(x, m),
                             function(k) unname(to_user[k]))
  x
}

#' The uncorrected ANOVA table with both effect sizes
#'
#' The sphericity correction is applied separately (see
#' \code{.rm_corrected_table()}), from epsilons that are available even where
#' car cannot supply them. \code{observed} is passed on: afex's anova() method
#' does not read the value stored on the fit, so leaving it out would report
#' generalised eta squared as though every factor were manipulated.
#' @noRd
.rm_base_table <- function(fit, observed, to_user) {
  got <- .collect_conditions(as.data.frame(stats::anova(
    fit, correction = "none", es = c("pes", "ges"), observed = observed,
    MSE = TRUE, intercept = FALSE, p_adjust_method = "none")))
  tab <- got$value
  if (inherits(tab, "error")) {
    .stopf("The ANOVA table could not be computed from the afex fit: %s",
           .rm_unmap_text(conditionMessage(tab), to_user))
  }
  out <- data.frame(term = .rm_unmap_terms(rownames(tab), to_user),
                    num_df = tab[["num Df"]],
                    den_df = tab[["den Df"]],
                    mse = tab[["MSE"]],
                    statistic = tab[["F"]],
                    partial_eta_sq = tab[["pes"]],
                    generalised_eta_sq = tab[["ges"]],
                    p_value = tab[["Pr(>F)"]],
                    term_int = rownames(tab),
                    stringsAsFactors = FALSE, row.names = NULL)
  attr(out, "note") <- .said_note(.rm_unmap_text(got$said, to_user),
                                  "afex's anova() method")
  out
}

#' Mauchly's test and the epsilon estimates for each within-subject term
#'
#' Read from car's summary of the multivariate model, after each term's error
#' matrix is rescaled to unit average variance: car judges singularity with an
#' absolute tolerance, so a response in small units would otherwise lose its
#' corrections although nothing is singular. Every statistic car derives here
#' is invariant to that rescaling, which is by a power of two and so exact.
#' Where the error matrix is genuinely singular (fewer error degrees of
#' freedom than contrasts), the epsilons are computed directly from it.
#'
#' @return list with \code{table} (internal term names; \code{mauchly_w},
#'   \code{p_value}, \code{gg_epsilon}, \code{hf_epsilon_raw} and
#'   \code{direct}, TRUE where car could not supply the epsilons),
#'   \code{reason} when there is no table, \code{error_df} and \code{said}.
#' @noRd
.rm_sphericity <- function(fit) {
  A <- fit$Anova
  if (!inherits(A, "Anova.mlm") || !isTRUE(A$repeated)) {
    return(list(table = NULL, reason = "unavailable", said = character(0)))
  }
  err_df <- A$error.df
  p <- vapply(A$P, ncol, integer(1))
  terms <- A$terms
  multi <- terms[p >= 2L]
  if (length(multi) == 0L) {
    return(list(table = NULL, reason = "two_levels", said = character(0)))
  }
  if (!is.finite(err_df) || err_df < 1) {
    return(list(table = NULL, reason = "no_error_df", said = character(0)))
  }

  contrast_cov <- function(t) {
    A$SSPE[[t]] %*% solve(crossprod(A$P[[t]]))
  }
  scaled <- A
  for (t in terms) {
    M <- contrast_cov(t)
    tr <- sum(diag(M))
    k2 <- if (is.finite(tr) && tr > 0) 2^(-round(log2(tr / ncol(M)))) else 1
    scaled$SSP[[t]] <- A$SSP[[t]] * k2
    scaled$SSPE[[t]] <- A$SSPE[[t]] * k2
    if (!is.null(scaled$singular)) {
      ev <- eigen(scaled$SSPE[[t]], symmetric = TRUE, only.values = TRUE)$values
      scaled$singular[[t]] <- !all(is.finite(ev)) ||
        sum(ev >= sqrt(.Machine$double.eps)) < ncol(scaled$SSPE[[t]])
    }
  }
  got <- .collect_conditions(summary(scaled, multivariate = FALSE))
  s <- if (inherits(got$value, "error")) NULL else got$value
  known <- "error SSP matrix|HF eps > 1 treated as 1"
  said <- .said_note(got$said[!grepl(known, got$said)], "car::Anova()")

  sp <- if (!is.null(s$sphericity.tests)) as.matrix(unclass(s$sphericity.tests)) else NULL
  eps <- if (!is.null(s$pval.adjustments)) as.matrix(unclass(s$pval.adjustments)) else NULL
  car_value <- function(m, t, col) {
    if (is.null(m) || !t %in% rownames(m) || !col %in% colnames(m)) return(NA_real_)
    as.numeric(m[t, col])
  }

  rows <- lapply(multi, function(t) {
    gg <- car_value(eps, t, "GG eps")
    hf <- car_value(eps, t, "HF eps")
    direct <- is.na(gg)
    zero <- FALSE
    if (direct) {
      # Greenhouse-Geisser from the covariance of the (orthonormalised)
      # contrasts, tr(S)^2 / (p tr(S^2)), which stays defined when S is
      # singular; Huynh-Feldt from it with Lecoutre's correction.
      M <- contrast_cov(t)
      pk <- ncol(M)
      tr1 <- sum(diag(M))
      tr2 <- sum(diag(M %*% M))
      # No within-subject error variance at all: nothing to correct, and the
      # degenerate-design note explains why no test is possible.
      zero <- !is.finite(tr1) || tr1 <= 0
      gg <- if (is.finite(tr1) && is.finite(tr2) && tr2 > 0) tr1^2 / (pk * tr2) else NA_real_
      den <- pk * (err_df - pk * gg)
      hf <- if (!is.na(gg) && is.finite(den) && den > 0) {
        ((err_df + 1) * pk * gg - 2) / den
      } else NA_real_
      if (!is.na(hf) && hf <= 0) hf <- NA_real_
    }
    data.frame(term = t,
               mauchly_w = car_value(sp, t, "Test statistic"),
               p_value = car_value(sp, t, "p-value"),
               gg_epsilon = gg, hf_epsilon_raw = hf,
               direct = direct & !zero, zero = zero, k = unname(p[t]),
               stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  row.names(tab) <- NULL
  list(table = tab, reason = NULL, error_df = err_df, said = said)
}

#' Notes describing how sphericity was assessed
#' @noRd
.rm_sphericity_notes <- function(sph, tab, correction) {
  if (is.null(tab)) {
    return(switch(sph$reason %||% "unavailable",
      two_levels = "Sphericity is not an issue for this design: every within-subject factor has two levels, where the assumption holds automatically.",
      no_error_df = "Sphericity could not be assessed: the design leaves no error degrees of freedom.",
      "Sphericity could not be assessed for this design."))
  }
  notes <- character(0)
  if (any(tab$direct)) {
    dir <- tab[tab$direct, , drop = FALSE]
    why <- if (all(sph$error_df < dir$k)) {
      sprintf("there are fewer error degrees of freedom (%d) than contrasts (%s)",
              as.integer(sph$error_df), paste(unique(dir$k), collapse = "/"))
    } else {
      "the within-subject contrasts are linearly dependent in these data"
    }
    notes <- c(notes, sprintf(
      "car could not estimate the sphericity corrections for %s: the error matrix of the within-subject contrasts is singular because %s, so Mauchly's test is undefined and reported as NA. The Greenhouse-Geisser and Huynh-Feldt epsilons were computed directly from the covariance matrix of the orthonormal within-subject contrasts, and %s.",
      paste(dir$term_user, collapse = ", "), why,
      if (correction == "none") "no correction was applied (correction = \"none\")" else
        sprintf("the %s correction in $anova uses them wherever they are defined",
                if (correction == "GG") "Greenhouse-Geisser" else "Huynh-Feldt")))
  }
  capped <- !is.na(tab$hf_epsilon_raw) & tab$hf_epsilon_raw > 1
  if (any(capped)) {
    notes <- c(notes, sprintf(
      "The Huynh-Feldt epsilon estimate exceeds 1 for %s (%s). Epsilon cannot exceed 1, so it is capped at 1 in $sphericity$hf_epsilon and wherever it is used (p_hf%s); the uncapped estimate is in $sphericity$hf_epsilon_raw.",
      paste(tab$term_user[capped], collapse = ", "),
      paste(format(signif(tab$hf_epsilon_raw[capped], 4)), collapse = ", "),
      if (correction == "HF") ", and the Huynh-Feldt correction in $anova" else ""))
  }
  undefined <- is.na(if (correction == "HF") tab$hf_epsilon_raw else tab$gg_epsilon) &
    !tab$zero
  if (correction != "none" && any(undefined)) {
    notes <- c(notes, sprintf(
      "The %s epsilon is undefined for %s with this few subjects, so the corrected p-value is reported as NA rather than uncorrected.",
      if (correction == "GG") "Greenhouse-Geisser" else "Huynh-Feldt",
      paste(tab$term_user[undefined], collapse = ", ")))
  }
  notes
}

#' Apply the sphericity correction to the uncorrected table
#'
#' As afex does: both degrees of freedom are multiplied by epsilon, the MSE is
#' the error sum of squares over the corrected error degrees of freedom, and
#' the p-value comes from the F distribution on the corrected degrees of
#' freedom. A term whose epsilon cannot be computed keeps its degrees of
#' freedom and gets an NA p-value, never an uncorrected one under a corrected
#' label.
#' @noRd
.rm_corrected_table <- function(base, sph_tab, correction) {
  out <- base[c("term", "num_df", "den_df", "mse", "statistic",
                "partial_eta_sq", "p_value", "term_int")]
  if (correction != "none" && !is.null(sph_tab)) {
    eps <- if (correction == "GG") sph_tab$gg_epsilon else sph_tab$hf_epsilon
    idx <- match(sph_tab$term, base$term_int)
    ok <- !is.na(eps)
    i <- idx[ok]
    e <- eps[ok]
    out$num_df[i] <- base$num_df[i] * e
    out$den_df[i] <- base$den_df[i] * e
    out$mse[i] <- base$mse[i] / e
    out$p_value[i] <- .rm_pf(base$statistic[i], out$num_df[i], out$den_df[i])
    out$p_value[idx[!ok]] <- NA_real_
  }
  row.names(out) <- NULL
  attr(out, "correction") <- correction
  attr(out, "statistic") <- "F"
  out
}

#' Notes for a design with no error degrees of freedom or no error variance
#' @noRd
.rm_degenerate_notes <- function(base) {
  notes <- character(0)
  no_df <- is.finite(base$den_df) & base$den_df < 1
  if (any(no_df)) {
    notes <- c(notes, sprintf(
      "No error degrees of freedom remain for %s: every between-subject cell holds a single subject, so the effects cannot be separated from the subjects' own variation. No F test is possible: F, the p-values and the MSE are NaN, partial and generalised eta squared are 1 by construction, and the marginal means have no standard errors. More subjects per between-subject cell are needed.",
      paste(base$term[no_df], collapse = ", ")))
  }
  zero <- !no_df & (is.na(base$statistic) | !is.finite(base$statistic) |
                      (!is.na(base$mse) & base$mse == 0))
  if (any(zero)) {
    notes <- c(notes, sprintf(
      "The error variance is zero for %s: every subject shows exactly the same pattern across the cells, so the F statistic is infinite or undefined and no test is possible. The pairwise comparisons have zero standard errors for the same reason.",
      paste(base$term[zero], collapse = ", ")))
  }
  notes
}

#' Within-subject residuals and fitted values, in the row order of the frame
#'
#' Each value minus its subject's mean and its cell's mean within its
#' between-subject group, plus that group's grand mean: the residuals of the
#' model the within-subject F tests are built on, which do not carry the
#' subjects' own level.
#' @return list with \code{residuals}, \code{fitted}, \code{standardised}
#'   (residuals divided by their within-subject cell's standard deviation)
#'   and \code{cell} (a factor labelled with the within-subject levels).
#' @noRd
.rm_within_residuals <- function(d, response, subject, within, between) {
  y <- d[[response]]
  grp <- if (length(between) > 0L) {
    do.call(paste, c(lapply(between, function(b) as.integer(d[[b]])),
                     list(sep = "\r")))
  } else rep("", nrow(d))
  wcodes <- lapply(within, function(w) as.integer(d[[w]]))
  wkey <- do.call(paste, c(wcodes, list(sep = "\r")))
  e <- y - stats::ave(y, d[[subject]]) - stats::ave(y, grp, wkey) +
    stats::ave(y, grp)
  # Residuals at the level of rounding error are zero: testing them would
  # test the arithmetic.
  scale <- max(abs(y))
  if (is.finite(scale) && all(abs(e) <= 1e-11 * scale)) e[] <- 0

  ord <- do.call(order, wcodes)
  keys <- unique(wkey[ord])
  labs <- do.call(paste, c(lapply(within, function(w) as.character(d[[w]])),
                           list(sep = " : ")))
  lab_levels <- labs[match(keys, wkey)]
  if (anyDuplicated(lab_levels)) lab_levels <- make.unique(lab_levels, sep = " #")
  cell <- factor(wkey, levels = keys, labels = lab_levels)
  cell_sd <- stats::ave(e, cell, FUN = function(z) rep(stats::sd(z), length(z)))
  z <- e / cell_sd
  z[!is.finite(z)] <- NA_real_
  list(residuals = e, fitted = y - e, standardised = z, cell = cell)
}

#' Give an emmeans grid built on the internal names the user's names and labels
#' @noRd
.rm_relabel_grid <- function(grid, map) {
  vars <- names(grid@levels)
  new <- lapply(vars, function(v) {
    lab <- map$levels[[v]]
    int <- grid@levels[[v]]
    if (is.null(lab)) return(int)
    lab[match(int, sprintf("L%0*d", nchar(as.character(length(lab))), seq_along(lab)))]
  })
  names(new) <- .rm_unmap_terms(vars, map$to_user)
  got <- .collect_conditions(stats::update(grid, levels = new))
  if (inherits(got$value, "error")) {
    return(list(grid = grid, note = sprintf(
      "The marginal means could not be relabelled with your factor names and levels (%s); they use the internal names listed in $internal_names.",
      conditionMessage(got$value))))
  }
  g <- got$value
  avg <- g@misc$avgd.over
  if (length(avg) > 0L) g@misc$avgd.over <- .rm_unmap_terms(avg, map$to_user)
  list(grid = g, note = character(0))
}

#' Say which factors the marginal means average over
#'
#' emmeans annotates its summary with the factors it averaged over, and the
#' tidy table drops the annotation. When one of those factors interacts
#' significantly with a factor that is shown, the averages can hide or reverse
#' the differences within each of its levels.
#' @noRd
.rm_averaged_note <- function(mesg, specs, within, between, anova_tab,
                              term_int, map) {
  line <- grep("^Results are averaged over the levels of:", mesg, value = TRUE)
  averaged <- setdiff(c(within, between), specs)
  if (length(line) == 0L && length(averaged) == 0L) return(character(0))
  txt <- if (length(line) > 0L) {
    sub("^Results are averaged over the levels of:",
        "The marginal means in $emmeans and the comparisons in $posthoc are averaged over the levels of:",
        line[1L])
  } else {
    sprintf("The marginal means in $emmeans and the comparisons in $posthoc are averaged over the levels of: %s",
            paste(averaged, collapse = ", "))
  }
  txt <- sub("\\.?\\s*$", ".", txt)
  spec_int <- unname(map$to_internal[specs])
  avg_int <- unname(map$to_internal[averaged])
  parts <- strsplit(term_int, ":", fixed = TRUE)
  hit <- vapply(parts, function(p) any(p %in% spec_int) && any(p %in% avg_int),
                logical(1)) & !is.na(anova_tab$p_value) & anova_tab$p_value < 0.05
  if (any(hit)) {
    txt <- paste(txt, sprintf(
      "The ANOVA table shows a significant interaction with a factor averaged over (%s), so these averages can hide or reverse the differences within its levels; set emm_specs = c(%s) to compare the cells.",
      paste(sprintf("%s, %s", anova_tab$term[hit],
                    .fmt_p(anova_tab$p_value[hit], digits = 3)), collapse = "; "),
      paste(sprintf("\"%s\"", unique(c(specs, averaged))), collapse = ", ")))
  }
  txt
}

#' Upper-tail F probability that is NA, not a warning, when undefined
#'
#' With no error variance an epsilon is 0/0, and \code{pf()} given NaN or
#' non-positive degrees of freedom warns "NaNs produced" to the console. The
#' $notes already explain the degenerate design; the p-value is simply NA.
#' @noRd
.rm_pf <- function(q, df1, df2) {
  n <- max(length(q), length(df1), length(df2))
  q <- rep_len(q, n); df1 <- rep_len(df1, n); df2 <- rep_len(df2, n)
  ok <- !is.na(q) & is.finite(df1) & is.finite(df2) & df1 > 0 & df2 > 0
  out <- rep(NA_real_, n)
  out[ok] <- stats::pf(q[ok], df1[ok], df2[ok], lower.tail = FALSE)
  out
}
