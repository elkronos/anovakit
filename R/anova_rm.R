#' Repeated Measures Analysis of Variance
#'
#' Fits a repeated measures ANOVA to data in long format using
#' \code{\link[afex]{aov_ez}}, reports Mauchly's test of sphericity with the
#' Greenhouse-Geisser and Huynh-Feldt corrections, generalised and partial eta
#' squared, residual diagnostics, and estimated marginal means.
#'
#' @details
#' \strong{Sphericity comes from the model that was fitted.} \pkg{afex} already
#' computes Mauchly's test and both epsilon corrections as part of
#' \code{aov_ez()}; they are read from there rather than recomputed by reshaping
#' the data and refitting a multivariate model.
#'
#' \strong{Unbalanced subjects.} A repeated measures ANOVA needs every subject
#' to appear in every within-subject cell. Subjects who do not are removed
#' before fitting, and both the count and the identifiers are reported in
#' \code{$subjects_dropped} and in \code{$notes}.
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
#' @param between Character vector or \code{NULL}. Between-subject factors.
#' @param factorize Logical. Convert the response to numeric through its labels.
#'   Default \code{TRUE}; with \code{FALSE} a non-numeric response is an error.
#'   The subject, within- and between-subject columns are always coerced to
#'   factors, whatever this is set to.
#' @param emm_specs Character vector or \code{NULL}. Factors to compute
#'   estimated marginal means over. Defaults to all within-subject factors.
#' @param adjust Character. Multiplicity adjustment for the pairwise
#'   comparisons, applied by \pkg{emmeans}. One of \code{"tukey"},
#'   \code{"sidak"}, \code{"scheffe"}, \code{"dunnettx"}, \code{"bonferroni"},
#'   \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"BH"},
#'   \code{"BY"}, \code{"fdr"} or \code{"none"}. Default \code{"tukey"}.
#' @param conf_level Numeric in (0, 1). Level for every interval returned.
#'   Default \code{0.95}.
#' @param correction Character. Sphericity correction reported in the ANOVA
#'   table: \code{"GG"} (default), \code{"HF"} or \code{"none"}.
#' @param posthoc Logical. Compute pairwise comparisons. There are
#'   \code{choose(k, 2)} of them, so this is worth turning off when the number
#'   of cells is large; \code{$notes} records that they were skipped. Default
#'   \code{TRUE}.
#' @param plots Logical. Build \pkg{ggplot2} objects. They are returned in
#'   \code{$plots}, never drawn. Default \code{TRUE}.
#' @param verbose Logical. Emit progress through \code{\link[base]{message}}.
#'   Default \code{FALSE}.
#' @param ... Further arguments passed to \code{\link[afex]{aov_ez}}, for
#'   example \code{fun_aggregate} or \code{observed}.
#'
#' @return An \code{\link{anovakit_fit}} object, with
#'   \code{$subjects_dropped} listing any subjects removed for incomplete
#'   within-subject cells, \code{$sphericity} holding Mauchly's test and the
#'   epsilon corrections, \code{$residuals} the within-cell residuals the
#'   diagnostics are built on, and \code{$n_removed_missing},
#'   \code{$n_removed_unbalanced} and \code{$n_removed_aggregated} breaking
#'   \code{$n_removed} into its three causes.
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
#'   # means over the within-subject factor
#'   print(fit$effect_sizes)
#'   print(fit$emmeans)
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
  .check_columns(data, response, "Response column")
  .check_columns(data, subject, "Subject column")
  .check_columns(data, within, "Within-subject factor(s)")
  if (!is.null(between)) .check_columns(data, between, "Between-subject factor(s)")
  .check_conf_level(conf_level)
  .check_flag(factorize, "factorize")
  .check_flag(posthoc, "posthoc")
  .check_flag(plots, "plots")
  .check_flag(verbose, "verbose")

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
  prep <- .prepare_frame(data, cols,
                         factors = c(subject, within, between),
                         keep_all = FALSE)
  d <- prep$data
  notes <- c(notes, prep$notes)

  for (w in within) {
    if (nlevels(d[[w]]) < 2L) {
      .stopf("Within-subject factor `%s` has %d level(s); at least 2 are required.",
             w, nlevels(d[[w]]))
    }
  }

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
  if (length(unique(d[[subject]])) < 2L) {
    .stopf("Fewer than two subjects remain with a complete set of within-subject cells.")
  }
  if (bal$duplicates > 0L) {
    notes <- c(notes, sprintf(
      "%d subject-by-cell combination(s) appear more than once. afex will aggregate them by their mean; check whether that is what you want.",
      bal$duplicates))
  }

  n_input <- nrow(data)
  n_dropped_pre <- n_input - nrow(d)

  ## Fit ---------------------------------------------------------------------
  .say(verbose, "Fitting the repeated measures model with afex::aov_ez().")
  got <- .collect_conditions(
    afex::aov_ez(id = subject, dv = response, data = d,
                 within = within, between = between, ...))
  fit <- got$value
  if (inherits(fit, "error")) {
    .stopf("afex::aov_ez() failed: %s", conditionMessage(fit))
  }
  # afex warns about aggregation and about contrast defaults. Those are worth
  # knowing and belong in $notes, not on the console.
  notes <- c(notes, .said_note(got$said, "afex::aov_ez()"))

  # $data_used must be the frame afex actually modelled, and the row arithmetic
  # must close over it: nrow($data_used) + $n_removed == nrow(input). When afex
  # aggregates duplicated subject-by-cell rows into their means, those rows
  # leave the analysed frame too, so they are counted here rather than left
  # unaccounted for between the input and the model.
  afex_long <- tryCatch(fit$data$long, error = function(e) NULL)
  data_used <- afex_long %||% d
  n_aggregated <- max(0L, nrow(d) - nrow(data_used))
  if (n_aggregated > 0L) {
    notes <- c(notes, sprintf(
      "afex aggregated %d row(s) down to %d subject-by-cell means before fitting; $data_used is the aggregated frame it modelled, and those %d row(s) are included in $n_removed.",
      nrow(d), nrow(data_used), n_aggregated))
  }
  n_dropped_total <- n_input - nrow(data_used)

  ## ANOVA table and sphericity ----------------------------------------------
  sph <- .sphericity_from_afex(fit)
  notes <- c(notes, sph$note)
  anova_tab <- .afex_anova_table(fit, correction)
  notes <- c(notes, attr(anova_tab, "note"))
  eff <- .afex_effect_sizes(fit)

  if (!is.null(sph$table) && any(sph$table$p_value < 0.05, na.rm = TRUE)) {
    bad <- sph$table$term[!is.na(sph$table$p_value) & sph$table$p_value < 0.05]
    notes <- c(notes, sprintf(
      "Mauchly's test rejects sphericity for %s. The reported table uses the %s correction.",
      paste(bad, collapse = ", "),
      switch(correction, GG = "Greenhouse-Geisser", HF = "Huynh-Feldt",
             none = "no")))
  }

  ## Residual diagnostics ----------------------------------------------------
  res <- .rm_residuals(fit)
  norm <- .check_normality(res$residuals)
  notes <- c(notes, norm$note)

  ## Marginal means ----------------------------------------------------------
  specs <- emm_specs %||% within
  .check_names(specs, "emm_specs")
  unknown <- setdiff(specs, c(within, between))
  if (length(unknown) > 0L) {
    .stopf("`emm_specs` must name within- or between-subject factors; %s did not.",
           paste(sprintf("`%s`", unknown), collapse = ", "))
  }
  emm <- .emmeans_grid(fit, specs, type = "link")
  notes <- c(notes, emm$note)
  emm_tab <- .emmeans_table(emm$grid, conf_level, protect = specs)
  notes <- c(notes, emm_tab$note)
  ph <- if (posthoc) {
    .emmeans_pairs(emm$grid, adjust = adjust, conf_level = conf_level)
  } else list(table = NULL,
              note = .no_posthoc_note(NROW(emm_tab$table)))
  notes <- c(notes, ph$note)

  ## Plots -------------------------------------------------------------------
  plot_list <- list()
  if (plots) {
    .say(verbose, "Building plots.")
    plot_list$residuals <- .plot_resid_fitted(
      res$fitted, res$residuals,
      title = "Residuals vs fitted (within-cell)")
    plot_list$qq <- .plot_qq(
      res$residuals,
      title = "Normal Q-Q plot of within-cell residuals")
    plot_list$index <- .plot_resid_index(res$residuals)
    if (!is.null(emm_tab$table)) {
      plot_list$emmeans <- .plot_emmeans(
        emm_tab$table, specs, estimate = "estimate",
        lower = "conf_low", upper = "conf_high", conf_level = conf_level,
        title = "Estimated marginal means",
        ylab = sprintf("Estimated %s", response))
    }
  }

  .new_fit(
    method       = "Repeated measures analysis of variance",
    call         = cl,
    model        = fit,
    anova        = anova_tab,
    effect_sizes = eff,
    emmeans      = emm_tab$table,
    emmeans_object = emm$grid,
    posthoc      = ph$table,
    assumptions  = list(sphericity = sph$table, normality = norm$test),
    plots        = plot_list,
    data_used    = data_used,
    n_removed    = n_dropped_total,
    conf_level   = conf_level,
    notes        = notes,
    extra        = list(subjects_dropped = bal$dropped,
                        n_removed_missing = prep$n_removed,
                        n_removed_unbalanced = n_dropped_pre - prep$n_removed,
                        n_removed_aggregated = n_aggregated,
                        sphericity = sph$table,
                        residuals = res$residuals)
  )
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

#' Drop subjects that do not appear in every within-subject cell
#' @noRd
.balance_subjects <- function(d, subject, within) {
  # A separator that cannot appear in a level. Pasting with sep = "" would
  # make subject "a" in cell "bc" indistinguishable from subject "ab" in
  # cell "c", miscounting the duplicates the note reports. Written as an
  # escape rather than a literal control byte, which is invisible in source.
  cell <- .cell_vector(d, within, sep = "\001")
  key <- paste(d[[subject]], cell, sep = "\001")
  duplicates <- sum(duplicated(key))
  n_cells <- nlevels(cell)
  complete_by_subject <- tapply(as.character(cell), d[[subject]],
                                function(z) length(unique(z)))
  keep_ids <- names(complete_by_subject)[complete_by_subject == n_cells]
  dropped <- setdiff(as.character(unique(d[[subject]])), keep_ids)
  out <- d[as.character(d[[subject]]) %in% keep_ids, , drop = FALSE]
  out[[subject]] <- droplevels(out[[subject]])
  row.names(out) <- NULL
  list(data = out, dropped = dropped, duplicates = duplicates)
}

#' Mauchly's test and the epsilon corrections, read from the afex fit
#' @noRd
.sphericity_from_afex <- function(fit) {
  s <- tryCatch(suppressWarnings(summary(fit$Anova)), error = function(e) e)
  if (inherits(s, "error") || is.null(s$sphericity.tests)) {
    return(list(table = NULL, note = .sphericity_absent_note(fit)))
  }
  # summary.Anova.mlm returns these as classed matrices; as.data.frame() on
  # them collapses to a single list column, so go through as.matrix() first.
  sp <- as.matrix(unclass(s$sphericity.tests))
  if (nrow(sp) == 0L || ncol(sp) < 2L) {
    return(list(table = NULL, note = .sphericity_absent_note(fit)))
  }
  out <- data.frame(term = rownames(sp),
                    mauchly_w = as.numeric(sp[, 1L]),
                    p_value = as.numeric(sp[, 2L]),
                    stringsAsFactors = FALSE, row.names = NULL)
  eps <- tryCatch(as.matrix(unclass(s$pval.adjustments)),
                  error = function(e) NULL)
  if (!is.null(eps) && nrow(eps) > 0L) {
    idx <- match(out$term, rownames(eps))
    cn <- colnames(eps)
    if ("GG eps" %in% cn) out$gg_epsilon <- as.numeric(eps[idx, "GG eps"])
    if ("HF eps" %in% cn) out$hf_epsilon <- as.numeric(eps[idx, "HF eps"])
    if ("Pr(>F[GG])" %in% cn) out$p_gg <- as.numeric(eps[idx, "Pr(>F[GG])"])
    if ("Pr(>F[HF])" %in% cn) out$p_hf <- as.numeric(eps[idx, "Pr(>F[HF])"])
  }
  list(table = out, note = character(0))
}

#' Tidy the afex ANOVA table at the requested sphericity correction
#' @noRd
.afex_anova_table <- function(fit, correction = "GG") {
  # anova.afex_aov() warns when a Huynh-Feldt epsilon above 1 is truncated.
  # Like every other third-party condition in the package, that belongs in
  # $notes rather than on the console.
  got <- .collect_conditions(
    as.data.frame(stats::anova(fit, correction = correction, es = "pes")))
  tab <- got$value
  said <- got$said
  if (inherits(tab, "error")) {
    tab <- as.data.frame(fit$anova_table)
    said <- character(0)
  }
  out <- data.frame(term = rownames(tab), stringsAsFactors = FALSE)
  ren <- c("num Df" = "num_df", "den Df" = "den_df", "MSE" = "mse",
           "F" = "statistic", "ges" = "generalised_eta_sq",
           "pes" = "partial_eta_sq", "Pr(>F)" = "p_value")
  for (nm in names(tab)) {
    target <- if (nm %in% names(ren)) unname(ren[nm]) else
      gsub("[^A-Za-z0-9]+", "_", tolower(nm))
    if (!target %in% names(out)) out[[target]] <- tab[[nm]]
  }
  attr(out, "correction") <- correction
  attr(out, "note") <- .said_note(said, "afex's anova() method")
  row.names(out) <- NULL
  out
}

#' Effect sizes from the afex table
#' @noRd
.afex_effect_sizes <- function(fit) {
  pes <- tryCatch(as.data.frame(stats::anova(fit, es = "pes")),
                  error = function(e) NULL)
  ges <- tryCatch(as.data.frame(stats::anova(fit, es = "ges")),
                  error = function(e) NULL)
  if (is.null(pes) && is.null(ges)) return(NULL)
  terms <- rownames(pes %||% ges)
  data.frame(
    term = terms,
    partial_eta_sq = if (!is.null(pes) && "pes" %in% names(pes)) pes$pes else NA_real_,
    generalised_eta_sq = if (!is.null(ges) && "ges" %in% names(ges)) ges$ges else NA_real_,
    stringsAsFactors = FALSE, row.names = NULL
  )
}

#' Within-cell residuals and fitted values for a repeated measures fit
#'
#' Taken from the multivariate model afex builds, flattened explicitly to a
#' vector rather than left as a matrix for a downstream function to mis-index.
#' @noRd
.rm_residuals <- function(fit) {
  lm_fit <- fit$lm
  if (is.null(lm_fit)) {
    return(list(residuals = numeric(0), fitted = numeric(0)))
  }
  r <- stats::residuals(lm_fit)
  f <- stats::fitted(lm_fit)
  list(residuals = as.numeric(r), fitted = as.numeric(f))
}

#' Explain why sphericity is absent, rather than saying only that it is
#' @noRd
.sphericity_absent_note <- function(fit) {
  n_lev <- tryCatch({
    idata <- fit$data$idata
    if (is.null(idata)) NULL else vapply(idata, function(z) nlevels(as.factor(z)),
                                         integer(1))
  }, error = function(e) NULL)
  if (!is.null(n_lev) && length(n_lev) > 0L && all(n_lev <= 2L)) {
    return("Sphericity is not an issue for this design: every within-subject factor has two levels, where the assumption holds automatically.")
  }
  "Sphericity could not be assessed for this design."
}
