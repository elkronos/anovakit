# Estimated marginal means -------------------------------------------------
#
# One wrapper around emmeans, so that the argument spelling is right in one
# place and a failure is reported rather than silently degraded.

#' Clear the emmeans options that change what a result means
#'
#' \code{emmeans::emm_options()} sets session-wide defaults for the adjustment,
#' sidedness, null value, degrees of freedom and scale of every summary. Left in
#' force, a setting made for some other analysis would silently change
#' \code{$emmeans} and \code{$posthoc} here while the \code{adjust} argument
#' still claimed otherwise. Only the computational limits (such as
#' \code{rg.limit}) are kept, since those are how a user lets a large grid
#' through.
#' @return the previous value of \code{getOption("emmeans")}, for restoring.
#' @noRd
.emm_options_guard <- function() {
  old <- getOption("emmeans")
  limits <- c("rg.limit", "lmer.df", "disable.pbkrtest", "pbkrtest.limit",
              "disable.lmerTest", "lmerTest.limit")
  keep <- old[intersect(names(old), limits)]
  options(emmeans = if (length(keep) > 0L) keep else NULL)
  old
}

#' Build an emmeans grid
#'
#' The covariance argument emmeans accepts is \code{vcov.} with a trailing dot,
#' and passing \code{vcov. = NULL} is an error rather than a no-op, so the
#' argument is only supplied when there is a matrix to supply.
#'
#' Two emmeans defaults are overridden. A covariate with only two distinct
#' values (a 0/1 indicator, say) is by default treated as a factor and
#' averaged over its two values, so the "adjusted" means would be read at
#' 0.5 rather than at the covariate's mean; \code{cov.keep = character(0)}
#' holds every covariate at its mean. And a quasi-likelihood GLM is by default
#' given asymptotic (z) inference, while its own coefficient table and omnibus
#' test use t and F; its residual degrees of freedom are passed so that all
#' three agree.
#'
#' @param model a fitted model.
#' @param specs character vector of factor names to marginalise over.
#' @param type \code{"response"} or \code{"link"}.
#' @param vcov_matrix optional covariance matrix (for example from
#'   \code{sandwich::vcovHC}).
#' @param ... further arguments passed to \code{emmeans::emmeans}.
#' @return list with \code{grid} (an \code{emmGrid} or \code{NULL}) and
#'   \code{note} (character, empty on success).
#' @noRd
.emmeans_grid <- function(model, specs, type = c("response", "link"),
                          vcov_matrix = NULL, ...) {
  type <- match.arg(type)
  old <- .emm_options_guard()
  on.exit(options(emmeans = old), add = TRUE)
  args <- list(object = model,
               specs = .formula(NULL, specs, quote_terms = TRUE))
  if (type == "response") args$type <- "response"
  if (!is.null(vcov_matrix)) args$vcov. <- vcov_matrix
  dots <- list(...)
  if (!"cov.keep" %in% names(dots)) args$cov.keep <- character(0)
  if (!"df" %in% names(dots) && inherits(model, "glm") &&
      .estimates_dispersion(model)) {
    dfr <- stats::df.residual(model)
    if (!is.null(dfr) && is.finite(dfr) && dfr >= 1) args$df <- dfr
  }
  args <- c(args, dots)

  # emmeans itself is quiet, but the methods it calls to build a reference
  # grid are not: summary.lm() warns about a perfect fit, for one. Capture
  # rather than suppress, so nothing is lost and nothing reaches the console.
  got <- .collect_conditions(do.call(emmeans::emmeans, args))
  grid <- got$value
  if (inherits(grid, "error")) {
    msg <- conditionMessage(grid)
    if (grepl("rg.limit", msg, fixed = TRUE)) {
      return(list(grid = NULL, note = paste(
        "Estimated marginal means were not computed: the reference grid has more cells than emmeans allows by default.",
        "Raise the limit with emmeans::emm_options(rg.limit = ...) before calling, or set posthoc = FALSE and plots = FALSE.")))
    }
    return(list(grid = NULL, note = sprintf(
      "Estimated marginal means could not be computed: %s", msg)))
  }
  list(grid = grid, note = .said_note(got$said, "emmeans"))
}

#' Standard column names used across the package
#'
#' emmeans names its estimate column \code{emmean}, \code{response},
#' \code{rate} or \code{ratio} depending on the link and the request, and its
#' interval columns \code{lower.CL}/\code{upper.CL} or
#' \code{asymp.LCL}/\code{asymp.UCL} depending on the degrees of freedom. Both
#' are normalised here so that callers see one set of names.
#'
#' A grouping column may already be called \code{se}, \code{df} or
#' \code{conf_low}. Renaming onto it would leave the group labels sitting under
#' a statistic's name -- so when a target is taken by a column that is not
#' itself being renamed, the statistic gets an \code{emm_} prefix instead, and
#' the collision is reported so nobody has to notice it from the column order.
#'
#' @return the data frame, with a \code{"renamed"} attribute naming any column
#'   that had to take a prefixed name.
#' @noRd
.normalise_emm_names <- function(df, protect = character(0)) {
  df <- as.data.frame(df)
  ren <- c(
    emmean = "estimate", response = "estimate", rate = "estimate",
    ratio = "estimate", prob = "estimate", odds.ratio = "estimate",
    SE = "se", df = "df",
    lower.CL = "conf_low", upper.CL = "conf_high",
    asymp.LCL = "conf_low", asymp.UCL = "conf_high",
    LCL = "conf_low", UCL = "conf_high",
    t.ratio = "statistic", z.ratio = "statistic",
    p.value = "p_value", null = "null"
  )
  nm <- names(df)
  # A grouping column may itself be called `df`, `se` or `estimate`. It is not
  # an emmeans statistic and must not be renamed or written over. emmeans puts
  # the grid's factor columns first, so protect the FIRST occurrence of each
  # protected name only: a second column of the same name is emmeans' own
  # statistic and still needs a distinct name.
  protect_idx <- vapply(intersect(protect, nm), function(p) match(p, nm),
                        integer(1))
  idx <- which(nm %in% names(ren) & !(seq_along(nm) %in% protect_idx))
  if (length(idx) == 0L) return(df)

  # Columns that are staying put: anything not being renamed, protected
  # grouping columns included.
  fixed <- nm[-idx]
  displaced <- character(0)
  new <- nm
  for (i in idx) {
    target <- unname(ren[nm[i]])
    if (target %in% fixed || target %in% new[-i]) {
      prefixed <- paste0("emm_", target)
      k <- 1L
      while (prefixed %in% c(fixed, new[-i])) {
        k <- k + 1L
        prefixed <- paste0("emm_", target, k)
      }
      displaced <- c(displaced, sprintf("%s -> %s", target, prefixed))
      new[i] <- prefixed
    } else {
      new[i] <- target
    }
  }
  names(df) <- new
  # emmeans' summary carries a subclass whose print method reads attributes that
  # no longer describe these columns. Once renamed, this is a plain data frame.
  df <- as.data.frame(unclass(df)[names(df)], stringsAsFactors = FALSE)
  names(df) <- new
  if (length(displaced) > 0L) attr(df, "renamed") <- displaced
  df
}

#' The note describing any emmeans column that had to be renamed
#' @noRd
.renamed_note <- function(df) {
  r <- attr(df, "renamed")
  if (is.null(r)) return(character(0))
  sprintf("A grouping column shares its name with a result column, so the result column was renamed (%s). The grouping column is unchanged.",
          paste(r, collapse = "; "))
}

#' The note explaining an empty \code{$posthoc}
#'
#' \code{posthoc = FALSE} is an efficiency switch, not a silent one. Without a
#' note the caller cannot tell a suppressed comparison from one that could not
#' be computed: both leave \code{$posthoc} \code{NULL}.
#'
#' @param n_cells number of estimable cells, or \code{NA} when unknown.
#' @param each how many times the full set of comparisons would have been made
#'   (one per response, in \code{\link{anova_manova}}).
#' @param extra anything else to append to the sentence.
#' @param n_comparisons the number of comparisons, when it is not
#'   \code{choose(n_cells, 2)} (per-factor comparisons, for instance).
#' @noRd
.no_posthoc_note <- function(n_cells = NA_integer_, each = 1L, extra = NULL,
                             n_comparisons = NULL) {
  n <- if (!is.null(n_comparisons)) {
    n_comparisons * each
  } else if (!is.na(n_cells) && n_cells >= 2L) {
    # choose() returns a double, and sprintf("%d", ) rejects anything above
    # .Machine$integer.max -- which choose(65537, 2) exceeds. Many groups is
    # the exact case this switch exists for, so it must not be the case that
    # breaks the message.
    choose(n_cells, 2L) * each
  } else {
    NA_real_
  }
  base <- if (!is.na(n) && n >= 1) {
    sprintf("Pairwise comparisons were not computed (posthoc = FALSE); there would have been %s.",
            format(n, scientific = FALSE, big.mark = ""))
  } else {
    "Pairwise comparisons were not computed (posthoc = FALSE)."
  }
  paste(c(base, extra), collapse = " ")
}

#' Summarise an emmeans grid into a tidy data frame
#'
#' @param grid an \code{emmGrid}.
#' @param conf_level confidence level.
#' @return list with \code{table} (data.frame or \code{NULL}), \code{note},
#'   and \code{mesg}, the annotations emmeans prints under a summary (which
#'   factors were averaged over, for instance) and which \code{as.data.frame()}
#'   discards.
#' @noRd
.emmeans_table <- function(grid, conf_level = 0.95, protect = character(0)) {
  if (is.null(grid)) return(list(table = NULL, note = character(0), mesg = character(0)))
  old <- .emm_options_guard()
  on.exit(options(emmeans = old), add = TRUE)
  s <- tryCatch(
    suppressWarnings(suppressMessages(
      summary(grid, infer = c(TRUE, FALSE), level = conf_level))),
    error = function(e) e)
  if (inherits(s, "error")) {
    return(list(table = NULL, mesg = character(0), note = sprintf(
      "Estimated marginal means could not be summarised: %s",
      conditionMessage(s))))
  }
  tab <- .normalise_emm_names(as.data.frame(s), protect = protect)
  list(table = tab, note = .renamed_note(tab),
       mesg = as.character(attr(s, "mesg")))
}

#' Pairwise contrasts from an emmeans grid
#'
#' Uses \code{emmeans::contrast()}. \code{emmeans::pairs} is not an exported
#' object, so calling it fails.
#'
#' The table carries both the unadjusted p-value (\code{p_value}) and the
#' multiplicity-adjusted one (\code{p_adjusted}), with the method in
#' \code{adjustment}, which is the layout \code{\link{anova_welch}} and
#' \code{\link{anova_kw}} use as well. The intervals are emmeans' adjusted
#' (simultaneous) intervals. emmeans cannot turn a step-down method such as
#' Holm or Benjamini-Hochberg into intervals and uses Bonferroni for them
#' instead; when that happens a note says so.
#'
#' @param grid an \code{emmGrid}.
#' @param adjust multiplicity adjustment passed to \code{emmeans::contrast}.
#' @param conf_level confidence level.
#' @param ratios \code{TRUE} for a log or logit link on the response scale,
#'   where emmeans back-transforms contrasts to ratios.
#' @param regrid \code{TRUE} to compare on the response scale for a link that
#'   emmeans cannot express as ratios (inverse, probit, square root, ...).
#'   Without it the contrasts would be differences on the link scale.
#' @noRd
#
# No `protect` here. A pairwise contrast table has no grouping columns -- the
# groups have been contrasted away into a `contrast` column -- so a name in
# `protect` can only match emmeans' own statistic, which then escapes
# normalisation and leaves the schema inconsistent between functions.
.emmeans_pairs <- function(grid, adjust = "tukey", conf_level = 0.95,
                           ratios = FALSE, regrid = FALSE) {
  if (is.null(grid)) return(list(table = NULL, note = character(0)))
  old <- .emm_options_guard()
  on.exit(options(emmeans = old), add = TRUE)
  if (regrid) {
    rg <- tryCatch(suppressWarnings(suppressMessages(emmeans::regrid(grid))),
                   error = function(e) e)
    if (inherits(rg, "error")) {
      return(list(table = NULL, note = sprintf(
        "Pairwise comparisons could not be computed on the response scale: %s",
        conditionMessage(rg))))
    }
    grid <- rg
  }
  ct <- tryCatch(
    suppressWarnings(suppressMessages(
      emmeans::contrast(grid, method = "pairwise", adjust = adjust))),
    error = function(e) e)
  if (inherits(ct, "error")) {
    return(list(table = NULL, note = sprintf(
      "Pairwise comparisons could not be computed: %s", conditionMessage(ct))))
  }
  s <- tryCatch(
    suppressWarnings(suppressMessages(
      summary(ct, infer = c(TRUE, TRUE), level = conf_level, adjust = adjust))),
    error = function(e) e)
  raw <- tryCatch(
    suppressWarnings(suppressMessages(
      summary(ct, infer = c(FALSE, TRUE), adjust = "none"))),
    error = function(e) e)
  if (inherits(s, "error") || inherits(raw, "error")) {
    err <- if (inherits(s, "error")) s else raw
    return(list(table = NULL, note = sprintf(
      "Pairwise comparisons could not be summarised: %s", conditionMessage(err))))
  }
  out <- .normalise_emm_names(as.data.frame(s))
  note <- .renamed_note(out)
  if (ratios && "estimate" %in% names(out) && !"ratio" %in% names(out)) {
    names(out)[names(out) == "estimate"] <- "ratio"
  }
  if ("p_value" %in% names(out)) {
    pos <- match("p_value", names(out))
    adj <- out$p_value
    out$p_value <- as.data.frame(raw)[["p.value"]]
    out <- data.frame(out[seq_len(pos)], p_adjusted = adj, adjustment = adjust,
                      out[-seq_len(pos)], stringsAsFactors = FALSE,
                      check.names = FALSE)
  }
  mesg <- as.character(attr(s, "mesg"))
  ci_line <- grep("^Conf-level adjustment:", mesg, value = TRUE)
  ci_method <- if (length(ci_line) > 0L) {
    sub("^Conf-level adjustment: ([^ ]+) method.*$", "\\1", ci_line[1L])
  } else {
    "none"
  }
  p_method <- if (identical(adjust, "none")) "none" else adjust
  if (!identical(ci_method, p_method) &&
      !(p_method == "tukey" && ci_method %in% c("tukey", "none"))) {
    note <- c(note, sprintf(
      "The comparison intervals use the %s adjustment while the p-values use %s: emmeans cannot turn %s into intervals.",
      ci_method, p_method, p_method))
  }
  list(table = out, note = note)
}

#' Marginal means and pairwise comparisons for the grouping factors of a model
#'
#' When the grouping factors enter a model additively, a comparison between
#' two cells of their grid is a sum of main-effect differences, and every
#' main-effect difference is repeated once for each level of the other
#' factors. Tukey-adjusting all of those cell pairs inflates the family (and
#' the running time grows with the fourth power of the number of cells) for
#' comparisons that are not separately informative. For an additive model the
#' marginal means and comparisons are therefore computed for each factor on
#' its own, averaged over the others, with one multiplicity family per factor.
#' When the model contains interactions among the grouping factors the cell
#' grid is used.
#'
#' Cells of the grid that contain no data are either not estimable (dropped,
#' with a note) or, under a model without the full interaction, estimated by
#' extrapolation from the other cells (kept, with a note naming them).
#'
#' @param model fitted model.
#' @param groups grouping factor names.
#' @param additive does the model contain no interaction among \code{groups}?
#' @param data the analysed frame, used to find the observed cells.
#' @param link the link of a GLM, which decides how comparisons are expressed.
#' @param max_comparisons skip the comparisons, with a note, above this many.
#' @return list with \code{grid} (an \code{emmGrid}, or a named list of them
#'   when per factor), \code{table}, \code{posthoc}, \code{notes} and
#'   \code{per_factor}.
#' @noRd
.emm_block <- function(model, groups, additive = FALSE, type = "response",
                       vcov_matrix = NULL, conf_level = 0.95, adjust = "tukey",
                       posthoc = TRUE, link = "identity", data = NULL,
                       max_comparisons = 5000L, n_each = 1L, ...) {
  notes <- character(0)
  per_factor <- isTRUE(additive) && length(groups) > 1L
  ratios <- type == "response" && link %in% c("log", "logit")
  regrid <- type == "response" && !link %in% c("identity", "log", "logit")

  if (per_factor) {
    grids <- stats::setNames(vector("list", length(groups)), groups)
    tabs <- list()
    for (g in groups) {
      e <- .emmeans_grid(model, g, type = type, vcov_matrix = vcov_matrix, ...)
      notes <- c(notes, e$note)
      grids[g] <- list(e$grid)
      t <- .emmeans_table(e$grid, conf_level, protect = groups)
      notes <- c(notes, t$note)
      if (!is.null(t$table)) tabs[[g]] <- t$table
    }
    table <- NULL
    term_col <- NULL
    if (length(tabs) > 0L) {
      stat_cols <- unique(unlist(lapply(tabs, function(t) setdiff(names(t), groups))))
      tabs <- lapply(names(tabs), function(g) {
        t <- tabs[[g]]
        for (h in groups) t[[h]] <- if (h == g) as.character(t[[h]]) else NA_character_
        for (sc in setdiff(stat_cols, names(t))) t[[sc]] <- NA
        t <- t[c(groups, stat_cols)]
        .cbind_label(t, "term", g)
      })
      table <- do.call(rbind, tabs)
      row.names(table) <- NULL
      term_col <- names(table)[1L]
      attr(table, "per_factor") <- TRUE
      attr(table, "term_col") <- term_col
    }
    notes <- c(notes, sprintf(
      "The grouping factors enter the model additively, so $emmeans and $posthoc are reported for each factor separately, averaged over the others (column `%s`), and each factor's comparisons are a separate multiplicity family. Set interaction = TRUE to compare cells instead.",
      term_col %||% "term"))
    k <- vapply(grids, function(gr) if (is.null(gr)) 0L else nrow(summary(gr)), integer(1))
    n_comp <- sum(choose(k, 2L))
    ph <- NULL
    if (!posthoc) {
      notes <- c(notes, .no_posthoc_note(n_comparisons = n_comp, each = n_each))
    } else if (n_comp * n_each > max_comparisons) {
      notes <- c(notes, .too_many_note(n_comp * n_each, max_comparisons))
    } else {
      phs <- list()
      for (g in groups) {
        if (is.null(grids[[g]]) || k[[g]] < 2L) next
        pr <- .emmeans_pairs(grids[[g]], adjust = adjust, conf_level = conf_level,
                             ratios = ratios, regrid = regrid)
        notes <- c(notes, pr$note)
        if (!is.null(pr$table)) phs[[g]] <- .cbind_label(pr$table, "term", g)
      }
      if (length(phs) > 0L) {
        ph <- do.call(rbind, unname(phs))
        row.names(ph) <- NULL
      }
    }
    return(list(grid = grids, table = table, posthoc = ph,
                notes = unique(notes), per_factor = TRUE))
  }

  e <- .emmeans_grid(model, groups, type = type, vcov_matrix = vcov_matrix, ...)
  notes <- c(notes, e$note)
  t <- .emmeans_table(e$grid, conf_level, protect = groups)
  notes <- c(notes, t$note)
  table <- t$table
  est_col <- if (!is.null(table)) .emm_col(table, "estimate", groups) else NULL

  # Cells of the grid that hold no data.
  if (!is.null(table) && !is.null(data) && length(groups) > 1L &&
      all(groups %in% names(data)) && all(groups %in% names(table))) {
    key_t <- do.call(paste, c(lapply(table[groups], as.character), sep = "\r"))
    key_d <- unique(do.call(paste, c(lapply(data[groups], as.character), sep = "\r")))
    empty <- !key_t %in% key_d
    if (any(empty)) {
      lab <- do.call(paste, c(lapply(table[empty, groups, drop = FALSE], as.character),
                              sep = " : "))
      est_na <- if (is.null(est_col)) rep(TRUE, nrow(table)) else is.na(table[[est_col]])
      nonest <- empty & est_na
      if (any(nonest)) {
        notes <- c(notes, sprintf(
          "%d level combination(s) contain no data and are not estimable under this model, so they are left out of $emmeans and $posthoc: %s.",
          sum(nonest), .abbrev(lab[nonest[empty]])))
      }
      if (any(empty & !nonest)) {
        notes <- c(notes, sprintf(
          "%d level combination(s) contain no data; their marginal means are extrapolated from the other cells under the fitted model rather than observed: %s.",
          sum(empty & !nonest), .abbrev(lab[!nonest[empty]])))
      }
      keep <- !nonest
      at <- attributes(table)[c("renamed")]
      table <- table[keep, , drop = FALSE]
      row.names(table) <- NULL
      if (!is.null(at$renamed)) attr(table, "renamed") <- at$renamed
    }
  }
  if (!is.null(table) && !is.null(est_col)) {
    # A non-estimable cell that is not empty cannot arise from data, but
    # emmeans can still mark one (a fully aliased contrast); drop it too.
    table <- table[!is.na(table[[est_col]]), , drop = FALSE]
    row.names(table) <- NULL
  }

  n_cells <- NROW(table)
  n_comp <- choose(n_cells, 2L)
  ph <- NULL
  if (!posthoc) {
    notes <- c(notes, .no_posthoc_note(n_cells, each = n_each))
  } else if (n_cells >= 2L && n_comp * n_each > max_comparisons) {
    notes <- c(notes, .too_many_note(n_comp * n_each, max_comparisons))
  } else if (n_cells >= 2L) {
    pr <- .emmeans_pairs(e$grid, adjust = adjust, conf_level = conf_level,
                         ratios = ratios, regrid = regrid)
    notes <- c(notes, pr$note)
    ph <- pr$table
    if (!is.null(ph)) {
      ec <- intersect(c("estimate", "ratio", "emm_estimate"), names(ph))[1L]
      if (!is.na(ec)) ph <- ph[!is.na(ph[[ec]]), , drop = FALSE]
      row.names(ph) <- NULL
    }
  }
  list(grid = e$grid, table = table, posthoc = ph, notes = unique(notes),
       per_factor = FALSE)
}

#' The note for comparisons skipped because there would be too many
#' @noRd
.too_many_note <- function(n, limit) {
  sprintf("Pairwise comparisons were skipped: there would be %s, more than the %s this function computes by default (the time and memory emmeans needs grow with the fourth power of the number of cells). Compare a subset with emmeans::contrast() on $emmeans_object, or set posthoc = FALSE.",
          format(n, scientific = FALSE, big.mark = ","),
          format(limit, big.mark = ","))
}

#' Multiplicity adjustment methods accepted by emmeans
#'
#' \code{"mvt"} is deliberately absent: it is computed by simulation, so it
#' returns a slightly different p-value on every call and consumes the caller's
#' random stream. \code{"tukey"} and \code{"sidak"} give deterministic answers
#' for the same problems.
#' @noRd
.adjust_choices <- function() {
  c("tukey", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
    "none", "sidak", "scheffe", "dunnettx")
}

#' Validate a multiplicity adjustment against a known set
#' @noRd
.check_adjust <- function(adjust, allowed, arg = "adjust") {
  if (!is.character(adjust) || length(adjust) != 1L || is.na(adjust)) {
    .stopf("`%s` must be a single character string.", arg)
  }
  if (!adjust %in% allowed) {
    .stopf("`%s` must be one of: %s; got \"%s\".",
           arg, paste(allowed, collapse = ", "), adjust)
  }
  invisible(adjust)
}

#' Robust covariance matrix, when requested and {sandwich} is installed
#'
#' Three situations make a sandwich covariance wrong rather than merely
#' different, and each falls back to the model-based covariance with a note:
#' fitted values on the boundary of the parameter space (separation in a
#' binary model, an all-zero group in a count model), where the sandwich
#' collapses towards zero and reports astronomically small p-values; and an
#' observation with leverage 1 (a single-observation cell), where HC2 to HC4
#' divide by zero. Frequency weights are a fourth, which the function cannot
#' detect: it says what the sandwich assumes instead.
#'
#' @param model a fitted model.
#' @param vcov_type \code{"model"} or an HC type such as \code{"HC0"}.
#' @return list with \code{matrix} (or \code{NULL}) and \code{note}.
#' @noRd
.robust_vcov <- function(model, vcov_type = "model") {
  if (identical(vcov_type, "model")) {
    return(list(matrix = NULL, note = character(0)))
  }
  if (!requireNamespace("sandwich", quietly = TRUE)) {
    return(list(matrix = NULL, note = sprintf(
      "vcov_type = \"%s\" requested but the {sandwich} package is not installed; model-based standard errors were used instead.",
      vcov_type)))
  }
  if (inherits(model, "glm")) {
    fam <- tryCatch(stats::family(model)$family, error = function(e) "")
    mu <- stats::fitted(model)
    at_boundary <- if (fam %in% c("binomial", "quasibinomial")) {
      any(mu < 1e-8 | mu > 1 - 1e-8)
    } else if (fam %in% c("poisson", "quasipoisson") || inherits(model, "negbin") ||
               grepl("^Negative Binomial", fam)) {
      any(mu < 1e-8)
    } else {
      FALSE
    }
    if (at_boundary) {
      return(list(matrix = NULL, note = sprintf(
        "Robust standard errors (%s) were not used: some fitted values are on the boundary (a group with no events, or all events), where the sandwich covariance collapses towards zero and reports spuriously small p-values. Model-based standard errors were used instead.",
        vcov_type)))
    }
  }
  got <- .collect_conditions(sandwich::vcovHC(model, type = vcov_type))
  V <- got$value
  if (inherits(V, "error")) {
    return(list(matrix = NULL, note = sprintf(
      "Robust covariance (%s) failed, model-based standard errors were used instead: %s",
      vcov_type, conditionMessage(V))))
  }
  if (!all(is.finite(V))) {
    h <- tryCatch(stats::hatvalues(model), error = function(e) numeric(0))
    return(list(matrix = NULL, note = sprintf(
      "Robust covariance (%s) is undefined here: %d observation(s) have leverage 1 (a group or cell with a single observation), and %s divides by 1 - leverage. Model-based standard errors were used instead; HC0 and HC1 do not have this problem.",
      vcov_type, sum(h > 1 - 1e-8), vcov_type)))
  }
  notes <- .said_note(got$said, "sandwich::vcovHC()")
  w <- tryCatch(stats::weights(model, type = "prior"), error = function(e) NULL)
  if (is.null(w) && !is.null(model$prior.weights)) w <- model$prior.weights
  if (!is.null(w) && any(w != 1)) {
    notes <- c(notes, "Robust standard errors treat each row as one independent unit. That is right when a row is one binomial observation of several trials, or when the weights are sampling weights; it is wrong when the weights count identical rows (frequency weights), which it treats as single units with inflated scores. For frequency-weighted data, expand to one row per observation or use vcov_type = \"model\".")
  }
  if (!inherits(model, "glm") ||
      identical(tryCatch(stats::family(model)$family, error = function(e) ""), "gaussian")) {
    notes <- c(notes, "Comparisons that use robust standard errors keep the residual degrees of freedom of the model; with very small groups of unequal variance they can be anti-conservative, and for a one-way design anova_welch() is better calibrated.")
  }
  list(matrix = V, note = notes)
}

#' Valid heteroscedasticity-consistent covariance types
#' @noRd
.vcov_choices <- function() c("model", "HC0", "HC1", "HC2", "HC3", "HC4")
