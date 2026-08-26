# Estimated marginal means -------------------------------------------------
#
# One wrapper around emmeans, so that the argument spelling is right in one
# place and a failure is reported rather than silently degraded.

#' Build an emmeans grid
#'
#' The covariance argument emmeans accepts is \code{vcov.} with a trailing dot,
#' and passing \code{vcov. = NULL} is an error rather than a no-op, so the
#' argument is only supplied when there is a matrix to supply.
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
  args <- list(object = model,
               specs = .formula(NULL, specs, quote_terms = TRUE))
  if (type == "response") args$type <- "response"
  if (!is.null(vcov_matrix)) args$vcov. <- vcov_matrix
  args <- c(args, list(...))

  # emmeans itself is quiet, but the methods it calls to build a reference
  # grid are not: summary.lm() warns about a perfect fit, for one. Capture
  # rather than suppress, so nothing is lost and nothing reaches the console.
  got <- .collect_conditions(do.call(emmeans::emmeans, args))
  grid <- got$value
  if (inherits(grid, "error")) {
    return(list(grid = NULL, note = sprintf(
      "Estimated marginal means could not be computed: %s",
      conditionMessage(grid))))
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
#' @noRd
.no_posthoc_note <- function(n_cells = NA_integer_, each = 1L, extra = NULL) {
  base <- if (!is.na(n_cells) && n_cells >= 2L) {
    # choose() returns a double, and sprintf("%d", ) rejects anything above
    # .Machine$integer.max -- which choose(65537, 2) exceeds. Many groups is
    # the exact case this switch exists for, so it must not be the case that
    # breaks the message.
    n <- choose(n_cells, 2L) * each
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
#' @return list with \code{table} (data.frame or \code{NULL}) and \code{note}.
#' @noRd
.emmeans_table <- function(grid, conf_level = 0.95, protect = character(0)) {
  if (is.null(grid)) return(list(table = NULL, note = character(0)))
  s <- tryCatch(
    suppressWarnings(suppressMessages(
      summary(grid, infer = c(TRUE, FALSE), level = conf_level))),
    error = function(e) e)
  if (inherits(s, "error")) {
    return(list(table = NULL, note = sprintf(
      "Estimated marginal means could not be summarised: %s",
      conditionMessage(s))))
  }
  tab <- .normalise_emm_names(as.data.frame(s), protect = protect)
  list(table = tab, note = .renamed_note(tab))
}

#' Pairwise contrasts from an emmeans grid
#'
#' Uses \code{emmeans::contrast()}. \code{emmeans::pairs} is not an exported
#' object, so calling it fails.
#'
#' @param grid an \code{emmGrid}.
#' @param adjust multiplicity adjustment passed to \code{emmeans::contrast}.
#' @param conf_level confidence level.
#' @param ratios \code{TRUE} for a log link, where contrasts are ratios.
#' @noRd
#
# No `protect` here. A pairwise contrast table has no grouping columns -- the
# groups have been contrasted away into a `contrast` column -- so a name in
# `protect` can only match emmeans' own statistic, which then escapes
# normalisation and leaves the schema inconsistent between functions.
.emmeans_pairs <- function(grid, adjust = "tukey", conf_level = 0.95,
                           ratios = FALSE) {
  if (is.null(grid)) return(list(table = NULL, note = character(0)))
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
      summary(ct, infer = c(TRUE, TRUE), level = conf_level))),
    error = function(e) e)
  if (inherits(s, "error")) {
    return(list(table = NULL, note = sprintf(
      "Pairwise comparisons could not be summarised: %s", conditionMessage(s))))
  }
  out <- .normalise_emm_names(as.data.frame(s))
  note <- .renamed_note(out)
  if (ratios && "estimate" %in% names(out) && !"ratio" %in% names(out)) {
    names(out)[names(out) == "estimate"] <- "ratio"
  }
  list(table = out, note = note)
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
  V <- tryCatch(sandwich::vcovHC(model, type = vcov_type),
                error = function(e) e)
  if (inherits(V, "error")) {
    return(list(matrix = NULL, note = sprintf(
      "Robust covariance (%s) failed, model-based standard errors were used instead: %s",
      vcov_type, conditionMessage(V))))
  }
  list(matrix = V, note = character(0))
}

#' Valid heteroscedasticity-consistent covariance types
#' @noRd
.vcov_choices <- function() c("model", "HC0", "HC1", "HC2", "HC3", "HC4")
