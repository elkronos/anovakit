# Shared input validation and data preparation ------------------------------
#
# Every exported function routes through these helpers so that a bad argument
# produces the same error text wherever it is supplied, and so that missing
# values, factor coercion and cell construction are handled once.

#' Stop with a formatted message and no call context
#' @param fmt sprintf-style format string.
#' @param ... values interpolated into \code{fmt}.
#' @noRd
.stopf <- function(fmt, ...) {
  stop(sprintf(fmt, ...), call. = FALSE)
}

#' Backtick-quote a column name so it can be used in a formula
#'
#' \code{stats::reformulate()} does not quote non-syntactic names, so a column
#' called \code{"my group"} produces a parse error. Quoting first makes every
#' column name usable.
#'
#' @param x character vector of column names.
#' @return \code{x} with each element wrapped in backticks.
#' @noRd
.bq <- function(x) {
  if (length(x) == 0L) return(character(0))
  bad <- grepl("`", x, fixed = TRUE)
  if (any(bad)) {
    # A backtick inside a name cannot be escaped in an R formula. Stripping it
    # would silently retarget the term at a different column.
    .stopf("Column name(s) containing a backtick cannot be used in a model formula: %s. Rename the column.",
           paste(x[bad], collapse = ", "))
  }
  paste0("`", x, "`")
}

#' Build a formula from column names, quoting them
#'
#' @param response single column name, or \code{NULL} for a one-sided formula.
#' @param terms character vector of term labels (already-composed terms such as
#'   \code{"a:b"} are passed through untouched when \code{quote_terms = FALSE}).
#' @param quote_terms should \code{terms} be backtick-quoted?
#' @param env environment attached to the formula.
#' @noRd
.formula <- function(response, terms, quote_terms = TRUE,
                     env = parent.frame()) {
  tl <- if (quote_terms) .bq(terms) else terms
  if (length(tl) == 0L) tl <- "1"
  stats::reformulate(tl,
                     response = if (is.null(response)) NULL else .bq(response),
                     env = env)
}

#' Coerce anything data.frame-like to a plain data.frame
#'
#' \code{data.table} and \code{tbl_df} objects change the meaning of
#' \code{x[, cols, drop = FALSE]}, which is how every function here subsets.
#' Converting once at the door removes a whole class of failure.
#'
#' @param data object to coerce.
#' @param arg argument name used in the error message.
#' @noRd
.as_df <- function(data, arg = "data") {
  if (missing(data) || is.null(data)) {
    .stopf("`%s` is required.", arg)
  }
  if (!is.data.frame(data)) {
    .stopf("`%s` must be a data.frame (or anything that inherits from one, such as a data.table or a tibble); got %s.",
           arg, paste(class(data), collapse = "/"))
  }
  if (!identical(class(data), "data.frame")) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }
  if (nrow(data) == 0L) {
    .stopf("`%s` has no rows.", arg)
  }
  dup <- unique(names(data)[duplicated(names(data))])
  if (length(dup) > 0L) {
    .stopf("`%s` has duplicated column name(s): %s. Referring to one by name would silently pick the first.",
           arg, paste(dup, collapse = ", "))
  }
  data
}

#' Validate that an argument is a single non-empty character string
#' @noRd
.check_name <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    .stopf("`%s` must be a single column name given as a character string.", arg)
  }
  invisible(x)
}

#' Validate that an argument is a character vector of column names
#' @noRd
.check_names <- function(x, arg, min_length = 1L) {
  if (!is.character(x) || length(x) < min_length || anyNA(x) || !all(nzchar(x))) {
    .stopf("`%s` must be a character vector of at least %d column name(s).",
           arg, min_length)
  }
  if (anyDuplicated(x)) {
    .stopf("`%s` contains duplicated column names: %s.",
           arg, paste(unique(x[duplicated(x)]), collapse = ", "))
  }
  invisible(x)
}

#' Validate that named columns exist in the data
#' @noRd
.check_columns <- function(data, cols, arg) {
  missing_cols <- setdiff(cols, names(data))
  if (length(missing_cols) > 0L) {
    .stopf("%s not found in `data`: %s. Available columns: %s.",
           arg, paste(missing_cols, collapse = ", "),
           paste(utils::head(names(data), 30L), collapse = ", "))
  }
  invisible(TRUE)
}

#' Validate a confidence level
#' @noRd
.check_conf_level <- function(conf_level) {
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    .stopf("`conf_level` must be a single number strictly between 0 and 1; got %s.",
           paste(format(conf_level), collapse = ", "))
  }
  invisible(conf_level)
}

#' Validate a length-one logical
#' @noRd
.check_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    .stopf("`%s` must be TRUE or FALSE.", arg)
  }
  invisible(x)
}

#' Validate that a column is numeric
#' @noRd
.check_numeric_col <- function(data, col, role = "Response") {
  if (!is.numeric(data[[col]])) {
    .stopf("%s `%s` must be numeric; it is %s.",
           role, col, paste(class(data[[col]]), collapse = "/"))
  }
  invisible(TRUE)
}

#' Generate a column name that does not collide with existing ones
#'
#' @param data data.frame whose names must be avoided.
#' @param candidate preferred name.
#' @noRd
.safe_name <- function(data, candidate) {
  nm <- candidate
  i <- 0L
  while (nm %in% names(data)) {
    i <- i + 1L
    nm <- paste0(candidate, ".", i)
  }
  nm
}

#' Subset to the analysed columns, coerce grouping variables, drop incomplete rows
#'
#' @param data a plain data.frame.
#' @param cols character vector of every column the analysis touches.
#' @param factors character vector of columns to coerce to factor and droplevels.
#' @param keep_all keep every column of \code{data} (\code{TRUE}) or only
#'   \code{cols} (\code{FALSE}). Complete-case filtering always uses \code{cols}
#'   alone, so an unrelated column full of \code{NA} never removes a row.
#' @return list with \code{data}, \code{n_removed} and \code{notes}.
#' @noRd
.prepare_frame <- function(data, cols, factors = character(0), keep_all = TRUE) {
  cols <- unique(cols)
  sub <- data[, cols, drop = FALSE]
  keep <- stats::complete.cases(sub)

  # Inf and -Inf survive complete.cases() but cannot be modelled. Treat them as
  # missing rather than letting them reach a variance or a link function.
  n_infinite <- 0L
  for (cn in cols) {
    v <- sub[[cn]]
    if (is.numeric(v)) {
      bad <- is.finite(v) == FALSE & !is.na(v)
      n_infinite <- n_infinite + sum(bad & keep)
      keep <- keep & !bad
    }
  }
  n_removed <- sum(!keep)
  out <- data[keep, if (keep_all) seq_along(data) else match(cols, names(data)),
              drop = FALSE]
  row.names(out) <- NULL

  notes <- character(0)
  if (n_removed > 0L) {
    notes <- c(notes, sprintf(
      "Dropped %d row(s) with missing%s values in: %s.",
      n_removed, if (n_infinite > 0L) " or infinite" else "",
      paste(cols, collapse = ", ")))
  }
  if (nrow(out) == 0L) {
    .stopf("No complete cases remain after removing rows with missing values in: %s.",
           paste(cols, collapse = ", "))
  }

  for (f in factors) {
    if (!is.factor(out[[f]])) {
      n_lev <- length(unique(out[[f]]))
      if (is.numeric(out[[f]]) && n_lev > 20L) {
        notes <- c(notes, sprintf(
          "Grouping variable `%s` is numeric with %d distinct values and was converted to a factor. If it is a covariate rather than a group, this is not the function you want.",
          f, n_lev))
      }
      if (is.character(out[[f]])) {
        # factor() orders character levels by the session's collation, so the
        # reference level of a character column depends on the locale. Sort in
        # the C locale so the same data gives the same reference everywhere.
        out[[f]] <- factor(out[[f]],
                           levels = sort(unique(out[[f]]), method = "radix"))
      } else {
        out[[f]] <- factor(out[[f]])
      }
    }
    out[[f]] <- droplevels(out[[f]])
  }

  list(data = out, n_removed = n_removed, notes = notes)
}

#' Check that grouping variables are usable
#'
#' @param data prepared data.frame.
#' @param groups character vector of grouping columns.
#' @param min_levels minimum number of levels each variable must have.
#' @param min_n minimum number of observations required in each populated cell.
#' @noRd
.check_groups <- function(data, groups, min_levels = 2L, min_n = 1L) {
  for (g in groups) {
    lv <- nlevels(droplevels(data[[g]]))
    if (lv < min_levels) {
      .stopf("Grouping variable `%s` has %d level(s) after removing missing values; at least %d are required.",
             g, lv, min_levels)
    }
  }
  cell <- .cell_vector(data, groups)
  tab <- table(cell)
  populated <- sum(tab > 0L)
  if (populated < 2L) {
    .stopf("Only %d group combination(s) contain data; at least 2 are required. Populated cell(s): %s.",
           populated, .abbrev(names(tab)[tab > 0L]))
  }
  if (min_n > 1L && any(tab[tab > 0L] < min_n)) {
    small <- names(tab)[tab > 0L & tab < min_n]
    .stopf("%d group combination(s) have fewer than %d observations: %s.",
           length(small), min_n, .abbrev(small))
  }
  invisible(TRUE)
}

#' Build the combined cell factor for a set of grouping variables
#'
#' Unlike a bare \code{interaction()} call this drops level combinations that
#' never occur, so downstream counts of "how many groups are there" match the
#' number of groups that actually contain data.
#'
#' @noRd
.cell_vector <- function(data, groups, sep = " : ") {
  if (length(groups) == 1L) {
    return(droplevels(as.factor(data[[groups]])))
  }
  f <- lapply(groups, function(g) droplevels(as.factor(data[[g]])))
  droplevels(do.call(interaction, c(f, list(drop = TRUE, sep = sep))))
}

#' Attach the combined cell factor under a non-colliding name
#'
#' @return list with \code{data} (with the new column) and \code{cell} (its name).
#' @noRd
.add_cell <- function(data, groups, sep = " : ", candidate = ".cell") {
  nm <- .safe_name(data, candidate)
  data[[nm]] <- .cell_vector(data, groups, sep = sep)
  list(data = data, cell = nm)
}

#' Summarise counts per populated cell
#' @noRd
.cell_counts <- function(data, groups, sep = " : ") {
  cell <- .cell_vector(data, groups, sep = sep)
  tab <- table(cell)
  data.frame(cell = names(tab), n = as.integer(tab),
             stringsAsFactors = FALSE, row.names = NULL)
}

#' Normalise a family argument the way stats::glm() accepts it
#'
#' \code{glm()} accepts a character string, a function or a family object.
#' Anything downstream that reads \code{family$family} needs the object.
#'
#' @noRd
.as_family <- function(family) {
  if (is.character(family)) {
    if (length(family) != 1L) {
      .stopf("`family` must be a single family name, a family function, or a family object.")
    }
    fn <- tryCatch(get(family, mode = "function", envir = parent.frame(2L)),
                   error = function(e) NULL)
    if (is.null(fn)) {
      .stopf("`family` = \"%s\" is not a known family. See ?family.", family)
    }
    family <- fn()
  }
  if (is.function(family)) family <- family()
  if (!inherits(family, "family")) {
    .stopf("`family` must be a single family name, a family function, or a family object; got %s.",
           paste(class(family), collapse = "/"))
  }
  family
}

#' Emit a message only when verbose is TRUE
#'
#' Uses message() so that suppressMessages() works, unlike cat().
#' @noRd
.say <- function(verbose, fmt, ...) {
  if (isTRUE(verbose)) message(sprintf(fmt, ...))
  invisible(NULL)
}

#' Check that a response is something the requested family can model
#'
#' \code{anova_glm()} accepts any family, so the response may legitimately be a
#' factor (binomial) or a matrix of counts. What it must never be is a type the
#' family cannot use at all, or one the plots will choke on later.
#' @noRd
.check_response_for_family <- function(data, response, family) {
  y <- data[[response]]
  fam <- family$family
  # A two-column successes/failures matrix is a legitimate binomial response
  # for stats::glm(), but the package works one column at a time --
  # .prepare_frame() subsets it as a vector -- so it fails downstream with
  # "(subscript) logical subscript too long". Refuse it here, where the message
  # can say what to do instead. A ONE-column matrix (what scale() returns, and
  # a common way to build a column) is fine: every other function in the
  # package accepts it, and .prepare_frame() drops the dimension.
  if (is.matrix(y) && ncol(y) > 1L) {
    .stopf("Response `%s` is a %d-column matrix. A successes/failures matrix response is not supported; supply the outcome as one column per row, with `weights` naming the trial counts if the data are aggregated.",
           response, ncol(y))
  }
  if (fam == "binomial") {
    ok <- is.numeric(y) || is.logical(y) || is.factor(y) || is.character(y)
    if (!ok) {
      .stopf("Response `%s` must be numeric, logical or a factor for a binomial family; it is %s.",
             response, paste(class(y), collapse = "/"))
    }
    return(invisible(TRUE))
  }
  if (!is.numeric(y) && !is.logical(y)) {
    .stopf("Response `%s` must be numeric for a %s family; it is %s. For a two-level factor use family = \"binomial\", or see anova_bin().",
           response, fam, paste(class(y), collapse = "/"))
  }
  if (fam == "poisson" && is.numeric(y) &&
      any(y[is.finite(y)] < 0)) {
    .stopf("Response `%s` contains negative values, which a Poisson family cannot model.",
           response)
  }
  invisible(TRUE)
}

#' Abbreviate a long vector of labels for an error message
#' @noRd
.abbrev <- function(x, n = 8L) {
  x <- as.character(x)
  if (length(x) <= n) return(paste(x, collapse = ", "))
  sprintf("%s and %d more", paste(utils::head(x, n), collapse = ", "),
          length(x) - n)
}

#' cbind a labelling column without colliding with what is already there
#' @noRd
.cbind_label <- function(df, name, value) {
  df <- as.data.frame(df)
  nm <- .safe_name(df, name)
  out <- cbind(stats::setNames(data.frame(value, stringsAsFactors = FALSE), nm),
               df)
  row.names(out) <- NULL
  out
}
