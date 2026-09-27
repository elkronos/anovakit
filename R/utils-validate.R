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
.prepare_frame <- function(data, cols, factors = character(0), keep_all = TRUE,
                           keep_infinite = character(0),
                           quiet_factors = character(0)) {
  cols <- unique(cols)
  sub <- data[, cols, drop = FALSE]
  keep <- stats::complete.cases(sub)

  # A factor can carry NA as an explicit level (factor(x, exclude = NULL), or
  # addNA()). complete.cases() sees a valid code and keeps the row, but every
  # model and test downstream treats the value as missing and drops it without
  # counting it. Treat it as missing here, where it is counted and reported.
  n_na_level <- 0L
  for (cn in cols) {
    v <- sub[[cn]]
    if (is.factor(v) && anyNA(levels(v))) {
      bad <- is.na(levels(v))[as.integer(v)] %in% TRUE
      n_na_level <- n_na_level + sum(bad & keep)
      keep <- keep & !bad
    }
  }

  # Inf and -Inf survive complete.cases() but cannot be modelled. Treat them as
  # missing rather than letting them reach a variance or a link function. A
  # rank-based analysis can use them, so a caller may exempt columns.
  n_infinite <- 0L
  for (cn in setdiff(cols, keep_infinite)) {
    v <- sub[[cn]]
    if (is.numeric(v)) {
      bad <- is.finite(v) == FALSE & !is.na(v)
      n_infinite <- n_infinite + sum(bad & keep)
      keep <- keep & !bad
    }
  }
  n_removed <- sum(!keep)
  if (keep_all) {
    # Columns without a name cannot be carried into a model frame: lm() fails
    # on them with "attempt to use zero-length variable name".
    named <- !is.na(names(data)) & nzchar(names(data))
    out_cols <- which(named | names(data) %in% cols)
  } else {
    out_cols <- match(cols, names(data))
  }
  out <- data[keep, out_cols, drop = FALSE]
  row.names(out) <- NULL

  notes <- character(0)
  if (n_removed > 0L) {
    what <- c("missing", if (n_infinite > 0L) "infinite")
    notes <- c(notes, sprintf(
      "Dropped %d row(s) with %s values in: %s.",
      n_removed, paste(what, collapse = " or "),
      paste(cols, collapse = ", ")))
    if (n_na_level > 0L) {
      notes <- c(notes, sprintf(
        "%d of them had NA as an explicit factor level, which is treated as missing. Recode it to a named level to keep those rows.",
        n_na_level))
    }
  }
  if (nrow(out) == 0L) {
    .stopf("No complete cases remain after removing rows with missing values in: %s.",
           paste(cols, collapse = ", "))
  }

  for (f in factors) {
    if (!is.factor(out[[f]])) {
      n_lev <- length(unique(out[[f]]))
      if (is.numeric(out[[f]]) && n_lev > 20L && !f %in% quiet_factors) {
        notes <- c(notes, sprintf(
          "Grouping variable `%s` is numeric with %d distinct values and was converted to a factor. If it is a covariate rather than a group, this is not the function you want.",
          f, n_lev))
      }
      out[[f]] <- .as_group_factor(out[[f]], f)
    }
    out[[f]] <- droplevels(out[[f]])
  }

  list(data = out, n_removed = n_removed, notes = notes)
}

#' Turn a grouping column into a factor whose levels are its distinct values
#'
#' \code{factor()} builds levels from \code{as.character()} of the values. Two
#' distinct values that print the same -- two instants in the repeated hour at
#' the end of daylight saving time, or two doubles equal to 15 significant
#' digits -- would then be merged into one group without a word. Levels are
#' built from the distinct values instead, and labelled with a representation
#' that keeps them apart.
#'
#' Character levels are sorted in the C locale, so the reference level of a
#' character column does not depend on the session's collation. Strings read
#' in a UTF-8 session are often marked with an "unknown" encoding, which the
#' radix sort refuses when they are not ASCII, so a UTF-8 copy is sorted and
#' the original strings are used as the levels.
#' @noRd
.as_group_factor <- function(x, name = "") {
  if (is.character(x)) {
    lev <- unique(x[!is.na(x)])
    lev <- lev[order(enc2utf8(lev), method = "radix")]
    return(factor(x, levels = lev))
  }
  if (is.logical(x)) return(factor(x, levels = c(FALSE, TRUE)[c(FALSE, TRUE) %in% x]))
  u <- sort(unique(x[!is.na(x)]))
  lab <- as.character(u)
  if (anyDuplicated(lab)) {
    lab <- if (inherits(u, "POSIXt")) {
      format(u, "%Y-%m-%d %H:%M:%OS6 %z")
    } else if (is.numeric(u)) {
      format(u, digits = 17L, trim = TRUE)
    } else {
      lab
    }
  }
  if (anyDuplicated(lab)) {
    .stopf("Grouping variable `%s` has distinct values that cannot be told apart once printed (%s). Convert it to character or factor yourself.",
           name, .abbrev(unique(lab[duplicated(lab)])))
  }
  factor(match(x, u), levels = seq_along(u), labels = lab)
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
      .stopf("Grouping variable `%s` has %d level(s) in the rows that are analysed (after dropping rows with missing values, or with zero weight); at least %d are required.",
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
  cell <- droplevels(do.call(interaction, c(f, list(drop = TRUE, sep = sep))))
  # interaction() labels a cell by pasting its levels, so a:"b : c" and
  # "a : b":c get the same label and would be merged into one group. Check the
  # labels against the combinations of level codes they are meant to represent.
  codes <- do.call(paste, c(lapply(f, as.integer), list(sep = "\r")))
  if (length(unique(codes)) != nlevels(cell)) {
    .stopf("The levels of %s produce the same label for different combinations when joined with \"%s\", so distinct groups would be merged. Rename the levels that contain \"%s\".",
           paste(sprintf("`%s`", groups), collapse = ", "), sep, trimws(sep))
  }
  cell
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
#' @param values \code{FALSE} to check only the column's type, which is done
#'   on the input; the levels and values are checked (the default) on the rows
#'   that are analysed, so that a level or value seen only in a dropped row
#'   neither triggers nor hides an error.
#' @noRd
.check_response_for_family <- function(data, response, family, values = TRUE) {
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
  if (fam %in% c("binomial", "quasibinomial")) {
    ok <- is.numeric(y) || is.logical(y) || is.factor(y) || is.character(y)
    if (!ok) {
      .stopf("Response `%s` must be numeric, logical, character or a factor for a %s family; it is %s.",
             response, fam, paste(class(y), collapse = "/"))
    }
    if (values && (is.factor(y) || is.character(y))) {
      lv <- if (is.factor(y)) levels(droplevels(y)) else
        sort(unique(y[!is.na(y)]), method = "radix")
      if (length(lv) != 2L) {
        .stopf("Response `%s` has %d observed level(s); a %s model of a factor response needs exactly two. Recode it to a two-level factor, or 0/1, first.",
               response, length(lv), fam)
      }
    }
    return(invisible(TRUE))
  }
  if (!is.numeric(y) && !is.logical(y)) {
    .stopf("Response `%s` must be numeric for a %s family; it is %s. For a two-level factor use family = \"binomial\", or see anova_bin().",
           response, fam, paste(class(y), collapse = "/"))
  }
  if (values && fam == "poisson" && is.numeric(y) &&
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

#' Refuse a column that plays two roles in one analysis
#'
#' A response that is also a grouping variable, or a subject identifier that is
#' also a between-subjects factor, produces an analysis that is either
#' meaningless or silently different from the one asked for. Every function
#' calls this once, with every role it takes.
#'
#' Two names are reserved: \code{car::Anova()} and \code{summary.manova()} label
#' the error row \code{Residuals} and the intercept \code{(Intercept)}, so a
#' column of either name makes their tables ambiguous.
#'
#' @param roles named list; each element is the column name(s) for one role,
#'   and its name is how that role is described in the message.
#' @noRd
.check_roles <- function(roles) {
  roles <- roles[!vapply(roles, is.null, logical(1))]
  seen <- character(0)
  owner <- character(0)
  for (r in names(roles)) {
    for (col in unique(roles[[r]])) {
      hit <- match(col, seen)
      if (!is.na(hit)) {
        .stopf("Column `%s` is given as both %s and %s; a column can play only one role.",
               col, owner[hit], r)
      }
      seen <- c(seen, col)
      owner <- c(owner, r)
    }
  }
  reserved <- intersect(seen, c("Residuals", "(Intercept)"))
  if (length(reserved) > 0L) {
    .stopf("Column name(s) %s are reserved: the analysis tables use them for the error and intercept rows. Rename the column.",
           paste(sprintf("`%s`", reserved), collapse = ", "))
  }
  invisible(TRUE)
}
