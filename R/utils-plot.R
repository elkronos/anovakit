# Shared plotting helpers ----------------------------------------------------
#
# Every plot in the package is built here, with tidy-evaluation aesthetics
# (.data[[ ]]) rather than the deprecated aes_string(), and returned rather than
# printed.

#' Rotate x labels only when they are long or numerous
#' @noRd
.gg_theme_auto <- function(labels) {
  needs_rotation <- length(labels) > 6L || any(nchar(as.character(labels)) > 8L)
  ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(face = "bold"),
      axis.text.x = if (needs_rotation) {
        ggplot2::element_text(angle = 45, hjust = 1)
      } else {
        ggplot2::element_text()
      },
      legend.position = "bottom"
    )
}

#' Q-Q plot of a numeric vector
#' @noRd
.plot_qq <- function(x, title = "Normal Q-Q plot of residuals",
                     subtitle = NULL) {
  df <- data.frame(value = as.numeric(x))
  df <- df[stats::complete.cases(df), , drop = FALSE]
  ggplot2::ggplot(df, ggplot2::aes(sample = .data[["value"]])) +
    ggplot2::stat_qq(alpha = 0.7) +
    ggplot2::stat_qq_line(linetype = "dashed") +
    ggplot2::labs(title = title, subtitle = subtitle,
                  x = "Theoretical quantiles", y = "Sample quantiles") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}

#' Residuals against fitted values
#' @noRd
.plot_resid_fitted <- function(fitted, residuals,
                               title = "Residuals vs fitted",
                               ylab = "Residuals") {
  df <- data.frame(fitted = as.numeric(fitted),
                   residuals = as.numeric(residuals))
  df <- df[stats::complete.cases(df), , drop = FALSE]
  ggplot2::ggplot(df, ggplot2::aes(x = .data[["fitted"]], y = .data[["residuals"]])) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::labs(title = title, x = "Fitted values", y = ylab) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}

#' Residuals against observation index
#' @noRd
.plot_resid_index <- function(residuals, title = "Residuals vs index") {
  df <- data.frame(index = seq_along(residuals),
                   residuals = as.numeric(residuals))
  ggplot2::ggplot(df, ggplot2::aes(x = .data[["index"]], y = .data[["residuals"]])) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::labs(title = title, x = "Observation", y = "Residuals") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}

#' Group means with confidence intervals
#'
#' @param stats data.frame with \code{group}, \code{mean}, \code{conf_low},
#'   \code{conf_high} and optionally \code{n}.
#' @noRd
.plot_means <- function(stats_df, response, conf_level, xlab = "Group",
                        title = NULL, ylab = NULL) {
  df <- stats_df
  if ("n" %in% names(df)) {
    df$label <- sprintf("%s\n(n = %d)", as.character(df$group), df$n)
  } else {
    df$label <- as.character(df$group)
  }
  df$label <- factor(df$label, levels = df$label[order(df$group)])
  ggplot2::ggplot(df, ggplot2::aes(x = .data[["label"]], y = .data[["mean"]],
                                   fill = .data[["label"]])) +
    ggplot2::geom_col(colour = "grey25", show.legend = FALSE) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data[["conf_low"]], ymax = .data[["conf_high"]]),
      width = 0.2) +
    ggplot2::labs(
      title = title %||% sprintf("Group means with %s%% confidence intervals",
                                 .pct(conf_level)),
      x = xlab, y = ylab %||% sprintf("Mean %s", response)) +
    .gg_theme_auto(as.character(df$group))
}

#' Boxplot of a response by cell, labelled with the cell's own count
#'
#' The counts shown are the counts of the cells actually tested, not the
#' marginal counts of the first grouping variable.
#' @noRd
.plot_box <- function(data, response, cell, title = NULL, xlab = "Group") {
  df <- data.frame(value = data[[response]],
                   cell  = droplevels(as.factor(data[[cell]])),
                   stringsAsFactors = FALSE)
  counts <- table(df$cell)
  med <- tapply(df$value, df$cell, stats::median, na.rm = TRUE)
  # Look levels up by position, not by name: a level that is the empty string
  # cannot be used as a name, and would come out as an NA label.
  ord <- order(med)
  lab_map <- sprintf("%s\n(n = %d)", levels(df$cell), as.integer(counts))
  df$label <- factor(lab_map[as.integer(df$cell)], levels = lab_map[ord])
  ggplot2::ggplot(df, ggplot2::aes(x = .data[["label"]], y = .data[["value"]],
                                   fill = .data[["label"]])) +
    ggplot2::geom_boxplot(outlier.alpha = 0.6, show.legend = FALSE) +
    ggplot2::labs(title = title %||% "Distribution of response by group",
                  subtitle = "Groups ordered by median",
                  x = xlab, y = response) +
    .gg_theme_auto(names(counts))
}

#' Estimated marginal means, with every factor in the grid mapped
#'
#' The first factor goes on the x axis; any remaining factors are mapped to
#' colour and to the line grouping, so a two-factor grid produces one line per
#' level of the second factor rather than a single line that doubles back.
#'
#' @param emm data.frame from \code{as.data.frame(emmeans_object)}.
#' @param factors names of the factor columns in \code{emm}.
#' @param estimate name of the estimate column. When a grouping column has
#'   displaced it the real column is \code{emm_estimate}, so both are tried:
#'   reading the wrong one would plot the factor codes and say nothing.
#' @param lower,upper names of the interval columns, or \code{NULL}.
#' @noRd
.plot_emmeans <- function(emm, factors, estimate = "estimate",
                          lower = "conf_low", upper = "conf_high",
                          conf_level = 0.95, title = NULL, ylab = NULL) {
  term_col <- attr(emm, "term_col")
  df <- as.data.frame(emm)
  factors <- intersect(factors, names(df))
  estimate <- .emm_col(df, estimate, factors)
  lower <- .emm_col(df, lower, factors)
  upper <- .emm_col(df, upper, factors)
  if (length(factors) == 0L || is.null(estimate)) return(NULL)
  has_ci <- !is.null(lower) && !is.null(upper)
  subtitle <- if (has_ci) {
    sprintf("Error bars are %s%% confidence intervals", .pct(conf_level))
  } else NULL

  # Marginal means of an additive model, one factor at a time: one panel per
  # factor, each showing that factor's levels.
  if (!is.null(term_col) && term_col %in% names(df)) {
    lvl <- rep(NA_character_, nrow(df))
    for (g in factors) {
      v <- as.character(df[[g]])
      lvl[!is.na(v)] <- v[!is.na(v)]
    }
    df$.level <- factor(lvl, levels = unique(lvl))
    df$.term <- factor(df[[term_col]], levels = unique(df[[term_col]]))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[".level"]], y = .data[[estimate]])) +
      ggplot2::geom_point(size = 2.5)
    if (has_ci) {
      p <- p + ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data[[lower]], ymax = .data[[upper]]), width = 0.15)
    }
    return(p +
      ggplot2::facet_wrap(ggplot2::vars(.data[[".term"]]), scales = "free_x") +
      ggplot2::labs(title = title %||% "Estimated marginal means",
                    subtitle = subtitle, x = NULL,
                    y = ylab %||% "Estimated marginal mean") +
      .gg_theme_auto(levels(df$.level)))
  }
  xvar <- factors[1L]
  df[[xvar]] <- as.factor(df[[xvar]])

  if (length(factors) > 1L) {
    df$.series <- droplevels(interaction(df[factors[-1L]], drop = TRUE, sep = " : "))
    p <- ggplot2::ggplot(df, ggplot2::aes(
      x = .data[[xvar]], y = .data[[estimate]],
      colour = .data[[".series"]], group = .data[[".series"]]))
    lab_colour <- paste(factors[-1L], collapse = " : ")
  } else {
    df$.series <- factor("all")
    p <- ggplot2::ggplot(df, ggplot2::aes(
      x = .data[[xvar]], y = .data[[estimate]], group = 1L))
    lab_colour <- NULL
  }

  p <- p +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2.5)
  if (has_ci) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data[[lower]], ymax = .data[[upper]]),
      width = 0.15)
  }
  p +
    ggplot2::labs(
      title = title %||% "Estimated marginal means",
      subtitle = subtitle,
      x = xvar, y = ylab %||% "Estimated marginal mean", colour = lab_colour) +
    .gg_theme_auto(levels(df[[xvar]]))
}

#' Stacked proportion plot for a binary response
#' @noRd
.plot_proportions <- function(prop_table, cell, response, success_level,
                              xlab = "Group") {
  df <- prop_table
  df[[cell]] <- as.factor(df[[cell]])
  df[[response]] <- as.factor(df[[response]])
  ggplot2::ggplot(df, ggplot2::aes(x = .data[[cell]], y = .data[["proportion"]],
                                   fill = .data[[response]])) +
    ggplot2::geom_col(colour = "grey25") +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.0f%%", 100 * .data[["proportion"]])),
      position = ggplot2::position_stack(vjust = 0.5), size = 3) +
    ggplot2::scale_y_continuous(labels = function(x) sprintf("%.0f%%", 100 * x)) +
    ggplot2::labs(
      title = sprintf("Proportion of %s by group", response),
      subtitle = sprintf("Modelled outcome: %s = %s", response, success_level),
      x = xlab, y = "Proportion", fill = response) +
    .gg_theme_auto(levels(df[[cell]]))
}

#' Scatterplot of response against covariate, by group
#' @noRd
.plot_covariate <- function(data, response, covariate, cell, conf_level = 0.95,
                            xlab = covariate) {
  data <- data[, unique(c(response, covariate, cell)), drop = FALSE]
  ggplot2::ggplot(data, ggplot2::aes(x = .data[[covariate]],
                                     y = .data[[response]],
                                     colour = .data[[cell]])) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                         level = conf_level) +
    ggplot2::labs(title = sprintf("%s against %s by group", response, covariate),
                  subtitle = sprintf("Bands are %s%% confidence intervals for each group's line",
                                     .pct(conf_level)),
                  x = xlab, y = response, colour = "Group") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   legend.position = "bottom")
}

#' Canonical discriminant scatterplot
#' @noRd
.plot_canonical <- function(scores, group, group_name = "Group") {
  df <- as.data.frame(scores)
  df$.group <- droplevels(as.factor(group))
  if (ncol(scores) == 1L) {
    return(
      ggplot2::ggplot(df, ggplot2::aes(x = .data[[".group"]], y = .data[["Can1"]],
                                       fill = .data[[".group"]])) +
        ggplot2::geom_boxplot(show.legend = FALSE) +
        ggplot2::labs(title = "Canonical discriminant scores",
                      subtitle = "Only one discriminant axis is estimable",
                      x = group_name, y = "Can1") +
        .gg_theme_auto(levels(df$.group))
    )
  }
  centroids <- stats::aggregate(df[, c("Can1", "Can2")],
                                by = list(.group = df$.group), FUN = mean)
  ggplot2::ggplot(df, ggplot2::aes(x = .data[["Can1"]], y = .data[["Can2"]],
                                   colour = .data[[".group"]])) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::geom_point(data = centroids, size = 5, shape = 4, stroke = 1.5) +
    ggplot2::labs(title = "Canonical discriminant analysis",
                  subtitle = "Crosses mark group centroids",
                  x = "Can1", y = "Can2", colour = group_name) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   legend.position = "bottom")
}

#' Null-coalescing operator
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Find the column holding a statistic, allowing for a displaced name
#'
#' When a grouping column is called \code{estimate} or \code{conf_low}, the
#' emmeans statistic of that name is parked under an \code{emm_} prefix. A
#' plot that asks for the bare name would silently draw the factor codes.
#'
#' @param df the table.
#' @param target the wanted statistic, or \code{NULL}.
#' @param protect column names that are grouping variables, not statistics.
#' @return the column name to use, or \code{NULL}.
#' @noRd
.emm_col <- function(df, target, protect = character(0)) {
  if (is.null(target)) return(NULL)
  prefixed <- paste0("emm_", target)
  # A grouping column may itself be called emm_estimate; it is never the
  # statistic, so only a prefixed name that is not a grouping column counts.
  if (prefixed %in% names(df) && !prefixed %in% protect) return(prefixed)
  if (target %in% names(df) && !target %in% protect) return(target)
  NULL
}
