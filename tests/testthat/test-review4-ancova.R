# Detailed regression tests: run locally and on CI (NOT_CRAN=true), skipped
# on CRAN to keep its check time short. The core suite runs everywhere.
skip_on_cran()

# Regression tests for the round-4 review items in the ANCOVA area. Each is
# checked against an independent reference (base R, car, emmeans or ggplot2's
# own layer data).

# output-F1: summary() prints the homogeneity-of-slopes test once ---------------

test_that("output-F1: summary() prints the homogeneity-of-slopes test once", {
  d <- withr::with_seed(1, data.frame(
    g = rep(c("a", "b", "c"), each = 10),
    y = stats::rnorm(30) + rep(c(0, 1, 5), each = 10),
    x = stats::rnorm(30)))
  f <- anova_ancova(d, "y", "g", "x", plots = FALSE)
  expect_identical(f$slopes_test, f$assumptions$slopes)
  # Reference: the nested-model F test of the covariate-by-group interaction
  ref <- stats::anova(stats::lm(y ~ g + x, data = d),
                      stats::lm(y ~ g * x, data = d))
  expect_equal(f$slopes_test$statistic, ref$F[2])
  out <- capture.output(summary(f))
  expect_identical(sum(grepl("additive vs covariate-by-group", out, fixed = TRUE)), 1L)
  # It is still printed (under the assumption checks), with the reference F
  expect_true(any(grepl(format(round(ref$F[2], 3), nsmall = 3), out, fixed = TRUE)))

  # anova_rm()'s sphericity table is still printed once as well
  if (requireNamespace("afex", quietly = TRUE)) {
    long <- withr::with_seed(2, data.frame(
      id = factor(rep(1:8, each = 3)), cond = factor(rep(c("a", "b", "c"), 8)),
      y = stats::rnorm(24)))
    frm <- anova_rm(long, "y", "id", within = "cond", plots = FALSE)
    expect_identical(frm$sphericity, frm$assumptions$sphericity)
    out_rm <- capture.output(summary(frm))
    expect_identical(sum(grepl("mauchly_w", out_rm, fixed = TRUE)), 1L)
  }
})

# ancova-F1: a saturated forced-interaction model prints nothing ----------------

test_that("ancova-F1: force_interaction = TRUE with no residual df is silent", {
  d <- data.frame(g = rep(c("a", "b", "c"), each = 2),
                  x = c(1, 2, 1.5, 3, 2, 2.5),
                  y = c(2, 3.1, 2.2, 4.8, 3.9, 4.1))
  # Reference: the interaction model leaves no residual degrees of freedom
  ref <- stats::lm(y ~ g * x, data = d)
  expect_identical(stats::df.residual(ref), 0L)
  expect_silent(f <- anova_ancova(d, "y", "g", "x", force_interaction = TRUE,
                                  plots = FALSE))
  # The means at the covariate mean are still the fitted lines' values there
  xm <- mean(d$x)
  expect_equal(f$emmeans$estimate,
               unname(stats::predict(ref, data.frame(g = c("a", "b", "c"), x = xm))))
  expect_true(all(is.na(f$emmeans$conf_low)))

  d2 <- data.frame(g = factor(c("a", "a", "b", "b")), x = c(1, 2, 1, 3),
                   y = c(1, 2, 2, 5))
  expect_silent(anova_ancova(d2, "y", "g", "x", plots = FALSE,
                             force_interaction = TRUE, vcov_type = "HC3"))
  # The helper itself: the cells it names are emmeans' non-estimable ones
  # (none here), with no warning on the way
  expect_silent(cells <- anovakit:::.nonest_cells(ref, "g"))
  expect_identical(cells, character(0))
})

# ancova-R3: non-estimable means are really left out of $emmeans -------------

test_that("ancova-R3: non-estimable group means are left out of $emmeans, as the note says", {
  d <- withr::with_seed(3, {
    g <- factor(rep(c("a", "b", "c"), each = 20))
    x <- ifelse(g == "a", stats::rnorm(60, 5, 1), ifelse(g == "b", 2, 8))
    data.frame(g = g, x = x, y = stats::rnorm(60))
  })
  ref <- as.data.frame(summary(suppressMessages(emmeans::emmeans(
    stats::lm(y ~ g * I(x - mean(x)), data = d), ~g))))
  est_ref <- ref[!is.na(ref$emmean), ]
  expect_identical(as.character(est_ref$g), "a")
  for (args in list(list(plots = FALSE), list(plots = TRUE),
                    list(plots = FALSE, vcov_type = "HC3"))) {
    f <- do.call(anova_ancova, c(list(d, "y", "g", "x", force_interaction = TRUE),
                                 args))
    expect_identical(nrow(f$emmeans), 1L)
    expect_false(anyNA(f$emmeans$g))
    expect_identical(as.character(f$emmeans$g), "a")
    expect_equal(f$emmeans$estimate, est_ref$emmean)
    expect_true(any(grepl("b, c are not estimable", f$notes, fixed = TRUE)))
    expect_true(any(grepl("left out of $emmeans", f$notes, fixed = TRUE)))
    expect_false(any(grepl("extrapolation", f$notes, fixed = TRUE)))
  }
  f <- anova_ancova(d, "y", "g", "x", force_interaction = TRUE, plots = FALSE)
  expect_equal(f$emmeans$se, est_ref$SE)
  expect_equal(f$emmeans$df, est_ref$df)
})

# ancova-manova-merge-F1: Levene leaves out cells the model fits exactly ------

test_that("ancova-manova-merge-F1: Levene's test leaves out exact-fit cells", {
  m <- mtcars
  m$cyl <- factor(m$cyl)
  m$am <- factor(m$am)
  f <- anova_ancova(m, "mpg", c("cyl", "am"), c("wt", "hp"),
                    force_interaction = TRUE, plots = FALSE)
  # Reference: car on the interaction model's residuals, without the cells of
  # one or two observations, those it fits exactly (residuals zero) and the
  # one (6 : 0, four observations, three coefficients) it leaves with a single
  # residual degree of freedom
  ref_m <- stats::lm(mpg ~ cyl * am * (wt + hp), data = m)
  r <- stats::residuals(ref_m)
  cell <- interaction(m$cyl, m$am, sep = " : ", drop = TRUE)
  n <- table(cell)
  zero <- tapply(abs(r), cell, max) < 1e-8
  expect_identical(sort(names(which(zero & n >= 3))), c("4 : 0", "6 : 1"))
  drop <- c(names(n)[n < 3 | zero], "6 : 0")
  keep <- !(as.character(cell) %in% drop)
  ref <- car::leveneTest(r[keep] ~ droplevels(cell[keep]))
  expect_equal(f$assumptions$levene$statistic, ref[["F value"]][1])
  expect_equal(f$assumptions$levene$p_value, ref[["Pr(>F)"]][1])
  expect_equal(f$assumptions$levene$df1, 1)
  expect_true(any(grepl("leaves out 2 group(s) that the fitted model fits exactly by construction (4 : 0, 6 : 1)",
                        f$notes, fixed = TRUE)))
  expect_true(any(grepl("leaves out 1 group(s) that the fitted model leaves with fewer than 2.5 residual degrees of freedom (6 : 0 with 1;",
                        f$notes, fixed = TRUE)))
  expect_false(any(grepl("Levene's test rejects", f$notes, fixed = TRUE)))

  # Two exact-fit cells and only one other: the test is skipped, with a reason
  d <- withr::with_seed(5, {
    g <- factor(rep(c("a", "b", "c"), c(3, 3, 12)))
    d <- data.frame(g = g, x1 = stats::rnorm(18), x2 = stats::rnorm(18))
    d$y <- d$x1 + stats::rnorm(18)
    d
  })
  f2 <- anova_ancova(d, "y", "g", c("x1", "x2"), force_interaction = TRUE,
                     plots = FALSE, posthoc = FALSE)
  r2 <- stats::residuals(stats::lm(y ~ g * (x1 + x2), data = d))
  expect_true(all(abs(r2[d$g %in% c("a", "b")]) < 1e-8))
  expect_null(f2$assumptions$levene)
  expect_true(any(grepl("Levene's test was skipped: it needs at least two groups with three or more observations that the model does not fit exactly, and only one group has that many",
                        f2$notes, fixed = TRUE)))
})

# ancova-manova-merge-F2: covariate plot with a column named .lwr or .upr -----

test_that("ancova-manova-merge-F2: the covariate plot builds with columns named .lwr or .upr", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  an <- MASS::anorexia
  for (case in c("response", "covariate")) {
    d <- an
    if (case == "response") {
      names(d)[names(d) == "Postwt"] <- ".upr"
      resp <- ".upr"; cv <- "Prewt"
    } else {
      names(d)[names(d) == "Prewt"] <- ".lwr"
      resp <- "Postwt"; cv <- ".lwr"
    }
    f <- anova_ancova(d, resp, "Treat", cv)
    expect_silent(b <- ggplot2::ggplot_build(f$plots$covariate))
    # Reference: the bands are lm()'s confidence intervals for the fitted model
    form <- if (is.null(f$simple_slopes)) stats::reformulate(c("Treat", cv), resp)
            else stats::reformulate(paste0("Treat * `", cv, "`"), resp)
    ref <- stats::lm(form, data = d)
    rib <- ggplot2::layer_data(f$plots$covariate, 2)
    lev <- levels(d$Treat)
    for (k in seq_along(lev)) {
      s <- rib[rib$group == k, ]
      nd <- data.frame(s$x, factor(lev[k], levels = lev))
      names(nd) <- c(cv, "Treat")
      p <- stats::predict(ref, nd, interval = "confidence")
      expect_equal(s$ymin, unname(p[, "lwr"]))
      expect_equal(s$ymax, unname(p[, "upr"]))
    }
  }
})

# ancova-R1: a cell with a constant response stays in Levene's test ----------

test_that("ancova-R1: a large cell whose residuals are zero because its response is constant stays in Levene's test", {
  d <- withr::with_seed(7, {
    g <- factor(rep(c("ctrl", "low", "high"), each = 12))
    x <- stats::rnorm(36, 10, 2)
    y <- ifelse(g == "high", 0, 5 + 0.5 * x + stats::rnorm(36, 0, 2))
    data.frame(g, x, y)
  })
  f <- anova_ancova(d, "y", "g", "x", force_interaction = TRUE, plots = FALSE)
  r <- stats::residuals(stats::lm(y ~ g * x, data = d))
  expect_true(all(abs(r[d$g == "high"]) < 1e-8))
  ref <- car::leveneTest(r ~ d$g)
  expect_equal(f$assumptions$levene$df1, 2)
  expect_equal(f$assumptions$levene$statistic, ref[["F value"]][1])
  expect_equal(f$assumptions$levene$p_value, ref[["Pr(>F)"]][1])
  expect_lt(f$assumptions$levene$p_value, 0.05)
  expect_true(any(grepl("Levene's test rejects", f$notes, fixed = TRUE)))
  expect_false(any(grepl("leaves out", f$notes, fixed = TRUE)))
})

# ancova-C1: cells with one residual degree of freedom are left out ----------

test_that("ancova-C1: Levene's test keeps its level when small cells have one residual df", {
  # Cells of four observations under y ~ g * (x1 + x2) have three coefficients
  # of their own, so one residual degree of freedom; left in, they make the
  # test reject null data about 20% of the time. A simulation under
  # homoscedastic data is the reference: the rejection rate must be near 5%.
  pv <- withr::with_seed(7, vapply(seq_len(120), function(i) {
    g <- factor(rep(c("a", "b", "c", "d"), c(4, 4, 15, 15)))
    d <- data.frame(g = g, x1 = stats::rnorm(38), x2 = stats::rnorm(38))
    d$y <- stats::rnorm(38)
    f <- anova_ancova(d, "y", "g", c("x1", "x2"), force_interaction = TRUE,
                      plots = FALSE, posthoc = FALSE)
    if (i == 1L) {
      # The statistic is car's on the two cells that are kept
      r <- stats::residuals(stats::lm(y ~ g * (x1 + x2), data = d))
      k <- d$g %in% c("c", "d")
      ref <- car::leveneTest(r[k] ~ droplevels(d$g[k]))
      expect_equal(f$assumptions$levene$p_value, ref[["Pr(>F)"]][1])
      expect_true(any(grepl("fewer than 2.5 residual degrees of freedom (a with 1, b with 1;",
                            f$notes, fixed = TRUE)))
    }
    f$assumptions$levene$p_value
  }, numeric(1)))
  # Binomial(120, 0.05) has over 99.9% of its mass below 15 rejections; before
  # the fix about 24 were expected
  expect_lt(sum(pv < 0.05), 15)
})
