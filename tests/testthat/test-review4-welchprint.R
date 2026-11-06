# Review round 4, area welchprint: anova_welch() omnibus test at extreme
# scales, its variance-ratio and bias-correction notes, the formatting of
# p-value columns in print()/summary(), and matrix or data-frame grouping
# columns in anova_welch() and anova_kw().

test_that("welchkw-F1: Welch's F is computed at any response scale, and no note calls a NaN 'reported'", {
  pg <- datasets::PlantGrowth
  ref <- stats::oneway.test(weight ~ group, data = pg, var.equal = FALSE)

  # Scaling by a power of two is exact, so the result must be identical.
  e <- pg
  e$weight <- e$weight * 2^-520
  fit <- anova_welch(e, "weight", "group", plots = FALSE)
  expect_identical(fit$anova$statistic, unname(ref$statistic))
  expect_identical(fit$anova$den_df, unname(ref$parameter[2L]))
  expect_identical(fit$anova$p_value, ref$p.value)

  # At 1e-155 the group variances are subnormal and oneway.test() itself
  # returns NaN; the multiplication is inexact, so compare with a tolerance.
  e <- pg
  e$weight <- e$weight * 1e-155
  expect_true(is.nan(suppressWarnings(
    stats::oneway.test(weight ~ group, data = e)$statistic)))
  fit <- anova_welch(e, "weight", "group", plots = FALSE)
  expect_equal(fit$anova$statistic, unname(ref$statistic), tolerance = 1e-8)
  expect_equal(fit$anova$den_df, unname(ref$parameter[2L]), tolerance = 1e-8)
  expect_equal(fit$anova$p_value, ref$p.value, tolerance = 1e-8)
  expect_length(fit$notes, 0L)

  # A noise pair under 5% of F: the "is reported" note must sit beside a
  # finite statistic, the one at an ordinary scale.
  d4 <- data.frame(g = rep(c("a", "b", "c"), each = 4),
                   y = c(1e9 + c(0, 1, 2, 3) * 1e-6, 1e9 + c(1, 2, 3, 4) * 1e-6,
                         1e9 + c(5, 6, 7, 8)))
  e <- d4
  e$y <- e$y * 2^-500
  ref4 <- stats::oneway.test(y ~ g, data = d4, var.equal = FALSE)
  fit <- anova_welch(e, "y", "g", plots = FALSE)
  expect_identical(fit$anova$statistic, unname(ref4$statistic))
  expect_true(any(grepl("so the omnibus test is reported", fit$notes)))
  e <- d4
  e$y <- e$y * 1e-150
  fit <- anova_welch(e, "y", "g", plots = FALSE)
  expect_true(is.finite(fit$anova$statistic))
  expect_equal(fit$anova$statistic, unname(ref4$statistic), tolerance = 0.01)
  expect_true(any(grepl("so the omnibus test is reported", fit$notes)))

  # Variances spanning more than double precision can weight together: NA
  # with a note, not a silent NaN.
  d <- data.frame(g = rep(c("a", "b", "c"), each = 4),
                  y = c(1e150 * (1 + c(0, 1, 2, 3) * 1e-5),
                        1e-150 * (1 + c(1, 2, 4, 3) * 1e-10),
                        1e-150 * (1 + c(5, 1, 3, 2) * 1e-10)))
  expect_true(is.nan(suppressWarnings(stats::oneway.test(y ~ g, d)$statistic)))
  fit <- anova_welch(d, "y", "g", plots = FALSE)
  expect_true(is.na(fit$anova$statistic) && !is.nan(fit$anova$statistic))
  expect_true(is.na(fit$anova$p_value))
  expect_true(any(grepl("could not be computed in double precision", fit$notes)))
  expect_false(any(grepl("is reported", fit$notes)))
})

test_that("welchkw-F2: the variance-ratio note prints a large ratio readably", {
  d <- data.frame(g = rep(c("a", "b", "c"), each = 3),
                  y = 1e6 + c(0, 1e-9, 2e-9, 5, 5 + 1e-9, 5 + 2e-9, 1, 3, 7))
  fit <- anova_welch(d, "y", "g", plots = FALSE)
  v <- tapply(d$y, d$g, stats::var)
  ratio <- max(v) / min(v)
  expect_gt(ratio, 1e15)
  note <- grep("^Largest group variance", fit$notes, value = TRUE)
  expect_length(note, 1L)
  shown <- as.numeric(sub("^Largest group variance is (\\S+) times.*", "\\1", note))
  expect_equal(shown, ratio, tolerance = 0.005)
  expect_false(grepl("[0-9]{8,}", note))

  # An ordinary ratio keeps three significant digits.
  set.seed(4)
  d <- data.frame(g = rep(c("a", "b"), each = 20),
                  y = c(stats::rnorm(20, 0, 1), stats::rnorm(20, 0, 3)))
  v <- tapply(d$y, d$g, stats::var)
  note <- grep("^Largest group variance",
               anova_welch(d, "y", "g", plots = FALSE)$notes, value = TRUE)
  expect_true(grepl(sprintf("is %s times", sprintf("%.3g", max(v) / min(v))),
                    note, fixed = TRUE))
})

test_that("welch-kw-print-inputs-F1: p-value columns print each value at `digits` significant digits, unpadded", {
  aq <- datasets::airquality
  aq$Month <- factor(aq$Month)
  fit <- anova_kw(aq, "Ozone", "Month", diagnostics = TRUE, plots = FALSE)
  p <- fit$posthoc$p_adjusted
  shown <- anovakit:::.round_df(fit$posthoc, 4L)$p_adjusted
  big <- !is.na(p) & p >= 1e-4
  expect_identical(shown[big], sprintf("%.4g", p[big]))
  pn <- fit$assumptions$normality$p_value
  shown_n <- anovakit:::.round_df(fit$assumptions$normality, 4L)$p_value
  ok <- !is.na(pn) & pn >= 1e-4
  expect_identical(shown_n[ok], sprintf("%.4g", pn[ok]))

  # A Bonferroni-capped 1 is "1", not "1.0000000", and below eps is "< eps".
  df <- data.frame(p_value = c(1, 0.9679927, 0.0496380123, 1e-20, NA))
  expect_identical(anovakit:::.round_df(df, 4L)$p_value,
                   c("1", "0.968", "0.04964",
                     format.pval(1e-20, digits = 4L), "NA"))
  out <- utils::capture.output(summary(
    anova_welch(datasets::PlantGrowth, "weight", "group",
                adjust = "bonferroni", plots = FALSE)))
  expect_false(any(grepl("1\\.0000", out)))
})

test_that("welch-kw-print-inputs-F2: the Hedges factor note agrees with the numbers it prints at 4 df", {
  d <- data.frame(g = rep(c("a", "b"), each = 3), y = c(1, 2, 3, 4, 5, 6))
  fit <- anova_welch(d, "y", "g", plots = FALSE)
  nu <- 4
  J <- gamma(nu / 2) / (sqrt(nu / 2) * gamma((nu - 1) / 2))
  expect_lt(J, 0.8)
  # Cohen's d on sqrt((1 + 1) / 2) = 1 is -3.
  expect_equal(fit$effect_sizes$hedges_g, -3 * J, tolerance = 1e-12)
  note <- grep("bias correction multiplied", fit$notes, value = TRUE)
  expect_length(note, 1L)
  expect_true(grepl(sprintf("(4 df, factor %.3f)", J), note, fixed = TRUE))
  expect_true(grepl("about 4 degrees of freedom or fewer", note, fixed = TRUE))
  expect_false(grepl("fewer than 4", note, fixed = TRUE))
  # Whatever the factor, it is never printed as 0.8 beside "less than 0.8".
  expect_false(grepl("factor 0.8)", note, fixed = TRUE))
})

test_that("welch-kw-print-inputs-F3: matrix and data-frame grouping columns are refused or flattened by name", {
  set.seed(2)
  d <- data.frame(y = stats::rnorm(24))
  g <- rep(c("a", "b", "c"), each = 8)
  d$g2 <- cbind(g, rep(c("x", "y"), 12))
  expect_error(anova_kw(d, "y", "g2", plots = FALSE),
               "Grouping column `g2` is a 2-column matrix; supply a single column.",
               fixed = TRUE)
  expect_error(anova_welch(d, "y", "g2", plots = FALSE),
               "Grouping column `g2` is a 2-column matrix; supply a single column.",
               fixed = TRUE)
  d$g2 <- NULL
  d$gd <- data.frame(a = g)
  plain <- data.frame(y = d$y, g = g)
  ref_w <- stats::oneway.test(y ~ g, data = plain, var.equal = FALSE)
  ref_k <- stats::kruskal.test(y ~ g, data = plain)
  fw <- anova_welch(d, "y", "gd", plots = FALSE)
  fk <- anova_kw(d, "y", "gd", plots = FALSE)
  expect_equal(fw$anova$statistic, unname(ref_w$statistic))
  expect_equal(fw$anova$p_value, ref_w$p.value)
  expect_equal(fk$anova$statistic, unname(ref_k$statistic))
  expect_equal(fk$anova$p_value, ref_k$p.value)
  d$gd2 <- data.frame(a = g, b = g)
  expect_error(anova_welch(d, "y", "gd2", plots = FALSE),
               "Grouping column `gd2` is a 2-column data frame", fixed = TRUE)
})

test_that("welchprint-R1: a finite Welch F at the response's own scale is never lost to rescaling", {
  # Group c's variance (~1e-301) is a normal double, so oneway.test() gives a
  # finite F unscaled; dividing the response by 2^16 would make it subnormal
  # and the weights overflow.
  d <- data.frame(g = rep(c("a", "b", "c"), each = 4),
                  y = c(1e5 + c(0, 1, 2, 4), 1e5 + c(3, 5, 6, 9),
                        1e-150 * c(1, 2, 4, 3)))
  ref <- stats::oneway.test(y ~ g, data = d, var.equal = FALSE)
  expect_true(is.finite(ref$statistic))
  fit <- anova_welch(d, "y", "g", plots = FALSE)
  expect_identical(fit$anova$statistic, unname(ref$statistic))
  expect_identical(fit$anova$den_df, unname(ref$parameter[2L]))
  expect_identical(fit$anova$p_value, ref$p.value)
  expect_false(any(grepl("could not be computed", fit$notes)))
  # Hand computation of the variance ratio from the unscaled variances.
  v <- tapply(d$y, d$g, stats::var)
  expect_true(any(grepl(sprintf("Largest group variance is %s times",
                                format(signif(max(v) / min(v), 3))),
                        fit$notes, fixed = TRUE)))

  set.seed(1)
  d <- data.frame(g = rep(c("a", "b", "c"), each = 6),
                  y = c(stats::rnorm(6, 1e5, 1), stats::rnorm(6, 1e5 + 3, 2),
                        stats::rnorm(6, 1e-150, 1e-151)))
  ref <- stats::oneway.test(y ~ g, data = d, var.equal = FALSE)
  fit <- anova_welch(d, "y", "g", plots = FALSE)
  expect_identical(fit$anova$statistic, unname(ref$statistic))
  expect_identical(fit$anova$p_value, ref$p.value)

  # Where only rescaling helps, the result equals oneway.test() on the
  # response multiplied by a power of two (exact), done by hand: here group
  # c's variance is subnormal and the means differ by 1e5.
  d <- data.frame(g = rep(c("a", "b", "c"), each = 4),
                  y = c(1 + c(0, 1, 2, 4), 2 + c(3, 5, 6, 9),
                        1e-160 * c(1, 2, 4, 3)))
  expect_false(is.finite(suppressWarnings(
    stats::oneway.test(y ~ g, data = d, var.equal = FALSE)$statistic)))
  e <- d
  e$y <- e$y * 2^200
  ref <- stats::oneway.test(y ~ g, data = e, var.equal = FALSE)
  expect_true(is.finite(ref$statistic))
  fit <- anova_welch(d, "y", "g", plots = FALSE)
  expect_equal(fit$anova$statistic, unname(ref$statistic), tolerance = 1e-12)
  expect_equal(fit$anova$den_df, unname(ref$parameter[2L]), tolerance = 1e-12)
  expect_equal(fit$anova$p_value, ref$p.value, tolerance = 1e-12)
  expect_false(any(grepl("could not be computed", fit$notes)))
})

test_that("welchprint-C1: every function refuses a multi-column matrix grouping column and flattens a one-column data frame", {
  set.seed(2)
  d <- data.frame(y = stats::rnorm(24), cnt = stats::rpois(24, 3),
                  b = stats::rbinom(24, 1, 0.5), y2 = stats::rnorm(24),
                  x = stats::rnorm(24))
  d$g2 <- cbind(rep(c("a", "b", "c"), each = 8), rep(c("x", "y"), 12))
  d$gd <- data.frame(a = rep(c("a", "b", "c"), each = 8))
  msg <- "Grouping column `g2` is a 2-column matrix; supply a single column."
  expect_error(anova_glm(d, "y", "g2", plots = FALSE), msg, fixed = TRUE)
  expect_error(anova_bin(d, "b", "g2", plots = FALSE), msg, fixed = TRUE)
  expect_error(anova_count(d, "cnt", "g2", plots = FALSE), msg, fixed = TRUE)
  expect_error(anova_manova(d, c("y", "y2"), "g2", plots = FALSE), msg,
               fixed = TRUE)
  expect_error(anova_ancova(d, "y", "g2", "x", plots = FALSE), msg,
               fixed = TRUE)

  # The one-column data frame gives the same test as the plain column, and
  # that agrees with car::Anova() on lm().
  d0 <- data.frame(y = d$y, x = d$x, gd = d$gd$a)
  fit <- anova_glm(d, "y", "gd", plots = FALSE)
  ref <- car::Anova(stats::lm(y ~ gd, data = d0), type = 2)
  expect_equal(fit$anova$p_value[1L], ref[["Pr(>F)"]][1L], tolerance = 1e-8)
  fit <- anova_ancova(d, "y", "gd", "x", plots = FALSE)
  expect_equal(fit$anova, anova_ancova(d0, "y", "gd", "x", plots = FALSE)$anova)
  ref <- car::Anova(stats::lm(y ~ gd + x, data = transform(d0, x = x - mean(x))),
                    type = 3)
  expect_equal(fit$anova$p_value[fit$anova$term == "gd"],
               ref["gd", "Pr(>F)"], tolerance = 1e-8)

  # anova_rm(): a within-subject factor given as a matrix column.
  r <- data.frame(y = stats::rnorm(24), s = factor(rep(1:8, 3)),
                  w = factor(rep(1:3, each = 8)))
  r$wm <- cbind(as.character(r$w), "z")
  expect_error(anova_rm(r, "y", "s", "wm", plots = FALSE),
               "Within-subject factor `wm` is a 2-column matrix", fixed = TRUE)
})
