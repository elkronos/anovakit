# Regression tests for the fifth review's anova_welch() findings
# (welchprint-F1, welch-print-groups-F1, welch-print-groups-F2). Each is
# checked against an independent reference: base R's oneway.test() and
# t.test() on the unscaled data (Welch's F and t do not depend on the
# response's units), or a hand computation of the variance ratio on the log
# scale.

test_that("welchprint-F1: a variance ratio that overflows still gets its note", {
  d <- data.frame(g = rep(c("a", "b", "c"), each = 4),
                  y = c(1e150 * c(1, 2, 4, 3), c(1, 3, 2, 5),
                        1e-150 * c(1, 2, 4, 3)))
  fit <- anova_welch(d, "y", "g", plots = FALSE)
  # Reference: the ratio by hand, as a power of ten from the logarithms.
  v <- tapply(d$y, d$g, var)
  expect_false(is.finite(max(v) / min(v)))
  log10_ratio <- log10(max(v)) - log10(min(v))
  expect_equal(log10_ratio, 600, tolerance = 1e-12)
  expect_identical(fit$assumptions$variance_ratio, Inf)
  note <- grep("Largest group variance is", fit$notes, value = TRUE)
  expect_length(note, 1L)
  expect_match(note, "Largest group variance is 1e+600 times the smallest",
               fixed = TRUE)
  # The test itself is unaffected.
  ref <- oneway.test(y ~ g, data = d)
  expect_equal(fit$anova$statistic, unname(ref$statistic), tolerance = 1e-10)

  # A finite ratio keeps its usual formatting.
  d2 <- data.frame(g = rep(c("a", "b"), each = 4),
                   y = c(c(1, 2, 4, 3), 10 * c(1, 3, 2, 5)))
  v2 <- tapply(d2$y, d2$g, var)
  fit2 <- anova_welch(d2, "y", "g", plots = FALSE)
  expect_true(any(grepl(sprintf("Largest group variance is %s times",
                                format(signif(max(v2) / min(v2), 3))),
                        fit2$notes, fixed = TRUE)))
})

test_that("welch-print-groups-F1: no finite but overflowed Welch F at 1e-155", {
  ref <- oneway.test(weight ~ feed, data = chickwts)
  for (s in c(1e-155, 1e-154, 1e-156)) {
    d <- chickwts
    d$weight <- d$weight * s
    # At 1e-155 each weight n / s^2 is finite but their sum is not, which is
    # what made oneway.test() return a finite, wrong F there.
    if (s == 1e-155) {
      w <- with(d, tapply(weight, feed, length) / tapply(weight, feed, var))
      expect_true(all(is.finite(w)))
      expect_false(is.finite(sum(w)))
    }
    fit <- anova_welch(d, "weight", "feed", plots = FALSE)
    expect_equal(fit$anova$statistic, unname(ref$statistic), tolerance = 1e-10)
    expect_equal(fit$anova$num_df, unname(ref$parameter[1L]), tolerance = 1e-10)
    expect_equal(fit$anova$den_df, unname(ref$parameter[2L]), tolerance = 1e-10)
    expect_equal(fit$anova$p_value, ref$p.value, tolerance = 1e-8)
    expect_length(fit$notes, 0L)
  }
})

test_that("welch-print-groups-F2: subnormal group variances are not taken for resolution noise", {
  ref <- oneway.test(weight ~ feed, data = chickwts)
  lv <- levels(chickwts$feed)
  for (s in c(1e-161, 1e-162, 1e-163, 2^-541)) {
    d <- chickwts
    d$weight <- d$weight * s
    expect_true(any(with(d, tapply(weight, feed, var)) < .Machine$double.xmin))
    fit <- anova_welch(d, "weight", "feed", plots = FALSE)
    expect_equal(fit$anova$statistic, unname(ref$statistic), tolerance = 1e-10)
    expect_equal(fit$anova$den_df, unname(ref$parameter[2L]), tolerance = 1e-10)
    expect_equal(fit$anova$p_value, ref$p.value, tolerance = 1e-8)
    expect_false(any(grepl("floating-point resolution", fit$notes)))
    # Pairwise comparisons: base t.test() on the unscaled data, with the
    # difference and interval in the scaled response's units.
    ph <- fit$posthoc
    expect_false(anyNA(ph$p_value))
    for (r in seq_len(nrow(ph))) {
      x <- chickwts$weight[chickwts$feed == ph$group1[r]]
      z <- chickwts$weight[chickwts$feed == ph$group2[r]]
      tt <- t.test(x, z)
      expect_equal(ph$statistic[r], unname(tt$statistic), tolerance = 1e-10)
      expect_equal(ph$df[r], unname(tt$parameter), tolerance = 1e-10)
      expect_equal(ph$p_value[r], tt$p.value, tolerance = 1e-8)
      expect_equal(ph$difference[r] / s, mean(x) - mean(z), tolerance = 1e-10)
      expect_equal(ph$conf_low[r] / s, tt$conf.int[1L], tolerance = 1e-10)
      expect_equal(ph$conf_high[r] / s, tt$conf.int[2L], tolerance = 1e-10)
    }
    # Standardised differences are unit-free: the same as on the raw data.
    raw <- anova_welch(chickwts, "weight", "feed", plots = FALSE)
    expect_equal(fit$effect_sizes$hedges_g, raw$effect_sizes$hedges_g,
                 tolerance = 1e-10)
    # Group summaries: base mean(), sd() and t.test() on the unscaled data.
    em <- fit$emmeans
    for (g in seq_along(lv)) {
      x <- chickwts$weight[chickwts$feed == lv[g]]
      ci <- t.test(x)$conf.int
      expect_equal(em$mean[g] / s, mean(x), tolerance = 1e-10)
      expect_equal(em$sd[g] / s, sd(x), tolerance = 1e-10)
      expect_equal(em$conf_low[g] / s, ci[1L], tolerance = 1e-10)
      expect_equal(em$conf_high[g] / s, ci[2L], tolerance = 1e-10)
    }
    expect_equal(fit$assumptions$variance_ratio,
                 raw$assumptions$variance_ratio, tolerance = 1e-12)
  }
})
