# Regression tests for the round-5 review items in the ancova area.

# ancova-levene-plot-F1: a cell with fractional residual df is not "exact" ----

test_that("ancova-levene-plot-F1: a cell with residual df 0.43 is not called exactly fitted", {
  d <- withr::with_seed(1, {
    d <- data.frame(g = rep(c("a", "b", "c"), c(3, 10, 10)))
    d$x1 <- c(-1, 0, 1, stats::rnorm(20, sd = 0.18))
    d$x2 <- c(0.5, -1, 0.5, stats::rnorm(20, sd = 0.18))
    d$y <- stats::rnorm(23)
    d
  })
  f <- anova_ancova(d, "y", "g", c("x1", "x2"), force_interaction = FALSE,
                    plots = FALSE)
  # Reference: lm()'s leverages and residuals for the additive model
  m <- stats::lm(y ~ g + x1 + x2, data = d)
  rdf <- tapply(1 - stats::hatvalues(m), d$g, sum)
  expect_equal(unname(round(rdf[["a"]], 2)), 0.43)
  expect_true(all(rdf[c("b", "c")] > 2.5))
  # Cell a's residuals are not zero, so it must not be called exactly fitted
  expect_gt(min(abs(stats::residuals(m)[d$g == "a"])), 0.01)
  lev <- grep("^Levene", f$notes, value = TRUE)
  expect_false(any(grepl("fits exactly", lev, fixed = TRUE)))
  expect_false(any(grepl("zero whatever the response", lev, fixed = TRUE)))
  expect_true(any(grepl("fewer than 2.5 residual degrees of freedom (a with 0.43;",
                        lev, fixed = TRUE)))
  # The test itself is car's on the two cells that are kept
  r <- stats::residuals(m)
  k <- d$g != "a"
  ref <- car::leveneTest(r[k] ~ factor(d$g[k]))
  expect_equal(f$assumptions$levene$statistic, ref[["F value"]][1])
  expect_equal(f$assumptions$levene$p_value, ref[["Pr(>F)"]][1])
})

# ancova-levene-plot-F2: cells with two residual df are left out -------------

test_that("ancova-levene-plot-F2: cells left with fewer than 2.5 residual df are left out of Levene's test", {
  # Three groups of four with a slope each: two residual df per cell, which
  # made the test reject homoscedastic data about 11% of the time at 5%.
  # Reference: the cells' residual df from lm()'s leverages.
  d <- withr::with_seed(11, {
    d <- data.frame(g = rep(letters[1:3], each = 4), x = stats::rnorm(12))
    d$y <- d$x + stats::rnorm(12)
    d
  })
  f <- anova_ancova(d, "y", "g", "x", force_interaction = TRUE,
                    plots = FALSE, posthoc = FALSE)
  m <- stats::lm(y ~ g * x, data = d)
  expect_equal(as.vector(tapply(1 - stats::hatvalues(m), d$g, sum)), c(2, 2, 2))
  expect_null(f$assumptions$levene)
  expect_true(any(grepl("Levene's test was skipped: it needs at least two groups with three or more observations that the model leaves with at least 2.5 residual degrees of freedom, and no group has that many",
                        f$notes, fixed = TRUE)))
  expect_true(any(grepl("(a with 2, b with 2, c with 2;", f$notes, fixed = TRUE)))
  expect_false(any(grepl("Levene's test rejects", f$notes, fixed = TRUE)))

  # A cell of four beside two large ones: the test runs on the large ones
  d2 <- withr::with_seed(3, {
    d <- data.frame(g = rep(c("a", "b", "c"), c(4, 15, 15)),
                    x = stats::rnorm(34))
    d$y <- d$x + stats::rnorm(34)
    d
  })
  f2 <- anova_ancova(d2, "y", "g", "x", force_interaction = TRUE,
                     plots = FALSE, posthoc = FALSE)
  r2 <- stats::residuals(stats::lm(y ~ g * x, data = d2))
  k <- d2$g != "a"
  ref <- car::leveneTest(r2[k] ~ factor(d2$g[k]))
  expect_equal(f2$assumptions$levene$df1, 1)
  expect_equal(f2$assumptions$levene$statistic, ref[["F value"]][1])
  expect_equal(f2$assumptions$levene$p_value, ref[["Pr(>F)"]][1])
  expect_true(any(grepl("(a with 2;", f2$notes, fixed = TRUE)))

  # Under the additive model a cell of three has at most two residual df
  # (each cell has an intercept of its own), so it is left out too
  d3 <- withr::with_seed(4, {
    d <- data.frame(g = rep(c("a", "b", "c"), c(3, 15, 15)),
                    x = stats::rnorm(33))
    d$y <- d$x + stats::rnorm(33)
    d
  })
  f3 <- anova_ancova(d3, "y", "g", "x", force_interaction = FALSE,
                     plots = FALSE, posthoc = FALSE)
  m3 <- stats::lm(y ~ g + x, data = d3)
  rdf3 <- tapply(1 - stats::hatvalues(m3), d3$g, sum)
  expect_lt(rdf3[["a"]], 2)
  r3 <- stats::residuals(m3)
  k3 <- d3$g != "a"
  ref3 <- car::leveneTest(r3[k3] ~ factor(d3$g[k3]))
  expect_equal(f3$assumptions$levene$statistic, ref3[["F value"]][1])
  expect_equal(f3$assumptions$levene$p_value, ref3[["Pr(>F)"]][1])
})
