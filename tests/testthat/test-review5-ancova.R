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
  expect_true(any(grepl(sprintf("fewer than 2.5 residual degrees of freedom (a with %s;",
                                signif(rdf[["a"]], 3)), lev, fixed = TRUE)))
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

# ancova-C1: cells of three are left out also under additive factors ---------

test_that("ancova-C1: cells of three observations are left out of Levene's test whatever their residual df", {
  # Additive A + B + x with three observations per cell: each cell keeps more
  # than 2.5 residual df (reference: lm()'s leverages), but the test on
  # groups of three cannot reject. Reference for that: car::leveneTest on
  # raw groups of three, with the SD tripled in 40% of them, essentially
  # never rejects.
  d <- withr::with_seed(5, {
    d <- expand.grid(A = letters[1:5], B = LETTERS[1:5], r = 1:3)
    d$x <- stats::rnorm(75)
    d$y <- d$x + stats::rnorm(75)
    d
  })
  m <- stats::lm(y ~ A + B + x, data = d)
  rdf <- tapply(1 - stats::hatvalues(m), interaction(d$A, d$B), sum)
  expect_true(all(rdf > 2.5))
  p_car <- withr::with_seed(6, replicate(200, {
    g <- factor(rep(1:25, each = 3))
    y <- stats::rnorm(75, sd = ifelse(as.integer(g) <= 10, 3, 1))
    car::leveneTest(y ~ g)[["Pr(>F)"]][1]
  }))
  expect_lt(mean(p_car < 0.05), 0.02)
  f <- anova_ancova(d, "y", c("A", "B"), "x", interaction = FALSE,
                    force_interaction = FALSE, plots = FALSE, posthoc = FALSE)
  expect_null(f$assumptions$levene)
  expect_true(any(grepl("Levene's test was skipped: it needs at least two groups with four or more observations, and no group has that many",
                        f$notes, fixed = TRUE)))
  expect_true(any(grepl("Levene's test leaves out 25 group(s) of exactly three observations",
                        f$notes, fixed = TRUE)))

  # Mixed sizes: the cells of three are dropped, the rest is car's test
  d2 <- withr::with_seed(7, {
    d <- expand.grid(A = letters[1:5], B = LETTERS[1:5], r = 1:4)
    d <- d[!(d$A == "a" & d$B == "A" & d$r > 3), ]
    d$x <- stats::rnorm(nrow(d))
    d$y <- d$x + stats::rnorm(nrow(d))
    d
  })
  f2 <- anova_ancova(d2, "y", c("A", "B"), "x", interaction = FALSE,
                     force_interaction = FALSE, plots = FALSE, posthoc = FALSE)
  m2 <- stats::lm(y ~ A + B + x, data = d2)
  cell <- interaction(d2$A, d2$B, sep = " : ")
  expect_gt(tapply(1 - stats::hatvalues(m2), cell, sum)[["a : A"]], 2.5)
  k <- cell != "a : A"
  ref <- car::leveneTest(stats::residuals(m2)[k] ~ droplevels(cell[k]))
  expect_equal(f2$assumptions$levene$df1, 23)
  expect_equal(f2$assumptions$levene$statistic, ref[["F value"]][1])
  expect_equal(f2$assumptions$levene$p_value, ref[["Pr(>F)"]][1])
  expect_true(any(grepl("leaves out 1 group(s) of exactly three observations (a : A)",
                        f2$notes, fixed = TRUE)))
})

# ancova-R1: the printed residual df stays below the threshold ---------------

test_that("ancova-R1: the Levene note prints a residual df below 2.5 as below 2.5 and does not misdescribe the cell", {
  d <- withr::with_seed(1937, {
    d <- data.frame(g = rep(c("a", "b", "c"), c(4, 12, 12)),
                    x = stats::rnorm(28))
    d$y <- d$x + stats::rnorm(28)
    d
  })
  rdf <- tapply(1 - stats::hatvalues(stats::lm(y ~ g + x, data = d)), d$g, sum)
  expect_lt(rdf[["a"]], 2.5)
  expect_gt(rdf[["a"]], 2.45)
  f <- anova_ancova(d, "y", "g", "x", force_interaction = FALSE,
                    plots = FALSE, posthoc = FALSE)
  lev <- grep("residual degrees of freedom \\(a with", f$notes, value = TRUE)
  expect_length(lev, 1L)
  shown <- as.numeric(sub(".*\\(a with ([0-9.]+);.*", "\\1", lev))
  expect_lt(shown, 2.5)
  expect_equal(shown, signif(rdf[["a"]], 3))
  expect_false(grepl("group of three", lev, fixed = TRUE))
  expect_false(grepl("at most two observations more", lev, fixed = TRUE))
  # Rounding never lifts a value below 2.5 to 2.5
  expect_identical(.rdf_below(2.4997), "2.4997")
  expect_identical(.rdf_below(2), "2")
})
