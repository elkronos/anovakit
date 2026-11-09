# Detailed regression tests: run locally and on CI (NOT_CRAN=true), skipped
# on CRAN to keep its check time short. The core suite runs everywhere.
skip_on_cran()

# The finite end of a diverged coefficient's interval is found by stepping
# away from the estimate and refitting with the coefficient fixed.

test_that("a zero-count rate ratio's finite end is found more than 32 units from the estimate", {
  skip_if_not_installed("MASS")
  # The estimate diverges to about 37 log units; the profile crosses the
  # cut-off about 34 units from it. References: profiles over the full
  # parameter vector (theta fixed at the fit's value for the negative
  # binomial model).
  dc <- withr::with_seed(3, {
    dc <- data.frame(g = factor(rep(c("a", "b", "c"), each = 8)),
                     h = factor(rep(c("u", "v"), 12)), e = stats::runif(24, 1, 5))
    dc$y <- c(rep(0, 8), stats::rnbinom(16, mu = 5 * dc$e[9:24], size = 1.5))
    dc
  })
  nb <- anova_count(dc, "y", c("g", "h"), offset = "e", model = "negbin",
                    plots = FALSE, posthoc = FALSE)$effect_sizes
  expect_equal(nb$conf_low[nb$term %in% c("gb", "gc")], c(35.76, 54.55), tolerance = 1e-3)
  expect_identical(nb$conf_high[nb$term %in% c("gb", "gc")], c(Inf, Inf))
  po <- anova_count(dc, "y", c("g", "h"), offset = "e", model = "poisson",
                    plots = FALSE, posthoc = FALSE)$effect_sizes
  expect_equal(po$conf_low[po$term %in% c("gb", "gc")], c(43.61, 66.45), tolerance = 1e-3)
})

test_that("a refit that stops far above the minimum does not end an unidentified interval", {
  # On some platforms glm.fit(), started from the data at a very large
  # offset, reports convergence hundreds of deviance units above the
  # minimum. That is simulated here, so the test does not depend on the
  # platform: gb:ss is not identified, and its interval stays (0, Inf).
  d <- withr::with_seed(4, {
    d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 20)),
                    s = factor(rep(c("n", "s"), 30)))
    d$c <- stats::rpois(60, 3)
    d$c[d$g == "b"] <- 0
    d
  })
  real_fit <- stats::glm.fit
  local_mocked_bindings(
    glm.fit = function(x, y, ..., offset = NULL) {
      f <- real_fit(x, y, ..., offset = offset)
      if (!is.null(offset) && max(abs(offset)) > 100) {
        f$deviance <- f$deviance + 500
        f$converged <- TRUE
      }
      f
    },
    .package = "stats")
  fc <- anova_count(d, "c", c("g", "s"), interaction = TRUE, plots = FALSE,
                    model = "poisson")
  rc <- fc$effect_sizes[fc$effect_sizes$term == "gb:ss", ]
  expect_identical(c(rc$conf_low, rc$conf_high), c(0, Inf))
  expect_true(any(grepl("gb:ss is open on both sides", fc$notes)))
})
