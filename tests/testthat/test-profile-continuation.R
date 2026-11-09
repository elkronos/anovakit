# Detailed regression tests: run locally and on CI (NOT_CRAN=true), skipped
# on CRAN to keep its check time short. The core suite runs everywhere.
skip_on_cran()

# A zero-count group drives its rate ratio to Inf; the finite end is where the
# profile deviance crosses the cutoff, which can lie more than 32 log units
# from the diverged estimate. References: profiles over the full parameter
# vector (theta fixed at the fit's value for the negative binomial model).
test_that("a zero-count rate ratio's finite end is found however far it lies from the estimate", {
  skip_if_not_installed("MASS")
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
