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

test_that("a crossing past a point no refit can reach is found by stepping back", {
  skip_if_not_installed("MASS")
  # The step from 16 to 32 units leaps past the crossing to where no refit
  # converges; the search steps back and solves between verified points.
  # Reference: full-parameter profile at the fit's theta.
  d <- withr::with_seed(3, {
    d <- expand.grid(r = 1:10, g = factor(letters[1:3]), s = factor(c("A", "B")))
    d$e <- stats::runif(nrow(d), 0.5, 2)
    d$y <- stats::rnbinom(nrow(d), mu = 30 * d$e, size = 2)
    d$y[d$g == "b" & d$s == "B"] <- 0
    d
  })
  r <- suppressWarnings(anova_count(d[, c("y", "g", "s")], "y", c("g", "s"),
                                    interaction = TRUE, model = "negbin",
                                    plots = FALSE, posthoc = FALSE))$effect_sizes
  i <- r$term == "gb:sB"
  expect_identical(r$conf_low[i], 0)
  expect_equal(r$conf_high[i], 0.00705, tolerance = 1e-3)
})

test_that("an end no refit could verify is not reported as open", {
  # gb:sB is bounded above (the two zero cells can only push it down); its
  # profile crosses the cut-off between 32 and 64 units from the estimate.
  d <- data.frame(g = factor(rep(letters[1:4], each = 4)),
                  s = factor(rep(c("A", "B"), each = 16)),
                  y = c(0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0,
                        0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1))
  r <- suppressWarnings(anova_glm(d, "y", c("g", "s"), family = binomial("cloglog"),
                                  interaction = TRUE, plots = FALSE, posthoc = FALSE))
  cf <- r$assumptions$coefficients
  expect_equal(cf$conf_high[cf$term == "gb:sB"], -2.3679, tolerance = 1e-4)
  expect_false(any(grepl("gb:sB[^.]*rule out no value", r$notes)))
})
