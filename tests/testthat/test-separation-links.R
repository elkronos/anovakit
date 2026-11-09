# Zero-event cells under links other than logit: the interval of each diverged
# coefficient is open on the side it diverged in, and its finite end is the
# point where the profile deviance crosses the chi-square cutoff.

sep_data <- function() {
  data.frame(g = factor(rep(c("a", "b", "c"), each = 10)),
             y = c(rep(0, 10), rep(0:1, 5), c(1, 1, 1, 1, 0, 1, 1, 1, 1, 0)))
}

# Reference: the end found by solving the profile deviance directly. Group a
# has no events, so with the intercept at b every other cell is fitted freely
# and the profile deviance is that of the a cell alone.
intercept_end <- function(link) {
  inv <- stats::binomial(link)$linkinv
  stats::uniroot(function(b) -2 * 10 * log(1 - inv(b)) - stats::qchisq(0.95, 1),
                 c(-8, 2), tol = 1e-10)$root
}

test_that("a probit fit with a zero-event group is recognised as separated", {
  f <- anova_glm(sep_data(), "y", "g", family = binomial("probit"), plots = FALSE,
                 posthoc = FALSE)
  expect_true(any(grepl("separation detected", f$notes)))
  cf <- f$assumptions$coefficients
  i <- cf$term == "(Intercept)"
  expect_identical(cf$conf_low[i], -Inf)
  expect_equal(cf$conf_high[i], intercept_end("probit"), tolerance = 1e-5)
  expect_identical(cf$conf_high[cf$term != "(Intercept)"], c(Inf, Inf))
  expect_true(all(is.finite(cf$conf_low[cf$term != "(Intercept)"])))
})

test_that("a cloglog fit with a zero-event group gets finite profile ends", {
  f <- anova_glm(sep_data(), "y", "g", family = binomial("cloglog"), plots = FALSE,
                 posthoc = FALSE)
  cf <- f$assumptions$coefficients
  i <- cf$term == "(Intercept)"
  expect_equal(cf$conf_high[i], intercept_end("cloglog"), tolerance = 1e-5)
  expect_true(all(is.finite(cf$conf_low[!i])))
  expect_identical(cf$conf_high[!i], c(Inf, Inf))
  expect_false(any(grepl("search failed", f$notes)))
})

test_that("under a link the direct profile does not suit, the finite end is NA, not open", {
  f <- anova_glm(sep_data(), "y", "g", family = binomial("cauchit"), plots = FALSE,
                 posthoc = FALSE)
  cf <- f$assumptions$coefficients
  expect_false(any(is.infinite(cf$conf_low) & is.infinite(cf$conf_high)))
  expect_true(any(grepl("NA where that search failed", f$notes)))
  expect_false(any(grepl("rule out no value", f$notes)))
})
