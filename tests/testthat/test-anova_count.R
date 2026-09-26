# anova_count(): Poisson, quasi-Poisson and negative binomial regression.
#
# Covers the model choice made from the Pearson dispersion, the exposure offset
# that has to sit inside the formula for emmeans to see it, and the marginal
# rates that follow from it. Negative binomial blocks skip when MASS is absent.

test_that("marginal means and post-hoc comparisons exist on the DEFAULT path", {
  d <- fx_counts()
  fit <- anova_count(d, "count", c("g1", "g2"), interaction = TRUE, plots = FALSE)

  expect_equal(nrow(fit$emmeans), 4L)
  expect_true(all(c("estimate", "se", "conf_low", "conf_high") %in%
                    names(fit$emmeans)))
  expect_equal(nrow(fit$posthoc), choose(4, 2))
  expect_true(inherits(fit$emmeans_object, "emmGrid"))
  expect_length(grep("could not be computed", fit$notes), 0L)
})

test_that("marginal rates match an independent emmeans call", {
  d <- fx_counts()
  fit <- anova_count(d, "count", c("g1", "g2"), interaction = TRUE, plots = FALSE)
  ref <- stats::glm(count ~ g1 * g2, data = d, family = stats::poisson())
  emm <- as.data.frame(suppressMessages(
    emmeans::emmeans(ref, ~ g1:g2, type = "response")))
  expect_equal(fit$emmeans$estimate, emm$rate)
  expect_equal(fit$emmeans$se, emm$SE)
})

test_that("vcov_type = 'HC0' changes the standard errors", {
  skip_if_not_installed("sandwich")
  d <- fx_counts()
  plain <- anova_count(d, "count", c("g1", "g2"), interaction = TRUE,
                       plots = FALSE)
  robust <- anova_count(d, "count", c("g1", "g2"), interaction = TRUE,
                        vcov_type = "HC0", plots = FALSE)

  expect_equal(plain$emmeans$estimate, robust$emmeans$estimate)
  expect_false(isTRUE(all.equal(plain$emmeans$se, robust$emmeans$se)))

  ref <- stats::glm(count ~ g1 * g2, data = d, family = stats::poisson())
  V <- sandwich::vcovHC(ref, type = "HC0")
  emm <- as.data.frame(suppressMessages(
    emmeans::emmeans(ref, ~ g1:g2, type = "response", vcov. = V)))
  expect_equal(robust$emmeans$se, emm$SE)
})

test_that("the negative binomial refit happens and matches MASS", {
  skip_if_not_installed("MASS")
  d <- fx_overdispersed()
  fit <- anova_count(d, "y", "g", plots = FALSE)

  expect_identical(fit$model_type, "negbin")
  expect_gt(fit$dispersion, 1.5)

  ref <- MASS::glm.nb(y ~ g, data = d)
  expect_equal(unname(stats::coef(fit$model)), unname(stats::coef(ref)),
               tolerance = 1e-6)
  expect_equal(fit$anova$p_value,
               as.data.frame(car::Anova(ref, test.statistic = "LR"))$`Pr(>Chisq)`,
               tolerance = 1e-6)
})

test_that("a negative binomial fit works with an exposure offset too", {
  skip_if_not_installed("MASS")
  d <- withr::with_seed(11, {
    dd <- data.frame(g = factor(rep(c("A", "B", "C"), each = 100)))
    dd$expo <- stats::runif(300, 0.5, 3)
    dd$y <- stats::rnbinom(300, mu = rep(c(5, 6, 7), each = 100) * dd$expo,
                           size = 0.7)
    dd
  })
  fit <- anova_count(d, "y", "g", offset = "expo", plots = FALSE)
  ref <- MASS::glm.nb(y ~ g + offset(log(expo)), data = d)

  expect_identical(fit$model_type, "negbin")
  expect_equal(unname(stats::coef(fit$model)), unname(stats::coef(ref)),
               tolerance = 1e-6)
  expect_equal(nrow(fit$emmeans), 3L)
})

test_that("the offset lives in the formula, so emmeans and terms can see it", {
  d <- fx_counts()
  fit <- anova_count(d, "count", "g1", offset = "hours",
                     model = "poisson", plots = FALSE)

  expect_false(is.null(attr(stats::terms(fit$model), "offset")))
  # One row per cell, not one row per observation.
  expect_equal(nrow(fit$emmeans), 2L)
  expect_true(any(grepl("rates at hours = 1", fit$notes)))

  ref <- stats::glm(count ~ g1 + offset(log(hours)), data = d,
                    family = stats::poisson())
  expect_equal(unname(stats::coef(fit$model)), unname(stats::coef(ref)))

  # The rates really are at one unit of exposure, as the note says, and not
  # at the mean log exposure emmeans would otherwise use
  emm <- as.data.frame(emmeans::emmeans(ref, ~ g1, type = "response", offset = 0))
  expect_equal(fit$emmeans$estimate, emm$rate)
  expect_equal(fit$emmeans$estimate, unname(exp(cumsum(stats::coef(ref)))))
  expect_equal(fit$emmeans$estimate, c(3.431959, 7.874468), tolerance = 1e-6)
})

test_that("the model choice can be forced and is always reported", {
  skip_if_not_installed("MASS")
  d <- fx_overdispersed()
  expect_identical(anova_count(d, "y", "g", model = "poisson",
                               plots = FALSE)$model_type, "poisson")
  expect_identical(anova_count(d, "y", "g", model = "quasipoisson",
                               plots = FALSE)$model_type, "quasipoisson")
  expect_identical(anova_count(d, "y", "g", model = "negbin",
                               plots = FALSE)$model_type, "negbin")

  auto <- anova_count(d, "y", "g", plots = FALSE)
  expect_true(any(grepl("Pearson dispersion is", auto$notes)))
})

test_that("overdispersion_threshold decides the automatic model choice", {
  skip_if_not_installed("MASS")
  d <- fx_overdispersed()
  phi <- anova_count(d, "y", "g", model = "poisson", plots = FALSE)$dispersion
  above <- anova_count(d, "y", "g", overdispersion_threshold = phi + 0.01,
                       plots = FALSE)
  below <- anova_count(d, "y", "g", overdispersion_threshold = phi - 0.01,
                       plots = FALSE)
  expect_identical(above$model_type, "poisson")
  expect_true(any(grepl(sprintf("at or below the threshold of %.2f", phi + 0.01),
                        above$notes)))
  expect_identical(below$model_type, "negbin")
  expect_true(any(grepl(sprintf("above the threshold of %.2f", phi - 0.01),
                        below$notes)))
})

test_that("a quasi-Poisson model is tested with F, as car::Anova does", {
  d <- fx_overdispersed()
  fit <- anova_count(d, "y", "g", model = "quasipoisson", plots = FALSE)
  ref <- car::Anova(stats::glm(y ~ g, family = stats::quasipoisson(), data = d),
                    type = 2, test.statistic = "F")
  expect_identical(attr(fit$anova, "statistic"), "F")
  expect_equal(fit$anova$statistic[1], ref$`F value`[1])
  expect_equal(fit$anova$p_value[1], ref$`Pr(>F)`[1])
  expect_true(any(grepl("has no likelihood, so the analysis of deviance uses an F test",
                        fit$notes)))
})

test_that("dispersion is the Pearson statistic", {
  d <- fx_counts()
  fit <- anova_count(d, "count", c("g1", "g2"), interaction = TRUE,
                     model = "poisson", plots = FALSE)
  ref <- stats::glm(count ~ g1 * g2, data = d, family = stats::poisson())
  expected <- sum(stats::residuals(ref, type = "pearson")^2) /
    stats::df.residual(ref)
  expect_equal(fit$dispersion, expected)
})

test_that("type and test_statistic are separate arguments", {
  d <- fx_counts()
  t2 <- anova_count(d, "count", c("g1", "g2"), interaction = TRUE, type = "II",
                    plots = FALSE)
  t3 <- anova_count(d, "count", c("g1", "g2"), interaction = TRUE, type = "III",
                    plots = FALSE)

  expect_true(all(unlist(t3$model$contrasts) == "contr.sum"))
  ref <- as.data.frame(car::Anova(
    withr::with_options(list(contrasts = c("contr.sum", "contr.poly")),
                        stats::glm(count ~ g1 * g2, data = d,
                                   family = stats::poisson())),
    type = 3, test.statistic = "LR"))
  expect_equal(t3$anova$statistic, ref$`LR Chisq`)
  expect_false(isTRUE(all.equal(t2$anova$statistic[1], t3$anova$statistic[1])))

  # test_statistic varies independently of type
  wald <- anova_count(d, "count", c("g1", "g2"), interaction = TRUE, type = "II",
                      model = "poisson", test_statistic = "Wald", plots = FALSE)
  wref <- as.data.frame(car::Anova(
    stats::glm(count ~ g1 * g2, data = d, family = stats::poisson()),
    type = 2, test.statistic = "Wald"))
  expect_identical(attr(wald$anova, "statistic"), "Wald chi-square")
  expect_equal(wald$anova$statistic, wref$Chisq)
  expect_equal(wald$anova$p_value, wref$`Pr(>Chisq)`)
})

test_that("a theta that collapses towards Poisson is reported as such", {
  skip_if_not_installed("MASS")
  # Converges (no th.warn) to a theta above 1000: essentially Poisson data
  d <- data.frame(g = factor(rep(c("a", "b"), each = 20)),
                  y = c(3, 3, 6, 9, 4, 2, 4, 7, 4, 2, 4, 7, 5, 3, 4, 5, 1, 3, 5, 9,
                        8, 5, 6, 6, 8, 5, 7, 10, 12, 6, 11, 9, 7, 12, 2, 8, 3, 5, 10, 6))
  ref <- MASS::glm.nb(y ~ g, data = d)
  expect_null(ref$th.warn)
  expect_gt(ref$theta, 1000)
  fit <- anova_count(d, "y", "g", model = "negbin", plots = FALSE)
  expect_true(any(grepl("very large \\(theta = 1\\.18e\\+03\\)", fit$notes)))
  expect_false(any(grepl("did not converge|theta = [0-9.]+ \\(SE", fit$notes)))
})

test_that("the model is additive by default even with three factors", {
  d <- withr::with_seed(4, {
    dd <- data.frame(a = factor(sample(letters[1:3], 300, TRUE)),
                     b = factor(sample(c("x", "y"), 300, TRUE)),
                     c3 = factor(sample(c("p", "q"), 300, TRUE)))
    dd$cnt <- stats::rpois(300, 6)
    dd
  })
  add <- anova_count(d, "cnt", c("a", "b", "c3"), plots = FALSE)
  full <- anova_count(d, "cnt", c("a", "b", "c3"), interaction = TRUE,
                      plots = FALSE)

  expect_setequal(attr(stats::terms(add$model), "term.labels"),
                  c("a", "b", "c3"))
  expect_length(attr(stats::terms(full$model), "term.labels"), 7L)
})

test_that("counts are validated", {
  d <- fx_counts()
  bad <- d; bad$count[1] <- -1
  expect_error(anova_count(bad, "count", "g1"), "cannot be negative")
  frac <- d; frac$count <- frac$count + 0.5
  expect_error(anova_count(frac, "count", "g1"), "whole numbers")
  chr <- d; chr$count <- as.character(chr$count)
  expect_error(anova_count(chr, "count", "g1"), "must be numeric")
  zero <- d; zero$hours[1] <- 0
  expect_error(anova_count(zero, "count", "g1", offset = "hours"),
               "strictly positive")
})

test_that("sparse and empty cells are reported", {
  d <- fx_nested()
  d$count <- withr::with_seed(112, stats::rpois(nrow(d), 4))
  fit <- anova_count(d, "count", c("g1", "g2"), plots = FALSE)
  expect_true(any(grepl("contain no observations", fit$notes)))
})

test_that("nothing is written to stdout", {
  d <- fx_counts()
  expect_silent_stdout(anova_count(d, "count", c("g1", "g2")))
  expect_silent_stdout(anova_count(d, "count", "g1", verbose = FALSE))
})
