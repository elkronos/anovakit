# anova_bin(): logistic regression for a binary response.
#
# Tables are checked against car::Anova and intervals against MASS's profile
# likelihood. The rest covers what makes a logistic fit go wrong quietly:
# separation, a response with one observed level, and the choice of which level
# counts as a success.

test_that("the deviance table reproduces car::Anova on the same model", {
  d <- fx_binary()
  fit <- anova_bin(d, "y", "g", plots = FALSE)
  ref <- as.data.frame(car::Anova(
    stats::glm(y ~ g, data = d, family = stats::binomial()),
    type = 2, test.statistic = "LR"))

  expect_equal(fit$anova$statistic, ref$`LR Chisq`)
  expect_equal(fit$anova$p_value, ref$`Pr(>Chisq)`)
})

test_that("odds ratios are the exponentiated coefficients", {
  d <- fx_binary()
  fit <- anova_bin(d, "y", "g", plots = FALSE)
  ref <- stats::glm(y ~ g, data = d, family = stats::binomial())
  expected <- exp(stats::coef(ref))[-1]

  expect_equal(fit$effect_sizes$odds_ratio, unname(expected))
  expect_false("(Intercept)" %in% fit$effect_sizes$term)
})

test_that("profile intervals are the default and match stats::confint", {
  d <- fx_binary()
  fit <- anova_bin(d, "y", "g", plots = FALSE)
  ref <- stats::glm(y ~ g, data = d, family = stats::binomial())
  ci <- suppressMessages(stats::confint(ref))[-1, , drop = FALSE]

  expect_identical(unique(fit$effect_sizes$ci_method), "profile")
  expect_equal(fit$effect_sizes$conf_low, unname(exp(ci[, 1])))
  expect_equal(fit$effect_sizes$conf_high, unname(exp(ci[, 2])))
})

test_that("wald intervals can be requested and differ from profile ones", {
  d <- fx_binary()
  prof <- anova_bin(d, "y", "g", ci_method = "profile", plots = FALSE)
  wald <- anova_bin(d, "y", "g", ci_method = "wald", plots = FALSE)
  expect_identical(unique(wald$effect_sizes$ci_method), "wald")
  expect_false(isTRUE(all.equal(prof$effect_sizes$conf_low,
                                wald$effect_sizes$conf_low)))
})

test_that("a response with only one observed level is rejected", {
  d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 20)),
                  y = rep(0, 60))
  expect_error(anova_bin(d, "y", "g"), "exactly two distinct values")

  d2 <- data.frame(g = factor(rep(c("a", "b"), each = 20)),
                   y = factor(rep("yes", 40), levels = c("no", "yes")))
  expect_error(anova_bin(d2, "y", "g"), "exactly two distinct values")
})

test_that("a non-binary numeric response is rejected", {
  d <- data.frame(g = factor(rep(c("a", "b"), each = 20)), y = rep(0:3, 10))
  expect_error(anova_bin(d, "y", "g"), "other than 0 and 1")
})

test_that("separation is detected and named", {
  d <- fx_separated()
  fit <- anova_bin(d, "y", "g", plots = FALSE)
  sep <- grep("separation", fit$notes, value = TRUE)

  expect_length(sep, 1L)
  expect_match(sep, "The affected coefficient\\(s\\): gc")
  # The offending odds ratio really is unusable, which is why it is flagged.
  expect_gt(fit$effect_sizes$odds_ratio[fit$effect_sizes$term == "gc"], 1e5)
})

test_that("a well-behaved fit is not flagged for separation", {
  d <- fx_binary()
  fit <- anova_bin(d, "y", "g", plots = FALSE)
  expect_length(grep("separation", fit$notes), 0L)
})

test_that("interactions between grouping variables can be fitted", {
  d <- withr::with_seed(2, {
    dd <- data.frame(g1 = factor(sample(c("a", "b"), 400, TRUE)),
                     g2 = factor(sample(c("x", "y"), 400, TRUE)))
    lp <- -1 + (dd$g1 == "b") + (dd$g2 == "y") +
      3 * (dd$g1 == "b") * (dd$g2 == "y")
    dd$y <- stats::rbinom(400, 1, stats::plogis(lp))
    dd
  })
  add <- anova_bin(d, "y", c("g1", "g2"), plots = FALSE)
  full <- anova_bin(d, "y", c("g1", "g2"), interaction = TRUE, plots = FALSE)

  expect_setequal(add$anova$term, c("g1", "g2"))
  expect_setequal(full$anova$term, c("g1", "g2", "g1:g2"))
  expect_lt(full$anova$p_value[full$anova$term == "g1:g2"], 0.001)
})

test_that("marginal probabilities and pairwise odds ratios are returned", {
  d <- fx_binary()
  fit <- anova_bin(d, "y", "g", plots = FALSE)

  expect_equal(nrow(fit$emmeans), 3L)
  expect_true(all(c("estimate", "conf_low", "conf_high") %in% names(fit$emmeans)))
  expect_true(all(fit$emmeans$estimate > 0 & fit$emmeans$estimate < 1))
  expect_equal(nrow(fit$posthoc), 3L)
  expect_true("ratio" %in% names(fit$posthoc))
  expect_true(inherits(fit$emmeans_object, "emmGrid"))
})

test_that("marginal probabilities match an independent emmeans call", {
  d <- fx_binary()
  fit <- anova_bin(d, "y", "g", plots = FALSE)
  ref <- stats::glm(y ~ g, data = d, family = stats::binomial())
  emm <- as.data.frame(suppressMessages(
    emmeans::emmeans(ref, ~ g, type = "response")))
  expect_equal(fit$emmeans$estimate, emm$prob)
})

test_that("success and reference levels can be chosen", {
  d <- fx_binary()
  default <- anova_bin(d, "y", "g", plots = FALSE)
  flipped <- anova_bin(d, "y", "g", success = "0", plots = FALSE)
  expect_equal(flipped$effect_sizes$odds_ratio,
               1 / default$effect_sizes$odds_ratio)

  relevelled <- anova_bin(d, "y", "g", reference = list(g = "c"), plots = FALSE)
  expect_setequal(relevelled$effect_sizes$term, c("ga", "gb"))
  expect_error(anova_bin(d, "y", "g", reference = list(g = "zzz")),
               "not found in")
  expect_error(anova_bin(d, "y", "g", reference = list(nope = "a")),
               "must be grouping variables")
})

test_that("robust standard errors change the answer when asked for", {
  skip_if_not_installed("sandwich")
  d <- fx_binary()
  plain <- anova_bin(d, "y", "g", ci_method = "wald", plots = FALSE)
  robust <- anova_bin(d, "y", "g", ci_method = "wald", vcov_type = "HC0",
                      plots = FALSE)
  expect_false(isTRUE(all.equal(plain$effect_sizes$se_log_or,
                                robust$effect_sizes$se_log_or)))
})

test_that("conf_level reaches the odds-ratio intervals", {
  d <- fx_binary()
  wide <- anova_bin(d, "y", "g", conf_level = 0.99, plots = FALSE)
  narrow <- anova_bin(d, "y", "g", conf_level = 0.80, plots = FALSE)
  expect_true(all((wide$effect_sizes$conf_high - wide$effect_sizes$conf_low) >
                    (narrow$effect_sizes$conf_high - narrow$effect_sizes$conf_low)))
})

test_that("reserved arguments cannot be smuggled through ...", {
  d <- fx_binary()
  # `family` is not an argument of anova_bin() at all, so R rejects it before
  # anything can act on it.
  expect_error(anova_bin(d, "y", "g", family = stats::binomial()),
               "unused argument")
})

test_that("nothing is written to stdout, including the plot", {
  d <- fx_binary()
  expect_silent_stdout(anova_bin(d, "y", "g"))
  fit <- anova_bin(d, "y", "g")
  expect_true(all(vapply(fit$plots, inherits, logical(1), "ggplot")))
})
