# anova_manova(): multivariate analysis of variance and covariance.
#
# The multivariate table is checked against stats::manova for all four
# statistics (a one-way design, where Type I, II and III coincide; the
# multi-factor Type II/III cases are in test-fixes-manova.R). Mardia's tests,
# Box's M and the canonical discriminant analysis are implemented here rather
# than delegated, so each is checked against an independent implementation.

test_that("the multivariate test reproduces stats::manova for every statistic", {
  d <- fx_multivariate()
  for (test in c("Pillai", "Wilks", "Hotelling-Lawley", "Roy")) {
    fit <- anova_manova(d, c("score1", "score2"), "g", test = test,
                        assumptions = FALSE, plots = FALSE)
    ref <- as.data.frame(summary(
      stats::manova(cbind(score1, score2) ~ g, data = d), test = test)$stats)

    expect_equal(fit$multivariate$statistic, ref[["g", test]],
                 info = test)
    expect_equal(fit$multivariate$p_value, ref[["g", "Pr(>F)"]], info = test)
    expect_equal(fit$multivariate$approx_f, ref[["g", "approx F"]], info = test)
  }
})

test_that("estimated marginal means come from emmeans, not from raw means", {
  d <- fx_multivariate()
  fit <- anova_manova(d, c("score1", "score2"), "g", covariates = "age",
                      assumptions = FALSE, plots = FALSE)
  ref <- as.data.frame(suppressMessages(
    emmeans::emmeans(stats::lm(score1 ~ age + g, data = d), ~ g)))
  got <- fit$emmeans[fit$emmeans$response == "score1", ]

  expect_equal(got$estimate, ref$emmean)
  expect_equal(got$se, ref$SE)
  expect_equal(got$conf_low, ref$lower.CL)
  expect_equal(got$conf_high, ref$upper.CL)

  raw <- unname(tapply(d$score1, d$g, mean))
  expect_false(isTRUE(all.equal(got$estimate, raw)))
  # The adjusted standard error is the model's, not a per-group one.
  expect_lt(max(got$se), min(unname(tapply(d$score1, d$g, function(v)
    stats::sd(v) / sqrt(length(v))))))
})

test_that("Box's M matches an independent implementation (rstatix::box_m)", {
  d <- fx_multivariate()
  fit <- anova_manova(d, c("score1", "score2"), "g", plots = FALSE)

  # rstatix::box_m(d[, c("score1", "score2")], d$g) on this fixture. rstatix
  # is not a dependency, so its values are pasted in rather than recomputed.
  expect_equal(fit$assumptions$box_m$statistic, 9.91179139028051, tolerance = 1e-10)
  expect_equal(fit$assumptions$box_m$p_value, 0.128416162751328, tolerance = 1e-10)
  expect_equal(fit$assumptions$box_m$df, 6)
})

test_that("Mardia's tests are computed and have the right degrees of freedom", {
  d <- fx_multivariate()
  fit <- anova_manova(d, c("score1", "score2"), "g", plots = FALSE)
  m <- fit$assumptions$mardia

  expect_equal(nrow(m), 2L)
  expect_setequal(m$test, c("Mardia skewness", "Mardia kurtosis"))
  p <- 2
  expect_equal(m$df[m$test == "Mardia skewness"], p * (p + 1) * (p + 2) / 6)
  expect_true(all(m$p_value >= 0 & m$p_value <= 1))
})

test_that("Mardia's skewness detects a strongly skewed response", {
  d <- withr::with_seed(31, {
    dd <- data.frame(g = factor(rep(c("a", "b"), each = 60)))
    dd$y1 <- stats::rexp(120, rate = 0.5)
    dd$y2 <- stats::rexp(120, rate = 0.5)
    dd
  })
  fit <- anova_manova(d, c("y1", "y2"), "g", plots = FALSE)
  skew_p <- fit$assumptions$mardia$p_value[
    fit$assumptions$mardia$test == "Mardia skewness"]
  expect_lt(skew_p, 0.01)
  expect_true(any(grepl("reject multivariate normality", fit$notes)))
})

test_that("canonical discriminant analysis is computed and is consistent", {
  d <- fx_multivariate()
  fit <- anova_manova(d, c("score1", "score2"), "g", plots = TRUE)

  expect_false(is.null(fit$canonical))
  expect_equal(nrow(fit$canonical), 2L)
  expect_true(all(fit$canonical$canonical_r >= 0 &
                    fit$canonical$canonical_r <= 1))
  expect_true(all(diff(fit$canonical$eigenvalue) <= 0))

  # The eigenvalues are those of solve(E) %*% H from the same model.
  s <- summary(stats::manova(cbind(score1, score2) ~ g, data = d))
  ev <- sort(Re(eigen(solve(s$SS$Residuals) %*% s$SS$g)$values),
             decreasing = TRUE)
  expect_equal(fit$canonical$eigenvalue, ev[seq_len(nrow(fit$canonical))])

  # Wilks' lambda is the product of 1/(1 + eigenvalue), over the fit's own
  # eigenvalues (every nonzero one is kept here).
  wilks <- as.data.frame(summary(
    stats::manova(cbind(score1, score2) ~ g, data = d),
    test = "Wilks")$stats)[["g", "Wilks"]]
  expect_equal(nrow(fit$canonical), 2L)
  expect_equal(prod(1 / (1 + fit$canonical$eigenvalue)), wilks)
  expect_equal(sum(fit$canonical$prop_variance), 1)
  expect_equal(fit$canonical$canonical_r,
               sqrt(fit$canonical$eigenvalue / (1 + fit$canonical$eigenvalue)))

  expect_s3_class(fit$plots$canonical, "ggplot")
})

test_that("effect sizes are per response, floored, and free of residual rows", {
  d <- fx_multivariate()
  fit <- anova_manova(d, c("score1", "score2"), "g", assumptions = FALSE,
                      plots = FALSE)

  expect_setequal(fit$effect_sizes$response, c("score1", "score2"))
  expect_false("Residuals" %in% fit$effect_sizes$term)
  expect_true("partial_omega_sq" %in% names(fit$effect_sizes))
  expect_true(all(fit$effect_sizes$partial_omega_sq >= 0))
  # Partial omega squared, df (MS - MSE) / (df MS + (N - df) MSE)
  a <- stats::anova(stats::lm(score2 ~ g, data = d))
  om <- a$Df[1] * (a$`Mean Sq`[1] - a$`Mean Sq`[2]) /
    (a$Df[1] * a$`Mean Sq`[1] + (nrow(d) - a$Df[1]) * a$`Mean Sq`[2])
  expect_equal(fit$effect_sizes$partial_omega_sq[fit$effect_sizes$response == "score2"],
               om)

  ref <- stats::anova(stats::lm(score1 ~ g, data = d))
  expected <- ref$`Sum Sq`[1] / sum(ref$`Sum Sq`)
  expect_equal(fit$effect_sizes$partial_eta_sq[
    fit$effect_sizes$response == "score1"], expected)
})

test_that("long column names do not break the univariate follow-ups", {
  d <- withr::with_seed(32, {
    dd <- data.frame(score1 = stats::rnorm(400), score2 = stats::rnorm(400))
    for (v in paste0("predictor_variable_number_", 1:5)) {
      dd[[v]] <- factor(sample(c("lo", "hi"), 400, replace = TRUE))
    }
    dd
  })
  gv <- grep("^predictor", names(d), value = TRUE)
  fit <- anova_manova(d, c("score1", "score2"), gv, assumptions = FALSE,
                      plots = FALSE)

  expect_setequal(fit$multivariate$term, gv)
  expect_setequal(attr(stats::terms(fit$univariate$score1$model),
                       "term.labels"), gv)
})

test_that("a single response falls back to a univariate analysis with a note", {
  d <- fx_multivariate()
  fit <- anova_manova(d, "score1", "g", plots = FALSE)

  expect_match(fit$method, "^Univariate")
  expect_true(any(grepl("Only one response", fit$notes)))
  expect_null(fit$multivariate)
  expect_null(fit$canonical)
  ref <- stats::anova(stats::lm(score1 ~ g, data = d))
  expect_equal(fit$anova$p_value[1], ref$`Pr(>F)`[1])
})

test_that("per-response pairwise comparisons are returned", {
  d <- fx_multivariate()
  fit <- anova_manova(d, c("score1", "score2"), "g", assumptions = FALSE,
                      plots = FALSE)
  expect_equal(nrow(fit$posthoc), 2L * choose(3, 2))
  expect_setequal(fit$posthoc$response, c("score1", "score2"))
})

test_that("input validation names the problem", {
  d <- fx_multivariate()
  expect_error(anova_manova(d, c("score1", "g"), "g"),
               "both a response and a grouping variable")
  expect_error(anova_manova(d, c("score1", "nope"), "g"), "not found in `data`")
  expect_error(anova_manova(d, c("score1", "score2"), "g", test = "Bogus"))
})

test_that("nothing is written to stdout", {
  d <- fx_multivariate()
  expect_silent_stdout(anova_manova(d, c("score1", "score2"), "g"))
})
