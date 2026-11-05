# Regression tests for the anova_manova() findings of the second review
# (R2-MANOVA-*). Each expectation is checked against an independent reference:
# stats (lm, anova.mlm, summary.manova), car or emmeans, or a hand count.

test_that("R2-MANOVA-01: a constant response is refused, not whitened into a spurious response", {
  d <- withr::with_seed(4, {
    d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 15)))
    d$y1 <- stats::rnorm(45) + (d$g == "b")
    d
  })
  d$y2 <- 5
  # The scenario: y2 has no variation at all, yet lm() leaves rounding residuals
  # that a purely relative check (residual SS <= 1e-12 * total SS = 0) misses.
  expect_equal(stats::var(d$y2), 0)
  expect_gt(sum(stats::residuals(stats::lm(y2 ~ g, d))^2), 0)
  expect_error(anova_manova(d, c("y1", "y2"), "g"),
               "fits `y2` exactly", fixed = TRUE)
  d$y2 <- 0.1
  expect_error(anova_manova(d, c("y1", "y2"), "g"),
               "fits `y2` exactly", fixed = TRUE)

  # A response with a large level and a small (but real) spread still gets its
  # test, equal to summary.manova's.
  d$y2 <- 1e3 + withr::with_seed(5, stats::rnorm(45))
  f <- anova_manova(d, c("y1", "y2"), "g", plots = FALSE)
  ref <- summary(stats::manova(cbind(y1, y2) ~ g, data = d), test = "Pillai")$stats
  expect_equal(f$multivariate$statistic[f$multivariate$term == "g"],
               ref["g", "Pillai"], tolerance = 1e-8)
  expect_equal(f$multivariate$p_value[f$multivariate$term == "g"],
               ref["g", "Pr(>F)"], tolerance = 1e-8)
})

test_that("R2-MANOVA-02: a Roy slopes test with min(df, p) > 1 is flagged as a lower-bound p-value", {
  d <- withr::with_seed(3, {
    n <- 60
    d <- data.frame(g = sample(c("a", "b", "c"), n, TRUE), x = stats::rnorm(n))
    d$y1 <- stats::rnorm(n) + (d$g == "b") + d$x * (1 + (d$g == "c"))
    d$y2 <- stats::rnorm(n) + d$x
    d$y3 <- stats::rnorm(n)
    d
  })
  f <- anova_manova(d, c("y1", "y2", "y3"), "g", covariates = "x",
                    test = "Roy", plots = FALSE)
  ref <- stats::anova(stats::lm(cbind(y1, y2, y3) ~ g + x, d),
                      stats::lm(cbind(y1, y2, y3) ~ g * x, d), test = "Roy")
  expect_equal(f$slopes_test$statistic, ref$Roy[2], tolerance = 1e-8)
  expect_equal(f$slopes_test$p_value, ref[["Pr(>F)"]][2], tolerance = 1e-8)
  # min(2 df, 3 responses) = 2 > 1: Roy's F is only an upper bound.
  expect_gt(min(abs(ref$Df[2]), 3), 1)
  expect_true(any(grepl("homogeneity-of-slopes test, Roy's largest root has only an upper bound",
                        f$notes, fixed = TRUE)))
  rej <- grep("Homogeneity of regression slopes is rejected", f$notes, value = TRUE)
  expect_length(rej, 1L)
  expect_match(rej, "a lower bound", fixed = TRUE)

  # Pillai's trace is calibrated, and gets no such caveat.
  fp <- anova_manova(d, c("y1", "y2", "y3"), "g", covariates = "x",
                     plots = FALSE)
  expect_false(any(grepl("upper bound", fp$notes, fixed = TRUE)))
})

test_that("R2-MANOVA-03: $posthoc labels the response in the same column as $emmeans", {
  d <- withr::with_seed(20, {
    d <- data.frame(response = sample(c("u", "v", "w"), 60, TRUE))
    d$y1 <- stats::rnorm(60)
    d$y2 <- stats::rnorm(60)
    d
  })
  f <- anova_manova(d, c("y1", "y2"), "response", plots = FALSE)
  expect_true(all(c("response.1", "response") %in% names(f$emmeans)))
  expect_setequal(unique(f$emmeans$response.1), c("y1", "y2"))
  expect_true("response.1" %in% names(f$posthoc))
  expect_false("response" %in% names(f$posthoc))
  expect_setequal(unique(f$posthoc$response.1), c("y1", "y2"))
  # The y1 comparisons are emmeans' own on the one-way model for y1.
  ref <- as.data.frame(emmeans::contrast(
    emmeans::emmeans(stats::lm(y1 ~ response, d), "response"),
    "pairwise", adjust = "tukey"))
  ph <- f$posthoc[f$posthoc$response.1 == "y1", ]
  expect_equal(ph$estimate, ref$estimate, tolerance = 1e-8)
})

test_that("R2-MANOVA-04: the skipped-comparisons note points to the per-response grids and the count spans responses", {
  d <- withr::with_seed(23, {
    n <- 800
    d <- data.frame(g = sample(sprintf("g%02d", 1:60), n, TRUE))
    for (i in 1:3) d[[paste0("y", i)]] <- stats::rnorm(n)
    d
  })
  f <- anova_manova(d, paste0("y", 1:3), "g", plots = FALSE)
  expect_null(f$posthoc)
  expect_null(f$emmeans_object)
  note <- grep("Pairwise comparisons were skipped", f$notes, value = TRUE)
  expect_length(note, 1L)
  # Hand count: choose(60, 2) per response, three responses.
  expect_match(note, format(3 * choose(60, 2), big.mark = ","), fixed = TRUE)
  expect_match(note, sprintf("%s for each of the 3 responses",
                             format(choose(60, 2), big.mark = ",")), fixed = TRUE)
  expect_match(note, "$univariate[[response]]$emmeans_object", fixed = TRUE)
  expect_false(grepl("on $emmeans_object", note, fixed = TRUE))
  # The place it points to holds a grid emmeans can contrast.
  expect_s4_class(f$univariate$y1$emmeans_object, "emmGrid")
})

test_that("R2-MANOVA-05: with one response the slopes note reports an exact F and no multivariate checks", {
  d <- withr::with_seed(14, {
    d <- data.frame(g = sample(c("a", "b"), 40, TRUE), x = stats::rnorm(40))
    d$y1 <- stats::rnorm(40)
    d
  })
  f <- anova_manova(d, "y1", "g", covariates = "x", plots = FALSE)
  ref <- stats::anova(stats::lm(y1 ~ g + x, d), stats::lm(y1 ~ g * x, d))
  expect_lt(ref[["Pr(>F)"]][2], 0.05)
  note <- grep("Homogeneity of regression slopes is rejected", f$notes, value = TRUE)
  expect_length(note, 1L)
  expect_match(note, sprintf("F(%g, %g) = %.4g, p = %.3g", ref$Df[2],
                             ref$Res.Df[2], ref$F[2], ref[["Pr(>F)"]][2]),
               fixed = TRUE)
  expect_false(grepl("approximate", note, fixed = TRUE))
  expect_false(grepl("Mardia|Box's M|car::Manova", note))
  expect_false(any(c("mardia", "box_m") %in% names(f$assumptions)))
})

test_that("R2-MANOVA-06: a constant covariate is refused (superseded by R3-48), not tested for slopes", {
  d <- withr::with_seed(4, {
    d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 15)), x = 1)
    d$y1 <- stats::rnorm(45) + (d$g == "b")
    d$y2 <- stats::rnorm(45)
    d
  })
  # Reference: lm() cannot estimate the covariate next to the intercept, and
  # anova.mlm compares two models of equal rank (0 degrees of freedom).
  expect_true(is.na(stats::coef(stats::lm(y1 ~ g + x, d))[["x"]]))
  ref <- stats::anova(stats::lm(cbind(y1, y2) ~ g + x, d),
                      stats::lm(cbind(y1, y2) ~ g * x, d))
  expect_equal(ref$Df[2], 0)
  expect_error(anova_manova(d, c("y1", "y2"), "g", covariates = "x", plots = FALSE),
               "Covariate `x` has the same value \\(1\\) in every analysed row, so there is nothing to adjust for")
})

test_that("R2-MANOVA-07: Type III tests on responses in small units need no fallback and match car", {
  d <- withr::with_seed(18, {
    n <- 60
    d <- data.frame(g = sample(c("a", "b", "c"), n, TRUE), x = stats::rnorm(n))
    d$y1 <- (stats::rnorm(n) + (d$g == "b") + d$x) * 1e-5
    d$y2 <- stats::rnorm(n) * 1e-5
    d
  })
  f <- anova_manova(d, c("y1", "y2"), "g", covariates = "x", type = "III",
                    plots = FALSE)
  expect_false(any(grepl("car::Manova() failed", f$notes, fixed = TRUE)))
  expect_false(any(grepl("deficient rank", f$notes, fixed = TRUE)))

  # car on the same data in ordinary units, under sum-to-zero contrasts.
  big <- d
  big$y1 <- big$y1 * 1e5
  big$y2 <- big$y2 * 1e5
  big$g <- factor(big$g)
  m <- withr::with_options(list(contrasts = c("contr.sum", "contr.poly")),
                           stats::lm(cbind(y1, y2) ~ x + g, big))
  M <- car::Manova(m, type = 3)
  pillai <- function(term) {
    ev <- Re(eigen(solve(M$SSPE, M$SSP[[term]]), only.values = TRUE)$values)
    sum(ev / (1 + ev))
  }
  mv <- f$multivariate
  expect_equal(mv$statistic[mv$term == "g"], pillai("g"), tolerance = 1e-8)
  expect_equal(mv$statistic[mv$term == "x"], pillai("x"), tolerance = 1e-8)
  # The SSCP matrices are scaled back to the data's own units.
  small <- withr::with_options(list(contrasts = c("contr.sum", "contr.poly")),
                               stats::lm(cbind(y1, y2) ~ x + factor(g), d))
  ssp <- .manova_ssp(small, "III")
  E <- crossprod(stats::residuals(small))
  expect_equal(unname(ssp$SSPE), unname(E), tolerance = 1e-8)
  expect_equal(unname(ssp$SSP[["factor(g)"]]), unname(M$SSP[["g"]]) * 1e-10,
               tolerance = 1e-8)
})
