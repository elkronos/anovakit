# Regression tests for the second review's ANCOVA findings. Each result is
# checked against an independent reference: emmeans on a model fitted here
# with lm(), base R's anova(), car::Anova(), or the function's own arguments.

test_that("R2-ANCOVA-01: non-estimable group means are dropped and named, not described as extrapolations", {
  # Cell path: x is constant in groups b and c, so under group-specific slopes
  # their means at the covariate mean are not estimable.
  set.seed(3)
  g <- factor(rep(c("a", "b", "c"), each = 20))
  x <- ifelse(g == "a", rnorm(60, 5, 1), ifelse(g == "b", 2, 8))
  d <- data.frame(g = g, x = x, y = rnorm(60))
  f <- anova_ancova(d, "y", "g", "x", force_interaction = TRUE, plots = FALSE)
  m <- stats::lm(y ~ g * I(x - mean(x)), data = d)
  ref <- suppressMessages(as.data.frame(summary(emmeans::emmeans(m, ~g))))
  estimable <- as.character(ref$g[!is.na(ref$emmean)])
  expect_identical(estimable, "a")
  expect_identical(as.character(f$emmeans$g), estimable)
  expect_equal(f$emmeans$estimate, ref$emmean[!is.na(ref$emmean)],
               tolerance = 1e-8)
  expect_null(f$posthoc)
  nonest <- grep("not estimable under this model", f$notes, value = TRUE)
  expect_length(nonest, 1L)
  expect_match(nonest, "b, c", fixed = TRUE)
  expect_true(any(grepl("fewer than two marginal means are estimable",
                        f$notes, fixed = TRUE)))
  # Nothing was extrapolated for b and c, so no note may say so.
  expect_false(any(grepl("extrapolations along the fitted slope", f$notes,
                         fixed = TRUE)))

  # Per-factor path: x constant within A = a2.
  set.seed(4)
  d <- data.frame(A = factor(rep(c("a1", "a2"), each = 40)),
                  B = factor(rep(c("b1", "b2"), 40)), x = rnorm(80))
  d$x[d$A == "a2"] <- 1.5
  d$y <- rnorm(80)
  f2 <- anova_ancova(d, "y", c("A", "B"), "x", interaction = FALSE,
                     force_interaction = TRUE, plots = FALSE)
  m2 <- stats::lm(y ~ A + B + I(x - mean(x)) * (A + B), data = d)
  refA <- suppressMessages(as.data.frame(summary(emmeans::emmeans(m2, ~A))))
  refB <- suppressMessages(as.data.frame(summary(emmeans::emmeans(m2, ~B))))
  expect_identical(as.character(refA$A[!is.na(refA$emmean)]), "a1")
  expect_true(all(is.na(refB$emmean)))
  expect_false(anyNA(f2$emmeans$estimate))
  expect_identical(f2$emmeans$A[f2$emmeans$term == "A"], "a1")
  expect_equal(f2$emmeans$estimate, refA$emmean[!is.na(refA$emmean)],
               tolerance = 1e-8)
  # b1 - b2 is estimable even though the B means are not; a1 - a2 is not.
  refBp <- suppressMessages(as.data.frame(summary(pairs(emmeans::emmeans(m2, ~B)))))
  expect_false(anyNA(f2$posthoc$estimate))
  expect_identical(f2$posthoc$contrast, "b1 - b2")
  expect_equal(f2$posthoc$estimate, refBp$estimate, tolerance = 1e-8)
  nonest2 <- grep("not estimable under this model", f2$notes, value = TRUE)
  expect_length(nonest2, 1L)
  expect_match(nonest2, "A = a2, B = b1, B = b2", fixed = TRUE)
  expect_false(any(grepl("extrapolations along the fitted slope", f2$notes,
                         fixed = TRUE)))

  # The extrapolation note still names an estimable group that is
  # extrapolated under the additive model (common slope).
  fa <- anova_ancova(d, "y", c("A", "B"), "x", interaction = FALSE,
                     force_interaction = FALSE, plots = FALSE)
  expect_true(any(grepl("a2 : b1 [1.5, 1.5]", fa$notes, fixed = TRUE)))
})

test_that("R2-ANCOVA-02: estimable incomplete crossed designs are analysed, not refused", {
  set.seed(15)
  cells <- data.frame(A = rep(1:10, each = 2), B = c(rbind(1:10, c(2:10, 1))))
  d <- cells[rep(1:20, each = 4), ]
  d$A <- factor(d$A)
  d$B <- factor(d$B)
  d$x <- rnorm(80)
  d$y <- d$x + rnorm(80)
  m <- stats::lm(y ~ x + A * B, data = d)
  expect_identical(m$rank, 21L)
  ref <- as.data.frame(suppressMessages(car::Anova(m, type = 2)))
  f <- anova_ancova(d, "y", c("A", "B"), "x", force_interaction = FALSE,
                    plots = FALSE)
  expect_s3_class(f, "anovakit_fit")
  expect_identical(stats::df.residual(f$model), stats::df.residual(m))
  for (tm in c("x", "A", "B", "A:B")) {
    expect_equal(f$anova$statistic[f$anova$term == tm],
                 ref[tm, "F value"], tolerance = 1e-6, info = tm)
    expect_equal(f$anova$df[f$anova$term == tm], ref[tm, "Df"], info = tm)
  }
  # 20 of the 100 cells are observed; the empty ones are left out.
  expect_identical(nrow(f$emmeans), 20L, info = paste(f$notes, collapse = " | "))
  # A design that really has too few rows is still refused, with the
  # observed-cell count as the reason.
  one <- d[!duplicated(d[c("A", "B")]), ]
  expect_error(anova_ancova(one, "y", c("A", "B"), "x", plots = FALSE),
               "21 estimable coefficients but there are only 20 rows")
})

test_that("R2-ANCOVA-03: a saturated slopes model gives NA, not NaN, and says why", {
  d <- data.frame(g = factor(c("a", "a", "b", "b")), x = c(1, 2, 1, 3),
                  y = c(1, 2, 2, 5))
  ref <- stats::anova(stats::lm(y ~ g + x, d), stats::lm(y ~ g * x, d))
  expect_identical(ref[["Res.Df"]][2L], 0)
  expect_true(is.nan(ref[["F"]][2L]))
  f <- anova_ancova(d, "y", "g", "x", plots = FALSE)
  expect_true(is.na(f$slopes_test$statistic) && !is.nan(f$slopes_test$statistic))
  expect_true(is.na(f$slopes_test$p_value) && !is.nan(f$slopes_test$p_value))
  expect_identical(f$slopes_test$df, ref[["Df"]][2L])
  expect_true(is.na(f$slopes_test$homogeneous))
  expect_true(any(grepl("no residual degrees of freedom", f$notes, fixed = TRUE)))
  expect_false(any(grepl("Check that the covariate varies", f$notes, fixed = TRUE)))
  # The same for a larger saturated design.
  set.seed(3)
  d <- data.frame(g = factor(rep(1:20, each = 3)), x = rnorm(60),
                  x2 = rnorm(60), x3 = rnorm(60))
  d$y <- d$x + rnorm(60)
  ref <- stats::anova(stats::lm(y ~ g + x + x2 + x3, d),
                      stats::lm(y ~ (x + x2 + x3) * g, d))
  expect_identical(ref[["Res.Df"]][2L], 0)
  st <- anova_ancova(d, "y", "g", c("x", "x2", "x3"), plots = FALSE)$slopes_test
  expect_false(is.nan(st$statistic) || is.nan(st$p_value))
  expect_true(is.na(st$p_value))
})

test_that("R2-ANCOVA-04: the Levene note does not claim robust $anova when the robust tests fell back", {
  skip_if_not_installed("sandwich")
  set.seed(2)
  d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 30)),
                  x = rnorm(90, 50, 10))
  d$y <- 0.5 * d$x + rnorm(90)
  d$x[d$g == "c"] <- 60
  f <- anova_ancova(d, "y", "g", "x", force_interaction = TRUE,
                    vcov_type = "HC3", plots = FALSE)
  # Independent reference: car refuses a user covariance for this aliased model.
  m <- stats::lm(y ~ g * I(x - mean(x)), data = d)
  expect_true(any(is.na(stats::coef(m))))
  expect_error(suppressMessages(
    car::Anova(m, type = 2, vcov. = sandwich::vcovHC(m, type = "HC3"))))
  expect_identical(attr(f$anova, "statistic"), "F")
  expect_false(any(grepl("\n", f$notes, fixed = TRUE)))
  expect_true(any(grepl("Robust (HC3) Wald F tests could not be computed",
                        f$notes, fixed = TRUE)))
  lev <- grep("^Levene's test rejects", f$notes, value = TRUE)
  if (length(lev) > 0L) {
    expect_false(grepl("already allows for this in $anova", lev, fixed = TRUE))
    expect_match(lev, "F tests in $anova", fixed = TRUE)
  }
  # Unit check of the wording, independent of whether Levene rejects here.
  d2 <- data.frame(g = factor(rep(c("a", "b"), each = 40)),
                   y = c(rnorm(40, 0, 1), rnorm(40, 0, 6)))
  mm <- stats::lm(y ~ g, data = d2)
  stopifnot(car::leveneTest(stats::residuals(mm) ~ d2$g)[["Pr(>F)"]][1L] < 0.05)
  lv <- .levene_on_model(mm, d2, "g", vcov_type = "HC3", robust_anova = FALSE)
  expect_false(any(grepl("already allows for this in $anova", lv$note, fixed = TRUE)))
  lv2 <- .levene_on_model(mm, d2, "g", vcov_type = "HC3", robust_anova = TRUE)
  expect_true(any(grepl("already allows for this in $anova", lv2$note, fixed = TRUE)))
})

test_that("R2-ANCOVA-05: the robust-SE note in an ANCOVA does not recommend anova_welch()", {
  skip_if_not_installed("sandwich")
  # anova_welch() has no covariate argument, so it cannot replace an ANCOVA.
  expect_false(any(grepl("covariate", names(formals(anova_welch)))))
  set.seed(9)
  d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 30)),
                  x = rnorm(90, 50, 10))
  d$y <- 0.5 * d$x + rnorm(90)
  f <- anova_ancova(d, "y", "g", "x", vcov_type = "HC3", plots = FALSE)
  expect_false(any(grepl("anova_welch", f$notes, fixed = TRUE)))
  expect_true(any(grepl("robust standard errors keep the residual degrees of freedom",
                        f$notes, fixed = TRUE)))
  # Where there is no covariate the hint is still given.
  expect_true(any(grepl("anova_welch",
                        .robust_vcov(stats::lm(y ~ g, d), "HC3")$note, fixed = TRUE)))
})
