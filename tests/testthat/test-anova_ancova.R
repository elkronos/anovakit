# anova_ancova(): analysis of covariance.
#
# Covers the two decisions the function makes on the caller's behalf -- whether
# to centre the covariates and whether to keep the covariate-by-group
# interaction -- and checks that both are reported in $notes and can be
# overridden. Tables are checked against car::Anova on the same model.

test_that("covariates are centred by default and this is recorded", {
  d <- fx_ancova()
  fit <- anova_ancova(d, "dv", "iv", "cov", plots = FALSE)

  expect_equal(unname(fit$covariate_means["cov"]), mean(d$cov))
  expect_equal(mean(fit$data_used$cov), 0)
  expect_true(any(grepl("mean-centred", fit$notes)))
})

test_that("centring is what makes the group row interpretable", {
  # Slopes genuinely differ, and the covariate is centred on 100.
  d <- fx_ancova(slope_shift = 0.6)
  centred <- anova_ancova(d, "dv", "iv", "cov", plots = FALSE)
  raw <- anova_ancova(d, "dv", "iv", "cov", center_covariates = FALSE,
                      plots = FALSE)

  expect_false(centred$slopes_test$homogeneous)

  p_centred <- centred$anova$p_value[centred$anova$term == "iv"]
  p_raw <- raw$anova$p_value[raw$anova$term == "iv"]
  expect_lt(p_centred, 1e-6)
  expect_gt(p_raw, 0.05)

  # The interaction row is unaffected by where the covariate's zero sits.
  expect_equal(centred$anova$sum_sq[grepl(":", centred$anova$term)],
               raw$anova$sum_sq[grepl(":", raw$anova$term)])
  expect_true(any(grepl("NOT centred", raw$notes)))
})

test_that("effect sizes carry no residual row and no value of exactly 0.5", {
  d <- fx_ancova(slope_shift = 0.6)
  fit <- anova_ancova(d, "dv", "iv", "cov", plots = FALSE)

  expect_false("Residuals" %in% fit$effect_sizes$term)
  expect_false("(Intercept)" %in% fit$effect_sizes$term)
  expect_false(any(abs(fit$effect_sizes$partial_eta_sq - 0.5) < 1e-12))
  expect_true(all(fit$effect_sizes$partial_eta_sq >= 0 &
                    fit$effect_sizes$partial_eta_sq <= 1))
  expect_true(all(fit$effect_sizes$omega_sq >= 0))
})

test_that("partial eta squared matches its definition", {
  d <- fx_ancova(slope_shift = 0)
  fit <- anova_ancova(d, "dv", "iv", "cov", force_interaction = FALSE,
                      plots = FALSE)
  ref <- as.data.frame(car::Anova(fit$model, type = 3))
  ss_err <- ref["Residuals", "Sum Sq"]
  expected <- ref["iv", "Sum Sq"] / (ref["iv", "Sum Sq"] + ss_err)
  expect_equal(fit$effect_sizes$partial_eta_sq[fit$effect_sizes$term == "iv"],
               expected)
})

test_that("the slopes test decides the model and force_interaction overrides it", {
  homog <- fx_ancova(slope_shift = 0)
  hetero <- fx_ancova(slope_shift = 0.6)

  expect_true(anova_ancova(homog, "dv", "iv", "cov",
                           plots = FALSE)$slopes_test$homogeneous)
  expect_false(anova_ancova(hetero, "dv", "iv", "cov",
                            plots = FALSE)$slopes_test$homogeneous)

  forced_on <- anova_ancova(homog, "dv", "iv", "cov", force_interaction = TRUE,
                            plots = FALSE)
  expect_true(any(grepl(":", forced_on$anova$term)))

  forced_off <- anova_ancova(hetero, "dv", "iv", "cov",
                             force_interaction = FALSE, plots = FALSE)
  expect_false(any(grepl(":", forced_off$anova$term)))
  expect_null(forced_off$simple_slopes)
})

test_that("simple slopes are returned when the interaction is kept", {
  d <- fx_ancova(slope_shift = 0.6)
  fit <- anova_ancova(d, "dv", "iv", "cov", plots = FALSE)

  expect_false(is.null(fit$simple_slopes))
  expect_equal(nrow(fit$simple_slopes), 2L)
  expect_true("slope" %in% names(fit$simple_slopes))

  ref <- as.data.frame(suppressMessages(
    emmeans::emtrends(fit$model, ~ iv, var = "cov")))
  expect_equal(fit$simple_slopes$slope, ref$cov.trend)
  expect_true(any(grepl("comparison at the covariate mean", fit$notes)))
})

test_that("adjusted means are model based, not raw group means", {
  d <- fx_ancova(slope_shift = 0)
  fit <- anova_ancova(d, "dv", "iv", "cov", force_interaction = FALSE,
                      plots = FALSE)
  ref <- as.data.frame(suppressMessages(emmeans::emmeans(fit$model, ~ iv)))

  expect_equal(fit$emmeans$estimate, ref$emmean)
  expect_equal(fit$emmeans$se, ref$SE)
  expect_false(isTRUE(all.equal(fit$emmeans$estimate,
                                unname(tapply(d$dv, d$iv, mean)))))
})

test_that("columns not used in the model survive into data_used", {
  d <- fx_ancova()
  d$subject_id <- paste0("s", seq_len(nrow(d)))
  d$unused_all_na <- NA_real_
  fit <- anova_ancova(d, "dv", "iv", "cov", plots = FALSE)

  expect_true(all(c("subject_id", "unused_all_na") %in% names(fit$data_used)))
  # A column full of NA that the model never touches must not delete every row.
  expect_identical(fit$n_removed, 0L)
  expect_equal(nrow(fit$data_used), nrow(d))
})

test_that("Levene's test stays aligned when rows are dropped", {
  d <- fx_ancova()
  d$cov[c(3, 9, 40)] <- NA
  fit <- anova_ancova(d, "dv", "iv", "cov", plots = FALSE)

  expect_identical(fit$n_removed, 3L)
  expect_false(is.null(fit$assumptions$levene))
  expect_equal(fit$assumptions$levene$df1 + fit$assumptions$levene$df2,
               nrow(fit$data_used) - 1)
})

test_that("input validation names the problem", {
  d <- fx_ancova()
  expect_error(anova_ancova(d, "dv", "iv", "not_there"), "not found in `data`")
  expect_error(anova_ancova(d, "dv", "iv", "iv"), "must be numeric")
  expect_error(anova_ancova(d, "dv", "iv", "cov", homogeneity_alpha = 0),
               "strictly between 0 and 1")
  expect_error(anova_ancova(d, "dv", c("iv", "cov"), "cov"),
               "cannot be both a group and a covariate")
})

test_that("nothing is written to stdout and no deprecation warnings fire", {
  d <- fx_ancova()
  expect_silent_stdout(anova_ancova(d, "dv", "iv", "cov"))
  expect_no_warning(anova_ancova(d, "dv", "iv", "cov"))
  fit <- anova_ancova(d, "dv", "iv", "cov")
  expect_named(fit$plots, c("residuals", "qq", "covariate", "emmeans"))
})
