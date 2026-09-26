# anova_kw(): Kruskal-Wallis with Dunn post-hoc comparisons.
#
# The omnibus test is checked against stats::kruskal.test. Dunn's test is
# implemented in this package rather than delegated, so it is checked against
# values from rstatix::dunn_test() and scikit_posthocs.posthoc_dunn(), including
# on data with heavy ties where the correction matters.

test_that("the omnibus test reproduces stats::kruskal.test exactly", {
  d <- fx_oneway()
  fit <- anova_kw(d, "value", "group", plots = FALSE)
  ref <- stats::kruskal.test(value ~ group, data = d)

  expect_equal(fit$anova$statistic, unname(ref$statistic))
  expect_equal(fit$anova$df, unname(ref$parameter))
  expect_equal(fit$anova$p_value, unname(ref$p.value))
})

test_that("Dunn's z statistics match rstatix and scikit_posthocs", {
  # Reference values, recorded because neither is a dependency:
  #   rstatix::dunn_test(fx_oneway(), value ~ group, p.adjust.method = "none")
  #   scikit_posthocs.posthoc_dunn(d, "value", "group") gives the same p.
  # rstatix reports group2 minus group1; this package group1 minus group2.
  d <- fx_oneway()
  fit <- anova_kw(d, "value", "group", adjust = "none", plots = FALSE)

  expect_equal(fit$posthoc$z, -c(4.16584609464, -2.47225445803, -6.63810055267),
               tolerance = 1e-9)
  expect_equal(fit$posthoc$p_value,
               c(3.10199929443e-05, 1.34263911514e-02, 3.17751084406e-11),
               tolerance = 1e-9)
  expect_equal(fit$posthoc$p_value, 2 * stats::pnorm(-abs(fit$posthoc$z)))
})

test_that("the tie correction is applied", {
  # Sixty draws from 1:5: heavily tied. References as above, with
  # p.adjust.method = "holm" for the adjusted values; scikit_posthocs agrees.
  d <- withr::with_seed(42, data.frame(
    g = factor(rep(c("a", "b", "c"), each = 20)),
    y = sample(1:5, 60, replace = TRUE)))
  fit <- anova_kw(d, "y", "g", adjust = "holm", plots = FALSE)

  expect_equal(fit$posthoc$z, -c(2.188227984013, 0.647156786761, -1.541071197252),
               tolerance = 1e-9)
  expect_equal(fit$posthoc$p_value,
               c(0.0286529996612, 0.5175304758187, 0.1232994589587),
               tolerance = 1e-9)
  expect_equal(fit$posthoc$p_adjusted,
               c(0.0859589989835, 0.5175304758187, 0.246598917917),
               tolerance = 1e-9)

  # Without the tie correction the first z would be smaller in magnitude.
  N <- nrow(d); rbar <- tapply(rank(d$y), d$g, mean)
  z_naive <- (rbar[["a"]] - rbar[["b"]]) / sqrt(N * (N + 1) / 12 * (2 / 20))
  expect_lt(abs(z_naive), abs(fit$posthoc$z[1]))
})

test_that("epsilon squared matches its definition", {
  d <- fx_oneway()
  fit <- anova_kw(d, "value", "group", plots = FALSE)
  H <- unname(stats::kruskal.test(value ~ group, data = d)$statistic)
  n <- nrow(d)
  expect_equal(fit$effect_sizes$estimate[1], H / ((n^2 - 1) / (n + 1)))
})

test_that("eta squared H matches rstatix::kruskal_effsize", {
  # rstatix::kruskal_effsize(d, y ~ g)$effsize on the tie-heavy data
  d <- withr::with_seed(42, data.frame(
    g = factor(rep(c("a", "b", "c"), each = 20)),
    y = sample(1:5, 60, replace = TRUE)))
  fit <- anova_kw(d, "y", "g", plots = FALSE)
  expect_identical(fit$effect_sizes$measure, c("epsilon_squared", "eta_squared_H"))
  expect_equal(fit$effect_sizes$estimate[2], 0.0535912754581, tolerance = 1e-10)
  expect_identical(sign(fit$posthoc$mean_rank_diff), sign(fit$posthoc$z))
})

test_that("six observations are enough: the optional diagnostic is optional", {
  d <- data.frame(g = rep(c("a", "b"), each = 3), y = c(1, 2, 3, 4, 5, 6))
  fit <- anova_kw(d, "y", "g", plots = FALSE)
  expect_equal(fit$anova$p_value,
               unname(stats::kruskal.test(y ~ factor(g), data = d)$p.value))
  expect_length(fit$assumptions, 0L)
})

test_that("diagnostics are off by default and available on request", {
  d <- fx_oneway()
  expect_length(anova_kw(d, "value", "group", plots = FALSE)$assumptions, 0L)

  with_diag <- anova_kw(d, "value", "group", diagnostics = TRUE, plots = FALSE)
  expect_true("normality" %in% names(with_diag$assumptions))
  expect_equal(nrow(with_diag$assumptions$normality), 3L)
  expect_true(any(grepl("contextual only", with_diag$notes)))
})

test_that("the adjustment is applied, recorded, and validated", {
  d <- fx_oneway()
  fit <- anova_kw(d, "value", "group", adjust = "bonferroni", plots = FALSE)
  expect_equal(fit$posthoc$p_adjusted,
               stats::p.adjust(fit$posthoc$p_value, "bonferroni"))
  expect_identical(unique(fit$posthoc$adjustment), "bonferroni")
  expect_error(anova_kw(d, "value", "group", adjust = "fdr_but_misspelled"))
})

test_that("boxplot labels carry the count of the cell that was tested", {
  d <- fx_twoway()
  fit <- anova_kw(d, "value", c("g1", "g2"))
  labels <- levels(fit$plots$box$data$label)

  # Four cells of 30, not two marginal groups of 60.
  expect_length(labels, 4L)
  expect_true(all(grepl("n = 30", labels)))
  expect_equal(sort(fit$emmeans$n), rep(30L, 4L))
})

test_that("a multi-factor call says what it can and cannot test", {
  d <- fx_twoway()
  fit <- anova_kw(d, "value", c("g1", "g2"), plots = FALSE)
  expect_true(any(grepl("cannot separate main effects", fit$notes)))
  expect_equal(fit$anova$df, 3)
})

test_that("user columns named N or .cell do not break the call", {
  d <- fx_oneway()
  d$N <- 999L
  d$.cell <- "mine"
  fit <- anova_kw(d, "value", "group", plots = FALSE)
  expect_equal(nrow(fit$emmeans), 3L)
})

test_that("group summaries include medians and mean ranks", {
  d <- fx_oneway()
  fit <- anova_kw(d, "value", "group", plots = FALSE)
  expect_true(all(c("median", "mean_rank") %in% names(fit$emmeans)))
  expect_equal(fit$emmeans$median[fit$emmeans$group == "A"],
               stats::median(d$value[d$group == "A"]))
  expect_equal(fit$emmeans$mean_rank,
               as.numeric(tapply(rank(d$value), d$group, mean)[
                 as.character(fit$emmeans$group)]))
})

test_that("nothing is written to stdout", {
  d <- fx_oneway()
  expect_silent_stdout(anova_kw(d, "value", "group"))
  expect_silent_stdout(anova_kw(d, "value", "group", diagnostics = TRUE))
})
