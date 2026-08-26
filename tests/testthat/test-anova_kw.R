# anova_kw(): Kruskal-Wallis with Dunn post-hoc comparisons.
#
# The omnibus test is checked against stats::kruskal.test. Dunn's test is
# implemented in this package rather than delegated, so it is checked against a
# hand computation of the tie-corrected statistic, including on data with heavy
# ties where the correction matters.

test_that("the omnibus test reproduces stats::kruskal.test exactly", {
  d <- fx_oneway()
  fit <- anova_kw(d, "value", "group", plots = FALSE)
  ref <- stats::kruskal.test(value ~ group, data = d)

  expect_equal(fit$anova$statistic, unname(ref$statistic))
  expect_equal(fit$anova$df, unname(ref$parameter))
  expect_equal(fit$anova$p_value, unname(ref$p.value))
})

test_that("Dunn's z statistics match the closed-form definition", {
  d <- fx_oneway()
  fit <- anova_kw(d, "value", "group", plots = FALSE)

  x <- d$value; g <- droplevels(factor(d$group))
  N <- length(x); r <- rank(x)
  ns <- table(g); rbar <- tapply(r, g, mean)
  ties <- table(x)
  sigma <- N * (N + 1) / 12 - sum(ties^3 - ties) / (12 * (N - 1))
  pairs <- utils::combn(levels(g), 2L, simplify = FALSE)
  z <- vapply(pairs, function(p) {
    (rbar[[p[1]]] - rbar[[p[2]]]) /
      sqrt(sigma * (1 / ns[[p[1]]] + 1 / ns[[p[2]]]))
  }, numeric(1))

  expect_equal(fit$posthoc$z, z)
  expect_equal(fit$posthoc$p_value, 2 * stats::pnorm(-abs(z)))
})

test_that("the tie correction is applied", {
  d <- withr::with_seed(42, data.frame(
    g = factor(rep(c("a", "b", "c"), each = 20)),
    y = sample(1:5, 60, replace = TRUE)))
  fit <- anova_kw(d, "y", "g", plots = FALSE)

  N <- nrow(d); r <- rank(d$y)
  ns <- table(d$g); rbar <- tapply(r, d$g, mean)
  ties <- table(d$y)
  sigma_tied <- N * (N + 1) / 12 - sum(ties^3 - ties) / (12 * (N - 1))
  sigma_naive <- N * (N + 1) / 12

  expect_lt(sigma_tied, sigma_naive)   # there really are ties here
  z <- (rbar[["a"]] - rbar[["b"]]) / sqrt(sigma_tied * (1 / ns[["a"]] + 1 / ns[["b"]]))
  expect_equal(fit$posthoc$z[1], z)
})

test_that("epsilon squared matches its definition", {
  d <- fx_oneway()
  fit <- anova_kw(d, "value", "group", plots = FALSE)
  H <- unname(stats::kruskal.test(value ~ group, data = d)$statistic)
  n <- nrow(d)
  expect_equal(fit$effect_sizes$estimate[1], H / ((n^2 - 1) / (n + 1)))
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
