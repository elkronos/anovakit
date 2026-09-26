# anova_welch(): Welch's one-way ANOVA with pairwise Welch t-tests.
#
# The omnibus test and every pairwise row are checked against stats::oneway.test
# and stats::t.test on the same data, so the wrapper cannot drift from the
# functions it wraps. The rest covers what the wrapper adds: per-group
# normality, the standardised mean differences (average-variance standardiser,
# Bonett's interval, the small-sample correction), the variance-ratio note, and
# combining several grouping variables into one cell factor.

test_that("the omnibus test reproduces stats::oneway.test exactly", {
  d <- fx_oneway()
  fit <- anova_welch(d, "value", "group", plots = FALSE)
  ref <- stats::oneway.test(value ~ group, data = d, var.equal = FALSE)

  expect_equal(fit$anova$statistic, unname(ref$statistic))
  expect_equal(fit$anova$num_df, unname(ref$parameter[1]))
  expect_equal(fit$anova$den_df, unname(ref$parameter[2]))
  expect_equal(fit$anova$p_value, unname(ref$p.value))
})

test_that("with two groups it reproduces the Welch t-test", {
  d <- fx_oneway(means = c(A = 5, B = 7), sds = c(A = 1, B = 2))
  fit <- anova_welch(d, "value", "group", plots = FALSE)
  ref <- stats::t.test(value ~ group, data = d, var.equal = FALSE)

  expect_equal(fit$anova$p_value, unname(ref$p.value))
  expect_equal(fit$posthoc$p_value, unname(ref$p.value))
  expect_equal(fit$posthoc$difference, unname(ref$estimate[1] - ref$estimate[2]))
  expect_equal(fit$posthoc$statistic, unname(ref$statistic))
  expect_equal(fit$posthoc$conf_low, ref$conf.int[1])
  expect_equal(fit$posthoc$conf_high, ref$conf.int[2])
})

test_that("group summaries match hand-computed means and t intervals", {
  d <- fx_oneway()
  fit <- anova_welch(d, "value", "group", conf_level = 0.9, plots = FALSE)
  a <- d$value[d$group == "A"]
  n <- length(a); se <- stats::sd(a) / sqrt(n)
  tcrit <- stats::qt(0.95, df = n - 1)

  row <- fit$emmeans[fit$emmeans$group == "A", ]
  expect_equal(row$mean, mean(a))
  expect_equal(row$sd, stats::sd(a))
  expect_equal(row$se, se)
  expect_equal(row$conf_low, mean(a) - tcrit * se)
  expect_equal(row$conf_high, mean(a) + tcrit * se)
})

test_that("Hedges' g matches the closed-form value and is named for what it is", {
  # Average-variance standardiser with the bias correction at its
  # Satterthwaite degrees of freedom (Delacre et al., 2021), and Bonett's
  # (2008) interval, written out from the published formulas.
  d <- fx_oneway(means = c(A = 0, B = 1), sds = c(A = 1, B = 2))
  fit <- anova_welch(d, "value", "group", plots = FALSE)
  x <- d$value[d$group == "A"]; y <- d$value[d$group == "B"]
  n1 <- length(x); n2 <- length(y); v1 <- stats::var(x); v2 <- stats::var(y)
  s_star <- sqrt((v1 + v2) / 2)
  d_star <- (mean(x) - mean(y)) / s_star
  df_star <- (n1 - 1) * (n2 - 1) * (v1 + v2)^2 / ((n2 - 1) * v1^2 + (n1 - 1) * v2^2)
  J <- gamma(df_star / 2) / (sqrt(df_star / 2) * gamma((df_star - 1) / 2))
  se <- sqrt(d_star^2 * (v1^2 / (n1 - 1) + v2^2 / (n2 - 1)) / (8 * s_star^4) +
               (v1 / (n1 - 1) + v2 / (n2 - 1)) / s_star^2)

  expect_true("hedges_g" %in% names(fit$effect_sizes))
  expect_false("cohen_d" %in% names(fit$effect_sizes))
  expect_equal(fit$effect_sizes$hedges_g, J * d_star)
  expect_equal(fit$effect_sizes$conf_low, d_star - stats::qnorm(0.975) * se)
  expect_equal(fit$effect_sizes$conf_high, d_star + stats::qnorm(0.975) * se)
  expect_identical(fit$effect_sizes$standardiser, "sqrt((s1^2 + s2^2) / 2)")

  uncorrected <- anova_welch(d, "value", "group", hedges_correction = FALSE,
                             plots = FALSE)
  expect_true("cohens_d" %in% names(uncorrected$effect_sizes))
  expect_equal(uncorrected$effect_sizes$cohens_d, d_star)
  skip_if_not_installed("effectsize")
  expect_equal(uncorrected$effect_sizes$cohens_d,
               effectsize::cohens_d(x, y, pooled_sd = FALSE)$Cohens_d)
})

test_that("conf_level reaches the pairwise intervals, not only the summaries", {
  d <- fx_oneway()
  wide <- anova_welch(d, "value", "group", conf_level = 0.99, plots = FALSE)
  narrow <- anova_welch(d, "value", "group", conf_level = 0.50, plots = FALSE)

  w <- wide$posthoc$conf_high - wide$posthoc$conf_low
  n <- narrow$posthoc$conf_high - narrow$posthoc$conf_low
  expect_true(all(w > n))

  we <- wide$effect_sizes$conf_high - wide$effect_sizes$conf_low
  ne <- narrow$effect_sizes$conf_high - narrow$effect_sizes$conf_low
  expect_true(all(we > ne))
})

test_that("the p-value adjustment is applied and is recorded", {
  d <- fx_oneway()
  fit <- anova_welch(d, "value", "group", adjust = "bonferroni", plots = FALSE)
  expect_equal(fit$posthoc$p_adjusted,
               stats::p.adjust(fit$posthoc$p_value, "bonferroni"))
  expect_identical(unique(fit$posthoc$adjustment), "bonferroni")
  expect_error(anova_welch(d, "value", "group", adjust = "nonsense"))
})

test_that("normality is assessed within groups, not on pooled residuals", {
  # Both groups are normal; their variances differ by a factor of 100.
  d <- fx_oneway(n_per = 40, means = c(A = 0, B = 0), sds = c(A = 1, B = 10))
  fit <- anova_welch(d, "value", "group", plots = FALSE)

  expect_equal(nrow(fit$assumptions$normality), 2L)
  expect_true(all(fit$assumptions$normality$p_value > 0.01))

  pooled <- stats::residuals(stats::lm(value ~ group, data = d))
  expect_lt(stats::shapiro.test(pooled)$p.value, 0.01)
})

test_that("a data.table input is accepted", {
  skip_if_not_installed("data.table")
  d <- fx_oneway()
  dt <- data.table::as.data.table(d)
  from_df <- anova_welch(d, "value", "group", plots = FALSE)
  from_dt <- anova_welch(dt, "value", "group", plots = FALSE)
  expect_equal(from_dt$anova, from_df$anova)
})

test_that("non-syntactic column names are handled", {
  d <- fx_oneway()
  names(d) <- c("my group", "y value")
  fit <- anova_welch(d, "y value", "my group", plots = FALSE)
  expect_fit_shape(fit)
  expect_equal(nrow(fit$emmeans), 3L)
})

test_that("empty cells are dropped rather than counted as groups", {
  d <- fx_nested()
  fit <- anova_welch(d, "value", c("g1", "g2"), plots = FALSE)

  expect_equal(nrow(fit$emmeans), 2L)
  expect_equal(nrow(fit$posthoc), 1L)
  expect_equal(fit$anova$p_value,
               unname(stats::t.test(value ~ g1, data = d)$p.value))
})

test_that("degenerate designs error with a message that names the problem", {
  d <- fx_oneway()
  d1 <- d[d$group == "A", ]
  expect_error(anova_welch(d1, "value", "group"), "at least 2 are required")

  d2 <- rbind(d[d$group != "C", ], d[d$group == "C", ][1, ])
  expect_error(anova_welch(d2, "value", "group"), "fewer than 2 observations")

  expect_error(anova_welch(d, "value", "not_a_column"), "not found in `data`")
  expect_error(anova_welch(d, "group", "group"), "must be numeric")
  expect_error(anova_welch(d, "value", "group", conf_level = 1.5),
               "strictly between 0 and 1")
})

test_that("missing values are dropped and counted", {
  d <- fx_oneway()
  d$value[c(2, 15, 40)] <- NA
  fit <- anova_welch(d, "value", "group", plots = FALSE)
  expect_identical(fit$n_removed, 3L)
  expect_equal(nrow(fit$data_used), nrow(d) - 3L)
  expect_true(any(grepl("Dropped 3 row", fit$notes)))
})

test_that("a user column called .cell does not collide with the internal one", {
  d <- fx_oneway()
  d$.cell <- "mine"
  fit <- anova_welch(d, "value", "group", plots = FALSE)
  expect_equal(nrow(fit$emmeans), 3L)
})

test_that("nothing is written to stdout and plots are returned, not drawn", {
  d <- fx_oneway()
  expect_silent_stdout(anova_welch(d, "value", "group"))
  fit <- anova_welch(d, "value", "group")
  expect_named(fit$plots, c("means", "box", "qq"))
  expect_true(all(vapply(fit$plots, inherits, logical(1), "ggplot")))
  expect_length(anova_welch(d, "value", "group", plots = FALSE)$plots, 0L)
})
