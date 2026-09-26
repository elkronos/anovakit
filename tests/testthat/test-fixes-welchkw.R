# Regression tests for the anova_welch() and anova_kw() findings of the
# adversarial review. Each is checked against an independent reference: base
# R (t.test, oneway.test, kruskal.test, lm), effectsize, a hand computation
# from the published formula, or constants from rstatix::dunn_test(),
# rstatix::kruskal_effsize(), scikit_posthocs.posthoc_dunn() and
# scipy.stats.kruskal() (neither rstatix nor Python is a dependency, so their
# results are recorded here with the call that produced them).

# Local fixtures -------------------------------------------------------------

wk_two <- function(n = c(10, 50), sds = c(4, 1), seed = 4) {
  withr::with_seed(seed, data.frame(
    g = rep(c("A", "B"), n),
    y = c(stats::rnorm(n[1], 0, sds[1]), stats::rnorm(n[2], 0, sds[2]))))
}

wk_three <- function(seed = 11) {
  withr::with_seed(seed, {
    d <- data.frame(g = rep(c("a", "b", "c", "d"), each = 15))
    d$y <- stats::rnorm(60) + 0.7 * (d$g == "b")
    d
  })
}

# The tie-heavy fixture: 60 draws from 1:5 in three groups of 20.
wk_ties <- function() {
  data.frame(
    g = factor(rep(c("a", "b", "c"), each = 20)),
    y = c(1, 5, 1, 1, 2, 4, 2, 2, 1, 4, 1, 5, 4, 2, 2, 3, 1, 1, 3, 4, 5, 5,
          5, 4, 2, 4, 3, 2, 1, 2, 3, 2, 4, 4, 2, 5, 4, 5, 4, 2, 2, 3, 1, 5,
          2, 2, 2, 4, 3, 5, 2, 2, 2, 5, 1, 1, 4, 5, 2, 1))
}

# Hand computation of the average-variance standardised difference, written
# in Delacre et al.'s (2021) and Bonett's (2008) own notation.
wk_smd_ref <- function(x, y, level = 0.95) {
  n1 <- length(x); n2 <- length(y)
  v1 <- stats::var(x); v2 <- stats::var(y)
  s_star <- sqrt((v1 + v2) / 2)
  d <- (mean(x) - mean(y)) / s_star
  df_star <- (n1 - 1) * (n2 - 1) * (v1 + v2)^2 / ((n2 - 1) * v1^2 + (n1 - 1) * v2^2)
  J <- gamma(df_star / 2) / (sqrt(df_star / 2) * gamma((df_star - 1) / 2))
  se <- sqrt(d^2 * (v1^2 / (n1 - 1) + v2^2 / (n2 - 1)) / (8 * s_star^4) +
               (v1 / (n1 - 1) + v2 / (n2 - 1)) / s_star^2)
  z <- stats::qnorm(1 - (1 - level) / 2)
  list(d = d, g = J * d, low = d - z * se, high = d + z * se)
}

# F052 / F053: standardised differences valid under unequal variances ---------

test_that("the standardised difference uses the average-variance SD and Bonett's interval (F052)", {
  d <- wk_two()
  x <- d$y[d$g == "A"]; y <- d$y[d$g == "B"]
  ref <- wk_smd_ref(x, y)

  fit_d <- anova_welch(d, "y", "g", hedges_correction = FALSE, plots = FALSE)
  expect_equal(fit_d$effect_sizes$cohens_d, ref$d)
  expect_equal(fit_d$effect_sizes$conf_low, ref$low)
  expect_equal(fit_d$effect_sizes$conf_high, ref$high)
  expect_identical(fit_d$effect_sizes$standardiser, "sqrt((s1^2 + s2^2) / 2)")

  fit_g <- anova_welch(d, "y", "g", plots = FALSE)
  expect_equal(fit_g$effect_sizes$hedges_g, ref$g)
  # the interval is for the population value, the same with or without J
  expect_equal(fit_g$effect_sizes$conf_low, ref$low)
  expect_equal(fit_g$effect_sizes$conf_high, ref$high)

  # Before the fix: g = 1.13 [0.43, 1.83] beside a Welch p of 0.147. Now the
  # effect-size interval agrees with the Welch test in the same row.
  tt <- stats::t.test(x, y)
  expect_gt(tt$p.value, 0.05)
  expect_lt(fit_g$effect_sizes$conf_low, 0)
  expect_gt(fit_g$effect_sizes$conf_high, 0)

  skip_if_not_installed("effectsize")
  es <- effectsize::cohens_d(x, y, pooled_sd = FALSE, ci = 0.95)
  expect_equal(fit_d$effect_sizes$cohens_d, es$Cohens_d)
  expect_lt(abs(fit_g$effect_sizes$hedges_g -
                  effectsize::hedges_g(x, y, pooled_sd = FALSE)$Hedges_g), 0.02)
})

test_that("the effect-size interval covers close to nominal under unequal variances (F052)", {
  withr::local_seed(1)
  smd <- anovakit:::.std_mean_diff
  cover <- function(n1, n2, sd1, sd2, delta, R = 1000) {
    s_bar <- sqrt((sd1^2 + sd2^2) / 2)
    mean(replicate(R, {
      e <- smd(stats::rnorm(n1, delta * s_bar, sd1), stats::rnorm(n2, 0, sd2),
               0.95, TRUE)
      e$conf_low <= delta && e$conf_high >= delta
    }))
  }
  # The pooled-SD interval covered 81% here and 67% with a variance ratio of 16.
  expect_gt(cover(10, 40, 2, 1, 0), 0.925)
  expect_gt(cover(10, 50, 4, 1, 0), 0.925)
  expect_gt(cover(10, 40, 2, 1, 0.8), 0.925)
  expect_lt(cover(40, 10, 4, 1, 0.8), 0.975)
})

test_that("small groups get an interval that does not undercover (F053)", {
  d <- data.frame(g = rep(c("A", "B"), each = 2), y = c(1, 2, 3.5, 5))
  fit <- anova_welch(d, "y", "g", plots = FALSE, hedges_correction = FALSE)
  ref <- wk_smd_ref(d$y[1:2], d$y[3:4])
  expect_equal(fit$effect_sizes$conf_low, ref$low)
  expect_equal(fit$effect_sizes$conf_high, ref$high)
  # Welch's CI includes 0 (p = 0.11); so must the standardised one now.
  expect_gt(fit$posthoc$p_value, 0.05)
  expect_lt(fit$effect_sizes$conf_low, 0)
  expect_gt(fit$effect_sizes$conf_high, 0)

  withr::local_seed(99)
  smd <- anovakit:::.std_mean_diff
  hit <- replicate(2000, {
    e <- smd(stats::rnorm(3, 1), stats::rnorm(3), 0.95, FALSE)
    e$conf_low <= 1 && e$conf_high >= 1
  })
  expect_gte(mean(hit), 0.94)       # the z interval on Var(d) gave 0.926
})

# F054: the Q-Q plot matches the per-group normality premise -----------------

test_that("the Welch Q-Q plot shows within-group standardised deviations (F054)", {
  d <- withr::with_seed(101, data.frame(
    group = factor(rep(c("A", "B"), each = 40)),
    value = c(stats::rnorm(40, 0, 1), stats::rnorm(40, 0, 10))))
  fit <- anova_welch(d, "value", "group")
  plotted <- fit$plots$qq$data$value
  ref <- unlist(lapply(split(d$value, d$group), function(v) as.numeric(scale(v))))
  expect_equal(sort(plotted), sort(unname(ref)))
  # not the pooled lm residuals, the mixture the docs reject
  pooled <- stats::residuals(stats::lm(value ~ group, data = d))
  expect_false(isTRUE(all.equal(sort(plotted), sort(unname(pooled)))))
  expect_gt(stats::shapiro.test(plotted)$p.value, 0.05)
  expect_match(fit$plots$qq$labels$subtitle, "group SD", fixed = TRUE)

  kw <- anova_kw(d, "value", "group", diagnostics = TRUE)
  expect_equal(sort(kw$plots$qq$data$value), sort(unname(ref)))
})

# F136: adjust is validated like everywhere else -----------------------------

test_that("`adjust` is matched exactly, with the documented defaults (F136)", {
  d <- wk_three()
  glm_msg <- tryCatch(anova_glm(d, "y", "g", adjust = "bon", plots = FALSE),
                      error = conditionMessage)
  expect_match(glm_msg, "`adjust` must be one of")
  for (f in list(anova_welch, anova_kw)) {
    expect_error(f(d, "y", "g", adjust = "bon", plots = FALSE),
                 "`adjust` must be one of: holm, hochberg, hommel, bonferroni, BH, BY, fdr, none; got \"bon\".",
                 fixed = TRUE)
    expect_error(f(d, "y", "g", adjust = NULL, plots = FALSE),
                 "`adjust` must be a single character string.", fixed = TRUE)
    expect_error(f(d, "y", "g", adjust = "tukey", plots = FALSE),
                 "`adjust` must be one of", fixed = TRUE)
    expect_error(f(d, "y", "g", adjust = c("holm", "BH"), plots = FALSE),
                 "`adjust` must be a single character string.", fixed = TRUE)
  }
  expect_identical(unique(anova_welch(d, "y", "g", plots = FALSE)$posthoc$adjustment),
                   "holm")
  expect_identical(unique(anova_kw(d, "y", "g", plots = FALSE)$posthoc$adjustment),
                   "BH")
  expect_identical(formals(anova_welch)$adjust, "holm")
  expect_identical(formals(anova_kw)$adjust, "BH")
})

test_that("a column cannot be both the response and a group, in either function", {
  d <- wk_three()
  d$y2 <- d$y
  expect_error(anova_welch(d, "y", c("g", "y")), "both the response and a grouping variable")
  expect_error(anova_kw(d, "y", c("g", "y")), "both the response and a grouping variable")
  names(d)[names(d) == "g"] <- "Residuals"
  expect_error(anova_kw(d, "y", "Residuals"), "reserved")
  expect_error(anova_welch(d, "y", "Residuals"), "reserved")
})

# F167: several grouping columns are labelled as combined cells ---------------

test_that("combined cells are labelled as cells, with a note, in both functions (F167)", {
  d <- withr::with_seed(4, {
    dd <- expand.grid(A = paste0("a", 1:3), B = paste0("b", 1:2), rep = 1:15)
    dd$y <- 1.0 * (dd$A == "a3") + stats::rnorm(nrow(dd))
    dd
  })
  w <- anova_welch(d, "y", c("A", "B"), plots = FALSE)
  k <- anova_kw(d, "y", c("A", "B"), plots = FALSE)
  expect_identical(w$anova$term, "A x B cells")
  expect_identical(k$anova$term, "A x B cells")
  expect_true(any(grepl("combined into 6 cells.*cannot separate main effects", w$notes)))
  expect_true(any(grepl("combined into 6 cells.*cannot separate main effects", k$notes)))

  cells <- interaction(d$A, d$B)
  expect_equal(w$anova$statistic,
               unname(stats::oneway.test(d$y ~ cells)$statistic))
  expect_equal(k$anova$statistic, unname(stats::kruskal.test(d$y, cells)$statistic))
  expect_false(any(grepl("combined", anova_welch(d, "y", "A", plots = FALSE)$notes)))
})

# F055: the posthoc = FALSE note names the right quantity ----------------------

test_that("the posthoc = FALSE note names Cohen's d when that is what would be computed (F055)", {
  d <- withr::with_seed(3, data.frame(g = rep(c("a", "b", "c"), each = 6),
                                      y = stats::rnorm(18)))
  off_d <- anova_welch(d, "y", "g", hedges_correction = FALSE, posthoc = FALSE,
                       plots = FALSE)$notes
  off_g <- anova_welch(d, "y", "g", posthoc = FALSE, plots = FALSE)$notes
  expect_true(any(grepl("Cohen's d is a per-pair quantity", off_d, fixed = TRUE)))
  expect_false(any(grepl("Hedges' g", off_d, fixed = TRUE)))
  expect_true(any(grepl("Hedges' g is a per-pair quantity", off_g, fixed = TRUE)))
})

# F001: the pairwise table's values, not only its shape -----------------------

test_that("difference, statistic, interval, variance ratio and labels are pinned (F001)", {
  d <- wk_three()
  fit <- anova_welch(d, "y", "g", adjust = "bonferroni", plots = FALSE)
  for (i in seq_len(nrow(fit$posthoc))) {
    row <- fit$posthoc[i, ]
    x <- d$y[d$g == row$group1]; y <- d$y[d$g == row$group2]
    tt <- stats::t.test(x, y, var.equal = FALSE)
    expect_equal(row$difference, mean(x) - mean(y))
    expect_equal(row$statistic, unname(tt$statistic))
    expect_equal(row$df, unname(tt$parameter))
    expect_true(row$conf_low <= row$difference && row$difference <= row$conf_high)
    ref <- wk_smd_ref(x, y)
    es <- fit$effect_sizes[i, ]
    expect_equal(es$hedges_g, ref$g)
    expect_equal(c(es$conf_low, es$conf_high), c(ref$low, ref$high))
    expect_identical(sign(es$hedges_g), sign(row$difference))
  }
  v <- tapply(d$y, d$g, stats::var)
  expect_equal(fit$assumptions$variance_ratio, max(v) / min(v))

  mag <- anovakit:::.effect_magnitude(c(0.19, 0.2, -0.49, 0.5, 0.79, 0.8, -3, NA))
  expect_identical(mag, c("negligible", "small", "small", "medium", "medium",
                          "large", "large", NA))
})

test_that("the variance-ratio note fires on the variance ratio, not the SD ratio (F001)", {
  d <- withr::with_seed(12, data.frame(
    g = rep(c("a", "b"), each = 40),
    y = c(stats::rnorm(40, 0, 1), stats::rnorm(40, 0, 3))))
  v <- tapply(d$y, d$g, stats::var)
  expect_gt(max(v) / min(v), 4)
  fit <- anova_welch(d, "y", "g", plots = FALSE)
  expect_true(any(grepl(sprintf("%.1f times the smallest", max(v) / min(v)),
                        fit$notes, fixed = TRUE)))
})

# F004: skipped Shapiro-Wilk groups are named, with the reason ---------------

test_that("a group whose Shapiro-Wilk test was skipped is named in $notes (F004)", {
  d <- withr::with_seed(5, data.frame(g = rep(c("A", "B", "C"), c(6000, 50, 30)),
                                      y = stats::rnorm(6080)))
  fit <- anova_welch(d, "y", "g", plots = FALSE)
  expect_true(is.na(fit$assumptions$normality$p_value[1]))
  expect_equal(fit$assumptions$normality$p_value[2],
               stats::shapiro.test(d$y[d$g == "B"])$p.value)
  note <- grep("Shapiro-Wilk was not run", fit$notes, value = TRUE)
  expect_length(note, 1L)
  expect_match(note, "group \"A\"", fixed = TRUE)
  expect_match(note, "5000", fixed = TRUE)
  expect_false(grepl("\"B\"", note, fixed = TRUE))
})

# F049: one p-value convention across the package ------------------------------

test_that("the posthoc tables carry raw p_value, p_adjusted and adjustment (F049)", {
  d <- wk_three()
  w <- anova_welch(d, "y", "g", adjust = "bonferroni", plots = FALSE)
  expect_named(w$posthoc, c("group1", "group2", "difference", "conf_low",
                            "conf_high", "statistic", "df", "p_value",
                            "p_adjusted", "adjustment"))
  tt <- stats::t.test(d$y[d$g == "a"], d$y[d$g == "b"])
  expect_equal(w$posthoc$p_value[1], tt$p.value)
  expect_equal(w$posthoc$p_adjusted[1], min(1, 6 * tt$p.value))
  # the interval is the unadjusted per-comparison one
  expect_equal(c(w$posthoc$conf_low[1], w$posthoc$conf_high[1]), tt$conf.int[1:2])

  k <- anova_kw(d, "y", "g", adjust = "bonferroni", plots = FALSE)
  expect_named(k$posthoc, c("group1", "group2", "mean_rank_diff", "z",
                            "p_value", "p_adjusted", "adjustment"))
  expect_equal(k$posthoc$p_adjusted, pmin(1, 6 * k$posthoc$p_value))

  g <- anova_glm(d, "y", "g", adjust = "bonferroni", plots = FALSE)
  pos <- match(c("p_value", "p_adjusted", "adjustment"), names(g$posthoc))
  expect_identical(diff(pos), c(1L, 1L))
  expect_identical(unique(g$posthoc$adjustment), "bonferroni")
})

# F007: an empty-string group level ------------------------------------------

test_that("a group level that is the empty string is analysed like any other (F007)", {
  d <- data.frame(g = rep(c("", "b", "c"), each = 5),
                  y = c(1, 4, 2, 7, 3, 8, 5, 6, 12, 9, 13, 10, 15, 11, 14))
  fit <- anova_kw(d, "y", "g", plots = FALSE)
  expect_equal(fit$anova$p_value, stats::kruskal.test(y ~ factor(g), data = d)$p.value)
  # rstatix::dunn_test(p.adjust.method = "none") with "" recoded to "zz":
  # b-c 1.626345597, b-zz -1.626345597, c-zz -3.252691193 (group2 - group1)
  ph <- fit$posthoc
  expect_equal(ph$z[ph$group1 == "" & ph$group2 == "b"], -1.626345597, tolerance = 1e-8)
  expect_equal(ph$z[ph$group1 == "" & ph$group2 == "c"], -3.252691193, tolerance = 1e-8)
  expect_equal(ph$z[ph$group1 == "b" & ph$group2 == "c"], -1.626345597, tolerance = 1e-8)
  expect_equal(ph$p_value[ph$group1 == "" & ph$group2 == "c"], 0.001143176597,
               tolerance = 1e-8)

  em <- fit$emmeans
  expect_equal(em$median[em$group == ""], stats::median(c(1, 4, 2, 7, 3)))
  expect_equal(em$mean_rank[em$group == ""], mean(rank(d$y)[d$g == ""]))
  expect_false(anyNA(em$median))
  expect_false(anyNA(em$mean_rank))

  w <- anova_welch(d, "y", "g", plots = FALSE)
  expect_equal(w$posthoc$difference[w$posthoc$group1 == "" & w$posthoc$group2 == "b"],
               mean(d$y[d$g == ""]) - mean(d$y[d$g == "b"]))

  # straight from a CSV with blank group cells
  tf <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("g,y", ",1", ",2", "a,3", "a,4", "b,5", "b,6"), tf)
  expect_s3_class(anova_kw(utils::read.csv(tf), "y", "g", plots = FALSE),
                  "anovakit_fit")
})

# F008: ties counted on the ranks ---------------------------------------------

test_that("ties are counted on the ranks, in the omnibus test and in Dunn's (F008)", {
  site1 <- c(1, 2, 3, 3, 3, 4, 7) / 10
  site2 <- rep(cumsum(rep(0.1, 3))[3], 6)   # prints as 0.3, is 0.30000000000000004
  d <- data.frame(g = rep(c("a", "b", "c"), c(7, 6, 7)),
                  y = c(site1, site2, c(3, 3, 3, 5, 6, 8, 9) / 10))
  fit <- anova_kw(d, "y", "g", plots = FALSE)
  # scipy.stats.kruskal (ties counted on the ranks): H = 3.050340136054415,
  # p = 0.21758404962837605; stats::kruskal.test(y) gives 3.681445.
  expect_equal(fit$anova$statistic, 3.050340136054415, tolerance = 1e-10)
  expect_equal(fit$anova$p_value, 0.21758404962837605, tolerance = 1e-10)
  # rstatix::dunn_test(p.adjust.method = "none") and scikit_posthocs agree:
  # p = 0.2118154570, 0.0947455742, 0.7212866827; rstatix statistic 1.2485889472
  expect_equal(fit$posthoc$p_value, c(0.2118154570, 0.0947455742, 0.7212866827),
               tolerance = 1e-8)
  expect_equal(fit$posthoc$z[1], -1.2485889472, tolerance = 1e-8)
  expect_true(any(grepl("ranked as distinct values", fit$notes)))
  # data without such values carry no such note
  expect_false(any(grepl("ranked as distinct", anova_kw(wk_ties(), "y", "g",
                                                          plots = FALSE)$notes)))
})

# F009: degenerate responses --------------------------------------------------

test_that("an all-identical response is refused, and eta squared H is NA when n = k (F009)", {
  flat <- data.frame(g = rep(c("a", "b", "c"), each = 5), y = 3)
  expect_error(anova_kw(flat, "y", "g", plots = FALSE),
               "every one of the 15 observations of `y` takes the same value")

  one <- data.frame(g = c("a", "b", "c"), y = 1:3)
  fit <- anova_kw(one, "y", "g", plots = FALSE)
  H <- unname(stats::kruskal.test(y ~ factor(g), data = one)$statistic)
  expect_equal(fit$effect_sizes$estimate[1], H / (3 - 1))
  expect_true(is.na(fit$effect_sizes$estimate[2]))
  expect_false(any(is.nan(fit$effect_sizes$estimate)))
  expect_true(any(grepl("eta_squared_H is undefined", fit$notes)))
})

# F058 / F093 / F111: the rank effect sizes ------------------------------------

test_that("the rank effect sizes match effectsize, rstatix and the rank ANOVA (F058, F111)", {
  d <- wk_ties()
  fit <- anova_kw(d, "y", "g", plots = FALSE)
  es <- stats::setNames(fit$effect_sizes$estimate, fit$effect_sizes$measure)
  # rstatix::kruskal_effsize(d, y ~ g)$effsize = 0.0535912754581 (eta2[H])
  expect_equal(unname(es["eta_squared_H"]), 0.0535912754581, tolerance = 1e-10)
  # the documented identities: R^2 and adjusted R^2 of an ANOVA on the ranks
  s <- summary(stats::lm(rank(y) ~ g, data = d))
  expect_equal(unname(es["epsilon_squared"]), s$r.squared)
  expect_equal(unname(es["eta_squared_H"]), s$adj.r.squared)
  expect_gte(es[["epsilon_squared"]], es[["eta_squared_H"]])

  skip_if_not_installed("effectsize")
  expect_equal(unname(es["epsilon_squared"]),
               effectsize::rank_epsilon_squared(y ~ g, data = d, ci = NULL)[[1]])
  expect_equal(unname(es["eta_squared_H"]),
               effectsize::rank_eta_squared(y ~ g, data = d, ci = NULL)[[1]])
})

test_that("a negative eta squared H is floored with a note giving its value (F093)", {
  d <- withr::with_seed(2, data.frame(g = factor(rep(c("A", "B", "C", "D"), each = 10)),
                                      y = stats::rnorm(40)))
  fit <- anova_kw(d, "y", "g", plots = FALSE, posthoc = FALSE)
  # rstatix::kruskal_effsize(d, y ~ g)$effsize = -0.0192276422764
  expect_identical(fit$effect_sizes$estimate[2], 0)
  note <- grep("eta_squared_H was negative", fit$notes, value = TRUE)
  expect_length(note, 1L)
  expect_match(note, "-0.0192", fixed = TRUE)
})

# F059 / F111: direction of Dunn's differences --------------------------------

test_that("mean_rank_diff and z are group1 minus group2 (F059, F111)", {
  d <- wk_ties()
  fit <- anova_kw(d, "y", "g", adjust = "holm", plots = FALSE)
  rbar <- tapply(rank(d$y), d$g, mean)
  expect_equal(fit$posthoc$mean_rank_diff,
               unname(c(rbar["a"] - rbar["b"], rbar["a"] - rbar["c"],
                        rbar["b"] - rbar["c"])))
  expect_identical(sign(fit$posthoc$mean_rank_diff), sign(fit$posthoc$z))
  # rstatix::dunn_test(d, y ~ g, p.adjust.method = "holm"): statistic
  # 2.1882279840, 0.6471567868, -1.5410711973 (group2 minus group1)
  expect_equal(fit$posthoc$z, -c(2.1882279840, 0.6471567868, -1.5410711973),
               tolerance = 1e-9)
})

# F057: infinite responses are ranks, not missing ------------------------------

test_that("infinite responses are kept for the rank test and reported honestly (F057)", {
  d <- data.frame(g = rep(c("a", "b"), each = 5), y = c(1:5, 6, 7, Inf, Inf, Inf))
  fit <- anova_kw(d, "y", "g", plots = FALSE)
  # scipy.stats.kruskal and stats::kruskal.test both keep them: p = 0.008207736107
  expect_equal(fit$anova$p_value, stats::kruskal.test(y ~ g, data = d)$p.value)
  expect_equal(fit$anova$p_value, 0.008207736107348873, tolerance = 1e-10)
  expect_identical(fit$n_removed, 0L)
  expect_identical(fit$emmeans$n, c(5L, 5L))
  b <- fit$emmeans[fit$emmeans$group == "b", ]
  expect_identical(b$mean, Inf)
  expect_true(is.na(b$sd) && is.na(b$conf_low) && is.na(b$conf_high))
  expect_equal(b$median, Inf)
  expect_equal(b$mean_rank, mean(rank(d$y)[6:10]))
  expect_equal(fit$emmeans$sd[1], stats::sd(1:5))
  expect_true(any(grepl("infinite value(s) of `y` were kept", fit$notes, fixed = TRUE)))

  both <- data.frame(g = rep(c("a", "b"), each = 4), y = c(1:4, -Inf, 5, 6, Inf))
  f2 <- anova_kw(both, "y", "g", plots = FALSE)
  expect_true(is.na(f2$emmeans$mean[2]))
  expect_false(any(is.nan(unlist(f2$emmeans[, -1]))))

  # a missing value is still dropped and counted
  d$y[1] <- NA
  expect_identical(anova_kw(d, "y", "g", plots = FALSE)$n_removed, 1L)
})

# F060: the stored htest names the user's columns --------------------------------

test_that("$model$data.name names the response and the grouping column (F060)", {
  d <- withr::with_seed(123, data.frame(group = rep(c("A", "B", "C"), each = 25),
                                        value = stats::rnorm(75)))
  expect_identical(anova_kw(d, "value", "group", plots = FALSE)$model$data.name,
                   "value by group")
  expect_false(grepl(".cell", anova_welch(d, "value", "group",
                                          plots = FALSE)$model$data.name,
                     fixed = TRUE))
})
