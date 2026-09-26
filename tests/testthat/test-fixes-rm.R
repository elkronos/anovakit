# anova_rm(): regression tests for the findings of the adversarial review.
#
# Every expectation is checked against an independent reference: afex, car,
# effectsize, emmeans, stats::lm(), stats::shapiro.test(), pingouin (values
# hard-coded from pingouin 0.5 rm_anova() / epsilon()), or a hand computation
# from the published formula.

rm_oneway <- function(n = 12, k = 3, seed = 1, shift = 0.3) {
  withr::with_seed(seed, {
    d <- expand.grid(id = factor(seq_len(n)),
                     time = factor(paste0("t", seq_len(k))))
    d$y <- stats::rnorm(n)[d$id] + shift * as.numeric(d$time) +
      stats::rnorm(nrow(d))
    d
  })
}

# Wide matrix of a one-way within design, subjects in rows
rm_wide <- function(d, y = "y", id = "id", time = "time") {
  w <- stats::reshape(d[c(id, time, y)], idvar = id, timevar = time,
                      direction = "wide")
  as.matrix(w[, -1L])
}

# Greenhouse-Geisser epsilon by hand from orthonormal Helmert contrasts
gg_by_hand <- function(W) {
  k <- ncol(W)
  C <- stats::contr.helmert(k)
  C <- sweep(C, 2L, sqrt(colSums(C^2)), "/")
  S <- stats::cov(W %*% C)
  sum(diag(S))^2 / ((k - 1) * sum(diag(S %*% S)))
}

# F018 ------------------------------------------------------------------------

test_that("generalised eta squared honours `observed` (F018)", {
  skip_if_not_installed("afex")
  skip_if_not_installed("effectsize")
  d <- withr::with_seed(7, {
    d <- expand.grid(id = factor(1:20), time = factor(c("t1", "t2", "t3")))
    d$age_grp <- factor(ifelse(as.integer(d$id) <= 10, "young", "old"))
    d$y <- stats::rnorm(20)[d$id] + 6 * (d$age_grp == "old") +
      as.numeric(d$time) + stats::rnorm(nrow(d))
    d
  })
  fit <- anova_rm(d, "y", subject = "id", within = "time", between = "age_grp",
                  observed = "age_grp", plots = FALSE, posthoc = FALSE)
  ref <- suppressMessages(afex::aov_ez("id", "y", d, within = "time",
                                       between = "age_grp",
                                       observed = "age_grp"))
  es <- effectsize::eta_squared(ref, generalized = "age_grp", verbose = FALSE)
  expect_equal(fit$effect_sizes$generalised_eta_sq,
               es$Eta2_generalized[match(fit$effect_sizes$term, es$Parameter)])
  expect_equal(fit$effect_sizes$generalised_eta_sq,
               as.data.frame(ref$anova_table)[["ges"]])
  expect_equal(round(fit$effect_sizes$generalised_eta_sq, 4),
               c(0.8034, 0.0551, 0.0041))
  expect_identical(attr(fit$effect_sizes, "observed"), "age_grp")

  # Without `observed` the values differ: every factor is manipulated
  manip <- anova_rm(d, "y", subject = "id", within = "time",
                    between = "age_grp", plots = FALSE, posthoc = FALSE)
  expect_gt(manip$effect_sizes$generalised_eta_sq[2], 0.2)

  expect_error(anova_rm(d, "y", "id", "time", between = "age_grp",
                        observed = "sex", plots = FALSE),
               "`observed` must name within- or between-subject factors")
})

# F019 ------------------------------------------------------------------------

test_that("with fewer subjects than levels the requested correction is applied (F019)", {
  skip_if_not_installed("afex")
  d <- withr::with_seed(17, {
    d <- expand.grid(id = factor(1:4), time = factor(paste0("t", 1:6)))
    d$y <- stats::rnorm(nrow(d)) + as.numeric(d$time) / 3
    d
  })
  gg <- anova_rm(d, "y", "id", "time", plots = FALSE)
  none <- anova_rm(d, "y", "id", "time", plots = FALSE, correction = "none")
  hf <- anova_rm(d, "y", "id", "time", plots = FALSE, correction = "HF")

  # pingouin.rm_anova(correction = True): eps = 0.343580, p-GG-corr = 0.037225
  # pingouin.epsilon(correction = "hf") = 0.759942
  eps_gg <- gg_by_hand(rm_wide(d))
  expect_equal(eps_gg, 0.3435801983, tolerance = 1e-8)
  expect_equal(gg$sphericity$gg_epsilon, eps_gg)
  expect_equal(gg$sphericity$hf_epsilon, 0.7599419288, tolerance = 1e-8)
  expect_true(is.na(gg$sphericity$mauchly_w))
  expect_equal(gg$anova$p_value, 0.037225, tolerance = 1e-4)
  expect_equal(gg$anova$p_value,
               stats::pf(gg$anova$statistic, 5 * eps_gg, 15 * eps_gg,
                         lower.tail = FALSE))
  expect_equal(gg$anova$num_df, 5 * eps_gg)
  expect_identical(attr(gg$anova, "correction"), "GG")
  expect_true(any(grepl("computed directly", gg$notes)))

  expect_equal(hf$anova$p_value,
               stats::pf(hf$anova$statistic, 5 * 0.7599419288, 15 * 0.7599419288,
                         lower.tail = FALSE), tolerance = 1e-7)
  # pingouin p-unc = 0.001599: the uncorrected p is labelled as such
  expect_equal(none$anova$p_value, 0.001599495, tolerance = 1e-6)
  expect_identical(attr(none$anova, "correction"), "none")
  expect_gt(gg$anova$p_value, 10 * none$anova$p_value)
})

test_that("a response in small units keeps its sphericity corrections (F019)", {
  skip_if_not_installed("afex")
  d <- rm_oneway(n = 12, k = 3, seed = 3, shift = 1 / 3)
  tiny <- d
  tiny$y <- tiny$y * 1e-10
  a <- anova_rm(d, "y", "id", "time", plots = FALSE)
  b <- anova_rm(tiny, "y", "id", "time", plots = FALSE)
  # pingouin.rm_anova(correction = True) on d: W = 0.8463, p-spher = 0.434132,
  # eps = 0.866776, p-GG-corr = 0.01343; epsilon(correction = "hf") = 1.0
  for (f in list(a, b)) {
    expect_equal(f$sphericity$mauchly_w, 0.8463, tolerance = 1e-4)
    expect_equal(f$sphericity$p_value, 0.434132, tolerance = 1e-5)
    expect_equal(f$sphericity$gg_epsilon, 0.866776, tolerance = 1e-5)
    expect_equal(f$sphericity$hf_epsilon, 1)
    expect_equal(f$anova$p_value, 0.01343, tolerance = 1e-3)
  }
  expect_equal(b$sphericity, a$sphericity)
  expect_equal(b$anova$p_value, a$anova$p_value)
  expect_equal(b$anova$mse / 1e-20, a$anova$mse)
  expect_false(any(grepl("could not", b$notes)))
})

# F021 ------------------------------------------------------------------------

test_that("arguments anova_rm() cannot honour are refused by name (F021)", {
  skip_if_not_installed("afex")
  d <- withr::with_seed(1, {
    d <- data.frame(id = rep(1:20, 3), time = rep(c("t1", "t2", "t3"), each = 20),
                    arm = rep(rep(c("p", "q"), 10), 3),
                    age = rep(stats::rnorm(20, 40, 5), 3))
    d$y <- stats::rnorm(60); d$w <- stats::runif(60)
    d
  })
  expect_error(anova_rm(d, "y", "id", "time", between = "arm",
                        interaction = FALSE, weights = "w", conf.level = 0.8,
                        vcov_type = "HC5", plots = FALSE),
               "no argument `interaction`, `weights`, `conf.level`, `vcov_type`")
  expect_error(anova_rm(d, "y", "id", "time", covariate = "age", plots = FALSE),
               "does not fit covariates \\(`covariate`\\)")
  expect_error(anova_rm(d, "y", "id", "time", id = "id", plots = FALSE),
               "`id` is set from `subject`")
  expect_error(anova_rm(d, "y", "id", "time", return = "nice", plots = FALSE),
               "`return` is set by anova_rm\\(\\) itself")
  # Every formal comes before `...`, so only a 14th positional argument lands there
  expect_error(anova_rm(d, "y", "id", "time", NULL, TRUE, NULL, "tukey", 0.95,
                        "GG", TRUE, FALSE, FALSE, median), "must be named")
  expect_error(anova_rm(d, "y", "id", "time", fun_aggregate = "median",
                        plots = FALSE), "must be a function")
  expect_error(anova_rm(d, "y", "id", "time", type = 1, plots = FALSE),
               "`type` must be 3")
  expect_error(anova_rm(d, "y", "id", "time", anova_table = list(es = "pes"),
                        plots = FALSE), "may hold only `p_adjust_method`")
  expect_error(anova_rm(d, "y", "id", "time", transformation = "sqrt",
                        plots = FALSE), "`transformation` is not supported")
  # A column in two roles
  expect_error(anova_rm(d, "y", "id", "time", between = "id", plots = FALSE),
               "given as both the subject identifier and a between-subject factor")
})

test_that("the whitelisted arguments are honoured (F021, F069)", {
  skip_if_not_installed("afex")
  d <- withr::with_seed(9, {
    d <- expand.grid(id = factor(1:14), time = factor(c("t1", "t2", "t3", "t4")))
    d$arm <- factor(ifelse(as.integer(d$id) <= 5, "ctl", "trt"))
    d$y <- stats::rnorm(14)[d$id] +
      as.numeric(d$time) * stats::rnorm(14, 1, 0.8)[d$id] +
      (d$arm == "trt") * as.numeric(d$time) * 0.3 + stats::rnorm(nrow(d))
    d
  })
  ref2 <- suppressMessages(afex::aov_ez("id", "y", d, within = "time",
                                        between = "arm", type = 2))
  t2 <- anova_rm(d, "y", "id", "time", between = "arm", type = 2, plots = FALSE)
  expect_equal(t2$anova$statistic, as.data.frame(ref2$anova_table)[["F"]])
  expect_identical(attr(t2$anova, "ss_type"), "II")
  expect_equal(attr(t2$model, "type"), 2)

  ref <- suppressMessages(afex::aov_ez("id", "y", d, within = "time",
                                       between = "arm", type = 3))
  adj <- anova_rm(d, "y", "id", "time", between = "arm", plots = FALSE,
                  anova_table = list(p_adjust_method = "bonferroni"))
  r <- as.data.frame(stats::anova(ref, correction = "GG",
                                  p_adjust_method = "bonferroni"))
  expect_equal(adj$anova$p_value, r[["Pr(>F)"]])
  expect_identical(attr(adj$anova, "p_adjust_method"), "bonferroni")
  expect_true(any(grepl("\\$sphericity are not adjusted", adj$notes)))
  plain <- anova_rm(d, "y", "id", "time", between = "arm", plots = FALSE)
  expect_equal(adj$sphericity$p_gg, plain$sphericity$p_gg)
})

# F022 / F068 -----------------------------------------------------------------

test_that("column names are never parsed: any name works as any role (F022)", {
  skip_if_not_installed("afex")
  d <- fx_repeated()
  base <- anova_rm(d, "score", "id", "time", between = "arm", plots = FALSE)
  odd <- c("participant id", "5id", "subject-id", "log(x)",
           "stop(\"evaluated\")", "a:b")
  for (nm in odd) {
    d2 <- d
    names(d2)[names(d2) == "id"] <- nm
    expect_silent(f <- anova_rm(d2, "score", nm, "time", between = "arm",
                                plots = FALSE))
    expect_equal(f$anova$statistic, base$anova$statistic, info = nm)
    expect_identical(names(f$data_used)[2], nm)
  }
  # ... and as the response, the within and the between factor
  d3 <- d
  names(d3)[match(c("score", "time", "arm"), names(d3))] <-
    c("my score", "time point", "arm (group)")
  f <- anova_rm(d3, "my score", "id", "time point", between = "arm (group)",
                emm_specs = c("time point", "arm (group)"))
  expect_identical(f$anova$term, c("arm (group)", "time point",
                                   "arm (group):time point"))
  expect_equal(f$anova$statistic, base$anova$statistic)
  expect_true(all(c("time point", "arm (group)") %in% names(f$emmeans)))
  expect_identical(f$sphericity$term, c("time point", "arm (group):time point"))
  expect_true(all(c("time point", "arm (group)") %in% names(f$plots$emmeans$data)))

  # User names that coincide with the internal ones are not confused with them
  d4 <- d
  names(d4)[match(c("score", "id", "time", "arm"), names(d4))] <-
    c("ID", "Y", "B1", "W1")
  f4 <- anova_rm(d4, "ID", "Y", "B1", between = "W1", plots = FALSE)
  expect_identical(f4$anova$term, c("W1", "B1", "W1:B1"))
  expect_equal(f4$anova$statistic, base$anova$statistic)
  expect_identical(f4$sphericity$term, c("B1", "W1:B1"))
  expect_true(any(grepl("averaged over the levels of: W1", f4$notes)))
})

test_that("level labels come back as given (F068)", {
  skip_if_not_installed("afex")
  d <- withr::with_seed(6, {
    d <- expand.grid(id = factor(1:10), week = c(0, 4, 12))
    d$y <- stats::rnorm(10)[d$id] + d$week / 4 + stats::rnorm(nrow(d))
    d
  })
  f <- anova_rm(d, "y", "id", within = "week", plots = FALSE)
  expect_identical(as.character(f$emmeans$week), c("0", "4", "12"))
  expect_identical(levels(f$data_used$week), c("0", "4", "12"))
  # emmeans prefixes numeric labels with the factor name, as it does for every
  # function in the package
  expect_identical(f$posthoc$contrast,
                   c("week0 - week4", "week0 - week12", "week4 - week12"))
  # One-way within design: the marginal means are the cell means
  expect_equal(f$emmeans$estimate, as.vector(tapply(d$y, d$week, mean)))

  # Labels whose pasted combinations collide in afex's wide format
  d3 <- withr::with_seed(2, {
    d3 <- expand.grid(id = factor(1:10), A = c("x_y", "x"), B = c("z", "y_z"),
                      stringsAsFactors = FALSE)
    d3$y <- stats::rnorm(nrow(d3))
    d3
  })
  f3 <- anova_rm(d3, "y", "id", within = c("A", "B"), plots = FALSE)
  recoded <- d3
  recoded$A <- factor(recoded$A, levels = c("x", "x_y"), labels = c("a1", "a2"))
  recoded$B <- factor(recoded$B, levels = c("y_z", "z"), labels = c("b1", "b2"))
  ref <- afex::aov_ez("id", "y", recoded, within = c("A", "B"))
  expect_equal(f3$anova$statistic, as.data.frame(ref$anova_table)[["F"]])
  expect_setequal(as.character(f3$emmeans$A), c("x", "x_y"))
  expect_setequal(as.character(f3$emmeans$B), c("z", "y_z"))
})

# F140 ------------------------------------------------------------------------

test_that("$residuals line up with $data_used row by row (F140)", {
  skip_if_not_installed("afex")
  d <- withr::with_seed(1, {
    d <- expand.grid(time = factor(c("t1", "t2", "t3")), id = factor(1:24))[, 2:1]
    d$score <- 10 + 2 * as.numeric(d$time) + stats::rnorm(nrow(d), 0, 1)
    d
  })
  d$score[d$id == 5 & d$time == "t1"] <- 30
  fit <- anova_rm(d, "score", subject = "id", within = "time")
  du <- fit$data_used
  expect_length(fit$residuals, nrow(du))
  # stats::lm reference: subject plus time, fitted to $data_used as it stands
  ref <- stats::residuals(stats::lm(score ~ id + time, data = du))
  expect_equal(fit$residuals, unname(ref))
  worst <- du[which.max(abs(fit$residuals)), ]
  expect_identical(as.character(worst$id), "5")
  expect_identical(as.character(worst$time), "t1")
  expect_equal(fit$plots$index$data$residuals, fit$residuals)

  # With a between-subject factor the subjects are nested in it
  g <- fx_repeated()
  fg <- anova_rm(g[sample(nrow(g)), ], "score", "id", "time", between = "arm",
                 plots = FALSE)
  refg <- stats::residuals(stats::lm(score ~ id + time * arm,
                                     data = fg$data_used))
  expect_equal(fg$residuals, unname(refg))
})

# F158 ------------------------------------------------------------------------

test_that("normality is tested per within-subject cell on within-subject residuals (F158)", {
  skip_if_not_installed("afex")
  # sleep: a two-level within design is the paired t-test, whose assumption is
  # the normality of the differences
  f <- anova_rm(sleep, "extra", "ID", "group")
  ref <- stats::shapiro.test(sleep$extra[11:20] - sleep$extra[1:10])
  expect_identical(f$assumptions$normality$cell, c("1", "2"))
  expect_equal(f$assumptions$normality$p_value, rep(ref$p.value, 2))
  expect_equal(f$assumptions$normality$statistic, rep(unname(ref$statistic), 2))

  # The Q-Q plot shows the same residuals, standardised within each cell
  z <- f$residuals / stats::ave(f$residuals, f$data_used$group,
                                FUN = stats::sd)
  expect_equal(sort(f$plots$qq$data$value), sort(z))

  # Multivariate normal data with very unequal time-point variances: the
  # per-cell test holds its level (the pooled test rejected every time)
  withr::local_seed(1)
  p <- replicate(40, {
    n <- 30; sds <- c(1, 3, 9); u <- stats::rnorm(n)
    y <- sapply(1:3, function(j) 10 + j + sds[j] * (0.5 * u + sqrt(0.75) * stats::rnorm(n)))
    dd <- data.frame(id = factor(rep(1:n, 3)), time = factor(rep(1:3, each = n)),
                     y = as.vector(y))
    anova_rm(dd, "y", "id", "time", plots = FALSE,
             posthoc = FALSE)$assumptions$normality$p_value
  })
  expect_lt(mean(p < 0.05), 0.12)
  expect_gt(mean(p < 0.05), 0)
})

# F020 ------------------------------------------------------------------------

test_that("an integer subject identifier is not mistaken for a covariate (F020)", {
  skip_if_not_installed("afex")
  d <- withr::with_seed(2, {
    d <- expand.grid(id = 1:24, time = c("t1", "t2", "t3"))
    d$y <- stats::rnorm(24)[d$id] + as.numeric(d$time) + stats::rnorm(nrow(d))
    d
  })
  fit <- anova_rm(d, "y", subject = "id", within = "time", plots = FALSE,
                  posthoc = FALSE)
  expect_false(any(grepl("covariate", fit$notes)))
  expect_identical(nlevels(fit$data_used$id), 24L)
})

# F066 ------------------------------------------------------------------------

test_that("the Huynh-Feldt epsilon is capped at 1 wherever it is reported (F066)", {
  skip_if_not_installed("afex")
  d <- rm_oneway(n = 12, k = 3, seed = 1)
  fit <- anova_rm(d, "y", "id", within = "time", plots = FALSE, posthoc = FALSE)
  ref <- afex::aov_ez("id", "y", d, within = "time")
  raw <- as.matrix(unclass(suppressWarnings(
    summary(ref$Anova, multivariate = FALSE))$pval.adjustments))["time", "HF eps"]
  expect_gt(raw, 1)                                  # 1.1756
  expect_equal(fit$sphericity$hf_epsilon_raw, unname(raw))
  expect_identical(fit$sphericity$hf_epsilon, 1)     # pingouin.epsilon: 1.0
  expect_equal(fit$sphericity$p_hf,
               stats::pf(fit$anova$statistic, 2, 22, lower.tail = FALSE))
  expect_true(any(grepl("capped at 1", fit$notes)))  # under the default GG
})

# F067 ------------------------------------------------------------------------

test_that("the aggregation note counts combinations after balancing and names the function (F067)", {
  skip_if_not_installed("afex")
  d <- withr::with_seed(3, {
    d <- expand.grid(id = factor(sprintf("s%02d", 1:10)),
                     time = factor(c("t1", "t2", "t3")))
    d$y <- stats::rnorm(10)[d$id] + as.numeric(d$time) + stats::rnorm(nrow(d))
    d
  })
  extra <- d[d$id == "s03" & d$time == "t2", ]
  extra2 <- extra
  extra$y <- extra$y + 5
  extra2$y <- extra2$y - 1
  d1 <- rbind(d, extra, extra2)
  f1 <- anova_rm(d1, "y", "id", within = "time", fun_aggregate = median,
                 plots = FALSE, posthoc = FALSE)
  note <- grep("more than one row", f1$notes, value = TRUE)
  expect_length(note, 1L)
  expect_match(note, "^1 subject-by-cell combination\\(s\\)")
  expect_match(note, "(2 extra row(s))", fixed = TRUE)
  expect_match(note, "`median`", fixed = TRUE)
  cell <- f1$data_used$id == "s03" & f1$data_used$time == "t2"
  expect_equal(f1$data_used$y[cell],
               stats::median(d1$y[d1$id == "s03" & d1$time == "t2"]))
  expect_identical(f1$n_removed_aggregated, 2L)
  expect_identical(nrow(f1$data_used) + f1$n_removed, nrow(d1))

  # The only duplicate belongs to a subject that is then dropped
  d2 <- d[!(d$id == "s05" & d$time == "t3"), ]
  d2 <- rbind(d2, d[d$id == "s05" & d$time == "t1", ])
  f2 <- anova_rm(d2, "y", "id", within = "time", plots = FALSE, posthoc = FALSE)
  expect_false(any(grepl("more than one row|aggregated", f2$notes)))
  expect_identical(f2$n_removed_aggregated, 0L)
  expect_identical(f2$subjects_dropped, "s05")
})

# F069 ------------------------------------------------------------------------

test_that("global afex options do not change the result (F069)", {
  skip_if_not_installed("afex")
  d <- withr::with_seed(9, {
    d <- expand.grid(id = factor(1:14), time = factor(c("t1", "t2", "t3", "t4")))
    d$arm <- factor(ifelse(as.integer(d$id) <= 5, "ctl", "trt"))
    d$y <- stats::rnorm(14)[d$id] +
      as.numeric(d$time) * stats::rnorm(14, 1, 0.8)[d$id] +
      (d$arm == "trt") * as.numeric(d$time) * 0.3 + stats::rnorm(nrow(d))
    d
  })
  a <- anova_rm(d, "y", "id", within = "time", between = "arm", plots = FALSE)
  old <- afex::afex_options()
  withr::defer(do.call(afex::afex_options, old))
  afex::afex_options(type = 2, emmeans_model = "univariate", include_aov = TRUE,
                     check_contrasts = FALSE, return_aov = "nice")
  b <- anova_rm(d, "y", "id", within = "time", between = "arm", plots = FALSE)
  do.call(afex::afex_options, old)

  expect_equal(b$anova, a$anova)
  expect_equal(b$posthoc, a$posthoc)
  expect_identical(attr(a$anova, "ss_type"), "III")
  expect_identical(attr(a$emmeans, "emmeans_model"), "multivariate")
  expect_identical(attr(a$posthoc, "emmeans_model"), "multivariate")
  expect_true(any(grepl("Type III", utils::capture.output(print(a)))))

  # References: afex's Type III table and emmeans' multivariate model
  ref <- suppressMessages(afex::aov_ez("id", "y", d, within = "time",
                                       between = "arm", type = 3))
  expect_equal(a$anova$statistic, as.data.frame(ref$anova_table)[["F"]])
  mv <- summary(emmeans::contrast(
    emmeans::emmeans(ref, ~ time, model = "multivariate"), "pairwise"))
  expect_equal(a$posthoc$se, mv$SE)
  expect_equal(a$posthoc$df, mv$df)
})

# F091 ------------------------------------------------------------------------

test_that("averaging over factors left out of emm_specs is noted (F091)", {
  skip_if_not_installed("afex")
  d <- withr::with_seed(7, {
    d <- expand.grid(id = factor(1:24), time = factor(c("t1", "t2", "t3")))
    d$arm <- factor(rep(rep(c("ctrl", "trt"), each = 12), 3))
    d$score <- 10 + ifelse(d$arm == "trt", 1, -1) * (as.numeric(d$time) - 2) * 3 +
      stats::rnorm(nrow(d), 0, 1)
    d
  })
  f <- anova_rm(d, "score", "id", within = "time", between = "arm", plots = FALSE)
  # emmeans' own annotation of the same grid
  ref <- suppressMessages(afex::aov_ez("id", "score", d, within = "time",
                                       between = "arm"))
  mesg <- attr(summary(emmeans::emmeans(ref, ~ time, model = "multivariate")), "mesg")
  expect_true(any(grepl("averaged over the levels of: arm", mesg)))
  note <- grep("averaged over the levels of: arm", f$notes, value = TRUE)
  expect_length(note, 1L)
  expect_match(note, "significant interaction")
  expect_match(note, "emm_specs = c(\"time\", \"arm\")", fixed = TRUE)

  full <- anova_rm(d, "score", "id", within = "time", between = "arm",
                   emm_specs = c("time", "arm"), plots = FALSE)
  expect_false(any(grepl("averaged over", full$notes)))
})

# F130 ------------------------------------------------------------------------

test_that("a design with no error df or no error variance says so (F130)", {
  skip_if_not_installed("afex")
  d <- withr::with_seed(19, {
    d <- expand.grid(id = factor(1:3), time = factor(c("t1", "t2", "t3")))
    d$arm <- factor(c("A", "B", "C"))[as.integer(d$id)]
    d$y <- stats::rnorm(9)
    d
  })
  f <- anova_rm(d, "y", "id", "time", between = "arm", plots = FALSE)
  expect_true(all(f$anova$den_df == 0))
  expect_true(any(grepl("No error degrees of freedom remain", f$notes)))

  d2 <- expand.grid(id = factor(1:6), time = factor(c("t1", "t2", "t3")))
  d2$y <- as.numeric(d2$id) + as.numeric(d2$time)       # y = subject + time
  f2 <- anova_rm(d2, "y", "id", "time", plots = FALSE)
  # stats::aov gives the same degenerate within-subject error stratum
  a <- summary(stats::aov(y ~ time + Error(id / time), data = d2))
  expect_equal(a[["Error: id:time"]][[1]]["Residuals", "Sum Sq"], 0,
               tolerance = 1e-12)
  expect_true(any(grepl("The error variance is zero for time", f2$notes)))
  expect_false(any(grepl("too few|this few subjects", f2$notes)))
})

# F162 ------------------------------------------------------------------------

test_that("dropped subjects are listed in a fixed order whatever the row order (F162)", {
  skip_if_not_installed("afex")
  d <- withr::with_seed(3, {
    d <- expand.grid(id = sprintf("s%02d", 1:30), time = c("t1", "t2", "t3"),
                     stringsAsFactors = FALSE)
    d$score <- stats::rnorm(nrow(d), 10 + as.integer(factor(d$time)))
    d
  })
  d <- d[!(d$id %in% sprintf("s%02d", seq(3, 25, by = 2)) & d$time == "t2"), ]
  p <- withr::with_seed(2, d[sample(nrow(d)), ])
  a <- anova_rm(d, "score", "id", "time", plots = FALSE)
  b <- anova_rm(p, "score", "id", "time", plots = FALSE)
  expected <- sort(sprintf("s%02d", seq(3, 25, by = 2)))
  expect_identical(a$subjects_dropped, expected)
  expect_identical(b$subjects_dropped, expected)
  expect_identical(grep("removed", a$notes, value = TRUE),
                   grep("removed", b$notes, value = TRUE))
  expect_equal(b$anova, a$anova)
})

# F023 ------------------------------------------------------------------------

test_that("the sphericity note reads correctly under correction = \"none\" (F023)", {
  skip_if_not_installed("afex")
  d <- withr::with_seed(42, {
    d <- expand.grid(id = factor(1:12), time = factor(c("t1", "t2", "t3", "t4")))
    d$y <- stats::rnorm(12)[d$id] +
      as.numeric(d$time) * stats::rnorm(12, 1, 1)[d$id] + stats::rnorm(nrow(d))
    d
  })
  f <- anova_rm(d, "y", "id", within = "time", correction = "none",
                plots = FALSE, posthoc = FALSE)
  ref <- afex::aov_ez("id", "y", d, within = "time")
  sp <- as.matrix(unclass(summary(ref$Anova, multivariate = FALSE)$sphericity.tests))
  expect_lt(sp["time", "p-value"], 0.05)             # car: sphericity rejected
  note <- grep("Mauchly", f$notes, value = TRUE)
  expect_identical(note, "Mauchly's test rejects sphericity for time; the reported table is uncorrected (correction = \"none\"), so its p-values for that term may be too small.")
  expect_true(any(grepl("no sphericity correction",
                        utils::capture.output(print(f)))))
  g <- anova_rm(d, "y", "id", within = "time", plots = FALSE, posthoc = FALSE)
  expect_true(any(grepl("uses the Greenhouse-Geisser correction", g$notes)))
})
