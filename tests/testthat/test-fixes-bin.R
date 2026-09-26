# Fixes to anova_bin() from the adversarial review. Every expectation is
# checked against something computed independently of the package: stats::glm
# fits, car, sandwich, statsmodels (numbers quoted), or a hand computation from
# cell proportions or the published formula.

# Local fixtures ---------------------------------------------------------------

fb_three <- function(seed = 105, n = 300) {
  withr::with_seed(seed, {
    g <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
    data.frame(g = g, y = stats::rbinom(n, 1, c(a = .2, b = .5, c = .8)[as.character(g)]))
  })
}

fb_aggregated <- function() {
  # 300 subjects, one row per group and outcome, the count as the weight
  data.frame(g = rep(c("a", "b", "c"), each = 2), y = rep(c(0, 1), 3),
             n = c(80, 20, 50, 50, 20, 80))
}

fb_sep <- function(seed = 13) {
  withr::with_seed(seed, {
    d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 20)))
    d$y <- c(stats::rbinom(20, 1, .3), stats::rbinom(20, 1, .6), rep(1, 20))
    d
  })
}

fb_sep_big <- function(seed = 106) {
  withr::with_seed(seed, data.frame(
    g = factor(c(rep("a", 100), rep("b", 100), rep("c", 12))),
    y = c(stats::rbinom(100, 1, 0.3), stats::rbinom(100, 1, 0.5), rep(1, 12))))
}

cell_odds <- function(y, g) {
  p <- tapply(as.numeric(y), g, mean)
  stats::setNames(as.vector(p / (1 - p)), names(p))
}

# F028 / F077 / F116: model statistics -----------------------------------------

test_that("model_stats use the prior weights: aggregated data match statsmodels", {
  fit <- anova_bin(fb_aggregated(), "y", "g", weights = "n", plots = FALSE)
  ms <- fit$model_stats
  # statsmodels GLM(Binomial, freq_weights = n): llf -169.39520276363208,
  # llnull -207.94415416798356, aic 344.79040552726417
  expect_equal(ms$mcfadden_r2, 0.18538127007509175, tolerance = 1e-8)
  expect_equal(ms$lr_vs_null$statistic, 77.09790280870294, tolerance = 1e-8)
  expect_identical(ms$lr_vs_null$df, 2L)
  expect_equal(ms$lr_vs_null$p_value,
               stats::pchisq(77.09790280870294, 2, lower.tail = FALSE),
               tolerance = 1e-6)
  expect_equal(ms$logLik, -169.39520276363208, tolerance = 1e-8)
  expect_equal(ms$AIC, 344.79040552726417, tolerance = 1e-8)
})

test_that("integer weights give the same statistics as the expanded rows", {
  d <- fb_three()
  d$w <- withr::with_seed(6, sample(1:3, nrow(d), replace = TRUE))
  fit <- anova_bin(d, "y", "g", weights = "w", plots = FALSE)
  long <- d[rep(seq_len(nrow(d)), d$w), ]
  m <- stats::glm(y ~ g, data = long, family = stats::binomial())
  m0 <- stats::glm(y ~ 1, data = long, family = stats::binomial())
  ll <- stats::logLik(m); ll0 <- stats::logLik(m0)
  expect_equal(fit$model_stats$mcfadden_r2, as.numeric(1 - ll / ll0))
  expect_equal(fit$model_stats$lr_vs_null$statistic, as.numeric(2 * (ll - ll0)))
  expect_equal(fit$model_stats$logLik, as.numeric(ll))
  expect_equal(fit$model_stats$AIC, stats::AIC(m))
  expect_gt(fit$model_stats$mcfadden_r2, 0)
})

test_that("unweighted model_stats match glm, logLik and anova(null, full)", {
  d <- fb_three()
  fit <- anova_bin(d, "y", "g", plots = FALSE)
  m <- stats::glm(y ~ g, data = d, family = stats::binomial())
  m0 <- stats::glm(y ~ 1, data = d, family = stats::binomial())
  lrt <- stats::anova(m0, m, test = "LRT")
  ms <- fit$model_stats
  expect_equal(ms$AIC, stats::AIC(m))
  expect_equal(ms$BIC, stats::BIC(m))
  expect_equal(ms$logLik, as.numeric(stats::logLik(m)))
  expect_equal(ms$mcfadden_r2, as.numeric(1 - stats::logLik(m) / stats::logLik(m0)))
  expect_equal(ms$lr_vs_null$statistic, lrt$Deviance[2])
  expect_equal(ms$lr_vs_null$df, lrt$Df[2])
  expect_equal(ms$lr_vs_null$p_value, lrt$`Pr(>Chi)`[2])
  expect_identical(ms$success_level, "1")
})

test_that("non-integer weights: deviance-based statistics, no rounded logLik", {
  d <- fb_three()
  d$w <- withr::with_seed(7, stats::runif(nrow(d), 0.1, 0.45))
  fit <- anova_bin(d, "y", "g", weights = "w", plots = FALSE)
  m <- suppressWarnings(stats::glm(y ~ g, data = d, family = stats::binomial(),
                                   weights = w))
  expect_equal(fit$model_stats$lr_vs_null$statistic, m$null.deviance - m$deviance)
  expect_equal(fit$model_stats$mcfadden_r2, 1 - m$deviance / m$null.deviance)
  # a one-factor model: model versus null is the test of the factor
  expect_equal(fit$model_stats$lr_vs_null$statistic,
               suppressWarnings(car::Anova(m, type = 2))$`LR Chisq`)
  expect_true(is.na(fit$model_stats$AIC))
  expect_true(is.na(fit$model_stats$BIC))
  expect_true(is.na(fit$model_stats$logLik))
  expect_true(any(grepl("not whole numbers", fit$notes)))
  expect_lt(fit$model_stats$mcfadden_r2, 1)
})

# F030: odds ratios against the reference level, whatever the coding ----------

test_that("odds ratios are level-versus-reference under Type III and sum contrasts", {
  d <- withr::with_seed(1, {
    dd <- data.frame(g = factor(sample(c("a", "b", "c"), 300, TRUE)))
    dd$y <- stats::rbinom(300, 1, c(a = .2, b = .5, c = .8)[as.character(dd$g)])
    dd
  })
  o <- cell_odds(d$y, d$g)

  f2 <- anova_bin(d, "y", "g", plots = FALSE)
  f3 <- anova_bin(d, "y", "g", type = "III", plots = FALSE)
  expect_equal(f3$effect_sizes$odds_ratio, unname(o[c("b", "c")] / o["a"]))
  expect_equal(f3$effect_sizes$odds_ratio, f2$effect_sizes$odds_ratio)
  expect_identical(f3$effect_sizes$comparison, c("b vs a", "c vs a"))
  expect_identical(f3$effect_sizes$factor, c("g", "g"))
  # profile intervals too: the same as for the treatment-coded glm
  ci <- suppressMessages(stats::confint(
    stats::glm(y ~ g, data = d, family = stats::binomial())))[-1, ]
  expect_equal(f3$effect_sizes$conf_low, unname(exp(ci[, 1])), tolerance = 1e-6)
  # ... while $anova is still the Type III table
  expect_identical(attr(f3$anova, "ss_type"), "III")

  f3c <- anova_bin(d, "y", "g", type = "III", reference = list(g = "c"),
                   plots = FALSE)
  expect_equal(f3c$effect_sizes$odds_ratio, unname(o[c("a", "b")] / o["c"]))
  expect_identical(f3c$effect_sizes$comparison, c("a vs c", "b vs c"))

  withr::local_options(contrasts = c("contr.sum", "contr.poly"))
  fs <- anova_bin(d, "y", "g", plots = FALSE)
  expect_equal(fs$effect_sizes$odds_ratio, unname(o[c("b", "c")] / o["a"]))
})

test_that("an ordered group gets level-versus-reference odds ratios, not trends", {
  d <- fb_three()
  d$g <- factor(d$g, levels = c("a", "b", "c"), ordered = TRUE)
  o <- cell_odds(d$y, d$g)
  fit <- anova_bin(d, "y", "g", plots = FALSE)
  expect_equal(fit$effect_sizes$odds_ratio, unname(o[c("b", "c")] / o["a"]))
  expect_identical(fit$effect_sizes$comparison, c("b vs a", "c vs a"))
  expect_false(any(grepl("\\.L|\\.Q", fit$effect_sizes$term)))
})

test_that("odds ratios in an interaction model say where they hold", {
  d <- withr::with_seed(2, {
    dd <- expand.grid(g = c("a", "b"), site = c("N", "S"), rep = 1:100)
    lp <- -1 + (dd$g == "b") + 0.5 * (dd$site == "S") +
      1.5 * (dd$g == "b") * (dd$site == "S")
    dd$y <- stats::rbinom(nrow(dd), 1, stats::plogis(lp))
    dd
  })
  o <- tapply(d$y, list(d$g, d$site), mean)
  o <- o / (1 - o)
  fit <- anova_bin(d, "y", c("g", "site"), interaction = TRUE, type = "III",
                   plots = FALSE)
  es <- fit$effect_sizes
  expect_identical(es$comparison,
                   c("b vs a at site = N", "S vs N at g = a", "(b vs a) x (S vs N)"))
  expect_equal(es$odds_ratio[1], o["b", "N"] / o["a", "N"])
  expect_equal(es$odds_ratio[2], o["a", "S"] / o["a", "N"])
  expect_equal(es$odds_ratio[3],
               (o["b", "S"] / o["a", "S"]) / (o["b", "N"] / o["a", "N"]))
})

# F078: reference on an ordered factor -------------------------------------------

test_that("a reference on an ordered factor drops the ordering instead of permuting it", {
  d <- withr::with_seed(105, {
    g <- factor(sample(c("low", "mid", "high"), 300, TRUE),
                levels = c("low", "mid", "high"), ordered = TRUE)
    data.frame(g = g, y = stats::rbinom(300, 1, c(low = .2, mid = .5, high = .8)[as.character(g)]))
  })
  o <- cell_odds(d$y, d$g)
  fit <- anova_bin(d, "y", "g", reference = list(g = "mid"), plots = FALSE)
  expect_false(is.ordered(fit$data_used$g))
  expect_identical(levels(fit$data_used$g), c("mid", "low", "high"))
  expect_equal(fit$effect_sizes$odds_ratio, unname(o[c("low", "high")] / o["mid"]))
  expect_identical(fit$effect_sizes$comparison, c("low vs mid", "high vs mid"))
  expect_true(any(grepl("ordered factor", fit$notes)))
  # the omnibus test is the same test of g
  ref <- car::Anova(stats::glm(y ~ g, data = d, family = stats::binomial()))
  expect_equal(fit$anova$statistic, ref$`LR Chisq`)
})

# F153: robust SEs at the boundary --------------------------------------------------

test_that("robust SEs under separation fall back to model-based ones", {
  skip_if_not_installed("sandwich")
  d <- withr::with_seed(3, {
    dd <- data.frame(g = factor(rep(c("a", "b", "c"), each = 20)))
    dd$y <- stats::rbinom(60, 1, c(a = .3, b = .5, c = .6)[as.character(dd$g)])
    dd$y[dd$g == "c"] <- 1
    dd
  })
  fit <- anova_bin(d, "y", "g", vcov_type = "HC3", plots = FALSE)
  m <- stats::glm(y ~ g, data = d, family = stats::binomial())
  se_model <- sqrt(diag(stats::vcov(m)))[-1]
  expect_equal(fit$effect_sizes$se_log_or, unname(se_model))
  gc <- fit$effect_sizes$term == "gc"
  expect_gt(fit$effect_sizes$p_value[gc], 0.9)     # as the separation note says
  expect_true(all(fit$posthoc$p_adjusted > 0.5))
  expect_true(any(grepl("boundary", fit$notes)))
  # the collapse the fallback avoids: the sandwich SE of gc is tiny
  expect_lt(sqrt(sandwich::vcovHC(m, type = "HC3")["gc", "gc"]), 1)
})

# F031: separation found whatever the coding and n ----------------------------------

test_that("separation is reported under Type III and names the odds ratio", {
  fit <- anova_bin(fb_sep_big(), "y", "g", type = "III", plots = FALSE)
  sep <- grep("separation", fit$notes, value = TRUE)
  expect_length(sep, 1L)
  expect_match(sep, "g = c", fixed = TRUE)
  expect_match(sep, "The affected coefficient(s): gc.", fixed = TRUE)

  two <- data.frame(g = factor(rep(c("a", "b"), c(50, 8))),
                    y = c(withr::with_seed(4, stats::rbinom(50, 1, 0.4)), rep(1, 8)))
  f2 <- anova_bin(two, "y", "g", type = "III", plots = FALSE)
  expect_match(grep("separation", f2$notes, value = TRUE),
               "The affected coefficient(s): gb.", fixed = TRUE)
})

test_that("separation is reported in a large Type II sample", {
  d <- withr::with_seed(1, data.frame(
    g = factor(c(rep("a", 20000), rep("b", 20000), rep("c", 5))),
    y = c(stats::rbinom(20000, 1, 0.9), stats::rbinom(20000, 1, 0.95), rep(1, 5))))
  fit <- anova_bin(d, "y", "g", posthoc = FALSE, plots = FALSE)
  m <- stats::glm(y ~ g, data = d, family = stats::binomial())
  # the textbook signature: a log-odds SE of about 88 from five observations
  expect_gt(stats::coef(summary(m))["gc", 2], 50)
  expect_true(any(grepl("separation.*g = c", fit$notes)))
  gc <- fit$effect_sizes$term == "gc"
  expect_identical(fit$effect_sizes$conf_high[gc], Inf)
})

# F116: the separation note, pinned exactly ------------------------------------------

test_that("the separation note names exactly the affected coefficients", {
  fit <- anova_bin(fb_sep_big(), "y", "g", plots = FALSE)
  sep <- grep("separation", fit$notes, value = TRUE)
  expect_length(sep, 1L)
  expect_match(sep, "The affected coefficient(s): gc. ", fixed = TRUE)

  # a separated BASELINE: its intercept runs away too, and is not listed
  rel <- anova_bin(fb_sep_big(), "y", "g", reference = list(g = "c"), plots = FALSE)
  sep <- grep("separation", rel$notes, value = TRUE)
  expect_match(sep, "The affected coefficient(s): ga, gb. ", fixed = TRUE)
  expect_false(grepl("Intercept", sep))
})

test_that("a rare but identified event rate is not separation", {
  d <- data.frame(g = factor(rep(c("a", "b"), each = 5000)),
                  y = 0)
  d$y[c(1:3, 5001:5010)] <- 1                    # 3 and 10 events in 5000
  fit <- anova_bin(d, "y", "g", posthoc = FALSE, plots = FALSE)
  m <- stats::glm(y ~ g, data = d, family = stats::binomial())
  expect_lt(min(stats::fitted(m)), 1e-3)          # small fitted probabilities
  expect_false(any(grepl("separation", fit$notes)))
  expect_true(all(is.finite(fit$effect_sizes$conf_high)))
})

# F045: profile intervals under separation ------------------------------------------

test_that("a separated odds ratio gets its true profile bound and an open end", {
  d <- fb_sep()
  fit <- anova_bin(d, "y", "g", plots = FALSE)
  es <- fit$effect_sizes
  gc <- es$term == "gc"

  # Hand computation from the definition: fix the log odds ratio at b, fit
  # the rest with b as an offset, keep b while the deviance rises by at most
  # qchisq(0.95, 1). The upper end never reaches the cut-off.
  m <- stats::glm(y ~ g, data = d, family = stats::binomial())
  gap <- function(b) {
    stats::deviance(stats::glm(y ~ I(g == "b"), data = d, family = stats::binomial(),
                               offset = b * (g == "c"))) -
      stats::deviance(m) - stats::qchisq(0.95, 1)
  }
  low <- exp(stats::uniroot(gap, c(0, 6), tol = 1e-10)$root)
  expect_equal(es$conf_low[gc], low, tolerance = 1e-5)
  expect_gt(es$conf_low[gc], 15)
  expect_identical(es$conf_high[gc], Inf)
  expect_identical(es$ci_method[gc], "profile")
  # the unaffected odds ratio keeps stats::confint()'s interval
  ci <- suppressWarnings(suppressMessages(stats::confint(m)))["gb", ]
  expect_equal(es$conf_low[es$term == "gb"], unname(exp(ci[1])))
  expect_true(any(grepl("open on that side", fit$notes)))
})

test_that("an odds ratio driven to zero is reported from 0 to its profile bound", {
  d <- fb_sep_big()
  fit <- anova_bin(d, "y", "g", reference = list(g = "c"), plots = FALSE)
  es <- fit$effect_sizes
  m <- stats::glm(y ~ g, data = transform(d, g = stats::relevel(g, "c")),
                  family = stats::binomial())
  gap <- function(b) {
    stats::deviance(stats::glm(y ~ I(g == "b"), data = d, family = stats::binomial(),
                               offset = b * (g == "a"))) -
      stats::deviance(m) - stats::qchisq(0.95, 1)
  }
  high <- exp(stats::uniroot(gap, c(-6, 0), tol = 1e-10)$root)
  expect_identical(es$conf_low[es$term == "ga"], 0)
  expect_equal(es$conf_high[es$term == "ga"], high, tolerance = 1e-5)
})

# F029 / F033 / F081: the proportion table --------------------------------------------

test_that("the proportion table uses the weights", {
  fit <- anova_bin(fb_aggregated(), "y", "g", weights = "n")
  p <- fit$assumptions$proportions
  ones <- p[p$y == "1", ]
  expect_equal(ones$n, c(20, 50, 80))
  expect_equal(ones$group_total, c(100, 100, 100))
  expect_equal(ones$proportion, c(0.2, 0.5, 0.8))
  expect_equal(ones$rows, c(1L, 1L, 1L))
  expect_equal(ones$group_rows, c(2L, 2L, 2L))
  # the same as the fitted probabilities of the saturated model
  expect_equal(ones$proportion, fit$emmeans$estimate, tolerance = 1e-6)
  labels <- ggplot2::ggplot_build(fit$plots$proportions)$data[[2]]$label
  expect_setequal(labels, c("20%", "50%", "80%"))

  d <- fb_three()
  d$w <- withr::with_seed(8, stats::runif(nrow(d), 0.5, 3))
  fw <- anova_bin(d, "y", "g", weights = "w", plots = FALSE)
  pw <- fw$assumptions$proportions
  pw1 <- pw[pw$y == "1", ]
  ref <- vapply(split(d, d$g), function(s) stats::weighted.mean(s$y, s$w), numeric(1))
  expect_equal(pw1$proportion, unname(ref))
  expect_equal(pw1$rows, as.integer(table(d$g, d$y)[, "1"]))
})

test_that("a response called n, proportion or group_total keeps its levels", {
  base <- fb_three()
  tab <- table(base$g, base$y)
  for (nm in c("n", "proportion", "group_total")) {
    d <- data.frame(g = base$g)
    d[[nm]] <- base$y
    fit <- anova_bin(d, nm, "g")
    p <- fit$assumptions$proportions
    expect_setequal(p[[nm]], c("0", "1"))
    stat_cols <- setdiff(names(p), c(".cell", nm))
    expect_length(stat_cols, 3L)
    cnt <- p[[stat_cols[1]]]
    prop <- p[[stat_cols[3]]]
    expect_equal(cnt[p[[nm]] == "1"], as.integer(tab[, "1"]))
    expect_equal(prop[p[[nm]] == "1"], unname(prop.table(tab, 1)[, "1"]))
    b <- ggplot2::ggplot_build(fit$plots$proportions)
    expect_setequal(b$data[[2]]$label,
                    sprintf("%.0f%%", 100 * as.vector(prop.table(tab, 1))))
    expect_identical(fit$plots$proportions$labels$fill, nm)
  }
})

test_that("an empty-string group gets its totals and proportions", {
  d <- withr::with_seed(105, {
    g <- sample(c("", "b", "c"), 300, TRUE)
    data.frame(g = g, y = stats::rbinom(300, 1, c(.2, .5, .8)[match(g, c("", "b", "c"))]))
  })
  p <- anova_bin(d, "y", "g", plots = FALSE)$assumptions$proportions
  tab <- table(d$g, d$y)
  blank <- which(rownames(tab) == "")         # "" cannot index by name
  e <- p[p$.cell == "", ]
  expect_equal(e$group_total, rep(as.integer(sum(tab[blank, ])), 2))
  expect_equal(e$proportion, as.vector(prop.table(tab, 1)[blank, ]))
  expect_false(anyNA(p$proportion))
})

# F079 / F035 / F080: the response and the arguments -------------------------------

test_that("an ordered binary response accepts a success level", {
  base <- fb_three()
  d <- data.frame(g = base$g,
                  resp = factor(ifelse(base$y == 1, "high", "low"),
                                levels = c("low", "high"), ordered = TRUE))
  fit <- anova_bin(d, "resp", "g", success = "low", plots = FALSE)
  expect_identical(fit$model_stats$success_level, "low")
  o <- cell_odds(d$resp == "low", d$g)
  expect_equal(fit$effect_sizes$odds_ratio, unname(o[c("b", "c")] / o["a"]))
  expect_equal(round(fit$effect_sizes$odds_ratio, 3), c(0.336, 0.110))
})

test_that("the success level of a character response does not depend on the locale", {
  d <- withr::with_seed(4, data.frame(g = rep(c("a", "b", "c"), each = 20),
                                      resp = sample(c("Yes", "no"), 60, TRUE)))
  # the C-locale (radix) order: upper case first, so "no" is the second level
  expected <- sort(c("Yes", "no"), method = "radix")[2]
  locs <- c("C", "C.UTF-8", "en_US.UTF-8")
  ok <- vapply(locs, function(l) {
    old <- Sys.getlocale("LC_COLLATE")
    on.exit(Sys.setlocale("LC_COLLATE", old))
    nzchar(suppressWarnings(Sys.setlocale("LC_COLLATE", l)))
  }, logical(1))
  locs <- locs[ok]
  differs <- vapply(locs, function(l) {
    withr::with_collate(l, sort(c("Yes", "no"))[2] != expected)
  }, logical(1))
  skip_if(!any(differs), "no installed locale collates differently from C")
  for (l in locs) {
    s <- withr::with_collate(l, anova_bin(d, "resp", "g", plots = FALSE)$model_stats$success_level)
    expect_identical(s, expected, info = l)
  }
})

test_that("success and reference must be single values, named once", {
  d <- fb_three()
  expect_error(anova_bin(d, "y", "g", success = c("0", "1")), "single")
  expect_error(anova_bin(d, "y", "g", success = NA), "single")
  expect_error(anova_bin(d, "y", "g", reference = list(g = c("b", "c"))), "single")
  expect_error(anova_bin(d, "y", "g", reference = list(g = character(0))), "single")
  expect_error(anova_bin(d, "y", "g", reference = list(g = "b", g = "c")),
               "at most once")
})

# F149: sparse data ------------------------------------------------------------------

test_that("sparse cells get a note on the liberal LR test", {
  mk <- function(y, n) data.frame(
    g = factor(rep(c("a", "b", "c"), each = n)),
    r = unlist(lapply(y, function(k) c(rep(1, k), rep(0, n - k)))))
  d <- mk(c(0, 3, 4), 20)
  fit <- anova_bin(d, "r", "g", plots = FALSE)
  note <- grep("Sparse data", fit$notes, value = TRUE)
  expect_length(note, 1L)
  expect_match(note, "liberal")
  expect_match(note, "Wald", fixed = TRUE)
  expect_match(note, "Rao", fixed = TRUE)
  # Expected counts under the null, as chisq.test() computes them
  ex <- suppressWarnings(stats::chisq.test(table(d$g, d$r)))$expected
  expect_lt(min(ex), 5)
  expect_match(note, sprintf("smallest %s", format(signif(min(ex), 3))), fixed = TRUE)

  # 10 expected events per group: no note
  ok <- anova_bin(mk(c(8, 10, 12), 50), "r", "g", plots = FALSE)
  ex_ok <- stats::chisq.test(table(mk(c(8, 10, 12), 50)$g, mk(c(8, 10, 12), 50)$r))$expected
  expect_gte(min(ex_ok), 5)
  expect_false(any(grepl("Sparse data", ok$notes)))
})

# F032: dispersion --------------------------------------------------------------------

test_that("no Pearson dispersion is reported for a 0/1 response", {
  d <- fb_three()
  fit <- anova_bin(d, "y", "g", plots = FALSE)
  expect_true(is.na(fit$assumptions$dispersion))
  expect_true(any(grepl("dispersion is NA", fit$notes)))
  # why: for a model saturated in the cells it is N / (N - k) whatever the data
  m <- stats::glm(y ~ g, data = d, family = stats::binomial())
  expect_equal(sum(stats::residuals(m, "pearson")^2) / stats::df.residual(m),
               nrow(d) / (nrow(d) - 3))
  agg <- anova_bin(fb_aggregated(), "y", "g", weights = "n", plots = FALSE)
  expect_true(is.na(agg$assumptions$dispersion))
})

# F006: roles ---------------------------------------------------------------------------

test_that("the response cannot also be a grouping variable", {
  d <- fb_three()
  d$bin <- d$y
  expect_error(anova_bin(d, "bin", c("g", "bin"), plots = FALSE),
               "both the response and a grouping variable")
})

test_that("nothing is printed", {
  expect_silent_stdout(anova_bin(fb_sep(), "y", "g", type = "III"))
  expect_silent(anova_bin(fb_aggregated(), "y", "g", weights = "n"))
})
