# Detailed regression tests: run locally and on CI (NOT_CRAN=true), skipped
# on CRAN to keep its check time short. The core suite runs everywhere.
skip_on_cran()

# Regression tests for the third review of anova_bin() and anova_count().
# Each is checked against an independent computation (base R, MASS, emmeans,
# sandwich or a hand computation), not against the package's own internals.

test_that("R3-03: anova_bin's BIC on a frequency table is the BIC of the expanded data", {
  agg <- data.frame(g = rep(c("a", "b", "c"), each = 2), y = rep(0:1, 3),
                    w = c(30, 20, 25, 25, 15, 35))
  ex <- agg[rep(1:6, agg$w), c("g", "y")]
  ref <- stats::BIC(stats::glm(y ~ g, stats::binomial(), ex))
  fa <- anova_bin(agg, "y", "g", weights = "w", plots = FALSE)
  expect_equal(fa$model_stats$BIC, ref, tolerance = 1e-10)
  expect_equal(fa$model_stats$BIC,
               anova_bin(ex, "y", "g", plots = FALSE)$model_stats$BIC,
               tolerance = 1e-10)
  expect_equal(fa$model_stats$BIC,
               anova_glm(agg, "y", "g", family = stats::binomial(), weights = "w",
                         plots = FALSE)$model_stats$BIC, tolerance = 1e-10)
  # AIC and the log-likelihood were already those of the expanded data
  expect_equal(fa$model_stats$AIC,
               stats::AIC(stats::glm(y ~ g, stats::binomial(), ex)), tolerance = 1e-10)
  # weights that are not frequencies (all 1) keep stats::BIC()
  ones <- transform(ex, w = 1)
  expect_equal(anova_bin(ones, "y", "g", weights = "w", plots = FALSE)$model_stats$BIC,
               ref, tolerance = 1e-10)
})

test_that("R3-09: a negative binomial fit at an enormous theta gets the Poisson-limit omnibus test", {
  skip_if_not_installed("MASS")
  d2 <- data.frame(y = c(3, 3, 3, 5, 5, 5, 8, 8, 8), g = rep(c("a", "b", "c"), each = 3))
  r <- anova_count(d2, "y", "g", model = "negbin", plots = FALSE)
  expect_identical(r$model_type, "negbin")
  expect_gt(r$model$theta, 1e6, label = paste(r$notes, collapse = " | "))
  ref <- stats::anova(stats::glm(y ~ g, stats::poisson(), d2), test = "Chisq")
  expect_equal(r$anova$statistic[1], ref$Deviance[2], tolerance = 1e-6)
  expect_equal(r$anova$p_value[1], ref$`Pr(>Chi)`[2], tolerance = 1e-6)
  # the same limit reached through a large but representable theta
  big <- stats::glm(y ~ g, MASS::negative.binomial(1e8), d2)
  big0 <- stats::glm(y ~ 1, MASS::negative.binomial(1e8), d2)
  expect_equal(r$anova$statistic[1], stats::deviance(big0) - stats::deviance(big),
               tolerance = 1e-5)
  expect_true(any(grepl("analysis of deviance is that of the Poisson fit", r$notes,
                        fixed = TRUE)))
  # the omnibus test and the Tukey comparisons now agree in direction
  expect_lt(r$anova$p_value[1], 0.05)
  expect_lt(min(r$posthoc$p_value), 0.05)

  # an ordinary theta keeps car's negative binomial test (theta held fixed)
  withr::local_seed(1)
  d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 20)))
  d$y <- stats::rnbinom(60, mu = rep(c(3, 4, 6), each = 20), size = 2)
  f <- anova_count(d, "y", "g", model = "negbin", plots = FALSE)
  nb0 <- stats::glm(y ~ 1, MASS::negative.binomial(f$model$theta), d)
  expect_equal(f$anova$statistic[1], stats::deviance(nb0) - stats::deviance(f$model),
               tolerance = 1e-6)
  expect_false(any(grepl("that of the Poisson fit", f$notes, fixed = TRUE)))
})

test_that("R3-12: a constant count response is a clear error, unless an offset makes rates differ", {
  d <- data.frame(y = 3, g = rep(c("a", "b", "c"), each = 3))
  for (m in c("negbin", "quasipoisson")) {
    expect_error(anova_count(d, "y", "g", model = m, plots = FALSE),
                 "takes the same value \\(3\\) in every one of its 9 rows")
  }
  # constant after rows with a missing group are dropped
  d2 <- rbind(d, data.frame(y = 7, g = NA))
  expect_error(anova_count(d2, "y", "g", model = "negbin", plots = FALSE),
               "takes the same value")
  # with an exposure the rates differ, and the Poisson fit is the reference
  d$t <- rep(c(1, 2, 4), each = 3)
  f <- anova_count(d, "y", "g", offset = "t", model = "poisson", plots = FALSE)
  ref <- stats::anova(stats::glm(y ~ g + offset(log(t)), stats::poisson(), d),
                      test = "Chisq")
  expect_equal(f$anova$statistic[1], ref$Deviance[2], tolerance = 1e-8)
})

test_that("bincount-C3/R3: a constant rate errors only where the model needs a dispersion", {
  skip_if_not_installed("MASS")
  d <- data.frame(y = 3, g = rep(c("a", "b", "c"), each = 3), e = 1)
  # a constant offset leaves the rate constant
  for (m in c("negbin", "quasipoisson")) {
    expect_error(anova_count(d, "y", "g", offset = "e", model = m, plots = FALSE),
                 "The rate `y` / `e` takes the same value \\(3\\)")
  }
  # so does a count proportional to its exposure
  d2 <- data.frame(g = rep(c("a", "b", "c"), each = 3), e = rep(1:3, 3))
  d2$y <- 2 * d2$e
  expect_error(anova_count(d2, "y", "g", offset = "e", model = "negbin", plots = FALSE),
               "The rate `y` / `e` takes the same value \\(2\\)")
  # the Poisson LR test is well defined (statistic 0, p = 1), as in anova_glm()
  # and base R
  ref <- stats::anova(stats::glm(y ~ g, stats::poisson(), d), test = "Chisq")
  for (m in c("auto", "poisson")) {
    f <- anova_count(d, "y", "g", model = m, plots = FALSE)
    expect_identical(f$model_type, "poisson")
    expect_equal(f$anova$statistic[1], ref$Deviance[2], tolerance = 1e-8)
    expect_equal(f$anova$p_value[1], 1)
  }
  expect_equal(anova_glm(d, "y", "g", family = "poisson", plots = FALSE)$anova$p_value[1],
               f$anova$p_value[1])
  fo <- anova_count(d2, "y", "g", offset = "e", model = "poisson", plots = FALSE)
  expect_equal(fo$anova$p_value[1], 1)
  # rates that do vary are fitted as before
  d2$y[1] <- 5
  f2 <- anova_count(d2, "y", "g", offset = "e", model = "quasipoisson", plots = FALSE)
  ref2 <- stats::anova(stats::glm(y ~ g + offset(log(e)), stats::quasipoisson(), d2),
                       test = "F")
  expect_equal(f2$anova$statistic[1], ref2$F[2], tolerance = 1e-5)
})

test_that("R3-13: with robust standard errors the overdispersion note does not call them too small", {
  skip_if_not_installed("sandwich")
  withr::local_seed(2)
  d <- data.frame(g = rep(c("a", "b", "c"), each = 12))
  d$cnt <- stats::rnbinom(36, mu = exp(1.5 + rep(c(0, .3, .6), each = 12)), size = 2)
  m <- stats::glm(cnt ~ g, stats::poisson(), d)
  V <- sandwich::vcovHC(m, type = "HC1")
  ref <- as.data.frame(summary(emmeans::emmeans(m, ~ g, vcov. = V, type = "response")))

  f1 <- anova_count(d, "cnt", "g", model = "poisson", vcov_type = "HC1", plots = FALSE)
  # the reported standard errors are the sandwich ones, not the model's
  expect_equal(f1$emmeans$se, ref$SE, tolerance = 1e-6)
  note <- grep("assumes a dispersion of 1", f1$notes, value = TRUE)
  expect_length(note, 1L)
  expect_false(grepl("too small, so its intervals are too narrow", note, fixed = TRUE))
  expect_match(note, "robust (HC1) covariance", fixed = TRUE)
  expect_match(note, "Only the omnibus likelihood-ratio test", fixed = TRUE)
  # the quoted size is the hand formula at the Pearson dispersion
  phi <- sum(stats::residuals(m, type = "pearson")^2) / stats::df.residual(m)
  expect_match(note, sprintf("roughly %.0f%%", 100 * stats::pchisq(
    stats::qchisq(0.95, 1) / phi, 1, lower.tail = FALSE)), fixed = TRUE)

  # a robust Wald omnibus test leaves nothing resting on the dispersion
  f2 <- anova_count(d, "cnt", "g", model = "poisson", vcov_type = "HC1",
                    test_statistic = "Wald", plots = FALSE)
  expect_equal(f2$anova$statistic[1],
               suppressMessages(car::Anova(m, test.statistic = "Wald", vcov. = V))$Chisq[1],
               tolerance = 1e-6)
  expect_true(any(grepl("None of the results reported here rests on that assumption",
                        f2$notes, fixed = TRUE)))
  expect_false(any(grepl("too small", f2$notes)))

  # anova_glm's note says the same
  g <- anova_glm(d, "cnt", "g", family = stats::poisson(), vcov_type = "HC1",
                 plots = FALSE)
  gn <- grep("^Pearson dispersion", g$notes, value = TRUE)
  expect_length(gn, 1L)
  expect_false(grepl("times too small", gn))
  expect_match(gn, "robust (HC1) covariance", fixed = TRUE)
  # without a robust covariance the model-based warning is unchanged
  f0 <- anova_count(d, "cnt", "g", model = "poisson", plots = FALSE)
  expect_true(any(grepl(sprintf("about %.2f times too small", sqrt(phi)), f0$notes,
                        fixed = TRUE)))
})

test_that("R3-14: the separation note describes the Wald omnibus test that $anova holds", {
  d <- data.frame(g = rep(c("a", "b", "c"), each = 12),
                  bin = c(rep(0:1, c(10, 2)), rep(0:1, c(8, 4)), rep(1, 12)),
                  cnt = c(rep(c(2, 3, 4), 4), rep(c(1, 3, 5), 4), rep(0, 12)))
  mb <- stats::glm(bin ~ g, stats::binomial(), d)
  mc <- stats::glm(cnt ~ g, stats::poisson(), d)

  b <- anova_bin(d, "bin", "g", test_statistic = "Wald", plots = FALSE)
  expect_equal(b$anova$statistic[1], car::Anova(mb, test.statistic = "Wald")$Chisq[1],
               tolerance = 1e-6)
  # the Wald omnibus collapses where the LR test does not
  expect_gt(b$anova$p_value[1], 0.5)
  expect_lt(stats::anova(mb, test = "Chisq")$`Pr(>Chi)`[2], 1e-4)
  nb <- grep("separation", b$notes, value = TRUE)
  expect_match(nb, "Wald omnibus test in $anova", fixed = TRUE)
  expect_match(nb, "test_statistic = \"LR\"", fixed = TRUE)

  k <- anova_count(d, "cnt", "g", model = "poisson", test_statistic = "Wald",
                   plots = FALSE)
  expect_equal(k$anova$statistic[1], car::Anova(mc, test.statistic = "Wald")$Chisq[1],
               tolerance = 1e-4)
  expect_gt(k$anova$p_value[1], 0.5)
  expect_lt(stats::anova(mc, test = "Chisq")$`Pr(>Chi)`[2], 1e-8)
  nk <- grep("No events were observed", k$notes, value = TRUE)
  expect_match(nk, "Wald omnibus test in $anova", fixed = TRUE)
  expect_false(grepl("remains usable", nk))

  g <- anova_glm(d, "cnt", "g", family = stats::poisson(), test_statistic = "Wald",
                 plots = FALSE)
  expect_match(grep("No events were observed", g$notes, value = TRUE),
               "Wald omnibus test in $anova", fixed = TRUE)

  # with the default likelihood-ratio test the note keeps its old advice
  kl <- anova_count(d, "cnt", "g", model = "poisson", plots = FALSE)
  expect_match(grep("No events were observed", kl$notes, value = TRUE),
               "likelihood-ratio omnibus test remains usable", fixed = TRUE)
  bl <- anova_bin(d, "bin", "g", plots = FALSE)
  expect_false(grepl("Wald omnibus", grep("separation", bl$notes, value = TRUE)))
})

test_that("R3-18: anova_bin and anova_count give the two-cell note anova_glm gives", {
  withr::local_seed(5)
  d <- data.frame(g = rep(c("a", "b"), each = 15))
  d$x <- stats::rnorm(30); d$y <- stats::rnorm(30); d$y2 <- stats::rnorm(30)
  d$bin <- stats::rbinom(30, 1, .5)
  d$cnt <- stats::rpois(30, 3)
  mb <- stats::glm(bin ~ g, stats::binomial(), d)
  b <- anova_bin(d, "bin", "g", plots = FALSE)
  # the omnibus LR p and the Wald comparison p differ
  expect_equal(b$anova$p_value[1], stats::anova(mb, test = "Chisq")$`Pr(>Chi)`[2],
               tolerance = 1e-8)
  expect_equal(b$posthoc$p_value[1], summary(mb)$coefficients["gb", 4],
               tolerance = 1e-6)
  expect_false(isTRUE(all.equal(b$anova$p_value[1], b$posthoc$p_value[1])))
  note <- grep("two cells", b$notes, value = TRUE)
  expect_length(note, 1L)
  expect_identical(note, grep("two cells", anova_glm(d, "bin", "g",
    family = stats::binomial(), plots = FALSE)$notes, value = TRUE))
  expect_match(note, "can differ from the omnibus likelihood-ratio test", fixed = TRUE)

  k <- anova_count(d, "cnt", "g", plots = FALSE)
  expect_match(grep("two cells", k$notes, value = TRUE), "can differ", fixed = TRUE)
  # with a Wald omnibus test the comparison is the same test: z^2 = chi-square
  kw <- anova_count(d, "cnt", "g", test_statistic = "Wald", plots = FALSE)
  mc <- stats::glm(cnt ~ g, stats::poisson(), d)
  expect_equal(kw$anova$statistic[1], summary(mc)$coefficients["gb", 3]^2,
               tolerance = 1e-6)
  expect_equal(kw$anova$p_value[1], kw$posthoc$p_value[1], tolerance = 1e-6)
  expect_match(grep("two cells", kw$notes, value = TRUE), "is the omnibus test",
               fixed = TRUE)
  # three groups: no two-cell note
  d3 <- rbind(d, transform(d[1:15, ], g = "c"))
  expect_false(any(grepl("two cells", anova_bin(d3, "bin", "g", plots = FALSE)$notes)))
})

test_that("R3-23/bincount-R2: profile intervals cost nothing per row and match aggregated data", {
  withr::local_seed(1)
  n <- 3e5
  d <- data.frame(g = factor(sample(c("a", "b", "c", "d"), n, TRUE)))
  d$b <- stats::rbinom(n, 1, .4)
  el <- system.time(f <- anova_bin(d, "b", "g", posthoc = FALSE, plots = FALSE))[3]
  expect_identical(unique(f$effect_sizes$ci_method), "profile")
  expect_lt(el, 30)
  # reference: stats::confint() on the same data as a frequency table, whose
  # likelihood is the same
  agg <- as.data.frame(table(g = d$g, b = d$b))
  ma <- stats::glm(b ~ g, stats::binomial(), agg, weights = Freq)
  ci <- exp(suppressMessages(stats::confint(ma)))[-1, ]
  expect_equal(f$effect_sizes$conf_low, unname(ci[, 1]), tolerance = 1e-5)
  expect_equal(f$effect_sizes$conf_high, unname(ci[, 2]), tolerance = 1e-5)

  # rare events: individual rows and the frequency table give the same profile
  # intervals, equal to stats::confint() on the individual rows
  withr::local_seed(3)
  r <- data.frame(g = factor(rep(c("a", "b", "c"), each = 25000)))
  r$y <- stats::rbinom(75000, 1, rep(c(.0004, .0012, .002), each = 25000))
  fr <- anova_bin(r, "y", "g", posthoc = FALSE, plots = FALSE)
  ra <- as.data.frame(table(g = r$g, y = r$y))
  fa <- anova_bin(ra, "y", "g", weights = "Freq", posthoc = FALSE, plots = FALSE)
  ref <- exp(suppressMessages(stats::confint(stats::glm(y ~ g, stats::binomial(), r))))[-1, ]
  expect_identical(unique(fr$effect_sizes$ci_method), "profile")
  expect_identical(unique(fa$effect_sizes$ci_method), "profile")
  expect_equal(fr$effect_sizes$conf_high, unname(ref[, 2]), tolerance = 1e-5)
  expect_equal(fr$effect_sizes$conf_low, unname(ref[, 1]), tolerance = 1e-5)
  expect_equal(fr$effect_sizes$conf_high, fa$effect_sizes$conf_high, tolerance = 1e-6)
  # here Wald and profile intervals differ by more than 2% at the upper end
  wald <- exp(stats::confint.default(stats::glm(y ~ g, stats::binomial(), r)))[-1, 2]
  expect_gt(max(abs(fr$effect_sizes$conf_high / unname(wald) - 1)), 0.02)
  expect_false(any(grepl("so the odds-ratio intervals are Wald", fr$notes, fixed = TRUE)))

  # many coefficients: the same fallback as anova_glm(), noted up front
  withr::local_seed(1)
  k <- 60
  dk <- data.frame(g = factor(rep(sprintf("L%02d", 1:k), each = 8)))
  dk$b <- stats::rbinom(nrow(dk), 1, .5)
  dk$b[seq(1, nrow(dk), by = 8)] <- 0
  dk$b[seq(2, nrow(dk), by = 8)] <- 1
  fk <- anova_bin(dk, "b", "g", posthoc = FALSE, plots = FALSE)
  mk <- stats::glm(b ~ g, stats::binomial(), dk)
  expect_identical(unique(fk$effect_sizes$ci_method), "wald")
  expect_equal(fk$effect_sizes$conf_high, unname(exp(stats::confint.default(mk)[-1, 2])),
               tolerance = 1e-8)
  expect_true(any(grepl("so the odds-ratio intervals are Wald intervals", fk$notes,
                        fixed = TRUE)))
  expect_true(any(grepl("so the coefficient intervals are Wald intervals",
    anova_glm(dk, "b", "g", family = stats::binomial(), posthoc = FALSE,
              plots = FALSE)$notes, fixed = TRUE)))

  # asked for by name, profiling is honoured
  small <- d[1:2000, ]
  ms <- stats::glm(b ~ g, stats::binomial(), small)
  fp <- anova_bin(small, "b", "g", ci_method = "profile", posthoc = FALSE, plots = FALSE)
  expect_identical(unique(fp$effect_sizes$ci_method), "profile")
  expect_equal(fp$effect_sizes$conf_low,
               unname(exp(suppressMessages(stats::confint(ms))[-1, 1])), tolerance = 1e-4)
  # a small problem keeps profile intervals by default
  expect_identical(unique(anova_bin(small, "b", "g", posthoc = FALSE,
                                    plots = FALSE)$effect_sizes$ci_method), "profile")
})

test_that("bincount-C1: one separated level among many profiles only that odds ratio", {
  withr::local_seed(1)
  k <- 200
  d <- data.frame(g = factor(rep(sprintf("L%03d", 1:k), each = 10)))
  d$b <- stats::rbinom(nrow(d), 1, .4)
  m <- stats::glm(b ~ g, stats::binomial(), d)
  se <- summary(m)$coefficients[, 2]
  sep <- setdiff(names(se)[se > 10], "(Intercept)")
  expect_identical(sep, "gL090")      # one level is all 0 or all 1
  el <- system.time(f <- anova_bin(d, "b", "g", posthoc = FALSE, plots = FALSE))[3]
  expect_lt(el, 30)
  es <- f$effect_sizes
  expect_identical(es$term[es$ci_method == "profile"], sep)
  expect_identical(sum(es$ci_method == "wald"), as.integer(k - 2))
  # the others are the Wald intervals of glm()
  w <- es$ci_method == "wald"
  expect_equal(es$conf_low[w], unname(exp(stats::confint.default(m)[es$term[w], 1])),
               tolerance = 1e-8)
  # the separated one is open on its diverged side, and its finite end is
  # where the profile deviance (by refitting glm() with the coefficient as an
  # offset) rises by qchisq(0.95, 1)
  j <- which(es$term == sep)
  expect_identical(es$conf_low[j], 0)
  bb <- log(es$conf_high[j])
  # Every other level has its own free coefficient, so with gL090 fixed at bb
  # the profile deviance is that of the L001 and L090 cells sharing an
  # intercept (L090 offset by bb) plus that of the rest, each level on its own.
  two <- d[d$g %in% c("L001", "L090"), ]
  two$off <- bb * (two$g == "L090")
  rest <- droplevels(d[!d$g %in% c("L001", "L090"), ])
  dev_b <- stats::glm(b ~ 1, stats::binomial(), two, offset = off)$deviance +
    stats::glm(b ~ g, stats::binomial(), rest)$deviance
  expect_equal(dev_b - m$deviance, stats::qchisq(0.95, 1), tolerance = 1e-4)
  note <- grep("has 200 coefficients", f$notes, value = TRUE)
  expect_length(note, 1L)
  expect_match(note, "except for any odds ratio that separation drives", fixed = TRUE)
  expect_false(any(grepl("Profile-likelihood intervals were computed for", f$notes)))
})

test_that("R3-42: the sparse-data check does not depend on the scale of non-frequency weights", {
  withr::local_seed(11)
  d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 100)))
  d$y <- stats::rbinom(300, 1, rep(c(.4, .5, .6), each = 100))
  d$w <- 0.01
  d$w2 <- stats::runif(300, 0.5, 1.5)
  d$w3 <- d$w2 / 300
  for (wc in c("w", "w2", "w3")) {
    f <- anova_bin(d, "y", "g", weights = wc, plots = FALSE)
    expect_false(any(grepl("Sparse data", f$notes)), info = wc)
  }
  # hand computation: expected events under the null, with the weights
  # rescaled to average 1, in small groups of 6 rows
  s <- d[c(1:6, 101:106, 201:206), ]
  s$w <- s$w2 / 50
  ws <- s$w / mean(s$w)
  pbar <- sum(ws * s$y) / sum(ws)
  tot <- tapply(ws, droplevels(s$g), sum)
  smallest <- min(pmin(tot * pbar, tot * (1 - pbar)))
  note <- grep("Sparse data", anova_bin(s, "y", "g", weights = "w", plots = FALSE)$notes,
               value = TRUE)
  expect_match(note, sprintf("smallest %s", format(signif(smallest, 3))), fixed = TRUE)

  # frequency weights are counts and are used as they are: the same note as
  # the expanded data
  agg <- data.frame(g = rep(c("a", "b"), each = 2), y = rep(0:1, 2), n = c(9, 2, 6, 3))
  ex <- agg[rep(1:4, agg$n), c("g", "y")]
  expect_identical(
    grep("Sparse data", anova_bin(agg, "y", "g", weights = "n", plots = FALSE)$notes,
         value = TRUE),
    grep("Sparse data", anova_bin(ex, "y", "g", plots = FALSE)$notes, value = TRUE))
})

test_that("R3-43: a theta that ran away on overdispersed data is not called a collapse to Poisson", {
  skip_if_not_installed("MASS")
  d <- data.frame(y = c(0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 1, 0),
                  g = rep(c("a", "b"), each = 6))
  f <- anova_count(d, "y", "g", model = "negbin", plots = FALSE)
  expect_gt(f$model$theta, 1000)
  expect_false(is.null(f$model$th.warn))   # theta stopped at a limit
  # reference: the likelihood is far higher at a small theta than where
  # glm.nb() stopped, so the Poisson limit is not the maximum
  ll_small <- stats::logLik(stats::glm(y ~ g, MASS::negative.binomial(0.1), d,
                                       control = stats::glm.control(maxit = 100)))
  expect_gt(as.numeric(ll_small), as.numeric(stats::logLik(f$model)) + 10)
  expect_false(any(grepl("collapsed towards it", f$notes)))
  expect_true(any(grepl("did not converge", f$notes)))
  # Poisson Pearson dispersion, by hand
  pm <- stats::glm(y ~ g, stats::poisson(), d)
  expect_gt(sum(stats::residuals(pm, type = "pearson")^2) / stats::df.residual(pm), 5)

  # equidispersed data still collapse, as the likelihood says they should
  withr::local_seed(1)
  e <- data.frame(g = factor(rep(c("a", "b", "c"), each = 30)))
  e$y <- stats::rpois(90, rep(c(3, 4, 5), each = 30))
  fe <- anova_count(e, "y", "g", model = "negbin", plots = FALSE)
  ll_pois <- stats::logLik(stats::glm(y ~ g, stats::poisson(), e))
  ll_nb1 <- stats::logLik(stats::glm(y ~ g, MASS::negative.binomial(1), e))
  expect_gt(as.numeric(ll_pois), as.numeric(ll_nb1))
  expect_true(any(grepl("collapsed towards it", fe$notes)))
})

test_that("R3-51/bincount-C4: a kept Poisson fit says what the dispersion it cannot rule out would cost", {
  withr::local_seed(8)
  d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 10)),
                  y = stats::rnbinom(30, mu = 4, size = 10))
  m <- stats::glm(y ~ g, stats::poisson(), d)
  df <- stats::df.residual(m)
  x2 <- sum(stats::residuals(m, type = "pearson")^2)
  p_over <- stats::pchisq(x2, df, lower.tail = FALSE)
  expect_gt(p_over, 0.05)
  f <- anova_count(d, "y", "g", plots = FALSE)
  expect_identical(f$model_type, "poisson")
  note <- grep("Poisson model was kept", f$notes, value = TRUE)
  expect_match(note, "cannot tell the dispersion from 1", fixed = TRUE)
  # the upper 95% confidence bound of the dispersion, and the sizes there, by hand
  hi <- x2 / stats::qchisq(0.05, df)
  expect_match(note, sprintf("dispersion of up to %.2f", hi), fixed = TRUE)
  size2 <- stats::pchisq(stats::qchisq(0.95, 2) / hi, 2, lower.tail = FALSE)
  size1 <- stats::pchisq(stats::qchisq(0.95, 1) / hi, 1, lower.tail = FALSE)
  expect_match(note, sprintf("one on 2 (the largest term here) roughly %.0f%%", 100 * size2),
               fixed = TRUE)
  expect_match(note, sprintf("1 degree of freedom would reject roughly %.0f%%", 100 * size1),
               fixed = TRUE)

  # a significant dispersion quotes the 2-df size at the estimate
  withr::local_seed(3)
  d3 <- data.frame(g = factor(rep(c("a", "b", "c"), each = 100)),
                   y = stats::rnbinom(300, mu = 4, size = 10))
  m3 <- stats::glm(y ~ g, stats::poisson(), d3)
  phi3 <- sum(stats::residuals(m3, type = "pearson")^2) / stats::df.residual(m3)
  n3 <- grep("Poisson model was kept", anova_count(d3, "y", "g", plots = FALSE)$notes,
             value = TRUE)
  expect_match(n3, sprintf("roughly %.0f%%", 100 * stats::pchisq(
    stats::qchisq(0.95, 2) / phi3, 2, lower.tail = FALSE)), fixed = TRUE)
  expect_false(grepl("cannot tell the dispersion from 1", n3, fixed = TRUE))

  # with many degrees of freedom the bound is close to the estimate: no caveat
  # when a nominal 5% test would reject under 7.5% at the bound
  withr::local_seed(4)
  d4 <- data.frame(g = factor(rep(c("a", "b", "c"), each = 1000)))
  d4$y <- stats::rpois(3000, 5)
  m4 <- stats::glm(y ~ g, stats::poisson(), d4)
  hi4 <- sum(stats::residuals(m4, type = "pearson")^2) / stats::qchisq(0.05, 2997)
  expect_lt(stats::pchisq(stats::qchisq(0.95, 2) / hi4, 2, lower.tail = FALSE), 0.075)
  f4 <- anova_count(d4, "y", "g", plots = FALSE, posthoc = FALSE)
  expect_false(any(grepl("inflated|cannot tell", f4$notes)))
})

test_that("bincount-C4: under mild overdispersion nearly every kept Poisson fit is caveated", {
  skip_on_cran()
  # 10 per group, true dispersion 1.4 (negative binomial, mu 4, size 10)
  withr::local_seed(7)
  g <- factor(rep(c("a", "b", "c"), each = 10))
  M <- t(replicate(200, {
    d <- data.frame(g = g, y = stats::rnbinom(30, mu = 4, size = 10))
    f <- anova_count(d, "y", "g", posthoc = FALSE, plots = FALSE)
    c(pois = f$model_type == "poisson", rej = f$anova$p_value[1] < .05,
      cav = any(grepl("inflated", f$notes)))
  }))
  kept <- M[, "pois"] == 1
  # nearly every kept Poisson fit carries the caveat (it was about half before)
  expect_gt(mean(M[kept, "cav"]), 0.9)
  # the uncaveated ones are too few to matter
  expect_lt(sum(kept & !M[, "cav"]), 20)
})

test_that("bincount-C2: with a robust covariance the F-test note does not claim a dispersion of 1", {
  skip_if_not_installed("sandwich")
  withr::local_seed(2)
  d <- data.frame(g = c("a", "b", "c")[rep(1:3, each = 12)])
  d$cnt <- stats::rnbinom(36, mu = exp(1.5 + rep(c(0, .3, .6), each = 12)), size = 2)
  m <- stats::glm(cnt ~ g, stats::poisson(), d)
  V <- sandwich::vcovHC(m, type = "HC1")
  ref <- as.data.frame(summary(emmeans::emmeans(m, ~ g, vcov. = V, type = "response")))
  fc <- anova_count(d, "cnt", "g", model = "poisson", vcov_type = "HC1",
                    test_statistic = "F", plots = FALSE)
  fg <- anova_glm(d, "cnt", "g", family = stats::poisson(), vcov_type = "HC1",
                  test_statistic = "F", plots = FALSE)
  # what is reported does use the sandwich covariance
  expect_equal(fc$emmeans$se, ref$SE, tolerance = 1e-6)
  expect_equal(fg$emmeans$se, ref$SE, tolerance = 1e-6)
  se_gc <- sqrt(V["gc", "gc"])
  expect_equal(log(fc$effect_sizes$conf_high[2]),
               stats::coef(m)[["gc"]] + stats::qnorm(0.975) * se_gc, tolerance = 1e-6)
  for (n in list(fc$notes, fg$notes)) {
    fn <- grep("test_statistic = \"F\"", n, value = TRUE)
    expect_length(fn, 1L)
    expect_false(grepl("still assume a dispersion of 1", fn, fixed = TRUE))
    expect_match(fn, "robust (HC1) covariance", fixed = TRUE)
  }
  # without a robust covariance the note still says they assume 1
  f0 <- anova_count(d, "cnt", "g", model = "poisson", test_statistic = "F", plots = FALSE)
  expect_true(any(grepl("still assume a dispersion of 1", f0$notes, fixed = TRUE)))
})

test_that("bincount-R1: a near-Poisson negative binomial fit is a collapse, not a runaway", {
  skip_if_not_installed("MASS")
  d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 15)),
                  y = c(2, 3, 1, 0, 3, 0, 2, 2, 0, 1, 2, 2, 1, 1, 2, 3, 0, 2, 1, 3, 3, 4,
                        2, 1, 1, 4, 1, 0, 3, 2, 2, 4, 0, 4, 6, 1, 2, 5, 2, 2, 6, 3, 6, 0, 5))
  f <- anova_count(d, "y", "g", model = "negbin", plots = FALSE)
  # reference: the Pearson dispersion is below 1 and the likelihood barely
  # moves between the Poisson limit and a large theta
  p <- stats::glm(y ~ g, stats::poisson(), d)
  expect_lt(sum(stats::residuals(p, type = "pearson")^2) / stats::df.residual(p), 1)
  gain <- as.numeric(stats::logLik(stats::glm(y ~ g, MASS::negative.binomial(1e3), d))) -
    as.numeric(stats::logLik(p))
  expect_lt(abs(gain), 1e-3)
  expect_true(any(grepl("collapsed towards it", f$notes)))
  expect_false(any(grepl("did not converge|overdispersed \\(the variance", f$notes)))
})

test_that("bincount-R5: the two-cell note names a t test for an estimated dispersion", {
  mt <- mtcars
  mt$am <- factor(mt$am)
  f <- anova_count(mt, "carb", "am", model = "quasipoisson", test_statistic = "Wald",
                   plots = FALSE)
  m <- stats::glm(carb ~ am, stats::quasipoisson(), mt)
  # the comparison is the t test of summary.glm()
  expect_equal(abs(f$posthoc$statistic), abs(summary(m)$coefficients["am1", 3]),
               tolerance = 1e-6)
  expect_equal(f$posthoc$p_value, summary(m)$coefficients["am1", 4], tolerance = 1e-6)
  note <- grep("two cells", f$notes, value = TRUE)
  expect_match(note, "It is a t test with the estimated dispersion and can differ from the omnibus Wald chi-square test.",
               fixed = TRUE)
  # a fixed-dispersion family keeps the z test wording
  p <- grep("two cells", anova_count(mt, "carb", "am", model = "poisson",
                                     plots = FALSE)$notes, value = TRUE)
  expect_match(p, "It is a Wald z test with the model covariance and can differ from the omnibus likelihood-ratio test.",
               fixed = TRUE)
})
