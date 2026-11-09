# Detailed regression tests: run locally and on CI (NOT_CRAN=true), skipped
# on CRAN to keep its check time short. The core suite runs everywhere.
skip_on_cran()

# Regression tests for the second review of anova_bin() and anova_count().
# Each is checked against an independent computation (base R, emmeans, a hand
# computation), not against the package's own internals.

test_that("R2-BINCOUNT-01: colliding coefficient names do not corrupt effect sizes or emmeans", {
  withr::local_seed(14)
  n <- 400
  d <- data.frame(a = sample(c("b0", "b1"), n, TRUE),
                  ab = sample(c("0", "1"), n, TRUE))
  d$y <- stats::rbinom(n, 1, stats::plogis(-0.5 + 1.5 * (d$a == "b1") -
                                             1 * (d$ab == "1")))
  d$c <- stats::rpois(n, exp(1 + 0.7 * (d$a == "b1") - 0.5 * (d$ab == "1")))
  # the collision itself, in base R
  expect_true(anyDuplicated(names(stats::coef(
    stats::glm(y ~ a + ab, stats::binomial(), d)))) > 0L)

  # reference: the same data with a column name that cannot collide
  d2 <- d
  names(d2)[names(d2) == "ab"] <- "zz"
  d2$a <- factor(d2$a); d2$zz <- factor(d2$zz)
  rb <- stats::glm(y ~ a + zz, stats::binomial(), d2)
  rc <- stats::glm(c ~ a + zz, stats::poisson(), d2)

  for (ty in c("II", "III")) {
    f <- anova_bin(d, "y", c("a", "ab"), type = ty, plots = FALSE)
    expect_equal(nrow(f$effect_sizes), 2L, info = ty)
    expect_equal(f$effect_sizes$odds_ratio, unname(exp(stats::coef(rb)[-1])),
                 tolerance = 1e-6, info = ty)
    expect_identical(f$effect_sizes$ci_method, rep("profile", 2L), info = ty)
    expect_equal(f$effect_sizes$conf_low,
                 unname(exp(suppressMessages(stats::confint(rb))[-1, 1])),
                 tolerance = 1e-4, info = ty)
    em_a <- as.data.frame(emmeans::emmeans(rb, "a", type = "response"))
    em_z <- as.data.frame(emmeans::emmeans(rb, "zz", type = "response"))
    expect_equal(f$emmeans$estimate, c(em_a$prob, em_z$prob), tolerance = 1e-6,
                 info = ty)
    expect_false(any(grepl("Profile-likelihood intervals failed", f$notes)))
    expect_true(any(grepl("brackets", f$notes)), info = ty)

    fc <- anova_count(d, "c", c("a", "ab"), model = "poisson", type = ty,
                      plots = FALSE)
    expect_equal(fc$effect_sizes$IRR, unname(exp(stats::coef(rc)[-1])),
                 tolerance = 1e-6, info = ty)
    se <- sqrt(diag(stats::vcov(rc)))[-1]
    expect_equal(fc$effect_sizes$conf_high,
                 unname(exp(stats::coef(rc)[-1] + stats::qnorm(0.975) * se)),
                 tolerance = 1e-6, info = ty)
    em_a <- as.data.frame(emmeans::emmeans(rc, "a", type = "response"))
    em_z <- as.data.frame(emmeans::emmeans(rc, "zz", type = "response"))
    expect_equal(fc$emmeans$estimate, c(em_a$rate, em_z$rate),
                 tolerance = 1e-6, info = ty)
    expect_true(any(grepl("brackets", fc$notes)), info = ty)
    # the tests are those of the renamed data
    expect_equal(fc$anova$statistic,
                 anova_count(d2, "c", c("a", "zz"), model = "poisson",
                             type = ty, plots = FALSE)$anova$statistic,
                 info = ty)
  }

  # with the interaction, against the renamed fit
  fi <- anova_count(d, "c", c("a", "ab"), model = "poisson", interaction = TRUE,
                    plots = FALSE)
  ri <- stats::glm(c ~ a * zz, stats::poisson(), d2)
  expect_equal(fi$effect_sizes$IRR, unname(exp(stats::coef(ri)[-1])),
               tolerance = 1e-6)
})

test_that("R2-BINCOUNT-02: weights that are not frequencies get a scale note and a scale-free dispersion", {
  withr::local_seed(11)
  d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 100)))
  d$y <- stats::rbinom(300, 1, rep(c(.4, .5, .6), each = 100))
  d$c <- stats::rpois(300, rep(c(3, 4, 5), each = 100))
  d$w <- 0.01
  d$w2 <- 100.5

  # anova_bin uses them as glm() does (documented), and says that it does
  f <- anova_bin(d, "y", "g", weights = "w", plots = FALSE)
  m <- suppressWarnings(stats::glm(y ~ g, stats::binomial(), d, weights = w))
  expect_equal(f$anova$statistic, m$null.deviance - m$deviance, tolerance = 1e-6)
  expect_true(any(grepl("not all whole numbers.*average 0.01", f$notes)))
  # weights averaging 1 need no such note
  d$w1 <- d$w * 100
  d$w1[1] <- 1.0001
  f1 <- anova_bin(d, "y", "g", weights = "w1", plots = FALSE)
  expect_false(any(grepl("rather than 1", f1$notes)))

  # anova_count's dispersion and model choice do not depend on the scale
  ref <- stats::glm(c ~ g, stats::poisson(), d)
  phi <- sum(stats::residuals(ref, type = "pearson")^2) / stats::df.residual(ref)
  for (w in c("w", "w2")) {
    fc <- anova_count(d, "c", "g", weights = w, plots = FALSE)
    expect_equal(fc$dispersion, phi, tolerance = 1e-8, info = w)
    expect_identical(fc$model_type, "poisson")
    expect_true(any(grepl("rescaled to average 1", fc$notes)), info = w)
  }
  # non-constant weights: Pearson statistic over (rows - p) times their mean
  d$w3 <- withr::with_seed(2, stats::runif(300, 50, 150))
  fc <- anova_count(d, "c", "g", weights = "w3", plots = FALSE)
  m3 <- stats::glm(c ~ g, stats::poisson(), d, weights = w3)
  expect_equal(fc$dispersion,
               sum(stats::residuals(m3, type = "pearson")^2) /
                 (stats::df.residual(m3) * mean(d$w3)), tolerance = 1e-8)
  expect_identical(fc$model_type, "poisson")
})

test_that("R2-BINCOUNT-03: a negative binomial fit that collapsed to Poisson is not called non-converged", {
  skip_if_not_installed("MASS")
  withr::local_seed(1)
  d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 30)))
  d$y <- stats::rpois(90, rep(c(3, 4, 5), each = 30))
  f <- anova_count(d, "y", "g", model = "negbin", plots = FALSE)
  pois <- stats::glm(y ~ g, stats::poisson(), d)
  # the fit is the Poisson maximum likelihood, and its regression converged
  expect_equal(unname(stats::coef(f$model)), unname(stats::coef(pois)),
               tolerance = 1e-4)
  expect_true(isTRUE(f$model$converged))
  expect_gt(f$model$theta, 1000)
  expect_true(any(grepl("collapsed towards it", f$notes)))
  expect_false(any(grepl("did not converge", f$notes)))
})

test_that("R2-BINCOUNT-04: profile intervals for many coefficients come with a note on their cost", {
  withr::local_seed(20)
  k <- 32
  d <- data.frame(g = factor(rep(sprintf("g%02d", 1:k), each = 12)))
  d$y <- stats::rbinom(nrow(d), 1, 0.5)
  d$y[seq(1, nrow(d), by = 12)] <- 0
  d$y[seq(2, nrow(d), by = 12)] <- 1
  f <- anova_bin(d, "y", "g", posthoc = FALSE, plots = FALSE)
  m <- stats::glm(y ~ g, stats::binomial(), d)
  ci <- suppressMessages(stats::confint(m))
  expect_equal(f$effect_sizes$conf_low, unname(exp(ci[-1, 1])), tolerance = 1e-4)
  expect_true(any(grepl(sprintf("computed for %d odds ratios", k - 1L), f$notes)))
  # few coefficients: no note
  f2 <- anova_bin(d[d$g %in% c("g01", "g02", "g03"), ], "y", "g",
                  posthoc = FALSE, plots = FALSE)
  expect_false(any(grepl("odds ratios. Each refits", f2$notes, fixed = TRUE)))
})

test_that("R2-BINCOUNT-05: cell counts under frequency weights match by position, not name", {
  withr::local_seed(6)
  d <- data.frame(g = rep(c("", "b", "c"), each = 20), y = stats::rpois(60, 4),
                  w = sample(1:3, 60, TRUE))
  cc <- anova_count(d, "y", "g", weights = "w", plots = FALSE)$assumptions$cell_counts
  ref <- vapply(c("", "b", "c"), function(l) sum(d$w[d$g == l]), numeric(1))
  expect_equal(cc$n, unname(ref))
  expect_equal(cc$rows, c(20L, 20L, 20L))
})

test_that("R2-BINCOUNT-06: a quasi-Poisson frequency table needs no expansion", {
  withr::local_seed(9)
  tab <- expand.grid(g = c("a", "b", "c"), y = 0:8)
  tab$n <- stats::rpois(nrow(tab), 4)
  tab <- tab[tab$n > 0, ]
  long <- tab[rep(seq_len(nrow(tab)), tab$n), c("g", "y")]
  b <- anova_count(tab, "y", "g", weights = "n", model = "quasipoisson",
                   plots = FALSE)
  ref <- stats::anova(stats::glm(y ~ g, stats::quasipoisson(), long), test = "F")
  expect_equal(b$anova$statistic[1], ref$F[2], tolerance = 1e-6)
  expect_equal(b$anova$p_value[1], ref$`Pr(>F)`[2], tolerance = 1e-6)
  # which is what the (corrected) Weights documentation now says
  expect_true(any(grepl("equal those of the data expanded", b$notes)))
})

test_that("R2-BINCOUNT-07: the method names the type of test actually used", {
  withr::local_seed(8)
  d <- expand.grid(g = c("a", "b", "c"), s = c("x", "y"), rep = 1:20)
  d <- d[!(d$g == "c" & d$s == "y"), ]
  d$y <- stats::rbinom(nrow(d), 1, 0.4)
  d$cnt <- stats::rpois(nrow(d), 3)
  # the empty cell aliases a coefficient, so there is no Type III test
  expect_true(anyNA(stats::coef(stats::glm(y ~ g * s, stats::binomial(), d))))
  fb <- anova_bin(d, "y", c("g", "s"), interaction = TRUE, type = "III",
                  plots = FALSE)
  fc <- anova_count(d, "cnt", c("g", "s"), interaction = TRUE, type = "III",
                    plots = FALSE)
  expect_match(fb$method, "Type II)", fixed = TRUE)
  expect_match(fc$method, "Type II)", fixed = TRUE)
  # without aliasing, Type III is kept and named
  full <- d[d$g != "c", ]
  expect_match(anova_bin(full, "y", c("g", "s"), interaction = TRUE,
                         type = "III", plots = FALSE)$method, "Type III)",
               fixed = TRUE)
  expect_match(anova_count(full, "cnt", c("g", "s"), interaction = TRUE,
                           type = "III", plots = FALSE)$method, "Type III)",
               fixed = TRUE)
})

test_that("R2-BINCOUNT-08: the observed-counts plot survives a failure of the marginal means", {
  skip_if_not_installed("ggplot2")
  old <- getOption("emmeans")
  withr::defer(options(emmeans = old))
  emmeans::emm_options(rg.limit = 1)
  withr::local_seed(13)
  d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 20)))
  d$y <- stats::rpois(60, 4)
  fc <- anova_count(d, "y", "g")
  expect_null(fc$emmeans)
  expect_identical(names(fc$plots), "observed")
  # the plot shows the raw counts: its data hold every row's count
  expect_equal(sort(fc$plots$observed$data$value), sort(d$y))
  expect_true(any(grepl("plot of estimated marginal means was skipped",
                        fc$notes)))
  expect_false(any(grepl("Plots were skipped", fc$notes)))
})

test_that("R2-BINCOUNT-09: a coefficient unidentified in both directions gets (0, Inf)", {
  withr::local_seed(4)
  d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 20)),
                  s = factor(rep(c("n", "s"), 30)))
  d$y <- stats::rbinom(60, 1, 0.4); d$y[d$g == "b"] <- 1
  d$c <- stats::rpois(60, 3); d$c[d$g == "b"] <- 0

  # reference: the profile deviance of gb:ss, by hand, is flat in both
  # directions far beyond the chi-square cut-off's reach
  prof <- function(fam, y) {
    m <- suppressWarnings(stats::glm(y ~ g * s, fam, data = d))
    X <- stats::model.matrix(m)
    j <- which(colnames(X) == "gb:ss")
    vapply(c(-32, -5, 5, 32), function(b) {
      f <- suppressWarnings(stats::glm.fit(X[, -j], y, offset = b * X[, j],
                                           family = fam,
                                           control = stats::glm.control(maxit = 100)))
      f$deviance - stats::deviance(m)
    }, numeric(1))
  }
  expect_true(all(prof(stats::poisson(), d$c) < stats::qchisq(0.95, 1)))
  # glm.fit() cannot follow a separated logistic fit, so for the binary case
  # by reasoning instead: the g = b rows (all events) add nothing to the
  # minimum deviance, so they are fitted perfectly by gb -> Inf whatever
  # value gb:ss takes, and no value of gb:ss is ruled out.
  full <- suppressWarnings(stats::glm(y ~ g * s, stats::binomial(), data = d))
  rest <- stats::glm(y ~ g * s, stats::binomial(),
                     data = droplevels(d[d$g != "b", ]))
  expect_equal(stats::deviance(full), stats::deviance(rest), tolerance = 1e-6)

  fb <- anova_bin(d, "y", c("g", "s"), interaction = TRUE, plots = FALSE)
  rb <- fb$effect_sizes[fb$effect_sizes$term == "gb:ss", ]
  expect_identical(c(rb$conf_low, rb$conf_high), c(0, Inf))
  fc <- anova_count(d, "c", c("g", "s"), interaction = TRUE, plots = FALSE,
                    model = "poisson")
  rc <- fc$effect_sizes[fc$effect_sizes$term == "gb:ss", ]
  expect_identical(c(rc$conf_low, rc$conf_high), c(0, Inf))
  expect_true(any(grepl("gb:ss is open on both sides", fb$notes)))
  expect_true(any(grepl("gb:ss is open on both sides", fc$notes)))
  expect_false(any(grepl("search failed", c(fb$notes, fc$notes))))
})
