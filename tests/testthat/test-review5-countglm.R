# Round-5 review fixes: count and GLM functions (exact-fit rule, frequency
# weights, weight-scaled overdispersion notes, saturated Wald tests).

test_that("countglm-F1: a large level with small real variation is not an exact fit", {
  d <- withr::with_seed(4, {
    d <- data.frame(g = gl(3, 10))
    d$t <- 1.7e9 + rep(c(0, 60, 120), each = 10) + stats::rnorm(30, 0, 300)
    d
  })
  # reference: lm() resolves the residuals (relative size ~1e-7) as real
  ref <- suppressWarnings(stats::anova(stats::lm(t ~ g, d)))
  f <- anova_glm(d, "t", "g", plots = FALSE)
  expect_equal(f$anova$statistic[1], ref[["F value"]][1], tolerance = 1e-8)
  expect_equal(f$anova$p_value[1], ref[["Pr(>F)"]][1], tolerance = 1e-8)
  expect_equal(f$anova$statistic[1], 2.1897, tolerance = 1e-4)
  # Gamma(log): car's F test on the same glm
  # (fitted to convergence, as anova_glm refits a deviance in small units)
  rg <- car::Anova(stats::glm(t ~ g, Gamma("log"), d,
                              control = stats::glm.control(epsilon = 1e-14, maxit = 100)),
                   test.statistic = "F")
  fg <- anova_glm(d, "t", "g", family = Gamma("log"), plots = FALSE)
  expect_equal(fg$anova$statistic[1], rg[["F value"]][1], tolerance = 1e-6)
  # pressures around 101325 Pa with an sd of 0.05
  dp <- withr::with_seed(5, data.frame(g = gl(3, 10),
    y = 101325 + rep(c(0, 0.03, 0.06), each = 10) + stats::rnorm(30, 0, 0.05)))
  rp <- suppressWarnings(stats::anova(stats::lm(y ~ g, dp)))
  fp <- anova_glm(dp, "y", "g", plots = FALSE)
  expect_equal(fp$anova$statistic[1], rp[["F value"]][1], tolerance = 1e-6)
  # a response constant within every group, one group at 0, is still exact
  dz <- data.frame(g = gl(3, 10), y = rep(c(0, 5, 1e7), each = 10))
  mz <- suppressWarnings(stats::glm(y ~ g, quasipoisson, dz))
  expect_lt(max(abs(mz$y - fitted(mz))), 1e-6)   # reference: glm's limit fit
  expect_error(anova_glm(dz, "y", "g", family = quasipoisson(), plots = FALSE),
               "fits response `y` exactly")
})

test_that("exactfit-dispersion-F1: a Gaussian response on a large baseline gets car's F test", {
  d <- withr::with_seed(4, {
    d <- data.frame(g = factor(rep(1:4, each = 8)))
    d$y <- 2.5e6 + rep(c(0, .2, .4, .1), each = 8) + stats::rnorm(32, sd = .3)
    d
  })
  ref <- car::Anova(stats::lm(y ~ g, d))
  f <- anova_glm(d, "y", "g", family = gaussian(), plots = FALSE)
  expect_equal(f$anova$statistic[1], ref[["F value"]][1], tolerance = 1e-8)
  expect_equal(f$anova$p_value[1], ref[["Pr(>F)"]][1], tolerance = 1e-8)
  expect_equal(f$anova$statistic[1], 6.2789, tolerance = 1e-4)
  # the likelihood-ratio test and the log link are not refused either
  rl <- car::Anova(stats::glm(y ~ g, gaussian, d), test.statistic = "LR")
  fl <- anova_glm(d, "y", "g", family = gaussian(), test_statistic = "LR", plots = FALSE)
  expect_equal(fl$anova$statistic[1], rl[["LR Chisq"]][1], tolerance = 1e-6)
  expect_no_error(anova_glm(d, "y", "g", family = gaussian("log"), plots = FALSE))
  # epoch seconds with an sd of 5
  d2 <- withr::with_seed(2, data.frame(g = gl(3, 10),
    y = 1.7e9 + rep(c(0, 3, 6), each = 10) + stats::rnorm(30, sd = 5)))
  r2 <- suppressWarnings(stats::anova(stats::lm(y ~ g, d2)))
  f2 <- anova_glm(d2, "y", "g", plots = FALSE)
  expect_equal(f2$anova$statistic[1], r2[["F value"]][1], tolerance = 1e-6)
})

test_that("countglm-F2: the pseudo-R2 note names the inverse.gaussian family grammatically", {
  d <- data.frame(y = c(1, 2, 4, 3, 5, 6), g = rep(c("a", "b", "c"), 2))
  f <- anova_glm(d, "y", "g", family = "inverse.gaussian", plots = FALSE)
  note <- grep("McFadden", f$notes, value = TRUE)
  expect_length(note, 1L)
  expect_match(note, "the log-likelihood of the inverse.gaussian family", fixed = TRUE)
  expect_false(grepl("of a inverse", note, fixed = TRUE))
  # reference: the value is NA and the deviance explained is reported
  m <- stats::glm(y ~ g, inverse.gaussian, transform(d, g = factor(g)))
  expect_true(is.finite(1 - m$deviance / m$null.deviance))
})

test_that("exactfit-dispersion-F2: frequency weights saturated in rows but not in observations are an exact fit", {
  d <- data.frame(g = factor(c("a", "b", "c")), y = c(3, 7, 12), w = c(5, 6, 7))
  ex <- d[rep(1:3, d$w), ]
  # reference: on the expanded data the dispersion is rounding error
  m <- stats::glm(y ~ g, quasipoisson, ex)
  expect_equal(stats::df.residual(m), 15L)
  expect_lt(summary(m)$dispersion, 1e-20)
  expect_error(anova_glm(ex, "y", "g", family = quasipoisson(), plots = FALSE),
               "fits response `y` exactly")
  # the frequency-weighted call gives the same outcome
  expect_error(anova_glm(d, "y", "g", family = quasipoisson(), weights = "w",
                         plots = FALSE), "fits response `y` exactly")
  expect_error(anova_count(d, "y", "g", model = "quasipoisson", weights = "w",
                           plots = FALSE), "fits response `y` exactly")
  expect_error(anova_count(d, "y", "g", model = "poisson", weights = "w",
                           test_statistic = "F", plots = FALSE),
               "fits response `y` exactly")
  # a fixed-dispersion LR test of the same table is defined (statistic as on
  # the expanded data) and is kept
  fp <- anova_glm(d, "y", "g", family = poisson(), weights = "w", plots = FALSE)
  rp <- car::Anova(stats::glm(y ~ g, poisson, ex))
  expect_equal(fp$anova$statistic[1], rp[["LR Chisq"]][1], tolerance = 1e-8)
  # frequency weights with residual variation still match the expanded data
  d2 <- data.frame(g = factor(rep(c("a", "b", "c"), each = 2)),
                   y = c(3, 5, 7, 9, 12, 10), w = c(2, 3, 4, 2, 3, 3))
  ex2 <- d2[rep(seq_len(6), d2$w), ]
  r2 <- car::Anova(stats::glm(y ~ g, quasipoisson, ex2), test.statistic = "F")
  f2 <- anova_glm(d2, "y", "g", family = quasipoisson(), weights = "w", plots = FALSE)
  expect_equal(f2$anova$statistic[1], r2[["F value"]][1], tolerance = 1e-8)
  k2 <- anova_count(d2, "y", "g", model = "quasipoisson", weights = "w", plots = FALSE)
  expect_equal(k2$anova$statistic[1], r2[["F value"]][1], tolerance = 1e-8)
})

test_that("exactfit-dispersion-F3: anova_count's inflation note describes the tests as weighted", {
  d <- withr::with_seed(3, {
    d <- data.frame(g = factor(rep(letters[1:6], each = 12)))
    d$w <- stats::runif(72, 2, 3)
    d$y <- stats::rpois(72, 8)
    d
  })
  k <- anova_count(d, "y", "g", model = "poisson", weights = "w", plots = FALSE)
  # reference: under a true dispersion of 1 the reported LR test, computed
  # with the weights as given, rejects far more than 10% of true nulls
  fit <- stats::glm(y ~ g, poisson, d, weights = w)
  sim <- withr::with_seed(11, replicate(300, {
    ys <- stats::rpois(72, 8)
    f1 <- stats::glm(ys ~ g, poisson, d, weights = w)
    f0 <- stats::glm(ys ~ 1, poisson, d, weights = w)
    stats::pchisq(f0$deviance - f1$deviance, 5, lower.tail = FALSE) < 0.05
  }))
  expect_gt(mean(sim), 0.35)
  # the note gives the factor the tests use: the dispersion bound times the
  # mean weight, and the sizes at that factor
  wm <- mean(d$w)
  x2 <- sum(stats::residuals(fit, type = "pearson")^2) / wm
  phi_hi <- x2 / stats::qchisq(0.05, 66)
  size5 <- stats::pchisq(stats::qchisq(0.95, 5) / (phi_hi * wm), 5, lower.tail = FALSE)
  note <- grep("Pearson dispersion is", k$notes, value = TRUE)
  expect_length(note, 1L)
  expect_match(note, sprintf("inflated by about %.2f", phi_hi * wm), fixed = TRUE)
  expect_match(note, sprintf("one on 5 (the largest term here) roughly %.0f%%", 100 * size5),
               fixed = TRUE)
  expect_gt(size5, 0.35)
  # InsectSprays: anova_count and anova_glm quote one dispersion and one size
  iw <- InsectSprays
  iw$w <- withr::with_seed(1, stats::runif(nrow(iw), 2, 3))
  mi <- stats::glm(count ~ spray, poisson, iw, weights = w)
  disp <- sum(stats::residuals(mi, type = "pearson")^2) / stats::df.residual(mi)
  kc <- anova_count(iw, "count", "spray", model = "poisson", weights = "w",
                    test_statistic = "F", plots = FALSE)
  kg <- anova_glm(iw, "count", "spray", family = poisson(), weights = "w", plots = FALSE)
  s1 <- sprintf("roughly %.0f%%", 100 * stats::pchisq(stats::qchisq(0.95, 1) / disp, 1,
                                                       lower.tail = FALSE))
  expect_true(any(grepl(sprintf("inflated by about %.2f", disp), kc$notes, fixed = TRUE)))
  expect_true(any(grepl(sprintf("Pearson residuals (%.2f)", disp), kc$notes, fixed = TRUE)))
  expect_true(any(grepl(s1, kc$notes, fixed = TRUE)))
  expect_true(any(grepl(sprintf("Pearson dispersion is %.2f", disp), kg$notes, fixed = TRUE)))
  expect_true(any(grepl(s1, kg$notes, fixed = TRUE)))
})

test_that("exactfit-dispersion-F4: a saturated Poisson fit keeps its Wald test", {
  d <- data.frame(g = factor(c("a", "b", "c")), y = c(3, 7, 12))
  m <- stats::glm(y ~ g, poisson, d)
  b <- stats::coef(m)[-1]
  V <- stats::vcov(m)[-1, -1]
  W <- drop(t(b) %*% solve(V) %*% b)   # reference: Wald chi-square by hand
  expect_equal(W, 4.939, tolerance = 1e-3)
  for (f in list(
    anova_glm(d, "y", "g", family = poisson(), test_statistic = "Wald", plots = FALSE),
    anova_count(d, "y", "g", model = "poisson", test_statistic = "Wald", plots = FALSE))) {
    expect_false(is.null(f$anova))
    expect_equal(f$anova$statistic[f$anova$term == "g"], W, tolerance = 1e-8)
    expect_equal(f$anova$p_value[f$anova$term == "g"],
                 stats::pchisq(W, 2, lower.tail = FALSE), tolerance = 1e-8)
    expect_false(any(grepl("could not be computed", f$notes)))
    expect_true(any(grepl("Wald tests are still valid", f$notes)))
  }
  # a saturated fit with an estimated dispersion still has NA inference
  q <- anova_glm(d, "y", "g", family = quasipoisson(), test_statistic = "Wald",
                 plots = FALSE)
  expect_true(all(is.na(q$anova$p_value[q$anova$term == "g"])) || is.null(q$anova))
})
