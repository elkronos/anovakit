# Round-5 review fixes: count and GLM functions (exact-fit rule, frequency
# weights, weight-scaled overdispersion notes, saturated Wald tests).


# Hand reference for the inflation of a weighted Poisson fit's tests when the
# counts' variance is phi * mu whatever the weights: the sandwich
# A^-1 B A^-1, A = X' diag(w mu) X, B = phi X' diag(w^2 mu) X, against A^-1,
# averaged over the eigenvalues for the term's coefficients.
hand_inflation <- function(fit, phi) {
  X <- stats::model.matrix(fit)[, -1, drop = FALSE]
  X <- cbind(1, X)
  mu <- stats::fitted(fit); w <- fit$prior.weights
  A <- crossprod(X * sqrt(w * mu)); B <- phi * crossprod(X * sqrt(w^2 * mu))
  Ai <- solve(A); V <- Ai %*% B %*% Ai
  J <- 2:ncol(X)
  lam <- eigen(solve(Ai[J, J], V[J, J]), only.values = TRUE)$values
  c(factor = mean(Re(lam)), c = sum(Re(lam)^2) / sum(Re(lam)),
    nu = sum(Re(lam))^2 / sum(Re(lam)^2))
}

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
  # the note gives the factor the tests use, from the sandwich covariance at
  # the dispersion of the weights rescaled to average 1 (hand computation),
  # and the size of the test of g, which matches the simulation
  wm <- mean(d$w)
  phi <- sum(stats::residuals(fit, type = "pearson")^2) / wm / 66
  hi <- hand_inflation(fit, phi)
  size5 <- stats::pchisq(stats::qchisq(0.95, 5) / hi[["c"]], hi[["nu"]], lower.tail = FALSE)
  note <- grep("Pearson dispersion is", k$notes, value = TRUE)
  expect_length(note, 1L)
  expect_false(grepl("little power", note, fixed = TRUE))
  expect_match(note, sprintf("are inflated by about %.2f", hi[["factor"]]), fixed = TRUE)
  expect_match(note, sprintf("a nominal 5%% test of g rejects roughly %.0f%%", 100 * size5),
               fixed = TRUE)
  expect_gt(size5, 0.35)
  # at the true dispersion (1) the sandwich size matches the simulated one
  h1 <- hand_inflation(fit, 1)
  expect_lt(abs(stats::pchisq(stats::qchisq(0.95, 5) / h1[["c"]], h1[["nu"]],
                              lower.tail = FALSE) - mean(sim)), 0.08)
  # InsectSprays: anova_count and anova_glm quote one inflation factor and
  # one size for the same weighted fit (hand computation)
  iw <- InsectSprays
  iw$w <- withr::with_seed(1, stats::runif(nrow(iw), 2, 3))
  mi <- stats::glm(count ~ spray, poisson, iw, weights = w)
  disp <- sum(stats::residuals(mi, type = "pearson")^2) / stats::df.residual(mi)
  hi <- hand_inflation(mi, disp / mean(iw$w))
  kc <- anova_count(iw, "count", "spray", model = "poisson", weights = "w",
                    test_statistic = "F", plots = FALSE)
  kg <- anova_glm(iw, "count", "spray", family = poisson(), weights = "w", plots = FALSE)
  s1 <- sprintf("roughly %.0f%%", 100 * stats::pchisq(stats::qchisq(0.95, 5) / hi[["c"]],
                                                       hi[["nu"]], lower.tail = FALSE))
  fac <- sprintf("inflated by about %.2f", hi[["factor"]])
  expect_true(any(grepl(fac, kc$notes, fixed = TRUE)))
  expect_true(any(grepl(sprintf("Pearson residuals (%.2f)", disp), kc$notes, fixed = TRUE)))
  expect_true(any(grepl(s1, kc$notes, fixed = TRUE)))
  # these weights are nearly equal (mean(w^2) / mean(w)^2 = 1.01), so their
  # unevenness moves no test and anova_glm keeps its overdispersion note: the
  # Pearson dispersion of the weighted fit is then the inflation of its tests
  expect_lt(mean(iw$w^2) / mean(iw$w)^2, 1.1)
  expect_true(any(grepl(sprintf("Pearson dispersion is %.2f", disp), kg$notes, fixed = TRUE)))
  expect_false(any(grepl("are unequal and not frequencies", kg$notes, fixed = TRUE)))
  expect_equal(hi[["factor"]], disp, tolerance = 0.05)
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

test_that("countglm-R1/C1: a nonzero response constant within groups on a large level is an exact fit", {
  # reference: the response is constant within every group (base R tapply),
  # so the residuals glm() leaves are rounding error (~1e-14 of the level)
  # however the deviance ratio comes out
  chk <- function(d, ...) {
    expect_true(all(tapply(d$y, d$g, function(v) diff(range(v))) == 0))
    expect_error(anova_glm(d, "y", "g", plots = FALSE, ...),
                 "fits response `y` exactly: the response is constant within every design cell")
  }
  d <- data.frame(g = gl(3, 5), y = 1e7 + rep(c(0, 10, 20), each = 5))
  # the deviance rule alone misses it: glm's deviance is not ~0 of the null
  gq <- suppressWarnings(stats::glm(y ~ g, stats::quasipoisson(), d))
  expect_gt(abs(gq$deviance / gq$null.deviance), 1e-12)
  chk(d, family = stats::quasipoisson())
  expect_error(anova_count(d, "y", "g", model = "quasipoisson", plots = FALSE),
               "fits response `y` exactly")
  chk(transform(d, y = 1.7e9 + rep(c(0, 10, 20), each = 5)), family = stats::quasipoisson())
  chk(transform(d, y = 1e8 + rep(c(0, 10, 20), each = 5)), family = stats::Gamma("log"))
  chk(data.frame(g = gl(3, 5), y = 1e10 + rep(c(0, 1, 2), each = 5)))
  chk(data.frame(g = gl(3, 4), y = rep(1.7e9 + c(0, .1, .2), each = 4)))
  chk(data.frame(g = gl(3, 4), y = rep(1.7e9 + c(0, .01, .02), each = 4)))
  chk(data.frame(g = gl(2, 5), y = rep(c(1e7, 1e7 + 1e-3), each = 5)))
  chk(data.frame(g = gl(2, 5), y = rep(c(1e7, 1e7 + 1), each = 5)), family = stats::quasipoisson())
  chk(data.frame(g = gl(2, 5), y = rep(c(1e7, 1e7 + 1e-2), each = 5)), family = stats::Gamma("log"))
  # real variation on the same level, however small against it, is tested:
  # 1e-3 noise on 1e10 (1e-13 of the level) matches lm()'s F to its rounding
  d6 <- withr::with_seed(1, data.frame(g = gl(3, 10),
    y = 1e10 + rep(c(0, 1, 2), each = 10) + stats::rnorm(30, 0, 1e-3)))
  ref <- suppressWarnings(stats::anova(stats::lm(y ~ g, d6)))
  f6 <- anova_glm(d6, "y", "g", plots = FALSE)
  expect_equal(f6$anova$statistic[1], ref[["F value"]][1], tolerance = 1e-2)
  # an additive two-way model whose cells are constant but not additive
  # leaves real residual variation, which is tested (lm reference)
  d2 <- expand.grid(a = gl(2, 1), b = gl(2, 1), r = 1:3)
  d2$y <- 1e7 + c(0, 1, 2, 7)[as.integer(d2$a) + 2L * (as.integer(d2$b) - 1L)]
  r2 <- suppressWarnings(stats::anova(stats::lm(y ~ a + b, d2)))
  f2 <- anova_glm(d2, "y", c("a", "b"), interaction = FALSE, type = "II", plots = FALSE)
  expect_equal(f2$anova$statistic[1:2], r2[["F value"]][1:2], tolerance = 1e-6)
})

test_that("countglm-C2: an all-boundary exact fit has no 'elsewhere' clause", {
  db <- data.frame(g = factor(c("a", "b", "c")), y = c(0, 1, 1), w = c(4, 5, 6))
  e <- tryCatch(anova_glm(db, "y", "g", family = stats::quasibinomial(),
                          weights = "w", plots = FALSE),
                error = conditionMessage)
  expect_match(e, "every response is 0 or 1, and every fitted value equals it", fixed = TRUE)
  expect_false(grepl("elsewhere", e, fixed = TRUE))
})

test_that("countglm-C3: messages name 'the inverse.gaussian family'", {
  dd <- data.frame(y = c(1, 2, 4, 3, 5, 6, 2, 3, 5), g = rep(c("a", "b", "c"), 3))
  n <- anova_glm(dd, "y", "g", family = "inverse.gaussian", test_statistic = "LR",
                 plots = FALSE)$notes
  expect_true(any(grepl("requested for the inverse.gaussian family", n, fixed = TRUE)))
  expect_false(any(grepl("a inverse", n, fixed = TRUE)))
  dd$c <- as.character(dd$y)
  expect_error(anova_glm(dd, "c", "g", family = "inverse.gaussian", plots = FALSE),
               "must be numeric for the inverse.gaussian family", fixed = TRUE)
})

test_that("countglm-R2: weights averaging 2.5 inflate the kept Poisson tests plainly, not hedged", {
  d <- withr::with_seed(2, {
    d <- data.frame(g = gl(4, 1000))
    d$y <- stats::rpois(4000, 6)
    d$w <- stats::runif(4000, 2, 3)
    d
  })
  r <- anova_count(d, "y", "g", plots = FALSE, weights = "w")
  expect_false(any(grepl("little power", r$notes, fixed = TRUE)))
  expect_false(any(grepl("cannot tell the dispersion from 1", r$notes, fixed = TRUE)))
  note <- grep("A Poisson model assumes a dispersion of 1", r$notes, value = TRUE)
  expect_length(note, 1L)
  # reference: the sandwich inflation by hand, at the dispersion of the
  # weights rescaled to average 1; about the Pearson dispersion of the glm as
  # weighted, since the weights vary little
  f <- stats::glm(y ~ g, stats::poisson(), d, weights = w)
  phi <- sum(stats::residuals(f, "pearson")^2) / mean(d$w) / stats::df.residual(f)
  hi <- hand_inflation(f, phi)
  expect_match(note, sprintf("are inflated by about %.2f", hi[["factor"]]), fixed = TRUE)
  sz3 <- stats::pchisq(stats::qchisq(0.95, 3) / hi[["c"]], hi[["nu"]], lower.tail = FALSE)
  expect_match(note, sprintf("a nominal 5%% test of g rejects roughly %.0f%%", 100 * sz3),
               fixed = TRUE)
  # simulation: true Poisson data, weights ~ U(2, 3): the LR test of g on 3
  # df rejects far more than 5% (about the quoted size)
  rej <- withr::with_seed(11, mean(replicate(300, {
    ds <- data.frame(g = gl(4, 50), w = stats::runif(200, 2, 3))
    ds$y <- stats::rpois(200, 6)
    m1 <- stats::glm(y ~ g, stats::poisson(), ds, weights = w)
    m0 <- stats::glm(y ~ 1, stats::poisson(), ds, weights = w)
    stats::pchisq(m0$deviance - m1$deviance, 3, lower.tail = FALSE) < 0.05
  })))
  expect_gt(rej, 0.25)
  expect_lt(abs(rej - sz3), 0.1)
  # weights barely above 1 on 27 df stay hedged: the inflation is then the
  # uncertain dispersion's, and the Pearson test has little power against
  # the dispersion at which the tests would reject 7.5% (hand computation)
  ds <- withr::with_seed(1, {
    ds <- data.frame(g = gl(3, 10))
    ds$y <- stats::rpois(30, 5)
    ds$w <- stats::runif(30, 1, 1.2)
    ds
  })
  rs <- anova_count(ds, "y", "g", plots = FALSE, weights = "w")
  expect_true(any(grepl("27 residual degrees of freedom it has little power", rs$notes,
                        fixed = TRUE)))
  phi_c <- .count_phi_mattering(2L, mean(ds$w))
  expect_equal(stats::pchisq(stats::qchisq(0.95, 2) / (phi_c * mean(ds$w)), 2,
                             lower.tail = FALSE), 0.075, tolerance = 1e-10)
  expect_lt(stats::pchisq(stats::qchisq(0.95, 27) / phi_c, 27, lower.tail = FALSE), 0.5)
})
