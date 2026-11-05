# Round-6 review fixes: deviance precision on a large level, the exact-fit
# guard, the inflation of tests by weights that are not frequencies, and
# small Welch and count notes.

test_that("a Gamma or quasi-Poisson test that rounding swamps is NA, with a Wald test still available", {
  g <- factor(rep(c("a", "b", "c"), each = 6))
  y <- withr::with_seed(42, 1e6 * (1 + 1e-9 * stats::rnorm(18))) + c(0, 1, 2)[g] * 1e-3
  d <- data.frame(g, y)
  f <- anova_glm(d, "y", "g", family = Gamma("inverse"), test_statistic = "F",
                 plots = FALSE, posthoc = FALSE)
  expect_true(is.na(f$anova$statistic[1]))
  expect_true(is.na(f$anova$p_value[1]))
  expect_true(any(grepl("varies too little relative to its level", f$notes, fixed = TRUE)))
  # counts of 1e12 differing by a few: car gives a negative F; now NA
  dn <- data.frame(g = g, y = 1e12 + c(0, 5, 10)[as.integer(g)] +
                     rep(c(-3, -1, 0, 0, 1, 3), 3))
  fq <- anova_glm(dn, "y", "g", family = quasipoisson(), plots = FALSE, posthoc = FALSE)
  expect_true(is.na(fq$anova$statistic[1]))
  fc <- anova_count(dn, "y", "g", model = "quasipoisson", plots = FALSE, posthoc = FALSE)
  expect_true(is.na(fc$anova$statistic[1]))
  # the Wald test does not use the deviance: b' V^-1 b by hand
  fw <- anova_glm(dn, "y", "g", family = quasipoisson(), test_statistic = "Wald",
                  plots = FALSE, posthoc = FALSE)
  m <- stats::glm(y ~ g, stats::quasipoisson(), dn)
  b <- stats::coef(m)[-1]
  W <- drop(t(b) %*% solve(stats::vcov(m)[-1, -1], b))
  expect_equal(fw$anova$statistic[1], W, tolerance = 1e-6)
  # Gaussian on the same counts resolves them. By hand: group means 1e12 + 0,
  # 5, 10 exactly; between mean square 6 * 50 / 2 = 150, within 60 / 15 = 4
  fg <- anova_glm(dn, "y", "g", plots = FALSE, posthoc = FALSE)
  expect_equal(fg$anova$statistic[1], 150 / 4, tolerance = 1e-6)
})

test_that("a Poisson likelihood-ratio test on a large level is kept: it is not divided by a dispersion", {
  dc <- withr::with_seed(8, data.frame(g = gl(2, 6),
    y = round(1e9 + stats::rnorm(12, sd = 20)) + rep(c(0, 30), each = 6)))
  f <- anova_count(dc, "y", "g", model = "poisson", plots = FALSE)
  m1 <- stats::glm(y ~ g, stats::poisson(), dc)
  m0 <- stats::glm(y ~ 1, stats::poisson(), dc)
  expect_equal(f$anova$statistic[1], m0$deviance - m1$deviance, tolerance = 1e-6)
  expect_false(any(grepl("varies too little", f$notes, fixed = TRUE)))
})

test_that("real variation beside a large constant group is analysed, not refused", {
  da <- withr::with_seed(1, data.frame(g = gl(3, 5),
    y = c(rep(1e7, 5), 1 + stats::rnorm(5, 0, 0.01), rep(5, 5))))
  ref <- suppressWarnings(stats::anova(stats::lm(y ~ g, da)))
  f <- anova_glm(da, "y", "g", plots = FALSE)
  expect_equal(f$anova$statistic[1], ref[["F value"]][1], tolerance = 1e-6)
  # counts of 1e8 varying by a few, beside a group at half that level
  dd <- data.frame(g = gl(3, 4), y = c(1e8 + c(0, 7, 3, 11), 1e8 + c(5, 2, 9, 1),
                                       5e7 + c(4, 0, 6, 2)))
  k <- anova_count(dd, "y", "g", model = "quasipoisson", plots = FALSE)
  rc <- car::Anova(stats::glm(y ~ g, stats::quasipoisson(), dd), test.statistic = "F")
  expect_equal(k$anova$statistic[1], rc[["F value"]][1], tolerance = 1e-3)
})

test_that("a response constant within every cell on a large level is still an exact fit", {
  dk <- data.frame(g = gl(3, 4), y = rep(1e7 + c(0, 10, 20), each = 4))
  expect_error(anova_glm(dk, "y", "g", family = quasipoisson(), plots = FALSE),
               "fits response `y` exactly")
  expect_error(anova_glm(dk, "y", "g", plots = FALSE), "fits response `y` exactly")
  # but whole-number counts of 1e12 that differ are not "the same value"
  d2 <- data.frame(g = gl(2, 6), y = 1e12 + c(0, 10, 20, 30, 40, 50, 5, 15, 25, 35, 45, 55))
  r <- tryCatch(anova_count(d2, "y", "g", model = "quasipoisson", plots = FALSE),
                error = function(e) conditionMessage(e))
  expect_false(is.character(r) && grepl("takes the same value", r, fixed = TRUE))
})

test_that("unequal weights that are not frequencies: the note's size matches simulation", {
  d <- expand.grid(a = gl(3, 1), b = gl(2, 1), r = 1:50)
  d$w <- ifelse(seq_len(300) %% 10 == 0, 9.1, 0.1)
  d$y <- withr::with_seed(11, stats::rpois(300, 5))
  k <- anova_count(d, "y", c("a", "b"), model = "poisson", weights = "w",
                   plots = FALSE, posthoc = FALSE)
  note <- grep("A Poisson model assumes a dispersion of 1", k$notes, value = TRUE)
  expect_length(note, 1L)
  q <- as.numeric(sub(".*test of a rejects roughly ([0-9]+)%.*", "\\1", note)) / 100
  # simulation: true Poisson counts, the same weights, car's type II LR test of a
  sim <- withr::with_seed(5, mean(replicate(300, {
    ys <- stats::rpois(300, 5)
    car::Anova(stats::glm(ys ~ a + b, stats::poisson(), d, weights = w))[1, 3] < 0.05
  })))
  expect_gt(sim, 0.3)
  expect_lt(abs(q - sim), 0.12)
  # anova_glm quotes the same size for the same fit
  g <- anova_glm(d, "y", c("a", "b"), family = poisson(), weights = "w",
                 plots = FALSE, posthoc = FALSE)
  expect_true(any(grepl(sprintf("test of a rejects roughly %.0f%%", 100 * q), g$notes,
                        fixed = TRUE)))
})

test_that("Welch's variance-ratio note fires above 4, not at exactly 4", {
  x <- c(8.72, 9.43, 8.78, 9.53, 9.38, 10.04, 9.09, 10.16)
  d <- data.frame(y = c(x, 2 * x), g = rep(c("a", "b"), each = 8))
  v <- tapply(d$y, d$g, stats::var)
  expect_equal(unname(v[["b"]] / v[["a"]]), 4)
  f <- anova_welch(d, "y", "g", plots = FALSE)
  expect_false(any(grepl("Largest group variance", f$notes, fixed = TRUE)))
  # an overflowed ratio with no omnibus F does not say Welch's test handles it
  d2 <- data.frame(g = rep(c("a", "b"), each = 4),
                   y = c(1e154 * c(1, 2, 4, 3), 1e-154 * c(1, 2, 4, 3)))
  f2 <- anova_welch(d2, "y", "g", plots = FALSE)
  if (is.na(f2$anova$statistic[1])) {
    expect_false(any(grepl("Welch's test handles this", f2$notes, fixed = TRUE)))
  }
})

test_that("a negative binomial on one row per group of frequency weights says why it fails", {
  skip_if_not_installed("MASS")
  dn <- data.frame(g = factor(c("a", "b", "c")), y = c(3, 7, 12), w = c(5, 6, 7))
  r <- tryCatch(anova_count(dn, "y", "g", model = "negbin", weights = "w", plots = FALSE),
                error = function(e) conditionMessage(e))
  if (is.character(r)) {
    expect_false(grepl("theta cannot be estimated when the model has no residual",
                       r, fixed = TRUE))
  } else {
    expect_s3_class(r, "anovakit_fit")
  }
})
