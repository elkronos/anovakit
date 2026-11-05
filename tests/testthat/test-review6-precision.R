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

test_that("binomial trial weights are inverse variances: no inflation note", {
  d <- withr::with_seed(7, {
    d <- data.frame(g = gl(4, 15))
    d$n <- sample(c(5, 10, 200), 60, replace = TRUE)
    d$y <- stats::rbinom(60, d$n, 0.3) / d$n
    d
  })
  f <- anova_glm(d, "y", "g", family = binomial(), weights = "n", plots = FALSE,
                 posthoc = FALSE)
  expect_false(any(grepl("sandwich covariance", f$notes, fixed = TRUE)))
})

test_that("Type II weight inflation of a main effect is that of the model without its interaction", {
  d <- data.frame(g = gl(4, 50), h = factor(rep(c("x", "y"), 100)))
  d$w <- ifelse(seq_len(200) %% 10 == 0, 9.1, 0.1)
  d$y <- withr::with_seed(3, stats::rpois(200, 5))
  m <- stats::glm(y ~ g * h, stats::poisson(), d, weights = w)
  m3 <- stats::glm(y ~ g * h, stats::poisson(), d, weights = w,
                   contrasts = list(g = "contr.sum", h = "contr.sum"))
  # at the true dispersion (1), the sandwich sizes of the tests of g
  w2 <- .weight_inflation(m, 1, "II")
  w3 <- .weight_inflation(m3, 1, "III")
  # simulation: true Poisson counts, car's LR tests of g with these weights
  sim <- withr::with_seed(9, rowMeans(replicate(200, {
    ys <- stats::rpois(200, 5)
    m <- stats::glm(ys ~ g * h, stats::poisson(), d, weights = w)
    m3 <- stats::glm(ys ~ g * h, stats::poisson(), d, weights = w,
                     contrasts = list(g = "contr.sum", h = "contr.sum"))
    c(car::Anova(m, type = 2)["g", 3] < 0.05, car::Anova(m3, type = 3)["g", 3] < 0.05)
  })))
  expect_lt(abs(w2$size[w2$term == "g"] - sim[1]), 0.1)
  expect_lt(abs(w3$size[w3$term == "g"] - sim[2]), 0.1)
  # Type II and III differ here (the heavy rows all have h = "y")
  expect_gt(w2$size[w2$term == "g"], 3 * w3$size[w3$term == "g"])
  f2 <- anova_count(d, "y", c("g", "h"), model = "poisson", weights = "w",
                    interaction = TRUE, plots = FALSE, posthoc = FALSE)
  expect_true(any(grepl("test of g rejects roughly", f2$notes, fixed = TRUE)))
})

test_that("a negative binomial fit with unequal non-frequency weights says its tests are inflated", {
  skip_if_not_installed("MASS")
  d <- data.frame(a = factor(rep(1:3, length.out = 300)),
                  b = factor(rep(1:2, each = 3, length.out = 300)))
  d$w <- ifelse(seq_len(300) %% 10 == 0, 9.1, 0.1)
  d$y <- withr::with_seed(2, { invisible(MASS::rnegbin(300, 5, 2)); stats::rpois(300, 5) })
  f <- anova_count(d, "y", c("a", "b"), model = "negbin", weights = "w",
                   plots = FALSE, posthoc = FALSE)
  note <- grep("negative binomial tests are inflated", f$notes, value = TRUE)
  expect_length(note, 1L)
  q <- as.numeric(sub(".*test of a rejects roughly ([0-9]+)%.*", "\\1", note)) / 100
  sim <- withr::with_seed(4, mean(replicate(200, {
    ys <- stats::rpois(300, 5)
    m <- suppressWarnings(MASS::glm.nb(ys ~ a + b, d, weights = w))
    car::Anova(m)["a", 3] < 0.05
  })))
  expect_lt(abs(q - sim), 0.12)
})

test_that("a Gaussian fit on a large level with many rows keeps its F test", {
  d <- withr::with_seed(3, data.frame(g = gl(2, 5000),
                                      y = round(1.7e9 + stats::rnorm(10000, sd = 0.5), 3)))
  f <- anova_glm(d, "y", "g", plots = FALSE, posthoc = FALSE)
  ref <- stats::anova(stats::lm(I(y - 1.7e9) ~ g, d))
  expect_equal(f$anova$statistic[1], ref[["F value"]][1], tolerance = 1e-6)
  expect_false(any(grepl("rounding", f$notes, fixed = TRUE)))
})

test_that("the exact-fit guard leaves real lack of fit and Gaussian zeros alone", {
  g <- factor(rep(c("a", "b", "c"), times = 3)); h <- factor(rep(c("x", "y", "z"), each = 3))
  d <- data.frame(g, h, y = 1e12 * (c(1, 2, 3)[g] + c(0, 1, 3)[h]) +
                    10 * c(1, -1, 0, -1, 1, 0, 0, 0, 0))
  ref <- car::Anova(stats::lm(y ~ g + h, d))
  f <- anova_glm(d, "y", c("g", "h"), plots = FALSE)
  expect_equal(f$anova$statistic[1:2], ref[["F value"]][1:2], tolerance = 1e-3)
  # counts of 1.5e12 that differ by 1 are not one value
  d2 <- data.frame(g = factor(rep(c("a", "b"), each = 4)),
                   y = 1.5e12 + c(0, 1, 0, 1, 1, 1, 1, 0))
  expect_s3_class(anova_glm(d2, "y", "g", plots = FALSE), "anovakit_fit")
  # a Gaussian group at 0 is an ordinary value
  e <- c(-1, 1, 1, -1)
  d3 <- data.frame(g = factor(rep(c("a", "b", "c"), each = 4)),
                   y = c(rep(0, 4), 1e9 + 0.05 * e, 2e9 + 0.05 * e))
  expect_s3_class(anova_glm(d3, "y", "g", plots = FALSE), "anovakit_fit")
})

test_that("a negative Poisson likelihood-ratio statistic within rounding is shown as 0", {
  g <- factor(rep(c("a", "b", "c"), each = 6))
  d <- data.frame(g = g, y = 1e12 + c(0, 5, 10)[as.integer(g)] + rep(c(-3, -1, 0, 0, 1, 3), 3))
  f <- anova_count(d, "y", "g", model = "poisson", plots = FALSE, posthoc = FALSE)
  expect_gte(f$anova$statistic[1], 0)
  expect_true(any(grepl("limited precision", f$notes, fixed = TRUE)))
})
