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
  # Where long double is no wider than double (ARM macOS) the variance of
  # values of 1e-154 underflows, and the function says so instead.
  f2 <- tryCatch(anova_welch(d2, "y", "g", plots = FALSE), error = function(e) e)
  if (inherits(f2, "error")) {
    expect_match(conditionMessage(f2), "cannot be represented in double precision",
                 fixed = TRUE)
  } else if (is.na(f2$anova$statistic[1])) {
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
  # The true LR statistic is about 3e-10; car's comes out as rounding
  # error of about +-0.004 (its sign depends on the platform). A negative
  # one is shown as 0, never below; either way p is 1 to within that error.
  expect_gte(f$anova$statistic[1], 0)
  expect_lt(f$anova$statistic[1], 0.01)
  expect_gt(f$anova$p_value[1], 0.99)
  # an ordinary null effect gets no precision note
  d2 <- data.frame(g = gl(2, 4), y = c(8, 12, 16, 12, 12, 16, 8, 12))
  f2 <- anova_count(d2, "y", "g", model = "poisson", plots = FALSE)
  expect_false(any(grepl("precision", f2$notes, fixed = TRUE)))
})

test_that("Gaussian responses constant within groups are refused, a zero group included", {
  d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 4)))
  for (v in list(c(0, 1, 3), c(0, 1e9, 2e9), c(1e9, 0, 2e9))) {
    d$y <- rep(v, each = 4)
    expect_error(anova_glm(d, "y", "g", plots = FALSE), "fits response `y` exactly")
  }
  # real variation of 1e-6 at level 1 beside a group at 1e8 is analysed (lm)
  d5 <- data.frame(g = gl(3, 5), y = c(rep(1e8, 5), 1 + c(-1, 1, 0, 1, -1) * 1e-6, rep(5, 5)))
  expect_s3_class(anova_glm(d5, "y", "g", plots = FALSE), "anovakit_fit")
})

test_that("an exactly multiplicative quasi-Poisson response near 1 is an exact fit", {
  d <- expand.grid(a = gl(3, 1), b = gl(3, 1))
  d$y <- with(d, 0.5 * c(1, 2, 3)[a] * c(1, 1.5, 2)[b])
  expect_error(anova_glm(d, "y", c("a", "b"), family = quasipoisson(), plots = FALSE),
               "fits response `y` exactly")
  d$y2 <- d$y * (1 + 1e-9 * c(1, -1, 0, 0, 1, -1, -1, 0, 1))
  expect_s3_class(anova_glm(d, "y2", c("a", "b"), family = quasipoisson(), plots = FALSE),
                  "anovakit_fit")
})

test_that("a saturated GLM keeps its deviance explained of 1", {
  sat <- data.frame(g = gl(3, 1), s = c(3, 7, 12))
  f <- anova_glm(sat, "s", "g", family = poisson(), plots = FALSE)
  expect_equal(f$effect_sizes$estimate[f$effect_sizes$measure == "deviance_explained"], 1)
})

test_that("binomial sampling weights and quasi-Poisson weights get the inflation note; trials do not", {
  d <- data.frame(g = gl(3, 100))
  d$w <- ifelse(seq_len(300) %% 10 == 0, 9.1, 0.1)
  d$y <- withr::with_seed(2, stats::rbinom(300, 1, 0.4))
  b <- anova_bin(d, "y", "g", weights = "w", plots = FALSE, posthoc = FALSE)
  note <- grep("are unequal and not frequencies", b$notes, value = TRUE)
  expect_length(note, 1L)
  q <- as.numeric(sub(".*test of g rejects roughly ([0-9]+)%.*", "\\1", note)) / 100
  sim <- withr::with_seed(6, mean(replicate(200, {
    yy <- stats::rbinom(300, 1, 0.4)
    suppressWarnings(car::Anova(stats::glm(yy ~ g, stats::binomial(), d, weights = w)))[1, 3] < 0.05
  })))
  expect_lt(abs(q - sim), 0.15)
  d$c <- withr::with_seed(3, stats::rpois(300, 5))
  gq <- anova_glm(d, "c", "g", family = quasipoisson(), weights = "w", plots = FALSE,
                  posthoc = FALSE)
  expect_true(any(grepl("so the F tests in $anova are inflated", gq$notes, fixed = TRUE)))
  d$nt <- withr::with_seed(4, sample(c(5, 10, 200), 300, TRUE))
  d$p <- withr::with_seed(5, stats::rbinom(300, d$nt, 0.3)) / d$nt
  t <- anova_glm(d, "p", "g", family = binomial(), weights = "nt", plots = FALSE, posthoc = FALSE)
  expect_false(any(grepl("sandwich", t$notes, fixed = TRUE)))
})

test_that("an F test is judged on the F scale: same verdict as the LR test on the same deviance", {
  g <- factor(rep(letters[1:6], each = 6))
  d <- withr::with_seed(3, data.frame(g = g, y = 5e7 + round(stats::rnorm(36, 0, 3))))
  fF <- anova_glm(d, "y", "g", family = Gamma("log"), plots = FALSE, posthoc = FALSE)
  fL <- anova_glm(d, "y", "g", family = Gamma("log"), test_statistic = "LR",
                  plots = FALSE, posthoc = FALSE)
  expect_true(is.na(fF$anova$statistic[1]))
  expect_true(is.na(fL$anova$statistic[1]))
})

test_that("exact multiplicative quasi-Poisson responses at small levels are refused", {
  d <- expand.grid(a = factor(1:2), b = factor(1:3), rep = 1:2)
  d$y <- c(1, 2)[d$a] * c(1, 3, 7)[d$b] * 1e-3
  expect_error(anova_glm(d, "y", c("a", "b"), family = quasipoisson(), plots = FALSE),
               "fits response `y` exactly")
  d2 <- expand.grid(a = gl(3, 1), b = gl(4, 1))
  d2$y <- 1e-6 * c(1, 2.5, 7)[d2$a] * c(1, 1.3, 0.4, 3)[d2$b]
  expect_error(anova_glm(d2, "y", c("a", "b"), family = quasipoisson(), plots = FALSE),
               "fits response `y` exactly")
})

test_that("Gaussian integer variation on a level of 1e14 is analysed (centred lm)", {
  d <- data.frame(g = gl(2, 3), y = 1e14 + c(0, 1, 2, 5, 6, 8))
  ref <- stats::anova(stats::lm(I(y - 1e14) ~ g, d))[["F value"]][1]
  expect_equal(anova_glm(d, "y", "g", plots = FALSE)$anova$statistic[1], ref, tolerance = 0.01)
})

test_that("equal group means keep a deviance explained of 0", {
  d <- data.frame(g = rep(c("ctrl", "trt"), each = 5), y = c(10, 14, 9, 15, 12, 11, 13, 12, 16, 8))
  f <- anova_glm(d, "y", "g", family = poisson(), plots = FALSE, posthoc = FALSE)
  m <- stats::glm(y ~ g, stats::poisson(), d)
  expect_equal(f$effect_sizes$estimate[f$effect_sizes$measure == "deviance_explained"],
               1 - stats::deviance(m) / m$null.deviance, tolerance = 1e-8)
})

test_that("weight notes describe the reported test: none beside a robust Wald test", {
  skip_if_not_installed("sandwich")
  aq <- stats::na.omit(airquality)
  aq$Month <- factor(aq$Month); aq$hot <- as.integer(aq$Temp > 80)
  aq$w <- withr::with_seed(5, stats::runif(nrow(aq), 0.5, 2))
  rw <- anova_glm(aq, "hot", "Month", family = binomial(), weights = "w",
                  test_statistic = "Wald", vcov_type = "HC3", plots = FALSE)
  expect_false(any(grepl("are unequal and not frequencies", rw$notes, fixed = TRUE)))
  rb <- anova_bin(aq, "hot", "Month", weights = "w", test_statistic = "Wald",
                  vcov_type = "HC3", plots = FALSE)
  expect_false(any(grepl("are unequal and not frequencies", rb$notes, fixed = TRUE)))
  d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 200)))
  d$w <- withr::with_seed(11, stats::rexp(600) + 0.2)
  d$y <- withr::with_seed(12, stats::rpois(600, 4))
  q <- anova_glm(d, "y", "g", family = quasipoisson(), weights = "w",
                 test_statistic = "Wald", vcov_type = "HC3", plots = FALSE)
  expect_false(any(grepl("are unequal and not frequencies", q$notes, fixed = TRUE)))
  # trial counts that are not whole numbers are still trials on proportions
  es <- esoph
  for (v in c("agegp", "alcgp", "tobgp")) es[[v]] <- factor(es[[v]], ordered = FALSE)
  es$n <- es$ncases + es$ncontrols; es$prop <- es$ncases / es$n; es$wn <- es$n + 0.5
  fe <- anova_glm(es, "prop", c("agegp", "alcgp", "tobgp"), family = binomial(),
                  weights = "wn", plots = FALSE, posthoc = FALSE)
  expect_false(any(grepl("sandwich", fe$notes, fixed = TRUE)))
})

test_that("an exactly additive Gaussian response with thousands of rows is still an exact fit", {
  d <- expand.grid(a = gl(5, 1), b = gl(7, 1), rep = 1:200)
  d$y <- withr::with_seed(3, 1 + round(stats::runif(5), 3)[d$a] / 7 +
                            round(stats::runif(7), 3)[d$b] * 0.37)
  expect_error(anova_glm(d, "y", c("a", "b"), plots = FALSE, posthoc = FALSE),
               "fits response `y` exactly")
})

test_that("barely unequal weights keep the overdispersion note; the weight note names the reported test", {
  skip_if_not_installed("MASS")
  d <- data.frame(g = gl(3, 50))
  d$y <- withr::with_seed(1, MASS::rnegbin(150, 5, 1.5))
  d$w <- withr::with_seed(2, 1 + stats::runif(150, 0, 0.05))
  f <- anova_glm(d, "y", "g", family = poisson(), weights = "w", plots = FALSE)
  expect_true(any(grepl("Pearson dispersion is", f$notes, fixed = TRUE)))
  expect_false(any(grepl("are unequal and not frequencies", f$notes, fixed = TRUE)))
  d2 <- data.frame(g = factor(rep(c("a", "b", "c"), each = 200)))
  d2$w <- withr::with_seed(11, stats::rexp(600) + 0.2)
  d2$c <- withr::with_seed(12, stats::rpois(600, 4))
  q <- anova_glm(d2, "c", "g", family = quasipoisson(), weights = "w",
                 test_statistic = "Wald", plots = FALSE)
  expect_true(any(grepl("so the Wald tests in $anova are inflated", q$notes, fixed = TRUE)))
})

test_that("glm.nb() failing on counts no more variable than Poisson gives the Poisson limit", {
  skip_if_not_installed("MASS")
  # As on ARM macOS, where theta.ml() can return NaN at the theta = Inf boundary.
  d <- data.frame(y = c(3, 3, 3, 5, 5, 5, 8, 8, 8), g = rep(c("a", "b", "c"), each = 3))
  real_fit <- .fit_glm
  local_mocked_bindings(.fit_glm = function(fml, d, family, ...) {
    if (is.null(family)) stop("missing value where TRUE/FALSE needed")
    real_fit(fml, d, family, ...)
  })
  r <- anova_count(d, "y", "g", model = "negbin", plots = FALSE)
  expect_identical(r$model_type, "negbin")
  ref <- stats::anova(stats::glm(y ~ g, stats::poisson(), d), test = "Chisq")
  expect_equal(r$anova$statistic[1], ref$Deviance[2], tolerance = 1e-6)
  expect_true(any(grepl("Poisson", r$notes, fixed = TRUE)))
})

test_that("emmeans gets anovakit's estimability basis, not estimability's svd of a wide R", {
  skip_if_not_installed("estimability")
  cells <- data.frame(A = rep(1:10, each = 2), B = c(rbind(1:10, c(2:10, 1))))
  d <- cells[rep(1:20, each = 4), ]
  d$A <- factor(d$A); d$B <- factor(d$B)
  d$x <- withr::with_seed(15, stats::rnorm(80)); d$y <- d$x + withr::with_seed(16, stats::rnorm(80))
  m <- stats::lm(y ~ x + A * B, data = d)
  expect_gt(ncol(m$qr$qr), nrow(m$qr$qr))   # 101 coefficients, 80 rows, rank 21
  m2 <- .emm_ready(m)
  nb <- estimability::nonest.basis(m2$qr)
  expect_identical(dim(nb), c(101L, 80L))
  # the same non-estimable space as base R's: the null space of X
  X <- stats::model.matrix(m)
  qx <- qr(t(X))
  ref <- qr.Q(qx, complete = TRUE)[, (qx$rank + 1L):ncol(X), drop = FALSE]
  proj <- nb %*% crossprod(nb, ref)
  expect_lt(max(abs(proj - ref)), 1e-8)
  expect_lt(max(abs(X %*% nb)), 1e-8)
  # and the same marginal means as emmeans' own route, where estimability's
  # method works (it fails on this model under R 4.6 with estimability 2.0)
  e2 <- summary(emmeans::emmeans(m2, ~ A:B))
  e1 <- tryCatch(summary(emmeans::emmeans(m, ~ A:B)), error = function(e) NULL)
  if (!is.null(e1)) {
    expect_equal(e2$emmean, e1$emmean, tolerance = 1e-10)
    expect_equal(e2$SE, e1$SE, tolerance = 1e-10)
  }
  skip_if_not(exists("nonest.basis.qr", envir = asNamespace("estimability"),
                     inherits = FALSE))
  # Under R 4.6 estimability 2.0's own method fails on this model ("infinite
  # or missing values in 'x'"); with it failing, anova_ancova() still gives
  # all 20 observed cells.
  local_mocked_bindings(
    nonest.basis.qr = function(x, ...) stop("infinite or missing values in 'x'"),
    .package = "estimability")
  expect_error(estimability::nonest.basis(m$qr), "infinite or missing")
  f <- anova_ancova(d, "y", c("A", "B"), "x", force_interaction = FALSE, plots = FALSE)
  expect_identical(nrow(f$emmeans), 20L, info = paste(f$notes, collapse = " | "))
})

test_that("zero repeated-measures error variance is recognised at the data's rounding level", {
  skip_if_not_installed("afex")
  d <- expand.grid(id = factor(1:5), time = factor(1:3))
  d$y <- as.numeric(d$id) + as.numeric(d$time)
  f <- anova_rm(d, "y", "id", "time", plots = FALSE)
  expect_true(is.na(f$anova$p_value[1]))
  # real variation of 1e-6 at a level of 5 is a test, not zero variance
  d$y2 <- d$y + withr::with_seed(2, stats::rnorm(15, sd = 1e-6))
  f2 <- anova_rm(d, "y2", "id", "time", plots = FALSE)
  expect_false(is.na(f2$anova$p_value[1]))
})
