# Regression tests for the third review's anova_welch() / anova_kw() findings
# (R3-25, R3-29, R3-35, R3-36, R3-46). Each is checked against an independent
# reference: base R (t.test, oneway.test, table), a hand computation from the
# published formula, or the noncentral t distribution function computed by a
# different numerical integral from the one the package uses.

# P(T <= t) for a noncentral t, integrating pnorm(t sqrt(V/df) - ncp) over
# V ~ chi-square(df) on the probability scale (V = qchisq(u, df)) with
# integrate(): neither the quadrature nor the conditioning of the package's
# .pnt().
r3wk_pnt <- function(t, df, ncp) {
  br <- c(0, 1e-12, 1e-8, 1e-5, 1e-3, 0.02, 0.1, 0.3, 0.5, 0.7, 0.9, 0.98,
          1 - 1e-3, 1 - 1e-5, 1 - 1e-8, 1)
  sum(vapply(seq_len(length(br) - 1L), function(k) {
    stats::integrate(function(u) {
      stats::pnorm(t * sqrt(stats::qchisq(u, df) / df) - ncp)
    }, br[k], br[k + 1L], rel.tol = 1e-10, abs.tol = 1e-14,
    subdivisions = 2000L)$value
  }, numeric(1)))
}

r3wk_two <- function(x, z) {
  data.frame(g = rep(c("a", "b"), c(length(x), length(z))), y = c(x, z))
}

test_that("R3-25: the Welch article describes the notes the function writes", {
  rmd <- testthat::test_path("..", "..", "vignettes", "articles", "welch.Rmd")
  skip_if_not(file.exists(rmd), "the article sources are not available")
  txt <- gsub("\\s+", " ", paste(readLines(rmd, warn = FALSE), collapse = " "))
  # The false claim: a note only "if you ask for a low conf_level".
  expect_false(grepl("note if you ask for a low `conf_level`", txt, fixed = TRUE))
  expect_true(grepl("fewer than 2 degrees of freedom", txt, fixed = TRUE))

  # What the article now says, at the default conf_level with a group of 30
  # beside one of two: the standardiser's df (hand computed) is below 2 and a
  # note says so.
  x <- c(0, 1); z <- 30 + seq(-0.1, 0.1, length.out = 30)
  f <- anova_welch(r3wk_two(x, z), "y", "g", plots = FALSE)
  expect_equal(f$conf_level, 0.95)
  nu <- (var(x) + var(z))^2 / (var(x)^2 / 1 + var(z)^2 / 29)
  expect_lt(nu, 2)
  expect_length(grep("fewer than 2 degrees of freedom", f$notes), 1L)
  expect_false(any(grepl("outside its", f$notes)))
})

test_that("R3-29: Welch and KW Shapiro-Wilk skip notes speak of observations, not residuals", {
  mw <- anova_welch(mtcars, "mpg", c("cyl", "am"), plots = FALSE)
  nt <- mw$assumptions$normality
  # Reference: the cell sizes in base R. Eight cylinders with a manual
  # gearbox has two cars.
  tab <- table(mtcars$cyl, mtcars$am)
  expect_equal(tab["8", "1"], 2L, ignore_attr = TRUE)
  row <- nt[nt$group == "8 : 1", ]
  expect_equal(row$n, 2)
  expect_equal(row$note,
               "Shapiro-Wilk skipped: 2 usable observation(s), at least 3 are required.")
  expect_false(any(grepl("residual", c(nt$note, mw$notes))))
  # Welch fits no model: there are no residuals to speak of.
  expect_s3_class(mw$model, "htest")

  kw <- anova_kw(mtcars, "mpg", c("cyl", "am"), plots = FALSE,
                 diagnostics = TRUE)
  expect_true(any(grepl("2 usable observation(s)", kw$notes, fixed = TRUE)))
  expect_false(any(grepl("residual", kw$assumptions$normality$note)))

  # Callers that do test residuals keep the word.
  nm <- anovakit:::.check_normality(c(1, 2))
  expect_match(nm$note, "2 usable residual(s)", fixed = TRUE)
  big <- anovakit:::.check_normality(seq_len(5001), unit = "observation")
  expect_match(big$note, "5001 observations exceeds", fixed = TRUE)
})

test_that("R3-35: an omnibus Welch F made of rounding noise is NA with a note", {
  d <- data.frame(g = rep(c("a", "b"), each = 3),
                  y = 1e6 + c(0, 1e-9, 2e-9, 5, 5 + 1e-9, 5 + 2e-9))
  # Reference: base R refuses these data, and with two groups F = t^2.
  expect_error(stats::t.test(y ~ g, d), "essentially constant")
  f <- anova_welch(d, "y", "g", plots = FALSE)
  expect_true(is.na(f$anova$statistic))
  expect_true(is.na(f$anova$den_df))
  expect_true(is.na(f$anova$p_value))
  expect_equal(f$anova$num_df, 1)
  expect_true(is.na(f$model$statistic))
  expect_true(is.na(f$model$p.value))
  expect_true(is.na(f$posthoc$statistic))
  note <- grep("omnibus test is left as NA", f$notes, value = TRUE)
  expect_length(note, 1L)
  expect_match(note, "a vs b")
  expect_match(note, "100% of Welch's F")
  expect_match(note, "Subtracting a constant")

  # Shifted by the constant the advice names, the data give oneway.test()'s F.
  d0 <- d; d0$y <- d0$y - 1e6
  f0 <- anova_welch(d0, "y", "g", plots = FALSE)
  ow <- stats::oneway.test(y ~ g, d0)
  expect_equal(f0$anova$statistic, unname(ow$statistic), tolerance = 1e-10)
  expect_false(any(grepl("floating-point resolution", f0$notes)))

  # Three groups, posthoc = FALSE: the noisy pair dominates the F all the same.
  d3 <- data.frame(g = rep(c("a", "b", "c"), each = 3),
                   y = 1e6 + c(0, 1e-9, 2e-9, 5, 5 + 1e-9, 5 + 2e-9, 3, 4, 6))
  f3 <- anova_welch(d3, "y", "g", plots = FALSE, posthoc = FALSE)
  expect_gt(unname(stats::oneway.test(y ~ g, d3)$statistic), 1e15)
  expect_true(is.na(f3$anova$statistic))
  expect_length(grep("omnibus test is left as NA", f3$notes), 1L)

  # A noisy pair that is a small part of F leaves the test as oneway.test()
  # gives it, and the note states the share, computed here by hand from the
  # definition of Welch's numerator, sum(w (m - mbar)^2).
  d4 <- data.frame(g = rep(c("a", "b", "c"), each = 4),
                   y = c(1e9 + c(0, 1, 2, 3) * 1e-6, 1e9 + c(1, 2, 3, 4) * 1e-6,
                         1e9 + c(5, 6, 7, 8)))
  f4 <- anova_welch(d4, "y", "g", plots = FALSE)
  expect_equal(f4$anova$statistic,
               unname(stats::oneway.test(y ~ g, d4)$statistic), tolerance = 1e-10)
  m <- tapply(d4$y, d4$g, mean); w <- 4 / tapply(d4$y, d4$g, var)
  num <- sum(w * (m - sum(w * m) / sum(w))^2)
  share <- w[["a"]] * w[["b"]] * (m[["a"]] - m[["b"]])^2 / sum(w) / num
  expect_lt(share, 0.05)
  expect_true(any(grepl(sprintf("a vs b, .* only %s%% of Welch's F",
                                format(signif(100 * share, 2))), f4$notes)))
})

test_that("R3-36: well-separated groups get exact intervals without integrate()'s cost", {
  # Exactness against the independent integral, in both regimes of .pnt():
  # over Z (|t| >= sqrt(2 df) or df < 100) and over U (large df).
  pnt <- anovakit:::.pnt
  for (a in list(c(80, 5, 70), c(80, 5, 110), c(500, 6, 300), c(-60, 3, -40),
                 c(200, 1.2, 150), c(40, 5000, 42), c(300, 4000, 310),
                 c(-200, 1e5, -201), c(5, 2, 35))) {
    expect_equal(pnt(a[1], a[2], a[3]), r3wk_pnt(a[1], a[2], a[3]),
                 tolerance = 1e-9, info = paste(a, collapse = ", "))
  }

  # The interval inverts that distribution function.
  set.seed(1)
  x <- 100 + stats::rnorm(5); z <- stats::rnorm(4)
  fit <- anova_welch(r3wk_two(x, z), "y", "g", plots = FALSE,
                     hedges_correction = FALSE)
  v1 <- var(x); v2 <- var(z)
  se <- sqrt(v1 / 5 + v2 / 4); s <- sqrt((v1 + v2) / 2)
  t <- (mean(x) - mean(z)) / se
  df <- se^4 / ((v1 / 5)^2 / 4 + (v2 / 4)^2 / 3)
  expect_gt(abs(t), 30)
  ncp <- unlist(fit$effect_sizes[, c("conf_low", "conf_high")]) * s / se
  expect_equal(r3wk_pnt(t, df, ncp[[1]]), 0.975, tolerance = 1e-7)
  expect_equal(r3wk_pnt(t, df, ncp[[2]]), 0.025, tolerance = 1e-6)

  # 150 well-separated groups (98% of pairs past |t| = 30) took 105 s with an
  # integrate() call per piece; about 6 s now, against about 3 s for null data.
  skip_on_cran()
  k <- 150L
  big <- withr::with_seed(1, data.frame(
    g = rep(sprintf("g%03d", seq_len(k)), each = 5),
    y = rep(seq_len(k), each = 5) * 10 + stats::rnorm(5 * k)))
  el <- system.time(f <- anova_welch(big, "y", "g", plots = FALSE))[["elapsed"]]
  expect_gt(mean(abs(f$posthoc$statistic) > 30), 0.9)
  expect_lt(el, 25)
})

test_that("R3-46: below 2 standardiser df the bias correction is not applied, with a note", {
  x <- c(0, 5); z <- 30 + seq(-0.5, 0.5, length.out = 30)
  fit <- anova_welch(r3wk_two(x, z), "y", "g", plots = FALSE)
  # Hand computation from the formulas in ?anova_welch.
  v1 <- var(x); v2 <- var(z)
  nu <- (v1 + v2)^2 / (v1^2 / 1 + v2^2 / 29)
  d <- (mean(x) - mean(z)) / sqrt((v1 + v2) / 2)
  J <- gamma(nu / 2) / (sqrt(nu / 2) * gamma((nu - 1) / 2))
  expect_lt(nu, 1.02)
  expect_lt(abs(J * d), 0.2)                  # the collapsed, "negligible" g
  expect_equal(fit$effect_sizes$hedges_g, d, tolerance = 1e-10)
  expect_equal(fit$effect_sizes$magnitude, "large")
  # The interval contains the estimate, so no outside-the-interval note.
  expect_lt(fit$effect_sizes$conf_low, d)
  expect_gt(fit$effect_sizes$conf_high, d)
  expect_false(any(grepl("outside its", fit$notes)))
  note <- grep("fewer than 2 degrees of freedom", fit$notes, value = TRUE)
  expect_length(note, 1L)
  expect_match(note, "a vs b")
  expect_match(note, "Cohen's d")
  # Without the correction nothing changes and nothing is said about it.
  fd <- anova_welch(r3wk_two(x, z), "y", "g", plots = FALSE,
                    hedges_correction = FALSE)
  expect_equal(fd$effect_sizes$cohens_d, d, tolerance = 1e-10)
  expect_false(any(grepl("fewer than 2 degrees of freedom", fd$notes)))

  # At and above nu = 2 the correction is applied as published.
  x2 <- c(0, 1, 2, 3, 4); z2 <- c(10, 11, 12, 13, 14)
  f2 <- anova_welch(r3wk_two(x2, z2), "y", "g", plots = FALSE)
  nu2 <- 8
  expect_equal((var(x2) + var(z2))^2 / (var(x2)^2 / 4 + var(z2)^2 / 4), nu2)
  J2 <- gamma(nu2 / 2) / (sqrt(nu2 / 2) * gamma((nu2 - 1) / 2))
  expect_equal(f2$effect_sizes$hedges_g,
               J2 * (mean(x2) - mean(z2)) / sqrt((var(x2) + var(z2)) / 2),
               tolerance = 1e-10)
  expect_false(any(grepl("fewer than 2 degrees of freedom", f2$notes)))
})
