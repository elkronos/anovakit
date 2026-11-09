# Detailed regression tests: run locally and on CI (NOT_CRAN=true), skipped
# on CRAN to keep its check time short. The core suite runs everywhere.
skip_on_cran()

# Regression tests for the second review's anova_welch() / anova_kw()
# findings (R2-WELCHKW-*). Each is checked against an independent reference:
# base R (t.test, oneway.test, kruskal.test), a hand computation from the
# published formula, or the noncentral t distribution function computed by a
# different numerical integral from the one the package uses.

# P(T <= t) for a noncentral t, integrating pnorm(t sqrt(V/df) - ncp) over
# V ~ chi-square(df) on the probability scale (V = qchisq(u, df)), which is
# exact where stats::pt() is (checked below) and does not share the package's
# conditioning on the normal variable.
r2wk_pnt <- function(t, df, ncp) {
  br <- c(0, 1e-12, 1e-8, 1e-5, 1e-3, 0.02, 0.1, 0.3, 0.5, 0.7, 0.9, 0.98,
          1 - 1e-3, 1 - 1e-5, 1 - 1e-8, 1)
  sum(vapply(seq_len(length(br) - 1L), function(k) {
    stats::integrate(function(u) {
      stats::pnorm(t * sqrt(stats::qchisq(u, df) / df) - ncp)
    }, br[k], br[k + 1L], rel.tol = 1e-10, abs.tol = 1e-14,
    subdivisions = 2000L)$value
  }, numeric(1)))
}

# The Welch noncentral-t interval for delta*, by the reference above.
r2wk_ref_ci <- function(x, z, cl = 0.95) {
  v1 <- var(x); v2 <- var(z); n1 <- length(x); n2 <- length(z)
  se <- sqrt(v1 / n1 + v2 / n2); s <- sqrt((v1 + v2) / 2)
  t <- (mean(x) - mean(z)) / se
  df <- se^4 / ((v1 / n1)^2 / (n1 - 1) + (v2 / n2)^2 / (n2 - 1))
  a <- (1 - cl) / 2
  r <- function(p) {
    stats::uniroot(function(n) r2wk_pnt(t, df, n) - p,
                   sort(c(0.001 * t, 3 * t)), tol = 1e-10)$root
  }
  c(r(1 - a), r(a)) * se / s
}

r2wk_two <- function(x, z) {
  data.frame(g = rep(c("a", "b"), c(length(x), length(z))), y = c(x, z))
}

test_that("R2-WELCHKW-01: effect-size intervals stay exact past |ncp| = 37.62", {
  # The reference agrees with pt() where pt() is exact.
  expect_equal(r2wk_pnt(3, 4, 2), stats::pt(3, 4, 2), tolerance = 1e-8)
  expect_equal(r2wk_pnt(-30, 4, -25), stats::pt(-30, 4, -25), tolerance = 1e-8)

  cases <- list(
    triplicates = list(c(10, 10.2, 10.1), c(20, 20.2, 20.1)),
    pairs_500   = list(c(0, 1), c(500, 502)),
    pairs_5000  = list(c(0, 1), c(5000, 5002)))
  for (nm in names(cases)) {
    x <- cases[[nm]][[1]]; z <- cases[[nm]][[2]]
    fit <- anova_welch(r2wk_two(x, z), "y", "g", plots = FALSE,
                       hedges_correction = FALSE)
    got <- unlist(fit$effect_sizes[1, c("conf_low", "conf_high")],
                  use.names = FALSE)
    expect_equal(got, r2wk_ref_ci(x, z), tolerance = 1e-6, label = nm)
  }
  # Values recorded in the finding (exact inversion): the upper bound no
  # longer sticks at -37.622 when the difference grows ten-fold.
  f500 <- anova_welch(r2wk_two(c(0, 1), c(500, 502)), "y", "g",
                      plots = FALSE, hedges_correction = FALSE)
  f5000 <- anova_welch(r2wk_two(c(0, 1), c(5000, 5002)), "y", "g",
                       plots = FALSE, hedges_correction = FALSE)
  expect_equal(f500$effect_sizes$conf_high, -40.088, tolerance = 1e-4)
  expect_equal(f5000$effect_sizes$conf_high, -400.579, tolerance = 1e-5)

  # Monte Carlo check of one bound: P(T <= t | ncp = bound) is the tail level.
  x <- c(10, 10.2, 10.1); z <- c(20, 20.2, 20.1)
  v1 <- var(x); v2 <- var(z); se <- sqrt(v1 / 3 + v2 / 3)
  t <- (mean(x) - mean(z)) / se; df <- se^4 / ((v1 / 3)^2 / 2 + (v2 / 3)^2 / 2)
  fit <- anova_welch(r2wk_two(x, z), "y", "g", plots = FALSE,
                     hedges_correction = FALSE)
  ncp_low <- fit$effect_sizes$conf_low * sqrt((v1 + v2) / 2) / se
  sim <- withr::with_seed(1, {
    n <- 4e5
    mean((stats::rnorm(n) + ncp_low) / sqrt(stats::rchisq(n, df) / df) <= t)
  })
  expect_equal(sim, 0.975, tolerance = 0.002)
})

test_that("R2-WELCHKW-02: a g collapsed by J(nu) near nu = 1 is not labelled 'negligible'", {
  d <- data.frame(g = rep(c("a", "b"), c(2, 30)),
                  y = c(0, 1, 30 + seq(-0.03, 0.03, length.out = 30)))
  fit <- anova_welch(d, "y", "g", plots = FALSE)
  x <- d$y[d$g == "a"]; z <- d$y[d$g == "b"]

  # Hand computation from the published formula (Delacre et al., 2021).
  v1 <- var(x); v2 <- var(z)
  nu <- (v1 + v2)^2 / (v1^2 / 1 + v2^2 / 29)
  J <- gamma(nu / 2) / (sqrt(nu / 2) * gamma((nu - 1) / 2))
  g <- J * (mean(x) - mean(z)) / sqrt((v1 + v2) / 2)
  expect_lt(nu, 1.01)
  expect_lt(abs(g), 0.2)                       # what would read "negligible"
  # Since R3-46 the correction is not applied below nu = 2: the reported
  # estimate is the uncorrected one, with a note.
  d_hand <- (mean(x) - mean(z)) / sqrt((v1 + v2) / 2)
  expect_equal(fit$effect_sizes$hedges_g, d_hand, tolerance = 1e-10)

  # Welch p (t.test) and the interval both rule out anything near zero.
  expect_lt(stats::t.test(x, z)$p.value, 0.05)
  expect_lt(fit$effect_sizes$conf_high, -0.8)
  expect_equal(fit$effect_sizes$magnitude, "large")

  expect_false(any(grepl("outside its 95% interval", fit$notes)))
  note <- grep("fewer than 2 degrees of freedom", fit$notes, value = TRUE)
  expect_length(note, 1L)
  expect_false(grepl("only with very small groups and a low conf_level", note))
})

test_that("R2-WELCHKW-03: a pair at floating-point resolution is NA with a note, not a crash", {
  d <- data.frame(g = rep(c("a", "b", "c"), each = 4),
                  y = c(1e9 + c(0, 1, 2, 3) * 1e-6, 1e9 + c(1, 2, 3, 4) * 1e-6,
                        1e9 + c(5, 6, 7, 8)))
  ya <- d$y[d$g == "a"]; yb <- d$y[d$g == "b"]; yc <- d$y[d$g == "c"]
  # Base R refuses this pair.
  expect_error(stats::t.test(ya, yb), "essentially constant")

  fit <- expect_no_error(anova_welch(d, "y", "g", plots = FALSE))
  ph <- fit$posthoc
  ab <- ph$group1 == "a" & ph$group2 == "b"
  expect_true(all(is.na(unlist(ph[ab, c("conf_low", "conf_high", "statistic",
                                        "df", "p_value", "p_adjusted")]))))
  expect_true(is.na(fit$effect_sizes[ab, "hedges_g"]))
  expect_true(any(grepl("a vs b, the spread within the groups is at the floating-point resolution",
                        fit$notes, fixed = TRUE)))

  # The resolvable pairs match t.test(); the omnibus matches oneway.test().
  for (pr in list(c("a", "c"), c("b", "c"))) {
    tt <- stats::t.test(d$y[d$g == pr[1]], d$y[d$g == pr[2]])
    row <- ph$group1 == pr[1] & ph$group2 == pr[2]
    expect_equal(ph$statistic[row], unname(tt$statistic), tolerance = 1e-8)
    expect_equal(ph$df[row], unname(tt$parameter), tolerance = 1e-8)
    expect_equal(ph$p_value[row], tt$p.value, tolerance = 1e-8)
  }
  ow <- stats::oneway.test(y ~ g, d)
  expect_equal(fit$anova$statistic, unname(ow$statistic), tolerance = 1e-10)
})

test_that("R2-WELCHKW-04: with hedges_correction = FALSE the note names Cohen's d", {
  d2 <- data.frame(g = rep(c("a", "b"), c(2, 2)), y = c(0, 1, 10, 10.5))
  fit <- anova_welch(d2, "y", "g", plots = FALSE, conf_level = 0.01,
                     hedges_correction = FALSE)
  x <- c(0, 1); z <- c(10, 10.5)
  d_hand <- (mean(x) - mean(z)) / sqrt((var(x) + var(z)) / 2)
  expect_equal(fit$effect_sizes$cohens_d, d_hand, tolerance = 1e-10)
  expect_equal(unlist(fit$effect_sizes[1, c("conf_low", "conf_high")],
                      use.names = FALSE),
               r2wk_ref_ci(x, z, 0.01), tolerance = 1e-6)
  expect_true(d_hand < fit$effect_sizes$conf_low)

  note <- grep("outside its 1% interval", fit$notes, value = TRUE)
  expect_length(note, 1L)
  expect_match(note, "Cohen's d")
  expect_false(grepl("Hedges", note))
  expect_false(grepl("bias-corrected|correction shrinks", note))
})

test_that("R2-WELCHKW-05: variance under- or overflow is not reported as constant data", {
  base <- data.frame(g = rep(c("a", "b"), each = 4), y = c(1:4, 3:6))
  for (scale in c(1e-200, 1e200)) {
    d <- base; d$y <- d$y * scale
    expect_gt(length(unique(d$y)), 2L)
    expect_error(anova_welch(d, "y", "g", plots = FALSE),
                 "cannot be represented in double precision.*Rescale")
    err <- tryCatch(anova_welch(d, "y", "g", plots = FALSE),
                    error = conditionMessage)
    expect_false(grepl("takes the same value", err))
  }
  # Rescaled, it runs and matches oneway.test (F = 4.8 on 1 and 6 df).
  fit <- anova_welch(base, "y", "g", plots = FALSE)
  ow <- stats::oneway.test(y ~ g, base)
  expect_equal(fit$anova$statistic, unname(ow$statistic), tolerance = 1e-10)
  expect_equal(fit$anova$statistic, 4.8, tolerance = 1e-10)
  # Truly constant groups keep the original message.
  const <- data.frame(g = rep(c("a", "b"), each = 3), y = c(1, 1, 1, 2, 3, 4))
  expect_error(anova_welch(const, "y", "g", plots = FALSE),
               "takes the same value")
})

test_that("R2-WELCHKW-06: a group lost to missing values is named in $notes", {
  d <- withr::with_seed(1, data.frame(
    g = factor(rep(c("ctrl", "low", "high"), each = 5),
               levels = c("ctrl", "low", "high")),
    y = stats::rnorm(15)))
  d$y[d$g == "high"] <- NA
  note <- "Level \"high\" of `g` had no complete rows, so it was dropped from the analysis."

  fw <- anova_welch(d, "y", "g", plots = FALSE)
  expect_true(note %in% fw$notes)
  ow <- stats::oneway.test(y ~ g, d)               # base R drops it silently
  expect_equal(fw$anova$statistic, unname(ow$statistic), tolerance = 1e-10)

  fk <- anova_kw(d, "y", "g", plots = FALSE)
  expect_true(note %in% fk$notes)
  kt <- stats::kruskal.test(y ~ g, d)
  expect_equal(fk$anova$statistic, unname(kt$statistic), tolerance = 1e-10)

  # No such note when every level keeps some rows.
  d2 <- d; d2$y[d2$g == "high"][1:2] <- c(0.5, -0.2)
  expect_false(any(grepl("had no complete rows",
                         anova_welch(d2, "y", "g", plots = FALSE)$notes)))
})

test_that("R2-WELCHKW-07: many-group pairwise comparisons are vectorised and match t.test", {
  k <- 30L
  big <- withr::with_seed(1, data.frame(
    g = rep(sprintf("g%02d", seq_len(k)), each = 5),
    y = stats::rnorm(5 * k, sd = rep(seq(0.5, 3, length.out = k), each = 5))))
  fit <- anova_welch(big, "y", "g", plots = FALSE)
  ph <- fit$posthoc
  expect_equal(nrow(ph), choose(k, 2))
  ref <- t(vapply(seq_len(nrow(ph)), function(r) {
    tt <- stats::t.test(big$y[big$g == ph$group1[r]],
                        big$y[big$g == ph$group2[r]])
    c(unname(tt$estimate[1] - tt$estimate[2]), tt$conf.int,
      unname(tt$statistic), unname(tt$parameter), tt$p.value)
  }, numeric(6)))
  got <- as.matrix(ph[, c("difference", "conf_low", "conf_high",
                          "statistic", "df", "p_value")])
  expect_equal(unname(got), ref, tolerance = 1e-10)
  expect_equal(ph$p_adjusted, stats::p.adjust(ref[, 6], "holm"),
               tolerance = 1e-10)

  # 150 groups (11,175 pairs) took about 14 s before; allow generous slack.
  skip_on_cran()
  k <- 150L
  huge <- withr::with_seed(1, data.frame(
    g = rep(sprintf("g%03d", seq_len(k)), each = 5),
    y = stats::rnorm(5 * k)))
  el <- system.time(anova_welch(huge, "y", "g", plots = FALSE))[["elapsed"]]
  expect_lt(el, 9)
})
