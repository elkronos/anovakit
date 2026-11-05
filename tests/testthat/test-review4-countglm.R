# Round-4 review fixes: count and GLM area ------------------------------------

test_that("glm-F1: a saturated GLM that estimates its dispersion reports NA inference, not 0 or NaN", {
  sat <- data.frame(y = c(1, 2, 4), g = c("a", "b", "c"))
  # reference: base R cannot form the F test either (no residual df)
  m <- stats::lm(y ~ g, sat)
  expect_identical(stats::df.residual(m), 0L)
  expect_false(is.finite(suppressWarnings(stats::anova(m))[["F value"]][1]))
  for (fam in list(gaussian(), Gamma("log"), quasipoisson(), inverse.gaussian())) {
    expect_silent(f <- anova_glm(sat, "y", "g", family = fam, plots = FALSE))
    a <- f$anova
    expect_true(is.na(a$statistic[1]) && !is.nan(a$statistic[1]))
    expect_true(is.na(a$p_value[1]) && !is.nan(a$p_value[1]))
    # the sums of squares are still those of the fit
    ref <- suppressWarnings(stats::glm(y ~ g, fam, transform(sat, g = factor(g))))
    expect_equal(a$sum_sq[1], ref$null.deviance - ref$deviance, tolerance = 1e-8)
    num <- unlist(Filter(is.double, c(f$posthoc, f$emmeans,
                                      f$assumptions$coefficients)))
    expect_false(any(is.nan(num)))
    expect_true(all(is.na(f$posthoc$p_value)))
    expect_true(any(grepl("saturated", f$notes)))
    expect_false(any(grepl("NaNs produced", f$notes)))
  }
  # a fixed-dispersion family keeps its valid tests
  fp <- anova_glm(sat, "y", "g", family = poisson(), plots = FALSE)
  rp <- car::Anova(stats::glm(y ~ g, poisson, transform(sat, g = factor(g))))
  expect_equal(fp$anova$statistic[1], rp[["LR Chisq"]][1], tolerance = 1e-8)
  # anova_count's quasi-Poisson fit is reported the same way
  fc <- anova_count(sat, "y", "g", model = "quasipoisson", plots = FALSE)
  expect_true(is.na(fc$anova$statistic[1]) && !is.nan(fc$anova$statistic[1]))
  expect_false(any(is.nan(fc$posthoc$se)))
  expect_false(any(grepl("NaNs produced", fc$notes)))
})

test_that("bincount-F1: anova_glm's overdispersion note gives the size for the largest term, as anova_count does", {
  d <- withr::with_seed(2, {
    dd <- data.frame(g = rep(c("a", "b", "c"), each = 12))
    dd$cnt <- stats::rnbinom(36, mu = exp(1.5 + rep(c(0, .3, .6), each = 12)), size = 2)
    dd
  })
  m <- stats::glm(cnt ~ g, poisson, d)
  phi <- sum(stats::residuals(m, type = "pearson")^2) / stats::df.residual(m)
  size <- function(k) stats::pchisq(stats::qchisq(0.95, k) / phi, k, lower.tail = FALSE)
  expect_equal(round(100 * c(size(1), size(2))), c(28, 40))
  clause <- sprintf("on 1 degree of freedom rejects roughly %.0f%% of true null hypotheses, and one on 2 (the largest term here) roughly %.0f%%",
                    100 * size(1), 100 * size(2))
  g <- anova_glm(d, "cnt", "g", family = poisson(), plots = FALSE)
  expect_true(any(grepl(clause, g$notes, fixed = TRUE)))
  skip_if_not_installed("sandwich")
  gr <- anova_glm(d, "cnt", "g", family = poisson(), vcov_type = "HC3", plots = FALSE)
  expect_true(any(grepl(clause, gr$notes, fixed = TRUE)))
  k <- anova_count(d, "cnt", "g", model = "poisson", plots = FALSE)
  expect_true(any(grepl(clause, k$notes, fixed = TRUE)))
})

test_that("count-glm-merge-F1: a constant count over varying exposures is analysed, not refused", {
  d <- withr::with_seed(2, data.frame(g = gl(3, 8), y = 5L,
                                      expo = round(stats::runif(24, 0.5, 4), 2)))
  m <- stats::glm(y ~ g + offset(log(expo)), quasipoisson, d)
  ref <- car::Anova(m, test.statistic = "F")
  expect_gt(m$deviance, 1)
  f <- anova_count(d, "y", "g", offset = "expo", model = "quasipoisson", plots = FALSE)
  expect_equal(f$anova$statistic[1], ref[["F value"]][1], tolerance = 1e-8)
  expect_equal(f$anova$p_value[1], ref[["Pr(>F)"]][1], tolerance = 1e-8)
  expect_equal(f$anova$statistic[1], 0.1217, tolerance = 1e-3)
  fp <- anova_count(d, "y", "g", offset = "expo", model = "poisson",
                    test_statistic = "F", plots = FALSE)
  rp <- car::Anova(stats::glm(y ~ g + offset(log(expo)), poisson, d),
                   test.statistic = "F")
  expect_equal(fp$anova$statistic[1], rp[["F value"]][1], tolerance = 1e-8)
  # a constant count over a constant exposure is still refused
  d$expo <- 2
  expect_error(anova_count(d, "y", "g", offset = "expo", model = "poisson",
                           test_statistic = "F", plots = FALSE),
               "takes the same value \\(5\\)")
})

test_that("count-glm-merge-F2: a response constant within every group at a boundary is an exact fit", {
  d <- data.frame(g = gl(2, 10), y = rep(0:1, each = 10))
  # reference: glm's fit matches the data to rounding error, and its F is huge
  m <- suppressWarnings(stats::glm(y ~ g, binomial, d))
  expect_lt(max(abs(m$y - fitted(m))), 1e-9)
  expect_gt(car::Anova(m, test.statistic = "F")[["F value"]][1], 1e10)
  expect_error(anova_glm(d, "y", "g", family = "binomial", test_statistic = "F",
                         plots = FALSE), "fits response `y` exactly")
  expect_error(anova_glm(d, "y", "g", family = "quasibinomial", plots = FALSE),
               "fits response `y` exactly")
  dc <- data.frame(g = gl(2, 10), y = rep(c(0L, 3L), each = 10))
  mc <- stats::glm(y ~ g, quasipoisson, dc)
  expect_lt(max(abs(mc$y - fitted(mc))), 1e-9)
  expect_error(anova_count(dc, "y", "g", model = "quasipoisson", plots = FALSE),
               "fits response `y` exactly")
  expect_error(anova_glm(dc, "y", "g", family = quasipoisson(), plots = FALSE),
               "fits response `y` exactly")
  # one group at the boundary, the other varying: a real F test
  dv <- data.frame(g = gl(2, 10), y = c(rep(0L, 10), 1:10))
  fv <- suppressWarnings(anova_count(dv, "y", "g", model = "quasipoisson", plots = FALSE))
  rv <- car::Anova(suppressWarnings(stats::glm(y ~ g, quasipoisson, dv)), test.statistic = "F")
  expect_equal(fv$anova$statistic[1], rv[["F value"]][1], tolerance = 1e-6)
  expect_lt(fv$anova$statistic[1], 1e6)
})

test_that("count-glm-merge-F3: robust F notes agree, do not repeat, and name $effect_sizes in anova_count", {
  skip_if_not_installed("sandwich")
  k <- anova_count(InsectSprays, "count", "spray", model = "poisson",
                   test_statistic = "F", vcov_type = "HC3", plots = FALSE)
  m <- stats::glm(count ~ spray, poisson, InsectSprays)
  expect_equal(k$anova$statistic[1], car::Anova(m, test.statistic = "F")[["F value"]][1],
               tolerance = 1e-8)
  # one note says which covariance the F test uses, and none contradicts it
  expect_length(grep("not the robust (HC3) covariance", k$notes, fixed = TRUE), 1L)
  expect_false(any(grepl("rests on the model-based variance", k$notes)))
  expect_false(any(grepl("coefficient table", k$notes)))
  expect_length(grep("$effect_sizes, $emmeans and $posthoc", k$notes, fixed = TRUE), 1L)
  g <- anova_glm(InsectSprays, "count", "spray", family = poisson(),
                 test_statistic = "F", vcov_type = "HC3", plots = FALSE)
  expect_length(grep("not the robust (HC3) covariance", g$notes, fixed = TRUE), 1L)
  expect_false(any(grepl("rests on the model-based variance", g$notes)))
  expect_length(grep("use the robust (HC3) covariance", g$notes, fixed = TRUE), 0L)
  # a quasi-Poisson F test with a robust covariance keeps the omnibus note
  q <- anova_count(InsectSprays, "count", "spray", model = "quasipoisson",
                   vcov_type = "HC3", plots = FALSE)
  qn <- grep("omnibus F test", q$notes, value = TRUE)
  expect_length(qn, 1L)
  expect_match(qn, "only $effect_sizes, the marginal means", fixed = TRUE)
  # the likelihood-ratio + robust notes name $effect_sizes too
  lr <- anova_count(InsectSprays, "count", "spray", model = "poisson",
                    vcov_type = "HC3", plots = FALSE)
  expect_false(any(grepl("coefficient table", lr$notes)))
  expect_true(any(grepl("used only by $effect_sizes", lr$notes, fixed = TRUE)))
  expect_true(any(grepl("$effect_sizes, the marginal means and the comparisons use the robust",
                        lr$notes, fixed = TRUE)))
})

test_that("bin-merge-F2: anova_bin's profile intervals do not fall back on the number of rows", {
  d <- withr::with_seed(2, {
    n <- 25000
    dd <- data.frame(g = factor(sample(letters[1:10], n, TRUE)))
    dd$y <- stats::rbinom(n, 1, 0.3)
    dd
  })
  # 10 coefficients x 25,000 rows exceeds the 200,000 the vignette once named
  f <- anova_bin(d, "y", "g", plots = FALSE, posthoc = FALSE)
  expect_true(all(f$effect_sizes$ci_method == "profile"))
  expect_false(any(grepl("Wald intervals", f$notes)))
  # reference: profile intervals of the same model on the data aggregated by cell
  agg <- stats::aggregate(cbind(s = y, n = 1) ~ g, d, sum)
  ma <- stats::glm(cbind(s, n - s) ~ g, binomial, agg)
  ci <- exp(suppressMessages(stats::confint(ma, "gb")))
  expect_equal(f$effect_sizes$conf_low[f$effect_sizes$term == "gb"], unname(ci[1]),
               tolerance = 1e-4)
  expect_equal(f$effect_sizes$conf_high[f$effect_sizes$term == "gb"], unname(ci[2]),
               tolerance = 1e-4)
  rmd <- test_path("..", "..", "vignettes", "articles", "binary.Rmd")
  skip_if_not(file.exists(rmd))
  txt <- paste(readLines(rmd, warn = FALSE), collapse = " ")
  expect_false(grepl("200,000", txt, fixed = TRUE))
  expect_true(grepl("number of rows does not add to the cost", txt, fixed = TRUE))
})

test_that("docs-consistency-F1: under separation the default still falls back to Wald, except for the diverged odds ratio", {
  d <- withr::with_seed(4, {
    dd <- data.frame(g = factor(rep(sprintf("L%02d", 1:55), each = 20)))
    dd$y <- stats::rbinom(nrow(dd), 1, 0.4)
    dd$y[dd$g == "L55"] <- 0
    dd
  })
  f <- anova_bin(d, "y", "g", plots = FALSE, posthoc = FALSE)
  expect_equal(as.vector(table(f$effect_sizes$ci_method)[c("profile", "wald")]), c(1L, 53L))
  expect_identical(f$effect_sizes$term[f$effect_sizes$ci_method == "profile"], "gL55")
  # reference: a Wald interval of a non-diverged odds ratio from the glm itself
  m <- stats::glm(y ~ g, binomial, d)
  w <- exp(stats::confint.default(m, "gL02"))
  expect_equal(f$effect_sizes$conf_low[f$effect_sizes$term == "gL02"], unname(w[1]),
               tolerance = 1e-6)
  rmd <- test_path("..", "..", "vignettes", "articles", "binary.Rmd")
  skip_if_not(file.exists(rmd))
  txt <- paste(readLines(rmd, warn = FALSE), collapse = " ")
  expect_false(grepl("unless separation was detected", txt, fixed = TRUE))
  expect_true(grepl("That holds under separation too", txt, fixed = TRUE))
})

test_that("docs-consistency-F2: anova_glm says when non-integer weights set the scale of a Poisson fit", {
  d <- withr::with_seed(1, data.frame(g = rep(c("a", "b", "c"), each = 10),
                                      y = stats::rpois(30, 5),
                                      w = stats::runif(30, 0.5, 2)))
  d$w10 <- d$w * 10
  f1 <- anova_glm(d, "y", "g", family = poisson(), weights = "w", plots = FALSE)
  f10 <- anova_glm(d, "y", "g", family = poisson(), weights = "w10", plots = FALSE)
  # reference: glm's own statistics and Pearson dispersion scale with the weights
  m1 <- stats::glm(y ~ g, poisson, d, weights = w)
  m10 <- stats::glm(y ~ g, poisson, d, weights = w10)
  expect_equal(f10$anova$statistic[1], car::Anova(m10)[["LR Chisq"]][1], tolerance = 1e-8)
  expect_equal(f10$anova$statistic[1] / f1$anova$statistic[1], 10, tolerance = 1e-6)
  disp <- function(m) sum(m$prior.weights * (m$y - fitted(m))^2 / fitted(m)) / m$df.residual
  expect_equal(f10$assumptions$dispersion, disp(m10), tolerance = 1e-8)
  expect_equal(f1$assumptions$dispersion, disp(m1), tolerance = 1e-8)
  for (f in list(f1, f10)) {
    expect_length(grep("so their scale matters", f$notes, fixed = TRUE), 1L)
    expect_true(any(grepl("Pearson dispersion in $assumptions$dispersion, and so the overdispersion check, scale with the weights",
                          f$notes, fixed = TRUE)))
  }
  # a family that estimates its dispersion absorbs the scale: no note
  q <- anova_glm(d, "y", "g", family = quasipoisson(), weights = "w10", plots = FALSE)
  expect_false(any(grepl("so their scale matters", q$notes, fixed = TRUE)))
  # whole-number weights are frequencies: no note
  d$k <- rep(1:3, 10)
  fk <- anova_glm(d, "y", "g", family = poisson(), weights = "k", plots = FALSE)
  expect_false(any(grepl("so their scale matters", fk$notes, fixed = TRUE)))
})

test_that("countglm-C1/R1: a constant large-valued group does not hide real variation in the others", {
  d <- withr::with_seed(1, data.frame(g = gl(2, 10),
                                      y = c(stats::rpois(10, 5), rep(1e7, 10))))
  m <- stats::glm(y ~ g, quasipoisson, d)
  # reference: real residual variation, and car's F test
  expect_gt(m$deviance, 1)
  ref <- car::Anova(m, test.statistic = "F")
  g <- anova_glm(d, "y", "g", family = "quasipoisson", plots = FALSE)
  expect_equal(g$anova$statistic[1], ref[["F value"]][1], tolerance = 1e-6)
  k <- anova_count(d, "y", "g", model = "quasipoisson", plots = FALSE)
  expect_equal(k$anova$statistic[1], ref[["F value"]][1], tolerance = 1e-6)
  expect_no_error(anova_count(d, "y", "g", model = "negbin", test_statistic = "F",
                              plots = FALSE))
  d3 <- withr::with_seed(3, data.frame(g = gl(3, 8),
    y = c(stats::runif(8, 1, 3), stats::runif(8, 2, 5), rep(5e6, 8))))
  d3$y[1:16] <- round(d3$y[1:16])
  m3 <- stats::glm(y ~ g, quasipoisson, d3)
  ref3 <- car::Anova(m3, test.statistic = "F")
  expect_equal(stats::df.residual(m3), 21L)
  k3 <- anova_count(d3, "y", "g", model = "quasipoisson", plots = FALSE)
  expect_equal(k3$anova$statistic[1], ref3[["F value"]][1], tolerance = 1e-6)
  expect_equal(k3$anova$p_value[1], ref3[["Pr(>F)"]][1], tolerance = 1e-6)
  g3 <- anova_glm(d3, "y", "g", family = quasipoisson(), plots = FALSE)
  expect_equal(g3$anova$statistic[1], ref3[["F value"]][1], tolerance = 1e-6)
  # a Gamma fit with real variation is not refused as exact (car's own refit
  # of this model fails, which the notes report, as before)
  d2 <- withr::with_seed(1, data.frame(g = gl(2, 10),
    y = c(stats::rgamma(10, 100, 100), rep(1e6, 10))))
  expect_no_error(anova_glm(d2, "y", "g", family = Gamma(link = "log"), plots = FALSE))
  # a large constant group next to a group constant at 0 is still exact
  dz <- data.frame(g = gl(2, 10), y = rep(c(0, 1e7), each = 10))
  expect_error(anova_glm(dz, "y", "g", family = quasipoisson(), plots = FALSE),
               "fits response `y` exactly")
})

test_that("countglm-C2: anova_glm's F note on a saturated Poisson fit says why F is NA", {
  sat <- data.frame(y = c(1L, 2L, 4L), g = c("a", "b", "c"))
  m <- stats::glm(y ~ g, poisson, transform(sat, g = factor(g)))
  expect_identical(stats::df.residual(m), 0L)   # reference: nothing to estimate from
  r <- anova_glm(sat, "y", "g", family = poisson(), test_statistic = "F", plots = FALSE)
  expect_true(is.na(r$anova$statistic[1]))
  expect_true(any(grepl("cannot be estimated here as there are no residual degrees of freedom (so the F test is NA)",
                        r$notes, fixed = TRUE)))
  expect_false(any(grepl("Pearson statistic (NA)", r$notes, fixed = TRUE)))
})

test_that("countglm-R2: MASS::negative.binomial(theta) in glm() absorbs the scale of the weights", {
  ins <- InsectSprays
  ins$w <- withr::with_seed(1, stats::runif(72, .5, 2))
  ins$w10 <- ins$w * 10
  fam <- MASS::negative.binomial(10)
  for (ts in c("LR", "Wald")) {
    a <- anova_glm(ins, "count", "spray", family = fam, weights = "w",
                   test_statistic = ts, plots = FALSE)
    b <- anova_glm(ins, "count", "spray", family = fam, weights = "w10",
                   test_statistic = ts, plots = FALSE)
    ref <- car::Anova(stats::glm(count ~ spray, fam, ins, weights = w10),
                      test.statistic = ts)
    expect_equal(b$anova$statistic[1], ref[[if (ts == "LR") "LR Chisq" else "Chisq"]][1],
                 tolerance = 1e-6)
    expect_equal(a$anova$statistic[1], b$anova$statistic[1], tolerance = 1e-6)
  }
})
