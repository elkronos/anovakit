# Regression tests for the third review of anova_glm() (and the parts of
# anova_bin() and anova_count() the same items touch). Each is checked against
# an independent reference: base R, car, emmeans, or a hand computation.

test_that("R3-01: a frequency-weighted 0/1 table fits like the expanded data under cloglog", {
  agg <- data.frame(g = rep(c("a", "b"), each = 2), y = c(0, 1, 0, 1),
                    w = c(60, 40, 45, 55))
  ex <- agg[rep(seq_len(nrow(agg)), agg$w), c("g", "y")]
  expect_silent(fa <- anova_glm(agg, "y", "g", family = binomial("cloglog"),
                                weights = "w", plots = FALSE))
  ref <- stats::glm(y ~ g, family = binomial("cloglog"), data = ex)
  ref0 <- stats::glm(y ~ 1, family = binomial("cloglog"), data = ex)
  expect_equal(fa$anova$statistic, ref0$deviance - ref$deviance, tolerance = 1e-6)
  expect_equal(fa$anova$statistic, 4.528568, tolerance = 1e-6)
  expect_equal(fa$model$deviance, ref$deviance, tolerance = 1e-8)
  expect_equal(unname(stats::coef(fa$model)), unname(stats::coef(ref)),
               tolerance = 1e-6)
  # the marginal probabilities are the observed proportions
  expect_equal(fa$emmeans$estimate, c(0.40, 0.55), tolerance = 1e-6)
  expect_equal(fa$model_stats$AIC, stats::AIC(ref), tolerance = 1e-8)
  expect_false(any(grepl("separation|did not converge", fa$notes)))

  # A two-factor table: car's refits of reduced models start well too.
  ex2 <- withr::with_seed(3, {
    dd <- data.frame(a = sample(c("x", "y", "z"), 300, TRUE),
                     b = sample(c("p", "q"), 300, TRUE))
    dd$y <- stats::rbinom(300, 1, 0.4)
    dd
  })
  agg2 <- stats::aggregate(w ~ a + b + y, transform(ex2, w = 1), sum)
  for (type in c("II", "III")) {
    f2 <- anova_glm(agg2, "y", c("a", "b"), family = binomial("cloglog"),
                    weights = "w", type = type, plots = FALSE)
    r2 <- if (type == "II") {
      car::Anova(stats::glm(y ~ a + b, binomial("cloglog"), ex2))
    } else {
      withr::with_options(list(contrasts = c("contr.sum", "contr.poly")),
        car::Anova(stats::glm(y ~ a + b, binomial("cloglog"), ex2), type = 3))
    }
    expect_equal(f2$anova$statistic, unname(r2[["LR Chisq"]]), tolerance = 1e-6)
  }
})

test_that("R3-02: binomial AIC and BIC are NA with non-integer weights, as in anova_bin()", {
  d <- withr::with_seed(2, {
    dd <- data.frame(g = sample(c("a", "b", "c"), 200, TRUE))
    dd$y <- stats::rbinom(200, 1, 0.4)
    dd$w <- stats::runif(200, 0.2, 0.45)
    dd
  })
  fg <- anova_glm(d, "y", "g", family = binomial, weights = "w", plots = FALSE)
  m <- suppressWarnings(stats::glm(y ~ g, binomial, d, weights = w))
  # R's log-likelihood rounds every weight here to 0: it describes no data
  expect_equal(as.numeric(stats::logLik(m)), 0)
  expect_true(is.na(fg$model_stats$AIC))
  expect_true(is.na(fg$model_stats$BIC))
  expect_true(any(grepl("AIC and BIC in $model_stats are NA", fg$notes,
                        fixed = TRUE)))
  fb <- suppressWarnings(anova_bin(d, "y", "g", weights = "w", plots = FALSE))
  expect_true(is.na(fb$model_stats$AIC))
  # Whole-number weights keep AIC and BIC
  d$w2 <- 1 + (seq_len(nrow(d)) %% 3)
  f2 <- anova_glm(d, "y", "g", family = binomial, weights = "w2", plots = FALSE)
  m2 <- stats::glm(y ~ g, binomial, d, weights = w2)
  expect_equal(f2$model_stats$AIC, stats::AIC(m2))
  expect_false(any(grepl("AIC and BIC", f2$notes)))
})

collide_data <- function() {
  withr::with_seed(14, {
    n <- 400
    d <- data.frame(a = sample(c("b0", "b1"), n, TRUE),
                    ab = sample(c("0", "1"), n, TRUE))
    d$y <- 1 + 1.5 * (d$a == "b1") - 1 * (d$ab == "1") + stats::rnorm(n)
    d$c <- stats::rpois(n, exp(1 + 0.7 * (d$a == "b1") - 0.5 * (d$ab == "1")))
    d$s <- stats::rbinom(n, 1, stats::plogis(-0.5 + 1.5 * (d$a == "b1") -
                                               1 * (d$ab == "1")))
    d
  })
}

test_that("R3-04: colliding coefficient names do not corrupt Gaussian marginal means", {
  d <- collide_data()
  f <- anova_glm(d, "y", c("a", "ab"), plots = FALSE)
  expect_false(anyDuplicated(names(stats::coef(f$model))) > 0L)
  expect_true(any(grepl("produce the same coefficient name", f$notes)))
  expect_false(any(grepl("unbounded", f$notes)))
  # reference: lm on renamed columns, predictions averaged over the grid
  d2 <- d; names(d2)[names(d2) == "ab"] <- "zz"
  m <- stats::lm(y ~ a + zz, d2)
  grid <- expand.grid(a = c("b0", "b1"), zz = c("0", "1"))
  grid$p <- stats::predict(m, grid)
  ref <- c(tapply(grid$p, grid$a, mean), tapply(grid$p, grid$zz, mean))
  expect_equal(f$emmeans$estimate, unname(ref), tolerance = 1e-10)
  expect_equal(nrow(f$assumptions$coefficients), 3L)
  # the saturated model's marginal means are the raw cell means
  fi <- anova_glm(d, "y", c("a", "ab"), interaction = TRUE, plots = FALSE)
  cm <- stats::aggregate(y ~ a + ab, d, mean)
  key <- paste(fi$emmeans$a, fi$emmeans$ab)
  expect_equal(fi$emmeans$estimate, cm$y[match(key, paste(cm$a, cm$ab))],
               tolerance = 1e-10)
  # the tests are unchanged by the relabelling
  expect_equal(f$anova$statistic, anova_glm(d2, "y", c("a", "zz"),
                                            plots = FALSE)$anova$statistic)
})

test_that("R3-06: estimates and tests do not depend on the units of the response", {
  d <- withr::with_seed(4, {
    dd <- data.frame(g = factor(rep(c("a", "b", "c"), each = 15)))
    dd$y <- stats::rgamma(45, 3, 1 / c(1, 2, 3)[dd$g])
    dd
  })
  d2 <- transform(d, y = y * 1e-6)
  f2 <- anova_glm(d2, "y", "g", family = gaussian(link = "log"), plots = FALSE,
                  ci_method = "wald")
  ref <- stats::glm(y ~ g, data = d2, family = gaussian(link = "log"),
                    control = stats::glm.control(epsilon = 1e-14, maxit = 100),
                    start = c(log(mean(d2$y)), 0, 0))
  ra <- car::Anova(ref, test.statistic = "F")
  expect_equal(f2$anova$statistic[1], ra[["F value"]][1], tolerance = 1e-6)
  expect_equal(f2$anova$statistic[1], 5.37681, tolerance = 1e-5)
  expect_equal(unname(stats::coef(f2$model))[-1], unname(stats::coef(ref))[-1],
               tolerance = 1e-6)
  f1 <- anova_glm(d, "y", "g", family = gaussian(link = "log"), plots = FALSE,
                  ci_method = "wald")
  expect_equal(f2$anova$p_value[1], f1$anova$p_value[1], tolerance = 1e-6)
  expect_equal(f2$posthoc$ratio, f1$posthoc$ratio, tolerance = 1e-6)
  # inverse Gaussian, log link, in very large units
  fi1 <- anova_glm(d, "y", "g", family = inverse.gaussian(link = "log"),
                   plots = FALSE, ci_method = "wald")
  fi2 <- anova_glm(transform(d, y = y * 1e12), "y", "g",
                   family = inverse.gaussian(link = "log"), plots = FALSE,
                   ci_method = "wald")
  expect_equal(fi2$anova$statistic[1], fi1$anova$statistic[1], tolerance = 1e-6)
  ri <- stats::glm(y ~ g, data = d, family = inverse.gaussian(link = "log"),
                   control = stats::glm.control(epsilon = 1e-14, maxit = 100))
  expect_equal(fi1$anova$statistic[1],
               car::Anova(ri, test.statistic = "F")[["F value"]][1],
               tolerance = 1e-6)
  # precision weights on a tiny scale
  d$w <- withr::with_seed(5, stats::runif(45, 0.5, 2))
  fw1 <- anova_glm(d, "y", "g", family = gaussian(link = "log"), weights = "w",
                   plots = FALSE, ci_method = "wald")
  fw2 <- anova_glm(transform(d, w = w * 1e-12), "y", "g",
                   family = gaussian(link = "log"), weights = "w",
                   plots = FALSE, ci_method = "wald")
  expect_equal(fw2$anova$statistic[1], fw1$anova$statistic[1], tolerance = 1e-6)
})

test_that("R3-10: a constant response stops quasi and Gamma fits instead of giving F = 1e12", {
  d <- data.frame(y = 3, g = rep(c("a", "b", "c"), each = 3))
  # reference: the estimated dispersion is rounding error
  m <- stats::glm(y ~ g, quasipoisson, d)
  expect_lt(summary(m)$dispersion, 1e-20)
  expect_error(anova_glm(d, "y", "g", family = quasipoisson(), plots = FALSE),
               "takes the same value \\(3\\) in every analysed row")
  expect_error(anova_count(d, "y", "g", model = "quasipoisson", plots = FALSE),
               "takes the same value")
  expect_error(anova_glm(data.frame(y = 2, g = rep(c("a", "b"), each = 3)), "y",
                         "g", family = Gamma(link = "log"), plots = FALSE),
               "fits it exactly")
  # constant within every group
  dw <- data.frame(y = rep(c(2, 4, 7), each = 3), g = rep(c("a", "b", "c"), each = 3))
  expect_error(anova_glm(dw, "y", "g", family = quasipoisson(), plots = FALSE),
               "fits response `y` exactly")
  # A fixed-dispersion fit of the same data has a defined LR test of 0
  fp <- anova_glm(d, "y", "g", family = poisson(), plots = FALSE)
  expect_equal(fp$anova$statistic, 0, tolerance = 1e-8)
  # data with residual variation are unaffected
  d$y <- c(3, 4, 2, 5, 6, 4, 9, 8, 10)
  f <- anova_glm(d, "y", "g", family = quasipoisson(), plots = FALSE)
  ref <- car::Anova(stats::glm(y ~ g, quasipoisson, d), test.statistic = "F")
  expect_equal(f$anova$statistic[1], ref[["F value"]][1])
})

test_that("R3-11: a Gaussian perfect fit stops with a clear error, not NaN or F = 1e31", {
  expect_error(anova_glm(data.frame(y = 5, g = rep(c("a", "b", "c"), each = 3)),
                         "y", "g", plots = FALSE),
               "takes the same value \\(5\\)")
  d2 <- data.frame(y = rep(c(1, 2, 3), each = 3), g = rep(c("a", "b", "c"), each = 3))
  # reference: base R itself calls this fit unreliable
  expect_warning(stats::anova(stats::lm(y ~ g, d2)), "perfect fit")
  expect_error(anova_glm(d2, "y", "g", plots = FALSE), "fits response `y` exactly")
  expect_error(anova_glm(d2, "y", "g", family = inverse.gaussian(), plots = FALSE),
               "exactly")
  # a saturated model keeps its NA-inference path
  sat <- data.frame(g = factor(c("a", "b", "c")), y = c(1, 2, 4))
  expect_silent(fs <- anova_glm(sat, "y", "g", plots = FALSE))
  expect_true(any(grepl("saturated", fs$notes)))
})

test_that("R3-16: an LR request on a quasi family becomes F; other estimated dispersions are labelled", {
  d <- withr::with_seed(12, {
    dd <- data.frame(g = rep(c("a", "b", "c"), each = 25))
    dd$cnt <- stats::rnbinom(75, mu = rep(c(3, 5, 8), each = 25), size = 1.5)
    dd
  })
  f <- anova_glm(d, "cnt", "g", family = quasipoisson(), test_statistic = "LR",
                 plots = FALSE)
  m <- stats::glm(cnt ~ g, quasipoisson, d)
  ref <- stats::anova(stats::glm(cnt ~ 1, quasipoisson, d), m, test = "F")
  expect_equal(f$anova$statistic[1], ref$F[2], tolerance = 1e-5)
  expect_equal(f$anova$p_value[1], ref[["Pr(>F)"]][2], tolerance = 1e-4)
  expect_identical(.omnibus_header(f$anova), "Omnibus test (Type II, F)")
  expect_true(any(grepl("has no likelihood, so the analysis of deviance uses an F test",
                        f$notes)))
  k <- anova_count(d, "cnt", "g", model = "quasipoisson", test_statistic = "LR",
                   plots = FALSE)
  expect_equal(f$anova$statistic[1], k$anova$statistic[1])
  # Gamma: kept as asked, but not called a likelihood-ratio test
  d$pos <- d$cnt + 0.5
  fg <- anova_glm(d, "pos", "g", family = Gamma("log"), test_statistic = "LR",
                  plots = FALSE)
  rg <- car::Anova(stats::glm(pos ~ g, Gamma("log"), d), test.statistic = "LR")
  expect_equal(fg$anova$statistic[1], rg[["LR Chisq"]][1])
  expect_no_match(.omnibus_header(fg$anova), "likelihood-ratio")
  expect_match(.omnibus_header(fg$anova), "scaled-deviance chi-square")
  expect_true(any(grepl("scaled-deviance chi-square rather than a likelihood-ratio test",
                        fg$notes)))
  # a fixed-dispersion LR test keeps its label
  fp <- anova_glm(d, "cnt", "g", family = poisson(), plots = FALSE)
  expect_match(.omnibus_header(fp$anova), "likelihood-ratio chi-square")
})

test_that("R3-21: rows where the profile interval and the Wald p-value disagree are noted", {
  d <- data.frame(g = factor(rep(c("a", "b"), each = 8)),
                  s = c(1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0))
  m <- stats::glm(s ~ g, binomial, d)
  prof <- exp(suppressMessages(stats::confint(m))[2, ])
  wald_p <- summary(m)$coefficients[2, 4]
  expect_gt(prof[1], 1)          # the profile interval excludes 1 ...
  expect_gt(wald_p, 0.05)        # ... while the Wald p-value is above 0.05
  b <- anova_bin(d, "s", "g", plots = FALSE)
  expect_equal(unname(c(b$effect_sizes$conf_low, b$effect_sizes$conf_high)),
               unname(prof), tolerance = 1e-6)
  expect_true(any(grepl("the interval of gb excludes the null odds ratio of 1 while its p-value is 0.05 or more",
                        b$notes, fixed = TRUE)))
  expect_false(any(grepl("opposite conclusions",
                         anova_bin(d, "s", "g", ci_method = "wald",
                                   plots = FALSE)$notes)))

  d2 <- withr::with_seed(17, {
    dd <- data.frame(g = factor(rep(c("a", "b"), c(4, 4))))
    dd$y <- stats::rgamma(8, 1, 1 / c(1, 3)[as.integer(dd$g)])
    dd
  })
  f <- anova_glm(d2, "y", "g", family = Gamma("log"), plots = FALSE,
                 posthoc = FALSE)
  co <- f$assumptions$coefficients
  m2 <- stats::glm(y ~ g, Gamma("log"), d2)
  expect_equal(co$p_value[2], summary(m2)$coefficients[2, 4])
  expect_lt(co$p_value[2], 0.05)
  expect_lt(co$conf_low[2], 0)
  expect_true(any(grepl("the interval of gb includes the null value 0 while its p-value is below 0.05",
                        f$notes, fixed = TRUE)))
  # agreement gives no note
  ok <- anova_glm(fx_overdispersed(), "y", "g", family = "poisson", plots = FALSE)
  expect_false(any(grepl("opposite conclusions", ok$notes)))
})

test_that("R3-22: Gaussian identity profile intervals are the t intervals, without profiling", {
  d <- withr::with_seed(1, data.frame(y = stats::rnorm(2000),
    g = factor(sample(sprintf("L%02d", 1:20), 2000, TRUE))))
  # Profiling must not be attempted: a failure would surface as a note.
  testthat::local_mocked_bindings(.profile_ci = function(...) stop("profiled"))
  f <- anova_glm(d, "y", "g", posthoc = FALSE, plots = FALSE)
  co <- f$assumptions$coefficients
  expect_equal(unname(as.matrix(co[, c("conf_low", "conf_high")])),
               unname(stats::confint(stats::lm(y ~ g, d))), tolerance = 1e-10)
  expect_true(all(co$ci_method == "profile"))
  expect_false(any(grepl("profiled|failed", f$notes)))
})

test_that("R3-22: default profile intervals give way to Wald when rows times coefficients squared is large", {
  d <- withr::with_seed(2, data.frame(
    c = stats::rpois(30000, 3),
    g = factor(sample(sprintf("L%02d", 1:20), 30000, TRUE))))
  testthat::local_mocked_bindings(.profile_ci = function(...) stop("profiled"))
  f <- anova_glm(d, "c", "g", family = poisson(), posthoc = FALSE, plots = FALSE)
  co <- f$assumptions$coefficients
  expect_true(all(co$ci_method == "wald"))
  m <- stats::glm(c ~ g, poisson, d)
  expect_equal(unname(as.matrix(co[, c("conf_low", "conf_high")])),
               unname(stats::confint.default(m)), tolerance = 1e-8)
  expect_true(any(grepl("20 coefficients and 30,000 rows", f$notes, fixed = TRUE)))
  # asked for by name, profiling is attempted
  g <- anova_glm(d[1:3000, ], "c", "g", family = poisson(), posthoc = FALSE,
                 plots = FALSE, ci_method = "profile")
  expect_true(any(grepl("Profile-likelihood intervals failed", g$notes)))
})

test_that("R3-27: a Type III request with an empty cell still reports effect sizes", {
  mt <- transform(mtcars, cyl = factor(cyl), gear = factor(gear))
  expect_equal(sum(mt$cyl == "8" & mt$gear == "4"), 0L)
  f <- anova_glm(mt, "mpg", c("cyl", "gear"), interaction = TRUE, type = "III",
                 plots = FALSE)
  expect_match(f$method, "Type II\\)")
  expect_equal(nrow(f$effect_sizes), 3L)
  ref <- suppressMessages(car::Anova(stats::lm(mpg ~ cyl * gear, mt), type = 2))
  ss <- ref[["Sum Sq"]]
  expect_equal(f$effect_sizes$partial_eta_sq, ss[1:3] / (ss[1:3] + ss[4]),
               tolerance = 1e-8)
  src <- testthat::test_path("..", "..", "R", "anova_glm.R")
  skip_if_not(file.exists(src))
  expect_false(any(grepl("a Type III table\\s*$", readLines(src))))
  expect_false(any(grepl("with an empty cell, say", readLines(src), fixed = TRUE)))
})

test_that("R3-33: colliding coefficient names do not corrupt Poisson or binomial results", {
  skip_if_not_installed("emmeans")
  d <- collide_data()
  d2 <- d; names(d2)[names(d2) == "ab"] <- "zz"
  d2$a <- factor(d2$a); d2$zz <- factor(d2$zz)
  f <- anova_glm(d, "c", c("a", "ab"), family = poisson(), plots = FALSE)
  ref <- stats::glm(c ~ a + zz, poisson, d2)
  ea <- summary(emmeans::emmeans(ref, ~ a, type = "response"))
  ez <- summary(emmeans::emmeans(ref, ~ zz, type = "response"))
  expect_equal(f$emmeans$estimate, c(ea$rate, ez$rate), tolerance = 1e-8)
  rz <- summary(emmeans::contrast(emmeans::emmeans(ref, ~ zz, type = "response"),
                                  "pairwise"))
  expect_equal(f$posthoc$ratio[f$posthoc$term == "ab"], rz$ratio, tolerance = 1e-8)
  expect_equal(nrow(f$assumptions$coefficients), 3L)
  expect_equal(f$assumptions$coefficients$estimate, unname(stats::coef(ref)),
               tolerance = 1e-8)
  fb <- anova_glm(d, "s", c("a", "ab"), family = binomial(), plots = FALSE)
  rb <- stats::glm(s ~ a + zz, binomial, d2)
  pb <- c(summary(emmeans::emmeans(rb, ~ a, type = "response"))$prob,
          summary(emmeans::emmeans(rb, ~ zz, type = "response"))$prob)
  expect_equal(fb$emmeans$estimate, pb, tolerance = 1e-8)
  expect_equal(fb$emmeans$estimate,
               anova_bin(d, "s", c("a", "ab"), plots = FALSE)$emmeans$estimate,
               tolerance = 1e-8)
})

test_that("R3-52: a quasi-likelihood F test on a Poisson or binomial fit says so in $method", {
  d <- withr::with_seed(2, {
    dd <- data.frame(g = factor(rep(c("a", "b", "c"), each = 20)))
    dd$cnt <- stats::rnbinom(60, mu = c(a = 5, b = 5, c = 7.5)[as.character(dd$g)],
                             size = 1.5)
    dd
  })
  k <- anova_count(d, "cnt", "g", model = "poisson", test_statistic = "F",
                   plots = FALSE)
  g <- anova_glm(d, "cnt", "g", family = poisson(), test_statistic = "F",
                 plots = FALSE)
  ref <- car::Anova(stats::glm(cnt ~ g, quasipoisson, d), test.statistic = "F")
  expect_equal(k$anova$statistic[1], ref[["F value"]][1])
  expect_equal(g$anova$statistic[1], ref[["F value"]][1])
  expect_match(k$method, "quasi-likelihood F test", fixed = TRUE)
  expect_match(g$method, "quasi-likelihood F test", fixed = TRUE)
  # the default tests are not relabelled
  expect_no_match(anova_count(d, "cnt", "g", model = "poisson", plots = FALSE)$method,
                  "quasi-likelihood")
  expect_no_match(anova_glm(d, "cnt", "g", family = poisson(), plots = FALSE)$method,
                  "quasi-likelihood")
  expect_no_match(anova_count(d, "cnt", "g", model = "quasipoisson",
                              plots = FALSE)$method, "quasi-likelihood F")
})

test_that("R3-59: the overdispersion note follows the Pearson test at the 5% level", {
  pois_data <- function(seed) withr::with_seed(seed, {
    dd <- data.frame(g = factor(rep(c("a", "b", "c"), each = 20)))
    dd$y <- stats::rnbinom(60, mu = rep(c(4, 6, 8), each = 20), size = 12)
    dd
  })
  pearson_p <- function(d) {
    m <- stats::glm(y ~ g, poisson, d)
    x2 <- sum(stats::residuals(m, type = "pearson")^2)
    stats::pchisq(x2, stats::df.residual(m), lower.tail = FALSE)
  }
  under <- pois_data(1)          # Pearson test p just under 0.05
  over <- pois_data(37)          # dispersion above 1, p just over 0.05
  expect_gt(pearson_p(under), 0.02)
  expect_lt(pearson_p(under), 0.05)
  expect_gt(pearson_p(over), 0.05)
  expect_lt(pearson_p(over), 0.1)
  fu <- anova_glm(under, "y", "g", family = "poisson", plots = FALSE)
  fo <- anova_glm(over, "y", "g", family = "poisson", plots = FALSE)
  expect_gt(fo$assumptions$dispersion, 1)
  expect_true(any(grepl("^Pearson dispersion is", fu$notes)))
  expect_true(any(grepl(sprintf("p = %.3f", pearson_p(under)), fu$notes,
                        fixed = TRUE)) ||
                any(grepl(.fmt_p(pearson_p(under)), fu$notes, fixed = TRUE)))
  expect_false(any(grepl("^Pearson dispersion is", fo$notes)))
})
