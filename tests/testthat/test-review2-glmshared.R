# Detailed regression tests: run locally and on CI (NOT_CRAN=true), skipped
# on CRAN to keep its check time short. The core suite runs everywhere.
skip_on_cran()

# Review 2, area glmshared: anova_glm() and the helpers it shares.
# Every expectation is checked against base R, car, emmeans, sandwich or a
# hand computation, not against the package's own code.

r2_bin_table <- function() {
  data.frame(g = factor(rep(c("a", "b", "c"), each = 2)), y = rep(c(0, 1), 3),
             w = c(30, 10, 20, 20, 12, 28))
}
r2_count_table <- function() {
  data.frame(g = factor(rep(c("a", "b", "c"), each = 4)),
             y = c(0, 1, 2, 3, 1, 2, 3, 5, 2, 4, 6, 7),
             w = c(3, 5, 2, 1, 4, 4, 2, 1, 1, 2, 3, 4))
}
r2_expand <- function(d) d[rep(seq_len(nrow(d)), d$w), , drop = FALSE]

test_that("R2-GLMSHARED-01: quasibinomial and NB(theta) frequency tables equal the expanded data", {
  skip_if_not_installed("MASS")
  d <- r2_bin_table()
  de <- r2_expand(d)
  f <- anova_glm(d, "y", "g", family = "quasibinomial", weights = "w", plots = FALSE)
  m <- stats::glm(y ~ g, family = stats::quasibinomial(), data = de)
  ref <- car::Anova(m, test.statistic = "F")
  expect_equal(f$anova$statistic[1], ref[["F value"]][1], tolerance = 1e-8)
  expect_equal(f$anova$p_value[1], ref[["Pr(>F)"]][1], tolerance = 1e-8)
  expect_equal(f$anova$df[2], stats::df.residual(m))
  ct <- stats::coef(summary(m))
  expect_equal(f$assumptions$coefficients$se, unname(ct[, "Std. Error"]), tolerance = 1e-4)
  em <- summary(emmeans::emmeans(m, "g", type = "response"))
  expect_equal(f$emmeans$df, rep(stats::df.residual(m), 3))
  expect_equal(f$emmeans$estimate, em$prob, tolerance = 1e-6)
  expect_equal(f$emmeans$se, em$SE, tolerance = 1e-4)
  expect_true(any(grepl("frequency weights", f$notes)))

  dp <- r2_count_table()
  dpe <- r2_expand(dp)
  fn <- anova_glm(dp, "y", "g", family = MASS::negative.binomial(2), weights = "w",
                  plots = FALSE)
  mn <- stats::glm(y ~ g, family = MASS::negative.binomial(2), data = dpe)
  refn <- car::Anova(mn, test.statistic = "F")
  expect_equal(fn$anova$statistic[1], refn[["F value"]][1], tolerance = 1e-8)
  expect_equal(fn$anova$p_value[1], refn[["Pr(>F)"]][1], tolerance = 1e-8)
  expect_equal(fn$assumptions$dispersion, summary(mn)$dispersion, tolerance = 1e-4)
  # a Gaussian fit keeps reading the weights as precision weights
  fg <- anova_glm(d, "y", "g", weights = "w", plots = FALSE)
  expect_equal(fg$anova$df[2], stats::df.residual(stats::lm(y ~ g, d, weights = w)))
})

test_that("R2-GLMSHARED-02: a factor named like an emmeans column keeps the others' statistics", {
  withr::local_seed(9)
  d <- data.frame(g = factor(rep(c("a", "b", "c"), 20)),
                  SE = factor(rep(c("u", "v"), each = 30)),
                  rate = factor(sample(c("p", "q"), 60, TRUE)))
  d$y <- stats::rpois(60, 3 + as.integer(d$g))
  f <- anova_glm(d, "y", c("g", "SE", "rate"), family = stats::poisson(), plots = FALSE)
  # reference: the same model with neutrally named factors
  dn <- d
  names(dn)[names(dn) == "SE"] <- "s1"
  names(dn)[names(dn) == "rate"] <- "r1"
  m <- stats::glm(y ~ g + s1 + r1, stats::poisson(), dn)
  for (h in c("g", "SE", "rate")) {
    hn <- c(g = "g", SE = "s1", rate = "r1")[[h]]
    ref <- summary(emmeans::emmeans(m, hn, type = "response"))
    got <- f$emmeans[f$emmeans$term == h, ]
    expect_equal(got$estimate, ref$rate, tolerance = 1e-6, info = h)
    expect_equal(got$se, ref$SE, tolerance = 1e-6, info = h)
    expect_identical(as.character(got[[h]]), as.character(ref[[hn]]), info = h)
  }

  d2 <- data.frame(a = factor(rep(c("a", "b", "c"), 20)),
                   df = factor(rep(c("u", "v"), each = 30)))
  d2$y <- stats::rnorm(60) + as.integer(d2$a)
  f2 <- anova_glm(d2, "y", c("a", "df"), plots = FALSE)
  ref2 <- summary(emmeans::emmeans(stats::lm(y ~ a + df, d2), "a"))
  rows <- f2$emmeans$term == "a"
  expect_equal(f2$emmeans$emm_df[rows], ref2$df)
  expect_equal(f2$emmeans$estimate[rows], ref2$emmean, tolerance = 1e-8)
  expect_false(anyNA(f2$emmeans$emm_df))
})

test_that("R2-GLMSHARED-03: anova() on the stored model of a Poisson frequency table works", {
  dp <- r2_count_table()
  f <- anova_glm(dp, "y", "g", family = "poisson", weights = "w", plots = FALSE)
  ref <- stats::anova(stats::glm(y ~ g, stats::poisson(), r2_expand(dp)), test = "Chisq")
  got <- stats::anova(f$model, test = "Chisq")
  expect_equal(got$Df[2], ref$Df[2])
  expect_equal(got$Deviance[2], ref$Deviance[2], tolerance = 1e-8)
  expect_equal(got[["Pr(>Chi)"]][2], ref[["Pr(>Chi)"]][2], tolerance = 1e-8)
  got2 <- stats::anova(stats::update(f$model, . ~ 1), f$model, test = "Chisq")
  expect_equal(got2$Df[2], 2)
  expect_equal(got2[["Pr(>Chi)"]][2], ref[["Pr(>Chi)"]][2], tolerance = 1e-8)
  # the package's own results stay those of the expanded data
  expect_equal(f$model_stats$df_residual, sum(dp$w) - 3)
})

test_that("R2-GLMSHARED-04: ci_method = 'wald' under separation gives Wald intervals", {
  withr::local_seed(5)
  d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 20)))
  d$y <- c(rep(0, 20), stats::rbinom(20, 1, .5), stats::rbinom(20, 1, .4))
  m <- suppressWarnings(stats::glm(y ~ g, stats::binomial(), d))
  fw <- anova_glm(d, "y", "g", family = stats::binomial(), ci_method = "wald",
                  plots = FALSE)
  ct <- fw$assumptions$coefficients
  ref <- stats::confint.default(m)
  expect_equal(ct$conf_low, unname(ref[, 1]), tolerance = 1e-6)
  expect_equal(ct$conf_high, unname(ref[, 2]), tolerance = 1e-6)
  expect_true(all(ct$ci_method == "wald"))
  # profiled rows say so
  fp <- anova_glm(d, "y", "g", family = stats::binomial(), plots = FALSE)
  cp <- fp$assumptions$coefficients
  expect_true(all(cp$ci_method == "profile"))
  expect_identical(cp$conf_high[cp$term == "gb"], Inf)
})

test_that("R2-GLMSHARED-05: a two-cell Wald comparison in a fixed-dispersion family is the omnibus test", {
  skip_if_not_installed("sandwich")
  withr::local_seed(22)
  d <- data.frame(g = factor(rep(c("a", "b"), each = 40)))
  d$y <- stats::rbinom(80, 1, c(.3, .6)[as.integer(d$g)])
  f <- anova_glm(d, "y", "g", family = stats::binomial(), test_statistic = "Wald",
                 vcov_type = "HC0", plots = FALSE)
  m <- stats::glm(y ~ g, stats::binomial(), d)
  V <- sandwich::vcovHC(m, type = "HC0")
  z <- stats::coef(m)[2] / sqrt(V[2, 2])
  ref <- suppressMessages(car::Anova(m, test.statistic = "Wald", vcov. = V))
  expect_equal(unname(z^2), ref$Chisq[1], tolerance = 1e-8)
  expect_equal(f$posthoc$p_value, ref[["Pr(>Chisq)"]][1], tolerance = 1e-8)
  note <- grep("two cells", f$notes, value = TRUE)
  expect_match(note, "is the omnibus test")
  # an estimated dispersion: t comparison against a chi-square omnibus test
  d$yg <- stats::rgamma(80, shape = 2, rate = c(1, 2)[as.integer(d$g)])
  fg <- anova_glm(d, "yg", "g", family = stats::Gamma(link = "log"),
                  test_statistic = "Wald", plots = FALSE)
  expect_match(grep("two cells", fg$notes, value = TRUE), "can differ")
})

test_that("R2-GLMSHARED-06: BIC of a Poisson frequency table counts observations", {
  dp <- r2_count_table()
  f <- anova_glm(dp, "y", "g", family = "poisson", weights = "w", plots = FALSE)
  m <- stats::glm(y ~ g, stats::poisson(), r2_expand(dp))
  expect_equal(f$model_stats$BIC, stats::BIC(m), tolerance = 1e-8)
  expect_equal(f$model_stats$AIC, stats::AIC(m), tolerance = 1e-8)
})

test_that("R2-GLMSHARED-07: whole-number weights are frequencies only in the documented families", {
  skip_if_not_installed("sandwich")
  d <- r2_bin_table()
  de <- r2_expand(d)
  # Gaussian (a linear probability model): precision weights, one row per unit
  fg <- anova_glm(d, "y", "g", weights = "w", vcov_type = "HC0", plots = FALSE)
  refg <- sqrt(diag(sandwich::vcovHC(stats::lm(y ~ g, d, weights = w), type = "HC0")))
  expect_equal(fg$assumptions$coefficients$se, unname(refg), tolerance = 1e-8)
  # binomial on a 0/1 response: frequency weights, as on the expanded data
  fb <- anova_glm(d, "y", "g", family = "binomial", weights = "w",
                  vcov_type = "HC0", plots = FALSE)
  refb <- sqrt(diag(sandwich::vcovHC(stats::glm(y ~ g, stats::binomial(), de),
                                     type = "HC0")))
  expect_equal(fb$assumptions$coefficients$se, unname(refb), tolerance = 1e-6)
})

test_that("R2-GLMSHARED-08: an unidentified intercept under sum-to-zero coding is (-Inf, Inf)", {
  withr::local_seed(5)
  d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 20)))
  d$y <- c(rep(0, 20), rep(1, 20), stats::rbinom(20, 1, .4))
  f <- anova_glm(d, "y", "g", family = stats::binomial(), type = "III", plots = FALSE)
  ct <- f$assumptions$coefficients
  int <- ct[ct$term == "(Intercept)", ]
  expect_identical(c(int$conf_low, int$conf_high), c(-Inf, Inf))
  expect_true(any(grepl("rule out no value", f$notes)))
  # Hand check: the intercept is the mean of the three cell logits. Fixed at
  # any b, the deviance still comes as close as we like to its infimum, that
  # of cell c alone: hold logit(c) at its estimate and send logit(a) to -M and
  # logit(b) to 3b - logit(c) + M.
  pc <- mean(d$y[d$g == "c"])
  dmin <- -2 * sum(stats::dbinom(d$y[d$g == "c"], 1, pc, log = TRUE))
  dev <- function(eta) -2 * sum(stats::dbinom(d$y, 1, stats::plogis(eta), log = TRUE))
  M <- 200
  for (b in c(-30, 0, 30)) {
    cell <- c(a = -M, b = 3 * b - stats::qlogis(pc) + M, c = stats::qlogis(pc))
    expect_equal(mean(cell), b)
    expect_lt(dev(cell[as.character(d$g)]) - dmin, 1e-6)
  }
})

test_that("R2-GLMSHARED-09: many coefficients get Wald intervals by default, profile on request", {
  withr::local_seed(17)
  d <- data.frame(g = factor(rep(sprintf("L%02d", 1:52), each = 8)))
  d$y <- stats::rpois(nrow(d), 4)
  f <- anova_glm(d, "y", "g", family = stats::poisson(), plots = FALSE, posthoc = FALSE)
  ct <- f$assumptions$coefficients
  ref <- stats::confint.default(stats::glm(y ~ g, stats::poisson(), d))
  expect_true(all(ct$ci_method == "wald"))
  expect_equal(ct$conf_low, unname(ref[, 1]), tolerance = 1e-6)
  expect_true(any(grepl("ci_method = \"profile\" explicitly", f$notes, fixed = TRUE)))
  fp <- anova_glm(d, "y", "g", family = stats::poisson(), plots = FALSE,
                  posthoc = FALSE, ci_method = "profile")
  expect_true(all(fp$assumptions$coefficients$ci_method == "profile"))
  # a small model keeps the profile default (stats::confint())
  small <- d[d$g %in% c("L01", "L02", "L03"), ]
  small$g <- droplevels(small$g)
  fs <- anova_glm(small, "y", "g", family = stats::poisson(), plots = FALSE)
  prof <- suppressMessages(stats::confint(stats::glm(y ~ g, stats::poisson(), small)))
  expect_equal(fs$assumptions$coefficients$conf_low, unname(prof[, 1]), tolerance = 1e-4)
})

test_that("R2-GLMSHARED-10: a call made through do.call() does not carry the data", {
  withr::local_seed(1)
  d <- data.frame(g = rep(c("a", "b", "c"), 100), y = stats::rnorm(300),
                  junk = stats::rnorm(300))
  f <- do.call(anova_glm, list(data = d, response = "y", groups = "g",
                               family = stats::gaussian(), plots = FALSE))
  expect_false(any(vapply(as.list(f$call), is.data.frame, logical(1))))
  expect_false(any(vapply(as.list(f$call), is.function, logical(1))))
  expect_identical(f$call[[1]], as.name("anova_glm"))
  expect_lt(as.numeric(utils::object.size(f$call)), as.numeric(utils::object.size(d)))
  out <- utils::capture.output(print(f))
  call_line <- grep("^Call:", out, value = TRUE)
  expect_lt(nchar(call_line), 200)
  # a direct call is kept as written
  g <- anova_glm(d, "y", "g", plots = FALSE)
  expect_identical(g$call, quote(anova_glm(data = d, response = "y", groups = "g", plots = FALSE)))
})
