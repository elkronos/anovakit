# Regression tests for the anova_glm() findings of the adversarial review.
# Each is checked against an independent reference: stats, car, emmeans,
# sandwich, effectsize, or a hand computation from the published formula.

# Local fixtures ----------------------------------------------------------------

glm_small_gauss <- function() {
  data.frame(g = factor(rep(c("a", "b", "c"), each = 2)),
             y = c(1.0, 2.0, 2.9, 4.1, 3.0, 5.2))
}

glm_gamma <- function(n_per = 4, seed = 4) {
  withr::with_seed(seed, {
    d <- data.frame(g = factor(rep(c("a", "b", "c"), each = n_per)))
    d$y <- stats::rgamma(3 * n_per, shape = 3,
                         rate = 3 / c(a = 1, b = 2, c = 3)[as.character(d$g)])
    d
  })
}

glm_all_ones <- function(seed = 3) {
  withr::with_seed(seed, {
    d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 20)))
    d$y <- stats::rbinom(60, 1, c(a = .3, b = .5, c = .6)[as.character(d$g)])
    d$y[d$g == "c"] <- 1
    d
  })
}

# F040 -------------------------------------------------------------------------

test_that("F040: Gaussian profile intervals equal confint(lm()), a t interval", {
  d <- glm_small_gauss()
  fit <- anova_glm(d, "y", "g", plots = FALSE)
  co <- fit$assumptions$coefficients
  ref <- stats::confint(stats::lm(y ~ g, data = d))
  expect_identical(unique(co$ci_method), "profile")
  expect_equal(co$conf_low, unname(ref[, 1]), tolerance = 1e-6)
  expect_equal(co$conf_high, unname(ref[, 2]), tolerance = 1e-6)
  # consistent with the t test in the same row: p > 0.05 and 0 inside
  expect_gt(co$p_value[co$term == "gc"], 0.05)
  expect_lt(co$conf_low[co$term == "gc"], 0)
})

test_that("F040: estimated-dispersion profiles are cut at qt(df.residual)", {
  d <- glm_gamma()
  fit <- anova_glm(d, "y", "g", family = stats::Gamma("log"), plots = FALSE)
  co <- fit$assumptions$coefficients
  m <- stats::glm(y ~ g, data = d, family = stats::Gamma("log"))
  phi <- summary(m)$dispersion
  X <- stats::model.matrix(m)
  crit <- stats::qt(0.975, stats::df.residual(m))
  # Hand computation: at each bound, the constrained fit's signed-root
  # deviance increase, scaled by the dispersion, equals the t cutoff.
  for (j in seq_len(ncol(X))) {
    for (b in c(co$conf_low[j], co$conf_high[j])) {
      mj <- stats::glm.fit(X[, -j, drop = FALSE], d$y, offset = X[, j] * b,
                           family = stats::Gamma("log"))
      expect_equal(sqrt((mj$deviance - stats::deviance(m)) / phi), crit,
                   tolerance = 1e-3)
    }
  }
  # and it is wider than the normal-cutoff interval confint() gives
  ref <- suppressMessages(stats::confint(m))
  expect_true(all(co$conf_high - co$conf_low > ref[, 2] - ref[, 1]))
})

test_that("F040: fixed-dispersion families keep confint()'s normal cutoffs", {
  d <- fx_binary(n = 150)
  fit <- anova_glm(d, "y", "g", family = "binomial", plots = FALSE)
  ref <- suppressMessages(stats::confint(
    stats::glm(y ~ g, family = stats::binomial(), data = d)))
  expect_equal(unname(as.matrix(fit$assumptions$coefficients[, c("conf_low", "conf_high")])),
               unname(ref), tolerance = 1e-6)
})

# F041 -------------------------------------------------------------------------

test_that("F041: inverse-link comparisons are response-scale differences, not 'ratio'", {
  d <- withr::with_seed(1, {
    dd <- data.frame(g = factor(rep(c("a", "b", "c"), each = 40)))
    dd$y <- stats::rgamma(120, shape = 5,
                          rate = 5 / c(a = 2, b = 4, c = 8)[as.character(dd$g)])
    dd
  })
  fit <- anova_glm(d, "y", "g", family = stats::Gamma(), plots = FALSE)
  expect_false("ratio" %in% names(fit$posthoc))
  # a one-way Gamma fit reproduces the group means
  mu <- tapply(d$y, d$g, mean)
  expect_equal(fit$posthoc$estimate[fit$posthoc$contrast == "a - b"],
               unname(mu["a"] - mu["b"]), tolerance = 1e-6)
  # emmeans' own response-scale comparison, for the standard errors
  m <- stats::glm(y ~ g, data = d, family = stats::Gamma())
  ref <- as.data.frame(emmeans::contrast(
    emmeans::regrid(emmeans::emmeans(m, "g", type = "response")),
    method = "pairwise"))
  expect_equal(fit$posthoc$se, ref$SE, tolerance = 1e-6)
  expect_true(any(grepl("inverse link", fit$notes)))
})

# F042 -------------------------------------------------------------------------

test_that("F042: McFadden's R squared is NA for a Gamma family, whatever the units", {
  d <- withr::with_seed(3, {
    dd <- data.frame(g = factor(rep(c("a", "b", "c"), each = 30)))
    dd$y <- stats::rgamma(90, shape = 20,
                          rate = 20 / c(a = 2, b = 3, c = 5)[as.character(dd$g)])
    dd
  })
  es <- lapply(c(1, 1e-3, 1e3), function(s) {
    f <- anova_glm(transform(d, y = y * s), "y", "g",
                   family = stats::Gamma("log"), plots = FALSE)
    expect_true(any(grepl("McFadden.*not reported", f$notes)))
    f$effect_sizes
  })
  m <- stats::glm(y ~ g, data = d, family = stats::Gamma("log"))
  for (e in es) {
    expect_true(is.na(e$estimate[e$measure == "mcfadden_r2"]))
    expect_equal(e$estimate[e$measure == "deviance_explained"],
                 1 - m$deviance / m$null.deviance)
  }
})

test_that("F042: binomial McFadden does not depend on aggregation", {
  agg <- withr::with_seed(16, {
    a <- data.frame(g = factor(rep(c("a", "b", "c"), each = 4)), n = 50)
    a$k <- stats::rbinom(12, 50, c(a = .3, b = .4, c = .5)[as.character(a$g)])
    a$p <- a$k / a$n
    a
  })
  long <- do.call(rbind, lapply(seq_len(12), function(i)
    data.frame(g = agg$g[i], y = rep(c(1, 0), c(agg$k[i], 50 - agg$k[i])))))
  fa <- anova_glm(agg, "p", "g", family = "binomial", weights = "n", plots = FALSE)
  fl <- anova_glm(long, "y", "g", family = "binomial", plots = FALSE)
  ref <- 1 - as.numeric(stats::logLik(stats::glm(y ~ g, binomial, data = long)) /
                          stats::logLik(stats::glm(y ~ 1, binomial, data = long)))
  mc <- function(f) f$effect_sizes$estimate[f$effect_sizes$measure == "mcfadden_r2"]
  expect_equal(mc(fa), ref, tolerance = 1e-8)
  expect_equal(mc(fl), ref, tolerance = 1e-8)
  expect_true(any(grepl("saturated model", fa$notes)))
})

test_that("F042: Poisson and negative binomial McFadden match logLik()", {
  skip_if_not_installed("MASS")
  d <- fx_overdispersed(n_per = 40)
  nbf <- MASS::negative.binomial(1.5)
  for (fam in list(stats::poisson(), nbf)) {
    fit <- anova_glm(d, "y", "g", family = fam, plots = FALSE, ci_method = "wald")
    full <- stats::glm(y ~ g, family = fam, data = d)
    null <- stats::glm(y ~ 1, family = fam, data = d)
    expect_equal(fit$effect_sizes$estimate[fit$effect_sizes$measure == "mcfadden_r2"],
                 as.numeric(1 - stats::logLik(full) / stats::logLik(null)),
                 tolerance = 1e-8)
  }
})

# F043 -------------------------------------------------------------------------

test_that("F043: nothing prints, and results do not depend on options(warn)", {
  d <- withr::with_seed(1, {
    dd <- data.frame(g = factor(rep(c("a", "b", "c"), each = 20)))
    dd$b <- stats::rbinom(60, 1, 0.5)
    dd$w <- stats::runif(60, 0.5, 2)
    dd
  })
  expect_silent(f1 <- anova_glm(d, "b", "g", family = "binomial", weights = "w",
                                plots = FALSE))
  mc <- f1$effect_sizes$estimate[f1$effect_sizes$measure == "mcfadden_r2"]
  # For a 0/1 response the Bernoulli McFadden is 1 - D / D0.
  m <- suppressWarnings(stats::glm(b ~ g, binomial, data = d, weights = w))
  expect_equal(mc, 1 - m$deviance / m$null.deviance, tolerance = 1e-8)
  withr::with_options(list(warn = 2), {
    f2 <- anova_glm(d, "b", "g", family = "binomial", weights = "w", plots = FALSE)
    expect_equal(f2$effect_sizes, f1$effect_sizes)
    sat <- data.frame(g = factor(c("a", "b", "c")), y = c(1, 2, 4))
    expect_silent(anova_glm(sat, "y", "g", vcov_type = "HC0", plots = FALSE))
    # a Poisson family with a non-integer response: no likelihood, no warnings
    dp <- data.frame(g = factor(rep(c("a", "b"), each = 6)),
                     y = c(1.5, 2, 2.5, 3, 1, 2, 4.5, 5, 6, 5.5, 4, 7))
    expect_silent(fp <- anova_glm(dp, "y", "g", family = "poisson", plots = FALSE))
  })
  expect_true(is.na(fp$effect_sizes$estimate[fp$effect_sizes$measure == "mcfadden_r2"]))
  expect_true(any(grepl("integer", fp$notes)))
})

# F044 / F161 -----------------------------------------------------------------

test_that("F044: partial omega squared ignores zero-weight rows (effectsize)", {
  skip_if_not_installed("effectsize")
  d <- withr::with_seed(10, {
    dd <- data.frame(g = factor(rep(c("a", "b", "c"), each = 10)))
    dd$y <- stats::rnorm(30, rep(c(0, .5, 1), each = 10))
    dd$w <- 1
    dd$w[c(1:4, 11:14, 21:24)] <- 0
    dd
  })
  fit <- anova_glm(d, "y", "g", weights = "w", plots = FALSE)
  ref <- suppressMessages(effectsize::omega_squared(
    car::Anova(stats::lm(y ~ g, data = d, weights = w)), partial = TRUE, ci = NULL))
  # one-way, so effectsize labels the partial measure plain "Omega2"
  expect_equal(fit$effect_sizes$partial_omega_sq, ref[[2]], tolerance = 1e-8)
  expect_identical(fit$n_removed, 12L)
})

test_that("F161: AIC and BIC with zero weights equal those of the deleted rows", {
  d <- withr::with_seed(21, {
    dd <- data.frame(g = rep(c("a", "b", "c"), each = 40))
    dd$y <- stats::rnorm(120, c(a = 0, b = 0.5, c = 1)[dd$g])
    dd$yb <- stats::rbinom(120, 1, c(a = .3, b = .5, c = .6)[dd$g])
    dd$w <- 1
    dd$w[sample(120, 30)] <- 0
    dd
  })
  del <- d[d$w > 0, ]
  mg <- stats::glm(y ~ g, data = del)
  mb <- stats::glm(yb ~ g, data = del, family = stats::binomial())
  fg <- anova_glm(d, "y", "g", weights = "w", plots = FALSE)
  fb <- anova_glm(d, "yb", "g", family = "binomial", weights = "w", plots = FALSE)
  expect_equal(fg$model_stats$AIC, stats::AIC(mg))
  expect_equal(fg$model_stats$BIC, stats::BIC(mg))
  expect_equal(fb$model_stats$BIC,
               -2 * as.numeric(stats::logLik(mb)) + log(90) * attr(stats::logLik(mb), "df"))
})

# F046 -------------------------------------------------------------------------

test_that("F046: the two-cell note claims equivalence only when t^2 = F", {
  d <- fx_oneway(means = c(A = 0, B = 1))
  fit <- anova_glm(d, "value", "group", plots = FALSE)
  ref <- stats::anova(stats::lm(value ~ group, data = d))
  expect_equal(fit$posthoc$statistic^2, ref$`F value`[1], tolerance = 1e-8)
  expect_true(any(grepl("is the omnibus test", fit$notes)))

  sep <- withr::with_seed(17, data.frame(
    g = factor(rep(c("a", "b"), c(40, 15))),
    y = c(stats::rbinom(40, 1, .3), rep(1, 15))))
  fb <- anova_glm(sep, "y", "g", family = "binomial", plots = FALSE)
  expect_false(any(grepl("is the omnibus test", fb$notes)))
  expect_true(any(grepl("can differ from the omnibus likelihood-ratio test", fb$notes)))

  skip_if_not_installed("sandwich")
  h <- withr::with_seed(8, {
    hh <- data.frame(g = factor(rep(c("a", "b"), c(40, 10))))
    hh$y <- stats::rnorm(50, c(a = 0, b = 1)[as.character(hh$g)],
                         c(a = .3, b = 3)[as.character(hh$g)])
    hh
  })
  fh <- anova_glm(h, "y", "g", vcov_type = "HC3", plots = FALSE)
  expect_false(any(grepl("is the omnibus test", fh$notes)))
  expect_true(any(grepl("robust \\(HC3\\) covariance and can differ", fh$notes)))
})

# F051 -------------------------------------------------------------------------

test_that("F051: a Gaussian fit whose ANOVA table fails never reports McFadden", {
  d <- withr::with_seed(11, {
    dd <- expand.grid(A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3"), rep = 1:6)
    dd <- dd[!(dd$A == "a3" & dd$B == "b3"), ]
    dd$y <- stats::rnorm(nrow(dd)) + (dd$A == "a2") + 0.5 * (dd$B == "b2")
    dd
  })
  # car refuses the same model, so there is no table to take sums of squares from
  m <- withr::with_options(list(contrasts = c("contr.sum", "contr.poly")),
                           stats::glm(y ~ A * B, data = d))
  expect_error(car::Anova(m, type = 3, test.statistic = "F"), "aliased")
  fit <- anova_glm(d, "y", c("A", "B"), interaction = TRUE, type = "III",
                   plots = FALSE, ci_method = "wald")
  expect_identical(nrow(fit$effect_sizes), 0L)
  expect_true("partial_eta_sq" %in% names(fit$effect_sizes))
  expect_false("measure" %in% names(fit$effect_sizes))
  expect_true(any(grepl("Partial eta and omega squared could not be computed", fit$notes)))
})

test_that("F051: a Gaussian LR or Wald test still gets partial eta squared", {
  d <- fx_oneway()
  fit <- anova_glm(d, "value", "group", test_statistic = "LR", plots = FALSE)
  ref <- stats::anova(stats::lm(value ~ group, data = d))
  expect_equal(fit$effect_sizes$partial_eta_sq,
               ref$`Sum Sq`[1] / sum(ref$`Sum Sq`))
})

# F089 / F135 ------------------------------------------------------------------

test_that("F089: negative.binomial(theta) defaults to the F test its t table implies", {
  skip_if_not_installed("MASS")
  d <- withr::with_seed(99, {
    g <- factor(rep(c("a", "b", "c"), each = 100))
    data.frame(g = g, cnt = stats::rnbinom(
      300, mu = c(a = 4, b = 5, c = 5.5)[as.character(g)], size = 1.5))
  })
  fit <- anova_glm(d, "cnt", "g", family = MASS::negative.binomial(1.5),
                   plots = FALSE, ci_method = "wald")
  m <- stats::glm(cnt ~ g, data = d, family = MASS::negative.binomial(1.5))
  ref <- car::Anova(m, test.statistic = "F")
  expect_identical(attr(fit$anova, "statistic"), "F")
  expect_equal(fit$anova$p_value[1], ref$`Pr(>F)`[1])
  expect_identical(unique(fit$assumptions$coefficients$distribution), "t")
})

test_that("F135: an invalid test_statistic is an error, not a fit without a table", {
  d <- fx_oneway()
  expect_error(anova_glm(d, "value", "group", test_statistic = "lr"),
               "`test_statistic` must be one of")
  expect_error(anova_glm(d, "value", "group", test_statistic = c("LR", "F")),
               "single character string")
})

# F155 -------------------------------------------------------------------------

test_that("F155: a group whose weights are all zero is not a group", {
  d <- withr::with_seed(5, data.frame(g = factor(rep(c("a", "b"), each = 20)),
                                      y = stats::rnorm(40)))
  d$w <- ifelse(d$g == "b", 0, 1)
  expect_error(anova_glm(d, "y", "g", weights = "w", plots = FALSE),
               "has 1 level")
  d3 <- rbind(d, withr::with_seed(6, data.frame(g = "c", y = stats::rnorm(20), w = 1)))
  fit <- anova_glm(d3, "y", "g", weights = "w", plots = FALSE)
  expect_identical(nrow(fit$data_used), 40L)
  expect_identical(fit$n_removed, 20L)
  expect_setequal(as.character(fit$emmeans$g), c("a", "c"))
  expect_equal(fit$anova$df[1], 1)
})

# F156 -------------------------------------------------------------------------

test_that("F156: separation and unbounded profile intervals are noted", {
  d <- glm_all_ones()
  fit <- anova_glm(d, "y", "g", family = "binomial", plots = FALSE)
  # stats::confint() is no reference here: depending on the platform its
  # spline returns NA or an extrapolated number for the side that never closes.
  # The finite end is checked against the deviance in the separation test below.
  # open above (Inf, not NA or a spline extrapolation), with a finite lower end
  expect_identical(fit$assumptions$coefficients$conf_high[
    fit$assumptions$coefficients$term == "gc"], Inf)
  expect_true(any(grepl("separation detected.*gc", fit$notes)))
  expect_true(any(grepl("estimate of gc has diverged", fit$notes)))

  qb <- anova_glm(d, "y", "g", family = "quasibinomial", plots = FALSE)
  expect_true(any(grepl("separation detected", qb$notes)))

  dp <- withr::with_seed(2, data.frame(g = factor(rep(c("a", "b", "c"), each = 10)),
                                       y = stats::rpois(30, 3)))
  dp$y[dp$g == "c"] <- 0
  fp <- anova_glm(dp, "y", "g", family = "poisson", plots = FALSE)
  expect_true(any(grepl("No events were observed in g = c", fp$notes)))
})

# F165 -------------------------------------------------------------------------

test_that("F165: when the marginal means fail, the notes say so and name the plot", {
  d <- fx_oneway()
  # emmeans refuses a grid with more rows than rg.limit
  fit <- withr::with_options(list(emmeans = list(rg.limit = 2)),
                             anova_glm(d, "value", "group"))
  expect_null(fit$emmeans)
  expect_null(fit$plots$emmeans)
  expect_false(any(grepl("fewer than two estimable cells", fit$notes)))
  expect_true(any(grepl("plot of estimated marginal means was skipped", fit$notes)))
  expect_true(any(grepl("Pairwise comparisons were not computed because", fit$notes)))
})

# F093 -------------------------------------------------------------------------

test_that("F093: a floored partial omega squared is noted", {
  d <- withr::with_seed(2, data.frame(g = factor(rep(c("A", "B", "C", "D"), each = 10)),
                                      y = stats::rnorm(40)))
  fit <- anova_glm(d, "y", "g", plots = FALSE, posthoc = FALSE)
  a <- stats::anova(stats::lm(y ~ g, data = d))
  raw <- (a$Df[1] * (a$`Mean Sq`[1] - a$`Mean Sq`[2])) /
    (a$Df[1] * a$`Mean Sq`[1] + (40 - a$Df[1]) * a$`Mean Sq`[2])
  expect_lt(raw, 0)
  expect_identical(fit$effect_sizes$partial_omega_sq, 0)
  expect_true(any(grepl("Partial omega squared was negative for g", fit$notes)))
})

# F047 / F088 ------------------------------------------------------------------

test_that("F047: a factor response with three levels is refused for binomial", {
  d <- withr::with_seed(13, data.frame(
    g = factor(rep(c("a", "b"), each = 60)),
    resp = factor(sample(c("low", "mid", "high"), 120, TRUE),
                  levels = c("low", "mid", "high"))))
  expect_error(anova_glm(d, "resp", "g", family = "binomial", plots = FALSE),
               "3 observed level")
  d2 <- droplevels(d[d$resp != "mid", ])
  fit <- anova_glm(d2, "resp", "g", family = "binomial", plots = FALSE)
  expect_true(any(grepl("equals \"high\"", fit$notes)))
  expect_equal(fit$emmeans$estimate,
               as.vector(tapply(d2$resp == "high", d2$g, mean)), tolerance = 1e-6)
})

test_that("F088: character responses are accepted, and factors for quasibinomial", {
  d <- withr::with_seed(15, data.frame(g = factor(rep(c("a", "b", "c"), each = 10)),
                                       y = stats::rpois(30, 3)))
  d$chr <- ifelse(d$y > 3, "hi", "lo")
  d$fac <- factor(d$chr)
  fc <- anova_glm(d, "chr", "g", family = "binomial", plots = FALSE)
  ref <- car::Anova(stats::glm(fac ~ g, binomial, data = d))
  expect_equal(fc$anova$statistic, ref$`LR Chisq`)
  fq <- anova_glm(d, "fac", "g", family = "quasibinomial", plots = FALSE)
  refq <- car::Anova(stats::glm(fac ~ g, quasibinomial, data = d),
                     test.statistic = "F")
  expect_equal(fq$anova$p_value[1], refq$`Pr(>F)`[1])
})

# F037 -------------------------------------------------------------------------

test_that("F037: quasi-binomial comparisons use t on the residual df", {
  d <- withr::with_seed(12, {
    dd <- data.frame(g = factor(rep(c("a", "b"), each = 12)))
    dd$y <- stats::rbinom(24, 1, ifelse(dd$g == "a", .3, .7))
    dd
  })
  fit <- anova_glm(d, "y", "g", family = "quasibinomial", plots = FALSE)
  ref <- stats::coef(summary(stats::glm(y ~ g, quasibinomial, data = d)))
  expect_equal(fit$posthoc$df, 22)
  expect_equal(fit$posthoc$p_value, ref["gb", 4], tolerance = 1e-8)
})

# F152 / F157 / F013 -------------------------------------------------------------

test_that("F152: robust SEs with zero weights equal those on the positive rows", {
  skip_if_not_installed("sandwich")
  d <- withr::with_seed(8, {
    dd <- data.frame(g = factor(rep(c("a", "b", "c"), each = 40)))
    dd$y <- stats::rnorm(120, c(a = 0, b = .4, c = .8)[as.character(dd$g)],
                         sd = c(a = 1, b = 1.5, c = 2)[as.character(dd$g)])
    dd$w <- rep(c(1, 0), 60)
    dd
  })
  fit <- anova_glm(d, "y", "g", weights = "w", vcov_type = "HC3", plots = FALSE)
  ref <- sqrt(diag(sandwich::vcovHC(stats::glm(y ~ g, data = d[d$w > 0, ]), "HC3")))
  expect_equal(fit$assumptions$coefficients$se, unname(ref), tolerance = 1e-8)
})

test_that("F157: frequency weights with robust SEs are flagged", {
  skip_if_not_installed("sandwich")
  d <- as.data.frame(UCBAdmissions)
  fit <- anova_glm(d, "Admit", "Gender", family = "binomial", weights = "Freq",
                   vcov_type = "HC0", plots = FALSE)
  expect_true(any(grepl("frequency weights", fit$notes)))
  # a binary response cannot be overdispersed: no Pearson-dispersion note
  expect_false(any(grepl("^Pearson dispersion is", fit$notes)))
  expect_true(is.na(fit$assumptions$dispersion))
})

test_that("F013: the omnibus test says which covariance it uses", {
  skip_if_not_installed("sandwich")
  d <- fx_oneway(sds = c(1, 2, 4))
  fw <- anova_glm(d, "value", "group", vcov_type = "HC3", test_statistic = "Wald",
                  plots = FALSE)
  m <- stats::glm(value ~ group, data = d)
  ref <- suppressMessages(
    car::Anova(m, test.statistic = "Wald", vcov. = sandwich::vcovHC(m, "HC3")))
  expect_equal(fw$anova$statistic, ref$Chisq)
  expect_true(any(grepl("omnibus Wald test uses the robust \\(HC3\\)", fw$notes)))
  expect_false(any(grepl("Coefficient covariances computed", fw$notes)))
  ff <- anova_glm(d, "value", "group", vcov_type = "HC3", plots = FALSE)
  expect_equal(ff$anova$statistic[1],
               stats::anova(stats::lm(value ~ group, data = d))$`F value`[1])
  expect_true(any(grepl("omnibus F test compares deviances", ff$notes)))
})

# F096 -------------------------------------------------------------------------

test_that("F096: the fit does not carry unanalysed input columns", {
  wide <- withr::with_seed(1, {
    w <- data.frame(g = factor(sample(letters[1:3], 2000, TRUE)), y = stats::rnorm(2000))
    for (j in 1:20) w[[paste0("junk", j)]] <- stats::rnorm(2000)
    w
  })
  fit <- anova_glm(wide, "y", "g", plots = FALSE)
  narrow <- anova_glm(wide[c("g", "y")], "y", "g", plots = FALSE)
  env <- environment(fit$model$formula)
  expect_false(exists("data", envir = env, inherits = FALSE))
  expect_identical(names(env$data_used), c("y", "g"))
  size <- function(x) length(serialize(x, NULL))
  expect_lt(size(fit), 1.05 * size(narrow))
  # the model can still be refitted from what it carries
  expect_equal(stats::coef(stats::update(fit$model, data = fit$data_used)),
               stats::coef(fit$model))
})

# Dispersion notes shared with anova_count (F146, F151, F154) ------------------

test_that("a Poisson frequency table gets the dispersion of the expanded data", {
  withr::local_seed(21)
  long <- data.frame(g = factor(rep(c("a", "b", "c"), each = 200)))
  long$y <- stats::rpois(600, rep(c(2, 3, 4), each = 200))
  tab <- stats::aggregate(list(n = rep(1, 600)), by = list(g = long$g, y = long$y), FUN = sum)
  f_long <- anova_glm(long, "y", "g", family = "poisson", plots = FALSE)
  f_tab <- anova_glm(tab, "y", "g", family = "poisson", weights = "n", plots = FALSE)
  ref <- sum(stats::residuals(f_long$model, type = "pearson")^2) /
    stats::df.residual(f_long$model)
  expect_equal(f_tab$assumptions$dispersion, ref, tolerance = 1e-8)
  expect_false(any(grepl("^Pearson dispersion is", f_tab$notes)))
  expect_true(any(grepl("frequency weights", f_tab$notes)))
})

test_that("mild Poisson overdispersion is reported with its consequence", {
  withr::local_seed(22)
  d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 150)))
  d$y <- stats::rnbinom(450, mu = 6, size = 12)
  fit <- anova_glm(d, "y", "g", family = "poisson", plots = FALSE)
  phi <- sum(stats::residuals(fit$model, type = "pearson")^2) / stats::df.residual(fit$model)
  expect_true(phi > 1.2 && phi < 1.5)
  size <- stats::pchisq(stats::qchisq(0.95, 1) / phi, 1, lower.tail = FALSE)
  note <- grep("^Pearson dispersion is", fit$notes, value = TRUE)
  expect_length(note, 1L)
  expect_match(note, sprintf("roughly %.0f%%", 100 * size), fixed = TRUE)
})

test_that("an F test on a Poisson fit says it is a quasi-likelihood test", {
  d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 20)),
                  y = rep(c(1, 3, 6, 2, 8, 0), 10))
  fit <- anova_glm(d, "y", "g", family = "poisson", test_statistic = "F",
                   plots = FALSE)
  ref <- car::Anova(stats::glm(y ~ g, data = d, family = poisson), test.statistic = "F")
  expect_equal(fit$anova$statistic[1], ref[["F value"]][1])
  expect_true(any(grepl("quasi-likelihood test", fit$notes)))
  expect_true(any(grepl(sprintf("%.3g", sum(stats::residuals(fit$model, type = "pearson")^2) /
                                   stats::df.residual(fit$model)), fit$notes, fixed = TRUE)))
})

# Separation, dispersion of 0/1 data (F045 and F032, glm parts) ----------------

test_that("a separated logistic coefficient has an open side and a true profile end", {
  withr::local_seed(106)
  d <- data.frame(g = factor(c(rep("a", 100), rep("b", 100), rep("c", 12))),
                  y = c(stats::rbinom(100, 1, 0.3), stats::rbinom(100, 1, 0.5), rep(1, 12)))
  fit <- anova_glm(d, "y", "g", family = "binomial", plots = FALSE)
  row <- fit$assumptions$coefficients[fit$assumptions$coefficients$term == "gc", ]
  expect_identical(row$conf_high, Inf)
  expect_true(is.finite(row$conf_low))
  # Independent check: fixing gc at the lower end raises the deviance by the
  # chi-square cut-off.
  X <- stats::model.matrix(~ g, d)
  full <- stats::glm(y ~ g, data = d, family = binomial)
  fixed <- stats::glm(d$y ~ 0 + X[, c("(Intercept)", "gb")],
                      offset = row$conf_low * X[, "gc"], family = binomial)
  expect_equal(stats::deviance(fixed) - stats::deviance(full),
               stats::qchisq(0.95, 1), tolerance = 1e-4)
  expect_true(any(grepl("diverged", fit$notes)))
})

test_that("a zero-count Poisson cell gets an interval open below", {
  withr::local_seed(6)
  d <- data.frame(g = factor(rep(c("A", "B", "C"), each = 20)))
  d$y <- c(stats::rpois(20, 3), stats::rpois(20, 5), rep(0, 20))
  fit <- anova_glm(d, "y", "g", family = "poisson", plots = FALSE)
  row <- fit$assumptions$coefficients[fit$assumptions$coefficients$term == "gC", ]
  expect_identical(row$conf_low, -Inf)
  expect_true(is.na(row$conf_high) || row$conf_high > row$estimate)
  expect_true(any(grepl("No events were observed in g = C", fit$notes)))
})

test_that("the dispersion of a 0/1 response is NA, with a reason", {
  d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 30)),
                  y = rep(c(0, 1, 1, 0, 1), 18))
  fit <- anova_glm(d, "y", "g", family = "binomial", plots = FALSE)
  expect_true(is.na(fit$assumptions$dispersion))
  expect_true(any(grepl("dispersion is reported as NA", fit$notes)))
})
