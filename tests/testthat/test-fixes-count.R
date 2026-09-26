# Regression tests for the anova_count() findings of the adversarial review.
#
# Each test would fail on the code before the fix and checks the package
# against an independent reference: stats::glm(), MASS::glm.nb(), car::Anova(),
# emmeans, sandwich, or a hand computation from the published formula.

# Local fixtures ---------------------------------------------------------------

fc_zero_group <- function() {
  withr::with_seed(6, {
    d <- data.frame(g = factor(rep(c("A", "B", "C"), each = 20)))
    d$y <- c(stats::rnbinom(20, mu = 3, size = 0.8),
             stats::rnbinom(20, mu = 5, size = 0.8), rep(0, 20))
    d
  })
}

fc_offset <- function() {
  withr::with_seed(1, {
    d <- data.frame(g = rep(c("A", "B", "C"), each = 40), stringsAsFactors = FALSE)
    d$hours <- stats::runif(120, 0.5, 4)
    d$w <- stats::runif(120, 0.5, 2)
    d$y <- stats::rpois(120, c(2, 4, 6)[match(d$g, c("A", "B", "C"))] * d$hours)
    d
  })
}

fc_small_nb <- function(n = 8, seed = 5) {
  withr::with_seed(seed, data.frame(
    g = factor(rep(c("a", "b", "c"), each = n)),
    y = stats::rnbinom(3 * n, mu = 3, size = 1)))
}

fc_freq <- function() {
  raw <- withr::with_seed(10, {
    d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 200)))
    d$cnt <- stats::rpois(600, c(a = 2, b = 2.5, c = 3)[as.character(d$g)])
    d
  })
  agg <- stats::aggregate(list(freq = rep(1, 600)),
                          by = list(g = raw$g, cnt = raw$cnt), FUN = sum)
  list(raw = raw, agg = agg)
}

pearson_dispersion <- function(m) {
  sum(stats::residuals(m, type = "pearson")^2) / stats::df.residual(m)
}

# F034: all-zero groups and all-zero responses --------------------------------

test_that("F034: an all-zero group is flagged on the default (negative binomial) path", {
  skip_if_not_installed("MASS")
  d <- fc_zero_group()
  fit <- anova_count(d, "y", "g", plots = FALSE)
  expect_identical(fit$model_type, "negbin")
  ref <- car::Anova(MASS::glm.nb(y ~ g, data = d), test.statistic = "LR")
  expect_equal(fit$anova$statistic, ref$`LR Chisq`, tolerance = 1e-6)
  expect_true(any(grepl("No events were observed in g = C", fit$notes)))
  expect_true(any(grepl("Wald p-value near 1", fit$notes)))
})

test_that("F034: a response with no events at all is refused", {
  z <- data.frame(g = factor(rep(c("A", "B"), each = 10)), y = 0)
  expect_error(anova_count(z, "y", "g", plots = FALSE), "zero in every row")
  expect_error(anova_count(z, "y", "g", model = "negbin", plots = FALSE),
               "zero in every row")
  # ... also when the only events sit in rows that are not analysed
  z$y[1] <- 4; z$w <- c(0, rep(1, 19))
  expect_error(anova_count(z, "y", "g", weights = "w", plots = FALSE),
               "zero in every row")
})

test_that("F034: an explicit model = 'negbin' that cannot be fitted is an error, not a substitution", {
  skip_if_not_installed("MASS")
  agg <- data.frame(region = c("North", "South", "East"),
                    cases = c(30, 52, 81), pyears = c(10000, 12000, 15000))
  # glm.nb() itself cannot fit a saturated model
  expect_error(suppressWarnings(MASS::glm.nb(
    cases ~ region + offset(log(pyears)), data = agg)))
  expect_error(anova_count(agg, "cases", "region", offset = "pyears",
                           model = "negbin", plots = FALSE),
               "negative binomial model could not be fitted")
})

# F036: offset models rebuild their marginal means from the analysed rows ----

test_that("F036: offset plus weights gives the marginal rates of the weighted model", {
  d <- fc_offset()
  fit <- anova_count(d, "y", "g", offset = "hours", weights = "w",
                     model = "poisson", plots = FALSE)
  d$g <- factor(d$g)
  ref <- stats::glm(y ~ g + offset(log(hours)), family = stats::poisson(),
                    data = d, weights = w)
  emm <- as.data.frame(emmeans::emmeans(ref, ~ g, type = "response", offset = 0))
  expect_equal(fit$emmeans$estimate, emm$rate, tolerance = 1e-8)
  expect_identical(nrow(fit$posthoc), 3L)
})

test_that("F036: a numeric-coded group with an offset keeps its marginal rates", {
  d <- fc_offset()
  d$site <- match(d$g, c("A", "B", "C"))
  fit <- anova_count(d, "y", "site", offset = "hours", model = "poisson",
                     plots = FALSE)
  expect_false(any(grepl("could not be computed", fit$notes)))
  ref <- stats::glm(y ~ factor(site) + offset(log(hours)),
                    family = stats::poisson(), data = d)
  rates <- exp(stats::coef(ref)[1] + c(0, stats::coef(ref)[-1]))
  expect_equal(fit$emmeans$estimate, unname(rates), tolerance = 1e-8)
  expect_equal(fit$emmeans$estimate, c(1.8146, 3.9091, 5.8380), tolerance = 1e-4)
})

test_that("F036: a row dropped for an infinite count does not come back as a group", {
  d <- fc_offset()
  d2 <- rbind(d[, c("g", "hours", "y")], data.frame(g = "D", hours = 1, y = Inf))
  fit <- anova_count(d2, "y", "g", offset = "hours", model = "poisson",
                     plots = FALSE)
  expect_identical(fit$n_removed, 1L)
  expect_identical(as.character(fit$emmeans$g), c("A", "B", "C"))
  expect_false(any(grepl("D", fit$posthoc$contrast)))
})

test_that("F036: a character group in a non-C collation keeps its rates with an offset", {
  old <- Sys.getlocale("LC_COLLATE")
  on.exit(suppressWarnings(Sys.setlocale("LC_COLLATE", old)), add = TRUE)
  ok <- suppressWarnings(Sys.setlocale("LC_COLLATE", "en_US.UTF-8"))
  skip_if(!nzchar(ok), "the en_US.UTF-8 collation is not available")
  d <- fc_offset()
  d$h <- ifelse(d$g == "A", "a", d$g)
  skip_if(identical(sort(c("a", "B")), c("B", "a")), "collation did not change")
  fit <- anova_count(d, "y", "h", offset = "hours", model = "poisson",
                     plots = FALSE)
  expect_false(any(grepl("could not be computed", fit$notes)))
  expect_identical(nrow(fit$emmeans), 3L)
  byg <- stats::setNames(fit$emmeans$estimate, as.character(fit$emmeans$h))
  expect_equal(unname(byg[c("a", "B", "C")]), c(1.8146, 3.9091, 5.8380),
               tolerance = 1e-4)
})

# F082, F147, F144, F085: the negative binomial notes --------------------------

test_that("F082: a small, imprecise theta is not called weak overdispersion", {
  skip_if_not_installed("MASS")
  d <- data.frame(g = factor(rep(c("A", "B"), each = 6)),
                  y = c(0, 4, 5, 0, 1, 5, 2, 0, 1, 6, 0, 21))
  fit <- anova_count(d, "y", "g", plots = FALSE)
  ref <- MASS::glm.nb(y ~ g, data = d)
  expect_identical(fit$model_type, "negbin")
  note <- grep("theta =", fit$notes, value = TRUE)
  expect_length(note, 1L)
  expect_match(note, sprintf("theta = %.3f (SE %.3f)", ref$theta, ref$SE.theta),
               fixed = TRUE)
  expect_match(note, "poorly determined")
  expect_false(grepl("weakly overdispersed|mildly", note))
})

test_that("F082: 'weakly overdispersed' is kept for a large, imprecise theta", {
  skip_if_not_installed("MASS")
  agg <- fc_freq()$agg
  fit <- anova_count(agg, "cnt", "g", weights = "freq", model = "negbin",
                     plots = FALSE)
  ref <- MASS::glm.nb(cnt ~ g, data = agg, weights = freq)
  expect_gt(ref$theta, 10)
  expect_lt(ref$theta / ref$SE.theta, 2)
  note <- grep("theta =", fit$notes, value = TRUE)
  expect_match(note, sprintf("theta = %.3f", ref$theta), fixed = TRUE)
  expect_match(note, "weakly overdispersed")
})

test_that("F147: the theta note states the small-sample anti-conservatism for the smallest group", {
  skip_if_not_installed("MASS")
  d <- fc_small_nb(n = 8)
  fit <- anova_count(d, "y", "g", plots = FALSE)
  expect_identical(fit$model_type, "negbin")
  note <- grep("theta =", fit$notes, value = TRUE)
  expect_match(note, sprintf("With %d observations in the smallest group",
                             min(table(d$g))), fixed = TRUE)
  expect_match(note, "markedly anti-conservative")
  expect_match(note, "10-15%", fixed = TRUE)
  expect_match(note, "intervals")
  expect_false(grepl("mildly", note))

  # With 100 per group the excess is small, and the note says so
  big <- anova_count(fx_overdispersed(), "y", "g", plots = FALSE)
  expect_match(grep("theta =", big$notes, value = TRUE), "mildly anti-conservative")
})

test_that("F144: non-convergence is read from the fit, in any session language", {
  skip_if_not_installed("MASS")
  d <- fx_overdispersed(n_per = 60, seed = 77)
  real <- MASS::glm.nb
  german <- "Grenze der Alternierungen erreicht"   # MASS's de translation
  notes <- testthat::with_mocked_bindings(
    anova_count(d, "y", "g", model = "negbin", plots = FALSE)$notes,
    glm.nb = function(...) {
      m <- real(...)
      warning(german)
      m$th.warn <- german
      m
    },
    .package = "MASS")
  expect_true(any(grepl("did not converge", notes)))
  expect_true(any(grepl(german, notes, fixed = TRUE)))
  expect_false(any(grepl("theta = [0-9.]+ \\(SE", notes)))
  # the translated warning is the same event, not a second note
  expect_false(any(grepl("MASS::glm.nb() reported", notes, fixed = TRUE)))

  # glm.nb() fits whose final IRLS step did not converge are reported too
  notes2 <- testthat::with_mocked_bindings(
    anova_count(d, "y", "g", model = "negbin", plots = FALSE)$notes,
    glm.nb = function(...) { m <- real(...); m$converged <- FALSE; m },
    .package = "MASS")
  expect_true(any(grepl("did not converge", notes2)))
})

test_that("F144: a real alternation-limit fit gets the non-convergence note", {
  skip_if_not_installed("MASS")
  d <- withr::with_seed(630, {
    k <- sample(2:3, 1); n <- sample(8:40, 1)
    g <- factor(rep(letters[1:k], each = n))
    mu <- exp(stats::rnorm(k, 1, 1))[as.integer(g)]
    expo <- exp(stats::rnorm(length(g), 0, stats::runif(1, 2, 4)))
    data.frame(g = g, expo = expo,
               y = stats::rnbinom(length(g), mu = mu * expo, size = stats::runif(1, 1, 20)))
  })
  ref <- suppressWarnings(MASS::glm.nb(y ~ g + offset(log(expo)), data = d))
  expect_false(is.null(ref$th.warn))
  fit <- anova_count(d, "y", "g", offset = "expo", plots = FALSE)
  expect_identical(fit$model_type, "negbin")
  expect_true(any(grepl("did not converge", fit$notes)))
})

test_that("F085: other glm.nb() warnings reach $notes, on both branches", {
  skip_if_not_installed("MASS")
  d <- fx_overdispersed(n_per = 60, seed = 77)
  real <- MASS::glm.nb
  notes <- testthat::with_mocked_bindings(
    anova_count(d, "y", "g", model = "negbin", plots = FALSE)$notes,
    glm.nb = function(...) {
      m <- real(...)
      warning("glm.fit: algorithm did not converge")
      m
    },
    .package = "MASS")
  expect_true(any(grepl("MASS::glm.nb() reported: glm.fit: algorithm did not converge",
                        notes, fixed = TRUE)))

  # model = "auto" falls back to quasi-Poisson when glm.nb() fails, and keeps
  # what it said on the way
  fb <- testthat::with_mocked_bindings(
    anova_count(d, "y", "g", plots = FALSE),
    glm.nb = function(...) { warning("theta went astray"); stop("no fit") },
    .package = "MASS")
  expect_identical(fb$model_type, "quasipoisson")
  expect_true(any(grepl("The negative binomial fit failed (no fit)", fb$notes,
                        fixed = TRUE)))
  expect_true(any(grepl("MASS::glm.nb() reported: theta went astray", fb$notes,
                        fixed = TRUE)))
})

# F146, F039, F084: the dispersion notes ---------------------------------------

test_that("F146: a Poisson fit kept below the threshold states what its dispersion costs", {
  d <- withr::with_seed(3, data.frame(
    g = factor(rep(c("a", "b", "c"), each = 100)),
    y = stats::rnbinom(300, mu = 4, size = 10)))
  fit <- anova_count(d, "y", "g", plots = FALSE)
  ref <- stats::glm(y ~ g, family = stats::poisson(), data = d)
  phi <- pearson_dispersion(ref)
  expect_gt(phi, 1); expect_lte(phi, 1.5)
  expect_identical(fit$model_type, "poisson")
  note <- grep("Poisson model was kept", fit$notes, value = TRUE)
  expect_length(note, 1L)
  # Poisson LR under dispersion phi is about phi * chisq: the size of a
  # nominal 5% test on 1 df is P(chisq_1 > 3.84 / phi)
  size <- stats::pchisq(stats::qchisq(0.95, 1) / phi, 1, lower.tail = FALSE)
  expect_match(note, sprintf("inflated by about that factor"), fixed = TRUE)
  expect_match(note, sprintf("at %.2f", phi), fixed = TRUE)
  expect_match(note, sprintf("roughly %.0f%%", 100 * size), fixed = TRUE)
  expect_match(note, sprintf("about %.2f times too small", sqrt(phi)), fixed = TRUE)
  expect_match(note, "quasipoisson")

  # below 1 there is nothing to warn about
  under <- withr::with_seed(4, data.frame(g = factor(rep(c("a", "b"), each = 50)),
                                          y = stats::rbinom(100, 8, 0.5)))
  uf <- anova_count(under, "y", "g", plots = FALSE)
  expect_lt(uf$dispersion, 1)
  expect_false(any(grepl("inflated", uf$notes)))
})

test_that("F039: an explicit Poisson model on overdispersed data gets the overdispersion note", {
  d <- withr::with_seed(2, {
    dd <- data.frame(A = sample(c("a1", "a2", "a3"), 200, TRUE))
    dd$cnt <- stats::rnbinom(200, mu = exp(1 + 0.4 * (dd$A == "a2")), size = 2)
    dd
  })
  fit <- anova_count(d, "cnt", "A", model = "poisson", plots = FALSE)
  ref <- stats::glm(cnt ~ A, family = stats::poisson(), data = d)
  expect_identical(fit$model_type, "poisson")
  expect_true(any(grepl(sprintf("Pearson dispersion is %.2f, above the threshold",
                                pearson_dispersion(ref)), fit$notes)))
  expect_true(any(grepl("overdispersed", fit$notes)))
})

test_that("F084: model = 'auto' with no residual df says overdispersion was not checked", {
  d <- data.frame(g = factor(c("A", "B", "C")), y = c(3, 7, 12))
  fit <- anova_count(d, "y", "g", plots = FALSE)
  expect_identical(stats::df.residual(stats::glm(y ~ g, stats::poisson(), d)), 0L)
  expect_identical(fit$model_type, "poisson")
  expect_true(is.na(fit$dispersion))
  expect_true(any(grepl("overdispersion could not be checked", fit$notes)))
})

# F151: frequency weights --------------------------------------------------------

test_that("F151: a frequency table and its expanded data make the same model choice", {
  fr <- fc_freq()
  raw <- anova_count(fr$raw, "cnt", "g", plots = FALSE)
  tab <- anova_count(fr$agg, "cnt", "g", weights = "freq", plots = FALSE)
  ref <- stats::glm(cnt ~ g, family = stats::poisson(), data = fr$raw)
  expect_equal(tab$dispersion, pearson_dispersion(ref))
  expect_equal(raw$dispersion, tab$dispersion)
  expect_identical(tab$model_type, "poisson")
  expect_equal(tab$anova$statistic,
               car::Anova(ref, test.statistic = "LR")$`LR Chisq`)
  expect_true(any(grepl("treated as frequency weights", tab$notes)))
  expect_true(any(grepl("600 observations", tab$notes)))

  # quasi-Poisson still counts rows, and says so
  q <- anova_count(fr$agg, "cnt", "g", weights = "freq", model = "quasipoisson",
                   plots = FALSE)
  expect_true(any(grepl("estimates its dispersion from the 27 rows", q$notes)))

  # non-integer weights are not frequencies
  w <- fr$agg; w$freq <- w$freq + 0.5
  nw <- anova_count(w, "cnt", "g", weights = "freq", model = "poisson",
                    plots = FALSE)
  expect_equal(nw$dispersion,
               pearson_dispersion(stats::glm(cnt ~ g, stats::poisson(), w,
                                             weights = freq)))
  expect_false(any(grepl("frequency weights", nw$notes)))
})

# F154: an F test on a fixed-dispersion model -----------------------------------

test_that("F154: test_statistic = 'F' on a Poisson or NB model says what the F test estimated", {
  d <- withr::with_seed(2, {
    dd <- data.frame(g = factor(rep(c("a", "b", "c"), each = 20)))
    dd$cnt <- stats::rnbinom(60, mu = c(a = 5, b = 5, c = 7.5)[as.character(dd$g)],
                             size = 1.5)
    dd
  })
  fit <- anova_count(d, "cnt", "g", model = "poisson", test_statistic = "F",
                     plots = FALSE)
  qref <- car::Anova(stats::glm(cnt ~ g, stats::quasipoisson(), d),
                     test.statistic = "F")
  expect_equal(fit$anova$p_value[1], qref$`Pr(>F)`[1])
  pref <- stats::glm(cnt ~ g, stats::poisson(), d)
  note <- grep("test_statistic = \"F\"", fit$notes, value = TRUE, fixed = TRUE)
  expect_length(note, 1L)
  expect_match(note, sprintf("(%.2f)", pearson_dispersion(pref)), fixed = TRUE)
  expect_match(note, "$effect_sizes and $posthoc still assume a dispersion of 1",
               fixed = TRUE)

  skip_if_not_installed("MASS")
  nb <- anova_count(d, "cnt", "g", model = "negbin", test_statistic = "F",
                    plots = FALSE)
  nref <- MASS::glm.nb(cnt ~ g, data = d)
  expect_true(any(grepl(sprintf("negative binomial model, so car::Anova() estimated a dispersion from the Pearson residuals (%.2f)",
                                pearson_dispersion(nref)), nb$notes, fixed = TRUE)))
  # a quasi-Poisson F test is what was asked for: no note
  qp <- anova_count(d, "cnt", "g", model = "quasipoisson", test_statistic = "F",
                    plots = FALSE)
  expect_false(any(grepl("test_statistic = \"F\" was requested", qp$notes, fixed = TRUE)))
})

# F038, F083: empty and sparse cells ------------------------------------------

test_that("F038: an empty cell of an additive model is not called inestimable", {
  d <- withr::with_seed(18, {
    dd <- expand.grid(a = c("p", "q"), b = c("x", "y", "z"), rep = 1:15)
    dd <- dd[!(dd$a == "q" & dd$b == "z"), ]
    dd$y <- stats::rpois(nrow(dd), 5)
    dd
  })
  fit <- anova_count(d, "y", c("a", "b"), model = "poisson", plots = FALSE)
  ref <- stats::glm(y ~ a + b, family = stats::poisson(), data = d)
  cell <- as.data.frame(emmeans::emmeans(ref, ~ a:b, type = "response"))
  expect_true(all(is.finite(cell$rate)))          # q:z is estimable
  note <- grep("contain no observations", fit$notes, value = TRUE)
  expect_length(note, 1L)
  expect_false(grepl("not estimable", note))
  expect_match(note, "additive model still estimates them")

  full <- anova_count(d, "y", c("a", "b"), interaction = TRUE, model = "poisson",
                      plots = FALSE)
  expect_true(is.na(stats::coef(stats::glm(y ~ a * b, stats::poisson(), d))["aq:bz"]))
  expect_true(any(grepl("contain no observations; contrasts involving them are not estimable",
                        full$notes)))
})

test_that("F083: aggregated rate data are judged on events, not rows", {
  agg <- data.frame(region = c("North", "South", "East"), cases = c(30, 52, 81),
                    pyears = c(10000, 12000, 15000))
  fit <- anova_count(agg, "cases", "region", offset = "pyears", plots = FALSE)
  m1 <- stats::glm(cases ~ region + offset(log(pyears)), agg, family = stats::poisson())
  m0 <- stats::glm(cases ~ offset(log(pyears)), agg, family = stats::poisson())
  expect_equal(fit$anova$statistic, stats::anova(m0, m1, test = "Chisq")$Deviance[2])
  expect_false(any(grepl("undefined|eta squared", fit$notes)))
  expect_false(any(grepl("fewer than 5", fit$notes)))
  expect_true(any(grepl("still valid", fit$notes)))

  # many rows but very few events is what makes a cell sparse
  d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 50)),
                  y = c(rep(c(3, 4), 25), rep(c(5, 6), 25), c(1, 2, rep(0, 48))))
  events <- tapply(d$y, d$g, sum)
  sp <- anova_count(d, "y", "g", model = "poisson", plots = FALSE)
  expect_true(any(grepl(sprintf("^%d group combination\\(s\\) hold fewer than 5 events",
                                sum(events < 5)), sp$notes)))
})

# F086: the offset's role -------------------------------------------------------

test_that("F086: the offset cannot also be the response or a grouping column", {
  d <- withr::with_seed(34, {
    dd <- data.frame(g = factor(rep(c("A", "B"), each = 30)),
                     hours = sample(1:3, 60, TRUE))
    dd$y <- stats::rpois(60, 3 * dd$hours) + 1
    dd
  })
  expect_error(anova_count(d, "y", c("g", "hours"), offset = "hours", plots = FALSE),
               "both a grouping variable and the exposure offset")
  expect_error(anova_count(d, "y", "g", offset = "y", plots = FALSE),
               "both the response and the exposure offset")
})

# F030 (IRR part): rate ratios against the reference level ---------------------

test_that("F030: IRRs are level-vs-reference ratios under type III and global sum contrasts", {
  d <- withr::with_seed(1, {
    dd <- data.frame(g = factor(sample(c("a", "b", "c"), 300, TRUE)))
    dd$cnt <- stats::rpois(300, c(a = 1, b = 2, c = 4)[as.character(dd$g)])
    dd
  })
  ref <- stats::glm(cnt ~ g, family = stats::poisson(), data = d,
                    contrasts = list(g = "contr.treatment"))
  rc <- stats::coef(summary(ref))[-1, ]
  obs <- as.vector(tapply(d$cnt, d$g, mean))         # hand check
  t3 <- anova_count(d, "cnt", "g", type = "III", model = "poisson", plots = FALSE)
  expect_true(all(unlist(t3$model$contrasts) == "contr.sum"))
  expect_identical(t3$effect_sizes$term, c("gb", "gc"))
  expect_identical(t3$effect_sizes$comparison, c("b vs a", "c vs a"))
  expect_identical(t3$effect_sizes$factor, c("g", "g"))
  expect_equal(t3$effect_sizes$IRR, unname(exp(rc[, 1])))
  expect_equal(t3$effect_sizes$IRR, obs[-1] / obs[1])
  expect_equal(t3$effect_sizes$conf_low,
               unname(exp(rc[, 1] - stats::qnorm(0.975) * rc[, 2])))
  expect_equal(t3$effect_sizes$p_value, unname(rc[, 4]))

  gl <- withr::with_options(list(contrasts = c("contr.sum", "contr.poly")),
                            anova_count(d, "cnt", "g", model = "poisson", plots = FALSE))
  expect_equal(gl$effect_sizes$IRR, unname(exp(rc[, 1])))

  # an ordered factor gets ratios against its first level, not polynomial terms
  d$dose <- factor(d$g, ordered = TRUE)
  od <- anova_count(d, "cnt", "dose", model = "poisson", plots = FALSE)
  expect_equal(od$effect_sizes$IRR, unname(exp(rc[, 1])))
  expect_identical(od$effect_sizes$comparison, c("b vs a", "c vs a"))
})

test_that("F030: NB, quasi-Poisson, robust and interaction IRRs match a treatment-coded refit", {
  skip_if_not_installed("MASS")
  d <- fc_small_nb(n = 30, seed = 9)
  nb <- anova_count(d, "y", "g", model = "negbin", type = "III", plots = FALSE)
  nref <- stats::coef(summary(MASS::glm.nb(y ~ g, data = d)))[-1, ]
  expect_equal(nb$effect_sizes$IRR, unname(exp(nref[, 1])), tolerance = 1e-6)
  expect_equal(nb$effect_sizes$p_value, unname(nref[, 4]), tolerance = 1e-5)

  qp <- anova_count(d, "y", "g", model = "quasipoisson", type = "III", plots = FALSE)
  qfit <- stats::glm(y ~ g, family = stats::quasipoisson(), data = d)
  qref <- stats::coef(summary(qfit))[-1, ]
  expect_equal(qp$effect_sizes$p_value, unname(qref[, 4]))
  expect_equal(qp$effect_sizes$conf_high,
               unname(exp(qref[, 1] + stats::qt(0.975, stats::df.residual(qfit)) * qref[, 2])))

  skip_if_not_installed("sandwich")
  hc <- anova_count(d, "y", "g", model = "poisson", type = "III",
                    vcov_type = "HC3", plots = FALSE)
  pfit <- stats::glm(y ~ g, family = stats::poisson(), data = d)
  se <- sqrt(diag(sandwich::vcovHC(pfit, type = "HC3")))[-1]
  expect_equal(hc$effect_sizes$conf_low,
               unname(exp(stats::coef(pfit)[-1] - stats::qnorm(0.975) * se)))

  e <- withr::with_seed(18, {
    dd <- expand.grid(a = c("p", "q"), b = c("x", "y", "z"), rep = 1:15)
    dd <- dd[!(dd$a == "q" & dd$b == "z"), ]
    dd$y <- stats::rpois(nrow(dd), 5)
    dd
  })
  ix <- anova_count(e, "y", c("a", "b"), interaction = TRUE, type = "III",
                    model = "poisson", plots = FALSE)
  iref <- stats::coef(stats::glm(y ~ a * b, family = stats::poisson(), data = e))[-1]
  expect_identical(ix$effect_sizes$term, names(iref))
  expect_equal(ix$effect_sizes$IRR, unname(exp(iref)))     # NA for aq:bz
  expect_identical(ix$effect_sizes$comparison[1], "q vs p at b = x")
  expect_identical(ix$effect_sizes$comparison[4], "(q vs p) x (y vs x)")
  expect_identical(ix$effect_sizes$factor[4], "a:b")
})
