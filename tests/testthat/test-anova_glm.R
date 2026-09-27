# anova_glm(): the general wrapper, for any family.
#
# Checks that a Gaussian fit reproduces stats::anova, that type = "III" fits
# under sum-to-zero contrasts rather than relabelling a Type I table, that the
# reference distribution follows whether the family estimates its dispersion,
# and that family may be given as an object, a function or a string.

test_that("a Gaussian fit reproduces stats::anova on the same model", {
  d <- fx_oneway()
  fit <- anova_glm(d, "value", "group", plots = FALSE)
  ref <- stats::anova(stats::lm(value ~ group, data = d))

  expect_equal(fit$anova$statistic[1], ref$`F value`[1])
  expect_equal(fit$anova$p_value[1], ref$`Pr(>F)`[1])
  expect_equal(fit$anova$sum_sq[1], ref$`Sum Sq`[1])
})

test_that("partial eta squared matches its definition and has no residual row", {
  d <- fx_oneway()
  fit <- anova_glm(d, "value", "group", plots = FALSE)
  ref <- stats::anova(stats::lm(value ~ group, data = d))
  expected <- ref$`Sum Sq`[1] / (ref$`Sum Sq`[1] + ref$`Sum Sq`[2])

  expect_equal(fit$effect_sizes$partial_eta_sq, expected)
  expect_false("Residuals" %in% fit$effect_sizes$term)
  expect_false("(Intercept)" %in% fit$effect_sizes$term)
  expect_false(any(abs(fit$effect_sizes$partial_eta_sq - 0.5) < 1e-12))
})

test_that("type = 'III' fits under sum-to-zero contrasts and matches car", {
  d <- fx_twoway()
  fit <- anova_glm(d, "value", c("g1", "g2"), interaction = TRUE,
                   type = "III", plots = FALSE)

  expect_true(all(unlist(fit$model$contrasts) == "contr.sum"))

  ref_model <- withr::with_options(
    list(contrasts = c("contr.sum", "contr.poly")),
    stats::glm(value ~ g1 * g2, data = d))
  ref <- as.data.frame(car::Anova(ref_model, type = 3, test.statistic = "F"))
  ref <- ref[rownames(ref) != "(Intercept)", ]

  expect_equal(fit$anova$sum_sq, ref$`Sum Sq`)
  expect_equal(fit$anova$p_value, ref$`Pr(>F)`)
})

test_that("Type II and Type III differ on an unbalanced design with interaction", {
  d <- withr::with_seed(3, {
    dd <- data.frame(A = factor(c(rep("a1", 80), rep("a2", 20))),
                     B = factor(c(rep(c("b1", "b2"), each = 40),
                                  rep(c("b1", "b2"), each = 10))))
    dd$y <- 0.05 * (dd$B == "b2") + 3 * (dd$A == "a2") * (dd$B == "b2") +
      stats::rnorm(100)
    dd
  })
  t2 <- anova_glm(d, "y", c("A", "B"), interaction = TRUE, type = "II",
                  plots = FALSE)
  t3 <- anova_glm(d, "y", c("A", "B"), interaction = TRUE, type = "III",
                  plots = FALSE)

  b2 <- t2$anova$p_value[t2$anova$term == "B"]
  b3 <- t3$anova$p_value[t3$anova$term == "B"]
  expect_false(isTRUE(all.equal(b2, b3)))
  # The interaction row is invariant to the sums-of-squares type.
  expect_equal(t2$anova$sum_sq[t2$anova$term == "A:B"],
               t3$anova$sum_sq[t3$anova$term == "A:B"])
})

test_that("family accepts a string, a function and an object", {
  d <- fx_binary()
  as_string <- anova_glm(d, "y", "g", family = "binomial", plots = FALSE)
  as_fun    <- anova_glm(d, "y", "g", family = stats::binomial, plots = FALSE)
  as_object <- anova_glm(d, "y", "g", family = stats::binomial(), plots = FALSE)

  expect_equal(as_string$anova, as_fun$anova)
  expect_equal(as_string$anova, as_object$anova)
  expect_error(anova_glm(d, "y", "g", family = "not_a_family"),
               "not a known family")
})

test_that("post-hoc comparisons are produced for factorial designs", {
  d <- fx_twoway()
  fit <- anova_glm(d, "value", c("g1", "g2"), interaction = TRUE, plots = FALSE)
  expect_false(is.null(fit$posthoc))
  expect_equal(nrow(fit$posthoc), choose(4, 2))
  expect_true(all(c("contrast", "estimate", "p_value") %in% names(fit$posthoc)))
})

test_that("a two-level factor gets its comparison and an explanatory note", {
  d <- fx_oneway(means = c(A = 0, B = 1))
  fit <- anova_glm(d, "value", "group", plots = FALSE)
  expect_equal(nrow(fit$posthoc), 1L)
  expect_true(any(grepl("two cells", fit$notes)))
})

test_that("interaction = FALSE really is additive and TRUE really is factorial", {
  d <- fx_twoway()
  add <- anova_glm(d, "value", c("g1", "g2"), plots = FALSE)
  full <- anova_glm(d, "value", c("g1", "g2"), interaction = TRUE, plots = FALSE)

  expect_setequal(add$anova$term, c("g1", "g2", "Residuals"))
  expect_setequal(full$anova$term, c("g1", "g2", "g1:g2", "Residuals"))
  expect_lt(full$anova$p_value[full$anova$term == "g1:g2"], 0.05)
})

test_that("reserved arguments cannot be smuggled through ...", {
  d <- fx_oneway()
  # There is no ... to smuggle through: the modelling call is built entirely
  # from named arguments, so anything else is rejected by R itself.
  expect_error(anova_glm(d, "value", "group", formula = value ~ group),
               "unused argument")
  expect_error(anova_glm(d, "value", "group",
                         contrasts = list(group = "contr.treatment")),
               "unused argument")
})

test_that("nothing is written to stdout and no deprecation warnings fire", {
  d <- fx_oneway()
  expect_silent_stdout(anova_glm(d, "value", "group"))
  expect_no_warning(anova_glm(d, "value", "group"))
  fit <- anova_glm(d, "value", "group")
  expect_true(all(vapply(fit$plots, inherits, logical(1), "ggplot")))
})

test_that("deviance explained is 1 - deviance / null deviance", {
  d <- fx_counts()
  fit <- anova_glm(d, "count", "g1", family = "poisson", plots = FALSE)
  m <- stats::glm(count ~ g1, family = stats::poisson(), data = d)
  expect_equal(fit$effect_sizes$estimate[fit$effect_sizes$measure == "deviance_explained"],
               1 - m$deviance / m$null.deviance)
})

test_that("a Gamma fit gets car's F test by default", {
  d <- fx_oneway()
  d$pos <- abs(d$value) + 1
  fit <- anova_glm(d, "pos", "group", family = stats::Gamma(link = "log"),
                   plots = FALSE)
  m <- stats::glm(pos ~ group, data = d, family = stats::Gamma(link = "log"))
  ref <- car::Anova(m, test.statistic = "F")
  expect_identical(attr(fit$anova, "statistic"), "F")
  expect_equal(fit$anova$statistic[1], ref$`F value`[1])
  expect_equal(fit$anova$p_value[1], ref$`Pr(>F)`[1])
})

test_that("a log link gives ratios of marginal means in a 'ratio' column", {
  d <- fx_counts()
  fit <- anova_glm(d, "count", "g1", family = "poisson", plots = FALSE)
  expect_true("ratio" %in% names(fit$posthoc))
  mu <- tapply(d$count, d$g1, mean)
  expect_equal(fit$posthoc$ratio, unname(mu["A"] / mu["B"]), tolerance = 1e-6)
})

test_that("robust standard errors reach the marginal means", {
  skip_if_not_installed("sandwich")
  d <- fx_oneway(sds = c(1, 2, 4))
  fit <- anova_glm(d, "value", "group", vcov_type = "HC3", plots = FALSE)
  m <- stats::glm(value ~ group, data = d)
  ref <- as.data.frame(emmeans::emmeans(m, "group",
                                        vcov. = sandwich::vcovHC(m, type = "HC3")))
  expect_equal(fit$emmeans$se, ref$SE, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(
    fit$emmeans$se, anova_glm(d, "value", "group", plots = FALSE)$emmeans$se)))
})

test_that("Poisson overdispersion is noted above a Pearson dispersion of 1.2", {
  d <- fx_overdispersed()
  fit <- anova_glm(d, "y", "g", family = "poisson", plots = FALSE)
  m <- stats::glm(y ~ g, family = stats::poisson(), data = d)
  phi <- sum(stats::residuals(m, type = "pearson")^2) / stats::df.residual(m)
  expect_equal(fit$assumptions$dispersion, phi)
  expect_gt(phi, 1.5)
  expect_lt(phi, 15)
  expect_true(any(grepl(sprintf("Pearson dispersion is %.2f", phi), fit$notes,
                        fixed = TRUE)))
  expect_true(any(grepl("family = \"quasipoisson\"", fit$notes, fixed = TRUE)))
  ok <- anova_glm(fx_counts(), "count", c("g1", "g2"), interaction = TRUE,
                  family = "poisson", plots = FALSE)
  expect_lt(ok$assumptions$dispersion, 1.2)
  expect_false(any(grepl("^Pearson dispersion is", ok$notes)))
})
