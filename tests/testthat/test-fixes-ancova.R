# anova_ancova(): regression tests for the adversarial-review findings.
#
# Every expectation is checked against an independent computation -- lm(),
# predict(), car, emmeans on a model fitted here, sandwich, per-cell OLS
# slopes, or a hand computation from the definition -- never against the
# package's own helpers.

# Fixtures ---------------------------------------------------------------------

anc_two_factor <- function(seed = 5) {
  withr::with_seed(seed, {
    d <- expand.grid(rep = 1:40, A = c("a1", "a2"), B = c("b1", "b2"))
    d$rep <- NULL
    d$x <- stats::rnorm(nrow(d), 20, 4)
    # True cell means at x = 20: 12, 12, 12, 16 -- an A:B interaction only
    d$y <- 2 + 0.5 * d$x + 4 * (d$A == "a2" & d$B == "b2") +
      stats::rnorm(nrow(d))
    d
  })
}

anc_crossover_slopes <- function(seed = 6) {
  withr::with_seed(seed, {
    d <- expand.grid(rep = 1:40, A = c("a1", "a2"), B = c("b1", "b2"))
    d$rep <- NULL
    d$x <- stats::rnorm(nrow(d), 20, 4)
    # Slope 1 in a1:b1 and a2:b2, 0 elsewhere: no slope difference by A or by
    # B on average, only by cell
    d$y <- 2 + ifelse((d$A == "a1") == (d$B == "b1"), 1, 0) * (d$x - 20) +
      stats::rnorm(nrow(d))
    d
  })
}

anc_hetero <- function(seed = 3) {
  withr::with_seed(seed, {
    n <- 120
    g <- factor(rep(c("a", "b", "c"), each = n / 3))
    x <- stats::rnorm(n, 50, 10)
    sdv <- c(a = 1, b = 3, c = 6)[as.character(g)]
    y <- 1 + c(a = 0, b = 0.5, c = 1.5)[as.character(g)] * x + as.integer(g) +
      stats::rnorm(n, 0, sdv)
    data.frame(y = y, g = g, x = x)
  })
}

cell_slope <- function(d, rows, response = "y", covariate = "x") {
  unname(stats::coef(stats::lm(stats::reformulate(covariate, response),
                               data = d[rows, ]))[2])
}

# F010: the group-by-group interaction and per-cell slopes ---------------------

test_that("F010: two grouping factors are crossed by default, as in a factorial ANCOVA", {
  d <- anc_two_factor()
  fit <- anova_ancova(d, "y", c("A", "B"), "x", plots = FALSE)

  expect_true("A:B" %in% fit$anova$term)
  ref_model <- withr::with_options(
    list(contrasts = c("contr.sum", "contr.poly")),
    stats::lm(y ~ x + A * B, data = d))
  ref <- as.data.frame(car::Anova(ref_model, type = 3))
  expect_equal(fit$anova$statistic[fit$anova$term == "A:B"], ref["A:B", "F value"])
  expect_equal(fit$anova$sum_sq[fit$anova$term == "Residuals"],
               sum(stats::residuals(ref_model)^2))

  # Adjusted cell means are the crossed model's predictions at the mean of x
  nd <- fit$emmeans[c("A", "B")]
  nd$x <- mean(d$x)
  expect_equal(fit$emmeans$estimate, unname(stats::predict(ref_model, nd)))
})

test_that("F010: the slopes test and simple slopes are per cell of the crossed design", {
  d <- anc_crossover_slopes()
  fit <- anova_ancova(d, "y", c("A", "B"), "x", plots = FALSE)

  ref <- stats::anova(stats::lm(y ~ x + A * B, d), stats::lm(y ~ x * A * B, d))
  expect_equal(fit$slopes_test$df, ref[["Df"]][2])
  expect_equal(fit$slopes_test$statistic, ref[["F"]][2])
  expect_false(fit$slopes_test$homogeneous)

  # With a slope per cell, each simple slope is that cell's own OLS slope
  ss <- fit$simple_slopes
  for (i in seq_len(nrow(ss))) {
    rows <- d$A == ss$A[i] & d$B == ss$B[i]
    expect_equal(ss$slope[i], cell_slope(d, rows), info = i)
  }
  expect_equal(ss$slope, c(1, 0, 0, 1), tolerance = 0.1)
})

test_that("F010: interaction = FALSE enters the groups additively and tests slopes by factor", {
  d <- anc_crossover_slopes()
  fit <- anova_ancova(d, "y", c("A", "B"), "x", interaction = FALSE,
                      plots = FALSE)
  expect_false("A:B" %in% fit$anova$term)
  ref <- stats::anova(stats::lm(y ~ x + A + B, d),
                      stats::lm(y ~ x + A + B + x:A + x:B, d))
  expect_equal(fit$slopes_test$statistic, ref[["F"]][2])
  # Additive groups: marginal means per factor, averaged over the other
  expect_true("term" %in% names(fit$emmeans))
  expect_setequal(unique(fit$emmeans$term), c("A", "B"))
  # interaction = 1 is the same model
  one <- anova_ancova(d, "y", c("A", "B"), "x", interaction = 1, plots = FALSE)
  expect_equal(one$anova, fit$anova)

  expect_error(anova_ancova(d, "y", c("A", "B"), "x", interaction = 3),
               "exceeds the number of grouping variables")
  expect_error(anova_ancova(d, "y", c("A", "B"), "x", interaction = "yes"),
               "must be TRUE, FALSE, or a whole number")
})

test_that("F010: an empty cell of the crossed design is left out, not invented", {
  d <- anc_two_factor()
  d <- d[!(d$A == "a2" & d$B == "b2"), ]
  fit <- anova_ancova(d, "y", c("A", "B"), "x", plots = FALSE)
  expect_equal(nrow(fit$emmeans), 3L)
  expect_false(any(fit$emmeans$A == "a2" & fit$emmeans$B == "b2"))
  expect_true(any(grepl("a2 : b2", fit$notes, fixed = TRUE)))
})

# F061: a two-valued covariate is held at its mean ------------------------------

test_that("F061: a 0/1 covariate is held at its mean, not averaged over 0 and 1", {
  d <- withr::with_seed(11, {
    n <- 200
    d <- data.frame(g = factor(rep(c("a", "b", "c"), length.out = n)))
    d$sex <- stats::rbinom(n, 1, 0.2)
    d$y <- 5 + 10 * d$sex + as.integer(d$g) + stats::rnorm(n)
    d
  })
  fit <- anova_ancova(d, "y", "g", "sex", force_interaction = FALSE,
                      plots = FALSE)
  m <- stats::lm(y ~ sex + g, data = d)
  ref <- stats::predict(m, data.frame(g = levels(d$g), sex = mean(d$sex)))
  expect_equal(fit$emmeans$estimate, unname(ref))

  # With the interaction, $posthoc compares the groups at the same point as the
  # Type III group row, so for two groups t^2 equals that row's F
  e <- withr::with_seed(12, {
    e <- data.frame(g = factor(rep(c("ctl", "trt"), each = 150)))
    e$smoker <- stats::rbinom(300, 1, 0.15)
    e$y <- 10 + 1 * (e$g == "trt") + 8 * e$smoker * (e$g == "trt") +
      stats::rnorm(300)
    e
  })
  fi <- anova_ancova(e, "y", "g", "smoker", force_interaction = TRUE,
                     plots = FALSE)
  mi <- stats::lm(y ~ smoker * g, data = e)
  p <- stats::predict(mi, data.frame(g = c("ctl", "trt"), smoker = mean(e$smoker)))
  expect_equal(fi$posthoc$estimate, unname(p[1] - p[2]))
  expect_equal(fi$posthoc$statistic^2,
               fi$anova$statistic[fi$anova$term == "g"])
})

# F011: $emmeans are at the covariate mean whether or not it is centred --------

test_that("F011: the uncentred note does not claim $emmeans sit at covariate = 0", {
  d <- withr::with_seed(7, {
    n <- 150
    d <- data.frame(grp = factor(rep(c("control", "treated"), each = n / 2)),
                    baseline = stats::rnorm(n, 100, 5))
    d$score <- 2 * d$baseline + ifelse(d$grp == "treated", 6, 0) +
      stats::rnorm(n, 0, 3)
    d
  })
  c1 <- anova_ancova(d, "score", "grp", "baseline", plots = FALSE)
  u1 <- anova_ancova(d, "score", "grp", "baseline", center_covariates = FALSE,
                     plots = FALSE)
  m <- stats::lm(score ~ baseline + grp, data = d)
  at_mean <- unname(stats::predict(
    m, data.frame(grp = levels(d$grp), baseline = mean(d$baseline))))
  expect_equal(u1$emmeans$estimate, at_mean)
  expect_equal(c1$emmeans$estimate, at_mean)

  note <- grep("NOT centred", u1$notes, value = TRUE)
  expect_length(note, 1L)
  expect_false(grepl("emmeans and the intercept are read at covariate = 0",
                     note, fixed = TRUE))
  expect_match(note, "$emmeans and $posthoc are still evaluated at the covariate mean",
               fixed = TRUE)
  expect_match(grep("mean-centred", c1$notes, value = TRUE),
               "$emmeans are evaluated at the covariate mean either way",
               fixed = TRUE)
})

# F012 / F013: vcov_type reaches $simple_slopes and $anova ----------------------

test_that("F012: robust standard errors reach the simple slopes", {
  skip_if_not_installed("sandwich")
  d <- anc_hetero()
  fit <- anova_ancova(d, "y", "g", "x", vcov_type = "HC3", plots = FALSE)
  expect_false(is.null(fit$simple_slopes))
  # HC3 is invariant to the parameterisation, so the per-group slope model
  # y ~ 0 + g + g:x gives the slopes' robust SEs directly
  cell_model <- stats::lm(y ~ 0 + g + g:x, data = d)
  V <- sandwich::vcovHC(cell_model, type = "HC3")
  slopes <- grep(":x$", names(stats::coef(cell_model)), value = TRUE)
  expect_equal(fit$simple_slopes$slope, unname(stats::coef(cell_model)[slopes]))
  expect_equal(fit$simple_slopes$se, unname(sqrt(diag(V)[slopes])))
  # ... and they differ from the model-based ones
  mb <- anova_ancova(d, "y", "g", "x", plots = FALSE)
  expect_false(isTRUE(all.equal(fit$simple_slopes$se, mb$simple_slopes$se)))
})

test_that("F013: robust standard errors reach the F tests, and the Levene advice is true", {
  skip_if_not_installed("sandwich")
  d <- withr::with_seed(7, {
    d <- data.frame(g = rep(c("a", "b", "c"), each = 30))
    d$cv <- stats::rnorm(90, 50, 10)
    d$y <- 0.3 * d$cv + (d$g == "b") +
      stats::rnorm(90, 0, ifelse(d$g == "c", 4, 1))
    d
  })
  a0 <- anova_ancova(d, "y", "g", "cv", force_interaction = FALSE,
                     plots = FALSE)
  a3 <- anova_ancova(d, "y", "g", "cv", force_interaction = FALSE,
                     vcov_type = "HC3", plots = FALSE)

  ref_model <- withr::with_options(
    list(contrasts = c("contr.sum", "contr.poly")),
    stats::lm(y ~ cv + g, data = d))
  # car announces which covariance it used through message()
  ref <- as.data.frame(suppressMessages(
    car::Anova(ref_model, type = 3, white.adjust = "hc3")))
  expect_equal(a3$anova$statistic[a3$anova$term == "g"], ref["g", "F"])
  expect_equal(a3$anova$p_value[a3$anova$term == "g"], ref["g", "Pr(>F)"])
  expect_equal(a3$anova$p_value[a3$anova$term == "g"], 0.006756, tolerance = 1e-3)
  expect_true(any(grepl("Wald F tests", a3$notes, fixed = TRUE)))

  # Type II too
  a3ii <- anova_ancova(d, "y", "g", "cv", force_interaction = FALSE,
                       vcov_type = "HC3", type = "II", plots = FALSE)
  ref2 <- as.data.frame(suppressMessages(car::Anova(
    stats::lm(y ~ cv + g, data = d), type = 2, white.adjust = "hc3")))
  expect_equal(a3ii$anova$statistic[a3ii$anova$term == "g"], ref2["g", "F"])

  # Effect sizes stay model-based: identical to the model-based fit, and to
  # the definition on car's ordinary Type III table
  expect_equal(a3$effect_sizes, a0$effect_sizes)
  ord <- as.data.frame(car::Anova(ref_model, type = 3))
  expect_equal(a3$effect_sizes$partial_eta_sq[a3$effect_sizes$term == "g"],
               ord["g", "Sum Sq"] / (ord["g", "Sum Sq"] + ord["Residuals", "Sum Sq"]))

  # The Levene note recommends vcov_type for the F tests, and following that
  # advice now changes them
  expect_true(any(grepl("Levene's test rejects", a0$notes) &
                    grepl("vcov_type", a0$notes)))
  expect_false(isTRUE(all.equal(a0$anova$p_value, a3$anova$p_value)))
  expect_false(any(grepl("consider vcov_type", a3$notes, fixed = TRUE)))
})

# F014: a grouping column called `slope` ----------------------------------------

test_that("F014: a grouping column named slope or <covariate>.trend keeps its labels", {
  d <- withr::with_seed(9, {
    n <- 90
    d <- data.frame(slope = factor(rep(c("a", "b", "c"), each = n / 3)))
    d$x <- stats::rnorm(n, 10, 2)
    d$y <- 1 + d$x * c(a = 1, b = 2, c = 3)[as.character(d$slope)] +
      stats::rnorm(n)
    d
  })
  fit <- anova_ancova(d, "y", "slope", "x", force_interaction = TRUE,
                      plots = FALSE)
  ss <- fit$simple_slopes
  expect_true(is.factor(ss$slope))
  expect_identical(as.character(ss$slope), c("a", "b", "c"))
  expect_true(is.numeric(ss$emm_slope))
  ref <- vapply(c("a", "b", "c"), function(l) cell_slope(d, d$slope == l),
                numeric(1))
  expect_equal(ss$emm_slope, unname(ref))
  expect_true(any(grepl("slope -> emm_slope", fit$notes, fixed = TRUE)))

  names(d)[names(d) == "slope"] <- "x.trend"
  fit2 <- anova_ancova(d, "y", "x.trend", "x", force_interaction = TRUE,
                       plots = FALSE)
  expect_identical(as.character(fit2$simple_slopes$x.trend), c("a", "b", "c"))
  expect_equal(fit2$simple_slopes$slope, unname(ref))
})

# F062: covariates that cannot adjust anything ----------------------------------

test_that("F062: a constant or group-determined covariate is refused by name", {
  d <- withr::with_seed(8, {
    n <- 60
    data.frame(g = factor(rep(c("a", "b", "c"), each = n / 3)),
               y = stats::rnorm(n), z = stats::rnorm(n))
  })
  d$k <- 5
  expect_error(anova_ancova(d, "y", "g", "k", plots = FALSE),
               "Covariate `k` has the same value")
  # A group-level attribute: constant within every group (sd within = 0)
  d$attr <- as.numeric(d$g) * 3
  expect_equal(as.vector(tapply(d$attr, d$g, stats::sd)), c(0, 0, 0))
  expect_error(anova_ancova(d, "y", "g", "attr", plots = FALSE),
               "`attr` are determined by the grouping variable")
  # An exact linear combination of another covariate and the groups
  d$w <- 2 * d$z + d$attr
  expect_error(anova_ancova(d, "y", "g", c("z", "w"), plots = FALSE),
               "exact linear combination")
  # A covariate constant within some groups only is fine for the additive fit
  d$partial <- ifelse(d$g == "a", d$z, 1)
  expect_s3_class(anova_ancova(d, "y", "g", "partial", plots = FALSE),
                  "anovakit_fit")
})

# F129: small-unit responses keep their table -----------------------------------

test_that("F129: a response in small units keeps its ANOVA table and effect sizes", {
  d <- withr::with_seed(24, {
    n <- 60
    d <- data.frame(g = rep(c("a", "b", "c"), each = 20), x = stats::rnorm(n))
    d$y <- (d$x + rep(c(0, 0.5, 1), each = 20) + stats::rnorm(n)) * 1e-5
    d
  })
  f <- anova_ancova(d, "y", "g", "x", force_interaction = FALSE, plots = FALSE)
  ref <- stats::anova(stats::lm(y ~ x + g, data = d))
  expect_false(is.null(f$anova))
  expect_equal(f$anova$statistic[f$anova$term == "g"], ref["g", "F value"])
  expect_equal(f$anova$sum_sq[f$anova$term == "g"], ref["g", "Sum Sq"])
  expect_equal(nrow(f$effect_sizes), 2L)
})

# F133: Levene with two observations per group ----------------------------------

test_that("F133: Levene's test is skipped, silently, when it is degenerate", {
  d <- withr::with_seed(40, {
    d <- data.frame(g = rep(c("a", "b", "c", "d"), each = 2), x = stats::rnorm(8))
    d$y <- d$x + stats::rnorm(8)
    d
  })
  expect_silent(f <- anova_ancova(d, "y", "g", "x", plots = FALSE))
  # By definition, the absolute deviations from each pair's median are equal
  r <- stats::residuals(stats::lm(y ~ x + g, data = d))
  z <- abs(r - stats::ave(r, d$g, FUN = stats::median))
  expect_equal(as.vector(tapply(z, d$g, stats::var)), rep(0, 4))
  expect_null(f$assumptions$levene)
  expect_true(any(grepl("Levene's test was skipped", f$notes, fixed = TRUE)))
  expect_false(any(grepl("Levene's test rejects", f$notes, fixed = TRUE)))
})

test_that("F112: Levene's statistic is car's on the model's residuals and cells", {
  d <- anc_hetero()
  d$y[c(5, 50)] <- NA
  fit <- anova_ancova(d, "y", "g", "x", force_interaction = TRUE, plots = FALSE)
  keep <- !is.na(d$y)
  ref_model <- stats::lm(y ~ x * g, data = d[keep, ])
  ref <- car::leveneTest(stats::residuals(ref_model) ~ d$g[keep])
  expect_equal(fit$assumptions$levene$statistic, ref[["F value"]][1])
  expect_equal(fit$assumptions$levene$p_value, ref[["Pr(>F)"]][1])
})

test_that("F112: homogeneity_alpha is the threshold the model choice uses", {
  d <- withr::with_seed(21, {
    n <- 60
    g <- factor(rep(c("a", "b"), each = n / 2))
    x <- stats::rnorm(n, 10, 2)
    data.frame(g = g, x = x,
               y = x * ifelse(g == "a", 1, 1.25) + stats::rnorm(n))
  })
  p <- stats::anova(stats::lm(y ~ x + g, d), stats::lm(y ~ x * g, d))[["Pr(>F)"]][2]
  expect_true(p > 0.001 && p < 0.9)
  above <- anova_ancova(d, "y", "g", "x", homogeneity_alpha = min(0.99, p * 1.5),
                        plots = FALSE)
  below <- anova_ancova(d, "y", "g", "x", homogeneity_alpha = p / 1.5,
                        plots = FALSE)
  expect_equal(above$slopes_test$p_value, p)
  expect_false(above$slopes_test$homogeneous)
  expect_false(is.null(above$simple_slopes))
  expect_true(below$slopes_test$homogeneous)
  expect_null(below$simple_slopes)
})

# F148: the pretest's cost is stated on the branch where it arises ---------------

test_that("F148: an additive fit with imbalanced covariate means carries a bias warning", {
  # Covariate means one pooled SD apart; slopes differ, but n is small enough
  # that the slopes test misses it
  mk <- function(shift, seed) withr::with_seed(seed, {
    n <- 12
    g <- factor(rep(c("a", "b", "c"), each = n))
    x <- stats::rnorm(3 * n, rep(c(-shift, 0, shift), each = n)) + 100
    data.frame(g = g, x = x,
               y = c(0.8, 1, 1.2)[as.integer(g)] * (x - 100) + stats::rnorm(3 * n))
  })
  imb <- mk(1.5, 31)
  fit <- anova_ancova(imb, "y", "g", "x", plots = FALSE)
  expect_true(fit$slopes_test$homogeneous)

  # Hand computation of the imbalance: covariate-mean range over the pooled
  # within-group SD
  m <- tapply(imb$x, imb$g, mean)
  sp <- sqrt(sum(tapply(imb$x, imb$g, function(v) sum((v - mean(v))^2))) /
               (nrow(imb) - 3))
  gap <- (max(m) - min(m)) / sp
  expect_gt(gap, 0.5)
  note <- grep("covariate means differ between groups", fit$notes, value = TRUE)
  expect_length(note, 1L)
  expect_match(note, sprintf("x %.2g pooled", gap), fixed = TRUE)
  expect_match(note, "force_interaction = TRUE", fixed = TRUE)
  expect_false(any(grepl("mildly optimistic", fit$notes, fixed = TRUE)))

  # Balanced covariate: no such note
  bal <- mk(0, 32)
  fb <- anova_ancova(bal, "y", "g", "x", force_interaction = FALSE,
                     plots = FALSE)
  mb <- tapply(bal$x, bal$g, mean)
  expect_lt(diff(range(mb)) / stats::sd(bal$x), 0.5)
  expect_false(any(grepl("covariate means differ", fb$notes, fixed = TRUE)))

  # When the pretest retains the interaction, the selection caveat is there,
  # without "mildly"
  het <- anova_ancova(anc_hetero(), "y", "g", "x", plots = FALSE)
  expect_false(het$slopes_test$homogeneous)
  expect_true(any(grepl("chosen by a test on these same data", het$notes,
                        fixed = TRUE)))
  expect_false(any(grepl("mildly", het$notes, fixed = TRUE)))
})

# F160: covariate names that collide with coefficient names ---------------------

test_that("F160: a covariate named like a contrast coefficient does not corrupt the results", {
  d <- withr::with_seed(1, {
    n <- 90
    d <- data.frame(dose = rep(c("high", "low", "mid"), each = 30),
                    dose1 = stats::rnorm(n, 50, 10))
    d$y <- 2 + c(high = 3, low = 0, mid = 1)[d$dose] + 0.05 * d$dose1 +
      stats::rnorm(n)
    d
  })
  ok <- d
  names(ok)[names(ok) == "dose1"] <- "baseline"
  ref_model <- stats::lm(y ~ baseline + dose, data = ok)
  ref <- unname(stats::predict(ref_model, data.frame(
    dose = c("high", "low", "mid"), baseline = mean(ok$baseline))))

  bad <- anova_ancova(d, "y", "dose", "dose1", force_interaction = FALSE,
                      plots = FALSE)
  expect_equal(bad$emmeans$estimate, ref)
  expect_equal(bad$posthoc$estimate,
               c(ref[1] - ref[2], ref[1] - ref[3], ref[2] - ref[3]))
  expect_equal(anyDuplicated(names(stats::coef(bad$model))), 0L)
  expect_true(any(grepl("put in brackets in $model (dose[1], dose[2]", bad$notes,
                        fixed = TRUE)))

  # Type II names the treatment coefficients dose<level>
  d2 <- d
  names(d2)[names(d2) == "dose1"] <- "doselow"
  bad2 <- anova_ancova(d2, "y", "dose", "doselow", type = "II",
                       force_interaction = FALSE, plots = FALSE)
  expect_equal(bad2$emmeans$estimate, ref)

  # Simple slopes: each group's own OLS slope
  s_bad <- anova_ancova(d, "y", "dose", "dose1", force_interaction = TRUE,
                        plots = FALSE)
  own <- vapply(c("high", "low", "mid"), function(l)
    cell_slope(d, d$dose == l, covariate = "dose1"), numeric(1))
  expect_equal(s_bad$simple_slopes$slope, unname(own))
})

# F015 / F050 / F097: the covariate plot and $data_used -------------------------

test_that("F015/F050: the covariate plot is on the original scale with bands at conf_level", {
  d <- withr::with_seed(104, {
    g <- factor(rep(c("A", "B"), each = 75))
    x <- stats::rnorm(150, 100, 5)
    data.frame(dv = 2 * x + stats::rnorm(150, 0, 3), iv = g, cov = x)
  })
  f <- anova_ancova(d, "dv", "iv", "cov", conf_level = 0.5)
  pts <- ggplot2::layer_data(f$plots$covariate, 1)
  expect_equal(range(pts$x), range(d$cov))

  s <- ggplot2::layer_data(f$plots$covariate, 2)
  sA <- s[s$group == 1, ]
  mA <- stats::lm(dv ~ cov, data = d[d$iv == "A", ])
  band <- stats::predict(mA, data.frame(cov = sA$x), interval = "confidence",
                         level = 0.5)
  expect_equal(sA$ymin, unname(band[, "lwr"]))
  expect_equal(sA$ymax, unname(band[, "upr"]))
  expect_match(f$plots$covariate$labels$subtitle, "50%", fixed = TRUE)
})

test_that("F097/F094: $data_used and the plot keep only the analysed columns", {
  d <- withr::with_seed(1, {
    n <- 60
    d <- data.frame(rn = seq_len(n), g = rep(c("a", "b", "c"), each = n / 3),
                    x = stats::rnorm(n))
    d$y <- d$x + stats::rnorm(n)
    for (j in 1:3) d[[paste0("junk", j)]] <- stats::rnorm(n)
    d
  })
  names(d)[1] <- ""                     # an unnamed column
  fit <- anova_ancova(d, "y", "g", "x")
  expect_setequal(names(fit$data_used), c("y", "g", "x", ".cell"))
  expect_setequal(names(fit$plots$covariate$data), c("y", "x", ".cell"))
  # The centred covariate in $data_used is the raw column minus its mean
  expect_equal(fit$data_used$x, d$x - mean(d$x))
  expect_equal(unname(fit$covariate_means), mean(d$x))
  # $covariate_means is returned uncentred too
  unc <- anova_ancova(d, "y", "g", "x", center_covariates = FALSE, plots = FALSE)
  expect_equal(unname(unc$covariate_means), mean(d$x))
  expect_equal(unc$data_used$x, d$x)
})

# F064: homogeneity_alpha validation ---------------------------------------------

test_that("F064: a missing homogeneity_alpha gets the validation message", {
  d <- withr::with_seed(1, data.frame(g = factor(rep(c("a", "b"), each = 30)),
                                      x = stats::rnorm(60), y = stats::rnorm(60)))
  for (bad in list(NA_real_, NaN, NA, "0.05", c(0.05, 0.1))) {
    expect_error(anova_ancova(d, "y", "g", "x", homogeneity_alpha = bad),
                 "`homogeneity_alpha` must be a single number strictly between 0 and 1",
                 fixed = TRUE)
  }
})

# F132: an exact fit ------------------------------------------------------------

test_that("F132: the slopes test on an exact fit is not computable, whatever the row order", {
  d <- withr::with_seed(9, {
    d <- data.frame(g = rep(c("a", "b", "c"), each = 10), x = stats::rnorm(30))
    d$y <- 2 + 3 * d$x + c(a = 0, b = 1, c = 2)[d$g]
    d
  })
  # The interaction coefficients are zero to rounding error
  ic <- stats::coef(stats::lm(y ~ x * g, data = d))
  expect_true(all(abs(ic[grepl(":", names(ic))]) < 1e-10))
  withr::local_seed(10)
  for (i in 1:4) {
    dd <- d[sample(nrow(d)), ]
    expect_silent(f <- anova_ancova(dd, "y", "g", "x", plots = FALSE))
    expect_true(is.na(f$slopes_test$p_value), info = i)
    expect_true(is.na(f$slopes_test$homogeneous), info = i)
    expect_null(f$simple_slopes)
    expect_true(any(grepl("fits the data exactly", f$notes, fixed = TRUE)),
                info = i)
    expect_false(any(grepl("Slopes differ", f$notes, fixed = TRUE)), info = i)
  }
})

# F159: extrapolated adjusted means --------------------------------------------

test_that("F159: a group observed away from the covariate mean is flagged", {
  d <- withr::with_seed(15, {
    g <- factor(rep(c("junior", "senior"), c(40, 60)))
    yrs <- ifelse(g == "junior", stats::runif(100, 1, 11), stats::runif(100, 15, 50))
    data.frame(g = g, yrs = yrs, pay = 50 + 1.2 * yrs + 5 * (g == "senior") +
                 stats::rnorm(100, 0, 4))
  })
  fit <- anova_ancova(d, "pay", "g", "yrs", plots = FALSE)
  m <- mean(d$yrs)
  rng <- range(d$yrs[d$g == "junior"])
  expect_true(m > rng[2])                      # hand check: outside junior
  expect_true(m > min(d$yrs[d$g == "senior"])) # ... but inside senior
  note <- grep("lies outside the range observed", fit$notes, value = TRUE)
  expect_length(note, 1L)
  expect_match(note, sprintf("junior [%.4g, %.4g]", rng[1], rng[2]), fixed = TRUE)
  expect_false(grepl("senior [", note, fixed = TRUE))

  # A randomised covariate raises no such note
  ok <- anova_ancova(anc_hetero(), "y", "g", "x", plots = FALSE)
  expect_false(any(grepl("lies outside the range observed", ok$notes)))
})

# Nothing prints on the new paths ----------------------------------------------

test_that("the new paths write nothing to the console", {
  skip_if_not_installed("sandwich")
  expect_silent(anova_ancova(anc_two_factor(), "y", c("A", "B"), "x"))
  expect_silent(anova_ancova(anc_hetero(), "y", "g", "x", vcov_type = "HC3"))
  d <- anc_hetero()
  names(d)[names(d) == "x"] <- "g1"
  expect_silent(anova_ancova(d, "y", "g", "g1", force_interaction = TRUE))
})
