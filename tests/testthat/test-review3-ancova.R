# anova_ancova(): regression tests for the third review round.
#
# Every expectation is checked against an independent computation -- lm(),
# predict(), hatvalues(), car, emmeans and sandwich on a model fitted here, or
# a hand computation -- never against the package's own helpers.

# R3-07: an uncentred covariate with a large offset ----------------------------

test_that("R3-07: an offset covariate is not called constant within groups", {
  d <- withr::with_seed(10, {
    d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 10)),
                    t = 1.7e9 + stats::runif(30, 0, 60))
    d$y <- as.integer(d$g) + 0.05 * (d$t - 1.7e9) + stats::rnorm(30)
    d
  })
  # The covariate varies by tens of seconds within every group ...
  expect_true(all(tapply(d$t, d$g, function(v) diff(range(v))) > 30))
  # ... but uncentred, lm() itself cannot separate it from the intercept
  expect_true(is.na(stats::coef(stats::lm(y ~ t + g, data = d))[["t"]]))

  err <- tryCatch(anova_ancova(d, "y", "g", "t", center_covariates = FALSE,
                               plots = FALSE), error = function(e) e)
  expect_s3_class(err, "error")
  expect_false(grepl("constant within every group", conditionMessage(err)))
  expect_match(conditionMessage(err), "vary too little relative to their distance from zero")
  expect_match(conditionMessage(err), "center_covariates = TRUE", fixed = TRUE)

  # The default (centred) fit matches car on the centred covariate
  f <- anova_ancova(d, "y", "g", "t", plots = FALSE)
  d$tc <- d$t - mean(d$t)
  ref <- as.data.frame(car::Anova(stats::lm(y ~ tc + g, data = d,
                                            contrasts = list(g = "contr.sum")),
                                  type = 3))
  expect_equal(f$anova$statistic[f$anova$term == "g"], ref["g", "F value"])
  expect_equal(f$anova$p_value[f$anova$term == "t"], ref["tc", "Pr(>F)"])

  # A covariate that really is constant within every group keeps its message,
  # uncentred and with an offset too
  d$lvl <- 1.7e9 + c(a = 1, b = 5, c = 9)[as.character(d$g)]
  expect_error(anova_ancova(d, "y", "g", "lvl", center_covariates = FALSE,
                            plots = FALSE),
               "constant within every group")
})

# R3-08: digits in the covariate notes ------------------------------------------

test_that("R3-08: the notes print covariate values with enough digits", {
  d <- withr::with_seed(14, {
    d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 10)))
    d$x <- 2020 + c(stats::runif(10, .10, .40), stats::runif(10, .5, .9),
                    stats::runif(10, .45, .95))
    d$y <- as.integer(d$g) + d$x - 2020 + stats::rnorm(30)
    d
  })
  f <- anova_ancova(d, "y", "g", "x", plots = FALSE)
  ext <- grep("lies outside the range observed in a ", f$notes, value = TRUE)
  expect_length(ext, 1L)
  nums <- as.numeric(regmatches(
    ext, regexec("x = ([0-9.]+) lies outside the range observed in a \\[([0-9.]+), ([0-9.]+)\\]",
                 ext))[[1]][-1])
  ra <- range(d$x[d$g == "a"])
  # Each printed value is the true one to well within the spread ...
  expect_equal(nums, c(mean(d$x), ra), tolerance = 1e-6)
  # ... so the printed note is self-consistent: the mean lies above the range
  expect_true(nums[1] > nums[3])
  cen <- grep("mean-centred before fitting", f$notes, value = TRUE)
  m <- as.numeric(sub(".*\\(x: ([0-9.]+)\\).*", "\\1", cen))
  expect_equal(m, mean(d$x), tolerance = 1e-6)
  expect_false(grepl("x: 2021)", cen, fixed = TRUE))

  unc <- anova_ancova(d, "y", "g", "x", center_covariates = FALSE, plots = FALSE)
  rng <- grep("Observed ranges", unc$notes, value = TRUE)
  got <- as.numeric(regmatches(rng, regexec("x \\[([0-9.]+), ([0-9.]+)\\]", rng))[[1]][-1])
  expect_equal(got, range(d$x), tolerance = 1e-6)
})

# R3-15: the Levene note when the robust covariance was refused -----------------

test_that("R3-15: the Levene note does not recommend a refused vcov_type", {
  skip_if_not_installed("sandwich")
  d <- withr::with_seed(4, {
    d <- data.frame(g = rep(c("a", "b", "c"), c(2, 20, 20)), x = stats::rnorm(42))
    d$y <- d$x + c(stats::rnorm(2), stats::rnorm(20, 0, 0.2), stats::rnorm(20, 0, 5))
    d
  })
  f <- anova_ancova(d, "y", "g", "x", force_interaction = TRUE,
                    vcov_type = "HC3", plots = FALSE)
  # Reference: the two rows of group a have leverage 1 in this model
  ref <- stats::lm(y ~ x * g, data = d)
  expect_equal(unname(stats::hatvalues(ref)[d$g == "a"]), c(1, 1))
  expect_identical(attr(f$anova, "statistic"), "F")
  lev <- grep("Levene's test rejects", f$notes, value = TRUE)
  expect_length(lev, 1L)
  expect_false(grepl("bases them on a heteroscedasticity-consistent", lev))
  expect_false(grepl("already allows for this", lev))
  expect_match(lev, "could not be used for this fit")

  # In the additive model $simple_slopes is NULL and the note does not name it
  h <- withr::with_seed(7, {
    h <- data.frame(g = rep(c("a", "b", "c"), each = 30))
    h$cv <- stats::rnorm(90, 50, 10)
    h$y <- 0.3 * h$cv + (h$g == "b") +
      stats::rnorm(90, 0, ifelse(h$g == "c", 4, 1))
    h
  })
  a0 <- anova_ancova(h, "y", "g", "cv", force_interaction = FALSE, plots = FALSE)
  expect_null(a0$simple_slopes)
  lev0 <- grep("Levene's test rejects", a0$notes, value = TRUE)
  expect_length(lev0, 1L)
  expect_lt(car::leveneTest(stats::residuals(stats::lm(y ~ cv + g, data = h)) ~
                              factor(h$g))[1, "Pr(>F)"], 0.05)
  expect_false(grepl("simple_slopes", lev0, fixed = TRUE))
  expect_match(lev0, "vcov_type = \"HC3\"", fixed = TRUE)
})

# R3-17: the covariate plot with a group of two -------------------------------

test_that("R3-17: drawing the covariate plot is silent with groups of two", {
  d <- data.frame(g = rep(c("a", "b", "c"), each = 2),
                  x = c(1, 2, 1.5, 3, 2, 2.5),
                  y = c(2, 3.1, 2.2, 4.8, 3.9, 4.1))
  f <- anova_ancova(d, "y", "g", "x")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(print(plot(f, "covariate")))
  # The bands come from the fitted (additive) model, which has residual df
  ref <- stats::lm(y ~ x + g, data = d)
  expect_equal(stats::df.residual(ref), 2L)
  rib <- ggplot2::layer_data(f$plots$covariate, 2)
  lev <- levels(factor(d$g))
  for (k in seq_along(lev)) {
    s <- rib[rib$group == k, ]
    p <- stats::predict(ref, data.frame(x = s$x, g = lev[k]),
                        interval = "confidence")
    expect_equal(s$ymin, unname(p[, "lwr"]))
    expect_equal(s$ymax, unname(p[, "upr"]))
  }

  # The per-group fallback draws a band only where a line has residual df
  p2 <- anovakit:::.plot_covariate(
    data.frame(y = c(d$y, 5, 6, 6.5), x = c(d$x, 3, 4, 5),
               g = c(d$g, "d", "d", "d")), "y", "x", "g")
  expect_silent(print(p2))
  # Layer 2: the band, for group d alone (three observations, x from 3 to 5),
  # equal to lm()'s interval for that group's own line
  bands <- ggplot2::layer_data(p2, 2)
  expect_equal(range(bands$x), c(3, 5))
  dd <- data.frame(y = c(5, 6, 6.5), x = c(3, 4, 5))
  pd <- stats::predict(stats::lm(y ~ x, data = dd), data.frame(x = bands$x),
                       interval = "confidence")
  expect_equal(bands$ymin, unname(pd[, "lwr"]))
  # Layer 3: lines without a band for the groups of two
  nob <- ggplot2::layer_data(p2, 3)
  expect_equal(length(unique(nob$group)), 3L)
  expect_true(is.null(nob$ymin) || all(is.na(nob$ymin)))
  expect_match(p2$labels$subtitle, "two observations have a line but no band")
})

# R3-34: Levene's test leaves out cells with one or two observations ------------

test_that("R3-34: Levene's test uses only cells with at least three observations", {
  sizes <- c(rep(1, 15), rep(2, 10), rep(6, 5))
  d <- withr::with_seed(2, {
    d <- data.frame(g = rep(sprintf("g%02d", seq_along(sizes)), sizes))
    d$x <- stats::rnorm(nrow(d))
    d$y <- stats::rnorm(nrow(d)) + d$x
    d
  })
  f <- anova_ancova(d, "y", "g", "x", plots = FALSE, posthoc = FALSE)
  # Reference: car on the fitted model's residuals, restricted to cells of 3+
  m <- if (is.null(f$simple_slopes)) stats::lm(y ~ x + g, data = d) else
    stats::lm(y ~ x * g, data = d)
  r <- stats::residuals(m)
  keep <- d$g %in% names(which(table(d$g) >= 3))
  ref <- car::leveneTest(r[keep] ~ factor(d$g[keep]))
  expect_equal(f$assumptions$levene$statistic, ref[["F value"]][1])
  expect_equal(f$assumptions$levene$p_value, ref[["Pr(>F)"]][1])
  expect_equal(f$assumptions$levene$df1, 4)
  expect_true(any(grepl("Levene's test uses only the 5 groups with at least three observations; 25 group(s)",
                        f$notes, fixed = TRUE)))

  # Fewer than two cells of three or more: skipped, with a reason
  d2 <- d[d$g %in% c(sprintf("g%02d", c(1:4, 16:19)), "g26"), ]
  f2 <- anova_ancova(d2, "y", "g", "x", plots = FALSE, posthoc = FALSE)
  expect_null(f2$assumptions$levene)
  expect_true(any(grepl("Levene's test was skipped: it needs at least two groups with three or more observations, and only one group has that many",
                        f2$notes, fixed = TRUE)))
})

# R3-47: a grouping column called `covariate` -----------------------------------

test_that("R3-47: a grouping column named covariate keeps its name, with a note", {
  d <- withr::with_seed(9, {
    n <- 90
    d <- data.frame(covariate = factor(rep(c("a", "b", "c"), each = n / 3)))
    d$x <- stats::rnorm(n, 10, 2)
    d$y <- 1 + d$x * c(a = 1, b = 2, c = 3)[as.character(d$covariate)] +
      stats::rnorm(n)
    d
  })
  fit <- anova_ancova(d, "y", "covariate", "x", plots = FALSE)
  ss <- fit$simple_slopes
  expect_identical(names(ss)[1:3], c("emm_covariate", "covariate", "slope"))
  expect_identical(ss$emm_covariate, rep("x", 3))
  expect_identical(as.character(ss$covariate), c("a", "b", "c"))
  ref <- as.data.frame(emmeans::emtrends(stats::lm(y ~ x * covariate, data = d),
                                         ~covariate, var = "x"))
  expect_equal(ss$slope, ref$x.trend)
  expect_true(any(grepl("covariate -> emm_covariate", fit$notes, fixed = TRUE)))
  # An ordinary grouping name keeps the plain label column and no note
  names(d)[1] <- "grp"
  plain <- anova_ancova(d, "y", "grp", "x", plots = FALSE)
  expect_identical(names(plain$simple_slopes)[1:2], c("covariate", "grp"))
  expect_false(any(grepl("covariate ->", plain$notes, fixed = TRUE)))
})

# R3-49: the post-selection caveat on the additive branch ----------------------

test_that("R3-49: the additive branch chosen by the test carries the selection caveat", {
  # Randomised groups (equal expected covariate means) with unequal slopes:
  # the additive model is chosen here although the slopes differ.
  d <- withr::with_seed(2, {
    g <- factor(rep(c("a", "b", "c"), each = 20))
    x <- stats::rnorm(60)
    data.frame(g = g, x = x,
               y = c(0.8, 1, 1.2)[as.integer(g)] * x + stats::rnorm(60))
  })
  ref_p <- stats::anova(stats::lm(y ~ x + g, data = d),
                        stats::lm(y ~ x * g, data = d))[2, "Pr(>F)"]
  expect_gt(ref_p, 0.05)
  f <- anova_ancova(d, "y", "g", "x", plots = FALSE)
  expect_null(f$simple_slopes)
  sel <- grep("Slopes are consistent with being equal", f$notes, value = TRUE)
  expect_length(sel, 1L)
  expect_match(sel, "do not allow for the model having been chosen", fixed = TRUE)
  expect_match(sel, "even when the groups were randomised", fixed = TRUE)
  expect_match(sel, "force_interaction = TRUE", fixed = TRUE)
  # With the model fixed in advance there was no selection to warn about
  ff <- anova_ancova(d, "y", "g", "x", force_interaction = FALSE, plots = FALSE)
  expect_false(any(grepl("do not allow for the model having been chosen",
                         ff$notes, fixed = TRUE)))
  # The help page no longer claims randomised groups are unaffected
  src <- testthat::test_path("..", "..", "R", "anova_ancova.R")
  skip_if_not(file.exists(src))
  txt <- paste(readLines(src), collapse = " ")
  expect_false(grepl("Randomised groups were #' not affected", txt, fixed = TRUE))
  expect_false(grepl("Randomised groups were not affected", txt, fixed = TRUE))
})

# R3-50: the covariate plot shows the fitted model ------------------------------

test_that("R3-50: the covariate plot draws the fitted model's lines", {
  d <- withr::with_seed(7, {
    n <- 150
    d <- data.frame(grp = factor(rep(c("control", "treated"), each = n / 2)),
                    baseline = stats::rnorm(n, 100, 5))
    d$score <- 2 * d$baseline + ifelse(d$grp == "treated", 6, 0) +
      stats::rnorm(n, 0, 3)
    d
  })
  fit <- anova_ancova(d, "score", "grp", "baseline")
  expect_null(fit$simple_slopes)                  # the additive model
  ln <- ggplot2::layer_data(fit$plots$covariate, 3)
  common <- stats::coef(stats::lm(score ~ baseline + grp, data = d))[["baseline"]]
  slope_of <- function(s) unname(stats::coef(stats::lm(y ~ x, data = s))[2])
  sl <- vapply(split(ln, ln$group), slope_of, numeric(1))
  expect_equal(unname(sl), rep(common, 2))
  # Each line spans its group's observed baseline, on the original scale
  expect_equal(range(ln$x[ln$group == 1]), range(d$baseline[d$grp == "control"]))
  # Its height is the additive model's prediction, raw covariate scale
  ref <- stats::lm(score ~ baseline + grp, data = d)
  s1 <- ln[ln$group == 2, ]
  expect_equal(s1$y, unname(stats::predict(ref, data.frame(baseline = s1$x,
                                                            grp = "treated"))))
  expect_match(fit$plots$covariate$labels$subtitle, "one common slope", fixed = TRUE)
  # Each group's own line is still there, dashed, as a check on parallelism
  own_ln <- ggplot2::layer_data(fit$plots$covariate, 4)
  expect_true(all(own_ln$linetype == "dashed"))
  own_add <- vapply(c("control", "treated"), function(l) unname(stats::coef(
    stats::lm(score ~ baseline, data = d[d$grp == l, ]))[2]), numeric(1))
  expect_equal(unname(vapply(split(own_ln, own_ln$group), slope_of, numeric(1))),
               unname(own_add))
  expect_false(isTRUE(all.equal(unname(own_add), rep(common, 2))))

  # With the interaction, each group's line is that group's own slope
  fi <- anova_ancova(d, "score", "grp", "baseline", force_interaction = TRUE)
  li <- ggplot2::layer_data(fi$plots$covariate, 3)
  own <- vapply(c("control", "treated"), function(l) unname(stats::coef(
    stats::lm(score ~ baseline, data = d[d$grp == l, ]))[2]), numeric(1))
  expect_equal(unname(vapply(split(li, li$group), slope_of, numeric(1))),
               unname(own))
  expect_length(fi$plots$covariate$layers, 3L)   # no dashed duplicate
})

# R3-58: vcov_type reaches $emmeans and $posthoc ------------------------------

test_that("R3-58: an HC covariance reaches $emmeans and $posthoc", {
  skip_if_not_installed("sandwich")
  d <- withr::with_seed(11, {
    d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 25)),
                    x = stats::rnorm(75, 10, 2))
    d$y <- 1 + 0.5 * d$x + stats::rnorm(75, 0, c(a = 0.5, b = 1, c = 4)[as.character(d$g)])
    d
  })
  f <- anova_ancova(d, "y", "g", "x", force_interaction = FALSE,
                    vcov_type = "HC3", plots = FALSE)
  d$xc <- d$x - mean(d$x)
  m <- stats::lm(y ~ xc + g, data = d, contrasts = list(g = "contr.sum"))
  V <- sandwich::vcovHC(m, type = "HC3")
  em <- emmeans::emmeans(m, ~g, vcov. = V)
  ref_emm <- as.data.frame(summary(em))
  expect_equal(f$emmeans$estimate, ref_emm$emmean)
  expect_equal(f$emmeans$se, ref_emm$SE)
  ref_ph <- as.data.frame(summary(emmeans::contrast(em, "pairwise",
                                                    adjust = "tukey")))
  expect_equal(f$posthoc$se, ref_ph$SE)
  expect_equal(f$posthoc$p_adjusted, ref_ph$p.value, tolerance = 1e-3)
  # ... and they differ from the model-based ones, so the check has teeth
  mb <- as.data.frame(summary(emmeans::emmeans(m, ~g)))
  expect_false(isTRUE(all.equal(f$emmeans$se, mb$SE)))
})
