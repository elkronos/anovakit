# Regression tests. Each block below pins down a defect that was found and
# fixed before release, so that it cannot come back unnoticed. The comments
# say what the wrong behaviour was, not just what the right one is.

test_that("a homogeneous cell of an ADDITIVE model is not called separation", {
  # The old cell-based rule flagged ~48% of clean additive two-factor fits.
  d <- withr::with_seed(73, {
    dd <- expand.grid(A = factor(c("a1", "a2", "a3")),
                      B = factor(c("b1", "b2", "b3")), rep = 1:25)
    lp <- -0.3 + 1.2 * (dd$A == "a2") + 0.5 * (dd$A == "a3") +
      0.9 * (dd$B == "b2") - 0.4 * (dd$B == "b3")
    dd$y <- stats::rbinom(nrow(dd), 1, stats::plogis(lp))
    dd$y[dd$A == "a1" & dd$B == "b3"] <- 0     # one homogeneous crossed cell
    dd
  })
  fit <- anova_bin(d, "y", c("A", "B"), interaction = FALSE, plots = FALSE)

  expect_length(grep("separation", fit$notes), 0L)
  # and the fit really is well identified, which is why flagging it was wrong
  expect_lt(max(abs(stats::coef(fit$model))), 5)
  expect_true(all(is.finite(fit$effect_sizes$conf_high)))
})

test_that("genuine separation is still detected, and names the coefficient", {
  fit <- anova_bin(fx_separated(), "y", "g", plots = FALSE)
  sep <- grep("separation", fit$notes, value = TRUE)
  expect_length(sep, 1L)
  expect_match(sep, "coefficient")
  expect_gt(fit$effect_sizes$odds_ratio[fit$effect_sizes$term == "gc"], 1e5)
})

test_that("coefficient p-values use t for families that estimate dispersion", {
  d <- fx_oneway()
  fit <- anova_glm(d, "value", "group", plots = FALSE)
  ref <- summary(fit$model)$coefficients

  expect_identical(unique(fit$assumptions$coefficients$distribution), "t")
  expect_equal(fit$assumptions$coefficients$p_value, unname(ref[, 4]))

  # Gamma: the old code used the normal and was out by orders of magnitude
  dg <- d; dg$pos <- abs(dg$value) + 1
  fg <- anova_glm(dg, "pos", "group", family = stats::Gamma(link = "log"),
                  plots = FALSE)
  expect_equal(fg$assumptions$coefficients$p_value,
               unname(summary(fg$model)$coefficients[, 4]))
})

test_that("coefficient p-values use z for families with fixed dispersion", {
  fit <- anova_glm(fx_binary(), "y", "g", family = "binomial", plots = FALSE)
  expect_identical(unique(fit$assumptions$coefficients$distribution), "z")
  expect_equal(fit$assumptions$coefficients$p_value,
               unname(summary(fit$model)$coefficients[, 4]))
})

test_that("partial omega squared matches effectsize and never exceeds 1", {
  skip_if_not_installed("effectsize")
  d <- fx_twoway()
  fit <- anova_glm(d, "value", c("g1", "g2"), interaction = TRUE, type = "II",
                   plots = FALSE)
  # effectsize needs an lm-based table; the model is the same either way.
  ref <- effectsize::omega_squared(
    car::Anova(stats::lm(value ~ g1 * g2, data = d), type = 2),
    partial = TRUE, ci = NULL)
  expect_equal(fit$effect_sizes$partial_omega_sq, ref$Omega2_partial,
               tolerance = 1e-8)

  # a near-aliased design used to give "191% of the variance"
  near <- withr::with_seed(93, {
    A <- factor(c(rep("a1", 30), rep("a2", 30)))
    B <- factor(c(rep("b1", 29), "b2", "b1", rep("b2", 29)))
    data.frame(A = A, B = B,
               y = 10 * (A == "a2") - 10 * (B == "b2") +
                 stats::rnorm(60, sd = 0.2))
  })
  f3 <- anova_glm(near, "y", c("A", "B"), type = "III", plots = FALSE)
  expect_true(all(f3$effect_sizes$partial_omega_sq >= 0 &
                    f3$effect_sizes$partial_omega_sq <= 1))
})

test_that("a grouping column named like a result column is not overwritten", {
  for (nm in c("se", "df", "conf_low", "estimate")) {
    d <- data.frame(g = rep(c("a", "b", "c"), each = 20),
                    y = withr::with_seed(1, stats::rnorm(60)))
    names(d)[1] <- nm
    fit <- anova_glm(d, "y", nm, plots = FALSE)

    # the grouping column still holds the group labels
    expect_setequal(as.character(fit$emmeans[[nm]]), c("a", "b", "c"))
    # and nothing is duplicated
    expect_false(anyDuplicated(names(fit$emmeans)) > 0L)
    # and the displacement is reported, not silent
    expect_true(any(grepl("shares its name", fit$notes)), info = nm)
  }
})

test_that("Box's M and the discriminant analysis are computed in the adjusted space", {
  # Covariate spread differs by group; the residual covariances are identical.
  d <- withr::with_seed(5, {
    n <- 180
    g <- factor(rep(c("a", "b", "c"), each = n / 3))
    age <- stats::rnorm(n, 40, ifelse(g == "c", 30, 3))
    data.frame(g = g, age = age,
               score1 = 5 + 0.5 * age + stats::rnorm(n),
               score2 = 3 + 0.3 * age + stats::rnorm(n))
  })
  fit <- anova_manova(d, c("score1", "score2"), "g", covariates = "age",
                      plots = FALSE)
  expect_gt(fit$assumptions$box_m$p_value, 0.05)
  expect_true(any(grepl("covariates partialled out", fit$notes)))
})

test_that("discriminant scores discriminate as well as their eigenvalue claims", {
  d <- withr::with_seed(9, {
    n <- 180
    g <- factor(rep(c("a", "b", "c"), each = n / 3))
    age <- stats::rnorm(n, 40, 15)
    data.frame(g = g, age = age,
               score1 = 2 * age + 3 * (g == "b") + stats::rnorm(n),
               score2 = 1 * age - 3 * (g == "b") + stats::rnorm(n))
  })
  fit <- anova_manova(d, c("score1", "score2"), "g", covariates = "age",
                      plots = TRUE)
  scores <- fit$plots$canonical$data$Can1
  gg <- fit$data_used$g
  between_within <- sum(table(gg) * (tapply(scores, gg, mean) - mean(scores))^2) /
    sum((scores - stats::ave(scores, gg))^2)

  # The plotted axis must reproduce the eigenvalue it is labelled with. Before
  # the fix it was 0.075 against a reported 4.75, a factor of 63.
  expect_equal(between_within, fit$canonical$eigenvalue[1], tolerance = 0.01)
})

test_that("the canonical analysis names the term it describes", {
  d <- withr::with_seed(17, {
    dd <- data.frame(A = factor(rep(c("a1", "a2"), each = 90)),
                     B = factor(rep(rep(c("b1", "b2", "b3"), each = 30), 2)))
    dd$y1 <- 2 * (dd$A == "a2") + stats::rnorm(180)
    dd$y2 <- 3 * (dd$B == "b3") + stats::rnorm(180)
    dd
  })
  fit <- anova_manova(d, c("y1", "y2"), c("A", "B"), plots = FALSE)
  expect_false(is.na(fit$canonical_term))
  expect_true(any(grepl("canonical discriminant analysis describes", fit$notes)))
})

test_that("anova_manova survives non-syntactic grouping column names", {
  d <- fx_multivariate()
  names(d)[names(d) == "g"] <- "my group"
  fit <- anova_manova(d, c("score1", "score2"), "my group",
                      assumptions = FALSE, plots = FALSE)
  expect_equal(nrow(fit$multivariate), 1L)
  expect_equal(nrow(fit$emmeans), 6L)

  for (nm in c("if", "TRUE", "2way", "a-b")) {
    d2 <- fx_multivariate()
    names(d2)[names(d2) == "g"] <- nm
    expect_s3_class(anova_manova(d2, c("score1", "score2"), nm,
                                 assumptions = FALSE, plots = FALSE),
                    "anovakit_fit")
  }
})

test_that("a column name containing a backtick is refused, not silently rewritten", {
  d <- fx_oneway()
  names(d) <- c("group", "y`1")
  expect_error(anova_welch(d, "y`1", "group"), "backtick")
})

test_that("duplicated column names are refused", {
  d <- data.frame(g = rep(c("a", "b"), each = 20),
                  y = withr::with_seed(2, stats::rnorm(40)),
                  z = withr::with_seed(3, stats::rnorm(40)))
  names(d)[3] <- "y"
  expect_error(anova_welch(d, "y", "g"), "duplicated column name")
})

test_that("Inf in the response is treated as missing, not passed to a test", {
  d <- fx_oneway()
  d$value[3] <- Inf
  fit <- anova_welch(d, "value", "group", plots = FALSE)
  expect_identical(fit$n_removed, 1L)
  expect_true(any(grepl("infinite", fit$notes)))
  expect_true(all(is.finite(fit$data_used$value)))

  expect_s3_class(anova_kw(d, "value", "group", diagnostics = TRUE,
                           plots = FALSE), "anovakit_fit")
})

test_that("a zero-variance group is refused rather than returning NaN", {
  d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 5)),
                  y = c(withr::with_seed(7, stats::rnorm(5)), rep(2, 5),
                        withr::with_seed(8, stats::rnorm(5, 3))))
  expect_error(anova_welch(d, "y", "g"), "positive variance in every group")
  # the rank-based route still works on the same data
  expect_s3_class(anova_kw(d, "y", "g", plots = FALSE), "anovakit_fit")
})

test_that("a factor response with a binomial family does not crash the plots", {
  d <- withr::with_seed(4, {
    dd <- data.frame(A = factor(sample(c("a1", "a2"), 200, replace = TRUE)))
    dd$fy <- factor(sample(c("yes", "no"), 200, replace = TRUE))
    dd
  })
  fit <- anova_glm(d, "fy", "A", family = "binomial")
  expect_s3_class(fit, "anovakit_fit")
  expect_false("box" %in% names(fit$plots))

  # and a factor response with a Gaussian family is refused with a useful message
  expect_error(anova_glm(d, "fy", "A"), "must be numeric for a gaussian family")
})

test_that("weights are named as a column and stay aligned when rows are dropped", {
  d <- fx_binary()
  d$w <- withr::with_seed(6, sample(1:3, nrow(d), replace = TRUE))
  d$y[1:5] <- NA

  fit <- anova_bin(d, "y", "g", weights = "w", plots = FALSE)
  expect_identical(fit$n_removed, 5L)
  expect_equal(unname(stats::coef(fit$model)),
               unname(stats::coef(stats::glm(
                 y ~ g, data = d[!is.na(d$y), ], family = stats::binomial(),
                 weights = d$w[!is.na(d$y)]))))

  expect_error(anova_bin(d, "y", "g", weights = "not_a_column"),
               "not found in `data`")
  # A predictor cannot also be the weights: caught by name, before the frame
  # is prepared and the column is coerced to a factor
  expect_error(anova_glm(fx_oneway(), "value", "group", weights = "group"),
               "both a grouping variable and the weights")
  expect_error(anova_glm(fx_oneway(), "value", "group", weights = "value"),
               "both the response and the weights")
  d2 <- fx_oneway(); d2$w <- as.character(seq_len(nrow(d2)))
  expect_error(anova_glm(d2, "value", "group", weights = "w"),
               "must be numeric")
})

test_that("an ordered grouping factor accepts a reference level", {
  d <- fx_binary()
  d$g <- factor(d$g, ordered = TRUE)
  fit <- anova_bin(d, "y", "g", reference = list(g = "b"), plots = FALSE)
  expect_s3_class(fit, "anovakit_fit")
  expect_identical(levels(fit$data_used$g)[1], "b")
})

test_that("the response may not also be a covariate or a group", {
  d <- fx_ancova()
  expect_error(anova_ancova(d, "dv", "iv", "dv"), "cannot also be a covariate")
  expect_error(anova_manova(d, c("dv", "cov"), "iv", covariates = "cov"),
               "cannot be both a response and a predictor")
  expect_error(anova_manova(d, "dv", "iv", covariates = "iv"),
               "cannot be both a group and a covariate")
})

test_that("the uncentred-covariate note only claims zero is outside when it is", {
  inside <- withr::with_seed(13, {
    dd <- data.frame(g = factor(rep(c("a", "b"), each = 60)),
                     cv = stats::rnorm(120, 1, 50))
    dd$y <- 2 * dd$cv + ifelse(dd$g == "b", 5, 0) + stats::rnorm(120)
    dd
  })
  n_in <- anova_ancova(inside, "y", "g", "cv", center_covariates = FALSE,
                       plots = FALSE)$notes
  expect_false(any(grepl("Zero lies outside", n_in)))

  n_out <- anova_ancova(fx_ancova(), "dv", "iv", "cov",
                        center_covariates = FALSE, plots = FALSE)$notes
  expect_true(any(grepl("Zero lies outside", n_out)))
})

test_that("simple slopes and covariate plots cover every covariate", {
  d <- fx_ancova(slope_shift = 0.6)
  d$cov2 <- withr::with_seed(14, stats::rnorm(nrow(d), 10, 2))
  fit <- anova_ancova(d, "dv", "iv", c("cov", "cov2"), force_interaction = TRUE)

  expect_setequal(unique(fit$simple_slopes$covariate), c("cov", "cov2"))
  expect_true(all(c("covariate_cov", "covariate_cov2") %in% names(fit$plots)))
})

test_that("posthoc = FALSE is available everywhere and is reported", {
  callers <- list(
    welch  = function() anova_welch(fx_oneway(), "value", "group",
                                    posthoc = FALSE, plots = FALSE),
    kw     = function() anova_kw(fx_oneway(), "value", "group",
                                 posthoc = FALSE, plots = FALSE),
    glm    = function() anova_glm(fx_oneway(), "value", "group",
                                  posthoc = FALSE, plots = FALSE),
    bin    = function() anova_bin(fx_binary(), "y", "g",
                                  posthoc = FALSE, plots = FALSE),
    count  = function() anova_count(fx_counts(), "count", "g1",
                                    posthoc = FALSE, plots = FALSE),
    ancova = function() anova_ancova(fx_ancova(), "dv", "iv", "cov",
                                     posthoc = FALSE, plots = FALSE),
    manova = function() anova_manova(fx_multivariate(), c("score1", "score2"),
                                     "g", posthoc = FALSE, assumptions = FALSE,
                                     plots = FALSE)
  )
  for (nm in names(callers)) {
    fit <- callers[[nm]]()
    expect_null(fit$posthoc, info = nm)
    expect_false(is.null(fit$anova), info = nm)
  }
})

test_that("the count dispersion of the fitted model is reported separately", {
  skip_if_not_installed("MASS")
  d <- fx_overdispersed()
  fit <- anova_count(d, "y", "g", plots = FALSE)

  expect_identical(fit$model_type, "negbin")
  expect_gt(fit$dispersion, 1.5)                    # the Poisson statistic
  expect_lt(fit$model_dispersion, 2)                # the model actually returned
  expect_equal(fit$model_dispersion,
               sum(stats::residuals(fit$model, type = "pearson")^2) /
                 stats::df.residual(fit$model))
})

test_that("an ill-determined negative binomial theta is reported as such", {
  skip_if_not_installed("MASS")
  # Exactly Poisson: theta has nothing to estimate and runs away.
  d <- withr::with_seed(41, data.frame(
    g = factor(rep(c("a", "b", "c"), each = 200)),
    y = stats::rpois(600, lambda = rep(c(6, 6.5, 7), each = 200))))
  fit <- suppressWarnings(
    anova_count(d, "y", "g", model = "negbin", plots = FALSE))

  expect_identical(fit$model_type, "negbin")
  expect_true(any(grepl("not meaningfully overdispersed|not well determined|did not converge",
                        fit$notes)))
  # and the note must not present a runaway theta as a three-decimal estimate
  runaway <- fit$model$theta > 1000 ||
    fit$model$theta / fit$model$SE.theta < 2
  expect_true(runaway)
})

test_that("anova_glm says when robust standard errors force Wald intervals", {
  skip_if_not_installed("sandwich")
  fit <- anova_glm(fx_oneway(), "value", "group", vcov_type = "HC3",
                   ci_method = "profile", plots = FALSE)
  expect_identical(unique(fit$assumptions$coefficients$ci_method), "wald")
  expect_true(any(grepl("Robust standard errors", fit$notes)))
})

test_that("anova_rm accounts for every row that left the analysis", {
  skip_if_not_installed("afex")
  d <- fx_repeated()
  d$score[1] <- NA                                   # 1 missing
  fit <- anova_rm(d, "score", subject = "id", within = "time",
                  between = "arm", plots = FALSE)

  # 1 row missing plus the other 2 rows of that now-incomplete subject
  expect_identical(fit$n_removed, 3L)
  expect_identical(fit$n_removed_missing, 1L)
  expect_identical(fit$n_removed_unbalanced, 2L)
  expect_equal(nrow(fit$data_used) + fit$n_removed, nrow(d))
})

test_that("a character grouping column gets the same reference level everywhere", {
  d <- data.frame(g = rep(c("apple", "Banana", "cherry"), each = 20),
                  y = withr::with_seed(15, stats::rnorm(60)))
  old <- Sys.getlocale("LC_COLLATE")
  on.exit(suppressWarnings(Sys.setlocale("LC_COLLATE", old)), add = TRUE)

  terms_in <- function(loc) {
    ok <- suppressWarnings(Sys.setlocale("LC_COLLATE", loc))
    if (!nzchar(ok)) return(NULL)
    levels(anova_glm(d, "y", "g", plots = FALSE)$data_used$g)
  }
  a <- terms_in("C")
  skip_if(is.null(a), "the C collation could not be set")
  expect_identical(a[1], "Banana")   # C collation: upper case sorts first
  # The decisive comparison needs a collation that differs from C's.
  b <- NULL
  for (loc in c("en_US.UTF-8", "en_GB.UTF-8", "de_DE.UTF-8", "fr_FR.UTF-8",
                "en_US.utf8", "English_United States.1252")) {
    b <- terms_in(loc)
    if (!is.null(b)) break
  }
  skip_if(is.null(b), "no locale with a non-C collation is installed")
  expect_identical(a, b)
})

test_that("the mvt adjustment, which is randomised, is not offered", {
  expect_error(anova_glm(fx_oneway(), "value", "group", adjust = "mvt"),
               "must be one of")
  seed_before <- if (exists(".Random.seed", .GlobalEnv)) {
    get(".Random.seed", .GlobalEnv)
  } else NULL
  invisible(anova_glm(fx_oneway(), "value", "group", plots = FALSE))
  if (!is.null(seed_before)) {
    expect_identical(get(".Random.seed", .GlobalEnv), seed_before)
  }
})

test_that("print honours its digits argument regardless of the global option", {
  fit <- anova_welch(fx_oneway(), "value", "group", plots = FALSE)
  old <- options(digits = 3)
  on.exit(options(old), add = TRUE)
  out <- utils::capture.output(print(fit, digits = 7))
  stat_line <- grep("^1 ", out, value = TRUE)[1]
  expect_false(is.na(stat_line))
  expect_match(stat_line, "[0-9]\\.[0-9]{4,}")
})

test_that("per-group normality reports the n the test actually used", {
  v <- c(1, 2, NA, 4, 5, 1, 2, 3, 4, 5)
  g <- factor(rep(c("a", "b"), each = 5))
  tab <- anovakit:::.normality_by_group(v, g)
  expect_identical(tab$n, c(4L, 5L))
})

# --- Validation, row accounting and weights -----------------------------------

test_that("Inf in a count response is dropped, not fed to the integer check", {
  d <- fx_counts()
  d$count[c(3L, 40L)] <- Inf
  fit <- anova_count(d, "count", "g1", plots = FALSE)
  expect_identical(fit$n_removed, 2L)
  expect_identical(nrow(fit$data_used) + fit$n_removed, nrow(d))
  expect_true(all(is.finite(fit$data_used$count)))

  # A genuinely non-integer count is still refused
  d2 <- fx_counts(); d2$count[1] <- 1.5
  expect_error(anova_count(d2, "count", "g1", plots = FALSE), "whole numbers")

  # ... and a column with no finite value at all says so
  d3 <- fx_counts(); d3$count <- NA_real_
  expect_error(anova_count(d3, "count", "g1", plots = FALSE), "no finite values")
})

test_that("anova_rm's row arithmetic closes when afex aggregates duplicates", {
  skip_if_not_installed("afex")
  d <- fx_repeated()
  dup <- d[d$id %in% levels(d$id)[1:2], ]                 # 2 subjects x 3 cells
  d2 <- rbind(d, dup)                                     # each appears twice
  fit <- suppressMessages(
    anova_rm(d2, "score", subject = "id", within = "time", between = "arm",
             plots = FALSE, fun_aggregate = mean))

  expect_identical(nrow(fit$data_used) + fit$n_removed, nrow(d2))
  expect_identical(fit$n_removed_aggregated, nrow(dup))
  expect_true(any(grepl("aggregated", fit$notes)))
  expect_true(any(grepl("included in \\$n_removed", fit$notes)))
})

test_that("a negative weight is caught even on a row that would be dropped", {
  d <- fx_oneway()
  d$w <- 1
  d$w[1] <- -1
  d$value[1] <- NA                       # the offending row leaves anyway
  expect_error(anova_glm(d, "value", "group", weights = "w", plots = FALSE),
               "negative value")
})

test_that("weights that leave nothing to fit on produce a named error", {
  d <- fx_oneway()
  d$w <- 0
  expect_error(anova_glm(d, "value", "group", weights = "w", plots = FALSE),
               "Every weight")

  d$w[1] <- 1
  expect_error(anova_glm(d, "value", "group", weights = "w", plots = FALSE),
               "non-zero weight")

  # A non-finite weight is missingness, not an error: the row is dropped
  d2 <- fx_oneway(); d2$w <- 1; d2$w[1] <- NA
  fit <- anova_glm(d2, "value", "group", weights = "w", plots = FALSE)
  expect_identical(fit$n_removed, 1L)
})

test_that("posthoc = FALSE explains the empty tables it leaves behind", {
  callers <- list(
    welch  = function() anova_welch(fx_oneway(), "value", "group",
                                    posthoc = FALSE, plots = FALSE),
    kw     = function() anova_kw(fx_oneway(), "value", "group",
                                 posthoc = FALSE, plots = FALSE),
    glm    = function() anova_glm(fx_oneway(), "value", "group",
                                  posthoc = FALSE, plots = FALSE),
    bin    = function() anova_bin(fx_binary(), "y", "g",
                                  posthoc = FALSE, plots = FALSE),
    count  = function() anova_count(fx_counts(), "count", "g1",
                                    posthoc = FALSE, plots = FALSE),
    ancova = function() anova_ancova(fx_ancova(), "dv", "iv", "cov",
                                     posthoc = FALSE, plots = FALSE),
    manova = function() anova_manova(fx_multivariate(), c("score1", "score2"),
                                     "g", posthoc = FALSE, assumptions = FALSE,
                                     plots = FALSE)
  )
  for (nm in names(callers)) {
    fit <- callers[[nm]]()
    expect_null(fit$posthoc, info = nm)
    expect_true(any(grepl("posthoc = FALSE", fit$notes, fixed = TRUE)),
                info = nm)
    # Exactly one note per fit, even when the function loops over responses
    expect_identical(sum(grepl("posthoc = FALSE", fit$notes, fixed = TRUE)), 1L,
                     info = nm)
  }

  # Welch's effect sizes go with its comparisons, and the note says so
  w <- callers$welch()
  expect_null(w$effect_sizes)
  expect_true(any(grepl("Hedges", w$notes)))

  # The others keep the effect sizes that come from the model itself
  expect_false(is.null(callers$kw()$effect_sizes))
  expect_false(is.null(callers$bin()$effect_sizes))
  expect_false(is.null(callers$count()$effect_sizes))
})

test_that("a data column called wts cannot displace the prior weights", {
  d <- fx_binary()
  d$w <- rep(c(1, 3), length.out = nrow(d))
  d$wts <- 0                            # a decoy that used to be the weights
  with_wts <- anova_bin(d, "y", "g", weights = "w", plots = FALSE)
  d2 <- d; d2$wts <- NULL
  without <- anova_bin(d2, "y", "g", weights = "w", plots = FALSE)
  expect_equal(stats::coef(with_wts$model), stats::coef(without$model))
  expect_equal(unname(stats::weights(with_wts$model, "prior")), d$w)
})

test_that("a rank-deficient response set is refused by name", {
  d <- fx_multivariate()
  d$score3 <- d$score1 + d$score2         # exactly collinear
  expect_error(
    anova_manova(d, c("score1", "score2", "score3"), "g", plots = FALSE),
    "collinear|could not be computed")
})

test_that("the canonical plot is coloured by the term it describes", {
  d <- fx_multivariate()
  d$site <- factor(rep(c("n", "s"), length.out = nrow(d)))
  fit <- anova_manova(d, c("score1", "score2"), c("g", "site"),
                      assumptions = FALSE)
  skip_if(is.null(fit$plots$canonical))
  expect_identical(fit$canonical_term, "g")
  # One colour per level of g, not per cell of g x site
  expect_identical(nlevels(droplevels(as.factor(fit$plots$canonical$data$.group))),
                   nlevels(d$g))
})

test_that("a grouping column named df displaces only emmeans' own df", {
  d <- fx_oneway()
  names(d)[names(d) == "group"] <- "df"
  fit <- anova_glm(d, "value", "df", plots = FALSE)
  expect_true(is.factor(fit$emmeans$df) || is.character(fit$emmeans$df))
  expect_setequal(as.character(fit$emmeans$df), levels(as.factor(d$df)))
  expect_true("emm_df" %in% names(fit$emmeans))
  expect_true(is.numeric(fit$emmeans$emm_df))
  expect_identical(anyDuplicated(names(fit$emmeans)), 0L)
})

test_that("emmeans plots find their columns when a grouping name was protected", {
  # The estimate itself is the hardest case: the plot's y aesthetic is the
  # displaced column, so resolving it by bare name would plot the grouping
  # factor against itself.
  for (nm in c("se", "estimate", "conf_low")) {
    d <- fx_oneway()
    names(d)[names(d) == "group"] <- nm
    fit <- anova_glm(d, "value", nm)
    expect_s3_class(fit$plots$emmeans, "ggplot")
    expect_true(paste0("emm_", nm) %in% names(fit$emmeans), info = nm)
    expect_true(is.numeric(fit$emmeans[[paste0("emm_", nm)]]), info = nm)
    expect_setequal(as.character(fit$emmeans[[nm]]), levels(as.factor(d[[nm]])))
    expect_true(any(grepl("renamed", fit$notes)), info = nm)
    # The plot's y values are the estimates, not the grouping labels
    py <- fit$plots$emmeans$data
    ycol <- anovakit:::.emm_col(py, "estimate", protect = nm)
    expect_true(is.numeric(py[[ycol]]), info = nm)
  }
})

# --- Reported counts, saturated models and captured conditions ----------------

test_that("the suppressed-comparison count is the number that would be made", {
  # anova_manova stacks $posthoc across every response, so the promise must
  # count them all, not just the first response's share
  d <- fx_multivariate()
  d$score3 <- withr::with_seed(31, stats::rnorm(nrow(d)))
  responses <- c("score1", "score2", "score3")
  off <- anova_manova(d, responses, "g", posthoc = FALSE, assumptions = FALSE,
                      plots = FALSE)
  on  <- anova_manova(d, responses, "g", posthoc = TRUE, assumptions = FALSE,
                      plots = FALSE)
  promised <- as.integer(sub(".*would have been ([0-9]+)\\..*", "\\1",
                             grep("would have been", off$notes, value = TRUE)))
  expect_identical(promised, nrow(on$posthoc))

  # ... and the same holds for the single-response functions
  for (f in list(
    function(p) anova_welch(fx_oneway(), "value", "group", posthoc = p,
                            plots = FALSE),
    function(p) anova_kw(fx_oneway(), "value", "group", posthoc = p,
                         plots = FALSE),
    function(p) anova_glm(fx_oneway(), "value", "group", posthoc = p,
                          plots = FALSE))) {
    note <- grep("would have been", f(FALSE)$notes, value = TRUE)
    promised <- as.integer(sub(".*would have been ([0-9]+)\\..*", "\\1", note))
    expect_identical(promised, nrow(f(TRUE)$posthoc))
  }
})

test_that("the suppressed-comparison count survives more cells than an integer", {
  # choose(65537, 2) exceeds .Machine$integer.max, and sprintf("%d") rejects it
  note <- anovakit:::.no_posthoc_note(65537L)
  expect_gt(choose(65537L, 2L), .Machine$integer.max)
  expect_match(note, "2147516416", fixed = TRUE)
  expect_false(grepl("e+", note, fixed = TRUE))
  expect_match(anovakit:::.no_posthoc_note(3L), "3.", fixed = TRUE)
  expect_match(anovakit:::.no_posthoc_note(NA_integer_), "posthoc = FALSE")
})

test_that("a saturated model reports NA inference instead of warning", {
  d <- fx_oneway()
  d$w <- 0
  d$w[c(1L, 30L)] <- 1                  # 2 positive weights, 3 parameters
  expect_silent(fit <- anova_glm(d, "value", "group", weights = "w",
                                 plots = FALSE, posthoc = FALSE))
  expect_true(any(grepl("saturated", fit$notes)))
  expect_true(any(grepl("undefined or degenerate", fit$notes)))
  expect_true(is.na(anovakit:::.coef_crit(fit$model, 0.95)))

  # A model with residual degrees of freedom is unaffected
  ok <- anova_glm(fx_oneway(), "value", "group", plots = FALSE)
  expect_true(is.finite(anovakit:::.coef_crit(ok$model, 0.95)))
  expect_false(any(grepl("saturated", ok$notes)))
})

test_that("third-party warnings and messages reach $notes, not the console", {
  # Non-integer prior weights make glm() warn about "non-integer #successes"
  d <- fx_binary()
  d$w <- withr::with_seed(41, stats::runif(nrow(d), 0.5, 2))
  expect_silent(fit <- anova_bin(d, "y", "g", weights = "w", plots = FALSE))
  expect_true(any(grepl("stats::glm\\(\\) reported", fit$notes)))
  expect_true(any(grepl("non-integer", fit$notes)))

  # An aliased design makes car::Anova() emit a message()
  ali <- fx_twoway()
  ali$copy <- ali$g1
  expect_silent(a <- anova_glm(ali, "value", c("g1", "copy"), plots = FALSE))
  expect_false(is.null(a$anova))
})

test_that("afex's own warnings are captured rather than printed", {
  skip_if_not_installed("afex")
  d <- fx_repeated()
  d2 <- rbind(d, d[d$id %in% levels(d$id)[1L], ])

  # Without fun_aggregate, afex warns that it is aggregating. The warning must
  # reach $notes and not the console.
  expect_silent(loud <- anova_rm(d2, "score", subject = "id", within = "time",
                                 between = "arm", plots = FALSE))
  said <- grep("aov_ez\\(\\) reported", loud$notes, value = TRUE)
  expect_length(said, 1L)
  expect_match(said, "aggregating data")
  expect_false(grepl("\n", said, fixed = TRUE))   # flattened to one line

  # With fun_aggregate given, afex says nothing, so neither does $notes
  expect_silent(quiet <- anova_rm(d2, "score", subject = "id", within = "time",
                                  between = "arm", plots = FALSE,
                                  fun_aggregate = mean))
  expect_false(any(grepl("aov_ez\\(\\) reported", quiet$notes)))

  # The package's own account of the aggregation is there either way
  for (f in list(loud, quiet)) {
    expect_true(any(grepl("afex aggregated", f$notes)))
    expect_identical(nrow(f$data_used) + f$n_removed, nrow(d2))
  }
})

test_that("the collector returns the value and never re-raises", {
  q <- anovakit:::.collect_conditions({
    warning("first"); message("second"); 42
  })
  expect_identical(q$value, 42)
  expect_identical(q$said, c("first", "second"))

  e <- anovakit:::.collect_conditions(stop("boom"))
  expect_s3_class(e$value, "error")
  expect_identical(conditionMessage(e$value), "boom")

  expect_identical(anovakit:::.said_note(character(0), "x"), character(0))
  expect_match(anovakit:::.said_note(c("a", "b"), "x()"), "x\\(\\) reported: a / b")
})

test_that("the repeated-measures duplicate count cannot be fooled by labels", {
  skip_if_not_installed("afex")
  d <- expand.grid(id = c("a", "ab"), time = c("bc", "c"),
                   stringsAsFactors = FALSE)
  d$score <- c(1, 2, 3, 4)
  d$id <- factor(d$id); d$time <- factor(d$time)
  # id "a" + cell "bc" and id "ab" + cell "c" both paste to "abc" with sep = ""
  bal <- anovakit:::.balance_subjects(d, "id", "time")
  expect_identical(bal$duplicates, 0L)
})

# --- Reference distributions, captured conditions and note truthfulness -------

test_that("a negative binomial fit uses z, as MASS does, not t", {
  skip_if_not_installed("MASS")
  d <- withr::with_seed(8, {
    n <- 24
    g <- factor(rep(c("a", "b"), each = n / 2))
    data.frame(g = g, cnt = c(stats::rnbinom(n / 2, size = 0.5, mu = 5),
                              stats::rnbinom(n / 2, size = 0.5, mu = 14)))
  })
  fit <- anova_count(d, "cnt", "g", model = "negbin", plots = FALSE,
                     posthoc = FALSE)
  ref <- MASS::glm.nb(cnt ~ g, data = d)
  rc <- stats::coef(summary(ref))

  # MASS::glm.nb fixes the dispersion at 1 and reports a z value; using a t
  # reference with df.residual made the p-value and the interval too wide
  expect_false(anovakit:::.estimates_dispersion(ref))
  expect_equal(anovakit:::.coef_crit(ref, 0.95), stats::qnorm(0.975))
  expect_equal(fit$effect_sizes$p_value[1], unname(rc["gb", 4]))
  expect_equal(fit$effect_sizes$conf_low[1],
               exp(rc["gb", 1] - stats::qnorm(0.975) * rc["gb", 2]))

  # ... while Gaussian and Gamma still use t, and binomial and Poisson still z
  ok <- anova_glm(fx_oneway(), "value", "group", plots = FALSE)
  expect_true(anovakit:::.estimates_dispersion(ok$model))
  b <- anova_bin(fx_binary(), "y", "g", plots = FALSE)
  expect_false(anovakit:::.estimates_dispersion(b$model))
})

test_that("routine third-party chatter does not become a note", {
  skip_if_not_installed("afex")
  # afex announces the contrasts it sets on every between-subjects fit. That
  # is the default this wrapper asks for, so it is not news.
  fit <- anova_rm(fx_repeated(), "score", subject = "id", within = "time",
                  between = "arm", plots = FALSE)
  expect_false(any(grepl("Contrasts set to", fit$notes)))
  expect_identical(anovakit:::.said_note(
    "Contrasts set to contr.sum for the following variables: arm", "afex"),
    character(0))

  # A real warning still gets through
  expect_match(anovakit:::.said_note("glm.fit: algorithm did not converge",
                                         "stats::glm()"),
               "did not converge")
})

test_that("no fitting or marginal-means call leaks a condition to the console", {
  # A vanishing offset makes glm() warn about numerically zero fitted rates
  z <- withr::with_seed(31, {
    g <- factor(rep(c("a", "b", "c"), each = 30))
    data.frame(g = g,
               cnt = c(stats::rpois(30, 6), stats::rpois(30, 9), rep(0L, 30)),
               expo = c(stats::runif(60, 0.9, 1.1), rep(1e-20, 30)))
  })
  expect_silent(cf <- anova_count(z, "cnt", "g", offset = "expo",
                                  plots = FALSE, posthoc = FALSE))
  expect_true(any(grepl("numerically 0", cf$notes)))

  # A perfect fit makes summary.lm() warn from inside emmeans' reference grid
  a <- withr::with_seed(1, {
    g <- factor(rep(c("a", "b", "c"), each = 20))
    x <- stats::rnorm(60, 10, 2)
    data.frame(g = g, x = x, y = 3 + 0.5 * x)
  })
  expect_silent(af <- anova_ancova(a, "y", "g", "x", plots = FALSE,
                                   posthoc = FALSE))
  expect_true(any(grepl("perfect fit", af$notes)))
})

test_that("the ancova model-choice note is emitted once and is true", {
  aliased <- withr::with_seed(3, {
    g <- factor(rep(c("a", "b", "c"), each = 20))
    data.frame(g = g, x = as.numeric(g), y = stats::rnorm(60))
  })
  ok <- fx_ancova()
  grid <- expand.grid(fi = list(NULL, TRUE, FALSE), cc = c(TRUE, FALSE),
                      dat = c("ok", "aliased"), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(grid))) {
    d <- if (grid$dat[[i]] == "ok") ok else aliased
    fit <- anova_ancova(d, if (grid$dat[[i]] == "ok") "dv" else "y",
                        if (grid$dat[[i]] == "ok") "iv" else "g",
                        if (grid$dat[[i]] == "ok") "cov" else "x",
                        force_interaction = grid$fi[[i]],
                        center_covariates = grid$cc[[i]], plots = FALSE)
    info <- sprintf("force_interaction=%s centre=%s data=%s",
                    if (is.null(grid$fi[[i]])) "NULL" else grid$fi[[i]],
                    grid$cc[[i]], grid$dat[[i]])
    # Exactly one note says which model was fitted
    chose <- grepl("additive model|covariate-by-group interaction", fit$notes)
    expect_identical(sum(chose), 1L, info = info)
    note <- fit$notes[chose]
    interaction_fitted <- any(grepl(":", attr(stats::terms(fit$model),
                                              "term.labels"), fixed = TRUE))
    expect_identical(grepl("interaction was (retained|fitted)", note),
                     interaction_fitted, info = info)
    # It never blames the slopes test for a choice the caller made
    if (!is.null(grid$fi[[i]])) {
      expect_match(note, "force_interaction", info = info)
    }
    # ... and it never says "gave could not be computed"
    expect_false(grepl("gave could not be computed", note), info = info)
    # The stated reference point matches the centring that was applied
    if (grepl("comparison at", note)) {
      expect_identical(grepl("at the covariate mean", note), grid$cc[[i]],
                       info = info)
    }
  }
})

test_that("McFadden's R squared is finite and correct for non-Gaussian fits", {
  d <- fx_binary()
  fit <- anova_glm(d, "y", "g", family = "binomial", plots = FALSE)
  mc <- fit$effect_sizes$estimate[fit$effect_sizes$measure == "mcfadden_r2"]
  ref <- stats::glm(y ~ g, family = stats::binomial(), data = fit$data_used)
  null <- stats::glm(y ~ 1, family = stats::binomial(), data = fit$data_used)
  expect_equal(mc, as.numeric(1 - stats::logLik(ref) / stats::logLik(null)))
  expect_true(is.finite(mc))

  # With prior weights, the null model must carry them too, or the reference
  # likelihood is not comparable. (anova_glm() has no offset argument, so the
  # null model never needs one.)
  cd <- fx_counts()
  cd$w <- withr::with_seed(5, sample(1:3, nrow(cd), replace = TRUE))
  cf <- anova_glm(cd, "count", "g1", family = "poisson", weights = "w",
                  plots = FALSE)
  full <- stats::glm(count ~ g1, family = stats::poisson(), data = cd, weights = w)
  null_w <- stats::glm(count ~ 1, family = stats::poisson(), data = cd, weights = w)
  expect_equal(cf$effect_sizes$estimate[cf$effect_sizes$measure == "mcfadden_r2"],
               as.numeric(1 - stats::logLik(full) / stats::logLik(null_w)))

  # A quasi family has no likelihood: NA, not an error
  qf <- anova_glm(fx_counts(), "count", "g1", family = "quasipoisson",
                  plots = FALSE)
  expect_true(is.na(qf$effect_sizes$estimate[
    qf$effect_sizes$measure == "mcfadden_r2"]))
})

test_that("a matrix response is refused by name, not by an internal error", {
  d <- fx_binary()
  d$Y <- cbind(succ = d$y, fail = 1L - d$y)
  expect_error(anova_glm(d, "Y", "g", family = "binomial", plots = FALSE),
               "2-column matrix")

  # A ONE-column matrix is what scale() returns and is a common way to build a
  # column. Every other function accepts it, so anova_glm must too.
  z <- fx_oneway()
  z$value <- scale(z$value)
  expect_true(is.matrix(z$value))
  fit <- anova_glm(z, "value", "group", plots = FALSE)
  ref <- stats::anova(stats::lm(as.numeric(value) ~ group,
                                data = transform(z, value = as.numeric(value))))
  expect_equal(fit$anova$p_value[1], ref[["Pr(>F)"]][1])
  for (f in list(
    function() anova_welch(z, "value", "group", plots = FALSE),
    function() anova_kw(z, "value", "group", plots = FALSE),
    function() anova_count(transform(fx_counts(), count = scale(count) * 0 +
                                       fx_counts()$count),
                           "count", "g1", plots = FALSE))) {
    expect_s3_class(f(), "anovakit_fit")
  }
})

# --- Reference distributions, captured conditions and conditional notes -------

test_that("only MASS's negbin fit fixes the dispersion; the family alone does not", {
  skip_if_not_installed("MASS")
  d <- withr::with_seed(99, {
    n <- 300
    g <- factor(rep(c("a", "b", "c"), each = n / 3))
    data.frame(g = g, cnt = stats::rnbinom(
      n, mu = c(a = 4, b = 8, c = 12)[as.character(g)], size = 1.5))
  })

  # MASS::glm.nb() estimates theta then fixes the dispersion at 1 -> z
  nb <- MASS::glm.nb(cnt ~ g, data = d)
  expect_s3_class(nb, "negbin")
  expect_false(anovakit:::.estimates_dispersion(nb))
  expect_equal(anovakit:::.coef_crit(nb, 0.95), stats::qnorm(0.975))

  # glm(family = negative.binomial(theta)) carries the SAME family string but
  # is an ordinary glm whose dispersion summary.glm() estimates -> t
  fx <- stats::glm(cnt ~ g, data = d, family = MASS::negative.binomial(2))
  # Both carry the same family string; only the class distinguishes them, so
  # matching on the string alone gave the fixed-theta fit a z reference
  expect_true(grepl("^Negative Binomial", stats::family(nb)$family))
  expect_true(grepl("^Negative Binomial", stats::family(fx)$family))
  expect_false(inherits(fx, "negbin"))
  expect_true(anovakit:::.estimates_dispersion(fx))

  ref <- stats::coef(summary(fx))
  got <- anova_glm(d, "cnt", "g", family = MASS::negative.binomial(2),
                   plots = FALSE, ci_method = "wald")$assumptions$coefficients
  expect_equal(got$p_value, unname(ref[, 4]))
  expect_identical(unique(got$distribution), "t")
})

test_that("aliasing is visible in the result, not only in a suppressed message", {
  # A covariate constant within one group makes cov:g partly inestimable
  d <- withr::with_seed(31, {
    g <- factor(rep(c("a", "b", "c"), each = 30))
    cov <- ifelse(g == "c", 7, stats::rnorm(60, 7, 2))
    data.frame(g = g, cov = cov,
               y = stats::rnorm(90, 2 + 0.4 * cov + as.integer(g), 1))
  })
  for (ty in c("II", "III")) {
    fit <- anova_ancova(d, "y", "g", "cov", force_interaction = TRUE,
                        type = ty, plots = FALSE)
    expect_true(any(is.na(stats::coef(fit$model))), info = ty)
    expect_true(any(grepl("could not be estimated", fit$notes)), info = ty)
  }

  # anova_manova reports it once, not once per response
  m <- fx_multivariate()
  m$dup <- m$g
  mf <- anova_manova(m, c("score1", "score2"), c("g", "dup"),
                     assumptions = FALSE, plots = FALSE)
  expect_lte(sum(grepl("could not be estimated", mf$notes)), 1L)
})

test_that("the ancova centring note claims an effect only where there is one", {
  d <- withr::with_seed(61, {
    n <- 120
    g <- factor(rep(c("ctl", "trt"), each = n / 2))
    base <- stats::rnorm(n, 100, 5)
    data.frame(g = g, base = base,
               score = 2 * base + 6 * (g == "trt") + stats::rnorm(n, 0, 3))
  })
  grp_row <- function(f) f$anova$p_value[f$anova$term == "g"]

  # Additive model: the group row is bit-identical either way, so the note
  # must not send the reader away from it
  # expect_equal, not expect_identical: centring changes the design matrix
  # from X to X - mean, so the QR decomposition agrees to about 13 significant
  # figures rather than bit-for-bit, and how many depends on the BLAS.
  add <- lapply(c(TRUE, FALSE), function(cc)
    anova_ancova(d, "score", "g", "base", force_interaction = FALSE,
                 center_covariates = cc, plots = FALSE))
  expect_equal(grp_row(add[[1]]), grp_row(add[[2]]))
  expect_match(grep("NOT centred", add[[2]]$notes, value = TRUE), "unaffected")
  expect_false(any(grepl("Read \\$emmeans rather than the ANOVA row",
                         add[[2]]$notes)))

  # Type II with the interaction: also invariant
  t2 <- lapply(c(TRUE, FALSE), function(cc)
    anova_ancova(d, "score", "g", "base", force_interaction = TRUE,
                 center_covariates = cc, type = "II", plots = FALSE))
  expect_equal(grp_row(t2[[1]]), grp_row(t2[[2]]))
  expect_match(grep("NOT centred", t2[[2]]$notes, value = TRUE), "unaffected")

  # Type III with the interaction: it genuinely moves, and the note says so
  t3 <- lapply(c(TRUE, FALSE), function(cc)
    anova_ancova(d, "score", "g", "base", force_interaction = TRUE,
                 center_covariates = cc, type = "III", plots = FALSE))
  expect_false(isTRUE(all.equal(grp_row(t3[[1]]), grp_row(t3[[2]]))))
  expect_match(grep("NOT centred", t3[[2]]$notes, value = TRUE),
               "read \\$emmeans instead")
  expect_match(grep("mean-centred", t3[[1]]$notes, value = TRUE),
               "Type III group row")
})

test_that("an untestable slopes assumption is NA, not a claim of heterogeneity", {
  d <- withr::with_seed(9, {
    n <- 120
    g <- factor(rep(c("ctl", "trt"), each = n / 2))
    base <- ifelse(g == "trt", 100, stats::rnorm(n / 2, 100, 5))
    data.frame(g = g, base = base,
               score = 2 * base + 6 * (g == "trt") + stats::rnorm(n, 0, 3))
  })
  fit <- anova_ancova(d, "score", "g", "base", plots = FALSE)
  expect_true(is.na(fit$slopes_test$homogeneous))
  expect_true(is.na(fit$slopes_test$p_value))
  expect_match(paste(fit$notes, collapse = " "), "could not be computed")

  # A computable test still gives TRUE/FALSE
  ok <- anova_ancova(fx_ancova(), "dv", "iv", "cov", plots = FALSE)
  expect_false(is.na(ok$slopes_test$homogeneous))
})

test_that("emtrends does not write to the console either", {
  sat <- data.frame(g = factor(rep(c("a", "b"), each = 3)),
                    base = c(1, 2, 3, 1, 2, 3),
                    score = c(2, 4, 6, 3, 5, 7))
  expect_silent(fit <- anova_ancova(sat, "score", "g", "base",
                                    force_interaction = TRUE, plots = FALSE))
  expect_true(any(grepl("perfect fit", fit$notes)))
})

test_that("a runaway negative binomial theta is reported as non-convergence", {
  skip_if_not_installed("MASS")
  d <- withr::with_seed(77, {
    n <- 200
    g <- factor(rep(c("a", "b"), each = n / 2))
    data.frame(g = g, y = stats::rnbinom(
      n, mu = c(a = 5, b = 8)[as.character(g)], size = 0.7))
  })
  real <- MASS::glm.nb
  notes <- testthat::with_mocked_bindings(
    anova_count(d, "y", "g", model = "negbin", plots = FALSE)$notes,
    glm.nb = function(...) {
      m <- real(...)
      warning("alternation limit reached")
      m
    },
    .package = "MASS")
  expect_true(any(grepl("did not converge", notes)))

  # Without the warning it is reported as an ordinary estimate
  clean <- anova_count(d, "y", "g", model = "negbin", plots = FALSE)$notes
  expect_false(any(grepl("did not converge", clean)))
  expect_true(any(grepl("dispersion parameter theta", clean)))
})
