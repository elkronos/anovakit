# anova_rm(): repeated measures ANOVA, fitted with afex::aov_ez().
#
# Every block skips when afex is absent. Covers the table and sphericity read
# off the fit, the removal of subjects with an incomplete within-subject
# design, and the row accounting that must close between the input and
# $data_used.

test_that("the ANOVA table reproduces afex::aov_ez under every correction", {
  skip_if_not_installed("afex")
  d <- fx_repeated()
  ref <- suppressMessages(afex::aov_ez(id = "id", dv = "score", data = d,
                                       within = "time", between = "arm"))
  ref_tab <- as.data.frame(ref$anova_table)          # afex's default, GG
  for (corr in c("GG", "HF", "none")) {
    fit <- anova_rm(d, "score", subject = "id", within = "time",
                    between = "arm", correction = corr, plots = FALSE)
    r <- suppressWarnings(as.data.frame(
      stats::anova(ref, correction = corr, es = "pes")))
    expect_identical(fit$anova$term, rownames(r), info = corr)
    expect_equal(fit$anova$num_df, r[["num Df"]], info = corr)
    expect_equal(fit$anova$den_df, r[["den Df"]], info = corr)
    expect_equal(fit$anova$mse, r[["MSE"]], info = corr)
    expect_equal(fit$anova$statistic, r[["F"]], info = corr)
    expect_equal(fit$anova$partial_eta_sq, r[["pes"]], info = corr)
    expect_equal(fit$anova$p_value, r[["Pr(>F)"]], info = corr)
    expect_identical(attr(fit$anova, "correction"), corr)
    if (corr == "GG") expect_equal(fit$anova$p_value, ref_tab[["Pr(>F)"]])
  }
  # The three corrections really differ on this fixture (HF eps < 1)
  ps <- vapply(c("GG", "HF", "none"), function(corr) {
    anova_rm(d, "score", "id", "time", between = "arm", correction = corr,
             plots = FALSE)$anova$p_value[2]
  }, numeric(1))
  expect_true(ps[["GG"]] > ps[["HF"]] && ps[["HF"]] > ps[["none"]])

  # Generalised eta squared is afex's (every factor manipulated here)
  fit <- anova_rm(d, "score", subject = "id", within = "time",
                  between = "arm", plots = FALSE)
  ges <- as.data.frame(stats::anova(ref, es = "ges"))
  expect_equal(fit$effect_sizes$generalised_eta_sq, ges[["ges"]])
  expect_false(isTRUE(all.equal(fit$effect_sizes$generalised_eta_sq,
                                fit$effect_sizes$partial_eta_sq)))
})

test_that("sphericity is reported, not silently dropped", {
  skip_if_not_installed("afex")
  d <- fx_repeated()
  fit <- anova_rm(d, "score", subject = "id", within = "time",
                  between = "arm", plots = FALSE)

  expect_false(is.null(fit$sphericity))
  expect_true(all(c("term", "mauchly_w", "p_value") %in% names(fit$sphericity)))
  expect_true(all(c("gg_epsilon", "hf_epsilon") %in% names(fit$sphericity)))

  ref <- suppressMessages(afex::aov_ez(id = "id", dv = "score", data = d,
                                       within = "time", between = "arm"))
  s <- suppressWarnings(summary(ref$Anova, multivariate = FALSE))
  sph <- as.matrix(unclass(s$sphericity.tests))
  expect_equal(fit$sphericity$mauchly_w, as.numeric(sph[, 1]))
  expect_equal(fit$sphericity$p_value, as.numeric(sph[, 2]))

  # The epsilons and corrected p-values, value by value
  eps <- as.matrix(unclass(s$pval.adjustments))[fit$sphericity$term, , drop = FALSE]
  expect_equal(fit$sphericity$gg_epsilon, unname(eps[, "GG eps"]))
  expect_equal(fit$sphericity$hf_epsilon, pmin(1, unname(eps[, "HF eps"])))
  expect_equal(fit$sphericity$hf_epsilon_raw, unname(eps[, "HF eps"]))
  expect_equal(fit$sphericity$p_gg, unname(eps[, "Pr(>F[GG])"]))
  expect_equal(fit$sphericity$p_hf, unname(eps[, "Pr(>F[HF])"]))
})

test_that("diagnostic plots are produced for a within-subjects design", {
  skip_if_not_installed("afex")
  d <- fx_repeated()
  fit <- anova_rm(d, "score", subject = "id", within = "time", between = "arm")

  expect_named(fit$plots, c("residuals", "qq", "index", "emmeans"))
  expect_true(all(vapply(fit$plots, inherits, logical(1), "ggplot")))
  expect_length(anova_rm(d, "score", subject = "id", within = "time",
                         plots = FALSE)$plots, 0L)
})

test_that("a factor response is read by its labels, not its level codes", {
  skip_if_not_installed("afex")
  d <- fx_repeated()
  d$score <- factor(round(d$score, 0))
  fit <- anova_rm(d, "score", subject = "id", within = "time", plots = FALSE)

  expect_equal(fit$data_used$score, as.numeric(as.character(d$score)))
  expect_true(any(grepl("via its labels", fit$notes)))

  bad <- d
  levels(bad$score)[1] <- "not a number"
  expect_error(anova_rm(bad, "score", subject = "id", within = "time"),
               "cannot be read as numbers")
})

test_that("subjects with an incomplete design are dropped and named", {
  skip_if_not_installed("afex")
  d <- fx_repeated()
  d <- d[!(d$id == "1" & d$time == "t2"), ]
  fit <- anova_rm(d, "score", subject = "id", within = "time",
                  between = "arm", plots = FALSE)

  expect_identical(fit$subjects_dropped, "1")
  expect_equal(length(unique(fit$data_used$id)), 23L)
  expect_true(any(grepl("were removed because", fit$notes)))
})

test_that("a large design does not fail on the Shapiro-Wilk size limit", {
  skip_if_not_installed("afex")
  # Normality is tested per within-subject cell, so the 5000 limit applies to
  # the number of subjects.
  d <- withr::with_seed(21, {
    dd <- expand.grid(id = factor(seq_len(5001)),
                      time = factor(c("t1", "t2")))
    dd$y <- stats::rnorm(nrow(dd))
    dd
  })
  fit <- anova_rm(d, "y", subject = "id", within = "time", plots = FALSE)
  expect_s3_class(fit, "anovakit_fit")
  expect_true(all(is.na(fit$assumptions$normality$p_value)))
  expect_true(all(grepl("exceeds the 5000", fit$assumptions$normality$note)))
  expect_true(any(grepl("Shapiro-Wilk skipped", fit$notes)))
})

test_that("the afex object is returned intact", {
  skip_if_not_installed("afex")
  d <- fx_repeated()
  fit <- anova_rm(d, "score", subject = "id", within = "time",
                  between = "arm", plots = FALSE)

  expect_s3_class(fit$model, "afex_aov")
  expect_setequal(names(fit$model$data), c("long", "wide", "idata"))
  # afex was given internal names; $internal_names says which is which
  nm <- fit$internal_names
  expect_identical(nm$name, c("score", "id", "time", "arm"))
  w <- nm$internal[nm$name == "time"]
  g <- suppressMessages(emmeans::emmeans(fit$model, stats::as.formula(paste("~", w)),
                                         model = "multivariate"))
  expect_true(inherits(g, "emmGrid"))
  expect_equal(summary(g)$emmean, fit$emmeans$estimate)
  # ... while the grid returned for further use speaks the user's names
  expect_identical(names(fit$emmeans_object@levels), "time")
  expect_identical(fit$emmeans_object@levels$time, c("t1", "t2", "t3"))
})

test_that("the marginal-means plot maps every factor in the grid", {
  skip_if_not_installed("afex")
  d <- fx_repeated()
  fit <- anova_rm(d, "score", subject = "id", within = "time", between = "arm",
                  emm_specs = c("time", "arm"))
  p <- fit$plots$emmeans

  expect_true(all(c("x", "y", "colour", "group") %in% names(p$mapping)))
  expect_equal(nrow(p$data), 6L)
  expect_true(all(c("conf_low", "conf_high") %in% names(p$data)))
})

test_that("emm_specs is validated", {
  skip_if_not_installed("afex")
  d <- fx_repeated()
  expect_error(anova_rm(d, "score", subject = "id", within = "time",
                        emm_specs = "not_a_factor"),
               "must name within- or between-subject factors")
})

test_that("a single-level within factor is rejected", {
  skip_if_not_installed("afex")
  d <- fx_repeated()
  d1 <- d[d$time == "t1", ]
  expect_error(anova_rm(d1, "score", subject = "id", within = "time"),
               "at least 2 are required")
})

test_that("nothing is written to stdout at the default verbosity", {
  skip_if_not_installed("afex")
  d <- fx_repeated()
  expect_silent_stdout(anova_rm(d, "score", subject = "id", within = "time",
                                between = "arm"))
})
