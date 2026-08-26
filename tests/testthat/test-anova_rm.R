# anova_rm(): repeated measures ANOVA, fitted with afex::aov_ez().
#
# Every block skips when afex is absent. Covers the table and sphericity read
# off the fit, the removal of subjects with an incomplete within-subject
# design, and the row accounting that must close between the input and
# $data_used.

test_that("the ANOVA table reproduces afex::aov_ez", {
  skip_if_not_installed("afex")
  d <- fx_repeated()
  fit <- anova_rm(d, "score", subject = "id", within = "time",
                  between = "arm", plots = FALSE)
  ref <- suppressMessages(afex::aov_ez(id = "id", dv = "score", data = d,
                                       within = "time", between = "arm"))
  ref_tab <- as.data.frame(ref$anova_table)

  expect_equal(fit$anova$statistic, ref_tab[["F"]])
  expect_equal(fit$anova$p_value, ref_tab[["Pr(>F)"]])
  expect_setequal(fit$anova$term, rownames(ref_tab))
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
  sph <- as.matrix(unclass(summary(ref$Anova)$sphericity.tests))
  expect_equal(fit$sphericity$mauchly_w, as.numeric(sph[, 1]))
  expect_equal(fit$sphericity$p_value, as.numeric(sph[, 2]))
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
  d <- withr::with_seed(21, {
    dd <- expand.grid(id = factor(seq_len(2000)),
                      time = factor(c("t1", "t2", "t3")))
    dd$y <- stats::rnorm(nrow(dd))
    dd
  })
  fit <- anova_rm(d, "y", subject = "id", within = "time", plots = FALSE)
  expect_s3_class(fit, "anovakit_fit")
  expect_null(fit$assumptions$normality)
  expect_true(any(grepl("Shapiro-Wilk skipped", fit$notes)))
})

test_that("the afex object is returned intact", {
  skip_if_not_installed("afex")
  d <- fx_repeated()
  fit <- anova_rm(d, "score", subject = "id", within = "time",
                  between = "arm", plots = FALSE)

  expect_s3_class(fit$model, "afex_aov")
  expect_setequal(names(fit$model$data), c("long", "wide", "idata"))
  expect_true(inherits(suppressMessages(emmeans::emmeans(fit$model, ~ time)),
                       "emmGrid"))
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
