# Detailed regression tests: run locally and on CI (NOT_CRAN=true), skipped
# on CRAN to keep its check time short. The core suite runs everywhere.
skip_on_cran()

# Regression tests for the second review of anova_rm() (R2-RM-*).

test_that("R2-RM-01: residuals survive a large constant offset in the response", {
  skip_if_not_installed("afex")
  set.seed(19)
  d <- expand.grid(id = factor(1:12), cond = factor(c("a", "b", "c")))
  d$t <- 1.7e9 + round(0.40 + 0.01 * as.numeric(d$cond) +
                         rnorm(12, 0, 0.05)[d$id] + rnorm(36, 0, 0.004), 3)
  f <- anova_rm(d, "t", "id", "cond", plots = FALSE)
  # Hand computation on the centred data: y - subject mean - cell mean + grand.
  y0 <- d$t - 1.7e9
  e_ref <- y0 - ave(y0, d$id) - ave(y0, d$cond) + mean(y0)
  expect_equal(f$residuals, e_ref, tolerance = 1e-4)
  expect_gt(max(abs(f$residuals)), 1e-3)
  ref_p <- tapply(e_ref, d$cond, function(z) stats::shapiro.test(z)$p.value)
  expect_equal(f$assumptions$normality$p_value, unname(as.numeric(ref_p)),
               tolerance = 1e-3)
  expect_false(any(grepl("Shapiro-Wilk skipped", f$notes)))

  # A genuinely error-free pattern on a large offset is still zeroed.
  d2 <- expand.grid(id = factor(1:5), time = factor(1:3))
  d2$y <- 1e9 + as.numeric(d2$id) + as.numeric(d2$time)
  f2 <- anova_rm(d2, "y", "id", "time", plots = FALSE)
  expect_true(all(f2$residuals == 0))
})

test_that("R2-RM-02: a between level emptied by subject removal is named in $notes", {
  skip_if_not_installed("afex")
  set.seed(16)
  d <- expand.grid(id = factor(1:12), time = factor(c("t1", "t2", "t3")))
  d$arm <- factor(rep(rep(c("a", "b", "c"), each = 4), 3))
  d$y <- rnorm(36)
  d <- d[!(d$arm == "c" & d$time == "t3"), ]
  f <- anova_rm(d, "y", "id", "time", "arm", plots = FALSE)
  # Reference: base aov on the complete subjects, which hold only arms a and b.
  dc <- droplevels(d[d$arm != "c", ])
  ref <- summary(stats::aov(y ~ arm * time + Error(id / time), data = dc))
  ref_df <- ref[["Error: id"]][[1]]["arm", "Df"]
  expect_equal(f$anova$num_df[f$anova$term == "arm"], ref_df)
  expect_true(any(grepl("Level(s) c of between-subject factor `arm`", f$notes,
                        fixed = TRUE)))
  expect_true(any(grepl("compares 2 of its 3 levels", f$notes, fixed = TRUE)))

  # No such note when every level keeps subjects.
  full <- expand.grid(id = factor(1:12), time = factor(c("t1", "t2", "t3")))
  full$arm <- factor(rep(rep(c("a", "b", "c"), each = 4), 3))
  full$y <- rnorm(36)
  h <- anova_rm(full, "y", "id", "time", "arm", plots = FALSE)
  expect_false(any(grepl("of between-subject factor", h$notes)))
})

test_that("R2-RM-03: zero error variance gives p = NA under every correction", {
  skip_if_not_installed("afex")
  d <- expand.grid(id = factor(1:5), time = factor(1:3))
  d$y <- as.numeric(d$id) + as.numeric(d$time)
  # Hand reference: the within-subject error sum of squares is exactly zero.
  e <- d$y - ave(d$y, d$id) - ave(d$y, d$time) + mean(d$y)
  expect_equal(sum(e^2), 0)
  for (cc in c("GG", "HF", "none")) {
    f <- anova_rm(d, "y", "id", "time", correction = cc, plots = FALSE)
    expect_true(is.na(f$anova$p_value[f$anova$term == "time"]), info = cc)
    expect_true(any(grepl("no test is possible", f$notes)), info = cc)
  }
  d2 <- expand.grid(id = factor(1:5), time = factor(1:2))
  d2$y <- as.numeric(d2$id) + as.numeric(d2$time)
  f2 <- anova_rm(d2, "y", "id", "time", plots = FALSE)
  expect_true(is.na(f2$anova$p_value))
  expect_true(any(grepl("no test is possible", f2$notes)))
})

test_that("R2-RM-04: summary prints the sphericity table once", {
  skip_if_not_installed("afex")
  set.seed(1)
  d <- expand.grid(id = factor(1:12), time = factor(c("t1", "t2", "t3")))
  d$y <- rnorm(36)
  f <- anova_rm(d, "y", "id", "time", plots = FALSE)
  out <- utils::capture.output(print(summary(f)))
  expect_length(grep("^ *[Ss]phericity$", out), 1L)
  # Reference for the content: car's Mauchly statistic on the wide data.
  wide <- stats::reshape(d, idvar = "id", timevar = "time", direction = "wide")
  mlm <- stats::lm(cbind(y.t1, y.t2, y.t3) ~ 1, data = wide)
  W <- stats::mauchly.test(mlm, X = ~1)$statistic
  expect_equal(f$sphericity$mauchly_w, unname(W), tolerance = 1e-6)
})

test_that("R2-RM-05: within-only type = 2 does not relay car's Type III note", {
  skip_if_not_installed("afex")
  set.seed(8)
  d <- expand.grid(id = factor(1:6), time = factor(1:3))
  d$y <- rnorm(18)
  f2 <- anova_rm(d, "y", "id", "time", type = 2, plots = FALSE)
  expect_identical(attr(f2$anova, "ss_type"), "II")
  expect_false(any(grepl("type-III tests substituted", f2$notes)))
  # Reference: the uncorrected within-subject F from base aov, identical for
  # Type II and III here.
  ref <- summary(stats::aov(y ~ time + Error(id / time), data = d))
  ref_F <- ref[["Error: id:time"]][[1]]["time", "F value"]
  expect_equal(f2$anova$statistic, ref_F, tolerance = 1e-8)
})

test_that("R2-RM-06: no error df gives correction 'none' with a note saying why", {
  skip_if_not_installed("afex")
  set.seed(12)
  d <- expand.grid(id = factor(1:2), time = factor(1:3))
  d$g <- factor(rep(c("a", "b"), 3))
  d$y <- rnorm(6)
  f <- anova_rm(d, "y", "id", "time", "g", plots = FALSE)
  # Hand reference: error df = subjects - between groups = 2 - 2 = 0.
  expect_equal(nlevels(d$id) - nlevels(d$g), 0L)
  expect_identical(attr(f$anova, "correction"), "none")
  expect_true(any(grepl("Sphericity could not be assessed", f$notes)))
  expect_true(all(is.na(f$anova$p_value)))
})
