# Regression tests for the round-3 review items in the manova area.

test_that("R3-40: Box's M is not computed when a group covariance matrix is singular up to rounding", {
  d <- withr::with_seed(4, data.frame(
    g = factor(rep(c("a", "b", "c"), each = 20)),
    y1 = stats::rnorm(60), y2 = stats::rnorm(60), y3 = stats::rnorm(60)))
  i <- d$g == "c"
  for (w in list(c(0.3, 0.7), c(1, 1), c(0.5, 0.25), c(2, -1))) {
    d$y3[i] <- w[1] * d$y1[i] + w[2] * d$y2[i]
    # Reference: the covariance of group c is singular (eigenvalue ratio at
    # rounding level), so Box's M is undefined there.
    ev <- eigen(stats::cov(d[i, c("y1", "y2", "y3")]), symmetric = TRUE)$values
    expect_lt(min(ev) / max(ev), 1e-12)
    f <- anova_manova(d, c("y1", "y2", "y3"), "g", posthoc = FALSE,
                      plots = FALSE)
    expect_true(is.na(f$assumptions$box_m$statistic))
    expect_true(is.na(f$assumptions$box_m$p_value))
    expect_true(any(grepl(
      "Box's M not computed: at least one covariance matrix is singular",
      f$notes, fixed = TRUE)))
    expect_false(any(grepl("Box's M rejects", f$notes, fixed = TRUE)))
  }

  # The tolerance is relative and scale-free: a response on a scale 1e8 times
  # larger is not mistaken for a singular one, and Box's M (scale invariant)
  # still equals rstatix::box_m(d[, c("score1", "score2")], d$g) on the
  # fixture, 9.91179139028051 on 6 df.
  m <- fx_multivariate()
  m$score2 <- m$score2 * 1e8
  fm <- anova_manova(m, c("score1", "score2"), "g", plots = FALSE)
  expect_equal(fm$assumptions$box_m$statistic, 9.91179139028051,
               tolerance = 1e-8)
  expect_equal(fm$assumptions$box_m$df, 6)
  m$score2 <- m$score2 * 1e-16
  fm <- anova_manova(m, c("score1", "score2"), "g", plots = FALSE)
  expect_equal(fm$assumptions$box_m$statistic, 9.91179139028051,
               tolerance = 1e-8)
})

test_that("R3-41: the Box's M size note speaks of cells and names the small one", {
  d <- expand.grid(A = factor(1:2), B = factor(1:3), r = 1:10)[, 1:2]
  d <- d[!(d$A == 1 & d$B == 3) | seq_len(nrow(d)) %in% c(5, 11), ]
  d <- withr::with_seed(6, {
    d$y1 <- stats::rnorm(nrow(d)); d$y2 <- stats::rnorm(nrow(d)); d
  })
  # Reference: the A x B table. No level of A or B is small; the cell
  # A = 1, B = 3 has 2 rows, no more than the 2 responses.
  tab <- table(d$A, d$B)
  expect_gt(min(table(d$A)), 2L)
  expect_gt(min(table(d$B)), 2L)
  expect_identical(min(tab), 2L)
  expect_identical(which(tab == 2L, arr.ind = TRUE)[1, ], c(row = 1L, col = 3L))

  for (inter in c(FALSE, TRUE)) {
    f <- anova_manova(d, c("y1", "y2"), c("A", "B"), interaction = inter,
                      posthoc = FALSE, plots = FALSE)
    note <- grep("^Box's M not computed", f$notes, value = TRUE)
    expect_identical(note,
      "Box's M not computed: every cell (combination of `A` and `B`) needs more than 2 observations, and cell (`A` = `1`, `B` = `3`) has 2.")
    expect_true(is.na(f$assumptions$box_m$statistic))
  }

  # One grouping variable: the unit is a group, named the same way.
  o <- withr::with_seed(2, data.frame(
    g = factor(rep(c("a", "b", "c"), c(10, 10, 2))),
    y1 = stats::rnorm(22), y2 = stats::rnorm(22)))
  expect_identical(min(table(o$g)), 2L)
  f <- anova_manova(o, c("y1", "y2"), "g", posthoc = FALSE, plots = FALSE)
  expect_identical(grep("^Box's M not computed", f$notes, value = TRUE),
    "Box's M not computed: every group needs more than 2 observations, and group `c` has 2.")
})

test_that("R3-48: anova_manova refuses constant and group-determined covariates", {
  d <- withr::with_seed(8, data.frame(
    g = factor(rep(c("a", "b", "c"), each = 20)),
    y1 = stats::rnorm(60), y2 = stats::rnorm(60)))
  d$k <- 5
  d$x <- as.numeric(d$g)
  d$x2 <- 3 * d$x - 1
  # Reference: lm() cannot estimate any of them next to g.
  expect_true(is.na(stats::coef(stats::lm(y1 ~ k + g, d))[["k"]]))
  expect_true(anyNA(stats::coef(stats::lm(y1 ~ x + g, d))))
  expect_identical(stats::sd(d$k), 0)

  expect_error(anova_manova(d, c("y1", "y2"), "g", covariates = "k",
                            plots = FALSE),
               "Covariate `k` has the same value \\(5\\) in every analysed row")
  expect_error(anova_manova(d, c("y1", "y2"), "g", covariates = "x",
                            plots = FALSE),
               "`x` are determined by the grouping variable")
  expect_error(anova_manova(d, c("y1", "y2"), "g", covariates = "x2",
                            type = "III", plots = FALSE),
               "`x2` are determined by the grouping variable")

  # Determined by one factor of an additive two-factor design.
  d$B <- factor(rep(c("u", "v"), 30))
  expect_error(anova_manova(d, c("y1", "y2"), c("g", "B"), covariates = "x",
                            plots = FALSE),
               "`x` are determined by the grouping variable")

  # Valid neighbours still run and match stats::manova (sequential with g
  # last equals Type II for g): a covariate with a huge offset and a tiny
  # spread, and one constant in a single group only.
  d$z <- 1e6 + withr::with_seed(9, stats::rnorm(60, sd = 1e-3))
  f <- anova_manova(d, c("y1", "y2"), "g", covariates = "z", plots = FALSE)
  zc <- d$z - mean(d$z)
  ref <- summary(stats::manova(cbind(y1, y2) ~ zc + g, data = d))
  expect_equal(f$multivariate$df[f$multivariate$term == "g"], 2)
  expect_equal(f$multivariate$p_value[f$multivariate$term == "g"],
               as.data.frame(ref$stats)[["g", "Pr(>F)"]])
  d$w <- withr::with_seed(10, stats::rnorm(60))
  d$w[d$g == "a"] <- 0
  f <- anova_manova(d, c("y1", "y2"), "g", covariates = "w", plots = FALSE)
  ref <- summary(stats::manova(cbind(y1, y2) ~ w + g, data = d))
  expect_equal(f$multivariate$p_value[f$multivariate$term == "g"],
               as.data.frame(ref$stats)[["g", "Pr(>F)"]])
})

test_that("R3-60: the tests guard package behaviour, not R's argument matching", {
  # The package's own Wilks statistic, not only a product of its
  # eigenvalues, against stats::manova.
  d <- fx_multivariate()
  fit <- anova_manova(d, c("score1", "score2"), "g", test = "Wilks",
                      plots = FALSE)
  ref <- as.data.frame(summary(stats::manova(cbind(score1, score2) ~ g, data = d),
                               test = "Wilks")$stats)
  expect_equal(fit$multivariate$statistic, ref[["g", "Wilks"]])
  expect_equal(fit$multivariate$approx_f, ref[["g", "approx F"]])
  expect_equal(fit$multivariate$p_value, ref[["g", "Pr(>F)"]])

  # anova_bin() and anova_glm() have no `...`, and the model they fit is the
  # one their documented arguments describe: binomial-logit for anova_bin(),
  # the requested family for anova_glm().
  expect_false("..." %in% names(formals(anova_bin)))
  expect_false("..." %in% names(formals(anova_glm)))
  b <- fx_binary()
  fb <- anova_bin(b, "y", "g", posthoc = FALSE, plots = FALSE)
  rb <- stats::glm(y ~ g, family = stats::binomial(), data = b)
  expect_identical(stats::family(fb$model)$family, "binomial")
  expect_identical(stats::family(fb$model)$link, "logit")
  expect_equal(unname(stats::fitted(fb$model)), unname(stats::fitted(rb)))
  o <- fx_oneway()
  fg <- anova_glm(o, "value", "group", family = stats::gaussian(link = "log"),
                  posthoc = FALSE, plots = FALSE)
  rg <- stats::glm(value ~ group, family = stats::gaussian(link = "log"),
                   data = o)
  expect_identical(stats::family(fg$model)$link, "log")
  expect_equal(unname(stats::fitted(fg$model)), unname(stats::fitted(rg)),
               tolerance = 1e-6)
})

test_that("manova-C1/R3/R4: the Box's M size note keeps cells, names and tied cells apart", {
  d <- expand.grid(A = factor(1:2), B = factor(1:3), r = 1:10)[, 1:2]
  d <- d[!(d$A == 1 & d$B == 3) | seq_len(nrow(d)) %in% c(5, 11), ]
  d <- withr::with_seed(6, {
    d$y1 <- stats::rnorm(nrow(d)); d$y2 <- stats::rnorm(nrow(d)); d
  })
  # Two cells tie for the smallest: the reference is the table of complete rows.
  d6 <- d
  d6$y1[d6$A == 2 & d6$B == 2][1:8] <- NA
  cc <- d6[stats::complete.cases(d6), ]
  tab <- table(cc$A, cc$B)
  expect_identical(min(tab), 2L)
  expect_identical(sum(tab == 2L), 2L)
  expect_identical(unname(tab["2", "2"]), 2L)
  expect_identical(unname(tab["1", "3"]), 2L)
  f <- anova_manova(d6, c("y1", "y2"), c("A", "B"), posthoc = FALSE, plots = FALSE)
  expect_identical(grep("^Box's M not computed", f$notes, value = TRUE),
    "Box's M not computed: every cell (combination of `A` and `B`) needs more than 2 observations, and 2 cells have 2: (`A` = `2`, `B` = `2`); (`A` = `1`, `B` = `3`).")

  # A comma inside a level and a space in a variable name stay inside backticks.
  d5 <- d
  levels(d5$A) <- c("lo, a", "hi")
  names(d5)[1] <- "my A"
  f <- anova_manova(d5, c("y1", "y2"), c("my A", "B"), posthoc = FALSE, plots = FALSE)
  expect_identical(grep("^Box's M not computed", f$notes, value = TRUE),
    "Box's M not computed: every cell (combination of `my A` and `B`) needs more than 2 observations, and cell (`my A` = `lo, a`, `B` = `3`) has 2.")

  # A comma inside a grouping-variable name is not rewritten to "and".
  m <- data.frame(cyl = factor(mtcars$cyl), `dose, mg` = factor(mtcars$am),
                  mpg = mtcars$mpg, hp = mtcars$hp, wt = mtcars$wt,
                  check.names = FALSE)
  tab <- table(m$cyl, m[["dose, mg"]])
  expect_identical(min(tab), 2L)
  expect_identical(sum(tab == 2L), 1L)
  expect_identical(unname(tab["8", "1"]), 2L)
  f <- anova_manova(m, c("mpg", "hp", "wt"), c("cyl", "dose, mg"), plots = FALSE)
  expect_identical(grep("^Box's M not computed", f$notes, value = TRUE),
    "Box's M not computed: every cell (combination of `cyl` and `dose, mg`) needs more than 3 observations, and cell (`cyl` = `8`, `dose, mg` = `1`) has 2.")

  # Three grouping variables are listed as "`a`, `b` and `c`".
  m3 <- data.frame(cyl = factor(mtcars$cyl), am = factor(mtcars$am),
                   vs = factor(mtcars$vs), mpg = mtcars$mpg, hp = mtcars$hp)
  f <- anova_manova(m3, c("mpg", "hp"), c("cyl", "am", "vs"),
                    interaction = FALSE, posthoc = FALSE, plots = FALSE)
  expect_match(grep("^Box's M not computed", f$notes, value = TRUE),
               "every cell (combination of `cyl`, `am` and `vs`)", fixed = TRUE)

  # More than three tied groups: three are named and the rest counted.
  o <- withr::with_seed(1, data.frame(
    g = factor(rep(letters[1:6], c(10, 10, 2, 2, 2, 2))),
    y1 = stats::rnorm(28), y2 = stats::rnorm(28)))
  expect_identical(names(which(table(o$g) == 2L)), c("c", "d", "e", "f"))
  f <- anova_manova(o, c("y1", "y2"), "g", plots = FALSE)
  expect_identical(grep("^Box's M not computed", f$notes, value = TRUE),
    "Box's M not computed: every group needs more than 2 observations, and 4 groups have 2: `c`; `d`; `e`; and 1 more.")

  # Three tied cells of two factors, separated by "; ".
  mt <- mtcars
  mt$gear <- factor(mt$gear)
  mt$carb <- factor(mt$carb)
  tab <- table(mt$gear, mt$carb)
  tab[tab == 0L] <- NA
  expect_identical(min(tab, na.rm = TRUE), 1L)
  expect_identical(unname(which(tab == 1L, arr.ind = TRUE)[, "row"]), c(3L, 3L, 3L))
  expect_identical(colnames(tab)[which(tab == 1L, arr.ind = TRUE)[, "col"]],
                   c("4", "6", "8"))
  f <- anova_manova(mt, c("mpg", "hp"), c("gear", "carb"), interaction = FALSE,
                    plots = FALSE)
  expect_identical(grep("^Box's M not computed", f$notes, value = TRUE),
    "Box's M not computed: every cell (combination of `gear` and `carb`) needs more than 2 observations, and 3 cells have 1: (`gear` = `5`, `carb` = `4`); (`gear` = `5`, `carb` = `6`); (`gear` = `5`, `carb` = `8`).")
})

test_that("manova-R5: a covariate constant up to rounding is refused, as stats::manova without it is the right fit", {
  d <- withr::with_seed(3, data.frame(
    g = factor(rep(c("a", "b", "c"), each = 15)),
    y1 = stats::rnorm(45), y2 = stats::rnorm(45)))
  d$x <- ifelse(seq_len(45) %% 2 == 0, 0.3, 0.1 + 0.2)
  # Two distinct doubles that differ only by rounding.
  expect_identical(length(unique(d$x)), 2L)
  expect_lt(diff(range(d$x)), 1e-15)
  msg <- "Covariate `x` has the same value \\(0.3\\) in every analysed row up to rounding error"
  expect_error(anova_manova(d, c("y1", "y2"), "g", covariates = "x", plots = FALSE),
               msg)
  expect_error(anova_ancova(d, "y1", "g", covariates = "x", plots = FALSE), msg)
  # Without the covariate the group test matches stats::manova (den Df 84).
  ref <- summary(stats::manova(cbind(y1, y2) ~ g, d))$stats
  f <- anova_manova(d, c("y1", "y2"), "g", plots = FALSE)
  mv <- f$multivariate
  row <- mv[mv$term == "g", ]
  expect_equal(row$den_df, ref["g", "den Df"])
  expect_equal(row$p_value, ref["g", "Pr(>F)"])
  # A covariate with a small but real spread is still accepted.
  d$x2 <- 1e6 + seq_len(45) * 1e-3
  expect_silent(anova_manova(d, c("y1", "y2"), "g", covariates = "x2", plots = FALSE))
})
