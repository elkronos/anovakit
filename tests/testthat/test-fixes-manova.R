# Regression tests for the anova_manova() findings of the adversarial review.
#
# Every expectation is checked against an independent reference: stats
# (summary.manova, anova.mlm, SSD), car, emmeans, MASS::lda, rstatix, or
# numpy/statsmodels values computed outside R and pasted in.

# --- Local fixtures -----------------------------------------------------------

mv_unbalanced <- function() {
  withr::with_seed(2, {
    d <- data.frame(A = factor(sample(c("a1", "a2"), 80, TRUE, prob = c(.8, .2))))
    d$B <- factor(ifelse(d$A == "a1",
                         sample(c("b1", "b2"), 80, TRUE, prob = c(.8, .2)),
                         sample(c("b1", "b2"), 80, TRUE, prob = c(.2, .8))))
    d$y1 <- 1.0 * (d$B == "b2") + stats::rnorm(80)
    d$y2 <- 0.5 * (d$B == "b2") + stats::rnorm(80)
    d$age <- stats::rnorm(80, 40, 8) + 3 * (d$B == "b2")
    d$y1 <- d$y1 + 0.05 * d$age
    d
  })
}

mv_balanced_2way <- function() {
  withr::with_seed(17, {
    dd <- data.frame(A = factor(rep(c("a1", "a2", "a3"), each = 60)),
                     B = factor(rep(rep(c("b1", "b2"), each = 30), 3)))
    dd$y1 <- 1 * (dd$A == "a2") + 4 * (dd$B == "b2") + stats::rnorm(180)
    dd$y2 <- 1 * (dd$A == "a3") + 4 * (dd$B == "b2") + stats::rnorm(180)
    dd
  })
}

mv_fx <- function() {
  withr::with_seed(110, {
    n <- 120
    d <- data.frame(g = factor(rep(c("a", "b", "c"), each = n / 3)),
                    age = stats::rnorm(n, 40, 8))
    d$score1 <- 5 + 2 * (d$g == "b") + 4 * (d$g == "c") + 0.1 * d$age +
      stats::rnorm(n)
    d$score2 <- 3 + 1 * (d$g == "b") + 2 * (d$g == "c") + 0.05 * d$age +
      stats::rnorm(n)
    d
  })
}

# Between-group over within-group sum of squares of a score.
bw_ratio <- function(s, g) {
  sum(table(g) * (tapply(s, g, mean) - mean(s))^2) / sum((s - stats::ave(s, g))^2)
}

stats_row <- function(sm, term, test = "Pillai") {
  st <- as.data.frame(sm$stats)
  c(statistic = st[[term, test]], p = st[[term, "Pr(>F)"]])
}

# --- F024: Type II / III multivariate tests -----------------------------------

test_that("F024: Type II tests do not depend on the order of groups and match the term-last sequential test", {
  d <- mv_unbalanced()
  ab <- anova_manova(d, c("y1", "y2"), c("A", "B"), plots = FALSE,
                     assumptions = FALSE)$multivariate
  ba <- anova_manova(d, c("y1", "y2"), c("B", "A"), plots = FALSE,
                     assumptions = FALSE)$multivariate
  expect_equal(ab[match(c("A", "B"), ab$term), -1L],
               ba[match(c("A", "B"), ba$term), -1L], ignore_attr = TRUE)

  # A Type II test of a term in an additive model is the sequential test with
  # that term entered last.
  a_last <- stats_row(summary(stats::manova(cbind(y1, y2) ~ B + A, data = d)), "A")
  b_last <- stats_row(summary(stats::manova(cbind(y1, y2) ~ A + B, data = d)), "B")
  expect_equal(ab$statistic[ab$term == "A"], a_last[["statistic"]])
  expect_equal(ab$p_value[ab$term == "A"], a_last[["p"]])
  expect_equal(ab$statistic[ab$term == "B"], b_last[["statistic"]])
  expect_equal(ab$p_value[ab$term == "B"], b_last[["p"]])
  expect_equal(ab$p_value[ab$term == "A"], 0.6267800, tolerance = 1e-6)
})

test_that("F024: the MANCOVA covariate row is adjusted for the groups", {
  d <- mv_unbalanced()
  fit <- anova_manova(d, c("y1", "y2"), "B", covariates = "age", plots = FALSE,
                      assumptions = FALSE)
  mv <- fit$multivariate
  age_last <- stats_row(summary(stats::manova(cbind(y1, y2) ~ B + age, data = d)), "age")
  b_last <- stats_row(summary(stats::manova(cbind(y1, y2) ~ age + B, data = d)), "B")
  expect_equal(mv$statistic[mv$term == "age"], age_last[["statistic"]])
  expect_equal(mv$p_value[mv$term == "age"], age_last[["p"]])
  expect_equal(mv$p_value[mv$term == "B"], b_last[["p"]])
  expect_equal(mv$p_value[mv$term == "age"], 0.02037211, tolerance = 1e-6)
})

test_that("F024: Type III tests match the model comparison that drops the term's columns", {
  d <- mv_unbalanced()
  fit <- anova_manova(d, c("y1", "y2"), c("A", "B"), interaction = TRUE,
                      type = "III", plots = FALSE, assumptions = FALSE)
  mv <- fit$multivariate
  expect_match(fit$method, "Type III")
  expect_identical(attr(mv, "ss_type"), "III")

  withr::local_options(contrasts = c("contr.sum", "contr.poly"))
  Y <- cbind(d$y1, d$y2)
  full <- stats::lm(Y ~ A * B, data = d)
  X <- stats::model.matrix(full)
  asg <- attr(X, "assign")
  labs <- attr(stats::terms(full), "term.labels")
  for (test in c("Pillai", "Wilks", "Hotelling-Lawley", "Roy")) {
    f <- anova_manova(d, c("y1", "y2"), c("A", "B"), interaction = TRUE,
                      type = "III", test = test, plots = FALSE,
                      assumptions = FALSE)$multivariate
    for (j in seq_along(labs)) {
      reduced <- stats::lm(Y ~ X[, asg != j] - 1)
      ref <- as.data.frame(stats::anova(full, reduced, test = test))[2L, ]
      expect_equal(f$statistic[f$term == labs[j]], ref[[test]],
                   info = paste(test, labs[j]))
      expect_equal(f$approx_f[f$term == labs[j]], ref[["approx F"]],
                   info = paste(test, labs[j]))
      expect_equal(f$p_value[f$term == labs[j]], ref[["Pr(>F)"]],
                   info = paste(test, labs[j]))
    }
  }
  expect_equal(mv$p_value[mv$term == "A"], 0.5871207, tolerance = 1e-6)
})

test_that("F024: Type II with an interaction uses the full model's error matrix", {
  d <- mv_unbalanced()
  mv <- anova_manova(d, c("y1", "y2"), c("A", "B"), interaction = TRUE,
                     plots = FALSE, assumptions = FALSE)$multivariate
  Y <- cbind(d$y1, d$y2)
  # anova.mlm tests every comparison against the largest model's residuals,
  # which is the Type II definition.
  ref <- as.data.frame(stats::anova(stats::lm(Y ~ A * B, d), stats::lm(Y ~ A + B, d),
                                    stats::lm(Y ~ B, d), test = "Pillai"))
  expect_equal(mv$statistic[mv$term == "A"], ref$Pillai[3])
  expect_equal(mv$p_value[mv$term == "A"], ref[["Pr(>F)"]][3])
  expect_equal(mv$p_value[mv$term == "A:B"], ref[["Pr(>F)"]][2])
})

test_that("F024: the canonical analysis uses the same (adjusted) hypothesis matrix as the test", {
  d <- mv_unbalanced()
  fit <- anova_manova(d, c("y1", "y2"), c("A", "B"), plots = FALSE)
  expect_identical(fit$canonical_term, "A")
  E <- stats::SSD(stats::lm(cbind(y1, y2) ~ A + B, data = d))$SSD
  H <- stats::SSD(stats::lm(cbind(y1, y2) ~ B, data = d))$SSD - E
  ev <- sort(Re(eigen(solve(E) %*% H)$values), decreasing = TRUE)
  expect_equal(fit$canonical$eigenvalue, ev[seq_len(nrow(fit$canonical))])
  expect_equal(fit$canonical$eigenvalue[1], 0.01236955, tolerance = 1e-5)
})

# --- F026 / F027: the marginal-means plot -------------------------------------

test_that("F026: a cell grid maps the other factors to colour and each panel is one response", {
  dd <- withr::with_seed(17, {
    dd <- data.frame(A = factor(rep(c("a1", "a2"), each = 90)),
                     B = factor(rep(rep(c("b1", "b2", "b3"), each = 30), 2)))
    dd$y1 <- 2 * (dd$A == "a2") + 3 * (dd$B == "b3") + stats::rnorm(180)
    dd$y2 <- 3 * (dd$B == "b3") + stats::rnorm(180)
    dd
  })
  p <- anova_manova(dd, c("y1", "y2"), c("A", "B"), interaction = TRUE,
                    assumptions = FALSE)$plots$emmeans
  pts <- ggplot2::layer_data(p, 2L)          # the points
  expect_identical(length(unique(pts$PANEL)), 2L)
  expect_identical(length(unique(pts$group)), 3L)   # one series per level of B
  # Points of different series are dodged apart, not stacked at one x
  expect_identical(anyDuplicated(pts[, c("PANEL", "x")]), 0L)
  ref <- as.data.frame(emmeans::emmeans(stats::lm(y1 ~ A * B, data = dd), ~ A * B))
  expect_equal(sort(pts$y[pts$PANEL == 1]), sort(ref$emmean))

  # Additive: one panel per response and factor, each factor's own means
  pa <- anova_manova(dd, c("y1", "y2"), c("A", "B"), assumptions = FALSE)$plots$emmeans
  lay <- ggplot2::ggplot_build(pa)$layout$layout
  expect_identical(nrow(lay), 4L)
  ref_a <- as.data.frame(emmeans::emmeans(stats::lm(y2 ~ A + B, data = dd), ~ A))
  pts <- ggplot2::layer_data(pa, 1L)
  panel <- lay$PANEL[lay[[".response"]] == "y2" & lay[[".term"]] == "A"]
  expect_equal(pts$y[pts$PANEL == panel], ref_a$emmean)
})

test_that("F027: the plot facets on the response even when a group is called response", {
  d <- withr::with_seed(110, {
    d <- data.frame(response = factor(rep(c("a", "b", "c"), each = 40)))
    d$score1 <- 5 + 2 * (d$response == "b") + 4 * (d$response == "c") + stats::rnorm(120)
    d$score2 <- 3 + 1 * (d$response == "b") + 2 * (d$response == "c") + stats::rnorm(120)
    d
  })
  fit <- anova_manova(d, c("score1", "score2"), "response", assumptions = FALSE)
  pts <- ggplot2::layer_data(fit$plots$emmeans, 2L)
  expect_identical(length(unique(pts$PANEL)), 2L)
  expect_identical(as.integer(table(pts$PANEL)), c(3L, 3L))
  ref <- as.data.frame(emmeans::emmeans(stats::lm(score2 ~ response, data = d), ~ response))
  expect_equal(pts$y[pts$PANEL == 2], ref$emmean)
})

# --- F070 / F071 / F072: canonical scores and adjustment -----------------------

test_that("F070: canonical scores and structure coefficients are adjusted for the other factors", {
  dd <- mv_balanced_2way()
  fit <- anova_manova(dd, c("y1", "y2"), c("A", "B"))
  expect_identical(fit$canonical_term, "A")
  sc <- fit$plots$canonical$data
  dfe <- stats::df.residual(fit$model)
  for (k in seq_len(nrow(fit$canonical))) {
    s <- sc[[paste0("Can", k)]]
    expect_equal(bw_ratio(s, dd$A), fit$canonical$eigenvalue[k], tolerance = 1e-8)
    expect_equal(sum((s - stats::ave(s, dd$A))^2) / dfe, 1, tolerance = 1e-8)
  }
  # Structure coefficients: diag(E)^-1/2 E v / sqrt(v'Ev), with E and H from stats
  E <- stats::SSD(stats::lm(cbind(y1, y2) ~ A + B, data = dd))$SSD
  H <- stats::SSD(stats::lm(cbind(y1, y2) ~ B, data = dd))$SSD - E
  V <- Re(eigen(solve(E) %*% H)$vectors)
  ref <- diag(1 / sqrt(diag(E))) %*% E %*% V %*% diag(1 / sqrt(diag(t(V) %*% E %*% V)))
  got <- as.matrix(fit$assumptions$structure_coefficients[, c("Can1", "Can2")])
  expect_equal(abs(unname(got)), abs(unname(ref)), tolerance = 1e-8)
  expect_equal(abs(got[, 1]), c(0.6309888, 0.8116607), tolerance = 1e-6)

  # With the interaction, the axes describe A:B and the plot shows A:B alone
  fi <- anova_manova(dd, c("y1", "y2"), c("A", "B"), interaction = TRUE)
  expect_identical(fi$canonical_term, "A:B")
  cell <- fi$plots$canonical$data$.group
  for (k in seq_len(nrow(fi$canonical))) {
    expect_equal(bw_ratio(fi$plots$canonical$data[[paste0("Can", k)]], cell),
                 fi$canonical$eigenvalue[k], tolerance = 1e-8)
  }
})

test_that("F071: mixed syntactic and non-syntactic names give the same adjusted analysis", {
  d <- withr::with_seed(9, {
    n <- 180
    g <- factor(rep(c("a", "b", "c"), each = n / 3))
    age <- stats::rnorm(n, 40, 15); bmi <- stats::rnorm(n, 25, 4)
    d <- data.frame(g = g, age = age, `body mass` = bmi, check.names = FALSE)
    d$score1 <- 2 * age + 0.5 * bmi + 3 * (g == "b") + stats::rnorm(n)
    d$score2 <- 1 * age - 0.5 * bmi - 3 * (g == "b") + stats::rnorm(n)
    d
  })
  mixed <- anova_manova(d, c("score1", "score2"), "g",
                        covariates = c("age", "body mass"))
  d2 <- d; names(d2)[3] <- "bodymass"
  plain <- anova_manova(d2, c("score1", "score2"), "g",
                        covariates = c("age", "bodymass"))
  s <- mixed$plots$canonical$data$Can1
  expect_equal(bw_ratio(s, d$g), mixed$canonical$eigenvalue[1], tolerance = 0.01)
  expect_equal(s, plain$plots$canonical$data$Can1)
  expect_equal(mixed$assumptions$box_m, plain$assumptions$box_m)
  expect_false(any(grepl("could not be", mixed$notes)))
  # Type II covariate-adjusted g test against the term-last sequential test
  ref <- stats_row(summary(stats::manova(cbind(score1, score2) ~ age + `body mass` + g,
                                          data = d)), "g")
  expect_equal(mixed$multivariate$p_value[mixed$multivariate$term == "g"], ref[["p"]])

  dd <- withr::with_seed(17, {
    dd <- data.frame(`treatment arm` = factor(rep(c("a1", "a2"), each = 90)),
                     B = factor(rep(rep(c("b1", "b2", "b3"), each = 30), 2)),
                     check.names = FALSE)
    dd$y1 <- 2 * (dd[["treatment arm"]] == "a2") + stats::rnorm(180)
    dd$y2 <- 3 * (dd$B == "b3") + stats::rnorm(180)
    dd
  })
  f <- anova_manova(dd, c("y1", "y2"), c("treatment arm", "B"),
                    interaction = TRUE, plots = FALSE)
  expect_identical(f$canonical_term, "treatment arm:B")
  expect_true(f$canonical_term %in% f$multivariate$term)
})

test_that("F072: collinear covariates leave Box's M, the structure coefficients and the plot intact", {
  d <- mv_fx()
  d$months <- d$age * 12
  r <- anova_manova(d, c("score1", "score2"), "g", covariates = c("age", "months"))
  one <- anova_manova(d, c("score1", "score2"), "g", covariates = "age")
  expect_false(is.null(r$assumptions$box_m))
  expect_equal(r$assumptions$box_m, one$assumptions$box_m)
  expect_true(all(is.finite(as.matrix(r$assumptions$structure_coefficients[, -1L]))))
  expect_true(all(is.finite(as.matrix(r$plots$canonical$data[, c("Can1", "Can2")]))))
  # The group test is unaffected by the aliasing
  ref <- stats_row(summary(stats::manova(cbind(score1, score2) ~ age + months + g,
                                          data = d)), "g")
  expect_equal(r$multivariate$p_value[r$multivariate$term == "g"], ref[["p"]])
  # The aliasing note does not call a covariate an empty cell
  expect_true(any(grepl("aliased", r$notes)))
  expect_false(any(grepl("cells are empty", r$notes)))
})

test_that("F072: an aliased Type III model falls back to Type II throughout, with a note", {
  d <- mv_fx()
  d$site <- factor(rep(c("n", "s"), length.out = nrow(d)))
  e <- d[!(d$g == "c" & d$site == "s"), ]
  fit <- anova_manova(e, c("score1", "score2"), c("g", "site"),
                      interaction = TRUE, type = "III", plots = FALSE)
  expect_match(fit$method, "Type II\\)")
  expect_true(any(grepl("Type III tests are not defined", fit$notes)))
  Y <- cbind(e$score1, e$score2)
  ref <- as.data.frame(stats::anova(stats::lm(Y ~ g * site, e), stats::lm(Y ~ g + site, e),
                                    stats::lm(Y ~ site, e), test = "Pillai"))
  expect_equal(fit$multivariate$p_value[fit$multivariate$term == "g"], ref[["Pr(>F)"]][3])
  expect_false(is.null(fit$univariate$score1$anova))
})

# --- F073 / F139: reserved names ------------------------------------------------

test_that("F073: a column called .Y is an ordinary covariate or group", {
  d <- mv_fx()
  d$.Y <- d$age
  printed <- utils::capture.output(expect_no_warning(
    a <- anova_manova(d, c("score1", "score2"), "g", covariates = ".Y",
                      plots = FALSE)))
  expect_identical(printed, character(0))
  b <- anova_manova(d, c("score1", "score2"), "g", covariates = "age", plots = FALSE)
  expect_equal(a$multivariate[, -1L], b$multivariate[, -1L], ignore_attr = TRUE)
  expect_identical(a$multivariate$term, c(".Y", "g"))

  d2 <- data.frame(.Y = d$g, score1 = d$score1, score2 = d$score2)
  f <- anova_manova(d2, c("score1", "score2"), ".Y", plots = FALSE)
  ref <- stats_row(summary(stats::manova(cbind(score1, score2) ~ .Y, data = d2)), ".Y")
  expect_equal(f$multivariate$p_value, ref[["p"]])
})

test_that("F139: a column named Residuals is refused up front", {
  d <- mv_fx()
  d$Residuals <- d$age
  expect_error(anova_manova(d, c("score1", "score2"), "g", covariates = "Residuals"),
               "reserved")
  names(d)[names(d) == "g"] <- "Residuals"
  d$Residuals.1 <- NULL
  expect_error(anova_manova(d[, c("Residuals", "score1", "score2")],
                            c("score1", "score2"), "Residuals"), "reserved")
})

# --- F074: Roy's bound ------------------------------------------------------------

test_that("F074: Roy's largest root is flagged as a lower bound when df > 1", {
  d <- withr::with_seed(1, data.frame(g = factor(rep(1:4, each = 10)),
                                      y1 = stats::rnorm(40), y2 = stats::rnorm(40),
                                      y3 = stats::rnorm(40)))
  f <- anova_manova(d, c("y1", "y2", "y3"), "g", test = "Roy", plots = FALSE)
  expect_true(any(grepl("Roy's largest root .*lower bound", f$notes)))
  ref <- stats_row(summary(stats::manova(cbind(y1, y2, y3) ~ g, data = d),
                           test = "Roy"), "g", "Roy")
  expect_equal(f$multivariate$statistic, ref[["statistic"]])
  expect_equal(f$multivariate$p_value, ref[["p"]])
  # With one hypothesis degree of freedom the F is exact and there is no note
  d2 <- d[d$g %in% c(1, 2), ]
  f2 <- anova_manova(d2, c("y1", "y2", "y3"), "g", test = "Roy", plots = FALSE)
  expect_false(any(grepl("Roy's largest root", f2$notes)))
  f3 <- anova_manova(d, c("y1", "y2", "y3"), "g", test = "Pillai", plots = FALSE)
  expect_false(any(grepl("Roy's largest root", f3$notes)))
})

# --- F137 / F138: MANCOVA ---------------------------------------------------------

test_that("F137: heterogeneous regression slopes are tested and reported", {
  d2 <- withr::with_seed(4, {
    d2 <- data.frame(g = factor(rep(c("a", "b", "c"), each = 30)), x = stats::rnorm(90))
    sl <- c(a = -2, b = 0, c = 2)[as.character(d2$g)]
    d2$y1 <- sl * d2$x + stats::rnorm(90, sd = .5)
    d2$y2 <- -sl * d2$x + stats::rnorm(90, sd = .5)
    d2
  })
  f <- anova_manova(d2, c("y1", "y2"), "g", covariates = "x", plots = FALSE)
  st <- f$slopes_test
  ref <- as.data.frame(stats::anova(stats::lm(cbind(y1, y2) ~ x + g, d2),
                                    stats::lm(cbind(y1, y2) ~ x * g, d2),
                                    test = "Pillai"))[2L, ]
  expect_equal(st$statistic, ref$Pillai)
  expect_equal(st$p_value, ref[["Pr(>F)"]])
  # statsmodels MANOVA, term x:C(g) of y1 + y2 ~ x * C(g)
  expect_equal(st$statistic, 0.939520877603, tolerance = 1e-9)
  expect_equal(st$approx_f, 37.2094801547, tolerance = 1e-9)
  expect_equal(c(st$num_df, st$den_df), c(4, 168))
  expect_false(st$homogeneous)
  expect_true(any(grepl("Homogeneity of regression slopes is rejected", f$notes)))

  # Homogeneous slopes: the test is reported and nothing fires
  ok <- anova_manova(mv_fx(), c("score1", "score2"), "g", covariates = "age",
                     plots = FALSE)
  expect_true(ok$slopes_test$homogeneous)
  expect_false(any(grepl("slopes is rejected", ok$notes)))
  expect_null(anova_manova(mv_fx(), c("score1", "score2"), "g",
                           plots = FALSE)$slopes_test)

  # summary() prints it
  out <- utils::capture.output(print(summary(f)))
  expect_true(any(grepl("Homogeneity of slopes", out)))
})

test_that("F138: a covariate with a large offset is kept and adjusted for", {
  d <- withr::with_seed(1, {
    n <- 90
    d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 30)),
                    ts = 1.7e9 + stats::rnorm(n, 0, 150))
    z <- (d$ts - 1.7e9) / 150
    d$y1 <- 2 * z + (d$g == "b") + stats::rnorm(n)
    d$y2 <- -z + (d$g == "c") + stats::rnorm(n)
    d
  })
  f <- anova_manova(d, c("y1", "y2"), "g", covariates = "ts", plots = FALSE)
  dc <- d
  dc$ts <- dc$ts - mean(dc$ts)
  ts_last <- stats_row(summary(stats::manova(cbind(y1, y2) ~ g + ts, data = dc)), "ts")
  g_last <- stats_row(summary(stats::manova(cbind(y1, y2) ~ ts + g, data = dc)), "g")
  mv <- f$multivariate
  expect_identical(mv$term, c("ts", "g"))
  expect_equal(mv$statistic[1], ts_last[["statistic"]])
  expect_equal(mv$p_value[2], g_last[["p"]])
  expect_equal(mv$p_value[2], 3.046035e-05, tolerance = 1e-5)
  expect_false(is.null(f$assumptions$box_m))
  expect_equal(unname(f$covariate_means), mean(d$ts))
  expect_true(any(grepl("mean-centred", f$notes)))
  expect_false(any(grepl("aliased", f$notes)))
})

# --- F166: which term the canonical analysis describes -------------------------

test_that("F166: the note says which term the canonical analysis describes and why", {
  d <- withr::with_seed(5, {
    d <- expand.grid(A = paste0("a", 1:3), B = paste0("b", 1:2),
                     C = paste0("c", 1:2), rep = 1:5)
    d$y1 <- 2 * ((d$B == "b2") == (d$C == "c2")) + stats::rnorm(nrow(d))
    d$y2 <- stats::rnorm(nrow(d))
    d
  })
  f <- anova_manova(d, c("y1", "y2"), c("A", "B", "C"), interaction = 2,
                    assumptions = FALSE)
  expect_identical(f$canonical_term, "A")
  expect_true(any(grepl("describes the `A` term only.*first grouping variable", f$notes)))
  f2 <- anova_manova(d, c("y1", "y2"), c("C", "B", "A"), interaction = 2,
                     assumptions = FALSE, plots = FALSE)
  expect_identical(f2$canonical_term, "C")
  f3 <- anova_manova(d, c("y1", "y2"), c("A", "B", "C"), interaction = TRUE,
                     assumptions = FALSE, plots = FALSE)
  expect_identical(f3$canonical_term, "A:B:C")
  expect_true(any(grepl("full interaction", f3$notes)))
  # The hypothesis matrix is the canonical term's own (Type II) one
  expect_identical(nlevels(droplevels(f$plots$canonical$data$.group)), 3L)
})

# --- F075 / F105 ------------------------------------------------------------------

test_that("F075: unprotected univariate follow-ups are flagged", {
  d <- withr::with_seed(3, {
    d <- data.frame(g = factor(rep(1:3, each = 20)),
                    matrix(stats::rnorm(60 * 6), 60, 6))
    d$X1 <- d$X1 + 0.8 * (d$g == 2)
    d
  })
  fit <- anova_manova(d, paste0("X", 1:6), "g", plots = FALSE, assumptions = FALSE)
  expect_gt(fit$multivariate$p_value, 0.05)
  expect_true(any(grepl("not protected by the multivariate test", fit$notes)))
  sig <- anova_manova(mv_fx(), c("score1", "score2"), "g", plots = FALSE)
  expect_false(any(grepl("not protected", sig$notes)))
})

test_that("F105: stacked tables have default row names", {
  fit <- anova_manova(mv_fx(), c("score1", "score2"), "g", plots = FALSE)
  for (nm in c("emmeans", "posthoc", "effect_sizes")) {
    expect_identical(rownames(fit[[nm]]), as.character(seq_len(nrow(fit[[nm]]))),
                     info = nm)
  }
})

# --- F128 / F131 / F025 / F108: Mardia and scale ----------------------------------

test_that("F108: Mardia's tests reproduce an independent computation on the residuals", {
  fit <- anova_manova(mv_fx(), c("score1", "score2"), "g", plots = FALSE)
  m <- fit$assumptions$mardia
  # numpy, from Mardia (1970), on residuals(lm(cbind(score1, score2) ~ g))
  expect_equal(m$statistic, c(1.318628772380, -1.207496117805), tolerance = 1e-9)
  expect_equal(m$p_value, c(0.858207051789, 0.227241139256), tolerance = 1e-9)
  expect_equal(m$df[1], 4)
  expect_false(any(grepl("reject multivariate normality", fit$notes)))
})

test_that("F128: Mardia's tests run in linear memory on a large sample", {
  Y <- withr::with_seed(5, matrix(stats::rnorm(2e5), 1e5, 2))
  res <- anovakit:::.mardia_test(Y)            # the n x n form would need 75 GB
  expect_false(is.null(res$table))
  # b2 against its direct definition, and b1 against the n^2 sum on a subsample
  Ys <- Y[1:300, ]
  Yc <- scale(Ys, scale = FALSE)
  D <- Yc %*% solve(crossprod(Yc) / 300) %*% t(Yc)
  small <- anovakit:::.mardia_test(Ys)$table
  expect_equal(small$statistic[1], 300 * (sum(D^3) / 300^2) / 6)
  expect_equal(small$statistic[2],
               (mean(diag(D)^2) - 8) / sqrt(8 * 2 * 4 / 300))
})

test_that("F131: Mardia's tests are skipped when the residual df barely exceed p", {
  for (s in 1:2) {
    d <- withr::with_seed(s, data.frame(
      g = rep(c("a", "b", "c"), each = 4),
      matrix(if (s == 1) stats::rnorm(108) else stats::rexp(108), 12)))
    f <- anova_manova(d, paste0("X", 1:9), "g", plots = FALSE)
    expect_null(f$assumptions$mardia)
    expect_true(any(grepl("Mardia's tests were not computed", f$notes)))
    expect_false(any(grepl("reject multivariate normality", f$notes)))
  }
})

test_that("F025: responses on very different scales give the same Mardia and canonical results", {
  skip_if_not_installed("MASS")
  d <- iris
  d$Sepal.Length <- d$Sepal.Length * 1e5
  d$Petal.Width <- d$Petal.Width * 1e-4
  f <- anova_manova(d, names(iris)[1:4], "Species", plots = FALSE)
  f0 <- anova_manova(iris, names(iris)[1:4], "Species", plots = FALSE)
  expect_false(is.null(f$canonical))
  expect_false(is.null(f$assumptions$mardia))
  expect_equal(f$canonical$eigenvalue, f0$canonical$eigenvalue, tolerance = 1e-8)
  expect_equal(f$assumptions$mardia, f0$assumptions$mardia, tolerance = 1e-8)
  expect_equal(f$assumptions$structure_coefficients,
               f0$assumptions$structure_coefficients, tolerance = 1e-8)
  # MASS::lda: eigenvalue = svd^2 * (k - 1) / (n - k)
  l <- MASS::lda(Species ~ ., data = iris)
  expect_equal(f$canonical$eigenvalue, l$svd^2 * 2 / 147, tolerance = 1e-8)
  expect_equal(f$canonical$eigenvalue, c(32.19193, 0.2853910), tolerance = 1e-6)
})

# --- F076 / F109: the canonical tables --------------------------------------------

test_that("F076: the structure coefficients name their responses", {
  fit <- anova_manova(iris, names(iris)[1:4], "Species", plots = FALSE)
  sc <- fit$assumptions$structure_coefficients
  expect_identical(sc$response, names(iris)[1:4])
  out <- utils::capture.output(print(summary(fit)))
  i <- grep("structure_coefficients", out)
  expect_true(any(grepl("Petal.Width", out[i + 1:6])))
})

test_that("F109: canonical axes have unit pooled within-group variance and match MASS::lda", {
  skip_if_not_installed("MASS")
  d <- mv_fx()
  fit <- anova_manova(d, c("score1", "score2"), "g")
  sc <- fit$plots$canonical$data
  S <- as.matrix(sc[, c("Can1", "Can2")])
  W <- apply(S, 2, function(s) s - stats::ave(s, d$g))
  expect_equal(unname(colSums(W^2)) / (nrow(d) - 3), c(1, 1))
  expect_equal(sum(fit$canonical$prop_variance), 1)
  expect_equal(fit$canonical$canonical_r,
               sqrt(fit$canonical$eigenvalue / (1 + fit$canonical$eigenvalue)))
  # Structure coefficients: pooled within-group correlations
  Yw <- apply(as.matrix(d[, c("score1", "score2")]), 2,
              function(y) y - stats::ave(y, d$g))
  expect_equal(unname(as.matrix(fit$assumptions$structure_coefficients[, -1L])),
               unname(stats::cor(Yw, W)))
  # MASS::lda scores have unit pooled within-group variance too: same up to sign
  l <- MASS::lda(g ~ score1 + score2, data = d)
  ref <- stats::predict(l)$x
  expect_equal(abs(unname(S)), abs(unname(ref)), tolerance = 1e-8)
  expect_equal(fit$canonical$eigenvalue, l$svd^2 * 2 / 117, tolerance = 1e-8)
})

# --- Wiring ------------------------------------------------------------------------

test_that("per-response means and comparisons come from .emm_block and emmeans", {
  d <- mv_fx()
  d$site <- factor(rep(c("n", "s"), length.out = nrow(d)))
  fit <- anova_manova(d, c("score1", "score2"), c("g", "site"), plots = FALSE,
                      assumptions = FALSE)
  # Additive factors: per-factor tables with a term column, under $univariate too
  expect_true(all(c("response", "term") %in% names(fit$emmeans)))
  ref <- as.data.frame(emmeans::emmeans(stats::lm(score2 ~ g + site, data = d), ~ site))
  got <- fit$emmeans[fit$emmeans$response == "score2" & fit$emmeans$term == "site", ]
  expect_equal(got$estimate, ref$emmean)
  expect_identical(fit$univariate$score2$emmeans$estimate,
                   fit$emmeans$estimate[fit$emmeans$response == "score2"])
  ph <- fit$posthoc[fit$posthoc$response == "score1" & fit$posthoc$term == "g", ]
  ref_ph <- as.data.frame(summary(emmeans::contrast(
    emmeans::emmeans(stats::lm(score1 ~ g + site, data = d), ~ g), "pairwise"),
    adjust = "tukey"))
  expect_equal(ph$p_adjusted, ref_ph$p.value)
  # The univariate fits refer to the analysed frame only
  expect_identical(ls(environment(stats::formula(fit$univariate$score1$model))),
                   "data_used")
  expect_setequal(names(fit$data_used), c("score1", "score2", "g", "site"))
  expect_s3_class(fit$model, "mlm")
})

# --- Findings shared with anova_ancova (F160, F061, F097) ------------------------

test_that("F160: a covariate named like a contrast coefficient does not corrupt the means", {
  d <- withr::with_seed(1, {
    n <- 90
    d <- data.frame(dose = rep(c("high", "low", "mid"), each = 30),
                    dose1 = stats::rnorm(n, 50, 10))
    d$y <- 2 + c(high = 3, low = 0, mid = 1)[d$dose] + 0.05 * d$dose1 +
      stats::rnorm(n)
    d$z <- 1 + c(high = 1, low = 0, mid = 2)[d$dose] - 0.03 * d$dose1 +
      stats::rnorm(n)
    d
  })
  ok <- d
  names(ok)[names(ok) == "dose1"] <- "baseline"
  ref <- function(resp) {
    m <- stats::lm(stats::reformulate(c("baseline", "dose"), resp), data = ok)
    unname(stats::predict(m, data.frame(dose = c("high", "low", "mid"),
                                        baseline = mean(ok$baseline))))
  }
  # Type III names the sum-to-zero coefficients dose1, dose2
  bad <- anova_manova(d, c("y", "z"), "dose", covariates = "dose1",
                      type = "III", plots = FALSE)
  expect_equal(bad$emmeans$estimate[bad$emmeans$response == "y"], ref("y"))
  expect_equal(bad$emmeans$estimate[bad$emmeans$response == "z"], ref("z"))
  r <- ref("y")
  expect_equal(bad$posthoc$estimate[bad$posthoc$response == "y"],
               c(r[1] - r[2], r[1] - r[3], r[2] - r[3]))
  expect_identical(anyDuplicated(names(stats::coef(bad$univariate$y$model))), 0L)
  expect_identical(anyDuplicated(rownames(stats::coef(bad$model))), 0L)
  expect_true(any(grepl("put in brackets", bad$notes, fixed = TRUE)))
  good <- anova_manova(ok, c("y", "z"), "dose", covariates = "baseline",
                       type = "III", plots = FALSE)
  expect_equal(bad$multivariate[, -1L], good$multivariate[, -1L], ignore_attr = TRUE)
  expect_equal(bad$canonical, good$canonical)

  # Type II names the treatment coefficients dose<level>
  d2 <- d
  names(d2)[names(d2) == "dose1"] <- "doselow"
  bad2 <- anova_manova(d2, c("y", "z"), "dose", covariates = "doselow",
                       plots = FALSE)
  expect_equal(bad2$emmeans$estimate[bad2$emmeans$response == "y"], ref("y"))
})

test_that("F061: a binary covariate is held at its mean, not at the midpoint of its values", {
  d <- withr::with_seed(8, {
    n <- 120
    d <- data.frame(g = factor(rep(c("a", "b", "c"), each = n / 3)),
                    smoker = stats::rbinom(n, 1, 0.2))
    d$y1 <- 1 + (d$g == "b") + 2 * d$smoker + stats::rnorm(n)
    d$y2 <- 2 - (d$g == "c") + 1 * d$smoker + stats::rnorm(n)
    d
  })
  expect_false(isTRUE(all.equal(mean(d$smoker), 0.5)))
  fit <- anova_manova(d, c("y1", "y2"), "g", covariates = "smoker", plots = FALSE)
  for (r in c("y1", "y2")) {
    m <- stats::lm(stats::reformulate(c("smoker", "g"), r), data = d)
    ref <- unname(stats::predict(m, data.frame(g = c("a", "b", "c"),
                                               smoker = mean(d$smoker))))
    expect_equal(fit$emmeans$estimate[fit$emmeans$response == r], ref, info = r)
  }
})

test_that("F097: $data_used holds the analysed columns only, and unrelated columns cannot break the call", {
  d <- mv_fx()
  d$notes <- NA                       # unrelated and entirely missing
  d$unused <- letters[seq_len(nrow(d)) %% 26 + 1L]
  d <- cbind(d, stats::setNames(data.frame(seq_len(nrow(d))), ""))  # unnamed column
  fit <- anova_manova(d, c("score1", "score2"), "g", covariates = "age",
                      plots = FALSE)
  expect_setequal(names(fit$data_used), c("score1", "score2", "g", "age"))
  expect_identical(nrow(fit$data_used) + fit$n_removed, nrow(d))
  expect_identical(fit$n_removed, 0L)
  ref <- stats_row(summary(stats::manova(cbind(score1, score2) ~ age + g,
                                          data = d[, c("score1", "score2", "g", "age")])), "g")
  expect_equal(fit$multivariate$p_value[fit$multivariate$term == "g"], ref[["p"]])
})
