# Regression tests for the issues a second, independent verification pass
# found in the first round of fixes. Each is checked against an independent
# computation: the expanded data, a direct root-find, or base R.

# Robust standard errors under frequency weights ------------------------------

ucb_frames <- function() {
  ucb <- as.data.frame(UCBAdmissions)
  ucb$y <- as.integer(ucb$Admit == "Admitted")
  ex <- ucb[rep(seq_len(nrow(ucb)), ucb$Freq), ]
  row.names(ex) <- NULL
  list(weighted = ucb, expanded = ex)
}

test_that("frequency-weighted robust SEs equal those of the expanded data, every HC type", {
  skip_if_not_installed("sandwich")
  u <- ucb_frames()
  mw <- stats::glm(y ~ Gender + Dept, data = u$weighted, family = binomial,
                   weights = Freq)
  me <- stats::glm(y ~ Gender + Dept, data = u$expanded, family = binomial)
  for (ty in c("HC0", "HC1", "HC2", "HC3", "HC4")) {
    expect_equal(.freq_sandwich(mw, ty), sandwich::vcovHC(me, type = ty),
                 tolerance = 1e-6, ignore_attr = TRUE, info = ty)
  }
})

test_that("anova_bin and anova_glm give the expanded data's robust results", {
  skip_if_not_installed("sandwich")
  u <- ucb_frames()
  fw <- anova_bin(u$weighted, "y", c("Gender", "Dept"), weights = "Freq",
                  vcov_type = "HC0", ci_method = "wald", plots = FALSE)
  fe <- anova_bin(u$expanded, "y", c("Gender", "Dept"), vcov_type = "HC0",
                  ci_method = "wald", plots = FALSE)
  expect_equal(fw$effect_sizes, fe$effect_sizes, tolerance = 1e-6)
  expect_equal(fw$emmeans, fe$emmeans, tolerance = 1e-6)
  expect_equal(fw$posthoc, fe$posthoc, tolerance = 1e-6)
  expect_true(any(grepl("4,526 observations the frequency weights represent",
                        fw$notes, fixed = TRUE)))
  # Row by row the sandwich overstated SE(GenderFemale) twelvefold.
  se_or <- log(fw$effect_sizes$conf_high[1] / fw$effect_sizes$odds_ratio[1]) /
    stats::qnorm(0.975)
  expect_lt(se_or, 0.1)

  gw <- anova_glm(u$weighted, "y", c("Gender", "Dept"), family = "binomial",
                  weights = "Freq", vcov_type = "HC3", plots = FALSE)
  ge <- anova_glm(u$expanded, "y", c("Gender", "Dept"), family = "binomial",
                  vcov_type = "HC3", plots = FALSE)
  expect_equal(gw$assumptions$coefficients, ge$assumptions$coefficients,
               tolerance = 1e-6)
})

test_that("a Poisson frequency table gives the robust results of the raw counts", {
  skip_if_not_installed("sandwich")
  tab <- as.data.frame(table(spray = InsectSprays$spray, count = InsectSprays$count))
  tab$count <- as.integer(as.character(tab$count))
  tab <- tab[tab$Freq > 0, ]
  for (mt in c("poisson", "quasipoisson")) {
    fw <- anova_count(tab, "count", "spray", model = mt, weights = "Freq",
                      vcov_type = "HC3", test_statistic = "Wald", plots = FALSE)
    fe <- anova_count(InsectSprays, "count", "spray", model = mt,
                      vcov_type = "HC3", test_statistic = "Wald", plots = FALSE)
    expect_equal(fw$anova, fe$anova, tolerance = 1e-6, info = mt)
    expect_equal(fw$emmeans, fe$emmeans, tolerance = 1e-6, info = mt)
  }
})

test_that("sampling weights keep each row as the unit of the sandwich", {
  skip_if_not_installed("sandwich")
  set.seed(31)
  d <- data.frame(g = gl(3, 20), y = rbinom(60, 1, 0.5), w = runif(60, 0.5, 2))
  fit <- suppressWarnings(anova_bin(d, "y", "g", weights = "w", vcov_type = "HC0",
                                    ci_method = "wald", plots = FALSE))
  m <- suppressWarnings(stats::glm(y ~ g, data = d, family = binomial, weights = w))
  V <- sandwich::vcovHC(m, type = "HC0")
  expect_equal(log(fit$effect_sizes$conf_high / fit$effect_sizes$odds_ratio) /
                 stats::qnorm(0.975), unname(sqrt(diag(V))[-1]), tolerance = 1e-6)
  expect_true(any(grepl("treat each row as one independent unit", fit$notes)))
})

test_that("a single-observation cell falls back to model-based SEs for every HC type", {
  skip_if_not_installed("sandwich")
  d <- data.frame(g = factor(c(rep("a", 8), rep("b", 8), "c")),
                  y = c(1:8, 3:10, 4))
  for (ty in c("HC0", "HC1", "HC3")) {
    fit <- anova_glm(d, "y", "g", vcov_type = ty, plots = FALSE)
    ref <- anova_glm(d, "y", "g", plots = FALSE)
    expect_equal(fit$emmeans$se, ref$emmeans$se, info = ty)
    expect_true(all(fit$emmeans$se > 0), info = ty)
    expect_true(any(grepl("leverage 1", fit$notes)), info = ty)
  }
})

# ANCOVA with many grouping levels ---------------------------------------------

test_that("the covariate check agrees with the model matrix and scales to large designs", {
  set.seed(32)
  d <- data.frame(a = factor(sample(letters[1:3], 40, TRUE)),
                  b = factor(sample(letters[1:4], 40, TRUE)),
                  x = rnorm(40))
  for (int in list(TRUE, FALSE)) {
    terms <- .group_terms(c("a", "b"), int)
    G <- stats::model.matrix(stats::reformulate(terms), data = d)
    cell <- interaction(d$a, d$b, drop = TRUE)
    x_cell <- stats::ave(d$x, cell)
    d$z <- x_cell
    determined <- qr(cbind(G, d$z))$rank == qr(G)$rank
    got <- tryCatch(.check_covariates_identified(d, "z", c("a", "b"), int),
                    error = function(e) FALSE)
    expect_identical(isTRUE(got), !determined, info = paste("interaction", int))
  }
  # Three 22-level factors: 10648 cells for 300 rows. Refused at once, where
  # it used to spend minutes in qr() and then in lm().
  big <- data.frame(f1 = factor(sample(1:22, 300, TRUE)),
                    f2 = factor(sample(1:22, 300, TRUE)),
                    f3 = factor(sample(1:22, 300, TRUE)),
                    cv = rnorm(300), y = rnorm(300))
  elapsed <- system.time(expect_error(
    anova_ancova(big, "y", c("f1", "f2", "f3"), "cv", plots = FALSE),
    "10,649 coefficients but there are only 300 rows"))[["elapsed"]]
  expect_lt(elapsed, 30)
  # Additively the same factors are estimable and fit quickly.
  fit <- anova_ancova(big, "y", c("f1", "f2", "f3"), "cv", interaction = FALSE,
                      posthoc = FALSE, plots = FALSE)
  expect_identical(length(stats::coef(fit$model)), 65L)
})

# Marginal means that cannot be computed ---------------------------------------

test_that("every function says when the marginal means, comparisons and plot are missing", {
  old <- getOption("emmeans")
  withr::defer(options(emmeans = old))
  emmeans::emm_options(rg.limit = 1)
  fits <- list(
    bin = anova_bin(fx_binary(), "y", "g"),
    glm = anova_glm(fx_oneway(), "value", "group"),
    ancova = anova_ancova(fx_ancova(), "dv", "iv", "cov"),
    manova = anova_manova(fx_multivariate(), c("score1", "score2"), "g",
                          assumptions = FALSE)
  )
  for (nm in names(fits)) {
    f <- fits[[nm]]
    expect_null(f$emmeans)
    expect_null(f$plots$emmeans)
    expect_true(any(grepl("plot of estimated marginal means was skipped",
                          f$notes, fixed = TRUE)), info = nm)
    expect_true(any(grepl("Pairwise comparisons were not computed because",
                          f$notes, fixed = TRUE)), info = nm)
  }
  add <- anova_glm(fx_twoway(), "value", c("g1", "g2"), plots = FALSE)
  expect_false(any(grepl("reported for each factor separately", add$notes)))
})

# Intervals of a diverged coefficient ------------------------------------------

zero_group <- function() {
  set.seed(1)
  data.frame(g = factor(rep(c("a", "b", "c"), each = 10)),
             y = c(stats::rpois(10, 3), stats::rpois(10, 4), rep(0, 10)))
}

test_that("an all-zero Poisson group gets a finite profile end and profiled neighbours", {
  d <- zero_group()
  fit <- anova_glm(d, "y", "g", family = "poisson", plots = FALSE)
  ct <- fit$assumptions$coefficients
  m <- stats::glm(y ~ g, data = d, family = poisson)
  X <- stats::model.matrix(m)
  dev <- function(b) stats::glm.fit(X[, -3], d$y, offset = b * X[, 3],
                                    family = stats::poisson())$deviance
  end <- stats::uniroot(function(b) dev(b) - stats::deviance(m) - stats::qchisq(0.95, 1),
                        c(-10, -1), tol = 1e-10)$root
  expect_equal(ct$conf_high[ct$term == "gc"], end, tolerance = 1e-5)
  expect_identical(ct$conf_low[ct$term == "gc"], -Inf)
  # The other coefficients are profiled, not replaced by Wald intervals.
  expect_identical(unique(ct$ci_method), "profile")
  ref <- suppressMessages(stats::confint(stats::glm(y ~ g, data = d[d$g != "c", ],
                                                    family = poisson)))
  expect_equal(ct$conf_low[1:2], unname(ref[, 1]), tolerance = 1e-4)
  expect_false(any(grepl("no value on that side can be ruled out", fit$notes)))
  expect_false(any(grepl("unbounded on at least one side (NA)", fit$notes, fixed = TRUE)))
})

test_that("a separated binomial coefficient has no stale 'NA' note", {
  d <- zero_group()
  d$y <- as.integer(d$y > 2)
  d$y[d$g == "a"][1:3] <- 0L
  fit <- anova_glm(d, "y", "g", family = "binomial", plots = FALSE)
  ct <- fit$assumptions$coefficients
  expect_true(is.finite(ct$conf_high[ct$term == "gc"]))
  expect_false(any(grepl("(NA)", fit$notes, fixed = TRUE)))
})

# Response levels are checked on the analysed rows -----------------------------

test_that("a binary factor response is checked after missing rows are dropped", {
  # "maybe" occurs only in a row dropped for its missing group.
  d <- data.frame(g = factor(c(rep("a", 6), rep("b", 5), NA)),
                  y = factor(c("no", "no", "yes", "no", "yes", "no",
                               "no", "yes", "no", "yes", "no", "maybe")))
  fit <- anova_glm(d, "y", "g", family = "binomial", plots = FALSE)
  expect_identical(fit$n_removed, 1L)
  expect_true(any(grepl("equals \"yes\"", fit$notes)))
  # Only one level is left once the missing row is gone.
  d2 <- data.frame(g = factor(c(rep("a", 5), rep("b", 5), NA)),
                   y = factor(c(rep("no", 10), "yes")))
  expect_error(anova_glm(d2, "y", "g", family = "binomial", plots = FALSE),
               "1 observed level")
})

# Result columns named like grouping columns ------------------------------------

test_that("per-factor marginal means keep every statistic when a group is named like one", {
  set.seed(33)
  base <- data.frame(y = rnorm(60), a = gl(3, 20), b = gl(2, 10, 60))
  for (nm in c("se", "df", "conf_low", "estimate")) {
    d <- base
    names(d)[names(d) == "b"] <- nm
    fit <- anova_glm(d, "y", c("a", nm), conf_level = 0.9)
    ref <- anova_glm(base, "y", c("a", "b"), conf_level = 0.9, plots = FALSE)
    stat <- .emm_col(fit$emmeans, "estimate", c("a", nm))
    expect_equal(fit$emmeans[[stat]], ref$emmeans$estimate, info = nm)
    se <- .emm_col(fit$emmeans, "se", c("a", nm))
    expect_false(anyNA(fit$emmeans[[se]]), info = nm)
    expect_s3_class(fit$plots$emmeans, "ggplot")
    expect_no_warning(ggplot2::ggplot_build(fit$plots$emmeans))
  }
  d <- base
  names(d)[2:3] <- c("estimate", "emm_estimate")
  fit <- anova_glm(d, "y", c("estimate", "emm_estimate"), interaction = TRUE)
  ref <- anova_glm(base, "y", c("a", "b"), interaction = TRUE, plots = FALSE)
  stat <- .emm_col(fit$emmeans, "estimate", c("estimate", "emm_estimate"))
  expect_identical(stat, "emm_estimate2")
  expect_equal(fit$emmeans[[stat]], ref$emmeans$estimate)
  expect_s3_class(fit$plots$emmeans, "ggplot")
})

# The stored call ----------------------------------------------------------------

test_that("update() and step() work on the stored model from any frame", {
  set.seed(34)
  d <- data.frame(g = gl(3, 10), h = gl(2, 5, 30), y = stats::rpois(30, 4))
  fit <- anova_glm(d, "y", c("g", "h"), family = "poisson", plots = FALSE)
  refit <- local(function(m) stats::update(m, . ~ 1))
  m0 <- refit(fit$model)
  expect_equal(stats::coef(m0),
               stats::coef(stats::glm(y ~ 1, data = d, family = poisson)))
  expect_s3_class(stats::step(fit$model, trace = 0), "glm")
  lr <- stats::anova(m0, fit$model, test = "LRT")
  expect_equal(lr$Deviance[2], stats::deviance(m0) - stats::deviance(fit$model))
})

# Printing p-values at the precision they were computed with ----------------------

test_that("Tukey p-values below ptukey's precision print as a bound", {
  df <- data.frame(contrast = c("a - b", "a - c", "b - c"),
                   p_value = c(1e-20, 1e-15, 0.02),
                   p_adjusted = c(3.963e-14, 2e-11, 0.05),
                   adjustment = "tukey")
  out <- .round_df(df, 4)
  expect_identical(out$p_adjusted[1:2], c("<1e-10", "<1e-10"))
  expect_match(out$p_value[1], "^<")
  expect_identical(out$p_value[2], "1e-15")
  df$adjustment <- "bonferroni"
  expect_identical(.round_df(df, 4)$p_adjusted[1], "3.963e-14")
  expect_identical(.round_df(data.frame(p_value = c(NA, 0.5)), 4)$p_value[1], "NA")
})
