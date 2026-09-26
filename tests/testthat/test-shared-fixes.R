# Regression tests for the shared helpers: input preparation, marginal means,
# robust standard errors, separation, printing. Each is checked against an
# independent computation rather than the package's own internals.

# Input preparation ---------------------------------------------------------

test_that("non-ASCII labels with an unknown encoding do not crash (UTF-8 session)", {
  skip_if_not(isTRUE(l10n_info()[["UTF-8"]]), "needs a UTF-8 session")
  lab <- c("Zürich", "Bern", "Genève")
  native <- enc2native(lab)
  Encoding(native) <- "unknown"
  d <- data.frame(site = rep(native, each = 4),
                  y = c(1.2, 2.1, 1.7, 1.5, 3.1, 2.9, 3.6, 3.3, 2.2, 2.8, 2.0, 2.4),
                  stringsAsFactors = FALSE)
  fit <- anova_welch(d, "y", "site", plots = FALSE)
  ref <- stats::oneway.test(y ~ site, data = d, var.equal = FALSE)
  expect_equal(fit$anova$statistic, unname(ref$statistic))
  # row order does not change the outcome
  fit2 <- anova_welch(d[c(5:12, 1:4), ], "y", "site", plots = FALSE)
  expect_equal(fit2$anova$statistic, fit$anova$statistic)
})

test_that("cells whose joined labels collide are refused, not merged", {
  d <- data.frame(g1 = c(rep("a", 6), rep("a : b", 6)),
                  g2 = c(rep("b : c", 6), rep("c", 6)),
                  y = c(1:6, 11:16))
  expect_error(anova_kw(d, "y", c("g1", "g2"), plots = FALSE),
               "same label for different combinations")
})

test_that("an explicit NA factor level is treated as missing and counted", {
  d <- data.frame(g = factor(c(rep("a", 6), rep("b", 6), NA, NA), exclude = NULL),
                  y = c(stats::rnorm(12), 5, 6))
  expect_true(anyNA(levels(d$g)))
  fit <- anova_welch(d, "y", "g", plots = FALSE)
  expect_identical(fit$n_removed, 2L)
  expect_identical(nrow(fit$data_used) + fit$n_removed, nrow(d))
  expect_false(anyNA(fit$emmeans$group))
  expect_true(any(grepl("explicit factor level", fit$notes)))
})

test_that("distinct instants that print alike stay distinct groups", {
  withr::local_timezone("America/New_York")
  t <- as.POSIXct(c(1636263000, 1636266600), origin = "1970-01-01")
  expect_identical(format(t[1]), format(t[2]))   # the repeated DST hour
  d <- data.frame(t = rep(t, each = 5), y = c(1:5, 11:15))
  fit <- anova_welch(d, "y", "t", plots = FALSE)
  expect_identical(nrow(fit$emmeans), 2L)
  expect_equal(unname(fit$anova$p_value),
               stats::oneway.test(y ~ factor(as.numeric(t)), data = d)$p.value)
})

test_that("a column cannot play two roles in any function", {
  d <- fx_binary()
  expect_error(anova_bin(d, "y", c("g", "y")), "both the response and a grouping variable")
  expect_error(anova_glm(fx_oneway(), "value", c("group", "value")), "both the response")
  cnt <- fx_counts()
  expect_error(anova_count(cnt, "count", "g1", offset = "g1"),
               "both a grouping variable and the exposure offset")
})

test_that("a column called Residuals is refused by name", {
  d <- fx_oneway(); names(d)[1] <- "Residuals"
  expect_error(anova_glm(d, "value", "Residuals", plots = FALSE), "reserved")
})

# Weights -------------------------------------------------------------------

test_that("zero-weight rows are dropped, so robust SEs match the subset fit", {
  skip_if_not_installed("sandwich")
  d <- fx_oneway(sds = c(1, 1.5, 2))
  d$w <- rep(c(1, 0), length.out = nrow(d))
  fit_w <- anova_glm(d, "value", "group", weights = "w", vcov_type = "HC3",
                     plots = FALSE)
  fit_s <- anova_glm(d[d$w > 0, ], "value", "group", vcov_type = "HC3",
                     plots = FALSE)
  ref <- sqrt(diag(sandwich::vcovHC(stats::glm(value ~ group, data = d[d$w > 0, ]),
                                    type = "HC3")))
  expect_equal(fit_w$assumptions$coefficients$se, unname(ref))
  expect_equal(fit_w$posthoc$se, fit_s$posthoc$se)
  expect_identical(fit_w$n_removed, sum(d$w == 0))
  expect_identical(nrow(fit_w$data_used) + fit_w$n_removed, nrow(d))
})

test_that("a group whose weights are all zero is reported, not analysed as empty", {
  d <- fx_oneway(n_per = 5)[1:10, ]
  d$group <- droplevels(d$group)
  d$w <- c(rep(1, 5), rep(0, 5))
  expect_error(anova_glm(d, "value", "group", weights = "w", plots = FALSE),
               "at least 2 are required")
})

# Marginal means and comparisons ---------------------------------------------

test_that("a two-valued covariate is held at its mean in the marginal means", {
  withr::local_seed(11)
  d <- data.frame(g = factor(rep(c("a", "b", "c"), length.out = 200)),
                  sex = stats::rbinom(200, 1, 0.2))
  d$y <- 5 + 10 * d$sex + as.integer(d$g) + stats::rnorm(200)
  m <- stats::lm(y ~ sex + g, data = d)
  e <- anovakit:::.emmeans_grid(m, "g", type = "link")
  est <- as.data.frame(summary(e$grid))$emmean
  ref <- stats::predict(m, data.frame(g = levels(d$g), sex = mean(d$sex)))
  expect_equal(est, unname(ref))
})

test_that("global emm_options cannot change the reported comparisons", {
  d <- fx_oneway()
  a <- anova_glm(d, "value", "group", plots = FALSE)
  old <- getOption("emmeans")
  withr::defer(options(emmeans = old))
  emmeans::emm_options(contrast = list(adjust = "none", side = ">"),
                       emmeans = list(level = 0.5))
  b <- anova_glm(d, "value", "group", plots = FALSE)
  expect_equal(b$posthoc$p_adjusted, a$posthoc$p_adjusted)
  expect_equal(b$emmeans$conf_low, a$emmeans$conf_low)
})

test_that("quasi-Poisson marginal comparisons use t, like the coefficient table", {
  withr::local_seed(3)
  d <- data.frame(g = factor(rep(c("a", "b"), each = 6)),
                  y = c(stats::rpois(6, 4), stats::rpois(6, 9)))
  fit <- anova_glm(d, "y", "g", family = "quasipoisson", plots = FALSE)
  ct <- stats::coef(summary(fit$model))
  expect_equal(fit$posthoc$p_value, unname(ct["gb", "Pr(>|t|)"]), tolerance = 1e-6)
  expect_equal(unique(fit$emmeans$df), stats::df.residual(fit$model))
})

test_that("comparisons under an inverse link are response-scale differences", {
  withr::local_seed(1)
  d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 40)))
  d$y <- stats::rgamma(120, shape = 5, rate = 5 / c(a = 2, b = 4, c = 8)[as.character(d$g)])
  fit <- anova_glm(d, "y", "g", family = stats::Gamma(), plots = FALSE)
  expect_false("ratio" %in% names(fit$posthoc))
  means <- stats::setNames(fit$emmeans$estimate, fit$emmeans$g)
  expect_equal(fit$posthoc$estimate[1], unname(means["a"] - means["b"]),
               tolerance = 1e-8)
})

test_that("log-link comparisons are still ratios", {
  d <- fx_counts()
  fit <- anova_glm(d, "count", "g1", family = "poisson", plots = FALSE)
  expect_true("ratio" %in% names(fit$posthoc))
  m <- stats::setNames(fit$emmeans$estimate, fit$emmeans$g1)
  expect_equal(fit$posthoc$ratio, unname(m["A"] / m["B"]))
})

test_that("an additive two-factor model compares each factor on its own", {
  withr::local_seed(4)
  d <- expand.grid(rep = 1:10, A = factor(paste0("a", 1:4)), B = factor(paste0("b", 1:3)))
  d$y <- as.integer(d$A) * 0.3 + stats::rnorm(nrow(d))
  fit <- anova_glm(d, "y", c("A", "B"))
  expect_true("term" %in% names(fit$posthoc))
  expect_identical(nrow(fit$posthoc), as.integer(choose(4, 2) + choose(3, 2)))
  ref <- as.data.frame(summary(emmeans::contrast(
    emmeans::emmeans(fit$model, ~A), "pairwise", adjust = "tukey")))
  expect_equal(fit$posthoc$p_adjusted[fit$posthoc$term == "A"], ref$p.value)
  expect_identical(nrow(fit$emmeans), 7L)
  expect_true(any(grepl("additively", fit$notes)))
  expect_true(inherits(fit$emmeans_object$A, "emmGrid"))
  p <- fit$plots$emmeans
  expect_s3_class(p, "ggplot")
  expect_no_warning(ggplot2::ggplot_build(p))
})

test_that("a full factorial with an empty cell leaves that cell out and says so", {
  d <- fx_twoway()
  d <- d[!(d$g1 == "b" & d$g2 == "y"), ]
  fit <- anova_glm(d, "value", c("g1", "g2"), interaction = TRUE, plots = FALSE)
  expect_identical(nrow(fit$emmeans), 3L)
  expect_false(anyNA(fit$emmeans$estimate))
  expect_false(anyNA(fit$posthoc$estimate))
  expect_true(any(grepl("contain no data and are not estimable", fit$notes)))
})

test_that("too many comparisons are skipped with a note rather than attempted", {
  withr::local_seed(5)
  d <- data.frame(g = factor(sprintf("g%03d", rep(1:101, each = 2))))
  d$y <- stats::rnorm(nrow(d))
  fit <- anova_glm(d, "y", "g", plots = FALSE)
  expect_null(fit$posthoc)
  expect_true(any(grepl("5,050", fit$notes)))
})

test_that("posthoc tables carry raw and adjusted p-values and the method", {
  fit <- anova_glm(fx_oneway(), "value", "group", adjust = "holm", plots = FALSE)
  expect_true(all(c("p_value", "p_adjusted", "adjustment") %in% names(fit$posthoc)))
  expect_equal(fit$posthoc$p_adjusted, stats::p.adjust(fit$posthoc$p_value, "holm"))
  expect_identical(unique(fit$posthoc$adjustment), "holm")
  expect_true(any(grepl("intervals use the bonferroni adjustment", fit$notes)))
})

# Robust standard errors ----------------------------------------------------

test_that("HC3 with a single-observation cell falls back instead of returning NaN", {
  skip_if_not_installed("sandwich")
  d <- fx_oneway(n_per = 6)
  d <- rbind(d, data.frame(group = "D", value = 5))
  expect_silent_stdout(
    fit <- anova_glm(d, "value", "group", vcov_type = "HC3", plots = FALSE))
  expect_false(anyNA(fit$posthoc$se))
  expect_true(any(grepl("leverage 1", fit$notes)))
})

test_that("robust SEs are not used when fitted values reach the boundary", {
  skip_if_not_installed("sandwich")
  withr::local_seed(3)
  d <- data.frame(g = factor(rep(c("a", "b", "c"), each = 20)))
  d$y <- stats::rbinom(60, 1, c(a = .3, b = .5, c = .6)[as.character(d$g)])
  d$y[d$g == "c"] <- 1
  fit <- anova_bin(d, "y", "g", vcov_type = "HC3", plots = FALSE)
  expect_true(any(grepl("boundary", fit$notes)))
  expect_true(all(fit$posthoc$p_adjusted > 1e-10))
})

test_that("robust SEs with non-unit weights say what they assume", {
  skip_if_not_installed("sandwich")
  d <- fx_binary()
  d$w <- rep(c(1, 2), length.out = nrow(d))
  fit <- anova_bin(d, "y", "g", weights = "w", vcov_type = "HC0", plots = FALSE)
  expect_true(any(grepl("frequency weights", fit$notes)))
})

# Separation ------------------------------------------------------------------

test_that("separation is detected whatever the contrasts", {
  d <- fx_separated()
  for (tp in c("II", "III")) {
    fit <- anova_bin(d, "y", "g", type = tp, plots = FALSE)
    expect_true(any(grepl("g = c", fit$notes[grepl("separation", fit$notes)])),
                info = tp)
  }
})

test_that("separation is detected in a large sample", {
  withr::local_seed(1)
  d <- data.frame(g = factor(c(rep("a", 20000), rep("b", 20000), rep("c", 5))),
                  y = c(stats::rbinom(20000, 1, 0.9), stats::rbinom(20000, 1, 0.95), rep(1, 5)))
  fit <- anova_bin(d, "y", "g", posthoc = FALSE, plots = FALSE)
  expect_true(any(grepl("separation.*g = c", fit$notes)))
})

test_that("a rare but well-determined event rate is not called separation", {
  withr::local_seed(2)
  d <- data.frame(g = factor(rep(c("a", "b"), each = 20000)))
  d$y <- stats::rbinom(40000, 1, ifelse(d$g == "a", 0.001, 0.003))
  fit <- anova_bin(d, "y", "g", posthoc = FALSE, plots = FALSE)
  expect_false(any(grepl("separation", fit$notes)))
})

test_that("an all-zero count group is flagged", {
  withr::local_seed(6)
  d <- data.frame(g = factor(rep(c("A", "B", "C"), each = 20)))
  d$y <- c(stats::rpois(20, 3), stats::rpois(20, 5), rep(0, 20))
  fit <- anova_count(d, "y", "g", model = "poisson", plots = FALSE)
  expect_true(any(grepl("No events were observed in g = C", fit$notes)))
})

# Model helpers -------------------------------------------------------------

test_that("term labels are not backtick-quoted", {
  d <- fx_oneway(); names(d)[1] <- "my group"
  fit <- anova_glm(d, "value", "my group", plots = FALSE)
  expect_identical(fit$anova$term[1], "my group")
  expect_identical(fit$effect_sizes$term[1], "my group")
})

test_that("a saturated Poisson model keeps valid tests, and says so", {
  d <- data.frame(g = factor(c("a", "b", "c")), y = c(5, 12, 30))
  fit <- anova_glm(d, "y", "g", family = "poisson", plots = FALSE)
  expect_true(any(grepl("still valid", fit$notes)))
  expect_false(any(grepl("undefined", fit$notes)))
  expect_equal(fit$anova$statistic,
               stats::anova(fit$model, test = "LRT")[2, "Deviance"])
})

test_that("an aliased covariate is not blamed on empty cells", {
  d <- fx_oneway()
  d$x <- as.numeric(d$group)
  m <- stats::lm(value ~ group + x, data = d)
  n <- anovakit:::.check_model_size(m)
  expect_true(any(grepl("aliased", n)))
  expect_false(any(grepl("because the corresponding cells are empty", n)))
})

test_that("small-unit responses keep their ANOVA table", {
  d <- fx_oneway()
  d$value <- d$value * 1e-6
  m <- stats::lm(value ~ group, data = d)
  av <- anovakit:::.car_anova(m, "II")
  expect_false(is.null(av$table))
  ref <- stats::anova(m)
  expect_equal(av$table$statistic[1], ref[["F value"]][1])
  expect_equal(av$table$sum_sq[1], ref[["Sum Sq"]][1])
})

test_that("the fitted model refits with update() on the analysed rows", {
  d <- fx_binary()
  d$w <- rep(c(1, 2), length.out = nrow(d))
  fit <- anova_bin(d, "y", "g", weights = "w", plots = FALSE)
  m2 <- stats::update(fit$model, data = fit$data_used)
  expect_equal(stats::coef(m2), stats::coef(fit$model))
  # the model's environment holds the analysed frame, not the caller's data
  env <- environment(stats::formula(fit$model))
  expect_setequal(setdiff(ls(env), "family_used"), "data_used")
})

test_that("an offset model rebuilds its marginal means from the analysed rows", {
  withr::local_seed(1)
  d <- data.frame(g = rep(c("A", "B", "C"), each = 40), stringsAsFactors = FALSE)
  d$hours <- stats::runif(120, 0.5, 4); d$w <- stats::runif(120, 0.5, 2)
  d$y <- stats::rpois(120, c(2, 4, 6)[match(d$g, c("A", "B", "C"))] * d$hours)
  a <- anova_count(d, "y", "g", offset = "hours", weights = "w", plots = FALSE)
  expect_identical(nrow(a$emmeans), 3L)
  d2 <- rbind(d[, c("g", "hours", "y")], data.frame(g = "D", hours = 1, y = Inf))
  b <- anova_count(d2, "y", "g", offset = "hours", plots = FALSE)
  expect_false("D" %in% b$emmeans$g)
})

test_that("binomial factor responses need exactly two levels and name the modelled one", {
  d <- data.frame(g = factor(rep(c("a", "b"), each = 30)),
                  y = factor(rep(c("lo", "mid", "hi"), 20)))
  expect_error(anova_glm(d, "y", "g", family = "binomial"), "exactly two")
  d$y <- factor(rep(c("no", "yes"), 30))
  fit <- anova_glm(d, "y", "g", family = "quasibinomial", plots = FALSE)
  expect_true(any(grepl("equals \"yes\"", fit$notes)))
  d$y <- as.character(d$y)
  fit2 <- anova_glm(d, "y", "g", family = "binomial", plots = FALSE)
  expect_equal(stats::coef(fit2$model),
               stats::coef(stats::glm(factor(y) ~ g, data = d, family = binomial)))
})

test_that("per-group normality says why a group was not tested", {
  tab <- anovakit:::.normality_by_group(c(1, 2, stats::rnorm(10)), rep(c("a", "b"), c(2, 10)))
  expect_true(grepl("at least 3", tab$note[tab$group == "a"]))
  expect_true(is.na(tab$note[tab$group == "b"]))
})

# Printing ----------------------------------------------------------------------

test_that("printing keeps counts whole and formats tiny p-values", {
  df <- data.frame(term = "g", df = 12345, statistic = 1.234567, p_value = 1e-300)
  r <- anovakit:::.round_df(df, 4)
  expect_identical(r$df, 12345)
  expect_match(r$p_value, "<")
  expect_equal(r$statistic, 1.235)
})

test_that("digits is validated", {
  fit <- anova_welch(fx_oneway(), "value", "group", plots = FALSE)
  expect_error(print(fit, digits = NA), "digits")
  expect_silent(utils::capture.output(print(fit, digits = 40)))
})

test_that("the omnibus header names the statistic", {
  fit <- anova_bin(fx_binary(), "y", "g", plots = FALSE)
  out <- utils::capture.output(print(fit))
  expect_true(any(grepl("Omnibus test (Type II, likelihood-ratio chi-square)", out,
                        fixed = TRUE)))
})

test_that("summary prints the method's own extra tables", {
  skip_if_not_installed("afex")
  fit <- anova_rm(fx_repeated(), "score", "id", "time", plots = FALSE)
  out <- utils::capture.output(print(summary(fit)))
  expect_true(any(grepl("^Sphericity", out)))
})

# Plot helpers ----------------------------------------------------------------

test_that("plot labels keep the confidence level exact", {
  expect_identical(anovakit:::.pct(0.9999), "99.99")
  expect_identical(anovakit:::.pct(0.95), "95")
})

test_that("the covariate band uses conf_level", {
  d <- fx_ancova()
  p <- anovakit:::.plot_covariate(d, "dv", "cov", "iv", conf_level = 0.5)
  expect_equal(p$layers[[2]]$stat_params$level, 0.5)
  b <- ggplot2::layer_data(p, 2)
  p95 <- anovakit:::.plot_covariate(d, "dv", "cov", "iv", conf_level = 0.95)
  b95 <- ggplot2::layer_data(p95, 2)
  expect_true(mean(b$ymax - b$ymin) < mean(b95$ymax - b95$ymin))
})

test_that("short labels are not rotated", {
  th <- anovakit:::.gg_theme_auto(c("a", "b", "c"))
  expect_true(is.null(th$axis.text.x$angle) || th$axis.text.x$angle == 0)
})

test_that("a grouping column called emm_estimate is not plotted against itself", {
  df <- data.frame(emm_estimate = c("x", "y"), estimate = c(1, 2),
                   conf_low = c(0, 1), conf_high = c(2, 3))
  expect_identical(anovakit:::.emm_col(df, "estimate", "emm_estimate"), "estimate")
})

test_that("an empty-string group gets a box-plot label, not NA", {
  d <- data.frame(g = rep(c("", "b", "c"), each = 6), y = c(1:6, 11:16, 21:26))
  p <- anovakit:::.plot_box(transform(d, g = factor(g)), "y", "g")
  labs <- levels(p$data$label)
  expect_false(anyNA(labs))
  expect_identical(length(labs), 3L)
  expect_true(any(startsWith(labs, "\n(n = 6)")))
})
