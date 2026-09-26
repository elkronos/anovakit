# Package-wide guarantees, checked across all eight functions against
# independent computations: conf_level reaches every interval, and every plot
# builds cleanly and draws the numbers its table reports.

lvl <- 0.8
z80 <- stats::qnorm(1 - (1 - lvl) / 2)

emm_ref <- function(grid, type = NULL) {
  s <- if (is.null(type)) summary(grid, level = lvl, infer = c(TRUE, FALSE)) else
    summary(grid, level = lvl, infer = c(TRUE, FALSE), type = type)
  s <- as.data.frame(s)
  lo <- intersect(c("lower.CL", "asymp.LCL"), names(s))[1]
  hi <- intersect(c("upper.CL", "asymp.UCL"), names(s))[1]
  cbind(s[[lo]], s[[hi]])
}

pairs_ref <- function(grid, type = NULL) {
  ct <- emmeans::contrast(grid, "pairwise", adjust = "tukey")
  s <- as.data.frame(if (is.null(type)) confint(ct, level = lvl) else
    confint(ct, level = lvl, type = type))
  lo <- intersect(c("lower.CL", "asymp.LCL"), names(s))[1]
  hi <- intersect(c("upper.CL", "asymp.UCL"), names(s))[1]
  cbind(s[[lo]], s[[hi]])
}

test_that("conf_level reaches the Welch comparisons and group means", {
  d <- fx_oneway()
  fit <- anova_welch(d, "value", "group", conf_level = lvl, plots = FALSE)
  ab <- d$group %in% c("A", "B")
  tt <- stats::t.test(value ~ group, data = droplevels(d[ab, ]), conf.level = lvl)
  expect_equal(unlist(fit$posthoc[1, c("conf_low", "conf_high")]),
               tt$conf.int[1:2], ignore_attr = TRUE)
  a <- d$value[d$group == "A"]
  expect_equal(fit$emmeans$conf_low[1],
               mean(a) - stats::qt(0.9, length(a) - 1) * stats::sd(a) / sqrt(length(a)))
})

test_that("conf_level reaches anova_glm's means, comparisons and coefficients", {
  d <- fx_oneway()
  fit <- anova_glm(d, "value", "group", conf_level = lvl, ci_method = "wald",
                   plots = FALSE)
  expect_equal(unname(as.matrix(fit$emmeans[, c("conf_low", "conf_high")])),
               emm_ref(fit$emmeans_object))
  expect_equal(unname(as.matrix(fit$posthoc[, c("conf_low", "conf_high")])),
               pairs_ref(fit$emmeans_object))
  m <- stats::lm(value ~ group, data = d)
  expect_equal(unname(as.matrix(fit$assumptions$coefficients[, c("conf_low", "conf_high")])),
               unname(stats::confint(m, level = lvl)))
})

test_that("conf_level reaches anova_bin's probabilities, odds ratios and comparisons", {
  d <- fx_binary()
  fit <- anova_bin(d, "y", "g", conf_level = lvl, ci_method = "wald", plots = FALSE)
  expect_equal(unname(as.matrix(fit$emmeans[, c("conf_low", "conf_high")])),
               emm_ref(fit$emmeans_object))
  expect_equal(unname(as.matrix(fit$posthoc[, c("conf_low", "conf_high")])),
               pairs_ref(fit$emmeans_object))
  m <- stats::glm(y ~ g, data = d, family = binomial)
  ct <- stats::coef(summary(m))[-1, ]
  expect_equal(fit$effect_sizes$conf_low, unname(exp(ct[, 1] - z80 * ct[, 2])))
  expect_equal(fit$effect_sizes$conf_high, unname(exp(ct[, 1] + z80 * ct[, 2])))
})

test_that("conf_level reaches anova_count's rates and rate ratios", {
  d <- fx_counts()
  fit <- anova_count(d, "count", "g1", model = "poisson", conf_level = lvl,
                     plots = FALSE)
  expect_equal(unname(as.matrix(fit$emmeans[, c("conf_low", "conf_high")])),
               emm_ref(fit$emmeans_object))
  expect_equal(unname(as.matrix(fit$posthoc[, c("conf_low", "conf_high")])),
               pairs_ref(fit$emmeans_object))
  m <- stats::glm(count ~ g1, data = d, family = poisson)
  ct <- stats::coef(summary(m))[-1, , drop = FALSE]
  expect_equal(fit$effect_sizes$conf_low, unname(exp(ct[, 1] - z80 * ct[, 2])))
})

test_that("conf_level reaches anova_ancova's means, comparisons and slopes", {
  d <- fx_ancova(slope_shift = 0.5)
  fit <- anova_ancova(d, "dv", "iv", "cov", force_interaction = TRUE,
                      conf_level = lvl, plots = FALSE)
  expect_equal(unname(as.matrix(fit$emmeans[, c("conf_low", "conf_high")])),
               emm_ref(fit$emmeans_object))
  expect_equal(unname(as.matrix(fit$posthoc[, c("conf_low", "conf_high")])),
               pairs_ref(fit$emmeans_object))
  dc <- transform(d, cov = cov - mean(cov))
  tr <- as.data.frame(confint(emmeans::emtrends(stats::lm(dv ~ cov * iv, data = dc),
                                                ~ iv, var = "cov"), level = lvl))
  expect_equal(fit$simple_slopes$conf_low, tr$lower.CL, tolerance = 1e-8)
})

test_that("conf_level reaches anova_manova's per-response means", {
  d <- fx_multivariate()
  fit <- anova_manova(d, c("score1", "score2"), "g", conf_level = lvl,
                      assumptions = FALSE, plots = FALSE)
  g1 <- fit$univariate$score1$emmeans_object
  rows <- fit$emmeans[fit$emmeans[[1]] == "score1", ]
  expect_equal(unname(as.matrix(rows[, c("conf_low", "conf_high")])), emm_ref(g1))
})

test_that("conf_level reaches anova_rm's means and comparisons", {
  skip_if_not_installed("afex")
  fit <- anova_rm(fx_repeated(), "score", "id", "time", conf_level = lvl,
                  plots = FALSE)
  expect_equal(unname(as.matrix(fit$emmeans[, c("conf_low", "conf_high")])),
               emm_ref(fit$emmeans_object))
  expect_equal(unname(as.matrix(fit$posthoc[, c("conf_low", "conf_high")])),
               pairs_ref(fit$emmeans_object))
})

# Plots -------------------------------------------------------------------------

build_all <- function(fit) {
  for (nm in names(fit$plots)) {
    p <- fit$plots[[nm]]
    expect_s3_class(p, "ggplot")
    expect_no_warning(ggplot2::ggplot_build(p), message = nm)
    f <- tempfile(fileext = ".png")
    expect_no_warning(ggplot2::ggsave(f, p, width = 4, height = 3, dpi = 40))
    unlink(f)
  }
}

errorbars_match <- function(p, lo, hi) {
  ld <- ggplot2::layer_data(p, which(vapply(p$layers, function(l)
    inherits(l$geom, "GeomErrorbar"), logical(1)))[1])
  expect_equal(sort(ld$ymin), sort(lo))
  expect_equal(sort(ld$ymax), sort(hi))
}

test_that("every plot of every function builds, and error bars match the tables", {
  fits <- list(
    welch  = anova_welch(fx_oneway(), "value", "group", conf_level = lvl),
    kw     = anova_kw(fx_oneway(), "value", "group", diagnostics = TRUE),
    glm    = anova_glm(fx_oneway(), "value", "group", conf_level = lvl),
    glm2   = anova_glm(fx_twoway(), "value", c("g1", "g2"), conf_level = lvl),
    bin    = anova_bin(fx_binary(), "y", "g", conf_level = lvl),
    count  = anova_count(fx_counts(), "count", c("g1", "g2"), interaction = TRUE,
                         conf_level = lvl),
    ancova = anova_ancova(fx_ancova(), "dv", "iv", "cov", conf_level = lvl),
    manova = anova_manova(fx_multivariate(), c("score1", "score2"), "g",
                          conf_level = lvl)
  )
  if (requireNamespace("afex", quietly = TRUE)) {
    fits$rm <- anova_rm(fx_repeated(), "score", "id", "time", conf_level = lvl)
  }
  for (nm in names(fits)) build_all(fits[[nm]])

  w <- fits$welch
  errorbars_match(w$plots$means, w$emmeans$conf_low, w$emmeans$conf_high)
  for (nm in c("glm", "glm2", "bin", "count", "ancova", "rm")) {
    f <- fits[[nm]]
    if (is.null(f)) next
    errorbars_match(f$plots$emmeans, f$emmeans$conf_low, f$emmeans$conf_high)
    expect_match(f$plots$emmeans$labels$subtitle, "80% confidence", info = nm)
  }
  # A two-factor interaction grid draws one series per level of the second
  # factor, not one line doubling back.
  cnt <- ggplot2::layer_data(fits$count$plots$emmeans, 1)
  expect_identical(length(unique(cnt$group)), 2L)
})

test_that("MANOVA's marginal-means plot separates the second factor", {
  d <- fx_multivariate()
  d$h <- factor(rep(c("x", "y"), length.out = nrow(d)))
  fit <- anova_manova(d, c("score1", "score2"), c("g", "h"), interaction = TRUE,
                      assumptions = FALSE)
  p <- fit$plots$emmeans
  expect_no_warning(ggplot2::ggplot_build(p))
  ld <- ggplot2::layer_data(p, 1)
  expect_gt(length(unique(ld$colour)), 1L)
})
