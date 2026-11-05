# Regression tests for the third review round, output area. Each test checks
# the package against an independent reference: base R, emmeans on an
# equivalent model, or a hand computation.

pkg_file <- function(...) {
  # The sources are only available when testing from a source tree.
  p <- testthat::test_path("..", "..", ...)
  if (!file.exists(p)) testthat::skip(sprintf("%s is not available", file.path(...)))
  paste(readLines(p, warn = FALSE), collapse = "\n")
}

r3_empty_additive <- function() {
  withr::with_seed(3, {
    e <- expand.grid(rep = 1:8, g = c("a", "b", "c"), h = c("x", "y"))
    e <- e[!(e$g == "c" & e$h == "y"), ]
    e$y <- stats::rnorm(nrow(e)); e$y2 <- stats::rnorm(nrow(e))
    e$x <- stats::rnorm(nrow(e))
    e$cnt <- stats::rpois(nrow(e), 4); e$bin <- stats::rbinom(nrow(e), 1, .5)
    e
  })
}

test_that("R3-05: ratio comparisons keep ' - ' inside level labels", {
  d <- withr::with_seed(3, {
    d <- data.frame(g = rep(c("a - b", "a / b", "c"), each = 40))
    d$b <- stats::rbinom(120, 1, rep(c(.2, .5, .8), each = 40))
    d$k <- stats::rpois(120, rep(c(2, 4, 8), each = 40))
    d
  })
  # Reference: emmeans' link-scale contrasts, whose labels name the real levels;
  # the ratio label swaps only the separator between the two sides.
  ref_lab <- c("(a - b) / (a / b)", "(a - b) / c", "(a / b) / c")
  for (fit in list(
    bin = anova_bin(d, "b", "g", plots = FALSE),
    count = anova_count(d, "k", "g", model = "poisson", plots = FALSE),
    glm = anova_glm(d, "k", "g", family = stats::poisson(), plots = FALSE))) {
    fam <- stats::family(fit$model)
    ref <- as.data.frame(emmeans::contrast(
      emmeans::emmeans(stats::glm(stats::reformulate("g", all.vars(stats::formula(fit$model))[1L]),
                                  family = fam, data = d), "g"), "pairwise"))
    expect_identical(as.character(ref$contrast),
                     c("(a - b) - (a / b)", "(a - b) - c", "(a / b) - c"))
    expect_identical(as.character(fit$posthoc$contrast), ref_lab)
    expect_false(anyDuplicated(as.character(fit$posthoc$contrast)) > 0L)
    ratio <- fit$posthoc[[intersect(c("ratio", "IRR", "odds_ratio"), names(fit$posthoc))[1L]]]
    expect_equal(ratio, exp(ref$estimate), tolerance = 1e-6)
  }
  # Ordinary labels are unchanged.
  f <- anova_count(warpbreaks, "breaks", "tension", plots = FALSE)
  expect_identical(as.character(f$posthoc$contrast), c("L / M", "L / H", "M / H"))
})

test_that("R3-19: every additive model notes the empty cell its marginal means average over", {
  e <- r3_empty_additive()
  fits <- list(
    count = anova_count(e, "cnt", c("g", "h"), plots = FALSE),
    glm = anova_glm(e, "y", c("g", "h"), plots = FALSE),
    bin = anova_bin(e, "bin", c("g", "h"), plots = FALSE),
    manova = anova_manova(e, c("y", "y2"), c("g", "h"), plots = FALSE))
  # Reference: the empty combination, from the data.
  tab <- table(e$g, e$h)
  expect_identical(sum(tab == 0), 1L)
  expect_identical(tab["c", "y"], 0L)
  for (nm in names(fits)) {
    note <- grep("contain no observations", fits[[nm]]$notes, value = TRUE)
    expect_length(note, 1L)
    expect_match(note, "1 of the 6 possible", fixed = TRUE, info = nm)
    expect_match(note, "c : y", fixed = TRUE, info = nm)
  }
  # The note is true: g = c's reported mean averages the observed c:x cell with
  # the extrapolated c:y cell of the additive model.
  gl <- fits$glm
  cells <- summary(emmeans::emmeans(stats::lm(y ~ g + h, e), ~ g:h))
  expect_equal(gl$emmeans$estimate[gl$emmeans$term == "g" & gl$emmeans$g == "c"],
               mean(cells$emmean[cells$g == "c"]), tolerance = 1e-8)
  # A full design gives no such note.
  full <- expand.grid(rep = 1:8, g = c("a", "b"), h = c("x", "y"))
  full$y <- withr::with_seed(4, stats::rnorm(nrow(full)))
  expect_false(any(grepl("contain no observations",
                         anova_glm(full, "y", c("g", "h"), plots = FALSE)$notes)))
})

test_that("R3-24: the counts article describes the profiled one-sided interval", {
  incidents <- withr::with_seed(1, data.frame(
    site = factor(rep(c("north", "south", "east"), each = 10),
                  levels = c("north", "south", "east")),
    events = c(stats::rpois(10, 2), stats::rpois(10, 4), stats::rpois(10, 0.05))))
  z <- anova_count(incidents, "events", "site", plots = FALSE)
  es <- z$effect_sizes[z$effect_sizes$comparison == "east vs north", ]
  # Reference: the likelihood-ratio profile bound for east's log ratio.
  full <- stats::glm(events ~ site, stats::poisson(), incidents)
  f <- function(b) {
    stats::deviance(stats::glm(events ~ I(site == "south"), stats::poisson(), incidents,
                               offset = b * (incidents$site == "east"))) -
      stats::deviance(full) - stats::qchisq(.95, 1)
  }
  expect_equal(es$conf_high, exp(stats::uniroot(f, c(-10, 0), tol = 1e-10)$root),
               tolerance = 1e-4)
  expect_equal(es$conf_low, 0)
  expect_true(is.finite(es$conf_high))
  txt <- gsub("\\s+", " ", pkg_file("vignettes", "articles", "counts.Rmd"))
  expect_no_match(txt, "with an interval from 0 to infinity", fixed = TRUE)
  expect_match(txt, "profiled and one-sided", fixed = TRUE)
})

test_that("R3-26: the 5000-comparison cap is documented for the emmeans-based functions only", {
  d <- withr::with_seed(1, {
    d <- data.frame(g = factor(rep(sprintf("g%03d", 1:102), each = 3)))
    d$y <- stats::rnorm(nrow(d)); d
  })
  expect_identical(choose(102, 2), 5151)
  w <- anova_welch(d, "y", "g", plots = FALSE)
  k <- anova_kw(d, "y", "g", plots = FALSE)
  g <- anova_glm(d, "y", "g", plots = FALSE)
  expect_identical(NROW(w$posthoc), 5151L)
  expect_identical(NROW(k$posthoc), 5151L)
  expect_null(g$posthoc)
  expect_true(any(grepl("5,151, more than the 5,000", g$notes, fixed = TRUE)))
  txt <- gsub("#'\\s*", "", pkg_file("R", "anovakit-package.R"))
  txt <- gsub("\\s+", " ", txt)
  expect_no_match(txt, "Above 5000 comparisons they are skipped by default", fixed = TRUE)
  expect_match(txt, "Welch's and Dunn's comparisons are always computed in full", fixed = TRUE)
})

test_that("R3-28: the plots article reads the warpbreaks cell means correctly", {
  f <- anova_count(warpbreaks, "breaks", c("tension", "wool"), interaction = TRUE,
                   plots = FALSE)
  # Reference: the raw cell means, which the saturated Poisson model reproduces.
  m <- tapply(warpbreaks$breaks, warpbreaks[c("tension", "wool")], mean)
  em <- f$emmeans
  est <- function(t, w) em$estimate[em$tension == t & em$wool == w]
  for (t in c("L", "M", "H")) for (w in c("A", "B")) {
    expect_equal(est(t, w), m[t, w], tolerance = 1e-6)
  }
  expect_gt(m["L", "A"], m["L", "B"])
  expect_lt(m["M", "A"], m["M", "B"])
  expect_gt(m["H", "A"], m["H", "B"])
  expect_gt(f$posthoc$p_adjusted[f$posthoc$contrast == "H A / H B"], 0.05)
  txt <- gsub("\\s+", " ", pkg_file("vignettes", "articles", "visuals.Rmd"))
  expect_no_match(txt, "no more often at the other two", fixed = TRUE)
  expect_match(txt, "less often at medium tension", fixed = TRUE)
})

test_that("R3-37: printed tables show no digits that are not in the data", {
  d <- withr::with_seed(1, data.frame(
    g = rep(c("a", "b", "c"), each = 10),
    y = stats::rnorm(30) + rep(c(0, 1, 5), each = 10), x = stats::rnorm(30)))
  fit <- anova_ancova(d, "y", "g", covariates = "x", plots = FALSE)
  out <- utils::capture.output(print(fit))
  expect_false(any(grepl("122.40000|23.94000|66.48000", out)))
  # Each printed number is the stored value at 4 significant digits, with no
  # padding (reference: signif() and a hand count of the digits shown).
  tab <- .round_df(fit$anova, 4L)
  for (col in c("sum_sq", "statistic")) {
    shown <- tab[[col]][!is.na(fit$anova[[col]])]
    stored <- fit$anova[[col]][!is.na(fit$anova[[col]])]
    expect_equal(as.numeric(shown), signif(stored, 4))
    sig <- nchar(gsub("^0+", "", gsub("[-.]", "", sub("e.*$", "", shown))))
    expect_true(all(sig <= 4L), info = col)
  }
  s <- utils::capture.output(summary(fit))
  ch <- fit$posthoc$conf_high
  for (v in ch) {
    expect_true(any(grepl(format(signif(v, 4), digits = 4), s, fixed = TRUE)))
  }
  expect_false(any(grepl("-3.66600", s, fixed = TRUE)))
})

test_that("R3-38: Welch and Kruskal-Wallis headers name their statistic", {
  d <- withr::with_seed(1, data.frame(
    g = rep(c("a", "b", "c"), each = 10),
    y = stats::rnorm(30) + rep(c(0, 1, 5), each = 10)))
  w <- anova_welch(d, "y", "g", plots = FALSE)
  k <- anova_kw(d, "y", "g", plots = FALSE)
  # Reference: the htest objects from stats.
  ow <- stats::oneway.test(y ~ g, d, var.equal = FALSE)
  kt <- stats::kruskal.test(y ~ factor(g), d)
  expect_equal(w$anova$statistic, unname(ow$statistic))
  expect_equal(k$anova$statistic, unname(kt$statistic))
  expect_identical(attr(w$anova, "statistic"), "Welch F")
  expect_identical(attr(k$anova, "statistic"), "Kruskal-Wallis chi-square")
  expect_true("Omnibus test (Welch F)" %in% utils::capture.output(print(w)))
  expect_true("Omnibus test (Kruskal-Wallis chi-square)" %in%
                utils::capture.output(print(k)))
})

test_that("R3-39: Scheffe comparisons with an empty cell use the estimable cells' rank", {
  d <- withr::with_seed(5, {
    d <- expand.grid(g = c("a", "b", "c"), h = c("u", "v"), r = 1:6)
    d$y <- stats::rnorm(nrow(d)) + as.integer(d$g)
    d$k <- stats::rpois(nrow(d), 2 + as.integer(d$g))
    d$b <- stats::rbinom(nrow(d), 1, c(.3, .5, .7)[as.integer(d$g)])
    d <- d[!(d$g == "c" & d$h == "v"), ]
    d$cell <- interaction(d$g, d$h, drop = TRUE, sep = " ")
    d
  })
  check <- function(fit, ref_model) {
    ref <- summary(emmeans::contrast(emmeans::emmeans(ref_model, "cell"), "pairwise",
                                     adjust = "scheffe"),
                   infer = TRUE, level = 0.9, adjust = "scheffe")
    expect_identical(nrow(fit$posthoc), 10L)
    expect_equal(fit$posthoc$p_adjusted, ref$p.value, tolerance = 1e-6)
    lo <- intersect(c("lower.CL", "asymp.LCL"), names(ref))
    expect_equal(fit$posthoc$conf_low, ref[[lo]], tolerance = 1e-6)
  }
  f <- anova_glm(d, "y", c("g", "h"), interaction = TRUE, adjust = "scheffe",
                 conf_level = 0.9, plots = FALSE)
  check(f, stats::lm(y ~ cell, d))
  # Hand computation: rank 4 = 5 estimable cells - 1.
  t1 <- f$posthoc$statistic[1]
  expect_equal(f$posthoc$p_adjusted[1],
               stats::pf(t1^2 / 4, 4, 25, lower.tail = FALSE), tolerance = 1e-8)
  ph <- f$posthoc[1, ]
  expect_equal((ph$estimate - ph$conf_low) / ph$se, sqrt(4 * stats::qf(.9, 4, 25)),
               tolerance = 1e-8)
  check(anova_glm(d, "y", c("g", "h"), interaction = TRUE, adjust = "scheffe",
                  conf_level = 0.9, plots = FALSE),
        stats::lm(y ~ cell, d))
  cnt <- anova_count(d, "k", c("g", "h"), interaction = TRUE, model = "poisson",
                     adjust = "scheffe", conf_level = 0.9, plots = FALSE)
  expect_equal(cnt$posthoc$p_adjusted,
               summary(emmeans::contrast(emmeans::emmeans(
                 stats::glm(k ~ cell, stats::poisson(), d), "cell"), "pairwise",
                 adjust = "scheffe"), adjust = "scheffe")$p.value, tolerance = 1e-6)
  bn <- anova_bin(d, "b", c("g", "h"), interaction = TRUE, adjust = "scheffe",
                  conf_level = 0.9, plots = FALSE)
  expect_equal(bn$posthoc$p_adjusted,
               summary(emmeans::contrast(emmeans::emmeans(
                 stats::glm(b ~ cell, stats::binomial(), d), "cell"), "pairwise",
                 adjust = "scheffe"), adjust = "scheffe")$p.value, tolerance = 1e-6)
})

test_that("R3-53: additive ANCOVA and every model function name an extrapolated empty cell", {
  d <- withr::with_seed(5, {
    d <- expand.grid(a = c("a1", "a2", "a3"), b = c("b1", "b2"), rep = 1:10)
    d <- d[!(d$a == "a3" & d$b == "b2"), ]
    d$y <- stats::rnorm(nrow(d)); d$x <- stats::rnorm(nrow(d))
    d$bin <- stats::rbinom(nrow(d), 1, .5); d$cnt <- stats::rpois(nrow(d), 4)
    d
  })
  fits <- list(
    glm = anova_glm(d, "y", c("a", "b"), plots = FALSE),
    bin = anova_bin(d, "bin", c("a", "b"), plots = FALSE),
    count = anova_count(d, "cnt", c("a", "b"), plots = FALSE),
    ancova = anova_ancova(d, "y", c("a", "b"), "x", interaction = FALSE, plots = FALSE),
    manova = anova_manova(d, c("y", "x"), c("a", "b"), plots = FALSE))
  for (nm in names(fits)) {
    expect_true(any(grepl("contain no observations: a3 : b2. The additive model",
                          fits[[nm]]$notes, fixed = TRUE)), info = nm)
  }
  # Reference: a3's additive marginal mean averages over the extrapolated cell,
  # so it differs from the mean of the only observed a3 cell.
  an <- fits$ancova
  ref <- summary(emmeans::emmeans(stats::lm(y ~ a + b + x, d), ~ a:b))
  rep_a3 <- an$emmeans$estimate[an$emmeans$term == "a" & an$emmeans$a == "a3"]
  expect_equal(rep_a3, mean(ref$emmean[ref$a == "a3"]), tolerance = 1e-8)
})

test_that("R3-54: summary() shows anova_bin's model statistics and model-versus-null test", {
  d <- withr::with_seed(7, {
    d <- data.frame(g = rep(c("a", "b", "c"), each = 30))
    d$y <- stats::rbinom(90, 1, c(.2, .5, .8)[as.integer(factor(d$g))])
    d
  })
  b <- anova_bin(d, "y", "g", plots = FALSE)
  # Reference: the likelihood-ratio test of the fit against the null model.
  m1 <- stats::glm(y ~ g, stats::binomial(), d)
  m0 <- stats::glm(y ~ 1, stats::binomial(), d)
  lr <- stats::anova(m0, m1, test = "LRT")
  expect_equal(b$model_stats$lr_vs_null$statistic, lr$Deviance[2])
  out <- utils::capture.output(summary(b))
  i <- grep("Model versus the null model", out, fixed = TRUE)
  expect_length(i, 1L)
  expect_match(out[i + 2L], format(signif(lr$Deviance[2], 4), digits = 4), fixed = TRUE)
  expect_match(out[i + 2L], format.pval(lr$`Pr(>Chi)`[2], digits = 4), fixed = TRUE)
  j <- grep("^Model fit$", out)
  expect_length(j, 1L)
  expect_match(out[j + 2L], format(signif(stats::AIC(m1), 4), digits = 4), fixed = TRUE)
})

test_that("R3-55: aliased coefficients are attributed to the empty cell or the collinear covariate", {
  e <- withr::with_seed(11, {
    e <- expand.grid(A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3"), rep = 1:6)
    e <- e[!(e$A == "a3" & e$B == "b3"), ]
    e$y <- stats::rnorm(nrow(e)); e$k <- stats::rpois(nrow(e), 3)
    e
  })
  # Reference: the empty cell from the data, and the aliased sum-to-zero
  # coefficient, which names a different cell.
  expect_identical(which(table(e$A, e$B) == 0, arr.ind = FALSE), 9L)
  ref <- stats::lm(y ~ A * B, e, contrasts = list(A = "contr.sum", B = "contr.sum"))
  expect_identical(names(stats::coef(ref))[is.na(stats::coef(ref))], "A2:B2")
  g3 <- anova_glm(e, "y", c("A", "B"), interaction = TRUE, type = "III",
                  ci_method = "wald", plots = FALSE)
  al <- grep("could not be estimated", g3$notes, value = TRUE)
  expect_length(al, 1L)
  expect_match(al, "A:B term", fixed = TRUE)
  expect_match(al, "contain no data: a3 : b3", fixed = TRUE)
  expect_no_match(al, "A2:B2", fixed = TRUE)
  expect_no_match(al, "collinear")

  cnt <- anova_count(e, "k", c("A", "B"), interaction = TRUE, plots = FALSE)
  expect_true(any(grepl("1 of the 9 possible group combination(s) contain no observations; contrasts involving them are not estimable. They are: a3 : b3.",
                        cnt$notes, fixed = TRUE)))

  # A covariate constant within one group: collinearity, not an empty cell.
  d <- withr::with_seed(10, {
    d <- data.frame(g = rep(c("a", "b", "c"), each = 10), x = stats::rnorm(30))
    d$y <- d$x + stats::rnorm(30)
    d$x3 <- ifelse(d$g == "b", 1, stats::rnorm(30))
    d
  })
  expect_true(is.na(stats::coef(stats::lm(y ~ x3 * g, d))["x3:gb"]))
  a <- anova_ancova(d, "y", "g", "x3", force_interaction = TRUE, plots = FALSE)
  al <- grep("could not be estimated", a$notes, value = TRUE)
  expect_length(al, 1L)
  expect_match(al, "collinear with other terms", fixed = TRUE)
  expect_match(al, "covariate(s) `x3`", fixed = TRUE)
  expect_no_match(al, "contain no data")
})
