# Regression tests for problems found while writing the website's walkthroughs.

test_that("the canonical plot's subtitle describes what is drawn", {
  skip_if_not_installed("MASS")
  sv <- MASS::survey
  mc <- anova_manova(sv, c("Wr.Hnd", "NW.Hnd"), "Sex", covariates = "Height",
                     posthoc = FALSE)
  p <- mc$plots$canonical
  expect_true(any(vapply(p$layers, function(l) inherits(l$geom, "GeomBoxplot"),
                         logical(1))))
  expect_no_match(p$labels$subtitle, "Crosses")
  expect_match(p$labels$subtitle, "Only one discriminant axis")
  expect_match(p$labels$subtitle, "`Sex` term")

  iris_fit <- anova_manova(iris, c("Sepal.Length", "Sepal.Width",
                                   "Petal.Length", "Petal.Width"), "Species",
                           posthoc = FALSE)
  expect_identical(iris_fit$plots$canonical$labels$subtitle,
                   "Crosses mark group centroids")
})

test_that("anova_rm records no correction when there is nothing to correct", {
  skip_if_not_installed("afex")
  fit <- anova_rm(datasets::sleep, "extra", subject = "ID", within = "group",
                  plots = FALSE)
  expect_null(fit$sphericity)
  expect_identical(attr(fit$anova, "correction"), "none")
  expect_match(.omnibus_header(fit$anova), "no sphericity correction")
  expect_equal(fit$anova$num_df, 1)
})

test_that("MANOVA's assumption notes do not recommend the statistic in use", {
  fit <- anova_manova(iris, c("Sepal.Length", "Sepal.Width", "Petal.Length",
                              "Petal.Width"), "Species", posthoc = FALSE,
                      plots = FALSE)
  box <- grep("Box's M rejects", fit$notes, value = TRUE)
  expect_length(box, 1L)
  expect_match(box, "equal group sizes")
  expect_match(box, "it is the statistic used here")
  expect_false(any(grepl("consider test", fit$notes)))
  wilks <- anova_manova(iris, c("Sepal.Length", "Sepal.Width", "Petal.Length",
                                "Petal.Width"), "Species", test = "Wilks",
                        posthoc = FALSE, plots = FALSE)
  expect_true(any(grepl("consider test = \"Pillai\"", wilks$notes, fixed = TRUE)))
})

test_that("the ANCOVA imbalance note reports the gap to two decimals", {
  skip_if_not_installed("MASS")
  bw <- MASS::birthwt
  bw$race <- factor(bw$race, labels = c("white", "black", "other"))
  fit <- anova_ancova(bw, "bwt", "race", "lwt", force_interaction = FALSE,
                      plots = FALSE)
  note <- grep("covariate means differ", fit$notes, value = TRUE)
  expect_length(note, 1L)
  expect_match(note, "lwt [0-9]+\\.[0-9]{2} pooled")
})

test_that("a single comparison under a step-down adjustment gets no mismatch note", {
  ucb <- as.data.frame(UCBAdmissions)
  fit <- anova_bin(ucb, "Admit", c("Gender", "Dept"), weights = "Freq",
                   adjust = "holm", plots = FALSE)
  mismatch <- grep("comparison intervals use the", fit$notes, value = TRUE)
  expect_false(any(grepl("none adjustment", mismatch)))
  one <- fit$posthoc[fit$posthoc$term == "Gender", ]
  expect_equal(one$p_adjusted, one$p_value)
})

test_that("anova_count counts observations, not rows, under frequency weights", {
  tab <- as.data.frame(table(spray = InsectSprays$spray,
                             count = InsectSprays$count))
  tab$count <- as.integer(as.character(tab$count))
  tab <- tab[tab$Freq > 0, ]
  fit <- anova_count(tab, "count", "spray", weights = "Freq")
  cc <- fit$assumptions$cell_counts
  expect_equal(cc$n, rep(12, 6))
  expect_identical(cc$rows, as.integer(table(droplevels(tab$spray))))
  ld <- ggplot2::layer_data(fit$plots$observed, 1)
  ref <- tapply(InsectSprays$count, InsectSprays$spray, stats::median)
  expect_equal(sort(ld$middle), sort(as.numeric(ref)))
})

test_that("with an offset the observed plot shows rates", {
  set.seed(41)
  d <- data.frame(g = gl(3, 20), hours = runif(60, 1, 10))
  d$y <- rpois(60, d$hours * c(1, 2, 3)[d$g])
  fit <- anova_count(d, "y", "g", offset = "hours", model = "poisson")
  p <- fit$plots$observed
  expect_match(p$labels$title, "Observed rates")
  ld <- ggplot2::layer_data(p, 1)
  ref <- tapply(d$y / d$hours, d$g, stats::median)
  expect_equal(sort(ld$middle), sort(as.numeric(ref)))
})

test_that("print() shows a long call without deparse's indentation", {
  fit <- anova_welch(warpbreaks, "breaks", c("wool", "tension"), plots = FALSE)
  out <- utils::capture.output(print(fit))
  call_line <- grep("^Call: ", out, value = TRUE)
  expect_length(call_line, 1L)
  expect_no_match(call_line, "  ")
})

test_that("the box plot draws infinite responses without warning", {
  worst <- data.frame(g = rep(c("a", "b"), each = 5),
                      y = c(1:5, 6, 7, Inf, Inf, Inf))
  fit <- anova_kw(worst, "y", "g")
  expect_no_warning(ggplot2::ggplot_build(fit$plots$box))
  ld <- ggplot2::layer_data(fit$plots$box, 1)
  expect_true(all(is.finite(ld$ymax)))
})

test_that("anova_glm's table has one column order whatever the statistic", {
  skip_if_not_installed("MASS")
  cols <- lapply(c("F", "LR", "Wald"), function(ts) {
    names(anova_glm(MASS::Cars93, "Price", "Type", family = Gamma(link = "log"),
                    test_statistic = ts, plots = FALSE)$anova)
  })
  expect_identical(cols[[1]], c("term", "sum_sq", "df", "statistic", "p_value"))
  expect_identical(cols[[2]], c("term", "df", "statistic", "p_value"))
  expect_identical(cols[[3]], cols[[2]])
})

test_that("marginal-means plots join points only across a second factor", {
  one <- anova_glm(chickwts, "weight", "feed")
  geoms <- vapply(one$plots$emmeans$layers, function(l) class(l$geom)[1],
                  character(1))
  expect_false("GeomLine" %in% geoms)
  tg <- transform(ToothGrowth, dose = factor(dose))
  two <- anova_glm(tg, "len", c("dose", "supp"), interaction = TRUE)
  geoms2 <- vapply(two$plots$emmeans$layers, function(l) class(l$geom)[1],
                   character(1))
  expect_true("GeomLine" %in% geoms2)
})

test_that("a Type III model with an empty cell falls back to Type II tests", {
  mt <- transform(mtcars, cyl = factor(cyl), gear = factor(gear))
  a <- anova_ancova(mt, "mpg", c("cyl", "gear"), "wt", plots = FALSE)
  expect_false(is.null(a$anova))
  expect_identical(attr(a$anova, "ss_type"), "II")
  expect_match(a$method, "Type II")
  expect_true(any(grepl("Type II tests were computed instead", a$notes)))
  g <- anova_glm(mt, "mpg", c("cyl", "gear"), interaction = TRUE, type = "III",
                 plots = FALSE)
  ref <- suppressMessages(car::Anova(stats::lm(mpg ~ cyl * gear, data = mt), type = 2))
  expect_equal(g$anova$statistic[1:3], unname(ref[["F value"]][1:3]))
})

test_that("a rate ratio against a cell with no events has a one-sided profile interval", {
  set.seed(10)
  traps <- data.frame(site = factor(rep(c("A", "B", "C"), each = 10)))
  traps$n <- rpois(30, c(A = 4, B = 2, C = 0.1)[as.character(traps$site)])
  expect_identical(sum(traps$n[traps$site == "C"]), 0L)
  fit <- anova_count(traps, "n", "site", model = "poisson", plots = FALSE)
  es <- fit$effect_sizes
  ref <- anova_glm(traps, "n", "site", family = "poisson", plots = FALSE)
  ref_ct <- ref$assumptions$coefficients
  expect_identical(es$conf_low[es$term == "siteC"], 0)
  expect_equal(es$conf_high[es$term == "siteC"],
               exp(ref_ct$conf_high[ref_ct$term == "siteC"]), tolerance = 1e-5)
  expect_true(is.finite(es$conf_high[es$term == "siteC"]))
  expect_true(any(grepl("profiling the likelihood directly", fit$notes)))
})

test_that("Shapiro-Wilk is skipped for a group with fewer than three distinct values", {
  w <- anova_welch(mtcars, "mpg", c("cyl", "am"), plots = FALSE)
  row <- w$assumptions$normality[w$assumptions$normality$group == "6 : 1", ]
  expect_true(is.na(row$p_value))
  expect_match(row$note, "distinct value")
})

test_that("anova_kw's diagnostics name the groups they could not test", {
  k <- anova_kw(mtcars, "mpg", c("cyl", "am"), diagnostics = TRUE, plots = FALSE)
  expect_true(any(grepl("Shapiro-Wilk was not run for group \"8 : 1\"", k$notes,
                        fixed = TRUE)))
})
