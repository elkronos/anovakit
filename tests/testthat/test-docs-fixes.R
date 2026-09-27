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
