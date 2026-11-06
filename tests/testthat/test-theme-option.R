# The anovakit.theme option sets the base theme of every plot.

theme_el <- function(p, el) ggplot2::calc_element(el, ggplot2::ggplot_build(p)$plot$theme)

test_that("unset, plots use theme_minimal with the package's adjustments", {
  withr::local_options(anovakit.theme = NULL)
  f <- anova_welch(PlantGrowth, "weight", "group")
  p <- plot(f, "means")
  expect_s3_class(theme_el(p, "panel.border"), "element_blank")   # theme_minimal
  expect_identical(theme_el(p, "plot.title")$face, "bold")
  expect_identical(theme_el(p, "legend.position"), "bottom")
})

test_that("a theme set in the option is the base of every plot, adjustments on top", {
  withr::local_options(anovakit.theme = ggplot2::theme_classic(base_size = 17))
  f <- anova_glm(chickwts, "weight", "feed")
  for (nm in names(f$plots)) {
    p <- f$plots[[nm]]
    expect_s3_class(theme_el(p, "axis.line"), "element_line")       # theme_classic
    expect_equal(theme_el(p, "text")$size, 17, info = nm)
    expect_identical(theme_el(p, "plot.title")$face, "bold", info = nm)
  }
  # long feed names still get rotated labels on the classic base
  expect_equal(theme_el(f$plots[["emmeans"]], "axis.text.x")$angle, 45)
  # MANOVA's own univariate plot uses the option too
  m <- anova_manova(iris, c("Sepal.Length", "Petal.Length"), "Species")
  expect_true(all(vapply(m$plots, function(p)
    inherits(theme_el(p, "axis.line"), "element_line"), logical(1))))
})

test_that("a function in the option is called when a plot is built", {
  withr::local_options(anovakit.theme = ggplot2::theme_get)
  old <- ggplot2::theme_set(ggplot2::theme_bw(base_size = 9))
  on.exit(ggplot2::theme_set(old), add = TRUE)
  p <- plot(anova_kw(PlantGrowth, "weight", "group"))
  expect_equal(theme_el(p, "text")$size, 9)
  expect_s3_class(theme_el(p, "panel.border"), "element_rect")      # theme_bw
})

test_that("the option changes no result", {
  withr::local_options(anovakit.theme = NULL)
  a <- anova_ancova(MASS::anorexia, "Postwt", "Treat", "Prewt")
  withr::local_options(anovakit.theme = ggplot2::theme_dark())
  b <- anova_ancova(MASS::anorexia, "Postwt", "Treat", "Prewt")
  for (tbl in c("anova", "effect_sizes", "emmeans", "posthoc", "notes")) {
    expect_identical(b[[tbl]], a[[tbl]], info = tbl)
  }
})

test_that("an invalid option is refused with a clear message", {
  withr::local_options(anovakit.theme = "classic")
  expect_error(anova_welch(PlantGrowth, "weight", "group"),
               "must be a ggplot2 theme.*class character")
  withr::local_options(anovakit.theme = function() 1)
  expect_error(anova_welch(PlantGrowth, "weight", "group"),
               "function returning an object of class numeric")
  # plots = FALSE never consults it
  expect_s3_class(anova_welch(PlantGrowth, "weight", "group", plots = FALSE),
                  "anovakit_fit")
})
