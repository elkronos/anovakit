test_that("an ordered first factor's marginal means are joined by a line", {
  has_line <- function(p) any(vapply(p$layers, function(l) inherits(l$geom, "GeomLine"), logical(1)))
  w <- warpbreaks
  f_un <- anova_count(w, "breaks", "tension", plots = TRUE)
  expect_false(has_line(f_un$plots$emmeans))
  w$tension <- factor(w$tension, ordered = TRUE)
  f_ord <- anova_count(w, "breaks", "tension", plots = TRUE)
  expect_true(has_line(f_ord$plots$emmeans))
  # Only the plot changes.
  expect_equal(f_ord$anova, f_un$anova)

  g <- anova_glm(w, "breaks", "tension", family = "poisson", plots = TRUE)
  expect_true(has_line(g$plots$emmeans))
  m <- anova_ancova(transform(mtcars, gear = factor(gear, ordered = TRUE)),
                    "mpg", "gear", "wt", plots = TRUE)
  expect_true(has_line(m$plots$emmeans))
})

test_that("the ANCOVA covariate plot's subtitle has one clause per line", {
  f <- anova_ancova(mtcars, "mpg", "cyl", "wt", plots = TRUE)
  sub <- f$plots$covariate$labels$subtitle
  expect_match(sub, ";\n", fixed = TRUE)
  expect_true(all(nchar(strsplit(sub, "\n", fixed = TRUE)[[1]]) < 110))
})
