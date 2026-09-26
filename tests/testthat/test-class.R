# The shared return class and its methods.
#
# The package's central promise is that all eight functions return the same
# shape, so these blocks run every one of them and check the components, the
# print/summary/plot methods, and the invariants that hold across all of them:
# nothing printed, no plot drawn, and the row accounting closing.

test_that("every function returns the same class with the same components", {
  d_num <- fx_oneway()
  d_bin <- fx_binary()
  d_cnt <- fx_counts()
  d_anc <- fx_ancova()
  d_mv  <- fx_multivariate()

  fits <- list(
    welch  = anova_welch(d_num, "value", "group", plots = FALSE),
    kw     = anova_kw(d_num, "value", "group", plots = FALSE),
    glm    = anova_glm(d_num, "value", "group", plots = FALSE),
    bin    = anova_bin(d_bin, "y", "g", plots = FALSE),
    count  = anova_count(d_cnt, "count", "g1", plots = FALSE),
    ancova = anova_ancova(d_anc, "dv", "iv", "cov", plots = FALSE),
    manova = anova_manova(d_mv, c("score1", "score2"), "g",
                          assumptions = FALSE, plots = FALSE)
  )
  if (requireNamespace("afex", quietly = TRUE)) {
    fits$rm <- anova_rm(fx_repeated(), "score", "id", "time", plots = FALSE)
  }
  for (nm in names(fits)) {
    expect_fit_shape(fits[[nm]])
    expect_true(is.data.frame(fits[[nm]]$anova), info = nm)
    expect_identical(fits[[nm]]$conf_level, 0.95, info = nm)
  }
})

test_that("print returns its input invisibly and mentions the method", {
  fit <- anova_welch(fx_oneway(), "value", "group", plots = FALSE)
  out <- utils::capture.output(res <- print(fit))
  expect_identical(res, fit)
  expect_true(any(grepl("Welch", out)))
  expect_true(any(grepl("Observations used", out)))
})

test_that("summary prints the sections the fit actually has", {
  fit <- anova_welch(fx_oneway(), "value", "group", plots = FALSE)
  res <- summary(fit)
  expect_s3_class(res, "summary.anovakit_fit")
  expect_identical(res$fit, fit)
  out <- utils::capture.output(print(res))
  expect_true(any(grepl("Assumption checks", out)))
  expect_true(any(grepl("Effect sizes", out)))
  # Welch fits no model, so its table is of group summaries, not EMMs
  expect_true(any(grepl("Group summaries", out)))
  expect_false(any(grepl("Estimated marginal means", out)))
  expect_true(any(grepl("Pairwise comparisons", out)))
})

test_that("plot selects by name and by index and rejects bad input", {
  fit <- anova_welch(fx_oneway(), "value", "group")
  expect_s3_class(plot(fit), "ggplot")
  expect_s3_class(plot(fit, "box"), "ggplot")
  expect_s3_class(plot(fit, 3L), "ggplot")
  expect_error(plot(fit, "not_a_plot"), "No plot called")
  expect_error(plot(fit, 99L), "whole number between 1 and")
  expect_error(plot(fit, NULL), "single plot name or index")
  expect_error(plot(fit, c(1, 2)), "single plot name or index")
})

test_that("plot on a fit built without plots explains itself", {
  fit <- anova_welch(fx_oneway(), "value", "group", plots = FALSE)
  expect_message(res <- plot(fit), "holds no plots")
  expect_null(res)
})

test_that("the notes vector is deduplicated and always character", {
  d <- fx_oneway()
  d$value[1] <- NA
  fit <- anova_welch(d, "value", "group", plots = FALSE)
  expect_type(fit$notes, "character")
  expect_identical(fit$notes, unique(fit$notes))
})

test_that("every function accepts a data.table and a tibble", {
  calls <- list(
    welch  = list(fx_oneway(), function(d) anova_welch(d, "value", "group", plots = FALSE)),
    kw     = list(fx_oneway(), function(d) anova_kw(d, "value", "group", plots = FALSE)),
    glm    = list(fx_oneway(), function(d) anova_glm(d, "value", "group", plots = FALSE)),
    bin    = list(fx_binary(), function(d) anova_bin(d, "y", "g", plots = FALSE)),
    count  = list(fx_counts(), function(d) anova_count(d, "count", "g1", plots = FALSE)),
    ancova = list(fx_ancova(), function(d) anova_ancova(d, "dv", "iv", "cov", plots = FALSE)),
    manova = list(fx_multivariate(), function(d) anova_manova(
      d, c("score1", "score2"), "g", assumptions = FALSE, plots = FALSE))
  )
  if (requireNamespace("afex", quietly = TRUE)) {
    calls$rm <- list(fx_repeated(), function(d) anova_rm(d, "score", "id", "time",
                                                         plots = FALSE))
  }
  for (nm in names(calls)) {
    d <- calls[[nm]][[1]]
    f <- calls[[nm]][[2]]
    ref <- f(d)
    if (requireNamespace("data.table", quietly = TRUE)) {
      got <- f(data.table::as.data.table(d))
      expect_equal(got$anova, ref$anova, info = paste(nm, "data.table"))
      expect_identical(got$n_removed, ref$n_removed, info = nm)
    }
    if (requireNamespace("tibble", quietly = TRUE)) {
      got <- f(tibble::as_tibble(d))
      expect_equal(got$anova, ref$anova, info = paste(nm, "tibble"))
      expect_identical(nrow(got$data_used) + got$n_removed, nrow(d), info = nm)
    }
  }
})

test_that("verbose = TRUE uses message(), never cat()", {
  d <- fx_oneway()
  expect_silent_stdout(suppressMessages(
    anova_welch(d, "value", "group", verbose = TRUE)))
  msgs <- testthat::capture_messages(
    anova_welch(d, "value", "group", verbose = TRUE))
  expect_gt(length(msgs), 1L)
  expect_match(msgs[1], "Preparing data")
  expect_silent_stdout(
    suppressMessages(anova_kw(d, "value", "group", verbose = TRUE)))
})

test_that("the shared validators reject the same things everywhere", {
  d <- fx_oneway()
  callers <- list(
    function(dd, ...) anova_welch(dd, "value", "group", ...),
    function(dd, ...) anova_kw(dd, "value", "group", ...),
    function(dd, ...) anova_glm(dd, "value", "group", ...)
  )
  for (f in callers) {
    expect_error(f(d, conf_level = 0), "strictly between 0 and 1")
    expect_error(f(d, conf_level = c(0.9, 0.95)), "strictly between 0 and 1")
    expect_error(f(d, plots = "yes"), "must be TRUE or FALSE")
    expect_error(f(data.frame()), "no rows")
    expect_error(f(list(a = 1)), "must be a data.frame")
  }
})
