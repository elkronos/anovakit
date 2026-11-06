# Third review, input handling: the note on dropped rows, integer64 and
# matrix columns, groups lost to zero weights, compact calls, and a character
# binomial response with an unknown encoding. Each test checks the result
# against an independent reference.

test_that("R3-20: the dropped-rows note names only the columns that held missing or infinite values", {
  d <- withr::with_seed(1, data.frame(g = rep(c("a", "b", "c"), each = 8),
                                      w = 1, y = stats::rnorm(24)))
  d$y[c(2, 9)] <- NA
  # Reference: the columns with a missing value, by base R.
  na_cols <- names(d)[colSums(is.na(d)) > 0]
  expect_identical(na_cols, "y")

  fw <- anova_welch(d, "y", "g", plots = FALSE)
  expect_true("Dropped 2 row(s) with missing values in: y." %in% fw$notes)
  fg <- anova_glm(d, "y", "g", weights = "w", plots = FALSE, ci_method = "wald")
  expect_true("Dropped 2 row(s) with missing values in: y." %in% fg$notes)
  expect_identical(fg$n_removed, sum(!stats::complete.cases(d)))

  # Only an infinite value: the note does not say "missing".
  d2 <- d
  d2$y[2] <- Inf
  d2$y[9] <- 1
  expect_identical(names(d2)[vapply(d2, function(v) any(is.infinite(v)), NA)], "y")
  f2 <- anova_glm(d2, "y", "g", plots = FALSE, ci_method = "wald")
  expect_true("Dropped 1 row(s) with infinite values in: y." %in% f2$notes)
  expect_false(any(grepl("missing", f2$notes)))

  # Both kinds, in different columns, listed in the order given.
  d3 <- d2
  d3$g[5] <- NA
  f3 <- anova_glm(d3, "y", "g", weights = "w", plots = FALSE, ci_method = "wald")
  expect_true("Dropped 2 row(s) with missing or infinite values in: y, g." %in% f3$notes)

  # An infinite value in a row already missing something else is still named.
  d5 <- d
  d5$w[2] <- Inf
  f5 <- anova_glm(d5, "y", "g", weights = "w", plots = FALSE, ci_method = "wald")
  expect_true("Dropped 2 row(s) with missing or infinite values in: y, w." %in% f5$notes)

  # NA as an explicit factor level counts as missing and names its column.
  d4 <- data.frame(g = factor(rep(c("a", "b", "c"), each = 8)), x = 1,
                   y = withr::with_seed(2, stats::rnorm(24)))
  d4$g[3] <- NA
  d4$g <- addNA(d4$g)
  f4 <- anova_welch(d4, "y", "g", plots = FALSE)
  expect_true("Dropped 1 row(s) with missing values in: g." %in% f4$notes)
})

test_that("R3-30: integer64 response, weights, covariate and offset are read as numbers", {
  skip_if_not_installed("bit64")
  d <- withr::with_seed(1, data.frame(g = rep(c("a", "b", "c"), each = 12),
                                      y = stats::rnorm(36), x = stats::rnorm(36),
                                      w = rep(1:3, 12)))
  d$y <- round((d$y + c(a = 0, b = 0.5, c = 1)[d$g]) * 100)
  d$x <- round(d$x * 100)
  D <- d
  D$y <- bit64::as.integer64(d$y)
  expect_true(is.numeric(D$y))

  # Response, against stats::anova(lm()) on the same data as double.
  ref <- stats::anova(stats::lm(y ~ g, d))
  f <- anova_glm(D, "y", "g", plots = FALSE, ci_method = "wald")
  expect_equal(f$anova$sum_sq, ref[["Sum Sq"]], tolerance = 1e-8)
  expect_equal(f$anova$p_value[1], ref[["Pr(>F)"]][1], tolerance = 1e-8)
  expect_identical(f$n_removed, 0L)
  expect_type(f$data_used$y, "double")

  # Weights.
  refw <- stats::anova(stats::lm(y ~ g, d, weights = w))
  fw <- anova_glm(transform(d, w = bit64::as.integer64(w)), "y", "g",
                  weights = "w", plots = FALSE, ci_method = "wald")
  expect_equal(fw$anova$sum_sq, refw[["Sum Sq"]], tolerance = 1e-8)
  expect_equal(fw$anova$statistic[1], refw[["F value"]][1], tolerance = 1e-8)

  # Covariate: the integer64 run equals the double run, which equals car.
  fa <- anova_ancova(transform(d, x = bit64::as.integer64(x)), "y", "g", "x",
                     plots = FALSE)
  fa_ref <- anova_ancova(d, "y", "g", "x", plots = FALSE)
  expect_equal(fa$anova, fa_ref$anova)
  car_ref <- car::Anova(stats::lm(y ~ x + g, d), type = 2)
  expect_equal(fa$anova$sum_sq[fa$anova$term == "g"], car_ref["g", "Sum Sq"],
               tolerance = 1e-8)

  # Welch and Kruskal-Wallis responses.
  fwe <- anova_welch(D, "y", "g", plots = FALSE)
  expect_equal(fwe$anova$statistic,
               unname(stats::oneway.test(y ~ g, d)$statistic), tolerance = 1e-8)
  fk <- anova_kw(D, "y", "g", plots = FALSE)
  expect_equal(fk$anova$statistic,
               unname(stats::kruskal.test(y ~ factor(g), d)$statistic),
               tolerance = 1e-8)

  # Count response and offset, against glm().
  dc <- withr::with_seed(4, data.frame(g = rep(c("a", "b", "c"), each = 10),
                                       e = rep(1:5, 6)))
  dc$y <- withr::with_seed(5, stats::rpois(30, 2 * dc$e))
  Dc <- transform(dc, y = bit64::as.integer64(y), e = bit64::as.integer64(e))
  fc <- anova_count(Dc, "y", "g", offset = "e", plots = FALSE)
  fc_ref <- anova_count(dc, "y", "g", offset = "e", plots = FALSE)
  expect_equal(fc$anova, fc_ref$anova)
  gref <- stats::glm(y ~ g + offset(log(e)), stats::poisson(), dc)
  lrt <- stats::anova(gref, test = "Chisq")
  expect_equal(fc_ref$anova$statistic[1], lrt$Deviance[2], tolerance = 1e-6)

  # A value a double cannot hold exactly is refused, without a console warning.
  big <- D
  big$y[1] <- bit64::as.integer64("9007199254740993")
  expect_error(anova_welch(big, "y", "g", plots = FALSE),
               "integer64 and holds values of 2\\^53 or more")
  expect_no_warning(try(anova_welch(big, "y", "g", plots = FALSE), silent = TRUE))
})

test_that("R3-31: a group whose weights are all zero is named and dropped from $data_used", {
  d <- withr::with_seed(2, data.frame(
    g = factor(rep(c("ctrl", "low", "high"), each = 10),
               levels = c("ctrl", "low", "high")),
    y = stats::rnorm(30), w = 1))
  d$w[d$g == "high"] <- 0
  lost <- "Level \"high\" of `g` had no rows with a non-zero weight in `w`, so it was dropped from the analysis."

  # Reference: the same analysis on the rows with a positive weight.
  sub <- droplevels(d[d$w > 0, ])
  ref <- stats::anova(stats::lm(y ~ g, sub))
  f <- anova_glm(d, "y", "g", weights = "w", plots = FALSE, ci_method = "wald")
  expect_true(lost %in% f$notes)
  expect_identical(levels(f$data_used$g), levels(sub$g))
  expect_equal(f$anova$sum_sq, ref[["Sum Sq"]], tolerance = 1e-8)
  expect_identical(f$n_removed, 10L)

  db <- d
  db$y <- withr::with_seed(3, stats::rbinom(30, 1, 0.5))
  fb <- anova_bin(db, "y", "g", weights = "w", plots = FALSE)
  expect_true(lost %in% fb$notes)
  expect_identical(levels(fb$data_used$g), c("ctrl", "low"))
  bref <- stats::anova(stats::glm(y ~ g, stats::binomial(),
                                  droplevels(db[db$w > 0, ])), test = "Chisq")
  expect_equal(fb$anova$statistic[1], bref$Deviance[2], tolerance = 1e-6)
  expect_equal(fb$anova$df[1], bref$Df[2])
  expect_equal(fb$anova$p_value[1], bref[["Pr(>Chi)"]][2], tolerance = 1e-6)

  dc <- d
  dc$y <- withr::with_seed(3, stats::rpois(30, 3))
  fc <- anova_count(dc, "y", "g", weights = "w", plots = FALSE)
  expect_true(lost %in% fc$notes)
  expect_identical(levels(fc$data_used$g), c("ctrl", "low"))
  gref <- stats::glm(y ~ g, stats::poisson(), droplevels(dc[dc$w > 0, ]))
  expect_equal(fc$anova$statistic[1],
               stats::anova(gref, test = "Chisq")$Deviance[2], tolerance = 1e-6)

  # With a positive weight left in every group, no level note.
  d2 <- d
  d2$w[d2$g == "high"][1] <- 1
  f2 <- anova_glm(d2, "y", "g", weights = "w", plots = FALSE, ci_method = "wald")
  expect_false(any(grepl("^Level", f2$notes)))
})

test_that("R3-32: a two-column matrix response or covariate is refused clearly in every function", {
  d <- withr::with_seed(3, data.frame(g = rep(c("a", "b", "c"), each = 10),
                                      x = stats::rnorm(30), id = rep(1:10, 3),
                                      t = rep(c("t1", "t2", "t3"), each = 10)))
  d$y <- withr::with_seed(4, cbind(stats::rnorm(30), stats::rnorm(30)))
  msg <- "Response `y` is a 2-column matrix; supply a single column."
  expect_error(anova_welch(d, "y", "g", plots = FALSE), msg, fixed = TRUE)
  expect_error(anova_kw(d, "y", "g", plots = FALSE), msg, fixed = TRUE)
  expect_error(anova_count(d, "y", "g", plots = FALSE), msg, fixed = TRUE)
  expect_error(anova_bin(d, "y", "g", plots = FALSE), msg, fixed = TRUE)
  expect_error(anova_ancova(d, "y", "g", "x", plots = FALSE), msg, fixed = TRUE)
  expect_error(anova_manova(d, c("y", "x"), "g", plots = FALSE), msg, fixed = TRUE)
  if (requireNamespace("afex", quietly = TRUE)) {
    expect_error(anova_rm(d, "y", "id", "t", plots = FALSE), msg, fixed = TRUE)
  }
  expect_error(anova_glm(d, "y", "g", plots = FALSE), "2-column matrix")

  d$y1 <- d$y[, 1]
  d$x2 <- cbind(d$x, d$x)
  expect_error(anova_ancova(d, "y1", "g", "x2", plots = FALSE),
               "Covariate `x2` is a 2-column matrix; supply a single column.",
               fixed = TRUE)
  d$w2 <- cbind(rep(1, 30), rep(1, 30))
  expect_error(anova_glm(d, "y1", "g", weights = "w2", plots = FALSE),
               "Weights column `w2` is a 2-column matrix", fixed = TRUE)

  # A one-column matrix (as scale() returns) is used as the vector it holds.
  d$ys <- scale(d$y1)
  fw <- anova_welch(d, "ys", "g", plots = FALSE)
  expect_equal(fw$anova$statistic,
               unname(stats::oneway.test(as.vector(ys) ~ g, d)$statistic))
  expect_null(dim(fw$data_used$ys))
})

test_that("R3-45: every function keeps a compact call when called through do.call()", {
  d <- withr::with_seed(5, data.frame(g = rep(c("a", "b", "c"), 100),
                                      y = stats::rpois(300, 3),
                                      x = stats::rnorm(300),
                                      id = rep(1:100, each = 3),
                                      t = rep(c("t1", "t2", "t3"), 100)))
  d$b <- as.integer(d$y > 3)
  calls <- list(
    anova_welch  = list(response = "y", groups = "g"),
    anova_kw     = list(response = "y", groups = "g"),
    anova_glm    = list(response = "y", groups = "g"),
    anova_count  = list(response = "y", groups = "g"),
    anova_bin    = list(response = "b", groups = "g"),
    anova_ancova = list(response = "x", groups = "g", covariates = "y"),
    anova_manova = list(responses = c("x", "y"), groups = "g"),
    anova_rm     = list(response = "x", subject = "id", within = "t")
  )
  if (!requireNamespace("afex", quietly = TRUE)) calls$anova_rm <- NULL
  for (fn in names(calls)) {
    args <- c(list(data = d), calls[[fn]], list(plots = FALSE))
    f <- do.call(getExportedValue("anovakit", fn), args)
    # Reference: the call as written by hand, with the data named by a symbol.
    expect_identical(f$call[[1]], as.name(fn), info = fn)
    expect_identical(f$call$data, as.name("<data frame: 300 rows>"), info = fn)
    expect_lt(as.numeric(utils::object.size(f$call)),
              as.numeric(utils::object.size(d)) / 5)
    out <- utils::capture.output(print(f))
    call_line <- out[startsWith(out, "Call: ")]
    expect_length(call_line, 1L)
    expect_lt(nchar(call_line), 300L)
  }
  # A call written by hand still shows the name of the data.
  f <- anova_welch(d, "y", "g", plots = FALSE)
  expect_identical(f$call, quote(anova_welch(data = d, response = "y",
                                             groups = "g", plots = FALSE)))
})

test_that("R3-56: a non-ASCII character binomial response with an unknown encoding fits", {
  skip_if_not(isTRUE(l10n_info()[["UTF-8"]]), "needs a UTF-8 session")
  r <- rep(c("\u00e4", "o", "o", "\u00e4", "o"), 6)
  Encoding(r) <- "unknown"
  # The radix sort the old check used refuses these strings.
  expect_error(sort(unique(r), method = "radix"), "encoding")
  d <- data.frame(g = rep(c("a", "b", "c"), each = 10), r = r,
                  stringsAsFactors = FALSE)
  # Row 1 is non-ASCII (the order in which the old check failed), row 2 is
  # ASCII (the order in which it happened to work).
  d$r[c(3, 11, 13, 21, 22, 23)] <- "o"
  Encoding(d$r) <- "unknown"
  expect_identical(d$r[1:2], c(r[1], "o"))
  ref <- stats::anova(stats::glm(I(enc2utf8(r) == "o") ~ g, stats::binomial(), d),
                      test = "Chisq")
  for (fam in list(stats::binomial(), stats::quasibinomial())) {
    for (ord in list(seq_len(30), c(2:30, 1))) {
      f <- anova_glm(d[ord, ], "r", "g", family = fam, plots = FALSE)
      expect_equal(f$data_used$r == "o", d$r[ord] == "o")
      if (fam$family == "binomial") {
        expect_equal(f$anova$p_value[1], ref[["Pr(>Chi)"]][2], tolerance = 1e-6)
      }
    }
  }
})

test_that("R3-30: integer64 columns in anova_rm and anova_manova match the double run and base R", {
  skip_if_not_installed("bit64")
  skip_if_not_installed("afex")
  d <- withr::with_seed(1, data.frame(id = factor(rep(1:10, 3)),
                                      t = factor(rep(c("t1", "t2", "t3"), each = 10)),
                                      g = rep(c("a", "b"), 15),
                                      y = round(stats::rnorm(30) * 100),
                                      x = round(stats::rnorm(30) * 100)))
  D <- transform(d, y = bit64::as.integer64(y), x = bit64::as.integer64(x))

  r64 <- anova_rm(D, "y", "id", "t", plots = FALSE)
  expect_equal(r64$anova, anova_rm(d, "y", "id", "t", plots = FALSE)$anova)
  # The F statistic (unchanged by a sphericity correction) against aov().
  aref <- summary(stats::aov(y ~ t + Error(id / t), d))[["Error: id:t"]][[1]]
  expect_equal(r64$anova$statistic[1], aref[["F value"]][1], tolerance = 1e-8)

  m64 <- anova_manova(D, c("y", "x"), "g", plots = FALSE)
  expect_equal(m64$anova, anova_manova(d, c("y", "x"), "g", plots = FALSE)$anova)
  mref <- summary(stats::manova(cbind(y, x) ~ g, d), test = "Pillai")$stats
  expect_equal(m64$anova$statistic[1], mref["g", "Pillai"], tolerance = 1e-8)
  expect_equal(m64$anova$p_value[1], mref["g", "Pr(>F)"], tolerance = 1e-8)
})

test_that("anova_glm unwraps a one-column data-frame response like the other functions", {
  ins <- datasets::InsectSprays
  ins$count <- data.frame(z = ins$count)
  f <- anova_glm(ins, "count", "spray", family = stats::poisson(), plots = FALSE)
  ref <- stats::anova(stats::glm(count ~ spray, stats::poisson(),
                                 datasets::InsectSprays), test = "Chisq")
  expect_equal(f$anova$statistic[1], ref$Deviance[2], tolerance = 1e-6)
  expect_null(dim(f$data_used$count))

  mt <- datasets::mtcars
  mt$am <- data.frame(a = mt$am)
  fb <- anova_glm(mt, "am", "cyl", family = stats::binomial(), plots = FALSE)
  bref <- stats::anova(stats::glm(am ~ factor(cyl), stats::binomial(),
                                  datasets::mtcars), test = "Chisq")
  expect_equal(fb$anova$statistic[1], bref$Deviance[2], tolerance = 1e-6)
  fbin <- anova_bin(mt, "am", "cyl", plots = FALSE)
  expect_equal(fbin$anova$statistic[1], bref$Deviance[2], tolerance = 1e-6)

  ins$count2 <- data.frame(a = 1:72, b = 1:72)
  expect_error(anova_glm(ins, "count2", "spray", plots = FALSE),
               "Response `count2` is a 2-column data frame; supply a single column.",
               fixed = TRUE)
  # A successes/failures matrix keeps its specific message.
  m <- datasets::mtcars
  m$am <- cbind(m$am, 1 - m$am)
  expect_error(anova_glm(m, "am", "cyl", family = stats::binomial(), plots = FALSE),
               "successes/failures matrix response is not supported", fixed = TRUE)
})
