# Deterministic fixtures.
#
# Every generator seeds itself from an argument with an explicit default, so a
# caller can vary the data on purpose. No generator resets the global stream as
# a side effect of being called with its default.

fx_oneway <- function(n_per = 25, means = c(A = 5, B = 7, C = 4),
                      sds = NULL, seed = 101) {
  if (is.null(sds)) sds <- rep(1, length(means))
  if (length(sds) != length(means)) {
    stop("`sds` must be the same length as `means`.", call. = FALSE)
  }
  withr::with_seed(seed, {
    g <- rep(names(means), each = n_per)
    y <- unlist(Map(function(m, s) stats::rnorm(n_per, m, s), means, sds),
                use.names = FALSE)
    data.frame(group = factor(g), value = y, stringsAsFactors = FALSE)
  })
}

fx_twoway <- function(n = 120, seed = 102) {
  withr::with_seed(seed, {
    d <- data.frame(
      g1 = factor(rep(c("a", "b"), each = n / 2)),
      g2 = factor(rep(rep(c("x", "y"), each = n / 4), 2))
    )
    d$value <- 10 + 2 * (d$g1 == "b") + 1 * (d$g2 == "y") +
      3 * (d$g1 == "b") * (d$g2 == "y") + stats::rnorm(n)
    d
  })
}

fx_nested <- function(n_per = 20, seed = 103) {
  # g1 and g2 are perfectly confounded: only a-x and b-y ever occur.
  withr::with_seed(seed, {
    data.frame(
      g1 = factor(rep(c("a", "b"), each = n_per)),
      g2 = factor(rep(c("x", "y"), each = n_per)),
      value = stats::rnorm(2 * n_per)
    )
  })
}

fx_ancova <- function(n = 150, covariate_mean = 100, slope_shift = 0,
                      seed = 104) {
  # slope_shift multiplies the covariate inside group B, so a non-zero value
  # makes the slopes genuinely differ.
  withr::with_seed(seed, {
    g <- factor(rep(c("A", "B"), each = n / 2))
    x <- stats::rnorm(n, covariate_mean, 5)
    y <- 2 * x + (g == "B") * slope_shift * x + stats::rnorm(n, 0, 3)
    data.frame(dv = y, iv = g, cov = x)
  })
}

fx_binary <- function(n = 300, probs = c(a = 0.2, b = 0.5, c = 0.8),
                      seed = 105) {
  withr::with_seed(seed, {
    g <- factor(sample(names(probs), n, replace = TRUE))
    data.frame(g = g, y = stats::rbinom(n, 1, probs[as.character(g)]))
  })
}

fx_separated <- function(seed = 106) {
  withr::with_seed(seed, {
    data.frame(
      g = factor(c(rep("a", 100), rep("b", 100), rep("c", 12))),
      y = c(stats::rbinom(100, 1, 0.3), stats::rbinom(100, 1, 0.5), rep(1, 12))
    )
  })
}

fx_counts <- function(n = 400, seed = 107) {
  withr::with_seed(seed, {
    d <- data.frame(
      g1 = factor(rep(c("A", "B"), each = n / 2)),
      g2 = factor(rep(rep(c("X", "Y"), each = n / 4), 2))
    )
    lambda <- with(d, ifelse(g1 == "A" & g2 == "X", 5,
                      ifelse(g1 == "A" & g2 == "Y", 10,
                      ifelse(g1 == "B" & g2 == "X", 15, 20))))
    d$count <- stats::rpois(n, lambda)
    d$hours <- stats::runif(n, 0.5, 4)
    d
  })
}

fx_overdispersed <- function(n_per = 100, theta = 0.7, seed = 108) {
  withr::with_seed(seed, {
    g <- factor(rep(c("A", "B", "C"), each = n_per))
    data.frame(g = g,
               y = stats::rnbinom(3 * n_per,
                                  mu = rep(c(5, 6, 7), each = n_per),
                                  size = theta))
  })
}

fx_repeated <- function(n_subj = 24, seed = 109) {
  withr::with_seed(seed, {
    d <- expand.grid(id = factor(seq_len(n_subj)),
                     time = factor(c("t1", "t2", "t3")))
    d$arm <- factor(rep(rep(c("ctrl", "trt"), each = n_subj / 2), 3))
    d$score <- 10 + 2 * as.numeric(d$time) + 1.5 * (d$arm == "trt") +
      stats::rnorm(nrow(d), 0, 2)
    d[order(d$id), ]
  })
}

fx_multivariate <- function(n = 120, seed = 110) {
  withr::with_seed(seed, {
    d <- data.frame(g = factor(rep(c("a", "b", "c"), each = n / 3)),
                    age = stats::rnorm(n, 40, 8))
    d$score1 <- 5 + 2 * (d$g == "b") + 4 * (d$g == "c") + 0.1 * d$age +
      stats::rnorm(n)
    d$score2 <- 3 + 1 * (d$g == "b") + 2 * (d$g == "c") + 0.05 * d$age +
      stats::rnorm(n)
    d
  })
}

#' Every component the shared class promises
expect_fit_shape <- function(fit) {
  testthat::expect_s3_class(fit, "anovatoolbox_fit")
  for (nm in c("method", "call", "anova", "assumptions", "plots",
               "data_used", "n_removed", "conf_level", "notes")) {
    testthat::expect_true(nm %in% names(fit),
                          info = paste("missing component:", nm))
  }
  testthat::expect_type(fit$method, "character")
  testthat::expect_type(fit$notes, "character")
  testthat::expect_true(is.data.frame(fit$data_used))
  testthat::expect_type(fit$n_removed, "integer")
  invisible(fit)
}

#' Assert that a call writes nothing to stdout
expect_silent_stdout <- function(expr) {
  out <- utils::capture.output(invisible(force(expr)))
  testthat::expect_identical(length(out), 0L)
}
