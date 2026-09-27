# The object returned by every anovakit function

All eight analysis functions return an object of class `anovakit_fit`.
It is a plain list, so `$` extraction works as usual, with `print`,
`summary` and `plot` methods for convenience.

## Components

- method:

  Character. The analysis that was run, as printed.

- call:

  The matched call, so a result can always say how it was made.

- model:

  The fitted model object: an `lm` from
  [`anova_ancova`](https://elkronos.github.io/anovakit/reference/anova_ancova.md),
  a `glm` (or a `negbin` from
  [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html)) from
  [`anova_bin`](https://elkronos.github.io/anovakit/reference/anova_bin.md),
  [`anova_count`](https://elkronos.github.io/anovakit/reference/anova_count.md)
  and
  [`anova_glm`](https://elkronos.github.io/anovakit/reference/anova_glm.md),
  the multivariate `mlm` from
  [`anova_manova`](https://elkronos.github.io/anovakit/reference/anova_manova.md),
  an `afex_aov` from
  [`anova_rm`](https://elkronos.github.io/anovakit/reference/anova_rm.md)
  (fitted on internal names, see its `$internal_names`), and the `htest`
  returned by [`oneway.test`](https://rdrr.io/r/stats/oneway.test.html)
  or [`kruskal.test`](https://rdrr.io/r/stats/kruskal.test.html) from
  [`anova_welch`](https://elkronos.github.io/anovakit/reference/anova_welch.md)
  and
  [`anova_kw`](https://elkronos.github.io/anovakit/reference/anova_kw.md),
  which fit no model. The call of a fitted `lm` or `glm` reaches the
  analysed rows from anywhere, so
  [`update()`](https://rdrr.io/r/stats/update.html),
  [`step()`](https://rdrr.io/r/stats/step.html) and
  [`lmtest::lrtest()`](https://rdrr.io/pkg/lmtest/man/lrtest.html) work
  on it directly: for instance `update(fit$model, . ~ 1)`.

- anova:

  Data frame. The omnibus test table. Its attributes record the type of
  sums of squares, the test statistic and (in `anova_rm`) the sphericity
  correction; [`print()`](https://rdrr.io/r/base/print.html) shows them.

- effect_sizes:

  Data frame, or `NULL`. What it holds depends on the method: partial
  eta squared with partial omega squared (`anova_ancova`,
  `anova_manova`, and `anova_glm` on a Gaussian family), a standardised
  mean difference on the average of the two group variances, with
  intervals (`anova_welch`), epsilon squared and eta squared for the
  rank statistic (`anova_kw`), partial and generalised eta squared
  (`anova_rm`), odds ratios (`anova_bin`) and incidence rate ratios
  (`anova_count`), each a level against its reference level whatever the
  contrasts, and the deviance explained (`anova_glm` on any other
  family) with McFadden's pseudo R squared for the binomial, Poisson and
  negative binomial families only. Both variance measures are *partial*:
  the classical omega squared is only a proportion of variance when the
  effect sums of squares partition the total, which Type II and Type III
  sums of squares do not.

- emmeans:

  Data frame, or `NULL`. Estimated marginal means with intervals at
  `conf_level` (group summaries for `anova_welch` and `anova_kw`). When
  grouping factors enter a model additively they are reported per
  factor, with a `term` column.

- emmeans_object:

  The `emmGrid`, a named list of them (one per factor) when the means
  are reported per factor, or `NULL`. Pass it to emmeans for contrasts
  the wrapper does not cover. It carries emmeans' own default confidence
  level rather than `conf_level`, so give `level =` when you summarise
  it. It is `NULL` for
  [`anova_welch`](https://elkronos.github.io/anovakit/reference/anova_welch.md)
  and
  [`anova_kw`](https://elkronos.github.io/anovakit/reference/anova_kw.md),
  which fit no model emmeans can use, and for
  [`anova_manova`](https://elkronos.github.io/anovakit/reference/anova_manova.md),
  where there is one grid per response under
  `$univariate[[response]]$emmeans_object`.

- posthoc:

  Data frame, or `NULL`. Pairwise comparisons, with the unadjusted
  p-value in `p_value`, the adjusted one in `p_adjusted` and the method
  in `adjustment`.

- assumptions:

  Named list of assumption checks. Contents vary by method; each element
  is a test object, a data frame or `NULL`. It is empty for
  [`anova_kw`](https://elkronos.github.io/anovakit/reference/anova_kw.md)
  unless `diagnostics = TRUE`, since the test assumes no distribution.

- plots:

  Named list of ggplot2 objects. Empty when `plots = FALSE`. Nothing is
  ever drawn as a side effect.

- data_used:

  Data frame. The rows and columns the model was fitted on: the analysed
  columns only, after missing and infinite values, zero weights and
  incomplete subjects are removed. Covariates are mean-centred there
  when `anova_ancova` or `anova_manova` centred them.

- n_removed:

  Integer. Input rows that are not in `data_used`: dropped for missing
  or infinite values (a factor level that is itself `NA` counts as
  missing) or a prior weight of zero, or, in
  [`anova_rm`](https://elkronos.github.io/anovakit/reference/anova_rm.md),
  for an incomplete within-subject design or because repeated
  subject-by-cell rows were aggregated (with `fun_aggregate`, the mean
  by default). `nrow(data_used) + n_removed` is always the number of
  rows given.

- conf_level:

  Numeric. The level used for every interval returned.

- notes:

  Character vector. Everything the function decided on your behalf,
  could not compute, or thinks you should know – including anything car,
  [`glm`](https://rdrr.io/r/stats/glm.html), emmeans, sandwich or afex
  said while the model was being fitted. Always read this.

## Components individual functions add

Each function returns everything above plus whatever its own method
produces. `names(fit)` lists them all.
[`anova_ancova`](https://elkronos.github.io/anovakit/reference/anova_ancova.md)
adds `$slopes_test`, `$simple_slopes`, `$covariate_means`,
`$model_additive` and `$model_interaction`;
[`anova_rm`](https://elkronos.github.io/anovakit/reference/anova_rm.md)
adds `$sphericity`, `$subjects_dropped`, the breakdown of `$n_removed`,
`$residuals` (within-subject residuals in the row order of `$data_used`)
and `$internal_names`;
[`anova_manova`](https://elkronos.github.io/anovakit/reference/anova_manova.md)
adds `$multivariate`, `$univariate`, `$canonical`, `$canonical_term`,
`$slopes_test` (with covariates), `$covariate_means` and `$test`;
[`anova_count`](https://elkronos.github.io/anovakit/reference/anova_count.md)
adds `$model_type`, `$dispersion` and `$model_dispersion`;
[`anova_bin`](https://elkronos.github.io/anovakit/reference/anova_bin.md)
adds `$model_stats`;
[`anova_glm`](https://elkronos.github.io/anovakit/reference/anova_glm.md)
adds `$model_stats` and `$family`. Each is documented on the function
that produces it.

## See also

[`print.anovakit_fit`](https://elkronos.github.io/anovakit/reference/print.anovakit_fit.md),
[`summary.anovakit_fit`](https://elkronos.github.io/anovakit/reference/summary.anovakit_fit.md),
[`plot.anovakit_fit`](https://elkronos.github.io/anovakit/reference/plot.anovakit_fit.md)
