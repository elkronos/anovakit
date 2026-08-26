# anovakit 0.1.0

First release. The package was previously distributed as loose scripts under the
name ANOVAtoolbox; it is now an installable package called anovakit, and the
repository has been renamed to match.

Eight analysis-of-variance workflows behind one argument convention and one
return class: `anova_welch()`, `anova_kw()`, `anova_ancova()`, `anova_rm()`,
`anova_manova()`, `anova_bin()`, `anova_count()` and `anova_glm()`. Each
validates its input, checks the assumptions that matter for its method, reports
effect sizes, computes estimated marginal means and pairwise comparisons, and
returns `ggplot2` objects rather than drawing them. See `vignette("anovakit")`
for a tour and the README for the design decisions behind it.

The package supersedes a collection of loose scripts distributed from the same
repository. If you used those, note the following before comparing results.

## The functions were renamed and now share one convention

`ancova_analysis()` is `anova_ancova()` and `manova_analysis()` is
`anova_manova()`. The grouping argument, previously `iv`, `group_var`,
`group_vars` or `group_vars_vec` depending on the file, is `groups` everywhere;
the response, previously `dv` or `response_var`, is `response`.
`anova_manova()` takes column names rather than a formula. Ten helpers that were
exported under generic names — `validate_vars()`, `check_normality()`,
`compute_emm()`, `calculate_effect_sizes()` and others — are now internal.

Every function returns an `anovakit_fit` object with `print()`, `summary()`
and `plot()` methods, in place of eight differently shaped lists.

## Results that changed

Numbers the old scripts produced for these analyses should not be relied on.

* **Count models.** Exposure offsets never reached the model, marginal means came
  back without standard errors, and the negative binomial refit failed silently
  on every call and fell back to Poisson — reporting p = 1.9e-07 on overdispersed
  data where the negative binomial gives p = 0.21.
* **Type III sums of squares.** The contrasts were set after the model was
  fitted, which does nothing, so the table was simple effects at the reference
  level under a Type III label. On an unbalanced factorial that reported
  p = 0.455 where the genuine Type III p-value was 4.2e-09.
* **ANCOVA.** The group effect was reported at covariate = 0 whenever the
  interaction was retained. Covariates are now mean-centred by default.
* **MANCOVA.** Marginal means fell back to raw group means with standard errors
  2.4 times too large.
* **Repeated measures.** Sphericity was never computed, no diagnostic plots were
  returned, and a factor response was converted through its level indices.
* **Welch.** Normality was tested on pooled residuals — which fails for perfectly
  normal groups whenever their variances differ, the exact case the method
  exists for — and the column labelled `cohen_d` held Hedges' g.
* **Effect sizes.** Omega squared used the non-partial formula on Type II and
  Type III tables, where it is not a proportion of variance; it is now partial
  omega squared, matching `effectsize::omega_squared(partial = TRUE)`. Every
  effect-size table also carried a spurious `Residuals` row at exactly 0.5.

## What the package now guarantees

* **Nothing prints.** No analysis function writes to the console and no plot is
  drawn as a side effect. Warnings and messages raised by `car`, `glm()`,
  `emmeans` and `afex` are captured into `$notes`.
* **`conf_level` sets every interval** in the tables a function returns.
* **`type = "III"` fits under sum-to-zero contrasts**, so the label matches the
  table.
* **Row accounting closes.** `nrow($data_used) + $n_removed` is always the number
  of rows supplied, including when rows leave for an incomplete repeated-measures
  design or afex aggregation.
* **`$notes` records every decision** made on the caller's behalf: the model
  chosen from a dispersion statistic or a slopes test, a centred covariate, a
  dropped subject, an unidentified odds ratio.
* **No escape hatches.** Only `anova_rm()` takes `...`, and it documents that the
  arguments go to `afex::aov_ez()`. Prior weights are named as a column, so they
  are subsetted with the data rather than misaligning when a row is dropped.

## New capabilities

* Mardia's tests of multivariate skewness and kurtosis, Box's M, and a canonical
  discriminant analysis with its plot, in `anova_manova()`. The previous README
  advertised these; the code contained none of them.
* Dunn's test with the standard tie correction, implemented directly.
* `interaction` on `anova_bin()`, `anova_count()`, `anova_glm()` and
  `anova_manova()`; `model` on `anova_count()` to choose the count family
  explicitly; `force_interaction` on `anova_ancova()` to fix the model in advance
  rather than select it with a test on the same data; per-group covariate slopes
  from `emmeans::emtrends()`.
* `posthoc = FALSE` on every function that computes comparisons, and
  heteroskedasticity-consistent standard errors through `vcov_type` that flow
  into the marginal means and the comparisons.

## Dependencies

Imports `car`, `emmeans`, `ggplot2` and base R. `afex`, `MASS` and `sandwich`
are Suggests behind `requireNamespace()` guards. The scripts previously needed
eighteen packages, including `dplyr`, `magrittr`, `rlang`, `multcomp`,
`nortest`, `lmtest`, `effectsize`, `effsize` and `dunn.test`.

## Tests

A `testthat` (edition 3) suite that checks results against `stats::oneway.test`,
`stats::kruskal.test`, `stats::manova`, `car::Anova`, `MASS::glm.nb`,
`afex::aov_ez` and hand-computed effect sizes, with a regression test for each
behaviour described above. The previous suite consisted of `cat()`-based
harnesses whose runners were commented out, so `R CMD check` ran no assertions.
