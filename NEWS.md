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

## Changes made after an adversarial review, before release

A multi-agent adversarial review reproduced 167 problems; all were fixed, each
with a regression test against an independent reference. The changes a user
will notice:

* **MANOVA** reports genuine Type II / Type III multivariate tests from
  `car::Manova()`; it previously reported sequential (Type I) tests whatever
  `type` said, so results depended on the order of `groups`. `$model` is now
  the multivariate `mlm`. With covariates, covariates are centred and a
  multivariate homogeneity-of-slopes test is reported in `$slopes_test`.
* **ANCOVA** gains `interaction`, defaulting to `TRUE`: with several grouping
  factors their interaction is fitted, and the slopes test and simple slopes
  follow the cells. `vcov_type` now reaches the F tests (robust Wald F) and
  the simple slopes. A 0/1 covariate is held at its mean in `$emmeans`.
* **Marginal means and comparisons.** For an additive multi-factor model
  (`interaction = FALSE`) they are reported per factor, with a `term` column,
  instead of for every cell of the grid. Empty cells are left out or flagged.
  Comparisons above 5000 are skipped with a note. `$posthoc` has
  `p_value` (unadjusted), `p_adjusted` and `adjustment` in every function.
  Non-log, non-logit links compare on the response scale instead of
  mislabelled link-scale "ratios". Quasi families use t. Global
  `emm_options()` no longer change results.
* **Effect sizes.** `anova_welch()` standardises by the average of the two
  group variances, with the noncentral-t interval of Welch's statistic (as
  `effectsize::cohens_d(pooled_sd = FALSE)` reports), so the interval is valid
  under unequal variances (a `standardiser` column says so). Odds ratios
  (`anova_bin()`) and rate ratios (`anova_count()`) are each a level against
  its reference level, whatever the contrasts, with `factor` and `comparison`
  columns. McFadden's R squared is reported only for binomial, Poisson and
  negative binomial families, and is invariant to aggregation.
* **Weights.** Zero-weight rows are dropped and counted in `$n_removed`.
  `anova_bin()`'s null model, model-versus-null test and proportion table use
  the weights. Whole-number weights above 1 are treated as frequencies for
  the dispersion in `anova_count()` and `anova_glm()`.
* **Intervals.** Profile intervals use t cut-offs when the dispersion is
  estimated. Under separation or a zero-count cell an interval is open on the
  side the estimate diverged in, and its finite end is profiled directly.
* **Robust standard errors** fall back to model-based ones, with a note, at
  the boundary (separation, an all-zero group) and at leverage 1.
* **Separation** is detected whatever the contrasts, and zero-count cells in
  count models are flagged.
* **Repeated measures.** Sphericity corrections are computed directly when
  afex cannot supply them, so an uncorrected p-value is never labelled "GG";
  `observed` reaches generalised eta squared; normality is tested per
  within-subject cell; `$residuals` follow `$data_used`; column names of any
  kind work; `...` accepts only `fun_aggregate`, `observed`, `type` and
  `anova_table = list(p_adjust_method = )`.
* **Count models.** An explicit `model = "negbin"` that cannot be fitted is an
  error rather than a silent switch to quasi-Poisson; the theta and
  overdispersion notes state the size of the problem they describe.
* **Kruskal-Wallis** keeps infinite values (ranks can use them), counts ties
  on the ranks, and `anova_welch()`/`anova_kw()` label a combined grouping
  "A x B cells".
* **Input.** Non-ASCII group labels no longer crash in UTF-8 sessions; an
  explicit `NA` factor level counts as missing; distinct values that print
  alike stay distinct groups; colliding cell labels, a column in two roles, and
  the names `Residuals` and `(Intercept)` are refused.
* **Printing.** `summary()` returns a `summary.anovakit_fit` object and shows
  the method's own tables; counts are not rounded; p-values below machine
  precision print as such; the omnibus heading names the test.
* **Packaging.** `ggplot2 (>= 3.4.0)` and `testthat (>= 3.1.8)` are declared.

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
