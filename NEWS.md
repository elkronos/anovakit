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
with a regression test against an independent reference, and a second,
independent pass re-verified every fix and closed the gaps it found. The
changes a user will notice:

* **MANOVA** reports genuine Type II / Type III multivariate tests from
  `car::Manova()`; it previously reported sequential (Type I) tests whatever
  `type` said, so results depended on the order of `groups`. `$model` is now
  the multivariate `mlm`. With covariates, covariates are centred and a
  multivariate homogeneity-of-slopes test is reported in `$slopes_test`.
* **ANCOVA** gains `interaction`, defaulting to `TRUE`: with several grouping
  factors their interaction is fitted, and the slopes test and simple slopes
  follow the cells. `vcov_type` now reaches the F tests (robust Wald F) and
  the simple slopes. A 0/1 covariate is held at its mean in `$emmeans`. A
  design with more coefficients than rows (say, the full factorial of three
  22-level factors on 300 rows) is refused at once, naming the remedy.
* **Marginal means and comparisons.** For an additive multi-factor model
  (`interaction = FALSE`) they are reported per factor, with a `term` column,
  instead of for every cell of the grid. Empty cells are left out or flagged.
  In the six model-based functions, comparisons above 5000 are skipped with a
  note. `$posthoc` has
  `p_value` (unadjusted), `p_adjusted` and `adjustment` in every function.
  Non-log, non-logit links compare on the response scale instead of
  mislabelled link-scale "ratios". Quasi families use t. Global
  `emm_options()` no longer change results. When the marginal means cannot
  be computed, every function says that the comparisons and the plot were
  skipped. A grouping column named like a result column (`se`, `estimate`,
  ...) no longer overwrites a statistic.
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
  the weights. Whole-number weights above 1 on a count or 0/1 response in a
  binomial, Poisson or negative binomial model are frequency weights: the dispersion (`anova_count()`, `anova_glm()`) and the
  robust covariance (`anova_bin()`, `anova_count()`, `anova_glm()`) are those
  of the data expanded to one row per observation. Before, robust standard
  errors of a frequency table were an order of magnitude too large.
* **Intervals.** Profile intervals use t cut-offs when the dispersion is
  estimated. Under separation or a zero-count cell an interval is open on the
  side the estimate diverged in, and its finite end is profiled directly, in
  every GLM family; the other coefficients keep their profile intervals.
* **Robust standard errors** fall back to model-based ones, with a note, at
  the boundary (separation, an all-zero group) and at leverage 1 (a
  single-observation cell, where HC0 and HC1 give a zero standard error).
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
  precision print as such, and Tukey, Dunnett and Sidak p-values below the
  precision of their computation print as a bound; the omnibus heading names
  the test.
* **The stored model** can be refitted directly: `update(fit$model, . ~ 1)`,
  `step()` and `lmtest::lrtest()` work on `fit$model` from any frame.
* **Packaging.** `ggplot2 (>= 3.4.0)` and `testthat (>= 3.1.8)` are declared.

## Plot theme option

* `options(anovakit.theme = ...)` sets the base theme of every plot the
  package draws: a complete ggplot2 theme, or a function returning one
  (`ggplot2::theme_get` follows `theme_set()`). Unset, it is
  `theme_minimal()`, as before. The package's adjustments (bold title, legend
  at the bottom, rotated long labels) are added on top, and no result changes.
  See `?anovakit` and the plots article.

## Changes made after two further reviews, before release

Two more adversarial reviews found about 100 further problems, mostly in edge
cases. Each is fixed and has a regression test.

* **Welch and Kruskal-Wallis.** Noncentral-t intervals for g are exact for
  very well-separated groups and faster. The Hedges correction is not applied
  to a pair whose standardiser has fewer than 2 degrees of freedom, and a note
  says so. When Welch's F is driven by rounding noise in near-constant groups,
  the omnibus test is NA with a note. A group lost entirely to missing values
  is named. Omnibus headers name their statistic.
* **ANCOVA.** Marginal means that are not estimable are left out and named,
  not reported as NA rows. The coefficient-count guard counts observed cells
  only. Levene's test leaves out cells of three or fewer observations, cells
  the model fits exactly (such as a cell with no more observations than
  coefficients under a slope per group) and cells it leaves with fewer than
  2.5 residual degrees of freedom; `$notes` names them. A covariate with a large offset is handled correctly. The covariate plot draws
  the fitted model's lines. The docs no longer say randomised groups are
  unaffected by the slopes pretest.
* **Repeated measures.** Residuals under a large offset, zero error variance
  and between-subject levels lost to incomplete subjects are handled, with
  notes.
* **MANOVA.** A constant or group-determined covariate is refused, as in
  `anova_ancova()`. Exact fits are detected when a response is constant. Box's
  M is NA with a note when a cell's covariance matrix is singular up to
  rounding, and its sample-size note refers to cells.
* **Binary, count and GLM models.** Grouping columns whose names and levels
  paste to the same coefficient (`a`/`b1` and `ab`/`1`) no longer corrupt the
  results. Frequency-weighted binomial fits start from the overall proportion,
  so a cloglog link no longer diverges. A constant or exactly fitted response
  with an estimated dispersion is an error instead of an absurd F. This
  includes a response that is constant within every group, at 0 or 1. A
  constant count with a varying offset is still analysed, because its rates
  differ. A response with a large level and small but real variation (epoch
  seconds, pressures in Pa) is analysed, not mistaken for an exact fit; where
  rounding swamps a Poisson, Gamma or negative binomial deviance test on such
  a level, the test is NA with a note, and a Wald test is still available. A
  saturated model that estimates its dispersion reports NA tests; a saturated
  Poisson or binomial model still gets its Wald chi-square table.
  `test_statistic = "LR"` on a quasi family gives F, and `$method` says
  "quasi-likelihood F test" when F is asked of a Poisson, binomial or negative
  binomial model.
* **Weights.** Frequency weights apply only to the binomial, Poisson and
  negative binomial families. With non-integer weights:
  * binomial AIC and BIC are NA;
  * in `anova_bin()` and `anova_count()`, the count dispersion and the
    sparse-data check do not depend on the weights' scale;
  * a note warns that the scale matters. The reported tests and the Pearson
    dispersion of `anova_glm()` scale with the weights, as in `glm()`.
  * With unequal weights that are not frequencies the tests are inflated by
    an amount that depends on where the heavy rows fall. `$notes` gives it,
    from the sandwich covariance, for Poisson, negative binomial,
    quasi-Poisson and binomial fits, and for F tests (`anova_glm()`,
    `anova_count()`, `anova_bin()`). A binomial's whole-number weights on
    proportions are its trials, which inflate nothing.

  A group whose weights are all zero is named and dropped.
* **Negative binomial.** A fit that collapses to Poisson is reported as such.
  A theta that runs away on overdispersed data is reported as not converged.
* **Intervals.** Separated coefficients profile to an open end, or to
  (-Inf, Inf) when the likelihood is flat. `$notes` names rows where a profile
  interval and its Wald p-value disagree. When left at the default,
  profiling falls back to Wald intervals for large models. This applies in
  `anova_bin()` too. Gaussian identity profile intervals are computed in
  closed form.
* **Marginal means and comparisons.** Every additive model notes the empty
  level combinations its marginal means average over. Ratio comparisons keep
  level labels that contain " - ". Scheffe adjustments with an empty cell use
  the rank of the estimable cells.
* **Input.** `bit64::integer64` columns are converted to double. A matrix
  column with more than one column is refused. The dropped-rows note names
  only the columns that held missing or infinite values. A non-ASCII
  character response no longer fails in binomial fits. `$call` stays compact
  under `do.call()`.
* **Printing.** Tables no longer pad values, p-values included, with false
  trailing zeros.
  `summary()` shows `$model_stats`.

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
