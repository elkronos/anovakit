# Welch's Analysis of Variance

Compares the mean of a numeric response across groups without assuming
equal variances, using
[`oneway.test`](https://rdrr.io/r/stats/oneway.test.html). Pairwise
Welch t-tests with a multiplicity adjustment, and standardised mean
differences, are also returned for any number of groups from two
upwards; with two groups the single comparison is the omnibus test.

## Usage

``` r
anova_welch(
  data,
  response,
  groups,
  conf_level = 0.95,
  adjust = "holm",
  hedges_correction = TRUE,
  posthoc = TRUE,
  plots = TRUE,
  verbose = FALSE
)
```

## Arguments

- data:

  A data frame, or anything inheriting from one, such as a `data.table`
  or a tibble.

- response:

  Character. Name of the numeric response column.

- groups:

  Character vector. One or more grouping columns.

- conf_level:

  Numeric in (0, 1). Level for every interval returned: group means,
  pairwise differences and standardised effect sizes. Default `0.95`.

- adjust:

  Character. Multiplicity adjustment for the pairwise comparisons'
  p-values, passed to
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html). One of `"holm"`,
  `"hochberg"`, `"hommel"`, `"bonferroni"`, `"BH"`, `"BY"`, `"fdr"` or
  `"none"`, spelled out in full. Default `"holm"`. Note that `"tukey"`
  is not available here, as it is on the functions whose comparisons
  come from emmeans: these are Welch t-tests, not linear contrasts on a
  common error term.

- hedges_correction:

  Logical. Apply the small-sample bias correction to the standardised
  mean differences (see Details). Default `TRUE`, in which case the
  effect size column is named `hedges_g`; otherwise `cohens_d`. Both use
  the average-variance standardiser.

- posthoc:

  Logical. Compute the pairwise comparisons. There are `choose(k, 2)` of
  them, so this is worth turning off when the number of groups is large.
  Default `TRUE`. The standardised mean differences are computed pair by
  pair alongside them, so `$effect_sizes` is empty when this is `FALSE`;
  `$notes` says so.

- plots:

  Logical. Build ggplot2 objects. They are returned in `$plots`, never
  drawn. Default `TRUE`.

- verbose:

  Logical. Emit progress through
  [`message`](https://rdrr.io/r/base/message.html). Default `FALSE`.

## Value

An
[`anovakit_fit`](https://elkronos.github.io/anovakit/reference/anovakit_fit.md)
object.

- `$anova`: the Welch test, one row. Its `term` is the grouping column,
  or `"A x B cells"` when several were combined.

- `$emmeans`: the group means with t intervals (raw group summaries; no
  model is fitted).

- `$posthoc`: one row per pair of groups, with `group1`, `group2`,
  `difference` (the mean of `group1` minus the mean of `group2`),
  `conf_low` and `conf_high` (the Welch interval for that difference at
  `conf_level`; an unadjusted, per-comparison interval, not a
  simultaneous one), `statistic` (Welch's t, same sign as `difference`),
  `df`, `p_value` (unadjusted), `p_adjusted` (adjusted by `adjust`) and
  `adjustment`. This is the same column layout as the functions built on
  emmeans, whose intervals are however adjusted.

- `$effect_sizes`: one row per pair, with `group1`, `group2`, `hedges_g`
  (or `cohens_d`), the same direction as `difference`, its interval
  `conf_low` and `conf_high`, `magnitude`, and `standardiser`, which
  names the standard deviation used.

- `$assumptions`: a per-group Shapiro-Wilk table in `normality`, with a
  `note` column giving the reason for any group that was not tested, and
  the ratio of the largest to the smallest group variance in
  `variance_ratio`.

- `$model`: the `htest` returned by
  [`oneway.test`](https://rdrr.io/r/stats/oneway.test.html); there is no
  fitted model object and `$emmeans_object` is `NULL`.

`$posthoc` and `$effect_sizes` are `NULL` when `posthoc = FALSE`.

## Details

Welch's test is the appropriate default for comparing means: it is
barely less powerful than the classical F test when variances are equal,
and it keeps its nominal error rate when they are not.

Because the method explicitly permits unequal variances, normality is
assessed *within each group* rather than on pooled residuals. Pooling
residuals from groups with different spreads produces a mixture
distribution that fails normality tests even when every group is
perfectly normal, which would argue against the very method being used.
For the same reason the Q-Q plot shows each observation's deviation from
its group mean divided by its group's standard deviation, so every group
contributes on its own scale. A group whose Shapiro-Wilk test could not
be run (fewer than 3 or more than 5000 observations) is named in
`$notes`, with the reason.

**Standardised mean differences.** The pairwise effect sizes do not
assume equal variances either. Each difference in means is divided by
the average-variance standard deviation of its two groups, \\s^\* =
\sqrt{(s_1^2 + s_2^2)/2}\\, rather than by the pooled standard
deviation, whose meaning and sampling variance both depend on equal
variances. The interval is for the population value \\\delta^\* =
(\mu_1 - \mu_2) / \sqrt{(\sigma_1^2 + \sigma_2^2)/2}\\. It inverts the
noncentral t distribution of Welch's statistic on Welch's degrees of
freedom and rescales the bounds on the noncentrality to the
average-variance standard deviation, which is the interval
`effectsize::cohens_d(pooled_sd = FALSE)` reports. In simulation
(variance ratios from 1/4 to 16 in either group, group sizes from 2 to
50, true differences of 0 and 0.8) a 95% interval covered 93-97%
whenever the smaller group had at least five observations, 93-98% for
equal groups of two or three, and 86-95% when a group of two or three
sat beside one ten or more times its size, where no interval does well:
with two observations a group's standard deviation is barely estimated.
Bonett's (2008) normal-theory interval for the same quantity was also
examined; it was wider than needed for equal small groups (99-100% at
two per group) and covered less in the unbalanced cases (81-94%).

With `hedges_correction = TRUE` the estimate is multiplied by the
small-sample bias correction \\J(\nu) = \Gamma(\nu/2) / (\sqrt{\nu/2}\\
\Gamma((\nu-1)/2))\\, evaluated at the Satterthwaite degrees of freedom
of the standardiser, \\\nu = (s_1^2 + s_2^2)^2 / (s_1^4/(n_1-1) +
s_2^4/(n_2-1))\\ (Delacre et al., 2021). When the two sample variances
are equal and the groups are the same size this is the familiar \\n_1 +
n_2 - 2\\; when the variances differ it removes the bias that the
pooled-variance correction leaves. The correction is a property of the
point estimate, so the interval, which is for \\\delta^\*\\ itself, is
the same with or without it. With two or three observations per group
and a low `conf_level`, the corrected estimate can fall outside that
interval; `$notes` says so when it does.

The standardiser is recorded in the `standardiser` column of
`$effect_sizes`. The magnitude labels use Cohen's conventional
thresholds of 0.2, 0.5 and 0.8 on the absolute value.

**Several grouping variables.** They are combined into a single cell
factor and unused level combinations are dropped, so the number of
groups reported is the number that actually contain data. The result is
one omnibus test across all populated cells, labelled `"A x B cells"` in
`$anova`: it cannot separate a main effect of one factor from a main
effect of another, and it cannot test an interaction. A note to that
effect is added to the returned object. For a factorial analysis use
[`anova_glm`](https://elkronos.github.io/anovakit/reference/anova_glm.md)
with `interaction = TRUE`.

## References

Bonett, D. G. (2008). Confidence intervals for standardized linear
contrasts of means. *Psychological Methods*, 13(2), 99-109.

Delacre, M., Lakens, D., Ley, C., Liu, L., & Leys, C. (2021). Why
Hedges' g\*s based on the non-pooled standard deviation should be
reported with Welch's t-test. *PsyArXiv*.

Welch, B. L. (1951). On the comparison of several mean values: an
alternative approach. *Biometrika*, 38(3/4), 330-336.

## See also

[`anova_kw`](https://elkronos.github.io/anovakit/reference/anova_kw.md)
for a rank-based alternative,
[`anova_ancova`](https://elkronos.github.io/anovakit/reference/anova_ancova.md)
to adjust for a covariate.

## Examples

``` r
set.seed(123)
d <- data.frame(
  group = rep(c("A", "B", "C"), each = 30),
  value = c(rnorm(30, 10, 2), rnorm(30, 12, 2.5), rnorm(30, 9, 1.5))
)
fit <- anova_welch(d, "value", "group")
fit
#> Welch's analysis of variance 
#> ----------------------------
#> Call: anova_welch(data = d, response = "value", groups = "group")
#> Observations used: 90
#> 
#> Omnibus test
#>    term statistic num_df den_df  p_value
#> 1 group     28.44      2  55.18 3.24e-09
#> 
#> Plots available: means, box, qq
#>   (use plot(x, which = "means"))
fit$posthoc
#>   group1 group2 difference     conf_low conf_high statistic       df
#> 1      A      B -2.5400534 -3.587216430 -1.492890 -4.855868 57.77758
#> 2      A      C  0.8691619  0.005291539  1.733032  2.020414 50.45176
#> 3      B      C  3.4092153  2.505770188  4.312660  7.584648 48.65250
#>        p_value   p_adjusted adjustment
#> 1 9.503682e-06 1.900736e-05       holm
#> 2 4.866484e-02 4.866484e-02       holm
#> 3 8.654470e-10 2.596341e-09       holm

# Standardised mean differences on the average-variance standard deviation,
# with an interval that does not assume equal variances
fit$effect_sizes
#>   group1 group2   hedges_g     conf_low  conf_high magnitude
#> 1      A      B -1.2374221 -1.804332292 -0.6941522     large
#> 2      A      C  0.5138686  0.003040726  1.0353311    medium
#> 3      B      C  1.9279766  1.314162066  2.5895992     large
#>              standardiser
#> 1 sqrt((s1^2 + s2^2) / 2)
#> 2 sqrt((s1^2 + s2^2) / 2)
#> 3 sqrt((s1^2 + s2^2) / 2)

# Two grouping variables are combined into one cell factor
d$site <- rep(c("north", "south"), 45)
anova_welch(d, "value", c("group", "site"), plots = FALSE)$anova
#>                 term statistic num_df   den_df     p_value
#> 1 group x site cells  11.41081      5 38.34597 8.68667e-07
```
