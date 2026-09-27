# Multivariate Analysis of Variance and Covariance

Compares two or more numeric responses jointly across groups, optionally
adjusting for covariates. Reports the multivariate test, follow-up
univariate analyses with effect sizes, covariate-adjusted estimated
marginal means, assumption checks (Mardia's multivariate normality
tests, Box's M and, with covariates, a multivariate test of homogeneous
regression slopes), and a canonical discriminant analysis with its plot.

## Usage

``` r
anova_manova(
  data,
  responses,
  groups,
  covariates = NULL,
  test = c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"),
  interaction = FALSE,
  type = c("II", "III"),
  conf_level = 0.95,
  adjust = "tukey",
  assumptions = TRUE,
  posthoc = TRUE,
  plots = TRUE,
  verbose = FALSE
)
```

## Arguments

- data:

  A data frame, or anything inheriting from one, such as a `data.table`
  or a tibble.

- responses:

  Character vector of at least two numeric response columns. A single
  response is accepted and falls back to a univariate analysis, with a
  note.

- groups:

  Character vector. One or more grouping columns.

- covariates:

  Character vector or `NULL`. Numeric covariates, which turn the
  analysis into a MANCOVA. They are mean-centred before fitting.

- test:

  Character. Multivariate statistic to report: `"Pillai"` (default),
  `"Wilks"`, `"Hotelling-Lawley"` or `"Roy"`. Pillai's trace is the most
  robust to departures from the assumptions. Roy's largest root has only
  an upper bound for its F statistic when the hypothesis has more than
  one degree of freedom, so its p-value is then a lower bound
  (anti-conservative).

- interaction:

  `FALSE` (additive, the default), `TRUE` (full factorial across the
  grouping variables), or a whole number giving the highest interaction
  order.

- type:

  Character. `"II"` (default) or `"III"`: the type of both the
  multivariate tests and the univariate follow-up tests. With `"III"`
  the multivariate model and the univariate ones are fitted under
  sum-to-zero contrasts.

- conf_level:

  Numeric in (0, 1). Level for every interval returned. Default `0.95`.

- adjust:

  Character. Multiplicity adjustment for the pairwise comparisons within
  each response. It does not adjust across responses. Default `"tukey"`.

- assumptions:

  Logical. Compute Mardia's tests and Box's M. Default `TRUE`. The
  homogeneity-of-slopes test is computed whenever there are covariates.

- posthoc:

  Logical. Compute pairwise comparisons within each response. Default
  `TRUE`.

- plots:

  Logical. Build ggplot2 objects. They are returned in `$plots`, never
  drawn. Default `TRUE`.

- verbose:

  Logical. Emit progress through
  [`message`](https://rdrr.io/r/base/message.html). Default `FALSE`.

## Value

An
[`anovakit_fit`](https://elkronos.github.io/anovakit/reference/anovakit_fit.md)
object. Its `$model` is the multivariate linear model (class `mlm`; pass
it to [`car::Manova()`](https://rdrr.io/pkg/car/man/Anova.html) for
further tests), and its `$anova` is the multivariate table. Extra
components: `$univariate` (a named list of per-response fits, each
holding that response's own `$model`, `$anova`, `$effect_sizes`,
`$emmeans`, `$emmeans_object` and `$posthoc`; the univariate p-values
are not adjusted across responses), `$multivariate` (the Type II or III
multivariate table, one row per term, with columns `term`, `df`,
`statistic`, `approx_f`, `num_df`, `den_df` and `p_value`), `$canonical`
(the discriminant axes), `$canonical_term` (the term they describe,
spelled as in `$multivariate$term`), `$slopes_test` (with covariates,
the multivariate test of homogeneous regression slopes; otherwise
`NULL`), `$covariate_means` (the means the covariates were centred at)
and `$test` (the multivariate statistic used).
`$assumptions$structure_coefficients` has a `response` column and one
column per axis. The top-level `$emmeans` and `$posthoc` stack every
response's table with a `response` column (`response.1` if a grouping
column is already called `response`); `$emmeans_object` is `NULL`,
because there is one grid per response – take them from
`$univariate[[r]]$emmeans_object` (a named list of grids, one per
factor, when the factors enter additively).

## Details

**The multivariate tests** are Type II or Type III tests, as `type`
says, computed by
[`car::Manova()`](https://rdrr.io/pkg/car/man/Anova.html) on the
multivariate linear model
`lm(cbind(<responses>) ~ <covariates> + <groups>)`. Type II tests each
term adjusted for every other term that does not contain it, so the
result does not depend on the order of `groups` or `covariates`; Type
III fits under sum-to-zero contrasts and tests each term adjusted for
all others. (`summary(stats::manova())` would give sequential,
order-dependent tests instead.) When the model has aliased coefficients
– an empty cell of a crossed design, or collinear covariates –
[`car::Manova()`](https://rdrr.io/pkg/car/man/Anova.html) refuses the
model, so the Type II tests are computed by the equivalent model
comparisons, and a note says so. Type III tests are not defined when
coefficients are aliased; Type II tests are then reported throughout,
with a note.

With `test = "Roy"` the F statistic for a term with more than one degree
of freedom (and more than one response) is an upper bound, so its
p-value is a lower bound: it is anti-conservative, and a note says so.

**Covariates** are mean-centred before fitting (the means are in
`$covariate_means` and in a note). Centring changes no test; it keeps a
covariate with a large offset and a small spread (a time stamp, say)
from being mistaken for a constant by the fitting routine, and it means
the adjusted marginal means are read at the covariate means. A MANCOVA
assumes that each covariate has the same slope in every group. That
assumption is tested by comparing the model with one in which each
covariate interacts with every grouping term, using the chosen
multivariate statistic; the result is in `$slopes_test`, and a note says
when it is rejected at the 0.05 level, since the common-slope adjusted
means and comparisons are then not interpretable as constant group
differences.

**Assumption checks.** Mardia's tests of multivariate skewness and
kurtosis (on the residuals of the multivariate model) and Box's M test
of equality of covariance matrices (across the cells of the grouping
variables, with the covariate effects removed) are computed from first
principles, so no additional package is required. Mardia's tests are
skipped, with a note, unless the residual degrees of freedom exceed the
number of responses by at least 10: when they equal it the statistics do
not depend on the data at all, and close to it they are dominated by the
design. Both are asymptotic tests, and Mardia's kurtosis test in
particular rejects too often when the sample is small relative to the
square of the number of responses. Box's M is notoriously sensitive to
non-normality: treat a small p-value as a prompt to look at the group
covariances rather than as a verdict.

**Univariate follow-ups.** One Type II or III ANOVA is fitted per
response and kept under `$univariate`. Their p-values are *not* adjusted
across responses and are reported whether or not the multivariate test
is significant; `adjust` applies only to the pairwise comparisons within
each response. A note says when a term's multivariate test is not
significant at the 0.05 level, since its univariate follow-ups are then
not protected by it.

**Canonical discriminant analysis.** The hypothesis and error sum of
squares and cross-products matrices of one term are taken from the same
[`car::Manova()`](https://rdrr.io/pkg/car/man/Anova.html) tests as
`$multivariate`, and the generalised eigenproblem is solved directly (on
responses scaled to unit error variance, so responses on very different
scales do not make it singular). `$canonical` reports the eigenvalues,
canonical correlations and the proportion of the term's between-group
variation on each axis; `$plots$canonical` plots the first two axes with
group centroids, and `$assumptions$structure_coefficients` gives the
pooled within-group correlations between each response and each axis,
computed from the error matrix. The analysis describes *one* term: the
full interaction of the grouping variables when it is in the model
(`interaction = TRUE`, or an order equal to the number of grouping
variables); otherwise the main effect of the *first* grouping variable,
`groups[1]` – list another variable first to describe it instead.
`$canonical_term` names the term, and a note says which term was used
and why. The scores are computed after removing the fitted effects of
every other term in the model (the covariates, and the other grouping
terms, under sum-to-zero coding), so that the plot shows the described
term alone. Axes are scaled to unit pooled within-group variance (the
error matrix divided by its degrees of freedom), so the spread of the
centroids on the plot means what it appears to mean; in a balanced
design the between-group to within-group ratio of the plotted scores
equals the eigenvalue exactly, and in an unbalanced one approximately.

**Estimated marginal means** come from emmeans applied to each
univariate model, so with covariates present they are adjusted means,
not raw group means. When several grouping variables enter additively
they are reported for each factor separately, averaged over the others,
with a `term` column.

Responses and covariates are named as columns rather than composed into
a formula. If you want a transformed response, create the column first.

## References

Mardia, K. V. (1970). Measures of multivariate skewness and kurtosis
with applications. *Biometrika*, 57(3), 519-530.

Box, G. E. P. (1949). A general distribution theory for a class of
likelihood criteria. *Biometrika*, 36(3/4), 317-346.

## See also

[`anova_ancova`](https://elkronos.github.io/anovakit/reference/anova_ancova.md)
for one response with covariates,
[`anova_glm`](https://elkronos.github.io/anovakit/reference/anova_glm.md)
for one response in any family.

## Examples

``` r
set.seed(1)
n <- 120
d <- data.frame(g = factor(rep(c("a", "b", "c"), each = n / 3)),
                age = rnorm(n, 40, 8))
d$score1 <- 5 + 2 * (d$g == "b") + 4 * (d$g == "c") + 0.1 * d$age + rnorm(n)
d$score2 <- 3 + 1 * (d$g == "b") + 2 * (d$g == "c") + 0.05 * d$age + rnorm(n)

fit <- anova_manova(d, c("score1", "score2"), "g")
fit
#> Multivariate analysis of variance (Pillai test, Type II) 
#> --------------------------------------------------------
#> Call: anova_manova(data = d, responses = c("score1", "score2"), groups = "g")
#> Observations used: 120
#> 
#> Omnibus test (Type II, Pillai's trace, approximate F)
#>   term df statistic approx_f num_df den_df   p_value
#> 1    g  2    0.7003    31.52      4    234 < 2.2e-16
#> 
#> Plots available: residuals_score1, qq_score1, residuals_score2, qq_score2, emmeans, canonical
#>   (use plot(x, which = "residuals_score1"))

# Assumption checks, computed directly rather than through a dependency
fit$assumptions$mardia
#>              test  statistic df   p_value
#> 1 Mardia skewness  1.2952272  4 0.8621848
#> 2 Mardia kurtosis -0.6280309 NA 0.5299837
fit$assumptions$box_m
#>   statistic df   p_value
#> 1  6.181818  6 0.4031341

# The discriminant axes, the term they describe, and how each response
# loads on them
fit$canonical
#>   axis   eigenvalue canonical_r prop_variance
#> 1 Can1 2.3344192069  0.83671841  9.999019e-01
#> 2 Can2 0.0002290369  0.01513223  9.810337e-05
fit$canonical_term
#> [1] "g"
fit$assumptions$structure_coefficients
#>   response      Can1       Can2
#> 1   score1 0.9284840 -0.3713724
#> 2   score2 0.5034467  0.8640263

# Follow-up tests, one per response, are kept separately as well as stacked
fit$univariate$score1$anova
#>        term   sum_sq  df statistic      p_value
#> 1         g 370.2427   2  117.7309 9.612167e-29
#> 2 Residuals 183.9721 117        NA           NA

# As a MANCOVA: the marginal means are adjusted for age, and the
# common-slope assumption is tested
mc <- anova_manova(d, c("score1", "score2"), "g",
                   covariates = "age", plots = FALSE)
mc$emmeans
#>   response g  estimate        se  df  conf_low conf_high
#> 1   score1 a  8.777619 0.1626659 116  8.455439  9.099800
#> 2   score1 b 11.246581 0.1626589 116 10.924415 11.568748
#> 3   score1 c 13.040960 0.1626567 116 12.718797 13.363122
#> 4   score2 a  5.096429 0.1553352 116  4.788768  5.404090
#> 5   score2 b  6.218451 0.1553284 116  5.910804  6.526098
#> 6   score2 c  7.101230 0.1553264 116  6.793587  7.408873
mc$slopes_test
#>                                   comparison df  statistic approx_f num_df
#> 1 common slopes vs covariate-by-group slopes  2 0.08020204 2.381249      4
#>   den_df    p_value homogeneous
#> 1    228 0.05244255        TRUE
```
