# Analysis of Covariance

Compares the mean of a numeric response across groups while adjusting
for one or more numeric covariates. Tests the homogeneity-of-slopes
assumption, checks residual normality and equality of variances, reports
Type II or Type III sums of squares with partial eta squared, and
returns covariate-adjusted estimated marginal means.

## Usage

``` r
anova_ancova(
  data,
  response,
  groups,
  covariates,
  center_covariates = TRUE,
  force_interaction = NULL,
  homogeneity_alpha = 0.05,
  interaction = TRUE,
  type = c("III", "II"),
  conf_level = 0.95,
  vcov_type = "model",
  adjust = "tukey",
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

- covariates:

  Character vector. One or more numeric covariate columns.

- center_covariates:

  Logical. Mean-centre the covariates before fitting. Default `TRUE`.
  See Details.

- force_interaction:

  Logical or `NULL`. `NULL` (default) lets the slopes test decide;
  `TRUE` always fits the covariate-by-group interaction; `FALSE` always
  fits the additive model. See Details for what the test-based choice
  costs.

- homogeneity_alpha:

  Numeric in (0, 1). Significance level for the slopes test. Default
  `0.05`.

- interaction:

  How the grouping variables combine: `TRUE` (the default, their full
  factorial), `FALSE` (additive), or a whole number giving the highest
  order of interaction among them. It has no effect with a single
  grouping variable. The covariate slopes are allowed to differ across
  every grouping term this includes. See Details.

- type:

  Character. `"III"` (default) or `"II"` sums of squares. The model is
  fitted under sum-to-zero contrasts when `"III"`. A model with aliased
  coefficients (an empty cell of the design, or collinear predictors)
  has no Type III tests; Type II tests are then computed instead, and
  `$notes`, the method and the table heading say so.

- conf_level:

  Numeric in (0, 1). Level for every interval returned, including the
  bands of the covariate plot. Default `0.95`.

- vcov_type:

  Character. `"model"` (default) or an HC type (`"HC0"` to `"HC4"`),
  which requires the sandwich package. An HC type reaches `$anova` (as
  Wald F tests), `$emmeans`, `$posthoc` and `$simple_slopes`; see
  Details.

- adjust:

  Character. Multiplicity adjustment for the pairwise comparisons,
  applied by emmeans. One of `"tukey"`, `"sidak"`, `"scheffe"`,
  `"dunnettx"`, `"bonferroni"`, `"holm"`, `"hochberg"`, `"hommel"`,
  `"BH"`, `"BY"`, `"fdr"` or `"none"`. Default `"tukey"`.

- posthoc:

  Logical. Compute pairwise comparisons. There are `choose(k, 2)` of
  them, so this is worth turning off when the number of cells is large;
  `$notes` records that they were skipped. Default `TRUE`.

- plots:

  Logical. Build ggplot2 objects. They are returned in `$plots`, never
  drawn. Default `TRUE`.

- verbose:

  Logical. Emit progress through
  [`message`](https://rdrr.io/r/base/message.html). Default `FALSE`.

## Value

An
[`anovakit_fit`](https://elkronos.github.io/anovakit/reference/anovakit_fit.md)
object, with five extra components: `$slopes_test` (the
homogeneity-of-slopes comparison, also in `$assumptions$slopes`);
`$simple_slopes` (the slope of each covariate in each cell of the
grouping structure, from
[`emtrends`](https://rvlenth.github.io/emmeans/reference/emtrends.html),
with intervals at `conf_level` and the covariance chosen by `vcov_type`;
`NULL` when the additive model is used); `$covariate_means` (the mean of
each covariate over the analysed rows, returned whether or not the
covariates were centred: they are the values subtracted when
`center_covariates = TRUE`, and the point at which `$emmeans` are
evaluated either way); and `$model_additive` and `$model_interaction`,
both candidate fits, so the one that was not chosen is still available.
`$assumptions` also holds `normality` (Shapiro-Wilk on the residuals)
and `levene` (Levene's test on the residuals across the cells).
`$data_used` holds the analysed columns and the cell factor, with the
covariates centred when `center_covariates = TRUE`. The covariate plot
shows each covariate on its original scale.

## Details

**Grouping structure.** With two or more grouping variables the default,
`interaction = TRUE`, fits their full factorial – `y ~ x + A * B` –
which is what a factorial ANCOVA is. `interaction = FALSE` enters them
additively and a whole number keeps interactions up to that order, as in
[`anova_glm`](https://elkronos.github.io/anovakit/reference/anova_glm.md).
The structure chosen is the group part of both candidate models below.
Under the full factorial `$emmeans` and `$posthoc` compare the cells,
and a cell with no data is left out with a note; with the grouping
variables entered additively they are reported for each factor
separately, averaged over the others.

**Covariates are centred.** Every covariate is mean-centred before the
model is fitted. That moves the intercept and, when the
covariate-by-group terms are in the model and `type = "III"`, the group
row of the ANOVA table: that row tests the group difference *at
covariate = 0*, which for an uncentred covariate may be a point nowhere
near the data, and the p-value can be anything. On the example below – a
baseline distributed around 100 – the same data gives p = 0.35 uncentred
against p = 5.5e-21 centred. In the additive model, and under Type II,
the group row is the same either way. `$emmeans`, `$posthoc` and
`$simple_slopes` do not depend on centring at all: they are always
evaluated with every covariate held at its mean, including a covariate
with only two values (a 0/1 indicator, say), which is held at its mean
rather than averaged over its two values. Set
`center_covariates = FALSE` only if the covariate's own zero is the
point you want the intercept (and the Type III group row) to refer to.
When the covariate mean lies outside the range a group was observed
over, that group's adjusted mean is an extrapolation along the fitted
slope, and `$notes` names the group.

A covariate that is constant, that is determined by the grouping
variables (constant within every group, such as a group-level
attribute), or that is an exact linear combination of the grouping
variables and the other covariates is refused with an error naming it:
its effect cannot be separated from the group effect, so there is
nothing to adjust for.

**Homogeneity of slopes.** The two candidate models are the grouping
structure with the covariates entered additively, and the same model
plus the product of every covariate with every grouping term in it, so
that under the full factorial each cell has its own slope. They are
compared by one nested F test, joint over all covariates: heterogeneity
in any one covariate retains the slope terms for all of them.
`homogeneity_alpha` is its level. The test is reported as not computable
(`NA`) when the covariate-by-group terms add nothing estimable, or when
the model without them already fits the data exactly. If the slopes
differ, the interaction model is used, and a note explains that the
group effect is then a comparison at one point on the covariate rather
than a constant difference, and that `$simple_slopes` is the thing to
read.

Choosing the model with a test on the same data has a cost, in both
directions. When the slopes really differ but the test misses the
difference, the adjusted comparisons of the additive model are biased by
the slope difference times the difference in covariate means. That
vanishes on average when the groups are randomised but not otherwise,
and the test often lacks the power to find a slope difference large
enough to matter. In simulations with three groups of 20, slopes 0.5, 1
and 1.5, and covariate means at -1, 0 and 1 within-group SDs, the
default pipeline rejected a true null of equal adjusted means for the
group row about 11–13% of the time at the 5% level (about 34% with 10
per group and means at -2, 0 and 2), against about 5% with
`force_interaction = TRUE`; always fitting the additive model
(`force_interaction = FALSE`) was worse still. Randomised groups were
not affected. So when the additive model is used and the covariate means
differ between groups by more than half a pooled within-group SD,
`$notes` says so. Conversely, when the test retains the interaction, the
p-values do not allow for that choice and can also be too small.
`force_interaction = TRUE` avoids the selection step: with centred
covariates and `type = "III"`, the group row is then a valid comparison
of the groups at the covariate mean whether or not the slopes differ,
which makes it the safer choice when the groups were not randomised.
`$notes` always reports which model was fitted and why.

**Robust standard errors.** With `vcov_type` other than `"model"`, the F
tests in `$anova` are Wald F tests (`car::Anova(vcov. = )`) built from
the heteroscedasticity-consistent covariance, on the model's residual
degrees of freedom, so that table has no sums of squares. `$emmeans`,
`$posthoc` and `$simple_slopes` use the same covariance. `$effect_sizes`
are still computed from the model-based sums of squares, which the
covariance does not change, and the homogeneity-of-slopes test stays the
model-based F test.

## See also

[`anova_welch`](https://elkronos.github.io/anovakit/reference/anova_welch.md)
when there is no covariate,
[`anova_manova`](https://elkronos.github.io/anovakit/reference/anova_manova.md)
for several responses at once.

## Examples

``` r
set.seed(7)
n <- 150
d <- data.frame(
  grp = factor(rep(c("control", "treated"), each = n / 2)),
  baseline = rnorm(n, mean = 100, sd = 5)
)
d$score <- 2 * d$baseline + ifelse(d$grp == "treated", 6, 0) + rnorm(n, 0, 3)

fit <- anova_ancova(d, "score", "grp", "baseline")
fit
#> Analysis of covariance (Type III) 
#> ---------------------------------
#> Call: anova_ancova(data = d, response = "score", groups = "grp", covariates = "baseline")
#> Observations used: 150
#> 
#> Omnibus test (Type III, F)
#>        term sum_sq  df statistic   p_value
#> 1  baseline  12400   1    1214.0 < 2.2e-16
#> 2       grp   1252   1     122.6 < 2.2e-16
#> 3 Residuals   1501 147        NA        NA
#> 
#> Notes
#>   - Covariate(s) mean-centred before fitting (baseline: 100.8), so the intercept refers to the covariate mean rather than to covariate = 0. $emmeans are evaluated at the covariate mean either way.
#>   - Slopes are consistent with being equal across groups (p = 0.6719), so the additive model was used. A non-significant test is not proof of equal slopes, only an absence of evidence against it.
#> 
#> Plots available: residuals, qq, covariate, emmeans
#>   (use plot(x, which = "residuals"))

# Group means adjusted to the covariate mean, not the raw group means
fit$emmeans
#>       grp estimate        se  df conf_low conf_high
#> 1 control 201.6520 0.3693537 147 200.9221  202.3819
#> 2 treated 207.4429 0.3693537 147 206.7129  208.1728
fit$effect_sizes
#>       term df    sum_sq partial_eta_sq partial_omega_sq
#> 1 baseline  1 12398.719      0.8920273        0.8899854
#> 2      grp  1  1252.047      0.4548248        0.4477946

# Why the additive model was chosen, and what centring did
fit$slopes_test
#>                                   comparison df statistic   p_value homogeneous
#> 1 additive vs covariate-by-group interaction  1 0.1801101 0.6719033        TRUE
fit$notes
#> [1] "Covariate(s) mean-centred before fitting (baseline: 100.8), so the intercept refers to the covariate mean rather than to covariate = 0. $emmeans are evaluated at the covariate mean either way."
#> [2] "Slopes are consistent with being equal across groups (p = 0.6719), so the additive model was used. A non-significant test is not proof of equal slopes, only an absence of evidence against it." 

# Fix the model in advance instead of letting the slopes test choose. With
# the interaction in the model, centring is what keeps the group row
# interpretable: uncentred, it tests the groups at baseline = 0.
int <- anova_ancova(d, "score", "grp", "baseline",
                    force_interaction = TRUE, plots = FALSE)
int$simple_slopes
#>   covariate     grp    slope         se  df conf_low conf_high statistic
#> 1  baseline control 1.994464 0.07795951 146 1.840389  2.148539  25.58333
#> 2  baseline treated 1.946220 0.08273512 146 1.782706  2.109733  23.52350
#>        p_value
#> 1 8.207346e-56
#> 2 1.599281e-51
int$anova$p_value[int$anova$term == "grp"]
#> [1] 5.501471e-21
unc <- anova_ancova(d, "score", "grp", "baseline", force_interaction = TRUE,
                    center_covariates = FALSE, plots = FALSE)
unc$anova$p_value[unc$anova$term == "grp"]
#> [1] 0.3544107

# Two grouping variables: their interaction is part of the model
d$site <- factor(rep(c("north", "south"), times = n / 2))
two <- anova_ancova(d, "score", c("grp", "site"), "baseline", plots = FALSE)
two$anova
#>        term       sum_sq  df    statistic      p_value
#> 1  baseline 12331.007064   1 1205.3293507 3.864746e-72
#> 2       grp  1247.927988   1  121.9822699 5.839292e-21
#> 3      site    13.657083   1    1.3349504 2.498268e-01
#> 4  grp:site     3.691966   1    0.3608817 5.489548e-01
#> 5 Residuals  1483.408683 145           NA           NA
```
