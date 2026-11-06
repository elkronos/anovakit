# Design choices

``` r

library(anovakit)
```

Every anovakit function follows the same rules about what it prints,
which settings can change a result, and what it decides for you. Each is
shown below on a small example; for the object itself see [Working with
results](https://elkronos.github.io/anovakit/articles/results.md).

## Nothing prints

No analysis function writes to the console, and no plot is drawn as a
side effect. Warnings and messages raised while fitting by car, emmeans,
sandwich, afex and base R’s model functions are captured and recorded in
`$notes` instead of reaching the console. Here a logistic regression
with separation in one group prints nothing, and the second of its notes
describes the problem:

``` r

sep <- data.frame(
  group = rep(c("a", "b", "c"), each = 10),
  y     = c(rep(0, 10), rep(c(0, 1), 5), rep(c(0, 1, 1, 1, 1), 2))
)
out <- capture.output(fit_sep <- anova_bin(sep, "y", "group", plots = FALSE))
length(out)
#> [1] 0
writeLines(strwrap(fit_sep$notes[2], 72))
#> Complete or quasi-complete separation detected: the fitted probability
#> is numerically 0 or 1 in group = a. The affected coefficient(s):
#> groupb, groupc. An affected odds ratio is not identified: it will be
#> enormous, its interval will be unbounded on one side, and its Wald
#> p-value will be near 1 no matter how strong the association is; the
#> same goes for that cell's marginal probability and the comparisons
#> involving it. With events this sparse the likelihood-ratio omnibus test
#> is also liberal. Consider a penalised fit such as logistf::logistf(),
#> or collapsing the offending level.
```

`verbose = TRUE` reports progress through
[`message()`](https://rdrr.io/r/base/message.html), so
[`suppressMessages()`](https://rdrr.io/r/base/message.html) silences it.
Plots are built and stored in `$plots`; they appear only when you print
one.

## `conf_level` sets every interval

One argument sets the level of every interval in the tables a function
returns: group and marginal means, pairwise differences, effect sizes,
odds and rate ratios, simple slopes and the error bars in the plots.

``` r

fit90 <- anova_glm(chickwts, "weight", "feed", conf_level = 0.90,
                   plots = FALSE)
fit90$conf_level
#> [1] 0.9
head(fit90$posthoc[, c("contrast", "estimate", "conf_low", "conf_high")], 3)
#>             contrast  estimate  conf_low conf_high
#> 1 casein - horsebean 163.38333 101.15257  225.6141
#> 2   casein - linseed 104.83333  45.49863  164.1680
#> 3  casein - meatmeal  46.67424 -13.99399  107.3425
```

It does not reach `$emmeans_object`. That is emmeans’ own reference
grid, and it carries emmeans’ default level of 0.95. Pass the level to
the function that summarises the grid, such as
[`confint()`](https://rdrr.io/r/stats/confint.html); given to
`contrast()` it has no effect on the intervals:

``` r

ctl <- emmeans::contrast(fit90$emmeans_object, "trt.vs.ctrl")
confint(ctl, level = fit90$conf_level)
#>  contrast           estimate   SE df lower.CL upper.CL
#>  horsebean - casein  -163.38 23.5 65   -217.3  -109.47
#>  linseed - casein    -104.83 22.4 65   -156.2   -53.43
#>  meatmeal - casein    -46.67 22.9 65    -99.2     5.88
#>  soybean - casein     -77.15 21.6 65   -126.7   -27.62
#>  sunflower - casein     5.33 22.4 65    -46.1    56.74
#> 
#> Degrees-of-freedom method: user-specified 
#> Confidence level used: 0.9 
#> Conf-level adjustment: dunnettx method for 5 estimates
```

## Type II by default, and what Type III means

[`anova_manova()`](https://elkronos.github.io/anovakit/reference/anova_manova.md),
[`anova_bin()`](https://elkronos.github.io/anovakit/reference/anova_bin.md),
[`anova_count()`](https://elkronos.github.io/anovakit/reference/anova_count.md)
and
[`anova_glm()`](https://elkronos.github.io/anovakit/reference/anova_glm.md)
report Type II tests by default. Each main effect is tested after the
other main effects, and the result does not depend on how the factors
are coded.
[`anova_ancova()`](https://elkronos.github.io/anovakit/reference/anova_ancova.md)
and
[`anova_rm()`](https://elkronos.github.io/anovakit/reference/anova_rm.md)
default to Type III; the
[ANCOVA](https://elkronos.github.io/anovakit/articles/ancova.md) and
[repeated
measures](https://elkronos.github.io/anovakit/articles/repeated-measures.md)
walkthroughs explain their tables. The type used is printed in the
heading and stored in `attr(fit$anova, "ss_type")`.

When you ask for `type = "III"`, the model is fitted under sum-to-zero
contrasts. In a model with an interaction, a Type III main-effect row
tests the main effect only under those contrasts. Under R’s default
treatment contrasts, `car::Anova(type = 3)` tests the simple effect of
each factor at the reference level of the other, and still labels the
table Type III.

The `warpbreaks` data with the first five rows removed make an
unbalanced two-way design:

``` r

wb <- warpbreaks[-(1:5), ]
t3 <- anova_glm(wb, "breaks", c("wool", "tension"), interaction = TRUE,
                type = "III", plots = FALSE)
unlist(t3$model$contrasts)
#>        wool     tension 
#> "contr.sum" "contr.sum"
t3$anova
#>           term   sum_sq df statistic      p_value
#> 1         wool  588.751  1  6.343522 0.0155817964
#> 2      tension 2002.165  2 10.786206 0.0001597943
#> 3 wool:tension 1121.446  2  6.041534 0.0048717187
#> 4    Residuals 3990.889 43        NA           NA
```

The same model fitted with treatment contrasts and passed to car gives a
different `wool` row, which is the effect of wool at `tension = "L"`
only:

``` r

car::Anova(lm(breaks ~ wool * tension, data = wb), type = 3)["wool", ]
#> Anova Table (Type III tests)
#> 
#> Response: breaks
#>      Sum Sq Df F value    Pr(>F)    
#> wool 1195.5  1  12.881 0.0008449 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

A model with aliased coefficients, such as one with an empty cell, has
no Type III tests. The function then computes Type II tests, and the
heading, `$method` and `$notes` say so.

## Covariates are centred

[`anova_ancova()`](https://elkronos.github.io/anovakit/reference/anova_ancova.md)
mean-centres every covariate before fitting
(`center_covariates = TRUE`), and
[`anova_manova()`](https://elkronos.github.io/anovakit/reference/anova_manova.md)
centres its covariates too. Centring always moves the intercept to the
covariate mean. When the model contains a group-by-covariate term and
`type = "III"`, it also moves the group row, because that row then tests
the groups at covariate = 0. On a baseline score distributed around 100,
zero is far from the data:

``` r

set.seed(2024)
n <- 240
trial <- data.frame(
  arm      = factor(rep(c("placebo", "low", "high"), each = n / 3)),
  baseline = rnorm(n, mean = 100, sd = 12)
)
effect <- c(placebo = 0, low = 3, high = 7)[as.character(trial$arm)]
trial$score <- 0.6 * trial$baseline + effect + rnorm(n, 0, 6)

centred <- anova_ancova(trial, "score", "arm", "baseline",
                        force_interaction = TRUE, plots = FALSE)
uncentred <- anova_ancova(trial, "score", "arm", "baseline",
                          force_interaction = TRUE,
                          center_covariates = FALSE, plots = FALSE)
data.frame(term        = centred$anova$term,
           p_centred   = signif(centred$anova$p_value, 3),
           p_uncentred = signif(uncentred$anova$p_value, 3))
#>           term p_centred p_uncentred
#> 1     baseline  3.80e-52    3.80e-52
#> 2          arm  6.53e-11    9.96e-01
#> 3 baseline:arm  6.91e-01    6.91e-01
#> 4    Residuals        NA          NA
all.equal(centred$emmeans, uncentred$emmeans)
#> [1] TRUE
```

On the same rows, the `arm` row has p = 6.5e-11 centred and p = 0.996
uncentred. The other rows and the marginal means are the same, because
`$emmeans`, `$posthoc` and `$simple_slopes` are always evaluated with
every covariate at its mean (a 0/1 covariate included). In the additive
model, and under Type II, centring does not change the group row. The
means subtracted are in `$covariate_means`.

## The ANCOVA model is chosen by a test unless you fix it

By default
[`anova_ancova()`](https://elkronos.github.io/anovakit/reference/anova_ancova.md)
chooses between the additive and the separate-slopes model with a
homogeneity-of-slopes test on the same data, and records the choice in
`$notes`. The p-values that follow do not allow for that selection. If
the test misses a real slope difference, the group row and `$posthoc`
reject more often than their nominal level: somewhat so when the groups
were randomised, and far more so when the groups’ covariate means
differ, which `$notes` points out. `force_interaction = TRUE` fixes the
separate-slopes model in advance and avoids the selection step. With
several grouping factors,
[`anova_ancova()`](https://elkronos.github.io/anovakit/reference/anova_ancova.md)
fits their interaction by default (`interaction = TRUE`), as a factorial
ANCOVA does. The [ANCOVA
walkthrough](https://elkronos.github.io/anovakit/articles/ancova.md)
covers both choices in detail.

## Comparisons follow the model

When grouping factors enter a model additively (`interaction = FALSE`,
the default except in
[`anova_ancova()`](https://elkronos.github.io/anovakit/reference/anova_ancova.md)),
`$emmeans` and `$posthoc` are reported for each factor on its own,
averaged over the others, with a `term` column and one multiplicity
family per factor. Comparing every cell would repeat each main-effect
difference once per level of the other factors, and adjust for
comparisons that carry no extra information.

``` r

additive <- anova_glm(warpbreaks, "breaks", c("wool", "tension"),
                      plots = FALSE)
additive$emmeans[, c("term", "wool", "tension", "estimate")]
#>      term wool tension estimate
#> 1    wool    A    <NA> 31.03704
#> 2    wool    B    <NA> 25.25926
#> 3 tension <NA>       L 36.38889
#> 4 tension <NA>       M 26.38889
#> 5 tension <NA>       H 21.66667
additive$posthoc[, c("term", "contrast", "p_adjusted")]
#>      term contrast  p_adjusted
#> 1    wool    A - B 0.073613669
#> 2 tension    L - M 0.033626219
#> 3 tension    L - H 0.001121788
#> 4 tension    M - H 0.447421021
```

When the model contains the interaction, the cells are compared:

``` r

cells <- anova_glm(warpbreaks, "breaks", c("wool", "tension"),
                   interaction = TRUE, plots = FALSE)
nrow(cells$emmeans)
#> [1] 6
nrow(cells$posthoc)
#> [1] 15
```

A level combination with no data is left out when the model cannot
estimate it, and named in `$notes` when the model estimates it by
extrapolation (as an additive model does). With no model to follow,
[`anova_welch()`](https://elkronos.github.io/anovakit/reference/anova_welch.md)
and
[`anova_kw()`](https://elkronos.github.io/anovakit/reference/anova_kw.md)
combine the grouping columns into one factor of the cells that contain
data and compare those.

## Every argument is named

No modelling function has a `...` argument except
[`anova_rm()`](https://elkronos.github.io/anovakit/reference/anova_rm.md),
and that one accepts only `fun_aggregate`, `observed`, `type` and
`anova_table = list(p_adjust_method = )`, passed to
[`afex::aov_ez()`](https://rdrr.io/pkg/afex/man/aov_car.html). Anything
else is refused by name rather than ignored:

``` r

rm_data <- data.frame(id = rep(1:6, each = 2), time = rep(c("pre", "post"), 6),
                      y = c(5, 7, 4, 6, 6, 9, 5, 5, 7, 8, 4, 7))
try(anova_rm(rm_data, "y", "id", "time", es = "pes"))
#> Registered S3 method overwritten by 'lme4':
#>   method           from
#>   na.action.merMod car
#> Error : anova_rm() has no argument `es`, so it would be ignored. The arguments passed on to afex::aov_ez() through `...` may only be: fun_aggregate, observed, type, anova_table.
```

Every argument that reaches a fitting call is named, so nothing can be
passed through that overrides a decision the function has made. A
`contrasts =` argument, for example, would defeat `type = "III"` while
the table still said Type III, so it is an unused argument:

``` r

try(anova_glm(warpbreaks, "breaks", "tension", contrasts = list()))
#> Error in anova_glm(warpbreaks, "breaks", "tension", contrasts = list()) : 
#>   unused argument (contrasts = list())
```

A column may play only one role (response, group, covariate, weights,
offset or subject):

``` r

try(anova_glm(warpbreaks, "breaks", "tension", weights = "breaks"))
#> Error : Column `breaks` is given as both the response and the weights; a column can play only one role.
```

Weights are named as a column too, so they stay aligned with the rows;
[Weights and robust standard
errors](https://elkronos.github.io/anovakit/articles/weights.html#zero-weights)
explains what happens to zero weights.

## Global settings do not change results

Options set for some other analysis in the same session do not change a
fit’s results.
[`afex::afex_options()`](https://rdrr.io/pkg/afex/man/afex_options.html)
cannot change an
[`anova_rm()`](https://elkronos.github.io/anovakit/reference/anova_rm.md)
fit, because the type, the contrasts and the emmeans model are all set
explicitly. The
[`emmeans::emm_options()`](https://rvlenth.github.io/emmeans/reference/emm_options.html)
settings that change a result (adjustment, sidedness, null, degrees of
freedom, scale) are cleared for the duration of the fit; only
computational limits such as `rg.limit` are kept, so you can still let a
large grid through. Here an `emm_options()` setting that turns off every
adjustment leaves the Tukey adjustment in place:

``` r

emmeans::emm_options(contrast = list(adjust = "none"))
tension <- anova_glm(warpbreaks, "breaks", "tension", plots = FALSE)
tension$posthoc[, c("contrast", "p_value", "p_adjusted", "adjustment")]
#>   contrast      p_value  p_adjusted adjustment
#> 1    L - M 0.0147169627 0.038459768      tukey
#> 2    L - H 0.0005008605 0.001431503      tukey
#> 3    M - H 0.2386143753 0.463083097      tukey
emmeans::emm_options(contrast = NULL)
```

The global `contrasts` option can change how a Type II model’s
coefficients are coded, but not its tests, marginal means or
comparisons; `type = "III"` always fits under sum-to-zero contrasts, and
odds and rate ratios always come from a treatment-coded fit.

The only option of its own that anovakit reads is `anovakit.theme`,
which sets the base theme of every plot and changes no result; the
[plots guide](https://elkronos.github.io/anovakit/articles/visuals.md)
shows how to use it.

## Assumption checks match the method

Each function checks the assumptions its own method makes.
[`anova_welch()`](https://elkronos.github.io/anovakit/reference/anova_welch.md)
tests normality within each group rather than on pooled residuals.
Pooling groups with different variances produces a mixture that can fail
a normality test when every group is normal, which would argue against
the method being used. For the same reason
[`anova_rm()`](https://elkronos.github.io/anovakit/reference/anova_rm.md)
tests the within-subject residuals in each within-subject cell.

``` r

welch <- anova_welch(chickwts, "weight", "feed", plots = FALSE)
welch$assumptions$normality[, c("group", "n", "statistic", "p_value")]
#>       group  n statistic   p_value
#> 1    casein 12 0.9166257 0.2591841
#> 2 horsebean 10 0.9375834 0.5264499
#> 3   linseed 12 0.9693128 0.9034734
#> 4  meatmeal 11 0.9791381 0.9611795
#> 5   soybean 14 0.9464029 0.5063768
#> 6 sunflower 12 0.9280884 0.3602904
```

[`anova_kw()`](https://elkronos.github.io/anovakit/reference/anova_kw.md)
assumes no distribution, so it runs no checks unless you ask with
`diagnostics = TRUE`:

``` r

length(anova_kw(chickwts, "weight", "feed", plots = FALSE)$assumptions)
#> [1] 0
names(anova_kw(chickwts, "weight", "feed", plots = FALSE,
               diagnostics = TRUE)$assumptions)
#> [1] "normality"
```

The [diagnostics
guide](https://elkronos.github.io/anovakit/articles/diagnostics.md)
lists every check by function and what to do when one fails.

## Robust standard errors

`vcov_type = "HC0"` to `"HC4"` on
[`anova_ancova()`](https://elkronos.github.io/anovakit/reference/anova_ancova.md),
[`anova_bin()`](https://elkronos.github.io/anovakit/reference/anova_bin.md),
[`anova_count()`](https://elkronos.github.io/anovakit/reference/anova_count.md)
and
[`anova_glm()`](https://elkronos.github.io/anovakit/reference/anova_glm.md)
changes the standard errors in the coefficient table and in the marginal
means and comparisons. In
[`anova_ancova()`](https://elkronos.github.io/anovakit/reference/anova_ancova.md)
it also reaches the F tests (as robust Wald F tests) and the simple
slopes, and in the other three the omnibus table when
`test_statistic = "Wald"`. Likelihood-ratio tests rest on the model.
Where fitted values sit on the boundary (separation, an all-zero group)
or a cell has one observation, the model-based covariance is used, with
a note. [Weights and robust standard
errors](https://elkronos.github.io/anovakit/articles/weights.md) gives
the rules, including how whole-number weights are treated as frequency
weights.

## Multiplicity adjustments depend on the method

[`anova_welch()`](https://elkronos.github.io/anovakit/reference/anova_welch.md)
and
[`anova_kw()`](https://elkronos.github.io/anovakit/reference/anova_kw.md)
compare groups directly and accept the
[`p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) methods,
defaulting to `"holm"` and `"BH"`. The other six take their comparisons
from emmeans, which also offers `"tukey"` (their default), `"sidak"`,
`"scheffe"` and `"dunnettx"`. Every `$posthoc` table reports the
unadjusted p-value in `p_value`, the adjusted one in `p_adjusted` and
the method in `adjustment`, as the table above shows.

The intervals in `$posthoc` from the emmeans-based functions are
simultaneous. emmeans cannot turn a step-down method such as Holm’s into
intervals and uses Bonferroni intervals instead; `$notes` says when that
happens.
[`anova_welch()`](https://elkronos.github.io/anovakit/reference/anova_welch.md)
reports a per-comparison interval for each pair.

## Notes record every decision

Whenever a function chooses something on your behalf, it says so in
`$notes`: the model chosen by a slopes test or a dispersion statistic, a
centred covariate, dropped rows or subjects, a fallback from Type III to
Type II, an odds ratio that is not identified. The ANCOVA above, fitted
without `force_interaction`, records two decisions:

``` r

writeLines(anova_ancova(trial, "score", "arm", "baseline", plots = FALSE)$notes)
#> Covariate(s) mean-centred before fitting (baseline: 100.3), so the intercept refers to the covariate mean rather than to covariate = 0. $emmeans are evaluated at the covariate mean either way.
#> Slopes are consistent with being equal across groups (p = 0.6908), so the additive model was used. A non-significant test is not proof of equal slopes, only an absence of evidence against it, and the p-values below do not allow for the model having been chosen by a test on these same data: if the slopes do differ, the group row and $posthoc can reject more often than their nominal level, somewhat so even when the groups were randomised and far more so when the groups' covariate means are unequal. Setting force_interaction = TRUE in advance avoids the selection step.
```

Missing optional packages are recorded the same way. Without MASS,
`anova_count(model = "auto")` fits a quasi-Poisson model where it would
have chosen a negative binomial, and without sandwich a `vcov_type`
request falls back to model-based standard errors; each records the
fallback in `$notes`.
[`anova_rm()`](https://elkronos.github.io/anovakit/reference/anova_rm.md)
without afex, and `anova_count(model = "negbin")` without MASS, stop
with an error instead.

[Working with
results](https://elkronos.github.io/anovakit/articles/results.html#the-notes)
shows how to read and search the notes, and how `$data_used` and
`$n_removed` account for every row you supplied.

## See also

- [Working with
  results](https://elkronos.github.io/anovakit/articles/results.md): the
  object every function returns.
- [Weights and robust standard
  errors](https://elkronos.github.io/anovakit/articles/weights.md).
- [Diagnostics](https://elkronos.github.io/anovakit/articles/diagnostics.md):
  the assumption checks and what to do when one fails.
- [Choose an
  analysis](https://elkronos.github.io/anovakit/articles/choosing.md):
  which function fits which design.
- [`?anovakit`](https://elkronos.github.io/anovakit/reference/anovakit-package.md)
  lists the shared conventions, and
  [`?anovakit_fit`](https://elkronos.github.io/anovakit/reference/anovakit_fit.md)
  the components of a fit.
