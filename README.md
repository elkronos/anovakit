# anovakit

<!-- badges: start -->
[![R-CMD-check](https://github.com/elkronos/anovakit/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/elkronos/anovakit/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Eight analysis-of-variance workflows that share one argument convention, one
return class, and one set of guarantees about what they will and will not do
quietly.

Every function names its columns as character strings, validates its input the
same way, drops incomplete rows and says how many, checks the assumptions that
matter for the method it implements, reports effect sizes, computes estimated
marginal means and pairwise comparisons, and returns `ggplot2` objects rather
than drawing them.

## Installation

```r
# install.packages("remotes")
remotes::install_github("elkronos/anovakit")
```

## Which function

| Function | Response | Use it when |
|---|---|---|
| `anova_welch()` | one numeric | Comparing means. The default choice: it does not assume equal variances. |
| `anova_kw()` | one numeric | No distributional assumption. Kruskal-Wallis with Dunn post-hoc. Convert an ordered factor with `as.integer()` first. |
| `anova_ancova()` | one numeric | Adjusting for one or more numeric covariates. |
| `anova_rm()` | one numeric | The same subjects measured repeatedly. |
| `anova_manova()` | two or more numeric | Several outcomes analysed jointly, with or without covariates. |
| `anova_bin()` | binary | Logistic regression, odds ratios, marginal probabilities. |
| `anova_count()` | counts | Poisson, quasi-Poisson or negative binomial, with exposure offsets. |
| `anova_glm()` | anything | Any other GLM family. |

Six take `data, response, groups`. `anova_manova()` takes `responses` (plural),
and `anova_rm()` takes `response` with `subject`, `within` and `between`.
Everything after that is shared.

## A first example

```r
library(anovakit)

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

fit$effect_sizes
#>   group1 group2   hedges_g    conf_low  conf_high magnitude
#> 1      A      B -1.2374969 -1.78385942 -0.6911345     large
#> 2      A      C  0.5148937  0.00698087  1.0228066    medium
#> 3      B      C  1.9329145  1.32538587  2.5404432     large

summary(fit)     # adds assumptions, effect sizes, marginal means, comparisons
plot(fit, "means")
```

## What every result contains

All eight functions return an object of class `anovakit_fit` — a plain
list, so `$` works as usual, with `print()`, `summary()` and `plot()` methods.

| Component | What is in it |
|---|---|
| `$method` | the analysis that was run |
| `$call` | the matched call |
| `$model` | the fitted model — `lm`, `glm`, `manova` or `afex_aov`; for `anova_welch()` and `anova_kw()`, which fit no model, the `htest` from `oneway.test()` or `kruskal.test()` |
| `$anova` | the omnibus table |
| `$effect_sizes` | Hedges' g, odds ratios, incidence rate ratios, partial eta and omega squared, epsilon squared, generalised eta squared or McFadden's pseudo R², depending on the method. Intervals accompany the standardised mean differences, odds ratios and rate ratios; the variance-explained measures are point estimates |
| `$emmeans` | estimated marginal means with intervals at `conf_level` |
| `$emmeans_object` | the `emmGrid`, for contrasts the wrapper does not cover. `NULL` for `anova_welch()` and `anova_kw()`, which fit no model emmeans can use; `anova_manova()` keeps one per response under `$univariate` |
| `$posthoc` | pairwise comparisons |
| `$assumptions` | the checks that matter for this method; empty for `anova_kw()` unless `diagnostics = TRUE` |
| `$plots` | named `ggplot2` objects; empty when `plots = FALSE` |
| `$data_used` | the rows and columns the model was fitted on |
| `$n_removed` | input rows that are not in `$data_used` — dropped for missing or infinite values, or in `anova_rm()` for an incomplete design or afex aggregation. `nrow($data_used) + $n_removed` is always the number of rows you supplied |
| `$conf_level` | the level used for every interval |
| `$notes` | everything the function decided on your behalf, could not compute, or thinks you should know |

Each function adds components of its own. `names(fit)` lists them, and each is
documented on the function that produces it.

| Function | What it adds |
|---|---|
| `anova_ancova()` | `$slopes_test` (the homogeneity-of-slopes comparison), `$simple_slopes` (per-group slopes, when the interaction is fitted), `$covariate_means`, and both candidate fits in `$model_additive` / `$model_interaction` |
| `anova_rm()` | `$sphericity` (Mauchly plus both epsilon corrections), `$subjects_dropped`, `$residuals`, and `$n_removed` broken out into `$n_removed_missing` / `$n_removed_unbalanced` / `$n_removed_aggregated` |
| `anova_manova()` | `$multivariate`, `$univariate` (one full fit per response), `$canonical`, `$canonical_term`, `$test` |
| `anova_count()` | `$model_type`, `$dispersion` (of the Poisson fit the choice was made on) and `$model_dispersion` (of the model returned) |
| `anova_bin()` | `$model_stats` — AIC, BIC, log-likelihood, McFadden's R², the level treated as a success, and `$lr_vs_null`, the model-versus-null test |
| `anova_glm()` | `$family` and `$model_stats` — AIC, BIC, the null and residual deviances, the residual degrees of freedom |

**Read `$notes`.** It is where the package tells you that it switched from
Poisson to negative binomial, that it centred your covariate, that two subjects
were dropped for an incomplete design, or that an odds ratio is not identified
because a cell is perfectly separated.

## Design decisions worth knowing about

**Nothing prints.** No analysis function writes to the console, and no plot is
drawn as a side effect. Warnings and messages raised by `car`, `glm()` and
`afex` are captured into `$notes` rather than escaping. `verbose = TRUE` emits
`message()`, which `suppressMessages()` can silence.

**`conf_level` means every interval.** Group means, pairwise differences,
standardised effect sizes, odds ratios, marginal means: one argument sets them
all in the tables the function returns. The exception is `$emmeans_object`,
which is emmeans' own grid and carries emmeans' default level — pass `level =`
yourself when you summarise it:

```r
d$covariate <- rnorm(90, 50, 10)
anc <- anova_ancova(d, "value", "group", "covariate", plots = FALSE)
emmeans::contrast(anc$emmeans_object, "trt.vs.ctrl", level = 0.90)
#>  contrast estimate    SE df t.ratio p.value
#>  B - A       2.541 0.472 86   5.379  <.0001
#>  C - A      -0.871 0.474 86  -1.835  0.1284
```

**Type III means Type III.** When you ask for `type = "III"` the model is
*fitted* under sum-to-zero contrasts. Setting the global contrast option after a
model exists does not change the contrasts stored on it, and `car::Anova(type = 3)`
would otherwise report simple effects at the reference level under the wrong
label.

```r
d$site <- rep(c("north", "south"), 45)
tw <- anova_glm(d, "value", c("group", "site"), interaction = TRUE, type = "III")
unlist(tw$model$contrasts)
#>     group      site
#> "contr.sum" "contr.sum"
```

**Covariates are centred.** `anova_ancova()` mean-centres covariates by default.
This matters most when the model contains a group-by-covariate term: that Type
III group row tests the groups *at covariate = 0*, which on a baseline
distributed around 100 is a point nowhere near the data. On the trial data in
the vignette that is p = 0.996 uncentred against p = 6.5e-11 centred, on the
same rows. In the additive model the group row is unaffected by centring; what
it buys there is that `$emmeans` and the intercept are read at the covariate
mean. `force_interaction` fixes the model in advance instead of letting a test
on the same data choose it.

**Empty cells are dropped.** When several grouping variables are combined, level
combinations that never occur are removed, so the number of groups reported is
the number that contain data.

**No escape hatches.** No modelling function has `...` except `anova_rm()`,
which forwards it to `afex::aov_ez()` and documents that it does. Every argument
that reaches a fitting call is named, so nothing can be passed through that
silently overrides a decision the function has already made — such as
`contrasts=`, which would defeat `type = "III"` while the output still said
"Type III". Prior weights are given as a column name, so they are subsetted with
the data rather than misaligning the moment a row is dropped:

```r
d$wt <- rep(c(1, 2), 45)
anova_glm(d, "value", "group", weights = "wt", plots = FALSE)$anova
```

**Assumption checks match the method.** `anova_welch()` tests normality within
each group rather than on pooled residuals, because pooling groups with
different variances produces a mixture that fails normality tests even when
every group is perfectly normal — which would argue against the very method
being used. `anova_kw()` assumes no distribution, so it runs no checks at all
unless you ask with `diagnostics = TRUE`.

**Robust standard errors are one argument.** `vcov_type = "HC0"` through
`"HC4"` on `anova_ancova()`, `anova_bin()`, `anova_count()` and `anova_glm()`
flows through to the marginal means and the comparisons, not only to the
coefficient table.

**Multiplicity adjustments differ by method.** `anova_welch()` and `anova_kw()`
compare directly and take the `p.adjust()` methods, defaulting to `"holm"` and
`"BH"`. The other six go through emmeans, which additionally offers `"tukey"`
(their default), `"sidak"`, `"scheffe"` and `"dunnettx"`.

## Arguments beyond the first three

Everything below is optional and every function that can support an argument
does, with the same name and the same meaning.

| Argument | Where | What it does |
|---|---|---|
| `conf_level` | all eight | The level for every interval. Default `0.95`. |
| `adjust` | all eight | Multiplicity adjustment for the comparisons. |
| `posthoc` | all eight | Skip the `choose(k, 2)` comparisons. |
| `plots` | all eight | Build the `ggplot2` objects. |
| `verbose` | all eight | Progress through `message()`. |
| `type` | ancova, manova, bin, count, glm | `"II"` or `"III"` sums of squares; `"III"` refits under sum-to-zero contrasts. |
| `interaction` | manova, bin, count, glm | `FALSE` (additive), `TRUE` (full factorial), or the highest order to include. |
| `weights` | bin, count, glm | Prior weights, named as a column. |
| `vcov_type` | ancova, bin, count, glm | `"HC0"`–`"HC4"` robust standard errors, which reach the marginal means too. |
| `ci_method` | bin, glm | `"profile"` (default) or `"wald"` coefficient intervals. |
| `test_statistic` | count, glm | The statistic `car::Anova()` uses. |
| `center_covariates`, `force_interaction`, `homogeneity_alpha` | ancova | Control the two decisions the function would otherwise make itself. |
| `covariates` | ancova, manova | Numeric covariates to adjust for. |
| `offset` | count | A strictly positive exposure column, entered as a log offset. |
| `model`, `overdispersion_threshold` | count | Choose the count family yourself, or move the threshold that chooses it. |
| `reference`, `success` | bin | Which group is the baseline, and which outcome level is a success. |
| `diagnostics` | kw | Optional normality diagnostics, off by default. |
| `hedges_correction` | welch | Hedges' g (default) or Cohen's d. |
| `subject`, `within`, `between`, `factorize`, `emm_specs`, `correction` | rm | The design, and which sphericity correction the table uses. |
| `test` | manova | `"Pillai"` (default), `"Wilks"`, `"Hotelling-Lawley"` or `"Roy"`. |
| `assumptions` | manova | Compute Mardia's tests and Box's M. |
| `family` | glm | Any `stats::family`, as an object, a function or a string. |

## Getting help

```r
vignette("anovakit")      # the whole package on one dataset
?anova_welch                  # any function
?anovakit_fit             # what the returned object holds
citation("anovakit")
```

Bug reports and feature requests:
<https://github.com/elkronos/anovakit/issues>.

## Dependencies

Imports `car`, `emmeans`, `ggplot2` and base R only.

`afex` (repeated measures), `MASS` (negative binomial) and `sandwich` (robust
standard errors) are Suggests. `anova_rm()` fails with a clear message without
`afex`, as does `anova_count(model = "negbin")` without `MASS`. Where a fallback
exists — `model = "auto"` without `MASS`, `vcov_type = "HC3"` without `sandwich`
— the function degrades to the nearest available method and records it in
`$notes`. The remaining Suggests (`knitr`, `rmarkdown`, `testthat`, `withr`,
`data.table`, `tibble`, `effectsize`) are for building and testing only.

Dunn's test, Mardia's tests of multivariate skewness and kurtosis, Box's M and
the canonical discriminant analysis are implemented directly rather than pulled
in as dependencies.

## Licence

GPL-3.
