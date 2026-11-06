# anovakit

<!-- badges: start -->
[![R-CMD-check](https://github.com/elkronos/anovakit/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/elkronos/anovakit/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/elkronos/anovakit/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/elkronos/anovakit/actions/workflows/pkgdown.yaml)
<!-- badges: end -->

anovakit has eight analysis-of-variance functions: Welch's ANOVA,
Kruskal-Wallis, ANCOVA, repeated measures, MANOVA, and analyses of deviance
for binary, count and other GLM outcomes. They take the same arguments and
return the same kind of object. Each one checks the assumptions its method
needs, reports effect sizes, marginal means and pairwise comparisons, builds
`ggplot2` figures without drawing them, and records in `$notes` every
decision it made for you.

## Installation

```r
# install.packages("remotes")
remotes::install_github("elkronos/anovakit", build_vignettes = TRUE)
```

`install_github()` skips vignettes unless you pass `build_vignettes = TRUE`,
so without it `vignette("anovakit")` is not available.

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
| `anova_glm()` | anything | A classical factorial ANOVA with main effects and interactions, Type II or III (the default Gaussian family), or any other GLM family. |

Six take `data, response, groups`, with columns named as character strings.
`anova_manova()` takes `responses` (plural), and `anova_rm()` takes `response`
with `subject`, `within` and `between`. [Choose an
analysis](https://elkronos.github.io/anovakit/articles/choosing.html) walks
through the choice.

## A first example

Fuel economy of 32 cars, by number of cylinders:

```r
library(anovakit)

fit <- anova_welch(mtcars, "mpg", "cyl")
fit
#> Welch's analysis of variance 
#> ----------------------------
#> Call: anova_welch(data = mtcars, response = "mpg", groups = "cyl")
#> Observations used: 32
#> 
#> Omnibus test (Welch F)
#>   term statistic num_df den_df   p_value
#> 1  cyl     31.62      2  18.03 1.271e-06
#> 
#> Notes
#>   - Largest group variance is 9.63 times the smallest. Welch's test handles this; the classical F test would not.
#> 
#> Plots available: means, box, qq
#>   (use plot(x, which = "means"))

fit$posthoc[, c("group1", "group2", "difference", "conf_low", "conf_high", "p_adjusted")]
#>   group1 group2 difference conf_low conf_high   p_adjusted
#> 1      4      6   6.920779 3.751376 10.090182 4.048495e-04
#> 2      4      8  11.563636 8.318518 14.808755 4.924043e-06
#> 3      6      8   4.642857 2.802925  6.482789 9.080709e-05

summary(fit)       # adds the assumption checks, effect sizes and group summaries
plot(fit, "means") # a ggplot2 object
```

The comparisons are Welch t tests with Holm-adjusted p-values. `cyl` is
numeric in `mtcars`; grouping columns are converted to factors.

## What you get

- The omnibus test in `$anova`, with the sums of squares, statistic and any
  correction named in the printed heading.
- Effect sizes in `$effect_sizes`, with intervals for standardised mean
  differences, odds ratios and rate ratios.
- Marginal means in `$emmeans` and pairwise comparisons in `$posthoc`, with
  both the unadjusted and the adjusted p-value.
- The assumption checks that fit the method in `$assumptions`.
- `ggplot2` objects in `$plots`, returned and never drawn.
- One `conf_level` for every interval, and no console output from any
  analysis function.
- `$notes`: dropped rows, centred covariates, a switch from Poisson to
  negative binomial, a perfectly separated cell.

## Documentation

The website, <https://elkronos.github.io/anovakit/>, has:

- [Get started](https://elkronos.github.io/anovakit/articles/anovakit.html):
  one analysis from start to finish, and the arguments every function shares.
- [Choose an analysis](https://elkronos.github.io/anovakit/articles/choosing.html):
  four questions that lead to one of the eight functions.
- [Walkthroughs and guides](https://elkronos.github.io/anovakit/articles/index.html):
  one worked analysis per function, and guides to results, diagnostics,
  plots, weights and the design choices behind the defaults.
- The [function reference](https://elkronos.github.io/anovakit/reference/index.html).

From R:

```r
vignette("anovakit")      # Get started
?anova_welch              # any function
?anovakit_fit             # what the returned object holds
citation("anovakit")
```

Report bugs and request features at
<https://github.com/elkronos/anovakit/issues>.

## Dependencies

Imports `car`, `emmeans`, `ggplot2` (3.4.0 or later) and base R.

`afex` (repeated measures), `MASS` (negative binomial) and `sandwich` (robust
standard errors) are Suggests. Without `afex`, `anova_rm()` stops with a
message, as does `anova_count(model = "negbin")` without `MASS`. Where a
fallback exists, the function uses it and says so in `$notes`:
`model = "auto"` without `MASS` fits quasi-Poisson where it would have chosen
the negative binomial, and `vcov_type` without `sandwich` uses model-based
standard errors. `bit64` is needed only to read `integer64` columns.

Dunn's test, Mardia's tests, Box's M and the canonical discriminant analysis
are implemented in the package rather than imported.

## Licence

GPL-3.
