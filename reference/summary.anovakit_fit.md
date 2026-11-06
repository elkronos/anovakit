# Summarise an anovakit result

Everything `print` shows, plus assumption checks, effect sizes, marginal
means, post-hoc comparisons and whatever further tables the method
produces (simple slopes, sphericity, the multivariate tests and
univariate follow-ups, the canonical axes, and the model fit statistics
in `$model_stats`, including `anova_bin`'s model-versus-null test).

## Usage

``` r
# S3 method for class 'anovakit_fit'
summary(object, digits = 4L, ...)

# S3 method for class 'summary.anovakit_fit'
print(x, digits = x$digits, ...)
```

## Arguments

- object:

  An
  [`anovakit_fit`](https://elkronos.github.io/anovakit/reference/anovakit_fit.md)
  object.

- digits:

  Number of significant digits for the printed tables. Default `4` for
  [`summary()`](https://rdrr.io/r/base/summary.html); the
  [`print()`](https://rdrr.io/r/base/print.html) method uses the value
  stored in the summary object.

- ...:

  Ignored.

- x:

  A `summary.anovakit_fit` object.

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `summary.anovakit_fit`, a list holding `fit` (the original
`object`) and `digits`. Its
[`print()`](https://rdrr.io/r/base/print.html) method prints the summary
and returns `x` invisibly.

## Examples

``` r
set.seed(1)
d <- data.frame(g = rep(c("a", "b", "c"), each = 20), y = rnorm(60))
summary(anova_welch(d, "y", "g", plots = FALSE))
#> Welch's analysis of variance 
#> ----------------------------
#> Call: anova_welch(data = d, response = "y", groups = "g", plots = FALSE)
#> Observations used: 60
#> 
#> Omnibus test (Welch F)
#>   term statistic num_df den_df p_value
#> 1    g     0.264      2  37.91  0.7694
#> 
#> Assumption checks
#> 
#>   normality
#>  group  n statistic p_value note
#>      a 20    0.9533  0.4195 <NA>
#>      b 20    0.9462  0.3127 <NA>
#>      c 20    0.9691  0.7352 <NA>
#> 
#>   variance_ratio
#>     1.272
#> 
#> Effect sizes
#>  group1 group2 hedges_g conf_low conf_high  magnitude            standardiser
#>       a      b   0.2163  -0.4025     0.841      small sqrt((s1^2 + s2^2) / 2)
#>       a      c  0.05873  -0.5604    0.6795 negligible sqrt((s1^2 + s2^2) / 2)
#>       b      c  -0.1693  -0.7926    0.4494 negligible sqrt((s1^2 + s2^2) / 2)
#> 
#> Group summaries (95% intervals)
#>  group  n      mean     sd     se conf_low conf_high
#>      a 20    0.1905 0.9133 0.2042  -0.2369    0.6179
#>      b 20 -0.006472 0.8714 0.1948  -0.4143    0.4013
#>      c 20    0.1388 0.8097 0.1811  -0.2402    0.5178
#> 
#> Pairwise comparisons
#>  group1 group2 difference conf_low conf_high statistic    df p_value p_adjusted
#>       a      b      0.197  -0.3744    0.7684    0.6979 37.92  0.4895          1
#>       a      c    0.05173   -0.501    0.6045    0.1895 37.46  0.8507          1
#>       b      c    -0.1453  -0.6838    0.3933   -0.5462  37.8  0.5882          1
#>  adjustment
#>        holm
#>        holm
#>        holm
```
