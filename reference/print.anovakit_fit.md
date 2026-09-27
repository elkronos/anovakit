# Print an anovakit result

Shows the method, the sample size, the omnibus table and any notes. Use
[`summary.anovakit_fit`](https://elkronos.github.io/anovakit/reference/summary.anovakit_fit.md)
for effect sizes, assumption checks and post-hoc comparisons.

## Usage

``` r
# S3 method for class 'anovakit_fit'
print(x, digits = 4L, ...)
```

## Arguments

- x:

  An
  [`anovakit_fit`](https://elkronos.github.io/anovakit/reference/anovakit_fit.md)
  object.

- digits:

  Number of significant digits for the printed tables. Default `4`.
  Unlike most printing in R this is not capped by `getOption("digits")`.

- ...:

  Ignored.

## Value

`x`, invisibly.

## Examples

``` r
set.seed(1)
d <- data.frame(g = rep(c("a", "b", "c"), each = 20), y = rnorm(60))
fit <- anova_welch(d, "y", "g", plots = FALSE)
print(fit)
#> Welch's analysis of variance 
#> ----------------------------
#> Call: anova_welch(data = d, response = "y", groups = "g", plots = FALSE)
#> Observations used: 60
#> 
#> Omnibus test
#>   term statistic num_df den_df p_value
#> 1    g     0.264      2  37.91  0.7694
```
