# Plot an anovakit result

Returns one of the ggplot2 objects the analysis built. The available
names are listed by [`print()`](https://rdrr.io/r/base/print.html) and
stored in `x$plots`.

## Usage

``` r
# S3 method for class 'anovakit_fit'
plot(x, which = 1L, ...)
```

## Arguments

- x:

  An
  [`anovakit_fit`](https://elkronos.github.io/anovakit/reference/anovakit_fit.md)
  object.

- which:

  Name or index of the plot. Defaults to the first one.

- ...:

  Ignored.

## Value

A ggplot2 object, invisibly returning `NULL` when the fit holds no
plots.

## Examples

``` r
set.seed(1)
d <- data.frame(g = rep(c("a", "b", "c"), each = 20), y = rnorm(60))
fit <- anova_welch(d, "y", "g")
p <- plot(fit, which = "means")
```
