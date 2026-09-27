# Analysis of Deviance for a Binary Response

Fits a logistic regression for a binary response across one or more
grouping variables and reports an analysis of deviance table, odds
ratios with confidence intervals, estimated marginal probabilities, and
pairwise comparisons on the odds-ratio scale.

## Usage

``` r
anova_bin(
  data,
  response,
  groups,
  success = NULL,
  reference = NULL,
  interaction = FALSE,
  type = c("II", "III"),
  test_statistic = c("LR", "Wald"),
  conf_level = 0.95,
  ci_method = c("profile", "wald"),
  vcov_type = "model",
  adjust = "tukey",
  posthoc = TRUE,
  plots = TRUE,
  weights = NULL,
  verbose = FALSE
)
```

## Arguments

- data:

  A data frame, or anything inheriting from one, such as a `data.table`
  or a tibble.

- response:

  Character. Name of the binary response column. May be a factor with
  two observed levels, a logical, a numeric 0/1 vector, or a character
  vector with two distinct values (whose levels are sorted in the C
  locale, so the default success level does not depend on the session).

- groups:

  Character vector. One or more grouping columns.

- success:

  Optional. A single value: the response level to model as the
  "success". By default the second level of the factor, which is what
  [`glm()`](https://rdrr.io/r/stats/glm.html) uses. Works for ordered
  factors too.

- reference:

  Optional named list of reference levels for the grouping factors, one
  level per factor, for example `list(site = "north")`. An ordered
  factor has no reference level (its contrasts are polynomial trends),
  so naming one here makes the function treat that factor as unordered,
  with a note; its levels are never reordered as a scale.

- interaction:

  `FALSE` (additive, the default), `TRUE` (full factorial), or a whole
  number giving the highest interaction order.

- type:

  Character. `"II"` (default) or `"III"` sums of squares. A model with
  aliased coefficients (an empty cell of the design, or collinear
  predictors) has no Type III tests; Type II tests are then computed
  instead, and `$notes`, the method and the table heading say so.

- test_statistic:

  Character. `"LR"` (default) or `"Wald"`, passed to
  [`Anova`](https://rdrr.io/pkg/car/man/Anova.html).

- conf_level:

  Numeric in (0, 1). Level for every interval returned. Default `0.95`.

- ci_method:

  Character. `"profile"` (default) for profile-likelihood intervals, or
  `"wald"`.

- vcov_type:

  Character. `"model"` (default) or an HC type (`"HC0"` to `"HC4"`) for
  robust standard errors, which requires the sandwich package. They
  reach the coefficient table, the marginal means and the comparisons;
  the omnibus table uses them only with `test_statistic = "Wald"`, since
  a likelihood-ratio test compares deviances. `$notes` says which
  applies.

- adjust:

  Character. Multiplicity adjustment for the pairwise comparisons,
  passed to
  [`contrast`](https://rvlenth.github.io/emmeans/reference/contrast.html).
  Default `"tukey"`.

- posthoc:

  Logical. Compute pairwise comparisons. There are `choose(k, 2)` of
  them, so this is worth turning off when the number of cells is large;
  `$notes` records that they were skipped. Default `TRUE`.

- plots:

  Logical. Build ggplot2 objects. They are returned in `$plots`, never
  drawn. Default `TRUE`.

- weights:

  Optional character. Name of a numeric column of prior weights. Given
  as a column name rather than a vector so that it is subsetted with the
  data: a vector supplied by the caller is evaluated against the
  original frame and silently misaligns as soon as one row is dropped
  for missing values. Whole-number weights (some above 1) are read as
  frequency weights, each row standing for that many identical
  observations, and the robust covariance (`vcov_type`) is then that of
  the data expanded to one row per observation; other weights are
  treated as sampling weights, with each row as one unit.

- verbose:

  Logical. Emit progress through
  [`message`](https://rdrr.io/r/base/message.html). Default `FALSE`.

## Value

An
[`anovakit_fit`](https://elkronos.github.io/anovakit/reference/anovakit_fit.md)
object. Besides the standard components:

- `$model_stats`:

  A list: `AIC`, `BIC` and `logLik`; `mcfadden_r2`, McFadden's pseudo R
  squared, `1 - deviance / null deviance`; `success_level`, the level
  modelled as the success; and `lr_vs_null`, the model-versus-null
  likelihood-ratio test (`statistic` = null deviance minus deviance,
  `df`, `p_value`). All of them use the prior weights. When any weight
  is not a whole number, `AIC`, `BIC` and `logLik` are `NA` with a note,
  because R's binomial log-likelihood rounds the weights; the
  deviance-based statistics are unaffected. `BIC` uses
  [`nobs()`](https://rdrr.io/r/stats/nobs.html), the number of rows with
  a non-zero weight.

- `$effect_sizes`:

  One row per coefficient of the reference-coded model (see Details):
  `term` (the coefficient), `factor` (the grouping variable, or
  `"g1:g2"` for an interaction), `comparison` (for example `"b vs a"`, a
  level against the reference level; in a model with interactions a
  lower-order odds ratio holds at the reference level of the factors it
  interacts with, `"b vs a at site = north"`, and an interaction row
  such as `"(b vs a) x (south vs north)"` is a ratio of odds ratios),
  `odds_ratio`, `conf_low`, `conf_high`, `se_log_or`, `statistic`,
  `p_value` and `ci_method`.

- `$assumptions`:

  `dispersion`, which is `NA`: the Pearson dispersion of a 0/1 response
  carries no information about overdispersion, so it is not reported (a
  note says why); and `proportions`, the observed proportion of each
  response level in each cell, with columns for the cell, the response
  level, `n` (the weighted count), `group_total` (the weighted cell
  total) and `proportion`; when `weights` is given, `rows` and
  `group_rows` give the number of data rows behind them. A column whose
  name the response already uses gets a suffix, e.g. `n.1`.

`$emmeans` holds the estimated marginal probabilities and `$posthoc` the
pairwise odds ratios.

## Details

**Interactions.** By default the model is additive. Set
`interaction = TRUE` for the full factorial, or to a whole number for
the highest order of interaction to include. An additive model cannot
detect an interaction between grouping variables, so if you have more
than one grouping variable this is a decision worth making deliberately.

**Type III sums of squares.** When `type = "III"` the model is *fitted*
under sum-to-zero contrasts, which is what makes Type III tests
meaningful. Setting the global contrast option after fitting has no
effect on an existing model.

**Odds ratios.** Under sum-to-zero or polynomial contrasts an
exponentiated coefficient is not an odds ratio against anything a reader
would recognise. `$effect_sizes` is therefore computed from a
reference-coded (treatment-coded) fit of the same model, with every
grouping factor treated as unordered and `reference` honoured, whatever
`type`, the global `contrasts` option or the class of the factor.
`$anova` still comes from the Type II or Type III fit, which describes
the same fitted probabilities.

**Separation.** Complete and quasi-complete separation are detected and
reported in `$notes`. Under separation, Wald odds ratios and their
intervals are meaningless even though
[`glm()`](https://rdrr.io/r/stats/glm.html) reports convergence, so
`ci_method = "profile"` is the default. The profile-likelihood interval
of an odds ratio that separation drives to infinity (or to zero) is open
on that side, and is reported with `conf_high = Inf` (or
`conf_low = 0`); its finite end is found by profiling the likelihood
directly, since
[`stats::confint()`](https://rdrr.io/r/stats/confint.html) cannot step
from a diverged estimate, and is `NA` when the likelihood rules out no
value on that side either.

**Sparse data.** When fewer than about 5 events or non-events are
expected in a cell under the null hypothesis, the likelihood-ratio
omnibus test is liberal (its false-positive rate can be well above the
nominal level) and `$notes` says so. `test_statistic = "Wald"` is not a
remedy; an exact test, or the score test, holds its level better.

**Weights.** The prior weights enter the fit, the model statistics and
the observed proportions alike. A row with weight `w` counts as `w`
identical observations, so data aggregated to one row per cell and
outcome, with the count as the weight, give the same results as the
individual rows (apart from BIC, which R bases on the number of rows).

## See also

[`anova_count`](https://elkronos.github.io/anovakit/reference/anova_count.md)
for counts,
[`anova_glm`](https://elkronos.github.io/anovakit/reference/anova_glm.md)
for other families.

## Examples

``` r
set.seed(1)
n <- 300
d <- data.frame(g = factor(sample(c("a", "b", "c"), n, replace = TRUE)))
d$y <- rbinom(n, 1, c(a = 0.2, b = 0.5, c = 0.8)[as.character(d$g)])
fit <- anova_bin(d, "y", "g")
fit
#> Analysis of deviance for a binary response (logistic regression) 
#> ----------------------------------------------------------------
#> Call: anova_bin(data = d, response = "y", groups = "g")
#> Observations used: 300
#> 
#> Omnibus test (Type II, likelihood-ratio chi-square)
#>   term df statistic   p_value
#> 1    g  2     70.57 4.738e-16
#> 
#> Notes
#>   - Modelling P(y = 1); the other level is the baseline.
#>   - $assumptions$dispersion is NA: the Pearson dispersion of a 0/1 response carries no information about overdispersion (in a model saturated in the cells it equals N / (N - k) whatever the data).
#> 
#> Plots available: proportions, emmeans
#>   (use plot(x, which = "proportions"))
fit$effect_sizes
#>   term factor comparison odds_ratio conf_low conf_high se_log_or statistic
#> 1   gb      g     b vs a   2.877857 1.620408   5.19644 0.2966147  3.563701
#> 2   gc      g     c vs a  14.483333 7.360953  30.08819 0.3578074  7.470496
#>        p_value ci_method
#> 1 3.656629e-04   profile
#> 2 7.989300e-14   profile

# Two factors, with their interaction
d$site <- factor(sample(c("north", "south"), n, replace = TRUE))
anova_bin(d, "y", c("g", "site"), interaction = TRUE, plots = FALSE)$anova
#>     term df statistic      p_value
#> 1      g  2 67.958757 1.749619e-15
#> 2   site  1  3.417930 6.449183e-02
#> 3 g:site  2  5.548029 6.241095e-02

# Aggregated data: one row per group and outcome, the count as the weight
agg <- data.frame(g = rep(c("a", "b", "c"), each = 2), y = rep(0:1, 3),
                  n = c(80, 20, 50, 50, 20, 80))
anova_bin(agg, "y", "g", weights = "n", plots = FALSE)$model_stats
#> $AIC
#> [1] 344.7904
#> 
#> $BIC
#> [1] 344.1657
#> 
#> $logLik
#> [1] -169.3952
#> 
#> $mcfadden_r2
#> [1] 0.1853813
#> 
#> $success_level
#> [1] "1"
#> 
#> $lr_vs_null
#>   statistic df      p_value
#> 1   77.0979  2 1.813022e-17
#> 
```
