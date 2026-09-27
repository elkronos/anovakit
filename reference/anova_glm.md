# Analysis of Deviance for a Generalised Linear Model

Fits a generalised linear model of a response on one or more grouping
variables in any [`family`](https://rdrr.io/r/stats/family.html), and
reports an analysis of deviance table, effect sizes, estimated marginal
means and pairwise comparisons.

## Usage

``` r
anova_glm(
  data,
  response,
  groups,
  family = stats::gaussian(),
  interaction = FALSE,
  type = c("II", "III"),
  test_statistic = NULL,
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

  Character. Name of the response column. For a binomial or
  quasi-binomial family it may be 0/1, logical, a proportion (with
  `weights` giving the numbers of trials), or a factor or character
  column with exactly two observed values; the model is then for the
  second level (in C-locale order for a character column), and `$notes`
  names it.

- groups:

  Character vector. One or more grouping columns.

- family:

  A [`family`](https://rdrr.io/r/stats/family.html) object, a family
  function, or the name of one as a character string. Default
  [`stats::gaussian()`](https://rdrr.io/r/stats/family.html).

- interaction:

  `FALSE` (additive, the default), `TRUE` (full factorial), or a whole
  number giving the highest interaction order.

- type:

  Character. `"II"` (default) or `"III"` sums of squares. A model with
  aliased coefficients (an empty cell of the design, or collinear
  predictors) has no Type III tests; Type II tests are then computed
  instead, and `$notes`, the method and the table heading say so.

- test_statistic:

  Character or `NULL`. One of `"LR"`, `"Wald"` or `"F"`, passed to
  [`Anova`](https://rdrr.io/pkg/car/man/Anova.html). When `NULL` (the
  default) it is chosen after fitting: an F test when the fitted model
  estimates its dispersion – the Gaussian, Gamma, inverse Gaussian,
  quasi, quasi-Poisson and quasi-binomial families, and
  `MASS::negative.binomial(theta)`, whose coefficient table also uses t
  – and a likelihood ratio test when the dispersion is fixed (binomial,
  Poisson).

- conf_level:

  Numeric in (0, 1). Level for every interval returned. Default `0.95`.

- ci_method:

  Character. `"profile"` (default) or `"wald"` for the coefficient
  intervals. Both use the reference distribution of the coefficient
  table's p-values. For a family with a fixed dispersion (binomial,
  Poisson) the profile interval is
  [`stats::confint()`](https://rdrr.io/r/stats/confint.html)'s, which
  cuts the signed-root deviance profile at normal quantiles. For a
  family that estimates its dispersion the same profile is cut at
  `qt(., df.residual)` instead, so that the interval agrees with the t
  test beside it; for the Gaussian identity model it equals
  `confint(lm())`. If profiling fails, Wald intervals are used, with a
  note. A bound the profile cannot reach (as under separation) is `NA`,
  also with a note.

- vcov_type:

  Character. `"model"` (default) or an HC type (`"HC0"` to `"HC4"`),
  which requires the sandwich package. The robust covariance is used by
  the coefficient table (with Wald intervals), the marginal means and
  the comparisons. The omnibus likelihood ratio and F tests are computed
  from the likelihood and cannot use it; with `test_statistic = "Wald"`
  the omnibus test uses it too. `$notes` says which applies.

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

- weights:

  Optional character. Name of a numeric column of prior weights. Given
  as a column name rather than a vector so that it is subsetted with the
  data: a vector supplied by the caller is evaluated against the
  original frame and silently misaligns as soon as one row is dropped
  for missing values. Rows with weight zero contribute nothing to the
  fit; they are dropped before fitting and counted in `$n_removed`.
  Whole-number weights (some above 1) on a count or 0/1 response are
  read as frequency weights, each row standing for that many identical
  observations: the robust covariance (`vcov_type`) and, for the Poisson
  families, the residual degrees of freedom are then those of the data
  expanded to one row per observation. Other weights (precision weights,
  or the trials behind a proportion) keep each row as one unit.

- verbose:

  Logical. Emit progress through
  [`message`](https://rdrr.io/r/base/message.html). Default `FALSE`.

## Value

An
[`anovakit_fit`](https://elkronos.github.io/anovakit/reference/anovakit_fit.md)
object, with `$family` (the family object actually used) and
`$model_stats` (AIC, BIC, the null and residual deviances, and the
residual degrees of freedom). `$anova` has the columns `term`, `sum_sq`
(F tests only), `df`, `statistic` and `p_value`, in that order whatever
`test_statistic`. In an F table of a non-Gaussian family, as car
computes it, a term's `sum_sq` is its change in deviance and the
`Residuals` row holds the Pearson chi-square, from which the dispersion
the F statistic divides by is estimated. `$assumptions` holds
`coefficients` – the coefficient table, with the intervals `ci_method`
selects – and `dispersion`. `$effect_sizes` depends on the family:

- Gaussian: partial eta squared and partial omega squared for each term,
  from the F table's sums of squares whatever `test_statistic` is.
  Partial omega squared is floored at 0, with a note. When the ANOVA
  table cannot be computed (a Type III table with an empty cell, say)
  the table is empty and a note says why.

- Any other family: `deviance_explained`,
  `1 - deviance / null deviance`, and `mcfadden_r2`, McFadden's pseudo R
  squared `1 - logLik(model) / logLik(intercept-only model)`. McFadden's
  measure is reported only for the binomial, Poisson and negative
  binomial families, whose log-likelihoods are probabilities; for a
  binomial response it is computed from the Bernoulli log-likelihood, so
  it is the same whether the data are given one row per trial or
  aggregated as proportions with trial weights (and for a 0/1 response
  it equals the deviance explained). For the Gamma, inverse Gaussian and
  quasi families it is `NA` with a note: a density's log-likelihood
  changes with the units of the response, and a quasi family has none.
  The deviance explained is measured against the saturated model, so for
  binomial data it does depend on the aggregation.

## Details

**Type III sums of squares.** When `type = "III"` the model is *fitted*
under sum-to-zero contrasts. This is the whole point: setting
`options(contrasts = )` after a model has been fitted does not change
the contrasts stored on it, so `car::Anova(type = 3)` would silently
report simple effects at the reference level under the label "Type III".

**Post-hoc.** Marginal means and comparisons come from emmeans and are
on the response scale. When the grouping variables enter the model
additively (the default with more than one), they are reported for each
factor on its own, averaged over the others; with an interaction they
are made across the cells of the grid. With a log or logit link a
comparison is a ratio (of means, or of odds); with any other
non-identity link (inverse, probit, square root, ...) it is a difference
between the marginal means on the response scale. Families that estimate
their dispersion use t with the residual degrees of freedom, as the
coefficient table does.

**Separation.** For binomial, quasi-binomial, Poisson, quasi-Poisson and
negative binomial families the fit is checked for groups whose fitted
probability is numerically 0 or 1 (or whose fitted rate is numerically
0), where coefficients, Wald tests and comparisons are not identified;
`$notes` says so when it happens.

For binary responses see
[`anova_bin`](https://elkronos.github.io/anovakit/reference/anova_bin.md),
and for counts see
[`anova_count`](https://elkronos.github.io/anovakit/reference/anova_count.md):
both add family-specific diagnostics and effect sizes that this general
wrapper does not.

## See also

[`anova_bin`](https://elkronos.github.io/anovakit/reference/anova_bin.md)
and
[`anova_count`](https://elkronos.github.io/anovakit/reference/anova_count.md),
which add family-specific diagnostics and effect sizes.

## Examples

``` r
set.seed(42)
d <- data.frame(
  A = factor(sample(c("a1", "a2"), 200, replace = TRUE)),
  B = factor(sample(c("b1", "b2"), 200, replace = TRUE))
)
d$y <- 3 + 2 * (d$A == "a2") + (d$B == "b2") +
  4 * (d$A == "a2") * (d$B == "b2") + rnorm(200)

# Type III tests, with the model fitted under sum-to-zero contrasts
fit <- anova_glm(d, "y", c("A", "B"), interaction = TRUE, type = "III")
fit$anova
#>        term   sum_sq  df statistic      p_value
#> 1         A 828.3928   1  924.0123 4.109856e-76
#> 2         B 378.7873   1  422.5099 8.435778e-51
#> 3       A:B 192.5110   1  214.7321 2.549481e-33
#> 4 Residuals 175.7173 196        NA           NA

# The coefficient table, with profile-likelihood intervals by default
fit$assumptions$coefficients
#>          term   estimate         se statistic       p_value distribution
#> 1 (Intercept)  5.4973240 0.06790639  80.95444 1.359261e-152            t
#> 2          A1 -2.0641894 0.06790639 -30.39757  4.109856e-76            t
#> 3          B1 -1.3958190 0.06790639 -20.55505  8.435778e-51            t
#> 4       A1:B1  0.9950827 0.06790639  14.65374  2.549481e-33            t
#>     conf_low conf_high ci_method
#> 1  5.3634030  5.631245   profile
#> 2 -2.1981104 -1.930268   profile
#> 3 -1.5297399 -1.261898   profile
#> 4  0.8611617  1.129004   profile

# A family may be given as a string, as stats::glm() allows. For a binary
# response anova_bin() is usually the better entry point.
d$hit <- rbinom(200, 1, plogis(-0.5 + 1.2 * (d$A == "a2")))
bin <- anova_glm(d, "hit", "A", family = "binomial", plots = FALSE)
bin$anova
#>   term df statistic      p_value
#> 1    A  1  22.84952 1.751931e-06
bin$effect_sizes
#>              measure   estimate
#> 1 deviance_explained 0.08270479
#> 2        mcfadden_r2 0.08270479

# Heteroskedasticity-consistent standard errors, which flow through to the
# marginal means and the comparisons
anova_glm(d, "y", "A", vcov_type = "HC3", plots = FALSE)$emmeans
#>    A estimate        se  df conf_low conf_high
#> 1 a1 3.365595 0.1004602 198 3.167486  3.563704
#> 2 a2 7.583053 0.2482866 198 7.093427  8.072679
```
