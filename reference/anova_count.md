# Analysis of Deviance for a Count Response

Fits a Poisson, quasi-Poisson or negative binomial regression for a
count response across one or more grouping variables, checks for
overdispersion, and reports an analysis of deviance table, incidence
rate ratios, estimated marginal rates and pairwise comparisons.

## Usage

``` r
anova_count(
  data,
  response,
  groups,
  offset = NULL,
  interaction = FALSE,
  model = c("auto", "poisson", "negbin", "quasipoisson"),
  overdispersion_threshold = 1.5,
  type = c("II", "III"),
  test_statistic = c("LR", "Wald", "F"),
  conf_level = 0.95,
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

  Character. Name of the count column. Must be non-negative whole
  numbers, not all zero.

- groups:

  Character vector. One or more grouping columns.

- offset:

  Optional character. Name of a strictly positive exposure column,
  entered as a log offset. It cannot also be the response, a grouping
  column or the weights.

- interaction:

  `FALSE` (additive, the default), `TRUE` (full factorial), or a whole
  number giving the highest interaction order.

- model:

  Character. `"auto"` (default), `"poisson"`, `"negbin"` or
  `"quasipoisson"`.

- overdispersion_threshold:

  Numeric. Pearson dispersion above which `model = "auto"` moves away
  from Poisson. Default `1.5`.

- type:

  Character. `"II"` (default) or `"III"` sums of squares. With `"III"`
  the model is fitted under sum-to-zero contrasts; this changes `$anova`
  only. A model with aliased coefficients (an empty cell of the design,
  or collinear predictors) has no Type III tests; Type II tests are then
  computed instead, and `$notes`, the method and the table heading say
  so.

- test_statistic:

  Character. `"LR"` (default), `"Wald"` or `"F"`, passed to
  [`Anova`](https://rdrr.io/pkg/car/man/Anova.html). `"F"` is the right
  choice for a quasi-Poisson model, which has no likelihood and is given
  it in place of `"LR"`. For a Poisson or negative binomial model, `"F"`
  makes [`car::Anova()`](https://rdrr.io/pkg/car/man/Anova.html)
  estimate a dispersion from the Pearson residuals, so the table becomes
  a quasi-likelihood F test while `$effect_sizes` and `$posthoc` still
  assume a dispersion of 1; `$notes` gives the dispersion the F test
  used. Use `model = "quasipoisson"` for quasi-likelihood inference
  throughout.

- conf_level:

  Numeric in (0, 1). Level for every interval returned. Default `0.95`.

- vcov_type:

  Character. `"model"` (default) or an HC type (`"HC0"` to `"HC4"`) for
  robust standard errors, which requires the sandwich package. They
  reach the coefficient table, the marginal means and the comparisons;
  the omnibus table uses them only with `test_statistic = "Wald"`, since
  a likelihood-ratio test compares deviances. `$notes` says which
  applies.

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

  Optional character. Name of a numeric column of prior weights, which
  multiply each row's log-likelihood as frequency weights would (see
  Details). With whole-number weights the robust covariance
  (`vcov_type`) is that of the data expanded to one row per observation.
  Given as a column name rather than a vector so that it is subsetted
  with the data: a vector supplied by the caller is evaluated against
  the original frame and silently misaligns as soon as one row is
  dropped for missing values.

- verbose:

  Logical. Emit progress through
  [`message`](https://rdrr.io/r/base/message.html). Default `FALSE`.

## Value

An
[`anovakit_fit`](https://elkronos.github.io/anovakit/reference/anovakit_fit.md)
object, with `$model_type` naming the model actually fitted,
`$dispersion` the Pearson dispersion of the *Poisson* fit (the statistic
the model choice was made on) and `$model_dispersion` that of the model
actually returned (both with the frequency-weight denominator when the
weights are frequencies). `$effect_sizes` holds incidence rate ratios
against the reference level, with columns `term` (the treatment-coded
coefficient), `factor` and `comparison` (what it compares), `IRR`,
`conf_low`, `conf_high` and `p_value`. `$emmeans` holds the estimated
marginal rates. `$assumptions` holds `poisson_dispersion` and
`model_dispersion` (as above) and `cell_counts`, the number of
observations in each cell (with frequency weights, `n` counts the
observations and `rows` the rows of the table). `$plots` holds
`emmeans`, the marginal rates with their intervals, and `observed`, box
plots of the observed counts by cell – of the observed rates, count
divided by exposure, when there is an `offset`.

## Details

**Overdispersion.** The Pearson dispersion statistic is computed from
the Poisson fit. With `model = "auto"` (the default) a dispersion above
`overdispersion_threshold` switches the model to negative binomial when
MASS is installed, and to quasi-Poisson otherwise (or when the negative
binomial fit fails). The threshold is a fixed cut-off, not a test, and a
Poisson model kept below it is still affected by a dispersion \\\phi\\
above 1: its likelihood-ratio and Wald statistics are inflated by about
\\\phi\\ and its intervals are too narrow, so at \\\phi = 1.5\\ a
nominal 5% test on 1 degree of freedom rejects about 11% of true null
hypotheses, however large the sample. `$notes` then gives the dispersion
and the size of the effect; use `model = "quasipoisson"` or `"negbin"`
when that matters. When the Poisson model has no residual degrees of
freedom the dispersion cannot be computed, and `$notes` says that
overdispersion was not checked. Which model was used, and why, is always
recorded in `$notes` and in `$model_type`. Set `model` explicitly to
take the decision out of the function's hands: an explicit
`model = "poisson"` is kept on overdispersed data (with a note), and an
explicit `model = "negbin"` that cannot be fitted is an error rather
than a substitution.

**Negative binomial inference.** The analysis of deviance treats the
dispersion parameter theta as known, and its tests and the Wald
intervals are asymptotic. With small groups they are anti-conservative,
and not only because theta is treated as known: with fewer than about 10
observations per group a nominal 5% test rejects roughly 10-15% of true
null hypotheses and 95% intervals cover about 90%; the excess fades by
about 30 per group. The note on theta states the size of the problem for
the smallest group in the data, whether theta is well determined, and
whether the fit converged. Convergence is read from the fitted object,
so the note is the same in every session language.

**Exposure.** An `offset` column enters the model as
`offset(log(exposure))` in the formula, not through
[`glm()`](https://rdrr.io/r/stats/glm.html)'s `offset` argument. That
matters: an offset supplied through the argument is invisible to
[`terms()`](https://rdrr.io/r/stats/terms.html), which means emmeans
cannot see it and `predict(newdata = )` silently recycles it. Estimated
marginal means are reported at an exposure of 1, so they are rates per
unit of exposure.

**Weights.** Prior weights multiply each row's contribution to the
likelihood, which is exactly what a frequency weight does: a row with
weight 3 counts as three identical observations, and a table with one
row per distinct count and a column of frequencies gives the same
estimates and likelihood-ratio tests as the expanded data. When every
weight is a whole number and some exceed 1, the weights are taken to be
frequencies and the Pearson dispersion (the one the model choice is made
on, reported in `$dispersion`) divides by the number of observations
they represent minus the number of parameters rather than by the number
of rows, which would inflate it by up to the ratio of the two. A
quasi-Poisson fit and an F test estimate their dispersion from the rows
all the same, so for those expand a frequency table to one row per
observation; `$notes` says so.

**Incidence rate ratios.** `$effect_sizes` is computed from a
treatment-coded parameterisation of the fitted model, whatever `type`
and the global `contrasts` option are, and with ordered factors treated
as unordered. Each ratio compares one level of a factor with its first
(reference) level, as the `factor` and `comparison` columns say (the
same layout as the odds ratios of
[`anova_bin`](https://elkronos.github.io/anovakit/reference/anova_bin.md));
in a model with interactions a main-effect ratio is taken at the
reference levels of the factors it interacts with, and an interaction
term, labelled `"(B vs A) x (Y vs X)"`, is a ratio of rate ratios. The
coefficients and their covariance (model-based or robust) are mapped
exactly onto that coding, which gives what a refit with
`contr.treatment` would give.

**Interactions.** The model is additive by default. Set
`interaction = TRUE` for the full factorial, or to a whole number for
the highest interaction order. With three or more grouping variables the
full factorial is often unestimable, so this is deliberately not the
default.

## See also

[`anova_bin`](https://elkronos.github.io/anovakit/reference/anova_bin.md)
for binary responses.

## Examples

``` r
set.seed(1)
n <- 400
d <- data.frame(
  g1 = factor(rep(c("A", "B"), each = n / 2)),
  g2 = factor(rep(rep(c("X", "Y"), each = n / 4), 2))
)
lambda <- with(d, ifelse(g1 == "A" & g2 == "X", 5,
                  ifelse(g1 == "A" & g2 == "Y", 10,
                  ifelse(g1 == "B" & g2 == "X", 15, 20))))
d$count <- rpois(n, lambda)
fit <- anova_count(d, "count", c("g1", "g2"), interaction = TRUE)
fit
#> Analysis of deviance for a count response (poisson regression) 
#> --------------------------------------------------------------
#> Call: anova_count(data = d, response = "count", groups = c("g1", "g2"), interaction = TRUE)
#> Observations used: 400
#> 
#> Omnibus test (Type II, likelihood-ratio chi-square)
#>    term df statistic   p_value
#> 1    g1  1    772.80 < 2.2e-16
#> 2    g2  1    232.80 < 2.2e-16
#> 3 g1:g2  1     23.96 9.855e-07
#> 
#> Notes
#>   - Pearson dispersion is 0.98, at or below the threshold of 1.50, so a Poisson model was kept.
#> 
#> Plots available: emmeans, observed
#>   (use plot(x, which = "emmeans"))
fit$emmeans
#>   g1 g2 estimate        se  df  conf_low conf_high
#> 1  A  X     5.10 0.2258314 Inf  4.676042  5.562397
#> 2  B  X    14.17 0.3764305 Inf 13.451088 14.927335
#> 3  A  Y     9.88 0.3143246 Inf  9.282749 10.515678
#> 4  B  Y    20.05 0.4477723 Inf 19.191313 20.947108
fit$posthoc
#>    contrast       IRR         se  df  conf_low conf_high null statistic
#> 1 A X / B X 0.3599153 0.01858534 Inf 0.3152005 0.4109735    1 -19.78939
#> 2 A X / A Y 0.5161943 0.02814524 Inf 0.4487240 0.5938095    1 -12.12798
#> 3 A X / B Y 0.2543641 0.01261484 Inf 0.2239357 0.2889271    1 -27.60411
#> 4 B X / A Y 1.4342105 0.05944385 Inf 1.2893467 1.5953505    1   8.70060
#> 5 B X / B Y 0.7067332 0.02452750 Inf 0.6464486 0.7726395    1 -10.00137
#> 6 A Y / B Y 0.4927681 0.01915403 Inf 0.4459379 0.5445162    1 -18.20714
#>         p_value   p_adjusted adjustment
#> 1  3.674528e-87 0.000000e+00      tukey
#> 2  7.508225e-34 0.000000e+00      tukey
#> 3 9.932148e-168 0.000000e+00      tukey
#> 4  3.301345e-18 3.297362e-14      tukey
#> 5  1.503059e-23 3.164136e-14      tukey
#> 6  4.529832e-74 0.000000e+00      tukey

# With an exposure offset, the marginal means are rates per unit exposure
d$hours <- runif(n, 0.5, 4)
anova_count(d, "count", "g1", offset = "hours", plots = FALSE)$emmeans
#>   g1  estimate        se  df conf_low conf_high
#> 1  A  4.309548 0.2183641 Inf 3.902128  4.759507
#> 2  B 10.077180 0.4719456 Inf 9.193367 11.045959
```
