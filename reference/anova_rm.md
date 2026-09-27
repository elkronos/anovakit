# Repeated Measures Analysis of Variance

Fits a repeated measures ANOVA to data in long format using
[`aov_ez`](https://rdrr.io/pkg/afex/man/aov_car.html), reports Mauchly's
test of sphericity with the Greenhouse-Geisser and Huynh-Feldt
corrections, generalised and partial eta squared, residual diagnostics,
and estimated marginal means.

## Usage

``` r
anova_rm(
  data,
  response,
  subject,
  within,
  between = NULL,
  factorize = TRUE,
  emm_specs = NULL,
  adjust = "tukey",
  conf_level = 0.95,
  correction = c("GG", "HF", "none"),
  posthoc = TRUE,
  plots = TRUE,
  verbose = FALSE,
  ...
)
```

## Arguments

- data:

  A data frame in long format: one row per subject per cell.

- response:

  Character. Name of the numeric response column.

- subject:

  Character. Name of the subject identifier column.

- within:

  Character vector. One or more within-subject factors.

- between:

  Character vector or `NULL`. Between-subject factors, each constant
  within a subject.

- factorize:

  Logical. Convert the response to numeric through its labels. Default
  `TRUE`; with `FALSE` a non-numeric response is an error. The subject,
  within- and between-subject columns are always coerced to factors,
  whatever this is set to.

- emm_specs:

  Character vector or `NULL`. Factors to compute estimated marginal
  means over. Defaults to all within-subject factors; when it leaves out
  a factor of the model, the means are averaged over that factor's
  levels and `$notes` says so.

- adjust:

  Character. Multiplicity adjustment for the pairwise comparisons,
  applied by emmeans. One of `"tukey"`, `"sidak"`, `"scheffe"`,
  `"dunnettx"`, `"bonferroni"`, `"holm"`, `"hochberg"`, `"hommel"`,
  `"BH"`, `"BY"`, `"fdr"` or `"none"`. Default `"tukey"`.

- conf_level:

  Numeric in (0, 1). Level for every interval returned. Default `0.95`.

- correction:

  Character. Sphericity correction applied to the ANOVA table: `"GG"`
  (default), `"HF"` or `"none"`. The correction actually applied is
  recorded in `attr(fit$anova, "correction")` and printed: `"none"` when
  every within-subject factor has two levels, where there is nothing to
  correct.

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

- ...:

  Only these arguments of
  [`aov_ez`](https://rdrr.io/pkg/afex/man/aov_car.html) are accepted;
  anything else is an error rather than silently ignored.

  `fun_aggregate`

  :   A function combining the rows of a subject that share a
      within-subject cell into one number. Default `mean`.

  `observed`

  :   Character. Factors that were observed (measured) rather than
      manipulated. Generalised eta squared depends on this (Olejnik and
      Algina, 2003).

  `type`

  :   Sums of squares: `3` (or `"III"`, the default) or `2` (`"II"`).

  `anova_table`

  :   A list holding only `p_adjust_method`, one of
      [`p.adjust.methods`](https://rdrr.io/r/stats/p.adjust.html), which
      adjusts the p-values of `$anova` across its terms. The `p_gg` and
      `p_hf` columns of `$sphericity` stay unadjusted.

  `covariate` is refused: `anova_rm()` does not fit covariates. `id`,
  `dv`, `return`, `include_aov` and `print.formula` are set by the
  function itself, and `transformation` is refused (transform the
  response column instead).

## Value

An
[`anovakit_fit`](https://elkronos.github.io/anovakit/reference/anovakit_fit.md)
object. `$anova` has the corrected degrees of freedom, the MSE, F,
partial eta squared and the p-value for each term; `$effect_sizes` has
partial and generalised eta squared, the latter honouring `observed`.
The function adds

- `$sphericity`:

  Mauchly's test (`mauchly_w`, `p_value`, `NA` where it is undefined),
  the Greenhouse-Geisser and Huynh-Feldt epsilons (`gg_epsilon`,
  `hf_epsilon`, the latter capped at 1), the p-values under each
  correction (`p_gg`, `p_hf`) and the uncapped Huynh-Feldt estimate
  (`hf_epsilon_raw`), one row per within-subject term with more than one
  degree of freedom.

- `$residuals`:

  The within-subject residuals (see Details), one per row of
  `$data_used` and in the same order, so `$residuals[i]` belongs to
  `$data_used[i, ]`.

- `$subjects_dropped`:

  The subjects removed for incomplete within-subject cells.

- `$internal_names`:

  A data frame with the `role`, your `name` and the `internal` name each
  column has in `$model`.

- `$n_removed_missing`, `$n_removed_unbalanced`,
  `$n_removed_aggregated`:

  `$n_removed` broken into its three causes.

`$data_used` has one row per subject and within-subject cell, in the
order the rows came in (a combined row takes the place of its first
row). `$assumptions$normality` has one Shapiro-Wilk test per
within-subject cell, with a `note` giving the reason when a cell could
not be tested.

## Details

**What is fixed rather than inherited.** The sums of squares are Type
III unless `type = 2` is given through `...`; the between-subject
factors are coded with sum-to-zero contrasts, which Type III needs; and
the marginal means and comparisons use afex's multivariate model, whose
standard errors and degrees of freedom come from each within-subject
cell's own variance rather than a pooled error term. None of these is
read from
[`afex::afex_options()`](https://rdrr.io/pkg/afex/man/afex_options.html),
so a global setting made for some other analysis cannot change the
result. The type is recorded in `attr(fit$anova, "ss_type")` (and
printed), the emmeans model in `attr(fit$emmeans, "emmeans_model")` and
`attr(fit$posthoc, "emmeans_model")`.

**Sphericity.** Mauchly's test and both epsilon estimates are read from
the multivariate model afex fitted (through
[`car::summary.Anova.mlm()`](https://rdrr.io/pkg/car/man/Anova.html)).
Each term's error matrix is first rescaled to unit average variance, so
a response measured in small units does not fall below car's *absolute*
singularity tolerance and lose its corrections. When the error matrix is
genuinely singular – fewer subjects than the term has contrasts – car
cannot give the corrections and Mauchly's test is undefined. The
Greenhouse-Geisser epsilon, \\tr(S)^2 / ((k-1) tr(S^2))\\, is still
defined, where \\S\\ is the covariance matrix of the orthonormal
within-subject contrasts; it and the Huynh-Feldt epsilon (with
Lecoutre's correction) are then computed directly, the requested
correction is applied with them, and `$notes` says so. The Huynh-Feldt
estimate can exceed 1; it is used, and reported in
`$sphericity$hf_epsilon`, capped at 1, with the uncapped value in
`$sphericity$hf_epsilon_raw`. The corrected degrees of freedom are the
uncorrected ones multiplied by epsilon, as in afex.

**Unbalanced subjects.** A repeated measures ANOVA needs every subject
to appear in every within-subject cell. Subjects who do not are removed
before fitting, and both the count and the identifiers (in the order of
the subject column's levels) are reported in `$subjects_dropped` and in
`$notes`.

**Repeated rows.** A subject with more than one row in the same
within-subject cell has those rows combined into one value before
fitting, with `fun_aggregate` (the mean unless another function is
given). The note says how many subject-by-cell combinations were
affected and how many extra rows were combined away; those rows are
counted in `$n_removed`.

**Column names and labels.** afex builds a model formula from the column
names and wide-format column names from the within-subject levels, so a
column name that is not a syntactic R name (`"subject id"`), or level
labels such as `0, 4, 12`, would be mangled or parsed as code. afex is
therefore given internal names (`Y`, `ID`, `W1`, `W2`, ..., `B1`, ...)
and level labels (`L01`, `L02`, ... in level order), and every table,
`$data_used`, `$emmeans_object` and the plots are mapped back to your
names and labels. `$model` is the afex fit itself and keeps the internal
names; `$internal_names` gives the correspondence.

**Normality.** The F tests for within-subject effects depend on the
within-subject errors: what is left of each value after its subject's
mean and its cell's mean (within each between-subject group) are
removed, \\y\_{ij} - \bar{y}\_{i\cdot} - \bar{y}\_{g(i)\cdot j} +
\bar{y}\_{g(i)\cdot\cdot}\\. These are `$residuals`, and Shapiro-Wilk is
run on them separately in each within-subject cell
(`$assumptions$normality`), as
[`anova_welch`](https://elkronos.github.io/anovakit/reference/anova_welch.md)
tests each group: the cells' variances usually differ, and pooling
values of different spread gives a mixture that fails a normality test
even when every cell is normal. With a single two-level within-subject
factor the residuals are plus and minus half of each subject's centred
difference, so both rows test the normality of the differences, which is
exactly the assumption of the equivalent paired t-test. The Q-Q plot
shows the same residuals standardised within each cell. Normality of the
subject means, which the between-subject tests rely on, is not tested.

**A factor response.** With `factorize = TRUE`, a response stored as a
factor is converted with `as.numeric(as.character(x))`, which recovers
the numbers. A bare
[`as.numeric()`](https://rdrr.io/r/base/numeric.html) on a factor
returns the level *indices*, so a response of 11, 17, 21 would silently
become 1, 2, 3.

Requires the afex package.

## References

Olejnik, S., & Algina, J. (2003). Generalized eta and omega squared
statistics: measures of effect size for some common research designs.
*Psychological Methods*, 8(4), 434-447.

Lecoutre, B. (1991). A correction for the epsilon-tilde approximate test
in repeated measures designs with two or more independent groups.
*Journal of Educational Statistics*, 16(4), 371-372.

## See also

[`anova_welch`](https://elkronos.github.io/anovakit/reference/anova_welch.md)
for independent groups,
[`aov_ez`](https://rdrr.io/pkg/afex/man/aov_car.html) for the fit
itself.

## Examples

``` r
if (requireNamespace("afex", quietly = TRUE)) {
  set.seed(1)
  d <- expand.grid(id = factor(1:24), time = factor(c("t1", "t2", "t3")))
  d$arm <- factor(rep(rep(c("ctrl", "trt"), each = 12), 3))
  d$score <- 10 + 2 * as.numeric(d$time) +
    1.5 * (d$arm == "trt") + rnorm(nrow(d), 0, 2)

  fit <- anova_rm(d, "score", subject = "id",
                  within = "time", between = "arm")
  print(fit)

  # Mauchly's test with both epsilon corrections, read off the fitted model
  print(fit$sphericity)

  # Generalised eta squared alongside partial eta squared, and the marginal
  # means over the within-subject factor (averaged over arm, as $notes says)
  print(fit$effect_sizes)
  print(fit$emmeans)

  # A factor that was measured rather than assigned (an age group, say) is
  # declared with `observed`, which changes generalised eta squared
  d$age_grp <- factor(rep(rep(c("younger", "older"), 12), 3))
  anova_rm(d, "score", subject = "id", within = "time",
           between = c("arm", "age_grp"), observed = "age_grp",
           plots = FALSE)$effect_sizes
}
#> Repeated measures analysis of variance 
#> --------------------------------------
#> Call: anova_rm(data = d, response = "score", subject = "id", within = "time", between = "arm")
#> Observations used: 72
#> 
#> Omnibus test (Type III, F, GG-corrected degrees of freedom)
#>       term num_df den_df   mse statistic partial_eta_sq   p_value
#> 1      arm  1.000  22.00 3.274   14.8000         0.4022 0.0008754
#> 2     time  1.994  43.87 3.672   29.7600         0.5750 7.013e-09
#> 3 arm:time  1.994  43.87 3.672    0.4833         0.0215 0.6194033
#> 
#> Notes
#>   - The Huynh-Feldt epsilon estimate exceeds 1 for time, arm:time (1.096, 1.096). Epsilon cannot exceed 1, so it is capped at 1 in $sphericity$hf_epsilon and wherever it is used (p_hf); the uncapped estimate is in $sphericity$hf_epsilon_raw.
#>   - The marginal means in $emmeans and the comparisons in $posthoc are averaged over the levels of: arm.
#> 
#> Plots available: residuals, qq, index, emmeans
#>   (use plot(x, which = "residuals"))
#>       term mauchly_w   p_value gg_epsilon hf_epsilon         p_gg         p_hf
#> 1     time 0.9969409 0.9683427  0.9969503          1 7.013415e-09 6.686437e-09
#> 2 arm:time 0.9969409 0.9683427  0.9969503          1 6.194033e-01 6.199669e-01
#>   hf_epsilon_raw
#> 1       1.096159
#> 2       1.096159
#>       term partial_eta_sq generalised_eta_sq
#> 1      arm      0.4021595         0.17208342
#> 2     time      0.5749716         0.48314858
#> 3 arm:time      0.0214969         0.01495398
#>   time estimate        se df conf_low conf_high
#> 1   t1 13.04973 0.3999464 22 12.22030  13.87917
#> 2   t2 14.80474 0.2902370 22 14.20282  15.40665
#> 3   t3 17.29008 0.4441658 22 16.36894  18.21123
#>               term partial_eta_sq generalised_eta_sq
#> 1              arm    0.410002321       0.1720834176
#> 2          age_grp    0.029364143       0.0090485174
#> 3      arm:age_grp    0.002793737       0.0008379477
#> 4             time    0.587583344       0.4831485847
#> 5         arm:time    0.022614367       0.0149539780
#> 6     age_grp:time    0.006697070       0.0044236914
#> 7 arm:age_grp:time    0.044381926       0.0304721791
```
