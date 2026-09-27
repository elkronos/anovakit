# Kruskal-Wallis Test with Dunn Post-hoc Comparisons

Compares the distribution of a numeric response across groups without
assuming normality, using
[`kruskal.test`](https://rdrr.io/r/stats/kruskal.test.html), followed by
Dunn's pairwise rank-sum comparisons with a multiplicity adjustment. An
ordered factor must be converted first, with
[`as.integer()`](https://rdrr.io/r/base/integer.html).

## Usage

``` r
anova_kw(
  data,
  response,
  groups,
  conf_level = 0.95,
  adjust = "BH",
  diagnostics = FALSE,
  posthoc = TRUE,
  plots = TRUE,
  verbose = FALSE
)
```

## Arguments

- data:

  A data frame, or anything inheriting from one, such as a `data.table`
  or a tibble.

- response:

  Character. Name of the numeric response column.

- groups:

  Character vector. One or more grouping columns.

- conf_level:

  Numeric in (0, 1). Used for the group summary intervals. Default
  `0.95`.

- adjust:

  Character. Multiplicity adjustment for Dunn's comparisons, passed to
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html). One of `"BH"`,
  `"holm"`, `"hochberg"`, `"hommel"`, `"bonferroni"`, `"BY"`, `"fdr"` or
  `"none"`, spelled out in full. Default `"BH"`. Note that `"tukey"` is
  not available here, as it is on the functions whose comparisons come
  from emmeans: Dunn's test is a rank comparison, not a linear contrast.

- diagnostics:

  Logical. Compute the optional normality diagnostics into
  `$assumptions$normality`, with a Q-Q plot. They are contextual only
  and are skipped by default, since Kruskal-Wallis makes no
  distributional assumption. Default `FALSE`.

- posthoc:

  Logical. Compute Dunn's pairwise comparisons. There are `choose(k, 2)`
  of them, so this is worth turning off when the number of groups is
  large. Default `TRUE`.

- plots:

  Logical. Build ggplot2 objects. They are returned in `$plots`, never
  drawn. Default `TRUE`.

- verbose:

  Logical. Emit progress through
  [`message`](https://rdrr.io/r/base/message.html). Default `FALSE`.

## Value

An
[`anovakit_fit`](https://elkronos.github.io/anovakit/reference/anovakit_fit.md)
object.

- `$anova`: the Kruskal-Wallis test, one row. Its `term` is the grouping
  column, or `"A x B cells"` when several were combined.

- `$posthoc`: Dunn's comparisons, one row per pair of groups, with
  `group1`, `group2`, `mean_rank_diff` (the mean rank of `group1` minus
  that of `group2`), `z` (the same difference over its standard error,
  so with the same sign), `p_value` (two-sided, unadjusted),
  `p_adjusted` (adjusted by `adjust`) and `adjustment`: the same p-value
  columns as every other function in the package. `rstatix::dunn_test()`
  reports `group2` minus `group1`, so its statistic has the opposite
  sign.

- `$effect_sizes`: `epsilon_squared` and `eta_squared_H`, defined in
  Details.

- `$emmeans`: one row per group with its size, median, mean rank (the
  quantity Dunn's comparisons are built from), mean, standard deviation
  and t interval for the mean. These are raw group summaries; no model
  is fitted.

- `$assumptions$normality`: present only when `diagnostics = TRUE`.

- `$model`: the `htest` from
  [`kruskal.test`](https://rdrr.io/r/stats/kruskal.test.html), computed
  on the ranks; its `data.name` names the response and the grouping
  column(s).

## Details

Dunn's test is computed directly from the rank sums, with the standard
tie correction, so the package does not depend on an external
implementation and the full range of
[`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) methods is
available.

**Ties.** Values are ranked exactly as stored, and ties are counted on
those ranks, both for the omnibus test and for Dunn's comparisons, so
the tie correction always matches the ranks used. Two values that differ
only through floating-point error (`0.1 + 0.2` and `0.3`, say) are
therefore distinct, not tied; `$notes` says so when the data contain
such values. Round the response first if they are meant to be equal.
([`kruskal.test`](https://rdrr.io/r/stats/kruskal.test.html) itself
counts ties on the printed values, so in that one situation its
statistic differs slightly from the one reported here, which equals
[`kruskal.test()`](https://rdrr.io/r/stats/kruskal.test.html) applied to
the ranks.)

**Infinite values** of the response are kept, not dropped: they rank as
the most extreme observations, which is how a rank analysis treats, for
example, non-completers coded as `Inf`. A group containing one has an
infinite (or, with both signs, undefined) mean and no standard
deviation, so those columns of `$emmeans` are reported as `Inf`, `-Inf`
or `NA`. Its mean rank is exact, as is its median unless that falls
between `-Inf` and `Inf` (then `NA`). The box plot cannot draw infinite
values. Missing values are dropped and counted as usual.

**Effect sizes.** With \\H\\ the tie-corrected statistic, \\n\\
observations and \\k\\ groups, `epsilon_squared` is \\H / (n - 1)\\ and
`eta_squared_H` is \\(H - k + 1) / (n - k)\\, the names used by Tomczak
and Tomczak (2014) and effectsize. The first is exactly the proportion
of rank variance between groups (the \\R^2\\ of a one-way analysis of
variance on the ranks), and the second is its bias-adjusted counterpart,
so `epsilon_squared` is never smaller than `eta_squared_H`: the reverse
of what the same names mean for a parametric analysis. A negative
`eta_squared_H` (less between-group variation than chance produces) is
reported as 0, with a note giving the unfloored value; with one
observation per group (\\n = k\\) it is undefined and reported as `NA`,
also with a note.

Normality is not required and is not tested. When `diagnostics = TRUE` a
Q-Q plot of within-group standardised deviations (each value minus its
group mean, divided by its group standard deviation) is produced for
context only; it plays no part in the inference.

When several grouping variables are supplied they are combined into a
single cell factor. The result is then one omnibus test across all
populated cells, labelled `"A x B cells"` in `$anova`: it cannot
separate a main effect of one factor from a main effect of another, and
it cannot test an interaction. A note to that effect is added to the
returned object. For a factorial rank-based analysis, consider the
Scheirer-Ray-Hare extension or an aligned rank transform.

A response in which every value is the same is refused: there is nothing
to rank.

## References

Cohen, B. H. (2008). *Explaining Psychological Statistics* (3rd ed.).
Wiley.

Dunn, O. J. (1964). Multiple comparisons using rank sums.
*Technometrics*, 6(3), 241-252.

Tomczak, M., & Tomczak, E. (2014). The need to report effect size
estimates revisited. An overview of some recommended measures of effect
size. *Trends in Sport Sciences*, 21(1), 19-25.

## See also

[`anova_welch`](https://elkronos.github.io/anovakit/reference/anova_welch.md)
for the parametric equivalent.

## Examples

``` r
set.seed(123)
d <- data.frame(
  group = rep(c("A", "B", "C"), each = 25),
  value = c(rnorm(25, 5), rnorm(25, 7, 1.5), rnorm(25, 4, 0.8))
)
fit <- anova_kw(d, "value", "group")
fit
#> Kruskal-Wallis rank sum test 
#> ----------------------------
#> Call: anova_kw(data = d, response = "value", groups = "group")
#> Observations used: 75
#> 
#> Omnibus test
#>    term statistic df   p_value
#> 1 group     45.71  2 1.187e-10
#> 
#> Plots available: box
#>   (use plot(x, which = "box"))
fit$posthoc
#>   group1 group2 mean_rank_diff         z      p_value   p_adjusted adjustment
#> 1      A      B         -24.56 -3.984158 6.771978e-05 1.015797e-04         BH
#> 2      A      C          16.88  2.738298 6.175816e-03 6.175816e-03         BH
#> 3      B      C          41.44  6.722456 1.786871e-11 5.360613e-11         BH

# Epsilon squared and eta squared for the rank statistic, and the group
# summary the comparisons are built from
fit$effect_sizes
#>           measure  estimate
#> 1 epsilon_squared 0.6176865
#> 2   eta_squared_H 0.6070667
fit$emmeans
#>   group  n   median mean_rank     mean        sd        se conf_low conf_high
#> 1     A 25 4.782025     35.44 4.966670 0.9467324 0.1893465 4.575878  5.357462
#> 2     B 25 6.907132     60.00 7.153206 1.3783100 0.2756620 6.584268  7.722145
#> 3     C 25 4.042403     18.56 4.008193 0.7779371 0.1555874 3.687076  4.329309

# Any p.adjust() method may be used; the default is "BH"
anova_kw(d, "value", "group", adjust = "bonferroni",
         plots = FALSE)$posthoc$p_adjusted
#> [1] 2.031593e-04 1.852745e-02 5.360613e-11

# Optional normality diagnostics, which play no part in the test itself
anova_kw(d, "value", "group", diagnostics = TRUE,
         plots = FALSE)$assumptions$normality
#>   group  n statistic   p_value note
#> 1     A 25 0.9674463 0.5812306 <NA>
#> 2     B 25 0.9768585 0.8166478 <NA>
#> 3     C 25 0.9891279 0.9927964 <NA>

# Small samples work: there is no minimum group size
small <- data.frame(g = rep(c("a", "b"), each = 3), y = c(1, 2, 3, 4, 5, 6))
anova_kw(small, "y", "g", plots = FALSE)$anova
#>   term statistic df    p_value
#> 1    g  3.857143  1 0.04953461

# Infinite values are kept and ranked as the most extreme outcomes
worst <- data.frame(g = rep(c("a", "b"), each = 5),
                    y = c(1:5, 6, 7, Inf, Inf, Inf))
anova_kw(worst, "y", "g", plots = FALSE)$anova
#>   term statistic df     p_value
#> 1    g  6.987578  1 0.008207736
```
