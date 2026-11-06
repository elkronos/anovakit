# Weights and robust standard errors

``` r

library(anovakit)
```

Three functions take prior weights and four take a robust covariance.

## Which functions take weights

| Function | `weights` | `vcov_type` |
|----|----|----|
| [`anova_bin()`](https://elkronos.github.io/anovakit/reference/anova_bin.md) | yes | yes |
| [`anova_count()`](https://elkronos.github.io/anovakit/reference/anova_count.md) | yes | yes |
| [`anova_glm()`](https://elkronos.github.io/anovakit/reference/anova_glm.md) | yes | yes |
| [`anova_ancova()`](https://elkronos.github.io/anovakit/reference/anova_ancova.md) | no | yes |

`weights` names a numeric column of the data, not a vector. The weights
are then subsetted with the rows, so they stay aligned when a row is
dropped for a missing value; a vector evaluated against the original
data would not. The column must be numeric and non-negative, and it
cannot also be the response, a grouping column or the offset. A missing
or infinite weight is treated like a missing value in any other analysed
column: the row is dropped and counted in `$n_removed`.

## How the weights are read

What a weight means depends on the weights and on the response. The
package decides which reading applies and says so in `$notes` when the
reading changes a result: the Pearson dispersion in
[`anova_count()`](https://elkronos.github.io/anovakit/reference/anova_count.md),
the robust covariance, and weights that are not whole numbers.

| The weights are | On | Read as | Each row is |
|----|----|----|----|
| whole numbers, some above 1 | a 0/1 response (binomial, quasi-binomial) or counts (Poisson, quasi-Poisson, negative binomial) | frequency weights | that many identical observations |
| whole numbers | a proportion (binomial) | numbers of trials | one binomial observation |
| not all whole numbers, or any weights in another family | any response | weights used as given, as [`glm()`](https://rdrr.io/r/stats/glm.html) uses them | one unit |

[`anova_bin()`](https://elkronos.github.io/anovakit/reference/anova_bin.md)
always has a 0/1 response and
[`anova_count()`](https://elkronos.github.io/anovakit/reference/anova_count.md)
always has counts, so in those two functions whole-number weights are
always frequencies.
[`anova_glm()`](https://elkronos.github.io/anovakit/reference/anova_glm.md)
also takes proportions, where whole-number weights are trials, and the
Gaussian, Gamma, inverse Gaussian and
[`quasi()`](https://rdrr.io/r/stats/family.html) families, where they
are precision weights. In both cases each row stays one unit.

## Frequency weights

A frequency weight says how many identical observations a row stands
for. `UCBAdmissions` is a table of 4,526 applications to graduate school
at Berkeley in 1973, by gender, department and outcome. As a data frame
it has one row per cell and the count in `Freq`:

``` r

ucb <- as.data.frame(UCBAdmissions)
head(ucb)
#>      Admit Gender Dept Freq
#> 1 Admitted   Male    A  512
#> 2 Rejected   Male    A  313
#> 3 Admitted Female    A   89
#> 4 Rejected Female    A   19
#> 5 Admitted   Male    B  353
#> 6 Rejected   Male    B  207
c(rows = nrow(ucb), applicants = sum(ucb$Freq))
#>       rows applicants 
#>         24       4526
```

For comparison, here is the same data with one row per applicant:

``` r

applicants <- ucb[rep(seq_len(nrow(ucb)), ucb$Freq),
                  c("Admit", "Gender", "Dept")]
row.names(applicants) <- NULL
nrow(applicants)
#> [1] 4526
```

The 24-row table with `weights = "Freq"` and the 4,526 rows without
weights give the same analysis:

``` r

from_table <- anova_bin(ucb, "Admit", c("Gender", "Dept"), weights = "Freq",
                        success = "Admitted", plots = FALSE)
from_rows <- anova_bin(applicants, "Admit", c("Gender", "Dept"),
                       success = "Admitted", plots = FALSE)
all.equal(from_table$anova, from_rows$anova)
#> [1] TRUE
data.frame(comparison = from_table$effect_sizes$comparison,
           or_table = from_table$effect_sizes$odds_ratio,
           or_rows = from_rows$effect_sizes$odds_ratio,
           se_table = from_table$effect_sizes$se_log_or,
           se_rows = from_rows$effect_sizes$se_log_or)
#>       comparison   or_table    or_rows   se_table    se_rows
#> 1 Female vs Male 1.10502735 1.10502735 0.08084646 0.08084644
#> 2         B vs A 0.95753028 0.95753028 0.10983890 0.10983890
#> 3         C vs A 0.28291804 0.28291804 0.10663289 0.10663288
#> 4         D vs A 0.27400567 0.27400567 0.10582342 0.10582342
#> 5         E vs A 0.17564230 0.17564230 0.12611350 0.12611349
#> 6         F vs A 0.03664494 0.03664494 0.16998176 0.16998115
c(table = from_table$model_stats$BIC, rows = from_rows$model_stats$BIC)
#>    table     rows 
#> 5246.412 5246.412
```

The odds ratios and standard errors agree to about six significant
figures; the rest is [`glm()`](https://rdrr.io/r/stats/glm.html)’s
convergence tolerance. BIC agrees too, because with frequency weights
its penalty uses the 4,526 observations the weights represent, not the
24 rows. The marginal means, the comparisons and the observed
proportions in `$assumptions$proportions` match in the same way.

Counts work the same way in
[`anova_count()`](https://elkronos.github.io/anovakit/reference/anova_count.md).
`InsectSprays` has 72 insect counts, 12 for each of six sprays.
Tabulated, it becomes one row per spray and count, with the frequency in
`Freq`:

``` r

spray_tab <- as.data.frame(table(spray = InsectSprays$spray,
                                 count = InsectSprays$count),
                           stringsAsFactors = FALSE)
spray_tab$count <- as.numeric(spray_tab$count)
c(rows = nrow(spray_tab), nonzero_rows = sum(spray_tab$Freq > 0),
  observations = sum(spray_tab$Freq))
#>         rows nonzero_rows observations 
#>          144           43           72
```

``` r

spray_rows <- anova_count(InsectSprays, "count", "spray", plots = FALSE)
spray_freq <- anova_count(spray_tab, "count", "spray", weights = "Freq",
                          plots = FALSE)
c(rows = spray_rows$model_type, table = spray_freq$model_type)
#>     rows    table 
#> "negbin" "negbin"
c(rows = spray_rows$dispersion, table = spray_freq$dispersion)
#>     rows    table 
#> 1.507713 1.507713
all.equal(spray_rows$anova, spray_freq$anova)
#> [1] TRUE
all.equal(spray_rows$effect_sizes, spray_freq$effect_sizes)
#> [1] TRUE
all.equal(spray_rows$emmeans, spray_freq$emmeans)
#> [1] TRUE
grep("frequency weights", spray_freq$notes, value = TRUE)
#> [1] "The weights in `Freq` are whole numbers, so they are treated as frequency weights: each row stands for that many observations. The Pearson dispersion divides by the 72 observations they represent minus the 6 parameters, not by the 43 rows, which would inflate it."
```

The Pearson dispersion, on which `model = "auto"` chooses between
Poisson and negative binomial, divides by the observations minus the
parameters. Divided by the rows of the table instead, it would be
larger, and could switch the model for no reason in the data. In
[`anova_glm()`](https://elkronos.github.io/anovakit/reference/anova_glm.md)
the residual degrees of freedom are counted in observations in the same
way, so a quasi-Poisson or quasi-binomial fit, an F test and a t
interval on a frequency table equal those of the expanded data. `$model`
is still the `glm` that
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) fitted to the rows,
and `$notes` says when its own
[`summary()`](https://rdrr.io/r/base/summary.html) would differ.

The [counts
walkthrough](https://elkronos.github.io/anovakit/articles/counts.md)
shows the `$assumptions$cell_counts` table and the `observed` plot for a
frequency table.

## Zero weights

`spray_tab` above still holds the combinations of spray and count that
never occurred, with a frequency of zero. A row with weight zero
contributes nothing to the fit, so it is dropped before fitting and
counted in `$n_removed`:

``` r

c(supplied = nrow(spray_tab), used = nrow(spray_freq$data_used),
  removed = spray_freq$n_removed)
#> supplied     used  removed 
#>      144       43      101
grep("zero", spray_freq$notes, value = TRUE)
#> [1] "Dropped 101 row(s) whose weight in `Freq` is zero: they contribute nothing to the fit, and keeping them would make row counts, robust standard errors and information criteria disagree with it."
```

Keeping those rows would make the row counts, the robust standard errors
and the information criteria disagree with the fit. As for every
function, `nrow($data_used) + $n_removed` is the number of rows
supplied.

A group whose weights are all zero is dropped from the analysis, and the
note names it:

``` r

no_f <- spray_tab
no_f$Freq[no_f$spray == "F"] <- 0
grep("Level", anova_count(no_f, "count", "spray", weights = "Freq",
                          plots = FALSE)$notes, value = TRUE)
#> [1] "Level \"F\" of `spray` had no rows with a non-zero weight in `Freq`, so it was dropped from the analysis."
```

When every weight is zero, or only one row has a non-zero weight, the
function stops with an error.

## Trials behind proportions

A binomial response can be given as a proportion, with the number of
trials behind each proportion as the weight. That form is for
[`anova_glm()`](https://elkronos.github.io/anovakit/reference/anova_glm.md):
[`anova_bin()`](https://elkronos.github.io/anovakit/reference/anova_bin.md)
takes a 0/1 response only. `esoph` is a case-control study of
oesophageal cancer, with the numbers of cases and controls in each of 88
combinations of age group, alcohol and tobacco consumption:

``` r

es <- esoph
es$agegp <- factor(es$agegp, ordered = FALSE)
es$alcgp <- factor(es$alcgp, ordered = FALSE)
es$trials <- es$ncases + es$ncontrols
es$prop <- es$ncases / es$trials
head(es[, c("agegp", "alcgp", "ncases", "trials", "prop")])
#>   agegp     alcgp ncases trials prop
#> 1 25-34 0-39g/day      0     40    0
#> 2 25-34 0-39g/day      0     10    0
#> 3 25-34 0-39g/day      0      6    0
#> 4 25-34 0-39g/day      0      5    0
#> 5 25-34     40-79      0     27    0
#> 6 25-34     40-79      0      7    0
```

(The factors are made unordered so that the coefficients are comparisons
with the first level rather than polynomial trends.) With one row per
subject instead, the same data have 975 rows of a 0/1 response:

``` r

per_subject <- rbind(
  data.frame(agegp = rep(es$agegp, es$ncases),
             alcgp = rep(es$alcgp, es$ncases), case = 1),
  data.frame(agegp = rep(es$agegp, es$ncontrols),
             alcgp = rep(es$alcgp, es$ncontrols), case = 0))
nrow(per_subject)
#> [1] 975
```

``` r

as_prop <- anova_glm(es, "prop", c("agegp", "alcgp"), family = binomial(),
                     weights = "trials", plots = FALSE)
as_rows <- anova_glm(per_subject, "case", c("agegp", "alcgp"),
                     family = binomial(), plots = FALSE)
all.equal(as_prop$anova, as_rows$anova)
#> [1] TRUE
data.frame(term = as_prop$emmeans$term,
           level = ifelse(is.na(as_prop$emmeans$agegp),
                          as.character(as_prop$emmeans$alcgp),
                          as.character(as_prop$emmeans$agegp)),
           prob_prop = as_prop$emmeans$estimate,
           prob_rows = as_rows$emmeans$estimate)
#>     term     level  prob_prop  prob_rows
#> 1  agegp     25-34 0.01253232 0.01253232
#> 2  agegp     35-44 0.06089869 0.06089869
#> 3  agegp     45-54 0.28070062 0.28070062
#> 4  agegp     55-64 0.39570900 0.39570900
#> 5  agegp     65-74 0.49748585 0.49748585
#> 6  agegp       75+ 0.51434499 0.51434499
#> 7  alcgp 0-39g/day 0.03978680 0.03978680
#> 8  alcgp     40-79 0.14813421 0.14813421
#> 9  alcgp    80-119 0.23568034 0.23568034
#> 10 alcgp      120+ 0.62161114 0.62161114
```

The likelihood-ratio tests and the marginal means are the same. The
deviance explained is not: it is measured against the saturated model,
which differs between the two forms, and a note says so. McFadden’s R
squared is computed from the Bernoulli log-likelihood and is the same
either way.

The trials are not frequencies, because each row is one binomial
observation rather than a count of identical ones. The Pearson
dispersion of the proportions carries information about overdispersion,
which a 0/1 response cannot show, so
[`anova_glm()`](https://elkronos.github.io/anovakit/reference/anova_glm.md)
checks it:

``` r

as_prop$assumptions$dispersion
#> [1] 1.487358
grep("Pearson dispersion", as_prop$notes, value = TRUE)
#> [1] "Pearson dispersion is 1.49 (Pearson test of overdispersion: p = 0.0032), where a binomial model assumes 1: its test statistics are inflated by about that factor, its standard errors are about 1.22 times too small, and a nominal 5% test on 1 degree of freedom rejects roughly 11% of true null hypotheses, and one on 5 (the largest term here) roughly 19%. Consider family = \"quasibinomial\", or a beta-binomial model."
```

A robust covariance also treats each of the 88 rows as one unit, so it
allows for the extra variation between rows that the dispersion of 1.49
measures. The per-subject form cannot see that variation, and its robust
standard errors stay close to the model-based ones:

``` r

prop_hc3 <- anova_glm(es, "prop", c("agegp", "alcgp"), family = binomial(),
                      weights = "trials", vcov_type = "HC3", plots = FALSE)
rows_hc3 <- anova_glm(per_subject, "case", c("agegp", "alcgp"),
                      family = binomial(), vcov_type = "HC3", plots = FALSE)
data.frame(term = as_prop$emmeans$term,
           se_model = as_prop$emmeans$se,
           se_hc3_prop = prop_hc3$emmeans$se,
           se_hc3_rows = rows_hc3$emmeans$se)
#>     term   se_model se_hc3_prop se_hc3_rows
#> 1  agegp 0.01262372  0.01385063 0.012971823
#> 2  agegp 0.02058721  0.02866157 0.020267469
#> 3  agegp 0.03880783  0.04939598 0.036329909
#> 4  agegp 0.03980602  0.06202781 0.043375963
#> 5  agegp 0.04952135  0.06751124 0.052921400
#> 6  agegp 0.09561464  0.12208635 0.093769792
#> 7  alcgp 0.01013718  0.01385172 0.009974735
#> 8  alcgp 0.02836986  0.03417081 0.028016447
#> 9  alcgp 0.04783171  0.05382355 0.046929813
#> 10 alcgp 0.08038267  0.11456357 0.093728591
```

Trials do not inflate the tests the way unequal sampling weights do
(below): they set the variance of each proportion, as precision weights
do.

## Sampling and other weights that are not whole numbers

Weights that are not whole numbers are used as given, as
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) uses prior weights:
each row’s contribution to the likelihood is multiplied by its weight,
as if it were that many observations. The fit, the tests, the intervals
and the marginal means all treat them that way. The robust covariance
treats each row as one unit instead, and two checks rescale the weights
to average 1 first: the sparse-data check of
[`anova_bin()`](https://elkronos.github.io/anovakit/reference/anova_bin.md)
and the Pearson dispersion on which
[`anova_count()`](https://elkronos.github.io/anovakit/reference/anova_count.md)
chooses its model.

Suppose the applicants above came from a survey that sampled men at a
higher rate than women, so that each woman’s response is weighted up and
each man’s down. The weights here are made up for illustration and
average 1:

``` r

applicants$w <- c(Male = 0.4, Female = 2.2)[as.character(applicants$Gender)]
applicants$w <- applicants$w / mean(applicants$w)
tapply(applicants$w, applicants$Gender, unique)
#>      Male    Female 
#> 0.3540501 1.9472758
```

### Their scale matters

Weights used as given have no natural scale, and the tests scale with
them. The same weights divided by 100 make the 4,526 applicants count as
about 45, and every likelihood-ratio statistic falls by the same factor:

``` r

sw <- anova_bin(applicants, "Admit", c("Gender", "Dept"), weights = "w",
                success = "Admitted", plots = FALSE)
applicants$w_small <- applicants$w / 100
sw_small <- anova_bin(applicants, "Admit", c("Gender", "Dept"),
                      weights = "w_small", success = "Admitted",
                      plots = FALSE)
data.frame(term = sw$anova$term, mean_1 = sw$anova$statistic,
           mean_0.01 = sw_small$anova$statistic)
#>     term     mean_1  mean_0.01
#> 1 Gender   6.100007 0.06100007
#> 2   Dept 684.795658 6.84795658
grep("average", sw_small$notes, value = TRUE)
#> [1] "The weights in `w_small` are not all whole numbers, so they are not read as frequencies, and they average 0.01 rather than 1. The fit, the tests and the intervals use them as given, as glm() uses prior weights, so their scale matters: multiplying every weight by a constant multiplies each likelihood-ratio and Wald statistic by about that constant. For sampling weights, rescale them to average 1, or use a design-based method such as survey::svyglm()."
#> [2] "The grouping factors enter the model additively, so $emmeans and $posthoc are reported for each factor separately, averaged over the others (column `term`), and each factor's comparisons are a separate multiplicity family. Set interaction = TRUE to compare cells instead."
```

[`anova_bin()`](https://elkronos.github.io/anovakit/reference/anova_bin.md),
[`anova_count()`](https://elkronos.github.io/anovakit/reference/anova_count.md)
and a Poisson or binomial
[`anova_glm()`](https://elkronos.github.io/anovakit/reference/anova_glm.md)
give that note whenever such weights do not average 1 (within about
10%). Rescale sampling weights to average 1, as above, or analyse them
with a design-based method such as `survey::svyglm()`.

In a binomial model, weights that are not whole numbers also make AIC,
BIC and the log-likelihood in `$model_stats` `NA`, with a note, because
R’s binomial log-likelihood rounds the weights. McFadden’s R squared and
the model-versus-null test use the deviances and are reported.
[`glm()`](https://rdrr.io/r/stats/glm.html)’s warning about “non-integer
\#successes” is recorded in `$notes` rather than printed.

In
[`anova_count()`](https://elkronos.github.io/anovakit/reference/anova_count.md)
the Pearson dispersion in `$dispersion` and `$model_dispersion` is
computed with the weights rescaled to average 1, so the choice between
Poisson and negative binomial does not depend on their scale. In
[`anova_glm()`](https://elkronos.github.io/anovakit/reference/anova_glm.md)
the dispersion is that of the model as weighted, so in a Poisson or
binomial fit it scales with the weights, and so does the overdispersion
check; the note on the weights’ scale says so. A family that estimates
its dispersion (Gaussian, Gamma, quasi families) absorbs the scale of
the weights.

### Unequal weights distort the tests

Weights that average 1 still distort the tests when they are unequal:
some terms are inflated, others deflated. A likelihood-ratio or Wald
test treats a woman’s row, of weight 1.95, as that many independent
observations, while it is one person. How much that moves a test depends
on where the heavy rows fall, so `$notes` computes the factor for each
term from the sandwich covariance:

``` r

sw$anova
#>     term df  statistic       p_value
#> 1 Gender  1   6.100007  1.351813e-02
#> 2   Dept  5 684.795658 9.518615e-146
grep("unequal", sw$notes, value = TRUE)
#> [1] "The weights in `w` are unequal and not frequencies, so the likelihood-ratio tests in $anova are inflated by between about 0.84 and 1.45, depending on the term (the weights are used as given and are not frequencies; computed from the sandwich covariance, assuming the weights are not inverse variances, which would not inflate them): a nominal 5% test of Dept rejects roughly 18% of true null hypotheses. The standard errors in the other tables are too small for the same reason; vcov_type = \"HC3\" gives ones that are not."
```

The calculation assumes the weights are not inverse variances; precision
weights, such as the trials behind a proportion, would not inflate the
tests. The note appears in
[`anova_bin()`](https://elkronos.github.io/anovakit/reference/anova_bin.md),
[`anova_count()`](https://elkronos.github.io/anovakit/reference/anova_count.md)
and
[`anova_glm()`](https://elkronos.github.io/anovakit/reference/anova_glm.md)
for Poisson, quasi-Poisson, negative binomial, binomial and
quasi-binomial fits (on a binomial, only for weights that are not whole
numbers on a 0/1 response). It never appears for Gaussian, Gamma or
inverse Gaussian fits, where weights are read as precision weights. It
is given only when some test would reject 7.5% or more of true null
hypotheses at a nominal 5%. Under the quasi families’ F tests the
estimated dispersion removes the average inflation but not the part due
to unequal weights, and the note says so.

The remedy the note names is a robust covariance, and an omnibus test
that uses it:

``` r

sw_hc3 <- anova_bin(applicants, "Admit", c("Gender", "Dept"), weights = "w",
                    success = "Admitted", vcov_type = "HC3",
                    test_statistic = "Wald", plots = FALSE)
sw_hc3$anova
#>     term df  statistic      p_value
#> 1 Gender  1   7.388409 6.564554e-03
#> 2   Dept  5 396.687392 1.535320e-83
data.frame(comparison = sw$effect_sizes$comparison,
           se_model = sw$effect_sizes$se_log_or,
           se_hc3 = sw_hc3$effect_sizes$se_log_or)
#>       comparison  se_model     se_hc3
#> 1 Female vs Male 0.1011743 0.09088179
#> 2         B vs A 0.1663584 0.14146922
#> 3         C vs A 0.1261355 0.13566444
#> 4         D vs A 0.1293345 0.13779505
#> 5         E vs A 0.1381565 0.15597408
#> 6         F vs A 0.1775987 0.21164049
grep("unequal", sw_hc3$notes, value = TRUE)
#> character(0)
```

Dept, the term the note names as the worst affected (its factor is the
larger one), falls from 685 to 397; Gender rises from 6.10 to 7.39.
Those are the directions of the note’s factors, but the sizes differ
from them, because the comparison also moves from a likelihood-ratio to
a Wald statistic. Most standard errors grow, and the inflation note is
gone. Weights alone do not describe strata or clusters; a design-based
analysis with survey accounts for those too.

## Robust standard errors with `vcov_type`

`vcov_type = "model"` (the default) uses the model’s own covariance.
`"HC0"` to `"HC4"` use a heteroscedasticity-consistent (sandwich)
covariance from the sandwich package, which does not rely on the
family’s variance function or on equal residual variances. HC0 is
White’s estimator; HC1 to HC4 refine it for small samples, HC1 by a
degrees-of-freedom factor and HC2 to HC4 by each row’s leverage. HC3,
the default of
[`sandwich::vcovHC()`](https://zeileis.codeberg.page/sandwich/reference/vcovHC.html)
and the type used throughout this guide, is the usual choice for small
or unbalanced groups.

The robust covariance does not reach every table:

| Result | [`anova_ancova()`](https://elkronos.github.io/anovakit/reference/anova_ancova.md) | [`anova_bin()`](https://elkronos.github.io/anovakit/reference/anova_bin.md) | [`anova_count()`](https://elkronos.github.io/anovakit/reference/anova_count.md) | [`anova_glm()`](https://elkronos.github.io/anovakit/reference/anova_glm.md) |
|----|----|----|----|----|
| Coefficients | no coefficient table | `$effect_sizes` (odds ratios) | `$effect_sizes` (rate ratios) | `$assumptions$coefficients` |
| `$emmeans` and `$posthoc` | yes | yes | yes | yes |
| Omnibus `$anova` | yes, as Wald F tests | with `test_statistic = "Wald"` | with `test_statistic = "Wald"` | with `test_statistic = "Wald"` |
| `$simple_slopes` | yes | no such table | no such table | no such table |
| Not reached | `$effect_sizes`, `$slopes_test` | `$anova` unless `test_statistic = "Wald"` | `$anova` unless `test_statistic = "Wald"`; `$dispersion` and the model choice | `$anova` unless `test_statistic = "Wald"`; `$effect_sizes` |

In
[`anova_bin()`](https://elkronos.github.io/anovakit/reference/anova_bin.md),
[`anova_count()`](https://elkronos.github.io/anovakit/reference/anova_count.md)
and
[`anova_glm()`](https://elkronos.github.io/anovakit/reference/anova_glm.md),
a likelihood-ratio or F test compares deviances, so it rests on the
model-based variance whatever `vcov_type` says; set
`test_statistic = "Wald"` for an omnibus test that uses the robust
covariance. `$notes` says which tests use it. Coefficient intervals in
[`anova_bin()`](https://elkronos.github.io/anovakit/reference/anova_bin.md)
and
[`anova_glm()`](https://elkronos.github.io/anovakit/reference/anova_glm.md)
become Wald intervals built from it, since a profile likelihood cannot
use a sandwich covariance.

On `warpbreaks` (numbers of warp breaks per loom, 9 looms in each of six
cells), the default likelihood-ratio table is unchanged by `vcov_type`,
while a Wald table uses it:

``` r

breaks_model <- anova_count(warpbreaks, "breaks", c("wool", "tension"),
                            interaction = TRUE, plots = FALSE)
breaks_lr <- anova_count(warpbreaks, "breaks", c("wool", "tension"),
                         interaction = TRUE, vcov_type = "HC3",
                         plots = FALSE)
breaks_wald <- anova_count(warpbreaks, "breaks", c("wool", "tension"),
                           interaction = TRUE, vcov_type = "HC3",
                           test_statistic = "Wald", plots = FALSE)
all.equal(breaks_model$anova, breaks_lr$anova)
#> [1] TRUE
breaks_wald$anova
#>           term df statistic      p_value
#> 1         wool  1  2.526565 0.1119436981
#> 2      tension  2 16.808766 0.0002238839
#> 3 wool:tension  2  6.765016 0.0339621730
data.frame(comparison = breaks_model$effect_sizes$comparison,
           low_model = breaks_model$effect_sizes$conf_low,
           high_model = breaks_model$effect_sizes$conf_high,
           low_hc3 = breaks_lr$effect_sizes$conf_low,
           high_hc3 = breaks_lr$effect_sizes$conf_high)
#>              comparison low_model high_model   low_hc3  high_hc3
#> 1 B vs A at tension = L 0.4651349  0.8625808 0.4369807 0.9181560
#> 2    M vs L at wool = A 0.3938590  0.7366784 0.3696585 0.7849066
#> 3    H vs L at wool = A 0.4032354  0.7532465 0.3679333 0.8255184
#> 4   (B vs A) x (M vs L) 1.2122753  2.9560520 1.1459874 3.1270404
#> 5   (B vs A) x (H vs L) 0.7667100  1.9009866 0.7290004 1.9993205
grep("omnibus", breaks_lr$notes, value = TRUE)
#> [1] "The omnibus likelihood-ratio test compares deviances, so it rests on the model-based variance; the robust (HC3) covariance is used only by $effect_sizes, the marginal means and the comparisons. For an omnibus test that uses it, set test_statistic = \"Wald\"."
```

In
[`anova_ancova()`](https://elkronos.github.io/anovakit/reference/anova_ancova.md)
the omnibus table itself changes: its F tests become Wald F tests built
from the robust covariance (through `car::Anova(vcov. = )`), on the
model’s residual degrees of freedom, so the table has no sums of
squares. Its `statistic` attribute, which
[`print()`](https://rdrr.io/r/base/print.html) shows in the heading,
says so. [`MASS::anorexia`](https://rdrr.io/pkg/MASS/man/anorexia.html)
has the weights of 72 young women before and after one of three
treatments:

``` r

an <- MASS::anorexia
an_hc3 <- anova_ancova(an, "Postwt", "Treat", "Prewt",
                       force_interaction = TRUE, vcov_type = "HC3",
                       plots = FALSE)
an_hc3$anova
#>          term df statistic     p_value
#> 1       Prewt  1  4.733147 0.033172900
#> 2       Treat  2  7.529746 0.001133501
#> 3 Prewt:Treat  2  3.211929 0.046650376
#> 4   Residuals 66        NA          NA
attr(an_hc3$anova, "statistic")
#> [1] "Wald F with HC3 covariance"
grep("Wald F", an_hc3$notes, value = TRUE)
#> [1] "vcov_type = \"HC3\": the F tests in $anova are Wald F tests built from that covariance (car::Anova(vcov. = )), on the model's residual degrees of freedom, so the table has no sums of squares. $effect_sizes are still computed from the model-based sums of squares, which the covariance does not change, and the homogeneity-of-slopes test in $slopes_test is the model-based F test."
```

`$effect_sizes` are still computed from the model-based sums of squares,
which the covariance does not change, and the homogeneity-of-slopes test
stays the model-based F test. In
[`anova_ancova()`](https://elkronos.github.io/anovakit/reference/anova_ancova.md)
and in a Gaussian
[`anova_glm()`](https://elkronos.github.io/anovakit/reference/anova_glm.md),
comparisons that use robust standard errors keep the model’s residual
degrees of freedom, and a note warns that with very small groups of
unequal variance they can be anti-conservative; for a one-way design
without covariates,
[`anova_welch()`](https://elkronos.github.io/anovakit/reference/anova_welch.md)
is better calibrated.

### Robust standard errors and weights

With frequency weights the sandwich is computed per observation, as it
would be on the data expanded to one row per observation.
[`sandwich::vcovHC()`](https://zeileis.codeberg.page/sandwich/reference/vcovHC.html)
applied to the weighted rows directly would count each row as one unit
with an inflated score, and overstate the standard errors:

``` r

table_hc3 <- anova_bin(ucb, "Admit", c("Gender", "Dept"), weights = "Freq",
                       success = "Admitted", vcov_type = "HC3",
                       plots = FALSE)
appl_hc3 <- anova_bin(applicants, "Admit", c("Gender", "Dept"),
                      success = "Admitted", vcov_type = "HC3",
                      plots = FALSE)
naive <- sqrt(diag(sandwich::vcovHC(from_table$model, type = "HC3")))[-1]
data.frame(comparison = table_hc3$effect_sizes$comparison,
           model_based = from_table$effect_sizes$se_log_or,
           hc3_table = table_hc3$effect_sizes$se_log_or,
           hc3_rows = appl_hc3$effect_sizes$se_log_or,
           hc3_row_by_row = unname(naive))
#>       comparison model_based  hc3_table   hc3_rows hc3_row_by_row
#> 1 Female vs Male  0.08084646 0.08045399 0.08045399       1.497407
#> 2         B vs A  0.10983890 0.11002565 0.11002565       3.751704
#> 3         C vs A  0.10663289 0.10631318 0.10631318       2.786957
#> 4         D vs A  0.10582342 0.10528448 0.10528448       2.821575
#> 5         E vs A  0.12611350 0.12648471 0.12648471       2.911143
#> 6         F vs A  0.16998176 0.16944692 0.16944691       2.843432
grep("frequency weights represent", table_hc3$notes, value = TRUE)
#> [1] "The robust covariance (HC3) treats each of the 4,526 observations the frequency weights represent as a unit, as it would on the data expanded to one row per observation; computed row by row it would count each row as one unit with an inflated score, and overstate the standard errors."
```

The HC3 standard errors from the table and from one row per applicant
are the same.
[`sandwich::vcovHC()`](https://zeileis.codeberg.page/sandwich/reference/vcovHC.html)
on the table is wrong by a factor of 17 to 34.

With any other weights (sampling weights, precision weights, trials) the
robust covariance treats each row as one independent unit, and `$notes`
says so. If such weights in fact count identical observations, expand
the data to one row per observation first.

## When the model-based covariance is used instead

In some fits a sandwich covariance is wrong, and the function uses the
model-based covariance instead, with a note saying why. That happens
when:

- the sandwich package is not installed;
- some fitted values are on the boundary: a group with no events or all
  events in a binomial model, a group of all-zero counts, or a separated
  coefficient, where the sandwich collapses towards zero and gives
  spuriously small p-values;
- a row has leverage 1, as a group or cell with a single observation
  does, so its residual is zero by construction and the sandwich would
  give a standard error of zero (HC0, HC1) or none at all (HC2 to HC4);
- the robust covariance fails or has non-finite entries.

A site where nothing happened puts the fitted rate on the boundary:

``` r

incidents <- data.frame(
  site = rep(c("north", "south", "east"), each = 8),
  events = c(2, 4, 3, 1, 5, 3, 2, 4,
             6, 3, 5, 7, 4, 6, 5, 8,
             0, 0, 0, 0, 0, 0, 0, 0))
grep("Robust", anova_count(incidents, "events", "site", vcov_type = "HC3",
                           plots = FALSE)$notes, value = TRUE)
#> [1] "Robust standard errors (HC3) were not used: some fitted values are on the boundary (a group with no events, or all events), where the sandwich covariance collapses towards zero and reports spuriously small p-values. Model-based standard errors were used instead."
```

A group with one observation has leverage 1:

``` r

small <- data.frame(g = c(rep("a", 6), rep("b", 6), "c"),
                    y = c(1:6, 3:8, 4))
grep("Robust", anova_glm(small, "y", "g", vcov_type = "HC3",
                         plots = FALSE)$notes, value = TRUE)
#> [1] "Robust standard errors (HC3) were not used: 1 row has leverage 1 (the model fits it exactly, as it does a group or cell with a single observation), so its residual is zero by construction and the sandwich would give the parameters it determines a standard error of zero (HC0, HC1) or none at all (HC2 to HC4). Model-based standard errors were used instead."
```

In
[`anova_ancova()`](https://elkronos.github.io/anovakit/reference/anova_ancova.md),
[`car::Anova()`](https://rdrr.io/pkg/car/man/Anova.html) refuses robust
Wald tests for a model with aliased coefficients (an empty cell of the
design, for example). `$anova` then keeps the model-based F tests, and
`$notes` says so, while `$emmeans` and `$posthoc` still use the robust
covariance.

## See also

- [Binary outcomes:
  `anova_bin()`](https://elkronos.github.io/anovakit/articles/binary.md),
  [Counts and rates:
  `anova_count()`](https://elkronos.github.io/anovakit/articles/counts.md)
  and [Any GLM family:
  `anova_glm()`](https://elkronos.github.io/anovakit/articles/glm.md)
  for weights in each function, including precision weights on averages
  in a Gamma model.
- [ANCOVA:
  `anova_ancova()`](https://elkronos.github.io/anovakit/articles/ancova.md)
  for robust tests alongside the homogeneity-of-slopes test.
- [Diagnostics](https://elkronos.github.io/anovakit/articles/diagnostics.md)
  for the overdispersion, separation and sparse-data checks mentioned
  here.
- [Working with
  results](https://elkronos.github.io/anovakit/articles/results.md) for
  `$notes`, `$data_used` and `$n_removed`.
