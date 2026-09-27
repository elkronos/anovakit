# Which function? Choosing an analysis

``` r

library(anovakit)
```

anovakit has eight analysis functions. They share one calling convention
and return the same kind of object, so moving between them costs
nothing. The choice between them is a choice about your *data*: what the
outcome is, how the observations are related, and what you need to
adjust for. This guide walks through that choice and links to a full
walkthrough of each function.

## Four questions

Answer these in order. The first answer that settles it wins.

**1. Is the outcome a yes/no, a count, or something other than a
continuous measurement?**

| Outcome | Function |
|----|----|
| Two categories (yes/no, success/failure, 0/1) | [`anova_bin()`](https://elkronos.github.io/anovakit/reference/anova_bin.md) |
| A count of events (0, 1, 2, …), possibly over differing exposure time | [`anova_count()`](https://elkronos.github.io/anovakit/reference/anova_count.md) |
| A proportion out of a known number of trials, a positive right-skewed measurement (costs, times), or anything else with a GLM family | [`anova_glm()`](https://elkronos.github.io/anovakit/reference/anova_glm.md) |
| A continuous measurement | go to question 2 |

**2. Are the same subjects measured more than once (before and after, at
several time points, under several conditions)?**

If yes, use
[`anova_rm()`](https://elkronos.github.io/anovakit/reference/anova_rm.md).
Treating repeated measurements as independent observations overstates
the sample size and makes every p-value too small.

**3. Are there several continuous outcomes that should be analysed
together?**

If yes, use
[`anova_manova()`](https://elkronos.github.io/anovakit/reference/anova_manova.md).
It tests the groups on all outcomes jointly and then follows up one
outcome at a time.

**4. Is there a numeric covariate to adjust for (a baseline score, age,
a pre-treatment measurement)?**

If yes, use
[`anova_ancova()`](https://elkronos.github.io/anovakit/reference/anova_ancova.md).
If no, use
[`anova_welch()`](https://elkronos.github.io/anovakit/reference/anova_welch.md),
or
[`anova_kw()`](https://elkronos.github.io/anovakit/reference/anova_kw.md)
when the outcome is ordinal or so skewed that a mean is not a useful
summary.

That covers most analyses. The rest of this guide explains each branch,
with the reasoning behind it.

## One continuous outcome, independent groups

### `anova_welch()`: the default for comparing means

Welch’s test compares group means without assuming the groups have equal
variances. When the variances are equal it gives almost the same answer
as the classical F test; when they are not, the classical test can be
badly wrong, especially when the smallest group has the largest
variance. There is little to lose and a lot to gain, so it is the
default choice.

``` r

fit <- anova_welch(PlantGrowth, "weight", "group", plots = FALSE)
fit$anova
#>    term statistic num_df   den_df    p_value
#> 1 group  5.180972      2 17.12842 0.01739282
```

The [Welch
walkthrough](https://elkronos.github.io/anovakit/articles/welch.md)
covers the standardised mean differences, the pairwise comparisons and
the per-group diagnostics.

### `anova_kw()`: when a mean is the wrong summary

The Kruskal-Wallis test compares the groups’ rank distributions. Use it
for ordinal outcomes (convert an ordered factor with
[`as.integer()`](https://rdrr.io/r/base/integer.html) first), for
heavy-tailed or very skewed data in small samples, or when outliers are
real values you cannot drop. It answers a different question from
Welch’s test — whether values in one group tend to be larger than in
another — so it is not a drop-in replacement when you care about means.

``` r

anova_kw(airquality, "Ozone", "Month", plots = FALSE)$anova
#>    term statistic df      p_value
#> 1 Month  29.26658  4 6.900714e-06
```

See the [Kruskal-Wallis
walkthrough](https://elkronos.github.io/anovakit/articles/kruskal-wallis.md).

### `anova_ancova()`: when a covariate explains part of the outcome

If a numeric variable measured before treatment (a baseline score, age,
weight) predicts the outcome, adjusting for it removes that part of the
variation and makes the group comparison more precise. In a randomised
study, it is also the right way to use a baseline measurement: better
than analysing change scores, and better than ignoring the baseline.

``` r

anc <- anova_ancova(MASS::anorexia, "Postwt", "Treat", "Prewt", plots = FALSE)
anc$anova
#>          term    sum_sq df statistic      p_value
#> 1       Prewt  503.4196  1 11.679513 0.0010866494
#> 2       Treat  785.1051  2  9.107357 0.0003214744
#> 3 Prewt:Treat  466.4783  2  5.411231 0.0066655907
#> 4   Residuals 2844.7843 66        NA           NA
```

The [ANCOVA
walkthrough](https://elkronos.github.io/anovakit/articles/ancova.md)
explains covariate centring, the homogeneity-of-slopes test and when to
fix the model in advance.

### Several grouping variables

Every function accepts several grouping columns.
[`anova_welch()`](https://elkronos.github.io/anovakit/reference/anova_welch.md)
and
[`anova_kw()`](https://elkronos.github.io/anovakit/reference/anova_kw.md)
combine them into one factor of the cells that contain data. The
model-based functions
([`anova_ancova()`](https://elkronos.github.io/anovakit/reference/anova_ancova.md),
[`anova_manova()`](https://elkronos.github.io/anovakit/reference/anova_manova.md),
[`anova_bin()`](https://elkronos.github.io/anovakit/reference/anova_bin.md),
[`anova_count()`](https://elkronos.github.io/anovakit/reference/anova_count.md)
and
[`anova_glm()`](https://elkronos.github.io/anovakit/reference/anova_glm.md))
fit them as factors, with or without their interaction (`interaction`),
and with Type II or Type III sums of squares (`type`). For a classical
two-way ANOVA of a continuous outcome, use
[`anova_glm()`](https://elkronos.github.io/anovakit/reference/anova_glm.md)
with its default Gaussian family:

``` r

tg <- transform(ToothGrowth, dose = factor(dose))
anova_glm(tg, "len", c("supp", "dose"), interaction = TRUE, plots = FALSE)$anova
#>        term   sum_sq df statistic      p_value
#> 1      supp  205.350  1 15.571979 2.311828e-04
#> 2      dose 2426.434  2 91.999965 4.046291e-18
#> 3 supp:dose  108.319  2  4.106991 2.186027e-02
#> 4 Residuals  712.106 54        NA           NA
```

See the [GLM
walkthrough](https://elkronos.github.io/anovakit/articles/glm.md).

## Repeated measurements: `anova_rm()`

When each subject contributes several measurements, the measurements are
correlated.
[`anova_rm()`](https://elkronos.github.io/anovakit/reference/anova_rm.md)
fits the repeated-measures ANOVA through **afex**, with within-subject
factors (`within`), between-subject factors (`between`), the
Greenhouse-Geisser or Huynh-Feldt correction for sphericity, and
generalised eta squared. Subjects missing a cell are removed and listed.
The [repeated-measures
walkthrough](https://elkronos.github.io/anovakit/articles/repeated-measures.md)
shows a mixed design.

A pre/post design with one follow-up measurement is a special case worth
knowing: if the groups were randomised, an ANCOVA of the post score
adjusting for the pre score usually has more power than the
repeated-measures analysis of the same data, and it answers the question
most studies ask (did the groups differ after treatment, given where
they started?).

## Several outcomes: `anova_manova()`

Use a MANOVA when the outcomes are conceptually one construct measured
several ways, or when you want one test that controls the error rate
across them before looking at each. It uses the correlations between the
outcomes, which separate ANOVAs ignore. With covariates it becomes a
MANCOVA.

``` r

anova_manova(iris, c("Sepal.Length", "Sepal.Width", "Petal.Length",
                     "Petal.Width"), "Species", plots = FALSE)$multivariate
#>      term df statistic approx_f num_df den_df      p_value
#> 1 Species  2  1.191899 53.46649      8    290 9.742163e-53
```

The [MANOVA
walkthrough](https://elkronos.github.io/anovakit/articles/manova.md)
covers the four test statistics, the canonical discriminant analysis,
and what to do when Box’s M rejects.

## Outcomes that are not continuous

### `anova_bin()`: binary outcomes

Logistic regression, reported as an analysis of deviance, with odds
ratios, marginal probabilities and pairwise comparisons. It detects
separation (a group in which every outcome is the same), which makes
ordinary intervals meaningless. Aggregated data (one row per group and
outcome, with a count) go in through `weights`. See the [binary-outcomes
walkthrough](https://elkronos.github.io/anovakit/articles/binary.md).

### `anova_count()`: counts and rates

Poisson regression when the counts are no more variable than a Poisson
distribution allows, and a negative binomial model when they are
overdispersed, chosen automatically from the Pearson dispersion unless
you choose yourself with `model`. When units were observed for different
lengths of time or at different sizes, give the exposure as `offset` and
the results become rates. See the [counts
walkthrough](https://elkronos.github.io/anovakit/articles/counts.md).

### `anova_glm()`: everything else

Any [`stats::family`](https://rdrr.io/r/stats/family.html): Gamma or
inverse Gaussian for positive, right-skewed measurements; binomial with
trial weights for proportions; quasi families when the variance does not
follow the family’s rule. With the default Gaussian family it is the
classical ANOVA with Type II or Type III sums of squares. See the [GLM
walkthrough](https://elkronos.github.io/anovakit/articles/glm.md).

## Common situations

| You have | Use | Not |
|----|----|----|
| Three or more groups, variances look different | [`anova_welch()`](https://elkronos.github.io/anovakit/reference/anova_welch.md) | a classical ANOVA, which assumes equal variances |
| Two groups | [`anova_welch()`](https://elkronos.github.io/anovakit/reference/anova_welch.md) (it is Welch’s t test) | a pooled t test |
| Ratings on a 1–5 scale | [`anova_kw()`](https://elkronos.github.io/anovakit/reference/anova_kw.md) | a comparison of means of ordinal codes |
| A baseline and a follow-up score, randomised groups | [`anova_ancova()`](https://elkronos.github.io/anovakit/reference/anova_ancova.md) with the baseline as covariate | change scores |
| The same people at three time points | [`anova_rm()`](https://elkronos.github.io/anovakit/reference/anova_rm.md) | an ANOVA that treats each measurement as a new person |
| Height, weight and waist size on the same people | [`anova_manova()`](https://elkronos.github.io/anovakit/reference/anova_manova.md) | three separate ANOVAs with no joint test |
| Passed / failed | [`anova_bin()`](https://elkronos.github.io/anovakit/reference/anova_bin.md) | a comparison of 0/1 means |
| Number of visits over different lengths of follow-up | [`anova_count()`](https://elkronos.github.io/anovakit/reference/anova_count.md) with `offset` | counts per person ignoring exposure |
| Successes out of a known number of trials | [`anova_glm()`](https://elkronos.github.io/anovakit/reference/anova_glm.md) with `family = "binomial"` and the trials as `weights` | percentages analysed as continuous |
| Costs or durations, positive and skewed | [`anova_glm()`](https://elkronos.github.io/anovakit/reference/anova_glm.md) with `family = Gamma(link = "log")` | a log transformation, when you want means on the original scale |

## Arguments that work the same everywhere

Once you have chosen, the rest is shared. Every function takes the data
first and names columns as character strings. These arguments mean the
same thing wherever they appear:

| Argument | Accepted by | What it does |
|----|----|----|
| `conf_level` | all eight | The level of every interval the function returns. |
| `adjust` | all eight | The multiplicity adjustment for the pairwise comparisons. |
| `posthoc` | all eight | `FALSE` skips the pairwise comparisons. |
| `plots` | all eight | `FALSE` skips building the figures, which is faster. |
| `verbose` | all eight | Progress messages. |
| `interaction` | ancova, manova, bin, count, glm | Additive, full factorial, or up to a given order. |
| `type` | ancova, manova, bin, count, glm | Type II or Type III tests. |
| `vcov_type` | ancova, bin, count, glm | Heteroscedasticity-consistent (sandwich) standard errors. |
| `weights` | bin, count, glm | Prior weights, named as a column. |

The object every function returns is described in [Working with
results](https://elkronos.github.io/anovakit/articles/results.md).
Before you report anything from a fit, read its `$notes`: it is where
the function records what it decided on your behalf and what it could
not compute.
