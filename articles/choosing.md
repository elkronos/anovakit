# Choose an analysis

``` r

library(anovakit)
```

The eight functions share one calling convention and return the same
kind of object, so the choice between them depends only on your data:
what the outcome is, how the observations are related, and what you need
to adjust for. Each function below links to a walkthrough that runs it
on real data.

## Four questions

Answer these in order and stop at the first one that settles it.

1.  Is the outcome binary, a count, or another non-continuous
    measurement?

    | Outcome | Function |
    |----|----|
    | Two categories (yes/no, success/failure, 0/1) | [`anova_bin()`](https://elkronos.github.io/anovakit/reference/anova_bin.md) |
    | A count of events (0, 1, 2, …), possibly over differing exposure time | [`anova_count()`](https://elkronos.github.io/anovakit/reference/anova_count.md) |
    | A proportion out of a known number of trials, a positive right-skewed measurement (costs, times), or anything else with a GLM family | [`anova_glm()`](https://elkronos.github.io/anovakit/reference/anova_glm.md) |
    | A continuous measurement | go to question 2 |

2.  Are the same subjects measured more than once (before and after, at
    several time points, under several conditions)? If yes, use
    [`anova_rm()`](https://elkronos.github.io/anovakit/reference/anova_rm.md).

3.  Are there several continuous outcomes to analyse together? If yes,
    use
    [`anova_manova()`](https://elkronos.github.io/anovakit/reference/anova_manova.md).

4.  Is there a numeric covariate to adjust for (a baseline score, age, a
    pre-treatment measurement)? If yes, use
    [`anova_ancova()`](https://elkronos.github.io/anovakit/reference/anova_ancova.md).
    If no, and the groups come from one grouping variable, use
    [`anova_welch()`](https://elkronos.github.io/anovakit/reference/anova_welch.md),
    or
    [`anova_kw()`](https://elkronos.github.io/anovakit/reference/anova_kw.md)
    when the outcome is ordinal or so skewed that a mean is not a useful
    summary. With two or more crossed factors whose main effects and
    interaction you want to test, use
    [`anova_glm()`](https://elkronos.github.io/anovakit/reference/anova_glm.md)
    with its default Gaussian family and `interaction = TRUE`, which
    gives the classical factorial ANOVA (see [Several grouping
    variables](#several-grouping-variables)).

## Common situations

| You have | Use | Not |
|----|----|----|
| Three or more groups, variances look different | [`anova_welch()`](https://elkronos.github.io/anovakit/reference/anova_welch.md) | a classical ANOVA, which assumes equal variances |
| Two groups | [`anova_welch()`](https://elkronos.github.io/anovakit/reference/anova_welch.md) (it is Welch’s t test) | a pooled t test |
| Ratings on a 1–5 scale | [`anova_kw()`](https://elkronos.github.io/anovakit/reference/anova_kw.md) | a comparison of means of ordinal codes |
| A baseline and a follow-up score, randomised groups | [`anova_ancova()`](https://elkronos.github.io/anovakit/reference/anova_ancova.md) with the baseline as covariate | change scores |
| The same people at three time points | [`anova_rm()`](https://elkronos.github.io/anovakit/reference/anova_rm.md) | an ANOVA that treats each measurement as a new person |
| Height, weight and waist size on the same people | [`anova_manova()`](https://elkronos.github.io/anovakit/reference/anova_manova.md) | three separate ANOVAs with no joint test |
| Two crossed factors and a continuous outcome | [`anova_glm()`](https://elkronos.github.io/anovakit/reference/anova_glm.md) with `interaction = TRUE` | [`anova_welch()`](https://elkronos.github.io/anovakit/reference/anova_welch.md), which cannot separate main effects from the interaction |
| Passed / failed | [`anova_bin()`](https://elkronos.github.io/anovakit/reference/anova_bin.md) | a comparison of 0/1 means |
| Number of visits over different lengths of follow-up | [`anova_count()`](https://elkronos.github.io/anovakit/reference/anova_count.md) with `offset` | counts per person ignoring exposure |
| Successes out of a known number of trials | [`anova_glm()`](https://elkronos.github.io/anovakit/reference/anova_glm.md) with `family = "binomial"` and the trials as `weights` | percentages analysed as continuous |
| Costs or durations, positive and skewed | [`anova_glm()`](https://elkronos.github.io/anovakit/reference/anova_glm.md) with `family = Gamma(link = "log")` | a log transformation, when you want means on the original scale |

## One outcome, independent groups

### `anova_welch()`

Welch’s test compares group means without assuming equal variances. When
the variances are equal it gives almost the same answer as the classical
F test. When they are not, the classical test can be badly wrong, most
of all when the smallest group has the largest variance. That makes it
the default for comparing means. The
[walkthrough](https://elkronos.github.io/anovakit/articles/welch.md)
runs it on data with unequal group variances.

### `anova_kw()`

The Kruskal-Wallis test compares the groups’ rank distributions, with
Dunn’s test for the pairwise comparisons. Use it for ordinal outcomes
(convert an ordered factor with
[`as.integer()`](https://rdrr.io/r/base/integer.html) first), for
heavy-tailed or very skewed data in small samples, or when outliers are
real values you cannot drop. It asks whether values in one group tend to
be larger than in another, which is a different question from whether
the means differ, and the [Kruskal-Wallis
walkthrough](https://elkronos.github.io/anovakit/articles/kruskal-wallis.md)
shows how to read it.

### `anova_ancova()`

When a numeric variable measured before treatment (a baseline score,
age, weight) predicts the outcome, adjusting for it removes that part of
the variation and makes the group comparison more precise. In a
randomised study it is the right way to use a baseline measurement:
better than analysing change scores, and better than ignoring the
baseline. The [ANCOVA
walkthrough](https://elkronos.github.io/anovakit/articles/ancova.md)
covers covariate centring, the homogeneity-of-slopes test and when to
fix the model in advance.

## Repeated measures and several outcomes

### `anova_rm()`

When each subject contributes several measurements, the measurements are
correlated. Treating them as independent observations overstates the
sample size and makes every p-value too small.
[`anova_rm()`](https://elkronos.github.io/anovakit/reference/anova_rm.md)
fits the repeated-measures ANOVA through `afex`, with within-subject
factors (`within`), between-subject factors (`between`), the
Greenhouse-Geisser or Huynh-Feldt correction for sphericity, and
generalised eta squared. Subjects missing a cell are removed and listed
in `$subjects_dropped`. The [repeated-measures
walkthrough](https://elkronos.github.io/anovakit/articles/repeated-measures.md)
runs a mixed design.

With one pre-treatment and one post-treatment measurement and randomised
groups, an ANCOVA of the post score adjusting for the pre score usually
has more power than the repeated-measures analysis of the same data. It
also answers the question most such studies ask: did the groups differ
after treatment, given where they started?

### `anova_manova()`

Use a MANOVA when the outcomes are one construct measured several ways,
or when you want one test that controls the error rate across them
before looking at each. It uses the correlations between the outcomes,
which separate ANOVAs ignore, and then follows up each outcome on its
own. With covariates it becomes a MANCOVA. The [MANOVA
walkthrough](https://elkronos.github.io/anovakit/articles/manova.md)
covers the four test statistics, the canonical discriminant analysis,
and what to do when Box’s M rejects.

## Binary, count and other outcomes

### `anova_bin()`

[`anova_bin()`](https://elkronos.github.io/anovakit/reference/anova_bin.md)
fits a logistic regression and reports it as an analysis of deviance,
with odds ratios, marginal probabilities and pairwise comparisons. It
detects separation (a group in which every outcome is the same), which
makes ordinary intervals meaningless, and says so in `$notes`.
Aggregated data, one row per group and outcome with a count, go in
through `weights`, as the [binary-outcomes
walkthrough](https://elkronos.github.io/anovakit/articles/binary.md)
shows.

### `anova_count()`

[`anova_count()`](https://elkronos.github.io/anovakit/reference/anova_count.md)
fits a Poisson model when the counts are no more variable than a Poisson
distribution allows, and a negative binomial model when the Pearson
dispersion is above `overdispersion_threshold` (1.5 by default). Set
`model` to choose the family yourself. When units were observed for
different lengths of time or at different sizes, give the exposure as
`offset` and the results become rates. The [counts
walkthrough](https://elkronos.github.io/anovakit/articles/counts.md)
models insurance claims per policyholder this way.

### `anova_glm()`

Any [`stats::family`](https://rdrr.io/r/stats/family.html): Gamma or
inverse Gaussian for positive, right-skewed measurements; binomial with
trial weights for proportions; quasi families when the variance does not
follow the family’s rule. With the default Gaussian family it is the
classical ANOVA with Type II or Type III sums of squares. The [GLM
walkthrough](https://elkronos.github.io/anovakit/articles/glm.md) fits a
Gamma model to car prices.

## Several grouping variables

Every function accepts several grouping columns.
[`anova_welch()`](https://elkronos.github.io/anovakit/reference/anova_welch.md)
and
[`anova_kw()`](https://elkronos.github.io/anovakit/reference/anova_kw.md)
combine them into one factor of the cells that contain data, so they
test the cells and cannot separate main effects from an interaction.
[`anova_ancova()`](https://elkronos.github.io/anovakit/reference/anova_ancova.md),
[`anova_manova()`](https://elkronos.github.io/anovakit/reference/anova_manova.md),
[`anova_bin()`](https://elkronos.github.io/anovakit/reference/anova_bin.md),
[`anova_count()`](https://elkronos.github.io/anovakit/reference/anova_count.md)
and
[`anova_glm()`](https://elkronos.github.io/anovakit/reference/anova_glm.md)
fit them as factors, with or without their interaction (`interaction`),
and with Type II or Type III tests (`type`). For a classical two-way
ANOVA of a continuous outcome, use
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

The [GLM
walkthrough](https://elkronos.github.io/anovakit/articles/glm.md) runs a
two-way ANOVA and explains the choice between Type II and Type III.

## See also

- [Get
  started](https://elkronos.github.io/anovakit/articles/anovakit.md):
  one analysis from start to finish, and the arguments every function
  shares (`conf_level`, `adjust`, `posthoc`, `plots`, `interaction`,
  `type`, `vcov_type`, `weights`) with their defaults.
- [Weights and robust standard
  errors](https://elkronos.github.io/anovakit/articles/weights.md):
  frequency, trial and sampling weights, and `vcov_type`.
- [Working with
  results](https://elkronos.github.io/anovakit/articles/results.md):
  what the returned object holds.
- [Diagnostics](https://elkronos.github.io/anovakit/articles/diagnostics.md):
  the assumption checks each function runs.
