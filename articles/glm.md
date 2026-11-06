# Any GLM family: analysis of deviance with anova_glm()

``` r

library(anovakit)
```

[`anova_glm()`](https://elkronos.github.io/anovakit/reference/anova_glm.md)
fits a generalised linear model of a response on one or more grouping
factors, in any `family`, and reports an analysis of deviance table,
effect sizes, marginal means and pairwise comparisons. With the default
Gaussian family it is an ordinary factorial analysis of variance. Other
families model responses whose spread changes with their mean.

Four functions cover particular cases better. For a one-way comparison
of means with unequal variances,
[`anova_welch()`](https://elkronos.github.io/anovakit/reference/anova_welch.md)
is better calibrated
([walkthrough](https://elkronos.github.io/anovakit/articles/welch.md)).
For 0/1 outcomes use
[`anova_bin()`](https://elkronos.github.io/anovakit/reference/anova_bin.md)
([walkthrough](https://elkronos.github.io/anovakit/articles/binary.md)),
and for counts
[`anova_count()`](https://elkronos.github.io/anovakit/reference/anova_count.md)
([walkthrough](https://elkronos.github.io/anovakit/articles/counts.md)):
both add diagnostics and effect sizes that
[`anova_glm()`](https://elkronos.github.io/anovakit/reference/anova_glm.md)
does not. To adjust for a numeric covariate use
[`anova_ancova()`](https://elkronos.github.io/anovakit/reference/anova_ancova.md)
([walkthrough](https://elkronos.github.io/anovakit/articles/ancova.md)).

This walkthrough follows two analyses side by side: a Gaussian two-way
ANOVA, and a Gamma model for prices.

## The data

`ToothGrowth` records the length of odontoblasts (the cells responsible
for tooth growth) in 60 guinea pigs. Each animal received vitamin C at
one of three doses (0.5, 1 or 2 mg/day) by one of two delivery methods:
orange juice (`OJ`) or ascorbic acid (`VC`). `dose` is stored as a
number; here it is a grouping variable with three levels, so make it a
factor.

``` r

tg <- ToothGrowth
tg$dose <- factor(tg$dose)
table(tg$supp, tg$dose)
#>     
#>      0.5  1  2
#>   OJ  10 10 10
#>   VC  10 10 10
round(tapply(tg$len, list(tg$supp, tg$dose), mean), 2)
#>      0.5     1     2
#> OJ 13.23 22.70 26.06
#> VC  7.98 16.77 26.14
```

This is a balanced 2 x 3 design with ten animals per cell. Besides
whether dose and supplement matter, the question is whether the effect
of the supplement depends on the dose. The cell means suggest it does:
orange juice is ahead at the two lower doses, but not at 2 mg/day.

[`MASS::Cars93`](https://rdrr.io/pkg/MASS/man/Cars93.html) lists 93 car
models on sale in the US in 1993. Their `Price` (in thousands of
dollars) is positive and right-skewed, and within each `Type` of car the
spread grows with the typical price:

``` r

cars <- MASS::Cars93[, c("Manufacturer", "Model", "Type", "Price")]
price <- split(cars$Price, cars$Type)
cs <- t(sapply(price, function(x) c(n = length(x), median = median(x),
                                    mean = mean(x), sd = sd(x),
                                    cv = sd(x) / mean(x))))
round(cs, 2)
#>          n median  mean    sd   cv
#> Compact 16  16.15 18.21  6.69 0.37
#> Large   11  20.90 24.30  6.34 0.26
#> Midsize 22  26.20 27.22 12.26 0.45
#> Small   21  10.00 10.17  1.95 0.19
#> Sporty  14  16.80 19.39  7.97 0.41
#> Van      9  19.10 19.10  1.88 0.10
```

A Gaussian model assumes the same standard deviation in every group, and
here it ranges from 1.88 for vans to 12.26 for midsize cars. A Gamma
family assumes instead that the standard deviation is proportional to
the mean, a constant coefficient of variation (`cv`). The `cv` column
fits that better, though not perfectly: vans are unusually uniform.

## Fitting the model

Columns are named as character strings. With `interaction = TRUE` the
model contains both factors and their interaction. Listing `dose` first
puts it on the x axis of the marginal-means plot.

``` r

fit <- anova_glm(tg, "len", c("dose", "supp"), interaction = TRUE)
fit
#> Analysis of deviance (gaussian family, identity link, Type II) 
#> --------------------------------------------------------------
#> Call: anova_glm(data = tg, response = "len", groups = c("dose", "supp"), interaction = TRUE)
#> Observations used: 60
#> 
#> Omnibus test (Type II, F)
#>        term sum_sq df statistic   p_value
#> 1      dose   2426  2        92 < 2.2e-16
#> 2      supp  205.4  1     15.57 0.0002312
#> 3 dose:supp  108.3  2     4.107   0.02186
#> 4 Residuals  712.1 54        NA        NA
#> 
#> Plots available: box, residuals, qq, emmeans
#>   (use plot(x, which = "box"))
```

The first line names the family, the link and the type of sums of
squares; a Gaussian family with the identity link is an ordinary linear
model. `Observations used: 60` says no rows were dropped. The heading
`Omnibus test (Type II, F)` names the type of test and the statistic,
which the table itself does not show. There is no `Notes` block because
nothing needed saying, and `Plots available` lists the four plots in
`$plots`. `summary(fit)` adds the assumption checks, the effect sizes,
the marginal means and the pairwise comparisons.

For the car prices, a Gamma family with a log link makes the group
effects multiplicative, so comparisons come out as ratios of means:

``` r

gam <- anova_glm(cars, "Price", "Type", family = Gamma(link = "log"))
gam
#> Analysis of deviance (Gamma family, log link, Type II) 
#> ------------------------------------------------------
#> Call: anova_glm(data = cars, response = "Price", groups = "Type", family = Gamma(link = "log"))
#> Observations used: 93
#> 
#> Omnibus test (Type II, F)
#>        term sum_sq df statistic   p_value
#> 1      Type  10.65  5     18.56 1.663e-12
#> 2 Residuals   9.98 87        NA        NA
#> 
#> Notes
#>   - McFadden's pseudo R squared is not reported (NA): the log-likelihood of the Gamma family is a log density, which shifts with the units of the response, so the ratio has no unit-free meaning. The deviance explained does not depend on the units and is reported instead.
#> 
#> Plots available: box, residuals, qq, emmeans
#>   (use plot(x, which = "box"))
```

`family` accepts a family object, as here, a family function, or a name
such as `"Gamma"`. The name alone gives the Gamma family’s default
inverse link, which is rarely what you want for group comparisons (see
[below](#the-link-decides-what-a-comparison-is)). The note explains why
one effect size is `NA`; [Effect sizes](#effect-sizes) has the details.

## The omnibus test

``` r

fit$anova
#>        term   sum_sq df statistic      p_value
#> 1      dose 2426.434  2 91.999965 4.046291e-18
#> 2      supp  205.350  1 15.571979 2.311828e-04
#> 3 dose:supp  108.319  2  4.106991 2.186027e-02
#> 4 Residuals  712.106 54        NA           NA
```

`term` names the effect, and `dose:supp` is the interaction. `sum_sq`
and `df` are each term’s sum of squares and degrees of freedom; the
`Residuals` row holds the error sum of squares. `statistic` is the F
statistic, the term’s mean square over the residual mean square, and
`p_value` its p-value.

The interaction is significant (F = 4.11, p = 0.022): the supplement
effect is not the same at every dose. A main-effect row averages over a
pattern that changes, so read the cell means and comparisons below
rather than the main effects alone.

### Type II or Type III

With `type = "II"` (the default) each main effect is tested after the
other main effect, but not after the interaction. With `type = "III"`
every term is tested after every other, interaction included. In a
balanced design like this one the two coincide:

``` r

fit3 <- anova_glm(tg, "len", c("dose", "supp"), interaction = TRUE,
                  type = "III", plots = FALSE)
fit3$anova
#>        term   sum_sq df statistic      p_value
#> 1      dose 2426.434  2 91.999965 4.046291e-18
#> 2      supp  205.350  1 15.571979 2.311828e-04
#> 3 dose:supp  108.319  2  4.106991 2.186027e-02
#> 4 Residuals  712.106 54        NA           NA
```

They differ when the cells are unequal. Suppose four animals had been
lost from the `VC`, 2 mg cell and two from the `OJ`, 0.5 mg cell:

``` r

lost <- c(which(tg$supp == "VC" & tg$dose == "2")[1:4],
          which(tg$supp == "OJ" & tg$dose == "0.5")[1:2])
tg_u <- tg[-lost, ]
table(tg_u$supp, tg_u$dose)
#>     
#>      0.5  1  2
#>   OJ   8 10 10
#>   VC  10 10  6
anova_glm(tg_u, "len", c("dose", "supp"), interaction = TRUE,
          plots = FALSE)$anova
#>        term     sum_sq df statistic      p_value
#> 1      dose 2173.01135  2 104.36382 3.330016e-18
#> 2      supp  154.03314  1  14.79559 3.529409e-04
#> 3 dose:supp   93.14518  2   4.47351 1.653811e-02
#> 4 Residuals  499.71600 48        NA           NA
anova_glm(tg_u, "len", c("dose", "supp"), interaction = TRUE,
          type = "III", plots = FALSE)$anova
#>        term     sum_sq df statistic      p_value
#> 1      dose 2249.45726  2 108.03531 1.692417e-18
#> 2      supp  125.31484  1  12.03706 1.111881e-03
#> 3 dose:supp   93.14518  2   4.47351 1.653811e-02
#> 4 Residuals  499.71600 48        NA           NA
```

The interaction row is the same in both, because it is adjusted for
everything either way; the main-effect rows differ. Type II has more
power to detect a main effect when there is no interaction. Type III
tests each main effect averaged over the levels of the other factor,
weighting each cell equally, which some analysts prefer when an
interaction is present because it does not depend on the cell sizes.

### Type III is fitted under sum-to-zero contrasts

A Type III test of `supp` removes the columns that code `supp` while
keeping the interaction, and what that tests depends on how the factors
are coded. Under R’s default treatment contrasts, the `supp` columns
describe the supplement effect at the reference dose only, so the “Type
III” row for `supp` tests the supplement at 0.5 mg/day. Under
sum-to-zero contrasts it tests the supplement effect averaged over the
doses, which is what a main effect means. So with `type = "III"`,
[`anova_glm()`](https://elkronos.github.io/anovakit/reference/anova_glm.md)
fits the model under `contr.sum`. Setting `options(contrasts = )` after
fitting would not help, because the contrasts are stored on the fitted
model. Here is what fitting with the defaults and asking
[`car::Anova()`](https://rdrr.io/pkg/car/man/Anova.html) for Type III
gives:

``` r

naive <- glm(len ~ dose * supp, data = tg)    # treatment contrasts
naive_tab <- car::Anova(naive, type = 3, test.statistic = "F")
naive_tab
#> Analysis of Deviance Table (Type III tests)
#> 
#> Response: len
#> Error estimate based on Pearson residuals 
#> 
#>           Sum Sq Df F values    Pr(>F)    
#> dose      885.26  2   33.565 3.363e-10 ***
#> supp      137.81  1   10.450  0.002092 ** 
#> dose:supp 108.32  2    4.107  0.021860 *  
#> Residuals 712.11 54                       
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
unlist(fit3$model$contrasts)
#>        dose        supp 
#> "contr.sum" "contr.sum"
```

The `supp` row has F = 10.45, the test of the supplement at the
reference dose, not the F = 15.57 of the averaged main effect. The
contrasts stored on `fit3$model` (`contr.sum`) are what make its Type
III row a test of the averaged effect.

### Choosing the test statistic

`test_statistic = NULL` (the default) chooses the test after fitting. A
family that estimates its dispersion, like the Gaussian and the Gamma,
gets an F test, and a family whose dispersion is fixed at 1 (binomial,
Poisson) gets a likelihood-ratio test. In the F table of a non-Gaussian
family the `sum_sq` column holds deviances: the `Type` row is the drop
in deviance when `Type` is added, and the `Residuals` row is the Pearson
chi-squared statistic, which divided by its degrees of freedom gives the
dispersion. F is the ratio of the two per degree of freedom:

``` r

a <- gam$anova
a
#>        term    sum_sq df statistic      p_value
#> 1      Type 10.645956  5  18.56058 1.662501e-12
#> 2 Residuals  9.980272 87        NA           NA
(a$sum_sq[1] / a$df[1]) / (a$sum_sq[2] / a$df[2])
#> [1] 18.56058
```

`"LR"` and `"Wald"` give the other two statistics. The columns come in
the same order whatever the statistic, and only an F table has `sum_sq`
and a `Residuals` row:

``` r

for (ts in c("LR", "Wald")) {
  g <- anova_glm(cars, "Price", "Type", family = Gamma(link = "log"),
                 test_statistic = ts, plots = FALSE)
  cat("\n", attr(g$anova, "statistic"), "\n", sep = "")
  print(g$anova)
}
#> 
#> scaled-deviance chi-square, estimated dispersion treated as known
#>   term df statistic      p_value
#> 1 Type  5   92.8029 1.730729e-18
#> 
#> Wald chi-square
#>   term df statistic      p_value
#> 1 Type  5  101.7011 2.314394e-20
```

All three agree that the car types differ. For a family that estimates
its dispersion, the “LR” statistic is not a true likelihood ratio. It is
the deviance drop divided by the estimated dispersion, referred to a
chi-squared distribution as though the dispersion were known, and the
heading calls it a scaled-deviance chi-square. It ignores the
uncertainty in the dispersion, which the F test allows for, so F is the
better default there. A quasi family has no likelihood, so
`test_statistic = "LR"` gives the F test, with a note. The Wald test is
built from the coefficients and their covariance, and is the only one of
the three that can use robust standard errors (see [Options](#options)).

## Effect sizes

For a Gaussian model `$effect_sizes` holds partial eta squared and
partial omega squared for each term, from the sums of squares of the F
table:

``` r

fit$effect_sizes
#>        term df   sum_sq partial_eta_sq partial_omega_sq
#> 1      dose  2 2426.434      0.7731092       0.75206604
#> 2      supp  1  205.350      0.2238254       0.19540824
#> 3 dose:supp  2  108.319      0.1320279       0.09384698
```

`partial_eta_sq` is the term’s sum of squares over itself plus the
residual sum of squares: the share of the variance not explained by the
other terms that this term explains. Dose accounts for 77% of it.
`partial_omega_sq` corrects eta squared for its upward bias in samples,
so it is a little smaller; a negative value is floored at 0, with a
note. Both are partial, so they do not add up to 1 across terms and are
not a partition of the total variance.

Any other family has no sums of squares, so the effect size is
`deviance_explained`, 1 minus the residual deviance over the null
deviance:

``` r

gam$effect_sizes
#>              measure estimate
#> 1 deviance_explained 0.542735
#> 2        mcfadden_r2       NA
```

Car type explains 54% of the deviance. McFadden’s pseudo R squared
(`mcfadden_r2`) is reported only for the binomial, Poisson and negative
binomial families. For a Gamma model it is `NA`, because the
log-likelihood of a density shifts with the units of the response, and
the note says so.

## Marginal means and comparisons

Because the ToothGrowth model contains the interaction, the marginal
means are for the six cells. For a Gaussian model with the identity link
they are the cell means, with intervals that use the pooled residual
variance:

``` r

fit$emmeans
#>   dose supp estimate       se df  conf_low conf_high
#> 1  0.5   OJ    13.23 1.148353 54 10.927691  15.53231
#> 2    1   OJ    22.70 1.148353 54 20.397691  25.00231
#> 3    2   OJ    26.06 1.148353 54 23.757691  28.36231
#> 4  0.5   VC     7.98 1.148353 54  5.677691  10.28231
#> 5    1   VC    16.77 1.148353 54 14.467691  19.07231
#> 6    2   VC    26.14 1.148353 54 23.837691  28.44231
```

`$posthoc` compares every pair of cells, 15 in all, with a Tukey
adjustment by default. The three comparisons of the supplements at the
same dose are:

``` r

same_dose <- c("dose0.5 OJ - dose0.5 VC", "dose1 OJ - dose1 VC",
               "dose2 OJ - dose2 VC")
ph <- fit$posthoc
print(ph[ph$contrast %in% same_dose, ], digits = 3)
#>                   contrast estimate   se df conf_low conf_high statistic
#> 3  dose0.5 OJ - dose0.5 VC     5.25 1.62 54    0.452     10.05    3.2327
#> 8      dose1 OJ - dose1 VC     5.93 1.62 54    1.132     10.73    3.6514
#> 12     dose2 OJ - dose2 VC    -0.08 1.62 54   -4.878      4.72   -0.0493
#>    p_value p_adjusted adjustment
#> 3  0.00209    0.02425      tukey
#> 8  0.00059    0.00739      tukey
#> 12 0.96089    1.00000      tukey
```

`estimate` is the first cell minus the second. `conf_low` and
`conf_high` are Tukey-adjusted (simultaneous) intervals. `p_value` is
unadjusted, and `p_adjusted` is adjusted across all 15 comparisons.
Orange juice gives longer odontoblasts than ascorbic acid at 0.5 and 1
mg/day, and the two are indistinguishable at 2 mg/day.

If only those three comparisons interest you, adjusting for all 15 is
more conservative than it needs to be. `$emmeans_object` is the emmeans
grid, so you can ask for exactly the family of comparisons you want. It
carries emmeans’ own default level, so give `level =` explicitly:

``` r

by_dose <- emmeans::contrast(fit$emmeans_object, method = "pairwise",
                             by = "dose")
summary(by_dose, infer = TRUE, level = fit$conf_level)
#> dose = 0.5:
#>  contrast estimate   SE df lower.CL upper.CL t.ratio p.value
#>  OJ - VC      5.25 1.62 54     1.99     8.51   3.233  0.0021
#> 
#> dose = 1:
#>  contrast estimate   SE df lower.CL upper.CL t.ratio p.value
#>  OJ - VC      5.93 1.62 54     2.67     9.19   3.651  0.0006
#> 
#> dose = 2:
#>  contrast estimate   SE df lower.CL upper.CL t.ratio p.value
#>  OJ - VC     -0.08 1.62 54    -3.34     3.18  -0.049  0.9609
#> 
#> Degrees-of-freedom method: user-specified 
#> Confidence level used: 0.95
```

With one comparison per dose there is nothing to adjust within each
family. [Working with
results](https://elkronos.github.io/anovakit/articles/results.md) covers
`$emmeans_object` further.

### Ratios on the log scale

The Gamma model’s marginal means are on the response scale, in thousands
of dollars. In a one-factor model they are the group means. Their
intervals are computed on the log scale and back-transformed, so they
are not symmetric about the estimate:

``` r

gam$emmeans
#>      Type estimate        se df  conf_low conf_high
#> 1 Compact 18.21250 1.5421305 87 15.391399  21.55068
#> 2   Large 24.30000 2.4815412 87 19.836025  29.76856
#> 3 Midsize 27.21818 1.9654379 87 23.579061  31.41895
#> 4   Small 10.16667 0.7514161 87  8.777667  11.77546
#> 5  Sporty 19.39286 1.7554535 87 16.199579  23.21560
#> 6     Van 19.10000 2.1563719 87 15.260827  23.90499
```

With a log link a pairwise comparison is a ratio of means, in a column
called `ratio`, tested against 1. These are the comparisons that survive
Tukey’s adjustment:

``` r

gph <- gam$posthoc
print(gph[gph$p_adjusted < 0.05,
          c("contrast", "ratio", "conf_low", "conf_high", "p_adjusted")],
      digits = 3)
#>             contrast ratio conf_low conf_high p_adjusted
#> 2  Compact / Midsize 0.669    0.484     0.925   6.57e-03
#> 3    Compact / Small 1.791    1.291     2.486   2.00e-05
#> 7      Large / Small 2.390    1.655     3.451   1.16e-08
#> 10   Midsize / Small 2.677    1.981     3.618   3.09e-10
#> 11  Midsize / Sporty 1.404    1.002     1.967   4.83e-02
#> 13    Small / Sporty 0.524    0.373     0.737   4.94e-06
#> 14       Small / Van 0.532    0.359     0.789   1.53e-04
```

Midsize cars, for example, average 2.68 times the price of small ones.

### The link decides what a comparison is

With the log link comparisons are ratios, and with the logit link they
are odds ratios. With any other non-identity link they are differences
between the means on the response scale, and `$notes` says so. The Gamma
family’s default link is the inverse:

``` r

inv <- anova_glm(cars, "Price", "Type", family = "Gamma", plots = FALSE)
inv$method
#> [1] "Analysis of deviance (Gamma family, inverse link, Type II)"
inv$notes[grepl("link", inv$notes)]
#> [1] "With the inverse link the pairwise comparisons are differences between the marginal means on the response scale; only log and logit links give ratios."
head(inv$posthoc[, c("contrast", "estimate", "conf_low", "conf_high")], 3)
#>            contrast  estimate   conf_low conf_high
#> 1   Compact - Large -6.087500 -14.601888  2.426888
#> 2 Compact - Midsize -9.005682 -16.286022 -1.725342
#> 3   Compact - Small  8.045833   3.046632 13.045035
```

### The coefficient table

The fitted coefficients are in `$assumptions$coefficients`:

``` r

fit$assumptions$coefficients
#>           term estimate       se  statistic      p_value distribution
#> 1  (Intercept)    13.23 1.148353 11.5208468 3.602548e-16            t
#> 2        dose1     9.47 1.624017  5.8312215 3.175641e-07            t
#> 3        dose2    12.83 1.624017  7.9001660 1.429712e-10            t
#> 4       suppVC    -5.25 1.624017 -3.2327258 2.092470e-03            t
#> 5 dose1:suppVC    -0.68 2.296706 -0.2960762 7.683076e-01            t
#> 6 dose2:suppVC     5.33 2.296706  2.3207148 2.410826e-02            t
#>     conf_low conf_high ci_method
#> 1 10.9276907 15.532309   profile
#> 2  6.2140429 12.725957   profile
#> 3  9.5740429 16.085957   profile
#> 4 -8.5059571 -1.994043   profile
#> 5 -5.2846186  3.924619   profile
#> 6  0.7253814  9.934619   profile
```

The table has t statistics (the dispersion is estimated) and
profile-likelihood intervals by default, which for a Gaussian model
equal `confint(lm())`. Under the default treatment contrasts the
intercept is the mean of the reference cell (0.5 mg/day, orange juice),
and the other coefficients are differences from it. `fit3`, fitted under
sum-to-zero contrasts, has different coefficients for the same fitted
values: the intercept is the average of the six cell means, `dose1` is
the average for 0.5 mg/day minus that grand mean, and `supp1` the
average for orange juice minus it.

``` r

fit3$assumptions$coefficients[, c("term", "estimate")]
#>          term   estimate
#> 1 (Intercept) 18.8133333
#> 2       dose1 -8.2083333
#> 3       dose2  0.9216667
#> 4       supp1  1.8500000
#> 5 dose1:supp1  0.7750000
#> 6 dose2:supp1  1.1150000
mean(fit$emmeans$estimate)
#> [1] 18.81333
```

The Gamma coefficients are on the log scale:

``` r

co <- gam$assumptions$coefficients
co
#>          term    estimate         se  statistic      p_value distribution
#> 1 (Intercept)  2.90210817 0.08467429 34.2737820 2.652780e-52            t
#> 2   TypeLarge  0.28836818 0.13265912  2.1737532 3.243991e-02            t
#> 3 TypeMidsize  0.40177703 0.11128382  3.6103814 5.103876e-04            t
#> 4   TypeSmall -0.58299378 0.11239391 -5.1870584 1.377488e-06            t
#> 5  TypeSporty  0.06279664 0.12395047  0.5066269 6.136982e-01            t
#> 6     TypeVan  0.04758016 0.14112382  0.3371519 7.368144e-01            t
#>     conf_low  conf_high ci_method
#> 1  2.7384042  3.0752703   profile
#> 2  0.0264585  0.5545836   profile
#> 3  0.1790690  0.6219030   profile
#> 4 -0.8077507 -0.3604911   profile
#> 5 -0.1832054  0.3101515   profile
#> 6 -0.2297083  0.3322378   profile
```

The intercept is the log of the mean price of the reference type,
`Compact`, and each other coefficient is the log of a ratio to it: small
cars cost exp(-0.583) = 0.56 times as much as compact ones on average.
The `distribution` column says the tests use t, because the dispersion
is estimated.

## Checking assumptions

A Gaussian model assumes independent observations, normally distributed
errors and the same error variance in every cell. `$assumptions` holds
the dispersion, and the plots in the next section cover normality.

### Dispersion

``` r

fit$assumptions$dispersion
#> [1] 13.18715
sqrt(fit$assumptions$dispersion)
#> [1] 3.631411
```

`dispersion` is the Pearson dispersion: the sum of squared Pearson
residuals divided by the residual degrees of freedom. For a Gaussian
model that is the residual mean square, so its square root, 3.63, is the
residual standard deviation in the units of the response. There is no
formal test of equal variances; compare the cell standard deviations:

``` r

cell_sd <- tapply(tg$len, list(tg$supp, tg$dose), sd)
round(cell_sd, 2)
#>     0.5    1    2
#> OJ 4.46 3.91 2.66
#> VC 2.75 2.52 4.80
```

They range from 2.52 to 4.80, so the largest cell variance is 3.64 times
the smallest. With ten animals in every cell, a balanced design
tolerates that well. If the variances differed badly, use robust
standard errors (`vcov_type`, under [Options](#options)) or, for a
one-way design,
[`anova_welch()`](https://elkronos.github.io/anovakit/reference/anova_welch.md).

For a Gamma family the dispersion is the squared coefficient of
variation:

``` r

gam$assumptions$dispersion
#> [1] 0.1147158
sqrt(gam$assumptions$dispersion)
#> [1] 0.3386972
```

The model puts the standard deviation of prices at about 34% of the mean
in every car type. The observed values in [the data](#the-data) range
from 10% (vans) to 45% (midsize cars), so the constant-CV assumption is
only roughly true; the `vcov_type` option below protects the comparisons
against that.

### Is Gamma the right family?

The vans show that the Gamma variance function is not exact. To check
that it still describes the prices better than the Gaussian, compare the
two fits by AIC. `$model_stats` holds the AIC, BIC, null and residual
deviances and residual degrees of freedom, and `$family` is the family
object used. The Gaussian and Gamma fits are full likelihoods for the
same response, so their AICs can be compared directly:

``` r

str(gam$model_stats)
#> List of 5
#>  $ AIC              : num 595
#>  $ BIC              : num 613
#>  $ null_deviance    : num 19.6
#>  $ residual_deviance: num 8.97
#>  $ df_residual      : int 87
gau <- anova_glm(cars, "Price", "Type", plots = FALSE)
c(gaussian = gau$model_stats$AIC, gamma = gam$model_stats$AIC)
#> gaussian    gamma 
#> 651.4658 594.8847
```

The Gamma model has the lower AIC by 57. The [diagnostics
guide](https://elkronos.github.io/anovakit/articles/diagnostics.md)
covers the checks across the package.

## Plots

`fit` has four plots. The `box` plot is the raw data behind the means:
the observed lengths in each cell, ordered by median and labelled with
the cell size.

``` r

plot(fit, "box")
```

![](glm_files/figure-html/plot-box-1.png)

In the `emmeans` plot the gap between the supplements closes at 2
mg/day. Dose, the first grouping factor, is on the x axis, and the
supplement sets the colour, with one line per supplement, so lines that
are not parallel are the interaction. The error bars are intervals for
each mean at `conf_level`.

``` r

plot(fit, "emmeans")
```

![](glm_files/figure-html/plot-emmeans-1.png)

The `residuals` and `qq` plots check the Gaussian assumptions:

``` r

plot(fit, "residuals")
plot(fit, "qq")
```

![](glm_files/figure-html/plot-resid-qq-1.png)![](glm_files/figure-html/plot-resid-qq-2.png)

On the left, the deviance residuals (for a Gaussian model, the ordinary
residuals) are plotted against the fitted values. They fall in vertical
stripes, one per cell mean. The two 2 mg/day cells have almost the same
mean (26.06 and 26.14), so their stripes overlap. Stripes of clearly
different heights mean unequal variances, and a funnel that widens with
the fitted value means a variance that grows with the mean, which
suggests another family. Here the heights vary somewhat, as the cell
standard deviations did, but there is no funnel.

On the right, the Q-Q plot compares the residuals with a normal
distribution; systematic curvature at the ends would mean skewness or
heavy tails. The points follow the line well. The lower tail is a little
shorter than a normal one, a mild departure for a balanced design with
60 observations.

If the variance of a Gamma model really is proportional to the mean
squared, its deviance residuals have roughly the same spread at every
fitted value. For the car prices the stripes differ in height:

``` r

plot(gam, "residuals")
```

![](glm_files/figure-html/gamma-residuals-1.png)

In a one-factor model the fitted values are the group means. Two stripes
sit side by side: vans (19.10) on the left and sporty cars (19.39) on
the right. The vans’ stripe is by far the shortest, because their prices
vary much less than the model assumes. The midsize stripe, at 27.22, is
the tallest.

The deviance residuals of a Gamma model are approximately normal when
the model fits. Here the middle follows the line and the upper tail
rises above it, so a few cars are dearer than a Gamma distribution with
this dispersion would predict.

``` r

plot(gam, "qq")
```

![](glm_files/figure-html/gamma-qq-1.png)

`gam$model` is an ordinary `glm`, so you can find those cars directly:

``` r

res <- residuals(gam$model, type = "deviance")
head(cars[order(res, decreasing = TRUE), ], 3)
#>     Manufacturer    Model    Type Price
#> 59 Mercedes-Benz     300E Midsize  61.9
#> 19     Chevrolet Corvette  Sporty  38.0
#> 48      Infiniti      Q45 Midsize  47.9
```

The Gamma `emmeans` plot shows the means on the response scale, with
their asymmetric intervals. `Type` has no natural order, so the points
are not joined by a line; lines appear only across a second grouping
factor, as in the ToothGrowth plot. The fourth plot, `plot(gam, "box")`,
shows the raw prices by type.

``` r

plot(gam, "emmeans")
```

![](glm_files/figure-html/gamma-emm-plot-1.png)

All plots are ggplot2 objects that you can change with `+`; see [Plots
and themes](https://elkronos.github.io/anovakit/articles/visuals.md).

## Options

### `ci_method`

`ci_method = "profile"` (the default) cuts the profile likelihood at the
t quantiles the tests beside it use. The p-values are Wald p-values,
though, so in a small sample a profile interval can exclude 0 while its
p-value is above 0.05, or the reverse; `$notes` names any such row.
`ci_method = "wald"` gives estimate ± t × SE, which always agrees with
the p-value:

``` r

wald <- anova_glm(cars, "Price", "Type", family = Gamma(link = "log"),
                  ci_method = "wald", plots = FALSE)$assumptions$coefficients
cbind(co[, c("term", "estimate")],
      profile = co[, c("conf_low", "conf_high")],
      wald = wald[, c("conf_low", "conf_high")])
#>          term    estimate profile.conf_low profile.conf_high wald.conf_low
#> 1 (Intercept)  2.90210817        2.7384042         3.0752703    2.73380885
#> 2   TypeLarge  0.28836818        0.0264585         0.5545836    0.02469382
#> 3 TypeMidsize  0.40177703        0.1790690         0.6219030    0.18058839
#> 4   TypeSmall -0.58299378       -0.8077507        -0.3604911   -0.80638884
#> 5  TypeSporty  0.06279664       -0.1832054         0.3101515   -0.18356834
#> 6     TypeVan  0.04758016       -0.2297083         0.3322378   -0.23291870
#>   wald.conf_high
#> 1      3.0704075
#> 2      0.5520425
#> 3      0.6229657
#> 4     -0.3595987
#> 5      0.3091616
#> 6      0.3280790
```

The two differ in the third decimal place here. Profile intervals follow
the shape of the likelihood, which matters more in small samples or near
a boundary; Wald intervals are symmetric on the log scale. Exponentiate
either to get an interval for the ratio.

### `vcov_type`

`vcov_type` takes an HC type (`"HC0"` to `"HC4"`) from the sandwich
package. The robust covariance does not rely on the family’s variance
function, so it protects the marginal means and comparisons when that
function is wrong, as the vans suggest it is here:

``` r

rob <- anova_glm(cars, "Price", "Type", family = Gamma(link = "log"),
                 vcov_type = "HC3", plots = FALSE)
cbind(rob$emmeans[, "Type", drop = FALSE],
      se_model = gam$emmeans$se, se_HC3 = rob$emmeans$se)
#>      Type  se_model    se_HC3
#> 1 Compact 1.5421305 1.7265476
#> 2   Large 2.4815412 2.0040958
#> 3 Midsize 1.9654379 2.6764077
#> 4   Small 0.7514161 0.4367684
#> 5  Sporty 1.7554535 2.2117883
#> 6     Van 2.1563719 0.6640312
grep("omnibus|Robust", rob$notes, value = TRUE)
#> [1] "The omnibus F test compares deviances scaled by a dispersion estimated from the Pearson residuals, so it uses the model's variance function, not the robust (HC3) covariance, which only the coefficient table, the marginal means and the comparisons use. For an omnibus test that uses it, set test_statistic = \"Wald\"."
#> [2] "Robust standard errors (HC3) were requested, so the coefficient intervals are Wald intervals built from them; profile-likelihood intervals cannot use a sandwich covariance."
```

The standard error for vans shrinks and the one for midsize cars grows:
each now reflects that type’s own spread. The robust covariance reaches
the coefficient table (with Wald intervals, since a profile likelihood
cannot use it), the marginal means and the comparisons. The omnibus F
test compares deviances and does not use it; `test_statistic = "Wald"`
gives an omnibus test that does:

``` r

anova_glm(cars, "Price", "Type", family = Gamma(link = "log"),
          vcov_type = "HC3", test_statistic = "Wald", plots = FALSE)$anova
#>   term df statistic      p_value
#> 1 Type  5  195.7549 2.298221e-40
```

[Weights and robust standard
errors](https://elkronos.github.io/anovakit/articles/weights.md) covers
`vcov_type` across the package, including when the function falls back
to the model-based covariance.

### `weights`

`weights` names a column of prior weights, and what a weight means
depends on the family and the response. Whole numbers on a 0/1 response
or on counts (binomial, Poisson and negative binomial families) are
frequency weights, and whole numbers on a binomial proportion are
numbers of trials. Any other weights, including all weights in the
Gaussian and Gamma families, are used as
[`glm()`](https://rdrr.io/r/stats/glm.html) uses them, with each row one
unit. [Weights and robust standard
errors](https://elkronos.github.io/anovakit/articles/weights.md) sets
out the rules; this section shows one use.

A response that is itself an average suits a Gamma model with weights.
The mean of m prices has variance φμ²/m, so its prior weight is m.
Suppose you had only the average price per manufacturer and type, with
the number of models behind each:

``` r

agg <- aggregate(Price ~ Manufacturer + Type, data = cars, FUN = mean)
agg$n_models <- aggregate(Price ~ Manufacturer + Type, data = cars,
                          FUN = length)$Price
table(agg$n_models)
#> 
#>  1  2 
#> 69 12
```

Weighting each average by its count recovers the marginal means of the
full data exactly; ignoring the counts does not:

``` r

w_fit <- anova_glm(agg, "Price", "Type", family = Gamma(link = "log"),
                   weights = "n_models", plots = FALSE)
u_fit <- anova_glm(agg, "Price", "Type", family = Gamma(link = "log"),
                   plots = FALSE)
data.frame(Type = gam$emmeans$Type, full_data = gam$emmeans$estimate,
           weighted = w_fit$emmeans$estimate,
           unweighted = u_fit$emmeans$estimate)
#>      Type full_data weighted unweighted
#> 1 Compact  18.21250 18.21250   18.60000
#> 2   Large  24.30000 24.30000   24.50500
#> 3 Midsize  27.21818 27.21818   27.31000
#> 4   Small  10.16667 10.16667   10.36875
#> 5  Sporty  19.39286 19.39286   19.16667
#> 6     Van  19.10000 19.10000   19.43125
```

The standard errors differ slightly from the full-data fit, because the
dispersion is now estimated from 81 averages rather than 93 cars. A row
with weight zero is dropped before fitting, counted in `$n_removed` and
recorded in `$notes`.

### `interaction`

The default, `interaction = FALSE`, fits the grouping factors
additively. The marginal means and comparisons are then reported for
each factor on its own, averaged over the others, with a `term` column
saying which, and each factor’s comparisons form their own multiplicity
family:

``` r

add <- anova_glm(tg, "len", c("dose", "supp"), plots = FALSE)
add$emmeans
#>   term dose supp estimate        se df  conf_low conf_high
#> 1 dose  0.5 <NA> 10.60500 0.8558752 56  8.890476  12.31952
#> 2 dose    1 <NA> 19.73500 0.8558752 56 18.020476  21.44952
#> 3 dose    2 <NA> 26.10000 0.8558752 56 24.385476  27.81452
#> 4 supp <NA>   OJ 20.66333 0.6988192 56 19.263430  22.06324
#> 5 supp <NA>   VC 16.96333 0.6988192 56 15.563430  18.36324
add$notes
#> [1] "The grouping factors enter the model additively, so $emmeans and $posthoc are reported for each factor separately, averaged over the others (column `term`), and each factor's comparisons are a separate multiplicity family. Set interaction = TRUE to compare cells instead."
```

`interaction = TRUE` fits the full factorial. With three or more factors
a whole number sets the highest order of interaction: `interaction = 2`
keeps the two-way interactions and drops the rest.

### `adjust`, `conf_level` and `posthoc`

`adjust` is applied by emmeans and accepts `"tukey"` (the default),
`"sidak"`, `"scheffe"`, `"dunnettx"`, `"bonferroni"`, `"holm"`,
`"hochberg"`, `"hommel"`, `"BH"`, `"BY"`, `"fdr"` or `"none"`. A
step-down method such as Holm cannot be turned into intervals, so
emmeans gives Bonferroni intervals beside Holm p-values, and `$notes`
says so:

``` r

holm <- anova_glm(tg, "len", c("dose", "supp"), interaction = TRUE,
                  adjust = "holm", plots = FALSE)
holm$notes
#> [1] "The comparison intervals use the bonferroni adjustment while the p-values use holm: emmeans cannot turn holm into intervals."
```

`conf_level` sets the level of every interval: the coefficient table,
the marginal means and the comparisons. `posthoc = FALSE` skips the
pairwise comparisons, of which there are `choose(k, 2)` for k cells
(counted per factor in an additive model); `$notes` records that they
were skipped and how many there would have been.

## What the notes say

``` r

fit$notes
#> character(0)
gam$notes
#> [1] "McFadden's pseudo R squared is not reported (NA): the log-likelihood of the Gamma family is a log density, which shifts with the units of the response, so the ratio has no unit-free meaning. The deviance explained does not depend on the units and is reported instead."
```

The ToothGrowth fit has no notes: no rows were dropped, the model is far
from saturated, nothing failed to compute, and the Gaussian family has
partial eta squared to report. The Gamma fit has one, explaining why
McFadden’s pseudo R squared is `NA`. Other notes in this walkthrough
record where the robust covariance is used, the per-factor reporting of
an additive model, Bonferroni intervals beside Holm p-values, and
differences rather than ratios under the inverse link.

## Reporting the result

These write-ups are built from the fits with inline R, so the numbers
cannot drift from the analysis:

> A two-way analysis of variance (Type II) found that odontoblast length
> depended on the dose of vitamin C, F(2, 54) = 92.00, p \< 0.001,
> partial η² = 0.77, and on the supplement, F(1, 54) = 15.57, p \<
> 0.001, with an interaction between them, F(2, 54) = 4.11, p = 0.022.
> Orange juice produced longer odontoblasts than ascorbic acid at 0.5
> mg/day (difference 5.25, Tukey-adjusted p = 0.024) and 1 mg/day (5.93,
> adjusted p = 0.007), but not at 2 mg/day (-0.08, adjusted p = 1.000).

> In a Gamma model with a log link, car prices differed by type, F(5,
> 87) = 18.56, p \< 0.001, with type explaining 54% of the deviance.
> Midsize cars cost 2.68 times as much as small cars on average (95% CI
> 1.98 to 3.62, Tukey-adjusted).

For a model on a non-identity link, name the family and link, and report
comparisons in the form they take on it: ratios for a log link.

## See also

- [`anova_bin()`](https://elkronos.github.io/anovakit/reference/anova_bin.md)
  ([walkthrough](https://elkronos.github.io/anovakit/articles/binary.md))
  and
  [`anova_count()`](https://elkronos.github.io/anovakit/reference/anova_count.md)
  ([walkthrough](https://elkronos.github.io/anovakit/articles/counts.md))
  for binary and count responses.
- [`anova_ancova()`](https://elkronos.github.io/anovakit/reference/anova_ancova.md)
  ([walkthrough](https://elkronos.github.io/anovakit/articles/ancova.md))
  to add a numeric covariate.
- [`anova_welch()`](https://elkronos.github.io/anovakit/reference/anova_welch.md)
  ([walkthrough](https://elkronos.github.io/anovakit/articles/welch.md))
  for a one-way comparison of means with unequal variances.
- [Weights and robust standard
  errors](https://elkronos.github.io/anovakit/articles/weights.md) for
  every kind of weight and for `vcov_type`.
- [Working with
  results](https://elkronos.github.io/anovakit/articles/results.md) for
  the fit object, `$emmeans_object` and reporting.
- [Diagnostics](https://elkronos.github.io/anovakit/articles/diagnostics.md)
  and [Plots and
  themes](https://elkronos.github.io/anovakit/articles/visuals.md).
