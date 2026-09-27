# anovakit: Analysis of Variance Workflows

A consistent set of wrappers around the analysis of variance designs
that applied work most often needs. Every function takes `data` first
and names its columns as character strings – six take `response` and
`groups`,
[`anova_manova`](https://elkronos.github.io/anovakit/reference/anova_manova.md)
takes `responses` (plural), and
[`anova_rm`](https://elkronos.github.io/anovakit/reference/anova_rm.md)
takes `response` with `subject`, `within` and `between`. All eight
validate their input the same way and return the same S3 class,
[`anovakit_fit`](https://elkronos.github.io/anovakit/reference/anovakit_fit.md).

## Choosing a function

- [`anova_welch`](https://elkronos.github.io/anovakit/reference/anova_welch.md):

  One numeric response, unequal variances permitted. The default choice
  for a one-way comparison of means.

- [`anova_kw`](https://elkronos.github.io/anovakit/reference/anova_kw.md):

  One numeric response, no distributional assumption. Kruskal-Wallis
  with Dunn post-hoc. Convert an ordered factor with
  [`as.integer()`](https://rdrr.io/r/base/integer.html) first.

- [`anova_ancova`](https://elkronos.github.io/anovakit/reference/anova_ancova.md):

  One numeric response, adjusting for one or more numeric covariates.

- [`anova_rm`](https://elkronos.github.io/anovakit/reference/anova_rm.md):

  One numeric response measured repeatedly on the same subjects.

- [`anova_manova`](https://elkronos.github.io/anovakit/reference/anova_manova.md):

  Two or more numeric responses analysed jointly, optionally adjusting
  for covariates.

- [`anova_bin`](https://elkronos.github.io/anovakit/reference/anova_bin.md):

  A binary response, via logistic regression.

- [`anova_count`](https://elkronos.github.io/anovakit/reference/anova_count.md):

  A count response, via Poisson or negative binomial regression.

- [`anova_glm`](https://elkronos.github.io/anovakit/reference/anova_glm.md):

  Any other generalised linear model family.

## Shared conventions

- `conf_level` sets every interval in the tables the function returns.
  The one thing it does not reach is `$emmeans_object`, which is
  emmeans' own grid and carries emmeans' default level; pass `level =`
  yourself when you summarise it.

- `plots = TRUE` builds ggplot2 objects and returns them; nothing is
  ever drawn as a side effect.

- `verbose = FALSE` by default; progress goes through
  [`message`](https://rdrr.io/r/base/message.html), never
  [`cat`](https://rdrr.io/r/base/cat.html). No analysis function writes
  to the console: warnings and messages raised by car,
  [`glm`](https://rdrr.io/r/stats/glm.html) and afex are captured into
  `$notes` instead.

- Rows with missing or infinite values in the analysed columns are
  dropped. The exception is the response of
  [`anova_kw`](https://elkronos.github.io/anovakit/reference/anova_kw.md):
  a rank test can use infinite values (as the best or worst possible
  outcome), so it keeps them. `$n_removed` counts the input rows that
  are not in `$data_used`, so `nrow($data_used) + $n_removed` is always
  the number of rows supplied.

- Grouping columns are coerced to factors; character columns get their
  levels in C-locale order, so the reference level does not depend on
  the session. A factor level that is itself `NA` counts as missing.

- [`anova_welch`](https://elkronos.github.io/anovakit/reference/anova_welch.md)
  and
  [`anova_kw`](https://elkronos.github.io/anovakit/reference/anova_kw.md)
  combine several grouping columns into one factor of the cells that
  contain data. The model-based functions fit the grouping columns as
  factors instead. When those factors enter the model additively,
  marginal means and pairwise comparisons are reported for each factor
  separately, averaged over the others. When the model contains their
  interaction, they are reported for the cells: a level combination with
  no data is left out when the model cannot estimate it, and named in
  `$notes` when the model estimates it by extrapolation from the other
  cells.

- `posthoc = FALSE` skips the pairwise comparisons and says in `$notes`
  how many there would have been. Above 5000 comparisons they are
  skipped by default, with a note.

- Prior weights, where a function accepts them, are named as a column
  rather than passed as a vector, so they are subsetted with the data.
  Rows with a weight of zero contribute nothing and are dropped (and
  counted in `$n_removed`). Whole-number weights on a count or 0/1
  response are frequency weights: robust standard errors and residual
  degrees of freedom are then those of the data expanded to one row per
  observation.

- `vcov_type`, where a function accepts it, asks for
  heteroscedasticity-consistent (sandwich) standard errors. Where the
  sandwich is degenerate – fitted values on the boundary (separation, a
  group with no events) or a cell with a single row – the model-based
  covariance is used instead, with a note.

- A column may play only one role (response, group, covariate, weights,
  offset, subject); the names `Residuals` and `(Intercept)` are
  reserved.

- Anything the function decided on your behalf, or could not compute, is
  recorded in `$notes`. Read it.

## Multiplicity adjustments

`adjust` takes different values depending on where the comparisons come
from.
[`anova_welch`](https://elkronos.github.io/anovakit/reference/anova_welch.md)
and
[`anova_kw`](https://elkronos.github.io/anovakit/reference/anova_kw.md)
compare directly and accept the
[`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) methods (`"holm"`,
`"hochberg"`, `"hommel"`, `"bonferroni"`, `"BH"`, `"BY"`, `"fdr"`,
`"none"`), defaulting to `"holm"` and `"BH"` respectively. The other six
go through emmeans and additionally accept `"tukey"` (their default),
`"sidak"`, `"scheffe"` and `"dunnettx"`.

Every `$posthoc` table reports the unadjusted p-value in `p_value`, the
adjusted one in `p_adjusted` and the method in `adjustment`. The
intervals of
[`anova_welch`](https://elkronos.github.io/anovakit/reference/anova_welch.md)
are per-comparison intervals; those of the emmeans-based functions are
adjusted (simultaneous) intervals, and a note says when emmeans had to
use a different adjustment for them (it cannot turn a step-down method
such as Holm into intervals, and uses Bonferroni instead). Tukey,
Dunnett and Sidak p-values are computed to a limited precision
([`ptukey()`](https://rdrr.io/r/stats/Tukey.html)'s upper tail is good
to about 1e-14), so [`print()`](https://rdrr.io/r/base/print.html) shows
those below 1e-10 (1e-12 for Sidak) as a bound; the stored values are
unchanged.

## See also

[`anovakit_fit`](https://elkronos.github.io/anovakit/reference/anovakit_fit.md)
for the returned object.

## Author

**Maintainer**: Justin Chase <jchase.msu@gmail.com>
