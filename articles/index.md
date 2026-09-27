# Articles

### Start here

- [Getting started with
  anovakit](https://elkronos.github.io/anovakit/articles/anovakit.md):

- [Which function? Choosing an
  analysis](https://elkronos.github.io/anovakit/articles/choosing.md):

  A decision guide: four questions about your data that lead to one of
  the eight functions, the situations each one is built for, and the
  mistakes the choice is meant to avoid.

- [Working with results: the anovakit_fit
  object](https://elkronos.github.io/anovakit/articles/results.md):

  What every function returns, how to print, extract and reuse it, how
  to go further with emmeans and the fitted model, and how to report a
  result.

### Walkthroughs

One article per function: the data, the fit, the output, effect sizes,
comparisons, diagnostics, plots, options and a model write-up.

- [Welch's ANOVA: comparing means when variances
  differ](https://elkronos.github.io/anovakit/articles/welch.md):

  A walkthrough of anova_welch(): the omnibus Welch test, pairwise Welch
  t tests, standardised mean differences that do not assume equal
  variances, and the three diagnostic plots.

- [Kruskal-Wallis: comparing groups by
  ranks](https://elkronos.github.io/anovakit/articles/kruskal-wallis.md):

  A walkthrough of anova_kw(): the Kruskal-Wallis test, epsilon squared,
  Dunn’s comparisons with the tie correction, missing and infinite
  values, and the optional diagnostics.

- [ANCOVA: comparing groups at the same value of a
  covariate](https://elkronos.github.io/anovakit/articles/ancova.md):

  A walkthrough of anova_ancova() on the anorexia treatment data:
  covariate-adjusted means, the homogeneity-of-slopes test, simple
  slopes, centring and robust tests.

- [Repeated measures ANOVA: the same subjects under several
  conditions](https://elkronos.github.io/anovakit/articles/repeated-measures.md):

  A walkthrough of anova_rm(): mixed designs, sphericity corrections,
  generalised eta squared, incomplete subjects and within-subject
  residual checks.

- [MANOVA: comparing groups on several responses at
  once](https://elkronos.github.io/anovakit/articles/manova.md):

  A walkthrough of anova_manova() on the iris measurements: the
  multivariate test, univariate follow-ups, canonical discriminant
  analysis, Box’s M and Mardia’s tests, and MANCOVA.

- [Binary outcomes: logistic regression with
  anova_bin()](https://elkronos.github.io/anovakit/articles/binary.md):

  A walkthrough of anova_bin() on the Berkeley admissions data: analysis
  of deviance, odds ratios, marginal probabilities, frequency weights
  and separation.

- [Counts and rates: analysis of deviance with
  anova_count()](https://elkronos.github.io/anovakit/articles/counts.md):

  A walkthrough of anova_count() on the warpbreaks and insurance-claims
  data: Poisson or negative binomial, rate ratios, marginal rates,
  exposure offsets, frequency weights and groups with no events.

- [Analysis of deviance: factorial designs and non-normal
  responses](https://elkronos.github.io/anovakit/articles/glm.md):

  A walkthrough of anova_glm(): a two-way Gaussian ANOVA with Type II
  and Type III tests, then a Gamma model for a positive, right-skewed
  response, with its coefficient table, dispersion, robust standard
  errors and weights.

### Guides

Topics that cut across every function.

- [Diagnostics: what anovakit checks and what to do about
  it](https://elkronos.github.io/anovakit/articles/diagnostics.md):

  Every assumption check the eight functions run, organised by what is
  being checked: how to read each one, a worked example that makes it
  fire, and the argument or function that helps when it fails.

- [Plots: every figure anovakit returns, and how to customise
  them](https://elkronos.github.io/anovakit/articles/visuals.md):

  A tour of every ggplot2 figure the eight functions build, grouped by
  kind, with what to look for in each; then how to change them, build
  your own from the result tables, combine and save them.
