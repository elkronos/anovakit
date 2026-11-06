# Articles

### Start here

- [Get started with
  anovakit](https://elkronos.github.io/anovakit/articles/anovakit.md):

  One analysis from start to finish (Welch’s ANOVA on mtcars, then the
  same question with a covariate in anova_ancova()), and the arguments
  every function shares.

- [Choose an
  analysis](https://elkronos.github.io/anovakit/articles/choosing.md):

  A decision guide: four questions about your data that lead to one of
  the eight functions, a table of common situations, and what each
  function is for.

### One outcome, independent groups

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

  A walkthrough of anova_ancova() on the anorexia treatment data: the
  homogeneity-of-slopes test, simple slopes, covariate centring,
  adjusted marginal means and robust standard errors (vcov_type).

### Repeated measures and several outcomes

- [Repeated measures ANOVA: the same subjects under several
  conditions](https://elkronos.github.io/anovakit/articles/repeated-measures.md):

  A walkthrough of anova_rm(): mixed designs, sphericity corrections,
  generalised eta squared, within-subject residual checks, incomplete
  subjects and repeated rows.

- [MANOVA: comparing groups on several responses at
  once](https://elkronos.github.io/anovakit/articles/manova.md):

  A walkthrough of anova_manova() on the iris measurements: the
  multivariate test, univariate follow-ups, canonical discriminant
  analysis, Box’s M and Mardia’s tests, MANCOVA and several grouping
  variables.

### Binary, count and other outcomes

- [Binary outcomes: logistic regression with
  anova_bin()](https://elkronos.github.io/anovakit/articles/binary.md):

  A walkthrough of anova_bin() on the Berkeley admissions data: analysis
  of deviance, odds ratios, marginal probabilities, separation,
  frequency weights and robust standard errors.

- [Counts and rates: analysis of deviance with
  anova_count()](https://elkronos.github.io/anovakit/articles/counts.md):

  A walkthrough of anova_count() on the warpbreaks and insurance-claims
  data: Poisson or negative binomial, rate ratios, marginal rates,
  exposure offsets, groups with no events and frequency weights.

- [Any GLM family: analysis of deviance with
  anova_glm()](https://elkronos.github.io/anovakit/articles/glm.md):

  A walkthrough of anova_glm(): a two-way Gaussian ANOVA with Type II
  and Type III tests, and a Gamma model with a log link for a positive,
  right-skewed response, with its coefficient table, dispersion, ratios
  of means and robust standard errors.

### Guides

Topics that cut across the functions.

- [Working with results: the anovakit_fit
  object](https://elkronos.github.io/anovakit/articles/results.md):

  What every function returns, how to print, extract and reuse it, how
  to go further with emmeans and the fitted model, and how to report a
  result.

- [Diagnostics: what anovakit checks and what to do about
  it](https://elkronos.github.io/anovakit/articles/diagnostics.md):

  Every assumption check the eight functions run, organised by what is
  being checked: how to read each one, a worked example that makes it
  fire, and the argument or function that helps when it fails. Also
  covers exact fits and the limits of numerical precision.

- [Plots and
  themes](https://elkronos.github.io/anovakit/articles/visuals.md):

  Every ggplot2 figure the eight functions build, grouped by kind, with
  what to look for in each; then how to change them, set a theme for all
  of them, build your own from the result tables, combine and save them.

- [Weights and robust standard
  errors](https://elkronos.github.io/anovakit/articles/weights.md):

  How anova_bin(), anova_count() and anova_glm() read prior weights as
  frequencies, trials or sampling weights, what happens to zero weights,
  and where vcov_type’s robust covariance reaches and when it is not
  used.

- [Design
  choices](https://elkronos.github.io/anovakit/articles/design.md):

  Why the defaults are what they are, which decisions a function makes
  for you and where it records them, and which settings can and cannot
  change a result.
