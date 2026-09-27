# Package index

## Analysis functions

Each takes a data frame and column names given as character strings, and
returns an `anovakit_fit`.

### Continuous outcomes

- [`anova_welch()`](https://elkronos.github.io/anovakit/reference/anova_welch.md)
  : Welch's Analysis of Variance
- [`anova_kw()`](https://elkronos.github.io/anovakit/reference/anova_kw.md)
  : Kruskal-Wallis Test with Dunn Post-hoc Comparisons
- [`anova_ancova()`](https://elkronos.github.io/anovakit/reference/anova_ancova.md)
  : Analysis of Covariance
- [`anova_rm()`](https://elkronos.github.io/anovakit/reference/anova_rm.md)
  : Repeated Measures Analysis of Variance
- [`anova_manova()`](https://elkronos.github.io/anovakit/reference/anova_manova.md)
  : Multivariate Analysis of Variance and Covariance

### Binary, count and other outcomes

- [`anova_bin()`](https://elkronos.github.io/anovakit/reference/anova_bin.md)
  : Analysis of Deviance for a Binary Response
- [`anova_count()`](https://elkronos.github.io/anovakit/reference/anova_count.md)
  : Analysis of Deviance for a Count Response
- [`anova_glm()`](https://elkronos.github.io/anovakit/reference/anova_glm.md)
  : Analysis of Deviance for a Generalised Linear Model

## The fitted object

What every function returns, and the methods for it.

- [`anovakit_fit`](https://elkronos.github.io/anovakit/reference/anovakit_fit.md)
  : The object returned by every anovakit function
- [`print(`*`<anovakit_fit>`*`)`](https://elkronos.github.io/anovakit/reference/print.anovakit_fit.md)
  : Print an anovakit result
- [`summary(`*`<anovakit_fit>`*`)`](https://elkronos.github.io/anovakit/reference/summary.anovakit_fit.md)
  [`print(`*`<summary.anovakit_fit>`*`)`](https://elkronos.github.io/anovakit/reference/summary.anovakit_fit.md)
  : Summarise an anovakit result
- [`plot(`*`<anovakit_fit>`*`)`](https://elkronos.github.io/anovakit/reference/plot.anovakit_fit.md)
  : Plot an anovakit result

## Package overview

- [`anovakit`](https://elkronos.github.io/anovakit/reference/anovakit-package.md)
  [`anovakit-package`](https://elkronos.github.io/anovakit/reference/anovakit-package.md)
  : anovakit: Analysis of Variance Workflows
