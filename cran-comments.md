# Submission comments

## Test environments

* local: macOS 26.5.2, aarch64-apple-darwin23, R 4.6.1
* Ubuntu 24.04, x86_64-pc-linux-gnu, R 4.3.3 (Ubuntu archive builds of every
  dependency)
* GitHub Actions (.github/workflows/R-CMD-check.yaml): macOS, Windows and
  Ubuntu on R release, Ubuntu on R devel and oldrel-1

Still to run before submission, since none can be run from the machines above:
win-builder (devel and release), the macOS builder, and R-hub
(`rhub::rhub_check()`). Record their results here.

## R CMD check results

0 errors | 0 warnings | the notes below.

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Justin Chase <jchase.msu@gmail.com>'
New submission
```

This is a first submission, so the note is expected.

```
* checking HTML version of manual ... NOTE
Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.
```

The local machine ships an old HTML Tidy. The manual builds without error and
this check is skipped rather than failed.

## Notes for the reviewer

`anova_rm()` requires the suggested package {afex}, so its example is guarded
with `requireNamespace()`. No example is wrapped in `\donttest{}`: every one
runs unconditionally.

Every use of a suggested package ({afex}, {MASS}, {sandwich}) is behind
`requireNamespace()` in R/ and `skip_if_not_installed()` in the tests. The suite
passes with all three absent.

Three of the methods cited in the Description are implemented directly rather
than delegated, so their DOIs point at the papers the code follows: Dunn's test
with the standard tie correction, Mardia's multivariate skewness and kurtosis,
and Box's M. Welch's test comes from `stats::oneway.test()`, and the
Greenhouse-Geisser and Huynh-Feldt corrections from {afex} and {car}.
