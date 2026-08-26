# Submission comments

## Test environments

* local: macOS 26.5.2, aarch64-apple-darwin23, R 4.6.1 — 0 errors | 0 warnings | 2 notes
* win-builder devel and release — not yet run
* macOS builder — not yet run
* R-hub: windows-x86_64-devel, ubuntu-gcc-release, fedora-clang-devel — not yet run

Fill in the remaining rows before submitting.

## R CMD check results

Two NOTEs locally, neither of which is a package issue.

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
runs unconditionally, the slowest takes 1.2 seconds, and the full set takes
under three.

Every use of a suggested package ({afex}, {MASS}, {sandwich}) is behind
`requireNamespace()` in R/ and `skip_if_not_installed()` in the tests. The suite
passes with all three absent.

The methods cited in the Description are implemented directly rather than
delegated, so the DOIs point at the papers the code follows: Dunn's test with
the standard tie correction, Mardia's multivariate skewness and kurtosis, and
Box's M.
