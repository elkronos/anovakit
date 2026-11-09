## Submission

This is a first submission.

## Test environments

* GitHub Actions: macOS (release), Windows (release), Ubuntu (devel, release,
  oldrel-1)
* local: macOS, aarch64-apple-darwin, R 4.6.1
* local: Ubuntu 24.04, x86_64-pc-linux-gnu, R 4.3.3; also with
  `_R_CHECK_DEPENDS_ONLY_=true`, so that no suggested package is available
* win-builder: R-devel and R-release
* macOS builder: R-release

## R CMD check results

0 errors | 0 warnings | 1 note

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Justin Chase <jchase.msu@gmail.com>'

New submission
```

If the note also lists possibly misspelled words in DESCRIPTION, these are the
surnames of the authors of the cited methods (Kruskal, Wallis, Dunn,
Greenhouse, Geisser, Huynh, Feldt, Mardia), each given with a DOI.

## Notes for the reviewer

`anova_rm()` needs the suggested package 'afex', so its example is guarded
with `requireNamespace()`. No example uses `\dontrun{}` or `\donttest{}`.

Every use of a suggested package in R/ is behind `requireNamespace()`, and the
tests that need 'afex', 'bit64', 'data.table', 'effectsize' or 'sandwich' skip
when it is absent. The suite passes with no suggested package installed.

On CRAN the tests run the core suite, which covers every function. A further
set of detailed regression tests is skipped there with `skip_on_cran()` to
keep the check time short; it runs on every push in continuous integration.

Three of the methods cited in the Description are implemented directly rather
than delegated, so their DOIs point at the papers the code follows: Dunn's test
with the standard tie correction, Mardia's multivariate skewness and kurtosis,
and Box's M. Welch's test comes from `stats::oneway.test()`, and the
Greenhouse-Geisser and Huynh-Feldt corrections from 'afex' and 'car', except
when the error matrix is singular and 'car' cannot give the epsilons: the
package then computes both itself, following Greenhouse and Geisser (1959) and
Huynh and Feldt (1976) with Lecoutre's (1991) correction.
