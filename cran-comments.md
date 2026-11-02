# CRAN submission comments

## Submission

`Stem` 1.0 was published on CRAN in 2012 and archived on 2022-06-13, as check
problems were not corrected. It is resubmitted here as version 2.0.0, a
thorough revision that corrects those problems and extends the package.

The maintainer changes: version 1.0 was maintained by Michela Cameletti
(`michela.cameletti@unibg.it`), who remains an author of version 2.0.0 and
agrees that Paolo Maranzano takes over as maintainer. The maintainer address is
`pmaranzano.ricercastatistica@gmail.com`, the same one used for the
maintainer's other CRAN packages.

## Response to the previous review

All the points raised in the previous review have been addressed.

* **"Please omit the redundant 'in R' from the end of your title."**
  The title is now `Spatio-Temporal Expectation-Maximization Models`.

* **"Please write references in the description of the DESCRIPTION file in the
  form authors (year) <doi:...> ..."**
  Every reference in `Description` now carries a `<doi:...>` link with no space
  after `doi:`, and the EM monograph is cited as `(2008, ISBN:9780470191613)`.

* **"Please write TRUE and FALSE instead of T and F."**
  `T` and `F` no longer appear as values or as names anywhere in `R/` or in the
  documentation. The two manual pages named in the review,
  `man/SCSTEM_Estimation.Rd` and `man/SCSTEM_Infocrit.Rd`, are regenerated from
  roxygen sources that use `TRUE`/`FALSE`.

* **"You write information messages to the console that cannot be easily
  suppressed."**
  All the `print()`/`cat()` diagnostics in `R/kalman.R`, `R/SCSTEM_Estimation.R`,
  `R/STEM_Bootstrap.fn.R`, `R/STEM_Estimation.R` and `R/SCSTEM_Bootstrap.R`
  have been replaced by `message()` calls guarded by a new `verbose` argument,
  which defaults to `FALSE`. The only remaining `cat()` calls are inside
  `print()` methods for the classes introduced by the package, where the review
  explicitly allows them.

* **"Please do not modify the .GlobalEnv."**
  The stray `.Random.seed` assignment has been removed. Where a seed is needed
  for reproducibility, an internal helper sets it and restores the previous RNG
  state on exit, so the user's workspace is left unchanged; this is covered by a
  regression test.

## Further changes made for this submission

* The dataset `pm10` is now stored as `data/pm10.rda` instead of a source file
  duplicated in `R/` and `data/`.
* The dependency on `SCDA`, which is not distributed on CRAN, has been removed;
  `dplyr` and `sf` are no longer needed either.
* A second dataset (`povalley`), five vignettes, a `testthat` suite and a
  conceptual map of the package have been added.

## Test environments

* local Windows 11, R 4.6.1 (`R CMD check --as-cran`, with the `\donttest`
  examples and the HTML math check)

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission; package was archived on CRAN. This is the resubmission of
  the archived package described above.

The words that a spell check may flag in `Description` are the names of the
cited authors (Fasso, Cameletti, Nicolis, McLachlan, Krishnan, Sugasawa,
Murakami, Cerqueti, Maranzano, Mattera) and technical terms (Kalman, Potts,
STEM, SC-STEM).

The examples marked `\donttest` fit spatio-temporal models and bootstraps on
real data and take about four minutes in all; the examples that are run by
default take a few seconds.
