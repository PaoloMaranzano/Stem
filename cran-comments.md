# CRAN submission comments

## Submission

`Stem` was archived from CRAN and is resubmitted here as version 2.0.0. The
maintainer address is `pmaranzano.ricercastatistica@gmail.com`, the same one
used for the maintainer's other CRAN packages.

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
  `man/SCSTEM_Estim.Rd` and `man/SCSTEM_Infocrit.Rd`, are regenerated from
  roxygen sources that use `TRUE`/`FALSE`.

* **"You write information messages to the console that cannot be easily
  suppressed."**
  All the `print()`/`cat()` diagnostics in `R/kalman.R`, `R/SCSTEM_Estim.R`,
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
* A second dataset (`povalley`), three vignettes, a `testthat` suite and a
  conceptual map of the package have been added.

## Test environments

* local Windows 11, R 4.5.1

## R CMD check results

See below; the check is run with `--as-cran`.

## Note on running the checks locally

The repository lives inside a OneDrive-synced folder. `R CMD check` creates
thousands of small files, and the sync client stalls the run — the check hangs
at "checking package dependencies" for as long as it is left there. Run the
build and the check from a directory that is not synced:

```r
# from any non-synced working directory
pkg <- "<path to the working copy of this repository>"
tar <- devtools::build(pkg, path = tempdir())
rcmdcheck::rcmdcheck(tar, args = c("--as-cran", "--no-manual"),
                     error_on = "warning")
```

Building the vignettes also needs pandoc on the path. If it is not, point R at
the copy that ships with RStudio:

```r
Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")
```
