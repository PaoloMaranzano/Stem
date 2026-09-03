# CRAN submission notes


## Maintainer-requested changes addressed

- Removed the redundant phrase "in R" from the package title.
- Added DOI/ISBN/URL-style links to references in `DESCRIPTION`.
- Replaced `T`/`F` defaults and comparisons with `TRUE`/`FALSE`.
- Replaced direct `cat()` diagnostics with suppressible `message()`/`warning()` calls.
- Removed the top-level `.Random.seed` assignment from `R/Stem-02-internal.R`.
- Added a package-specific vignette under `vignettes/STEM.Rmd`.
- Added an updated function map under `inst/extdata/`.

## To check locally

```r
devtools::document()
devtools::check()
rcmdcheck::rcmdcheck(args = c("--as-cran"), error_on = "warning")
```
