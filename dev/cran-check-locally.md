# Running the CRAN checks locally

Notes for the maintainer, not for CRAN (the comments for CRAN are in
`cran-comments.md`).

The repository lives inside a OneDrive-synced folder. `R CMD check` creates
thousands of small files, and the sync client stalls the run - the check hangs
at "checking package dependencies" for as long as it is left there. Run the
build and the check from a directory that is not synced:

```
cd <a folder that is not synced>
R CMD build <path to the working copy of this repository>
R CMD check --as-cran Stem_<version>.tar.gz
```

The check needs, in the library of the version of R that runs it:

- the packages of `Imports` and `Suggests` (`testthat` included: without it
  the tests are not run);
- `V8`, for the check of the math in the HTML manual (without it the check
  reports a NOTE that it skipped it);
- pandoc on the path, to build the vignettes and to check `README.md` and
  `NEWS.md`. If it is not there, use the copy that ships with RStudio:

```
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:$PATH"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
```

On R 4.6.1 (Windows 11) the build takes about two minutes and the check about
seven: four for the `\donttest` examples, one for the tests, one for the
vignettes. `spelling::spell_check_package()` and `urlchecker::url_check()`
check the spelling of the documentation and its URLs.
