## ===========================================================================
## Latent Spatial Regimes in Spatio-Temporal State-Space Models: Evidence on
## Local Price Leadership in Italian Fuel Markets
##
## Author:  Paolo Maranzano (University of Milano-Bicocca; Fondazione Eni Enrico
##          Mattei)
## Script:  AppFuels_01_RunSetting1.R, the run of setting 1 of the application
##          (Section 5 of the paper, Supplement C)
##
## What it does: fits the 22 SC-STEM models of setting 1 (eleven cities, two
## fuels of Y), in which Y is a station of a major brand and X the nearest
## independent station.
## It starts from the weekly data of data/, the dataset of the paper, and first
## builds the data of the setting (fu_setting_data() in AppFuels_00_Setup.R):
## the relative prices, their missing weeks completed for the lags by a pooled
## STEM model of every city and fuel, the pairs (Y, X) and the pulses. Then,
## for every model (fu_run_model() in AppFuels_00_Setup.R):
##   1. the design: the response (the relative prices of the stations Y over
##      the window, with their missing weeks) and the covariates (an
##      intercept, the four lags of each of the four blocks, the pulses of
##      the excise duties);
##   2. the regimes: SC-STEM fits over the grid of the number of regimes K
##      (1 to min(5, floor(n/5)) for n stations) and of the spatial penalty
##      phi, and the two-step rule on the BIC, which selects (K, phi) and the
##      partition; the fit with K = 1 is the pooled STEM model;
##   3. the refit-with-clustering bootstrap (B = 100) of the selected fit and
##      of the pooled model;
##   4. the tests: the bootstrap Wald tests of the lags of X, regime by regime
##      and in all the regimes at once, for both fits; the station-by-station
##      least-squares benchmark (F and HC1 tests) and its panel combination
##      (Dumitrescu-Hurlin);
##   5. the residual diagnostics of the selected fit;
##   6. the check of the local maxima of the selected fit;
##   7. the compact results of the model, read by AppFuels_04_Analysis.R.
## No table or figure of the paper is written here: AppFuels_04_Analysis.R
## draws them all from the results.
##
## Outputs, in <out>/setting1/<city>_<fuel>/ (one folder per model; <out> is
## output/runs/application by default):
##   model.csv              one row: K_max, the selected K and phi, the sizes of
##                          the regimes, the usable bootstrap draws, the
##                          stability of the partition, the minutes taken
##   grid.csv               the information criteria of every (K, phi) of the grid
##   regimes.csv            the stations Y, their neighbor X and their regime
##   estimates.csv          the estimates of the selected fit and of the pooled
##                          model, with bootstrap standard errors and intervals
##   draws.csv.gz           the aligned bootstrap draws of the regression
##                          coefficients of the selected fit
##   wald.csv               the bootstrap Wald tests (both fits)
##   pump_tests.csv         the station-by-station F and HC1 tests
##   dumitrescu_hurlin.csv  their panel combination
##   diagnostics.csv        the residual diagnostics, by regime
##   optima.csv, optima_coefficients.csv   the check of the local maxima
## and <out>/setting1/models.csv, the rows of model.csv of the setting; and,
## in <cache>/completion, the completed series of every city and fuel (a cache
## shared by the three settings). The results of the run of the paper are in
## results/application/setting1.
## ===========================================================================
##
## Run from a terminal, in the folder of this script,
##
##     Rscript AppFuels_01_RunSetting1.R
##     Rscript AppFuels_01_RunSetting1.R --cities=ME,FI --fuels=g --cores=4
##
## or open it in RStudio and press Source. Options (command line):
##   --data=<folder>     the weekly data (default: ../data)
##   --out=<folder>      the compact results (default:
##                       ../output/runs/application); AppFuels_04_Analysis.R
##                       reads the results shipped in results/application
##                       unless told otherwise (--results=)
##   --cache=<folder>    the completed series and the intermediate fits
##                       (default: ../output/runs/application/cache)
##   --cities=RM,MI,...  --fuels=g,d   a subset of the models
##   --cores=1           R processes the draws of a bootstrap are shared among;
##                       the draws do not depend on it
##
## Every model is written as soon as it is done, and a model whose results are
## already there is skipped: a stopped run resumes where it stopped. The seed
## of a model depends on its city and fuel only, so the results do not depend
## on which models are run together, nor on --cores, and the three settings can
## run at the same time on different machines.
##
## Requirements: R (>= 4.2) with data.table and geodist, and Stem at the release
## named in AppFuels_00_Setup.R, 2.0.0 (installed from GitHub when another
## version is installed). Time: on a first run, the completion of the series
## takes about
## an hour (most of it for Rome) and serves the three settings; in the run of
## the paper the 22 models took about 73 hours, with the bootstraps on seven
## cores (columns minutes_* of models.csv), 24 of them for Rome, gasoline; the
## grids of the largest cities (Rome, Milan, Turin) need some GB of memory.
## ===========================================================================

suppressPackageStartupMessages(library(data.table))
invisible(Sys.setlocale("LC_TIME", "C"))
HERE <- local({
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  f <- if (length(m)) {
    normalizePath(sub("^--file=", "", m[1]), winslash = "/")
  } else {
    NULL
  }
  if (is.null(f)) {
    for (i in rev(seq_len(sys.nframe()))) {
      of <- sys.frame(i)$ofile
      if (!is.null(of)) {
        f <- normalizePath(of, winslash = "/")
        break
      }
    }
  }
  if (is.null(f)) {
    normalizePath(getwd(), winslash = "/")
  } else {
    dirname(f)
  }
})
source(file.path(HERE, "AppFuels_00_Setup.R"))

SETTING <- 1L
OPT <- fu_options(list(
  data   = file.path(HERE, "..", "data"),
  out    = file.path(HERE, "..", "output", "runs", "application"),
  cache  = file.path(HERE, "..", "output", "runs", "application", "cache"),
  cities = names(FU_CITY),
  fuels  = names(FU_FUEL),
  cores  = 1L
))
SHA7 <- fu_require_stem()
options(width = 160)


## ===========================================================================
## The data of the setting, from the weekly data: relative prices, completed
## series (cached in <cache>/completion), pairs (Y, X), pulses
## ===========================================================================
CASES <- fu_setting_data(SETTING, OPT$data, file.path(OPT$cache, "completion"), OPT$cities)


## ===========================================================================
## The 22 models of the setting, the smaller first, and the table of the
## models (models.csv)
## ===========================================================================
fu_run_setting(SETTING, CASES, OPT, SHA7)
