## ===========================================================================
## APPLICATION, SCRIPT 2 OF 3: THE SC-STEM MODELS, THEIR BOOTSTRAP AND THE
## TESTS
##
## For every setting, city and fuel of Y (66 models), on the data written by
## 01-data-and-eda.R:
##
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
##   4. the tests: the bootstrap Wald tests of the lags of X (and of the other
##      hypotheses of fu_tests()), regime by regime and in all the regimes at
##      once, for both fits; the station-by-station least-squares benchmark
##      (F and HC1 tests) and its panel combination (Dumitrescu-Hurlin);
##   5. the residual diagnostics of the selected fit;
##   6. the check of the local maxima of the selected fit (fu_check_optima());
##   7. the compact results of the model, read by 03-results.R.
##
## Run from a terminal, in the folder of this script,
##
##     Rscript 02-models.R
##     Rscript 02-models.R --settings=2 --cities=ME,FI --fuels=g --cores=4
##
## or open it in RStudio and press Source. Options (command line):
##   --cases=<file>      the data of the settings (default:
##                       ../output/application/data/cases.rds)
##   --out=<folder>      the compact results (default:
##                       ../output/application/results); 03-results.R reads
##                       the results shipped with the package unless told
##                       otherwise (--results=)
##   --cache=<folder>    the intermediate fits (default:
##                       ../output/application/cache)
##   --settings=1,2,3    --cities=RM,MI,...   --fuels=g,d   a subset of the models
##   --cores=1           R processes the draws of a bootstrap are shared among;
##                       the draws do not depend on it
##
## Every model is written as soon as it is done, and a model whose results are
## already there is skipped: a stopped run resumes where it stopped. Within a
## model, the grid, the bootstraps and the check of the local maxima are cached
## (<cache>), so a model stopped halfway resumes from its last completed step.
## The seed of a model depends on its city and fuel only (fu_seed()), so the
## results do not depend on which models are run together, nor on --cores.
##
## Outputs, in <out>/setting<s>/<city>_<fuel>/ (one folder per model):
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
## and <out>/setting<s>/models.csv, the rows of model.csv of the setting.
##
## Requirements: Stem at the commit pinned in functions.R. Time: the whole run
## took about 70 hours on a desktop computer with the bootstraps on seven
## cores; the largest models (Rome, Milan and Turin in settings 1 and 3) take
## several hours each, and their grid needs some GB of memory. The models of
## setting 2 take minutes to an hour each.
## ===========================================================================

suppressPackageStartupMessages(library(data.table))
invisible(Sys.setlocale("LC_TIME", "C"))
HERE <- local({
  a <- commandArgs(trailingOnly = FALSE); m <- grep("^--file=", a, value = TRUE)
  f <- if (length(m)) normalizePath(sub("^--file=", "", m[1]), winslash = "/") else NULL
  if (is.null(f)) for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) { f <- normalizePath(of, winslash = "/"); break }
  }
  if (is.null(f)) normalizePath(getwd(), winslash = "/") else dirname(f)
})
source(file.path(HERE, "functions.R"))

OPT <- fu_options(list(
  cases    = file.path(HERE, "..", "output", "application", "data", "cases.rds"),
  out      = file.path(HERE, "..", "output", "application", "results"),
  cache    = file.path(HERE, "..", "output", "application", "cache"),
  settings = FU_SETTINGS$setting,
  cities   = names(FU_CITY),
  fuels    = names(FU_FUEL),
  cores    = 1L
))
if (!file.exists(OPT$cases))
  stop("no data of the settings: run 01-data-and-eda.R first (", OPT$cases, ")", call. = FALSE)
CASES <- readRDS(OPT$cases)$cases
SHA7 <- fu_require_stem()
options(width = 160)


## ===========================================================================
## One model
## ===========================================================================
run_model <- function(s, cc, fu) {
  id <- sprintf("%s_%s", cc, FU_FUEL[[fu]])
  od <- file.path(OPT$out, sprintf("setting%d", s), id)
  if (file.exists(file.path(od, "model.csv"))) { message("done already: setting ", s, ", ", id); return(invisible()) }
  cache <- file.path(OPT$cache, sprintf("setting%d", s))
  seed <- fu_seed(cc, fu)

  ## 1. The design
  dfc <- CASES[[as.character(s)]]
  des <- fu_design(dfc[dfc$city == cc, , drop = FALSE], fu)
  kmax <- fu_kmax(des$d)
  message(sprintf("\n== setting %d, %s, %s: %d stations Y, %d weeks, K_max = %d", s, FU_CITY[[cc]],
                  FU_FUEL[[fu]], des$d, des$Tn, kmax))

  ## 2. The regimes: the grid of (K, phi) and the two-step rule on the BIC.
  ## Only the table of the grid, the selected fit and the pooled fit are kept
  ## (the grid of a large city takes some GB).
  g <- fu_cached(cache, paste0(id, "_grid"), {
    ic <- Stem::SCSTEM_Infocrit(fu_model(des), K_grid = seq_len(kmax), phi_grid = FU_PHI,
                                knn = min(FU_KNN, des$d - 1L), distance = "geo", min_cluster_size = FU_M,
                                seed = seed, verbose = TRUE, control = Stem::STEM_control())
    sel <- Stem::SCSTEM_Select(ic, band = FU_BAND, criterion = "BIC")
    list(table = ic$table, K = sel$K_selected, phi = sel$phi_selected, fit = sel$fit,
         pooled = ic$fits[[which(ic$table$K == 1)[1]]])
  })
  K <- g$K
  message(sprintf("   K = %d, phi = %s, regimes of %s stations", K, format(g$phi),
                  paste(tabulate(g$fit$group, K), collapse = "/")))

  ## 3. The refit-with-clustering bootstrap: each draw simulates the data from
  ## the fitted model, regime by regime, and re-estimates SC-STEM at the
  ## selected (K, phi), partition included; the regimes of a draw are matched
  ## to the estimated ones (Stem::SCSTEM_BootInference()). The pooled model is
  ## bootstrapped in the same way, with the same seed.
  boot <- function(fit, name) fu_cached(cache, paste0(id, "_", name), {
    bt <- Stem::SCSTEM_Bootstrap(fit, B = FU_B, seed = seed, verbose = TRUE, cores = OPT$cores)
    list(B_valid = bt$B_valid, inf = Stem::SCSTEM_BootInference(bt, level = FU_LEVEL))
  })
  b0 <- boot(g$fit, "boot")
  bP <- if (K == 1L) b0 else boot(g$pooled, "boot_pooled")

  ## 4. The tests
  wald <- rbind(fu_wald(g$fit$phi_hat, g$fit$group, b0$inf, des, "selected"),
                if (K > 1L) fu_wald(g$pooled$phi_hat, g$pooled$group, bP$inf, des, "pooled"))
  pt <- fu_pump_tests(des)
  dh <- fu_dh(pt)

  ## 5. The residual diagnostics
  dg <- fu_diagnostics(g$fit, des)

  ## 6. The check of the local maxima (NULL for a pooled model)
  opt <- fu_cached(cache, paste0(id, "_optima"), fu_check_optima(g$fit, des, b0$inf, seed))

  ## 7. The compact results
  dir.create(od, recursive = TRUE, showWarnings = FALSE)
  w <- function(d, f) if (!is.null(d)) utils::write.csv(cbind(setting = s, model = id, d), file.path(od, f),
                                                         row.names = FALSE)
  w(g$table[, !vapply(g$table, is.list, TRUE)], "grid.csv")
  w(data.frame(des$meta, regime = g$fit$group, theta_bound = g$fit$theta_bound[g$fit$group]), "regimes.csv")
  w(rbind(fu_estimates(b0$inf, des$names, "selected"),
          if (K > 1L) fu_estimates(bP$inf, des$names, "pooled")), "estimates.csv")
  data.table::fwrite(cbind(setting = s, model = id, fu_draws(b0$inf, des$names)), file.path(od, "draws.csv.gz"))
  w(wald, "wald.csv"); w(pt, "pump_tests.csv"); w(dh, "dumitrescu_hurlin.csv"); w(dg, "diagnostics.csv")
  if (!is.null(opt)) { w(opt$regimes, "optima.csv"); w(opt$coefficients, "optima_coefficients.csv") }
  brands <- des$meta$y_brand
  row <- data.frame(
    setting = s, model = id, city = cc, fuel = FU_FUEL[[fu]], pumps = des$d, weeks = des$Tn,
    K_max = kmax, K = K, phi = g$phi, BIC = g$table$BIC[g$table$K == K & g$table$phi == g$phi][1],
    sizes = paste(tabulate(g$fit$group, K), collapse = "/"),
    B_valid = b0$B_valid, B_valid_pooled = bP$B_valid,
    ari_regimes_brands = if (length(unique(brands)) > 1 && K > 1) fu_ari(g$fit$group, brands) else NA_real_,
    ari_draws_median = stats::median(b0$inf$stability$ARI),
    minutes_grid = fu_minutes(g), minutes_boot = fu_minutes(b0),
    minutes_boot_pooled = if (K > 1L) fu_minutes(bP) else 0, minutes_optima = fu_minutes(opt),
    stem = SHA7, stringsAsFactors = FALSE)
  utils::write.csv(row, file.path(od, "model.csv"), row.names = FALSE)   # written last: the model is done
  message(sprintf("   done: X leads Y (8 lags, 5%%) in %d of %d regimes", sum(wald$fit == "selected" &
    wald$block == "nb_both" & wald$regime != "all" & wald$p_chisq < 0.05), K))
  invisible()
}


## ===========================================================================
## Every model: the smaller ones first
## ===========================================================================
todo <- do.call(rbind, lapply(OPT$settings, function(s) {
  d <- CASES[[as.character(s)]]
  n <- tapply(d$y_site, d$city, function(v) length(unique(v)))
  expand.grid(setting = s, city = intersect(names(FU_CITY), intersect(OPT$cities, names(n))),
              fuel = OPT$fuels, stringsAsFactors = FALSE)
}))
todo$n <- vapply(seq_len(nrow(todo)), function(i) {
  d <- CASES[[as.character(todo$setting[i])]]; length(unique(d$y_site[d$city == todo$city[i]])) }, 1L)
todo <- todo[order(todo$n, todo$setting), ]
message(sprintf("%d models; Stem %s; bootstrap on %d process(es)", nrow(todo), SHA7, OPT$cores))
for (i in seq_len(nrow(todo))) {
  r <- tryCatch(run_model(todo$setting[i], todo$city[i], todo$fuel[i]), error = function(e) e)
  if (inherits(r, "error"))
    message("   FAILED: setting ", todo$setting[i], ", ", todo$city[i], ", ", todo$fuel[i], ": ", conditionMessage(r))
}

## the table of the models of every setting, from the files of all its models
for (s in unique(todo$setting)) {
  f <- Sys.glob(file.path(OPT$out, sprintf("setting%d", s), "*", "model.csv"))
  if (length(f)) fu_write(do.call(rbind, lapply(f, utils::read.csv, stringsAsFactors = FALSE)),
                          file.path(OPT$out, sprintf("setting%d", s), "models.csv"))
}
message("done: results in ", OPT$out)
