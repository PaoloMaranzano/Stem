## ===========================================================================
## SC-STEM fuel application: MAIN SCRIPT, SETTING 2
##
##   Y  an independent pump (pompa bianca)
##   X  its neighbour among the pumps of the major brands (Eni, Api-Ip, Esso, Q8,
##      Tamoil)
##
## Question: do the major-brand pumps lead the independent pumps near them, on
## the same fuel or across fuels, and does the relation differ across spatial
## regimes of the city? This is also the case of the SMALL SAMPLES: 6 to 71
## independents per city; with fewer than 2m = 10 (Florence, Venice) the model
## is the pooled STEM, elsewhere K_max = min(3, floor(n / m)). Read with the
## number of regimes chosen, the stability of the partition under the
## bootstrap (B_valid, the bootstrap groups) and the ranges at a limit
## (regimes.csv, theta_bound). From the same fits: the pass-through of the excise
## duties by regime (the coefficients of the fiscal pulses).
##
## The model, for pump Y in regime g, fuel f (gasoline or diesel) and week t:
##
##   y*_t = a_g + sum_l b_gl y*_{t-l}        own lags, fuel f
##              + sum_l c_gl y*o_{t-l}       own lags, the other fuel
##              + sum_l d_gl x*_{t-l}        the neighbour X, fuel f
##              + sum_l e_gl x*o_{t-l}       the neighbour X, the other fuel
##              + pulses of the excise duties
##              + latent AR(1) of the regime + spatially correlated error
##
## with y* the weekly relative prices (cents per litre; Sunday price minus the
## national mean, minus the pump mean). The response keeps its missing weeks,
## which the Kalman filter handles; the lags come from the completed series
## (fuels-data.R, stage 3).
##
## The model is estimated twice on the same regimes. Without penalty
## (lambda = 0): the regimes are selected and every test is made on this fit.
## With a ridge on the four blocks of lags (lambda*): the lag coefficients are
## estimated with less variance, and the response is completed; the intercept
## and the pulses are never penalized. Why this order: the ridge would make the
## regimes look more alike, and it is not needed to find them (collinearity
## spoils the single coefficients, not the fit), so the regimes are selected
## without it and the penalty is tuned on the regimes found; the ridge
## corrected for its bias is the estimate at lambda = 0 itself, so the
## inference is made there (the proposition on the ridge in Section 2 of the
## paper).
##
## Four stages, run in this order (--stage=), each on the output of the one
## before it:
##
##   lags    the pre-analysis of the lags, pump by pump: ACF, PACF, cross-
##           correlations with the other three series, lag plots. It led to
##           the lags LAGS below.                              (no model)
##   grid    for each city and fuel, at lambda = 0: the grid of SC-STEM fits
##           over (K, phi) and the two-step rule of SCSTEM_Select() on the
##           BIC, which gives K*, phi* and the partition.
##   lambda  on that partition, held fixed: one fit for every ridge penalty in
##           `lambdas`; lambda* minimizes the AIC (`lambda_ic`).
##   final   the refit-with-clustering bootstrap of the model at lambda = 0,
##           of the model at lambda* and of the pooled model; the tests, all at
##           lambda = 0 (bootstrap Wald by regime, likelihood ratio on the
##           partition, pump by pump F and HC1, Dumitrescu-Hurlin); the lag
##           coefficients at lambda* beside those at lambda = 0; the completed
##           response at both; the diagnostics, the pulses, the maps.
##   all     grid, lambda and final one after the other, in one process.
##   summary the tables of grid, lambda and final, rebuilt from the files of
##           every model of the case in the folder; nothing fitted.
##   dry     the plan: the models and the fits of every stage; nothing fitted.
##
##     Rscript fuels-main-setting2.R --stage=lags
##     Rscript fuels-main-setting2.R --stage=grid
##     Rscript fuels-main-setting2.R --stage=lambda
##     Rscript fuels-main-setting2.R --stage=final
##     Rscript fuels-main-setting2.R --stage=all
##     Rscript fuels-main-setting2.R --stage=grid --job=1/3      one of three processes
##     Rscript fuels-main-setting2.R --stage=grid --cities=VE    one city only
##     Rscript fuels-main-setting2.R --stage=grid --xdef=rs      X = the mean within r*
##     Rscript fuels-main-setting2.R --stage=final --boot_cores=8  every bootstrap on 8 processes
##
## Several machines. Each machine runs the cities given to it, with as many
## processes for every bootstrap as it has cores, for instance
##
##     Rscript fuels-main-setting2.R --stage=all --cities=RM --boot_cores=7            machine 1
##     Rscript fuels-main-setting2.R --stage=all --cities=MI,TO,NA --boot_cores=3      machine 2
##     Rscript fuels-main-setting2.R --stage=all --cities=RM --fuels=d --boot_cores=8  one fuel only
##
## Every model writes its own files (grid/, lambda/, final/<model>/ and the
## cache), and the seed of a model depends on its city and fuel only, so how
## the cities are shared among the machines does not change the results. Once
## every machine has finished, copy the folder <case>/<xdef>/ of each into one
## (the files of the models have different names; the tables, which do not,
## are rebuilt) and run
##
##     Rscript fuels-main-setting2.R --stage=summary
##
## which rebuilds grid/selected.csv, lambda/selected.csv and final/summary.csv
## for all the cities. Never run the same city and fuel on two machines at once.
##
## Input:  <root>/setting2-Yindependent-Xmajor/data.RData   (fuels-data.R)
## Output: <root>/setting2-Yindependent-Xmajor/<xdef>/     lags/, grid/, lambda/, final/
##
## This is part of the REPLICATION MATERIAL of the paper, not of the Stem
## package.
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
source(file.path(HERE, "fuels-functions.R"))


## ===========================================================================
## THE SETUP. This is the part meant to be edited.
## ===========================================================================
SETTING <- 2L
CFG <- fu_config(list(
  root        = file.path(HERE, "fuels"),
  stage       = "lags",            # "lags", "grid", "lambda", "final", "all", "summary" or "dry"
  ## The neighbour X: "nn" the nearest pump of the X set; "rs" the mean of the X
  ## set within r* (the nearest when r* holds none). "nn" is the main analysis
  ## in every setting, "rs" a check here, where r* often holds several X and
  ## their mean is a different object from the nearest pump.
  xdef        = "nn",
  cities      = "",                # empty: every city of the case
  fuels       = c("g", "d"),       # the fuel of Y: g gasoline, d diesel
  m           = 5L,                # the smallest regime: K_max = min(k_cap, floor(n / m))
  k_cap       = 3L,
  knn         = 5L,                # neighbours of the graph of the Potts penalty
  phi_grid    = c(0, 0.025, 0.05, 0.1, 0.2, 0.5, 1),
  band        = c(0.025, 0.2),     # the band of the two-step rule
  ## the ridge penalties of stage "lambda" and the criterion choosing among
  ## them ("AIC"; "BIC" and "KIC" are also recorded)
  lambdas     = c(0, 0.001, 0.003, 0.01, 0.02, 0.03, 0.05, 0.1, 0.3, 1),
  lambda_ic   = "AIC",
  B           = 100L,              # bootstrap draws
  ## the R processes every bootstrap runs on (SCSTEM_Bootstrap(cores =)): 1, the
  ## default, runs the draws one after the other; the draws do not depend on it
  boot_cores  = 1L,
  boot_pooled = TRUE,              # bootstrap the pooled model too, when K* > 1
  level       = 0.95,
  seed        = 20261005L,
  job         = "1/1"              # "i/N": this process runs its share of the models
))

## The lags of every model, the same in the three settings, the eleven cities
## and the two fuels: the last four weeks of each of the four blocks
## (FU_BLOCKS: own, own_other, nb, nb_other), 16 lags, as in a VAR(4). From the
## pre-analysis of stage "lags" (2026-10-09), medians over the cities of the
## three settings: the PACF of the own price is beyond the band for 100% of the
## pumps at lag 1, 56-72% at lag 2, 26-41% at lag 3, 14-22% at lag 4 and 0-17%
## at lags 5 to 8; the order of the AR chosen by the AIC pump by pump has
## median 3-4 and quartiles 2 and 5-6. The own dynamics reach lag 4, with a weak
## tail to lag 8. Lag 52 is not needed: its PACF is beyond the band for at most
## 1% of the pumps, the ACF at lag 52 (14-30%) being the persistence of the
## series, not an annual season. The other three blocks take the same four
## lags, the convention of the Granger tests. Lags may be non-contiguous;
## integer(0) leaves a block out. The 52 weeks of 2021 are a presample: the lags
## cost no week of the window, and every model is fitted on the weeks of the
## window (2022-01 to 2026-06).
LAGS <- list(own = 1:4, own_other = 1:4, nb = 1:4, nb_other = 1:4)

## The ridge penalty of stage "final": NULL takes, city by city and fuel by
## fuel, lambda* of stage "lambda"; a number imposes that penalty on every
## model (it is then fitted on the partition of stage "grid"), for instance
##   FIT_LAMBDA <- 0.01
FIT_LAMBDA <- NULL


## ===========================================================================
## The data of the case, and the folders of its outputs
## ===========================================================================
CASE <- fu_case_dir(normalizePath(CFG$root[1], winslash = "/", mustWork = FALSE), SETTING)
if (!file.exists(file.path(CASE, "data.RData")))
  stop("no data for this case: run fuels-data.R first (", file.path(CASE, "data.RData"), ")", call. = FALSE)
load(file.path(CASE, "data.RData"))                    # fuels_data
XDEF <- match.arg(CFG$xdef[1], c("nn", "rs"))
OUT <- file.path(CASE, XDEF)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
CITIES <- if (any(nzchar(CFG$cities))) CFG$cities else unique(fuels_data$city)
STAGE <- CFG$stage[1]
if (!STAGE %in% c("lags", "grid", "lambda", "final", "all", "summary", "dry"))
  stop("unknown stage: ", STAGE, "; the stages are lags, grid, lambda, final, all, summary and dry",
       call. = FALSE)
options(width = 180)
message(sprintf("setting %d (%s), X = %s, stage %s, %d cities", SETTING, FU_CASES$folder[SETTING],
                XDEF, STAGE, length(CITIES)))
city_data <- function(cc) fuels_data[fuels_data$city == cc, , drop = FALSE]


## ===========================================================================
## Stage "lags": the pre-analysis of the lags, pump by pump
## ===========================================================================
## For each city and fuel: the median and the quartiles of the ACF, the PACF
## and the cross-correlations at lags 1 to 8 and 52, the share of pumps beyond
## the band 2 / sqrt(T), and one page of plots. Read lags_summary.csv first: a
## lag carries information when its share clearly exceeds the 5% expected by
## chance, in most cities. The lags LAGS were chosen on the ACF and the PACF;
## the cross-correlations of the persistent price series exceed the band at
## almost every lag and were not used.
if (STAGE == "lags") {
  LDIR <- file.path(OUT, "lags"); dir.create(LDIR, showWarnings = FALSE)
  stats_all <- list()
  for (fu in CFG$fuels) {
    grDevices::cairo_pdf(file.path(LDIR, sprintf("lags_%s.pdf", FU_FUEL[[fu]])), width = 16, height = 6.5,
                         onefile = TRUE)
    for (cc in CITIES) {
      dfc <- city_data(cc)
      st <- fu_lag_stats(dfc, fu, XDEF)
      stats_all[[length(stats_all) + 1]] <- st
      fu_lag_plots(dfc, fu, XDEF, st)
      message("  lags: ", FU_CITY[[cc]], ", ", FU_FUEL[[fu]])
    }
    invisible(grDevices::dev.off())
  }
  st <- do.call(rbind, stats_all)
  utils::write.csv(st, file.path(LDIR, "lags_by_city.csv"), row.names = FALSE)
  ## across the cities: the median share of pumps beyond the band, by measure,
  ## fuel and lag
  sm <- stats::aggregate(share_sig ~ measure + fuel + lag, data = st, FUN = stats::median)
  sm$share_sig <- round(100 * sm$share_sig, 1)
  sm <- stats::reshape(sm, idvar = c("measure", "fuel"), timevar = "lag", direction = "wide")
  names(sm) <- sub("^share_sig[.]", "lag_", names(sm))
  cat("\nSHARE OF PUMPS BEYOND THE BAND 2/sqrt(T), %, median over the cities\n")
  print(sm, row.names = FALSE)
  utils::write.csv(sm, file.path(LDIR, "lags_summary.csv"), row.names = FALSE)
  quit(save = "no", status = 0)
}


## ===========================================================================
## Stem, at the pinned commit, and the cache of the fits
## ===========================================================================
if (!STAGE %in% c("dry", "summary")) {
  SHA7 <- fu_require_stem()
  suppressPackageStartupMessages(library(Stem))
} else SHA7 <- substr(sub("^.*@", "", FU_STEM_REF), 1L, 7L)
CACHE <- file.path(OUT, "cache", SHA7)
## The lags enter the names of the cached fits, so that a fit made with other
## lags is never reused: one group per block, its runs of consecutive lags
## written "atb" (a single lag "a"), so that LAGS = 1:4 in every block gives
## "L1t4-1t4-1t4-1t4" and c(1, 2, 52) would give "1t2.52". Every cached value
## also carries the settings it was made with (the grid of phi, the
## penalties, the partition selected before it, B): when they no longer match,
## fu_cached() runs the stage again.
lag_runs <- function(l) {
  if (!length(l)) return("0")
  l <- sort(unique(l)); r <- split(l, cumsum(c(1, diff(l) != 1)))
  paste(vapply(r, function(x) if (length(x) == 1L) as.character(x) else paste0(x[1], "t", x[length(x)]), ""),
        collapse = ".")
}
LAG_KEY <- paste0("L", paste(vapply(FU_BLOCKS, function(b) lag_runs(LAGS[[b]]), ""), collapse = "-"))
with_settings <- function(val, s) { attr(val, "fu_settings") <- s; val }
same_settings <- function(s) function(val) identical(attr(val, "fu_settings"), s)


## ===========================================================================
## The models: one per city and fuel
## ===========================================================================
model_list <- function(cities) {
  out <- list()
  for (cc in cities) for (fu in CFG$fuels) {
    n <- length(unique(city_data(cc)$y_site))
    out[[length(out) + 1]] <- list(id = sprintf("%s_%s", cc, FU_FUEL[[fu]]), city = cc, fuel = fu, n = n,
                                   kmax = fu_kmax(n, CFG$m[1], CFG$k_cap[1]))
  }
  out
}
## the models this process runs (--cities, --fuels, --job), and every model of
## the case, from whose files the tables of the stages are rebuilt
MODELS <- model_list(CITIES)
MODELS_ALL <- model_list(unique(fuels_data$city))
sizes <- vapply(MODELS, `[[`, 1, "n")
MINE <- fu_my_share(sizes, CFG$job[1])
## the seed of a model depends on its city and fuel only, so that a run on some
## of the cities reproduces the full run
seed_of <- function(M) CFG$seed[1] + 10L * match(M$city, names(FU_CITY)) + match(M$fuel, names(FU_FUEL))
design_of <- function(M) fu_design(city_data(M$city), M$fuel, XDEF, LAGS)
cache_file <- function(name) file.path(CACHE, paste0(name, ".rds"))
## a stage that needs the output of the one before it stops instead of running
## that one silently
need <- function(name, stage) if (!file.exists(cache_file(name)))
  stop("no output of stage \"", stage, "\" for ", sub("_L.*$", "", name), ": run --stage=", stage, " first",
       call. = FALSE)

## Stage "grid", one model: the grid of (K, phi) at lambda = 0 and the
## two-step rule on the BIC.
grid_name <- function(M) sprintf("%s_%s_grid", M$id, LAG_KEY)
grid_of <- function(M, des) {
  s <- list(K_grid = seq_len(M$kmax), phi_grid = CFG$phi_grid, m = CFG$m[1], knn = CFG$knn[1], seed = seed_of(M))
  ic <- fu_cached(CACHE, grid_name(M), with_settings(
    Stem::SCSTEM_Infocrit(fu_model(des), K_grid = seq_len(M$kmax), phi_grid = CFG$phi_grid,
                          knn = min(CFG$knn[1], des$d - 1L), distance = "geo",
                          min_cluster_size = CFG$m[1], seed = seed_of(M), verbose = TRUE,
                          alpha = 0, lambda = 0, control = Stem::STEM_control()), s),
    valid = same_settings(s))
  list(ic = ic, sel = Stem::SCSTEM_Select(ic, band = CFG$band, criterion = "BIC"))
}

## Stage "lambda", one model: the fits along lambda on the partition of the
## grid, lambda* by the criterion lambda_ic.
lambda_name <- function(M, fixed = FALSE) sprintf("%s_%s_lambda%s", M$id, LAG_KEY, if (fixed) "_fixed" else "")
lambda_of <- function(M, des, sel, lambdas = CFG$lambdas, fixed = FALSE) {
  s <- list(lambdas = lambdas, ic = CFG$lambda_ic[1], K = sel$K_selected, phi = sel$phi_selected,
            partition = sel$fit$group)
  fu_cached(CACHE, lambda_name(M, fixed),
    with_settings(fu_lambda(des, sel, lambdas, CFG, seed_of(M), CFG$lambda_ic[1]), s), valid = same_settings(s))
}


## ===========================================================================
## Stage "dry": the plan
## ===========================================================================
plan <- do.call(rbind, lapply(seq_along(MODELS), function(k) with(MODELS[[k]], data.frame(
  model = id, city = FU_CITY[[city]], fuel = FU_FUEL[[fuel]], pumps = n, K_max = kmax,
  grid_fits = 1L + (kmax - 1L) * length(CFG$phi_grid),
  lambda_fits = length(CFG$lambdas),
  boot_refits_max = CFG$B[1] * (2L + (kmax > 1L && isTRUE(CFG$boot_pooled[1]))),
  this_job = k %in% MINE))))
if (STAGE == "dry") {
  cat(sprintf("\nTHE PLAN: %d models (%d in this job), X = %s, B = %d\n", nrow(plan), length(MINE), XDEF, CFG$B[1]))
  print(plan[order(-plan$pumps), ], row.names = FALSE)
  cat(sprintf(paste0("stage grid: %d SC-STEM fits; stage lambda: %d fits on a fixed partition; ",
                     "stage final: up to %d bootstrap refits\n"),
              sum(plan$grid_fits), sum(plan$lambda_fits), sum(plan$boot_refits_max)))
  utils::write.csv(plan, file.path(OUT, "plan.csv"), row.names = FALSE)
  quit(save = "no", status = 0)
}


## ===========================================================================
## The tables of the stages
## ===========================================================================
## Every model writes its own files. The table of a stage is rebuilt from the
## files of every model of the case present in this folder, whichever process
## or machine ran it: once the folders of several machines are put together,
## --stage=summary rebuilds all of them without fitting anything.
GDIR <- file.path(OUT, "grid"); LDIR <- file.path(OUT, "lambda"); FDIR <- file.path(OUT, "final")
for (dd in c(GDIR, LDIR, FDIR)) dir.create(dd, showWarnings = FALSE)
stack_rows <- function(path_of, keep = function(r) TRUE)
  do.call(rbind, lapply(MODELS_ALL, function(M) {
    f <- path_of(M)
    if (file.exists(f)) { r <- readRDS(f); if (keep(r)) r }
  }))
table_grid <- function() {
  tab <- stack_rows(function(M) file.path(GDIR, paste0(M$id, "_selected.rds")))
  if (is.null(tab)) return(invisible(NULL))
  utils::write.csv(tab, file.path(GDIR, "selected.csv"), row.names = FALSE)
  cat(sprintf("\nSELECTED AT lambda = 0 (two-step rule, BIC): %d of %d models\n", nrow(tab), length(MODELS_ALL)))
  print(tab, row.names = FALSE)
}
table_lambda <- function() {
  tab <- stack_rows(function(M) file.path(LDIR, paste0(M$id, "_selected.rds")))
  if (is.null(tab)) return(invisible(NULL))
  utils::write.csv(tab, file.path(LDIR, "selected.csv"), row.names = FALSE)
  cat(sprintf("\nlambda* BY THE %s, ON THE PARTITION OF STAGE grid: %d of %d models\n", CFG$lambda_ic[1],
              nrow(tab), length(MODELS_ALL)))
  print(tab, row.names = FALSE)
}
table_final <- function() {
  ## the models finished with this Stem only
  tab <- stack_rows(function(M) file.path(FDIR, M$id, "summary_row.rds"), function(r) identical(r$stem, SHA7))
  if (is.null(tab)) return(invisible(NULL))
  utils::write.csv(tab, file.path(FDIR, "summary.csv"), row.names = FALSE)
  cat(sprintf("\nSUMMARY: %d of %d models finished\n", nrow(tab), length(MODELS_ALL)))
  print(tab, row.names = FALSE)
}


## ===========================================================================
## Stage "grid": the regimes, at lambda = 0
## ===========================================================================
run_grid <- function() {
  for (k in MINE) {
    M <- MODELS[[k]]
    message(sprintf("\n== grid: %s, %s (%d pumps, K_max = %d)", FU_CITY[[M$city]], FU_FUEL[[M$fuel]], M$n, M$kmax))
    res <- tryCatch({
      des <- design_of(M); g <- grid_of(M, des)
      utils::write.csv(g$ic$table, file.path(GDIR, paste0(M$id, "_grid.csv")), row.names = FALSE)
      row <- data.frame(model = M$id, city = FU_CITY[[M$city]], fuel = FU_FUEL[[M$fuel]], pumps = M$n,
                        weeks = des$Tn, K_max = M$kmax, K = g$sel$K_selected, phi = g$sel$phi_selected,
                        BIC = g$sel$selected_row$BIC[1],
                        sizes = paste(tabulate(g$sel$fit$group), collapse = "/"),
                        min_grid = fu_minutes(g$ic), stringsAsFactors = FALSE)
      saveRDS(row, file.path(GDIR, paste0(M$id, "_selected.rds")))
      message(sprintf("   K* = %d, phi* = %s, regimes of %s pumps", row$K, format(row$phi), row$sizes))
    }, error = function(e) e)
    if (inherits(res, "error")) message("   FAILED: ", conditionMessage(res))
  }
  table_grid()
}


## ===========================================================================
## Stage "lambda": the ridge penalty, on the partition selected at lambda = 0
## ===========================================================================
run_lambda <- function() {
  for (k in MINE) {
    M <- MODELS[[k]]
    message(sprintf("\n== lambda: %s, %s (%d pumps)", FU_CITY[[M$city]], FU_FUEL[[M$fuel]], M$n))
    res <- tryCatch({
      need(grid_name(M), "grid")
      des <- design_of(M); g <- grid_of(M, des)
      lam <- lambda_of(M, des, g$sel)
      utils::write.csv(lam$table, file.path(LDIR, paste0(M$id, "_lambda.csv")), row.names = FALSE)
      row <- data.frame(model = M$id, city = FU_CITY[[M$city]], fuel = FU_FUEL[[M$fuel]], pumps = M$n,
                        K = g$sel$K_selected, phi = g$sel$phi_selected, criterion = CFG$lambda_ic[1],
                        lambda = lam$lambda,
                        df_lambda0 = lam$table$df[lam$table$lambda == 0][1],
                        df_lambda_star = lam$table$df[lam$table$selected][1],
                        min_lambda = fu_minutes(lam), stringsAsFactors = FALSE)
      saveRDS(row, file.path(LDIR, paste0(M$id, "_selected.rds")))
      message(sprintf("   lambda* = %s (%s)", format(row$lambda), row$criterion))
    }, error = function(e) e)
    if (inherits(res, "error")) message("   FAILED: ", conditionMessage(res))
  }
  table_lambda()
}


## ===========================================================================
## Stage "final": bootstrap, tests and tables of one model
## ===========================================================================
final_one <- function(M) {
  message(sprintf("\n== final: %s, %s (%d pumps, K_max = %d)", FU_CITY[[M$city]], FU_FUEL[[M$fuel]], M$n, M$kmax))
  need(grid_name(M), "grid")
  if (is.null(FIT_LAMBDA)) need(lambda_name(M), "lambda")
  od <- file.path(FDIR, M$id); dir.create(od, showWarnings = FALSE)
  des <- design_of(M)

  ## 1. the regimes at lambda = 0 (stage grid), the pooled model, lambda*
  ##    (stage lambda, or FIT_LAMBDA) and the fit at lambda* on those regimes
  g <- grid_of(M, des); sel <- g$sel; K <- sel$K_selected
  pooled <- g$ic$fits[[which(g$ic$table$K == 1)[1]]]
  lam <- if (is.null(FIT_LAMBDA)) lambda_of(M, des, sel) else lambda_of(M, des, sel, FIT_LAMBDA, fixed = TRUE)
  fitL <- lam$fit                                      # NULL when lambda* = 0

  ## 2. the refit-with-clustering bootstrap: of the model at lambda = 0 (for
  ##    the inference), of the model at lambda* (the spread of the ridge
  ##    estimates) and of the pooled model (for the comparison). The fit at
  ##    lambda* was made on a fixed partition (max_iter = 0, kept in its
  ##    control); its draws are given the control of the fit at lambda = 0, so
  ##    that they re-estimate the partition at (K*, phi*) as those of the
  ##    model at lambda = 0 do.
  boot <- function(fit, name, s, ...) fu_cached(CACHE, sprintf("%s_%s_%s", M$id, LAG_KEY, name), with_settings({
    bt <- Stem::SCSTEM_Bootstrap(fit, B = CFG$B[1], seed = seed_of(M), verbose = TRUE,
                                 cores = CFG$boot_cores[1], ...)
    list(info = bt$info, B_valid = bt$B_valid, inf = Stem::SCSTEM_BootInference(bt, level = CFG$level[1]),
         groups = bt$groups)
  }, s), valid = same_settings(s))
  s0 <- list(B = CFG$B[1], level = CFG$level[1], K = K, phi = sel$phi_selected, partition = sel$fit$group)
  b0 <- boot(sel$fit, "boot_lambda0", s0)
  bL <- if (is.null(fitL)) NULL else
    boot(fitL, "boot_lambdastar", c(s0, lambda = lam$lambda), control = sel$fit$input_args$control)
  bP <- if (K == 1L) b0 else if (isTRUE(CFG$boot_pooled[1]))
    boot(pooled, "boot_pooled", list(B = CFG$B[1], level = CFG$level[1], K = 1L)) else NULL

  ## 3. the tests, all at lambda = 0: bootstrap Wald regime by regime, for the
  ##    selected and the pooled model; likelihood ratio on the selected
  ##    partition; pump by pump, F and HC1, with the Dumitrescu-Hurlin
  ##    combination
  tests <- fu_tests(des)
  wald <- function(fit, b, label) if (!is.null(b)) do.call(rbind, lapply(names(tests), function(h)
    cbind(model = label, fu_wald_boot(fit, b$inf, des, tests[[h]], h))))
  wd <- rbind(wald(sel$fit, b0, sprintf("selected, K = %d", K)),
              if (K > 1L) wald(pooled, bP, "pooled"))
  sP <- list(partition = sel$fit$group, K = K)
  lr <- fu_cached(CACHE, sprintf("%s_%s_lr", M$id, LAG_KEY),
                  with_settings(fu_lr_fixed(des, sel$fit$group, K, CFG, seed_of(M)), sP), valid = same_settings(sP))
  pt <- fu_pump_tests(des); dh <- fu_dh(pt)

  ## 4. the ridge: the lag coefficients at lambda* beside those at lambda = 0,
  ##    and the completed response of both fits
  rt <- fu_ridge_table(b0$inf, if (is.null(bL)) NULL else bL$inf, des)
  saveRDS(list(weeks = des$weeks, sites = des$meta$y_site, missing = is.na(des$z), lambda = lam$lambda,
               lambda0 = Stem::SCSTEM_Complete(sel$fit),
               lambda_star = if (is.null(fitL)) NULL else Stem::SCSTEM_Complete(fitL)),
          file.path(od, "completed_response.rds"))

  ## 5. diagnostics, the pass-through of the excise duties, regimes and brands
  dg <- fu_diagnostics(sel$fit, des)
  pul <- rbind(cbind(model = sprintf("selected, K = %d", K), fu_pulse_table(b0$inf, des)),
               if (K > 1L && !is.null(bP)) cbind(model = "pooled", fu_pulse_table(bP$inf, des)))
  br <- fu_brand_table(sel$fit, des)

  ## 6. the outputs of the model
  est_table <- function(inf) { e <- as.data.frame(inf$summary)
    num <- suppressWarnings(as.integer(sub("^beta", "", e$parameter)))
    e$covariate <- ifelse(is.na(num), e$parameter, des$names[num]); e }
  utils::write.csv(g$ic$table, file.path(od, "grid.csv"), row.names = FALSE)
  utils::write.csv(lam$table, file.path(od, "lambda.csv"), row.names = FALSE)
  utils::write.csv(wd, file.path(od, "wald_bootstrap.csv"), row.names = FALSE)
  utils::write.csv(lr, file.path(od, "lr_fixed_partition.csv"), row.names = FALSE)
  utils::write.csv(pt, file.path(od, "pump_tests.csv"), row.names = FALSE)
  utils::write.csv(dh, file.path(od, "dumitrescu_hurlin.csv"), row.names = FALSE)
  utils::write.csv(dg, file.path(od, "diagnostics.csv"), row.names = FALSE)
  utils::write.csv(pul, file.path(od, "excise_pulses.csv"), row.names = FALSE)
  utils::write.csv(br$table, file.path(od, "regimes_by_brand.csv"), row.names = FALSE)
  utils::write.csv(est_table(b0$inf), file.path(od, "estimates_lambda0.csv"), row.names = FALSE)
  if (!is.null(bL)) {
    utils::write.csv(est_table(bL$inf), file.path(od, "estimates_lambdastar.csv"), row.names = FALSE)
    utils::write.csv(rt, file.path(od, "ridge_vs_lambda0.csv"), row.names = FALSE)
  }
  utils::write.csv(data.frame(des$meta, regime = sel$fit$group,
                              theta_bound = sel$fit$theta_bound[sel$fit$group]),
                   file.path(od, "regimes.csv"), row.names = FALSE)
  fu_map(des, sel$fit, file.path(od, "map.pdf"),
         sprintf("%s, %s (Y): K = %d, phi = %s", FU_CITY[[M$city]], FU_FUEL[[M$fuel]], K,
                 format(sel$phi_selected)))

  ## 7. one row of the summary of the case
  p_of <- function(h, model) { z <- wd[wd$block == h & startsWith(wd$model, model), ]
    if (!nrow(z)) NA_character_ else paste(format.pval(z$p_value, digits = 2), collapse = "; ") }
  row <- data.frame(model = M$id, city = FU_CITY[[M$city]], fuel = FU_FUEL[[M$fuel]], pumps = M$n,
                    weeks = des$Tn, K_max = M$kmax, K = K, phi = sel$phi_selected,
                    lambda = lam$lambda, lambda_from = if (is.null(FIT_LAMBDA)) CFG$lambda_ic[1] else "fixed",
                    B_valid = b0$B_valid, B_valid_ridge = if (is.null(bL)) NA_integer_ else bL$B_valid,
                    p_nb_both_regimes = p_of("nb_both", "selected"), p_nb_both_pooled = p_of("nb_both", "pooled"),
                    p_nb_regimes = p_of("nb", "selected"), p_nb_other_regimes = p_of("nb_other", "selected"),
                    p_own_other_regimes = p_of("own_other", "selected"),
                    pumps_reject_nb_both_F_pct = dh$reject_F_pct[dh$block == "nb_both"],
                    pumps_reject_nb_both_hc1_pct = dh$reject_hc1_pct[dh$block == "nb_both"],
                    ridge_se_ratio_median = if (is.null(rt)) NA_real_ else stats::median(rt$se_ratio),
                    ridge_shift_abs_median = if (is.null(rt)) NA_real_ else stats::median(abs(rt$shift)),
                    ari_regimes_brands = br$ari,
                    min_grid = fu_minutes(g$ic), min_lambda = fu_minutes(lam), min_boot = fu_minutes(b0),
                    stem = SHA7, stringsAsFactors = FALSE)
  saveRDS(row, file.path(od, "summary_row.rds"))
  message(sprintf("   K* = %d, phi* = %s, lambda* = %s | neighbour -> Y, p by regime: %s", K,
                  format(sel$phi_selected), format(lam$lambda), row$p_nb_both_regimes))
  invisible(row)
}

run_final <- function() {
  for (k in MINE) {
    res <- tryCatch(final_one(MODELS[[k]]), error = function(e) e)
    if (inherits(res, "error")) message("   FAILED: ", conditionMessage(res))
  }
  table_final()
}


## ===========================================================================
## Run the stage
## ===========================================================================
if (STAGE == "summary") { table_grid(); table_lambda(); table_final() }
if (STAGE %in% c("grid", "all")) run_grid()
if (STAGE %in% c("lambda", "all")) run_lambda()
if (STAGE %in% c("final", "all")) run_final()
cat("\noutputs in", OUT, "\n")
