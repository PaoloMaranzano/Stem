## ===========================================================================
## SC-STEM fuel application: MAIN SCRIPT, SETTING 3
##
##   Y  a pump of a major brand (Eni, Api-Ip, Esso, Q8, Tamoil)
##   X  its neighbour among the pumps of the OTHER major brands
##
## Question: is there price leadership among the major brands at the local
## scale, on the same fuel or across fuels? In the exploration Eni and Tamoil
## predicted their neighbours more than they were predicted by them, Esso, Q8
## and Api-Ip the reverse: is that heterogeneity spatial, regimes of the city,
## or a property of the brands? regimes_by_brand.csv and its Adjusted Rand
## Index (about 0 when the regimes are unrelated to the brands) answer it,
## together with the penalty phi selected. From the same fits: the
## pass-through of the excise
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
## national mean, minus the pump mean). The four blocks of lags carry a ridge
## penalty; the intercept and the pulses do not. The response keeps its missing
## weeks, which the Kalman filter handles; the lags come from the completed
## series (fuels-data.R, stage 3).
##
## Three stages, run in this order (--stage=):
##
##   lags    the pre-analysis of the lags, pump by pump: ACF, PACF, cross-
##           correlations with the other three series, lag plots. Read its
##           output, then edit LAG_CANDIDATES below.          (no model)
##   select  for each city and fuel, the pooled STEM model for every candidate
##           set of lags and every ridge penalty in `lambdas`: the BIC chooses
##           both.                                             (pooled fits)
##   fit     for each city and fuel, at the chosen lags and penalty: the grid of
##           SC-STEM fits over (K, phi), the two-step selection, the
##           refit-with-clustering bootstrap of the selected model and of the
##           pooled one, the tests, the diagnostics, the maps.
##   dry     the plan of stage "fit": the models and their cost, nothing fitted.
##
##     Rscript fuels-main-setting3.R --stage=lags
##     Rscript fuels-main-setting3.R --stage=select
##     Rscript fuels-main-setting3.R --stage=fit
##     Rscript fuels-main-setting3.R --stage=fit --job=1/3      one of three processes
##     Rscript fuels-main-setting3.R --stage=fit --xdef=rs      X = the mean within r*
##
## Input:  <root>/setting3-Ymajor-Xothermajor/data.RData   (fuels-data.R)
## Output: <root>/setting3-Ymajor-Xothermajor/<xdef>/     lags/, select/, fit/
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
SETTING <- 3L
CFG <- fu_config(list(
  root        = file.path(HERE, "fuels"),
  stage       = "lags",            # "lags", "select", "fit" or "dry"
  ## The neighbour X: "nn" the nearest pump of the X set; "rs" the mean of the X
  ## set within r* (the nearest when r* holds none). Proposal: "nn" as the main
  ## analysis in every setting; "rs" as a check here, where r* holds several X
  ## for most Y and their mean is a different object from the nearest pump.
  xdef        = "nn",
  cities      = "",                # empty: every city of the case
  fuels       = c("g", "d"),       # the fuel of Y: g gasoline, d diesel
  lambdas     = c(0, 0.03, 0.1, 0.3, 1),   # ridge penalties tried in stage "select"
  m           = 5L,                # the smallest regime: K_max = min(k_cap, floor(n / m))
  k_cap       = 3L,
  knn         = 5L,                # neighbours of the graph of the Potts penalty
  phi_grid    = c(0, 0.025, 0.05, 0.1, 0.2, 0.5, 1),
  band        = c(0.025, 0.2),     # the band of the two-step rule
  B           = 100L,              # bootstrap draws
  boot_pooled = TRUE,              # bootstrap the pooled model too, when K > 1 is selected
  level       = 0.95,
  seed        = 20261005L,
  job         = "1/1"              # "i/N": this process runs its share of the models
))

## The candidate sets of lags of stage "select", one list per candidate with
## the lags of each of the four blocks (FU_BLOCKS: own, own_other, nb,
## nb_other). Lags may be non-contiguous; integer(0) leaves a block out. EDIT
## after reading the output of stage "lags". The 52 weeks of 2021 are a
## presample: lags up to 52 cost no week of the window, and every candidate is
## fitted on the weeks of the window (2022-01 to 2026-06), so their BIC are
## comparable.
LAG_CANDIDATES <- list(
  short = list(own = 1:2, own_other = 1:2, nb = 1:2, nb_other = 1:2),
  own4  = list(own = 1:4, own_other = 1:2, nb = 1:2, nb_other = 1:2),
  all4  = list(own = 1:4, own_other = 1:4, nb = 1:4, nb_other = 1:4)
)

## The lags and the penalty of stage "fit": NULL takes, city by city and fuel
## by fuel, those chosen by the BIC in stage "select"; otherwise they are
## imposed on every model, for instance
##   FIT_FIXED <- list(lags = LAG_CANDIDATES$own4, lambda = 0.1)
FIT_FIXED <- NULL


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
## lag worth a candidate is one whose share clearly exceeds the 5% expected by
## chance, in most cities.
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
## Stem, at the pinned commit (stages "select" and "fit")
## ===========================================================================
if (STAGE != "dry") {
  SHA7 <- fu_require_stem()
  suppressPackageStartupMessages(library(Stem))
} else SHA7 <- substr(sub("^.*@", "", FU_STEM_REF), 1L, 7L)
CACHE <- file.path(OUT, "cache", SHA7)


## ===========================================================================
## Stage "select": the lags and the ridge penalty, by the BIC of the pooled model
## ===========================================================================
if (STAGE == "select") {
  SDIR <- file.path(OUT, "select"); dir.create(SDIR, showWarnings = FALSE)
  units <- expand.grid(city = CITIES, fuel = CFG$fuels, stringsAsFactors = FALSE)
  sizes <- vapply(units$city, function(cc) length(unique(city_data(cc)$y_site)), 1L)
  res <- list()
  for (k in fu_my_share(sizes, CFG$job[1])) {
    cc <- units$city[k]; fu <- units$fuel[k]
    message(sprintf("\n== select: %s, %s (%d pumps)", FU_CITY[[cc]], FU_FUEL[[fu]], sizes[k]))
    res[[length(res) + 1]] <- fu_cached(CACHE, sprintf("select_%s_%s", cc, fu),
      fu_select(city_data(cc), fu, XDEF, LAG_CANDIDATES, CFG$lambdas))
  }
  ## every process rewrites the table with all the units finished so far
  done <- lapply(seq_len(nrow(units)), function(k) {
    f <- file.path(CACHE, sprintf("select_%s_%s.rds", units$city[k], units$fuel[k]))
    if (file.exists(f)) readRDS(f)
  })
  tab <- do.call(rbind, done)
  utils::write.csv(tab, file.path(SDIR, "select_all.csv"), row.names = FALSE)
  chosen <- tab[tab$selected, c("city", "fuel", "candidate", "lambda", "BIC", "df", "t_start")]
  cat("\nCHOSEN BY THE BIC: lags and ridge penalty, by city and fuel\n"); print(chosen, row.names = FALSE)
  utils::write.csv(chosen, file.path(SDIR, "selected.csv"), row.names = FALSE)
  quit(save = "no", status = 0)
}


## ===========================================================================
## Stage "fit": the plan
## ===========================================================================
## The lags and the penalty of every model: imposed (FIT_FIXED) or chosen in
## stage "select".
spec_of <- function(cc, fu) {
  if (!is.null(FIT_FIXED)) return(list(lags = FIT_FIXED$lags, lambda = FIT_FIXED$lambda, from = "fixed"))
  f <- file.path(OUT, "select", "selected.csv")
  if (!file.exists(f)) stop("no selected lags: run --stage=select first, or set FIT_FIXED", call. = FALSE)
  s <- utils::read.csv(f, stringsAsFactors = FALSE)
  s <- s[s$city == cc & s$fuel == FU_FUEL[[fu]], ]
  if (!nrow(s)) stop("no selected lags for ", cc, ", ", FU_FUEL[[fu]], call. = FALSE)
  list(lags = LAG_CANDIDATES[[s$candidate[1]]], lambda = s$lambda[1], from = s$candidate[1])
}
MODELS <- list()
for (cc in CITIES) for (fu in CFG$fuels) {
  n <- length(unique(city_data(cc)$y_site))
  MODELS[[length(MODELS) + 1]] <- list(id = sprintf("%s_%s", cc, FU_FUEL[[fu]]), city = cc, fuel = fu, n = n,
                                       kmax = fu_kmax(n, CFG$m[1], CFG$k_cap[1]))
}
sizes <- vapply(MODELS, `[[`, 1, "n")
MINE <- fu_my_share(sizes, CFG$job[1])
plan <- do.call(rbind, lapply(seq_along(MODELS), function(k) with(MODELS[[k]], data.frame(
  model = id, city = FU_CITY[[city]], fuel = FU_FUEL[[fuel]], pumps = n, K_max = kmax,
  grid_fits = 1L + (kmax - 1L) * length(CFG$phi_grid),
  boot_refits = CFG$B[1] * (1L + (kmax > 1L && isTRUE(CFG$boot_pooled[1]))),
  this_job = k %in% MINE))))
cat(sprintf("\nTHE PLAN: %d models (%d in this job), B = %d, X = %s\n", nrow(plan), length(MINE), CFG$B[1], XDEF))
print(plan[order(-plan$pumps), ], row.names = FALSE)
cat(sprintf("total: %d grid fits, up to %d bootstrap refits\n", sum(plan$grid_fits), sum(plan$boot_refits)))
FDIR <- file.path(OUT, "fit"); dir.create(FDIR, showWarnings = FALSE)
utils::write.csv(plan, file.path(FDIR, "plan.csv"), row.names = FALSE)
if (STAGE == "dry") quit(save = "no", status = 0)
if (STAGE != "fit") stop("unknown stage: ", STAGE, call. = FALSE)


## ===========================================================================
## Stage "fit": one model, from the grid to the tables
## ===========================================================================
fit_one <- function(M, seed) {
  message(sprintf("\n== %s: %d pumps, K_max = %d", M$id, M$n, M$kmax))
  od <- file.path(FDIR, M$id); dir.create(od, showWarnings = FALSE)
  sp <- spec_of(M$city, M$fuel)
  des <- fu_design(city_data(M$city), M$fuel, XDEF, sp$lags)
  mod <- fu_model(des)
  tag <- sprintf("%s_lambda%s_%s", M$id, format(sp$lambda), sp$from)

  ## 1. the grid of (K, phi) and the two-step selection, the pooled model
  ##    competing; the ridge on the four blocks of lags
  ic <- fu_cached(CACHE, paste0(tag, "_grid"),
    Stem::SCSTEM_Infocrit(mod, K_grid = seq_len(M$kmax), phi_grid = CFG$phi_grid,
                          knn = min(CFG$knn[1], des$d - 1L), distance = "geo",
                          min_cluster_size = CFG$m[1], seed = seed, verbose = TRUE,
                          alpha = 0, lambda = sp$lambda, penalize = des$penalize,
                          control = Stem::STEM_control()))
  sel <- Stem::SCSTEM_Select(ic, band = CFG$band, criterion = "BIC")
  pooled <- ic$fits[[which(ic$table$K == 1)[1]]]
  K <- sel$K_selected

  ## 2. the refit-with-clustering bootstrap of the selected model, and of the
  ##    pooled one for the comparison
  boot <- function(fit, name) fu_cached(CACHE, paste0(tag, "_", name), {
    bt <- Stem::SCSTEM_Bootstrap(fit, B = CFG$B[1], seed = seed, verbose = TRUE)
    list(info = bt$info, B_valid = bt$B_valid, inf = Stem::SCSTEM_BootInference(bt, level = CFG$level[1]),
         groups = bt$groups)
  })
  b_sel <- boot(sel$fit, sprintf("boot_K%d", K))
  b_pool <- if (K == 1L) b_sel else if (isTRUE(CFG$boot_pooled[1])) boot(pooled, "boot_K1") else NULL

  ## 3. the tests: bootstrap Wald regime by regime, for the selected and the
  ##    pooled model; likelihood ratio on the selected partition (no ridge);
  ##    pump by pump, F and HC1, with the Dumitrescu-Hurlin combination
  tests <- fu_tests(des)
  wald <- function(fit, b, label) if (!is.null(b)) do.call(rbind, lapply(names(tests), function(h)
    cbind(model = label, fu_wald_boot(fit, b$inf, des, tests[[h]], h))))
  wd <- rbind(wald(sel$fit, b_sel, sprintf("selected, K = %d", K)),
              if (K > 1L) wald(pooled, b_pool, "pooled"))
  lr <- fu_cached(CACHE, paste0(tag, "_lr"), fu_lr_fixed(des, sel$fit$group, K, CFG, seed))
  pt <- fu_pump_tests(des); dh <- fu_dh(pt)

  ## 4. diagnostics, the pass-through of the excise duties, regimes and brands
  dg <- fu_diagnostics(sel$fit, des)
  pul <- rbind(cbind(model = sprintf("selected, K = %d", K), fu_pulse_table(b_sel$inf, des)),
               if (K > 1L && !is.null(b_pool)) cbind(model = "pooled", fu_pulse_table(b_pool$inf, des)))
  br <- fu_brand_table(sel$fit, des)

  ## 5. the outputs of the model
  utils::write.csv(ic$table, file.path(od, "grid.csv"), row.names = FALSE)
  utils::write.csv(wd, file.path(od, "wald_bootstrap.csv"), row.names = FALSE)
  utils::write.csv(lr, file.path(od, "lr_fixed_partition.csv"), row.names = FALSE)
  utils::write.csv(pt, file.path(od, "pump_tests.csv"), row.names = FALSE)
  utils::write.csv(dh, file.path(od, "dumitrescu_hurlin.csv"), row.names = FALSE)
  utils::write.csv(dg, file.path(od, "diagnostics.csv"), row.names = FALSE)
  utils::write.csv(pul, file.path(od, "excise_pulses.csv"), row.names = FALSE)
  utils::write.csv(br$table, file.path(od, "regimes_by_brand.csv"), row.names = FALSE)
  est <- as.data.frame(b_sel$inf$summary)
  num <- suppressWarnings(as.integer(sub("^beta", "", est$parameter)))
  est$covariate <- ifelse(is.na(num), est$parameter, des$names[num])
  utils::write.csv(est, file.path(od, "estimates.csv"), row.names = FALSE)
  utils::write.csv(data.frame(des$meta, regime = sel$fit$group,
                              theta_bound = sel$fit$theta_bound[sel$fit$group]),
                   file.path(od, "regimes.csv"), row.names = FALSE)
  fu_map(des, sel$fit, file.path(od, "map.pdf"),
         sprintf("%s, %s (Y): K = %d, phi = %s, lambda = %s", FU_CITY[[M$city]], FU_FUEL[[M$fuel]], K,
                 format(sel$phi_selected), format(sp$lambda)))

  ## 6. one row of the summary of the case
  p_of <- function(h, model) { z <- wd[wd$block == h & startsWith(wd$model, model), ]
    if (!nrow(z)) NA_character_ else paste(format.pval(z$p_value, digits = 2), collapse = "; ") }
  row <- data.frame(model = M$id, city = FU_CITY[[M$city]], fuel = FU_FUEL[[M$fuel]], pumps = M$n,
                    weeks = des$Tn, lags = sp$from, lambda = sp$lambda, K_max = M$kmax, K = K,
                    phi = sel$phi_selected, B_valid = b_sel$B_valid,
                    p_nb_both_regimes = p_of("nb_both", "selected"), p_nb_both_pooled = p_of("nb_both", "pooled"),
                    p_nb_regimes = p_of("nb", "selected"), p_nb_other_regimes = p_of("nb_other", "selected"),
                    p_own_other_regimes = p_of("own_other", "selected"),
                    pumps_reject_nb_both_F_pct = dh$reject_F_pct[dh$block == "nb_both"],
                    pumps_reject_nb_both_hc1_pct = dh$reject_hc1_pct[dh$block == "nb_both"],
                    ari_regimes_brands = br$ari,
                    min_grid = fu_minutes(ic), min_boot = fu_minutes(b_sel), stem = SHA7,
                    stringsAsFactors = FALSE)
  saveRDS(row, file.path(od, "summary_row.rds"))
  message(sprintf("   K = %d, phi = %s | neighbour -> Y, p by regime: %s", K, format(sel$phi_selected),
                  row$p_nb_both_regimes))
  invisible(row)
}

for (k in MINE) {
  res <- tryCatch(fit_one(MODELS[[k]], seed = CFG$seed[1] + k), error = function(e) e)
  if (inherits(res, "error")) message("   FAILED: ", conditionMessage(res))
}

## the summary of every model finished so far with this Stem
rows <- lapply(MODELS, function(M) {
  f <- file.path(FDIR, M$id, "summary_row.rds")
  if (file.exists(f)) { r <- readRDS(f); if (identical(r$stem, SHA7)) r }
})
summ <- do.call(rbind, rows)
if (!is.null(summ) && nrow(summ)) {
  utils::write.csv(summ, file.path(FDIR, "summary.csv"), row.names = FALSE)
  cat(sprintf("\nSUMMARY: %d of %d models finished\n", nrow(summ), length(MODELS)))
  print(summ, row.names = FALSE)
}
cat("\noutputs in", OUT, "\n")
