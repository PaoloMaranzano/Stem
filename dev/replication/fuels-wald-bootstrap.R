## ===========================================================================
## SC-STEM fuel application: THE BOOTSTRAP WALD TEST
##
## The test of leadership is post-estimation inference on the selected fit of
## a model: the hypotheses are the blocks of lags of fu_tests() (in a regime,
## the lags of the neighbor have zero coefficients). The test propagates the
## uncertainty of the estimated partition, because the bootstrap re-estimates
## SC-STEM at (K*, phi*) on every sample and matches the regimes of the sample
## to the observed ones by the majority rule.
##
## Row "all" is the hypothesis in all the regimes at once, on q K degrees of
## freedom: the coefficients of the K regimes stacked, with their joint
## bootstrap covariance.
##
## THE WALD TEST, for one model and one hypothesis (q lags):
##   1. on the data: the coefficients b_k of the q lags in every regime k of
##      the selected fit (stage "grid", K* and phi*);
##   2. the bootstrap: that of stage "final" of the main scripts
##      (Stem::SCSTEM_Bootstrap: B samples drawn regime by regime from the
##      fitted model, with the observed covariates, coordinates and missing
##      values, SC-STEM re-estimated on each; Stem::SCSTEM_BootInference: the
##      regimes matched by the majority rule), read from the cache of the main
##      scripts; computed here, with the same seed, when it is not there or B
##      differs;
##   3. V_k, the covariance of the draws of b_k, and W_k = b_k' V_k^-1 b_k;
##   4. p_chisq against the chi-square on q degrees of freedom; p_boot against
##      the draws W*_kb = (b*_kb - b_k)' V_k^-1 (b*_kb - b_k): centered at the
##      estimate, a draw reproduces the deviation of b_k from its true value,
##      which is what W_k measures when the null holds.
##
## Only exported functions of Stem are used (SCSTEM_Select, SCSTEM_Bootstrap,
## SCSTEM_BootInference): the test is specific to the application and lives in
## the replication material, not in the package.
##
## Edit THE SETUP and press Source in RStudio, or from a terminal
##
##     Rscript fuels-wald-bootstrap.R --settings=1,2,3 --cores=4
##
## (command-line options override THE SETUP: --settings, --cities, --fuels,
## --blocks, --B, --cores, --root, --out). Requirements: the stages grid and
## final of the main scripts, whose cached grids, bootstraps and summaries are
## read; Stem at the pinned commit. A bootstrap computed here is saved, when
## complete, in <out>/<case>/nn/wald-boot/. Outputs:
## <out>/<case>/nn/final/<model>/test_wald.csv and
## <out>/<case>/nn/final/test_wald_summary.csv.
##
## This is part of the REPLICATION MATERIAL of the paper, not of the Stem
## package.
## ===========================================================================


## ===========================================================================
## THE SETUP. This is the part meant to be edited.
## ===========================================================================
RUN <- list(
  ## the settings: 1, 2, 3, or several
  settings = 1:3,
  ## the cities: "all", or codes among RM MI NA TO PA BA CT BO ME FI VE
  cities   = "all",
  ## the fuel of Y: "g" gasoline, "d" diesel
  fuels    = c("g", "d"),
  ## the hypotheses (blocks of lags of fu_tests()): "nb_both" (the lags of the
  ## neighbor, both fuels: the leadership test of the paper), "nb", "nb_other",
  ## "own_other"
  blocks   = "nb_both",
  ## bootstrap draws (stage "final" of the main scripts uses 100)
  B        = 100L,
  ## R processes the draws of a model are shared among
  cores    = 4L,
  ## the folder with the data and the outputs of the main scripts: "" for
  ## fuels/ beside the scripts
  root     = "",
  ## where this script writes: "" for the same folder as root
  out      = ""
)


## ===========================================================================
## Machinery: nothing below needs editing
## ===========================================================================
suppressPackageStartupMessages(library(data.table))
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
RUN <- fu_bt_setup(RUN, HERE)
SHA7 <- fu_require_stem()

## --- the statistic --------------------------------------------------------------
## The Wald statistic of the coefficients est, with the covariance of their
## bootstrap draws D (one row per draw), and its draws centered at est.
wd_stat <- function(est, D) {
  if (nrow(D) <= length(est)) return(list(stat = NA_real_, null = numeric(0), se = rep(NA_real_, length(est))))
  V <- stats::cov(D)
  Vi <- tryCatch(solve(V), error = function(e) NULL)
  if (is.null(Vi)) return(list(stat = NA_real_, null = numeric(0), se = sqrt(diag(V))))
  C <- sweep(D, 2, est)
  list(stat = drop(crossprod(est, Vi %*% est)), null = rowSums((C %*% Vi) * C), se = sqrt(diag(V)))
}
## the bootstrap of the fitted model: the one of stage "final" (same seed, same
## B), or computed here
wd_bootstrap <- function(tk, fit, run) {
  f_main <- file.path(tk$case, "nn", "cache", SHA7, sprintf("%s_%s_boot_lambda0.rds", tk$model, FU_BT$lag_key))
  if (file.exists(f_main)) {
    b0 <- readRDS(f_main); s <- attr(b0, "fu_settings")
    if (identical(as.integer(s$B), as.integer(run$B)) && identical(as.integer(s$partition), as.integer(fit$group)))
      return(b0$inf)
  }
  f_state <- file.path(tk$out_case, "nn", "wald-boot", sprintf("%s_B%d.rds", tk$model, run$B))
  if (file.exists(f_state)) return(readRDS(f_state))
  dir.create(dirname(f_state), recursive = TRUE, showWarnings = FALSE)
  bt <- Stem::SCSTEM_Bootstrap(fit, B = run$B, seed = fu_bt_seed(tk), verbose = TRUE, cores = run$cores)
  inf <- Stem::SCSTEM_BootInference(bt, level = FU_BT$level)
  saveRDS(inf, f_state)
  inf
}

## --- one model, every hypothesis --------------------------------------------------
wd_model <- function(tk, run) {
  sel <- fu_bt_selected(tk, SHA7); des <- sel$des; fit <- sel$fit; K <- tk$K
  inf <- wd_bootstrap(tk, fit, run)
  A <- inf$aligned
  ari <- stats::median(inf$stability$ARI, na.rm = TRUE)
  tests <- fu_tests(des)
  fmt <- function(x) paste(sprintf("%.4f", x), collapse = "; ")
  do.call(rbind, lapply(intersect(run$blocks, names(tests)), function(blk) {
    par <- paste0("beta", tests[[blk]]); ip <- match(par, dimnames(A)[[3]]); q <- length(ip)
    est <- lapply(seq_len(K), function(g) unname(fit$phi_hat[g, par]))
    draws <- lapply(seq_len(K), function(g) matrix(A[, g, ip], ncol = q))
    res <- lapply(seq_len(K), function(g) wd_stat(est[[g]], draws[[g]][stats::complete.cases(draws[[g]]), , drop = FALSE]))
    D_all <- do.call(cbind, draws)
    res[[K + 1L]] <- wd_stat(unlist(est), D_all[stats::complete.cases(D_all), , drop = FALSE])
    fu_bt_table(tk, "Wald", blk, stat = vapply(res, `[[`, 1, "stat"), df = c(rep(q, K), K * q),
                pumps = tabulate(fit$group, K), null = lapply(res, `[[`, "null"), B = run$B, ari = ari,
                extra = data.frame(coef = c(vapply(est, fmt, ""), NA),
                                   se_boot = c(vapply(res[seq_len(K)], function(r) fmt(r$se), ""), NA)))
  }))
}

## --- the run ------------------------------------------------------------------
tasks <- fu_bt_tasks(RUN)
cat(sprintf("bootstrap Wald test: %d models; blocks %s; B = %d on %d processes; Stem %s\n  outputs in %s\n\n",
            nrow(tasks), paste(RUN$blocks, collapse = ", "), RUN$B, RUN$cores, SHA7, RUN$out))
fu_bt_run(tasks, function(tk) wd_model(tk, RUN), "test_wald.csv", "bootstrap Wald test")
