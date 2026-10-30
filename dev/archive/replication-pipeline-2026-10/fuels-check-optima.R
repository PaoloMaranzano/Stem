## ===========================================================================
## SC-STEM fuel application: CHECK OF THE LOCAL MAXIMA OF THE SELECTED FITS
##
## For every model of the three settings, the partition selected in stage
## "grid" (two-step rule on the BIC) is refitted regime by regime from least
## squares, the start of a fit with no alternation behind it, with the EM
## algorithm stopping only when the parameters and the log-likelihood have
## both settled (em_stop = "all", at most 2000 iterations); the log-likelihood
## of every regime is set beside that of the selected fit, whose regimes start
## from the last iteration of the alternation. Two starts reaching the same
## value do not prove a global maximum; two different values show that the
## likelihood of the regime has more than one. A pooled model (K* = 1) has one
## start only and is listed without a refit.
##
##     Rscript fuels-check-optima.R --root=<the fuels folder> --workers=6
##
## Options: --root (default fuels/ beside the scripts), --settings (default
## 1,2,3), --workers (default 4). Phase 1 reads the cached grids one at a time
## (those of Rome take several GB of memory) and keeps the selected fits;
## phase 2 refits them on --workers processes. Every model is saved as soon
## as it is done, in <root>/checks/optima-parts/, and a model already saved is
## skipped, so the check can be stopped and resumed. Output:
## <root>/checks/optima_by_regime.csv, optima_coefficients.csv (the
## coefficients at the two maxima, with the bootstrap standard errors of the
## selected fit) and optima_wald.csv (the bootstrap Wald tests at both).
## Neither the package nor the cached stages of the main scripts are changed.
##
## The settings below are those of the main scripts (fuels-main-setting*.R):
## keep them equal.
##
## This is part of the REPLICATION MATERIAL of the paper, not of the Stem
## package.
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

args <- commandArgs(trailingOnly = TRUE)
opt <- function(k, d) { a <- grep(paste0("^--", k, "="), args, value = TRUE)
  if (length(a)) sub("^--[^=]+=", "", a[1]) else d }
ROOT <- normalizePath(opt("root", file.path(HERE, "fuels")), winslash = "/", mustWork = TRUE)
SETTINGS <- as.integer(strsplit(opt("settings", "1,2,3"), ",")[[1]])
WORKERS <- as.integer(opt("workers", "4"))

## the settings of the main scripts
LAGS <- list(own = 1:4, own_other = 1:4, nb = 1:4, nb_other = 1:4)
LAG_KEY <- "L1t4-1t4-1t4-1t4"
M_MIN <- 5L; KNN <- 5L; BAND <- c(0.025, 0.2); SEED <- 20261005L
seed_of <- function(city, fuel) SEED + 10L * match(city, names(FU_CITY)) + match(fuel, names(FU_FUEL))

SHA7 <- fu_require_stem()
CHK <- file.path(ROOT, "checks"); PARTS <- file.path(CHK, "optima-parts")
dir.create(PARTS, recursive = TRUE, showWarnings = FALSE)

## the models, from the summary of stage "final" of every setting
tasks <- do.call(rbind, lapply(SETTINGS, function(s) {
  case <- fu_case_dir(ROOT, s)
  sm <- utils::read.csv(file.path(case, "nn", "final", "summary.csv"), stringsAsFactors = FALSE)
  data.frame(setting = s, case = case, model = sm$model, city = sub("_.*$", "", sm$model),
             fuel = names(FU_FUEL)[match(sm$fuel, FU_FUEL)], pumps = sm$pumps, K = sm$K, phi = sm$phi,
             stringsAsFactors = FALSE)
}))
tasks <- tasks[order(-tasks$pumps), ]
tasks$part <- file.path(PARTS, sprintf("s%d_%s_v2.rds", tasks$setting, tasks$model))
tasks$sel <- file.path(PARTS, sprintf("s%d_%s_selected.rds", tasks$setting, tasks$model))
cat(sprintf("%d models; %d already checked\n", nrow(tasks), sum(file.exists(tasks$part))))

## Phase 1: the selected fit of every model, one cached grid at a time
for (i in which(!file.exists(tasks$part) & !file.exists(tasks$sel))) {
  tk <- tasks[i, ]
  ic <- readRDS(file.path(tk$case, "nn", "cache", SHA7, sprintf("%s_%s_grid.rds", tk$model, LAG_KEY)))
  sel <- Stem::SCSTEM_Select(ic, band = BAND, criterion = "BIC")
  if (sel$K_selected != tk$K || abs(sel$phi_selected - tk$phi) > 1e-9)
    stop("the grid of ", tk$model, " (setting ", tk$setting, ") does not give the K* and phi* of the summary",
         call. = FALSE)
  f <- sel$fit
  saveRDS(list(group = f$group, loglik_g = f$loglik_g, em_converged = f$em_converged, phi_hat = f$phi_hat),
          tk$sel)
  rm(ic, sel, f); invisible(gc())
  cat(sprintf("%s  selected fit of setting %d, %s\n", format(Sys.time(), "%H:%M"), tk$setting, tk$model))
}

## Phase 2: the refit from least squares, on WORKERS processes. Besides the
## log-likelihoods, every regime gives the coefficients at the two maxima, with
## the bootstrap standard errors of the selected fit (stage "final"), and the
## bootstrap Wald tests of fu_tests() at both, on the covariance of the same
## bootstrap draws.
check_one <- function(tk) {
  if (file.exists(tk$part)) return(tk$model)
  suppressPackageStartupMessages(library(data.table))
  source(file.path(HERE, "fuels-functions.R"))
  f <- readRDS(tk$sel); K <- tk$K
  e <- new.env(); load(file.path(tk$case, "data.RData"), envir = e)
  des <- fu_design(e$fuels_data[e$fuels_data$city == tk$city, , drop = FALSE], tk$fuel, "nn", LAGS)
  t0 <- Sys.time()
  r <- if (K == 1L) NULL else Stem::SCSTEM_Estimation(
    fu_model(des), K = K, phi_penalty = 0, init_partition = f$group, max_iter = 0, distance = "geo",
    knn = min(KNN, des$d - 1L), min_cluster_size = M_MIN, lambda = 0, seed = seed_of(tk$city, tk$fuel),
    control = Stem::STEM_control(em_stop = "all", em_maxit = 2000L), verbose = FALSE)
  par <- function(x, p) if (is.null(x)) rep(NA_real_, K) else unname(x$phi_hat[, p])
  out <- data.frame(setting = tk$setting, model = tk$model, pumps = tk$pumps, K = K, phi = tk$phi,
                    regime = seq_len(K), size = tabulate(f$group, K),
                    loglik_selected = f$loglik_g, loglik_ls = if (is.null(r)) NA_real_ else r$loglik_g,
                    converged_selected = f$em_converged, converged_ls = if (is.null(r)) NA else r$em_converged,
                    theta_selected = unname(f$phi_hat[, "theta"]), theta_ls = par(r, "theta"),
                    G_selected = unname(f$phi_hat[, "G"]), G_ls = par(r, "G"),
                    sigma2eps_selected = unname(f$phi_hat[, "sigma2eps"]), sigma2eps_ls = par(r, "sigma2eps"),
                    minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
                    stringsAsFactors = FALSE)
  out$diff <- out$loglik_ls - out$loglik_selected
  coef <- wald <- NULL
  if (!is.null(r)) {
    inf <- readRDS(file.path(tk$case, "nn", "cache", SHA7_W,
                             sprintf("%s_%s_boot_lambda0.rds", tk$model, LAG_KEY_W)))$inf
    sm <- as.data.frame(inf$summary)
    blk <- rep("other", length(des$names))
    for (b in names(des$blocks)) blk[des$blocks[[b]]] <- b
    coef <- do.call(rbind, lapply(seq_len(K), function(g) {
      p <- paste0("beta", seq_along(des$names))
      s <- sm[sm$cluster == g, ]
      data.frame(setting = tk$setting, model = tk$model, regime = g, size = sum(f$group == g),
                 diff_loglik = out$diff[g], covariate = des$names, block = blk,
                 beta_selected = unname(f$phi_hat[g, p]), beta_ls = unname(r$phi_hat[g, p]),
                 se_boot = s$se[match(p, s$parameter)], stringsAsFactors = FALSE)
    }))
    tests <- fu_tests(des)
    wd <- function(fit, which) do.call(rbind, lapply(names(tests), function(h)
      cbind(fit = which, fu_wald_boot(fit, inf, des, tests[[h]], h))))
    a <- wd(list(phi_hat = f$phi_hat, group = f$group), "selected")
    b <- wd(list(phi_hat = r$phi_hat, group = f$group), "ls")
    wald <- data.frame(setting = tk$setting, model = tk$model, block = a$block, regime = a$regime,
                       pumps_regime = a$pumps, diff_loglik = out$diff[a$regime],
                       W_selected = a$W, p_selected = a$p_value, W_ls = b$W, p_ls = b$p_value,
                       stringsAsFactors = FALSE)
  }
  saveRDS(list(regimes = out, coef = coef, wald = wald), tk$part)
  tk$model
}
todo <- which(!file.exists(tasks$part))
cat(sprintf("\n%s  refitting %d models on %d processes\n", format(Sys.time(), "%H:%M"), length(todo), WORKERS))
if (length(todo)) {
  SHA7_W <- SHA7; LAG_KEY_W <- LAG_KEY
  cl <- parallel::makePSOCKcluster(min(WORKERS, length(todo)))
  parallel::clusterExport(cl, c("HERE", "LAGS", "KNN", "M_MIN", "SEED", "seed_of", "check_one", "SHA7_W",
                                "LAG_KEY_W"))
  res <- parallel::parLapplyLB(cl, split(tasks[todo, ], seq_along(todo)), function(tk)
    tryCatch(check_one(tk), error = function(e) paste("FAILED", tk$model, conditionMessage(e))))
  parallel::stopCluster(cl)
  bad <- grep("^FAILED", unlist(res), value = TRUE)
  if (length(bad)) cat(bad, sep = "\n")
}

## The tables
parts <- lapply(tasks$part[file.exists(tasks$part)], readRDS)
tab <- do.call(rbind, lapply(parts, `[[`, "regimes"))
tab <- tab[order(tab$setting, -tab$pumps, tab$model, tab$regime), ]
coef <- do.call(rbind, lapply(parts, `[[`, "coef"))
wald <- do.call(rbind, lapply(parts, `[[`, "wald"))
coef$z_diff <- (coef$beta_ls - coef$beta_selected) / coef$se_boot
utils::write.csv(tab, file.path(CHK, "optima_by_regime.csv"), row.names = FALSE)
utils::write.csv(coef, file.path(CHK, "optima_coefficients.csv"), row.names = FALSE)
utils::write.csv(wald, file.path(CHK, "optima_wald.csv"), row.names = FALSE)
cat(sprintf("\n%s  %d regimes of %d models in %s\n", format(Sys.time(), "%H:%M"), nrow(tab),
            length(unique(paste(tab$setting, tab$model))), CHK))
d <- tab[is.finite(tab$diff), ]
cat(sprintf("regimes refitted: %d; |difference| > 1: %d; > 10: %d; least squares higher by > 10: %d\n",
            nrow(d), sum(abs(d$diff) > 1), sum(abs(d$diff) > 10), sum(d$diff > 10)))
