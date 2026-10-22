## ===========================================================================
## SC-STEM fuel application: THE RESULTS OF THE REPLICATION PACKAGE
##
## Writes the compact results of the 66 models shipped with the replication
## package (replication/results/application in the repository) from the
## outputs of the original runs (fuels-main-setting*.R, stage "final";
## fuels-check-optima.R), in the format that replication/application/
## 02-models.R writes, with its functions. A full run of 02-models.R gives the
## same files: this script only spares the days of computing.
##
## What is read, for every model:
##   <root>/<case>/nn/final/summary.csv, <model>/{grid,regimes,pump_tests,
##     dumitrescu_hurlin,diagnostics,estimates_lambda0}.csv
##   <root>/<case>/nn/cache/09e9cdb/<model>_<lags>_boot_lambda0.rds and
##     _boot_pooled.rds (the bootstraps), _grid.rds (only the pooled fit, whose
##     estimates are needed at full precision: the bootstrap summary rounds
##     them to five digits)
##   <root>/checks/optima-parts/s<s>_<model>_selected.rds (the selected fit:
##     partition and estimates at full precision) and _v2.rds (the check of
##     the local maxima)
## The Wald tests are computed again with fu_wald() and checked against those
## of the original runs (wald_bootstrap.csv, test_wald.csv).
##
##     Rscript fuels-make-replication-results.R --root=<the fuels folder> --out=<replication/results/application>
##
## The estimates of the pooled fits are saved in <out>/../pooled_phi.rds while
## the grids are read (some GB each), so that a second run skips the reading;
## the file is deleted at the end.
##
## This is part of the development material, not of the replication package.
## ===========================================================================
suppressPackageStartupMessages(library(data.table))
HERE <- local({
  a <- commandArgs(trailingOnly = FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]), winslash = "/")) else getwd()
})
source(file.path(HERE, "..", "..", "replication", "application", "functions.R"))
args <- commandArgs(trailingOnly = TRUE)
opt <- function(k, d) { a <- grep(paste0("^--", k, "="), args, value = TRUE); if (length(a)) sub("^--[^=]+=", "", a[1]) else d }
ROOT <- normalizePath(opt("root", file.path(HERE, "fuels")), winslash = "/", mustWork = TRUE)
OUT <- normalizePath(opt("out", file.path(HERE, "..", "..", "replication", "results", "application")),
                     winslash = "/", mustWork = FALSE)
SETTINGS <- as.integer(strsplit(opt("settings", "1,2,3"), ",")[[1]])
SHA7 <- "09e9cdb"; KEY <- "L1t4-1t4-1t4-1t4"
FOLDER <- c("setting1-Ymajor-Xindependent", "setting2-Yindependent-Xmajor", "setting3-Ymajor-Xothermajor")
F_POOLED <- file.path(dirname(OUT), "pooled_phi.rds")
POOLED <- if (file.exists(F_POOLED)) readRDS(F_POOLED) else list()
worst <- 0

for (s in SETTINGS) {
  case <- file.path(ROOT, FOLDER[s]); fin <- file.path(case, "nn", "final"); cache <- file.path(case, "nn", "cache", SHA7)
  sm <- utils::read.csv(file.path(fin, "summary.csv"), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(sm))) {
    id <- sm$model[i]; cc <- sub("_.*$", "", id); fu <- names(FU_FUEL)[match(sm$fuel[i], FU_FUEL)]
    K <- sm$K[i]; fd <- file.path(fin, id); od <- file.path(OUT, sprintf("setting%d", s), id)
    dir.create(od, recursive = TRUE, showWarnings = FALSE)
    w <- function(d, f) if (!is.null(d)) utils::write.csv(cbind(setting = s, model = id, d), file.path(od, f), row.names = FALSE)

    ## the names of the covariates, and the blocks of the design (the same in
    ## every model: when the file of a model is missing, that of another model
    ## of the setting; the number of coefficients is checked below)
    f0 <- file.path(fd, "estimates_lambda0.csv")
    if (!file.exists(f0)) f0 <- Sys.glob(file.path(fin, "*", "estimates_lambda0.csv"))[1]
    e0 <- utils::read.csv(f0, stringsAsFactors = FALSE)
    e0 <- unique(e0[grepl("^beta", e0$parameter), c("parameter", "covariate")])
    nm <- e0$covariate[order(as.integer(sub("^beta", "", e0$parameter)))]
    blocks <- lapply(stats::setNames(FU_BLOCKS, FU_BLOCKS), function(b) grep(paste0("^", b, "_l"), nm))
    blocks$pulses <- grep("^pulse_", nm)
    des <- list(names = nm, blocks = blocks)

    ## the selected fit, its bootstrap, the pooled fit and its bootstrap
    sel <- readRDS(file.path(ROOT, "checks", "optima-parts", sprintf("s%d_%s_selected.rds", s, id)))
    reg <- utils::read.csv(file.path(fd, "regimes.csv"), stringsAsFactors = FALSE, colClasses = c(y_site = "character", nn_site = "character"))
    stopifnot(identical(as.integer(reg$regime), as.integer(sel$group)), nrow(sel$phi_hat) == K)
    b0 <- readRDS(file.path(cache, sprintf("%s_%s_boot_lambda0.rds", id, KEY)))
    bP <- if (K == 1L) b0 else readRDS(file.path(cache, sprintf("%s_%s_boot_pooled.rds", id, KEY)))
    stopifnot(length(grep("^beta", dimnames(b0$inf$aligned)[[3]])) == length(nm),
              length(grep("^beta", colnames(sel$phi_hat))) == length(nm))
    key <- paste(s, id)
    if (K > 1L && is.null(POOLED[[key]])) {
      message(sprintf("%s  reading the grid of setting %d, %s", format(Sys.time(), "%H:%M"), s, id))
      ic <- readRDS(file.path(cache, sprintf("%s_%s_grid.rds", id, KEY)))
      POOLED[[key]] <- ic$fits[[which(ic$table$K == 1)[1]]]$phi_hat
      rm(ic); invisible(gc()); saveRDS(POOLED, F_POOLED)
    }
    pooled_phi <- if (K == 1L) sel$phi_hat else POOLED[[key]]

    ## the Wald tests, checked against the original ones
    wald <- rbind(fu_wald(sel$phi_hat, sel$group, b0$inf, des, "selected"),
                  if (K > 1L) fu_wald(pooled_phi, rep(1L, length(sel$group)), bP$inf, des, "pooled"))
    old <- utils::read.csv(file.path(fd, "wald_bootstrap.csv"), stringsAsFactors = FALSE)
    old$fit <- ifelse(startsWith(old$model, "selected"), "selected", "pooled")
    old$regime <- as.character(old$regime)
    m <- merge(wald, old[, c("fit", "block", "regime", "W")], by = c("fit", "block", "regime"))
    stopifnot(nrow(m) == nrow(old))
    tw <- utils::read.csv(file.path(fd, "test_wald.csv"), stringsAsFactors = FALSE)
    m2 <- merge(wald[wald$fit == "selected", ], tw[, c("block", "regime", "statistic", "p_boot")], by = c("block", "regime"),
                suffixes = c("", "_old"))
    d1 <- max(abs(m$statistic / m$W - 1)); d2 <- max(abs(m2$statistic / m2$statistic_old - 1))
    d3 <- max(abs(m2$p_boot - m2$p_boot_old))
    worst <- max(worst, d1, d2, d3)
    if (d1 > 1e-8 || d2 > 1e-8 || d3 > 1e-12) stop("the Wald tests of setting ", s, ", ", id, " differ from the original ones")

    ## the files of 02-models.R
    g <- utils::read.csv(file.path(fd, "grid.csv"), stringsAsFactors = FALSE)
    w(g, "grid.csv")
    names(reg)[match(c("nn_site", "nn_brand", "nn_dist_km"), names(reg))] <- c("x_site", "x_brand", "x_dist_km")
    w(reg[, c("y_site", "y_brand", "y_group", "y_type", "y_municipality", "y_lon", "y_lat", "x_site", "x_brand",
              "x_dist_km", "regime", "theta_bound")], "regimes.csv")
    w(rbind(fu_estimates(b0$inf, nm, "selected"), if (K > 1L) fu_estimates(bP$inf, nm, "pooled")), "estimates.csv")
    data.table::fwrite(cbind(setting = s, model = id, fu_draws(b0$inf, nm)), file.path(od, "draws.csv.gz"))
    w(wald, "wald.csv")
    for (f in c("pump_tests.csv", "dumitrescu_hurlin.csv", "diagnostics.csv"))
      w(utils::read.csv(file.path(fd, f), stringsAsFactors = FALSE, colClasses = if (f == "pump_tests.csv") c(y_site = "character") else NA),
        f)
    pt <- file.path(ROOT, "checks", "optima-parts", sprintf("s%d_%s_v2.rds", s, id))
    opt_min <- NA_real_
    if (K > 1L) {
      p <- readRDS(pt); r <- p$regimes; wb <- p$wald[p$wald$block == "nb_both", ]
      wb <- wb[match(r$regime, wb$regime), ]
      w(data.frame(regime = r$regime, size = r$size, loglik_selected = r$loglik_selected, loglik_ls = r$loglik_ls,
                   diff = r$diff, converged_selected = r$converged_selected, converged_ls = r$converged_ls,
                   theta_selected = r$theta_selected, theta_ls = r$theta_ls, G_selected = r$G_selected, G_ls = r$G_ls,
                   sigma2eps_selected = r$sigma2eps_selected, sigma2eps_ls = r$sigma2eps_ls,
                   W_nb_both_selected = wb$W_selected, p_nb_both_selected = wb$p_selected,
                   W_nb_both_ls = wb$W_ls, p_nb_both_ls = wb$p_ls), "optima.csv")
      co <- p$coef
      w(data.frame(regime = co$regime, covariate = co$covariate, block = ifelse(co$block == "other", "intercept", co$block),
                   beta_selected = co$beta_selected, beta_ls = co$beta_ls, se_boot = co$se_boot), "optima_coefficients.csv")
      opt_min <- r$minutes[1]
    }
    brands <- reg$y_brand
    row <- data.frame(
      setting = s, model = id, city = cc, fuel = FU_FUEL[[fu]], pumps = sm$pumps[i], weeks = sm$weeks[i],
      K_max = sm$K_max[i], K = K, phi = sm$phi[i], BIC = g$BIC[g$K == K & g$phi == sm$phi[i]][1],
      sizes = paste(tabulate(sel$group, K), collapse = "/"), B_valid = b0$B_valid, B_valid_pooled = bP$B_valid,
      ari_regimes_brands = if (length(unique(brands)) > 1 && K > 1) fu_ari(sel$group, brands) else NA_real_,
      ari_draws_median = stats::median(b0$inf$stability$ARI),
      minutes_grid = sm$min_grid[i], minutes_boot = sm$min_boot[i],
      minutes_boot_pooled = if (K > 1L) attr(bP, "fu_minutes") else 0, minutes_optima = opt_min,
      stem = SHA7, stringsAsFactors = FALSE)
    stopifnot(isTRUE(all.equal(as.numeric(row$ari_regimes_brands), as.numeric(sm$ari_regimes_brands[i]))))
    utils::write.csv(row, file.path(od, "model.csv"), row.names = FALSE)
    message(sprintf("setting %d, %s: done", s, id))
  }
  f <- Sys.glob(file.path(OUT, sprintf("setting%d", s), "*", "model.csv"))
  utils::write.csv(do.call(rbind, lapply(f, utils::read.csv, stringsAsFactors = FALSE)),
                   file.path(OUT, sprintf("setting%d", s), "models.csv"), row.names = FALSE)
}
if (identical(sort(SETTINGS), 1:3)) unlink(F_POOLED)
cat(sprintf("largest relative difference of the Wald statistics from the original ones: %.2g\n", worst))
