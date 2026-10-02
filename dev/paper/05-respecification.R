## ---------------------------------------------------------------------------
## An alternative specification of the application.
##
## The reference specification regresses PM2.5 on the PM10 measured at the same
## station on the same day. That is close to regressing a part on the whole:
## PM2.5 is by definition the fine fraction of PM10, the two correlate at 0.91
## in these data, and their median ratio is 0.68. The consequence is visible in
## the fit -- the coefficient of PM10 comes out at 0.58 to 0.60 in every
## estimated regime, and altitude, whose marginal R2 is zero because the network
## was deliberately restricted to stations below 250 m, explains nothing. With
## 83% of the variance of the response accounted for by an almost accounting
## identity, there is very little structure left for spatial regimes to carry,
## and the cross-validation of 04-cv.R duly finds none.
##
## This script fits the model that the data can actually support: PM2.5 is
## explained by its own latent spatio-temporal structure, with only exogenous
## regressors. The design matrix keeps the intercept and the altitude and adds
## annual and semiannual harmonics of the day of the year, which are exogenous
## by construction; everything else -- level, persistence, spatial range and the
## split of the residual variance -- is left to the model, which is where the
## regimes can differ.
##
##   Rscript dev/paper/05-respecification.R [--cores 4] [--folds 5]
## ---------------------------------------------------------------------------

PKG <- "C:/Users/paulm/OneDrive/Documenti/GitHub/Stem"
Sys.setenv(STEM_CV_NORUN = "1")
source(file.path("dev", "paper", "04-cv.R"))   # definitions only: folds, run_fold
Sys.unsetenv("STEM_CV_NORUN")

args <- commandArgs(trailingOnly = TRUE)
getopt <- function(flag, default) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1L]
}
CORES <- as.integer(getopt("--cores", "4"))
NFOLD <- as.integer(getopt("--folds", "5"))

## ---- the alternative design matrix -----------------------------------------
doy <- as.integer(format(povalley$dates, "%j"))
harm <- cbind(sin(2 * pi * doy / 365.25), cos(2 * pi * doy / 365.25),
              sin(4 * pi * doy / 365.25), cos(4 * pi * doy / 365.25))

## stacked by station, as the package expects
X2 <- do.call(rbind, lapply(seq_len(d), function(i) {
  cbind(1, povalley$altitude[i], harm)
}))
colnames(X2) <- c("intercept", "altitude", "sin1", "cos1", "sin2", "cos2")
stopifnot(nrow(X2) == Tn * d)

r2 <- ncol(X2)
ols2 <- stats::lm.fit(x = X2, y = as.vector(z0))
cat("R2 of the alternative design:",
    round(1 - stats::var(ols2$residuals) / stats::var(as.vector(z0)), 4), "\n")

phi2 <- list(beta = matrix(ols2$coefficients, ncol = 1),
             sigma2eps = 0.5 * stats::var(ols2$residuals),
             sigma2omega = 0.4 * stats::var(ols2$residuals),
             theta = 1 / 60000,
             G = matrix(0.7, 1, 1),
             Sigmaeta = matrix(0.1 * stats::var(ols2$residuals), 1, 1),
             m0 = as.matrix(0), C0 = as.matrix(1))

## The helpers of 04-cv.R close over X0, phi0 and r, so they are rebound here:
## the fold machinery is identical, only the design changes.
X0   <- X2
phi0 <- phi2
r    <- r2

## ---- the grid ---------------------------------------------------------------
mod2 <- STEM_Model(z = z0, covariates = X0, coordinates = coords,
                   phi = phi0, A = matrix(1, d, 1))

grid2 <- if (file.exists(file.path(CACHE, "grid_respec.rds"))) {
  readRDS(file.path(CACHE, "grid_respec.rds"))
} else {
  message("grid over k and phi, alternative specification")
  g <- SCSTEM_Infocrit(mod2, K_grid = 1:4, phi_grid = c(0, 0.25, 0.5, 1, 2),
                       knn = KNN, distance = "geo", precision = 0.1,
                       precision_full_dataset = 0.01, max_iter = 8,
                       seed = SEED, verbose = TRUE)
  saveRDS(g, file.path(CACHE, "grid_respec.rds"))
  g
}
cat("\n=== grid, alternative specification ===\n")
print(grid2$table[, c("k", "phi", "k_eff", "admissible", "loglik", "BIC", "min_size")])

## ---- cross-validation -------------------------------------------------------
schemes <- c("random", "LKTO", "LKLO", "LKLHTO")
jobs <- list()
for (sch in schemes) {
  fl <- make_folds(sch, NFOLD)
  for (j in seq_along(fl)) for (ci in seq_along(CONFIGS)) {
    jobs[[length(jobs) + 1L]] <- list(scheme = sch, fold = j, cfg = ci, f = fl[[j]])
  }
}
message("cross-validation on the alternative specification: ", length(jobs), " jobs")

cl <- parallel::makeCluster(min(CORES, parallel::detectCores() - 1L))
on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterExport(cl, c("PKG", "z0", "X0", "coords", "d", "Tn", "r",
                              "phi0", "cov_rows", "krige_dropped", "run_fold",
                              "CONFIGS", "KNN", "SEED"), envir = environment())
invisible(parallel::clusterEvalQ(cl, {
  suppressMessages(pkgload::load_all(PKG, quiet = TRUE))
}))

t0 <- Sys.time()
out <- parallel::parLapply(cl, jobs, function(jb) {
  rr <- run_fold(jb$f, CONFIGS[[jb$cfg]])
  data.frame(scheme = jb$scheme, fold = jb$fold,
             k = CONFIGS[[jb$cfg]]$k, phi = CONFIGS[[jb$cfg]]$phi,
             rmse = if (is.null(rr)) NA_real_ else unname(rr["rmse"]),
             mae  = if (is.null(rr)) NA_real_ else unname(rr["mae"]),
             n    = if (is.null(rr)) NA_real_ else unname(rr["n"]))
})
cv2 <- do.call(rbind, out)
message("done in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
saveRDS(cv2, file.path(CACHE, "cv_respec.rds"))

cat("\n=== mean RMSE by scheme and configuration (alternative specification) ===\n")
m <- tapply(cv2$rmse, list(paste0("k=", cv2$k), cv2$scheme), mean, na.rm = TRUE)
print(round(m, 4))
cat("\n=== relative to the pooled model (= 100) ===\n")
print(round(100 * sweep(m, 2, m["k=1", ], "/"), 1))

cat("\n=== paired comparison against pooled, fold by fold ===\n")
for (sch in colnames(m)) {
  a <- cv2[cv2$scheme == sch, ]
  p <- a$rmse[a$k == 1]
  for (kk in 2:4) {
    q <- a$rmse[a$k == kk]
    if (all(is.na(q)) || all(is.na(p))) next
    tt <- stats::t.test(q, p, paired = TRUE)
    cat(sprintf("%-7s k=%d vs k=1: mean diff %+7.4f  (p = %.3f)%s\n",
                sch, kk, mean(q - p, na.rm = TRUE), tt$p.value,
                ifelse(tt$p.value < 0.05,
                       ifelse(mean(q - p, na.rm = TRUE) < 0, "  BETTER", "  WORSE"), "")))
  }
}
