## ---------------------------------------------------------------------------
## How do the information criteria behave as T and d grow?
##
## The question has a different answer under the null and under the alternative,
## and the asymmetry is the point.
##
##   Under the ALTERNATIVE, regimes really exist and the log-likelihood gain
##   from fitting them is proportional to the number of observations they
##   describe, hence O(dT). The parameter penalty is K log(dT) = O(log dT).
##   Every criterion therefore selects the true k once dT is large enough.
##
##   Under the NULL, there is nothing to find, and the gain is whatever the
##   search over partitions can extract by chance. For a FIXED partition the
##   likelihood ratio is O_p(1) in T, so the gain does not grow with T at all;
##   maximising it over the roughly 2^(d-1) binary partitions inflates it to
##   O(d), by the usual extreme-value argument for the maximum of many
##   correlated chi-squares. So the comparison is
##
##       gain ~ c * d        against        penalty ~ K log(dT).
##
##   Lengthening the series therefore HELPS the criteria -- the penalty grows,
##   the null gain does not -- while widening the network HURTS them, because
##   the null gain grows linearly in d and the penalty only logarithmically.
##
## This script measures both, on the generator of 06-dgp.R.
##
##   Rscript dev/paper/08-ic-scaling.R [--reps 30] [--cores 15]
## ---------------------------------------------------------------------------

PKG <- "C:/Users/paulm/OneDrive/Documenti/GitHub/Stem"
suppressMessages(pkgload::load_all(PKG, quiet = TRUE))
source(file.path("dev", "paper", "06-dgp.R"))

CACHE <- file.path("dev", "paper", "cache")
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
getopt <- function(f, d) { i <- match(f, args)
  if (is.na(i) || i == length(args)) d else args[i + 1L] }
REPS  <- as.integer(getopt("--reps", "30"))
CORES <- as.integer(getopt("--cores", "15"))

D_GRID <- c(18L, 36L, 54L)
T_GRID <- c(90L, 250L, 700L)
K_MAX  <- 4L

## log Stirling number of the second kind, for the extended BIC
logS2 <- function(d, k) {
  j <- 0:k
  keep <- (k - j) > 0
  terms <- lchoose(k, j)[keep] + d * log((k - j)[keep])
  sgn <- ((-1)^j)[keep]
  m <- max(terms)
  m + log(sum(sgn * exp(terms - m))) - lfactorial(k)
}

## synthetic coordinates on the extent of the real basin, so that the geometry
## is realistic while d can be varied freely
make_coords <- function(d, seed) {
  set.seed(seed)
  cbind(lon = stats::runif(d, 7.5, 13.5), lat = stats::runif(d, 44.7, 46.1))
}

one_rep <- function(r, dd, TT, hypothesis) {
  coords <- make_coords(dd, seed = 900000 + dd)
  labels <- as.integer(coords[, 1] > stats::median(coords[, 1])) + 1L

  psi <- if (hypothesis == "null") dgp_psi("S0") else dgp_psi("S1", level = 2L)
  x   <- dgp_covariate(coords, TT, seed = 7000 + r)
  z   <- dgp_simulate(labels, psi, x, coords, rho = 1,
                      seed = 10000 * (hypothesis == "alt") + 1000 * dd + r)
  X   <- dgp_design(x)

  phi0 <- list(beta = matrix(c(mean(z), 0), 2, 1),
               sigma2eps = 0.6 * stats::var(as.vector(z)),
               sigma2omega = 0.3 * stats::var(as.vector(z)),
               theta = 1 / 60000, G = matrix(0.7, 1, 1),
               Sigmaeta = matrix(0.1 * stats::var(as.vector(z)), 1, 1),
               m0 = as.matrix(0), C0 = as.matrix(1))
  mod <- STEM_Model(z = z, covariates = X, coordinates = coords,
                    phi = phi0, K = matrix(1, dd, 1))

  n <- dd * TT
  out <- list()
  for (kk in seq_len(K_MAX)) {
    f <- try(SCSTEM_Estimation(mod, k = kk, phi_penalty = 0, knn = min(5L, dd - 1L),
                          distance = "geo", precision = 0.1,
                          precision_full_dataset = 0.05, max_iter = 6,
                          seed = 1000 + r, verbose = FALSE), silent = TRUE)
    if (inherits(f, "try-error")) next
    ll <- unname(f$info_crit[["loglik"]])
    np <- unname(f$info_crit[["npar"]])
    if (!is.finite(ll)) next
    out[[length(out) + 1L]] <- data.frame(
      rep = r, d = dd, T = TT, hypothesis = hypothesis, k = kk,
      loglik = ll, npar = np,
      AIC  = -2 * ll + 2 * np,
      AICc = -2 * ll + 2 * np + 2 * np * (np + 1) / max(n - np - 1, 1),
      BIC  = -2 * ll + log(n) * np,
      EBIC = -2 * ll + log(n) * np + 2 * (if (kk > 1) logS2(dd, kk) else 0),
      ari  = Stem:::scstem_ari(labels, f$group))
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

cells <- expand.grid(d = D_GRID, T = T_GRID,
                     hypothesis = c("null", "alt"), stringsAsFactors = FALSE)
message("IC scaling study: ", nrow(cells), " cells x ", REPS, " replications")

cl <- parallel::makeCluster(min(CORES, parallel::detectCores() - 1L))
on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterExport(cl, c("PKG", "K_MAX", "logS2", "make_coords", "one_rep",
                              "dgp_base", "dgp_psi", "dgp_covariate",
                              "dgp_design", "dgp_simulate", "dgp_common_field",
                              "sigma_eta_of"), envir = environment())
invisible(parallel::clusterEvalQ(cl, {
  suppressMessages(pkgload::load_all(PKG, quiet = TRUE))
}))

res <- list()
for (i in seq_len(nrow(cells))) {
  cc <- cells[i, ]
  tag <- sprintf("d%d_T%d_%s", cc$d, cc$T, cc$hypothesis)
  f <- file.path(CACHE, paste0("icscale_", tag, ".rds"))
  if (file.exists(f)) { message("[cache] ", tag); res[[tag]] <- readRDS(f); next }
  message("[run  ] ", tag)
  t0 <- Sys.time()
  parallel::clusterExport(cl, c("cc"), envir = environment())
  out <- parallel::parLapply(cl, seq_len(REPS), function(r)
    one_rep(r, cc$d, cc$T, cc$hypothesis))
  out <- do.call(rbind, out[!vapply(out, is.null, logical(1))])
  message("   done in ",
          round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
  saveRDS(out, f)
  res[[tag]] <- out
}

all <- do.call(rbind, res)
saveRDS(all, file.path(CACHE, "icscale.rds"))

## ---- what each criterion selects -------------------------------------------
sel <- function(df, crit) {
  do.call(rbind, lapply(split(df, list(df$rep, df$d, df$T, df$hypothesis),
                             drop = TRUE), function(s) {
    data.frame(d = s$d[1], T = s$T[1], hypothesis = s$hypothesis[1],
               crit = crit, k = s$k[which.min(s[[crit]])])
  }))
}
picks <- do.call(rbind, lapply(c("AIC", "AICc", "BIC", "EBIC"),
                              function(cr) sel(all, cr)))
saveRDS(picks, file.path(CACHE, "icscale_picks.rds"))

cat("\n=== mean selected k, by criterion ===\n")
for (h in c("null", "alt")) {
  cat("\n--- ", h, " ---\n", sep = "")
  s <- picks[picks$hypothesis == h, ]
  print(round(tapply(s$k, list(paste0("d=", s$d, " T=", s$T), s$crit), mean), 2))
}

cat("\n=== share choosing k = 1 under the null / k = 2 under the alternative ===\n")
for (h in c("null", "alt")) {
  target <- if (h == "null") 1L else 2L
  s <- picks[picks$hypothesis == h, ]
  cat("\n--- ", h, ", target k = ", target, " ---\n", sep = "")
  print(round(100 * tapply(s$k == target,
                           list(paste0("d=", s$d, " T=", s$T), s$crit), mean), 1))
}

cat("\n=== the null gain: 2 * (loglik at k=2 minus loglik at k=1) ===\n")
nl <- all[all$hypothesis == "null" & all$k <= 2, ]
g <- do.call(rbind, lapply(split(nl, list(nl$rep, nl$d, nl$T), drop = TRUE),
                           function(s) {
  if (nrow(s) < 2) return(NULL)
  data.frame(d = s$d[1], T = s$T[1],
             gain = 2 * (s$loglik[s$k == 2] - s$loglik[s$k == 1]))
}))
cat("mean, by d and T (theory says it grows with d, not with T):\n")
print(round(tapply(g$gain, list(paste0("d=", g$d), paste0("T=", g$T)), mean), 1))
cat("\nthe BIC penalty it must beat, 9 log(dT):\n")
pen <- outer(D_GRID, T_GRID, function(a, b) 9 * log(a * b))
dimnames(pen) <- list(paste0("d=", D_GRID), paste0("T=", T_GRID))
print(round(pen, 1))
