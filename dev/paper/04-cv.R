## ---------------------------------------------------------------------------
## Superseded by Stem::SCSTEM_CV() (2026-09-28), which implements the same four
## schemes as an optional validation of chosen models, not as a selection rule;
## (k, phi) are chosen by SCSTEM_Select(). The argument below, that information
## criteria cannot arbitrate, compared the BIC penalty with the likelihood gain
## under the alternative instead of the null, and is withdrawn. The script is kept
## because 05-respecification.R sources its definitions.
##
## Spatio-temporal cross-validation for the choice of k and phi.
##
## Information criteria cannot arbitrate here. With n = dT observations the BIC
## penalty for an extra regime is 9 log(dT) ~ 100, while the likelihood gain
## between neighbouring values of k is in the hundreds or thousands: the penalty
## is two orders of magnitude too small. The extended BIC of Chen and Chen,
## which charges the size of the model space, addresses exactly the right
## pathology -- the partition is chosen by searching over S(d,k) ~ 1e20
## candidates -- but its correction is 2 gamma log S(d,k), of order d log k,
## whereas the likelihood gain is of order dT. The ratio is O(1/T), so with
## T = 1826 it cannot bite either. Out-of-sample validation is therefore not a
## fallback but the only criterion with the right scaling.
##
## Four blocking schemes, following the taxonomy of Otto, Fasso and Maranzano
## (2024) and the implementation of CAST (Meyer et al.):
##
##   random   leave random cells out           optimistic: neighbours in both
##                                             space and time stay in training
##   LKLO     leave K locations out            spatial prediction at unmonitored
##                                             sites
##   LKTO     leave K time blocks out          temporal prediction
##   LKLHTO   leave K locations and H times    both bands removed from training,
##            out                              scored on their intersection
##
##   Rscript dev/paper/04-cv.R [--folds 5] [--cores 8]
## ---------------------------------------------------------------------------

PKG <- "C:/Users/paulm/OneDrive/Documenti/GitHub/Stem"
suppressMessages(pkgload::load_all(PKG, quiet = TRUE))

CACHE <- file.path("dev", "paper", "cache")
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
getopt <- function(flag, default) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1L]
}
NFOLD <- as.integer(getopt("--folds", "5"))
CORES <- as.integer(getopt("--cores", "8"))
SEED  <- 20260904
KNN   <- 5

## the configurations put to the test: the pooled model and the clustered ones
## at a moderate penalty
CONFIGS <- list(list(k = 1L, phi = 0),
                list(k = 2L, phi = 0.5),
                list(k = 3L, phi = 0.5),
                list(k = 4L, phi = 0.5))

data(povalley)
z0 <- povalley$z; X0 <- povalley$covariates; coords <- povalley$coords
d  <- ncol(z0); Tn <- nrow(z0); r <- ncol(X0)

ols <- stats::lm.fit(x = as.matrix(X0), y = as.vector(z0))
phi0 <- list(beta = matrix(ols$coefficients, ncol = 1),
             sigma2eps = 0.5 * stats::var(ols$residuals),
             sigma2omega = 0.4 * stats::var(ols$residuals),
             theta = 1 / 60000,
             G = matrix(0.7, 1, 1),
             Sigmaeta = matrix(0.1 * stats::var(ols$residuals), 1, 1),
             m0 = as.matrix(0), C0 = as.matrix(1))

## covariate rows for a subset of locations, in the package's stacking order
cov_rows <- function(idx) as.vector(outer(seq_len(Tn), (idx - 1L) * Tn, "+"))

## ---------------------------------------------------------------------------
## Folds. Each fold is a list with `blank`, the cells removed from training, and
## `score`, the cells on which the prediction is evaluated. The two coincide
## except for LKLHTO, where whole bands are removed from training but only their
## intersection is scored, which is what makes that scheme strict.
## ---------------------------------------------------------------------------
make_folds <- function(scheme, K, seed = SEED) {
  set.seed(seed)
  folds <- vector("list", K)

  if (scheme == "random") {
    grp <- matrix(sample(rep_len(seq_len(K), Tn * d)), Tn, d)
    for (j in seq_len(K)) {
      m <- grp == j
      folds[[j]] <- list(blank = m, score = m, drop_loc = integer(0))
    }

  } else if (scheme == "LKLO") {
    lg <- sample(rep_len(seq_len(K), d))
    for (j in seq_len(K)) {
      m <- matrix(FALSE, Tn, d); m[, lg == j] <- TRUE
      folds[[j]] <- list(blank = m, score = m, drop_loc = which(lg == j))
    }

  } else if (scheme == "LKTO") {
    ## contiguous blocks in time: random time points would leave the immediate
    ## neighbours of every held-out day in the training set
    tb <- cut(seq_len(Tn), breaks = K, labels = FALSE)
    for (j in seq_len(K)) {
      m <- matrix(FALSE, Tn, d); m[tb == j, ] <- TRUE
      folds[[j]] <- list(blank = m, score = m, drop_loc = integer(0))
    }

  } else if (scheme == "LKLHTO") {
    lg <- sample(rep_len(seq_len(K), d))
    tb <- cut(seq_len(Tn), breaks = K, labels = FALSE)
    for (j in seq_len(K)) {
      band <- matrix(FALSE, Tn, d)
      band[, lg == j] <- TRUE          # the whole location band
      band[tb == j, ] <- TRUE          # the whole time band
      sc <- matrix(FALSE, Tn, d)
      sc[tb == j, lg == j] <- TRUE     # scored only on the intersection
      folds[[j]] <- list(blank = band, score = sc, drop_loc = which(lg == j))
    }
  } else {
    stop("unknown scheme: ", scheme)
  }
  folds
}

## ---------------------------------------------------------------------------
## Prediction at locations that were removed from the fit altogether.
##
## The held-out site inherits the regime of its nearest retained station -- the
## label of an unmonitored site is predicted by spatial contiguity, which is
## exactly what the Potts penalty encodes -- and is then predicted by the
## conditional mean given the retained stations of that regime,
##
##   zhat(s0,t) = x(s0,t)'beta_k + K(s0) yhat_t + w' Sigma_k^{-1} r_t ,
##
## with w the covariance between s0 and those stations. Sigma_k^{-1} does not
## depend on t, so this is one solve and a matrix product over all time points.
## ---------------------------------------------------------------------------
krige_dropped <- function(fitobj, keep_idx, drop_idx, group_keep) {
  Xall <- Stem:::changedimension_covariates(X0, d, r, Tn)   # d x r x Tn
  pred <- matrix(NA_real_, Tn, length(drop_idx))
  dmall <- geodist::geodist(coords, measure = "geodesic")

  ## nearest retained station of each dropped one, and hence its regime
  nearest <- vapply(drop_idx, function(i) keep_idx[which.min(dmall[i, keep_idx])],
                    integer(1))
  reg_of_drop <- group_keep[match(nearest, keep_idx)]

  for (g in unique(reg_of_drop)) {
    f <- fitobj[[g]]
    if (is.null(f)) next
    members <- keep_idx[group_keep == g]
    ph  <- f$estimates$phi.hat
    ysm <- as.numeric(f$estimates$y.smoothed)
    b   <- matrix(as.numeric(ph$beta), ncol = 1)

    dmg <- dmall[members, members, drop = FALSE]
    Sig <- ph$sigma2eps * diag(length(members)) +
           ph$sigma2omega * exp(-ph$theta * dmg)
    Sinv <- solve(Sig)

    ## residuals of the members at the fitted parameters
    Rt <- matrix(NA_real_, Tn, length(members))
    for (tt in seq_len(Tn)) {
      Rt[tt, ] <- z0[tt, members] -
        as.numeric(Xall[members, , tt] %*% b) - ysm[tt]
    }
    Rt[is.na(Rt)] <- 0     # a member blanked by the fold contributes nothing

    for (jj in which(reg_of_drop == g)) {
      i <- drop_idx[jj]
      w <- ph$sigma2omega * exp(-ph$theta * dmall[i, members])
      corr <- as.numeric(Rt %*% (Sinv %*% matrix(w, ncol = 1)))
      pred[, jj] <- as.numeric(t(Xall[i, , ]) %*% b) + ysm + corr
    }
  }
  pred
}

## ---------------------------------------------------------------------------
## One fold, one configuration
## ---------------------------------------------------------------------------
run_fold <- function(fold, cfg) {
  z <- z0
  z[fold$blank] <- NA_real_

  keep_idx <- setdiff(seq_len(d), fold$drop_loc)
  drop_idx <- fold$drop_loc

  ## drop the held-out locations from the model: a location observed nowhere
  ## cannot be assigned to a regime
  zk <- z[, keep_idx, drop = FALSE]
  Xk <- X0[cov_rows(keep_idx), , drop = FALSE]
  ck <- coords[keep_idx, , drop = FALSE]

  ## every retained location must keep at least one observation
  if (any(apply(zk, 2, function(cc) all(is.na(cc))))) return(NULL)

  mod <- STEM_Model(z = zk, covariates = Xk, coordinates = ck,
                    phi = phi0, A = matrix(1, length(keep_idx), 1))

  zhat_keep <- matrix(NA_real_, Tn, length(keep_idx))
  fits <- NULL; grp <- NULL

  if (cfg$k == 1L) {
    f <- try(STEM_Estimation(mod, precision = 0.01, max.iter = 40,
                             distance = "geo"), silent = TRUE)
    if (inherits(f, "try-error")) return(NULL)
    zhat_keep <- STEM_Complete(f, distance = "geo")
    fits <- list(f); grp <- rep(1L, length(keep_idx))
  } else {
    f <- try(SCSTEM_Estimation(mod, K = cfg$k, phi_penalty = cfg$phi, knn = KNN,
                          distance = "geo", precision = 0.1,
                          precision_full_dataset = 0.01, max_iter = 8,
                          seed = SEED, verbose = FALSE), silent = TRUE)
    if (inherits(f, "try-error")) return(NULL)
    zhat_keep <- SCSTEM_Complete(f)
    fits <- f$fit_list; grp <- f$group
  }

  zhat <- matrix(NA_real_, Tn, d)
  zhat[, keep_idx] <- zhat_keep
  if (length(drop_idx)) {
    zhat[, drop_idx] <- krige_dropped(fits, keep_idx, drop_idx, grp)
  }

  ## the fitted objects of a fold are large and are of no further use: releasing
  ## them here keeps a worker that runs many folds from accumulating them
  rm(fits, zhat_keep); invisible(gc(verbose = FALSE))

  sc <- fold$score & !is.na(z0)
  err <- zhat[sc] - z0[sc]
  err <- err[is.finite(err)]
  if (!length(err)) return(NULL)
  c(rmse = sqrt(mean(err^2)), mae = mean(abs(err)), n = length(err))
}

## ---------------------------------------------------------------------------
## Drive
## ---------------------------------------------------------------------------
## Sourcing the file with STEM_CV_NORUN set defines the pieces without running
## the whole exercise, which is what the smoke test does.
if (nzchar(Sys.getenv("STEM_CV_NORUN"))) {
  message("STEM_CV_NORUN set: definitions loaded, nothing run.")
} else {

schemes <- c("random", "LKTO", "LKLO", "LKLHTO")
jobs <- list()
for (sch in schemes) {
  fl <- make_folds(sch, NFOLD)
  for (j in seq_along(fl)) for (ci in seq_along(CONFIGS)) {
    jobs[[length(jobs) + 1L]] <- list(scheme = sch, fold = j, cfg = ci,
                                      f = fl[[j]])
  }
}
message("cross-validation: ", length(jobs), " jobs (",
        length(schemes), " schemes x ", NFOLD, " folds x ",
        length(CONFIGS), " configurations)")

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
  r <- run_fold(jb$f, CONFIGS[[jb$cfg]])
  data.frame(scheme = jb$scheme, fold = jb$fold,
             k = CONFIGS[[jb$cfg]]$k, phi = CONFIGS[[jb$cfg]]$phi,
             rmse = if (is.null(r)) NA_real_ else unname(r["rmse"]),
             mae  = if (is.null(r)) NA_real_ else unname(r["mae"]),
             n    = if (is.null(r)) NA_real_ else unname(r["n"]))
})
res <- do.call(rbind, out)
message("done in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
        " min")

saveRDS(res, file.path(CACHE, "cv.rds"))

cat("\n=== mean RMSE by scheme and configuration ===\n")
print(round(tapply(res$rmse, list(paste0("k=", res$k), res$scheme),
                   mean, na.rm = TRUE), 4))
cat("\n=== failed folds ===\n")
print(table(res$scheme[is.na(res$rmse)]))

}
