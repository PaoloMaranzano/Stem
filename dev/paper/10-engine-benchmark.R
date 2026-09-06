## ---------------------------------------------------------------------------
## Benchmark and numerical check of the estimation engine against the last
## unoptimised state of the code, commit 7d0b094.
##
## Two changes are being measured together:
##   the Woodbury form of the Kalman filter, which replaces a d x d inversion at
##   every time point by p x p algebra against a cached factorisation of
##   Sigma_e, and
##   the matrix form of the M-step accumulations, which replaces three loops
##   building lists of n matrices by crossprods.
##
## The check that matters is not the timing but the agreement: both changes are
## algebraic identities, so the estimates must coincide to floating point.
##
##   Rscript dev/paper/10-engine-benchmark.R [reference commit]
## ---------------------------------------------------------------------------

local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  h <- if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
  source(file.path(h, "00-setup.R"), chdir = TRUE)
}, envir = globalenv())

PKG  <- stem_root()
REF  <- (function(a) if (length(a)) a[1] else "7d0b094")(commandArgs(TRUE))
TMP  <- tempfile("stem-ref-"); dir.create(TMP)

load_env <- function(dir) {
  e <- new.env(parent = globalenv())
  for (f in list.files(dir, pattern = "[.]R$", full.names = TRUE))
    sys.source(f, envir = e, keep.source = FALSE)
  e
}

cur <- load_env(file.path(PKG, "R"))
for (f in system(sprintf('git -C "%s" ls-tree --name-only %s R/', PKG, REF),
                 intern = TRUE)) {
  writeLines(system(sprintf('git -C "%s" show %s:%s', PKG, REF, f), intern = TRUE),
             file.path(TMP, basename(f)))
}
ref <- load_env(TMP)

## a synthetic network of the requested size, generated from the model itself
make <- function(env, dn, TT, seed = 1) {
  set.seed(seed)
  co <- cbind(stats::runif(dn, 7.5, 13.5), stats::runif(dn, 44.7, 46.1))
  dm <- geodist::geodist(co, measure = "geodesic")
  L  <- chol(18 * diag(dn) + 5 * exp(-dm / 120000))
  y  <- as.numeric(stats::arima.sim(list(ar = 0.9), TT)) * 2
  X  <- cbind(1, rep(stats::runif(dn, 0, 250), each = TT),
              stats::rnorm(TT * dn, 30, 12))
  Xa <- env$changedimension_covariates(X, dn, 3L, TT)
  b  <- matrix(c(2.3, 0.004, 0.58), 3, 1)
  z  <- matrix(NA_real_, TT, dn)
  for (tt in seq_len(TT))
    z[tt, ] <- as.numeric(Xa[, , tt] %*% b) + y[tt] +
               as.numeric(crossprod(L, stats::rnorm(dn)))
  ols <- stats::lm.fit(x = X, y = as.vector(z))
  s2  <- stats::var(ols$residuals)
  env$STEM_Model(z = z, covariates = X, coordinates = co,
    phi = list(beta = matrix(ols$coefficients, ncol = 1),
               sigma2eps = 0.6 * s2, sigma2omega = 0.3 * s2, theta = 1/60000,
               G = matrix(0.7, 1, 1), Sigmaeta = matrix(0.1 * s2, 1, 1),
               m0 = as.matrix(0), C0 = as.matrix(1)),
    K = matrix(1, dn, 1))
}

fit <- function(env, dn, TT, iters) {
  m <- make(env, dn, TT)
  tm <- system.time(f <- env$STEM_Estimation(m, precision = 0.01,
         max.iter = iters, distance = "geo", cov.spat = env$Sigmastar.exp))
  list(t = tm[["elapsed"]],
       s = c(unlist(f$estimates$phi.hat), loglik = f$estimates$loglik))
}

cat(sprintf("%4s %5s | %10s %10s %8s | %s\n",
            "d", "T", "before (s)", "after (s)", "speedup", "max rel. diff"))
for (cfg in list(c(36, 365), c(60, 365), c(100, 365), c(200, 365), c(400, 200))) {
  dn <- cfg[1]; TT <- cfg[2]
  a <- try(fit(ref, dn, TT, 10L), silent = TRUE)
  b <- fit(cur, dn, TT, 10L)
  if (inherits(a, "try-error")) {
    ## the unoptimised code takes log(dmvnorm(.)) rather than the log-density
    ## itself, and that underflows to -Inf once the network is large enough, so
    ## there is nothing to compare against
    cat(sprintf("%4d %5d | %10s %10.1f %8s | %s\n", dn, TT,
                "FAILED", b$t, "-", "reference underflows to -Inf"))
    next
  }
  rel <- max(abs(a$s - b$s) / pmax(abs(a$s), 1e-12), na.rm = TRUE)
  cat(sprintf("%4d %5d | %10.1f %10.1f %7.1fx | %.2e\n",
              dn, TT, a$t, b$t, a$t / b$t, rel))
}
