## ---------------------------------------------------------------------------
## Where the time goes inside one EM fit, by sampling profiler.
##
##   Rscript dev/paper/11-engine-profile.R [d] [T] [iterations]
##
## The point is to see which of the remaining O(d^3) or O(T d) steps still
## matters at a network size we actually intend to fit, rather than guessing
## from the source. Rprof samples the call stack, so the number to read is
## "self" time -- the time spent inside a function rather than in what it calls.
## ---------------------------------------------------------------------------

local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  h <- if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
  source(file.path(h, "00-setup.R"), chdir = TRUE)
}, envir = globalenv())

a   <- commandArgs(TRUE)
dn  <- if (length(a) >= 1) as.integer(a[1]) else 200L
TT  <- if (length(a) >= 2) as.integer(a[2]) else 365L
it  <- if (length(a) >= 3) as.integer(a[3]) else 10L

stem_load()

set.seed(1)
co <- cbind(stats::runif(dn, 7.5, 13.5), stats::runif(dn, 44.7, 46.1))
dm <- geodist::geodist(co, measure = "geodesic")
L  <- chol(18 * diag(dn) + 5 * exp(-dm / 120000))
y  <- as.numeric(stats::arima.sim(list(ar = 0.9), TT)) * 2
X  <- cbind(1, rep(stats::runif(dn, 0, 250), each = TT),
            stats::rnorm(TT * dn, 30, 12))
Xa <- changedimension_covariates(X, dn, 3L, TT)
b  <- matrix(c(2.3, 0.004, 0.58), 3, 1)
z  <- matrix(NA_real_, TT, dn)
for (tt in seq_len(TT))
  z[tt, ] <- as.numeric(Xa[, , tt] %*% b) + y[tt] +
             as.numeric(crossprod(L, stats::rnorm(dn)))

ols <- stats::lm.fit(x = X, y = as.vector(z))
s2  <- stats::var(ols$residuals)
mod <- STEM_Model(z = z, covariates = X, coordinates = co,
  phi = list(beta = matrix(ols$coefficients, ncol = 1),
             sigma2eps = 0.6 * s2, sigma2omega = 0.3 * s2, theta = 1/60000,
             G = matrix(0.7, 1, 1), Sigmaeta = matrix(0.1 * s2, 1, 1),
             m0 = as.matrix(0), C0 = as.matrix(1)),
  K = matrix(1, dn, 1))

prof <- tempfile(fileext = ".out")
utils::Rprof(prof, interval = 0.005, line.profiling = FALSE)
fit <- STEM_Estimation(mod, precision = 0.01, max.iter = it,
                       distance = "geo", cov.spat = Sigmastar.exp)
utils::Rprof(NULL)

s <- summaryRprof(prof)
cat(sprintf("d = %d, T = %d, %d EM iterations, %.1f s total\n\n",
            dn, TT, it, s$sampling.time))
cat("self time, top 20\n")
print(utils::head(s$by.self, 20))
