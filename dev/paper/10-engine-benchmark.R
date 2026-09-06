setwd("C:/Users/paulm/OneDrive/Documenti/GitHub/Stem")
suppressMessages(pkgload::load_all(".", quiet = TRUE))

mk <- function(dn, TT, seed = 1) {
  set.seed(seed)
  co <- cbind(stats::runif(dn, 7.5, 13.5), stats::runif(dn, 44.7, 46.1))
  dm <- geodist::geodist(co, measure = "geodesic")
  Sig <- 18 * diag(dn) + 5 * exp(-dm / 120000)
  L <- chol(Sig)
  y <- as.numeric(stats::arima.sim(list(ar = 0.9), TT)) * 2
  X <- cbind(1, rep(stats::runif(dn, 0, 250), each = TT),
             stats::rnorm(TT * dn, 30, 12))
  Xa <- Stem:::changedimension_covariates(X, dn, 3, TT)
  b <- matrix(c(2.3, 0.004, 0.58), 3, 1)
  z <- matrix(NA_real_, TT, dn)
  for (tt in seq_len(TT))
    z[tt, ] <- as.numeric(Xa[, , tt] %*% b) + y[tt] +
               as.numeric(crossprod(L, stats::rnorm(dn)))
  ols <- stats::lm.fit(x = X, y = as.vector(z)); s2 <- stats::var(ols$residuals)
  STEM_Model(z = z, covariates = X, coordinates = co,
             phi = list(beta = matrix(ols$coefficients, ncol = 1),
                        sigma2eps = 0.6*s2, sigma2omega = 0.3*s2, theta = 1/60000,
                        G = matrix(0.7,1,1), Sigmaeta = matrix(0.1*s2,1,1),
                        m0 = as.matrix(0), C0 = as.matrix(1)),
             K = matrix(1, dn, 1))
}

cat(sprintf("%4s %5s | %9s %9s %8s | %s\n",
            "d", "T", "R (s)", "fast (s)", "speedup", "max rel. difference"))
for (cfg in list(c(36,365), c(60,365), c(100,365), c(200,365), c(400,200))) {
  dn <- cfg[1]; TT <- cfg[2]
  mod <- mk(dn, TT)
  t1 <- system.time(a <- try(STEM_Estimation(mod, precision=0.01, max.iter=10,
                       distance="geo", engine="R"), silent=TRUE))[["elapsed"]]
  t2 <- system.time(b <- try(STEM_Estimation(mod, precision=0.01, max.iter=10,
                       distance="geo", engine="fast"), silent=TRUE))[["elapsed"]]
  if (inherits(a,"try-error") || inherits(b,"try-error")) {
    cat(sprintf("%4d %5d | %9.1f %9.1f %8s | %s\n", dn, TT, t1, t2, "-",
                if (inherits(a,"try-error")) "R engine failed" else "fast engine failed"))
    next
  }
  sa <- c(unlist(a$estimates$phi.hat), loglik = a$estimates$loglik)
  sb <- c(unlist(b$estimates$phi.hat), loglik = b$estimates$loglik)
  rel <- max(abs(sa - sb) / pmax(abs(sa), 1e-12), na.rm = TRUE)
  cat(sprintf("%4d %5d | %9.1f %9.1f %7.1fx | %.2e\n", dn, TT, t1, t2, t1/t2, rel))
}
