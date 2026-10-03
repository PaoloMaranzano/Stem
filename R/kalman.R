#' @keywords internal
#' @noRd

### ---------------------------------------------------------------------------
### The Kalman filter and smoother of the STEM model
### ---------------------------------------------------------------------------
### A STEM fit is the state-space model
###
###   z_t = X_t beta + A' y_t + e_t,      e_t ~ N(0, sigma2omega Sigma*),
###   y_t = G y_{t-1} + eta_t,            eta_t ~ N(0, Sigmaeta),
###   y_0 ~ N(m0, C0),
###
### with Sigma* = R(theta) + b I and b = sigma2eps / sigma2omega. Internally the
### parameters travel as a list `phi` with elements A (p x d, the transposed
### loading matrix, fixed), sigma2omega, logtheta, logb, beta (r x 1), G,
### Sigmaeta, m0 (p x 1) and C0 (fixed). The data travel as a list `dat`, built
### once per fit by stem_em_data(): z (n x d), covariates (d x r x n), dist
### (d x d), obs (stem_obs_index(z)), n, d, r, p and cov.spat, the
### correlation function.
###
### kalman() runs, at the parameters `phi`, the three recursions every
### estimation algorithm of the package rests on:
###
###   the Kalman filter (filtering(), filterstep()): the filtered moments
###     y_{t|t}, P_{t|t} and, as its by-product, the log-likelihood by the
###     prediction-error decomposition;
###   the fixed-interval smoother (smoothing()): the smoothed moments
###     y_{t|n}, P_{t|n}, and those of the initial state, y_{0|n} and P_{0|n};
###   the lag-one covariance smoother of Shumway and Stoffer (1982):
###     P_{t,t-1|n} = Cov(y_t, y_{t-1} | z).
###
### The E-step (estep.R) forms the expected sufficient statistics from them, and
### STEM_Estimation() calls kalman() once more at the estimates, for the
### log-likelihood and the smoothed states it reports.
###
###   regularization  a ridge added to the matrices inverted here (the gains of
###                   the smoother at times 0 and n); 0 by default, see
###                   STEM_Estimation()
`kalman` <- function(phi, dat, regularization = 0) {
  n <- dat$n
  p <- dat$p
  G <- phi$G
  W <- phi$Sigmaeta
  C0 <- phi$C0

  ### the state-space model at phi
  Sigmastar <- dat$cov.spat(d = dat$d, logb = phi$logb, logtheta = phi$logtheta, dist = dat$dist)
  ss <- list(z = dat$z, Fmat = phi$A, Gmat = G,
             Vmat = phi$sigma2omega * Sigmastar, Wmat = W,
             m0 = t(phi$m0), C0 = C0, XXX = dat$covariates, beta = phi$beta,
             flag.cov = TRUE, n = n, p = p, d = dat$d)

  ### the Kalman filter and the fixed-interval smoother
  filt <- filtering(ss)
  smo <- smoothing(filt)
  Cf <- filt$C                         # filtered variances, a list of n
  Msm <- as.matrix(smo$m)              # smoothed means, n x p
  Csm <- smo$C                         # smoothed variances, a list of n

  ### The smoother gains B_t = P_{t|t} G' (G P_{t|t} G' + W)^{-1}, t = 1..n,
  ### and B_0 of the initial state; with them, the smoothed moments of y_0.
  R1 <- G %*% C0 %*% t(G) + W
  B0 <- C0 %*% t(G) %*% solve(diag(regularization, p) + R1)
  B <- vector("list", n)
  for (tt in seq_len(n - 1L)) {
    B[[tt]] <- Cf[[tt]] %*% t(G) %*% solve(G %*% Cf[[tt]] %*% t(G) + W)
  }
  B[[n]] <- Cf[[n]] %*% t(G) %*% solve(diag(regularization, p) + G %*% Cf[[n]] %*% t(G) + W)
  y0 <- t(ss$m0) + B0 %*% (matrix(Msm[1L, ], ncol = 1) - G %*% t(ss$m0))
  P0 <- C0 + B0 %*% (Csm[[1L]] - R1) %*% t(B0)

  ### The lag-one covariances, by the backward recursion of Shumway and Stoffer
  ### (1982). It starts at time n from the gain actually used there, built on
  ### the observed rows only (Durbin and Koopman 2012, Sect. 4.10); with
  ### nothing observed at time n that gain is zero.
  CCC <- vector("list", n)
  oi_n <- dat$obs$idx[[n]]
  Rn <- G %*% Cf[[n - 1L]] %*% t(G) + W
  if (length(oi_n) == 0L) {
    CCC[[n]] <- G %*% Cf[[n - 1L]]
  } else {
    Fn <- phi$A[, oi_n, drop = FALSE]
    Qn <- t(Fn) %*% Rn %*% Fn + ss$Vmat[oi_n, oi_n, drop = FALSE]
    An <- Rn %*% Fn %*% solve(diag(regularization, nrow(Qn)) + Qn)
    CCC[[n]] <- (diag(p) - An %*% t(Fn)) %*% G %*% Cf[[n - 1L]]
  }
  for (tt in n:2) {
    Bb <- if (tt >= 3L) B[[tt - 2L]] else B0
    CCC[[tt - 1L]] <- Cf[[tt - 1L]] %*% t(Bb) +
      B[[tt - 1L]] %*% (CCC[[tt]] - G %*% Cf[[tt - 1L]]) %*% t(Bb)
  }

  list(loglik = filt$loglik, ss = ss, Sigmastar = Sigmastar,
       filtered = list(m = filt$m, C = Cf),
       smoothed = list(m = smo$m, C = Csm, lag1 = CCC, y0 = y0, P0 = P0))
}
