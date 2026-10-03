#' @keywords internal
#' @noRd

### ---------------------------------------------------------------------------
### The E-step, shared by every estimation algorithm
### ---------------------------------------------------------------------------
### The model of a STEM fit is the state-space model
###
###   z_t = X_t beta + A' y_t + e_t,      e_t ~ N(0, sigma2omega Sigma*),
###   y_t = G y_{t-1} + eta_t,            eta_t ~ N(0, Sigmaeta),
###   y_0 ~ N(m0, C0),
###
### with Sigma* = R(theta) + b I, b = sigma2eps / sigma2omega. Internally the
### parameters travel as a list `phi` with elements A (p x d, the transposed
### loading matrix, fixed), sigma2omega, logtheta, logb, beta (r x 1), G,
### Sigmaeta, m0 (p x 1) and C0 (fixed).
###
### The data travel as a list `dat`, built once per fit by stem_em_data():
### z (n x d), covariates (d x r x n), dist (d x d), obs (stem_obs_index(z)),
### n, d, r, p and cov.spat, the correlation function.
###
### Given the parameters, the Kalman filter and smoother return the moments of
### the latent states given all the data. From them the E-step forms everything
### the M-step of either algorithm needs:
###
###   S11 = sum_t E(y_t y_t' | z),  S00 = sum_t E(y_{t-1} y_{t-1}' | z),
###   S10 = sum_t E(y_t y_{t-1}' | z),
###   BB  = sum_t E(e_t e_t' | z),  e_t = z_t - X_t beta - A' y_t,
###   zhat, the observations completed by their conditional expectations,
###   y0 and P0, the smoothed mean and variance of the initial state,
###
### and the log-likelihood at the parameters, the by-product of the filter.

### The filter and the smoother at the parameters `phi`: the log-likelihood,
### the filtered and the smoothed moments. On its own, this is what a fit
### reports at its estimates.
`stem_filter_smooth` <- function(phi, dat) {
  Sigmastar <- dat$cov.spat(d = dat$d, logb = phi$logb, logtheta = phi$logtheta, dist = dat$dist)
  ss <- list(z = dat$z, Fmat = phi$A, Gmat = phi$G,
             Vmat = phi$sigma2omega * Sigmastar, Wmat = phi$Sigmaeta,
             m0 = t(phi$m0), C0 = phi$C0, XXX = dat$covariates, beta = phi$beta,
             flag.cov = TRUE, n = dat$n, p = dat$p, d = dat$d)
  filt <- filtering(ss)
  list(ss = ss, Sigmastar = Sigmastar, filt = filt, smo = smoothing(filt),
       loglik = filt$loglik)
}

### The E-step at the parameters `phi`.
###   regularization  a ridge added to the matrices inverted here (the gains at
###                   time 0 and n and the blocks of the missing values); 0 by
###                   default, see STEM_Estimation()
`stem_estep` <- function(phi, dat, regularization = 0) {

  fs  <- stem_filter_smooth(phi, dat)
  n   <- dat$n
  d   <- dat$d
  p   <- dat$p
  G   <- phi$G
  W   <- phi$Sigmaeta
  C0  <- phi$C0
  Cf  <- fs$filt$C                     # filtered variances, a list of n
  Msm <- as.matrix(fs$smo$m)           # smoothed means, n x p
  Csm <- fs$smo$C                      # smoothed variances, a list of n

  ### the smoother gains B_t = C_t G' (G C_t G' + W)^{-1}, t = 1..n, and the
  ### gain B_0 of the initial state; with them, the smoothed moments of y_0
  R1 <- G %*% C0 %*% t(G) + W
  B0 <- C0 %*% t(G) %*% solve(diag(regularization, p) + R1)
  B <- vector("list", n)
  for (tt in seq_len(n - 1L)) {
    B[[tt]] <- Cf[[tt]] %*% t(G) %*% solve(G %*% Cf[[tt]] %*% t(G) + W)
  }
  B[[n]] <- Cf[[n]] %*% t(G) %*% solve(diag(regularization, p) + G %*% Cf[[n]] %*% t(G) + W)
  y0 <- t(fs$ss$m0) + B0 %*% (matrix(Msm[1L, ], ncol = 1) - G %*% t(fs$ss$m0))
  P0 <- C0 + B0 %*% (Csm[[1L]] - R1) %*% t(B0)

  ### The lag-one covariances Cov(y_t, y_{t-1} | z), by the backward recursion
  ### of Shumway and Stoffer (1982). It starts at time n from the gain actually
  ### used there, built on the observed rows only (Durbin and Koopman 2012,
  ### Sect. 4.10); with nothing observed at time n that gain is zero.
  CCC <- vector("list", n)
  oi_n <- dat$obs$idx[[n]]
  Rn <- G %*% Cf[[n - 1L]] %*% t(G) + W
  if (length(oi_n) == 0L) {
    CCC[[n]] <- G %*% Cf[[n - 1L]]
  } else {
    Fn <- phi$A[, oi_n, drop = FALSE]
    Qn <- t(Fn) %*% Rn %*% Fn + fs$ss$Vmat[oi_n, oi_n, drop = FALSE]
    An <- Rn %*% Fn %*% solve(diag(regularization, nrow(Qn)) + Qn)
    CCC[[n]] <- (diag(p) - An %*% t(Fn)) %*% G %*% Cf[[n - 1L]]
  }
  for (tt in n:2) {
    Bb <- if (tt >= 3L) B[[tt - 2L]] else B0
    CCC[[tt - 1L]] <- Cf[[tt - 1L]] %*% t(Bb) +
      B[[tt - 1L]] %*% (CCC[[tt]] - G %*% Cf[[tt - 1L]]) %*% t(Bb)
  }

  ### The three sums of outer products of the smoothed states, each a
  ### crossprod of the matrix of the smoothed means plus a sum of variances.
  lag  <- seq_len(n - 1L)
  S11 <- crossprod(Msm) + Reduce(`+`, Csm)
  S00 <- y0 %*% t(y0) + P0 + crossprod(Msm[lag, , drop = FALSE]) + Reduce(`+`, Csm[lag])
  S10 <- matrix(Msm[1L, ], ncol = 1) %*% t(y0) +
    crossprod(Msm[-1L, , drop = FALSE], Msm[lag, , drop = FALSE]) + Reduce(`+`, CCC)

  ### The second moment of the measurement error, BB = sum_t E(e_t e_t' | z).
  ### At a complete time point it is Z C^s_t Z' + r_t r_t', with r_t the
  ### residual at the smoothed state, and over the complete time points the sum
  ### collapses to Z (sum_t C^s_t) Z' + R'R. With missing entries the quantity
  ### has to be completed, because the EM algorithm maximizes the expected
  ### complete-data log-likelihood: the missing block of the residual is
  ### predicted from the observed one by Gaussian conditioning on Sigma_e (the
  ### algebra of kriging at a fixed time point), and its conditional variance
  ### Omega_t is added back. Writing S_t for the matrix that lifts a residual
  ### on the observed rows to the full vector (identity on the observed rows,
  ### the conditioning matrix on the missing ones),
  ###
  ###   BB_t = S_t Zo C^s_t Zo' S_t' + Omega_t + (S_t r_o)(S_t r_o)' .
  ###
  ### The conditional mean of a missing element is NOT the signal alone, which
  ### would be exact only for a diagonal Sigma_e. The completed observations
  ### zhat are what the update of beta of the EM algorithm uses.
  Zmat <- t(phi$A)                                      # d x p
  Xb <- matrix(0, d, n)
  bvec <- as.numeric(phi$beta)
  for (j in seq_len(dat$r)) Xb <- Xb + dat$covariates[, j, ] * bvec[j]
  sig <- Zmat %*% t(Msm)                                # d x n signal
  mu  <- Xb + sig
  Res <- t(t(as.matrix(dat$z)) - mu)                    # n x d residuals
  zhat <- matrix(NA_real_, n, d)

  cidx <- which(dat$obs$complete)
  if (length(cidx)) {
    BB <- Zmat %*% Reduce(`+`, Csm[cidx]) %*% t(Zmat) + crossprod(Res[cidx, , drop = FALSE])
    zhat[cidx, ] <- as.matrix(dat$z)[cidx, , drop = FALSE]
  } else {
    BB <- matrix(0, d, d)
  }
  blocks <- stem_blocks_cache(fs$ss$Vmat, d, regularization)
  for (tt in which(!dat$obs$complete)) {
    oi <- dat$obs$idx[[tt]]
    bl <- blocks(dat$obs$key[tt], oi)
    if (length(oi) == 0L) {
      BB <- BB + bl$Omega
      zhat[tt, ] <- mu[, tt]
    } else {
      SZo  <- bl$S %*% Zmat[oi, , drop = FALSE]
      rhat <- bl$S %*% matrix(Res[tt, oi], ncol = 1)
      BB <- BB + SZo %*% Csm[[tt]] %*% t(SZo) + bl$Omega + tcrossprod(rhat)
      zhat[tt, ] <- mu[, tt] + rhat
    }
  }

  list(loglik = fs$loglik, Sigmastar = fs$Sigmastar, y.smoothed = fs$smo$m,
       Msm = Msm, y0 = y0, P0 = P0, S11 = S11, S10 = S10, S00 = S00,
       BB = BB, zhat = zhat, sig = sig)
}
