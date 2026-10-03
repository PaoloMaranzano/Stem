#' @keywords internal
#' @noRd

### ---------------------------------------------------------------------------
### The E-step, shared by every estimation algorithm
### ---------------------------------------------------------------------------
### kalman() (R/kalman.R) gives, at the current parameters, the filtered and
### smoothed moments of the latent states and the log-likelihood. From them the
### E-step forms the expected sufficient statistics that the M-step of either
### algorithm needs:
###
###   S11 = sum_t E(y_t y_t' | z),  S00 = sum_t E(y_{t-1} y_{t-1}' | z),
###   S10 = sum_t E(y_t y_{t-1}' | z),
###   BB  = sum_t E(e_t e_t' | z),  e_t = z_t - X_t beta - A' y_t,
###   zhat, the observations completed by their conditional expectations,
###
### together with the smoothed moments of the initial state, y0 and P0.
###
###   regularization  a ridge added to the matrices inverted (the gains of the
###                   smoother at times 0 and n and the blocks of the missing
###                   values); 0 by default, see STEM_Estimation()
`stem_estep` <- function(phi, dat, regularization = 0) {

  ks  <- kalman(phi, dat, regularization)
  n   <- dat$n
  d   <- dat$d
  Msm <- as.matrix(ks$smoothed$m)      # smoothed means, n x p
  Csm <- ks$smoothed$C                 # smoothed variances, a list of n
  y0  <- ks$smoothed$y0
  P0  <- ks$smoothed$P0

  ### The three sums of outer products of the smoothed states, each a
  ### crossprod of the matrix of the smoothed means plus a sum of variances.
  lag <- seq_len(n - 1L)
  S11 <- crossprod(Msm) + Reduce(`+`, Csm)
  S00 <- y0 %*% t(y0) + P0 + crossprod(Msm[lag, , drop = FALSE]) + Reduce(`+`, Csm[lag])
  S10 <- matrix(Msm[1L, ], ncol = 1) %*% t(y0) +
    crossprod(Msm[-1L, , drop = FALSE], Msm[lag, , drop = FALSE]) + Reduce(`+`, ks$smoothed$lag1)

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
  blocks <- stem_blocks_cache(ks$ss$Vmat, d, regularization)
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

  list(loglik = ks$loglik, Sigmastar = ks$Sigmastar, y.smoothed = ks$smoothed$m,
       Msm = Msm, y0 = y0, P0 = P0, S11 = S11, S10 = S10, S00 = S00,
       BB = BB, zhat = zhat, sig = sig)
}
