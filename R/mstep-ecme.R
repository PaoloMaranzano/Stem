#' @keywords internal
#' @noRd

### ---------------------------------------------------------------------------
### The M-step of the ECME algorithm
### ---------------------------------------------------------------------------
### ECME, "expectation / conditional maximization either" (Liu and Rubin
### 1994), lets some of the conditional maximizations of the M-step maximize
### the observed log-likelihood itself instead of the expected complete-data
### log-likelihood Q. Here the latent process, sigma2omega and the spatial
### parameters are updated on Q exactly as in the EM algorithm, and then beta
### and m0 are updated on the observed log-likelihood, given all the others.
### The step on the likelihood comes last, after the steps on Q, which is what
### keeps the algorithm monotone (Liu and Rubin 1994).
###
### Why it helps. Given the other parameters, the observed log-likelihood
### depends on (beta, m0) only through the mean of the data,
### E(z_t) = X_t beta + A' G^t m0, so it is a quadratic in them and its
### maximizer is a generalized least-squares estimator under the covariance the
### model implies, latent process included. The EM algorithm instead updates
### beta and m0 given the smoothed latent states, which absorb part of any error
### in the level: along the direction in which the intercept and m0 move
### together it corrects only a small part of the error at each iteration. The
### step below takes that direction in one go.
`stem_mstep_ecme` <- function(est, phi, dat, opt) {
  lat <- stem_update_latent(est, phi, dat, opt)
  s2w <- stem_update_sigma2omega(est, dat, opt)
  spa <- stem_update_spatial(est, phi, dat, s2w, opt)

  phi$sigma2omega <- s2w
  phi$logtheta <- spa$logtheta
  phi$logb <- spa$logb
  phi$G <- lat$G
  phi$Sigmaeta <- lat$Sigmaeta

  mean <- stem_update_mean_observed(phi, dat, opt)
  phi$beta <- matrix(mean$beta, ncol = 1)
  phi$m0 <- mean$m0
  list(phi = phi, n_iter_NR = spa$n_iter,
       beta_df = mean$df, lambda_ref = mean$lambda_ref, lambda_eff = mean$lambda_eff)
}

### beta and m0 maximizing the observed log-likelihood given the other
### parameters of `phi`. Without a latent process m0 plays no role and is held.
### With a penalty on beta, m0 is profiled out and the quadratic left in beta
### goes to stem_beta_update(), as in the EM algorithm.
`stem_update_mean_observed` <- function(phi, dat, opt) {
  r <- dat$r
  p <- dat$p
  ib <- seq_len(r)
  q <- stem_gls_mean(phi, dat, with_m0 = opt$latent)
  if (opt$latent) {
    im <- r + seq_len(p)
    Mmm <- q$M[im, im, drop = FALSE]
    Hb <- solve(Mmm, q$M[im, ib, drop = FALSE])          # Mmm^{-1} M_m,beta
    hv <- solve(Mmm, q$v[im])                            # Mmm^{-1} v_m
    M <- q$M[ib, ib, drop = FALSE] - q$M[ib, im, drop = FALSE] %*% Hb
    v <- q$v[ib] - q$M[ib, im, drop = FALSE] %*% hv
  } else {
    M <- q$M
    v <- q$v
  }
  bet <- stem_beta_update(M = M, v = v, alpha = opt$alpha, lambda = opt$lambda,
                          w = opt$pen_w, beta0 = phi$beta, lambda_scale = opt$lambda_scale,
                          ridge_reg = opt$regularization)
  m0 <- if (opt$latent) phi$m0 + matrix(hv - Hb %*% bet$beta, p, 1) else phi$m0
  list(beta = bet$beta, m0 = m0, df = bet$df,
       lambda_ref = bet$lambda_ref, lambda_eff = bet$lambda_eff)
}

### The quadratic of the observed log-likelihood in delta = (beta, dm0), dm0 a
### change of m0 from its current value, by the augmented Kalman filter (de Jong
### 1991; Durbin and Koopman 2012, Ch. 6). The filter of the current parameters
### is run on every column of
###
###   [ z_t, X_t, A' G^t ]          (d x (1 + r + p); the last block with_m0 only),
###
### the first with the current m0 as initial mean, the others with zero. The
### gains do not depend on the data, so every column shares them, and the
### filter is linear in its data and in its initial mean: the innovations of
### z_t at delta are e_t - E_t delta, with [e_t, E_t] the innovations of the
### columns. Hence, with F_t the innovation covariance,
###
###   loglik(delta) = c - 1/2 sum_t (e_t - E_t delta)' F_t^{-1} (e_t - E_t delta),
###
### maximized at M^{-1} v with M = sum_t E_t' F_t^{-1} E_t and
### v = sum_t E_t' F_t^{-1} e_t, which is what is returned. F_t^{-1} is applied
### through the Woodbury identity and missing rows are dropped, as in
### filtering() and filterstep().
`stem_gls_mean` <- function(phi, dat, with_m0 = TRUE) {
  n <- dat$n; p <- dat$p; r <- dat$r
  G <- phi$G
  W <- phi$Sigmaeta
  V <- phi$sigma2omega *
    dat$cov.spat(d = dat$d, logb = phi$logb, logtheta = phi$logtheta, dist = dat$dist)
  Z <- t(phi$A)
  z <- dat$zmat
  obs <- dat$obs

  ### constants of the measurement equation, by missingness pattern
  store <- new.env(parent = emptyenv())
  pattern <- function(key, oi) {
    kk <- if (is.na(key)) "*complete*" else key
    hit <- store[[kk]]
    if (is.null(hit)) {
      Zo <- Z[oi, , drop = FALSE]
      Lc <- chol(V[oi, oi, drop = FALSE])
      VinvZ <- backsolve(Lc, forwardsolve(t(Lc), Zo))
      hit <- list(Zo = Zo, Lc = Lc, VinvZ = VinvZ, FU = crossprod(Zo, VinvZ))
      assign(kk, hit, envir = store)
    }
    hit
  }

  ncol <- 1L + r + if (with_m0) p else 0L
  a_f <- matrix(0, p, ncol)              # filtered means of the columns
  a_f[, 1] <- phi$m0
  P_f <- phi$C0                           # filtered variance, common to all
  Gt <- diag(p)
  S <- matrix(0, ncol, ncol)
  for (tt in seq_len(n)) {
    Gt <- G %*% Gt
    a <- G %*% a_f
    R <- G %*% P_f %*% t(G) + W
    oi <- obs$idx[[tt]]
    if (length(oi) == 0L) {
      a_f <- a
      P_f <- R
      next
    }
    cst <- pattern(obs$key[tt], oi)
    Y <- cbind(z[tt, oi], matrix(dat$covariates[oi, , tt], nrow = length(oi)),
               if (with_m0) cst$Zo %*% Gt)
    E <- Y - cst$Zo %*% a
    Mm <- solve(solve(R) + cst$FU)
    Ue <- crossprod(cst$VinvZ, E)                       # Z' V^{-1} E
    a_f <- a + R %*% (Ue - cst$FU %*% Mm %*% Ue)
    P_f <- R - R %*% (cst$FU - cst$FU %*% Mm %*% cst$FU) %*% R
    Le <- forwardsolve(t(cst$Lc), E)
    S <- S + crossprod(Le) - crossprod(Ue, Mm %*% Ue)   # E' F^{-1} E
  }
  list(M = S[-1, -1, drop = FALSE], v = S[-1, 1])
}
