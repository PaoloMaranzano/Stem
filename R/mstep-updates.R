#' @keywords internal
#' @noRd

### ---------------------------------------------------------------------------
### The conditional updates of the M-step
### ---------------------------------------------------------------------------
### The M-step is not a joint maximization: it is a sequence of conditional
### maximizations of the expected complete-data log-likelihood Q (the ECM
### algorithm of Meng and Rubin 1993), each over some parameters given the
### others. The updates below are its building blocks; mstep-em.R and
### mstep-ecme.R put them together. Every one takes the output `est` of
### stem_estep(), the current parameters `phi`, the data `dat` and the settings
### `opt` of the fit (see stem_em_settings()).

### The latent process: G and Sigmaeta from the smoothed moments, and m0 as the
### smoothed mean of the initial state. With a diagonal Sigmaeta its update
### uses the current G; with a diagonal G its update uses the new Sigmaeta.
### Without a latent process (latent = FALSE) nothing in the data speaks about
### the three, which are held at their current values.
`stem_update_latent` <- function(est, phi, dat, opt) {
  if (!opt$latent) return(list(G = phi$G, Sigmaeta = phi$Sigmaeta, m0 = phi$m0))
  n <- dat$n
  p <- dat$p
  S11 <- est$S11; S10 <- est$S10; S00 <- est$S00

  if (opt$Sigmaetadiag) {
    G <- phi$G
    num <- S11 - S10 %*% t(G) - G %*% t(S10) + G %*% S00 %*% t(G)
    Sigmaeta <- diag(diag(num) / n, p)
  } else {
    Sigmaeta <- (S11 - S10 %*% solve(S00) %*% t(S10)) / n
  }

  if (opt$Gdiag) {
    Wi <- solve(diag(opt$regularization, p) + Sigmaeta)
    G <- diag(diag(Wi %*% S10) / diag(Wi %*% S00), p)
  } else {
    G <- S10 %*% solve(S00)
  }

  list(G = G, Sigmaeta = Sigmaeta, m0 = est$y0)
}

### The scale of the measurement error, given the spatial correlation of the
### current parameters: sigma2omega = tr(Sigma*^{-1} BB) / (n d). The divisor is
### the complete-data count, because BB already carries the conditional
### variance of whatever was not observed.
`stem_update_sigma2omega` <- function(est, dat, opt) {
  D <- solve(diag(opt$regularization, dat$d) + est$Sigmastar, est$BB)
  sum(diag(D)) / (dat$n * dat$d)
}

### The regression coefficients given the smoothed latent states: the GLS of
### the completed observations net of the latent signal on the covariates,
### under Sigma_e = sigma2omega Sigma*. The two accumulations
###
###   M = sum_t X_t' Sigma_e^{-1} X_t,   v = sum_t X_t' Sigma_e^{-1} (zhat_t - A' y_t)
###
### are crossprods of the design blocks stacked by period. stem_beta_update()
### returns their maximizer, ordinary or penalized (R/stem-penalty.R).
`stem_update_beta` <- function(est, phi, dat, sigma2omega, opt) {
  d <- dat$d; n <- dat$n; r <- dat$r
  Sei  <- solve(diag(opt$regularization, d) + sigma2omega * est$Sigmastar)
  SX   <- Sei %*% matrix(dat$covariates, nrow = d)                 # d x (r n)
  Xbig <- matrix(aperm(dat$covariates, c(1, 3, 2)), d * n, r)
  SXbig <- matrix(aperm(array(SX, c(d, r, n)), c(1, 3, 2)), d * n, r)
  M <- crossprod(Xbig, SXbig)
  v <- crossprod(Xbig, as.vector(Sei %*% (t(est$zhat) - est$sig)))
  if (det(M) < 1e-7) {
    warning("Error in beta estimation! The matrix can not be inverted!!!!", call. = FALSE)
  }
  stem_beta_update(M = M, v = v, alpha = opt$alpha, lambda = opt$lambda,
                   w = opt$pen_w, beta0 = phi$beta, lambda_scale = opt$lambda_scale,
                   ridge_reg = opt$regularization)
}

### The spatial parameters (log theta, log b) given sigma2omega, by
### Newton-Raphson on the part of Q that depends on them,
###
###   n log|sigma2omega Sigma*| + tr((sigma2omega Sigma*)^{-1} BB) .
###
### When the Hessian is not negative definite enough, a grid search over the
### two parameters, from 1/100 to 10 times their current values, moves the
### point before the step is tried again. A step that is singular or not
### finite ends the iterations at the last valid point. With spatial = FALSE
### the correlation is the identity and nothing is estimated.
`stem_update_spatial` <- function(est, phi, dat, sigma2omega, opt) {
  logtheta <- phi$logtheta
  logb <- phi$logb
  if (!opt$spatial) return(list(logtheta = logtheta, logb = logb, n_iter = 0L))

  n <- dat$n; d <- dat$d; dist <- dat$dist; BB <- est$BB
  cov.spat <- dat$cov.spat
  takes_E <- "E" %in% names(formals(cov.spat))
  Qpart <- function(lt, lb) {
    Q_function_addendo1(sigma2omega = sigma2omega, n = n,
                        Sigmastar = cov.spat(d = d, logb = lb, logtheta = lt, dist = dist), B = BB)
  }

  it <- 0L
  converged <- FALSE
  while (!converged && it < opt$nr_maxit) {
    it <- it + 1L

    hess_ok <- FALSE
    attempt <- 0L
    while (!hess_ok && attempt < opt$nr_hess_maxit) {
      attempt <- attempt + 1L
      ### the kernel exp(-theta h), shared by the correlation and its two
      ### derivatives in log(theta); the derivatives in log(b) are multiples of
      ### the identity, carried as single numbers (see stem_xprod())
      Ker <- exp(-exp(logtheta) * dist)
      cs_args <- list(logb = logb, d = d, logtheta = logtheta, dist = dist)
      if (takes_E) cs_args$E <- Ker
      X   <- do.call(cov.spat, cs_args)
      Xi  <- solve(X)
      XiB <- Xi %*% BB
      d1t <- d1_Sigmastar_logtheta.exp(logtheta = logtheta, dist = dist, E = Ker)
      d2t <- d2_Sigmastar_logtheta.exp(logtheta = logtheta, dist = dist, E = Ker)
      d1b <- d1_Sigmastar_logb.exp(logb = logb, d = d)
      d2b <- d2_Sigmastar_logb.exp(logb = logb, d = d)
      Pt  <- Xi %*% d1t
      Pb  <- d1b * Xi

      g_t  <- d1_Q(n = n, X = X, d1_X = d1t, sigma2omega = sigma2omega, B = BB,
                   Xi = Xi, XiB = XiB, P = Pt)
      g_b  <- d1_Q(n = n, X = X, d1_X = d1b, sigma2omega = sigma2omega, B = BB,
                   Xi = Xi, XiB = XiB, P = Pb)
      h_tt <- d2_Q(n = n, X = X, d1_X = d1t, d2_X = d2t, sigma2omega = sigma2omega, B = BB,
                   Xi = Xi, XiB = XiB, P = Pt, P2 = Xi %*% d2t)
      h_bb <- d2_Q(n = n, X = X, d1_X = d1b, d2_X = d2b, sigma2omega = sigma2omega, B = BB,
                   Xi = Xi, XiB = XiB, P = Pb, P2 = Pb)
      h_tb <- d12_Q(n = n, X = X, d1_X_theta = d1t, d1_X_logb = d1b, sigma2omega = sigma2omega,
                    B = BB, Xi = Xi, XiB = XiB, Pt = Pt, Pb = Pb)
      H <- matrix(c(h_tt, h_tb, h_tb, h_bb), 2, 2)

      hess_ok <- isTRUE(det(H) > 1e-3)
      if (!all(is.finite(H))) break
      if (!hess_ok) {
        ### the grid search: the point of the grid where the objective is
        ### lowest becomes the new current point
        lt_grid <- log(seq(0.01, 10, length = 10) * max(exp(logtheta), .Machine$double.xmin))
        lb_grid <- log(seq(0.01, 10, length = 10) * max(exp(logb), .Machine$double.xmin))
        QQ <- outer(seq_along(lt_grid), seq_along(lb_grid),
                    Vectorize(function(i, j) Qpart(lt_grid[i], lb_grid[j])))
        col <- which.min(apply(QQ, 2, min))
        logtheta <- lt_grid[apply(QQ, 2, which.min)[col]]
        logb <- lb_grid[col]
      }
    }

    old <- c(logtheta, logb)
    delta <- try(solve(H + diag(opt$regularization, 2), c(g_t, g_b)), silent = TRUE)
    if (inherits(delta, "try-error") || !all(is.finite(delta))) break
    new <- old - as.numeric(delta)
    if (!all(is.finite(new))) break
    converged <- isTRUE(sqrt(sum((new - old)^2)) / max(sqrt(sum(old^2)), .Machine$double.eps) <
                          opt$nr_tol)
    logtheta <- new[1]
    logb <- new[2]
    if (isTRUE(opt$verbose)) message("*** NR Algorithm - iteration n. ", it)
  }

  list(logtheta = logtheta, logb = logb, n_iter = it)
}
