#' @keywords internal
#' @noRd

### ---------------------------------------------------------------------------
### The estimation algorithms: the wrapper
### ---------------------------------------------------------------------------
### Four algorithms, chosen by `algorithm` in STEM_control(), all on the
### state-space form of the model: every iteration runs the Kalman filter and
### smoother of kalman.R inside its E-step.
###
###   "EM"            the EM algorithm: E-step (estep.R), M-step (mstep-em.R);
###   "ECME"          the same E-step, with the M-step of mstep-ecme.R, in which
###                   beta and m0 maximize the observed log-likelihood;
###   "SQUAREM"       the EM iterations accelerated by SQUAREM (squarem.R);
###   "SQUAREM-ECME"  the ECME iterations accelerated by SQUAREM.
###
### The four reach the same maximum; they differ in how many iterations they
### take and in how close to it a given stopping rule leaves them. An iteration
### is a "map": the parameters in, the parameters out, and the log-likelihood at
### the parameters in. The drivers below run a map until the stopping rule.
###
### SQUAREM needs an objective that the iterations increase, to safeguard its
### extrapolation. Without a penalty it is the log-likelihood. With a ridge on
### beta (lambda > 0, alpha = 0) it is the penalized log-likelihood,
###
###   loglik(psi) - lambda/2 sum_j w_j s_j(psi)^2 beta_j^2 ,
###
### with the scale of the penalty, s_j^2 = [sum_t X_t' Sigma_e^{-1} X_t]_jj,
### measured at the point psi itself (stem_ridge_merit()): the scale moves with
### Sigma_e, so there is no fixed objective, but along the plain iterations this
### one decreases by at most 1e-5, and the accelerated iterations reach the
### fixed point of the plain ones (to 1e-7 on dynamic designs with collinear
### lags, in about half the time; CHANGELOG, 2026-10-05). With the lasso or the
### elastic net (alpha > 0) the fits run the plain iterations of their map.

`stem_em_fit` <- function(phi, dat, opt, ctl) {
  map <- switch(ctl$algorithm,
                "EM" = , "SQUAREM" = stem_map_em,
                "ECME" = , "SQUAREM-ECME" = stem_map_ecme)
  penalized <- opt$lambda > 0
  accelerate <- ctl$algorithm %in% c("SQUAREM", "SQUAREM-ECME") && (!penalized || opt$alpha <= 0)
  if (accelerate && penalized) {
    stem_iterate_squarem(stem_ridge_merit(map), phi, dat, opt, ctl)
  } else if (accelerate) {
    stem_iterate_squarem(map, phi, dat, opt, ctl)
  } else {
    stem_iterate_plain(map, phi, dat, opt, ctl)
  }
}

### A map whose value also carries the penalized log-likelihood of a ridge at
### its input, `merit`, which the safeguard of SQUAREM compares. The scale of
### the penalty is the GLS one of the M-step (stem_beta_update(), standardize =
### TRUE), here at the input parameters.
`stem_ridge_merit` <- function(map) {
  function(phi, dat, opt) {
    out <- map(phi, dat, opt)
    Sstar <- dat$cov.spat(d = dat$d, logb = phi$logb, logtheta = phi$logtheta, dist = dat$dist)
    Sei <- tryCatch(chol2inv(chol(phi$sigma2omega * Sstar)), error = function(e) NULL)
    if (is.null(Sei)) return(out)
    Xm <- matrix(dat$covariates, nrow = dat$d)                   # d x (r n): r fastest
    s2 <- rowSums(matrix(colSums(Xm * (Sei %*% Xm)), nrow = dat$r))
    out$merit <- out$loglik - opt$lambda / 2 * sum(opt$pen_w * s2 * as.numeric(phi$beta)^2)
    out
  }
}

### One iteration of each algorithm.
`stem_map_em` <- function(phi, dat, opt) {
  est <- stem_estep(phi, dat, opt$regularization)
  out <- stem_mstep_em(est, phi, dat, opt)
  out$loglik <- est$loglik
  out
}

`stem_map_ecme` <- function(phi, dat, opt) {
  est <- stem_estep(phi, dat, opt$regularization)
  out <- stem_mstep_ecme(est, phi, dat, opt)
  out$loglik <- est$loglik
  out
}

### The plain iterations of a map, until the stopping rule or em_maxit.
### Returns the last iterate, the number of iterations, the Newton-Raphson
### iterations of each, the stopping criteria at the end and the last value of
### the map (which carries what the update of beta reported).
`stem_iterate_plain` <- function(map, phi, dat, opt, ctl) {
  iterNR <- integer(0)
  crit <- stem_em_check(NULL)
  ll_old <- NA_real_
  last <- NULL
  for (i in seq_len(ctl$em_maxit)) {
    if (isTRUE(opt$verbose)) message("**** EM Algorithm - iteration n. ", i)
    out <- map(phi, dat, opt)
    ### An iteration that returns non-finite parameters (which the
    ### Newton-Raphson step can produce on small or nearly collinear sets of
    ### locations) is discarded and the last valid iterate returned, so that a
    ### clusterwise algorithm calling this on many candidate subsets degrades
    ### gracefully instead of failing.
    if (!stem_par_finite(out$phi) || !is.finite(out$loglik)) {
      if (i == 1L) {
        stop("The EM algorithm produced non-finite parameters at the first iteration: check the starting values in 'phi' and the conditioning of the data.", call. = FALSE)
      }
      warning("The EM algorithm produced non-finite parameters at iteration ", i,
              "; the last valid iterate is returned.", call. = FALSE)
      break
    }
    iterNR[i] <- out$n_iter_NR
    if (i > 1L) crit <- stem_em_check(out$phi, phi, out$loglik, ll_old, ctl)
    phi <- out$phi
    ll_old <- out$loglik
    last <- out
    if (crit$done) break
  }
  list(phi = phi, iter = length(iterNR), iterNR = iterNR, crit = crit, last = last)
}

### The stopping rule. Two criteria: the largest relative change of a free
### parameter, with a floor of 1e-3 on the denominator for parameters near
### zero, and the absolute change of the log-likelihood. The loading matrix and
### C0, which are not estimated, are not free. With em_stop = "any" either
### criterion stops the iterations, with "all" both are needed. Called with
### NULL, it returns the state before any comparison is possible.
`stem_em_free` <- function(phi) {
  c(phi$sigma2omega, phi$logtheta, phi$logb, as.numeric(phi$beta), as.numeric(phi$G),
    as.numeric(phi$Sigmaeta), as.numeric(phi$m0))
}

`stem_em_check` <- function(phi_new, phi_old = NULL, ll_new = NA, ll_old = NA, ctl = NULL) {
  if (is.null(phi_new)) {
    return(list(max_rel_par = Inf, delta_loglik = Inf, conv_par = FALSE, conv_log = FALSE,
                done = FALSE))
  }
  a <- stem_em_free(phi_new)
  b <- stem_em_free(phi_old)
  max_rel <- max(abs(a - b) / pmax(abs(b), 1e-3))
  dll <- abs(ll_new - ll_old)
  conv_par <- isTRUE(max_rel < ctl$em_tol_par)
  conv_log <- isTRUE(dll < ctl$em_tol_loglik)
  done <- if (identical(ctl$em_stop, "all")) conv_par && conv_log else conv_par || conv_log
  list(max_rel_par = max_rel, delta_loglik = dll, conv_par = conv_par, conv_log = conv_log,
       done = done)
}

`stem_par_finite` <- function(phi) all(is.finite(stem_em_free(phi)))

### The data of a fit, in the form the E- and M-steps read: built once per fit.
`stem_em_data` <- function(StemModel, distance, cov.spat) {
  d <- StemModel$data$d
  n <- StemModel$data$n
  r <- StemModel$data$r
  coordinates <- StemModel$data$coordinates
  dist <- if (distance == "geo") {
    as.matrix(geodist::geodist(coordinates, measure = "geodesic"))
  } else {
    as.matrix(stats::dist(coordinates, diag = TRUE))
  }
  z <- StemModel$data$z
  list(z = stats::ts(z), zmat = as.matrix(z),
       covariates = changedimension_covariates(StemModel$data$covariates, d = d, r = r, n = n),
       coordinates = coordinates, dist = dist, obs = stem_obs_index(z),
       n = n, d = d, r = r, p = StemModel$skeleton$p, cov.spat = cov.spat)
}
