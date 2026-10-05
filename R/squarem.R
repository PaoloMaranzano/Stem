#' @keywords internal
#' @noRd

### ---------------------------------------------------------------------------
### SQUAREM: acceleration of the EM (or ECME) iterations
### ---------------------------------------------------------------------------
### The EM algorithm is a fixed-point iteration psi <- F(psi), which converges
### linearly: along its slowest direction the distance from the maximum shrinks
### by a factor lambda close to one at every iteration. SQUAREM (Roland and
### Varadhan 2005; Varadhan and Roland 2008, scheme S3) extrapolates from two
### steps of the iteration,
###
###   r = F(psi) - psi,   w = F(F(psi)) - 2 F(psi) + psi,   alpha = ||r|| / ||w||,
###   psi' = psi + 2 alpha r + alpha^2 w ,
###
### and then takes one ordinary step from psi', which stabilizes the
### extrapolation. For a linear iteration in one direction, r = (lambda - 1) e
### and w = (lambda - 1)^2 e with e = psi - psi_hat, so alpha = 1/(1 - lambda)
### and psi' is the maximum itself; in several directions the extrapolation
### removes most of the slowest. With alpha = 1, psi' = F(F(psi)): two ordinary
### steps.
###
### The steplength and its safeguard follow the R package SQUAREM (Du and
### Varadhan 2020). alpha is capped by step_max, which starts at 1 and grows by
### a factor 4 every time the cap binds. The extrapolation is tried once per
### cycle and accepted only if the log-likelihood at psi' is not below that at
### psi; otherwise the cycle ends at F(F(psi)), whose increase the EM algorithm
### guarantees, and step_max shrinks. Every call of the iteration evaluates the
### log-likelihood at its input, so the safeguard costs nothing. With a ridge
### the map also returns the penalized log-likelihood, `merit`
### (stem_ridge_merit() in em-fit.R), and the safeguard compares that. The
### extrapolation is done on a scale on which every value is admissible: the
### logarithm of the variances (the Cholesky factor of a full Sigmaeta) and of
### the range.
###
### The stopping rule is the one of the ordinary iterations (stem_em_check()),
### applied at the start of every cycle to the step F(psi) - psi; the number of
### iterations reported is the number of calls of the iteration.

### The free parameters as one unconstrained vector, and back.
`stem_par_vec` <- function(phi, opt) {
  S <- phi$Sigmaeta
  s <- if (opt$Sigmaetadiag) {
    log(diag(S))
  } else {
    L <- chol(S)
    c(log(diag(L)), L[upper.tri(L)])
  }
  c(as.numeric(phi$beta), as.numeric(phi$m0), as.numeric(phi$G), s,
    log(phi$sigma2omega), phi$logb, phi$logtheta)
}

`stem_par_unvec` <- function(v, phi, opt) {
  r <- length(phi$beta)
  p <- nrow(phi$G)
  k <- 0L
  take <- function(m) {
    x <- v[k + seq_len(m)]
    k <<- k + m
    x
  }
  phi$beta <- matrix(take(r), ncol = 1)
  phi$m0 <- matrix(take(p), p, 1)
  phi$G <- matrix(take(p * p), p, p)
  if (opt$Sigmaetadiag) {
    phi$Sigmaeta <- diag(exp(take(p)), p)
  } else {
    L <- matrix(0, p, p)
    diag(L) <- exp(take(p))
    L[upper.tri(L)] <- take(p * (p - 1L) / 2L)
    phi$Sigmaeta <- crossprod(L)
  }
  phi$sigma2omega <- exp(take(1L))
  phi$logb <- take(1L)
  ### the range stays within the limits the distances can identify
  phi$logtheta <- min(max(take(1L), opt$logtheta_lim[1]), opt$logtheta_lim[2])
  phi
}

### What the safeguard compares: the penalized log-likelihood of a ridge when
### the map returns it, the log-likelihood otherwise.
`stem_merit` <- function(o) if (is.null(o$merit)) o$loglik else o$merit

### The accelerated iterations. Same arguments and value as stem_iterate_plain().
`stem_iterate_squarem` <- function(map, phi, dat, opt, ctl) {
  maxit <- ctl$em_maxit
  calls <- 0L
  iterNR <- integer(0)
  crit <- stem_em_check(NULL)
  ll_old <- NA_real_
  last <- NULL
  step_max <- 1
  grow <- 4

  ### one call of the iteration: its value, or NULL when it fails or returns
  ### anything that is not finite
  step <- function(ph) {
    calls <<- calls + 1L
    if (isTRUE(opt$verbose)) message("**** EM Algorithm (SQUAREM) - iteration n. ", calls)
    out <- tryCatch(map(ph, dat, opt), error = function(e) NULL)
    if (is.null(out) || !stem_par_finite(out$phi) || !is.finite(out$loglik)) return(NULL)
    iterNR[calls] <<- out$n_iter_NR
    out
  }
  failed <- function() {
    if (calls == 1L) {
      stop("The EM algorithm produced non-finite parameters at the first iteration: check the starting values in 'phi' and the conditioning of the data.", call. = FALSE)
    }
    warning("The EM algorithm produced non-finite parameters at iteration ", calls,
            "; the last valid iterate is returned.", call. = FALSE)
  }

  cycle <- 0L
  repeat {
    cycle <- cycle + 1L
    phi0 <- phi
    o1 <- step(phi0)
    if (is.null(o1)) { failed(); break }
    if (cycle > 1L) crit <- stem_em_check(o1$phi, phi0, o1$loglik, ll_old, ctl)
    ll_old <- o1$loglik
    phi <- o1$phi
    last <- o1
    if (crit$done || calls >= maxit) break

    o2 <- step(o1$phi)
    if (is.null(o2)) { failed(); break }
    phi <- o2$phi
    last <- o2
    if (calls >= maxit) break

    ### the extrapolation, on the unconstrained scale
    v <- tryCatch(list(stem_par_vec(phi0, opt), stem_par_vec(o1$phi, opt),
                       stem_par_vec(o2$phi, opt)), error = function(e) NULL)
    if (is.null(v) || !all(is.finite(unlist(v)))) next
    r <- v[[2]] - v[[1]]
    w <- v[[3]] - 2 * v[[2]] + v[[1]]
    alpha <- sqrt(sum(r^2) / sum(w^2))
    if (!is.finite(alpha)) alpha <- 1
    alpha <- max(1, min(step_max, alpha))
    if (abs(alpha - 1) > 0.01) {
      o3 <- step(stem_par_unvec(v[[1]] + 2 * alpha * r + alpha^2 * w, phi0, opt))
      if (!is.null(o3) && stem_merit(o3) >= stem_merit(o1) - 1e-6) {
        phi <- o3$phi
        last <- o3
      } else {
        if (alpha == step_max) step_max <- max(1, step_max / grow)
        alpha <- 1
      }
    }
    if (alpha == step_max) step_max <- grow * step_max
    if (calls >= maxit) break
  }

  list(phi = phi, iter = calls, iterNR = iterNR, crit = crit, last = last)
}
