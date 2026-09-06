#' @keywords internal
#' @noRd
#'
### ---------------------------------------------------------------------------
### A faster forward pass, algebraically identical to filtering()
### ---------------------------------------------------------------------------
### The straightforward implementation forms and inverts the d x d matrix
###
###     Q_t = Z R_t Z' + V ,        Z = t(Fmat),  V = Sigma_e ,
###
### at every time point, which costs O(T d^3) and dominates everything else once
### the network has more than a few dozen locations. But V does not depend on t,
### and Z R_t Z' has rank p, with p typically 1. The Woodbury identity therefore
### applies:
###
###     Q_t^{-1} = V^{-1} - U (R_t^{-1} + Z'V^{-1}Z)^{-1} U' ,   U = V^{-1} Z ,
###
### where U and Z'V^{-1}Z are computed ONCE per EM iteration, and the matrix
### determinant lemma gives
###
###     log|Q_t| = log|V| + log|R_t| + log|R_t^{-1} + Z'V^{-1}Z| .
###
### None of the quantities the recursion needs -- the gain applied to the
### innovation, the updated state variance, the quadratic form of the
### log-likelihood -- requires Q_t itself. They reduce to p x p algebra plus one
### triangular solve against the cached Cholesky factor of V, so the per-step
### cost falls from O(d^3) to O(d^2).
###
### Two further differences from the straightforward version, both of which
### matter at scale. The log-density is evaluated directly rather than through
### mvtnorm::dmvnorm, which re-checks the symmetry of Q with all.equal at every
### step -- on a 36-station network that check alone is about a seventh of the
### total running time. And the density is never formed on its natural scale, so
### it cannot underflow.
###
### Missing values are handled as in filterstep(), by restricting the
### measurement equation to the observed rows (Durbin and Koopman 2012,
### Sect. 4.10). The constants above then depend on the missingness pattern
### rather than being global, so they are computed once per DISTINCT pattern and
### cached; in a monitoring network the number of distinct patterns is far
### smaller than the number of time points.
###
### Arguments
###   ss   the state-space list assembled by kalman(), unchanged
###
### Value
###   ss with m, C and loglik filled in, exactly as filtering() returns them.

`filtering_fast` <- function(ss) {

  n <- ss$n; p <- ss$p; d <- ss$d
  Gm <- ss$Gmat; Wm <- ss$Wmat
  Zfull <- t(ss$Fmat)                       # d x p
  V     <- ss$Vmat                          # d x d, constant in t
  flag  <- isTRUE(ss$flag.cov)
  beta  <- ss$beta
  zmat  <- as.matrix(ss$z)

  obs <- stem_obs_index(zmat)

  ### Constants of the measurement equation, cached by missingness pattern.
  store <- new.env(parent = emptyenv())
  pattern_const <- function(key, oi) {
    kk <- if (is.na(key)) "*full*" else key
    hit <- store[[kk]]
    if (!is.null(hit)) return(hit)
    Zo <- Zfull[oi, , drop = FALSE]
    Vo <- V[oi, oi, drop = FALSE]
    Lc <- chol(Vo)                                  # Vo = t(Lc) %*% Lc
    VinvZ <- backsolve(Lc, forwardsolve(t(Lc), Zo))
    hit <- list(Lc = Lc, VinvZ = VinvZ,
                FU = crossprod(Zo, VinvZ),
                logdetV = 2 * sum(log(diag(Lc))),
                Zo = Zo, nobs = length(oi))
    assign(kk, hit, envir = store)
    hit
  }

  m <- matrix(NA_real_, n, p)
  C <- vector("list", n)
  loglik <- 0
  mprev <- t(ss$m0)          # 1 x p
  Cprev <- ss$C0

  for (tt in seq_len(n)) {

    a  <- Gm %*% t(mprev)                       # p x 1
    Rt <- Gm %*% Cprev %*% t(Gm) + Wm           # p x p

    oi <- obs$idx[[tt]]
    if (length(oi) == 0L) {
      ### nothing observed: no update and no likelihood contribution
      m[tt, ] <- a
      C[[tt]] <- Rt
      mprev <- matrix(a, nrow = 1); Cprev <- Rt
      next
    }

    cn <- pattern_const(obs$key[tt], oi)
    Ri <- solve(Rt)
    Mm <- solve(Ri + cn$FU)                     # p x p
    ZQiZ <- cn$FU - cn$FU %*% Mm %*% cn$FU      # Z' Q^{-1} Z

    zt <- matrix(zmat[tt, oi], ncol = 1)
    mu <- cn$Zo %*% a
    ## XXX is d x r x n; subsetting the third margin would leave a 3-d array,
    ## so the design block is rebuilt as a matrix explicitly
    if (flag) {
      Xt <- matrix(ss$XXX[oi, , tt], nrow = length(oi))
      mu <- mu + Xt %*% beta
    }
    e <- zt - mu

    Ue   <- crossprod(cn$VinvZ, e)              # Z' V^{-1} e, p x 1
    ZQie <- Ue - cn$FU %*% Mm %*% Ue            # Z' Q^{-1} e
    mt   <- a + Rt %*% ZQie
    Ct   <- Rt - Rt %*% ZQiZ %*% Rt

    Le   <- forwardsolve(t(cn$Lc), e)
    quad <- sum(Le^2) - as.numeric(crossprod(Ue, Mm %*% Ue))
    ldQ  <- cn$logdetV + log(det(Rt)) + log(det(Ri + cn$FU))
    loglik <- loglik - 0.5 * (cn$nobs * log(2 * pi) + ldQ + quad)

    m[tt, ] <- mt
    C[[tt]] <- Ct
    mprev <- matrix(mt, nrow = 1); Cprev <- Ct
  }

  if (stats::is.ts(ss$z)) {
    ss$m <- stats::ts(m, stats::start(ss$z), end = stats::end(ss$z),
                      frequency = stats::frequency(ss$z))
  } else {
    ss$m <- m
  }
  ss$C <- C
  ss$loglik <- loglik
  ss
}
