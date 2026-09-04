#' @keywords internal
#' @noRd
#'
### ---------------------------------------------------------------------------
### Missing observations in the response
### ---------------------------------------------------------------------------
### The treatment follows Durbin and Koopman (2012, 2nd ed.), Sections 2.7 and
### 4.10. When only some of the elements of the observation vector z_t are
### observed, Section 4.10 (p. 111) prescribes a selection matrix W_t whose rows
### are a subset of the rows of the identity, and replaces the measurement
### equation by
###
###     z*_t = Z*_t y_t + e*_t,     e*_t ~ N(0, H*_t),
###
### with z*_t = W_t z_t, Z*_t = W_t Z_t and H*_t = W_t H_t W_t'. The Kalman
### filter and the smoother then proceed exactly as in the complete case, on an
### observation vector whose dimension varies over time. When every element is
### missing the update is skipped altogether, which Durbin and Koopman obtain by
### setting Z_t = 0; the same device leaves the backward smoothing recursions
### unchanged, which is why smoothing() and smootherstep() need no modification
### here: they read the filtered moments and the transition, never the data.
###
### What does need care is the M-step. The EM algorithm maximizes the expected
### COMPLETE-data log-likelihood, so the sufficient statistics have to be
### completed: the missing entries of z_t enter through their conditional
### expectation given everything observed, and their conditional variance is
### added back as a correction term. The two helpers below supply the pieces.
###
### Note that the conditional expectation of a missing element is NOT simply the
### signal Z_t yhat_t. That would be exact only for a diagonal H_t. Here H_t is
### the spatial covariance Sigma_e, which couples the locations, so the missing
### block of the measurement error is predicted from the observed block by the
### usual Gaussian conditioning -- the same algebra as kriging, at a fixed time
### point.


### Observed-row index of every time point.
###
### Returns a list of integer vectors, one per time point, giving the rows of
### the observation vector that are actually observed, together with a key
### identifying the missingness pattern so that the expensive per-pattern
### quantities can be cached. Patterns repeat heavily in practice -- a station
### is typically down for a stretch of consecutive days -- so the number of
### distinct patterns is far smaller than the number of time points.
###
### Arguments
###   z   n x d numeric matrix of observations, possibly with NA
`stem_obs_index` <- function(z) {
  z <- as.matrix(z)
  n <- nrow(z)
  ok <- !is.na(z)
  idx <- vector("list", n)
  key <- character(n)
  complete <- rep(TRUE, n)
  for (tt in seq_len(n)) {
    o <- which(ok[tt, ])
    idx[[tt]] <- o
    complete[tt] <- length(o) == ncol(z)
    ### The key must distinguish a complete time point from one where nothing
    ### is observed: pasting an empty index would give the same empty string
    ### for both, and the fully missing pattern would then be served the
    ### complete-data shortcut from the cache.
    key[tt] <- if (complete[tt]) NA_character_ else paste0("p", paste(o, collapse = "."))
  }
  list(idx = idx, key = key, complete = complete,
       any_missing = !all(complete),
       n_obs = vapply(idx, length, integer(1)))
}


### Conditional blocks of the measurement error for one missingness pattern.
###
### For a pattern with observed rows `oi` and missing rows `mi`, and a
### measurement covariance Sigma (the full d x d Sigma_e), returns
###
###   P     = Sigma[mi, oi] Sigma[oi, oi]^{-1}
###   Omega = Sigma[mi, mi] - P Sigma[oi, mi]
###
### so that, conditionally on the latent state and on the observed part of the
### measurement error e_obs,
###
###   e_mis | e_obs ~ N( P e_obs , Omega ).
###
### `S` is the d x length(oi) matrix that lifts a residual defined on the
### observed rows to the full vector: identity on the observed rows, P on the
### missing ones. With nothing observed, S has no columns and Omega is the whole
### of Sigma, which is the correct degenerate case.
###
### Arguments
###   Sigma  d x d measurement covariance
###   oi     integer vector of observed rows
###   d      dimension of the observation vector
###   regularization  small ridge added before inversion, as elsewhere in the
###                   package
`stem_missing_blocks` <- function(Sigma, oi, d, regularization = 0) {
  mi <- setdiff(seq_len(d), oi)

  if (length(mi) == 0L) {
    S <- diag(1, d)
    return(list(oi = oi, mi = mi, P = NULL,
                Omega = matrix(0, d, d), S = S))
  }

  Omega_full <- matrix(0, d, d)

  if (length(oi) == 0L) {
    ### nothing observed: the measurement error keeps its whole variance
    Omega_full[mi, mi] <- Sigma
    return(list(oi = oi, mi = mi, P = NULL,
                Omega = Omega_full, S = matrix(0, d, 0)))
  }

  Soo <- Sigma[oi, oi, drop = FALSE]
  Smo <- Sigma[mi, oi, drop = FALSE]
  Som <- Sigma[oi, mi, drop = FALSE]
  Smm <- Sigma[mi, mi, drop = FALSE]

  Soo_inv <- try(solve(diag(regularization, nrow(Soo)) + Soo), silent = TRUE)
  if (inherits(Soo_inv, "try-error")) {
    Soo_inv <- MASS::ginv(Soo)
  }

  P <- Smo %*% Soo_inv
  Omega_full[mi, mi] <- Smm - P %*% Som

  S <- matrix(0, d, length(oi))
  S[cbind(oi, seq_along(oi))] <- 1
  S[mi, ] <- P

  list(oi = oi, mi = mi, P = P, Omega = Omega_full, S = S)
}


### Cache of the per-pattern blocks.
###
### Sigma_e is recomputed at every EM iteration, so the cache lives for one
### iteration only and is keyed by the missingness pattern.
`stem_blocks_cache` <- function(Sigma, d, regularization = 0) {
  store <- new.env(parent = emptyenv())
  function(key, oi) {
    ### Callers branch on the complete case before reaching here, so a missing
    ### key would be a programming error rather than a data condition.
    if (is.na(key)) {
      stop("stem_blocks_cache() called on a complete time point")
    }
    hit <- store[[key]]
    if (is.null(hit)) {
      hit <- stem_missing_blocks(Sigma, oi, d, regularization)
      hit$complete <- FALSE
      assign(key, hit, envir = store)
    }
    hit
  }
}
