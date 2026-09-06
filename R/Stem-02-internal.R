#' @keywords internal
#' @noRd

### ---------------------------------------------------------------------------
### Structured products against the inverse covariance
### ---------------------------------------------------------------------------
### The derivatives of the objective with respect to the two covariance
### parameters are traces of products of X^{-1} with dSigma/dparameter. For
### log(theta) that derivative is a full d x d matrix and the product is a
### genuine matrix product. For log(b) it is exp(logb) times the identity, so
### the product is a SCALING and the trace is a multiple of tr(X^{-1}): forming
### the dense diagonal matrix and multiplying by it costs O(d^3) to compute
### something that costs O(d^2).
###
### These two helpers let d1_Q(), d2_Q() and d12_Q() take the derivative either
### as a d x d matrix or as a single number standing for that multiple of the
### identity, and pick the cheap path in the second case.

`stem_xprod` <- function(Xi, A) if (length(A) == 1L) A * Xi else Xi %*% A

`stem_xtrace` <- function(Xi, A) {
  if (length(A) == 1L) A * sum(diag(Xi)) else sum(Xi * t(A))
}

NULL
