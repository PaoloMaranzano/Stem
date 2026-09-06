#' @keywords internal
#' @noRd

### First derivative of the profiled objective with respect to one of the two
### covariance parameters. The inverse of X is taken once instead of twice, the
### product X^{-1}B is supplied by the caller when it is shared between the
### derivatives, and both traces use tr(AB) = sum(A * t(B)), which is O(d^2).
###
### `d1_X` is either a d x d matrix or a single number standing for that
### multiple of the identity, which is what the derivative with respect to
### log(b) is; in the second case the product X^{-1} d1_X is a scaling rather
### than a matrix product. `P` is that product, which the caller supplies when
### the same derivative enters more than one of d1_Q(), d2_Q() and d12_Q().
###
### Superseded:
###   d1_Q <- n*sum(diag(solve(X) %*% d1_X)) -
###           (1/sigma2omega)*sum(diag(solve(X) %*% d1_X %*% solve(X) %*% B))

`d1_Q` <-
function(n, X, d1_X, sigma2omega, B, Xi = NULL, XiB = NULL, P = NULL){
	if (is.null(Xi))  Xi  <- solve(X)
	if (is.null(XiB)) XiB <- Xi %*% B
	if (is.null(P))   P   <- stem_xprod(Xi, d1_X)
	n * stem_xtrace(Xi, d1_X) - (1/sigma2omega) * sum(P * t(XiB))
}
