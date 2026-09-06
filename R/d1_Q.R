#' @keywords internal
#' @noRd

### First derivative of the profiled objective with respect to one of the two
### covariance parameters. The inverse of X is taken once instead of twice, the
### product X^{-1}B is supplied by the caller when it is shared between the
### derivatives, and both traces use tr(AB) = sum(A * t(B)), which is O(d^2).
###
### Superseded:
###   d1_Q <- n*sum(diag(solve(X) %*% d1_X)) -
###           (1/sigma2omega)*sum(diag(solve(X) %*% d1_X %*% solve(X) %*% B))

`d1_Q` <-
function(n, X, d1_X, sigma2omega, B, Xi = NULL, XiB = NULL){
	if (is.null(Xi))  Xi  <- solve(X)
	if (is.null(XiB)) XiB <- Xi %*% B
	P <- Xi %*% d1_X
	n * sum(Xi * t(d1_X)) - (1/sigma2omega) * sum(P * t(XiB))
}
