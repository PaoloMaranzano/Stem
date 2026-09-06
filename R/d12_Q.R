#' @keywords internal
#' @noRd

### Mixed second derivative, with the same treatment. `Pt` and `Pb` are
### X^{-1} d1_X_theta and X^{-1} d1_X_logb, shared with d1_Q() and d2_Q(); the
### second is a scaling rather than a matrix product, since the derivative with
### respect to log(b) is a multiple of the identity -- see d1_Q().
###
### Superseded:
###   d12_Q <- -n*sum(diag(solve(X) %*% d1_X_theta %*% solve(X) %*% d1_X_logb)) +
###            (2/sigma2omega)*sum(diag(solve(X) %*% d1_X_theta %*% solve(X) %*%
###                                     d1_X_logb %*% solve(X) %*% B))

`d12_Q` <-
function(n, X, d1_X_theta, d1_X_logb, sigma2omega, B, Xi = NULL, XiB = NULL,
         Pt = NULL, Pb = NULL){
	if (is.null(Xi))  Xi  <- solve(X)
	if (is.null(XiB)) XiB <- Xi %*% B
	if (is.null(Pt))  Pt  <- stem_xprod(Xi, d1_X_theta)
	if (is.null(Pb))  Pb  <- stem_xprod(Xi, d1_X_logb)
	-n * sum(Pt * t(Pb)) + (2/sigma2omega) * sum((Pt %*% Pb) * t(XiB))
}
