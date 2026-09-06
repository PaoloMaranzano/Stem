#' @keywords internal
#' @noRd

### Mixed second derivative, with the same treatment.
###
### Superseded:
###   d12_Q <- -n*sum(diag(solve(X) %*% d1_X_theta %*% solve(X) %*% d1_X_logb)) +
###            (2/sigma2omega)*sum(diag(solve(X) %*% d1_X_theta %*% solve(X) %*%
###                                     d1_X_logb %*% solve(X) %*% B))

`d12_Q` <-
function(n, X, d1_X_theta, d1_X_logb, sigma2omega, B, Xi = NULL, XiB = NULL){
	if (is.null(Xi))  Xi  <- solve(X)
	if (is.null(XiB)) XiB <- Xi %*% B
	Pt <- Xi %*% d1_X_theta
	Pb <- Xi %*% d1_X_logb
	-n * sum(Pt * t(Pb)) + (2/sigma2omega) * sum((Pt %*% Pb) * t(XiB))
}
