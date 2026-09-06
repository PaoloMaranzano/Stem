#' @keywords internal
#' @noRd

### Second derivative. The superseded expression called solve(X) SEVEN times and
### formed a dozen d x d products only to read their diagonals; at a few hundred
### locations that single line dominated the whole estimation. Here the inverse
### is taken once, the two products X^{-1} d1_X and X^{-1} B are formed once,
### and every trace is evaluated as sum(A * t(B)).
###
### Superseded:
###   d2_Q <- n*sum(diag(solve(X) %*% d2_X)) -
###           n*sum(diag(solve(X) %*% d1_X %*% solve(X) %*% d1_X)) -
###           (1/sigma2omega)*sum(diag(solve(X) %*% d2_X %*% solve(X) %*% B)) +
###           (2/sigma2omega)*sum(diag(solve(X) %*% d1_X %*% solve(X) %*% d1_X %*%
###                                    solve(X) %*% B))

`d2_Q` <-
function(n, X, d1_X, d2_X, sigma2omega, B, Xi = NULL, XiB = NULL){
	if (is.null(Xi))  Xi  <- solve(X)
	if (is.null(XiB)) XiB <- Xi %*% B
	P  <- Xi %*% d1_X            # X^{-1} d1_X
	P2 <- Xi %*% d2_X            # X^{-1} d2_X
	n * sum(Xi * t(d2_X)) -
	  n * sum(P * t(P)) -
	  (1/sigma2omega) * sum(P2 * t(XiB)) +
	  (2/sigma2omega) * sum((P %*% P) * t(XiB))
}
