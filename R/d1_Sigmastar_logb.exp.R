#' @keywords internal
#' @noRd

### Derivative of the covariance with respect to log(b). It is exp(logb) times
### the identity, and it is returned as that single number rather than as a
### dense d x d matrix: d1_Q(), d2_Q() and d12_Q() recognise a scalar as a
### multiple of the identity and replace the matrix products by scalings, which
### is O(d^2) instead of O(d^3). At d = 400 the dense form allocated 160,000
### doubles to carry 400 of them, and then multiplied by them.
###
### Superseded:
###   d1 = diag(exp(logb), d)

`d1_Sigmastar_logb.exp` <-
function(logb, d) {
		exp(logb)
}
