#' @keywords internal
#' @noRd

### First derivative of the exponential correlation with respect to log(theta).
###
### `E` is the exponential kernel exp(-theta h), shared with Sigmastar.exp() and
### with the second derivative; supplying it avoids recomputing a d x d
### elementwise exp() three times per Newton-Raphson step. The product
### theta * h is likewise formed once instead of twice.
###
### Superseded:
###   d1 = exp(-exp(logtheta) * dist) * (-exp(logtheta) * dist)
###   diag(d1) = 0

`d1_Sigmastar_logtheta.exp` <-
function(logtheta, dist, E = NULL) {
		td <- exp(logtheta) * dist
		if (is.null(E)) E <- exp(-td)
		d1 <- E * (-td)
		diag(d1) <- 0
		d1
}
