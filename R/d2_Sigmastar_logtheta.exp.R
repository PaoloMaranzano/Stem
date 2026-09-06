#' @keywords internal
#' @noRd

### Second derivative of the exponential correlation with respect to
### log(theta). See d1_Sigmastar_logtheta.exp() for `E`: the kernel and the
### product theta * h are shared with the first derivative and with
### Sigmastar.exp(), which used to evaluate them four times over between them.
###
### Superseded:
###   d2 = exp(-exp(logtheta) * dist) * (-exp(logtheta) * dist) *
###        (1 - exp(logtheta) * dist)
###   diag(d2) = 0

`d2_Sigmastar_logtheta.exp` <-
function(logtheta, dist, E = NULL) {
		td <- exp(logtheta) * dist
		if (is.null(E)) E <- exp(-td)
		d2 <- E * (-td) * (1 - td)
		diag(d2) <- 0
		d2
}
