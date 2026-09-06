#' @keywords internal
#' @noRd

### Second derivative of the covariance with respect to log(b). It coincides
### with the first, and like it is returned as the multiplier of the identity
### rather than as a dense matrix. See d1_Sigmastar_logb.exp().
###
### Superseded:
###   d2 = diag(exp(logb), d)

`d2_Sigmastar_logb.exp` <-
function(logb, d) {
		exp(logb)
}
