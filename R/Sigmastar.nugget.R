#' @keywords internal
#' @noRd

### The measurement correlation of a model with NO spatial structure: the
### identity, up to the same exp(logb) the exponential form adds to its
### diagonal, so that Sigma_e = sigma2omega * (1 + exp(logb)) I is a single
### variance. It is what STEM_Estimation() substitutes for the exponential
### covariance when `spatial = FALSE`.
###
### Only the SUM sigma2eps + sigma2omega is then identified: their split is
### whatever the starting value of logb made it, and the range has no meaning
### at all, which is why the Newton-Raphson step is skipped in that case rather
### than left to wander. The signature matches Sigmastar.exp() so that the two
### are interchangeable wherever cov.spat is passed.

`Sigmastar.nugget` <-
function(d, logb, logtheta, dist, E = NULL) {
	diag(1 + exp(logb), d)
}
