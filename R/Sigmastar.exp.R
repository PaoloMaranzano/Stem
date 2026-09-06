#' @keywords internal
#' @noRd

### The exponential spatial correlation, exp(-theta h), with exp(logb) added to
### the diagonal.
###
### `E` is the exponential kernel itself. Within one Newton-Raphson step the
### same kernel is needed here and by both derivatives with respect to
### log(theta), so the caller computes it once and passes it in; at a few
### hundred locations each evaluation is a d x d elementwise exp(), which is not
### free. When it is not supplied the function computes it, so the signature
### remains the one the rest of the package uses.
###
### The diagonal is updated in place rather than adding diag(exp(logb), d),
### which would allocate a second dense d x d matrix to hold d numbers.
###
### Superseded:
###   Sigmastar.exp <- function(d, logb, logtheta, dist)
###     diag(exp(logb), d) + exp(-exp(logtheta) * dist)

`Sigmastar.exp` <-
function(d, logb, logtheta, dist, E = NULL) {
	if (is.null(E)) E <- exp(-exp(logtheta) * dist)
	diag(E) <- diag(E) + exp(logb)
	E
}
