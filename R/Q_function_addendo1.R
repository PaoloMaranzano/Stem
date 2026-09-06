#' @keywords internal
#' @noRd

### tr(A B) is the sum of the elementwise product of A with the transpose of B,
### so it costs O(d^2) and does not need the O(d^3) product to be formed. The
### inverse is taken once and can be supplied by the caller, which reuses it
### across the several derivative evaluations of one Newton-Raphson step.
###
### Superseded:
###   Q_addendo1 <- n*log(det(sigma2omega * Sigmastar)) +
###                 sum(diag(solve(sigma2omega * Sigmastar) %*% B))

`Q_function_addendo1` <-
function(sigma2omega, n, Sigmastar, B, Si = NULL){
        S  <- sigma2omega * Sigmastar
        if (is.null(Si)) Si <- solve(S)
        n * determinant(S, logarithm = TRUE)$modulus[1] + sum(Si * t(B))
}
