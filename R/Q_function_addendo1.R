#' @keywords internal
#' @noRd

### tr(A B) is the sum of the elementwise product of A with the transpose of B,
### so it costs O(d^2) and does not need the O(d^3) product to be formed.
###
### The log determinant and the inverse come from ONE Cholesky factorisation
### rather than from two separate ones. The superseded version called solve(),
### which factorises S, and then determinant(), which factorises it again: two
### O(d^3) decompositions of the same positive definite matrix where one does.
### The determinant is read off the diagonal of the factor, and chol2inv()
### inverts from it. The fallback covers the case of a matrix that has lost
### positive definiteness numerically, which the grid search can produce at
### extreme values of the covariance parameters.
###
### `Si` is the inverse when the caller already holds it.
###
### Superseded:
###   Q_addendo1 <- n*log(det(sigma2omega * Sigmastar)) +
###                 sum(diag(solve(sigma2omega * Sigmastar) %*% B))
### and then
###   S  <- sigma2omega * Sigmastar
###   if (is.null(Si)) Si <- solve(S)
###   n * determinant(S, logarithm = TRUE)$modulus[1] + sum(Si * t(B))

`Q_function_addendo1` <-
function(sigma2omega, n, Sigmastar, B, Si = NULL){
        S <- sigma2omega * Sigmastar
        if (!is.null(Si)) {
                return(n * determinant(S, logarithm = TRUE)$modulus[1] +
                       sum(Si * t(B)))
        }
        R <- tryCatch(chol(S), error = function(e) NULL)
        if (is.null(R)) {
                return(n * determinant(S, logarithm = TRUE)$modulus[1] +
                       sum(solve(S) * t(B)))
        }
        2 * n * sum(log(diag(R))) + sum(chol2inv(R) * t(B))
}
