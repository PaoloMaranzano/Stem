#' @keywords internal
#' @noRd

### ---------------------------------------------------------------------------
### One step of the Kalman filter, in Woodbury form
### ---------------------------------------------------------------------------
### The step used to form and invert the d x d matrix
###
###     Q_t = Z R_t Z' + V ,      Z = t(Fmat),   V = Sigma_e ,
###
### at every time point. That costs O(d^3) per step and O(T d^3) over the whole
### pass, and it dominates everything else as soon as the network has more than
### a few dozen locations.
###
### Two facts make it avoidable. First, V does not depend on t: it is the
### measurement covariance, rebuilt once per EM iteration. Second, Z R_t Z' has
### rank p, the dimension of the latent state, which is typically 1. The
### Woodbury identity therefore gives
###
###     Q_t^{-1} = V^{-1} - U (R_t^{-1} + Z'V^{-1}Z)^{-1} U' ,     U = V^{-1} Z,
###
### and the matrix determinant lemma gives
###
###     log|Q_t| = log|V| + log|R_t| + log|R_t^{-1} + Z'V^{-1}Z| ,
###
### where U, Z'V^{-1}Z, the Cholesky factor of V and log|V| are all constants of
### the pass: they are computed once, in filtering(), and handed to this
### function in `cst`. What is left at each step is p x p algebra plus a single
### triangular solve against the cached factor, so the cost per step falls from
### O(d^3) to O(d^2). The recursion never needs Q_t itself: the gain enters only
### through Z'Q^{-1}Z and Z'Q^{-1}e, both p x p or p x 1.
###
### The log-density is also evaluated directly rather than through
### mvtnorm::dmvnorm. That is not only faster -- dmvnorm re-checks the symmetry
### of Q with all.equal at every step, which on a 36-station network was about a
### seventh of the total running time -- but necessary: the Gaussian density on
### d observations is of order exp(-d), so on a few hundred locations it falls
### below the smallest representable double and the previous
### log(dmvnorm(...)) returned -Inf, breaking the EM algorithm outright.
###
### Missing values are handled as before, by restricting the measurement
### equation to the rows actually observed (Durbin and Koopman 2012,
### Sect. 4.10); the constants above then depend on the missingness pattern, and
### filtering() caches them per distinct pattern.
###
### The superseded implementation, kept for reference:
###
###   `filterstep` <- function(z,Fmat,Gmat,Vt,Wt,mx,Cx,XXXcov,betacov,flag) {
###       a <- Gmat %*% t(mx)
###       R <- Gmat %*% Cx %*% t(Gmat) + Wt
###       z  <- as.matrix(z)
###       oi <- which(!is.na(z[, 1]))
###       if (length(oi) == 0L) return(list(m = a, C = R, loglikterm = 0, n_obs = 0L))
###       if (length(oi) < nrow(z)) {
###         z    <- z[oi, , drop = FALSE]
###         Fmat <- Fmat[, oi, drop = FALSE]
###         Vt   <- Vt[oi, oi, drop = FALSE]
###         if (!is.null(XXXcov)) XXXcov <- as.matrix(XXXcov)[oi, , drop = FALSE]
###       }
###       f <- t(Fmat) %*% a
###       Q <- t(Fmat) %*% R %*% Fmat + Vt          # <- the d x d matrix
###       if(flag==FALSE) e <- z - f else e <- z - (XXXcov %*% betacov) - f
###       A <- R %*% Fmat %*% solve(Q)              # <- the O(d^3) inversion
###       m <- a + A%*%e
###       C <- R - A%*%Q%*%t(A)
###       if (length(z)>1 && flag==TRUE)
###         loglikterm <- mvtnorm::dmvnorm(as.numeric(z),
###                         as.numeric(f+(XXXcov %*% betacov)), Q, log=TRUE)
###       ... (the three remaining branches, identical in structure)
###       list(m=m,C=C,loglikterm=loglikterm,n_obs=length(oi))
###   }
###
### Arguments
###   cst      constants of the current missingness pattern, from filtering():
###            Lc (Cholesky of the observed block of V), VinvZ, FU = Z'V^{-1}Z,
###            logdetV, Zo (the observed rows of Z) and nobs
###   zobs     the observed entries of z at this time point, a column vector
###   Xobs     the corresponding rows of the design matrix, or NULL
###   betacov  the regression coefficients
###   Gmat, Wt the transition matrix and the innovation covariance
###   mx, Cx   the filtered mean (1 x p) and variance at the previous step
###   flag     whether the regression term is present

`filterstep` <-
function(cst, zobs, Xobs, betacov, Gmat, Wt, mx, Cx, flag) {

    a  <- Gmat %*% t(mx)                              # p x 1
    R  <- Gmat %*% Cx %*% t(Gmat) + Wt                # p x p

    ### nothing observed: no update and no likelihood contribution, which is the
    ### Z_t = 0 device of Durbin and Koopman (2012), Sect. 4.10
    if (is.null(cst)) {
      return(list(m = a, C = R, loglikterm = 0, n_obs = 0L))
    }

    Ri <- solve(R)
    Mm <- solve(Ri + cst$FU)                          # p x p
    ZQiZ <- cst$FU - cst$FU %*% Mm %*% cst$FU         # Z' Q^{-1} Z

    mu <- cst$Zo %*% a
    if (isTRUE(flag)) mu <- mu + Xobs %*% betacov
    e <- zobs - mu

    Ue   <- crossprod(cst$VinvZ, e)                   # Z' V^{-1} e
    ZQie <- Ue - cst$FU %*% Mm %*% Ue                 # Z' Q^{-1} e
    m    <- a + R %*% ZQie
    C    <- R - R %*% ZQiZ %*% R

    ### quadratic form and log determinant, without ever forming Q
    Le   <- forwardsolve(t(cst$Lc), e)
    quad <- sum(Le^2) - as.numeric(crossprod(Ue, Mm %*% Ue))
    ldQ  <- cst$logdetV + log(det(R)) + log(det(Ri + cst$FU))
    loglikterm <- -0.5 * (cst$nobs * log(2 * pi) + ldQ + quad)

    list(m = m, C = C, loglikterm = loglikterm, n_obs = cst$nobs)
  }
