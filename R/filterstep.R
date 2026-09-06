#' @keywords internal
#' @noRd

### One step of the Kalman filter.
###
### Missing values in z are handled as prescribed by Durbin and Koopman (2012,
### 2nd ed.), Section 4.10, p. 111: the observation equation is restricted to
### the rows actually observed, that is z, the loading matrix and the
### measurement covariance are pre-multiplied by the selection matrix W_t whose
### rows are a subset of the rows of the identity. The recursion is then the
### standard one on an observation vector of reduced dimension. When nothing is
### observed the update is skipped entirely -- Durbin and Koopman obtain this by
### setting Z_t = 0 -- so that the filtered moments equal the predicted ones and
### the time point contributes nothing to the likelihood.
###
### Note that `Fmat` is p x d here, so t(Fmat) is the d x p loading matrix Z_t
### of Durbin and Koopman, and selecting observed rows of Z_t means selecting
### observed COLUMNS of Fmat.

`filterstep` <-
function(z,Fmat,Gmat,Vt,Wt,mx,Cx,XXXcov,betacov,flag)
  {
    a <- Gmat %*% t(mx)
    R <- Gmat %*% Cx %*% t(Gmat) + Wt

    ### selection matrix of Durbin and Koopman (2012), Sect. 4.10
    z  <- as.matrix(z)
    oi <- which(!is.na(z[, 1]))

    if (length(oi) == 0L) {
      ### every element missing: no update, no likelihood contribution
      return(list(m = a, C = R, loglikterm = 0, n_obs = 0L))
    }

    if (length(oi) < nrow(z)) {
      z       <- z[oi, , drop = FALSE]
      Fmat    <- Fmat[, oi, drop = FALSE]
      Vt      <- Vt[oi, oi, drop = FALSE]
      if (!is.null(XXXcov)) XXXcov <- as.matrix(XXXcov)[oi, , drop = FALSE]
    }

    f <- t(Fmat) %*% a
    Q <- t(Fmat) %*% R %*% Fmat + Vt

    if(flag==FALSE)        e <- z - f
    else                 e <- z - (XXXcov %*% betacov) - f

    A <- R %*% Fmat %*% solve(Q)
    m <- a + A%*%e
    C <- R - A%*%Q%*%t(A)

    ### The density has to be evaluated on the log scale directly. Taking the
    ### logarithm of the density itself underflows: with d observations the
    ### Gaussian density is of order exp(-d), and for a network of a few hundred
    ### locations it falls below the smallest representable double, so that
    ### log(dmvnorm(...)) returns -Inf and the EM algorithm breaks. The switch
    ### costs nothing and makes large networks estimable.
    if (length(z)>1  &&  flag==TRUE)  loglikterm <-mvtnorm::dmvnorm(as.numeric(z),as.numeric(f+(XXXcov %*% betacov)),Q,log=TRUE)
    if (length(z)>1  &&  flag==FALSE)  loglikterm <-mvtnorm::dmvnorm(as.numeric(z),as.numeric(f),Q,log=TRUE)
    ### The mean of the one-step-ahead predictive distribution is
    ### mu_t = X_t beta + K y_t^{t-1} (Fasso and Cameletti 2010, Eq. 9), so the
    ### regression term ADDS to f. The univariate branch subtracted it, which
    ### made the log-likelihood wrong whenever a fit had a single location --
    ### rare for a pooled model, reachable for a cluster-wise one, and now also
    ### reachable at a time point where a single location is observed.
    if (length(z)==1 &&  flag==TRUE)  loglikterm <- -0.5*( log(2*pi) + log(Q) + (z - (f+XXXcov %*% betacov)) ^2/Q)
    if (length(z)==1 &&  flag==FALSE)  loglikterm <- -0.5*( log(2*pi) + log(Q) + (z-f)^2/Q)

  list(m=m,C=C,loglikterm=loglikterm,n_obs=length(oi))
  }
