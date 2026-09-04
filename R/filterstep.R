#' @keywords internal
#' @noRd

`filterstep` <-
function(z,Fmat,Gmat,Vt,Wt,mx,Cx,XXXcov,betacov,flag)
  {
    a <- Gmat %*% t(mx)
    R <- Gmat %*% Cx %*% t(Gmat) + Wt
    f <- t(Fmat) %*% a
    Q <- t(Fmat) %*% R %*% Fmat + Vt

    if(flag==FALSE)        e <- z - f
    else                 e <- z - (XXXcov %*% betacov) - f

    A <- R %*% Fmat %*% solve(Q)
    m <- a + A%*%e
    C <- R - A%*%Q%*%t(A)

    if (length(z)>1  &&  flag==TRUE)  loglikterm <-log(mvtnorm::dmvnorm(as.numeric(z),as.numeric(f+(XXXcov %*% betacov)),Q))
    if (length(z)>1  &&  flag==FALSE)  loglikterm <-log(mvtnorm::dmvnorm(as.numeric(z),as.numeric(f),Q))
    ### The mean of the one-step-ahead predictive distribution is
    ### mu_t = X_t beta + K y_t^{t-1} (Fasso and Cameletti 2010, Eq. 9), so the
    ### regression term ADDS to f. The univariate branch subtracted it, which
    ### made the log-likelihood wrong whenever a fit had a single location --
    ### rare for a pooled model, reachable for a cluster-wise one.
    if (length(z)==1 &&  flag==TRUE)  loglikterm <- -0.5*( log(2*pi) + log(Q) + (z - (f+XXXcov %*% betacov)) ^2/Q)
    if (length(z)==1 &&  flag==FALSE)  loglikterm <- -0.5*( log(2*pi) + log(Q) + (z-f)^2/Q)

  list(m=m,C=C,loglikterm=loglikterm)
  }

