#' @keywords internal
#' @noRd

`smoothing` <-
  function(ss) {
    m <- ss$m
    C <- ss$C
    m0<- ss$m0
    C0<- ss$C0
    d <- ss$d
    nobs <- ss$n

    for (tt in (nobs-1):1)    {

      if (ss$p == 1)
        nextstep <-
          smootherstep.uni(
            m[tt,],
            C[[tt]],
            ss$Gmat,
            ss$Wmat,
            m[tt+1,],
            C[[tt+1]])

      else
        nextstep <-
          smootherstep(
            matrix(m[tt,],nrow=1),
            C[[tt]],
            ss$Gmat,
            ss$Wmat,
            matrix(m[tt+1,],nrow=1),
            C[[tt+1]])

      m[tt,]  <- nextstep$ms
      C[[tt]] <- nextstep$Cs
    }

    ### The smoother used to build, row by row inside the loop, the n x d matrix
    ### mu with mu[t,] = F' m_t, and return it as ss$mu. Nothing in the package
    ### ever read it: kalman() forms its own fitted values from the smoothed
    ### states, and the quantity is not part of any returned object. It is no
    ### longer computed. At n = 730 and d = 400 it allocated and threw away 2.3
    ### megabytes on every EM iteration, hundreds of times over an SC-STEM grid.
    ###
    ### Superseded:
    ###   mu <- matrix(NA, nobs, d)
    ###   mu[nobs,] <- t(ss$Fmat) %*% m[nobs,]        # before the loop
    ###   mu[tt,]   <- t(ss$Fmat) %*% m[tt,]          # inside the loop
    ###   ss$mu     <- mu

    ss$m0 <- m0
    ss$C0 <- C0
    if (stats::is.ts(ss$Z)) {
      ss$m <- stats::ts(m,start = stats::start(ss$z),
                 end=stats::end(ss$z),
                 frequency=stats::frequency(ss$z))
    } else {
      ss$m <- m
    }
    ss$C <- C
    ss
  }

