#' @keywords internal
#' @noRd

### ---------------------------------------------------------------------------
### The forward pass
### ---------------------------------------------------------------------------
### The loop itself is unchanged in structure. What changed is that the
### constants of the measurement equation are now hoisted out of it: the
### Cholesky factor of Sigma_e, the matrix U = Sigma_e^{-1} Z, the p x p product
### Z' Sigma_e^{-1} Z and log|Sigma_e| do not depend on the time point, so they
### are computed once here and handed to filterstep(), which then completes each
### step by the Woodbury identity in O(d^2) instead of O(d^3). See the header of
### R/filterstep.R for the algebra.
###
### With missing values the constants depend on which rows are observed, so they
### are computed once per DISTINCT missingness pattern and cached. In a
### monitoring network a station is typically out of service for a stretch of
### consecutive days, so the number of distinct patterns is far smaller than the
### number of time points and the cache almost always hits.
###
### The superseded version simply called filterstep() with the full Vmat at
### every step:
###
###   `filtering` <- function(ss) {
###     m <- matrix(NA,ss$n,ss$p); C <- vector("list",ss$n)
###     firststep <- filterstep(z = matrix(ss$z[1,]), Fmat = ss$Fmat,
###                             Gmat = ss$Gmat, Vt = ss$Vmat, Wt = ss$Wmat,
###                             mx = ss$m0, Cx = ss$C0, XXXcov = ss$XXX[,,1],
###                             betacov = ss$beta, flag = ss$flag.cov)
###     m[1,] <- firststep$m; C[[1]] <- firststep$C
###     loglik <- firststep$loglikterm
###     for (tt in 2:ss$n) {
###       nextstep <- filterstep(z = matrix(ss$z[tt,]), Fmat = ss$Fmat,
###                              Gmat = ss$Gmat, Vt = ss$Vmat, Wt = ss$Wmat,
###                              mx = matrix(m[tt-1,],nrow=1), Cx = C[[tt-1]],
###                              XXXcov = ss$XXX[,,tt], betacov = ss$beta,
###                              flag = ss$flag.cov)
###       m[tt,] <- nextstep$m; C[[tt]] <- nextstep$C
###       loglik <- loglik + nextstep$loglikterm
###     }
###     ss$m <- m; ss$C <- C; ss$loglik <- loglik
###     ss
###   }

`filtering` <-
function(ss) {

	Zfull <- t(ss$Fmat)                       # d x p loading matrix
	V     <- ss$Vmat                          # d x d, constant in t
	flag  <- isTRUE(ss$flag.cov)
	zmat  <- as.matrix(ss$z)

	obs <- stem_obs_index(zmat)

	### constants of the measurement equation, cached by missingness pattern
	store <- new.env(parent = emptyenv())
	pattern_const <- function(key, oi) {
		if (length(oi) == 0L) return(NULL)
		kk <- if (is.na(key)) "*complete*" else key
		hit <- store[[kk]]
		if (!is.null(hit)) return(hit)
		Zo <- Zfull[oi, , drop = FALSE]
		Vo <- V[oi, oi, drop = FALSE]
		Lc <- chol(Vo)                                   # Vo = t(Lc) %*% Lc
		VinvZ <- backsolve(Lc, forwardsolve(t(Lc), Zo))  # Vo^{-1} Zo
		hit <- list(Lc = Lc, VinvZ = VinvZ,
		            FU = crossprod(Zo, VinvZ),
		            logdetV = 2 * sum(log(diag(Lc))),
		            Zo = Zo, nobs = length(oi), oi = oi)
		assign(kk, hit, envir = store)
		hit
	}

	m <- matrix(NA_real_, ss$n, ss$p)
	C <- vector("list", ss$n)
	loglik <- 0
	mprev <- ss$m0                      # 1 x p
	Cprev <- ss$C0

	for (tt in seq_len(ss$n)) {
		oi  <- obs$idx[[tt]]
		cst <- pattern_const(obs$key[tt], oi)
		Xo  <- if (flag && length(oi)) matrix(ss$XXX[oi, , tt], nrow = length(oi)) else NULL

		step <- filterstep(cst    = cst,
		                   zobs   = if (length(oi)) matrix(zmat[tt, oi], ncol = 1) else NULL,
		                   Xobs   = Xo,
		                   betacov = ss$beta,
		                   Gmat   = ss$Gmat,
		                   Wt     = ss$Wmat,
		                   mx     = mprev,
		                   Cx     = Cprev,
		                   flag   = flag)

		m[tt, ]  <- step$m
		C[[tt]]  <- step$C
		loglik   <- loglik + step$loglikterm
		mprev    <- matrix(step$m, nrow = 1)
		Cprev    <- step$C
	}

	if (stats::is.ts(ss$z)) {
		ss$m <- stats::ts(m, stats::start(ss$z), end = stats::end(ss$z),
		                  frequency = stats::frequency(ss$z))
	} else {
		ss$m <- m
	}
	ss$C <- C
	ss$loglik <- loglik

	ss
}
