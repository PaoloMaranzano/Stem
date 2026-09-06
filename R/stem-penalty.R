#' @keywords internal
#' @noRd

### ---------------------------------------------------------------------------
### Regularized update of the regression coefficients
### ---------------------------------------------------------------------------
### The M-step for beta maximizes, over beta alone,
###
###     -1/2 beta' M beta + beta' v ,
###
###     M = sum_t X_t' Sigma_e^{-1} X_t ,   v = sum_t X_t' Sigma_e^{-1} v_t ,
###
### with v_t the completed residual of period t net of the latent contribution.
### Its unpenalized maximizer is M^{-1} v, which is what the package computes
### when no penalty is asked for.
###
### An elastic net on the coefficients subtracts, in the glmnet parameterization,
###
###     lambda { alpha ||D beta||_1 + (1 - alpha)/2 * beta' D beta } ,
###
### with D the diagonal indicator of the penalized coordinates. alpha = 0 is
### ridge, alpha = 1 is the lasso, and anything between is the elastic net.
###
### WHY THE PENALTY DOES NOT DISTURB THE EM. The E-step is untouched: the
### penalty is a function of beta alone and does not involve the latent states,
### so the usual Jensen argument applies verbatim to Q(psi) - pen(beta) and the
### PENALIZED observed-data log-likelihood increases at every iteration. What
### changes is only which point of the M-step is returned.
###
### RIDGE has a closed form. Differentiating,
###
###     v - M beta - lambda (1 - alpha) D beta = 0
###       =>  beta = (M + lambda (1 - alpha) D)^{-1} v ,
###
### one extra term in a matrix that is already formed, at no cost worth naming.
###
### THE LASSO AND THE ELASTIC NET have no closed form, but the objective is a
### quadratic plus a SEPARABLE penalty, which is exactly the setting in which
### coordinate descent applies: cycling over the coordinates, each update is a
### soft-thresholding,
###
###     beta_j <- S( v_j - sum_{k != j} M_jk beta_k , lambda alpha w_j )
###               / ( M_jj + lambda (1 - alpha) w_j ) ,
###
###     S(x, t) = sign(x) (|x| - t)_+ ,
###
### with w_j = 1 for a penalized coordinate and 0 otherwise -- an unpenalized
### coordinate, the intercept in particular, therefore keeps its ordinary least
### squares update given the others. The objective is convex and the penalty
### separable, so the cycle converges to the exact maximizer of the M-step: this
### is a genuine EM, not merely a generalized one.
###
### THE METRIC MATTERS. M is a GLS cross-product, not a Euclidean one: it
### carries Sigma_e^{-1}. Standardizing the columns of X in the usual Euclidean
### sense is therefore NOT the right normalization here. The scaling that makes
### lambda comparable across covariates, and across regimes in SC-STEM, is the
### one that puts the diagonal of M at one, and that is what `standardize`
### does -- internally, so that the coefficients returned are always on the
### original scale of the covariates.
###
### THE INTERCEPT IS NEVER PENALIZED BY DEFAULT, and not merely by convention:
### the model already carries a latent process whose initial mean m0 absorbs the
### level, so the intercept and m0 are only weakly identified apart. Shrinking
### the intercept would move that identification about arbitrarily.

### The soft-thresholding operator
`stem_soft` <- function(x, t) sign(x) * pmax(abs(x) - t, 0)

### Which coordinates are penalized. The default excludes the intercept, found
### as the column of the design that is constant.
`stem_penalized_index` <- function(r, penalize = NULL, intercept = 1L) {
	if (is.null(penalize)) {
		w <- rep(1, r)
		if (!is.na(intercept) && intercept >= 1L && intercept <= r) w[intercept] <- 0
		return(w)
	}
	if (is.logical(penalize)) return(as.numeric(penalize))
	w <- rep(0, r); w[penalize] <- 1; w
}

### The update itself.
###
###   M, v         the two accumulations of the M-step
###   alpha        0 ridge, 1 lasso, in between elastic net
###   lambda       the overall penalty strength; 0 returns the ordinary update
###   w            the 0/1 weights of stem_penalized_index()
###   beta0        the current coefficients, used to start the coordinate descent
###   standardize  put the diagonal of M at one before penalizing
###   ridge_reg    the small ridge the package already adds for conditioning
`stem_beta_update` <- function(M, v, alpha = 0, lambda = 0, w = NULL,
                               beta0 = NULL, standardize = TRUE,
                               ridge_reg = 0, tol = 1e-9, maxit = 1000L) {

	r <- nrow(M)
	if (is.null(w)) w <- stem_penalized_index(r)
	v <- as.numeric(v)

	### no penalty: the estimator the package has always computed
	if (lambda <= 0) {
		return(list(beta = solve(diag(ridge_reg, r) + M, v), df = r, iter = 0L))
	}

	### the scaling that makes lambda comparable across covariates. s_j is the
	### GLS norm of column j; beta is solved for on the scaled design and
	### returned on the original one.
	s <- if (isTRUE(standardize)) sqrt(pmax(diag(M), .Machine$double.eps)) else rep(1, r)
	Ms <- M / tcrossprod(s)
	vs <- v / s
	b  <- if (is.null(beta0)) rep(0, r) else as.numeric(beta0) * s

	l1 <- lambda * alpha * w
	l2 <- lambda * (1 - alpha) * w + ridge_reg

	### ridge: closed form
	if (alpha <= 0) {
		A  <- Ms + diag(l2, r)
		bs <- solve(A, vs)
		df <- sum(diag(solve(A, Ms)))
		return(list(beta = bs / s, df = df, iter = 0L))
	}

	### lasso and elastic net: cyclic coordinate descent
	dg <- diag(Ms)
	it <- 0L
	repeat {
		it <- it + 1L
		bo <- b
		for (j in seq_len(r)) {
			rj <- vs[j] - sum(Ms[j, ] * b) + dg[j] * b[j]
			b[j] <- if (l1[j] > 0) stem_soft(rj, l1[j]) / (dg[j] + l2[j])
			        else rj / (dg[j] + l2[j])
		}
		if (max(abs(b - bo)) < tol || it >= maxit) break
	}

	### Effective degrees of freedom. For the lasso the number of active
	### coefficients is unbiased for the degrees of freedom (Zou, Hastie and
	### Tibshirani 2007); for the elastic net the corresponding quantity is the
	### trace of the ridge projection restricted to the active set.
	act <- which(abs(b) > 0 | w == 0)
	df <- if (!length(act)) 0 else if (alpha >= 1) length(act) else {
		MA <- Ms[act, act, drop = FALSE]
		sum(diag(solve(MA + diag(l2[act], length(act)), MA)))
	}

	list(beta = b / s, df = df, iter = it)
}
