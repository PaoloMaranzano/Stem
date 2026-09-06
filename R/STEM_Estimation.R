#' ML Estimation
#'
#'
#' @description The function \code{STEM_Estimation} computes the maximum likelihood estimates of the unknown parameters of a hierarchical spatio-temporal model of class \dQuote{STEM_Model}. The estimates are obtained using Kalman filtering and EM algorithm.
#'
#' @param StemModel an object of class \dQuote{STEM_Model} given as output by the \code{\link{STEM_Model}} function.
#' @param precision  a small positive number used for the algorithm convergence. Default is equal to 0.01. See \code{DETAILS} below.
#' @param max.iter maximum number of iterations for the EM algorithm. Default is equal to 50.
#' @param flag.Gdiag logical, indicating whether the transition matrix \eqn{G} is diagonal.
#' @param flag.Sigmaetadiag logical, indicating whether the variance-covariance matrix of the state equation \eqn{\Sigma_\eta} is diagonal.
#' @param cov.spat type of spatial covariance function. For the moment only the \emph{exponential} function is implemented.
#' @param distance character, indicating the type of distance. 'euclidean' compute euclidean distance while 'geo' compute the geodedic distance. use 'geo' only if the coordinates format is Longitude, Latitude. Default is 'euclidean'.
#' @param regularization a small positive number to be added to the digonal of the matrices matrices that need to be inverted . Default is set to 0.01
#' @param verbose logical. If TRUE, the progress of the EM and Newton-Raphson iterations is reported through message(). Default is FALSE.
#'
#'
#' @return The function returns an object of class \dQuote{STEM_Model} which is a list given by:
#' \itemize{
#' \item{skeleton} As the \code{skeleton} component of the \code{StemModel} object given in input.
#' \item{data} As the \code{data} component of the \code{StemModel} object given in input.
#' \item{estimates}A list of four objects: \code{phi.hat}, \code{y.smoothed}, \code{loglik}, \code{convergence.par} here described.
#' \code{phi.hat} is a list with the parameter ML estimates (\code{sigma2omega}, \code{beta}, \code{G}, \code{Sigmaeta}, \code{m0}, \code{C0}, \code{theta}, \code{sigma2eps}).
#' \code{y.smoothed} is a \code{ts} object (\eqn{n} by \eqn{p}) which is the output of the Kalman filtering procedure. \code{loglik} is the log-likehood value.
#' \code{convergence.par} is a list of 4 objects with some information about the convergence of the algorithm: \code{conv.log} and \code{conv.par} are logical values
#' for the two convergence criteria described above; \code{iterEM} is the number of iterations for the EM algorithm and \code{iterNR} is the number of
#' Newton-Raphson iterations for each EM algorithm iteration.
#' }
#'
#'
#' @details This function estimates the vector parameter \code{phi} of the hierarchical spatio-temporal model of class \dQuote{STEM_Model} using Kalman filtering and EM algorithm.
#'
#' @section Computational form of the filter:
#' The forward pass does not build the \eqn{d \times d} matrix
#' \eqn{Q_t = Z R_t Z' + \Sigma_e} that the prediction step nominally requires,
#' and does not invert it. Two facts make that avoidable. The measurement
#' covariance \eqn{\Sigma_e} does not depend on \eqn{t}: it is rebuilt once per
#' EM iteration. And \eqn{Z R_t Z'} has rank \eqn{p}, the dimension of the
#' latent state, which is usually one. The Woodbury identity therefore gives
#' \deqn{Q_t^{-1} = \Sigma_e^{-1} - U (R_t^{-1} + Z' \Sigma_e^{-1} Z)^{-1} U',
#'       \qquad U = \Sigma_e^{-1} Z,}
#' and the matrix determinant lemma gives
#' \deqn{\log|Q_t| = \log|\Sigma_e| + \log|R_t| +
#'       \log|R_t^{-1} + Z' \Sigma_e^{-1} Z|.}
#' Here \eqn{U}, \eqn{Z' \Sigma_e^{-1} Z}, the Cholesky factor of
#' \eqn{\Sigma_e} and \eqn{\log|\Sigma_e|} are constants of the pass and are
#' computed once. Every quantity the recursion needs -- the gain applied to the
#' innovation, the updated state variance, the quadratic form of the
#' log-likelihood -- then reduces to \eqn{p} by \eqn{p} algebra plus one
#' triangular solve, so the cost per time point falls from \eqn{O(d^3)} to
#' \eqn{O(d^2)} and the cost of a pass from \eqn{O(T d^3)} to
#' \eqn{O(d^3 + T d^2)}. The two forms are algebraically identical and agree to
#' floating point. The identity holds for any \eqn{p}, and needs only that
#' \eqn{\Sigma_e} and \eqn{R_t} be invertible, which positive definiteness
#' already guarantees; the gain is largest when \eqn{p} is small relative to
#' \eqn{d}, which is the regime the model is written for.
#'
#' With missing values the constants depend on which rows are observed, so they
#' are computed once per distinct missingness pattern and cached.
#'
#' The log-density is evaluated directly rather than as the logarithm of the
#' density. On \eqn{d} observations the Gaussian density is of order
#' \eqn{e^{-d}}, so on a network of a few hundred locations it falls below the
#' smallest representable double and taking its logarithm afterwards returns
#' \code{-Inf}.
#'
#' @section Computational form of the M-step:
#' The sums the M-step accumulates over time are written as matrix products
#' rather than as loops over \eqn{t}. Because the loading matrix does not depend
#' on time, the second moment of the measurement error collapses to
#' \deqn{\sum_t \{ Z C^s_t Z' + r_t r_t' \} =
#'       Z (\sum_t C^s_t) Z' + R'R,}
#' with \eqn{R} the \eqn{T} by \eqn{d} matrix of residuals, so one cross-product
#' replaces \eqn{T} outer products and the \eqn{T} matrices of size \eqn{d} by
#' \eqn{d} that used to be held at once. The three sums of outer products of the
#' smoothed states are cross-products of the \eqn{T} by \eqn{p} matrix of those
#' states, and the two accumulations entering the update of \eqn{\beta} are
#' cross-products of the design blocks stacked by period. Finally, every trace
#' of a matrix product is evaluated as \eqn{tr(AB) = \sum_{ij} A_{ij} B_{ji}},
#' which costs \eqn{O(d^2)} instead of forming the product.
#' The algorithm details and formulas are given in Fasso' and Cameletti (2007, 2009). Note that some parameters (\code{beta}, \code{sigma2omega},
#'   \code{G}, \code{Sigmaeta} and \code{m0}) are updated using closed form solutions while \code{theta} and \code{sigma2epsilon} using the Newton-Raphson algorithm.
#'
#'   For initializing the algorithm the values contained in \code{StemModel$skeleton$phi} are used as initial values.
#'   The algorithm converges when the following convergence criteria (named in the output as \code{conv.par} and \code{conv.log} respectively) are jointly met
#'
#'   \deqn{\frac{\left\|\phi^{\left(i+1\right)}-\phi^{\left(i\right)}
#'   \right\|}{\left\|\phi^{\left(i\right)}\right\|} < \pi }{||\phi^{i+1}-\phi^{i}
#'   || / ||\phi^{i}|| < \pi }
#'
#'   \deqn{\frac{\left\|\log L\left( \phi^{\left(i+1\right)}\right)-\log
#'   L\left( \phi^{\left(i\right)}\right)\right\|}{\left\|\log
#'   L\left( \phi^{\left(i\right)}\right)\right\|}<\pi}{||log L(\phi^{i+1}-log
#'   L(\phi^{i})|| / ||log L(\phi^{i})||<\pi}
#'
#' where \eqn{\pi} is given by the \code{precision} option and \eqn{i} is the number of iteration. The use of these relative criteria instead of some other absolute ones makes it possible to
#' correct for the different parameter scales.
#'
#'
#'
#' @author Michela Cameletti  < michela.cameletti@unibg.it >
#'
#' @references
#' The model estimated by this function, its EM algorithm and the
#' parametric bootstrap are those of the following three companion works.
#'
#' Fasso, A., Cameletti, M., Nicolis, O. (2007) \emph{Air quality monitoring
#' using heterogeneous networks}. Environmetrics, 18, 245--264.
#' \doi{10.1002/env.837}
#'
#' Fasso, A., Cameletti, M. (2007) \emph{A general spatio-temporal model for
#' environmental data}. GRASPA Technical Report n. 27, The Italian Group of
#' Environmental Statistics. The reference that introduces this package.
#'
#' Fasso, A., Cameletti, M. (2010) \emph{A unified statistical approach for
#' simulation, modeling, analysis and mapping of environmental data}.
#' Simulation, 86, 139--153. \doi{10.1177/0037549709102150}
#'
#' Further background:
#'
#' Amisigo, B.A., Van De Giesen, N.C. (2005) \emph{Using a spatio-temporal
#' dynamic state-space model with the EM algorithm to patch gaps in daily
#' riverflow series}. Hydrology and Earth System Sciences, 9, 209--224.
#'
#' McLachlan, G.J., Krishnan, T. (2008) \emph{The EM Algorithm and
#' Extensions}, 2nd edition. Wiley, New York.
#'
#' Shumway, R.H., Stoffer, D.S. (2006) \emph{Time Series Analysis and Its
#' Applications: with R Examples}. Springer, New York.
#'
#' Xu, K., Wikle, C.K. (2007) \emph{Estimation of parameterized
#' spatio-temporal dynamic models}. Journal of Statistical Inference and
#' Planning, 137, 567--588.
#'
#'
#' @examples
#'
#' #load the data
#' data(pm10)
#'
#' #extract the data
#' coordinates <- pm10$coords
#' covariates <- pm10$covariates
#' z <- pm10$z
#'
#' #build the parameter list
#' #(the phi list is used for the algorithm starting values)
#' phi <- list(beta=matrix(c(3.65,0.046,-0.904),3,1),
#'             sigma2eps=0.1,
#'             sigma2omega=0.2,
#'             theta=0.01,
#'             G=matrix(0.77,1,1),
#'             Sigmaeta=matrix(0.3,1,1),
#'             m0=as.matrix(0),
#'             C0=as.matrix(1))
#'
#' K <-matrix(1,ncol(z),1)
#'
#' mod1 <- STEM_Model(z=z,covariates=covariates,
#'                    coordinates=coordinates,phi=phi,K=K)
#' class(mod1)
#'
#' #mod1 is given as output by the STEM_Model function
#' mod1.est <- STEM_Estimation(mod1)
#' phi.estimates <- unlist(mod1.est$estimates$phi.hat)
#'
#'
#' @seealso See Also \code{\link{STEM_Model}} and \code{\link{pm10}}
#'
#' @keywords models spatial
#'
#'
#' @export


STEM_Estimation <-
function(StemModel, precision=0.01, max.iter=50,flag.Gdiag=TRUE,flag.Sigmaetadiag=TRUE,cov.spat=Sigmastar.exp,distance="euclidean",regularization=0.01, verbose = FALSE)
{

z 		=  StemModel$data$z

p 		= StemModel$skeleton$p
n 		= StemModel$data$n
d 		= StemModel$data$d
r 		= StemModel$data$r

covariates 	= StemModel$data$covariates
covariates  = changedimension_covariates(covariates,d=d,r=r,n=n)
coordinates = StemModel$data$coordinates

phi_start 	= StemModel$skeleton$phi
phi_start$K = t(StemModel$skeleton$K)
n_par     	= length(unlist(phi_start))

phi_start$logb 		= log(phi_start$sigma2eps/phi_start$sigma2omega)
phi_start$logtheta	= log(phi_start$theta)


####################
###EM algorithm while loop
####################
#output matrices (two columns are added for number of iteration and -2loglik)
parameters_mat   = matrix(0,nrow=max.iter,ncol=n_par+2)
#colnames(parameters_mat) = c("n_iter","-2loglik",rep("k",dd*pp),"sigma2omega", "theta","logb",rep("beta",6),rep("G",pp*pp),rep("Sigmaeta",pp*pp),rep("m0",pp),rep("C0",pp*pp))
distance_mat    	= matrix(0,nrow=max.iter,ncol=1)
distancelog_mat	= c()
Q_mat       	= matrix(0,nrow=max.iter,ncol=2)
iterNR 		= c()

converged_EM_1 	= FALSE
converged_EM_2 	= FALSE
n_iter_EM     	= 1
step_last     	= NULL

	if (isTRUE(verbose)) message("**** EM Algorithm - iteration n. ", n_iter_EM)
while ((!converged_EM_1 | !converged_EM_2) && n_iter_EM < max.iter){
	step = kalman(	z            		= z,
			coordinates  	= coordinates,
         		p           		= p,
			n			= n,
			d			= d,
			r			= r,
			phi_j        		= phi_start,
			max.iter     	= max.iter,
			precision       	= precision,
			covariates   	= covariates,
			Gdiag        		= flag.Gdiag,
			Sigmaetadiag 	= flag.Sigmaetadiag,
			cov.spat		= cov.spat,
			distance = distance,
			regularization=regularization,
			verbose = verbose
	)

	iterNR[n_iter_EM] 	= step$n_iter_NR
	step$phi$loglik 		= -2*step$phi$loglik
	par           		= t(matrix(unlist(step$phi)))

	### Robustness: an EM step that returns non-finite parameters (which the
	### inner Newton-Raphson can produce on small or nearly collinear subsets of
	### locations) used to propagate NaN into the convergence tests and abort
	### the fit with "missing value where TRUE/FALSE needed". The iteration is
	### now discarded and the last valid iterate is returned instead, so that a
	### clusterwise algorithm calling this function on many candidate subsets
	### degrades gracefully rather than failing.
	if(!all(is.finite(par))) {
		if(n_iter_EM == 1) {
			stop("The EM algorithm produced non-finite parameters at the first iteration: check the starting values in 'phi' and the conditioning of the data.", call. = FALSE)
		}
		warning("The EM algorithm produced non-finite parameters at iteration ", n_iter_EM,
			"; the last valid iterate is returned.", call. = FALSE)
		break
	}

	parameters_mat[n_iter_EM,] = cbind(n_iter_EM, par)
	step_last = step

	if(n_iter_EM==1) {
  		prev_lik = 0
  		prev_par = rep(0,n_par)
   	} else {
  		prev_lik = (parameters_mat[n_iter_EM-1, 2])          #second column for -2loglik
  		prev_par = parameters_mat[n_iter_EM-1, -c(1,2)]   #no n_iter e -2loglik
  	}

  	if(n_iter_EM==1 | n_iter_EM==2) {
  		media_lik = 1
  	} else {
  		media_lik = mean(c(step$phi$loglik, prev_lik))
  	}

	dist_rel_num = sqrt(t(parameters_mat[n_iter_EM,-c(1,2)] - unlist(prev_par)) %*% (parameters_mat[n_iter_EM,-c(1,2)] - unlist(prev_par)))
	### Floor on the denominators: both relative criteria are undefined at the
	### first iteration, where the reference vector is exactly zero.
	dist_rel_den = max(sqrt(t(unlist(prev_par)) %*% unlist(prev_par)), .Machine$double.eps)
	dist_rel = dist_rel_num / dist_rel_den
	distance_mat[n_iter_EM,] = dist_rel

	diff_rel_loglik = abs(step$phi$loglik - prev_lik) / max(abs(prev_lik), .Machine$double.eps)
	distancelog_mat[n_iter_EM] = diff_rel_loglik

	###Check the convergence!
	### isTRUE() so that a non-finite criterion counts as "not converged"
	### instead of turning the loop condition into NA.
	converged_EM_1 = isTRUE(as.logical(diff_rel_loglik < precision))
	converged_EM_2 = isTRUE(as.logical(dist_rel < precision))

	Q_mat[n_iter_EM,] = c(unlist(step$Q_prev),unlist(step$Q_new))

	###Updating the iteration number and the parameter vector
	phi_start 	= step$phi[-1]
	n_iter_EM 	= n_iter_EM + 1
	if (isTRUE(verbose)) message("**** EM Algorithm - iteration n. ", n_iter_EM)

} # here the while loop ends

phi_start$theta = exp(phi_start$logtheta)
phi_start$sigma2eps = exp(phi_start$logb) * phi_start$sigma2omega
phi_start = phi_start[-which(names(phi_start) == "logtheta")]
phi_start = phi_start[-which(names(phi_start) == "K")]
phi_start = phi_start[-which(names(phi_start) == "logb")]
phi.estimated= phi_start

StemModel$estimates$phi.hat = phi.estimated
### the last VALID Kalman step, which differs from the last attempted one when
### the loop broke out on a non-finite iterate
StemModel$estimates$y.smoothed = step_last$m.smoother
StemModel$estimates$loglik = (parameters_mat[(n_iter_EM-1),2])*(-2)
convergence.par 			= list(conv.log = converged_EM_1,
						conv.par = converged_EM_2,
						iterEM   = n_iter_EM-1,
						iterNR   = iterNR)
StemModel$estimates$convergence.par = convergence.par
return (StemModel)
}

