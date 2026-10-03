#' ML Estimation
#'
#'
#' @description The function \code{STEM_Estimation} computes the maximum likelihood estimates of the unknown parameters of a hierarchical spatio-temporal model of class \dQuote{STEM_Model}. The estimates are obtained using Kalman filtering and the EM algorithm, by default accelerated by SQUAREM (see \code{algorithm} in \code{\link{STEM_control}}).
#'
#' This is one of the two estimation engines of the package. The entry point is
#' \code{\link{STEM_Fit}}, which calls this function when \code{K = 1} and
#' \code{\link{SCSTEM_Estimation}} when \code{K > 1}, and which takes the same
#' arguments and returns the same object. Call \code{STEM_Estimation} directly
#' only to bypass the dispatch; it is the historical interface of the package
#' and is kept unchanged.
#'
#' @param StemModel an object of class \dQuote{STEM_Model} given as output by the \code{\link{STEM_Model}} function.
#' @param precision optional; when given, it replaces \code{em_tol_par} of
#'   \code{control}, the tolerance of the relative change of every free
#'   parameter. Kept for the code written before \code{control} existed.
#' @param max.iter optional; when given, it replaces \code{em_maxit} of
#'   \code{control}, the maximum number of EM iterations.
#' @param control the computational settings, an object returned by
#'   \code{\link{STEM_control}} or a list of some of its settings. Default
#'   \code{NULL}, the defaults of \code{STEM_control()}. This function uses
#'   \code{algorithm}, \code{em_tol_par}, \code{em_tol_loglik},
#'   \code{em_maxit}, \code{em_stop}, \code{nr_maxit} and \code{nr_hess_maxit}.
#' @param flag.Gdiag logical, indicating whether the transition matrix \eqn{G} is diagonal.
#' @param flag.Sigmaetadiag logical, indicating whether the variance-covariance matrix of the state equation \eqn{\Sigma_\eta} is diagonal.
#' @param cov.spat type of spatial covariance function. For the moment only the \emph{exponential} function is implemented.
#' @param distance character, indicating the type of distance. 'euclidean' compute euclidean distance while 'geo' compute the geodedic distance. use 'geo' only if the coordinates format is Longitude, Latitude. Default is 'euclidean'.
#' @param regularization a non-negative number added to the diagonal of the
#'   matrices the EM algorithm inverts (the smoother gains at the first and last
#'   time point, \eqn{\Sigma^*} in the update of \eqn{\sigma^2_\omega},
#'   \eqn{\Sigma_e} in the update of \eqn{\beta}, the blocks of the missing
#'   values and the Hessian of the Newton step). Default 0. With a positive value
#'   the iterations no longer converge to the maximum of the likelihood the
#'   filter computes; versions before 2.0.0 used 0.01, which biased the variance
#'   parameters by several percent.
#' @param verbose logical. If TRUE, the progress of the EM and Newton-Raphson iterations is reported through message(). Default is FALSE.
#' @param alpha the elastic-net mixing parameter for the regression coefficients, in \eqn{[0,1]}: \code{0} is ridge, \code{1} is the lasso, anything in between is the elastic net. Ignored when \code{lambda} is zero. Default is 0.
#' @param lambda the strength of the penalty on the regression coefficients. \code{0}, the default, gives the ordinary maximum likelihood estimator. See \code{\link{STEM_Fit}} for the interface that selects it.
#' @param penalize which regression coefficients the penalty acts on: a logical vector of length \eqn{r}, an index vector, or \code{NULL} (the default) for every coefficient except the intercept.
#' @param lambda_scale how \code{lambda} is measured. \dQuote{relative}, the default, puts the L1 part of the penalty on the scale of the largest partial gradient, so that \code{lambda} is free of the units of the response and \code{lambda * alpha >= 1} zeroes every penalized coefficient; \dQuote{absolute} applies it to the scaled normal equations directly. The L2 part is unaffected: it multiplies a curvature whose diagonal is already one, and is scale free either way.
#' @param latent logical. \code{FALSE} switches the latent process off, by setting the loading matrix to zero and holding \eqn{G}, \eqn{\Sigma_\eta} and \eqn{m_0} at their input values, which they are unidentified without it. Default is \code{TRUE}.
#' @param spatial logical. \code{FALSE} replaces the exponential correlation by the identity, so that \eqn{\Sigma_e} is a single variance and the Newton-Raphson step is skipped. Only the sum \eqn{\sigma^2_\varepsilon + \sigma^2_\omega} is then identified, and the range has no meaning. Default is \code{TRUE}.
#'
#'
#' @return The function returns an object of class \dQuote{STEM_Model} which is a list given by:
#' \itemize{
#' \item{skeleton} As the \code{skeleton} component of the \code{StemModel} object given in input.
#' \item{data} As the \code{data} component of the \code{StemModel} object given in input.
#' \item{estimates}A list of four objects: \code{phi.hat}, \code{y.smoothed}, \code{loglik}, \code{convergence.par} here described.
#' \code{phi.hat} is a list with the parameter ML estimates (\code{sigma2omega}, \code{beta}, \code{G}, \code{Sigmaeta}, \code{m0}, \code{C0}, \code{theta}, \code{sigma2eps}).
#' \code{y.smoothed} is a \code{ts} object (\eqn{n} by \eqn{p}) with the smoothed latent states at the estimates. \code{loglik} is the log-likelihood at the estimates.
#' \code{convergence.par} is a list with some information about the convergence of the algorithm: \code{conv.log} and \code{conv.par} are logical values
#' for the two convergence criteria described below, and \code{converged} says whether the stopping rule was met (either criterion with
#' \code{em_stop = "any"}, both with \code{em_stop = "all"}) rather than the limit on the iterations reached; \code{max.rel.par} and \code{delta.loglik} are the values
#' of the two criteria at the last check; \code{iterEM} is the number of iterations (EM or ECME steps, those of the extrapolations of SQUAREM included), \code{iterNR} the number of
#' Newton-Raphson iterations within each, \code{algorithm} the algorithm used and \code{control} the settings.
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
#' @section Algorithms:
#' Four algorithms are available, through \code{algorithm} in
#' \code{\link{STEM_control}}. They reach the same maximum of the likelihood
#' and differ in the number of iterations they take and in how close to the
#' maximum a given stopping rule leaves them.
#' \describe{
#'   \item{\code{"EM"}}{the EM algorithm: the E-step is the Kalman smoother, the
#'     M-step the sequence of conditional updates above.}
#'   \item{\code{"ECME"}}{the ECME algorithm (Liu and Rubin 1994): the same, but
#'     \eqn{\beta} and \eqn{m_0} maximize the observed log-likelihood given the
#'     other parameters, a generalized least-squares step computed by the Kalman
#'     filter run on the columns of the covariates. It removes the slow
#'     direction in which the intercept and \eqn{m_0} move together, at about
#'     25\% more time per iteration.}
#'   \item{\code{"SQUAREM"}}{(the default) the EM iterations accelerated by
#'     SQUAREM (Varadhan and Roland 2008): two EM steps are extrapolated along
#'     the direction of slowest convergence, and the extrapolation is kept only
#'     if it does not lower the log-likelihood. The steplength follows the R
#'     package SQUAREM (Du and Varadhan 2020).}
#'   \item{\code{"SQUAREM-ECME"}}{the ECME iterations accelerated by SQUAREM.}
#' }
#' With a penalty on \eqn{\beta} (\code{lambda > 0}) the safeguard of SQUAREM
#' has no fixed objective to rely on, since the scale of the penalty is
#' re-measured at every iteration, and the fit runs the plain iterations of EM
#' or ECME.
#'
#' @section Starting values and stopping rule:
#'   For initializing the algorithm the values contained in \code{StemModel$skeleton$phi} are used as initial values.
#'   The algorithm stops when either of the following two criteria (named in the output as \code{conv.par} and \code{conv.log} respectively) is met, or, with \code{em_stop = "all"} in \code{control}, when both are
#'
#'   \deqn{\max_l \frac{|\psi_l^{(i)} - \psi_l^{(i-1)}|}{\max(|\psi_l^{(i-1)}|, 10^{-3})} < \texttt{em\_tol\_par}}{max_l |psi_l(i) - psi_l(i-1)| / max(|psi_l(i-1)|, 1e-3) < em_tol_par}
#'
#'   \deqn{|\log L(\psi^{(i)}) - \log L(\psi^{(i-1)})| < \texttt{em\_tol\_loglik}}{|log L(psi(i)) - log L(psi(i-1))| < em_tol_loglik}
#'
#' or after \code{em_maxit} iterations, with the settings taken from
#' \code{control} (see \code{\link{STEM_control}}). The first criterion is the
#' relative change of every free parameter taken one at a time, on the scale on
#' which the algorithm updates it (the range and the ratio of the two variances
#' on the log scale); the loading matrix and \eqn{C_0}, which are not
#' estimated, are excluded. The second is the absolute change of the
#' log-likelihood, which, unlike a relative change, does not loosen as the
#' log-likelihood grows with the size of the data. Earlier versions of the
#' package stopped when the relative change of the whole parameter vector,
#' fixed loadings included, and the relative change of the log-likelihood were
#' both below \code{precision}; with large log-likelihoods that rule stopped
#' while several units could still be gained. With SQUAREM the criteria are
#' checked at the start of every cycle, on the EM step taken from there. The
#' estimates returned are the last iterate, and the log-likelihood and the
#' smoothed states are computed at them.
#'
#'
#'
#' @author Michela Cameletti  < michela.cameletti@unibg.it>
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
#' The algorithms:
#'
#' Du, Y., Varadhan, R. (2020) \emph{SQUAREM: An R package for off-the-shelf
#' acceleration of EM, MM and other EM-like monotone algorithms}. Journal of
#' Statistical Software, 92(7), 1--41. \doi{10.18637/jss.v092.i07}
#'
#' Liu, C., Rubin, D.B. (1994) \emph{The ECME algorithm: a simple extension of
#' EM and ECM with faster monotone convergence}. Biometrika, 81, 633--648.
#' \doi{10.1093/biomet/81.4.633}
#'
#' Varadhan, R., Roland, C. (2008) \emph{Simple and globally convergent methods
#' for accelerating the convergence of any EM algorithm}. Scandinavian Journal
#' of Statistics, 35, 335--353. \doi{10.1111/j.1467-9469.2007.00585.x}
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
#' A <-matrix(1,ncol(z),1)
#'
#' mod1 <- STEM_Model(z=z,covariates=covariates,
#'                    coordinates=coordinates,phi=phi,A=A)
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
function(StemModel, precision = NULL, max.iter = NULL, flag.Gdiag = TRUE, flag.Sigmaetadiag = TRUE,
         cov.spat = Sigmastar.exp, distance = "euclidean", regularization = 0, verbose = FALSE,
         alpha = 0, lambda = 0, penalize = NULL, lambda_scale = c("relative", "absolute"),
         latent = TRUE, spatial = TRUE, control = NULL)
{
  lambda_scale <- match.arg(lambda_scale)

  ### The computational settings: those of `control`, with the two historical
  ### arguments, when given, in place of the corresponding settings.
  control <- stem_control_resolve(control)
  if (!is.null(precision)) control$em_tol_par <- stem_control_resolve(list(em_tol_par = precision))$em_tol_par
  if (!is.null(max.iter)) control$em_maxit <- stem_control_resolve(list(em_maxit = max.iter))$em_maxit

  ### With no spatial correlation the exponential function is replaced by the
  ### identity, so that Sigma_e is a single variance and the Newton-Raphson step,
  ### which would be estimating a range that has nothing to estimate, is skipped.
  ### Only the SUM sigma2eps + sigma2omega is then identified: the split reported
  ### is whatever the starting value made it.
  if (!spatial && missing(cov.spat)) cov.spat <- Sigmastar.nugget

  dat <- stem_em_data(StemModel, distance = distance, cov.spat = cov.spat)
  opt <- list(Gdiag = flag.Gdiag, Sigmaetadiag = flag.Sigmaetadiag,
              regularization = regularization, latent = latent, spatial = spatial,
              alpha = alpha, lambda = lambda, pen_w = stem_penalized_index(dat$r, penalize),
              lambda_scale = lambda_scale, nr_maxit = control$nr_maxit,
              nr_hess_maxit = control$nr_hess_maxit, nr_tol = control$em_tol_par,
              verbose = verbose)

  ### The starting values, on the scale of the algorithm: the range and the
  ### ratio of the two variances on the log scale. Switching off the latent
  ### process sets the loading matrix to zero: the state then contributes
  ### nothing, and G, Sigmaeta and m0 are held at their starting values.
  sk <- StemModel$skeleton$phi
  phi <- list(A = if (latent) t(StemModel$skeleton$A) else matrix(0, dat$p, dat$d),
              sigma2omega = sk$sigma2omega,
              logtheta = log(sk$theta),
              logb = log(sk$sigma2eps / sk$sigma2omega),
              beta = sk$beta, G = sk$G, Sigmaeta = sk$Sigmaeta,
              m0 = matrix(sk$m0, dat$p, 1), C0 = sk$C0)

  fit <- stem_em_fit(phi, dat, opt, control)

  ### the Kalman filter and smoother at the estimates: the log-likelihood and
  ### the smoothed latent states reported
  fin <- kalman(fit$phi, dat, regularization)
  ph <- fit$phi
  StemModel$estimates$phi.hat <- list(sigma2omega = ph$sigma2omega, beta = ph$beta, G = ph$G,
                                      Sigmaeta = ph$Sigmaeta, m0 = ph$m0, C0 = ph$C0,
                                      theta = exp(ph$logtheta),
                                      sigma2eps = exp(ph$logb) * ph$sigma2omega)
  StemModel$estimates$y.smoothed <- fin$smoothed$m
  StemModel$estimates$loglik <- fin$loglik
  StemModel$estimates$convergence.par <- list(conv.log = fit$crit$conv_log,
                                              conv.par = fit$crit$conv_par,
                                              converged = fit$crit$done,
                                              max.rel.par = fit$crit$max_rel_par,
                                              delta.loglik = fit$crit$delta_loglik,
                                              iterEM = fit$iter,
                                              iterNR = fit$iterNR,
                                              algorithm = control$algorithm,
                                              control = control)
  ### The penalty in force and the effective number of regression coefficients
  ### it leaves. With lambda = 0 the second is simply r, and everything
  ### downstream -- the information criteria in particular -- behaves as it
  ### always has.
  StemModel$estimates$penalty <- list(alpha = alpha, lambda = lambda,
                                      lambda_scale = lambda_scale, penalize = penalize,
                                      beta.df = fit$last$beta_df,
                                      lambda.ref = fit$last$lambda_ref,
                                      lambda.eff = fit$last$lambda_eff)
  StemModel$estimates$scope <- list(latent = latent, spatial = spatial)
  StemModel
}
