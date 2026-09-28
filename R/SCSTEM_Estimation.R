#' Spatially-clustered STEM estimation
#'
#' @description
#' \code{SCSTEM_Estimation} fits a spatially-clustered spatio-temporal
#' expectation-maximization (SC-STEM) model. The \eqn{d} monitoring locations
#' are partitioned into \eqn{k} latent spatial regimes, and a separate
#' \dQuote{STEM_Model} is estimated within each regime, so that regression
#' coefficients, variance components and latent temporal dynamics are all
#' cluster-specific. The partition and the parameters are estimated jointly by
#' alternating optimization of a Potts-penalized log-likelihood.
#'
#' This is one of the two estimation engines of the package. The entry point is
#' \code{\link{STEM_Fit}}, which calls this function when \code{k > 1} and
#' \code{\link{STEM_Estimation}} when \code{k = 1}, and which takes the same
#' arguments and returns the same object. Call \code{SCSTEM_Estimation} directly only
#' to bypass the dispatch.
#'
#' @details
#' \strong{Model.} Conditionally on location \eqn{i} belonging to cluster
#' \eqn{k}, the SC-STEM model is the cluster-specific STEM model
#' \deqn{z_{it} = x_{it}' \beta_k + K_i y^{(k)}_t + e_{it}, \qquad
#'       y^{(k)}_t = G_k y^{(k)}_{t-1} + \eta^{(k)}_t,}
#' with \eqn{e_t \sim N(0, \Sigma_{e,k})},
#' \eqn{\Sigma_{e,k} = \sigma^2_{\epsilon k} I + \sigma^2_{\omega k} C(h;\theta_k)}
#' and \eqn{\eta^{(k)}_t \sim N(0, \Sigma_{\eta k})}. Setting \eqn{k = 1}
#' recovers the pooled \code{\link{STEM_Estimation}} fit.
#'
#' \strong{Objective.} Labels \eqn{k_1,\ldots,k_d} and parameters are estimated
#' by maximizing the penalized log-likelihood
#' \deqn{Q = \sum_{i=1}^{d} \ell_i(k_i) +
#'           \phi\, c \sum_{(i,j) \in E} I(k_i = k_j),}
#' where \eqn{\ell_i(k)} is the log-likelihood contribution of location \eqn{i}
#' under the parameters of cluster \eqn{k}, \eqn{E} is the edge set of the
#' symmetrized \code{knn} graph with each unordered pair counted once, and
#' \eqn{c > 0} is the scale factor discussed
#' below. This is the Potts-type penalty introduced for spatially-clustered
#' regression by Sugasawa and Murakami (2021) and carried over to
#' spatially-clustered spatial autoregressive models by Cerqueti, Maranzano and
#' Mattera (2025); here it is applied to the spatio-temporal likelihood of a
#' STEM model.
#'
#' \strong{Scale of the penalty.} In cross-sectional spatially-clustered models
#' each unit contributes a single observation to the likelihood, so that
#' \eqn{\phi} of order one balances fit against spatial cohesion. Here each
#' location contributes \eqn{T} observations, and the log-likelihood differences
#' between clusters grow with \eqn{T} and with the scale of the response. The
#' factor \eqn{c} restores a comparable interpretation:
#' \itemize{
#'   \item \code{phi_scale = "auto"} (default) sets \eqn{c} to the median across
#'     locations of the spread \eqn{\max_k \ell_{ik} - \min_k \ell_{ik}},
#'     divided by the average number of neighbors. With this normalization
#'     \eqn{\phi = 1} is the point at which full agreement with the
#'     neighborhood is worth about as much as the typical gain from picking the
#'     best-fitting cluster. The factor is computed once, at the first sweep, and
#'     is returned in \code{phi_multiplier}. Since a penalized fit starts from
#'     the unpenalized solution (see below), the first sweep scores the
#'     locations with the parameters of the unpenalized fit: the factor is
#'     measured on the fitted clusters, does not depend on \code{init_method},
#'     and is of the order of the whole gain of the right regime over the
#'     wrong ones, so the informative range of \eqn{\phi} is small, a grid such
#'     as \eqn{\phi \in \{0, 0.025, 0.05, 0.1, 0.2\}}.
#'   \item \code{phi_scale = "per-observation"} sets \eqn{c = T}, which makes
#'     \eqn{\phi} invariant to the length of the series and directly comparable
#'     with the cross-sectional literature, but leaves it dependent on the scale
#'     of the response.
#'   \item \code{phi_scale = "raw"} sets \eqn{c = 1}, penalizing on the
#'     untransformed likelihood scale.
#' }
#' The effective penalty actually applied is always reported in
#' \code{phi_effective}.
#'
#' \strong{Penalized fits start from the unpenalized one.} With
#' \code{phi_penalty > 0} the fit with \code{phi_penalty = 0} at the same
#' \eqn{k} is run first, and its solution is the starting partition. A penalty
#' that is strong from the first sweep freezes whatever partition it is given,
#' so a start unrelated to the regimes would stay where it is; from the
#' unpenalized solution the penalty only has to decide how much spatial
#' smoothing that solution can afford. \code{\link{SCSTEM_Infocrit}} fits
#' \eqn{\phi = 0} once per \eqn{k} and passes its partition on, so a single fit
#' and the corresponding member of a grid coincide.
#'
#' \strong{Algorithm.} The two steps are iterated until convergence:
#' \enumerate{
#'   \item \emph{Parameter update given labels.} A STEM model is estimated by
#'     \code{\link{STEM_Estimation}} on each cluster. A cluster that falls below
#'     the minimum admissible size keeps its last valid parameters
#'     (\dQuote{stale freeze}) instead of aborting the algorithm; a cluster that
#'     never obtained a valid fit is excluded from the assignment step.
#'   \item \emph{Label update given parameters.} With
#'     \code{label_update = "ICM"} (the default) locations are visited
#'     sequentially and each label maximizes its own penalized contribution
#'     given the current labels of all the others, in the spirit of the
#'     Iterated Conditional Modes algorithm of Besag (1986). For fixed
#'     parameters this sweep cannot decrease \eqn{Q}, which rules out the label
#'     cycling that a simultaneous update can produce.
#'     \code{label_update = "simultaneous"} reproduces the joint update of all
#'     labels used in earlier versions of the package and in Sugasawa and
#'     Murakami (2021).
#' }
#' The loop stops when the partition is unchanged, when the improvement of
#' \eqn{Q} falls below \code{abs_tol}/\code{rel_tol}, when a previously visited
#' partition reappears, or when \code{max_iter} is reached.
#'
#' \strong{On monotonicity.} The label step is monotone at fixed parameters, but
#' the alternation as a whole is \emph{not} guaranteed to increase \eqn{Q}
#' monotonically, and \code{obj_trace} may well show a decrease. The reason is
#' structural rather than numerical. In the STEM measurement equation the
#' spatial covariance \eqn{\Sigma_{e,k}} couples the locations, so the exact
#' marginal likelihood does not factorize across them and the assignment score
#' \eqn{\ell_{ik}} has to be a conditional pseudo-likelihood, whereas the
#' parameter step maximizes the exact within-cluster likelihood through the EM
#' algorithm. The two objectives agree on what a good partition looks like but
#' are not the same function, so a parameter update can lower \eqn{Q} while
#' raising the exact likelihood. This is a property of the spatio-temporal
#' specification, and it is why the best partition visited along the iterations
#' is always the one returned.
#'
#' \strong{Degeneracy and the minimum-size constraint.} Within-cluster
#' homogeneity is exactly what the assignment step seeks, so the cluster-wise
#' variance components \eqn{\sigma^2_{\epsilon k} + \sigma^2_{\omega k}} shrink
#' as locations are reallocated. Left unconstrained, the cluster with the
#' smallest residual variance then attracts every location and the partition
#' collapses -- the clusterwise counterpart of the degenerate-likelihood problem
#' of Gaussian mixtures. With \code{enforce_min_size = TRUE}
#' (the default) a location may leave its cluster only if that cluster stays at
#' or above \code{min_cluster_size}, which keeps every visited configuration
#' admissible while preserving the monotonicity of the ICM sweep within the
#' feasible set. The constraint matters most when the number of locations is
#' small relative to \eqn{k}; it can be lifted with
#' \code{enforce_min_size = FALSE}.
#'
#' The constraint has a cost of its own: a location sitting in a cluster that is
#' exactly at the minimum size can never move, so with small networks a large
#' share of the partition may stay frozen at its initial value. Each sweep is
#' therefore followed by a \emph{size-preserving swap pass}
#' (\code{swap_pass = TRUE}, the default): the labels of two locations in
#' different clusters are exchanged whenever the exchange strictly increases
#' \eqn{Q}. A swap leaves every cluster size unchanged, so feasibility holds by
#' construction, and it is accepted only on strict improvement, so the
#' monotonicity of the algorithm is preserved. The number of accepted swaps is
#' reported in \code{obj_trace}.
#'
#' \strong{Final refit and information criteria.} On convergence the
#' cluster-wise STEM models are re-estimated once on the final partition. All
#' reported coefficients, variance components and information criteria come
#' from this refit and are based on the \emph{exact} cluster-wise
#' log-likelihoods returned by \code{\link{STEM_Estimation}}, not on the
#' pseudo-likelihood used to rank clusters during the assignment step. The
#' number of free parameters is \eqn{k_{eff} (r + 3 + 3p)} for a diagonal
#' specification with \eqn{r} covariates and latent dimension \eqn{p}, where
#' \eqn{k_{eff}} counts the clusters that could actually be re-estimated.
#'
#' \strong{Statistical features and scope.} The assumptions are inherited from
#' the STEM model and determine which datasets the family applies to. The
#' response is \emph{Gaussian} and \emph{univariate}: \code{z} is a \eqn{T} by
#' \eqn{d} matrix of one variable measured at \eqn{d} sites, so several
#' pollutants modeled jointly is a different specification. What is multivariate
#' is the \emph{latent state}, of dimension \eqn{p \ge 1}, loaded onto the
#' locations by the known matrix \eqn{K}, which is not estimated and is common
#' across regimes. The latent dynamics is a VAR(1), with \eqn{G} and
#' \eqn{\Sigma_\eta} diagonal by default; \eqn{m_0} is estimated and \eqn{C_0}
#' is held fixed. The spatial correlation function is exponential, hence
#' isotropic and stationary \emph{within} a regime, the partition itself being
#' the only source of non-stationarity across the domain. Locations are
#' point-referenced and fixed over time, and the time index is discrete and
#' regularly spaced.
#'
#' \strong{Missing values} are supported in the \emph{response}, and handled as
#' prescribed by Durbin and Koopman (2012, 2nd ed.), Sections 2.7 and 4.10: at
#' each time point the measurement equation is restricted to the locations
#' actually observed, through a selection matrix whose rows are a subset of the
#' rows of the identity, and a time point at which nothing is observed
#' contributes no update and no likelihood term. Because the EM algorithm
#' maximizes the expected complete-data log-likelihood, the M-step completes the
#' sufficient statistics: a missing value enters through its conditional
#' expectation given everything observed, and its conditional variance is added
#' back as a correction. Note that this conditional expectation is not the
#' signal alone, since the spatial covariance couples the locations, so the
#' missing block of the measurement error is predicted from the observed one by
#' the same algebra as kriging at a fixed time point. In the assignment step,
#' each location is scored on the time points at which it was observed.
#'
#' Missing values are \emph{not} supported in the covariates or the coordinates:
#' the design matrix enters the closed-form M-step updates directly and the
#' coordinates enter the distance matrix, so gaps there would require a
#' stochastic E-step. Impute them before fitting. Every location must retain at
#' least one observation.
#'
#' @param StemModel an object of class \dQuote{STEM_Model} given as output by
#'   the \code{\link{STEM_Model}} function.
#' @param k integer, the number of spatial clusters. \code{k = 1} returns the
#'   pooled STEM fit. Default is 3.
#' @param phi_penalty non-negative number, the weight of the Potts spatial
#'   penalty. \code{phi_penalty = 0} gives non-spatial clusterwise STEM.
#'   Default is 0.05.
#' @param phi_scale character, one of \code{"auto"} (default),
#'   \code{"per-observation"} or \code{"raw"}, setting the scale factor of the
#'   spatial penalty. See \code{Details}.
#' @param knn integer, the number of nearest neighbors used to build the
#'   spatial penalty graph. The graph is symmetrized, and the neighbors are
#'   measured with the metric given by \code{distance}, so that the penalty and
#'   the covariance see the same geometry. Default is 5.
#' @param distance character, \code{"euclidean"} for Euclidean distance or
#'   \code{"geo"} for geodesic distance. It governs both the covariance of the
#'   measurement error and the nearest neighbors of the penalty graph. Use
#'   \code{"geo"} only when the coordinates are longitude/latitude. Default is
#'   \code{"geo"}.
#' @param init_method character, one of \code{"departures"} (default),
#'   \code{"kmeans"} or \code{"coordinates"}, the starting partition of an
#'   unpenalized fit. \code{"departures"} fits the pooled model and runs
#'   k-means on how each location departs from it: the mean of its residual
#'   from the pooled signal, the slopes of that residual on its covariates, and
#'   the lag-one autocorrelation and log-variance of what the slopes leave.
#'   \code{"kmeans"} runs k-means on the location-wise covariate means, which
#'   carry no information on the regimes when the covariates are exogenous to
#'   them; \code{"coordinates"} clusters the spatial coordinates, which is also
#'   the fallback for intercept-only models. In every case the features are
#'   compressed by principal components, with multiple restarts and a
#'   minimum-cluster-size filter.
#' @param init_partition optional integer vector of length \eqn{d} giving a
#'   starting partition, overriding \code{init_method}. Use it to supply an
#'   externally computed initialization; for instance the AMKM partition used
#'   by versions of the package before 2.0.0 can be reproduced by passing
#'   \code{SCDA::SC_AMKM(...)$df$cluster}. The partition is repaired if it
#'   violates \code{min_cluster_size}. With \code{phi_penalty > 0} and no
#'   \code{init_partition}, the fit starts from the solution of the
#'   unpenalized fit at the same \eqn{k}, itself started from
#'   \code{init_method}; a given \code{init_partition} is taken as that start,
#'   which is how \code{\link{SCSTEM_Infocrit}} passes the partition of its
#'   unpenalized fit.
#' @param label_update character, \code{"ICM"} (default) for the sequential
#'   Iterated Conditional Modes sweep, or \code{"simultaneous"} for the joint
#'   update of all labels.
#' @param precision small positive number, the convergence tolerance of the EM
#'   algorithm in each cluster-wise fit. Default is 0.1.
#' @param precision_full_dataset small positive number, the convergence
#'   tolerance of the EM algorithm for the pooled fit used to initialize the
#'   procedure. Default is 0.01.
#' @param regularization small positive number added to the diagonal of the
#'   matrices that have to be inverted. Default is 0.01.
#' @param max_iter integer, the maximum number of alternating iterations.
#'   Default is 10.
#' @param abs_tol,rel_tol absolute and relative tolerances on the improvement
#'   of the penalized objective. Defaults are 1e-5 and 1e-6.
#' @param min_cluster_size integer or \code{NULL}. Minimum number of locations
#'   required to estimate a cluster-wise model. When \code{NULL} (default) it is
#'   set to \code{ncov + 2}.
#' @param enforce_min_size logical. If \code{TRUE} (default) the assignment
#'   step never lets a cluster fall below \code{min_cluster_size}. See
#'   \code{Details}. Set to \code{FALSE} to reproduce the unconstrained sweep.
#' @param swap_pass logical. If \code{TRUE} (default) each sweep is followed
#'   by a size-preserving swap pass that exchanges the labels of pairs of
#'   locations whenever this increases the objective. See \code{Details}.
#' @param share2conv number in \eqn{[0,1)}. Optional early-stopping rule kept
#'   for backward compatibility: the loop also stops when the share of
#'   locations changing cluster falls below this value. Set to 0 (default) to
#'   rely only on the objective-based criteria.
#' @param seed integer or \code{NULL}, seed used for the initialization step so
#'   that the fit is reproducible. Default is 123456789.
#' @param alpha the elastic-net mixing parameter for the regression
#'   coefficients, in \eqn{[0,1]}: \code{0} is ridge, \code{1} the lasso,
#'   anything between the elastic net. The penalty acts within each regime.
#'   Ignored when \code{lambda} is zero. Default is 0. See
#'   \code{\link{STEM_Fit}}.
#' @param lambda the strength of the penalty on the regression coefficients.
#'   \code{0}, the default, gives the unpenalized estimator, and the information
#'   criteria then count the nominal number of coefficients as before; with
#'   \code{lambda > 0} they count the effective number the penalty leaves.
#' @param penalize which coefficients the penalty acts on: a logical vector of
#'   length \eqn{r}, an index vector, or \code{NULL} (the default) for every
#'   coefficient except the intercept.
#' @param lambda_scale how \code{lambda} is measured, \dQuote{relative} (the
#'   default) or \dQuote{absolute}. See \code{\link{STEM_Estimation}}.
#' @param lambda_by how the penalty is distributed over the regimes.
#'   \dQuote{common}, the default, gives every regime the same \code{lambda},
#'   which under the relative scale already means the same proportional
#'   shrinkage, since the reference is computed inside each regime.
#'   \dQuote{size} sets \eqn{\lambda_g = \lambda \bar n / n_g} with
#'   \eqn{\bar n = d/k}, shrinking a regime of half the average size twice as
#'   hard. Either way \code{lambda} stays ONE hyperparameter: genuinely
#'   cluster-specific \eqn{(\alpha_g, \lambda_g)} is a different model and is
#'   not offered here.
#' @param latent logical, passed to \code{\link{STEM_Estimation}} within each
#'   regime. Default is \code{TRUE}.
#' @param spatial logical, passed to \code{\link{STEM_Estimation}} within each
#'   regime. Default is \code{TRUE}.
#' @param verbose logical. If \code{TRUE}, progress information is emitted via
#'   \code{message()}. Default is \code{FALSE}.
#'
#' @return An object of class \dQuote{SCSTEM_Estimation}, a list with components:
#' \itemize{
#'   \item \code{phi_hat}: \eqn{k} by \eqn{npar} matrix of cluster-wise
#'     parameter estimates from the final refit.
#'   \item \code{group}: integer vector of length \eqn{d} with the estimated
#'     cluster label of each location.
#'   \item \code{df}: data frame with the coordinates and the estimated labels.
#'   \item \code{fit_list}: list of the \dQuote{STEM_Model} objects returned by
#'     \code{\link{STEM_Estimation}} on the final partition.
#'   \item \code{idx_g}: list of the location indices of each cluster.
#'   \item \code{info_crit}: named vector with the total log-likelihood, the
#'     number of free parameters and AIC, BIC and KIC.
#'   \item \code{loglik_g}: cluster-wise exact log-likelihoods.
#'   \item \code{final_refit}: logical vector flagging the clusters that could
#'     be re-estimated on the final partition.
#'   \item \code{obj_trace}: data frame tracing the penalized objective, the
#'     number of label changes and the cluster sizes along the iterations.
#'   \item \code{convergence}: character describing the exit route.
#'   \item \code{penalized_obj}: value of the penalized objective at the exit.
#'   \item \code{input_args}: the arguments used for the fit, needed by
#'     \code{\link{SCSTEM_Bootstrap}} and \code{\link{SCSTEM_Select}}.
#' }
#'
#' @author Paolo Maranzano \email{pmaranzano.ricercastatistica@gmail.com},
#'   Francesco Caccia, Michela Cameletti
#'
#' @references
#' Besag, J. (1986) \emph{On the statistical analysis of dirty pictures}.
#' Journal of the Royal Statistical Society, Series B, 48, 259--302.
#'
#' Cerqueti, R., Maranzano, P., Mattera, R. (2025) \emph{Spatially-clustered
#' spatial autoregressive models with application to agricultural market
#' concentration in Europe}. Journal of Agricultural, Biological and
#' Environmental Statistics. \doi{10.1007/s13253-025-00685-7}
#'
#' Fasso, A., Cameletti, M., Nicolis, O. (2007) \emph{Air quality monitoring
#' using heterogeneous networks}. Environmetrics, 18, 245--264.
#' \doi{10.1002/env.837}
#'
#' Fasso, A., Cameletti, M. (2010) \emph{A unified statistical approach for
#' simulation, modeling, analysis and mapping of environmental data}.
#' Simulation, 86, 139--153. \doi{10.1177/0037549709102150}
#'
#' Sugasawa, S., Murakami, D. (2021) \emph{Spatially clustered regression}.
#' Spatial Statistics, 44, 100525. \doi{10.1016/j.spasta.2021.100525}
#'
#' @examples
#' # Daily PM2.5 at 36 background stations of the Po Valley. The first 180
#' # days are used to keep the example fast; the covariates are the intercept,
#' # the station altitude and the daily PM10 concentration.
#' data(povalley)
#'
#' Tn <- 180L
#' Tfull <- nrow(povalley$z)
#' d <- ncol(povalley$z)
#' # rows of the stacked covariate matrix belonging to the first Tn days
#' keep <- as.vector(outer(seq_len(Tn), (seq_len(d) - 1L) * Tfull, '+'))
#'
#' phi <- list(beta = matrix(c(1.25, -0.00003, 0.64), 3, 1),
#'             sigma2eps = 18.66,
#'             sigma2omega = 1e-06,
#'             theta = 2e-06,
#'             G = matrix(0.59, 1, 1),
#'             Sigmaeta = matrix(4.25, 1, 1),
#'             m0 = as.matrix(0),
#'             C0 = as.matrix(1))
#'
#' mod <- STEM_Model(z = povalley$z[seq_len(Tn), ],
#'                   covariates = povalley$covariates[keep, ],
#'                   coordinates = povalley$coords,
#'                   phi = phi, K = matrix(1, d, 1))
#'
#' \donttest{
#' # three spatial regimes with a moderate spatial penalty
#' fit <- SCSTEM_Estimation(mod, k = 3, phi_penalty = 0.05, distance = 'geo')
#' fit
#'
#' # the estimated regimes on the map
#' plot(povalley$coords, col = fit$group, pch = 19,
#'      xlab = 'Longitude', ylab = 'Latitude')
#' }
#'
#' @seealso \code{\link{STEM_Model}}, \code{\link{STEM_Estimation}},
#'   \code{\link{SCSTEM_Infocrit}}, \code{\link{SCSTEM_Select}},
#'   \code{\link{SCSTEM_Bootstrap}} and \code{\link{pm10}}
#'
#' @keywords models spatial
#'
#' @export
SCSTEM_Estimation <- function(StemModel,
                         k = 3,
                         phi_penalty = 0.05,
                         phi_scale = c("auto", "per-observation", "raw"),
                         knn = 5,
                         distance = c("geo", "euclidean"),
                         init_method = c("departures", "kmeans", "coordinates"),
                         init_partition = NULL,
                         label_update = c("ICM", "simultaneous"),
                         precision = 0.1,
                         precision_full_dataset = 0.01,
                         regularization = 0.01,
                         max_iter = 10,
                         abs_tol = 1e-5,
                         rel_tol = 1e-6,
                         min_cluster_size = NULL,
                         enforce_min_size = TRUE,
                         swap_pass = TRUE,
                         share2conv = 0,
                         seed = 123456789,
                         alpha = 0,
                         lambda = 0,
                         penalize = NULL,
                         lambda_scale = c("relative", "absolute"),
                         lambda_by = c("common", "size"),
                         latent = TRUE,
                         spatial = TRUE,
                         verbose = FALSE) {

  ### the arguments as given, for the unpenalized fit that starts a penalized one
  call_args <- as.list(environment())

  ##############################
  ########## Checks ###########
  ##############################

  if (!inherits(StemModel, "STEM_Model")) {
    stop("'StemModel' must be an object of class 'STEM_Model'.", call. = FALSE)
  }
  phi_scale <- match.arg(phi_scale)
  lambda_scale <- match.arg(lambda_scale)
  lambda_by <- match.arg(lambda_by)

  ### How the penalty is distributed over the regimes.
  ###
  ### "common" gives every regime the same lambda. Under the default relative
  ### scale that already means the same PROPORTIONAL shrinkage everywhere,
  ### because the reference against which lambda is measured is computed inside
  ### each regime and therefore carries its own information: a regime with half
  ### the locations has a proportionally smaller reference, so the same lambda
  ### buys the same fraction of the path. This is the pooled hyperparameter of
  ### the manuscript, and one number to select rather than k.
  ###
  ### "size" departs from that on purpose, giving lambda_g = lambda * nbar / n_g
  ### with nbar = d/k: a regime with half the average number of locations is
  ### shrunk twice as hard. The argument for it is that small regimes carry
  ### noisier coefficients than the proportional rule alone accounts for. It
  ### remains ONE hyperparameter; genuinely cluster-specific (alpha_g, lambda_g)
  ### is a different model and is not offered here.
  lambda_of <- function(n_g) {
    if (lambda <= 0 || lambda_by == "common") return(lambda)
    nbar <- length(StemModel$data$z[1, ]) / k
    lambda * nbar / max(n_g, 1)
  }
  distance <- match.arg(distance)
  init_method <- match.arg(init_method)
  label_update <- match.arg(label_update)

  if (length(k) != 1L || is.na(k) || k < 1 || k != round(k)) {
    stop("'k' must be a single positive integer.", call. = FALSE)
  }
  k <- as.integer(k)
  if (length(phi_penalty) != 1L || is.na(phi_penalty) || phi_penalty < 0) {
    stop("'phi_penalty' must be a single non-negative number.", call. = FALSE)
  }
  if (share2conv < 0 || share2conv >= 1) {
    stop("'share2conv' must lie in [0, 1).", call. = FALSE)
  }

  ##############################
  ########## Setup ############
  ##############################

  z <- StemModel$data$z
  coordinates <- StemModel$data$coordinates
  covariates <- StemModel$data$covariates
  phi0 <- StemModel$skeleton$phi
  Kmat <- StemModel$skeleton$K
  pdim <- StemModel$skeleton$p

  d <- ncol(z)
  Tobs <- nrow(z)
  ncov <- ncol(covariates)
  Nobs <- d * Tobs

  if (is.null(min_cluster_size)) min_cluster_size <- ncov + 2L
  min_cluster_size <- max(2L, as.integer(min_cluster_size))

  if (k > 1 && d < k * min_cluster_size) {
    stop("Too few locations (", d, ") for k = ", k,
         " clusters of at least ", min_cluster_size, " locations each.", call. = FALSE)
  }

  ### Penalty multiplier: see the "Scale of the penalty" paragraph in Details.
  ### Under "auto" the multiplier is data-driven and is computed once, at the
  ### first assignment step, from the spread of the location-wise
  ### log-likelihood contributions across clusters.
  pen_mult <- switch(phi_scale,
                     "per-observation" = Tobs,
                     "raw" = 1,
                     "auto" = NA_real_)
  phi_eff <- phi_penalty * pen_mult

  npar_g <- scstem_npar(ncov = ncov, pdim = pdim)

  if (k == 1L) {
    ### Pooled model: a single STEM fit on the whole network
    if (isTRUE(verbose)) message("Pooled STEM fit (k = 1) ...")
    pooled <- STEM_Estimation(StemModel, precision = precision_full_dataset,
                              distance = distance, regularization = regularization,
                              verbose = FALSE, alpha = alpha, lambda = lambda_of(d),
                              penalize = penalize, lambda_scale = lambda_scale,
                              latent = latent, spatial = spatial)
    par_names <- names(unlist(pooled$estimates$phi.hat))
    loglik <- as.numeric(pooled$estimates$loglik)
    ### effective parameter count under a penalty; see the k > 1 branch
    df_p <- pooled$estimates$penalty$beta.df
    npar_1 <- if (lambda > 0 && !is.null(df_p) && is.finite(df_p))
      npar_g - ncov + df_p else npar_g
    info <- c(loglik = loglik, k = npar_1,
              AIC = -2 * loglik + 2 * npar_1,
              BIC = -2 * loglik + log(Nobs) * npar_1,
              KIC = -2 * loglik + 3 * npar_1)
    out <- list(
      phi_hat = matrix(unlist(pooled$estimates$phi.hat), nrow = 1,
                       dimnames = list("cluster 1", par_names)),
      group = rep(1L, d),
      df = data.frame(coordinates, cluster = 1L),
      fit_list = list(pooled),
      idx_g = list(seq_len(d)),
      info_crit = info,
      loglik_g = loglik,
      final_refit = TRUE,
      obj_trace = data.frame(iter = integer(0), objective = numeric(0),
                             label_changes = integer(0)),
      convergence = "Pooled model (k = 1): no clustering performed",
      penalized_obj = NA_real_,
      input_args = list(StemModel = StemModel, k = 1L, phi_penalty = phi_penalty,
                        phi_scale = phi_scale, knn = knn, distance = distance,
                        init_method = init_method, label_update = label_update,
                        precision = precision,
                        precision_full_dataset = precision_full_dataset,
                        regularization = regularization, max_iter = max_iter,
                        abs_tol = abs_tol, rel_tol = rel_tol,
                        min_cluster_size = min_cluster_size,
                        enforce_min_size = enforce_min_size,
                        swap_pass = swap_pass,
                        share2conv = share2conv, seed = seed,
                        Tobs = Tobs, d = d, ncov = ncov, pdim = pdim,
                        npar_g = npar_g, Nobs = Nobs,
                        alpha = alpha, lambda = lambda, penalize = penalize,
                        lambda_scale = lambda_scale, lambda_by = lambda_by,
                        latent = latent, spatial = spatial)
    )
    class(out) <- c("SCSTEM_Estimation", "list")
    return(out)
  }

  ### Spatial penalty graph (symmetrized knn)
  nbinfo <- scstem_neighbors(coordinates, knn = knn, distance = distance)
  nb <- nbinfo$nb
  W <- nbinfo$W

  ### A penalized fit starts from the solution of the unpenalized fit at the
  ### same k. A penalty that is strong from the first sweep freezes whatever
  ### partition it is given, so a start unrelated to the regimes would stay
  ### where it is; from the unpenalized solution the penalty only has to decide
  ### how much spatial smoothing that solution can afford. The same start also
  ### fixes the scale of the penalty: the first sweep scores the locations with
  ### the parameters of the unpenalized fit, so the factor of phi_scale = "auto"
  ### is the one of the unpenalized fit and does not depend on init_method.
  ### A given init_partition is taken as that start (SCSTEM_Infocrit() passes
  ### the partition of its unpenalized fit).
  if (phi_penalty > 0 && is.null(init_partition)) {
    if (isTRUE(verbose)) message("* Fitting phi = 0 first: its solution is the start ...")
    args0 <- call_args
    args0$phi_penalty <- 0
    args0$verbose <- FALSE
    fit0 <- tryCatch(suppressWarnings(do.call(SCSTEM_Estimation, args0)),
                     error = function(e) NULL)
    if (!is.null(fit0) && length(unique(fit0$group)) == k) {
      init_partition <- fit0$group
    } else {
      warning("SCSTEM_Estimation: the unpenalized fit at k = ", k, " failed or ",
              "lost a cluster; the penalized fit starts from init_method.",
              call. = FALSE)
    }
  }

  ### Initial partition. The seed is applied through scstem_with_seed(), which
  ### restores the RNG stream on exit so that the user's workspace is left
  ### untouched (CRAN policy on .GlobalEnv).
  Xmeans <- scstem_covariate_means(covariates, d = d, Tobs = Tobs)
  if (!is.null(init_partition)) {
    ### user-supplied starting partition (for instance one obtained with an
    ### external clustering routine); it is validated and, if needed, repaired
    ### so that it satisfies the minimum-size requirement
    labels <- as.integer(as.factor(init_partition))
    if (length(labels) != d) {
      stop("'init_partition' must have one label per location (length ", d, ").",
           call. = FALSE)
    }
    if (length(unique(labels)) != k) {
      stop("'init_partition' defines ", length(unique(labels)),
           " groups, but k = ", k, " was requested.", call. = FALSE)
    }
    if (isTRUE(enforce_min_size) &&
        any(tabulate(labels, nbins = k) < min_cluster_size)) {
      labels <- scstem_repair_partition(labels, feat = Xmeans, k = k,
                                        min_size = min_cluster_size)
    }
  } else {
    ### "departures": k-means on how each location departs from the pooled fit
    ### (see scstem_departures()); should the pooled fit fail, the covariate
    ### means are used instead
    feat <- Xmeans
    method0 <- init_method
    if (init_method == "departures") {
      method0 <- "kmeans"
      pooled0 <- tryCatch(
        STEM_Estimation(StemModel, precision = precision_full_dataset,
                        distance = distance, regularization = regularization,
                        verbose = FALSE, alpha = alpha, lambda = lambda,
                        penalize = penalize, lambda_scale = lambda_scale,
                        latent = latent, spatial = spatial),
        error = function(e) NULL)
      if (is.null(pooled0)) {
        warning("SCSTEM_Estimation: the pooled fit of init_method = \"departures\" ",
                "failed; the partition is initialized on the covariate means.",
                call. = FALSE)
      } else {
        feat <- scstem_departures(z, covariates, STEM_Signal(pooled0), Tobs)
      }
    }
    labels <- scstem_with_seed(
      seed,
      scstem_init(Xmeans = feat, coords = coordinates, k = k,
                  method = method0, min_size = min_cluster_size)
    )
  }
  if (length(unique(labels)) < k) {
    stop("The initialization returned fewer than k = ", k,
         " non-empty clusters. Try a smaller k or a different init_method.",
         call. = FALSE)
  }

  ##########################################
  ########## Alternating algorithm #########
  ##########################################

  beta_g <- matrix(NA_real_, nrow = k, ncol = ncov)
  s2eps_g <- s2omega_g <- rep(NA_real_, k)
  ysm_g <- vector("list", k)
  fit <- vector("list", k)
  has_valid <- rep(FALSE, k)
  stale_warned <- FALSE

  obj_prev <- -Inf
  best_obj <- -Inf
  best_labels <- labels
  label_history <- character(0)
  obj_trace <- data.frame(iter = integer(0), objective = numeric(0), swaps = integer(0),
                          label_changes = integer(0), min_cluster = integer(0))
  convergence <- "Maximum number of iterations reached"

  for (it in seq_len(max_iter)) {

    labels_prev <- labels

    ### ---------------------------------------------------------------
    ### Step 1: cluster-wise parameter update, given the labels
    ### ---------------------------------------------------------------
    for (g in seq_len(k)) {
      idx <- which(labels == g)
      if (length(idx) >= min_cluster_size) {
        mod_g <- try(
          STEM_Model(z = z[, idx, drop = FALSE],
                     covariates = covariates[scstem_rows(idx, Tobs), , drop = FALSE],
                     coordinates = coordinates[idx, , drop = FALSE],
                     phi = phi0,
                     K = Kmat[idx, , drop = FALSE]),
          silent = TRUE)
        fit_g <- if (inherits(mod_g, "try-error")) mod_g else try(
          STEM_Estimation(mod_g, precision = precision, distance = distance,
                          regularization = regularization, verbose = FALSE,
                          alpha = alpha, lambda = lambda_of(length(idx)),
                          penalize = penalize, lambda_scale = lambda_scale,
                          latent = latent, spatial = spatial),
          silent = TRUE)

        if (!inherits(fit_g, "try-error") && !is.null(fit_g$estimates$phi.hat)) {
          fit[[g]] <- fit_g
          beta_g[g, ] <- as.numeric(fit_g$estimates$phi.hat$beta)
          s2eps_g[g] <- as.numeric(fit_g$estimates$phi.hat$sigma2eps)
          s2omega_g[g] <- as.numeric(fit_g$estimates$phi.hat$sigma2omega)
          ysm_g[[g]] <- as.matrix(fit_g$estimates$y.smoothed)
          has_valid[g] <- TRUE
        }
      } else if (has_valid[g]) {
        ### stale freeze: keep the last valid parameters
        if (!stale_warned && isTRUE(verbose)) {
          message("* Cluster ", g, " fell below the minimum size (", min_cluster_size,
                  ") at iteration ", it, ": its parameters are kept frozen.")
          stale_warned <- TRUE
        }
      }
    }

    if (!any(has_valid)) {
      stop("No cluster could be estimated: try a smaller k, a larger ",
           "min_cluster_size or looser convergence settings.", call. = FALSE)
    }

    ### ---------------------------------------------------------------
    ### Step 2: label update, given the parameters
    ### ---------------------------------------------------------------
    LL <- matrix(-Inf, nrow = d, ncol = k)
    for (g in seq_len(k)) {
      if (!has_valid[g]) next
      for (i in seq_len(d)) {
        LL[i, g] <- scstem_loglike_i(
          z_i = z[, i],
          X_i = covariates[scstem_rows(i, Tobs), , drop = FALSE],
          beta = beta_g[g, ],
          ysm = ysm_g[[g]],
          K_i = Kmat[i, , drop = FALSE],
          sigma2eps = s2eps_g[g],
          sigma2omega = s2omega_g[g]
        )
      }
    }
    LL[!is.finite(LL)] <- -Inf

    ### Data-driven penalty scale, fixed once at the first sweep so that the
    ### objective stays comparable along the iterations. The multiplier is the
    ### median across locations of the spread of the log-likelihood
    ### contributions across clusters, divided by the average number of
    ### neighbors: phi_penalty = 1 is then the point at which full agreement
    ### with the neighborhood is worth as much as the typical gain from
    ### picking the best-fitting cluster.
    if (phi_scale == "auto" && !is.finite(phi_eff)) {
      rng <- apply(LL, 1, function(r) {
        r <- r[is.finite(r)]
        if (length(r) < 2) NA_real_ else max(r) - min(r)
      })
      mean_nb <- mean(vapply(nb, length, integer(1)))
      pen_mult <- stats::median(rng, na.rm = TRUE) / max(mean_nb, 1)
      if (!is.finite(pen_mult) || pen_mult <= 0) pen_mult <- 1
      phi_eff <- phi_penalty * pen_mult
      if (isTRUE(verbose)) {
        message("* Penalty scale (auto): multiplier = ", signif(pen_mult, 4),
                " ; effective phi = ", signif(phi_eff, 4))
      }
    }

    ### Clusters that never obtained a valid fit cannot receive locations.
    LL[, !has_valid] <- -Inf

    if (label_update == "ICM") {
      ### Sequential (Gauss-Seidel) sweep: each label maximizes its own
      ### penalized contribution given the CURRENT labels of all the others,
      ### already-updated neighbors included.
      ###
      ### The sweep is constrained: a location may leave its cluster only if
      ### that cluster would stay at or above min_cluster_size. Without this
      ### constraint the cluster with the smallest residual variance attracts
      ### every location -- the clusterwise analogue of the degenerate-likelihood
      ### problem of Gaussian mixtures, and of the boundary solutions discussed
      ### and the partition
      ### collapses. Restricting the moves keeps every configuration admissible
      ### and preserves the monotonicity of the sweep within the feasible set.
      sizes <- tabulate(labels, nbins = k)
      for (i in seq_len(d)) {
        nbi <- nb[[i]]
        penvec <- if (length(nbi)) tabulate(labels[nbi], nbins = k) else rep(0, k)
        qd <- LL[i, ] + phi_eff * penvec
        gi <- labels[i]
        if (isTRUE(enforce_min_size) && sizes[gi] <= min_cluster_size) {
          ### the current cluster cannot afford to lose this location
          next
        }
        gnew <- which.max(qd)
        if (gnew != gi) {
          sizes[gi] <- sizes[gi] - 1L
          sizes[gnew] <- sizes[gnew] + 1L
          labels[i] <- gnew
        }
      }
    } else {
      ### Joint update of all labels, with the penalty evaluated at the labels
      ### of the previous iteration (pre-2.0.0 behavior). The joint update
      ### offers no way to impose the size constraint move by move, so an
      ### inadmissible configuration is repaired afterwards.
      Ind <- matrix(0, nrow = d, ncol = k)
      Ind[cbind(seq_len(d), labels_prev)] <- 1
      Pen <- W %*% Ind
      labels <- apply(LL + phi_eff * as.matrix(Pen), 1, which.max)
      if (isTRUE(enforce_min_size) && any(tabulate(labels, nbins = k) < min_cluster_size)) {
        labels <- scstem_repair_partition(labels, feat = Xmeans, k = k,
                                          min_size = min_cluster_size)
      }
    }

    ### ---------------------------------------------------------------
    ### Size-preserving swap pass
    ### ---------------------------------------------------------------
    ### The constrained sweep above cannot move a location out of a cluster
    ### sitting at the minimum size, which with small networks can freeze a
    ### large share of the partition at its initial value. Exchanging the
    ### labels of two locations leaves all cluster sizes unchanged, so it is
    ### always feasible, and it is accepted only when it strictly increases
    ### the penalized objective.
    n_swap <- 0L
    if (isTRUE(swap_pass) && isTRUE(enforce_min_size)) {
      sw <- scstem_swap_pass(labels, LL, phi_eff, nb)
      labels <- sw$labels
      n_swap <- sw$nswap
    }

    ### ---------------------------------------------------------------
    ### Penalized objective at the current (parameters, labels)
    ### ---------------------------------------------------------------
    obj <- sum(LL[cbind(seq_len(d), labels)]) +
      phi_eff * scstem_potts_pairs(labels, nb)
    n_changes <- sum(labels != labels_prev)
    obj_trace <- rbind(obj_trace,
                       data.frame(iter = it, objective = obj,
                                  label_changes = n_changes,
                                  swaps = n_swap,
                                  min_cluster = min(tabulate(labels, nbins = k))))
    if (is.finite(obj) && obj > best_obj) {
      best_obj <- obj
      best_labels <- labels
    }

    if (isTRUE(verbose)) {
      message("* Iteration ", it, ": penalized objective = ", round(obj, 4),
              " ; label changes = ", n_changes,
              " ; smallest cluster = ", min(tabulate(labels, nbins = k)))
    }

    ### ---------------------------------------------------------------
    ### Exit conditions
    ### ---------------------------------------------------------------
    if (n_changes == 0) {
      convergence <- "Convergence reached (labels stable)"
      break
    }
    abs_imp <- abs(obj - obj_prev)
    rel_imp <- abs_imp / (abs(obj_prev) + .Machine$double.eps)
    if (it > 1 && (abs_imp < abs_tol || rel_imp < rel_tol)) {
      convergence <- "Convergence reached (objective stable)"
      break
    }
    if (share2conv > 0 && (n_changes / d) <= share2conv) {
      convergence <- paste0("Convergence reached (fewer than ", share2conv * 100,
                            "% of the locations changed cluster)")
      break
    }
    lab_hash <- paste(labels, collapse = ",")
    if (lab_hash %in% label_history) {
      labels <- best_labels
      convergence <- "Exited on label cycle: best visited partition returned"
      warning("SCSTEM_Estimation: label cycle detected; the best visited partition was returned.",
              call. = FALSE)
      break
    }
    label_history <- c(label_history, lab_hash)
    obj_prev <- obj
  }

  ### Always return the best partition visited.
  ###
  ### The ICM sweep is monotone for FIXED parameters, but the alternation as a
  ### whole is not guaranteed to increase Q. The reason is structural: the
  ### assignment score is the conditional pseudo-likelihood of
  ### scstem_loglike_i(), whereas the parameter step maximizes the EXACT
  ### within-cluster likelihood through the EM algorithm. The two objectives
  ### agree on what a good partition looks like but are not the same function,
  ### so a parameter update can lower Q even while it raises the exact
  ### likelihood. The spatial covariance across locations is what makes an exact
  ### per-location decomposition unavailable in the first place.
  ### Keeping the best visited partition makes the returned solution
  ### well defined regardless of the path taken.
  final_obj <- if (nrow(obj_trace)) obj_trace$objective[nrow(obj_trace)] else -Inf
  if (is.finite(best_obj) && best_obj > final_obj + 1e-8) {
    labels <- best_labels
    convergence <- paste0(convergence, "; best visited partition returned")
  }

  ##########################################
  ########## Final refit ###################
  ##########################################

  loglik_g <- rep(NA_real_, k)
  final_refit <- rep(FALSE, k)
  idx_g <- vector("list", k)
  fit_final <- vector("list", k)
  par_list <- vector("list", k)

  for (g in seq_len(k)) {
    idx_g[[g]] <- which(labels == g)
    if (length(idx_g[[g]]) < min_cluster_size) next
    mod_g <- try(
      STEM_Model(z = z[, idx_g[[g]], drop = FALSE],
                 covariates = covariates[scstem_rows(idx_g[[g]], Tobs), , drop = FALSE],
                 coordinates = coordinates[idx_g[[g]], , drop = FALSE],
                 phi = phi0,
                 K = Kmat[idx_g[[g]], , drop = FALSE]),
      silent = TRUE)
    fit_g <- if (inherits(mod_g, "try-error")) mod_g else try(
      STEM_Estimation(mod_g, precision = precision, distance = distance,
                      regularization = regularization, verbose = FALSE,
                      alpha = alpha, lambda = lambda_of(length(idx_g[[g]])),
                      penalize = penalize, lambda_scale = lambda_scale,
                      latent = latent, spatial = spatial),
      silent = TRUE)
    if (!inherits(fit_g, "try-error") && !is.null(fit_g$estimates$phi.hat)) {
      fit_final[[g]] <- fit_g
      par_list[[g]] <- unlist(fit_g$estimates$phi.hat)
      loglik_g[g] <- as.numeric(fit_g$estimates$loglik)
      final_refit[g] <- TRUE
    }
  }

  if (!any(final_refit)) {
    stop("No cluster could be re-estimated on the final partition. ",
         "Try a smaller k, a larger min_cluster_size, or a stronger phi_penalty.",
         call. = FALSE)
  }

  par_names <- names(par_list[[which(final_refit)[1]]])
  phi_hat <- matrix(NA_real_, nrow = k, ncol = length(par_names),
                    dimnames = list(paste("cluster", seq_len(k)), par_names))
  for (g in which(final_refit)) phi_hat[g, ] <- par_list[[g]][par_names]

  if (any(!final_refit)) {
    warning("SCSTEM_Estimation: cluster(s) ",
            paste(which(!final_refit), collapse = ", "),
            " could not be re-estimated on the final partition (size below ",
            min_cluster_size, " or failed fit). Their estimates are NA.",
            call. = FALSE)
  }

  ### Exact total log-likelihood and information criteria
  ###
  ### With a penalty on the regression coefficients the number of free
  ### parameters is no longer the nominal one: what enters the criteria is the
  ### EFFECTIVE number of coefficients the penalty leaves, which each regime
  ### reports as estimates$penalty$beta.df. For ridge that is
  ### tr(M (M + lambda D)^{-1}), for the lasso the number of active
  ### coefficients, and for the elastic net the corresponding trace on the
  ### active set. With lambda = 0 every regime returns r and this reduces to
  ### k_eff * npar_g, the count the package has always used.
  k_eff <- sum(final_refit)
  loglik_tot <- sum(loglik_g[final_refit])
  k_par <- if (lambda > 0) {
    sum(vapply(which(final_refit), function(g) {
      df_g <- fit_final[[g]]$estimates$penalty$beta.df
      if (is.null(df_g) || !is.finite(df_g)) ncov else df_g
    }, numeric(1))) + k_eff * (npar_g - ncov)
  } else {
    k_eff * npar_g
  }
  info <- c(loglik = loglik_tot, k = k_par,
            AIC = -2 * loglik_tot + 2 * k_par,
            BIC = -2 * loglik_tot + log(Nobs) * k_par,
            KIC = -2 * loglik_tot + 3 * k_par)

  if (isTRUE(verbose)) {
    message("SC-STEM estimation ended (", convergence, "); k_eff = ", k_eff,
            " ; BIC = ", round(info[["BIC"]], 3))
  }

  out <- list(
    phi_hat = phi_hat,
    group = as.integer(labels),
    df = data.frame(coordinates, cluster = as.integer(labels)),
    fit_list = fit_final,
    idx_g = idx_g,
    info_crit = info,
    loglik_g = loglik_g,
    final_refit = final_refit,
    obj_trace = obj_trace,
    convergence = convergence,
    penalized_obj = if (is.finite(best_obj)) best_obj else NA_real_,
    best_objective = if (is.finite(best_obj)) best_obj else NA_real_,
    last_objective = if (nrow(obj_trace)) obj_trace$objective[nrow(obj_trace)] else NA_real_,
    phi_effective = phi_eff,
    phi_multiplier = pen_mult,
    input_args = list(StemModel = StemModel, k = k, phi_penalty = phi_penalty,
                      phi_scale = phi_scale, knn = knn, distance = distance,
                      init_method = init_method, label_update = label_update,
                      precision = precision,
                      precision_full_dataset = precision_full_dataset,
                      regularization = regularization, max_iter = max_iter,
                      abs_tol = abs_tol, rel_tol = rel_tol,
                      min_cluster_size = min_cluster_size,
                      enforce_min_size = enforce_min_size,
                      swap_pass = swap_pass,
                      share2conv = share2conv, seed = seed,
                      Tobs = Tobs, d = d, ncov = ncov, pdim = pdim,
                      npar_g = npar_g, Nobs = Nobs,
                      alpha = alpha, lambda = lambda, penalize = penalize,
                      lambda_scale = lambda_scale, lambda_by = lambda_by,
                      latent = latent, spatial = spatial)
  )
  class(out) <- c("SCSTEM_Estimation", "list")
  out
}


#' Print method for SC-STEM fits
#'
#' @param x an object of class \dQuote{SCSTEM_Estimation}.
#' @param digits integer, number of significant digits. Default is 4.
#' @param ... further arguments, currently ignored.
#'
#' @return \code{x}, invisibly. Called for its side effect of printing a
#'   compact summary of the fit.
#'
#' @export
print.SCSTEM_Estimation <- function(x, digits = 4, ...) {
  cat("Spatially-clustered STEM model\n")
  cat("  clusters requested : ", x$input_args$k, "\n", sep = "")
  cat("  clusters estimated : ", sum(x$final_refit), "\n", sep = "")
  cat("  spatial penalty    : phi = ", x$input_args$phi_penalty,
      " (", x$input_args$phi_scale, ", knn = ", x$input_args$knn, ")\n", sep = "")
  cat("  label update       : ", x$input_args$label_update, "\n", sep = "")
  cat("  convergence        : ", x$convergence, "\n", sep = "")
  cat("  cluster sizes      : ",
      paste(as.integer(table(factor(x$group, levels = seq_len(x$input_args$k)))),
            collapse = ", "), "\n", sep = "")
  cat("\nInformation criteria\n")
  print(round(x$info_crit, digits))
  cat("\nCluster-wise estimates\n")
  print(round(x$phi_hat, digits))
  invisible(x)
}
