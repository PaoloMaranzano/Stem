#' Fit a STEM or SC-STEM model, with or without a penalty on the coefficients
#'
#' @description
#' One entry point for the four estimators the package provides. Which one runs
#' is decided by two arguments and nothing else:
#'
#' \tabular{lll}{
#'   \strong{k} \tab \strong{lambda} \tab \strong{what is fitted} \cr
#'   \code{1}   \tab \code{0}        \tab the pooled STEM model, \code{\link{STEM_Estimation}} \cr
#'   \code{> 1} \tab \code{0}        \tab the spatially-clustered model, \code{\link{SCSTEM_Estim}} \cr
#'   \code{1}   \tab \code{> 0}      \tab the pooled model with an elastic net on \eqn{\beta} \cr
#'   \code{> 1} \tab \code{> 0}      \tab the clustered model with an elastic net on \eqn{\beta} within each regime
#' }
#'
#' The defaults \code{k = 1}, \code{alpha = 0}, \code{lambda = 0} therefore
#' reproduce the classical estimator exactly, to floating point.
#'
#' @details
#' \strong{The penalty.} In the parameterization of \code{glmnet}, the objective
#' subtracts from the expected complete-data log-likelihood
#' \deqn{\lambda \left\{ \alpha \|D\beta\|_1 +
#'       \frac{1-\alpha}{2}\, \beta' D \beta \right\},}
#' with \eqn{D} the diagonal indicator of the penalized coordinates. So
#' \code{alpha = 0} is ridge, \code{alpha = 1} is the lasso, and anything in
#' between is the elastic net. Only the regression coefficients are penalized:
#' the variance components, the range, the transition matrix and the initial
#' state are estimated by maximum likelihood as before.
#'
#' \strong{Only the M-step changes.} The penalty is a function of \eqn{\beta}
#' alone and does not involve the latent states, so the E-step is untouched and
#' the usual argument still gives an EM algorithm that increases the
#' \emph{penalized} observed-data log-likelihood at every iteration. For ridge
#' the M-step keeps its closed form,
#' \eqn{\hat\beta = (M + \lambda(1-\alpha)D)^{-1}v}; for the lasso and the
#' elastic net it is solved by cyclic coordinate descent, which converges to the
#' exact maximizer because the objective is convex and the penalty separable.
#'
#' \strong{The intercept is not penalized} unless \code{penalize} says
#' otherwise, and not merely by convention: the model already carries a latent
#' process whose initial mean absorbs the level, so the intercept and \eqn{m_0}
#' are only weakly identified apart.
#'
#' \strong{Scale.} The accumulation \eqn{M = \sum_t X_t' \Sigma_e^{-1} X_t} is a
#' generalized least squares cross-product, so the penalty acts in that metric
#' and not in the Euclidean one. The columns of the design are therefore scaled
#' internally to put the diagonal of \eqn{M} at one, which is what makes a given
#' \code{lambda} mean the same thing for every covariate and, in the clustered
#' model, for every regime. Coefficients are returned on the original scale.
#'
#' \strong{Choosing lambda and alpha.} They are hyperparameters like \code{k}
#' and \code{phi_penalty} and are chosen the same way: over a grid, by an
#' information criterion computed with the \emph{effective} number of
#' coefficients the penalty leaves, or by spatio-temporal cross-validation.
#' \code{\link{SCSTEM_Infocrit}} accepts \code{alpha} and \code{lambda} through
#' its \code{...} and returns criteria already corrected for the effective
#' degrees of freedom. Inference after selection is not valid naively: if
#' \code{lambda} is chosen from the data, the parametric bootstrap has to repeat
#' the selection on every draw, exactly as \code{\link{SCSTEM_Bootstrap}}
#' repeats the clustering.
#'
#' The derivation -- the penalized EM and its monotonicity, the M-step in both
#' forms, the generalized least squares metric, the effective degrees of freedom
#' and what is still open -- is kept as a standalone document in the repository,
#' at \code{dev/regularization/stem-elastic-net.tex}, since the material is a
#' study of its own rather than documentation of the software.
#'
#' @param StemModel an object of class \dQuote{STEM_Model}, from
#'   \code{\link{STEM_Model}}.
#' @param k number of spatial regimes. \code{1}, the default, fits the pooled
#'   model; \code{k > 1} fits the spatially-clustered model.
#' @param alpha the elastic-net mixing parameter, in \eqn{[0,1]}. \code{0} is
#'   ridge, \code{1} is the lasso, in between is the elastic net. Ignored when
#'   \code{lambda} is zero. Default is 0.
#' @param lambda the strength of the penalty on the regression coefficients.
#'   \code{0}, the default, gives the unpenalized maximum likelihood estimator.
#' @param penalize which coefficients the penalty acts on: a logical vector of
#'   length \eqn{r}, an index vector, or \code{NULL} (the default) for every
#'   coefficient except the intercept.
#' @param phi_penalty the strength of the Potts penalty on the partition, passed
#'   to \code{\link{SCSTEM_Estim}}. Ignored when \code{k = 1}.
#' @param distance \dQuote{geo} or \dQuote{euclidean}. Default is
#'   \dQuote{euclidean} for the pooled model and \dQuote{geo} for the clustered
#'   one, which are the defaults of the two functions being called.
#' @param verbose logical, passed on.
#' @param ... further arguments passed to \code{\link{STEM_Estimation}} when
#'   \code{k = 1} and to \code{\link{SCSTEM_Estim}} when \code{k > 1}.
#'
#' @return The object the underlying function returns: of class
#'   \dQuote{STEM_Model} when \code{k = 1}, of class \dQuote{SCSTEM_Estim} when
#'   \code{k > 1}. In both cases the penalty in force and the effective number
#'   of coefficients are recorded, under \code{estimates$penalty} and inside
#'   each regime's fit respectively.
#'
#' @seealso \code{\link{STEM_Estimation}}, \code{\link{SCSTEM_Estim}},
#'   \code{\link{SCSTEM_Infocrit}}, \code{\link{SCSTEM_Select}}
#'
#' @references
#' Zou, H., Hastie, T. (2005) \emph{Regularization and variable selection via
#' the elastic net}. Journal of the Royal Statistical Society B, 67, 301--320.
#' \doi{10.1111/j.1467-9868.2005.00503.x}
#'
#' Zou, H., Hastie, T., Tibshirani, R. (2007) \emph{On the degrees of freedom of
#' the lasso}. The Annals of Statistics, 35, 2173--2192.
#' \doi{10.1214/009053607000000127}
#'
#' Otto, P., Fasso, A., Maranzano, P. (2024) \emph{A review of regularised
#' estimation methods and cross-validation in spatiotemporal statistics}.
#' Statistics Surveys, 18, 299--340. \doi{10.1214/24-SS150}
#'
#' @examples
#' \donttest{
#' data(povalley)
#' Tn <- 120L
#' d  <- nrow(povalley$coords)
#' z  <- povalley$z[seq_len(Tn), ]
#' X  <- cbind(1, as.vector(povalley$altitude_std[seq_len(Tn), ]))
#' ols <- stats::lm.fit(x = X, y = as.vector(z))
#' s2  <- stats::var(ols$residuals)
#' mod <- STEM_Model(z = z, covariates = X, coordinates = povalley$coords,
#'   phi = list(beta = matrix(ols$coefficients, ncol = 1),
#'              sigma2eps = 0.7 * s2, sigma2omega = 0.3 * s2,
#'              theta = 1 / 100000, G = matrix(0.8, 1, 1),
#'              Sigmaeta = matrix(0.2 * s2, 1, 1),
#'              m0 = as.matrix(0), C0 = as.matrix(1)),
#'   K = matrix(1, d, 1))
#'
#' ## the classical estimator
#' fit0 <- STEM_Fit(mod, distance = "geo", max.iter = 5)
#'
#' ## the same model with a ridge on the coefficients
#' fitr <- STEM_Fit(mod, alpha = 0, lambda = 1, distance = "geo", max.iter = 5)
#' }
#'
#' @export
STEM_Fit <- function(StemModel, k = 1, alpha = 0, lambda = 0, penalize = NULL,
                     phi_penalty = 1, distance = NULL, verbose = FALSE, ...) {

  if (!inherits(StemModel, "STEM_Model")) {
    stop("'StemModel' must be an object of class 'STEM_Model'.", call. = FALSE)
  }
  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || k < 1L) {
    stop("'k' must be a single integer of at least one.", call. = FALSE)
  }
  if (length(alpha) != 1L || is.na(alpha) || alpha < 0 || alpha > 1) {
    stop("'alpha' must be a single number in [0, 1].", call. = FALSE)
  }
  if (length(lambda) != 1L || is.na(lambda) || lambda < 0) {
    stop("'lambda' must be a single non-negative number.", call. = FALSE)
  }

  if (k == 1L) {
    args <- list(StemModel = StemModel, alpha = alpha, lambda = lambda,
                 penalize = penalize, verbose = verbose, ...)
    if (!is.null(distance)) args$distance <- distance
    return(do.call(STEM_Estimation, args))
  }

  args <- list(StemModel = StemModel, k = k, phi_penalty = phi_penalty,
               alpha = alpha, lambda = lambda, penalize = penalize,
               verbose = verbose, ...)
  if (!is.null(distance)) args$distance <- distance
  do.call(SCSTEM_Estim, args)
}
