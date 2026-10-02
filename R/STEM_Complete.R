#' Completed response of a fitted STEM model
#'
#' @description
#' \code{STEM_Complete} returns the conditional expectation of the response given
#' the observed data, \eqn{\mathrm{E}[z \mid \mathcal{Z}]}, for a fitted STEM
#' model. It equals the observed value wherever the response was observed and
#' the model-based prediction wherever it was not, so it is at once the vector
#' of fitted values and the interpolator of the missing entries.
#'
#' @details
#' The prediction at a missing entry is \emph{not} the signal
#' \eqn{x_{it}'\hat\beta + A_i \hat{y}_t} alone. That would be exact only if the
#' measurement covariance were diagonal. Under the STEM specification
#' \eqn{\Sigma_e = \sigma^2_\epsilon I + \sigma^2_\omega C_\theta} couples the
#' locations, so the missing block of the measurement error is predicted from
#' the block observed at the same time point by Gaussian conditioning:
#' \deqn{\hat{z}_{M,t} = x_{M,t}'\hat\beta + A_M \hat{y}_t
#'       + P_t\left(z_{O,t} - x_{O,t}'\hat\beta - A_O \hat{y}_t\right),
#'       \qquad
#'       P_t = \left[\Sigma_e\right]_{MO}\left(\left[\Sigma_e\right]_{OO}\right)^{-1},}
#' which is the same algebra as kriging, applied at a fixed time point. This is
#' the quantity the EM algorithm uses to complete its sufficient statistics; see
#' Durbin and Koopman (2012, Sections 2.7 and 4.10) for the state-space
#' treatment of missing observations.
#'
#' At a time point where nothing was observed the correction vanishes and the
#' prediction is the signal, with the smoothed state carrying the information of
#' the neighboring time points.
#'
#' The function is the natural building block of a cross-validation exercise:
#' blanking a set of cells, refitting, and reading this matrix at those cells
#' gives out-of-sample predictions that never see the blanked values.
#'
#' @param StemModel an object returned by \code{\link{STEM_Estimation}}.
#' @param distance character, \code{"euclidean"} or \code{"geo"}, the distance
#'   used when the model was estimated. Default is \code{"euclidean"}.
#' @param cov.spat the spatial covariance function used when the model was
#'   estimated. Default is the exponential one.
#'
#' @return A \eqn{T} by \eqn{d} numeric matrix.
#'
#' @author Paolo Maranzano \email{pmaranzano.ricercastatistica@gmail.com}
#'
#' @references
#' Durbin, J., Koopman, S.J. (2012) \emph{Time Series Analysis by State Space
#' Methods}, 2nd edition. Oxford University Press, Oxford.
#'
#' @seealso \code{\link{STEM_Estimation}}, \code{\link{SCSTEM_Complete}}
#'
#' @examples
#' data(pm10)
#' phi <- list(beta = matrix(c(3.65, 0.046, -0.904), 3, 1),
#'             sigma2eps = 0.1, sigma2omega = 0.2, theta = 0.01,
#'             G = matrix(0.77, 1, 1), Sigmaeta = matrix(0.3, 1, 1),
#'             m0 = as.matrix(0), C0 = as.matrix(1))
#' z <- pm10$z
#' z[5, 2] <- NA                     # blank one cell
#' mod <- STEM_Model(z = z, covariates = pm10$covariates,
#'                   coordinates = pm10$coords, phi = phi,
#'                   A = matrix(1, ncol(z), 1))
#' \donttest{
#' fit <- STEM_Estimation(mod, precision = 0.05, max.iter = 5)
#' zhat <- STEM_Complete(fit)
#' zhat[5, 2]                        # the model-based prediction of the blank
#' }
#'
#' @keywords models spatial
#'
#' @export
STEM_Complete <- function(StemModel, distance = "euclidean",
                        cov.spat = Sigmastar.exp) {

  if (is.null(StemModel$estimates) || is.null(StemModel$estimates$phi.hat)) {
    stop("'StemModel' must be a fitted model, as returned by STEM_Estimation().",
         call. = FALSE)
  }

  z   <- as.matrix(StemModel$data$z)
  d   <- StemModel$data$d
  n   <- StemModel$data$n
  r   <- StemModel$data$r
  phi <- StemModel$estimates$phi.hat
  ysm <- as.matrix(StemModel$estimates$y.smoothed)
  Amat <- as.matrix(StemModel$skeleton$A)                    # d x p

  XX <- changedimension_covariates(StemModel$data$covariates, d, r, n)

  dm <- if (distance == "geo") {
    as.matrix(geodist::geodist(StemModel$data$coordinates, measure = "geodesic"))
  } else {
    as.matrix(stats::dist(StemModel$data$coordinates, diag = TRUE))
  }

  logb     <- log(phi$sigma2eps / phi$sigma2omega)
  logtheta <- log(phi$theta)
  Sigma_e  <- phi$sigma2omega *
    cov.spat(d = d, logb = logb, logtheta = logtheta, dist = dm)

  obs_ix <- stem_obs_index(z)
  blocks <- stem_blocks_cache(Sigma_e, d, regularization = 0)

  beta <- matrix(as.numeric(phi$beta), ncol = 1)
  out  <- matrix(NA_real_, n, d)

  for (tt in seq_len(n)) {
    signal <- XX[, , tt] %*% beta + Amat %*% matrix(ysm[tt, ], ncol = 1)
    oi <- obs_ix$idx[[tt]]

    if (obs_ix$complete[tt]) {
      out[tt, ] <- z[tt, ]
    } else if (length(oi) == 0L) {
      out[tt, ] <- as.numeric(signal)
    } else {
      bl <- blocks(obs_ix$key[tt], oi)
      r_o <- matrix(z[tt, oi], ncol = 1) - signal[oi, , drop = FALSE]
      out[tt, ] <- as.numeric(signal) + as.numeric(bl$S %*% r_o)
    }
  }
  dimnames(out) <- dimnames(z)
  out
}


#' Completed response of a fitted SC-STEM model
#'
#' @description
#' \code{SCSTEM_Complete} returns the conditional expectation of the response
#' given the observed data for a spatially-clustered STEM fit, by applying
#' \code{\link{STEM_Complete}} within each estimated regime and reassembling the
#' columns in their original order.
#'
#' @details
#' Because the regimes are fitted independently on disjoint sets of locations,
#' the prediction of a missing entry borrows only from the locations of its own
#' regime. This is a property of the model rather than of the implementation:
#' the measurement covariance \eqn{\Sigma_{e,k}} is defined within a regime, and
#' two locations in different regimes are uncorrelated by construction.
#'
#' Locations belonging to a regime that could not be re-estimated are returned
#' unchanged, that is with their missing entries still missing.
#'
#' @param SCSTEM an object of class \dQuote{SCSTEM_Estimation} returned by
#'   \code{\link{SCSTEM_Estimation}}.
#'
#' @return A \eqn{T} by \eqn{d} numeric matrix.
#'
#' @author Paolo Maranzano \email{pmaranzano.ricercastatistica@gmail.com}
#'
#' @seealso \code{\link{STEM_Complete}}, \code{\link{SCSTEM_Estimation}}
#'
#' @keywords models spatial
#'
#' @export
SCSTEM_Complete <- function(SCSTEM) {

  if (!inherits(SCSTEM, "SCSTEM_Estimation")) {
    stop("'SCSTEM' must be an object of class 'SCSTEM_Estimation' returned by SCSTEM_Estimation().",
         call. = FALSE)
  }

  dist <- SCSTEM$input_args$distance
  if (is.null(dist)) dist <- "euclidean"

  d <- length(SCSTEM$group)
  n <- nrow(as.matrix(SCSTEM$fit_list[[which(SCSTEM$final_refit)[1]]]$data$z))
  out <- matrix(NA_real_, n, d)

  for (g in which(SCSTEM$final_refit)) {
    idx <- SCSTEM$idx_g[[g]]
    if (!length(idx)) next
    out[, idx] <- STEM_Complete(SCSTEM$fit_list[[g]], distance = dist)
  }
  out
}
