#' Fitted values of a STEM model
#'
#' @description
#' \code{STEM_Signal} returns the fitted values of a STEM model in the sense the
#' word usually carries: the conditional mean of the response given the latent
#' path,
#' \deqn{\hat\mu_{ti} = x_{ti}' \hat\beta + K_i \hat y_t ,}
#' with \eqn{\hat y_t} the smoothed state. It is the systematic part of the
#' measurement equation --- the regression surface plus the latent process ---
#' and it excludes the measurement error, which is what makes it an estimate of
#' something rather than a copy of the data.
#'
#' @details
#' \strong{This is not \code{\link{STEM_Complete}}, and the difference matters.}
#' \code{STEM_Complete()} returns \eqn{E[z \mid z_{obs}]}: it fills the gaps of
#' an incomplete record and, wherever the response \emph{was} observed, it
#' returns the observation itself, because the conditional expectation of
#' something you have seen is that thing. On a complete record
#' \code{STEM_Complete()} is therefore the data, and scoring it against anything
#' measures nothing. Use it to impute; use \code{STEM_Signal()} to fit.
#'
#' In state-space language \eqn{K_i \hat y_t} is the signal, and the function is
#' named for it; the regression term is added because in this model the
#' systematic part is the sum of the two.
#'
#' \strong{What to score it against.} Against the response \eqn{z} the residual
#' contains the measurement error, which no model can remove, so an in-sample
#' comparison between two specifications is compressed by a constant that
#' belongs to neither. In a simulation the honest target is the true conditional
#' mean; on real data it is the response at held-out cells, with a blocking
#' scheme that respects both dependencies --- see Otto, Fasso and Maranzano
#' (2024).
#'
#' @param StemModel a fitted model of class \dQuote{STEM_Model}, as returned by
#'   \code{\link{STEM_Estimation}} or \code{\link{STEM_Fit}}.
#'
#' @return A \eqn{T} by \eqn{d} matrix of fitted values, with the dimnames of
#'   the response.
#'
#' @seealso \code{\link{STEM_Complete}} for the imputation of missing responses,
#'   \code{\link{STEM_Kriging}} for prediction at unobserved locations,
#'   \code{\link{SCSTEM_Signal}} for the clustered model.
#'
#' @references
#' Durbin, J., Koopman, S. J. (2012) \emph{Time Series Analysis by State Space
#' Methods}, 2nd edition. Oxford University Press.
#'
#' Otto, P., Fasso, A., Maranzano, P. (2024) \emph{A review of regularised
#' estimation methods and cross-validation in spatiotemporal statistics}.
#' Statistics Surveys, 18, 299--340. \doi{10.1214/24-SS150}
#'
#' @export
STEM_Signal <- function(StemModel) {

  if (is.null(StemModel$estimates) || is.null(StemModel$estimates$phi.hat)) {
    stop("'StemModel' must be a fitted model, as returned by STEM_Estimation().",
         call. = FALSE)
  }

  d   <- StemModel$data$d
  n   <- StemModel$data$n
  r   <- StemModel$data$r
  phi <- StemModel$estimates$phi.hat
  ysm <- as.matrix(StemModel$estimates$y.smoothed)          # n x p
  Kmat <- as.matrix(StemModel$skeleton$K)                   # d x p

  XX   <- changedimension_covariates(StemModel$data$covariates, d, r, n)
  beta <- matrix(as.numeric(phi$beta), ncol = 1)

  out <- matrix(NA_real_, n, d)
  for (tt in seq_len(n)) {
    out[tt, ] <- as.numeric(XX[, , tt] %*% beta +
                            Kmat %*% matrix(ysm[tt, ], ncol = 1))
  }
  dimnames(out) <- dimnames(as.matrix(StemModel$data$z))
  out
}


#' Fitted values of a fitted SC-STEM model
#'
#' @description
#' \code{SCSTEM_Signal} returns the fitted values of a spatially-clustered STEM
#' model: \code{\link{STEM_Signal}} applied within each estimated regime and
#' reassembled onto the original ordering of the locations. Every column is
#' computed from the parameters and the latent path of the regime its location
#' was assigned to.
#'
#' @details
#' See \code{\link{STEM_Signal}} for what the quantity is and for why it is not
#' \code{\link{SCSTEM_Complete}}. Regimes that the final refit could not
#' estimate contribute columns of \code{NA} rather than a value borrowed from
#' elsewhere.
#'
#' @param SCSTEM an object of class \dQuote{SCSTEM_Estimation}, as returned by
#'   \code{\link{SCSTEM_Estimation}} or by \code{\link{STEM_Fit}} with
#'   \code{k > 1}.
#'
#' @return A \eqn{T} by \eqn{d} matrix of fitted values.
#'
#' @seealso \code{\link{STEM_Signal}}, \code{\link{SCSTEM_Complete}},
#'   \code{\link{SCSTEM_Estimation}}
#'
#' @export
SCSTEM_Signal <- function(SCSTEM) {

  if (!inherits(SCSTEM, "SCSTEM_Estimation")) {
    stop("'SCSTEM' must be an object of class 'SCSTEM_Estimation' returned by SCSTEM_Estimation().",
         call. = FALSE)
  }

  d <- length(SCSTEM$group)
  n <- nrow(as.matrix(SCSTEM$fit_list[[which(!vapply(SCSTEM$fit_list, is.null,
                                                     logical(1)))[1]]]$data$z))
  out <- matrix(NA_real_, n, d)

  for (g in seq_along(SCSTEM$fit_list)) {
    fit <- SCSTEM$fit_list[[g]]
    if (is.null(fit) || is.null(fit$estimates$phi.hat)) next
    idx <- if (!is.null(SCSTEM$idx_g)) SCSTEM$idx_g[[g]] else which(SCSTEM$group == g)
    if (!length(idx)) next
    out[, idx] <- STEM_Signal(fit)
  }
  out
}
