#' Blocked cross-validation of the predictive accuracy of STEM and SC-STEM fits
#'
#' @description
#' \code{SCSTEM_CV} measures how well one or more fitted models predict data
#' they have not seen, under the form of spatio-temporal extrapolation an
#' application needs, and ranks them. It is a validation tool, optional and
#' separate from the choice of the hyperparameters, which
#' \code{\link{SCSTEM_Select}} makes: it can be applied at any time to the
#' selected fit, to the pooled model, or to models with different covariates,
#' numbers of regimes or penalties, all evaluated on the same folds.
#'
#' @details
#' Each fold removes a set of cells from the data, refits every model on what
#' is left with its own settings, and predicts the removed cells. Four blocking
#' schemes are available (Otto, Fasso and Maranzano, 2024):
#' \describe{
#'   \item{\code{"LKLO"}}{leave locations out: whole locations are removed.
#'     This is prediction at unmonitored sites, and the scheme on which a
#'     partition is tested: a removed location has no label of its own and
#'     inherits the regime of its nearest retained location.}
#'   \item{\code{"LKTO"}}{leave times out: contiguous blocks of time points are
#'     removed at every location, which is gap filling and reconstruction in
#'     time.}
#'   \item{\code{"LKLHTO"}}{leave locations and times out: a band of locations
#'     and a band of time points are both removed from the fit, and only their
#'     intersection is scored, so that neither the location nor the period of a
#'     scored cell contributes anything to the fit. The strictest scheme.}
#'   \item{\code{"random"}}{cells removed at random. Every removed cell keeps its
#'     neighbors in space and time, so this scheme measures interpolation and is
#'     the most optimistic.}
#' }
#'
#' A removed cell of a retained location is predicted by
#' \code{\link{SCSTEM_Complete}}, the conditional mean given the observed data
#' of its regime. A removed location is predicted by the same conditional mean,
#' computed with the parameters of the regime it inherits and given the
#' retained locations of that regime observed at the same time point:
#' \deqn{\hat{z}_{0t} = x_{0t}'\hat\beta_k + K_0 \hat{y}_{k,t}
#'   + \Sigma_{0O}\Sigma_{OO}^{-1}(z_{Ot} - x_{Ot}'\hat\beta_k - K_O\hat{y}_{k,t}),}
#' which is kriging within the regime. For the pooled model, fitted by
#' \code{SCSTEM_Estimation(..., k = 1)}, this is ordinary kriging.
#'
#' The refits reproduce the estimation of each model, with the settings stored
#' in its \code{input_args} (number of regimes, penalty, graph, initialization,
#' penalty on the coefficients), overridable through \code{...}. The models
#' must share the response, the locations and the time points; their
#' covariates may differ. The folds depend only on \code{seed}, the scheme and
#' the dimensions of the data, so separate calls with the same \code{seed}
#' use the same folds.
#'
#' @param SCSTEM an object of class \dQuote{SCSTEM_Estimation}, or a named list
#'   of such objects to be compared on the same folds.
#' @param scheme character vector, one or more of \code{"LKLO"} (default),
#'   \code{"LKTO"}, \code{"LKLHTO"} and \code{"random"}.
#' @param folds integer, the number of folds of each scheme. Default is 5.
#' @param seed optional integer seed for the allocation of locations and cells
#'   to the folds. The random number stream of the session is restored on exit.
#' @param verbose logical. If \code{TRUE}, the progress is reported via
#'   \code{message()}. Default is \code{FALSE}.
#' @param ... arguments of \code{\link{SCSTEM_Estimation}} overriding those of
#'   the original fits in every refit, for instance \code{max_iter}.
#'
#' @return An object of class \dQuote{SCSTEM_CV}, a list with
#' \itemize{
#'   \item \code{summary}: data frame with, for every scheme and model, the
#'     root mean squared error and the mean absolute error over all the scored
#'     cells, the number of scored cells, the number of folds whose refit
#'     succeeded, and the rank of the model by root mean squared error within
#'     the scheme.
#'   \item \code{folds}: data frame with the same errors fold by fold, and the
#'     error message of the folds whose refit failed.
#'   \item \code{models}: data frame with the number of regimes and the penalty
#'     of every model.
#'   \item \code{scheme}, \code{n_folds}, \code{seed}: the settings used.
#' }
#'
#' @author Paolo Maranzano \email{pmaranzano.ricercastatistica@gmail.com}
#'
#' @references
#' Otto, P., Fasso, A., Maranzano, P. (2024) \emph{A review of regularised
#' estimation methods and cross-validation in spatiotemporal statistics}.
#' Statistics Surveys, 18, 299--340. \doi{10.1214/24-SS150}
#'
#' @seealso \code{\link{SCSTEM_Select}}, \code{\link{SCSTEM_Complete}},
#'   \code{\link{STEM_Complete}}
#'
#' @examples
#' data(povalley)
#'
#' Tn <- 180L
#' Tfull <- nrow(povalley$z)
#' d <- ncol(povalley$z)
#' keep <- as.vector(outer(seq_len(Tn), (seq_len(d) - 1L) * Tfull, '+'))
#'
#' phi <- list(beta = matrix(c(1.25, -0.00003, 0.64), 3, 1),
#'             sigma2eps = 18.66, sigma2omega = 1e-06, theta = 2e-06,
#'             G = matrix(0.59, 1, 1), Sigmaeta = matrix(4.25, 1, 1),
#'             m0 = as.matrix(0), C0 = as.matrix(1))
#'
#' mod <- STEM_Model(z = povalley$z[seq_len(Tn), ],
#'                   covariates = povalley$covariates[keep, ],
#'                   coordinates = povalley$coords,
#'                   phi = phi, K = matrix(1, d, 1))
#'
#' \donttest{
#' pooled <- SCSTEM_Estimation(mod, k = 1, distance = 'geo')
#' fit2 <- SCSTEM_Estimation(mod, k = 2, phi_penalty = 0.05, distance = 'geo')
#' cv <- SCSTEM_CV(list(pooled = pooled, k2 = fit2), scheme = "LKLO",
#'                 folds = 4, seed = 1)
#' cv
#' }
#'
#' @keywords models spatial
#'
#' @export
SCSTEM_CV <- function(SCSTEM, scheme = "LKLO", folds = 5L, seed = NULL,
                      verbose = FALSE, ...) {

  fits <- if (inherits(SCSTEM, "SCSTEM_Estimation")) list(model = SCSTEM) else SCSTEM
  if (!is.list(fits) || !length(fits) ||
      !all(vapply(fits, inherits, logical(1), what = "SCSTEM_Estimation"))) {
    stop("'SCSTEM' must be an object of class 'SCSTEM_Estimation' returned by ",
         "SCSTEM_Estimation(), or a list of such objects.", call. = FALSE)
  }
  if (is.null(names(fits)) || any(!nzchar(names(fits)))) {
    names(fits) <- paste0("model", seq_along(fits))
  }
  if (anyDuplicated(names(fits))) stop("The names of the models must be distinct.", call. = FALSE)

  schemes <- c("LKLO", "LKTO", "LKLHTO", "random")
  if (!length(scheme) || !all(scheme %in% schemes)) {
    stop("'scheme' must contain one or more of ",
         paste(dQuote(schemes, FALSE), collapse = ", "), ".", call. = FALSE)
  }
  scheme <- unique(scheme)

  ### The models must share the response, the locations and the time points,
  ### or the folds would not be the same for all of them
  base1 <- fits[[1]]$input_args$StemModel
  z0 <- as.matrix(base1$data$z)
  coords0 <- as.matrix(base1$data$coordinates)
  for (m in fits[-1]) {
    b <- m$input_args$StemModel
    if (!identical(dim(as.matrix(b$data$z)), dim(z0)) ||
        !isTRUE(all.equal(as.matrix(b$data$z), z0, check.attributes = FALSE)) ||
        !isTRUE(all.equal(as.matrix(b$data$coordinates), coords0, check.attributes = FALSE))) {
      stop("The models must share the response, the locations and the time points; ",
           "only their covariates and settings may differ.", call. = FALSE)
    }
  }
  Tobs <- nrow(z0)
  d <- ncol(z0)

  if (length(folds) != 1L || is.na(folds) || folds < 2 || folds != round(folds)) {
    stop("'folds' must be a single integer of at least 2.", call. = FALSE)
  }
  folds <- as.integer(folds)
  if (any(scheme %in% c("LKLO", "LKLHTO")) && folds > d) {
    stop("'folds' cannot exceed the number of locations.", call. = FALSE)
  }
  if (any(scheme %in% c("LKTO", "LKLHTO")) && folds > Tobs) {
    stop("'folds' cannot exceed the number of time points.", call. = FALSE)
  }

  ### Folds: `blank` holds the cells removed from the fit, `score` those on
  ### which the prediction is evaluated, `drop` the locations removed whole.
  ### Each scheme draws its folds under its own seed, before any fit, so that
  ### they depend on nothing else.
  make_folds <- function(sch) {
    tb <- cut(seq_len(Tobs), breaks = folds, labels = FALSE)
    lg <- sample(rep_len(seq_len(folds), d))
    cells <- matrix(sample(rep_len(seq_len(folds), Tobs * d)), Tobs, d)
    lapply(seq_len(folds), function(j) {
      m <- matrix(FALSE, Tobs, d)
      switch(sch,
        LKLO = {
          m[, lg == j] <- TRUE
          list(blank = m, score = m, drop = which(lg == j))
        },
        LKTO = {
          m[tb == j, ] <- TRUE
          list(blank = m, score = m, drop = integer(0))
        },
        LKLHTO = {
          m[, lg == j] <- TRUE
          m[tb == j, ] <- TRUE
          sc <- matrix(FALSE, Tobs, d)
          sc[tb == j, lg == j] <- TRUE
          list(blank = m, score = sc, drop = which(lg == j))
        },
        random = {
          m <- cells == j
          list(blank = m, score = m, drop = integer(0))
        })
    })
  }
  fold_sets <- lapply(scheme, function(sch) {
    scstem_with_seed(if (is.null(seed)) NULL else seed + match(sch, schemes),
                     make_folds(sch))
  })
  names(fold_sets) <- scheme

  ### Everything a model needs to be refitted on a fold and to predict it
  prepare <- function(fit) {
    args <- fit$input_args
    base <- args$StemModel
    ncov <- ncol(base$data$covariates)
    refit_args <- args[setdiff(names(args), "StemModel")]
    ### a fit that carries `control` passes its settings through it only, so
    ### that a control given in `...` is not overridden by the historical fields
    if (!is.null(args$control)) {
      refit_args[c("precision", "precision_full_dataset", "max_iter", "abs_tol", "rel_tol")] <- NULL
    }
    refit_args$verbose <- FALSE
    user_args <- list(...)
    refit_args[names(user_args)] <- user_args
    refit_args <- refit_args[names(refit_args) %in% names(formals(SCSTEM_Estimation))]
    refit_args <- refit_args[!vapply(refit_args, is.null, logical(1))]
    dist_type <- if (is.null(args$distance)) "euclidean" else args$distance
    dm <- if (dist_type == "geo") {
      as.matrix(geodist::geodist(base$data$coordinates, measure = "geodesic"))
    } else {
      as.matrix(stats::dist(base$data$coordinates, diag = TRUE))
    }
    list(base = base, refit_args = refit_args, dm = dm,
         XX = changedimension_covariates(base$data$covariates, d, ncov, Tobs),
         K = as.matrix(base$skeleton$K))
  }

  ### One fold, one model: refit on the retained cells and locations, predict
  ### the removed ones, return the errors on the scored cells
  run_fold <- function(fo, pr) {
    z <- z0
    z[fo$blank] <- NA_real_
    keep <- setdiff(seq_len(d), fo$drop)
    zk <- z[, keep, drop = FALSE]
    if (any(colSums(!is.na(zk)) == 0L)) {
      stop("a retained location has no observation left in this fold")
    }
    rows <- as.vector(outer(seq_len(Tobs), (keep - 1L) * Tobs, "+"))
    mod <- STEM_Model(z = zk,
                      covariates = pr$base$data$covariates[rows, , drop = FALSE],
                      coordinates = pr$base$data$coordinates[keep, , drop = FALSE],
                      phi = pr$base$skeleton$phi,
                      K = pr$K[keep, , drop = FALSE])
    fit <- suppressWarnings(do.call(SCSTEM_Estimation, c(list(StemModel = mod), pr$refit_args)))
    zhat <- matrix(NA_real_, Tobs, d)
    zhat[, keep] <- SCSTEM_Complete(fit)
    if (length(fo$drop)) {
      zhat[, fo$drop] <- scstem_predict_new(fit, keep, fo$drop, zk, pr$dm, pr$XX, pr$K)
    }
    sc <- fo$score & !is.na(z0)
    err <- (zhat - z0)[sc]
    err <- err[is.finite(err)]
    c(sse = sum(err^2), sae = sum(abs(err)), n = length(err))
  }

  preps <- lapply(fits, prepare)
  rows <- list()
  for (sch in scheme) {
    for (j in seq_len(folds)) {
      for (mo in names(fits)) {
        if (isTRUE(verbose)) message("* ", sch, ", fold ", j, " / ", folds, ", ", mo)
        e <- tryCatch(run_fold(fold_sets[[sch]][[j]], preps[[mo]]),
                      error = function(err) conditionMessage(err))
        ok <- is.numeric(e) && e[["n"]] > 0
        rows[[length(rows) + 1L]] <- data.frame(
          scheme = sch, fold = j, model = mo,
          rmse = if (ok) sqrt(e[["sse"]] / e[["n"]]) else NA_real_,
          mae = if (ok) e[["sae"]] / e[["n"]] else NA_real_,
          n = if (ok) e[["n"]] else NA_real_,
          sse = if (ok) e[["sse"]] else NA_real_,
          sae = if (ok) e[["sae"]] else NA_real_,
          message = if (is.character(e)) e else NA_character_,
          stringsAsFactors = FALSE)
      }
    }
  }
  res <- do.call(rbind, rows)

  ### Errors over all the scored cells of a scheme, and the rank of each model
  summary <- do.call(rbind, lapply(scheme, function(sch) {
    s <- do.call(rbind, lapply(names(fits), function(mo) {
      x <- res[res$scheme == sch & res$model == mo & is.finite(res$rmse), , drop = FALSE]
      data.frame(scheme = sch, model = mo,
                 rmse = if (nrow(x)) sqrt(sum(x$sse) / sum(x$n)) else NA_real_,
                 mae = if (nrow(x)) sum(x$sae) / sum(x$n) else NA_real_,
                 n = sum(x$n), folds_ok = nrow(x), stringsAsFactors = FALSE)
    }))
    s$rank <- rank(s$rmse, na.last = "keep", ties.method = "min")
    s[order(s$rank), , drop = FALSE]
  }))
  rownames(summary) <- NULL

  models <- data.frame(model = names(fits),
                       k = vapply(fits, function(f) as.integer(f$input_args$k), integer(1)),
                       phi = vapply(fits, function(f) if (f$input_args$k == 1L) NA_real_ else
                         as.numeric(f$input_args$phi_penalty), numeric(1)),
                       stringsAsFactors = FALSE)
  rownames(models) <- NULL

  out <- list(summary = summary,
              folds = res[, c("scheme", "fold", "model", "rmse", "mae", "n", "message")],
              models = models, scheme = scheme, n_folds = folds, seed = seed)
  class(out) <- c("SCSTEM_CV", "list")
  out
}


### ---------------------------------------------------------------------------
### Prediction at locations removed from a fit
### ---------------------------------------------------------------------------
### Each new location inherits the regime of its nearest retained location and
### is predicted by the conditional mean given the retained locations of that
### regime observed at the same time point: STEM_Complete() on the regime
### augmented by the new locations, whose responses are all missing.
###
### Arguments
###   fit      SCSTEM_Estimation object fitted on the locations `keep`
###   keep     indices of the retained locations in the whole network
###   new      indices of the new locations in the whole network
###   z_train  T x length(keep) response the fit was estimated on
###   dm       distance matrix of the whole network
###   XX       d x r x T covariate array of the whole network
###   K        d x p loading matrix of the whole network
`scstem_predict_new` <- function(fit, keep, new, z_train, dm, XX, K) {
  Tobs <- nrow(z_train)
  pred <- matrix(NA_real_, Tobs, length(new))
  nearest <- keep[apply(dm[new, keep, drop = FALSE], 1L, which.min)]
  reg <- fit$group[match(nearest, keep)]
  for (g in intersect(unique(reg), which(fit$final_refit))) {
    f <- fit$fit_list[[g]]
    ph <- f$estimates$phi.hat
    ysm <- as.matrix(f$estimates$y.smoothed)
    beta <- matrix(as.numeric(ph$beta), ncol = 1)
    members <- keep[fit$group == g]
    add <- new[reg == g]
    A <- c(members, add)
    Sigma <- ph$sigma2omega *
      Sigmastar.exp(length(A), log(ph$sigma2eps / ph$sigma2omega),
                    log(ph$theta), dm[A, A, drop = FALSE])
    zA <- cbind(z_train[, match(members, keep), drop = FALSE],
                matrix(NA_real_, Tobs, length(add)))
    obs <- stem_obs_index(zA)
    blocks <- stem_blocks_cache(Sigma, length(A), regularization = 0)
    at_new <- length(members) + seq_along(add)
    out <- matrix(NA_real_, Tobs, length(add))
    for (tt in seq_len(Tobs)) {
      signal <- matrix(XX[A, , tt], nrow = length(A)) %*% beta +
        K[A, , drop = FALSE] %*% matrix(ysm[tt, ], ncol = 1)
      oi <- obs$idx[[tt]]
      if (!length(oi)) {
        out[tt, ] <- signal[at_new]
      } else {
        bl <- blocks(obs$key[tt], oi)
        r_o <- matrix(zA[tt, oi], ncol = 1) - signal[oi, , drop = FALSE]
        out[tt, ] <- (as.numeric(signal) + as.numeric(bl$S %*% r_o))[at_new]
      }
    }
    pred[, match(add, new)] <- out
  }
  pred
}


#' Print method for the blocked cross-validation of STEM and SC-STEM fits
#'
#' @param x an object of class \dQuote{SCSTEM_CV}.
#' @param digits integer, number of significant digits. Default is 3.
#' @param ... further arguments, currently ignored.
#'
#' @return \code{x}, invisibly. Called for its side effect of printing the
#'   errors and the ranking of the models under each scheme.
#'
#' @export
print.SCSTEM_CV <- function(x, digits = 3, ...) {
  cat("Blocked cross-validation of STEM and SC-STEM fits\n")
  cat("  folds per scheme : ", x$n_folds, "\n", sep = "")
  cat("  models           : ",
      paste0(x$models$model, " (",
             ifelse(x$models$k == 1L, "pooled",
                    paste0("k = ", x$models$k, ", phi = ", x$models$phi)), ")",
             collapse = ", "), "\n\n", sep = "")
  s <- x$summary
  s$rmse <- signif(s$rmse, digits)
  s$mae <- signif(s$mae, digits)
  print(s, row.names = FALSE)
  bad <- sum(!is.na(x$folds$message))
  if (bad) cat("\n", bad, " refit(s) failed; see $folds$message.\n", sep = "")
  invisible(x)
}
