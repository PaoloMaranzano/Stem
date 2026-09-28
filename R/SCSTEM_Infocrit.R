#' Information criteria over a grid of SC-STEM configurations
#'
#' @description
#' \code{SCSTEM_Infocrit} fits the spatially-clustered STEM model over a grid of
#' numbers of clusters \eqn{k} and spatial penalties \eqn{\phi}, and returns the
#' likelihood-based information criteria of every configuration together with
#' the diagnostics needed to decide which of them are admissible. It is the
#' input of the two-step tuning rule implemented in
#' \code{\link{SCSTEM_Select}}.
#'
#' @details
#' For each pair \eqn{(k, \phi)} the function calls \code{\link{SCSTEM_Estimation}}
#' and records the exact total log-likelihood of the final refit, the number of
#' free parameters \eqn{k_{eff}(r + 3 + 3p)}, and
#' \deqn{AIC = -2\ell + 2 K, \qquad BIC = -2\ell + \log(dT) K, \qquad
#'       KIC = -2\ell + 3 K,}
#' where \eqn{d} is the number of locations and \eqn{T} the number of time
#' points, so that \eqn{dT} is the number of observations entering the
#' likelihood.
#'
#' \strong{Admissibility.} A configuration is admissible when all its \eqn{k}
#' clusters could be re-estimated on the final partition, that is when
#' \code{k_eff == k}. A configuration with collapsed clusters has fewer
#' effective regimes than its nominal \eqn{k}, its likelihood is not comparable
#' with the others, and no parametric bootstrap can be generated from it: the
#' \code{admissible} column flags these cases and
#' \code{\link{SCSTEM_Select}} discards them.
#'
#' \strong{Caveats.} Two warnings apply to any likelihood-based comparison in
#' this setting. First, the partition is itself optimized on the data, so the
#' maximized likelihood retains an optimism bias that a parameter count of the
#' form \eqn{k(r + 3 + 3p)} does not fully correct; in-sample criteria therefore
#' tend to favor small \eqn{\phi} and large \eqn{k}. Second, the cluster-wise
#' variance components shrink as the assignment step pursues within-cluster
#' homogeneity, which inflates the likelihood of the configurations with the
#' weakest spatial regularization. Both effects are strongest at
#' \eqn{\phi \approx 0}, which is why \code{\link{SCSTEM_Select}} chooses
#' \eqn{k} within a moderate band of penalties.
#'
#' @param StemModel an object of class \dQuote{STEM_Model} given as output by
#'   the \code{\link{STEM_Model}} function.
#' @param k_grid integer vector of candidate numbers of clusters. Default is
#'   \code{1:3}. The pooled model \eqn{k = 1} is always a useful reference.
#' @param phi_grid numeric vector of candidate non-negative spatial penalties.
#'   Default is \code{c(0, 0.025, 0.05, 0.1, 0.2, 0.5, 1)}: the moderate values
#'   up to 0.2 form the band in which \code{\link{SCSTEM_Select}} chooses
#'   \eqn{k}, and the strong values 0.5 and 1 compete only in the choice of
#'   \eqn{\phi}. At every \eqn{k} the fit at \eqn{\phi = 0} starts the penalized
#'   ones (see \code{\link{SCSTEM_Estimation}}), so the grid should contain 0.
#' @param mink,maxk deprecated scalars kept for backward compatibility with
#'   versions of the package prior to 2.0.0. When supplied and \code{k_grid} is
#'   missing, the grid is set to \code{mink:maxk}.
#' @param verbose logical. If \code{TRUE}, the progress over the grid is
#'   reported via \code{message()}. Default is \code{FALSE}.
#' @param ... further arguments passed to \code{\link{SCSTEM_Estimation}}, such as
#'   \code{distance}, \code{knn}, \code{phi_scale}, \code{label_update},
#'   \code{precision} or \code{min_cluster_size}.
#'
#' @return An object of class \dQuote{SCSTEM_Infocrit}: a list with
#' \itemize{
#'   \item \code{table}: a data frame with one row per configuration and columns
#'     \code{k}, \code{phi}, \code{k_eff}, \code{admissible}, \code{loglik},
#'     \code{npar}, \code{AIC}, \code{BIC}, \code{KIC}, \code{min_size},
#'     \code{convergence}.
#'   \item \code{fits}: the list of \dQuote{SCSTEM_Estimation} objects, named
#'     \code{"k=<k>, phi=<phi>"}, so that the selected configuration can be used
#'     without refitting.
#'   \item \code{groups}: matrix of the estimated partitions, one column per
#'     configuration, which \code{\link{SCSTEM_Select}} compares with the
#'     selected one.
#'   \item \code{failed}: a data frame listing the configurations whose fit
#'     raised an error, with the error message.
#' }
#'
#' @author Paolo Maranzano \email{pmaranzano.ricercastatistica@gmail.com},
#'   Francesco Caccia, Michela Cameletti
#'
#' @references
#' Cerqueti, R., Maranzano, P., Mattera, R. (2025) \emph{Spatially-clustered
#' spatial autoregressive models with application to agricultural market
#' concentration in Europe}. Journal of Agricultural, Biological and
#' Environmental Statistics. \doi{10.1007/s13253-025-00685-7}
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
#' ic <- SCSTEM_Infocrit(mod, k_grid = 1:3, phi_grid = c(0, 0.05),
#'                       distance = 'geo')
#' ic
#' }
#'
#' @seealso \code{\link{SCSTEM_Estimation}} and \code{\link{SCSTEM_Select}}
#'
#' @keywords models spatial
#'
#' @export
SCSTEM_Infocrit <- function(StemModel,
                            k_grid = 1:3,
                            phi_grid = c(0, 0.025, 0.05, 0.1, 0.2, 0.5, 1),
                            mink = NULL,
                            maxk = NULL,
                            verbose = FALSE,
                            ...) {

  if (!inherits(StemModel, "STEM_Model")) {
    stop("'StemModel' must be an object of class 'STEM_Model'.", call. = FALSE)
  }

  ### Backward compatibility with the pre-2.0.0 signature
  ### SCSTEM_Infocrit(StemModel, mink, maxk, ...)
  if (missing(k_grid) && !is.null(mink) && !is.null(maxk)) {
    k_grid <- seq.int(from = as.integer(mink), to = as.integer(maxk))
  }

  k_grid <- sort(unique(as.integer(k_grid)))
  phi_grid <- sort(unique(as.numeric(phi_grid)))
  if (any(k_grid < 1)) stop("'k_grid' must contain positive integers.", call. = FALSE)
  if (any(phi_grid < 0)) stop("'phi_grid' must contain non-negative values.", call. = FALSE)

  d <- ncol(StemModel$data$z)

  ### The pooled model does not depend on phi: it is fitted once and reused.
  grid <- expand.grid(k = k_grid, phi = phi_grid, KEEP.OUT.ATTRS = FALSE)
  grid <- grid[!(grid$k == 1 & grid$phi != phi_grid[1]), , drop = FALSE]
  grid <- grid[order(grid$k, grid$phi), , drop = FALSE]
  rownames(grid) <- NULL
  nconf <- nrow(grid)

  res <- data.frame(
    k = integer(0), phi = numeric(0), k_eff = integer(0),
    admissible = logical(0), loglik = numeric(0), npar = numeric(0),
    AIC = numeric(0), BIC = numeric(0), KIC = numeric(0),
    min_size = integer(0), convergence = character(0),
    stringsAsFactors = FALSE
  )
  fits <- list()
  groups <- matrix(NA_integer_, nrow = d, ncol = 0)
  gnames <- character(0)
  failed <- data.frame(k = integer(0), phi = numeric(0),
                       message = character(0), stringsAsFactors = FALSE)

  ### A penalized fit starts from the solution of the unpenalized fit at the
  ### same k (see SCSTEM_Estimation()). The grid visits phi = 0 first at every
  ### k and hands its partition on, so it is not fitted again for every phi.
  start_k <- list()

  ### Under init_method = "departures" (the default) an unpenalized fit starts
  ### from the k-means of the departures of the locations from the pooled fit.
  ### The grid visits k = 1 first, so the departures are computed once from its
  ### pooled fit and every unpenalized fit at k > 1 starts from their k-means,
  ### exactly as it would on its own, without fitting the pooled model again.
  dots <- list(...)
  dep_start <- is.null(dots$init_partition) &&
    (is.null(dots$init_method) || identical(dots$init_method[1], "departures"))
  dep_feat <- NULL
  seed0 <- if ("seed" %in% names(dots)) dots$seed else formals(SCSTEM_Estimation)$seed
  ms0 <- if (is.null(dots$min_cluster_size)) ncol(StemModel$data$covariates) + 2L else
    dots$min_cluster_size
  ms0 <- max(2L, as.integer(ms0))

  for (j in seq_len(nconf)) {

    kk <- grid$k[j]
    pp <- grid$phi[j]
    tag <- paste0("k=", kk, ", phi=", pp)

    if (isTRUE(verbose)) {
      message("[", j, "/", nconf, "] fitting ", tag, " ...")
    }

    args_j <- c(list(StemModel = StemModel, k = kk, phi_penalty = pp,
                     verbose = FALSE), list(...))
    if (pp > 0 && !is.null(start_k[[as.character(kk)]])) {
      args_j$init_partition <- start_k[[as.character(kk)]]
    }
    if (kk > 1 && pp == 0 && dep_start && !is.null(dep_feat)) {
      args_j$init_partition <- scstem_with_seed(
        seed0,
        scstem_init(Xmeans = dep_feat, coords = StemModel$data$coordinates,
                    k = kk, method = "kmeans", min_size = ms0))
    }
    fit <- tryCatch(
      suppressWarnings(do.call(SCSTEM_Estimation, args_j)),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      failed <- rbind(failed, data.frame(k = kk, phi = pp,
                                         message = conditionMessage(fit),
                                         stringsAsFactors = FALSE))
      next
    }

    k_eff <- sum(fit$final_refit)
    sizes <- tabulate(fit$group, nbins = kk)
    res <- rbind(res, data.frame(
      k = kk, phi = pp, k_eff = k_eff,
      admissible = (k_eff == kk) && all(sizes > 0),
      loglik = unname(fit$info_crit[["loglik"]]),
      npar = unname(fit$info_crit[["k"]]),
      AIC = unname(fit$info_crit[["AIC"]]),
      BIC = unname(fit$info_crit[["BIC"]]),
      KIC = unname(fit$info_crit[["KIC"]]),
      min_size = min(sizes),
      convergence = fit$convergence,
      stringsAsFactors = FALSE
    ))
    if (kk > 1 && pp == 0 && length(unique(fit$group)) == kk) {
      start_k[[as.character(kk)]] <- fit$group
    }
    if (kk == 1 && dep_start && is.null(dep_feat) && length(fit$fit_list) == 1L) {
      dep_feat <- tryCatch(
        scstem_departures(StemModel$data$z, StemModel$data$covariates,
                          STEM_Signal(fit$fit_list[[1]]), nrow(StemModel$data$z)),
        error = function(e) NULL)
    }
    fits[[tag]] <- fit
    groups <- cbind(groups, fit$group)
    gnames <- c(gnames, tag)
  }

  if (ncol(groups)) colnames(groups) <- gnames
  rownames(res) <- NULL

  if (isTRUE(verbose) && nrow(failed)) {
    message("* ", nrow(failed), " configuration(s) failed and were dropped.")
  }

  out <- list(table = res, fits = fits, groups = groups, failed = failed,
              k_grid = k_grid, phi_grid = phi_grid)
  class(out) <- c("SCSTEM_Infocrit", "list")
  out
}


#' Print method for SC-STEM information criteria
#'
#' @param x an object of class \dQuote{SCSTEM_Infocrit}.
#' @param digits integer, number of significant digits. Default is 2.
#' @param ... further arguments, currently ignored.
#'
#' @return \code{x}, invisibly. Called for its side effect of printing the grid
#'   of information criteria.
#'
#' @export
print.SCSTEM_Infocrit <- function(x, digits = 2, ...) {
  tab <- x$table
  num <- c("loglik", "AIC", "BIC", "KIC")
  tab[num] <- lapply(tab[num], round, digits)
  cat("SC-STEM information criteria over ", nrow(tab), " configuration(s)\n\n", sep = "")
  print(tab, row.names = FALSE)
  if (nrow(x$failed)) {
    cat("\n", nrow(x$failed), " configuration(s) failed:\n", sep = "")
    print(x$failed, row.names = FALSE)
  }
  invisible(x)
}
