#' Refit-with-clustering parametric bootstrap for SC-STEM models
#'
#' @description
#' \code{SCSTEM_Bootstrap} implements the parametric bootstrap for
#' spatially-clustered STEM models. Bootstrap datasets are generated cluster by
#' cluster from the fitted model and the \emph{entire} SC-STEM procedure --
#' including the endogenous partitioning -- is re-estimated on every bootstrap
#' sample, so that the uncertainty of the partition is propagated to all
#' bootstrap statistics.
#'
#' @details
#' The distinguishing feature of the procedure is that it re-runs the clustering
#' at every draw, rather than conditioning on the estimated partition. The
#' algorithm adapts the refit-with-clustering bootstrap of Maranzano, Mattera
#' and Sugasawa (2026+) to the spatio-temporal setting:
#' \enumerate{
#'   \item Fit the SC-STEM model on the observed data at fixed \eqn{K} and
#'     \eqn{\phi}, obtaining the partition \eqn{\hat{P}} and the cluster-wise
#'     parameter estimates \eqn{\hat{\phi}_k}.
#'   \item For \eqn{b = 1,\ldots,B} and for each cluster \eqn{k}, simulate a
#'     spatio-temporal dataset \eqn{z^{*(b)}} from the STEM model of cluster
#'     \eqn{k} evaluated at \eqn{\hat{\phi}_k}, using the observed covariates,
#'     coordinates and loading matrix of the locations in that cluster. This is
#'     \code{\link{STEM_Simulation}} applied cluster by cluster, so that the
#'     within-cluster spatial covariance and the latent temporal dynamics of the
#'     fitted model are both reproduced.
#'   \item Restore the original ordering of the locations, so that the simulated
#'     data are aligned with the rows of the spatial penalty graph and every
#'     location keeps its own neighbors.
#'   \item Re-run \code{\link{SCSTEM_Estimation}} on \eqn{z^{*(b)}} with the same
#'     \eqn{K}, \eqn{\phi} and algorithmic settings, and store the cluster-wise
#'     estimates, the refit partition and the convergence diagnostics.
#' }
#' Because each bootstrap sample is generated cluster by cluster from the
#' observed covariates and coordinates, the spatial structure of the original
#' data is preserved: the fixed-effects surface inherits the contiguity of the
#' estimated regimes while every location retains its position in the
#' neighborhood graph.
#'
#' \strong{Failed refits.} A bootstrap sample whose refit raises an error, or
#' whose refit collapses a cluster, is recorded but excluded from the
#' inferential summaries computed by \code{\link{SCSTEM_BootInference}}. The
#' share of usable draws is reported.
#'
#' \strong{Label alignment.} The clusters of a refit carry arbitrary labels, so
#' they are relocated onto the clusters of the original fit by the majority
#' rule before any statistic is computed. The alignment is applied by
#' \code{\link{SCSTEM_BootInference}}, and the raw draws are returned here
#' unaligned, together with the refit partitions, so that any alternative
#' matching can be applied afterwards.
#'
#' @param SCSTEM an object of class \dQuote{SCSTEM_Estimation} returned by
#'   \code{\link{SCSTEM_Estimation}}.
#' @param B integer, the number of bootstrap replicates. Default is 100.
#' @param seed integer or \code{NULL}. When supplied, the draws are exactly
#'   reproducible; the RNG stream is restored on exit. Default is \code{NULL}.
#' @param verbose logical. If \code{TRUE}, progress is reported via
#'   \code{message()}. Default is \code{FALSE}.
#' @param ... further arguments passed to \code{\link{SCSTEM_Estimation}} for the
#'   refits, overriding the settings of the original fit.
#'
#' @return An object of class \dQuote{SCSTEM_Bootstrap}, a list with
#' \itemize{
#'   \item \code{draws}: a \eqn{B} by \eqn{K} by \eqn{npar} array of the
#'     cluster-wise parameter estimates of every refit, with \code{NA} for
#'     failed draws.
#'   \item \code{groups}: a \eqn{d} by \eqn{B} matrix of the refit partitions.
#'   \item \code{loglik}: matrix of the cluster-wise log-likelihoods of the
#'     refits.
#'   \item \code{info}: a data frame with one row per draw reporting whether the
#'     refit succeeded, how many clusters it recovered, its total
#'     log-likelihood, the number of iterations and the convergence message.
#'   \item \code{original}: the \dQuote{SCSTEM_Estimation} object given in input.
#'   \item \code{B}, \code{B_valid}: requested and usable number of draws.
#' }
#'
#' @author Paolo Maranzano \email{pmaranzano.ricercastatistica@gmail.com},
#'   Francesco Caccia, Michela Cameletti
#'
#' @references
#' Chatterjee, S., Lahiri, P., Li, H. (2008) \emph{Parametric bootstrap
#' approximation to the distribution of EBLUP and related prediction intervals
#' in linear mixed models}. Annals of Statistics, 36, 1221--1245.
#' \doi{10.1214/07-AOS512}
#'
#' Davison, A.C., Hinkley, D.V. (1997) \emph{Bootstrap Methods and their
#' Application}. Cambridge University Press.
#'
#' Hall, P., Maiti, T. (2006) \emph{On parametric bootstrap methods for small
#' area prediction}. Journal of the Royal Statistical Society, Series B, 68,
#' 221--238. \doi{10.1111/j.1467-9868.2006.00541.x}
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
#'                   phi = phi, A = matrix(1, d, 1))
#'
#' \donttest{
#' fit <- SCSTEM_Estimation(mod, K = 2, phi_penalty = 0.05, distance = 'geo')
#' boot <- SCSTEM_Bootstrap(fit, B = 20, seed = 1)
#' boot
#' }
#'
#' @seealso \code{\link{SCSTEM_Estimation}}, \code{\link{SCSTEM_BootInference}} and
#'   \code{\link{STEM_Bootstrap}}
#'
#' @keywords models spatial
#'
#' @export
SCSTEM_Bootstrap <- function(SCSTEM, B = 100, seed = NULL, verbose = FALSE, ...) {

  if (!inherits(SCSTEM, "SCSTEM_Estimation")) {
    stop("'SCSTEM' must be an object of class 'SCSTEM_Estimation' returned by SCSTEM_Estimation().",
         call. = FALSE)
  }
  if (length(B) != 1L || is.na(B) || B < 1 || B != round(B)) {
    stop("'B' must be a single positive integer.", call. = FALSE)
  }
  B <- as.integer(B)

  args <- SCSTEM$input_args
  base_model <- args$StemModel
  K <- args$K
  d <- args$d
  Tobs <- args$Tobs

  if (!any(SCSTEM$final_refit)) {
    stop("The fitted model has no re-estimated cluster: no bootstrap sample can be generated.",
         call. = FALSE)
  }
  if (any(!SCSTEM$final_refit)) {
    warning("Cluster(s) ", paste(which(!SCSTEM$final_refit), collapse = ", "),
            " were not re-estimated on the final partition: the bootstrap data ",
            "for their locations cannot be generated and the draws will be ",
            "restricted to the remaining clusters.", call. = FALSE)
  }

  par_names <- colnames(SCSTEM$phi_hat)
  npar <- length(par_names)

  draws <- array(NA_real_, dim = c(B, K, npar),
                 dimnames = list(paste0("b", seq_len(B)),
                                 paste("cluster", seq_len(K)), par_names))
  groups <- matrix(NA_integer_, nrow = d, ncol = B,
                   dimnames = list(NULL, paste0("b", seq_len(B))))
  loglik <- matrix(NA_real_, nrow = B, ncol = K,
                   dimnames = list(paste0("b", seq_len(B)),
                                   paste("cluster", seq_len(K))))
  info <- data.frame(draw = seq_len(B), ok = FALSE, K_eff = NA_integer_,
                     loglik = NA_real_, iter = NA_integer_,
                     convergence = NA_character_, message = NA_character_,
                     stringsAsFactors = FALSE)

  ### Arguments of the refits: the settings of the original fit, overridable
  ### through `...`. The refit must reproduce the estimation actually performed,
  ### so the penalty is passed on the scale that was effectively applied.
  refit_args <- list(K = K, phi_penalty = args$phi_penalty,
                     phi_scale = args$phi_scale, knn = args$knn,
                     distance = args$distance, init_method = args$init_method,
                     label_update = args$label_update,
                     regularization = args$regularization,
                     ### the computational settings of the original fit; a fit
                     ### made before `control` existed carries them one by one
                     control = args$control,
                     min_cluster_size = args$min_cluster_size,
                     enforce_min_size = args$enforce_min_size,
                     swap_pass = args$swap_pass, share2conv = args$share2conv,
                     seed = args$seed, verbose = FALSE,
                     ### The settings that shape the estimator itself. Before
                     ### they were carried over, a fit with a ridge on the
                     ### coefficients was resampled WITHOUT it, so the bootstrap
                     ### described a different estimator from the one reported.
                     alpha = args$alpha,
                     lambda = args$lambda, penalize = args$penalize,
                     lambda_scale = args$lambda_scale,
                     lambda_by = args$lambda_by, latent = args$latent,
                     spatial = args$spatial)
  if (is.null(args$control)) {
    refit_args[c("precision", "precision_full_dataset", "max_iter", "abs_tol", "rel_tol")] <-
      args[c("precision", "precision_full_dataset", "max_iter", "abs_tol", "rel_tol")]
  }
  user_args <- list(...)
  refit_args[names(user_args)] <- user_args
  ### every bootstrap sample is new data: no refit may be shared across them
  refit_args$refit_cache <- NULL
  ### Keep only what SCSTEM_Estimation() actually accepts. The settings are carried
  ### over from a stored list, so without this filter a change in the signature
  ### of SCSTEM_Estimation() would make every refit fail at run time with an
  ### "unused argument" error instead of being caught at build time.
  refit_args <- refit_args[names(refit_args) %in% names(formals(SCSTEM_Estimation))]
  refit_args <- refit_args[!vapply(refit_args, is.null, logical(1))]

  ### Cluster-wise generating models: the fitted STEM model of each cluster with
  ### the skeleton set at its own ML estimates, so that STEM_Simulation() draws
  ### from the estimated model rather than from the starting values.
  gen <- vector("list", K)
  for (g in which(SCSTEM$final_refit)) {
    gen[[g]] <- SCSTEM$fit_list[[g]]
    gen[[g]]$skeleton$phi <- gen[[g]]$estimates$phi.hat
  }

  na_pattern <- is.na(base_model$data$z)

  one_draw <- function(b) {

    z_star <- matrix(NA_real_, nrow = Tobs, ncol = d)
    for (g in which(SCSTEM$final_refit)) {
      idx <- SCSTEM$idx_g[[g]]
      if (!length(idx)) next
      sim_g <- STEM_Simulation(gen[[g]], distance = args$distance)
      z_star[, idx] <- as.matrix(sim_g)
    }
    ### locations of a cluster that could not be re-estimated keep their
    ### observed series, so that the refit still sees a complete network
    miss <- which(apply(z_star, 2, function(cc) all(is.na(cc))))
    if (length(miss)) z_star[, miss] <- base_model$data$z[, miss]

    ### The replicate has to reproduce the observed design, missing values
    ### included: a bootstrap sample with a complete response would understate
    ### the uncertainty of a fit obtained from an incomplete one.
    if (any(na_pattern)) z_star[na_pattern] <- NA_real_

    mod_star <- STEM_Model(z = z_star,
                           covariates = base_model$data$covariates,
                           coordinates = base_model$data$coordinates,
                           phi = base_model$skeleton$phi,
                           A = base_model$skeleton$A)

    do.call(SCSTEM_Estimation, c(list(StemModel = mod_star), refit_args))
  }

  ### The whole loop is evaluated under the requested seed. scstem_with_seed()
  ### restores the RNG stream on exit, so the caller's workspace is untouched.
  res <- scstem_with_seed(seed, lapply(seq_len(B), function(b) {
    if (isTRUE(verbose) && (b %% 10L == 0L || b == 1L)) {
      message("* bootstrap draw ", b, " / ", B)
    }
    fit_b <- tryCatch(suppressWarnings(one_draw(b)), error = function(e) e)
    if (inherits(fit_b, "error")) {
      return(list(ok = FALSE, message = conditionMessage(fit_b)))
    }
    list(ok = all(fit_b$final_refit),
         message = NA_character_,
         group = fit_b$group,
         phi_hat = fit_b$phi_hat,
         loglik_g = fit_b$loglik_g,
         K_eff = sum(fit_b$final_refit),
         loglik = unname(fit_b$info_crit[["loglik"]]),
         iter = if (nrow(fit_b$obj_trace)) max(fit_b$obj_trace$iter) else NA_integer_,
         convergence = fit_b$convergence)
  }))

  for (b in seq_len(B)) {
    rb <- res[[b]]
    info$message[b] <- rb$message
    if (is.null(rb$group)) next
    groups[, b] <- rb$group
    draws[b, , ] <- rb$phi_hat
    loglik[b, ] <- rb$loglik_g
    info$ok[b] <- rb$ok
    info$K_eff[b] <- rb$K_eff
    info$loglik[b] <- rb$loglik
    info$iter[b] <- rb$iter
    info$convergence[b] <- rb$convergence
  }

  B_valid <- sum(info$ok)
  if (isTRUE(verbose)) {
    message("Bootstrap completed: ", B_valid, " / ", B, " usable draws.")
  }
  if (B_valid == 0L) {
    warning("No bootstrap draw produced a fully admissible refit; the inferential ",
            "summaries will not be computable.", call. = FALSE)
  }

  out <- list(draws = draws, groups = groups, loglik = loglik, info = info,
              original = SCSTEM, B = B, B_valid = B_valid,
              refit_args = refit_args)
  class(out) <- c("SCSTEM_Bootstrap", "list")
  out
}


#' Print method for the SC-STEM parametric bootstrap
#'
#' @param x an object of class \dQuote{SCSTEM_Bootstrap}.
#' @param ... further arguments, currently ignored.
#'
#' @return \code{x}, invisibly. Called for its side effect of printing a
#'   compact summary of the bootstrap run.
#'
#' @export
print.SCSTEM_Bootstrap <- function(x, ...) {
  cat("Refit-with-clustering parametric bootstrap for SC-STEM\n")
  cat("  replicates requested : ", x$B, "\n", sep = "")
  cat("  usable replicates    : ", x$B_valid,
      " (", round(100 * x$B_valid / x$B, 1), "%)\n", sep = "")
  cat("  clusters             : ", x$original$input_args$K, "\n", sep = "")
  cat("  spatial penalty      : phi = ", x$original$input_args$phi_penalty, "\n", sep = "")
  nfail <- sum(!is.na(x$info$message))
  if (nfail) cat("  refits raising an error: ", nfail, "\n", sep = "")
  invisible(x)
}
