#' Inference from the SC-STEM parametric bootstrap
#'
#' @description
#' \code{SCSTEM_BootInference} post-processes the draws produced by
#' \code{\link{SCSTEM_Bootstrap}}. It aligns the partition of every refit onto
#' the clusters of the original fit, and returns standard errors, several
#' families of confidence intervals for the cluster-wise parameters, and
#' percentile tests for the differences between clusters.
#'
#' @details
#' \strong{Label alignment.} The clusters of a refit carry arbitrary labels, so
#' before any statistic is computed each refit cluster is mapped to the cluster
#' of the original fit with which it shares the largest number of locations,
#' resolving conflicts greedily by decreasing overlap (majority rule). Draws
#' whose refit did not recover all \eqn{k} clusters are discarded.
#'
#' \strong{Confidence intervals.} From the same set of aligned draws the
#' function computes, for every cluster-wise parameter, the intervals of
#' Davison and Hinkley (1997):
#' \describe{
#'   \item{\code{normal}}{\eqn{\hat\theta - \hat{b} \pm z_{1-\alpha/2}
#'     \, se^*}, with \eqn{\hat{b}} the bootstrap estimate of the bias.}
#'   \item{\code{basic}}{\eqn{(2\hat\theta - q^*_{1-\alpha/2},\;
#'     2\hat\theta - q^*_{\alpha/2})}.}
#'   \item{\code{percentile}}{\eqn{(q^*_{\alpha/2},\; q^*_{1-\alpha/2})}.}
#'   \item{\code{bc}}{the bias-corrected percentile interval, which shifts the
#'     percentile levels by the median bias of the bootstrap distribution.}
#' }
#' Studentized (bootstrap-\eqn{t}) intervals are not reported: they would
#' require an asymptotic standard error for every cluster-wise parameter at
#' every draw, which the EM estimation of a STEM model does not provide -- the
#' bootstrap is precisely the way standard errors are obtained here.
#'
#' \strong{Tests between clusters.} For every parameter and every pair of
#' clusters, the equal-tailed percentile test on the bootstrap distribution of
#' the difference gives a two-sided p-value
#' \eqn{2\min\{\Pr^*(\Delta \le 0),\, \Pr^*(\Delta \ge 0)\}}, truncated at one.
#' It answers the question of whether two estimated regimes really differ in a
#' given coefficient.
#'
#' \strong{Interpretation.} Because the partition is re-estimated at every draw,
#' these intervals incorporate the uncertainty of the clustering and are
#' therefore wider -- and more honest -- than intervals computed conditionally
#' on the estimated partition.
#'
#' @param SCSTEMboot an object of class \dQuote{SCSTEM_Bootstrap} returned by
#'   \code{\link{SCSTEM_Bootstrap}}.
#' @param level confidence level of the intervals. Default is 0.95.
#' @param parameters optional character vector selecting the parameters to
#'   report. Defaults to all the columns of \code{phi_hat}.
#' @param digits integer, rounding applied to the returned tables. Default is 5.
#'
#' @return An object of class \dQuote{SCSTEM_BootInference}, a list with
#' \itemize{
#'   \item \code{summary}: data frame with one row per cluster and parameter,
#'     reporting the original estimate, the bootstrap mean, bias, standard
#'     error, \eqn{t} statistic, normal-approximation p-value and the four
#'     families of confidence intervals.
#'   \item \code{tests}: data frame of the pairwise percentile tests between
#'     clusters.
#'   \item \code{aligned}: the aligned bootstrap draws, as a
#'     \eqn{B_{valid}} by \eqn{k} by \eqn{npar} array, so that any further
#'     statistic can be computed by the user.
#'   \item \code{coclustering}: \eqn{d} by \eqn{d} matrix of the bootstrap
#'     co-clustering frequencies, that is the share of usable draws in which
#'     two locations are assigned to the same regime. It is a direct measure of
#'     the stability of the estimated partition.
#'   \item \code{stability}: data frame with the Adjusted Rand Index between
#'     every refit partition and the original one.
#'   \item \code{level}, \code{B_used}: settings and number of draws used.
#' }
#'
#' @author Paolo Maranzano \email{pmaranzano.ricercastatistica@gmail.com}
#'
#' @references
#' Davison, A.C., Hinkley, D.V. (1997) \emph{Bootstrap Methods and their
#' Application}. Cambridge University Press.
#'
#' Efron, B., Tibshirani, R.J. (1993) \emph{An Introduction to the Bootstrap}.
#' Chapman and Hall, New York.
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
#' fit <- SCSTEM_Estimation(mod, k = 2, phi_penalty = 0.05, distance = 'geo')
#' boot <- SCSTEM_Bootstrap(fit, B = 20, seed = 1)
#' inf <- SCSTEM_BootInference(boot)
#' inf
#' inf$summary
#' }
#'
#' @seealso \code{\link{SCSTEM_Bootstrap}} and \code{\link{SCSTEM_Estimation}}
#'
#' @keywords models spatial
#'
#' @export
SCSTEM_BootInference <- function(SCSTEMboot, level = 0.95,
                                 parameters = NULL, digits = 5) {

  if (!inherits(SCSTEMboot, "SCSTEM_Bootstrap")) {
    stop("'SCSTEMboot' must be an object of class 'SCSTEM_Bootstrap' returned by SCSTEM_Bootstrap().",
         call. = FALSE)
  }
  if (length(level) != 1L || is.na(level) || level <= 0 || level >= 1) {
    stop("'level' must be a single number strictly between 0 and 1.", call. = FALSE)
  }

  orig <- SCSTEMboot$original
  k <- orig$input_args$k
  d <- orig$input_args$d
  par_names <- dimnames(SCSTEMboot$draws)[[3]]
  if (is.null(parameters)) parameters <- par_names
  parameters <- intersect(parameters, par_names)
  if (!length(parameters)) {
    stop("None of the requested parameters is available in the bootstrap draws.",
         call. = FALSE)
  }

  ok <- which(SCSTEMboot$info$ok)
  if (!length(ok)) {
    stop("No usable bootstrap draw: every refit failed or collapsed a cluster. ",
         "Consider increasing B, relaxing min_cluster_size, or strengthening the ",
         "spatial penalty.", call. = FALSE)
  }
  Bu <- length(ok)

  ##################################################
  ########## Label alignment (majority rule) #######
  ##################################################

  aligned <- array(NA_real_, dim = c(Bu, k, length(parameters)),
                   dimnames = list(paste0("b", ok), paste("cluster", seq_len(k)),
                                   parameters))
  ari <- rep(NA_real_, Bu)
  coclust <- matrix(0, d, d)

  for (m in seq_along(ok)) {
    b <- ok[m]
    grp_b <- SCSTEMboot$groups[, b]
    map <- scstem_align_labels(reference = orig$group, refit = grp_b, K = k)
    for (g in seq_len(k)) {
      target <- map[g]
      if (is.na(target)) next
      aligned[m, target, ] <- SCSTEMboot$draws[b, g, parameters]
    }
    ari[m] <- scstem_ari(grp_b, orig$group)
    coclust <- coclust + outer(grp_b, grp_b, "==")
  }
  coclust <- coclust / Bu
  dimnames(coclust) <- list(rownames(orig$df), rownames(orig$df))

  ##################################################
  ########## Confidence intervals ##################
  ##################################################

  alpha <- 1 - level
  zq <- stats::qnorm(1 - alpha / 2)
  lo_name <- paste0("lower_", round(level * 100), "_")
  up_name <- paste0("upper_", round(level * 100), "_")

  rows <- list()
  for (g in seq_len(k)) {
    for (p in parameters) {
      x <- aligned[, g, p]
      x <- x[is.finite(x)]
      est <- orig$phi_hat[g, p]
      if (length(x) < 2L || !is.finite(est)) {
        rows[[length(rows) + 1L]] <- data.frame(
          cluster = g, parameter = p, estimate = est,
          boot_mean = NA_real_, bias = NA_real_, se = NA_real_,
          t = NA_real_, p_value = NA_real_,
          normal_lo = NA_real_, normal_up = NA_real_,
          basic_lo = NA_real_, basic_up = NA_real_,
          perc_lo = NA_real_, perc_up = NA_real_,
          bc_lo = NA_real_, bc_up = NA_real_,
          n_draws = length(x), stringsAsFactors = FALSE)
        next
      }

      bmean <- mean(x)
      bias <- bmean - est
      se <- stats::sd(x)
      qs <- stats::quantile(x, probs = c(alpha / 2, 1 - alpha / 2), names = FALSE)

      ### bias-corrected percentile interval: the percentile levels are shifted
      ### by the median bias z0 of the bootstrap distribution
      prop <- mean(x < est)
      prop <- min(max(prop, 1 / (2 * length(x))), 1 - 1 / (2 * length(x)))
      z0 <- stats::qnorm(prop)
      a1 <- stats::pnorm(2 * z0 - zq)
      a2 <- stats::pnorm(2 * z0 + zq)
      bcq <- stats::quantile(x, probs = c(a1, a2), names = FALSE)

      tstat <- if (se > 0) est / se else NA_real_
      pval <- if (is.finite(tstat)) 2 * (1 - stats::pnorm(abs(tstat))) else NA_real_

      rows[[length(rows) + 1L]] <- data.frame(
        cluster = g, parameter = p, estimate = est,
        boot_mean = bmean, bias = bias, se = se,
        t = tstat, p_value = pval,
        normal_lo = est - bias - zq * se, normal_up = est - bias + zq * se,
        basic_lo = 2 * est - qs[2], basic_up = 2 * est - qs[1],
        perc_lo = qs[1], perc_up = qs[2],
        bc_lo = bcq[1], bc_up = bcq[2],
        n_draws = length(x), stringsAsFactors = FALSE)
    }
  }
  summ <- do.call(rbind, rows)

  summ$signif <- ifelse(is.na(summ$p_value), "",
                 ifelse(summ$p_value < 0.001, "***",
                 ifelse(summ$p_value < 0.01, "**",
                 ifelse(summ$p_value < 0.05, "*",
                 ifelse(summ$p_value < 0.1, ".", "")))))

  num <- vapply(summ, is.numeric, logical(1))
  num["cluster"] <- FALSE
  num["n_draws"] <- FALSE
  summ[num] <- lapply(summ[num], round, digits)
  names(summ)[names(summ) == "normal_lo"] <- paste0("normal_", lo_name, "CI")
  names(summ)[names(summ) == "normal_up"] <- paste0("normal_", up_name, "CI")
  names(summ)[names(summ) == "basic_lo"] <- paste0("basic_", lo_name, "CI")
  names(summ)[names(summ) == "basic_up"] <- paste0("basic_", up_name, "CI")
  names(summ)[names(summ) == "perc_lo"] <- paste0("perc_", lo_name, "CI")
  names(summ)[names(summ) == "perc_up"] <- paste0("perc_", up_name, "CI")
  names(summ)[names(summ) == "bc_lo"] <- paste0("bc_", lo_name, "CI")
  names(summ)[names(summ) == "bc_up"] <- paste0("bc_", up_name, "CI")
  rownames(summ) <- NULL

  ##################################################
  ########## Pairwise tests between clusters #######
  ##################################################

  tests <- NULL
  if (k >= 2L) {
    trows <- list()
    for (g in seq_len(k - 1L)) {
      for (h in (g + 1L):k) {
        for (p in parameters) {
          dlt <- aligned[, g, p] - aligned[, h, p]
          dlt <- dlt[is.finite(dlt)]
          if (length(dlt) < 2L) next
          pl <- mean(dlt <= 0)
          pu <- mean(dlt >= 0)
          pv <- min(1, 2 * min(pl, pu))
          qd <- stats::quantile(dlt, probs = c(alpha / 2, 1 - alpha / 2), names = FALSE)
          trows[[length(trows) + 1L]] <- data.frame(
            parameter = p, cluster_a = g, cluster_b = h,
            diff_estimate = orig$phi_hat[g, p] - orig$phi_hat[h, p],
            diff_boot_mean = mean(dlt),
            lower = qd[1], upper = qd[2], p_value = pv,
            n_draws = length(dlt), stringsAsFactors = FALSE)
        }
      }
    }
    if (length(trows)) {
      tests <- do.call(rbind, trows)
      tests$signif <- ifelse(tests$p_value < 0.001, "***",
                      ifelse(tests$p_value < 0.01, "**",
                      ifelse(tests$p_value < 0.05, "*",
                      ifelse(tests$p_value < 0.1, ".", ""))))
      tnum <- vapply(tests, is.numeric, logical(1))
      tnum["cluster_a"] <- FALSE
      tnum["cluster_b"] <- FALSE
      tnum["n_draws"] <- FALSE
      tests[tnum] <- lapply(tests[tnum], round, digits)
      names(tests)[names(tests) == "lower"] <- paste0(lo_name, "CI")
      names(tests)[names(tests) == "upper"] <- paste0(up_name, "CI")
      rownames(tests) <- NULL
    }
  }

  out <- list(summary = summ, tests = tests, aligned = aligned,
              coclustering = coclust,
              stability = data.frame(draw = ok, ARI = round(ari, digits),
                                     stringsAsFactors = FALSE),
              level = level, B_used = Bu)
  class(out) <- c("SCSTEM_BootInference", "list")
  out
}


#' Print method for the SC-STEM bootstrap inference
#'
#' @param x an object of class \dQuote{SCSTEM_BootInference}.
#' @param ... further arguments, currently ignored.
#'
#' @return \code{x}, invisibly. Called for its side effect of printing the
#'   cluster-wise summary table.
#'
#' @export
print.SCSTEM_BootInference <- function(x, ...) {
  cat("SC-STEM bootstrap inference (", x$B_used, " usable draws, level = ",
      x$level, ")\n\n", sep = "")
  keep <- c("cluster", "parameter", "estimate", "se", "t", "p_value", "signif")
  keep <- intersect(keep, names(x$summary))
  print(x$summary[, keep, drop = FALSE], row.names = FALSE)
  cat("\nMedian ARI of the refit partitions against the original: ",
      round(stats::median(x$stability$ARI, na.rm = TRUE), 3), "\n", sep = "")
  cat("Use $summary for the confidence intervals, $tests for the between-cluster\n",
      "comparisons and $coclustering for the stability of the partition.\n", sep = "")
  invisible(x)
}
