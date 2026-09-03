#' Two-step selection of the SC-STEM hyperparameters
#'
#' @description
#' \code{SCSTEM_Select} applies the two-step tuning rule for the number of
#' clusters \eqn{k} and the spatial penalty \eqn{\phi} of a spatially-clustered
#' STEM model, given the grid of fits produced by
#' \code{\link{SCSTEM_Infocrit}}.
#'
#' @details
#' In-sample information criteria cannot be used naively here: the partition is
#' optimized on the data, so the maximized likelihood keeps an optimism bias
#' that the parameter count does not fully correct, and the cluster-wise
#' variance components shrink as the assignment step pursues within-cluster
#' homogeneity. Both effects are strongest at \eqn{\phi \approx 0}, where
#' in-sample criteria systematically overselect \eqn{k}. The rule therefore
#' restricts the criterion to a \emph{moderate-penalty band} \eqn{\Phi_M} and
#' proceeds in two steps:
#'
#' \describe{
#'   \item{(S1) Number of clusters, by BIC within the band.}{For each
#'     \eqn{\phi \in \Phi_M}, record the BIC-minimizing \eqn{k} among the
#'     admissible configurations; select the modal winner \eqn{\hat{k}} across
#'     the band, resolving ties towards the smaller \eqn{k}. Excluding
#'     \eqn{\phi \approx 0} removes the region where the optimism bias is
#'     largest. Under-selection is the harmful direction, while residual
#'     over-selection is comparatively benign for prediction, because
#'     supernumerary clusters are either small or near-duplicates of existing
#'     regimes.}
#'   \item{(S2) Spatial penalty, at the onset of the stability plateau.}{At
#'     \eqn{k = \hat{k}}, compute for each \eqn{\phi \in \Phi_M} the stability
#'     index \eqn{S(\phi)}, the average Adjusted Rand Index between the
#'     partition estimated at \eqn{\phi} and those estimated at the adjacent
#'     grid values. Since a stronger penalty makes the partition more rigid,
#'     \eqn{S(\phi)} increases with \eqn{\phi} almost by construction and its
#'     maximizer would be biased towards over-smoothing. The rule therefore
#'     selects the \emph{smallest} penalty on the plateau,
#'     \eqn{\hat{\phi} = \min\{\phi \in \Phi_M : S(\phi) \ge \max_{\Phi_M} S -
#'     \delta\}}, that is the minimal amount of spatial forcing that already
#'     delivers a reproducible partition -- a clustering-stability argument in
#'     the spirit of von Luxburg (2010).}
#' }
#'
#' Only admissible configurations, in the sense of
#' \code{\link{SCSTEM_Infocrit}}, enter the rule, and the pooled model
#' \eqn{k = 1} is always retained as the reference against which the selected
#' configuration must be justified.
#'
#' @param infocrit an object of class \dQuote{SCSTEM_Infocrit} returned by
#'   \code{\link{SCSTEM_Infocrit}}.
#' @param band numeric vector of length two giving the moderate-penalty band
#'   \eqn{\Phi_M}. Default is \code{c(0.25, 1)}. Grid values falling inside the
#'   closed interval are used. If the band contains fewer than two grid values,
#'   the whole grid is used and a warning is issued.
#' @param criterion character, the information criterion used in step (S1). One
#'   of \code{"BIC"} (default), \code{"AIC"} or \code{"KIC"}.
#' @param delta small positive tolerance defining the stability plateau in step
#'   (S2). Default is 0.05.
#'
#' @return An object of class \dQuote{SCSTEM_Select}, a list with
#' \itemize{
#'   \item \code{k_selected}, \code{phi_selected}: the selected
#'     hyperparameters.
#'   \item \code{fit}: the corresponding \dQuote{SCSTEM_Estim} object, taken
#'     from \code{infocrit$fits} without refitting.
#'   \item \code{step1}: data frame with the criterion-minimizing \eqn{k} at
#'     each penalty in the band.
#'   \item \code{step2}: data frame with the stability index \eqn{S(\phi)} at
#'     \eqn{k = \hat{k}} and the plateau threshold.
#'   \item \code{ari_to_selected}: Adjusted Rand Index between every admissible
#'     partition on the grid and the selected one.
#'   \item \code{reference}: the pooled \eqn{k = 1} row of the criteria table,
#'     when available.
#'   \item \code{band}, \code{criterion}, \code{delta}: the settings used.
#' }
#'
#' @author Paolo Maranzano \email{pmaranzano.ricercastatistica@gmail.com}
#'
#' @references
#' von Luxburg, U. (2010) \emph{Clustering stability: an overview}. Foundations
#' and Trends in Machine Learning, 2, 235--274. \doi{10.1561/2200000008}
#'
#' Cerqueti, R., Maranzano, P., Mattera, R. (2025) \emph{Spatially-clustered
#' spatial autoregressive models with application to agricultural market
#' concentration in Europe}. Journal of Agricultural, Biological and
#' Environmental Statistics. \doi{10.1007/s13253-025-00685-7}
#'
#' @examples
#' \donttest{
#' data(pm10)
#'
#' phi <- list(beta = matrix(c(3.65, 0.046, -0.904), 3, 1),
#'             sigma2eps = 0.1, sigma2omega = 0.2, theta = 0.01,
#'             G = matrix(0.77, 1, 1), Sigmaeta = matrix(0.3, 1, 1),
#'             m0 = as.matrix(0), C0 = as.matrix(1))
#'
#' mod1 <- STEM_Model(z = pm10$z, covariates = pm10$covariates,
#'                    coordinates = pm10$coords * 1000, phi = phi,
#'                    K = matrix(1, ncol(pm10$z), 1))
#'
#' ic <- SCSTEM_Infocrit(mod1, k_grid = 1:3,
#'                       phi_grid = c(0, 0.25, 0.5, 0.75, 1),
#'                       distance = "euclidean")
#' sel <- SCSTEM_Select(ic)
#' sel
#' }
#'
#' @seealso \code{\link{SCSTEM_Infocrit}} and \code{\link{SCSTEM_Estim}}
#'
#' @keywords models spatial
#'
#' @export
SCSTEM_Select <- function(infocrit,
                          band = c(0.25, 1),
                          criterion = c("BIC", "AIC", "KIC"),
                          delta = 0.05) {

  if (!inherits(infocrit, "SCSTEM_Infocrit")) {
    stop("'infocrit' must be an object of class 'SCSTEM_Infocrit' returned by SCSTEM_Infocrit().",
         call. = FALSE)
  }
  criterion <- match.arg(criterion)
  if (length(band) != 2L || anyNA(band) || band[1] > band[2]) {
    stop("'band' must be a numeric vector of length two with band[1] <= band[2].",
         call. = FALSE)
  }
  if (length(delta) != 1L || is.na(delta) || delta < 0) {
    stop("'delta' must be a single non-negative number.", call. = FALSE)
  }

  tab <- infocrit$table
  reference <- tab[tab$k == 1, , drop = FALSE]

  ### Only admissible, multi-cluster configurations enter the rule
  adm <- tab[tab$admissible & tab$k > 1, , drop = FALSE]
  if (!nrow(adm)) {
    stop("No admissible configuration with k > 1 is available: every fit on the grid ",
         "collapsed at least one cluster. Consider a smaller k_grid, a smaller ",
         "min_cluster_size, or a stronger phi_penalty.", call. = FALSE)
  }

  ##############################################
  ########## (S1) number of clusters ###########
  ##############################################

  phis <- sort(unique(adm$phi))
  in_band <- phis[phis >= band[1] & phis <= band[2]]
  if (length(in_band) < 2L) {
    warning("The moderate-penalty band [", band[1], ", ", band[2], "] contains ",
            length(in_band), " admissible grid value(s); the whole grid is used instead.",
            call. = FALSE)
    in_band <- phis
  }

  step1 <- do.call(rbind, lapply(in_band, function(p) {
    sub <- adm[adm$phi == p, , drop = FALSE]
    if (!nrow(sub)) return(NULL)
    ### ties resolved towards the smaller k
    sub <- sub[order(sub[[criterion]], sub$k), , drop = FALSE]
    data.frame(phi = p, k_best = sub$k[1], value = sub[[criterion]][1],
               stringsAsFactors = FALSE)
  }))

  votes <- table(step1$k_best)
  top <- as.integer(names(votes)[votes == max(votes)])
  k_sel <- min(top)

  ##############################################
  ########## (S2) spatial penalty ##############
  ##############################################

  sub_k <- adm[adm$k == k_sel & adm$phi %in% in_band, , drop = FALSE]
  sub_k <- sub_k[order(sub_k$phi), , drop = FALSE]
  phis_k <- sub_k$phi

  tag_of <- function(kk, pp) paste0("k=", kk, ", phi=", pp)
  grp <- function(pp) infocrit$groups[, tag_of(k_sel, pp)]

  if (length(phis_k) == 1L) {
    phi_sel <- phis_k
    step2 <- data.frame(phi = phis_k, stability = NA_real_,
                        threshold = NA_real_, on_plateau = TRUE,
                        stringsAsFactors = FALSE)
  } else {
    stab <- vapply(seq_along(phis_k), function(j) {
      neigh <- c(j - 1L, j + 1L)
      neigh <- neigh[neigh >= 1L & neigh <= length(phis_k)]
      mean(vapply(neigh, function(m) scstem_ari(grp(phis_k[j]), grp(phis_k[m])),
                  numeric(1)), na.rm = TRUE)
    }, numeric(1))
    thr <- max(stab, na.rm = TRUE) - delta
    on_plateau <- is.finite(stab) & stab >= thr
    phi_sel <- if (any(on_plateau)) phis_k[which(on_plateau)[1]] else phis_k[which.max(stab)]
    step2 <- data.frame(phi = phis_k, stability = stab, threshold = thr,
                        on_plateau = on_plateau, stringsAsFactors = FALSE)
  }

  ##############################################
  ########## Output ############################
  ##############################################

  tag_sel <- tag_of(k_sel, phi_sel)
  fit_sel <- infocrit$fits[[tag_sel]]

  sel_group <- infocrit$groups[, tag_sel]
  ari_tab <- do.call(rbind, lapply(colnames(infocrit$groups), function(tg) {
    data.frame(configuration = tg,
               ARI = scstem_ari(infocrit$groups[, tg], sel_group),
               stringsAsFactors = FALSE)
  }))

  out <- list(
    k_selected = k_sel,
    phi_selected = phi_sel,
    fit = fit_sel,
    step1 = step1,
    step2 = step2,
    ari_to_selected = ari_tab,
    reference = reference,
    selected_row = tab[tab$k == k_sel & tab$phi == phi_sel, , drop = FALSE],
    band = band,
    criterion = criterion,
    delta = delta
  )
  class(out) <- c("SCSTEM_Select", "list")
  out
}


#' Print method for the SC-STEM selection rule
#'
#' @param x an object of class \dQuote{SCSTEM_Select}.
#' @param digits integer, number of significant digits. Default is 3.
#' @param ... further arguments, currently ignored.
#'
#' @return \code{x}, invisibly. Called for its side effect of printing the two
#'   steps of the selection rule and the selected configuration.
#'
#' @export
print.SCSTEM_Select <- function(x, digits = 3, ...) {
  cat("Two-step selection of the SC-STEM hyperparameters\n")
  cat("  moderate-penalty band : [", x$band[1], ", ", x$band[2], "]\n", sep = "")
  cat("  criterion (S1)        : ", x$criterion, "\n", sep = "")
  cat("  plateau tolerance (S2): delta = ", x$delta, "\n\n", sep = "")

  cat("(S1) criterion-minimizing k within the band\n")
  s1 <- x$step1
  s1$value <- round(s1$value, digits)
  print(s1, row.names = FALSE)

  cat("\n(S2) stability of the partition at k = ", x$k_selected, "\n", sep = "")
  s2 <- x$step2
  s2$stability <- round(s2$stability, digits)
  s2$threshold <- round(s2$threshold, digits)
  print(s2, row.names = FALSE)

  cat("\nSelected configuration: k = ", x$k_selected,
      " , phi = ", x$phi_selected, "\n", sep = "")
  if (nrow(x$selected_row)) {
    sr <- x$selected_row
    cat("  loglik = ", round(sr$loglik, digits),
        " | ", x$criterion, " = ", round(sr[[x$criterion]], digits),
        " | cluster sizes >= ", sr$min_size, "\n", sep = "")
  }
  if (nrow(x$reference)) {
    cat("  pooled reference (k = 1): ", x$criterion, " = ",
        round(x$reference[[x$criterion]][1], digits), "\n", sep = "")
  }
  invisible(x)
}
