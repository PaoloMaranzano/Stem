#' Two-step selection of the SC-STEM hyperparameters
#'
#' @description
#' \code{SCSTEM_Select} applies the two-step tuning rule for the number of
#' clusters \eqn{K} and the spatial penalty \eqn{\phi} of a spatially-clustered
#' STEM model, given the grid of fits produced by
#' \code{\link{SCSTEM_Infocrit}}.
#'
#' @details
#' In-sample information criteria cannot be used naively here: the partition is
#' optimized on the data, so the maximized likelihood keeps an optimism bias
#' that the parameter count does not fully correct, and the cluster-wise
#' variance components shrink as the assignment step pursues within-cluster
#' homogeneity. Both effects are strongest at \eqn{\phi \approx 0}, where
#' in-sample criteria systematically overselect \eqn{K}. The rule therefore
#' chooses \eqn{K} within a \emph{moderate-penalty band} \eqn{\Phi_M} and
#' proceeds in two steps:
#'
#' \describe{
#'   \item{(S1) Number of clusters, by BIC within the band.}{For each
#'     \eqn{\phi \in \Phi_M}, record the BIC-minimizing \eqn{K} among the
#'     admissible configurations at that \eqn{\phi} and the pooled model
#'     \eqn{K = 1}, whose criterion does not depend on \eqn{\phi}; select the
#'     modal winner \eqn{\hat{K}} across the band, resolving ties towards the
#'     smaller \eqn{K}. If \eqn{\hat{K} = 1} the data do not support a
#'     partition: the pooled fit is returned and (S2) does not apply. Excluding
#'     \eqn{\phi \approx 0} removes the region where the optimism bias is
#'     largest. Under-selection is the harmful direction, while residual
#'     over-selection is comparatively benign for prediction, because
#'     supernumerary clusters are either small or near-duplicates of existing
#'     regimes.}
#'   \item{(S2) Spatial penalty, by the same criterion at \eqn{\hat{K}}.}{Among
#'     the admissible fits with \eqn{K = \hat{K}} on the whole grid, select the
#'     one with the smallest criterion, resolving ties towards the smaller
#'     \eqn{\phi}. At a given \eqn{K} the parameter count does not depend on
#'     \eqn{\phi} (up to the effective degrees of freedom of a penalty on the
#'     coefficients), so the step compares the likelihoods reached by the
#'     partitions estimated at different penalties. Every penalized fit starts
#'     from the unpenalized one at the same \eqn{K} (see
#'     \code{\link{SCSTEM_Estimation}}): a penalty is selected only when it has
#'     moved the partition to one with a higher likelihood than the
#'     unpenalized solution, and the optimism bias, which favors
#'     \eqn{\phi = 0}, makes the step conservative. Penalties outside the band
#'     can therefore be kept on the grid: they are selected only when the
#'     likelihood supports them. This is the rule of Cerqueti, Maranzano and
#'     Mattera (2025).}
#' }
#'
#' Only admissible configurations, in the sense of
#' \code{\link{SCSTEM_Infocrit}}, enter the rule. The pooled model enters it
#' whenever \code{K_grid} contains 1, as it should: it is the answer when no
#' partition improves on it, and it is selected when no configuration with
#' \eqn{K > 1} is admissible.
#'
#' @param infocrit an object of class \dQuote{SCSTEM_Infocrit} returned by
#'   \code{\link{SCSTEM_Infocrit}}.
#' @param band numeric vector of length two giving the moderate-penalty band
#'   \eqn{\Phi_M} of step (S1). Default is \code{c(0.025, 0.2)}, the moderate
#'   non-zero values of the default grid of \code{\link{SCSTEM_Infocrit}}; its
#'   strong values 0.5 and 1 lie outside the band and enter step (S2) only.
#'   Grid values falling inside the closed interval are used. If the band
#'   contains fewer than two grid values, the whole grid is used and a warning
#'   is issued.
#' @param criterion character, the information criterion used in both steps.
#'   One of \code{"BIC"} (default), \code{"AIC"} or \code{"KIC"}.
#'
#' @return An object of class \dQuote{SCSTEM_Select}, a list with
#' \itemize{
#'   \item \code{K_selected}, \code{phi_selected}: the selected
#'     hyperparameters. When \eqn{\hat{K} = 1}, \code{phi_selected} is the grid
#'     value at which the pooled fit is stored, which has no effect on it.
#'   \item \code{fit}: the corresponding fitted object, taken from
#'     \code{infocrit$fits} without refitting.
#'   \item \code{step1}: data frame with the criterion-minimizing \eqn{K} at
#'     each penalty in the band.
#'   \item \code{step2}: data frame with the criterion of every admissible fit
#'     at \eqn{K = \hat{K}}, one row per penalty; \code{NULL} when
#'     \eqn{\hat{K} = 1}.
#'   \item \code{ari_to_selected}: Adjusted Rand Index between every admissible
#'     partition on the grid and the selected one.
#'   \item \code{reference}: the pooled \eqn{K = 1} row of the criteria table,
#'     when available.
#'   \item \code{band}, \code{criterion}: the settings used.
#' }
#'
#' @author Paolo Maranzano \email{pmaranzano.ricercastatistica@gmail.com}
#'
#' @references
#' Cerqueti, R., Maranzano, P., Mattera, R. (2025) \emph{Spatially-clustered
#' spatial autoregressive models with application to agricultural market
#' concentration in Europe}. Journal of Agricultural, Biological and
#' Environmental Statistics. \doi{10.1007/s13253-025-00685-7}
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
#' ic <- SCSTEM_Infocrit(mod, K_grid = 1:3, distance = 'geo')
#' sel <- SCSTEM_Select(ic)
#' sel
#' }
#'
#' @seealso \code{\link{SCSTEM_Infocrit}} and \code{\link{SCSTEM_Estimation}}
#'
#' @keywords models spatial
#'
#' @export
SCSTEM_Select <- function(infocrit,
                          band = c(0.025, 0.2),
                          criterion = c("BIC", "AIC", "KIC")) {

  if (!inherits(infocrit, "SCSTEM_Infocrit")) {
    stop("'infocrit' must be an object of class 'SCSTEM_Infocrit' returned by SCSTEM_Infocrit().",
         call. = FALSE)
  }
  criterion <- match.arg(criterion)
  if (length(band) != 2L || anyNA(band) || band[1] > band[2]) {
    stop("'band' must be a numeric vector of length two with band[1] <= band[2].",
         call. = FALSE)
  }

  tab <- infocrit$table
  reference <- tab[tab$K == 1, , drop = FALSE]

  ### Only admissible configurations enter the rule. The pooled model is fitted
  ### once, at the first grid value of phi, and competes at every phi of the
  ### band, since its criterion does not depend on the penalty.
  adm <- tab[tab$admissible & tab$K > 1, , drop = FALSE]
  pooled <- tab[tab$admissible & tab$K == 1 & is.finite(tab[[criterion]]), , drop = FALSE]
  pooled <- pooled[seq_len(min(1L, nrow(pooled))), , drop = FALSE]
  if (!nrow(adm) && !nrow(pooled)) {
    stop("No admissible configuration is available: every fit on the grid ",
         "collapsed at least one cluster and the grid holds no pooled fit. ",
         "Include K = 1 in K_grid.", call. = FALSE)
  }

  ##############################################
  ########## (S1) number of clusters ###########
  ##############################################

  phis <- sort(unique(adm$phi))
  in_band <- phis[phis >= band[1] & phis <= band[2]]
  if (length(phis) && length(in_band) < 2L) {
    warning("The moderate-penalty band [", band[1], ", ", band[2], "] contains ",
            length(in_band), " admissible grid value(s); the whole grid is used instead.",
            call. = FALSE)
    in_band <- phis
  }

  step1 <- do.call(rbind, lapply(in_band, function(p) {
    sub <- rbind(adm[adm$phi == p, , drop = FALSE], pooled)
    ### ties resolved towards the smaller K
    sub <- sub[order(sub[[criterion]], sub$K), , drop = FALSE]
    data.frame(phi = p, K_best = sub$K[1], value = sub[[criterion]][1],
               stringsAsFactors = FALSE)
  }))
  if (is.null(step1)) {
    ### no admissible partition at all: the pooled model is the answer
    step1 <- data.frame(phi = pooled$phi, K_best = 1L, value = pooled[[criterion]],
                        stringsAsFactors = FALSE)
  }

  votes <- table(step1$K_best)
  top <- as.integer(names(votes)[votes == max(votes)])
  K_sel <- min(top)

  ##############################################
  ########## (S2) spatial penalty ##############
  ##############################################

  tag_of <- function(kk, pp) paste0("K=", kk, ", phi=", pp)

  if (K_sel == 1L) {
    ### one regime: there is no partition and no penalty to choose
    phi_sel <- pooled$phi[1]
    step2 <- NULL
  } else {
    ### every admissible fit at the selected K, on the whole grid; ties
    ### resolved towards the smaller phi
    sub_k <- adm[adm$K == K_sel & is.finite(adm[[criterion]]), , drop = FALSE]
    sub_k <- sub_k[order(sub_k$phi), , drop = FALSE]
    best <- order(sub_k[[criterion]], sub_k$phi)[1]
    phi_sel <- sub_k$phi[best]
    step2 <- data.frame(phi = sub_k$phi, value = sub_k[[criterion]],
                        selected = seq_len(nrow(sub_k)) == best,
                        stringsAsFactors = FALSE)
  }

  ##############################################
  ########## Output ############################
  ##############################################

  tag_sel <- tag_of(K_sel, phi_sel)
  fit_sel <- infocrit$fits[[tag_sel]]

  sel_group <- infocrit$groups[, tag_sel]
  ari_tab <- do.call(rbind, lapply(colnames(infocrit$groups), function(tg) {
    data.frame(configuration = tg,
               ARI = scstem_ari(infocrit$groups[, tg], sel_group),
               stringsAsFactors = FALSE)
  }))

  out <- list(
    K_selected = K_sel,
    phi_selected = phi_sel,
    fit = fit_sel,
    step1 = step1,
    step2 = step2,
    ari_to_selected = ari_tab,
    reference = reference,
    selected_row = tab[tab$K == K_sel & tab$phi == phi_sel, , drop = FALSE],
    band = band,
    criterion = criterion
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
  cat("  criterion             : ", x$criterion, "\n\n", sep = "")

  cat("(S1) criterion-minimizing K within the band\n")
  s1 <- x$step1
  s1$value <- round(s1$value, digits)
  print(s1, row.names = FALSE)

  if (is.null(x$step2)) {
    cat("\n(S2) not applicable: no partition improves on the pooled model\n")
    cat("\nSelected configuration: K = 1 (the pooled model)\n")
  } else {
    cat("\n(S2) criterion over the penalties at K = ", x$K_selected, "\n", sep = "")
    s2 <- x$step2
    s2$value <- round(s2$value, digits)
    print(s2, row.names = FALSE)

    cat("\nSelected configuration: K = ", x$K_selected,
        " , phi = ", x$phi_selected, "\n", sep = "")
  }
  if (nrow(x$selected_row)) {
    sr <- x$selected_row
    cat("  loglik = ", round(sr$loglik, digits),
        " | ", x$criterion, " = ", round(sr[[x$criterion]], digits),
        " | cluster sizes >= ", sr$min_size, "\n", sep = "")
  }
  if (nrow(x$reference)) {
    cat("  pooled reference (K = 1): ", x$criterion, " = ",
        round(x$reference[[x$criterion]][1], digits), "\n", sep = "")
  }
  invisible(x)
}
