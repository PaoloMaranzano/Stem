#' Computational settings of the STEM and SC-STEM estimators
#'
#' @description
#' \code{STEM_control} collects the settings that govern the numerical side of
#' the estimation -- which algorithm maximizes the likelihood, when it stops,
#' how many iterations it may take, how long the alternation of SC-STEM runs
#' and how hard the Newton-Raphson step of the spatial parameters tries -- in
#' one object, in the spirit of the \code{control} argument of
#' \code{\link[stats]{optim}}. The defaults are meant for an analysis; an
#' expert user may change them, for instance to trade accuracy for speed in a
#' simulation study.
#'
#' @details
#' \strong{The algorithm.} \code{"EM"}, \code{"ECME"}, \code{"SQUAREM"} (the
#' default) or \code{"SQUAREM-ECME"}; see the section "Algorithms" of
#' \code{\link{STEM_Estimation}}. All four maximize the same likelihood. The EM
#' algorithm converges linearly and slowly on the regimes of SC-STEM, where the
#' latent process and the spatial field share most of their variation; ECME
#' updates the regression coefficients and \eqn{m_0} on the observed
#' likelihood, which makes their estimates accurate at the stop; SQUAREM
#' extrapolates the iterations and cuts their number.
#'
#' \strong{When the EM algorithm stops.} Two criteria; with
#' \code{em_stop = "any"} (the default) the algorithm stops as soon as either is
#' met, with \code{em_stop = "all"} only when both are met at the same
#' iteration, and in either case after \code{em_maxit} iterations:
#' \deqn{\max_l \frac{|\psi_l^{(i)} - \psi_l^{(i-1)}|}{\max(|\psi_l^{(i-1)}|,
#'   10^{-3})} < \texttt{tol\_par}, \qquad
#'   |\ell(\psi^{(i)}) - \ell(\psi^{(i-1)})| < \texttt{tol\_loglik}.}{
#'   max_l |psi_l(i) - psi_l(i-1)| / max(|psi_l(i-1)|, 1e-3) < tol_par  and
#'   |loglik(i) - loglik(i-1)| < tol_loglik.}
#' The first is the relative change of every free parameter taken one at a time
#' (the criterion of D-STEM, Wang, Finazzi and Fasso, 2021), so that a small
#' parameter still moving -- the persistence of the latent process, for
#' instance -- is not hidden by the large ones, as it is in the relative change
#' of the whole parameter vector. The floor \eqn{10^{-3}} keeps the criterion
#' finite for parameters near zero. The loading matrix and the initial variance
#' \eqn{C_0}, which the EM algorithm does not estimate, are excluded. The second
#' criterion is the absolute change of the log-likelihood: a relative change of
#' the log-likelihood depends on its scale, which grows with the size of the
#' data, while the differences that matter -- between fits, in the information
#' criteria -- are absolute. D-STEM stops when the relative change of the
#' parameters or of the log-likelihood is below \eqn{10^{-4}}: with
#' log-likelihoods of the order of \eqn{10^4} the second is met first, while
#' several units can still be gained. In a simulation study of SC-STEM,
#' tolerances of \eqn{10^{-3}} on both criteria left the estimates within 1 to
#' 2\% of those at \eqn{10^{-4}} and the selected models unchanged, while
#' \eqn{10^{-2}} moved the range and the variances by tens of percent.
#'
#' \strong{Session defaults.} Every estimation function takes
#' \code{control = NULL} by default, which means the defaults of
#' \code{STEM_control()}, or, when the option is set,
#' \code{options(Stem.control = list(...))}: a list of some of the settings,
#' which then apply to every fit of the session that does not pass its own
#' \code{control}.
#'
#' \strong{Two sets of EM settings.} The pooled model and the final refit of
#' SC-STEM -- the fits whose estimates and log-likelihoods are reported and
#' enter the information criteria -- use \code{em_tol_par},
#' \code{em_tol_loglik} and \code{em_maxit}. Inside the alternation of SC-STEM
#' the EM algorithm only has to rank the regimes, and every regime starts from
#' its estimates of the previous iteration, so looser settings suffice there:
#' \code{alt_em_tol_par}, \code{alt_em_tol_loglik} and \code{alt_em_maxit}.
#'
#' @param algorithm the algorithm that maximizes the likelihood:
#'   \code{"SQUAREM"} (default), \code{"EM"}, \code{"ECME"} or
#'   \code{"SQUAREM-ECME"}. It applies to the pooled fit and to the final
#'   refits of SC-STEM; inside the alternation, whose fits stop after a few
#'   iterations, the iterations of EM or ECME are not accelerated.
#' @param em_tol_par tolerance of the relative change of every free parameter,
#'   for the pooled fit and the final refit. Default \code{1e-3}.
#' @param em_tol_loglik tolerance of the absolute change of the
#'   log-likelihood, for the pooled fit and the final refit. Default
#'   \code{1e-3}.
#' @param em_maxit maximum number of EM iterations of the pooled fit and of
#'   every regime in the final refit. Default \code{500}.
#' @param em_stop how the two criteria combine, for every EM run:
#'   \code{"any"} (default), stop as soon as either is met; \code{"all"}, stop
#'   only when both are met.
#' @param alt_em_tol_par,alt_em_tol_loglik,alt_em_maxit the same three settings
#'   for the EM algorithm run inside the alternation of SC-STEM. Defaults
#'   \code{1e-2}, \code{1} and \code{50}.
#' @param alt_maxit maximum number of iterations of the alternation of
#'   SC-STEM. Default \code{10}. Zero refits the regimes of a given starting
#'   partition without moving it.
#' @param alt_abs_tol,alt_rel_tol absolute and relative tolerance on the
#'   penalized objective, below which the alternation stops. Defaults
#'   \code{1e-5} and \code{1e-6}.
#' @param nr_maxit maximum number of Newton-Raphson iterations of the spatial
#'   parameters within one EM iteration. Default \code{50}.
#'
#' @return An object of class \dQuote{STEM_control}: a named list of the
#'   settings.
#'
#' @references
#' Wang, Y., Finazzi, F., Fasso, A. (2021) \emph{D-STEM v2: A Software for
#' Modeling Functional Spatio-Temporal Data}. Journal of Statistical Software,
#' 99(10), 1--29. \doi{10.18637/jss.v099.i10}
#'
#' @seealso \code{\link{STEM_Estimation}}, \code{\link{SCSTEM_Estimation}}
#'
#' @examples
#' STEM_control()
#' ## the plain EM algorithm, with a tighter stopping rule
#' STEM_control(algorithm = "EM", em_tol_par = 1e-4, em_tol_loglik = 1e-4)
#'
#' @export
STEM_control <- function(algorithm = c("SQUAREM", "EM", "ECME", "SQUAREM-ECME"),
                         em_tol_par = 1e-3, em_tol_loglik = 1e-3, em_maxit = 500L,
                         em_stop = c("any", "all"),
                         alt_em_tol_par = 1e-2, alt_em_tol_loglik = 1, alt_em_maxit = 50L,
                         alt_maxit = 10L, alt_abs_tol = 1e-5, alt_rel_tol = 1e-6,
                         nr_maxit = 50L) {
  pos <- function(x, nm) {
    if (length(x) != 1L || !is.numeric(x) || is.na(x) || x <= 0) {
      stop("'", nm, "' must be a single positive number.", call. = FALSE)
    }
    x
  }
  int <- function(x, nm) as.integer(pos(round(x), nm))
  nonneg <- function(x, nm) {
    if (length(x) != 1L || !is.numeric(x) || is.na(x) || x < 0) {
      stop("'", nm, "' must be a single non-negative number.", call. = FALSE)
    }
    as.integer(round(x))
  }
  algorithm <- match.arg(algorithm)
  em_stop <- match.arg(em_stop)
  out <- list(algorithm = algorithm,
              em_tol_par = pos(em_tol_par, "em_tol_par"),
              em_tol_loglik = pos(em_tol_loglik, "em_tol_loglik"),
              em_maxit = int(em_maxit, "em_maxit"),
              em_stop = em_stop,
              alt_em_tol_par = pos(alt_em_tol_par, "alt_em_tol_par"),
              alt_em_tol_loglik = pos(alt_em_tol_loglik, "alt_em_tol_loglik"),
              alt_em_maxit = int(alt_em_maxit, "alt_em_maxit"),
              alt_maxit = nonneg(alt_maxit, "alt_maxit"),
              alt_abs_tol = pos(alt_abs_tol, "alt_abs_tol"),
              alt_rel_tol = pos(alt_rel_tol, "alt_rel_tol"),
              nr_maxit = int(nr_maxit, "nr_maxit"))
  class(out) <- c("STEM_control", "list")
  out
}

## A control given as NULL, as an object of STEM_control(), or as a plain named
## list of some of its settings (the others at their defaults). NULL means the
## session defaults: those of options(Stem.control = list(...)) if set, the
## defaults of STEM_control() otherwise.
stem_control_resolve <- function(control) {
  if (is.null(control)) control <- getOption("Stem.control")
  if (is.null(control)) return(STEM_control())
  if (inherits(control, "STEM_control")) return(control)
  if (!is.list(control)) {
    stop("'control' must be NULL, a list, or the output of STEM_control().", call. = FALSE)
  }
  bad <- setdiff(names(control), names(formals(STEM_control)))
  if (length(bad)) {
    stop("unknown setting(s) in 'control': ", paste(bad, collapse = ", "),
         "; see ?STEM_control.", call. = FALSE)
  }
  do.call(STEM_control, control)
}

## The EM settings a fit runs with: those of the final fits, or those of the
## alternation, written into the em_* fields that STEM_Estimation() reads.
## Inside the alternation the iterations are never accelerated: its fits stop
## after a few iterations, where the cycles of SQUAREM cost more than they
## save, and its extrapolations from points that are far from converged were
## seen to move the partitions the alternation reaches.
stem_control_em <- function(control, stage = c("final", "alternation")) {
  stage <- match.arg(stage)
  if (stage == "final") return(control)
  control$em_tol_par <- control$alt_em_tol_par
  control$em_tol_loglik <- control$alt_em_tol_loglik
  control$em_maxit <- control$alt_em_maxit
  control$algorithm <- switch(control$algorithm, "SQUAREM" = "EM", "SQUAREM-ECME" = "ECME",
                              control$algorithm)
  control
}

#' @export
print.STEM_control <- function(x, ...) {
  cat("STEM estimation settings\n")
  cat("  algorithm                       : ", x$algorithm, "\n", sep = "")
  cat("  the iterations stop when        : ",
      if (x$em_stop == "any") "either criterion is met" else "both criteria are met",
      " (em_stop = \"", x$em_stop, "\")\n", sep = "")
  cat("  final fits (pooled, final refit): tol_par = ", x$em_tol_par,
      ", tol_loglik = ", x$em_tol_loglik, ", maxit = ", x$em_maxit, "\n", sep = "")
  cat("  EM inside the alternation       : tol_par = ", x$alt_em_tol_par,
      ", tol_loglik = ", x$alt_em_tol_loglik, ", maxit = ", x$alt_em_maxit, "\n", sep = "")
  cat("  alternation                     : maxit = ", x$alt_maxit,
      ", abs_tol = ", x$alt_abs_tol, ", rel_tol = ", x$alt_rel_tol, "\n", sep = "")
  cat("  Newton-Raphson (spatial)        : maxit = ", x$nr_maxit, "\n", sep = "")
  invisible(x)
}
