#' Computational settings of the STEM and SC-STEM estimators
#'
#' @description
#' \code{STEM_control} collects the settings that govern the numerical side of
#' the estimation -- when the EM algorithm stops, how many iterations it may
#' take, how long the alternation of SC-STEM runs and how hard the
#' Newton-Raphson step of the spatial parameters tries -- in one object, in the
#' spirit of the \code{control} argument of \code{\link[stats]{optim}}. The
#' defaults are meant for an analysis; an expert user may change them, for
#' instance to trade accuracy for speed in a simulation study.
#'
#' @details
#' \strong{When the EM algorithm stops.} Two criteria, both required at the
#' same iteration (or the iteration limit):
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
#' the log-likelihood depends on its scale, and with log-likelihoods of the
#' order of \eqn{10^4} a relative tolerance of \eqn{10^{-3}} stops the algorithm
#' while several units can still be gained. On a regime with a persistent
#' latent process the EM algorithm is slow along the ridge between the
#' intercept and the level of the latent process; the defaults reach the
#' maximum within a few thousandths of a log-likelihood unit.
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
#' @param em_tol_par tolerance of the relative change of every free parameter,
#'   for the pooled fit and the final refit. Default \code{1e-4}.
#' @param em_tol_loglik tolerance of the absolute change of the
#'   log-likelihood, for the pooled fit and the final refit. Default
#'   \code{1e-3}.
#' @param em_maxit maximum number of EM iterations of the pooled fit and of
#'   every regime in the final refit. Default \code{500}.
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
#' @param nr_hess_maxit maximum number of attempts at a negative-definite
#'   Hessian within one Newton-Raphson iteration. Default \code{30}.
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
#' ## a faster, looser setting, for a simulation study
#' STEM_control(em_tol_par = 1e-3, em_tol_loglik = 1e-2, em_maxit = 200)
#'
#' @export
STEM_control <- function(em_tol_par = 1e-4, em_tol_loglik = 1e-3, em_maxit = 500L,
                         alt_em_tol_par = 1e-2, alt_em_tol_loglik = 1, alt_em_maxit = 50L,
                         alt_maxit = 10L, alt_abs_tol = 1e-5, alt_rel_tol = 1e-6,
                         nr_maxit = 50L, nr_hess_maxit = 30L) {
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
  out <- list(em_tol_par = pos(em_tol_par, "em_tol_par"),
              em_tol_loglik = pos(em_tol_loglik, "em_tol_loglik"),
              em_maxit = int(em_maxit, "em_maxit"),
              alt_em_tol_par = pos(alt_em_tol_par, "alt_em_tol_par"),
              alt_em_tol_loglik = pos(alt_em_tol_loglik, "alt_em_tol_loglik"),
              alt_em_maxit = int(alt_em_maxit, "alt_em_maxit"),
              alt_maxit = nonneg(alt_maxit, "alt_maxit"),
              alt_abs_tol = pos(alt_abs_tol, "alt_abs_tol"),
              alt_rel_tol = pos(alt_rel_tol, "alt_rel_tol"),
              nr_maxit = int(nr_maxit, "nr_maxit"),
              nr_hess_maxit = int(nr_hess_maxit, "nr_hess_maxit"))
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
stem_control_em <- function(control, stage = c("final", "alternation")) {
  stage <- match.arg(stage)
  if (stage == "final") return(control)
  control$em_tol_par <- control$alt_em_tol_par
  control$em_tol_loglik <- control$alt_em_tol_loglik
  control$em_maxit <- control$alt_em_maxit
  control
}

#' @export
print.STEM_control <- function(x, ...) {
  cat("STEM estimation settings\n")
  cat("  final fits (pooled, final refit): tol_par = ", x$em_tol_par,
      ", tol_loglik = ", x$em_tol_loglik, ", maxit = ", x$em_maxit, "\n", sep = "")
  cat("  EM inside the alternation       : tol_par = ", x$alt_em_tol_par,
      ", tol_loglik = ", x$alt_em_tol_loglik, ", maxit = ", x$alt_em_maxit, "\n", sep = "")
  cat("  alternation                     : maxit = ", x$alt_maxit,
      ", abs_tol = ", x$alt_abs_tol, ", rel_tol = ", x$alt_rel_tol, "\n", sep = "")
  cat("  Newton-Raphson (spatial)        : maxit = ", x$nr_maxit,
      ", Hessian attempts = ", x$nr_hess_maxit, "\n", sep = "")
  invisible(x)
}
