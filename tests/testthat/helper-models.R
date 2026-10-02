### Small, fast fixtures shared by the tests: a subset of the povalley network
### large enough for a three-regime fit but small enough to run in seconds.

po_subset <- function(Tn = 90L, d = NULL) {
  utils::data("povalley", package = "Stem", envir = environment())
  povalley <- get("povalley", envir = environment())
  Tfull <- nrow(povalley$z)
  if (is.null(d)) d <- ncol(povalley$z)
  idx <- seq_len(d)
  keep <- as.vector(outer(seq_len(Tn), (idx - 1L) * Tfull, "+"))
  list(z = povalley$z[seq_len(Tn), idx, drop = FALSE],
       covariates = povalley$covariates[keep, , drop = FALSE],
       coordinates = povalley$coords[idx, , drop = FALSE],
       Tn = Tn, d = d)
}

po_phi <- function() {
  list(beta = matrix(c(1.25, -0.00003, 0.64), 3, 1),
       sigma2eps = 18.66,
       sigma2omega = 1e-06,
       theta = 2e-06,
       G = matrix(0.59, 1, 1),
       Sigmaeta = matrix(4.25, 1, 1),
       m0 = as.matrix(0),
       C0 = as.matrix(1))
}

po_model <- function(Tn = 90L, d = NULL) {
  s <- po_subset(Tn = Tn, d = d)
  STEM_Model(z = s$z, covariates = s$covariates, coordinates = s$coordinates,
             phi = po_phi(), K = matrix(1, s$d, 1))
}

### The tests check behaviour, not accuracy: a loose, fast setting of the EM
### algorithm for the whole session, which the tests that set their own
### tolerance override.
options(Stem.control = list(em_tol_par = 0.05, em_tol_loglik = 1, em_maxit = 30,
                            alt_em_tol_par = 0.1, alt_em_tol_loglik = 10, alt_em_maxit = 10))
