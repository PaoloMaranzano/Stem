#' @keywords internal
#' @noRd

### ---------------------------------------------------------------------------
### The M-step of the EM algorithm
### ---------------------------------------------------------------------------
### The conditional updates of R/mstep-updates.R, in the order of Fasso and
### Cameletti (2007): the latent process (G, Sigmaeta and m0), sigma2omega,
### beta given the smoothed latent states, and the spatial parameters by
### Newton-Raphson. Every update maximizes the expected complete-data
### log-likelihood Q of the E-step `est`.
###
### Returns the new parameters, the number of Newton-Raphson iterations and
### what the update of beta reports about the penalty.
`stem_mstep_em` <- function(est, phi, dat, opt) {
  lat <- stem_update_latent(est, phi, dat, opt)
  s2w <- stem_update_sigma2omega(est, dat, opt)
  bet <- stem_update_beta(est, phi, dat, s2w, opt)
  spa <- stem_update_spatial(est, phi, dat, s2w, opt)

  phi$sigma2omega <- s2w
  phi$logtheta <- spa$logtheta
  phi$logb <- spa$logb
  phi$beta <- matrix(bet$beta, ncol = 1)
  phi$G <- lat$G
  phi$Sigmaeta <- lat$Sigmaeta
  phi$m0 <- lat$m0
  list(phi = phi, n_iter_NR = spa$n_iter,
       beta_df = bet$df, lambda_ref = bet$lambda_ref, lambda_eff = bet$lambda_eff)
}
