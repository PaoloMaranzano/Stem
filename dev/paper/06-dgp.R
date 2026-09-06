## ---------------------------------------------------------------------------
## Data-generating processes for the simulation study. Definitions only: this
## file is sourced by 07-simulation.R and by the calibration script, and runs
## nothing of its own.
##
## WHY THIS REPLACES THE FIRST DESIGN
##
## The first version generated each regime by calling STEM_Simulation() on its
## own sub-network, which is what the SC-STEM model literally says: regime k has
## its own latent process y^(k), and two locations in different regimes are
## uncorrelated. That is faithful to the model and it is also useless as an
## experiment, because it makes the partition identifiable from the correlation
## structure alone, whatever the parameters. Generating with IDENTICAL
## parameters in the two regimes and one latent path per regime recovers the
## true partition with ARI = 1 in every replication; generating the same data
## with one shared latent path gives ARI = 0.14. The recovery measured in that
## design was therefore the recovery of the latent split, not of any difference
## in Psi_k, and a design that separates the regimes only in the dynamics or
## only in the spatial covariance was testing nothing of the sort.
##
## THE FIX
##
## The coupling between the latent processes becomes an explicit factor of the
## design. The innovations of the k processes are drawn with cross-correlation
## rho:
##
##   rho = 1  the regimes share their dynamics, so the partition can be found
##            ONLY through the parameters -- the honest, hard case
##   rho = 0  independent processes: the model-consistent, easy case, kept as a
##            reference so that the size of the confound can be quantified
##
## With equal transition coefficients the correlation of the latent processes is
## exactly rho; with different ones it is
## rho * sqrt((1-G1^2)(1-G2^2)) / (1 - G1 G2), which is below one however large
## rho is. That is not a defect of the design but a fact about the object: two
## AR(1) processes with different persistence cannot be perfectly correlated, so
## a scenario that separates the regimes on their dynamics necessarily leaks a
## little information through the latent paths. The design measures that leak
## instead of hiding it.
##
## The measurement error is drawn as ONE global spatial field whenever the
## covariance parameters are common to the regimes, and block by regime only
## when they are not -- in which case the block structure is intrinsic to the
## scenario rather than an artefact of the generator.
##
## Effect sizes are expressed in interpretable units: the coefficient contrast
## in residual standard deviations, the variance contrast as a nugget share, the
## range in kilometers.
## ---------------------------------------------------------------------------

## ---------------------------------------------------------------------------
## Baseline parameters, from the pooled fit on the real network
## ---------------------------------------------------------------------------
dgp_base <- function() {
  list(beta        = c(2.34, 0.60),   # intercept and one standardised covariate
       sigma2eps   = 21.4,            # nugget
       sigma2omega = 4.19,            # partial sill
       theta       = 8.1e-06,         # 1/theta = 123 km
       G           = 0.90,            # persistence of the latent process
       var_y       = 6.0)             # stationary variance of the latent process
}

## The innovation variance that holds the stationary variance of an AR(1) at
## var_y, so that a scenario separating the regimes on G separates them on
## persistence alone and not on the amplitude of the latent signal.
sigma_eta_of <- function(G, var_y) var_y * (1 - G^2)

## ---------------------------------------------------------------------------
## A standardised, exogenous covariate with realistic structure: an AR(1) in
## time whose innovations are a spatially correlated field. A covariate that was
## independent across stations would make a coefficient contrast trivially
## visible; one that was perfectly common would make it indistinguishable from
## the latent process. This sits in between, and the two knobs are explicit.
## ---------------------------------------------------------------------------
dgp_covariate <- function(coords, TN, a_time = 0.7, range_km = 100, seed = 1) {
  d <- nrow(coords)
  dm <- geodist::geodist(coords, measure = "geodesic")
  Cx <- exp(-dm / (range_km * 1000))
  L <- chol(Cx + diag(1e-8, d))

  set.seed(seed)
  x <- matrix(NA_real_, TN, d)
  w <- as.numeric(crossprod(L, stats::rnorm(d)))
  x[1, ] <- w
  for (tt in 2:TN) {
    w <- as.numeric(crossprod(L, stats::rnorm(d)))
    x[tt, ] <- a_time * x[tt - 1, ] + sqrt(1 - a_time^2) * w
  }
  (x - mean(x)) / stats::sd(x)
}

## the covariate stacked by station, as STEM_Model() expects
dgp_design <- function(x) {
  TN <- nrow(x); d <- ncol(x)
  cbind(intercept = 1, xcov = as.vector(x))
}

## ---------------------------------------------------------------------------
## Parameter sets of the two regimes
##
## `scenario` says which component separates them, `level` how far apart they
## are. Regime 1 always carries the baseline.
## ---------------------------------------------------------------------------
dgp_psi <- function(scenario, level = 2L, base = dgp_base()) {

  A <- base; B <- base
  res_sd <- sqrt(base$sigma2eps + base$sigma2omega)

  ## coefficient contrast, in residual standard deviations. The covariate is
  ## standardised, so a contrast of delta*res_sd in beta moves the conditional
  ## mean by delta residual standard deviations per unit of the covariate.
  delta_beta <- c(0.25, 0.5, 1.0)[level]
  ## persistence of the second regime
  G2         <- c(0.80, 0.60, 0.30)[level]
  ## nugget share of the second regime, against 0.836 in the baseline
  share2     <- c(0.70, 0.50, 0.30)[level]
  ## range of the second regime, in km, against 123 km in the baseline
  range2     <- c(60, 30, 15)[level]

  if (scenario %in% c("S1", "S4")) {
    B$beta <- c(base$beta[1], base$beta[2] + delta_beta * res_sd)
  }
  if (scenario %in% c("S2", "S4")) {
    B$G <- G2
  }
  if (scenario %in% c("S3", "S4")) {
    tot <- base$sigma2eps + base$sigma2omega
    B$sigma2eps   <- tot * share2
    B$sigma2omega <- tot * (1 - share2)
    B$theta       <- 1 / (range2 * 1000)
  }
  list(A, B)
}

## TRUE when the two regimes share every parameter of the measurement
## covariance, in which case the error field can be drawn globally
dgp_common_field <- function(psi) {
  isTRUE(all.equal(psi[[1]][c("sigma2eps", "sigma2omega", "theta")],
                   psi[[2]][c("sigma2eps", "sigma2omega", "theta")]))
}

## ---------------------------------------------------------------------------
## The generator
##
##   labels  regime of each location
##   psi     list of parameter sets, one per regime
##   rho     cross-correlation of the innovations of the latent processes
##   force_block  draw the measurement error block by regime even when the
##                covariance parameters are common: this reproduces the first
##                design and is what makes the confound measurable
## ---------------------------------------------------------------------------
dgp_simulate <- function(labels, psi, x, coords, rho = 1, seed = 1,
                         force_block = FALSE) {

  TN <- nrow(x); d <- ncol(x)
  k  <- length(psi)
  dm <- geodist::geodist(coords, measure = "geodesic")

  set.seed(seed)

  ## ---- latent processes, with innovations correlated at rho ----------------
  Gs   <- vapply(psi, function(p) p$G, numeric(1))
  vy   <- psi[[1]]$var_y
  sds  <- sqrt(vapply(Gs, sigma_eta_of, numeric(1), var_y = vy))
  R    <- matrix(rho, k, k); diag(R) <- 1
  Seta <- diag(sds, k) %*% R %*% diag(sds, k)
  Leta <- chol(Seta + diag(1e-10, k))

  y <- matrix(NA_real_, TN, k)
  y[1, ] <- sqrt(vy) * as.numeric(crossprod(chol(R + diag(1e-10, k)),
                                            stats::rnorm(k)))
  for (tt in 2:TN) {
    eta <- as.numeric(crossprod(Leta, stats::rnorm(k)))
    y[tt, ] <- Gs * y[tt - 1, ] + eta
  }

  ## ---- measurement error ---------------------------------------------------
  e <- matrix(NA_real_, TN, d)
  if (dgp_common_field(psi) && !force_block) {
    ## one global spatial field: the partition leaves no trace in the residual
    ## covariance, which is what makes S0 a genuine null
    p1 <- psi[[1]]
    Sig <- p1$sigma2eps * diag(d) + p1$sigma2omega * exp(-p1$theta * dm)
    L <- chol(Sig)
    for (tt in seq_len(TN)) e[tt, ] <- as.numeric(crossprod(L, stats::rnorm(d)))
  } else {
    for (g in seq_len(k)) {
      idx <- which(labels == g)
      p <- psi[[g]]
      Sig <- p$sigma2eps * diag(length(idx)) +
             p$sigma2omega * exp(-p$theta * dm[idx, idx, drop = FALSE])
      L <- chol(Sig)
      for (tt in seq_len(TN)) {
        e[tt, idx] <- as.numeric(crossprod(L, stats::rnorm(length(idx))))
      }
    }
  }

  ## ---- assemble ------------------------------------------------------------
  z <- matrix(NA_real_, TN, d)
  for (i in seq_len(d)) {
    g <- labels[i]
    z[, i] <- psi[[g]]$beta[1] + psi[[g]]$beta[2] * x[, i] + y[, g] + e[, i]
  }
  attr(z, "latent") <- y
  z
}

## ---------------------------------------------------------------------------
## The scenarios
##
## rho = 1 throughout except for the reference cell S5, which reproduces the
## model-consistent generator of the first design so that the contribution of
## the latent split can be read off directly.
## ---------------------------------------------------------------------------
dgp_scenarios <- function() {
  rbind(
    data.frame(id = "S0", scenario = "S0", level = 2L, rho = 1, force_block = FALSE,
               label = "null: no regime at all"),
    data.frame(id = "S1a", scenario = "S1", level = 1L, rho = 1, force_block = FALSE,
               label = "coefficients, 0.25 residual sd"),
    data.frame(id = "S1b", scenario = "S1", level = 2L, rho = 1, force_block = FALSE,
               label = "coefficients, 0.50 residual sd"),
    data.frame(id = "S1c", scenario = "S1", level = 3L, rho = 1, force_block = FALSE,
               label = "coefficients, 1.00 residual sd"),
    data.frame(id = "S2a", scenario = "S2", level = 1L, rho = 1, force_block = FALSE,
               label = "persistence, G = 0.90 vs 0.80"),
    data.frame(id = "S2b", scenario = "S2", level = 3L, rho = 1, force_block = FALSE,
               label = "persistence, G = 0.90 vs 0.30"),
    data.frame(id = "S3a", scenario = "S3", level = 1L, rho = 1, force_block = FALSE,
               label = "covariance, nugget share .84 vs .70, range 123 vs 60 km"),
    data.frame(id = "S3b", scenario = "S3", level = 3L, rho = 1, force_block = FALSE,
               label = "covariance, nugget share .84 vs .30, range 123 vs 15 km"),
    data.frame(id = "S4", scenario = "S4", level = 2L, rho = 1, force_block = FALSE,
               label = "all three contrasts, middle level"),
    data.frame(id = "S5", scenario = "S4", level = 2L, rho = 0, force_block = TRUE,
               label = "reference: the model-consistent generator"),
    stringsAsFactors = FALSE
  )
}
