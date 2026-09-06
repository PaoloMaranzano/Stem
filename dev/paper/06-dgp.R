## ---------------------------------------------------------------------------
## Data-generating processes for the simulation study. Definitions only: this
## file is sourced by 07-simulation.R and runs nothing of its own.
##
## THE GEOMETRY: THE OVERLAP DESIGN
##
## The spatial configuration follows Morelli, Maranzano and Otto (2026), Spatial
## Statistics 73, 100960, Section 4. There the K = 4 cluster centres sit at the
## corners of a square of half-side d,
##
##   mu_sp = ((d,d), (-d,d), (d,-d), (-d,-d)),   Sigma_sp = nu_sp * I_2,
##
## so that d alone controls how much the clusters overlap in space: at d = 0 the
## four Gaussians coincide and the partition has no spatial signature at all; as
## d grows the clusters separate. Reducing the design to K = 2 and K = 3 we keep
## the NEAREST-NEIGHBOUR centre distance at 2d rather than the radius, so that a
## given d means the same degree of overlap whatever K:
##
##   K = 2   centres (-d, 0) and (d, 0)
##   K = 3   equilateral triangle of side 2d, i.e. circumradius 2d/sqrt(3)
##
## The standardised separation is then 2d / sqrt(nu_sp) for every K: with
## nu_sp = 0.4 it runs from 0 at d = 0 to 3.16 standard deviations at d = 1.
##
## The abstract plane is mapped onto a geographic box centred on the Po Valley,
## one abstract unit being UNIT_KM kilometres. This is not cosmetic: it is what
## makes the covariance parameters interpretable. With UNIT_KM = 100 and
## nu_sp = 0.4 a cluster has a standard deviation of 63 km and, at d = 1, the
## centres are 200 km apart, against a baseline correlation range of 123 km --
## the regime in which a monitoring network of the Po Valley actually sits.
##
## WHY THE GENERATOR IS NOT STEM_Simulation() CALLED REGIME BY REGIME
##
## The first version of this file generated each regime by calling
## STEM_Simulation() on its own sub-network, which is what the SC-STEM model
## literally says: regime k has its own latent process y^(k), and two locations
## in different regimes are uncorrelated. That is faithful to the model and it
## is also useless as an experiment, because it makes the partition identifiable
## from the correlation structure alone, whatever the parameters. Generating
## with IDENTICAL parameters in the two regimes and one latent path per regime
## recovers the true partition with ARI = 1 in every replication; generating the
## same data with one shared latent path gives ARI = 0.14. The recovery measured
## in that design was the recovery of the latent split, not of any difference in
## Psi_k, and a design that separates the regimes only in the dynamics or only
## in the spatial covariance was testing nothing of the sort.
##
## THE FIX
##
## The coupling between the latent processes is an explicit factor of the
## design. The innovations of the K processes are drawn with cross-correlation
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
## Effect sizes are in interpretable units: the coefficient contrast in residual
## standard deviations, the variance contrast as a nugget share, the range in
## kilometres.
## ---------------------------------------------------------------------------

NU_SP   <- 0.4          # variance of each coordinate within a cluster
UNIT_KM <- 100          # kilometres per abstract unit of the overlap design
LON0    <- 9.5          # centre of the geographic box, Po Valley
LAT0    <- 45.5

## ---------------------------------------------------------------------------
## Geometry
## ---------------------------------------------------------------------------

## centres of the K clusters at overlap d, nearest-neighbour distance 2d
dgp_centres <- function(K, d) {
  if (K == 1L) return(cbind(0, 0))
  if (K == 2L) return(cbind(c(-d, d), c(0, 0)))
  if (K == 3L) {
    r <- 2 * d / sqrt(3); a <- c(90, 210, 330) * pi / 180
    return(cbind(r * cos(a), r * sin(a)))
  }
  if (K == 4L) return(cbind(c(d, -d, d, -d), c(d, d, -d, -d)))
  stop("K must be 1, 2, 3 or 4")
}

## Labels. `balanced = TRUE` gives exactly equal regime sizes, the remainder
## spread over the first regimes; `balanced = FALSE` draws them with equal
## probabilities, as in the source design. The balanced version is the default
## because SC-STEM carries a minimum regime size, and an unlucky multinomial
## draw at n = 20 produces a regime the model cannot fit -- which would confound
## the recovery of the partition with the feasibility of the fit.
dgp_labels <- function(n, K, balanced = TRUE) {
  if (balanced) {
    g <- rep(seq_len(K), length.out = n)
    return(sort(g))
  }
  sample.int(K, n, replace = TRUE)
}

## Locations: the abstract cloud of the overlap design, mapped to longitude and
## latitude so that distances are kilometres and the covariance parameters keep
## their meaning.
dgp_locations <- function(n, K, d, nu_sp = NU_SP, balanced = TRUE, seed = 1) {
  set.seed(seed)
  mu <- dgp_centres(K, d)
  g  <- dgp_labels(n, K, balanced)
  xy <- cbind(mu[g, 1] + stats::rnorm(n, sd = sqrt(nu_sp)),
              mu[g, 2] + stats::rnorm(n, sd = sqrt(nu_sp)))
  coords <- cbind(
    lon = LON0 + xy[, 1] * UNIT_KM / (111.320 * cos(LAT0 * pi / 180)),
    lat = LAT0 + xy[, 2] * UNIT_KM / 110.574)
  list(coords = coords, labels = g, xy = xy, mu = mu)
}

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
  cbind(intercept = 1, xcov = as.vector(x))
}

## ---------------------------------------------------------------------------
## Parameter sets of the K regimes
##
## `scenario` says which component separates them, `level` how far apart the two
## EXTREME regimes are. Regime 1 always carries the baseline and regime K the
## full contrast; with K = 3 the middle regime sits halfway, so that K = 2
## reproduces exactly the earlier design and K = 3 extends it without changing
## the meaning of `level`.
## ---------------------------------------------------------------------------
dgp_psi <- function(scenario, level = 2L, K = 2L, base = dgp_base()) {

  res_sd <- sqrt(base$sigma2eps + base$sigma2omega)
  tot    <- base$sigma2eps + base$sigma2omega

  ## the contrast of the extreme regime
  delta_beta <- c(0.25, 0.5, 1.0)[level]   # in residual standard deviations
  G_last     <- c(0.80, 0.60, 0.30)[level] # against 0.90 in the baseline
  share_last <- c(0.70, 0.50, 0.30)[level] # nugget share, against 0.836
  range_last <- c(60, 30, 15)[level]       # km, against 123 km

  ## regime g sits at fraction w of the way from the baseline to the extreme
  w_of <- function(g) if (K == 1L) 0 else (g - 1) / (K - 1)

  lapply(seq_len(K), function(g) {
    w <- w_of(g)
    p <- base
    if (scenario %in% c("S1", "S4")) {
      p$beta <- c(base$beta[1], base$beta[2] + w * delta_beta * res_sd)
    }
    if (scenario %in% c("S2", "S4")) {
      p$G <- base$G + w * (G_last - base$G)
    }
    if (scenario %in% c("S3", "S4")) {
      share <- (base$sigma2eps / tot) + w * (share_last - base$sigma2eps / tot)
      rng   <- (1 / base$theta / 1000) + w * (range_last - 1 / base$theta / 1000)
      p$sigma2eps   <- tot * share
      p$sigma2omega <- tot * (1 - share)
      p$theta       <- 1 / (rng * 1000)
    }
    p
  })
}

## TRUE when all regimes share every parameter of the measurement covariance, in
## which case the error field can be drawn globally
dgp_common_field <- function(psi) {
  key <- function(p) unlist(p[c("sigma2eps", "sigma2omega", "theta")])
  all(vapply(psi[-1], function(p) isTRUE(all.equal(key(p), key(psi[[1]]))),
             logical(1)))
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
## The scenarios: what separates the regimes, and by how much
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

## ---------------------------------------------------------------------------
## The dimensions of the design
##
## T is read as a real observation window: 60 and 120 are five and ten years of
## monthly data, 365 and 730 one and two years of daily data. n spans the sizes
## a regional network actually takes: Northern Italy carries about 260 air
## quality stations, so 20 and 40 are a small sub-network, 60 and 100 a regional
## one, 200 and 400 a national or multi-regional one.
## ---------------------------------------------------------------------------
dgp_dims <- function() {
  list(TN = c(60L, 120L, 365L, 730L),
       n  = c(20L, 40L, 60L, 100L, 200L, 400L),
       K  = c(2L, 3L),
       d  = c(0, 1/3, 2/3, 1))
}

## ---------------------------------------------------------------------------
## One complete data set of the design, ready for STEM_Model()
## ---------------------------------------------------------------------------
dgp_draw <- function(n, TN, K, d, scenario_row, rep = 1L, balanced = TRUE) {
  seed <- 1000L * rep + 1L
  loc  <- dgp_locations(n, K, d, balanced = balanced, seed = seed)
  x    <- dgp_covariate(loc$coords, TN, seed = seed + 1L)
  psi  <- dgp_psi(scenario_row$scenario, scenario_row$level, K = K)
  z    <- dgp_simulate(loc$labels, psi, x, loc$coords,
                       rho = scenario_row$rho, seed = seed + 2L,
                       force_block = scenario_row$force_block)
  list(z = z, covariates = dgp_design(x), coordinates = loc$coords,
       labels = loc$labels, psi = psi)
}
