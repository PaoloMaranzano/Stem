## ---------------------------------------------------------------------------
## Data-generating processes for the simulation study. Definitions only: this
## file is sourced by 07-simulation.R and runs nothing of its own.
##
## THE GEOMETRY: THE OVERLAP DESIGN
##
## The spatial configuration follows Morelli, Maranzano and Otto (2026), Spatial
## Statistics 73, 100960, Section 4. THE OVERLAP PARAMETER IS THEIR d, WHICH IS
## CALLED omega THROUGHOUT THIS DESIGN: in the package "d" already means the
## number of locations, and carrying two meanings for one letter through the
## simulation code and the paper was a defect waiting to happen.
##
## In the source design the K = 4 cluster centres sit at the corners of a square
## of half-side omega,
##
##   mu_sp = ((w,w), (-w,w), (w,-w), (-w,-w)),   Sigma_sp = nu_sp * I_2,
##
## so that omega alone controls how much the clusters overlap in space: at
## omega = 0 the Gaussians coincide and the partition has no spatial signature
## at all; as omega grows the clusters separate. Reducing the design to K = 2
## and K = 3 we keep the NEAREST-NEIGHBOUR centre distance at 2*omega rather
## than the radius, so that a given omega means the same degree of overlap
## whatever K:
##
##   K = 2   centres (-omega, 0) and (omega, 0)
##   K = 3   equilateral triangle of side 2*omega, circumradius 2*omega/sqrt(3)
##   K = 4   the square of the source paper, radius omega*sqrt(2)
##
## The standardised separation is then 2*omega / sqrt(nu_sp) for every K: with
## nu_sp = 0.4 it runs from 0 at omega = 0 to 3.16 standard deviations at
## omega = 1.
##
## The abstract plane is mapped onto a geographic box centred on the Po Valley,
## one abstract unit being UNIT_KM kilometres. This is not cosmetic: it is what
## makes the covariance parameters interpretable. With UNIT_KM = 100 and
## nu_sp = 0.4 a cluster has a standard deviation of 63 km and, at omega = 1,
## the centres are 200 km apart, against a baseline correlation range of 123 km
## -- the regime in which a monitoring network of the Po Valley actually sits.
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

## ---------------------------------------------------------------------------
## THE GENERATOR, IN FULL
##
## Write K for the number of regimes, n for the number of locations, T for the
## number of time points, g_i in {1,...,K} for the regime of location i and
## I_g = {i : g_i = g} for its index set.
##
## (1) Regime sizes.  n_1,...,n_K from dgp_sizes(): equal up to the remainder
##     under "balanced", proportional to 1:2:...:K under "unbalanced", in both
##     cases with every n_g >= N_MIN. Labels are assigned in blocks, so the
##     partition is fixed within a cell and only the coordinates are random.
##
## (2) Coordinates.  In the abstract plane of the overlap design,
##
##        s_i | g_i = g  ~  N_2( mu_g , nu_sp I_2 ) ,        independent over i,
##
##     with mu_g = dgp_centres(K, omega): centres at nearest-neighbour distance
##     2*omega, so 2*omega / sqrt(nu_sp) is the standardised separation and is
##     the same for every K. For K = 1, nu_sp is replaced by
##     nu_sp + Var(mu | K = 2, omega). The plane is mapped affinely to longitude
##     and latitude at UNIT_KM kilometres per unit, centred on (LON0, LAT0), and
##     h_ij is the geodesic distance between s_i and s_j.
##
## (3) Covariate.  One standardised, exogenous covariate, an AR(1) in time whose
##     innovations are a spatially correlated field:
##
##        x_1 = w_1 ,   x_t = a x_{t-1} + sqrt(1 - a^2) w_t ,
##        w_t ~ N_n(0, C_x) ,   (C_x)_ij = exp(-h_ij / rho_x) ,
##
##     with a = 0.7 and rho_x = 100 km, then centred and scaled over all nT
##     values. A covariate independent across stations would make a coefficient
##     contrast trivially visible; a perfectly common one would make it
##     indistinguishable from the latent process. This sits in between.
##
## (4) Latent processes.  One AR(1) per regime, with cross-correlated
##     innovations:
##
##        y_t = diag(G_1,...,G_K) y_{t-1} + eta_t ,   eta_t ~ N_K(0, Sigma_eta),
##        Sigma_eta = D R_rho D ,  D = diag(sigma_eta,1 , ... , sigma_eta,K) ,
##        (R_rho)_gh = rho for g != h and 1 for g = h ,
##        sigma^2_eta,g = vbar_y (1 - G_g^2) ,   y_1 ~ N_K(0, vbar_y R_rho) .
##
##     The innovation variance is tied to G so that every regime has the same
##     stationary variance vbar_y: a scenario that separates the regimes on G
##     then separates them on persistence alone, not on the amplitude of the
##     signal.
##
## (5) Measurement error.  When all regimes share (sigma^2_eps, sigma^2_omega,
##     theta) and force_block is FALSE, ONE global field:
##
##        e_t ~ N_n(0, Sigma) ,  Sigma = sigma^2_eps I_n + sigma^2_omega exp(-theta h),
##
##     independent over t. Otherwise one field per regime, independent across
##     regimes:
##
##        e_{t,I_g} ~ N_{n_g}(0, Sigma_g) ,
##        Sigma_g = sigma^2_eps,g I_{n_g} + sigma^2_omega,g exp(-theta_g h_{I_g,I_g}) .
##
## (6) Response.
##
##        z_ti = beta_0,g_i + beta_1,g_i x_ti + y_t,g_i + e_ti .
##
##     This is the STEM measurement equation with loading matrix K_g = 1_{n_g}
##     within each regime, which is what SCSTEM_Estimation() fits.
##
## PSEUDOCODE
##
##   input  n, T, K, omega, balance, scenario s, level l, rho, replication r
##   1  n_1..n_K  <- sizes(n, K, balance)          ; g <- labels(n_1..n_K)
##   2  mu        <- centres(K, omega)
##   3  for i in 1..n:  s_i <- mu_{g_i} + N_2(0, nu_sp I)      ; map to lon/lat
##   4  h         <- geodesic distances between the s_i
##   5  C_x       <- exp(-h / rho_x)  ;  L_x <- chol(C_x)
##      x_1 <- L_x' N(0, I) ;  for t in 2..T:  x_t <- a x_{t-1} + sqrt(1-a^2) L_x' N(0,I)
##      x   <- (x - mean(x)) / sd(x)
##   6  Psi_1..Psi_K <- parameters(s, l, K)        ; regime g at fraction (g-1)/(K-1)
##   7  Sigma_eta <- D R_rho D  ;  L_eta <- chol(Sigma_eta)
##      y_1 <- sqrt(vbar_y) chol(R_rho)' N(0,I)
##      for t in 2..T:  y_t <- diag(G) y_{t-1} + L_eta' N(0, I)
##   8  if the regimes share the covariance and not force_block:
##         L <- chol(Sigma) ;  for t in 1..T:  e_t <- L' N(0, I)
##      else for g in 1..K:
##         L_g <- chol(Sigma_g) ;  for t in 1..T:  e_{t,I_g} <- L_g' N(0, I)
##   9  for t in 1..T, i in 1..n:
##         z_ti <- beta_0,g_i + beta_1,g_i x_ti + y_t,g_i + e_ti
##   output z (T x n), X = [1, vec(x)], coordinates, g, Psi_1..Psi_K
##
## Seeds are derived from the replication index, so a cell is reproducible on
## its own and the same replication uses the same geometry across scenarios.
## ---------------------------------------------------------------------------

NU_SP   <- 0.4          # variance of each coordinate within a cluster
UNIT_KM <- 100          # kilometres per abstract unit of the overlap design
LON0    <- 9.5          # centre of the geographic box, Po Valley
LAT0    <- 45.5
## The smallest regime the design allows. The package itself refuses fewer than
## r + 2 units, but the binding constraint is spatial rather than parametric: a
## regime estimates a range from its own pairwise distances, and six locations
## already give fifteen of them. Below that the range is not identified in any
## practical sense and a failure to recover the partition would be a failure to
## fit, not a failure to separate.
N_MIN   <- 6L
## The smallest ratio of largest to smallest regime that still counts as an
## unbalanced design. Below it the correction for N_MIN has flattened the
## allocation back to nearly equal sizes, and the cell would be a duplicate of
## the balanced one under a different name.
IMB_MIN <- 1.5

## ---------------------------------------------------------------------------
## Geometry
## ---------------------------------------------------------------------------

## centres of the K clusters at overlap omega, nearest-neighbour distance 2*omega
dgp_centres <- function(K, omega) {
  if (K == 1L) return(cbind(0, 0))
  if (K == 2L) return(cbind(c(-omega, omega), c(0, 0)))
  if (K == 3L) {
    r <- 2 * omega / sqrt(3); a <- c(90, 210, 330) * pi / 180
    return(cbind(r * cos(a), r * sin(a)))
  }
  if (K == 4L) return(cbind(c(omega, -omega, omega, -omega), c(omega, omega, -omega, -omega)))
  stop("K must be 1, 2, 3 or 4")
}

## Mean per-coordinate variance of the K centres. With K = 1 the design has no
## between-cluster spread, so the WITHIN-cluster dispersion is inflated by this
## amount: the pooled configuration then covers the same area as the clustered
## one at the same overlap, and the selection rule is not handed a free
## geometric cue for telling k = 1 from k > 1.
dgp_centre_var <- function(K, omega) {
  mu <- dgp_centres(K, omega)
  if (K == 1L) return(0)
  mean(apply(mu, 2, function(v) mean((v - mean(v))^2)))
}

## Regime sizes.
##
##   "balanced"    exactly equal, the remainder spread over the first regimes
##   "unbalanced"  sizes proportional to 1 : 2 : ... : K, so the largest regime
##                 is K times the smallest, then corrected so that no regime
##                 falls below n_min, the excess being taken from the largest
##
## The correction is what makes the unbalanced case usable: SC-STEM fits a full
## STEM model inside every regime, so a regime with fewer units than the r + 6
## parameters it has to estimate is not a hard case, it is an infeasible one,
## and a design that produced those would confound the recovery of the partition
## with the feasibility of the fit. `dgp_feasible()` says whether a cell admits
## the requested imbalance at all.
dgp_sizes <- function(n, K, balance = c("balanced", "unbalanced"),
                      n_min = N_MIN) {
  balance <- match.arg(balance)
  if (K == 1L) return(n)
  if (balance == "balanced") {
    s <- rep(n %/% K, K)
    if (n %% K) s[seq_len(n %% K)] <- s[seq_len(n %% K)] + 1L
    return(as.integer(s))
  }
  w <- seq_len(K) / sum(seq_len(K))
  s <- pmax(1L, as.integer(round(n * w)))
  s[K] <- n - sum(s[-K])
  short <- pmax(0L, n_min - s)
  if (any(short)) {
    s <- pmax(s, n_min)
    s[K] <- n - sum(s[-K])
  }
  as.integer(s)
}

## TRUE when the cell can carry K regimes of at least n_min units each AND, for
## an unbalanced design, when the imbalance survives that constraint
dgp_feasible <- function(n, K, balance = "balanced", n_min = N_MIN,
                         imb_min = IMB_MIN) {
  s <- dgp_sizes(n, K, balance, n_min)
  ok <- all(s >= n_min) && sum(s) == n
  if (ok && K > 1L && balance == "unbalanced") ok <- max(s) / min(s) >= imb_min
  ok
}

dgp_labels <- function(n, K, balance = "balanced", n_min = N_MIN) {
  rep(seq_len(K), times = dgp_sizes(n, K, balance, n_min))
}

## Locations: the abstract cloud of the overlap design, mapped to longitude and
## latitude so that distances are kilometres and the covariance parameters keep
## their meaning.
##
##   s_i | g_i = g  ~  N_2( mu_g , nu_sp I_2 )
##
## with mu_g = dgp_centres(K, omega). For K = 1 the dispersion is nu_sp + the
## between-centre variance the design would have had at the same overlap, so
## that the pooled case is not simply a smaller map.
dgp_locations <- function(n, K, omega, nu_sp = NU_SP, balance = "balanced",
                          seed = 1, k_ref = 2L) {
  set.seed(seed)
  mu <- dgp_centres(K, omega)
  g  <- dgp_labels(n, K, balance)
  sd_i <- sqrt(if (K == 1L) nu_sp + dgp_centre_var(k_ref, omega) else nu_sp)
  xy <- cbind(mu[g, 1] + stats::rnorm(n, sd = sd_i),
              mu[g, 2] + stats::rnorm(n, sd = sd_i))
  coords <- cbind(
    lon = LON0 + xy[, 1] * UNIT_KM / (111.320 * cos(LAT0 * pi / 180)),
    lat = LAT0 + xy[, 2] * UNIT_KM / 110.574)
  list(coords = coords, labels = g, xy = xy, mu = mu, sizes = tabulate(g, K))
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
  list(TN      = c(60L, 120L, 365L),
       n       = c(20L, 40L, 60L, 100L, 200L, 400L),
       K       = c(1L, 2L, 3L),
       omega   = c(0, 1/3, 2/3, 1),
       balance = c("balanced", "unbalanced"))
}

## ---------------------------------------------------------------------------
## One complete data set of the design, ready for STEM_Model()
## ---------------------------------------------------------------------------
dgp_draw <- function(n, TN, K, omega, scenario_row, rep = 1L,
                     balance = "balanced") {
  seed <- 1000L * rep + 1L
  loc  <- dgp_locations(n, K, omega, balance = balance, seed = seed)
  x    <- dgp_covariate(loc$coords, TN, seed = seed + 1L)
  psi  <- dgp_psi(scenario_row$scenario, scenario_row$level, K = K)
  z    <- dgp_simulate(loc$labels, psi, x, loc$coords,
                       rho = scenario_row$rho, seed = seed + 2L,
                       force_block = scenario_row$force_block)
  list(z = z, covariates = dgp_design(x), coordinates = loc$coords,
       labels = loc$labels, psi = psi, sizes = loc$sizes,
       latent = attr(z, "latent"), x = x)
}
