### ---------------------------------------------------------------------------
### Internal helpers shared by the spatially-clustered STEM (SC-STEM) routines.
###
### The algorithmic devices implemented here -- a Potts-penalized objective
### evaluated on undirected neighbor pairs, an ICM (Iterated Conditional Modes)
### label update, information criteria computed on the final refit, and a
### refit-with-clustering parametric bootstrap -- are adapted from the
### spatially-clustered small area estimation work of Maranzano, Mattera and
### Sugasawa (2026+). The statistical model they serve here is different: STEM
### is a hierarchical dynamic model for point-referenced spatio-temporal data
### with a latent state and a spatially correlated error, so each device had to
### be reworked for that setting rather than transferred.
###
### None of these functions is exported.
### ---------------------------------------------------------------------------


### ---------------------------------------------------------------------------
### Neighborhood structure
### ---------------------------------------------------------------------------
### The Potts penalty is defined on an UNDIRECTED graph, so the
### k-nearest-neighbor graph produced by spdep::knearneigh() -- which is
### asymmetric by construction -- is symmetrized before use: j is a neighbor of
### i whenever i is among the k nearest of j or vice versa. Without this step
### the same pair (i,j) would contribute to the assignment score of one unit but
### not of the other, and the sequential ICM sweep would not be maximizing a
### well-defined objective.
###
### The graph has to be built with the SAME metric as the covariance, otherwise
### the two spatial ingredients of the model disagree with one another. On
### geographic coordinates, neighbors computed on raw degrees are the neighbors
### of a planar metric in which a degree of longitude and a degree of latitude
### count the same: at 45 degrees of latitude a degree of longitude is about
### 78 km against 111 km, so such a graph systematically prefers north-south
### neighbors to east-west ones. With distance = "geo" the neighbors are
### measured on the sphere, as Sigma_e is.
`scstem_neighbors` <- function(coordinates, knn = 5,
                               distance = c("geo", "euclidean")) {

  distance <- match.arg(distance)
  d <- nrow(coordinates)
  if (knn < 1 || knn >= d) {
    stop("'knn' must be a positive integer strictly smaller than the number of locations.",
         call. = FALSE)
  }

  ### the package takes the first column for the longitude and the second for
  ### the latitude throughout, as geodist::geodist() does; naming them says so
  ### to spdep as well, which would otherwise assume it out loud
  cc <- as.matrix(coordinates)[, 1:2, drop = FALSE]
  colnames(cc) <- c("lon", "lat")

  nb <- spdep::knn2nb(spdep::knearneigh(cc, k = knn,
                                        longlat = identical(distance, "geo")))
  W <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)
  W <- pmax(W, t(W))
  diag(W) <- 0

  nb_list <- lapply(seq_len(d), function(i) which(W[i, ] > 0))

  list(W = W, nb = nb_list, knn = knn, distance = distance)
}


### ---------------------------------------------------------------------------
### Potts term: number of concordant neighbor pairs, each counted once
### ---------------------------------------------------------------------------
`scstem_potts_pairs` <- function(labels, nb) {
  s <- 0
  for (i in seq_along(labels)) {
    nbi <- nb[[i]]
    nbi <- nbi[nbi > i]
    if (length(nbi)) s <- s + sum(labels[nbi] == labels[i])
  }
  s
}


### ---------------------------------------------------------------------------
### Number of free parameters of a single cluster-wise STEM fit
### ---------------------------------------------------------------------------
### Per cluster: beta (ncov), sigma2eps, sigma2omega, theta (3), plus the state
### equation blocks G, Sigmaeta and m0. C0 is held fixed and is not counted.
`scstem_npar` <- function(ncov, pdim = 1, Gdiag = TRUE, Sigmaetadiag = TRUE) {
  nG <- if (isTRUE(Gdiag)) pdim else pdim^2
  nS <- if (isTRUE(Sigmaetadiag)) pdim else pdim * (pdim + 1) / 2
  ncov + 3 + nG + nS + pdim
}


### ---------------------------------------------------------------------------
### Per-location log-likelihood contribution used by the assignment step
### ---------------------------------------------------------------------------
### ALERT (methodological). Under the STEM measurement equation the observations
### of different locations at the same time point are spatially correlated
### through Sigma_e = sigma2eps * I + sigma2omega * C(h; theta), so the exact
### marginal log-likelihood does NOT factorize across locations and no exact
### per-location contribution exists. The assignment step therefore uses a
### PSEUDO-LIKELIHOOD: conditionally on the smoothed latent state path of
### cluster k, location i contributes
###
###   l_ik = sum_t log N( z_it ; x_it beta_k + K_i yhat_t(k) ,
###                       sigma2eps_k + sigma2omega_k )
###
### which is the exact conditional density of the series of location i given the
### latent state, with the marginal error variance
### sigma2eps_k + sigma2omega_k = diag(Sigma_e_k) (the exponential correlation
### function equals 1 at distance zero). It plays the role that the per-unit
### density plays in clusterwise regression, and it is used ONLY to rank
### clusters in the label update: all reported quantities -- coefficients, variance
### components, information criteria -- come from the exact cluster-wise
### likelihoods returned by STEM_Estimation() on the final partition.
###
### Arguments
###   z_i          T x 1 numeric, observations of location i
###   X_i          T x ncov numeric, covariates of location i (time-ordered)
###   beta         ncov x 1 numeric, cluster-wise regression coefficients
###   ysm          T x pdim numeric, smoothed latent state of the cluster
###   K_i          1 x pdim numeric, loading row of location i
###   sigma2eps    scalar, cluster-wise nugget variance
###   sigma2omega  scalar, cluster-wise spatial variance
`scstem_loglike_i` <- function(z_i, X_i, beta, ysm, K_i, sigma2eps, sigma2omega) {

  v <- as.numeric(sigma2eps) + as.numeric(sigma2omega)
  if (!is.finite(v) || v <= 0) return(-Inf)

  ysm <- as.matrix(ysm)
  K_i <- matrix(as.numeric(K_i), nrow = 1)

  fit_t <- as.numeric(as.matrix(X_i) %*% matrix(as.numeric(beta), ncol = 1)) +
    as.numeric(ysm %*% t(K_i))

  ### Only the time points at which this location was observed contribute, and
  ### the count of terms is its own T_i. The missingness pattern of a location
  ### does not depend on the cluster it is being scored against, so the scores
  ### stay comparable across clusters, which is all the assignment step needs.
  keep <- !is.na(z_i)
  Tobs <- sum(keep)
  if (Tobs == 0L) return(-Inf)

  res <- z_i[keep] - fit_t[keep]
  out <- -0.5 * Tobs * log(2 * pi * v) - 0.5 * sum(res^2) / v
  if (!is.finite(out)) return(-Inf)
  out
}


### ---------------------------------------------------------------------------
### Initial partition
### ---------------------------------------------------------------------------
### The starting partition comes from k-means on
### the location-wise summaries of the COVARIATES only (compressed by PCA at 90%
### of cumulative variance), with multiple external restarts and a
### minimum-cluster-size admissibility filter. Initializing on the covariates
### leaves spatial contiguity entirely to the Potts penalty, so that phi can be
### read as the price of spatial coherence rather than as a constraint built
### into the starting point. Intercept-only designs fall back to the
### coordinates.
###
### Versions of the package before 2.0.0 initialized the partition with
### SCDA::SC_AMKM(). That dependency has been dropped, because SCDA is not
### distributed on CRAN and a hard dependency on it would make this package
### unpublishable. Any external initialization -- AMKM included -- can still be
### used by passing it to SCSTEM_Estimation() through the `init_partition` argument.
###
### Arguments
###   Xmeans     d x ncov numeric, per-location averages of the covariates
###   coords     d x 2 numeric, spatial coordinates
###   k          number of clusters
###   method     "kmeans" (default) or "coordinates"
###   min_size   minimum admissible cluster size
`scstem_init` <- function(Xmeans, coords, k, method = c("kmeans", "coordinates"),
                          min_size = 2L, nstart_ext = 50L, nstart_int = 25L) {

  method <- match.arg(method)
  d <- nrow(coords)
  if (k == 1) return(rep(1L, d))

  ### Feature space for the k-means starts
  feat <- try({
    if (method == "coordinates" || is.null(Xmeans) || ncol(as.matrix(Xmeans)) == 0) {
      scale(coords)
    } else {
      Xm <- as.matrix(Xmeans)
      ### drop constant columns (typically the intercept), which carry no
      ### information and make prcomp(scale. = TRUE) fail
      keep <- apply(Xm, 2, function(x) stats::sd(x, na.rm = TRUE) > 0)
      if (!any(keep)) {
        scale(coords)
      } else {
        pca <- stats::prcomp(Xm[, keep, drop = FALSE], center = TRUE, scale. = TRUE)
        cum <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
        ncomp <- which(cum >= 0.90)[1]
        as.matrix(pca$x[, seq_len(ncomp), drop = FALSE])
      }
    }
  }, silent = TRUE)
  if (inherits(feat, "try-error")) feat <- scale(coords)

  WSS <- rep(Inf, nstart_ext)
  CL <- vector("list", nstart_ext)
  valid <- rep(FALSE, nstart_ext)
  for (m in seq_len(nstart_ext)) {
    CL[[m]] <- tryCatch(
      stats::kmeans(x = feat, centers = k, nstart = nstart_int, iter.max = 100),
      error = function(e) NULL
    )
    if (!is.null(CL[[m]])) {
      tab <- table(CL[[m]]$cluster)
      valid[m] <- length(tab) == k && all(tab >= min_size)
      WSS[m] <- CL[[m]]$tot.withinss
    }
  }

  if (any(valid)) {
    as.integer(CL[[which.min(ifelse(valid, WSS, Inf))]]$cluster)
  } else if (any(is.finite(WSS))) {
    ### No restart met the admissibility filter (typical when a couple of
    ### outlying locations dominate a skewed covariate): repair the best
    ### solution instead of returning a partition that would collapse at the
    ### first cluster-wise fit.
    lab <- as.integer(CL[[which.min(WSS)]]$cluster)
    lab <- scstem_repair_partition(lab, feat = feat, k = k, min_size = min_size)
    lab
  } else {
    stop("The initialization step failed: k-means could not produce any partition.",
         call. = FALSE)
  }
}


### ---------------------------------------------------------------------------
### Adjusted Rand Index (used by the stability step of the tuning rule)
### ---------------------------------------------------------------------------
`scstem_ari` <- function(x, y) {
  x <- as.integer(as.factor(x))
  y <- as.integer(as.factor(y))
  if (length(x) != length(y)) stop("Partitions of different length.", call. = FALSE)
  tab <- table(x, y)
  n <- sum(tab)
  if (n < 2) return(NA_real_)
  choose2 <- function(v) sum(v * (v - 1) / 2)
  sum_ij <- choose2(as.vector(tab))
  sum_i <- choose2(rowSums(tab))
  sum_j <- choose2(colSums(tab))
  tot <- n * (n - 1) / 2
  expected <- sum_i * sum_j / tot
  maxi <- (sum_i + sum_j) / 2
  if (isTRUE(all.equal(maxi, expected))) return(1)
  (sum_ij - expected) / (maxi - expected)
}


### ---------------------------------------------------------------------------
### Per-location averages of the covariates
### ---------------------------------------------------------------------------
### The covariate matrix of a STEM_Model stacks the (T x ncov) blocks of the d
### locations by row. This helper returns the d x ncov matrix of location-wise
### time averages used by the initialization step.
`scstem_covariate_means` <- function(covariates, d, Tobs) {
  X <- as.matrix(covariates)
  out <- matrix(NA_real_, nrow = d, ncol = ncol(X))
  for (i in seq_len(d)) {
    rows <- ((i - 1) * Tobs + 1):(i * Tobs)
    out[i, ] <- colMeans(X[rows, , drop = FALSE], na.rm = TRUE)
  }
  colnames(out) <- colnames(X)
  out
}


### ---------------------------------------------------------------------------
### Per-location departures from the pooled fit (init_method = "departures")
### ---------------------------------------------------------------------------
### The residual of each location from the signal of the pooled fit,
### r_i = z_i - muhat_i, summarized by four kinds of feature: its mean (how the
### level of the location departs from the pooled one), the slopes of r_i on the
### location's covariates (how its response to them departs), and the lag-one
### autocorrelation and the log-variance of what that regression leaves (how its
### dynamics and its noise depart). Locations of one regime share their
### departures and locations of different regimes do not, which the covariate
### means used by init_method = "kmeans" cannot see when the covariates are
### exogenous to the regimes.
###
### Arguments
###   z           T x d response, NA where missing
###   covariates  (d T) x ncov covariates, stacked by location
###   signal      T x d signal of the pooled fit (STEM_Signal())
###   Tobs        number of time points
### Returns a d x m matrix; a feature that cannot be computed at a location
### (too few observations) takes the median across locations.
`scstem_departures` <- function(z, covariates, signal, Tobs) {
  d <- ncol(z)
  X <- as.matrix(covariates)
  ### the slopes are taken on the covariates that vary; the intercept is the
  ### column of the regression below
  vary <- apply(X, 2, function(v) isTRUE(stats::sd(v, na.rm = TRUE) > 0))
  Xs <- X[, vary, drop = FALSE]
  m <- 1L + ncol(Xs) + 2L
  out <- matrix(NA_real_, nrow = d, ncol = m)
  for (i in seq_len(d)) {
    r <- z[, i] - signal[, i]
    Xi <- Xs[scstem_rows(i, Tobs), , drop = FALSE]
    ok <- !is.na(r) & stats::complete.cases(Xi)
    if (sum(ok) < ncol(Xs) + 4L) next
    fit <- stats::lm.fit(cbind(1, Xi[ok, , drop = FALSE]), r[ok])
    b <- fit$coefficients[-1]
    b[is.na(b)] <- 0
    u <- rep(NA_real_, length(r))
    u[ok] <- fit$residuals
    lag_ok <- !is.na(u[-1]) & !is.na(u[-length(u)])
    ar1 <- if (sum(lag_ok) >= 3L)
      suppressWarnings(stats::cor(u[-1][lag_ok], u[-length(u)][lag_ok])) else NA_real_
    out[i, ] <- c(mean(r[ok]), b, ar1,
                  log(max(stats::var(u, na.rm = TRUE), .Machine$double.eps)))
  }
  for (j in seq_len(m)) {
    miss <- !is.finite(out[, j])
    if (any(miss)) out[miss, j] <- stats::median(out[!miss, j])
  }
  colnames(out) <- c("level", paste0("slope_", seq_len(ncol(Xs))), "ar1", "logvar")
  out
}


### ---------------------------------------------------------------------------
### Row indices of a set of locations in the stacked covariate matrix
### ---------------------------------------------------------------------------
`scstem_rows` <- function(idx, Tobs) {
  unlist(lapply(idx, function(i) ((i - 1) * Tobs + 1):(i * Tobs)))
}


### ---------------------------------------------------------------------------
### Majority-rule alignment of a refit partition onto a reference partition
### ---------------------------------------------------------------------------
### Used by the bootstrap post-processing: every cluster of the refit is mapped
### to the reference cluster with which it shares the largest number of units,
### resolving conflicts greedily by decreasing overlap. Returns an integer
### vector of length K_refit giving the reference label of each refit cluster
### (NA when no reference cluster is left to claim).
`scstem_align_labels` <- function(reference, refit, K) {
  tab <- matrix(0L, nrow = K, ncol = K)
  for (g in seq_len(K)) {
    for (h in seq_len(K)) {
      tab[g, h] <- sum(refit == g & reference == h)
    }
  }
  map <- rep(NA_integer_, K)
  taken <- rep(FALSE, K)
  ord <- order(tab, decreasing = TRUE)
  for (pos in ord) {
    g <- ((pos - 1) %% K) + 1
    h <- ((pos - 1) %/% K) + 1
    if (is.na(map[g]) && !taken[h] && tab[g, h] > 0) {
      map[g] <- h
      taken[h] <- TRUE
    }
  }
  ### clusters with no overlap at all receive whatever reference label is left
  left <- which(is.na(map))
  free <- which(!taken)
  if (length(left) && length(free)) {
    map[left[seq_along(free)]] <- free
  }
  map
}


### ---------------------------------------------------------------------------
### Reproducible evaluation without leaving the RNG state modified
### ---------------------------------------------------------------------------
### CRAN policies forbid packages from modifying the user's global environment.
### This helper therefore evaluates expr under the requested seed and then puts
### the RNG stream back exactly as it was found -- including the case in which
### no .Random.seed existed before the call, in which case the object created by
### seeding is removed again. When seed is NULL nothing is touched at all and
### the user is expected to call set.seed() beforehand.
`scstem_with_seed` <- function(seed, expr) {
  if (is.null(seed)) return(force(expr))
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  set.seed(seed)
  force(expr)
}


### ---------------------------------------------------------------------------
### Enforce a minimum cluster size on a candidate partition
### ---------------------------------------------------------------------------
### k-means on strongly skewed covariates routinely isolates a couple of
### outlying locations, so that no restart satisfies the minimum-size
### admissibility filter. Rather than falling back on an inadmissible partition
### -- which collapses at the very first cluster-wise fit -- the candidate is
### repaired: while some cluster is short of min_size, the units closest (in the
### feature space used for the initialization) to the centroid of the most
### deficient cluster are moved into it, taken from the clusters that can
### afford to lose them.
`scstem_repair_partition` <- function(labels, feat, k, min_size) {

  labels <- as.integer(labels)
  feat <- as.matrix(feat)
  n <- length(labels)
  if (k * min_size > n) {
    stop("A partition into ", k, " clusters of at least ", min_size,
         " units each is impossible with ", n, " locations.", call. = FALSE)
  }

  centroid <- function(g) {
    idx <- which(labels == g)
    if (!length(idx)) return(rep(NA_real_, ncol(feat)))
    colMeans(feat[idx, , drop = FALSE])
  }

  guard <- 0L
  repeat {
    sizes <- tabulate(labels, nbins = k)
    short <- which(sizes < min_size)
    if (!length(short)) break
    guard <- guard + 1L
    if (guard > n * k) break

    g <- short[which.min(sizes[short])]
    cen <- centroid(g)
    ### an empty cluster has no centroid: seed it with the unit that is
    ### farthest from its own cluster centroid
    if (anyNA(cen)) {
      donors <- which(sizes[labels] > min_size)
      if (!length(donors)) donors <- which(sizes[labels] > 1)
      if (!length(donors)) break
      dd <- vapply(donors, function(i) {
        sum((feat[i, ] - centroid(labels[i]))^2)
      }, numeric(1))
      labels[donors[which.max(dd)]] <- g
      next
    }

    donors <- which(labels != g & sizes[labels] > min_size)
    if (!length(donors)) donors <- which(labels != g & sizes[labels] > 1)
    if (!length(donors)) break

    dd <- vapply(donors, function(i) sum((feat[i, ] - cen)^2), numeric(1))
    labels[donors[which.min(dd)]] <- g
  }

  labels
}


### ---------------------------------------------------------------------------
### Penalty contribution of the pairs incident to two locations
### ---------------------------------------------------------------------------
### Counts, once each, the concordant neighbor pairs that involve i or j. Used
### to evaluate the exact change of the Potts term produced by swapping the
### labels of i and j, without recomputing the whole quadratic form.
`scstem_pen_local` <- function(labels, i, j, nb) {
  ni <- nb[[i]]; ni <- ni[ni != j]
  nj <- nb[[j]]; nj <- nj[nj != i]
  s <- 0
  if (length(ni)) s <- s + sum(labels[ni] == labels[i])
  if (length(nj)) s <- s + sum(labels[nj] == labels[j])
  if (j %in% nb[[i]]) s <- s + as.integer(labels[i] == labels[j])
  s
}


### ---------------------------------------------------------------------------
### Size-preserving swap pass
### ---------------------------------------------------------------------------
### The constrained ICM sweep cannot move a location out of a cluster that sits
### at the minimum admissible size, so with few locations a large share of the
### network can stay frozen at the initial partition. A swap exchanges the
### labels of two locations in different clusters: it leaves every cluster size
### unchanged -- hence feasibility is preserved by construction -- and it is
### accepted only when it strictly increases the penalized objective, so the
### monotonicity of the alternating algorithm is preserved as well.
###
### The pass is greedy: candidate pairs are scanned and every improving swap is
### applied immediately, repeating until no improving swap is left or max_pass
### sweeps have been performed.
###
### Arguments
###   labels   current partition
###   LL       d x k matrix of log-likelihood contributions
###   phi_eff  effective penalty
###   nb       neighbor list
###   tol      minimum improvement required to accept a swap
###   max_pass maximum number of full scans
###
### Value: a list with the updated labels and the number of accepted swaps.
### The scan is made affordable by an EXACT pruning bound. The penalty term of
### a swap can gain at most phi_eff * (|nb_i| + |nb_j| + 1), because that is the
### largest number of concordant pairs the two locations can be involved in.
### A pair whose likelihood delta already falls below minus that bound can never
### improve the objective and is discarded without evaluating the Potts term.
### The likelihood deltas are computed cluster pair by cluster pair with an
### outer sum, so the only loop left runs over the surviving candidates, which
### are examined in decreasing order of likelihood gain.
`scstem_swap_pass` <- function(labels, LL, phi_eff, nb, tol = 1e-8, max_pass = 5L) {

  d <- length(labels)
  k <- ncol(LL)
  nswap <- 0L
  nbsize <- vapply(nb, length, integer(1))

  for (pass in seq_len(max_pass)) {

    improved <- FALSE

    for (a in seq_len(k - 1L)) {
      for (b in (a + 1L):k) {

        Ia <- which(labels == a)
        Ib <- which(labels == b)
        if (!length(Ia) || !length(Ib)) next

        ### likelihood gain of moving i from a to b, and j from b to a
        gi <- LL[Ia, b] - LL[Ia, a]
        gj <- LL[Ib, a] - LL[Ib, b]
        dll <- outer(gi, gj, "+")

        ### exact upper bound on the penalty gain of each candidate pair
        bound <- phi_eff * outer(nbsize[Ia], nbsize[Ib], "+") + phi_eff
        cand <- which(is.finite(dll) & (dll + bound > tol), arr.ind = TRUE)
        if (!nrow(cand)) next

        ### examine the most promising candidates first
        cand <- cand[order(dll[cand], decreasing = TRUE), , drop = FALSE]

        for (r in seq_len(nrow(cand))) {
          i <- Ia[cand[r, 1L]]
          j <- Ib[cand[r, 2L]]
          ### a previous accepted swap in this scan may have moved i or j
          if (labels[i] != a || labels[j] != b) next
          dll_ij <- LL[i, b] + LL[j, a] - LL[i, a] - LL[j, b]
          if (!is.finite(dll_ij)) next
          pen_before <- scstem_pen_local(labels, i, j, nb)
          labels[i] <- b; labels[j] <- a
          pen_after <- scstem_pen_local(labels, i, j, nb)
          if (dll_ij + phi_eff * (pen_after - pen_before) > tol) {
            nswap <- nswap + 1L
            improved <- TRUE
          } else {
            labels[i] <- a; labels[j] <- b
          }
        }
      }
    }

    if (!improved) break
  }

  list(labels = labels, nswap = nswap)
}
