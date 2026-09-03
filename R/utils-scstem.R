#' @keywords internal
#' @noRd
#'
#' Internal helpers shared by the spatially-clustered STEM (SC-STEM) routines.
#'
#' The design mirrors the spatially-clustered Fay-Herriot (SC-FH) toolbox of
#' Maranzano, Mattera and Sugasawa (2026+), so that the two model families share
#' the same algorithmic conventions: a Potts-penalised objective evaluated on
#' undirected neighbour pairs, an ICM (Iterated Conditional Modes) label update,
#' information criteria computed on the final refit, and a
#' refit-with-clustering parametric bootstrap.
NULL


### ---------------------------------------------------------------------------
### Neighbourhood structure
### ---------------------------------------------------------------------------
### The Potts penalty is defined on an UNDIRECTED graph, so the
### k-nearest-neighbour graph produced by spdep::knearneigh() -- which is
### asymmetric by construction -- is symmetrised before use: j is a neighbour of
### i whenever i is among the k nearest of j or vice versa. Without this step
### the same pair (i,j) would contribute to the assignment score of one unit but
### not of the other, and the sequential ICM sweep would not be maximising a
### well-defined objective.
`scstem_neighbours` <- function(coordinates, knn = 5) {

  d <- nrow(coordinates)
  if (knn < 1 || knn >= d) {
    stop("'knn' must be a positive integer strictly smaller than the number of locations.",
         call. = FALSE)
  }

  nb <- spdep::knn2nb(spdep::knearneigh(as.matrix(coordinates), k = knn))
  W <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)
  W <- pmax(W, t(W))
  diag(W) <- 0

  nb_list <- lapply(seq_len(d), function(i) which(W[i, ] > 0))

  list(W = W, nb = nb_list, knn = knn)
}


### ---------------------------------------------------------------------------
### Potts term: number of concordant neighbour pairs, each counted once
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
### marginal log-likelihood does NOT factorise across locations and no exact
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
### function equals 1 at distance zero). This is the direct STEM analogue of
### FHloglike_i() in the SC-FH toolbox, and it is used ONLY to rank clusters in
### the label update: all reported quantities -- coefficients, variance
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

  Tobs <- length(z_i)
  v <- as.numeric(sigma2eps) + as.numeric(sigma2omega)
  if (!is.finite(v) || v <= 0) return(-Inf)

  ysm <- as.matrix(ysm)
  K_i <- matrix(as.numeric(K_i), nrow = 1)

  fit_t <- as.numeric(as.matrix(X_i) %*% matrix(as.numeric(beta), ncol = 1)) +
    as.numeric(ysm %*% t(K_i))

  res <- z_i - fit_t
  out <- -0.5 * Tobs * log(2 * pi * v) - 0.5 * sum(res^2) / v
  if (!is.finite(out)) return(-Inf)
  out
}


### ---------------------------------------------------------------------------
### Initial partition
### ---------------------------------------------------------------------------
### Following the SC-FH strategy, the starting partition comes from k-means on
### the location-wise summaries of the COVARIATES only (compressed by PCA at 90%
### of cumulative variance), with multiple external restarts and a
### minimum-cluster-size admissibility filter. Initialising on the covariates
### leaves spatial contiguity entirely to the Potts penalty, so that phi can be
### read as the price of spatial coherence rather than as a constraint built
### into the starting point. Intercept-only designs fall back to the
### coordinates. The "AMKM" method reproduces the pre-2.0.0 behaviour and
### requires the (non-CRAN) SCDA package.
###
### Arguments
###   Xmeans     d x ncov numeric, per-location averages of the covariates
###   coords     d x 2 numeric, spatial coordinates
###   k          number of clusters
###   method     "kmeans" (default), "AMKM" or "coordinates"
###   min_size   minimum admissible cluster size
###   crs        CRS passed to SCDA::SC_AMKM when method = "AMKM"
`scstem_init` <- function(Xmeans, coords, k, method = c("kmeans", "AMKM", "coordinates"),
                          min_size = 2L, crs = 4326, nstart_ext = 50L, nstart_int = 25L) {

  method <- match.arg(method)
  d <- nrow(coords)
  if (k == 1) return(rep(1L, d))

  if (method == "AMKM") {
    if (!requireNamespace("SCDA", quietly = TRUE)) {
      stop("init_method = 'AMKM' requires the 'SCDA' package, which is not available on CRAN. ",
           "Install it from its development repository, or use init_method = 'kmeans'.",
           call. = FALSE)
    }
    dati <- cbind(as.data.frame(coords), as.data.frame(Xmeans))
    cl <- SCDA::SC_AMKM(Data_sf = sf::st_as_sf(dati, coords = c(1, 2)),
                        Method = "AMKM", MinNc = k, MaxNc = k, IndexCol = 0, CRS = crs)
    return(as.integer(cl$df$cluster))
  }

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
    stop("The initialisation step failed: k-means could not produce any partition.",
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
### time averages used by the initialisation step.
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
### feature space used for the initialisation) to the centroid of the most
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
