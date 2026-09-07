## ---------------------------------------------------------------------------
## Simulation study of the SC-STEM paper (Section 4 of the manuscript).
##
## The designs are run on the real geometry of the povalley network, so that the
## irregular spacing, the varying local density and the elongated shape of the
## basin -- all of which affect both the nearest-neighbor graph that carries the
## penalty and the distance matrix that enters the spatial covariance -- are the
## ones of the application.
##
## Four designs, organized by WHICH component of the parameter set separates the
## regimes. D1 is the direct spatio-temporal counterpart of the cross-sectional
## literature; D2, D3 and D4 have no such counterpart.
##
##   D1  separation in the regression coefficients
##   D2  separation in the latent dynamics (G, Sigmaeta)
##   D3  separation in the spatial covariance (nugget-to-sill ratio, range)
##   D4  joint separation
##
## Each in a clear- and a poor-separation version. The true partition splits the
## basin either by longitude or by altitude.
##
##   Rscript dev/paper/02-simulation.R [--reps N] [--cores N] [--designs D1,D2]
## ---------------------------------------------------------------------------

PKG <- "C:/Users/paulm/OneDrive/Documenti/GitHub/Stem"
suppressMessages(pkgload::load_all(PKG, quiet = TRUE))

CACHE <- file.path("dev", "paper", "cache")
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)

## ---- command line ----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
getopt <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
REPS    <- as.integer(getopt("--reps", "200"))
CORES   <- as.integer(getopt("--cores", "16"))
WHICH   <- strsplit(getopt("--designs", "D1,D2,D3,D4"), ",")[[1]]
TN      <- as.integer(getopt("--T", "365"))
PHI_GRID <- c(0, 0.25, 0.5, 1)
KNN     <- 5

## ---- the geometry and the covariates ---------------------------------------
data(povalley)
Tfull  <- nrow(povalley$z)
d      <- ncol(povalley$z)
coords <- povalley$coords
alt    <- povalley$altitude

## the first TN days of the real covariates, stacked by station as the package
## expects
keep  <- as.vector(outer(seq_len(TN), (seq_len(d) - 1L) * Tfull, "+"))
Xsim  <- povalley$covariates[keep, , drop = FALSE]
Ksim  <- matrix(1, d, 1)

## ---- true partitions --------------------------------------------------------
part_longitude <- as.integer(coords[, 1] > stats::median(coords[, 1])) + 1L
part_altitude  <- as.integer(alt > stats::median(alt)) + 1L

## ---- true parameters --------------------------------------------------------
## The baseline is the pooled fit on the real data, so that the simulated series
## look like the ones the model is meant for.
base <- list(beta = matrix(c(2.34, 0.004, 0.58), 3, 1),
             sigma2eps = 21.4, sigma2omega = 4.19,
             theta = 8.1e-06,
             G = matrix(0.97, 1, 1), Sigmaeta = matrix(0.63, 1, 1),
             m0 = as.matrix(0), C0 = as.matrix(1))

## `sep` is 1 for the clear-separation version and shrinks the contrast towards
## zero for the poor-separation one, so that the two versions of a design differ
## only in how far apart the regimes are.
make_truth <- function(design, sep) {
  A <- base; B <- base
  mix <- function(x, y) x + sep * (y - x)          # from A towards the contrast

  if (design %in% c("D1", "D4")) {
    ## the response to PM10 halves, and the intercept moves with it
    B$beta <- matrix(c(mix(2.34, 5.00), 0.004, mix(0.58, 0.30)), 3, 1)
  }
  if (design %in% c("D2", "D4")) {
    ## a much less persistent latent process, with innovations rescaled so that
    ## the stationary variance of the state stays comparable
    gB <- mix(0.97, 0.60)
    B$G <- matrix(gB, 1, 1)
    B$Sigmaeta <- matrix(base$Sigmaeta[1, 1] * (1 - gB^2) / (1 - 0.97^2), 1, 1)
  }
  if (design %in% c("D3", "D4")) {
    ## the same total residual variance, split differently between nugget and
    ## partial sill, and a much shorter range
    tot <- base$sigma2eps + base$sigma2omega
    shareB <- mix(base$sigma2eps / tot, 0.40)
    B$sigma2eps   <- tot * shareB
    B$sigma2omega <- tot * (1 - shareB)
    B$theta       <- mix(base$theta, 3e-05)
  }
  list(A, B)
}

## ---- one replication --------------------------------------------------------
## Data are generated regime by regime, as SCSTEM_Bootstrap() does: a STEM model
## is built on the stations of each regime with that regime's parameters, and
## the simulated columns are put back in the original order so that every
## station keeps its place in the neighborhood graph.
simulate_clustered <- function(truth, labels, seed) {
  set.seed(seed)
  z <- matrix(NA_real_, TN, d)
  for (g in sort(unique(labels))) {
    idx  <- which(labels == g)
    kp   <- as.vector(outer(seq_len(TN), (idx - 1L) * TN, "+"))
    modg <- STEM_Model(z = matrix(0, TN, length(idx)),
                       covariates = Xsim[kp, , drop = FALSE],
                       coordinates = coords[idx, , drop = FALSE],
                       phi = truth[[g]], K = matrix(1, length(idx), 1))
    z[, idx] <- as.matrix(STEM_Simulation(modg, distance = "geo"))
  }
  z
}

one_rep <- function(r, design, sep, labels, do_select) {
  truth <- make_truth(design, sep)
  z <- simulate_clustered(truth, labels, seed = 10000L * match(design, c("D1","D2","D3","D4")) +
                                              1000L * round(10 * sep) + r)
  mod <- STEM_Model(z = z, covariates = Xsim, coordinates = coords,
                    phi = base, K = Ksim)

  out <- list(rep = r, design = design, sep = sep)

  ## pooled benchmark
  pooled <- try(STEM_Estimation(mod, precision = 0.05, max.iter = 25,
                                distance = "geo"), silent = TRUE)
  out$pooled_ok <- !inherits(pooled, "try-error")
  if (out$pooled_ok) out$pooled_loglik <- pooled$estimates$loglik

  ## SC-STEM at the true k, over the penalty grid
  rows <- list()
  for (ph in PHI_GRID) {
    fit <- try(SCSTEM_Estimation(mod, k = 2L, phi_penalty = ph, knn = KNN,
                            distance = "geo", precision = 0.1,
                            precision_full_dataset = 0.05, max_iter = 8,
                            seed = 1000L + r, verbose = FALSE), silent = TRUE)
    if (inherits(fit, "try-error")) {
      rows[[length(rows) + 1L]] <- data.frame(phi = ph, ari = NA_real_,
                                              loglik = NA_real_, ok = FALSE)
      next
    }
    rows[[length(rows) + 1L]] <- data.frame(
      phi = ph,
      ari = Stem:::scstem_ari(labels, fit$group),
      loglik = unname(fit$info_crit[["loglik"]]),
      ok = TRUE)
  }
  out$grid <- do.call(rbind, rows)

  ## model selection, only on the sub-study that asks for it
  if (isTRUE(do_select)) {
    ic <- try(SCSTEM_Infocrit(mod, k_grid = 1:4, phi_grid = PHI_GRID,
                              knn = KNN, distance = "geo",
                              precision = 0.1, precision_full_dataset = 0.05,
                              max_iter = 6, seed = 1000L + r), silent = TRUE)
    if (!inherits(ic, "try-error")) {
      s <- try(SCSTEM_Select(ic, band = c(0.25, 1), criterion = "BIC"),
               silent = TRUE)
      if (!inherits(s, "try-error")) {
        out$k_selected   <- s$k_selected
        out$phi_selected <- s$phi_selected
        out$ari_selected <- Stem:::scstem_ari(labels, s$fit$group)
      }
    }
  }
  out
}

## ---- the scenarios ----------------------------------------------------------
scen <- expand.grid(design = c("D1", "D2", "D3", "D4"),
                    sep    = c(1, 0.4),
                    stringsAsFactors = FALSE)
scen$level <- ifelse(scen$sep == 1, "clear", "poor")
## D1 and D4 separate on the covariate response, which the altitude split makes
## most interpretable; D2 and D3 separate on dynamics and covariance, which the
## longitudinal split of the basin matches (ventilation improves eastwards)
scen$partition <- ifelse(scen$design %in% c("D1", "D4"), "altitude", "longitude")
scen <- scen[scen$design %in% WHICH, , drop = FALSE]

message("Designs: ", paste(unique(scen$design), collapse = ", "),
        " | replications: ", REPS, " | T = ", TN, " | cores: ", CORES)

cl <- parallel::makeCluster(min(CORES, parallel::detectCores() - 1L))
on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterExport(cl, c("PKG", "TN", "d", "coords", "alt", "Xsim", "Ksim",
                              "base", "make_truth", "simulate_clustered",
                              "one_rep", "PHI_GRID", "KNN"),
                        envir = environment())
invisible(parallel::clusterEvalQ(cl, {
  suppressMessages(pkgload::load_all(PKG, quiet = TRUE))
}))

for (i in seq_len(nrow(scen))) {
  s <- scen[i, ]
  tag <- paste0("sim_", s$design, "_", s$level)
  f <- file.path(CACHE, paste0(tag, ".rds"))
  if (file.exists(f)) { message("[cache] ", tag); next }

  labels <- if (s$partition == "altitude") part_altitude else part_longitude
  ## the workers see only what is exported, so the scenario is unpacked into
  ## plain variables rather than passed through the row of the data frame
  des  <- s$design
  sp   <- s$sep
  half <- max(1L, REPS %/% 2L)
  parallel::clusterExport(cl, c("labels", "des", "sp", "half"),
                          envir = environment())

  message("[run  ] ", tag, " (", s$partition, " split, sizes ",
          paste(table(labels), collapse = "/"), ")")
  t0 <- Sys.time()
  ## the selection sub-study is run on half the replications, since it fits the
  ## whole grid at every one of them
  res <- parallel::parLapply(cl, seq_len(REPS), function(r) {
    one_rep(r, design = des, sep = sp, labels = labels, do_select = (r <= half))
  })
  message("[done ] ", tag, " in ",
          round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
  saveRDS(list(scenario = s, labels = labels, res = res), f)
}

message("\nAll requested designs complete.")
