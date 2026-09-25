## ---------------------------------------------------------------------------
## Validation of the refit-with-clustering parametric bootstrap.
##
## The point estimator is one thing; the interval built around it is another,
## and the second does not follow from the first. This experiment asks whether
## the bootstrap of SCSTEM_Bootstrap() / SCSTEM_BootInference() is calibrated:
##
##   (a) COVERAGE. Does a nominal 95 per cent interval for a regime-specific
##       parameter contain the true value 95 per cent of the time? Reported for
##       all four interval families the package returns -- normal, basic,
##       percentile, bias-corrected -- at two nominal levels.
##   (b) STANDARD ERRORS. Is the average bootstrap standard error the same size
##       as the Monte Carlo standard deviation of the estimator it claims to
##       describe? Their ratio is the sharpest single diagnostic: it separates
##       an interval that is the wrong width from one that is the wrong shape.
##   (c) THE PARTITION. How often does a bootstrap refit recover the partition
##       of the fit it was generated from? The bootstrap re-runs the clustering
##       on every draw, so its draws carry label uncertainty; if the refits
##       scatter, the intervals widen for a reason that is real and should be
##       visible rather than hidden.
##
## Each Monte Carlo replication costs one fit plus B refits of the WHOLE
## procedure, so this is the expensive experiment of the study and is run on a
## deliberately small set of cells, at the TRUE number of regimes and a fixed
## penalty. That is the same choice the companion SC-SAE study makes: coverage
## is asked of the estimator, not of the selection rule, and mixing the two
## would leave a failure unattributable.
##
##   Rscript dev/paper/14-bootstrap-coverage.R --nrep=200 --B=200
##   Rscript dev/paper/14-bootstrap-coverage.R --scenario=S1b,S4 --omega=0,2/3
##   Rscript dev/paper/14-bootstrap-coverage.R --nrep=2 --B=20 --tag=smoke
##
## Results are appended replication by replication, so the run resumes.
## ---------------------------------------------------------------------------

local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  here <- if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
  source(file.path(here, "00-setup.R"), chdir = TRUE)
  source(file.path(here, "06-dgp.R"), chdir = TRUE)
}, envir = globalenv())

stem_load()
scen <- dgp_scenarios()

CFG <- stem_config(list(
  n        = 60L,
  TN       = 365L,
  K        = 2L,
  omega    = c(0, 2/3),
  balance  = "balanced",
  scenario = c("S1b", "S4"),
  phi      = 0.5,       # the penalty the fit and every refit are run at
  nrep     = 200L,      # Monte Carlo replications
  B        = 200L,      # bootstrap refits per replication
  levels   = c(0.90, 0.95),
  alpha    = 0,
  lambda   = 0,
  tag      = "coverage",
  seed0    = 5000L
))
stem_print_config(CFG)

CSV <- file.path(stem_cache_dir(), sprintf("14-coverage-%s.csv", CFG$tag))
CSA <- file.path(stem_cache_dir(), sprintf("14-coverage-%s-stability.csv", CFG$tag))

cells <- expand.grid(n = CFG$n, TN = CFG$TN, K = CFG$K, omega = CFG$omega,
                     id = CFG$scenario, balance = CFG$balance,
                     stringsAsFactors = FALSE)
cells <- cells[mapply(dgp_feasible, cells$n, cells$K, cells$balance), ]
rownames(cells) <- NULL
cat(sprintf("%d cells x %d replications x %d refits, writing to\n  %s\n\n",
            nrow(cells), CFG$nrep, CFG$B, CSV))

## ---------------------------------------------------------------------------
## One replication: fit, bootstrap, and record the interval against the truth
## ---------------------------------------------------------------------------
run_one <- function(cell, rep) {

  row <- scen[scen$id == cell$id, , drop = FALSE]
  dat <- dgp_draw(cell$n, cell$TN, cell$K, cell$omega, row,
                  rep = CFG$seed0[1] + rep, balance = cell$balance)

  ols <- stats::lm.fit(x = dat$covariates, y = as.vector(dat$z))
  s2  <- stats::var(ols$residuals)
  mod <- STEM_Model(
    z = dat$z, covariates = dat$covariates, coordinates = dat$coordinates,
    phi = list(beta = matrix(ols$coefficients, ncol = 1),
               sigma2eps = 0.7 * s2, sigma2omega = 0.3 * s2,
               theta = 1 / 100000, G = matrix(0.8, 1, 1),
               Sigmaeta = matrix(0.2 * s2, 1, 1),
               m0 = as.matrix(0), C0 = as.matrix(1)),
    K = matrix(1, cell$n, 1))

  t0 <- proc.time()[["elapsed"]]
  fit <- SCSTEM_Estimation(mod, k = cell$K, phi_penalty = CFG$phi[1],
                           distance = "geo", verbose = FALSE,
                           alpha = CFG$alpha[1], lambda = CFG$lambda[1])
  boot <- SCSTEM_Bootstrap(fit, B = CFG$B[1], seed = CFG$seed0[1] + rep,
                           verbose = FALSE)
  secs <- proc.time()[["elapsed"]] - t0

  ## The bootstrap aligns every refit onto the ORIGINAL fit. Aligning the
  ## original fit onto the TRUTH is this script's job, and without it the
  ## comparison would be between a regime and whichever true regime happens to
  ## carry the same arbitrary label.
  map <- scstem_align_labels(reference = dat$labels, refit = fit$group,
                             K = cell$K)
  tru <- dgp_truth(dat$psi)

  key <- data.frame(n = cell$n, TN = cell$TN, K = cell$K, omega = cell$omega,
                    id = cell$id, balance = cell$balance, rep = rep,
                    stringsAsFactors = FALSE)

  out <- list()
  for (lev in CFG$levels) {
    inf <- try(SCSTEM_BootInference(boot, level = lev, digits = 12),
               silent = TRUE)
    if (inherits(inf, "try-error")) next
    s <- inf$summary
    for (i in seq_len(nrow(tru))) {
      g_fit <- which(map == tru$regime[i])
      p     <- tru$parameter[i]
      j <- if (length(g_fit)) which(s$cluster == g_fit[1] & s$parameter == p) else integer(0)
      if (!length(j)) next
      r <- s[j[1], ]
      out[[length(out) + 1L]] <- cbind(key, data.frame(
        level = lev, regime = tru$regime[i], parameter = p,
        truth = tru$truth[i], estimate = r$estimate, se = r$se,
        normal_lo = r$normal_lo, normal_up = r$normal_up,
        basic_lo  = r$basic_lo,  basic_up  = r$basic_up,
        perc_lo   = r$perc_lo,   perc_up   = r$perc_up,
        bc_lo     = r$bc_lo,     bc_up     = r$bc_up,
        n_draws = r$n_draws, stringsAsFactors = FALSE))
    }
  }
  cover <- if (length(out)) do.call(rbind, out) else NULL

  ## how far the refits wander from the partition they were generated from
  stab <- cbind(key, data.frame(
    B_used = if (is.null(boot$B_used)) CFG$B[1] else boot$B_used,
    ari_mean = mean(inf$stability$ARI, na.rm = TRUE),
    ari_min  = min(inf$stability$ARI, na.rm = TRUE),
    ari_share1 = mean(inf$stability$ARI >= 0.999, na.rm = TRUE),
    secs = secs, stringsAsFactors = FALSE))

  list(coverage = cover, stability = stab)
}

## ---------------------------------------------------------------------------
append_csv <- function(df, path) {
  if (is.null(df)) return(invisible(NULL))
  utils::write.table(df, path, sep = ",", row.names = FALSE,
                     col.names = !file.exists(path), append = file.exists(path))
}
done <- if (file.exists(CSA)) utils::read.csv(CSA, stringsAsFactors = FALSE) else NULL
key  <- function(x) paste(x$n, x$TN, x$K, round(x$omega, 6), x$id, x$rep, sep = "|")

cat(sprintf("%4s %5s %2s %5s %-4s %5s | %8s %8s %8s\n",
            "n", "T", "K", "omega", "scen", "rep", "secs", "ariBoot", "share1"))

for (i in seq_len(nrow(cells))) {
  for (r in seq_len(CFG$nrep)) {
    cell <- cells[i, , drop = FALSE]; cell$rep <- r
    if (!is.null(done) && key(cell) %in% key(done)) next
    out <- try(run_one(cell, r), silent = TRUE)
    if (inherits(out, "try-error")) {
      cat(sprintf("%4d %5d %2d %5.2f %-4s %5d | FAILED: %s\n", cell$n, cell$TN,
                  cell$K, cell$omega, cell$id, r,
                  conditionMessage(attr(out, "condition"))))
      next
    }
    append_csv(out$coverage, CSV)
    append_csv(out$stability, CSA)
    s <- out$stability
    cat(sprintf("%4d %5d %2d %5.2f %-4s %5d | %8.1f %8.3f %8.3f\n",
                s$n, s$TN, s$K, s$omega, s$id, s$rep, s$secs,
                s$ari_mean, s$ari_share1))
    utils::flush.console()
  }
}
