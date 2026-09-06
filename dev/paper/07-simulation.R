## ---------------------------------------------------------------------------
## Monte Carlo driver for the simulation study.
##
## One replication is the whole procedure a user would run: build the model
## object, fit the (k, phi) grid with SCSTEM_Infocrit(), and apply the two-step
## rule with SCSTEM_Select(). What is recorded is the selected (k, phi), the
## Adjusted Rand Index of the selected partition against the truth, the
## log-likelihood, and the wall time.
##
##   Rscript dev/paper/07-simulation.R timing   one replication per cell, to
##                                              size the Monte Carlo
##   Rscript dev/paper/07-simulation.R full     NREP replications
##
## Results are appended to a CSV cell by cell, so a run that is interrupted
## keeps everything it has already produced, and a run that is resumed skips
## the cells already in the file.
##
## THE COST OF A CELL is driven by n and T, and only very weakly by K, by the
## overlap d and by the scenario: those change how many EM and ICM iterations
## are needed, not the cost of one. The timing design is therefore split in two
## blocks -- the full n by T cross at one representative (K, d, scenario), and
## a block that varies K, d and the scenario at a fixed, middling (n, T) to
## confirm that the second group of factors does not move the cost.
## ---------------------------------------------------------------------------

HERE <- "C:/Users/paulm/OneDrive/Documenti/GitHub/Stem"
OUT  <- file.path(HERE, "dev", "paper", "cache")
mode <- (function(a) if (length(a)) a[1] else "timing")(commandArgs(TRUE))

source(file.path(HERE, "dev", "paper", "06-dgp.R"))
suppressMessages(pkgload::load_all(HERE, quiet = TRUE))

K_GRID   <- 1:4
PHI_GRID <- c(0, 0.5, 1)
NREP     <- if (mode == "timing") 1L else 100L
CSV      <- file.path(OUT, sprintf("07-simulation-%s.csv", mode))

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------------------
## The cells
## ---------------------------------------------------------------------------
dims <- dgp_dims()
scen <- dgp_scenarios()

cells_timing <- rbind(
  ## block A: the n by T cross, at K = 2, moderate overlap, all contrasts
  expand.grid(n = dims$n, TN = dims$TN, K = 2L, d = 2/3, id = "S4",
              block = "A", stringsAsFactors = FALSE),
  ## block B: K, overlap and scenario at a middling (n, T)
  expand.grid(n = 60L, TN = 365L, K = dims$K, d = c(0, 2/3),
              id = c("S0", "S1b", "S3b", "S4"),
              block = "B", stringsAsFactors = FALSE)
)
## cheapest first, so that an interrupted run still covers the design
cells_timing <- cells_timing[order(cells_timing$n * cells_timing$TN), ]

cells_full <- expand.grid(n = dims$n, TN = dims$TN, K = dims$K, d = dims$d,
                          id = scen$id, block = "full",
                          stringsAsFactors = FALSE)

cells <- if (mode == "timing") cells_timing else cells_full
rownames(cells) <- NULL

## ---------------------------------------------------------------------------
## One replication
## ---------------------------------------------------------------------------
run_one <- function(cell, rep) {

  row <- scen[scen$id == cell$id, , drop = FALSE]
  dat <- dgp_draw(cell$n, cell$TN, cell$K, cell$d, row, rep = rep)

  ## starting values from the pooled OLS fit, as a user would
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
  ic <- SCSTEM_Infocrit(mod, k_grid = K_GRID, phi_grid = PHI_GRID,
                        distance = "geo", verbose = FALSE)
  sel <- SCSTEM_Select(ic)
  secs <- proc.time()[["elapsed"]] - t0

  gsel <- if (is.null(sel$fit)) NULL else sel$fit$group
  ari  <- if (is.null(gsel)) NA_real_ else scstem_ari(gsel, dat$labels)

  data.frame(block = cell$block, n = cell$n, TN = cell$TN, K = cell$K,
             d = cell$d, id = cell$id, rep = rep,
             k_hat = sel$k_selected, phi_hat = sel$phi_selected, ari = ari,
             nconf = nrow(ic$table), nfail = nrow(ic$failed),
             secs = secs, stringsAsFactors = FALSE)
}

## ---------------------------------------------------------------------------
## The loop, appending as it goes
## ---------------------------------------------------------------------------
done <- if (file.exists(CSV)) utils::read.csv(CSV, stringsAsFactors = FALSE) else NULL
key  <- function(x) paste(x$block, x$n, x$TN, x$K, x$d, x$id, x$rep, sep = "|")

cat(sprintf("%-5s %4s %5s %2s %5s %-4s %4s | %8s %6s %6s %6s\n",
            "block", "n", "T", "K", "d", "scen", "rep",
            "secs", "k_hat", "phi", "ARI"))

for (i in seq_len(nrow(cells))) {
  for (r in seq_len(NREP)) {
    cell <- cells[i, , drop = FALSE]
    cell$rep <- r
    if (!is.null(done) && key(cell) %in% key(done)) next

    out <- tryCatch(run_one(cell, r), error = function(e) {
      cbind(cell[c("block", "n", "TN", "K", "d", "id", "rep")],
            k_hat = NA_integer_, phi_hat = NA_real_, ari = NA_real_,
            nconf = NA_integer_, nfail = NA_integer_, secs = NA_real_)
    })

    utils::write.table(out, CSV, sep = ",", row.names = FALSE,
                       col.names = !file.exists(CSV), append = file.exists(CSV))
    cat(sprintf("%-5s %4d %5d %2d %5.2f %-4s %4d | %8.1f %6s %6s %6s\n",
                out$block, out$n, out$TN, out$K, out$d, out$id, out$rep,
                out$secs, format(out$k_hat), format(round(out$phi_hat, 2)),
                format(round(out$ari, 3))))
    utils::flush.console()
  }
}
