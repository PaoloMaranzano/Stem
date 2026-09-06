## ---------------------------------------------------------------------------
## Monte Carlo driver for the simulation study.
##
## One replication is the whole procedure a user would run: build the model
## object, fit the (k, phi) grid with SCSTEM_Infocrit(), and apply the two-step
## rule with SCSTEM_Select(). What is recorded is the selected (k, phi), the
## Adjusted Rand Index of the selected partition against the truth, and the wall
## time.
##
## EVERY FACTOR OF THE DESIGN IS A COMMAND-LINE OPTION, so a run can be narrowed
## to the margins that will actually be discussed. The values below are the
## defaults; each is a comma-separated list on the command line:
##
##   Rscript dev/paper/07-simulation.R                       the whole design
##   Rscript dev/paper/07-simulation.R --balance=balanced    balanced cells only
##   Rscript dev/paper/07-simulation.R --K=1,2 --omega=0,1 --nrep=25
##   Rscript dev/paper/07-simulation.R --n=40,100 --TN=120,365 --tag=core
##   Rscript dev/paper/07-simulation.R --scenario=S0,S4 --nrep=1 --tag=timing
##
## Fractions are allowed, so --omega=0,1/3,2/3,1 works. The script runs from any
## working directory: it finds the repository from its own location, and the
## output directory can be moved with the STEM_CACHE environment variable. That
## is what makes it usable on a virtual machine.
##
## Results are appended to a CSV cell by cell, so a run that is interrupted
## keeps everything it has already produced, and a run that is resumed skips the
## cells already in the file. `--tag` names the file, so several designs can sit
## side by side.
##
## THE COST OF A CELL is driven by n and T, and only weakly by K, by the overlap
## and by the scenario: those change how many EM and ICM iterations are needed,
## not the cost of one.
## ---------------------------------------------------------------------------

local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  here <- if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
  source(file.path(here, "00-setup.R"), chdir = TRUE)
  source(file.path(here, "06-dgp.R"), chdir = TRUE)
}, envir = globalenv())

stem_load()

dims <- dgp_dims()
scen <- dgp_scenarios()

CFG <- stem_config(list(
  ## the design
  n        = dims$n,
  TN       = dims$TN,
  K        = dims$K,
  omega    = dims$omega,
  balance  = dims$balance,
  scenario = scen$id,
  nrep     = 100L,
  ## what the estimator searches over
  k_grid   = 1:4,
  phi_grid = c(0, 0.5, 1),
  ## regularisation of the regression coefficients; 0, 0 is the unpenalised
  ## estimator, which is what the paper reports
  alpha    = 0,
  lambda   = 0,
  ## bookkeeping
  tag      = "full",
  seed0    = 1000L
))
stem_print_config(CFG)

CSV <- file.path(stem_cache_dir(), sprintf("07-simulation-%s.csv", CFG$tag))

## ---------------------------------------------------------------------------
## The cells
##
## K = 1 has a single regime, so the overlap, the balance and the scenario are
## all vacuous there: it enters once per (n, T), with the dispersion inflated so
## that the map is the same size as at K = 2 -- see dgp_locations(). Cells whose
## imbalance cannot be realised with regimes of at least N_MIN units are dropped
## rather than silently rebalanced.
## ---------------------------------------------------------------------------
cells <- NULL
if (1L %in% CFG$K) {
  cells <- expand.grid(n = CFG$n, TN = CFG$TN, K = 1L, omega = CFG$omega[1],
                       id = CFG$scenario[1], balance = CFG$balance[1],
                       stringsAsFactors = FALSE)
}
Kmulti <- setdiff(CFG$K, 1L)
if (length(Kmulti)) {
  cells <- rbind(cells,
    expand.grid(n = CFG$n, TN = CFG$TN, K = Kmulti, omega = CFG$omega,
                id = CFG$scenario, balance = CFG$balance,
                stringsAsFactors = FALSE))
}
cells <- cells[mapply(dgp_feasible, cells$n, cells$K, cells$balance), ]
## cheapest first, so that an interrupted run still covers the design
cells <- cells[order(cells$n * cells$TN), ]
rownames(cells) <- NULL

cat(sprintf("%d cells x %d replications, writing to\n  %s\n\n",
            nrow(cells), CFG$nrep, CSV))

## ---------------------------------------------------------------------------
## One replication
## ---------------------------------------------------------------------------
run_one <- function(cell, rep) {

  row <- scen[scen$id == cell$id, , drop = FALSE]
  dat <- dgp_draw(cell$n, cell$TN, cell$K, cell$omega, row, rep = rep,
                  balance = cell$balance)

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
  ic <- SCSTEM_Infocrit(mod, k_grid = CFG$k_grid, phi_grid = CFG$phi_grid,
                        distance = "geo", verbose = FALSE,
                        alpha = CFG$alpha[1], lambda = CFG$lambda[1])
  sel <- SCSTEM_Select(ic)
  secs <- proc.time()[["elapsed"]] - t0

  gsel <- if (is.null(sel$fit)) NULL else sel$fit$group
  ari  <- if (is.null(gsel)) NA_real_ else scstem_ari(gsel, dat$labels)

  data.frame(n = cell$n, TN = cell$TN, K = cell$K, omega = cell$omega,
             id = cell$id, balance = cell$balance,
             alpha = CFG$alpha[1], lambda = CFG$lambda[1], rep = rep,
             k_hat = sel$k_selected, phi_hat = sel$phi_selected, ari = ari,
             nconf = nrow(ic$table), nfail = nrow(ic$failed),
             secs = secs, stringsAsFactors = FALSE)
}

## ---------------------------------------------------------------------------
## The loop, appending as it goes
## ---------------------------------------------------------------------------
done <- if (file.exists(CSV)) utils::read.csv(CSV, stringsAsFactors = FALSE) else NULL
key  <- function(x) paste(x$n, x$TN, x$K, round(x$omega, 6), x$id, x$balance,
                          x$rep, sep = "|")

cat(sprintf("%4s %5s %2s %5s %-4s %-10s %5s | %8s %6s %6s %6s\n",
            "n", "T", "K", "omega", "scen", "balance", "rep",
            "secs", "k_hat", "phi", "ARI"))

for (i in seq_len(nrow(cells))) {
  for (r in seq_len(CFG$nrep)) {
    cell <- cells[i, , drop = FALSE]
    cell$rep <- r
    if (!is.null(done) && key(cell) %in% key(done)) next

    out <- tryCatch(run_one(cell, r), error = function(e) {
      cbind(cell[c("n", "TN", "K", "omega", "id", "balance", "rep")],
            alpha = CFG$alpha[1], lambda = CFG$lambda[1],
            k_hat = NA_integer_, phi_hat = NA_real_, ari = NA_real_,
            nconf = NA_integer_, nfail = NA_integer_, secs = NA_real_)
    })

    utils::write.table(out, CSV, sep = ",", row.names = FALSE,
                       col.names = !file.exists(CSV), append = file.exists(CSV))
    cat(sprintf("%4d %5d %2d %5.2f %-4s %-10s %5d | %8.1f %6s %6s %6s\n",
                out$n, out$TN, out$K, out$omega, out$id, out$balance,
                out$rep, out$secs, format(out$k_hat),
                format(round(out$phi_hat, 2)), format(round(out$ari, 3))))
    utils::flush.console()
  }
}
