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
  ## what the estimator searches over. phi_ref is the penalty at which parameter
  ## recovery is read off, and has to be a point of phi_grid.
  k_grid   = 1:4,
  phi_grid = c(0, 0.5, 1),
  phi_ref  = 0.5,
  ## regularisation of the regression coefficients; 0, 0 is the unpenalised
  ## estimator, which is what the paper reports
  alpha    = 0,
  lambda   = 0,
  ## bookkeeping
  tag      = "full",
  seed0    = 1000L
))
stem_print_config(CFG)

CSV     <- file.path(stem_cache_dir(), sprintf("07-simulation-%s.csv", CFG$tag))
CSV_PAR <- file.path(stem_cache_dir(), sprintf("07-simulation-%s-params.csv", CFG$tag))

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

  ## ------------------------------------------------------------------------
  ## Two questions, two fits, deliberately separated.
  ##
  ## SELECTION asks whether the rule finds the truth, and is read off the fit
  ## the rule chose. RECOVERY asks whether the estimator gets the parameters
  ## right when it is told the truth, and is read off the fit at the TRUE
  ## number of regimes. Measuring recovery on the selected fit would confound
  ## the two: a poor estimate would be indistinguishable from a poor selection.
  ## Both fits are already in the grid, so neither costs an extra run.
  ## ------------------------------------------------------------------------
  grab <- function(kk, pp) {
    j <- which(ic$table$k == kk & abs(ic$table$phi - pp) < 1e-8)
    if (!length(j)) NULL else ic$fits[[j[1]]]
  }
  fit_sel  <- sel$fit
  ## with one regime the penalty is vacuous, and the grid holds k = 1 only at
  ## the first value of phi
  fit_true <- grab(cell$K, if (cell$K == 1L) CFG$phi_grid[1] else CFG$phi_ref[1])
  fit_pool <- grab(1L, CFG$phi_grid[1])

  ## Clustering accuracy. The ARI is invariant to label switching; the share of
  ## correctly assigned locations is not, so the labels are first relocated onto
  ## the truth by the majority rule, exactly as the bootstrap does.
  acc <- function(fit) {
    if (is.null(fit) || is.null(fit$group)) return(c(ari = NA_real_, share = NA_real_))
    g  <- fit$group
    kk <- max(max(g), cell$K)
    map <- scstem_align_labels(reference = dat$labels, refit = g, K = kk)
    c(ari   = scstem_ari(g, dat$labels),
      share = mean(map[g] == dat$labels, na.rm = TRUE))
  }
  a_sel  <- acc(fit_sel)
  a_true <- acc(fit_true)

  ## Predictive accuracy against the CONDITIONAL MEAN, not against z: the noise
  ## is irreducible, and scoring against z would compress every comparison
  ## towards one.
  ##
  ## NOTE. This does NOT use SCSTEM_Fitted(), and the distinction matters.
  ## SCSTEM_Fitted() returns E[z | observed], which fills the gaps and therefore
  ## returns z itself wherever z was observed -- on complete data it is the data,
  ## and scoring it against anything would be meaningless. What is wanted here is
  ## the SIGNAL the model fits,
  ##
  ##     muhat_ti = x_ti' betahat_g + K_i yhat_t^(g) ,
  ##
  ## which is a different object. The expression below is specific to the design
  ## of this study -- intercept plus one covariate, p = 1, K = 1 -- and is
  ## written out rather than taken from the package, which has no function for
  ## it yet.
  signal <- function(fit) {
    if (is.null(fit) || is.null(fit$fit_list)) return(NULL)
    mh <- matrix(NA_real_, cell$TN, cell$n)
    for (g in seq_along(fit$fit_list)) {
      f <- fit$fit_list[[g]]
      if (is.null(f) || is.null(f$estimates$phi.hat)) next
      idx <- which(fit$group == g)
      if (!length(idx)) next
      b <- as.numeric(f$estimates$phi.hat$beta)
      y <- as.numeric(f$estimates$y.smoothed)
      for (i in idx) mh[, i] <- b[1] + b[2] * dat$x[, i] + y
    }
    mh
  }
  rmse <- function(fit) {
    mh <- signal(fit)
    if (is.null(mh) || all(is.na(mh))) return(NA_real_)
    sqrt(mean((mh - dat$mu)^2, na.rm = TRUE))
  }
  r_true <- rmse(fit_true)
  r_pool <- rmse(fit_pool)

  key <- data.frame(n = cell$n, TN = cell$TN, K = cell$K, omega = cell$omega,
                    id = cell$id, balance = cell$balance,
                    alpha = CFG$alpha[1], lambda = CFG$lambda[1], rep = rep,
                    stringsAsFactors = FALSE)

  summ <- cbind(key, data.frame(
    k_hat = sel$k_selected, phi_hat = sel$phi_selected,
    k_correct = as.integer(sel$k_selected == cell$K),
    ari_sel = a_sel[["ari"]],   share_sel = a_sel[["share"]],
    ari_true = a_true[["ari"]], share_true = a_true[["share"]],
    rmse_true = r_true, rmse_pooled = r_pool, rmse_ratio = r_true / r_pool,
    nconf = nrow(ic$table), nfail = nrow(ic$failed),
    secs = secs, stringsAsFactors = FALSE))

  ## ------------------------------------------------------------------------
  ## Parameter recovery, regime by regime, at the true number of regimes and
  ## after relocating the estimated labels onto the true ones. What is stored is
  ## the estimate beside the truth, one row each: bias, RMSE and coverage are
  ## Monte Carlo summaries of this file, not quantities a single replication
  ## could compute.
  ## ------------------------------------------------------------------------
  tru <- dgp_truth(dat$psi)
  est <- rep(NA_real_, nrow(tru))
  if (!is.null(fit_true) && !is.null(fit_true$phi_hat)) {
    ph  <- fit_true$phi_hat
    map <- scstem_align_labels(reference = dat$labels,
                               refit = fit_true$group, K = cell$K)
    est <- vapply(seq_len(nrow(tru)), function(i) {
      g_fit <- which(map == tru$regime[i])   # the fitted regime playing that role
      p     <- tru$parameter[i]
      if (!length(g_fit) || !(p %in% colnames(ph))) NA_real_ else ph[g_fit[1], p]
    }, numeric(1))
  }
  params <- cbind(key[rep(1L, nrow(tru)), ], tru, estimate = est)
  rownames(params) <- NULL

  list(summary = summ, params = params)
}

## ---------------------------------------------------------------------------
## The loop, appending as it goes
## ---------------------------------------------------------------------------
## Two files, because the two have different shapes: one row per replication for
## the selection and accuracy measures, one row per replication, regime and
## parameter for the recovery. Both are appended cell by cell, so an interrupted
## run keeps what it has and a resumed one skips it.
done <- if (file.exists(CSV)) utils::read.csv(CSV, stringsAsFactors = FALSE) else NULL
key  <- function(x) paste(x$n, x$TN, x$K, round(x$omega, 6), x$id, x$balance,
                          x$rep, sep = "|")

append_csv <- function(df, path) {
  utils::write.table(df, path, sep = ",", row.names = FALSE,
                     col.names = !file.exists(path), append = file.exists(path))
}

cat(sprintf("%4s %5s %2s %5s %-4s %-10s %5s | %7s %5s %6s %6s %6s\n",
            "n", "T", "K", "omega", "scen", "balance", "rep",
            "secs", "k_hat", "ARI", "share", "rmseR"))

for (i in seq_len(nrow(cells))) {
  for (r in seq_len(CFG$nrep)) {
    cell <- cells[i, , drop = FALSE]
    cell$rep <- r
    if (!is.null(done) && key(cell) %in% key(done)) next

    out <- tryCatch(run_one(cell, r), error = function(e) {
      k0 <- cbind(cell[c("n", "TN", "K", "omega", "id", "balance", "rep")],
                  alpha = CFG$alpha[1], lambda = CFG$lambda[1])
      list(summary = cbind(k0, data.frame(
             k_hat = NA_integer_, phi_hat = NA_real_, k_correct = NA_integer_,
             ari_sel = NA_real_, share_sel = NA_real_,
             ari_true = NA_real_, share_true = NA_real_,
             rmse_true = NA_real_, rmse_pooled = NA_real_, rmse_ratio = NA_real_,
             nconf = NA_integer_, nfail = NA_integer_, secs = NA_real_,
             stringsAsFactors = FALSE)),
           params = NULL)
    })

    append_csv(out$summary, CSV)
    if (!is.null(out$params)) append_csv(out$params, CSV_PAR)

    s <- out$summary
    cat(sprintf("%4d %5d %2d %5.2f %-4s %-10s %5d | %7.1f %5s %6s %6s %6s\n",
                s$n, s$TN, s$K, s$omega, s$id, s$balance, s$rep, s$secs,
                format(s$k_hat), format(round(s$ari_true, 3)),
                format(round(s$share_true, 3)), format(round(s$rmse_ratio, 3))))
    utils::flush.console()
  }
}
