## ===========================================================================
## SC-STEM fuel application: THE BOOTSTRAP LIKELIHOOD-RATIO TEST
##
## Not part of the paper or of the replication material (decision of
## 2026-10-20): the results are kept for a possible request of a referee.
##
## The two bootstrap tests of leadership, this one and fuels-wald-bootstrap.R,
## are post-estimation inference on the selected fit of a model: the same
## hypotheses (the blocks of lags of fu_tests(): in a regime, the lags of the
## neighbor have zero coefficients), the same regimes, the same tables. Both
## propagate the uncertainty of the estimated partition, by re-estimating
## SC-STEM at (K*, phi*) on every bootstrap sample and matching the regimes of
## the sample to the observed ones by the majority rule. They differ in the
## statistic and in what the bootstrap is asked for:
##
##                 Wald                              LR
##   statistic     W_k = b_k' V_k^-1 b_k             LR_k = 2 (l_full,k - l_restr,k)
##   draws from    the fitted model                  the restricted model (the null)
##   the bootstrap the covariance V_k of the         the distribution of LR_k
##   gives         coefficients and, centered at      under the null
##                 the estimate, the distribution
##                 of W_k under the null
##   p_chisq       W_k against chi-square(q)         LR_k against chi-square(q)
##   p_boot        share of the W*_k >= W_k          share of the LR*_k >= LR_k
##
## Row "all" is the hypothesis in all the regimes at once, on q K degrees of
## freedom: the coefficients of the K regimes stacked, with their joint
## bootstrap covariance (Wald); the sum of the K statistics (LR).
##
## THE LR TEST, for one model and one hypothesis (q lags):
##   1. on the data: in every regime k of the selected fit (stage "grid", K*
##      and phi*), the model without the q lags is fitted on the regime's
##      pumps, and LR_k = 2 (l_full,k - l_restr,k);
##   2. the bootstrap: B samples drawn regime by regime from the restricted
##      fits of step 1 (Stem::STEM_Simulation, the generator of
##      Stem::SCSTEM_Bootstrap, with the observed covariates, coordinates and
##      missing values), SC-STEM re-estimated on each at (K*, phi*) with the
##      settings of the selected fit, the regimes matched by the majority rule
##      of Stem::SCSTEM_BootInference. The null holds in all the regimes at
##      once, so that one bootstrap serves the K + 1 tests of a hypothesis;
##   3. on every sample, LR*_kb as in step 1 on the regimes of the sample;
##   4. p_chisq against the chi-square on q degrees of freedom; p_boot against
##      the draws LR*_kb.
## The likelihood of a regime may have several local maxima
## (fuels-check-optima.R), and a statistic built on two maximized likelihoods
## must not depend on where the EM algorithm stops. In every regime, on the
## data and on every sample: the restricted model is fitted from two starts
## (the full estimates without the q lags, and least squares) and the better
## kept; the full model is refitted from the restricted solution with the q
## lags at zero, where its likelihood equals the restricted one, so that the
## EM algorithm can only climb from there; the full likelihood is the larger of
## that refit and of the regime's fit within SC-STEM. LR is never negative.
##
## Only exported functions of Stem are used (SCSTEM_Select, SCSTEM_Estimation,
## STEM_Model, STEM_Estimation, STEM_Simulation): the test is specific to the
## application and lives in the replication material, not in the package.
##
## Edit THE SETUP and press Source in RStudio, or from a terminal
##
##     Rscript fuels-lr-bootstrap.R --settings=1,2,3 --cores=4
##
## (command-line options override THE SETUP: --settings, --cities, --fuels,
## --blocks, --B, --cores, --root, --out). Requirements: the stages grid and
## final of the main scripts, whose cached grids and summaries are read; Stem
## at the pinned commit. The bootstrap of a model is saved after every round
## of draws in <out>/<case>/nn/lr-boot/, and a stopped run resumes where it
## was. Outputs: <out>/<case>/nn/final/<model>/test_lr.csv and
## <out>/<case>/nn/final/test_lr_summary.csv.
##
## This is part of the REPLICATION MATERIAL of the paper, not of the Stem
## package.
## ===========================================================================


## ===========================================================================
## THE SETUP. This is the part meant to be edited.
## ===========================================================================
RUN <- list(
  ## the settings: 1, 2, 3, or several
  settings = 1:3,
  ## the cities: "all", or codes among RM MI NA TO PA BA CT BO ME FI VE
  cities   = "all",
  ## the fuel of Y: "g" gasoline, "d" diesel
  fuels    = c("g", "d"),
  ## the hypotheses (blocks of lags of fu_tests()): "nb_both" (the lags of the
  ## neighbor, both fuels: the leadership test of the paper), "nb", "nb_other",
  ## "own_other"
  blocks   = "nb_both",
  ## bootstrap draws (stage "final" of the main scripts uses 100)
  B        = 100L,
  ## R processes the draws of a model are shared among
  cores    = 4L,
  ## the folder with the data and the outputs of the main scripts: "" for
  ## fuels/ beside the scripts
  root     = "",
  ## where this script writes: "" for the same folder as root
  out      = ""
)


## ===========================================================================
## Machinery: nothing below needs editing
## ===========================================================================
suppressPackageStartupMessages(library(data.table))
HERE <- local({
  a <- commandArgs(trailingOnly = FALSE); m <- grep("^--file=", a, value = TRUE)
  f <- if (length(m)) normalizePath(sub("^--file=", "", m[1]), winslash = "/") else NULL
  if (is.null(f)) for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) { f <- normalizePath(of, winslash = "/"); break }
  }
  if (is.null(f)) normalizePath(getwd(), winslash = "/") else dirname(f)
})
## kept outside the replication material (dev/lr-bootstrap): the functions
## and the fuels folder are those of dev/replication
if (!file.exists(file.path(HERE, "fuels-functions.R"))) HERE <- file.path(dirname(HERE), "replication")
source(file.path(HERE, "fuels-functions.R"))
RUN <- fu_bt_setup(RUN, HERE)
SHA7 <- fu_require_stem()

## --- the statistic --------------------------------------------------------------
## rows of the pumps idx in the stacked covariates (pump by pump, weeks within)
lr_rows <- function(idx, Tn) unlist(lapply(idx, function(i) ((i - 1L) * Tn + 1L):(i * Tn)))
## a starting value for STEM_Model(): the estimates ph where usable, the
## coefficients beta, and phi0 for the rest (as the starts of SCSTEM_Estimation)
lr_phi <- function(ph, phi0, beta) {
  p <- nrow(phi0$G); out <- phi0
  beta <- as.numeric(beta)
  if (length(beta) == length(phi0$beta) && all(is.finite(beta))) out$beta <- matrix(beta, ncol = 1)
  take <- function(x, ref, shape) if (is.null(x) || length(x) != length(ref) || any(!is.finite(x))) ref else shape(x)
  if (!is.null(ph)) {
    out$G <- take(ph$G, phi0$G, function(x) matrix(x, p, p))
    out$Sigmaeta <- take(ph$Sigmaeta, phi0$Sigmaeta, function(x) matrix(x, p, p))
    out$m0 <- take(ph$m0, phi0$m0, function(x) matrix(x, p, 1))
    for (nm in c("sigma2eps", "sigma2omega", "theta")) {
      x <- ph[[nm]]
      if (length(x) == 1L && is.finite(x) && x > 0) out[[nm]] <- as.numeric(x)
    }
  }
  out
}
## least-squares coefficients of a regime
lr_ols <- function(z, X) {
  y <- as.vector(z); ok <- is.finite(y)
  b <- tryCatch(stats::lm.fit(X[ok, , drop = FALSE], y[ok])$coefficients, error = function(e) NULL)
  if (is.null(b)) rep(NA_real_, ncol(X)) else { b[!is.finite(b)] <- 0; b }
}
## one STEM fit of a regime from a start; when it fails, NULL with the reason
## as its attribute "error"
lr_fit <- function(z, X, coords, start, ctl) {
  fail <- function(m) structure(list(), error = m)
  mod <- tryCatch(Stem::STEM_Model(z = z, covariates = X, coordinates = coords, phi = start,
                                   A = matrix(1, ncol(z), 1)), error = function(e) fail(conditionMessage(e)))
  if (!inherits(mod, "STEM_Model")) return(fail(paste("STEM_Model:", attr(mod, "error"))))
  f <- tryCatch(suppressWarnings(Stem::STEM_Estimation(mod, control = ctl, distance = "geo")),
                error = function(e) fail(conditionMessage(e)))
  if (!is.null(attr(f, "error"))) return(fail(paste("STEM_Estimation:", attr(f, "error"))))
  if (is.null(f$estimates$phi.hat) || !is.finite(f$estimates$loglik)) return(fail("no finite log-likelihood"))
  f
}
lr_ll <- function(f) if (!length(f) || !is.null(attr(f, "error"))) NA_real_ else as.numeric(f$estimates$loglik)
## The statistic of one regime, with the nested starts. ph_full, ll_full: the
## regime's fit within SC-STEM; keep: the columns of the restricted model;
## phi0_full, phi0_rest: the starting values of the two models (the skeleton
## of fu_model()).
lr_regime <- function(z, X, coords, keep, ph_full, ll_full, phi0_full, phi0_rest, ctl) {
  Xr <- X[, keep, drop = FALSE]
  fa <- lr_fit(z, Xr, coords, lr_phi(ph_full, phi0_rest, as.numeric(ph_full$beta)[keep]), ctl)
  fb <- lr_fit(z, Xr, coords, lr_phi(NULL, phi0_rest, lr_ols(z, Xr)), ctl)
  la <- lr_ll(fa); lb <- lr_ll(fb)
  if (!is.finite(la) && !is.finite(lb))
    return(list(stat = NA_real_, ll_full = ll_full, ll_full_nested = NA_real_, ll_rest = NA_real_, phi_rest = NULL,
                error = paste(unique(c(attr(fa, "error"), attr(fb, "error"))), collapse = "; ")))
  rest <- if (!is.finite(lb) || (is.finite(la) && la >= lb)) fa else fb
  ph_r <- rest$estimates$phi.hat
  b_full <- numeric(ncol(X)); b_full[keep] <- as.numeric(ph_r$beta)
  fc <- lr_fit(z, X, coords, lr_phi(ph_r, phi0_full, b_full), ctl)
  l_full <- max(ll_full, lr_ll(fc), na.rm = TRUE)
  list(stat = 2 * (l_full - lr_ll(rest)), ll_full = ll_full, ll_full_nested = lr_ll(fc), ll_rest = lr_ll(rest),
       phi_rest = lr_phi(ph_r, phi0_rest, as.numeric(ph_r$beta)))
}
## every regime of an SC-STEM fit (fit_list, idx_g, loglik_g) on data z
lr_all_regimes <- function(fit, z, X, coords, Tn, keep, phi0_full, phi0_rest, ctl) {
  lapply(seq_along(fit$idx_g), function(g) {
    idx <- fit$idx_g[[g]]
    lr_regime(z[, idx, drop = FALSE], X[lr_rows(idx, Tn), , drop = FALSE], coords[idx, , drop = FALSE],
              keep, fit$fit_list[[g]]$estimates$phi.hat, fit$loglik_g[g], phi0_full, phi0_rest, ctl)
  })
}

## --- the bootstrap ----------------------------------------------------------------
## the majority rule of Stem::SCSTEM_BootInference (scstem_align_labels): the
## observed label of every regime of a refit
lr_align <- function(reference, refit, K) {
  tab <- matrix(0L, K, K)
  for (g in seq_len(K)) for (h in seq_len(K)) tab[g, h] <- sum(refit == g & reference == h)
  map <- rep(NA_integer_, K); taken <- rep(FALSE, K)
  for (pos in order(tab, decreasing = TRUE)) {
    g <- ((pos - 1L) %% K) + 1L; h <- ((pos - 1L) %/% K) + 1L
    if (is.na(map[g]) && !taken[h] && tab[g, h] > 0) { map[g] <- h; taken[h] <- TRUE }
  }
  left <- which(is.na(map)); free <- which(!taken)
  if (length(left) && length(free)) map[left[seq_along(free)]] <- free
  map
}
## the adjusted Rand index, without packages
lr_ari <- function(a, b) {
  tab <- table(a, b); n <- sum(tab); c2 <- function(x) sum(x * (x - 1) / 2)
  s <- c2(tab); sa <- c2(rowSums(tab)); sb <- c2(colSums(tab)); e <- sa * sb / (n * (n - 1) / 2)
  if (((sa + sb) / 2 - e) == 0) 1 else (s - e) / ((sa + sb) / 2 - e)
}
## one draw: a sample from the null model, SC-STEM on it, the statistics of
## its regimes, matched to the observed ones
lr_draw <- function(seed, ctx) {
  set.seed(seed)
  z <- matrix(NA_real_, ctx$Tn, ctx$d)
  for (g in seq_along(ctx$idx_g)) {
    idx <- ctx$idx_g[[g]]
    gen <- Stem::STEM_Model(z = ctx$z[, idx, drop = FALSE],
                            covariates = ctx$X[lr_rows(idx, ctx$Tn), ctx$keep, drop = FALSE],
                            coordinates = ctx$coords[idx, , drop = FALSE], phi = ctx$phi_null[[g]],
                            A = matrix(1, length(idx), 1))
    z[, idx] <- as.matrix(Stem::STEM_Simulation(gen, distance = "geo"))
  }
  z[ctx$na] <- NA_real_
  mod <- Stem::STEM_Model(z = z, covariates = ctx$X, coordinates = ctx$coords, phi = ctx$phi0_full,
                          A = matrix(1, ctx$d, 1))
  fb <- tryCatch(suppressWarnings(do.call(Stem::SCSTEM_Estimation, c(list(StemModel = mod), ctx$refit_args))),
                 error = function(e) e)
  if (inherits(fb, "error")) return(list(ok = FALSE, message = conditionMessage(fb)))
  if (!all(fb$final_refit)) return(list(ok = FALSE, message = "a regime was not re-estimated"))
  r <- lr_all_regimes(fb, z, ctx$X, ctx$coords, ctx$Tn, ctx$keep, ctx$phi0_full, ctx$phi0_rest, ctx$ctl)
  K <- length(ctx$idx_g)
  stat <- rep(NA_real_, K); stat[lr_align(ctx$group, fb$group, K)] <- vapply(r, `[[`, 1, "stat")
  list(ok = all(is.finite(stat)), stat = stat, message = NA_character_, ari = lr_ari(ctx$group, fb$group))
}
## a draw on a worker, where ctx has been exported: a function of the global
## environment, so that sending it does not send the fit of the model along
lr_task <- function(s) lr_draw(s, ctx)
## The bootstrap of a hypothesis: the statistic on the data and the draws,
## saved after every round of 2 x cores draws and resumed from the saved ones.
lr_bootstrap <- function(tk, des, fit, blk, run) {
  K <- tk$K; a <- fit$input_args
  refit_args <- list(K = a$K, phi_penalty = a$phi_penalty, phi_scale = a$phi_scale, knn = a$knn,
                     distance = a$distance, init_method = a$init_method, label_update = a$label_update,
                     regularization = a$regularization, control = a$control,
                     min_cluster_size = a$min_cluster_size, enforce_min_size = a$enforce_min_size,
                     swap_pass = a$swap_pass, share2conv = a$share2conv, seed = a$seed, verbose = FALSE,
                     alpha = a$alpha, lambda = a$lambda, penalize = a$penalize, lambda_scale = a$lambda_scale,
                     lambda_by = a$lambda_by, latent = a$latent, spatial = a$spatial)
  refit_args <- refit_args[names(refit_args) %in% names(formals(Stem::SCSTEM_Estimation))]
  refit_args <- refit_args[!vapply(refit_args, is.null, logical(1))]
  ctl <- if (is.null(a$control)) Stem::STEM_control() else a$control
  phi0_full <- a$StemModel$skeleton$phi
  keep <- setdiff(seq_along(des$names), fu_tests(des)[[blk]])
  phi0_rest <- fu_model(fu_subset_design(des, keep))$skeleton$phi
  f_state <- file.path(tk$out_case, "nn", "lr-boot", sprintf("%s_%s.rds", tk$model, blk))
  dir.create(dirname(f_state), recursive = TRUE, showWarnings = FALSE)
  st <- if (file.exists(f_state)) readRDS(f_state) else NULL
  if (is.null(st)) {
    obs <- lr_all_regimes(fit, des$z, des$covariates, des$coords, des$Tn, keep, phi0_full, phi0_rest, ctl)
    bad <- which(vapply(obs, function(o) is.null(o$phi_rest), TRUE))
    if (length(bad))
      stop("the restricted model of regime ", bad[1], " could not be fitted on the data: ", obs[[bad[1]]]$error,
           call. = FALSE)
    st <- list(obs = obs, draws = list())
    saveRDS(st, f_state)
  }
  seeds <- local({ set.seed(fu_bt_seed(tk) + 7919L * match(blk, FU_BT$blocks))
    sample.int(.Machine$integer.max, 10000L) })[seq_len(run$B)]
  todo <- setdiff(seq_len(run$B), as.integer(names(st$draws)))
  if (length(todo)) {
    ctx <- list(z = des$z, X = des$covariates, coords = des$coords, Tn = des$Tn, d = des$d, na = is.na(des$z),
                keep = keep, idx_g = fit$idx_g, group = fit$group, phi_null = lapply(st$obs, `[[`, "phi_rest"),
                phi0_full = phi0_full, phi0_rest = phi0_rest, ctl = ctl, refit_args = refit_args)
    nc <- min(run$cores, length(todo)); t0 <- Sys.time()
    cl <- NULL
    if (nc > 1L) {
      cl <- parallel::makePSOCKcluster(nc)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterExport(cl, c("lr_rows", "lr_phi", "lr_ols", "lr_fit", "lr_ll", "lr_regime",
                                    "lr_all_regimes", "lr_align", "lr_ari", "lr_draw", "lr_task"))
      parallel::clusterExport(cl, "ctx", envir = environment())
    }
    for (ch in split(todo, ceiling(seq_along(todo) / (2L * nc)))) {
      res <- if (is.null(cl)) lapply(seeds[ch], lr_draw, ctx = ctx) else parallel::clusterApplyLB(cl, seeds[ch], lr_task)
      names(res) <- as.character(ch)
      st$draws[names(res)] <- res
      saveRDS(st, f_state)
      cat(sprintf("  %s  %s, %s: %d / %d draws (%d usable), %.0f min\n", format(Sys.time(), "%H:%M"), tk$model,
                  blk, length(st$draws), run$B, sum(vapply(st$draws, function(x) isTRUE(x$ok), TRUE)),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
      utils::flush.console()
    }
  }
  st
}

## --- one model, every hypothesis --------------------------------------------------
lr_model <- function(tk, run) {
  sel <- fu_bt_selected(tk, SHA7); des <- sel$des; fit <- sel$fit; K <- tk$K
  tests <- fu_tests(des)
  do.call(rbind, lapply(intersect(run$blocks, names(tests)), function(blk) {
    st <- lr_bootstrap(tk, des, fit, blk, run); q <- length(tests[[blk]])
    ok <- Filter(function(x) isTRUE(x$ok), st$draws[as.character(seq_len(run$B))])
    S <- if (length(ok)) do.call(rbind, lapply(ok, `[[`, "stat")) else matrix(NA_real_, 0, K)
    S <- cbind(S, rowSums(S))
    stat <- vapply(st$obs, `[[`, 1, "stat")
    fu_bt_table(tk, "LR", blk, stat = c(stat, sum(stat)), df = c(rep(q, K), K * q),
                pumps = tabulate(fit$group, K), null = lapply(seq_len(K + 1L), function(j) S[, j]), B = run$B,
                ari = if (length(ok)) stats::median(vapply(ok, `[[`, 1, "ari")) else NA_real_,
                extra = data.frame(ll_full_scstem = c(vapply(st$obs, `[[`, 1, "ll_full"), NA),
                                   ll_full_nested = c(vapply(st$obs, `[[`, 1, "ll_full_nested"), NA),
                                   ll_restricted = c(vapply(st$obs, `[[`, 1, "ll_rest"), NA)))
  }))
}

## --- the run ------------------------------------------------------------------
tasks <- fu_bt_tasks(RUN)
cat(sprintf("bootstrap LR test: %d models; blocks %s; B = %d on %d processes; Stem %s\n  outputs in %s\n\n",
            nrow(tasks), paste(RUN$blocks, collapse = ", "), RUN$B, RUN$cores, SHA7, RUN$out))
fu_bt_run(tasks, function(tk) lr_model(tk, RUN), "test_lr.csv", "bootstrap LR test")
