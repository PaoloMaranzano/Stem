## ===========================================================================
## SC-STEM: numerical verification of the theory of the paper.
##
## Every statement of Section 2 that can be checked by simulation is checked
## here, against its own prediction, with a verdict. Nothing is tuned to agree:
## the predictions are the formulas of the paper, evaluated at the parameters
## the data were generated with.
##
##     Rscript check-theory.R                        all parts
##     Rscript check-theory.R --parts=lemma,cond     some of them
##     Rscript check-theory.R --reps=100             fewer replications
##
## Parts
##   lemma      the expected gap between the exact within-regime log-density and
##              the pseudo-likelihood is -(T/2) log|R|
##   estimated  the same gap with estimated parameters and smoothed states
##   cond       the conditional score exceeds the marginal one by Delta_i(g) in
##              expectation when i belongs to g, and falls below it when not
##   decision   which score assigns locations best, with the parameters known
##   algorithm  the full procedure under the three scores
##   label      the label step: ascent, finite termination, and the two
##              thresholds of the penalty
##   ic         the in-sample criteria under the null as the network grows
##   ic-stem    the same mechanism in SC-STEM itself
##   tables     the LaTeX tables of the paper, from the CSVs of the parts above
##
## Needs an installed Stem and run-simulations.R beside it, from which the
## generator is taken. Writes to <here>/theory/.
##
## This is part of the REPLICATION MATERIAL of the paper, not of the package.
## ===========================================================================

here <- local({
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]), winslash = "/"))
  else normalizePath(getwd(), winslash = "/")
})
SIM_DEFINE_ONLY <- TRUE
source(file.path(here, "run-simulations.R"), local = globalenv())

cond_scores <- utils::getFromNamespace("scstem_cond_scores", "Stem")
neighbors   <- function(...) suppressWarnings(
  utils::getFromNamespace("scstem_neighbors", "Stem")(...))
ari <- utils::getFromNamespace("scstem_ari", "Stem")

TH <- sim_config(list(
  parts = c("lemma", "estimated", "cond", "decision", "algorithm", "label", "ic",
            "ic-stem", "tables"),
  reps  = 400L,
  alg_reps = 10L,
  out   = file.path(here, "theory")
))
OUT <- TH$out[1]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
OUT <- normalizePath(OUT, winslash = "/")
cat("theory checks, writing to", OUT, "\n\n")

verdicts <- list()
verdict <- function(part, claim, ok, detail) {
  verdicts[[length(verdicts) + 1L]] <<- data.frame(
    part = part, claim = claim, verdict = if (ok) "CONSISTENT" else "NOT CONSISTENT",
    detail = detail, stringsAsFactors = FALSE)
  cat(sprintf("  [%s] %s -- %s\n", if (ok) "ok" else "!!", claim, detail))
}
save_csv <- function(df, name) {
  utils::write.csv(df, file.path(OUT, name), row.names = FALSE)
}

## ---------------------------------------------------------------------------
## Helpers
## ---------------------------------------------------------------------------
## a within-regime covariance on real geography, from the design generator
cov_of <- function(coords, s2e, s2o, theta) {
  dm <- as.matrix(geodist::geodist(coords, measure = "geodesic"))
  S <- s2o * exp(-theta * dm); diag(S) <- s2e + s2o
  list(S = S, dm = dm)
}
rmvn <- function(Tn, S) {
  U <- chol(S)
  matrix(stats::rnorm(Tn * ncol(S)), Tn) %*% U
}
ldmvn_rows <- function(R, S) {                  # sum over rows of log N(r_t; 0, S)
  U <- chol(S)
  Z <- backsolve(U, t(R), transpose = TRUE)
  -0.5 * nrow(R) * ncol(R) * log(2 * pi) - nrow(R) * sum(log(diag(U))) - 0.5 * sum(Z^2)
}
ld_indep <- function(R, s2) {
  -0.5 * length(R) * log(2 * pi * s2) - 0.5 * sum(R^2) / s2
}
base <- dgp_base()
tot  <- base$sigma2eps + base$sigma2omega
## (nugget share, range in km): from the baseline to a strongly correlated field
COVS <- data.frame(
  label = c("baseline", "low nugget", "short range", "strong"),
  share = c(base$sigma2eps / tot, 0.30, 0.30, 0.05),
  range = c(123, 123, 15, 200))


## ===========================================================================
## LEMMA. E[l_g - l~_g] = -(T/2) log|R_g| = T KL(N(0,Sigma) || N(0,D))
## ===========================================================================
if ("lemma" %in% TH$parts) {
  cat("LEMMA: the expected gap of the pseudo-likelihood\n")
  set.seed(101)
  Tn <- 200L
  res <- NULL
  for (ng in c(5L, 15L, 40L)) {
    loc <- dgp_locations(3L * ng, 3L, 1, seed = 7)
    co  <- loc$coords[loc$labels == 1L, , drop = FALSE][seq_len(ng), , drop = FALSE]
    for (j in seq_len(nrow(COVS))) {
      cv <- cov_of(co, tot * COVS$share[j], tot * (1 - COVS$share[j]),
                   1 / (COVS$range[j] * 1000))
      S <- cv$S; s2 <- S[1, 1]
      Rm <- stats::cov2cor(S)
      pred <- -0.5 * Tn * as.numeric(determinant(Rm, logarithm = TRUE)$modulus)
      A <- diag(ng) - Rm                                   # I - D^{-1} Sigma
      sd_pred <- sqrt(0.5 * Tn * sum(A * t(A)))            # sd of the gap
      gap <- replicate(TH$reps, {
        R <- rmvn(Tn, S)
        ldmvn_rows(R, S) - ld_indep(R, s2)
      })
      res <- rbind(res, data.frame(n_g = ng, covariance = COVS$label[j],
        predicted = pred, observed = mean(gap), mc_se = stats::sd(gap) / sqrt(length(gap)),
        sd_gap = stats::sd(gap), sd_predicted = sd_pred))
    }
  }
  res$z <- (res$observed - res$predicted) / res$mc_se
  print(res, row.names = FALSE, digits = 4)
  save_csv(res, "lemma.csv")
  verdict("lemma", "expected gap equals -(T/2) log|R|", all(abs(res$z) < 3.5),
          sprintf("max |z| = %.2f over %d settings, %d replications each",
                  max(abs(res$z)), nrow(res), TH$reps))
  verdict("lemma", "the gap is nonnegative and zero only without correlation",
          all(res$predicted >= 0),
          sprintf("predicted gap from %.1f to %.1f", min(res$predicted), max(res$predicted)))
  verdict("lemma", "the Monte Carlo spread matches its formula",
          all(abs(res$sd_gap / res$sd_predicted - 1) < 0.15),
          sprintf("sd ratio from %.2f to %.2f", min(res$sd_gap / res$sd_predicted),
                  max(res$sd_gap / res$sd_predicted)))
  cat("\n")
}


## ===========================================================================
## THE SAME GAP UNDER ESTIMATION: parameters estimated, latent path smoothed
## ===========================================================================
if ("estimated" %in% TH$parts) {
  cat("ESTIMATED: the gap with estimated parameters and smoothed states\n")
  scen <- dgp_scenarios()
  res <- NULL
  for (rp in seq_len(min(TH$reps, 30L))) {
    dat <- dgp_draw(15L, 120L, 1L, 0, scen[scen$id == "S0", ], rep = 900L + rp)
    mod <- sim_model(dat, 15L)
    f <- try(Stem::STEM_Estimation(mod, distance = "geo", precision = 0.01), silent = TRUE)
    if (inherits(f, "try-error")) next
    ph <- f$estimates$phi.hat
    co <- dat$coordinates
    cv <- cov_of(co, ph$sigma2eps, ph$sigma2omega, ph$theta)
    xb <- matrix(as.numeric(dat$covariates %*% as.numeric(ph$beta)), nrow = 120L)
    R  <- dat$z - xb - as.numeric(f$estimates$y.smoothed)
    Rm <- stats::cov2cor(cv$S)
    res <- rbind(res, data.frame(rep = rp,
      gap = ldmvn_rows(R, cv$S) - ld_indep(R, cv$S[1, 1]),
      predicted = -60 * as.numeric(determinant(Rm, logarithm = TRUE)$modulus)))
  }
  res$ratio <- res$gap / res$predicted
  save_csv(res, "lemma-estimated.csv")
  ## fits whose estimated range is zero predict no gap at all; their ratio is a
  ## division by zero and says nothing
  ok <- res$predicted > 1e-6
  print(summary(res$ratio[ok]))
  ## The Lemma holds GIVEN the latent path. With the smoothed state plugged in,
  ## the common part of the spatial error is absorbed by the latent process, the
  ## residuals keep little of the correlation the fitted covariance implies, and
  ## the realized gap is far below -(T/2) log|R-hat|. This is what the paper
  ## states after the Lemma, and what is checked here.
  verdict("estimated", "with the smoothed state plugged in, the realized gap is far below the Lemma's",
          stats::median(res$ratio[ok]) < 0.5,
          sprintf("median observed/predicted = %.2f over %d fits (%d with no fitted correlation left out)",
                  stats::median(res$ratio[ok]), sum(ok), sum(!ok)))
  cat("\n")
}


## ===========================================================================
## THE CONDITIONAL SCORE: E[lc_i(g) - l_i(g)] = +Delta_i(g) if i in g,
##   T { -1/2 log(1-q) - q (1+a) / (2 (1-q)) } < 0 if i is independent of g
## ===========================================================================
if ("cond" %in% TH$parts) {
  cat("COND: the conditional score against the marginal one\n")
  set.seed(202)
  Tn <- 100L; ng <- 15L
  loc <- dgp_locations(60L, 3L, 1, seed = 11)
  res <- NULL
  for (j in seq_len(nrow(COVS))) {
    s2e <- tot * COVS$share[j]; s2o <- tot * (1 - COVS$share[j])
    th  <- 1 / (COVS$range[j] * 1000)
    idx_g <- which(loc$labels == 1L)[seq_len(ng + 1L)]
    ## the target inside the regime: its most central member
    co_g <- loc$coords[idx_g, , drop = FALSE]
    cv   <- cov_of(co_g, s2e, s2o, th)
    tgt  <- which.min(rowSums(cv$dm))
    J    <- setdiff(seq_len(ng + 1L), tgt)
    Rm   <- stats::cov2cor(cv$S)
    q_in <- as.numeric(t(Rm[J, tgt]) %*% solve(Rm[J, J], Rm[J, tgt]))
    ## the target outside: the location of regime 2 nearest to the regime
    out_i <- which(loc$labels == 2L)
    dmix <- as.matrix(geodist::geodist(loc$coords[out_i, , drop = FALSE], co_g[J, , drop = FALSE],
                                       measure = "geodesic"))
    o <- out_i[which.min(apply(dmix, 1, min))]
    co_all <- rbind(co_g[J, , drop = FALSE], loc$coords[o, , drop = FALSE])
    cv2 <- cov_of(co_all, s2e, s2o, th)
    rho_o <- stats::cov2cor(cv2$S)[seq_len(ng), ng + 1L]
    q_out <- as.numeric(t(rho_o) %*% solve(Rm[J, J], rho_o))

    for (case in c("member", "outsider", "outsider, variance x2")) {
      a <- if (case == "outsider, variance x2") 2 else 1
      d_obs <- replicate(TH$reps, {
        if (case == "member") {
          E <- rmvn(Tn, cv$S)
          Z <- E[, c(J, tgt)]
          dmm <- cv$dm[c(J, tgt), c(J, tgt)]
        } else {
          Z <- cbind(rmvn(Tn, cv$S[J, J]), stats::rnorm(Tn, sd = sqrt(a * (s2e + s2o))))
          dmm <- cv2$dm
        }
        cs <- cond_scores(z = Z, covariates = matrix(0, (ng + 1L) * Tn, 1), Tobs = Tn,
                          beta = 0, ysm = matrix(0, Tn, 1), Kmat = matrix(0, ng + 1L, 1),
                          sigma2eps = s2e, sigma2omega = s2o, theta = th,
                          members = seq_len(ng), dm = dmm)
        marg <- ld_indep(Z[, ng + 1L], s2e + s2o)
        cs$cond[ng + 1L] - marg
      })
      q <- if (case == "member") q_in else q_out
      pred <- if (case == "member") -0.5 * Tn * log(1 - q) else
        Tn * (-0.5 * log(1 - q) - q * (1 + a) / (2 * (1 - q)))
      res <- rbind(res, data.frame(covariance = COVS$label[j], case = case, q = q,
        predicted = pred, observed = mean(d_obs), mc_se = stats::sd(d_obs) / sqrt(length(d_obs))))
    }
  }
  res$z <- (res$observed - res$predicted) / res$mc_se
  print(res, row.names = FALSE, digits = 4)
  save_csv(res, "conditional-score.csv")
  verdict("cond", "member: the conditional score gains Delta_i(g) in expectation",
          all(abs(res$z[res$case == "member"]) < 3.5),
          sprintf("max |z| = %.2f", max(abs(res$z[res$case == "member"]))))
  verdict("cond", "outsider: it loses the predicted amount, with the predicted sign",
          all(abs(res$z[res$case != "member"]) < 3.5) &&
            all(res$observed[res$case != "member" & res$q > 0.01] < 0),
          sprintf("max |z| = %.2f", max(abs(res$z[res$case != "member"]))))
  cat("\n")
}


## ===========================================================================
## THE DECISION, with the parameters and the latent path known
##
## Each location is scored against each regime, whose members are the true
## ones (the location itself excluded), and assigned to the best. No penalty,
## no iteration: this isolates what each score knows.
## ===========================================================================
if ("decision" %in% TH$parts) {
  cat("DECISION: which score assigns best, parameters known\n")
  set.seed(303)
  n <- 60L; Tn <- 60L; K <- 3L
  res_sd <- sqrt(tot)
  cases <- list(
    shared = list(beta1 = base$beta[2] + c(0, 0.125, 0.25) * res_sd,
                  share = rep(base$sigma2eps / tot, 3), range = rep(123, 3), field = "regime"),
    covariance = list(beta1 = rep(base$beta[2], 3), share = c(base$sigma2eps / tot, 0.57, 0.30),
                      range = c(123, 69, 15), field = "regime"),
    global = list(beta1 = base$beta[2] + c(0, 0.125, 0.25) * res_sd,
                  share = rep(base$sigma2eps / tot, 3), range = rep(123, 3), field = "global"))
  res <- NULL
  for (cn in names(cases)) for (w in SIM_LEVELS$omega) {
    cs0 <- cases[[cn]]
    err <- matrix(0, 0, 4)
    for (rp in seq_len(max(20L, TH$reps %/% 8L))) {
      loc <- dgp_locations(n, K, w, seed = 5000L + rp)
      g <- loc$labels
      dmall <- as.matrix(geodist::geodist(loc$coords, measure = "geodesic"))
      x <- matrix(stats::rnorm(n * Tn), Tn, n)
      E <- matrix(0, Tn, n)
      if (cs0$field == "global") {
        S <- tot * (1 - cs0$share[1]) * exp(-dmall / (cs0$range[1] * 1000))
        diag(S) <- tot
        E <- rmvn(Tn, S)
      } else for (h in seq_len(K)) {
        ii <- which(g == h)
        S <- tot * (1 - cs0$share[h]) * exp(-dmall[ii, ii] / (cs0$range[h] * 1000))
        diag(S) <- tot
        E[, ii] <- rmvn(Tn, S)
      }
      Z <- x * matrix(cs0$beta1[g], Tn, n, byrow = TRUE) + E
      ## scores of every location against every regime, three ways
      M <- C <- A <- matrix(NA_real_, n, K)
      for (h in seq_len(K)) {
        s2e <- tot * cs0$share[h]; s2o <- tot * (1 - cs0$share[h])
        th  <- 1 / (cs0$range[h] * 1000)
        R <- Z - x * cs0$beta1[h]
        M[, h] <- apply(R, 2, function(r) ld_indep(r, tot))
        sc <- cond_scores(z = R, covariates = matrix(0, n * Tn, 1), Tobs = Tn, beta = 0,
                          ysm = matrix(0, Tn, 1), Kmat = matrix(0, n, 1),
                          sigma2eps = s2e, sigma2omega = s2o, theta = th,
                          members = which(g == h), dm = dmall)
        C[, h] <- sc$cond
        A[, h] <- M[, h] + sc$delta
      }
      ## a location is on a boundary if its nearest neighbour is in another regime
      diag(dmall) <- Inf
      bnd <- g[apply(dmall, 1, which.min)] != g
      err <- rbind(err, cbind(max.col(M) != g, max.col(A) != g, max.col(C) != g, bnd))
    }
    for (b in c(FALSE, TRUE)) {
      sel <- err[, 4] == b
      res <- rbind(res, data.frame(case = cn, omega = w, boundary = b, n_loc = sum(sel),
        marginal = mean(err[sel, 1]), corrected = mean(err[sel, 2]),
        conditional = mean(err[sel, 3])))
    }
  }
  print(res, row.names = FALSE, digits = 3)
  save_csv(res, "decision.csv")
  modl <- res[res$case != "global", ]
  verdict("decision", "under the model, the conditional score never assigns worse",
          all(modl$conditional <= modl$marginal + 0.02),
          sprintf("largest excess of conditional over marginal: %.3f",
                  max(modl$conditional - modl$marginal)))
  cov_rows <- res[res$case == "covariance", ]
  verdict("decision", "when regimes differ only in covariance the marginal score is blind",
          all(cov_rows$marginal > 0.5),
          sprintf("marginal error %.2f-%.2f, conditional %.2f-%.2f (chance = 0.67)",
                  min(cov_rows$marginal), max(cov_rows$marginal),
                  min(cov_rows$conditional), max(cov_rows$conditional)))
  ## the corrected score adds a reward for proximity to the members of a
  ## regime: it cannot help where proximity says nothing about membership
  ## pooled over boundary and interior locations, weighted by their numbers
  wmean <- function(x, w) sum(x * w) / sum(w)
  by_w <- do.call(rbind, lapply(split(res[res$case == "shared", ], res$omega[res$case == "shared"]),
    function(s) data.frame(omega = s$omega[1], marginal = wmean(s$marginal, s$n_loc),
                           corrected = wmean(s$corrected, s$n_loc))))
  bd <- res[res$case == "shared" & res$boundary, ]
  verdict("decision", "the corrected score is a proximity prior: no help without separation, harm at boundaries",
          by_w$corrected[by_w$omega == 0] >= by_w$marginal[by_w$omega == 0] - 0.01 &&
            all(bd$corrected > bd$marginal),
          sprintf("overall error, corrected minus marginal, by omega: %s; at boundaries: %s",
                  paste(sprintf("%.2f: %+.3f", by_w$omega, by_w$corrected - by_w$marginal), collapse = ", "),
                  paste(sprintf("%+.3f", bd$corrected - bd$marginal), collapse = ", ")))
  gl <- res[res$case == "global", ]
  verdict("decision", "with one field across regimes the conditional score loses where proximity misleads",
          all(gl$conditional[gl$boundary] > gl$marginal[gl$boundary]),
          sprintf("boundary error, conditional minus marginal: %s",
                  paste(sprintf("%+.3f", gl$conditional[gl$boundary] - gl$marginal[gl$boundary]),
                        collapse = ", ")))
  cat("\n")
}


## ===========================================================================
## THE FULL PROCEDURE under the three scores
## ===========================================================================
if ("algorithm" %in% TH$parts) {
  cat("ALGORITHM: the full procedure under the three scores\n")
  scen <- dgp_scenarios()
  cells <- rbind(
    data.frame(id = "S4", omega = SIM_LEVELS$omega, K = 3L),
    data.frame(id = c("S3b", "S1a", "S0"), omega = 0.70, K = 3L))
  n <- 60L; Tn <- 120L
  res <- NULL; traces <- NULL
  for (ci in seq_len(nrow(cells))) for (rp in seq_len(TH$alg_reps)) {
    cl <- cells[ci, ]
    dat <- dgp_draw(n, Tn, cl$K, cl$omega, scen[scen$id == cl$id, ], rep = 7000L + rp)
    mod <- sim_model(dat, n)
    W <- neighbors(dat$coordinates, knn = 5, distance = "geo")$W
    ij <- which(upper.tri(W) & W > 0, arr.ind = TRUE)
    for (sc in c("marginal", "conditional", "corrected")) for (ph in c(0, 0.5)) {
      t0 <- proc.time()[["elapsed"]]
      f <- try(Stem::SCSTEM_Estimation(mod, k = cl$K, phi_penalty = ph, distance = "geo",
                                       knn = 5, precision = 0.05, min_cluster_size = N_MIN,
                                       score = sc, verbose = FALSE), silent = TRUE)
      secs <- proc.time()[["elapsed"]] - t0
      if (inherits(f, "try-error")) next
      conc <- mean(f$group[ij[, 1]] == f$group[ij[, 2]])
      ## the same share under a random partition with the same sizes
      conc0 <- mean(replicate(50, { p <- sample(f$group); mean(p[ij[, 1]] == p[ij[, 2]]) }))
      res <- rbind(res, data.frame(id = cl$id, omega = cl$omega, rep = rp, score = sc,
        phi = ph, ari = ari(f$group, dat$labels), concordance = conc,
        concordance_random = conc0, secs = secs))
      tr <- f$obj_trace
      if (nrow(tr) > 1)
        traces <- rbind(traces, data.frame(id = cl$id, omega = cl$omega, rep = rp,
          score = sc, phi = ph, iter = tr$iter, objective = tr$objective,
          objective_before = tr$objective_before))
    }
  }
  save_csv(res, "algorithm.csv"); save_csv(traces, "algorithm-traces.csv")
  agg <- stats::aggregate(cbind(ari, concordance, concordance_random, secs) ~ id + omega + score + phi,
                          data = res, FUN = mean)
  print(agg[order(agg$id, agg$omega, agg$phi, agg$score), ], row.names = FALSE, digits = 3)

  ## label step: never a decrease; parameter step: decreases do occur
  lab_ok <- all(traces$objective >= traces$objective_before - 1e-6)
  sp <- split(traces, list(traces$id, traces$omega, traces$rep, traces$score, traces$phi), drop = TRUE)
  par_steps <- unlist(lapply(sp, function(s) {
    s <- s[order(s$iter), ]
    if (nrow(s) < 2) return(NULL)
    s$objective_before[-1] - s$objective[-nrow(s)]
  }))
  verdict("algorithm", "the label step never lowers the objective", lab_ok,
          sprintf("%d label steps checked", nrow(traces)))
  verdict("algorithm", "the parameter step can lower it: the alternation is not monotone",
          any(par_steps < -1e-6),
          sprintf("%d of %d parameter steps lowered the objective",
                  sum(par_steps < -1e-6), length(par_steps)))
  s0 <- agg[agg$id == "S0", ]
  verdict("algorithm", "S0: no partition is recovered, whatever the score",
          all(abs(s0$ari) < 0.15),
          sprintf("mean ARI from %.3f to %.3f", min(s0$ari), max(s0$ari)))
  verdict("algorithm", "S0: the penalty makes the returned partition compact",
          all(s0$concordance[s0$phi > 0] > s0$concordance_random[s0$phi > 0]),
          sprintf("concordant edges %.2f-%.2f against %.2f at random",
                  min(s0$concordance[s0$phi > 0]), max(s0$concordance[s0$phi > 0]),
                  mean(s0$concordance_random)))
  cat("\n")
}


## ===========================================================================
## THE LABEL STEP at fixed scores: ascent, finite termination, thresholds
##
## Proposition: at a fixed point of the label step with effective penalty
## phi_eff, (i) a location whose best-versus-second-best score gap exceeds
## phi_eff times its degree carries its data-best label, and (ii) a location
## whose score spread is below phi_eff carries a label of maximal count among its
## neighbours. Checked on score matrices produced by the package.
## ===========================================================================
if ("label" %in% TH$parts) {
  cat("LABEL: the label step at fixed scores\n")
  scen <- dgp_scenarios()
  icm <- function(LL, nb, phi_eff, lab) {           # the rule of the package, unconstrained
    Q <- function(l) sum(LL[cbind(seq_along(l), l)]) +
      phi_eff * sum(vapply(seq_along(l), function(i) sum(l[nb[[i]]] == l[i]), 0)) / 2
    q <- Q(lab); sweeps <- 0L; mono <- TRUE
    repeat {
      sweeps <- sweeps + 1L; changed <- FALSE
      for (i in seq_along(lab)) {
        v <- LL[i, ] + phi_eff * tabulate(lab[nb[[i]]], ncol(LL))
        b <- which.max(v)
        if (v[b] > v[lab[i]] + 1e-12) { lab[i] <- b; changed <- TRUE }
      }
      q2 <- Q(lab); if (q2 < q - 1e-9) mono <- FALSE; q <- q2
      if (!changed || sweeps > 1000L) break
    }
    list(lab = lab, sweeps = sweeps, mono = mono)
  }
  res <- NULL
  for (rp in 1:4) {
    dat <- dgp_draw(100L, 120L, 3L, 0.70, scen[scen$id == "S4", ], rep = 8000L + rp)
    f <- Stem::SCSTEM_Estimation(sim_model(dat, 100L), k = 3, phi_penalty = 0.5,
                                 distance = "geo", knn = 5, precision = 0.05,
                                 min_cluster_size = N_MIN, verbose = FALSE)
    LL <- f$score_last
    nb <- neighbors(dat$coordinates, knn = 5, distance = "geo")$nb
    deg <- vapply(nb, length, 0L)
    srt <- t(apply(LL, 1, sort, decreasing = TRUE))
    gap <- srt[, 1] - srt[, 2]; spread <- srt[, 1] - srt[, ncol(srt)]
    cscale <- stats::median(spread) / mean(deg)
    best <- max.col(LL)
    for (ph in c(0, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16)) {
      pe <- ph * cscale
      r <- icm(LL, nb, pe, best)
      cnt <- t(vapply(seq_along(nb), function(i) tabulate(r$lab[nb[[i]]], ncol(LL)), numeric(ncol(LL))))
      is_major <- cnt[cbind(seq_along(nb), r$lab)] == apply(cnt, 1, max)
      data_dom <- gap > pe * deg
      nb_dom   <- spread < pe
      res <- rbind(res, data.frame(rep = rp, phi = ph, sweeps = r$sweeps, monotone = r$mono,
        follows_data = mean(r$lab == best), bound_data = mean(data_dom),
        viol_data = sum(data_dom & r$lab != best),
        follows_majority = mean(is_major), bound_majority = mean(nb_dom),
        viol_majority = sum(nb_dom & !is_major)))
    }
  }
  agg <- stats::aggregate(cbind(sweeps, follows_data, bound_data, follows_majority, bound_majority)
                          ~ phi, data = res, FUN = mean)
  print(agg, row.names = FALSE, digits = 3)
  save_csv(res, "label-step.csv")
  verdict("label", "every sweep raises the objective and the step terminates",
          all(res$monotone) && all(res$sweeps < 1000L),
          sprintf("at most %d sweeps to a fixed point", max(res$sweeps)))
  verdict("label", "(i) a large enough data gap always wins", sum(res$viol_data) == 0,
          sprintf("%d violations", sum(res$viol_data)))
  verdict("label", "(ii) a small enough spread always follows the neighbours",
          sum(res$viol_majority) == 0, sprintf("%d violations", sum(res$viol_majority)))
  cat("\n")
}


## ===========================================================================
## THE IN-SAMPLE CRITERIA UNDER THE NULL
##
## Canonical case: d independent locations, T periods, known variance. The
## best two-group split of the location means gains, in expectation, d/pi in
## log-likelihood, while BIC charges (1/2) log(dT) for the extra mean: the
## criterion must eventually split a homogeneous network.
## ===========================================================================
if ("ic" %in% TH$parts) {
  cat("IC: the in-sample criteria under the null\n")
  set.seed(404)
  best_split_gain <- function(m, Tn) {          # exact optimum of the 2-means split
    s <- sort(m); d <- length(s); cs <- cumsum(s); tot_s <- cs[d]
    k <- seq_len(d - 1L)
    bss <- cs[k]^2 / k + (tot_s - cs[k])^2 / (d - k) - tot_s^2 / d
    0.5 * Tn * max(bss)                          # sigma^2 = 1
  }
  res <- NULL
  for (Tn in c(60L, 365L)) for (d in c(10L, 20L, 50L, 100L, 200L, 500L, 1000L)) {
    G <- replicate(TH$reps, best_split_gain(stats::rnorm(d, sd = 1 / sqrt(Tn)), Tn))
    res <- rbind(res, data.frame(T = Tn, d = d, gain_per_location = mean(G) / d,
      predicted = 1 / pi, se = stats::sd(G) / d / sqrt(length(G)),
      p_split = mean(2 * G > log(d * Tn)),
      p_split_threshold = as.numeric(d / pi > 0.5 * log(d * Tn))))
  }
  print(res, row.names = FALSE, digits = 3)
  save_csv(res, "ic-canonical.csv")
  big <- res[res$d >= 200, ]
  verdict("ic", "the null gain grows linearly in d, at rate 1/pi",
          all(abs(big$gain_per_location - 1 / pi) < 0.03),
          sprintf("gain per location at d >= 200: %.3f-%.3f against 1/pi = %.3f",
                  min(big$gain_per_location), max(big$gain_per_location), 1 / pi))
  verdict("ic", "BIC splits a homogeneous network with probability tending to one",
          all(res$p_split[res$d >= 200] > 0.95),
          sprintf("P(split) at d = 1000: %.3f (T = 60), %.3f (T = 365)",
                  res$p_split[res$d == 1000 & res$T == 60], res$p_split[res$d == 1000 & res$T == 365]))
  cat("\n")
}

## ===========================================================================
## THE SAME MECHANISM IN SC-STEM: the null gain of k = 2 over k = 1, as the
## network grows, on data with no regimes at all
## ===========================================================================
if ("ic-stem" %in% TH$parts) {
  cat("IC-STEM: the null gain of a second regime in SC-STEM\n")
  scen <- dgp_scenarios()
  res <- NULL
  for (n in c(20L, 50L, 100L, 200L)) for (rp in seq_len(5L)) {
    dat <- dgp_draw(n, 60L, 1L, 0, scen[scen$id == "S0", ], rep = 9000L + rp)
    ic <- try(Stem::SCSTEM_Infocrit(sim_model(dat, n), k_grid = 1:2, phi_grid = 0,
                                    distance = "geo", knn = 5, verbose = FALSE,
                                    min_cluster_size = N_MIN), silent = TRUE)
    if (inherits(ic, "try-error")) next
    tb <- ic$table
    l1 <- tb$loglik[tb$k == 1]; l2 <- tb$loglik[tb$k == 2]
    res <- rbind(res, data.frame(n = n, rep = rp, gain = l2 - l1,
      gain_per_location = (l2 - l1) / n,
      bic_prefers_2 = tb$BIC[tb$k == 2] < tb$BIC[tb$k == 1]))
  }
  agg <- stats::aggregate(cbind(gain, gain_per_location, bic_prefers_2) ~ n, data = res, FUN = mean)
  print(agg, row.names = FALSE, digits = 3)
  save_csv(res, "ic-stem.csv")
  fit <- stats::lm(gain ~ n, data = res)
  verdict("ic-stem", "in SC-STEM too the null gain grows linearly in the network size",
          stats::coef(fit)[["n"]] > 0 && summary(fit)$r.squared > 0.5,
          sprintf("slope %.2f per location, R^2 = %.2f; BIC prefers k = 2 in %.0f%% of fits at n = %d",
                  stats::coef(fit)[["n"]], summary(fit)$r.squared,
                  100 * agg$bic_prefers_2[nrow(agg)], agg$n[nrow(agg)]))
  cat("\n")
}

## ===========================================================================
## TABLES for the paper, from the CSVs of the parts above (--parts=tables)
## ===========================================================================
if ("tables" %in% TH$parts) {
  cat("TABLES\n")
  rd <- function(f) { p <- file.path(OUT, f); if (file.exists(p)) utils::read.csv(p) else NULL }
  f2 <- function(x, d = 2) formatC(x, format = "f", digits = d)
  tex <- function(lines, name) { writeLines(lines, file.path(OUT, name)); message("wrote ", name) }

  dec <- rd("decision.csv")
  if (!is.null(dec)) {
    agg <- do.call(rbind, lapply(split(dec, list(dec$case, dec$omega), drop = TRUE), function(s)
      data.frame(case = s$case[1], omega = s$omega[1],
                 marginal = sum(s$marginal * s$n_loc) / sum(s$n_loc),
                 corrected = sum(s$corrected * s$n_loc) / sum(s$n_loc),
                 conditional = sum(s$conditional * s$n_loc) / sum(s$n_loc),
                 m_b = s$marginal[s$boundary], c_b = s$corrected[s$boundary],
                 k_b = s$conditional[s$boundary])))
    agg <- agg[order(match(agg$case, c("shared", "covariance", "global")), agg$omega), ]
    tex(c("\\begin{tabular}{llrrrrrr}", "\\toprule",
          "& & \\multicolumn{3}{c}{all locations} & \\multicolumn{3}{c}{boundary locations} \\\\",
          "\\cmidrule(lr){3-5}\\cmidrule(l){6-8}",
          "regimes & $\\omega$ & marginal & $+\\Delta_i$ & conditional & marginal & $+\\Delta_i$ & conditional \\\\",
          "\\midrule",
          sprintf("%s & %s & %s & %s & %s & %s & %s & %s \\\\",
                  ifelse(duplicated(agg$case), "", paste0("\\emph{", agg$case, "}")),
                  f2(agg$omega), f2(agg$marginal), f2(agg$corrected), f2(agg$conditional),
                  f2(agg$m_b), f2(agg$c_b), f2(agg$k_b)),
          "\\bottomrule", "\\end{tabular}"), "tab_theory_decision.tex")
  }

  alg <- rd("algorithm.csv")
  if (!is.null(alg)) {
    a <- stats::aggregate(cbind(ari, concordance, concordance_random) ~ id + omega + phi + score,
                          data = alg, FUN = mean)
    w <- stats::reshape(a[, c("id", "omega", "phi", "score", "ari")], direction = "wide",
                        idvar = c("id", "omega", "phi"), timevar = "score")
    cc <- stats::aggregate(concordance ~ id + omega + phi, data = alg[alg$score == "marginal", ], FUN = mean)
    cr <- stats::aggregate(concordance_random ~ id + omega + phi, data = alg, FUN = mean)
    w <- merge(merge(w, cc), cr)
    ord <- c("S4", "S1a", "S3b", "S0")
    w <- w[order(match(w$id, ord), w$omega, w$phi), ]
    tex(c("\\begin{tabular}{lrrrrrrr}", "\\toprule",
          "& & & \\multicolumn{3}{c}{ARI at the true $K$} & \\multicolumn{2}{c}{concordant edges} \\\\",
          "\\cmidrule(lr){4-6}\\cmidrule(l){7-8}",
          "scenario & $\\omega$ & $\\phi$ & marginal & $+\\Delta_i$ & conditional & marginal & random \\\\",
          "\\midrule",
          sprintf("%s & %s & %s & %s & %s & %s & %s & %s \\\\", w$id, f2(w$omega), f2(w$phi, 1),
                  f2(w$ari.marginal), f2(w$ari.corrected), f2(w$ari.conditional),
                  f2(w$concordance), f2(w$concordance_random)),
          "\\bottomrule", "\\end{tabular}"), "tab_theory_algorithm.tex")
  }

  lem <- rd("lemma.csv")
  if (!is.null(lem)) tex(c("\\begin{tabular}{rlrrrrr}", "\\toprule",
    "$d_h$ & covariance & $-\\frac{T}{2}\\log|R_h|$ & observed & s.e. & sd, formula & sd, observed \\\\",
    "\\midrule",
    sprintf("%d & %s & %s & %s & %s & %s & %s \\\\", lem$n_g, lem$covariance, f2(lem$predicted, 1),
            f2(lem$observed, 1), f2(lem$mc_se), f2(lem$sd_predicted, 1), f2(lem$sd_gap, 1)),
    "\\bottomrule", "\\end{tabular}"), "tab_theory_lemma.tex")

  cs <- rd("conditional-score.csv")
  if (!is.null(cs)) tex(c("\\begin{tabular}{llrrrr}", "\\toprule",
    "covariance & location & $q_{i,h}$ & predicted & observed & s.e. \\\\", "\\midrule",
    sprintf("%s & %s & %s & %s & %s & %s \\\\", cs$covariance, cs$case, f2(cs$q, 3),
            f2(cs$predicted, 1), f2(cs$observed, 1), f2(cs$mc_se)),
    "\\bottomrule", "\\end{tabular}"), "tab_theory_cond.tex")

  lb <- rd("label-step.csv")
  if (!is.null(lb)) {
    a <- stats::aggregate(cbind(sweeps, follows_data, bound_data, follows_majority, bound_majority)
                          ~ phi, data = lb, FUN = mean)
    v <- stats::aggregate(cbind(viol_data, viol_majority) ~ phi, data = lb, FUN = sum)
    a$viol <- v$viol_data + v$viol_majority
    tex(c("\\begin{tabular}{rrrrrrr}", "\\toprule",
      "$\\phi$ & sweeps & follow the data & bound (i) & follow the neighbours & bound (ii) & violations \\\\",
      "\\midrule",
      sprintf("%s & %s & %s & %s & %s & %s & %d \\\\", f2(a$phi), f2(a$sweeps, 1), f2(a$follows_data),
              f2(a$bound_data), f2(a$follows_majority), f2(a$bound_majority),
              as.integer(a$viol)),
      "\\bottomrule", "\\end{tabular}"), "tab_theory_label.tex")
  }

  icc <- rd("ic-canonical.csv")
  if (!is.null(icc)) tex(c("\\begin{tabular}{rrrrr}", "\\toprule",
    "$T$ & $d$ & $G_d/d$ & s.e. & $\\Pr(\\text{BIC splits})$ \\\\", "\\midrule",
    sprintf("%d & %d & %s & %s & %s \\\\", icc$T, icc$d, f2(icc$gain_per_location, 3), f2(icc$se, 3),
            f2(icc$p_split)),
    "\\bottomrule", "\\end{tabular}"), "tab_theory_ic.tex")

  ics <- rd("ic-stem.csv")
  if (!is.null(ics)) {
    a <- stats::aggregate(cbind(gain, gain_per_location, bic_prefers_2) ~ n, data = ics, FUN = mean)
    tex(c("\\begin{tabular}{rrrr}", "\\toprule",
      "$d$ & gain of $k=2$ over $k=1$ & per location & $\\Pr(\\text{BIC prefers } k = 2)$ \\\\", "\\midrule",
      sprintf("%d & %s & %s & %s \\\\", a$n, f2(a$gain, 1), f2(a$gain_per_location, 3),
              f2(a$bic_prefers_2)),
      "\\bottomrule", "\\end{tabular}"), "tab_theory_icstem.tex")
  }
  cat("\n")
}

## ---------------------------------------------------------------------------
if (!length(verdicts)) quit(save = "no")
V <- do.call(rbind, verdicts)
save_csv(V, "verdicts.csv")
cat("VERDICTS\n")
print(V[, c("part", "verdict", "claim")], row.names = FALSE, right = FALSE)
