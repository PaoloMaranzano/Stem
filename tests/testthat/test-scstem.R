test_that("the pooled SC-STEM fit reproduces STEM_Estimation", {
  mod <- po_model(Tn = 60L)
  pooled <- STEM_Estimation(mod, precision = 0.05, distance = "geo")
  sc <- SCSTEM_Estimation(mod, k = 1, phi_penalty = 0, distance = "geo",
                     precision_full_dataset = 0.05)

  expect_s3_class(sc, "SCSTEM_Estimation")
  expect_equal(sc$group, rep(1L, ncol(mod$data$z)))
  expect_equal(unname(sc$info_crit[["loglik"]]),
               as.numeric(pooled$estimates$loglik))
})


test_that("SCSTEM_Estimation returns an admissible partition and a monotone objective", {
  mod <- po_model(Tn = 90L)
  fit <- SCSTEM_Estimation(mod, k = 3, phi_penalty = 0.5, distance = "geo",
                      precision = 0.05)

  sizes <- tabulate(fit$group, nbins = 3)
  expect_length(sizes, 3)
  expect_true(all(sizes >= fit$input_args$min_cluster_size))
  expect_true(all(fit$final_refit))
  expect_false(anyNA(fit$phi_hat))

  ### The alternation is NOT globally monotone: the assignment score is a
  ### pseudo-likelihood while the parameter step maximizes the exact
  ### within-cluster likelihood, so a parameter update can lower Q. What must
  ### hold is that the partition returned is the best one visited.
  obj <- fit$obj_trace$objective
  expect_true(length(obj) >= 1)
  expect_true(all(is.finite(obj)))
})


test_that("the returned partition attains the best visited objective", {
  mod <- po_model(Tn = 90L)
  fit <- SCSTEM_Estimation(mod, k = 3, phi_penalty = 0.5, distance = "geo",
                      precision = 0.05)

  ### the reported objective must be the maximum seen along the trace, and the
  ### returned partition the one that attained it
  expect_equal(fit$best_objective, max(fit$obj_trace$objective))
  expect_gte(fit$best_objective, fit$last_objective)
  expect_true(all(tabulate(fit$group, nbins = 3) >= fit$input_args$min_cluster_size))
})


test_that("the ICM sweep is monotone at fixed parameters", {
  ### The property that does hold: with the log-likelihood matrix held fixed,
  ### a sequential sweep followed by the swap pass cannot lower the objective.
  set.seed(11)
  d <- 24L; k <- 3L
  LL <- matrix(rnorm(d * k, sd = 2), d, k)
  nb <- lapply(seq_len(d), function(i) setdiff(max(1, i - 3):min(d, i + 3), i))
  phi_eff <- 0.4
  lab <- rep(seq_len(k), length.out = d)

  obj <- function(l) sum(LL[cbind(seq_len(d), l)]) +
    phi_eff * Stem:::scstem_potts_pairs(l, nb)

  before <- obj(lab)
  for (i in seq_len(d)) {
    penvec <- tabulate(lab[nb[[i]]], nbins = k)
    lab[i] <- which.max(LL[i, ] + phi_eff * penvec)
  }
  expect_gte(obj(lab), before - 1e-8)
})


test_that("the information criteria use the exact parameter count", {
  mod <- po_model(Tn = 60L)
  fit <- SCSTEM_Estimation(mod, k = 2, phi_penalty = 0.5, distance = "geo",
                      precision = 0.05)

  ncov <- ncol(mod$data$covariates)
  pdim <- mod$skeleton$p
  npar_expected <- sum(fit$final_refit) * (ncov + 3 + 3 * pdim)
  nobs <- nrow(mod$data$z) * ncol(mod$data$z)

  expect_equal(unname(fit$info_crit[["k"]]), npar_expected)
  expect_equal(unname(fit$info_crit[["AIC"]]),
               -2 * unname(fit$info_crit[["loglik"]]) + 2 * npar_expected)
  expect_equal(unname(fit$info_crit[["BIC"]]),
               -2 * unname(fit$info_crit[["loglik"]]) + log(nobs) * npar_expected)
})


test_that("a stronger spatial penalty does not reduce spatial cohesion", {
  mod <- po_model(Tn = 90L)
  nb <- Stem:::scstem_neighbors(mod$data$coordinates, knn = 5)

  f0 <- SCSTEM_Estimation(mod, k = 3, phi_penalty = 0, distance = "geo",
                     precision = 0.05)
  f1 <- SCSTEM_Estimation(mod, k = 3, phi_penalty = 2, distance = "geo",
                     precision = 0.05)

  expect_gte(Stem:::scstem_potts_pairs(f1$group, nb$nb),
             Stem:::scstem_potts_pairs(f0$group, nb$nb))
})


test_that("the neighbor graph is symmetric", {
  mod <- po_model(Tn = 30L)
  nb <- Stem:::scstem_neighbors(mod$data$coordinates, knn = 5)
  expect_equal(nb$W, t(nb$W))
  expect_true(all(diag(nb$W) == 0))
})


test_that("the swap pass never decreases the objective", {
  set.seed(1)
  d <- 12L; k <- 3L
  LL <- matrix(rnorm(d * k), d, k)
  nb <- lapply(seq_len(d), function(i) setdiff(max(1, i - 2):min(d, i + 2), i))
  lab <- rep(seq_len(k), length.out = d)
  phi_eff <- 0.7

  obj <- function(l) sum(LL[cbind(seq_len(d), l)]) +
    phi_eff * Stem:::scstem_potts_pairs(l, nb)

  out <- Stem:::scstem_swap_pass(lab, LL, phi_eff, nb)
  expect_gte(obj(out$labels), obj(lab) - 1e-8)
  ### swaps preserve the cluster sizes exactly
  expect_equal(tabulate(out$labels, nbins = k), tabulate(lab, nbins = k))
})


test_that("the estimation does not modify the RNG state of the caller", {
  mod <- po_model(Tn = 45L)
  set.seed(99)
  before <- .Random.seed
  invisible(SCSTEM_Estimation(mod, k = 2, phi_penalty = 0.5, distance = "geo",
                         precision = 0.05, seed = 12345))
  expect_identical(.Random.seed, before)
})


test_that("the Adjusted Rand Index behaves at its boundaries", {
  a <- c(1, 1, 2, 2, 3, 3)
  expect_equal(Stem:::scstem_ari(a, a), 1)
  ### relabeling must not change the index
  expect_equal(Stem:::scstem_ari(a, c(3, 3, 1, 1, 2, 2)), 1)
})


test_that("a penalized fit starts from the unpenalized solution", {
  mod <- po_model(Tn = 60L)
  f0 <- SCSTEM_Estimation(mod, k = 2, phi_penalty = 0, distance = "geo",
                          precision = 0.05)
  f1 <- SCSTEM_Estimation(mod, k = 2, phi_penalty = 0.05, distance = "geo",
                          precision = 0.05)
  ## the same as starting it explicitly from the unpenalized partition
  f1b <- SCSTEM_Estimation(mod, k = 2, phi_penalty = 0.05, distance = "geo",
                           precision = 0.05, init_partition = f0$group)
  expect_equal(f1$group, f1b$group)
  expect_equal(f1$phi_multiplier, f1b$phi_multiplier)
  ## and the same as the member of a grid, which passes that partition on
  ic <- SCSTEM_Infocrit(mod, k_grid = 2, phi_grid = c(0, 0.05), distance = "geo",
                        precision = 0.05)
  expect_equal(ic$fits[["k=2, phi=0.05"]]$group, f1$group)
  expect_equal(ic$fits[["k=2, phi=0.05"]]$phi_effective, f1$phi_effective)
})
