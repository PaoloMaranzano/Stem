test_that("STEM_control validates its settings and merges partial lists", {
  expect_s3_class(STEM_control(), "STEM_control")
  expect_error(STEM_control(em_tol_par = -1), "positive")
  expect_error(STEM_control(alt_maxit = -1), "non-negative")
  expect_identical(STEM_control(alt_maxit = 0)$alt_maxit, 0L)
  ctl <- Stem:::stem_control_resolve(list(em_maxit = 7))
  expect_identical(ctl$em_maxit, 7L)
  expect_identical(ctl$em_tol_par, STEM_control()$em_tol_par)
  expect_error(Stem:::stem_control_resolve(list(nonsense = 1)), "unknown setting")
  old <- options(Stem.control = list(em_maxit = 11)); on.exit(options(old))
  expect_identical(Stem:::stem_control_resolve(NULL)$em_maxit, 11L)
})

test_that("the EM stopping rule ignores the fixed loadings and reports both criteria", {
  mod <- po_model(Tn = 60L)
  fit <- STEM_Estimation(mod, distance = "geo",
                         control = list(em_tol_par = 1e-3, em_tol_loglik = 1e-2, em_maxit = 200))
  cp <- fit$estimates$convergence.par
  expect_true(all(c("conv.par", "conv.log", "converged", "max.rel.par", "delta.loglik") %in% names(cp)))
  expect_true(cp$iterEM <= 200)
  ## each criterion reported as met is below its tolerance
  if (isTRUE(cp$conv.par)) expect_lt(cp$max.rel.par, 1e-3)
  if (isTRUE(cp$conv.log)) expect_lt(cp$delta.loglik, 1e-2)
  ## the historical arguments still work, as overrides of control
  fit2 <- STEM_Estimation(mod, distance = "geo", precision = 0.5, max.iter = 3)
  expect_lte(fit2$estimates$convergence.par$iterEM, 3)
})

test_that("em_stop combines the two criteria as asked", {
  expect_identical(STEM_control()$em_stop, "any")
  expect_error(STEM_control(em_stop = "both"))
  mod <- po_model(Tn = 60L)
  ctl <- list(em_tol_par = 1e-3, em_tol_loglik = 1e-2, em_maxit = 200)
  any <- STEM_Estimation(mod, distance = "geo", control = c(ctl, em_stop = "any"))$estimates$convergence.par
  all <- STEM_Estimation(mod, distance = "geo", control = c(ctl, em_stop = "all"))$estimates$convergence.par
  if (isTRUE(any$converged)) expect_true(any$conv.par || any$conv.log)
  if (isTRUE(all$converged)) expect_true(all$conv.par && all$conv.log)
  expect_lte(any$iterEM, all$iterEM)
})

test_that("SC-STEM records its settings and starts the regimes from their own data", {
  mod <- po_model(Tn = 60L)
  fit <- SCSTEM_Estimation(mod, k = 2, phi_penalty = 0, distance = "geo", seed = 1)
  expect_s3_class(fit$input_args$control, "STEM_control")
  expect_length(fit$em_converged, 2L)
  ## a fixed partition refitted without moving it
  refit <- SCSTEM_Estimation(mod, k = 2, phi_penalty = 0, distance = "geo",
                             init_partition = fit$group, max_iter = 0)
  expect_identical(refit$group, fit$group)
  ## the starting values of a regime: least squares on its own locations
  s <- po_subset(Tn = 60L)
  idx <- which(fit$group == 1L)
  rows <- as.vector(outer(seq_len(60L), (idx - 1L) * 60L, "+"))
  ph <- Stem:::scstem_phi_ols(po_phi(), s$z[, idx], s$covariates[rows, ])
  ols <- stats::lm.fit(s$covariates[rows, ], as.vector(s$z[, idx]))$coefficients
  expect_equal(as.numeric(ph$beta), as.numeric(ols))
})
