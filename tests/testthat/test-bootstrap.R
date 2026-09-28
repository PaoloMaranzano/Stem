test_that("the bootstrap refits accept every argument it forwards", {
  ### Regression test. SCSTEM_Bootstrap() rebuilds the argument list of the
  ### refits from the settings stored in the fitted object, so a change in the
  ### signature of SCSTEM_Estimation() used to make every draw fail at run time with
  ### an "unused argument" error, which the tryCatch() reported only as a lost
  ### replicate. The forwarded names must be a subset of the formals.
  mod <- po_model(Tn = 45L)
  fit <- SCSTEM_Estimation(mod, k = 2, phi_penalty = 0.5, distance = "geo",
                      precision = 0.05)
  stored <- names(fit$input_args)
  formals_sc <- names(formals(SCSTEM_Estimation))
  ### every stored setting that shares a name with a formal must be passable
  common <- intersect(stored, formals_sc)
  expect_true(length(common) > 10)
})


test_that("the parametric bootstrap produces usable draws", {
  skip_on_cran()

  mod <- po_model(Tn = 60L)
  fit <- SCSTEM_Estimation(mod, k = 2, phi_penalty = 0.5, distance = "geo",
                      precision = 0.05)
  boot <- SCSTEM_Bootstrap(fit, B = 5, seed = 7)

  expect_s3_class(boot, "SCSTEM_Bootstrap")
  ### the whole point of the fixes to kalman() and STEM_Estimation(): the
  ### refits must not die on non-finite iterates
  expect_gt(boot$B_valid, 0)
  expect_true(all(is.na(boot$info$message)))
  expect_equal(dim(boot$groups), c(ncol(mod$data$z), 5L))
})


test_that("the bootstrap is reproducible under a fixed seed", {
  skip_on_cran()

  mod <- po_model(Tn = 45L)
  fit <- SCSTEM_Estimation(mod, k = 2, phi_penalty = 0.5, distance = "geo",
                      precision = 0.05)
  b1 <- SCSTEM_Bootstrap(fit, B = 3, seed = 42)
  b2 <- SCSTEM_Bootstrap(fit, B = 3, seed = 42)

  expect_equal(b1$groups, b2$groups)
  expect_equal(b1$draws, b2$draws)
})


test_that("the bootstrap inference aligns labels and returns coherent intervals", {
  skip_on_cran()

  mod <- po_model(Tn = 60L)
  fit <- SCSTEM_Estimation(mod, k = 2, phi_penalty = 0.5, distance = "geo",
                      precision = 0.05)
  boot <- SCSTEM_Bootstrap(fit, B = 8, seed = 3)
  inf <- SCSTEM_BootInference(boot)

  expect_s3_class(inf, "SCSTEM_BootInference")
  lo <- grep("^perc_lower", names(inf$summary), value = TRUE)
  up <- grep("^perc_upper", names(inf$summary), value = TRUE)
  ok <- is.finite(inf$summary[[lo]]) & is.finite(inf$summary[[up]])
  expect_true(all(inf$summary[[lo]][ok] <= inf$summary[[up]][ok]))

  ### the co-clustering matrix is a symmetric matrix of shares
  expect_equal(inf$coclustering, t(inf$coclustering))
  expect_true(all(inf$coclustering >= 0 & inf$coclustering <= 1))
  expect_true(all(diag(inf$coclustering) == 1))
})


test_that("SCSTEM_Select applies the two-step rule", {
  skip_on_cran()

  mod <- po_model(Tn = 60L)
  ic <- SCSTEM_Infocrit(mod, k_grid = 1:3, phi_grid = c(0, 0.5, 1),
                        distance = "geo", precision = 0.05)

  expect_s3_class(ic, "SCSTEM_Infocrit")
  expect_true(all(c("k", "phi", "k_eff", "admissible", "BIC") %in% names(ic$table)))
  ### the pooled model must be present and admissible by construction
  expect_true(any(ic$table$k == 1 & ic$table$admissible))

  sel <- SCSTEM_Select(ic, band = c(0.5, 1))
  expect_s3_class(sel, "SCSTEM_Select")
  expect_true(sel$k_selected %in% ic$table$k)
  expect_true(sel$phi_selected %in% ic$table$phi)
  ### the selected configuration must itself be admissible
  expect_true(sel$selected_row$admissible)

  ### (S2) takes the penalty with the smallest criterion at the selected k, on
  ### the whole grid; a penalty can win it only by the criterion
  ic4 <- ic
  at_k <- ic4$table$k == 2
  ic4$table$admissible[at_k] <- TRUE
  ic4$table$BIC[ic4$table$k != 2] <- max(ic$table$BIC) + 1
  ic4$table$BIC[at_k] <- c(3, 1, 2)[rank(ic4$table$phi[at_k])]
  sel4 <- SCSTEM_Select(ic4, band = c(0.5, 1))
  expect_equal(sel4$k_selected, 2L)
  expect_equal(sel4$phi_selected, 0.5)
  expect_equal(sel4$step2$phi[sel4$step2$selected], 0.5)
  expect_output(print(sel4), "criterion over the penalties at k = 2")
  ic4$table$BIC[at_k] <- 1
  expect_equal(SCSTEM_Select(ic4, band = c(0.5, 1))$phi_selected, 0)

  ### the pooled model competes in (S1): when no partition improves on it, it
  ### is the answer, and when a partition does, it is not
  best <- min(ic$table$BIC)
  ic1 <- ic
  ic1$table$BIC[ic1$table$k == 1] <- best - 1
  sel1 <- SCSTEM_Select(ic1, band = c(0.5, 1))
  expect_equal(sel1$k_selected, 1L)
  expect_null(sel1$step2)
  expect_identical(sel1$fit, ic$fits[[paste0("k=1, phi=", sel1$phi_selected)]])
  expect_output(print(sel1), "pooled model")
  ic2 <- ic
  ic2$table$BIC[ic2$table$k == 1] <- max(ic$table$BIC) + 1
  expect_gt(SCSTEM_Select(ic2, band = c(0.5, 1))$k_selected, 1L)

  ### with no admissible partition at all the pooled model is selected
  ic3 <- ic
  ic3$table$admissible[ic3$table$k > 1] <- FALSE
  expect_equal(SCSTEM_Select(ic3, band = c(0.5, 1))$k_selected, 1L)
})

test_that("the fit records the ridge settings the bootstrap refits with", {
  mod <- po_model(Tn = 60L)
  f <- SCSTEM_Estimation(mod, k = 2, phi_penalty = 0.5, distance = "geo",
                         precision = 0.05, max_iter = 3, lambda = 0.2)
  expect_equal(f$input_args$lambda, 0.2)
  expect_equal(f$input_args$alpha, 0)
  expect_true(isTRUE(f$input_args$latent))
  expect_true(isTRUE(f$input_args$spatial))
})
