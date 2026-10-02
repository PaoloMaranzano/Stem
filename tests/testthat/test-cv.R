test_that("a removed location is predicted by kriging within its regime", {
  s <- po_subset(Tn = 60L)
  d <- s$d
  drop <- c(3L, 10L, 25L)
  keep <- setdiff(seq_len(d), drop)
  rows <- function(idx) as.vector(outer(seq_len(s$Tn), (idx - 1L) * s$Tn, "+"))
  mk <- function(idx) STEM_Model(z = s$z[, idx], covariates = s$covariates[rows(idx), , drop = FALSE],
                                 coordinates = s$coordinates[idx, , drop = FALSE],
                                 phi = po_phi(), A = matrix(1, length(idx), 1))
  fit <- SCSTEM_Estimation(mk(keep), K = 1, distance = "geo", precision = 0.05)

  ### the same conditional mean through STEM_Complete(): the fit carried over to
  ### the whole network, with the removed columns missing
  mfull <- mk(seq_len(d))
  full <- fit$fit_list[[1]]
  full$data <- mfull$data
  full$data$z[, drop] <- NA
  full$skeleton$A <- mfull$skeleton$A
  ref <- STEM_Complete(full, distance = "geo")[, drop]

  dm <- as.matrix(geodist::geodist(s$coordinates, measure = "geodesic"))
  XX <- changedimension_covariates(mfull$data$covariates, d, ncol(mfull$data$covariates), s$Tn)
  got <- scstem_predict_new(fit, keep, drop, s$z[, keep], dm, XX, mfull$skeleton$A)
  expect_equal(got, unname(ref), tolerance = 1e-8)
})

test_that("SCSTEM_CV evaluates and ranks fits on the same folds", {
  skip_on_cran()

  s <- po_subset(Tn = 60L)
  full <- SCSTEM_Estimation(po_model(Tn = 60L), K = 1, distance = "geo", precision = 0.05)
  ### a model with fewer covariates on the same data
  phi2 <- po_phi()
  phi2$beta <- phi2$beta[1:2, , drop = FALSE]
  mod2 <- STEM_Model(z = s$z, covariates = s$covariates[, 1:2, drop = FALSE],
                     coordinates = s$coordinates, phi = phi2, A = matrix(1, s$d, 1))
  reduced <- SCSTEM_Estimation(mod2, K = 1, distance = "geo", precision = 0.05)

  cv <- SCSTEM_CV(list(full = full, reduced = reduced), scheme = c("LKLO", "LKLHTO"),
                  folds = 2, seed = 1)
  expect_s3_class(cv, "SCSTEM_CV")
  expect_equal(nrow(cv$summary), 4L)
  expect_true(all(is.na(cv$folds$message)))
  expect_true(all(cv$summary$n > 0))
  for (sch in c("LKLO", "LKLHTO")) {
    expect_setequal(cv$summary$rank[cv$summary$scheme == sch], 1:2)
  }
  expect_output(print(cv), "Blocked cross-validation")

  ### the folds depend only on the seed and the scheme
  cv1 <- SCSTEM_CV(full, scheme = "LKLO", folds = 2, seed = 1)
  expect_equal(cv1$folds$rmse,
               cv$folds$rmse[cv$folds$scheme == "LKLO" & cv$folds$model == "full"])

  ### models on different data cannot be compared
  mod3 <- STEM_Model(z = s$z + 1, covariates = s$covariates, coordinates = s$coordinates,
                     phi = po_phi(), A = matrix(1, s$d, 1))
  other <- SCSTEM_Estimation(mod3, K = 1, distance = "geo", precision = 0.05)
  expect_error(SCSTEM_CV(list(a = full, b = other), folds = 2, seed = 1),
               "share the response")
})
