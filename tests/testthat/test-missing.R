### Missing values in the response, following Durbin and Koopman (2012, 2nd
### ed.), Sections 2.7 and 4.10.

test_that("missing values are accepted in z but not in covariates or coordinates", {
  s <- po_subset(Tn = 40L, d = 8L)

  z_na <- s$z
  z_na[3, 2] <- NA
  expect_silent(
    STEM_Model(z = z_na, covariates = s$covariates, coordinates = s$coordinates,
               phi = po_phi(), K = matrix(1, s$d, 1))
  )

  cov_na <- s$covariates
  cov_na[5, 2] <- NA
  expect_error(
    STEM_Model(z = s$z, covariates = cov_na, coordinates = s$coordinates,
               phi = po_phi(), K = matrix(1, s$d, 1)),
    "coordinates, covariates"
  )

  coord_na <- s$coordinates
  coord_na[1, 1] <- NA
  expect_error(
    STEM_Model(z = s$z, covariates = s$covariates, coordinates = coord_na,
               phi = po_phi(), K = matrix(1, s$d, 1)),
    "coordinates, covariates"
  )
})

test_that("a location with no observation at all is rejected", {
  s <- po_subset(Tn = 40L, d = 8L)
  z_na <- s$z
  z_na[, 4] <- NA
  expect_error(
    STEM_Model(z = z_na, covariates = s$covariates, coordinates = s$coordinates,
               phi = po_phi(), K = matrix(1, s$d, 1)),
    "no observed value"
  )
})

test_that("complete data are unaffected by the missing-value machinery", {
  ### stem_obs_index() must report a complete panel, and the filter must take
  ### the unrestricted path. This is the guard against a regression in the
  ### complete-data results, which have to stay bit-for-bit what they were.
  s <- po_subset(Tn = 30L, d = 6L)
  ix <- Stem:::stem_obs_index(s$z)
  expect_false(ix$any_missing)
  expect_true(all(ix$complete))
  expect_true(all(is.na(ix$key)))
  expect_equal(unique(ix$n_obs), 6L)
})

test_that("the observed index and the pattern key separate the two extremes", {
  z <- matrix(1:12, nrow = 4, ncol = 3)
  z[2, 2] <- NA          # partially observed
  z[3, ]  <- NA          # nothing observed
  ix <- Stem:::stem_obs_index(z)

  expect_true(ix$any_missing)
  expect_equal(ix$n_obs, c(3L, 2L, 0L, 3L))
  expect_equal(ix$idx[[2]], c(1L, 3L))
  expect_equal(ix$idx[[3]], integer(0))
  ### a complete row and a fully missing row must never share a cache key
  expect_true(is.na(ix$key[1]))
  expect_false(is.na(ix$key[3]))
  expect_false(identical(ix$key[2], ix$key[3]))
})

test_that("the conditional blocks reproduce Gaussian conditioning", {
  set.seed(1)
  d <- 5
  A <- matrix(stats::rnorm(d * d), d, d)
  Sigma <- crossprod(A) + diag(d)
  oi <- c(1L, 2L, 4L)
  mi <- c(3L, 5L)

  bl <- Stem:::stem_missing_blocks(Sigma, oi, d)

  expect_equal(bl$mi, mi)
  expect_equal(bl$P, Sigma[mi, oi] %*% solve(Sigma[oi, oi]))
  expect_equal(bl$Omega[mi, mi],
               Sigma[mi, mi] - bl$P %*% Sigma[oi, mi])
  ### the observed block of Omega carries no residual uncertainty
  expect_equal(bl$Omega[oi, oi], matrix(0, length(oi), length(oi)))
  ### S is the identity on the observed rows and P on the missing ones
  expect_equal(bl$S[oi, ], diag(1, length(oi)))
  expect_equal(bl$S[mi, ], bl$P)
})

test_that("with nothing observed the error keeps its whole variance", {
  set.seed(2)
  d <- 4
  A <- matrix(stats::rnorm(d * d), d, d)
  Sigma <- crossprod(A) + diag(d)

  bl <- Stem:::stem_missing_blocks(Sigma, integer(0), d)
  expect_equal(bl$Omega, Sigma)
  expect_equal(ncol(bl$S), 0L)
})

test_that("with nothing missing the blocks are the identity and zero", {
  set.seed(3)
  d <- 4
  A <- matrix(stats::rnorm(d * d), d, d)
  Sigma <- crossprod(A) + diag(d)

  bl <- Stem:::stem_missing_blocks(Sigma, seq_len(d), d)
  expect_equal(bl$S, diag(1, d))
  expect_equal(bl$Omega, matrix(0, d, d))
  expect_equal(bl$mi, integer(0))
})

test_that("estimation runs with scattered gaps and with an empty time point", {
  s <- po_subset(Tn = 45L, d = 8L)

  set.seed(4321)
  z_na <- s$z
  z_na[sample(length(z_na), round(0.10 * length(z_na)))] <- NA
  z_na[7, ] <- NA                      # a time point with nothing observed
  ### keep every location with at least one observation
  for (j in seq_len(ncol(z_na))) if (all(is.na(z_na[, j]))) z_na[1, j] <- s$z[1, j]

  mod <- STEM_Model(z = z_na, covariates = s$covariates,
                    coordinates = s$coordinates,
                    phi = po_phi(), K = matrix(1, s$d, 1))
  fit <- STEM_Estimation(mod, precision = 0.5, max.iter = 3)

  p <- fit$estimates$phi.hat
  expect_true(all(is.finite(unlist(p[c("sigma2eps", "sigma2omega", "theta")]))))
  expect_true(all(is.finite(as.numeric(p$beta))))
  expect_true(is.finite(utils::tail(fit$estimates$loglik, 1)))
})

test_that("the assignment score uses only the observed time points", {
  set.seed(5)
  Tn <- 20
  X  <- cbind(1, stats::rnorm(Tn))
  b  <- matrix(c(0.5, 1.5), 2, 1)
  ysm <- matrix(stats::rnorm(Tn), Tn, 1)
  Ki <- 1
  z  <- as.numeric(X %*% b) + as.numeric(ysm) + stats::rnorm(Tn, sd = 0.3)

  full <- Stem:::scstem_loglike_i(z, X, b, ysm, Ki, 0.05, 0.04)

  ### blanking a time point must give exactly the score of the shortened series
  z_na <- z; z_na[5] <- NA
  gapped <- Stem:::scstem_loglike_i(z_na, X, b, ysm, Ki, 0.05, 0.04)
  manual <- Stem:::scstem_loglike_i(z[-5], X[-5, , drop = FALSE], b,
                                    ysm[-5, , drop = FALSE], Ki, 0.05, 0.04)
  expect_equal(gapped, manual)
  expect_false(isTRUE(all.equal(gapped, full)))

  ### a location observed nowhere cannot be scored
  expect_equal(Stem:::scstem_loglike_i(rep(NA_real_, Tn), X, b, ysm, Ki, 0.05, 0.04),
               -Inf)
})

test_that("SC-STEM estimation tolerates gaps in the response", {
  s <- po_subset(Tn = 45L, d = 14L)
  set.seed(99)
  z_na <- s$z
  z_na[sample(length(z_na), round(0.08 * length(z_na)))] <- NA
  for (j in seq_len(ncol(z_na))) if (all(is.na(z_na[, j]))) z_na[1, j] <- s$z[1, j]

  mod <- STEM_Model(z = z_na, covariates = s$covariates,
                    coordinates = s$coordinates,
                    phi = po_phi(), K = matrix(1, s$d, 1))

  fit <- SCSTEM_Estim(mod, k = 2, phi_penalty = 0.5, knn = 3,
                      precision = 0.5, precision_full_dataset = 0.5,
                      max_iter = 2, seed = 1)

  expect_s3_class(fit, "SCSTEM_Estim")
  expect_length(fit$group, s$d)
  expect_true(all(fit$group %in% 1:2))
})

test_that("STEM_Fitted reproduces the observations and predicts the blanks", {
  s <- po_subset(Tn = 40L, d = 8L)

  set.seed(77)
  cells <- sample(length(s$z), 25)
  truth <- s$z[cells]
  z_na <- s$z; z_na[cells] <- NA

  mod <- STEM_Model(z = z_na, covariates = s$covariates,
                    coordinates = s$coordinates,
                    phi = po_phi(), K = matrix(1, s$d, 1))
  fit <- STEM_Estimation(mod, precision = 0.5, max.iter = 3)
  zhat <- STEM_Fitted(fit)

  expect_equal(dim(zhat), dim(s$z))
  ### wherever the response was observed the completion returns it unchanged
  obs <- !is.na(z_na)
  expect_equal(zhat[obs], z_na[obs])
  ### and the blanks are filled with finite predictions
  expect_true(all(is.finite(zhat[cells])))
  ### which beat the station means, the natural naive alternative
  naive <- matrix(colMeans(z_na, na.rm = TRUE), nrow(z_na), ncol(z_na),
                  byrow = TRUE)[cells]
  expect_lt(sqrt(mean((zhat[cells] - truth)^2)),
            sqrt(mean((naive - truth)^2)))
})

test_that("a fully missing time point falls back on the signal", {
  s <- po_subset(Tn = 40L, d = 8L)
  z_na <- s$z
  z_na[11, ] <- NA

  mod <- STEM_Model(z = z_na, covariates = s$covariates,
                    coordinates = s$coordinates,
                    phi = po_phi(), K = matrix(1, s$d, 1))
  fit <- STEM_Estimation(mod, precision = 0.5, max.iter = 3)
  zhat <- STEM_Fitted(fit)

  ### with nothing observed there is no spatial correction, so the prediction
  ### is exactly the signal x'beta + K yhat
  XX <- Stem:::changedimension_covariates(s$covariates, s$d, ncol(s$covariates),
                                          nrow(z_na))
  b <- matrix(as.numeric(fit$estimates$phi.hat$beta), ncol = 1)
  ysm <- as.numeric(fit$estimates$y.smoothed)
  signal <- as.numeric(XX[, , 11] %*% b) + ysm[11]
  expect_equal(as.numeric(zhat[11, ]), signal)
})
