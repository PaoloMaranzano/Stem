### The four estimation algorithms (R/em-fit.R): the same maximum, the step of
### ECME against an explicit computation, and the settings.

### A small, well-identified model: an AR(1) latent process common to the
### locations, a spatial field and a covariate. (The povalley fixture of
### helper-models.R runs to the boundary of the parameter space, where no
### algorithm converges.)
sim_small <- function(d = 10L, n = 50L, seed = 1L) {
  set.seed(seed)
  co <- matrix(stats::runif(2 * d), d)
  R <- exp(-2 * as.matrix(stats::dist(co)))
  L <- t(chol(0.3 * R + diag(0.2, d)))
  y <- stats::filter(stats::rnorm(n, sd = sqrt(0.3)), 0.7, method = "recursive")
  x <- matrix(stats::rnorm(n * d), n, d)
  z <- 1 + 0.5 * x + matrix(y, n, d) + t(L %*% matrix(stats::rnorm(d * n), d, n))
  STEM_Model(z = z, covariates = cbind(1, as.vector(x)), coordinates = co,
             phi = list(beta = matrix(c(0.8, 0.4), 2, 1), sigma2eps = 0.3, sigma2omega = 0.3,
                        theta = 1, G = matrix(0.5), Sigmaeta = matrix(0.5), m0 = as.matrix(0),
                        C0 = as.matrix(1)),
             A = matrix(1, d, 1))
}

test_that("the algorithm is a validated setting", {
  expect_identical(STEM_control()$algorithm, "SQUAREM")
  expect_identical(STEM_control(algorithm = "ECME")$algorithm, "ECME")
  expect_error(STEM_control(algorithm = "Newton"), "should be one of")
})

test_that("EM, ECME, SQUAREM and SQUAREM-ECME reach the same maximum", {
  mod <- sim_small()
  tight <- list(em_tol_par = 1e-7, em_tol_loglik = 1e-7, em_stop = "all", em_maxit = 2000)
  fits <- lapply(c("EM", "ECME", "SQUAREM", "SQUAREM-ECME"), function(a)
    STEM_Estimation(mod, control = c(tight, algorithm = a)))
  expect_true(all(vapply(fits, function(f) f$estimates$convergence.par$converged, TRUE)))
  ll <- vapply(fits, function(f) f$estimates$loglik, 1)
  expect_lt(max(ll) - min(ll), 1e-6)
  th <- vapply(fits, function(f) f$estimates$phi.hat$theta, 1)
  expect_lt(max(abs(th - th[1])) / th[1], 1e-4)
  ### the two accelerated algorithms take far fewer iterations than the plain ones
  it <- vapply(fits, function(f) f$estimates$convergence.par$iterEM, 1)
  expect_lt(it[3], it[1] / 2)
  expect_lt(it[4], it[2] / 2)
  expect_identical(fits[[4]]$estimates$convergence.par$algorithm, "SQUAREM-ECME")
})

test_that("the reported log-likelihood is that of the reported estimates", {
  mod <- sim_small(d = 8L, n = 40L)
  f <- STEM_Estimation(mod)
  ph <- f$estimates$phi.hat
  dat <- Stem:::stem_em_data(f, "euclidean", Sigmastar.exp)
  phi <- list(A = t(f$skeleton$A), sigma2omega = ph$sigma2omega, logtheta = log(ph$theta),
              logb = log(ph$sigma2eps / ph$sigma2omega), beta = ph$beta, G = ph$G,
              Sigmaeta = ph$Sigmaeta, m0 = ph$m0, C0 = ph$C0)
  expect_equal(Stem:::kalman(phi, dat)$loglik, f$estimates$loglik)
})

test_that("the ECME step is the GLS of beta and m0 under the covariance of the model", {
  ### six locations, 25 periods, scattered gaps and an empty time point: the
  ### augmented filter against the explicit covariance of the stacked data
  set.seed(11)
  d <- 6; n <- 25; r <- 2
  co <- matrix(stats::runif(2 * d), d)
  X <- array(0, c(d, r, n)); X[, 1, ] <- 1; X[, 2, ] <- stats::rnorm(d * n)
  z <- matrix(stats::rnorm(n * d), n, d)
  z[3, 2] <- NA; z[7, ] <- NA; z[10, c(1, 4)] <- NA
  dist <- as.matrix(stats::dist(co))
  phi <- list(A = matrix(1, 1, d), sigma2omega = 0.5, logtheta = log(2), logb = log(0.4),
              beta = matrix(0, r, 1), G = matrix(0.9), Sigmaeta = matrix(0.3), m0 = matrix(0.2),
              C0 = matrix(1))
  dat <- list(zmat = z, covariates = X, dist = dist, obs = Stem:::stem_obs_index(z),
              n = n, d = d, r = r, p = 1L, cov.spat = Sigmastar.exp)
  q <- Stem:::stem_gls_mean(phi, dat)
  got <- solve(q$M, q$v)

  V <- 0.5 * Sigmastar.exp(d = d, logb = log(0.4), logtheta = log(2), dist = dist)
  vt <- numeric(n); v <- 1
  for (t in 1:n) { v <- 0.81 * v + 0.3; vt[t] <- v }
  Py <- matrix(0, n, n)
  for (s in 1:n) for (t in s:n) Py[s, t] <- Py[t, s] <- 0.9^(t - s) * vt[s]
  Om <- kronecker(Py, matrix(1, d, d)) + kronecker(diag(n), V)
  Xs <- do.call(rbind, lapply(1:n, function(t) cbind(X[, , t], 0.9^t)))
  y <- as.vector(t(z)) - rep(0.9^(1:n), each = d) * 0.2      # data net of the current m0
  ok <- !is.na(y)
  Oi <- solve(Om[ok, ok])
  want <- solve(t(Xs[ok, ]) %*% Oi %*% Xs[ok, ], t(Xs[ok, ]) %*% Oi %*% y[ok])
  expect_equal(as.numeric(got), as.numeric(want), tolerance = 1e-10)
})

test_that("a penalized fit runs the plain iterations under SQUAREM", {
  mod <- sim_small(d = 8L, n = 40L)
  ctl <- list(em_tol_par = 1e-6, em_tol_loglik = 1e-6, em_maxit = 40)
  a <- STEM_Estimation(mod, alpha = 0, lambda = 0.2, control = c(ctl, algorithm = "EM"))
  b <- STEM_Estimation(mod, alpha = 0, lambda = 0.2, control = c(ctl, algorithm = "SQUAREM"))
  expect_identical(a$estimates$convergence.par$iterEM, b$estimates$convergence.par$iterEM)
  expect_equal(a$estimates$phi.hat, b$estimates$phi.hat)
})
