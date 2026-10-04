### The quantities the package reports, against independent dense
### computations: the log-likelihood, the smoothed states, the maximum, and the
### number of observations of the BIC. Until 2.0.0 the package returned four
### times the log-likelihood, a defect that checks of internal consistency
### alone could not see.

### A small model with a few gaps, and its covariance written out in full; its
### maximum is interior (a draw whose Sigmaeta runs to zero would not do)
dense_case <- function(seed = 3L, d = 7L, n = 30L) {
  set.seed(seed)
  co <- matrix(stats::runif(2 * d), d)
  h <- as.matrix(stats::dist(co))
  L <- t(chol(0.3 * exp(-2 * h) + diag(0.2, d)))
  y <- stats::filter(stats::rnorm(n, sd = sqrt(0.3)), 0.7, method = "recursive")
  x <- matrix(stats::rnorm(n * d), n, d)
  z <- 1 + 0.5 * x + matrix(y, n, d) + t(L %*% matrix(stats::rnorm(d * n), d, n))
  z[sample(length(z), 12)] <- NA
  X <- cbind(1, as.vector(x))
  mod <- STEM_Model(z = z, covariates = X, coordinates = co,
                    phi = list(beta = matrix(c(0.8, 0.4), 2, 1), sigma2eps = 0.3, sigma2omega = 0.3,
                               theta = 1, G = matrix(0.5), Sigmaeta = matrix(0.5), m0 = as.matrix(0),
                               C0 = as.matrix(1)),
                    A = matrix(1, d, 1))
  list(mod = mod, z = z, X = X, h = h, d = d, n = n)
}

### the log-likelihood of the observed values and E[y_t | z], from the joint
### covariance of the observations stacked by time
dense_fit <- function(cs, ph) {
  d <- cs$d; n <- cs$n
  G <- ph$G[1, 1]; Q <- ph$Sigmaeta[1, 1]
  Se <- ph$sigma2eps * diag(d) + ph$sigma2omega * exp(-ph$theta * cs$h)
  V <- numeric(n); v <- ph$C0[1, 1]
  for (t in seq_len(n)) { v <- G^2 * v + Q; V[t] <- v }
  Cy <- outer(seq_len(n), seq_len(n), function(t, s) G^abs(t - s) * V[pmin(t, s)])
  Sig <- kronecker(Cy, matrix(1, d, d)) + kronecker(diag(n), Se)
  mu_y <- G^seq_len(n) * ph$m0[1, 1]
  xb <- vapply(seq_len(n), function(t)
    as.vector(cs$X[(seq_len(d) - 1L) * n + t, , drop = FALSE] %*% ph$beta), numeric(d))
  mu <- as.vector(xb + matrix(mu_y, d, n, byrow = TRUE))
  zz <- as.vector(t(cs$z)); ok <- !is.na(zz)
  U <- chol(Sig[ok, ok]); r <- zz[ok] - mu[ok]
  ll <- -0.5 * (sum(ok) * log(2 * pi) + 2 * sum(log(diag(U))) +
                  sum(backsolve(U, r, transpose = TRUE)^2))
  Cyz <- kronecker(Cy, matrix(1, 1, d))[, ok, drop = FALSE]
  list(loglik = ll, ey = mu_y + as.vector(Cyz %*% chol2inv(U) %*% r))
}

test_that("the log-likelihood is the Gaussian log-likelihood of the observed values", {
  cs <- dense_case()
  f <- STEM_Estimation(cs$mod)
  ex <- dense_fit(cs, f$estimates$phi.hat)
  expect_equal(f$estimates$loglik, ex$loglik, tolerance = 1e-10)
})

test_that("the smoothed states are the conditional means of the states", {
  cs <- dense_case()
  f <- STEM_Estimation(cs$mod)
  ex <- dense_fit(cs, f$estimates$phi.hat)
  expect_equal(as.numeric(f$estimates$y.smoothed), ex$ey, tolerance = 1e-10)
})

test_that("the estimates are a maximum of the exact log-likelihood", {
  cs <- dense_case()
  f <- STEM_Estimation(cs$mod, control = list(em_tol_par = 1e-9, em_tol_loglik = 1e-9,
                                              em_stop = "all", em_maxit = 5000))
  ph <- f$estimates$phi.hat
  l0 <- dense_fit(cs, ph)$loglik
  ### no coordinate move of 1% raises it
  moves <- list(
    function(p, s) { p$beta[1] <- p$beta[1] * s; p }, function(p, s) { p$beta[2] <- p$beta[2] * s; p },
    function(p, s) { p$sigma2eps <- p$sigma2eps * s; p }, function(p, s) { p$sigma2omega <- p$sigma2omega * s; p },
    function(p, s) { p$theta <- p$theta * s; p }, function(p, s) { p$G <- p$G * s; p },
    function(p, s) { p$Sigmaeta <- p$Sigmaeta * s; p }, function(p, s) { p$m0 <- p$m0 * s; p })
  gain <- vapply(moves, function(mv) max(dense_fit(cs, mv(ph, 1.01))$loglik,
                                         dense_fit(cs, mv(ph, 0.99))$loglik) - l0, 1)
  expect_true(all(gain < 1e-8))
})

test_that("the BIC counts the observed values", {
  cs <- dense_case()
  fit <- SCSTEM_Estimation(cs$mod, K = 1, distance = "euclidean")
  ic <- fit$info_crit
  expect_equal(unname(ic[["BIC"]]),
               unname(-2 * ic[["loglik"]] + log(sum(!is.na(cs$z))) * ic[["df"]]))
})
