### Regularization of the regression coefficients.
###
### The property that matters most is the first one: with lambda = 0 nothing
### about the estimator may change, bit for bit. The rest checks the two closed
### forms against independent algebra and the iterative one against its own
### optimality conditions.

make_model <- function(dn = 24L, TT = 90L, seed = 1L) {
  set.seed(seed)
  co <- cbind(stats::runif(dn, 7.5, 13.5), stats::runif(dn, 44.7, 46.1))
  dm <- geodist::geodist(co, measure = "geodesic")
  L  <- chol(18 * diag(dn) + 5 * exp(-dm / 120000))
  y  <- as.numeric(stats::arima.sim(list(ar = 0.9), TT)) * 2
  X  <- cbind(1, rep(stats::runif(dn, 0, 250), each = TT),
              stats::rnorm(TT * dn, 30, 12))
  Xa <- changedimension_covariates(X, dn, 3L, TT)
  b  <- matrix(c(2.3, 0.004, 0.58), 3, 1)
  z  <- matrix(NA_real_, TT, dn)
  for (tt in seq_len(TT))
    z[tt, ] <- as.numeric(Xa[, , tt] %*% b) + y[tt] +
               as.numeric(crossprod(L, stats::rnorm(dn)))
  ols <- stats::lm.fit(x = X, y = as.vector(z))
  s2  <- stats::var(ols$residuals)
  STEM_Model(z = z, covariates = X, coordinates = co,
    phi = list(beta = matrix(ols$coefficients, ncol = 1),
               sigma2eps = 0.6 * s2, sigma2omega = 0.3 * s2, theta = 1 / 60000,
               G = matrix(0.7, 1, 1), Sigmaeta = matrix(0.1 * s2, 1, 1),
               m0 = as.matrix(0), C0 = as.matrix(1)),
    K = matrix(1, dn, 1))
}

test_that("lambda = 0 reproduces the unpenalized estimator exactly", {
  a <- STEM_Estimation(make_model(), precision = 0.01, max.iter = 6,
                       distance = "geo")
  b <- STEM_Fit(make_model(), k = 1, alpha = 0, lambda = 0,
                precision = 0.01, max.iter = 6, distance = "geo")
  expect_identical(unlist(a$estimates$phi.hat), unlist(b$estimates$phi.hat))
  expect_identical(a$estimates$loglik, b$estimates$loglik)
  expect_equal(b$estimates$penalty$beta.df, ncol(a$data$covariates))
})

test_that("the ridge update matches an independent closed form", {
  set.seed(3); r <- 6L
  A <- matrix(stats::rnorm(40 * r), 40, r)
  M <- crossprod(A); v <- as.numeric(crossprod(A, stats::rnorm(40)))
  w <- stem_penalized_index(r)
  expect_identical(w, c(0, 1, 1, 1, 1, 1))

  lam <- 2.5
  s <- sqrt(diag(M))
  ref <- solve(M / tcrossprod(s) + diag(lam * w, r), v / s) / s
  got <- stem_beta_update(M, v, alpha = 0, lambda = lam, w = w)$beta
  expect_equal(got, ref, tolerance = 1e-12)

  ## and with no penalty it is the ordinary solve
  expect_equal(stem_beta_update(M, v, lambda = 0, ridge_reg = 0.01)$beta,
               as.numeric(solve(diag(0.01, r) + M, v)), tolerance = 1e-12)
})

test_that("the lasso and elastic-net updates satisfy their optimality conditions", {
  set.seed(5); r <- 6L
  A <- matrix(stats::rnorm(40 * r), 40, r)
  M <- crossprod(A); v <- as.numeric(crossprod(A, stats::rnorm(40)))
  w <- stem_penalized_index(r)
  s <- sqrt(diag(M)); Ms <- M / tcrossprod(s); vs <- v / s

  for (al in c(1, 0.5)) {
    lam <- 1.2
    b  <- stem_beta_update(M, v, alpha = al, lambda = lam, w = w)$beta
    bs <- b * s
    g  <- as.numeric(vs - Ms %*% bs - lam * (1 - al) * w * bs)
    kkt <- vapply(seq_len(r), function(j) {
      if (w[j] == 0) abs(g[j])
      else if (abs(bs[j]) > 1e-10) abs(g[j] - lam * al * sign(bs[j]))
      else max(0, abs(g[j]) - lam * al)
    }, numeric(1))
    expect_lt(max(kkt), 1e-8)
  }
})

test_that("the lasso shrinks coefficients to exactly zero and the count follows", {
  fit <- STEM_Fit(make_model(), k = 1, alpha = 1, lambda = 5,
                  precision = 0.01, max.iter = 6, distance = "geo")
  b <- as.numeric(fit$estimates$phi.hat$beta)
  expect_true(any(b[-1] == 0))
  ## the intercept is never penalized, so it survives
  expect_true(b[1] != 0)
  expect_equal(fit$estimates$penalty$beta.df, sum(b != 0))
  expect_lt(fit$estimates$penalty$beta.df, 3)
})

test_that("STEM_Fit dispatches on k and validates its arguments", {
  m <- make_model()
  expect_s3_class(STEM_Fit(m, k = 1, precision = 0.05, max.iter = 3,
                           distance = "geo"), "STEM_Model")
  expect_error(STEM_Fit(m, alpha = 2), "alpha")
  expect_error(STEM_Fit(m, lambda = -1), "lambda")
  expect_error(STEM_Fit(m, k = 0), "k")
  expect_error(STEM_Fit(list()), "STEM_Model")
})
