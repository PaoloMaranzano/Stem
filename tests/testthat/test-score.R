## The assignment scores of the label step. scstem_loglike_i() is the marginal
## score; scstem_cond_scores() returns the conditional score and the
## expected-gap term. These tests pin the second against a brute-force
## computation from the joint Gaussian density, which is its definition.

cond_toy <- function(d = 7L, Tn = 25L, seed = 1L, na_share = 0) {
  set.seed(seed)
  coords <- cbind(stats::runif(d), stats::runif(d))
  dm <- as.matrix(stats::dist(coords, diag = TRUE))
  X <- cbind(1, stats::rnorm(d * Tn))
  beta <- c(0.5, -1)
  ysm <- matrix(stats::rnorm(Tn), ncol = 1)
  Kmat <- matrix(1, d, 1)
  s2e <- 0.4; s2o <- 1.1; theta <- 2.5
  z <- matrix(as.numeric(X %*% beta), nrow = Tn) +
    as.numeric(ysm) + matrix(stats::rnorm(d * Tn), nrow = Tn)
  if (na_share > 0) z[sample(length(z), round(na_share * length(z)))] <- NA
  list(z = z, X = X, beta = beta, ysm = ysm, Kmat = Kmat, s2e = s2e, s2o = s2o,
       theta = theta, dm = dm, d = d, Tn = Tn)
}

## the definition: log p(r_i | r_J) = log N(r_{J u i}) - log N(r_J), period by
## period, J being the members observed at that period other than i
brute_cond <- function(tt, i, members) {
  R <- tt$z - matrix(as.numeric(tt$X %*% tt$beta), nrow = tt$Tn) - as.numeric(tt$ysm)
  Sig <- tt$s2o * exp(-tt$theta * tt$dm); diag(Sig) <- tt$s2e + tt$s2o
  ldn <- function(x, S) {
    if (!length(x)) return(0)
    U <- chol(S)
    -0.5 * length(x) * log(2 * pi) - sum(log(diag(U))) -
      0.5 * sum(backsolve(U, x, transpose = TRUE)^2)
  }
  out <- 0
  for (t in seq_len(tt$Tn)) {
    if (is.na(R[t, i])) next
    J <- setdiff(members, i)
    J <- J[!is.na(R[t, J])]
    out <- out + ldn(R[t, c(J, i)], Sig[c(J, i), c(J, i), drop = FALSE]) -
      ldn(R[t, J], Sig[J, J, drop = FALSE])
  }
  out
}

test_that("the conditional score is the conditional Gaussian log-density", {
  cond_scores <- getFromNamespace("scstem_cond_scores", "Stem")
  for (na in c(0, 0.15)) {
    tt <- cond_toy(na_share = na)
    members <- c(1L, 2L, 4L, 6L)
    cs <- cond_scores(tt$z, tt$X, tt$Tn, tt$beta, tt$ysm, tt$Kmat, tt$s2e, tt$s2o,
                      tt$theta, members, tt$dm, spatial = TRUE)
    for (i in seq_len(tt$d)) {
      expect_equal(cs$cond[i], brute_cond(tt, i, members), tolerance = 1e-9)
    }
    expect_true(all(cs$delta >= 0))
  }
})

test_that("without spatial correlation the conditional score is the marginal one", {
  cond_scores <- getFromNamespace("scstem_cond_scores", "Stem")
  loglike_i <- getFromNamespace("scstem_loglike_i", "Stem")
  tt <- cond_toy(na_share = 0.1)
  cs <- cond_scores(tt$z, tt$X, tt$Tn, tt$beta, tt$ysm, tt$Kmat, tt$s2e, tt$s2o,
                    tt$theta, c(1L, 3L, 5L), tt$dm, spatial = FALSE)
  marg <- vapply(seq_len(tt$d), function(i) {
    rows <- ((i - 1) * tt$Tn + 1):(i * tt$Tn)
    loglike_i(tt$z[, i], tt$X[rows, , drop = FALSE], tt$beta, tt$ysm,
              tt$Kmat[i, , drop = FALSE], tt$s2e, tt$s2o)
  }, numeric(1))
  expect_equal(cs$cond, marg, tolerance = 1e-10)
  expect_equal(cs$delta, rep(0, tt$d))
})

test_that("the expected gap is the Schur-complement formula", {
  cond_scores <- getFromNamespace("scstem_cond_scores", "Stem")
  tt <- cond_toy()
  members <- c(2L, 3L, 5L, 7L)
  cs <- cond_scores(tt$z, tt$X, tt$Tn, tt$beta, tt$ysm, tt$Kmat, tt$s2e, tt$s2o,
                    tt$theta, members, tt$dm, spatial = TRUE)
  Sig <- tt$s2o * exp(-tt$theta * tt$dm); diag(Sig) <- tt$s2e + tt$s2o
  Rm <- stats::cov2cor(Sig)
  for (i in seq_len(tt$d)) {
    J <- setdiff(members, i)
    rho <- Rm[J, i]
    q <- as.numeric(t(rho) %*% solve(Rm[J, J], rho))
    expect_equal(cs$delta[i], -0.5 * tt$Tn * log(1 - q), tolerance = 1e-9)
  }
})

test_that("the default score is unchanged and the alternatives run", {
  mod <- po_model(Tn = 60L)
  a <- SCSTEM_Estimation(mod, k = 2, phi_penalty = 0.5, distance = "geo",
                         precision = 0.05, max_iter = 4)
  b <- SCSTEM_Estimation(mod, k = 2, phi_penalty = 0.5, distance = "geo",
                         precision = 0.05, max_iter = 4, score = "marginal")
  expect_identical(a$group, b$group)
  expect_equal(a$input_args$score, "marginal")
  for (s in c("conditional", "corrected")) {
    f <- SCSTEM_Estimation(mod, k = 2, phi_penalty = 0.5, distance = "geo",
                           precision = 0.05, max_iter = 4, score = s)
    expect_true(all(tabulate(f$group, nbins = 2) >= f$input_args$min_cluster_size))
    expect_equal(f$input_args$score, s)
  }
})

test_that("the fit records the settings the bootstrap must reproduce", {
  mod <- po_model(Tn = 60L)
  f <- SCSTEM_Estimation(mod, k = 2, phi_penalty = 0.5, distance = "geo",
                         precision = 0.05, max_iter = 3, lambda = 0.2)
  expect_equal(f$input_args$lambda, 0.2)
  expect_equal(f$input_args$alpha, 0)
  expect_true(isTRUE(f$input_args$latent))
  expect_true(isTRUE(f$input_args$spatial))
})

test_that("the label step never lowers the objective, whatever the score", {
  mod <- po_model(Tn = 60L)
  for (s in c("marginal", "conditional")) {
    f <- SCSTEM_Estimation(mod, k = 3, phi_penalty = 0.5, distance = "geo",
                           precision = 0.05, max_iter = 5, score = s)
    tr <- f$obj_trace
    expect_true(all(tr$objective >= tr$objective_before - 1e-8))
    expect_equal(dim(f$score_last), c(ncol(mod$data$z), 3L))
  }
})
