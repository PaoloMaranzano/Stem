# Stem 2.0.0

First release of the package after its archival on CRAN, and the first to
include the spatially-clustered STEM (SC-STEM) model family. The classical STEM
workflow of version 1.0 is unchanged in its statistical content; what changes is
its numerical robustness, the whole SC-STEM layer, and the packaging.

## Breaking changes: names aligned with the paper and with SCDA

* The loading matrix is `A` (it was `K`): `STEM_Model(..., A = )`,
  `skeleton$A`, `STEM_Kriging(..., A.newlocations = )`. Scripts written for
  version 1.0 must replace `K =` by `A =`.
* The number of regimes is `K`, with index `k = 1, ..., K`, as in the paper and
  as in the SCDA package: `STEM_Fit(K = )`, `SCSTEM_Estimation(K = )`,
  `SCSTEM_Infocrit(K_grid = )`, `SCSTEM_Select()$K_selected`, and the columns
  `K`, `K_eff` of the grid. The transition matrix stays `G`.
* The number of free parameters is `df`: `info_crit["df"]` (it was `"k"`) and
  the column `df` of the grid (it was `npar`).
* The deprecated arguments `mink` and `maxk` of `SCSTEM_Infocrit()` are removed.
* There are no aliases: an old name stops with a message naming the new one.
  This matters for `k`, which R would otherwise complete silently to `knn`.

## New features

### Regularized regression coefficients

* `STEM_Fit()` is a single entry point for the pooled and the clustered model,
  with or without a penalty on the regression coefficients. `k = 1` fits the
  pooled STEM model and `k > 1` the spatially-clustered one; `lambda = 0`, the
  default, gives the ordinary maximum likelihood estimator and reproduces
  `STEM_Estimation()` exactly.
* `alpha` and `lambda` follow the `glmnet` parameterization: `alpha = 0` is
  ridge, `alpha = 1` the lasso, anything between the elastic net. Only the
  regression coefficients are penalized; the intercept is excluded by default.
  In the clustered model the penalty acts within each regime, with one
  `(alpha, lambda)` shared by all of them.
* `lambda` is dimensionless. Its L1 part is measured against the largest partial
  gradient, so `lambda` in `(0, 1]` traverses the whole lasso path and means the
  same thing at any error variance and in any regime; its L2 part is left alone,
  being already unit-free. `lambda_scale = "absolute"` recovers the raw
  convention.
* `lambda_by = "size"` spreads the penalty over the regimes in proportion to
  `1/n_g` instead of equally, shrinking the smaller regimes more. It remains one
  hyperparameter.
* `latent = FALSE` switches the latent process off and `spatial = FALSE`
  replaces the spatial correlation by the identity. Together with
  `regularization = 0` they reduce the model exactly to penalized linear
  regression, and with `k > 1` to clusterwise penalized regression.
* The penalty leaves the E-step untouched, so the algorithm remains an EM on the
  penalized likelihood: the M-step keeps a closed form under a ridge and is
  solved exactly by coordinate descent otherwise. The design is scaled
  internally in the generalized least squares metric the model works in, so a
  given `lambda` means the same thing for every covariate and every regime.
* The information criteria count the effective number of coefficients the
  penalty leaves rather than the nominal one.
* The theory is written up in `dev/regularization/stem-elastic-net.tex` in the
  source repository.

### Spatially-clustered STEM models

* `SCSTEM_Estimation()` fits an SC-STEM model: the monitoring locations are
  partitioned into `k` latent spatial regimes and a separate STEM model is
  estimated within each of them, so that regression coefficients, variance
  components and latent temporal dynamics are all cluster-specific. Labels and
  parameters are estimated jointly by maximizing a Potts-penalized
  log-likelihood.
* `SCSTEM_Infocrit()` now fits a full `(k, phi)` grid and returns exact
  log-likelihoods, AIC, BIC and KIC, an admissibility flag, the cluster sizes
  and the estimated partitions. The pre-2.0.0 `mink`/`maxk` calling convention
  still works.
* `SCSTEM_Select()` implements a two-step rule for choosing the
  hyperparameters: **(S1)** the modal BIC-minimizing `k` inside a
  moderate-penalty band, ties resolved towards the smaller `k`; **(S2)** the
  `phi` with the smallest BIC at that `k` over the whole grid, ties resolved
  towards the smaller `phi`. Only admissible configurations enter the rule.
  The pooled `k = 1` model competes in (S1) and is returned when no partition
  improves on it.
* `SCSTEM_Bootstrap()` is now a **refit-with-clustering** parametric bootstrap:
  data are generated cluster by cluster from the fitted model and the *entire*
  procedure, endogenous partitioning included, is re-estimated on every draw, so
  that the uncertainty of the partition is propagated. It replaces the previous
  bootstrap, which conditioned on the estimated partition and therefore
  understated the uncertainty. The refits use the same `alpha`, `lambda`,
  `penalize`, `lambda_scale`, `lambda_by`, `latent` and `spatial` as the
  original fit, so a fit with a ridge on the coefficients is resampled with it.
* `SCSTEM_BootInference()` aligns every refit onto the original clusters by the
  majority rule and returns bootstrap standard errors; normal, basic,
  percentile and bias-corrected confidence intervals; pairwise percentile tests
  for the differences between clusters; the co-clustering matrix; and the ARI of
  each refit against the original partition.
* `SCSTEM_CV()` is an optional validation of predictive accuracy, separate
  from the choice of `(k, phi)`: blocked cross-validation leaving out
  locations (`"LKLO"`), time blocks (`"LKTO"`), both (`"LKLHTO"`) or random
  cells, of one fit or of a named list of fits on the same data (the pooled
  model, the selected one, models with different covariates), all refitted on
  the same folds and ranked by root mean squared error. A removed location
  inherits the regime of its nearest retained location and is kriged within it.
* New classes `SCSTEM_Estimation`, `SCSTEM_Infocrit`, `SCSTEM_Select`,
  `SCSTEM_Bootstrap`, `SCSTEM_BootInference` and `SCSTEM_CV`, each with a
  `print()` method.

### Algorithmic changes in the SC-STEM assignment step

* The label update is now **ICM** (Iterated Conditional Modes, Besag 1986) by
  default: locations are visited sequentially and the Potts penalty is
  recomputed on the fly, so the sweep cannot decrease the objective and cannot
  cycle. The previous simultaneous update is still available via
  `label_update = "simultaneous"`.
* The penalized objective
  `Q = sum_i l_{i,k_i} + phi * c * #{concordant neighbor pairs}` is now
  computed explicitly and traced along the iterations (`obj_trace`).
  Convergence is declared on label stability, on the improvement of `Q`, on
  cycle detection (the best visited partition is returned) or at `max_iter`.
* **Minimum-size constraint** (`enforce_min_size = TRUE`). Within-cluster
  homogeneity is what the assignment step seeks, so the cluster-wise variance
  components shrink and, left unconstrained, the cluster with the smallest
  residual variance attracts every location. On the `pm10` example the
  unconstrained sweep collapses to a single cluster even at `phi = 0`. A
  location may now leave its cluster only if that cluster stays at or above
  `min_cluster_size`.
* **Size-preserving swap pass** (`swap_pass = TRUE`). Since the constraint above
  freezes any location sitting in a minimum-size cluster, each sweep is followed
  by a pass that exchanges the labels of two locations whenever this strictly
  increases `Q`. Swaps leave cluster sizes unchanged, so feasibility and
  monotonicity both hold. An exact pruning bound keeps the scan affordable.
* The `knn` graph is **symmetrized**: the Potts penalty is defined on an
  undirected graph, whereas `spdep::knearneigh()` returns an asymmetric one.
* The `knn` graph is built with the **same metric as the covariance**. The
  nearest neighbors used to be taken on raw coordinates, so on longitude and
  latitude they were the neighbors of a planar metric in which a degree of
  longitude and a degree of latitude count the same; at 45 degrees of latitude
  the first is about 78 km and the second about 111 km, so the graph preferred
  north-south neighbors while `Sigma_e` was measured on the sphere. The
  `distance` argument now governs both.
* New `phi_scale` argument. Each location contributes `T` observations to the
  likelihood, so a penalty calibrated for cross-sectional models is not
  transferable. The default `"auto"` normalizes the penalty by the median spread
  of the location-wise log-likelihood contributions; `"per-observation"` and
  `"raw"` are also available. The penalty actually applied is reported in
  `phi_effective`.
* The partition of an unpenalized fit is now initialized, by default, on the
  **departures of each location from the pooled model**
  (`init_method = "departures"`): the mean of its residual from the pooled
  signal, the slopes of that residual on its covariates, and the lag-one
  autocorrelation and log-variance of what the slopes leave. The covariate
  means, the previous default, carry no information on the regimes when the
  covariates are exogenous to them; they remain available as
  `init_method = "kmeans"`. `SCSTEM_Infocrit()` computes the departures once,
  from its fit at `k = 1`.
* A fit with `phi_penalty > 0` **starts from the solution of the unpenalized
  fit** at the same `k` (`SCSTEM_Estimation()` fits `phi = 0` first;
  `SCSTEM_Infocrit()` passes on the partition of the `phi = 0` fit it already
  has). A penalty that is strong from the first sweep froze whatever partition
  it was given, so a start unrelated to the regimes stayed where it was. The
  automatic scale, computed at the first sweep, is then the one of the
  unpenalized fit and no longer depends on the initialization. Measured there
  it is several times larger than before, so the defaults are recalibrated:
  `phi_penalty = 0.05`, `phi_grid = c(0, 0.025, 0.05, 0.1, 0.2, 0.5, 1)` and
  the band `c(0.025, 0.2)` of `SCSTEM_Select()`; the strong values 0.5 and 1
  lie outside the band and compete only in the choice of `phi`.
* On convergence the cluster-wise models are **re-estimated once** on the final
  partition, and every reported quantity comes from that refit.

### Corrections to the SC-STEM assignment score

The score used to allocate locations to clusters before 2.0.0 contained four
defects, all fixed:

* the quadratic form of the state equation was divided by `m0` instead of
  `Sigmaeta`, which flips its sign whenever `m0 < 0`;
* `sigma2omega` and `sigma2eps` were taken from the pooled fit rather than from
  the candidate cluster, which made the corresponding term constant across
  clusters and turned the `nugget_var` switch into a no-op;
* the predictor used the lagged smoothed state `y_{t-1}` instead of the
  contemporaneous `K_i y_t` of the measurement equation;
* the score was divided by the number of locations, which left `phi` without a
  comparable scale.

The contribution is now an explicit conditional pseudo-likelihood, documented as
such, and used *only* to rank clusters: coefficients, variance components and
information criteria all come from the exact cluster-wise likelihoods.

### Information criteria

* AIC, BIC and KIC are computed on the exact total log-likelihood of the final
  refit, with `k_eff * (ncov + 3 + 3p)` free parameters and `n = d * T`
  observations. The previous implementation used the *number of clusters* as the
  number of parameters and the number of locations as the sample size.

### Convergence of the EM algorithm and computational settings

* New `STEM_control()`: the computational settings of every estimation
  function in one object, in the spirit of `optim(control = )`. All the
  estimation functions take `control`; `options(Stem.control = list(...))` sets
  them for a whole session.
* Four estimation algorithms, chosen by `STEM_control(algorithm = )`: the EM
  algorithm; the ECME algorithm (Liu and Rubin 1994), in which the regression
  coefficients and `m0` are updated on the observed likelihood by a
  generalized least-squares step computed with the Kalman filter; and both
  accelerated by SQUAREM (Varadhan and Roland 2008). They reach the same
  maximum. The default is `"SQUAREM"`: in the simulation study it fitted the
  regimes in about 40% of the time of the EM algorithm, closer to the maximum,
  and ran the whole grid of `SCSTEM_Infocrit()` in about half the time on the
  large cells (20% more on the small ones), with the same selections. The
  algorithm applies to the pooled fit and to the final refits; inside the
  alternation of SC-STEM, whose fits stop after a few iterations, the
  iterations are not accelerated. Penalized fits (`lambda > 0`) run the plain
  iterations.
* `regularization` defaults to 0 (it was 0.01) in `STEM_Estimation()`,
  `SCSTEM_Estimation()` and `STEM_Bootstrap()`. The constant was added to every
  matrix the EM algorithm inverts but not to the filter that computes the
  likelihood, so the iterations converged to a point that is not the maximum:
  the variance parameters were off by 3 to 13% on the regimes of the
  simulation study, and the variance of the latent innovations of a pooled fit
  by 84%, whatever the tolerance.
* The estimates, the log-likelihood and the smoothed states returned by
  `STEM_Estimation()` refer to the same parameters: the log-likelihood and the
  smoothed states used to be those of the last-but-one iteration.
* The EM code is split by role: the E-step (`R/estep.R`), the conditional
  updates of the M-step (`R/mstep-updates.R`), the M-steps of the EM and ECME
  algorithms (`R/mstep-em.R`, `R/mstep-ecme.R`), SQUAREM (`R/squarem.R`) and
  the wrapper that runs the algorithm chosen (`R/em-fit.R`). With
  `algorithm = "EM"` and `regularization = 0.01` the iterations reproduce those
  of the previous code to 1e-13. The smoother gain of the initial state now
  uses the transpose of `G`, which matters only for a non-symmetric `G`.
* New stopping rule of the EM algorithm, on two criteria: the largest relative
  change of a free parameter, taken one at a time (as in D-STEM v2), and the
  absolute change of the log-likelihood, with a tolerance of `1e-3` on both by
  default. With `em_stop = "any"` (the default, as in D-STEM v2) the algorithm
  stops when either criterion is met, with `em_stop = "all"` only when both are;
  in either case after `em_maxit` iterations (default 500). The previous
  rule stopped when the relative change of the whole parameter vector and the
  relative change of the log-likelihood fell below `precision`. That vector
  contained the fixed loadings, so the rule loosened with the number of
  locations, and with log-likelihoods of the order of `1e4` a relative tolerance
  stopped the algorithm while several units could still be gained.
* SC-STEM: the first time a regime is fitted it starts from the least-squares
  coefficients of its own locations instead of the pooled starting values; every
  later fit of the regime, the final refit included, starts from its previous
  estimates. Started from the pooled coefficients, a regime with a persistent
  latent process stayed near them. A penalized fit, which starts from the
  partition of the unpenalized fit, starts its regimes from the final refit of
  that partition when it is at hand (always within `SCSTEM_Infocrit()`).
* SC-STEM: the final refit runs with the same settings as the pooled fit at
  `k = 1`, so that the log-likelihoods compared across `k` are computed to the
  same accuracy; the fit reports in `em_converged` whether the EM algorithm of
  each regime met its stopping rule.
* SC-STEM: fits of the same data that end at the same partition share their
  final refit. `SCSTEM_Infocrit()` shares the refits across its grid, and a
  penalized fit that ends at the partition of the unpenalized fit it starts
  from reuses that refit. The refit depends on the partition and the settings,
  not on `phi_penalty`; two refits of one partition from different starts
  stopped at different points of the likelihood, and that difference decided
  between penalties with the same partition. The fit reports
  `refit_reused`; the argument `refit_cache` carries the shared refits and is
  not a setting to tune.
* The historical arguments (`precision`, `max.iter`, `precision_full_dataset`,
  `max_iter`, `abs_tol`, `rel_tol`) still work, as overrides of `control`.
* The limits of the Newton-Raphson step of the spatial parameters are settings
  of `control` (`nr_maxit`, `nr_hess_maxit`).

## Numerical robustness of the STEM core

These were pre-existing defects, harmless for a single pooled fit but fatal for
a clusterwise algorithm, which refits the model on hundreds of different subsets
of locations. Both aborted the fit with
`missing value where TRUE/FALSE needed`:

* `kalman()`: the Newton-Raphson step is validated before use; a singular or
  non-finite Hessian, or a non-finite step, now exits the loop keeping the last
  valid iterate instead of propagating `NaN` into the convergence test. The
  relative criterion has a floor on its denominator, the Hessian condition is
  evaluated through `isTRUE()`, and the fallback grid search is guarded against
  an underflowed scale that produced `log(0)`.
* `STEM_Estimation()`: a non-finite EM iterate is discarded with a warning and
  the last valid one is returned; both relative convergence criteria have a
  floor on the denominator and are wrapped in `isTRUE()`.

On the `pm10` example these fixes take the parametric bootstrap from 4 usable
draws out of 12 to 12 out of 12.

* `STEM_Simulation()` called `mvrnorm()` unqualified, relying on a `NAMESPACE`
  import; it now calls `MASS::mvrnorm()`.

## Larger networks

* The Kalman filter no longer forms or inverts the `d x d` predictive
  covariance. Since the measurement covariance does not depend on time and the
  state contribution has rank `p`, the Woodbury identity and the matrix
  determinant lemma reduce each step to `p x p` algebra against a Cholesky
  factor computed once per pass, so a forward pass costs `O(d^3 + T d^2)`
  instead of `O(T d^3)`. On 200 stations over 365 days an EM iteration is about
  ten times faster; the estimates are unchanged to floating point, the two forms
  being algebraic identities.
* The log-likelihood term is evaluated as a log density rather than as the
  logarithm of a density. The Gaussian density on `d` observations is of order
  `exp(-d)`, so past roughly 300 locations it underflowed to zero and the fit
  aborted with a non-finite log-likelihood. Networks of that size can now be
  fitted at all.
* The M-step accumulations are matrix products rather than loops over time, and
  no longer hold `T` matrices of size `d x d` at once. The derivative helpers
  receive the inverse of the scaled covariance from the caller instead of
  recomputing it up to seven times each, and every trace of a matrix product is
  evaluated without forming the product.
* The Newton-Raphson step that updates the two covariance parameters shares
  what it used to rebuild: the exponential kernel is evaluated once per
  iteration rather than three times, the derivative with respect to `log b` is
  carried as the multiplier of the identity that it is rather than as a dense
  matrix, the products against the inverse covariance are formed once instead
  of three times each, and the log determinant and the inverse now come from a
  single Cholesky factorisation. The `d x d` matrix products per inner
  iteration go from twelve to five.
* `mvtnorm` moves from `Imports` to `Suggests`: no function in the package uses
  it any more, though the test suite still does, as an independent reference
  implementation of the likelihood.

## CRAN compliance

Addresses the review comments received on the previous submission:

* the redundant "in R" has been removed from the package title;
* the references in `DESCRIPTION` now carry `<doi:...>`, `<https:...>` and
  ISBN links in the required format;
* `T` and `F` have been replaced by `TRUE` and `FALSE` throughout the code and
  the documentation;
* all `print()`/`cat()` diagnostics have been replaced by `message()` behind a
  `verbose = FALSE` argument, threaded through `kalman()`,
  `STEM_Estimation()`, `STEM_Bootstrap()` and the SC-STEM routines. The only
  remaining `cat()` calls are inside `print()` methods, where they belong;
* the package no longer writes to `.GlobalEnv`. Where a seed is needed for
  reproducibility, the RNG stream is saved and restored on exit.

Further packaging work required before submission:

* `pm10` moved from a 250 KB source file duplicated in `R/` and `data/` to a
  proper `data/pm10.rda` (25 KB), with `R/pm10.R` reduced to documentation;
* the dependency on **SCDA**, which is not distributed on CRAN, has been
  removed. The AMKM initialization it provided is replaced by an internal
  k-means initialization on the PCA-compressed covariate means, with multiple
  restarts, a minimum-cluster-size filter and a repair step. Any external
  initialization, AMKM included, can still be supplied through the new
  `init_partition` argument;
* `dplyr` and `sf` were likewise dropped, as nothing in the package needs them
  any more; every external call is written as `package::function()` and
  `NAMESPACE` no longer carries `importFrom` directives;
* `URL` and `BugReports` fields added, along with `LICENSE`, `NEWS.md` and a
  `testthat` suite.

## Documentation

* New conceptual map of the package under `inst/extdata/`, showing how the STEM
  and SC-STEM functions call one another, and reachable from R with
  `system.file("extdata", package = "Stem")`.
* Vignette rewritten around the full workflow, from `STEM_Model()` to the
  SC-STEM selection rule and bootstrap.


# Stem 1.0

* Original CRAN release by Michela Cameletti, subsequently archived.
