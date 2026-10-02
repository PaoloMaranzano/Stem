# Stem <img src="man/figures/logo.png" align="right" height="132" alt="" />

<!-- badges: start -->
[![License: GPL v2+](https://img.shields.io/badge/License-GPL%20(%3E%3D%202)-blue.svg)](https://www.gnu.org/licenses/gpl-2.0)
[![R >= 4.1](https://img.shields.io/badge/R-%3E%3D%204.1-blue.svg)](https://cran.r-project.org/)
<!-- badges: end -->

**Spatio-temporal expectation-maximization models, and their spatially-clustered
extension.**

`Stem` estimates hierarchical spatio-temporal models by maximum likelihood,
using Kalman filtering and smoothing inside an EM algorithm. Version 2.0.0 adds
the **spatially-clustered STEM (SC-STEM)** family, in which regression
coefficients, variance components and latent dynamics are allowed to differ
across spatial regimes that are estimated from the data rather than imposed, and
an optional **ridge, lasso or elastic-net penalty** on the regression
coefficients for the collinear designs environmental data usually produce.

## The models

### STEM: one process, three sources of variation

A pollutant is observed at $d$ locations over $T$ days. The model separates what
is measured from what is real, and splits the real part into a regression
surface, a dynamic component shared by the whole network, and small-scale
spatial noise. In matrix form, with $Z_t$ the $d$-vector of observations at
time $t$:

```math
Z_t = X_t \beta + A y_t + e_t
```

```math
y_t = G y_{t-1} + \eta_t
```

```math
\Sigma_e = \sigma^2_\varepsilon I + \sigma^2_\omega \exp(-\theta h)
```

with $e_t \sim N(0, \Sigma_e)$, $\eta_t \sim N(0, \Sigma_\eta)$ and
$y_0 \sim N(m_0, C_0)$. The three error components are mutually independent and
white over time.

| Symbol | Role | How to read it |
|---|---|---|
| $\beta$ | Regression coefficients on the covariates $X_t$ | The part of the concentration explained by observable drivers: altitude, emissions, a co-pollutant. Constant over space and time. |
| $y_t$ | Latent temporal process, dimension $p$ | The unobserved "regional level" that moves the whole network together day by day - weather, seasonality, anything not in $X_t$. |
| $A$ | $d \times p$ loading matrix | How strongly each location feels the shared process. Usually $A = 1$, meaning one common level; it can also carry EOF loadings. |
| $G$ | Transition matrix of the state equation | The **persistence** of the latent process. With $p = 1$ it is a scalar autoregressive coefficient: near 1 means long memory, near 0 means the level is renewed each day. |
| $\Sigma_\eta$ | Innovation variance of the state equation | How much genuinely new information enters the latent process at each step. |
| $\sigma^2_\varepsilon$ | Measurement error variance | Instrumental noise, independent across locations. Geostatistically it is the **nugget**: the discontinuity of the covariance at distance zero. |
| $\sigma^2_\omega$ | Variance of the small-scale spatial component | The **partial sill**: local spatial structure that the regression and the latent process do not capture. |
| $\theta$ | Range parameter of the exponential covariance | How fast spatial correlation decays with distance $h$. Its reciprocal $1/\theta$ is the characteriztic length beyond which two locations are effectively uncorrelated. |
| $\Sigma_e$ | Observation covariance | Nugget plus spatially correlated component. Its diagonal is $\sigma^2_\varepsilon + \sigma^2_\omega$, the total variance of a single observation. |
| $m_0, C_0$ | Initial state distribution | Starting values for the Kalman recursions; $C_0$ is held fixed. |

Estimation is by maximum likelihood, through Kalman filtering and smoothing
inside an EM algorithm: $\beta$, $\sigma^2_\omega$, $G$, $\Sigma_\eta$ and $m_0$
have closed-form updates, while $\theta$ and $\sigma^2_\varepsilon$ need a
Newton-Raphson step.

### SC-STEM: spatial regimes

The STEM model imposes **one** $\beta$, **one** $G$, **one** covariance on the
entire domain. That is a strong assumption. Tobler's first law of geography says
near things are more related than distant ones, and $\Sigma_e$ encodes exactly
that. But geography also obeys a **second law - spatial heterogeneity**
(Goodchild 2004): geographic variation is not uniform, and it is not only the
*values* that change across space, it is the *relationships* themselves. The
response of PM2.5 to altitude in the middle of an alluvial plain need not be the
response near the Alpine foothills.

SC-STEM addresses this by introducing **spatial regimes**. Each location $i$ is
assigned to one of $K$ latent regimes, and conditionally on belonging to regime
$k$ it follows its own STEM model:

```math
z_{it} = x_{it}^{\top} \beta_k + A_i y_t^{(k)} + e_{it}
```

```math
y_t^{(k)} = G_k y_{t-1}^{(k)} + \eta_t^{(k)}
```

```math
\Sigma_{e,k} = \sigma^2_{\varepsilon k} I + \sigma^2_{\omega k} \exp(-\theta_k h)
```

So **every regime carries a full parameter set of its own**,

```math
\Psi_k = ( \beta_k, \sigma^2_{\varepsilon k}, \sigma^2_{\omega k}, \theta_k, G_k, \Sigma_{\eta k}, m_{0k} )
```

and the model estimates $K$ such sets together with the partition itself. Not
only the regression coefficients vary: so do the nugget, the spatial range and
the persistence of the latent dynamics, so a regime can differ from another in
*how* it behaves in time and space, not merely in level.

The partition is not imposed. Labels $k_1, \ldots, k_d$ and parameters are
estimated jointly by maximizing a Potts-penalized log-likelihood,

```math
Q = \sum_{i=1}^{d} \ell_i(k_i) + \phi \, c \sum_{(i,j) \in E} I(k_i = k_j)
```

where $\ell_i(k)$ is the log-likelihood contribution of location $i$ under the
parameters of regime $k$, $E$ is the edge set of the symmetrized
$k$-nearest-neighbor graph with each unordered pair counted once,
$I(\cdot)$ the indicator function, and $c$ a scale factor documented in
`?SCSTEM_Estimation`.

The penalty is what makes the regimes *spatial*. The first term rewards fit and
would happily scatter the labels; the second rewards neighboring locations
sharing a label. The hyperparameter $\phi \ge 0$ arbitrates between them:
$\phi = 0$ gives ordinary clusterwise regression with no spatial structure,
large $\phi$ gives contiguous and rigid regimes. Setting $K = 1$ returns the
pooled model, which stays available as the reference against which any clustered
fit must justify itself.

### Ridge, lasso and elastic net on the coefficients

Environmental covariates are collinear by construction: meteorological drivers
are measured on overlapping supports and chemical species share sources. The
variance of $\hat\beta$ is governed by the inverse of

```math
M = \sum_{t} X_t^{\top} \Sigma_e^{-1} X_t
```

and it is that inverse a penalty stabilizes. The argument bites hardest in
SC-STEM, where each regime fits a complete model on a *subset* of the network:
collinearity that a 200-station network absorbs, a regime of twelve does not,
and the smallest regime that can still be fitted is in practice what caps the
number of regimes worth entertaining.

`STEM_Fit()` adds an elastic net on the regression coefficients, in the
parameterization of `glmnet`:

```math
\ell_\pi(\psi) = \ell(\psi) - \lambda \left\{ \alpha \lVert D\beta \rVert_1 + \frac{1-\alpha}{2} \beta^{\top} D \beta \right\}
```

with $D$ the diagonal indicator of the penalized coordinates. So $\alpha = 0$ is
**ridge**, $\alpha = 1$ the **lasso**, and anything in between the **elastic
net**. Only $\beta$ is penalized: the variance components, the range, the
transition matrix and the initial state keep their maximum likelihood values,
because they describe the error process rather than the mean.

**The E-step does not change.** The penalty is a function of $\beta$ alone and
does not involve the latent states, so it passes through the conditional
expectation untouched and the usual argument still gives an EM algorithm that
increases the *penalized* likelihood at every iteration. The Kalman filter and
smoother -- the expensive part -- are not touched at all; only the point
returned by the M-step differs. Under a ridge that point keeps a closed form,

```math
\hat\beta = (M + \lambda D)^{-1} v , \qquad v = \sum_t X_t^{\top} \Sigma_e^{-1} v_t
```

one added diagonal on a matrix already assembled; under a lasso or an elastic
net it is found by cyclic coordinate descent, which converges to the *exact*
maximizer because the objective is convex and the penalty separable.

Three points are specific to this model and worth knowing before choosing
$\lambda$.

- **The metric is not Euclidean.** $M$ is a *generalized* least squares
  cross-product: it carries an estimated $\Sigma_e^{-1}$. Standardizing the
  columns of $X$ the usual way is therefore not what makes $\lambda$ comparable
  across covariates here. The scaling that does is the one putting the diagonal
  of $M$ at one, and the package applies it internally, returning coefficients
  on the original scale. In SC-STEM this is also what makes a single $\lambda$
  mean the same thing in regimes of different size.
- **Two scales have to be removed, and only one of them from each half.** After
  the scaling above, the $L_2$ part is *already* free of the units of the
  response: it multiplies a curvature whose diagonal is one, so it shrinks by
  $1/(1+\lambda)$ whatever the units. The $L_1$ part is not, because a threshold
  has to be compared with the gradient. `lambda_scale = "relative"`, the
  default, therefore rescales the $L_1$ part only, by the largest partial
  gradient: $\lambda$ is then unit-free for every $\alpha$, and every penalized
  coefficient is exactly zero once $\lambda\alpha \ge 1$, so $\lambda \in (0,1]$
  traverses the whole lasso path as in `glmnet`. Rescaling both halves would
  break the ridge to fix the lasso.
- **$\alpha$ and $\lambda$ are pooled**, one pair for the whole partition, and
  under the scalings above that already means the same *proportional* shrinkage
  in every regime: the reference is computed inside each regime, so a regime
  with half the locations has a proportionally smaller reference and the same
  $\lambda$ buys the same fraction of the path. `lambda_by = "size"` departs
  from that on purpose, setting $\lambda_g = \lambda \bar n / n_g$ so that a
  regime of half the average size is shrunk twice as hard. Either way $\lambda$
  stays **one** hyperparameter; genuinely cluster-specific
  $(\alpha_g, \lambda_g)$ is a different model and is not offered here.
- **The intercept is not penalized**, and not merely by convention: the model
  already carries a latent process whose initial mean $m_0$ absorbs the level,
  so shrinking $\beta_0$ would not shrink "the level" but move it into $m_0$ at
  a rate depending on $G$.
- **The criteria count effective parameters.** With a penalty in force the
  nominal $r$ overstates the flexibility of the fit, so AIC, BIC and KIC use
  $\mathrm{tr}( M (M + \lambda D)^{-1} )$ for a ridge, the number of active
  coefficients for a lasso, and the corresponding trace on the active set for an
  elastic net. At $\lambda = 0$ this reduces to $r$ and nothing changes.

`STEM_Fit()` is the single entry point for all four estimators, and which one
runs is decided by two arguments and nothing else:

| `K` | `lambda` | what is fitted |
|---|---|---|
| `1` | `0` | the pooled STEM model |
| `> 1` | `0` | the spatially-clustered model |
| `1` | `> 0` | the pooled model with an elastic net on $\beta$ |
| `> 1` | `> 0` | the clustered model, penalized within each regime |

The defaults `K = 1`, `alpha = 0`, `lambda = 0` reproduce `STEM_Estimation()`
bit for bit, which the test suite asserts on the whole parameter vector and on
the log-likelihood.

```r
fit0 <- STEM_Fit(mod)                                  # classical STEM
fitr <- STEM_Fit(mod, alpha = 0,   lambda = 0.3)       # ridge
fitl <- STEM_Fit(mod, alpha = 1,   lambda = 0.3)       # lasso
fite <- STEM_Fit(mod, alpha = 0.5, lambda = 0.3)       # elastic net
fitc <- STEM_Fit(mod, K = 3, phi_penalty = 0.05,       # SC-STEM, ridge per regime
                 alpha = 0, lambda = 0.3)
```

**Penalized linear regression is the degenerate case.** Two switches take the
model down to it: `latent = FALSE` sets the loading matrix to zero, so the state
contributes nothing and the parameters describing it, now unidentified, are held
where they started; `spatial = FALSE` replaces the exponential correlation by
the identity, so $\Sigma_e = \sigma^2 I$ and the Newton-Raphson step is skipped.
What is left is $z_{ti} = x_{ti}'\beta + e_{ti}$ with $e \sim N(0, \sigma^2 I)$,
estimated by penalized least squares.

```r
fit <- STEM_Fit(mod, alpha = 0.5, lambda = 0.3,
                latent = FALSE, spatial = FALSE, regularization = 0)
```

Set `regularization = 0` as well: the small ridge the package adds for
conditioning is otherwise the only thing separating the two, and on a collinear
design it is not negligible. With that, the agreement with an elastic net
computed directly on $X'X$ and $X'y$ is between $10^{-14}$ and $10^{-11}$ for
the ridge, the lasso, the elastic net and the unpenalized case alike. With
`K > 1` the same switches give clusterwise penalized regression, the partition
still estimated.

$\lambda$ and $\alpha$ are hyperparameters like $K$ and $\phi$ and are chosen
the same way: by an information criterion computed with the effective degrees of
freedom, or -- the honest route -- by spatio-temporal cross-validation with
blocking that respects both dependencies (Otto, Fasso and Maranzano 2024).
`SCSTEM_Infocrit()` accepts `alpha` and `lambda` and returns criteria already
corrected for the effective degrees of freedom, so a grid over all four
hyperparameters is a loop over $(\alpha, \lambda)$ around the $(K, \phi)$ grid
it already traverses; the two-step rule of `SCSTEM_Select()` still arbitrates
$(K, \phi)$ only. One caveat: if $\lambda$ is selected from the data, the parametric
bootstrap must repeat the selection on every draw, exactly as
`SCSTEM_Bootstrap()` repeats the clustering; and with $\alpha > 0$ the estimator
is not smooth, so intervals for a coefficient at the boundary do not have their
usual coverage interpretation.

The derivation -- the penalized EM and its monotonicity, both forms of the
M-step, the GLS metric, the degrees of freedom and what is still open -- is a
standalone document at
[`dev/regularization/stem-elastic-net.tex`](dev/regularization/stem-elastic-net.tex).

### Statistical features and scope

What the model class covers, and what it does not. Every entry below reflects
the current implementation, not the model on paper.

| Feature | Supported | Notes |
|---|---|---|
| Distribution of the response | Gaussian only | The EM closed forms and the Kalman recursions both rely on normality. Non-Gaussian responses have to be transformed first. |
| Response dimension | univariate | `z` is a $T \times d$ matrix of **one** variable measured at $d$ sites. Several pollutants modeled jointly is a different specification. |
| Latent state | multivariate | $p \ge 1$ latent processes, default $p = 1$. This is what is multivariate in the model. |
| Loading matrix `A` | known, user-supplied | $d \times p$, not estimated, and common across regimes in SC-STEM. |
| Latent dynamics | VAR(1) | `G` and `Sigmaeta` are diagonal by default; both can be made full via `flag.Gdiag` and `flag.Sigmaetadiag`. |
| Initial condition | `m0` estimated, `C0` fixed | |
| Spatial covariance | exponential only | $\sigma^2_\epsilon I + \sigma^2_\omega \exp(-\theta h)$: isotropic and stationary *within* a regime. Only this function has the analytical derivatives the Newton-Raphson step needs. |
| Non-stationarity in space | across regimes only | Each regime carries its own $\theta_k$, $\sigma^2_{\epsilon k}$, $\sigma^2_{\omega k}$, so the partition is itself a coarse form of non-stationarity. |
| Spatial support | point-referenced | Distances `"euclidean"` or `"geo"`. Areal data only through centroids, and only when the units are small relative to the distances between them. |
| Time | discrete, regularly spaced | The state equation links consecutive time points; irregular spacing would need a continuous-time formulation. |
| Panel | may be unbalanced | Gaps in the response are allowed; see the note below. |
| **Missing values in `z`** | **supported** | Handled by the EM algorithm following Durbin and Koopman (2012), Sect. 2.7 and 4.10. Every location must keep at least one observation. |
| **Missing values in covariates or coordinates** | **not supported** | The design matrix enters the closed-form M-step and the coordinates the distance matrix. Impute before fitting. |
| Network composition | fixed over time | A station entering or leaving the window is represented by marking the unobserved periods as `NA` in the response. |
| Regime sizes (SC-STEM only) | bounded below | Each regime needs enough locations for its own fit, which caps the number of regimes that can be entertained on a given network. |
| Regularization of `beta` | ridge, lasso, elastic net | Through `STEM_Fit(alpha, lambda)`, within each regime in SC-STEM. Only the regression coefficients are penalized; the variance components and the range stay at their maximum likelihood values. `lambda = 0`, the default, reproduces the classical estimator exactly. |

> **On missing values in the response.** The treatment follows Durbin and
> Koopman (2012, 2nd ed.), Sections 2.7 and 4.10. At each time point the
> measurement equation is restricted to the locations actually observed, through
> a selection matrix whose rows are a subset of the rows of the identity; a time
> point at which nothing is observed contributes no update and no likelihood
> term, which is their `Z_t = 0` device. The backward smoothing recursions need
> no change at all, since they read the filtered moments and the transition,
> never the data.
>
> The M-step is where the work is. EM maximizes the *expected complete-data*
> log-likelihood, so the sufficient statistics are completed rather than
> truncated: a missing value enters through its conditional expectation given
> everything observed, and its conditional variance is added back as a
> correction. That conditional expectation is *not* the signal alone -- which
> would be exact only for a diagonal covariance -- because the spatial covariance
> couples the locations, so the missing block of the measurement error is
> predicted from the observed one by the same algebra as kriging at a fixed time
> point. The divisor of the variance update stays the complete-data count, for
> the same reason. Kriging and both bootstraps follow: the predictor conditions
> on the observed sub-vector, and a bootstrap replicate reproduces the observed
> design, gaps included.
>
> On complete data every result is bit-for-bit what it was before the feature
> was added; this is asserted by the test suite.

## Installation

```r
# install.packages("remotes")
remotes::install_github("PaoloMaranzano/Stem")
```

## What is in the package

Two functions are all a first session needs: `STEM_Model()` to build the object,
`STEM_Fit()` to estimate it. Everything else is either a layer above (selection,
uncertainty, prediction) or an engine below.

**Build and fit**

| Function | Purpose |
|---|---|
| `STEM_Model()` | build the model object from data, coordinates and starting values |
| `STEM_Fit()` | **the entry point.** Data and hyperparameters in, one fitted model out. `K` chooses pooled or clustered, `(alpha, lambda)` unpenalized or penalized |

**Choosing the hyperparameters**

| Function | Purpose |
|---|---|
| `SCSTEM_Infocrit()` | fit a grid of `(K, phi)` and return the criteria, with the effective degrees of freedom when a penalty is in force |
| `SCSTEM_Select()` | the two-step rule: `K` by the modal criterion in a band, `phi` by the smallest criterion at that `K` |

**After the fit**

| Function | Purpose |
|---|---|
| `STEM_Signal()`, `SCSTEM_Signal()` | **fitted values**: the conditional mean the model estimates, regression plus latent process, without the measurement error |
| `STEM_Complete()`, `SCSTEM_Complete()` | gap-filling: `E[z given observed]`, which returns the data wherever the data exist |
| `STEM_Kriging()` | spatial prediction at unobserved locations |
| `STEM_Simulation()` | simulate from a fitted or specified model |
| `STEM_Bootstrap()` | parametric bootstrap for the pooled model |
| `SCSTEM_Bootstrap()` | refit-with-clustering bootstrap: the partition is re-estimated on every draw |
| `SCSTEM_BootInference()` | standard errors, intervals and between-regime tests |
| `SCSTEM_CV()` | optional validation: blocked cross-validation (locations, times, both, random cells) of one or more fits on the same folds, ranked by predictive error; it does not choose `(k, phi)` |

**The two engines**, which `STEM_Fit()` dispatches to. Call them directly only
if you want to bypass the dispatch; the arguments and the return values are the
same either way, and the historical API is preserved.

| Function | Purpose |
|---|---|
| `STEM_Estimation()` | the pooled fit: EM with Kalman filtering and smoothing |
| `SCSTEM_Estimation()` | the clustered fit: `STEM_Estimation()` per regime, alternated with an ICM sweep on the labels |

Two datasets ship with the package: `pm10` (22 stations, 366 days, the original
example) and `povalley` (36 background stations of the Po Valley, daily PM2.5
over 2019--2023, with altitude and PM10 as covariates).

## A short example

```r
library(Stem)
data(povalley)

Tn <- 365L
Tfull <- nrow(povalley[["z"]]); d <- ncol(povalley[["z"]])
keep <- as.vector(outer(seq_len(Tn), (seq_len(d) - 1L) * Tfull, "+"))

phi <- list(beta = matrix(c(1.25, -0.00003, 0.64), 3, 1),
            sigma2eps = 18.66, sigma2omega = 1e-06, theta = 2e-06,
            G = matrix(0.59, 1, 1), Sigmaeta = matrix(4.25, 1, 1),
            m0 = as.matrix(0), C0 = as.matrix(1))

mod <- STEM_Model(z = povalley[["z"]][seq_len(Tn), ],
                  covariates = povalley[["covariates"]][keep, ],
                  coordinates = povalley[["coords"]],
                  phi = phi, A = matrix(1, d, 1))

# one entry point for every estimator: K chooses pooled or clustered,
# lambda chooses penalized or not
fit <- STEM_Fit(mod, distance = "geo")                 # classical STEM

# explore the grid and let the two-step rule choose K and phi
ic  <- SCSTEM_Infocrit(mod, K_grid = 1:4, distance = "geo")
sel <- SCSTEM_Select(ic)
sel

# the same clustered model with a ridge on the coefficients of every regime
fitr <- STEM_Fit(mod, K = sel[["K_selected"]], phi_penalty = sel[["phi_selected"]],
                 alpha = 0, lambda = 2, distance = "geo")

# uncertainty, with the partition re-estimated at every draw
boot <- SCSTEM_Bootstrap(sel[["fit"]], B = 200, seed = 1)
SCSTEM_BootInference(boot)
```

## Design notes

Six points make the implementation behave sensibly in practice. The first five
concern SC-STEM and are documented in detail in `?SCSTEM_Estimation`; the sixth is
what makes a large network fittable at all.

**An unpenalized fit starts from the departures from the pooled model.** The
default `init_method = "departures"` fits the pooled model and clusters the
locations on how they depart from it: the mean of the residual, its slopes on
the covariates, and the autocorrelation and variance of what is left. Locations
of one regime share their departures; the covariate means, the earlier default
and still an option, say nothing about the regimes when the covariates are
exogenous to them.

**A penalized fit starts from the unpenalized one.** A fit with $\phi > 0$
starts from the partition of the fit with $\phi = 0$ at the same $K$, and the
automatic scale of the penalty is measured there, on regimes that are already
fitted. A penalty that is strong from the first sweep would freeze whatever
partition it is given, and a scale measured on the starting partition would
make the same $\phi$ mean different things for different initializations. On a
grid the rule costs nothing: `SCSTEM_Infocrit()` passes the partition of its
$\phi = 0$ fit on. Measured on fitted regimes the scale is of the order of the
whole gain of the right regime over the wrong ones, so the useful values of
$\phi$ are small: the default grid is
$\phi \in \{0, 0.025, 0.05, 0.1, 0.2, 0.5, 1\}$. `SCSTEM_Select()` chooses $K$
by the modal BIC over the moderate values $[0.025, 0.2]$, the pooled model
competing, and then $\phi$ by the smallest BIC at that $K$ over the whole grid:
since every penalized fit starts from the unpenalized one, a penalty is chosen
only when it leads to a partition with a higher likelihood, and the strong
values $0.5$ and $1$ are kept on the grid for data on which they do.

**Labels are updated sequentially (ICM).** Each location maximizes its own
penalized contribution given the current labels of all the others, so a sweep
cannot decrease the objective at fixed parameters and cannot cycle - unlike a
simultaneous update, which remains available as an option.

**Degeneracy is controlled.** Within-cluster homogeneity is exactly what the
assignment step pursues, so the regime-specific variances shrink and, left
unconstrained, the regime with the smallest residual variance attracts every
location. A minimum-size constraint prevents the collapse, and a size-preserving
swap pass restores mobility without breaking feasibility.

**Uncertainty includes the partition.** `SCSTEM_Bootstrap()` re-runs the whole
procedure, clustering included, on every draw, so the reported intervals are not
conditional on a partition that is itself estimated.

**The filter never inverts a $d \times d$ matrix.** The predictive covariance
$Q_t = A P_t A' + \Sigma_\varepsilon$ has $\Sigma_\varepsilon$ constant in $t$
and a rank-$p$ update on top of it, so the Woodbury identity and the matrix
determinant lemma reduce each step to $p \times p$ algebra against a Cholesky
factor of $\Sigma_\varepsilon$ that is computed once per pass. The cost of a
forward pass is $O(d^3 + T d^2)$ instead of $O(T d^3)$, and the log density is
evaluated in closed form rather than as the logarithm of a density, which on a
few hundred locations underflows to `-Inf`. On a network of 200 stations over
365 days one EM iteration is about ten times faster than the direct form, and a
400-station network, which the direct form cannot fit at all, takes about a
second per iteration. The two forms are algebraic identities of each other and
agree to about $10^{-13}$ in relative terms.

## Documentation

```r
vignette("getting-started", package = "Stem")      # the classical STEM workflow
vignette("SCSTEM", package = "Stem")               # spatially-clustered STEM
vignette("function-map", package = "Stem")         # map of the package
vignette("computational-notes", package = "Stem")  # how the engine is computed
```

The last one explains why a forward pass costs $O(d^3 + T d^2)$ rather than
$O(T d^3)$, what that means, and under which assumptions -- the Woodbury
identity, the determinant lemma, the log-density, and the matrix form of the
M-step. The theory behind the regularized estimator is kept outside the package,
at [`dev/regularization/stem-elastic-net.tex`](dev/regularization/stem-elastic-net.tex).

A reading guide mapping each equation of the papers onto the function that
implements it ships with the package as a standalone document:

```r
browseURL(system.file("extdata", "STEM_model_verification.html", package = "Stem"))
```

The development record - the defects fixed in this release and the reasoning
behind the algorithmic changes - is kept in the repository, outside the built
package, at [`dev/code-changes-report.html`](dev/code-changes-report.html). The
user-facing summary is [`NEWS.md`](NEWS.md), and the dated development log is
[`CHANGELOG.md`](CHANGELOG.md).

## Theoretical references

The STEM model implemented here comes from three companion works by Fasso and
Cameletti, which play different roles.

**Fasso, A., Cameletti, M. and Nicolis, O. (2007).** Air quality monitoring
using heterogeneous networks. *Environmetrics*, 18, 245-264.
<https://doi.org/10.1002/env.837>
> Introduces the *geostatistical dynamical calibration* (GDC) model, the most
> general of the three: it adds instrument calibration components - an additive
> bias `A(t)` and a multiplicative bias `B(t)` - so that a network of
> heterogeneous instruments (gravimetric and TEOM monitors) can be modeled
> jointly, with the loading matrix obtained by empirical orthogonal functions.
> The model estimated by this package is the special case with no calibration
> bias.

**Fasso, A. and Cameletti, M. (2007).** A general spatio-temporal model for
environmental data. *GRASPA Technical Report* n. 27.
> The direct theoretical reference for this package, which it announces by name.
> It states the three-stage hierarchy, the scaled spatial covariance
> `Gamma(h) = 1 + gamma` at `h = 0` and `C_theta(h)` otherwise with
> `gamma = sigma2eps/sigma2omega`, the choice of estimating `log(gamma)` rather
> than `sigma2eps` for positive-definiteness, the EM algorithm with closed-form
> M-steps for `beta`, `sigma2omega`, `G`, `Sigma_eta` and `m0` and
> Newton-Raphson for the spatial covariance parameters, and the spatio-temporal
> parametric bootstrap.

**Fasso, A. and Cameletti, M. (2010).** A unified statistical approach for
simulation, modeling, analysis and mapping of environmental data. *Simulation*,
86, 139-153. <https://doi.org/10.1177/0037549709102150>
> The most complete published statement: the same model and EM algorithm, plus
> the kriging predictor used by `STEM_Kriging()`, the bootstrap of
> `STEM_Bootstrap()`, and a sensitivity analysis of the model components. Its
> Equations (12)-(18) are what the estimation code implements.

### Spatial heterogeneity

**Goodchild, M. F. (2004).** The validity and usefulness of laws in geographic
information science and geography. *Annals of the Association of American
Geographers*, 94(2), 300-303.
<https://doi.org/10.1111/j.1467-8306.2004.09402008.x>
> Articulates spatial heterogeneity as a second law of geography, alongside
> Tobler's first law on spatial dependence. The motivation for letting the
> regression relationship itself vary across space.

**Zhu, A.-X. and Turner, M. (2022).** How is the Third Law of Geography
different? *Annals of GIS*, 28(1), 57-67.
<https://doi.org/10.1080/19475683.2022.2026467>
> Situates spatial dependence, spatial heterogeneity and geographic similarity
> with respect to one another, and clarifies what each principle does and does
> not claim.

### The spatially-clustered apparatus

These works supply the algorithmic devices that SC-STEM adapts; the statistical
model remains the STEM one above.

**Besag, J. (1986).** On the statistical analysis of dirty pictures.
*JRSS-B*, 48, 259-302.
> The Iterated Conditional Modes algorithm used for the label update.

**Sugasawa, S. and Murakami, D. (2021).** Spatially clustered regression.
*Spatial Statistics*, 44, 100525.
<https://doi.org/10.1016/j.spasta.2021.100525>
> The Potts-type spatial penalty on the partition.

**Cerqueti, R., Maranzano, P. and Mattera, R. (2025).** Spatially-clustered
spatial autoregressive models with application to agricultural market
concentration in Europe. *JABES*.
<https://doi.org/10.1007/s13253-025-00685-7>
> The same penalty carried over to spatial econometric models, and the
> hyperparameter-selection practice this package follows.

**Maranzano, P., Mattera, R. and Sugasawa, S. (2026).** Small area estimation
under spatial regimes: spatially clustered Fay-Herriot models for agricultural
indicators. *arXiv:2608.13638*. <https://arxiv.org/abs/2608.13638>
> A different class of model - area-level small area estimation with known
> sampling variances - in which the same apparatus is developed: the ICM label
> update, the information criteria on the final refit, the two-step rule for the
> number of regimes and the penalty, and the refit-with-clustering parametric
> bootstrap. SC-STEM adapts those devices to the point-referenced
> spatio-temporal setting.

### Regularization

**Zou, H. and Hastie, T. (2005).** Regularization and variable selection via the
elastic net. *JRSS-B*, 67, 301-320.
<https://doi.org/10.1111/j.1467-9868.2005.00503.x>
> The `(alpha, lambda)` parameterization `STEM_Fit()` follows, and the argument
> for the mixed penalty when the covariates are correlated in groups - which is
> what environmental drivers are.

**Zou, H., Hastie, T. and Tibshirani, R. (2007).** On the degrees of freedom of
the lasso. *The Annals of Statistics*, 35, 2173-2192.
<https://doi.org/10.1214/009053607000000127>
> Why the number of active coefficients is the right count to put into an
> information criterion, which is what the package uses once a penalty is in
> force.

**Friedman, J., Hastie, T. and Tibshirani, R. (2010).** Regularization paths for
generalized linear models via coordinate descent. *Journal of Statistical
Software*, 33, 1-22. <https://doi.org/10.18637/jss.v033.i01>
> The coordinate-descent scheme the M-step uses when `alpha > 0`. The objective
> here is the expected complete-data log-likelihood rather than a residual sum
> of squares, but it is a quadratic with a separable penalty and the update is
> the same soft-thresholding.

**Otto, P., Fasso, A. and Maranzano, P. (2024).** A review of regularised
estimation methods and cross-validation in spatiotemporal statistics.
*Statistics Surveys*, 18, 299-340. <https://doi.org/10.1214/24-SS150>
> Regularization and cross-validation in exactly this setting, including why
> random K-fold is anti-conservative when the residuals are correlated in space
> and time, and which blocking schemes to use instead when selecting `lambda`.

## Authors

| Author | Role | ORCID | GitHub |
|---|---|---|---|
| Michela Cameletti | aut, cph | [0000-0002-6502-7779](https://orcid.org/0000-0002-6502-7779) | [@michelacameletti](https://github.com/michelacameletti) |
| Francesco Caccia | aut, cph | | |
| Paolo Maranzano | aut, cre, cph | [0000-0002-9228-2759](https://orcid.org/0000-0002-9228-2759) | [@PaoloMaranzano](https://github.com/PaoloMaranzano) |

Licensed under GPL (>= 2).
