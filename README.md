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

**STEM** separates what is measured from what is real, and splits the real part
into a regression surface, a dynamic component shared by the whole network, and
small-scale spatial noise. With $z_t$ the $d$-vector of observations at time
$t$:

```math
z_t = X_t \beta + A y_t + e_t, \qquad y_t = G y_{t-1} + \eta_t, \qquad
\Sigma_e = \sigma^2_\varepsilon I + \sigma^2_\omega \exp(-\theta h),
```

with $e_t \sim N(0, \Sigma_e)$, $\eta_t \sim N(0, \Sigma_\eta)$ and
$y_0 \sim N(m_0, C_0)$: a linear Gaussian state-space model.

**SC-STEM** assigns each location to one of $K$ spatial regimes, estimated from
the data, and gives every regime its own STEM model and its own full parameter
set $\Psi_k = (\beta_k, \sigma^2_{\varepsilon k}, \sigma^2_{\omega k}, \theta_k, G_k, \Sigma_{\eta k}, m_{0k})$.
Labels and parameters are estimated together by maximizing a Potts-penalized
log-likelihood,

```math
Q = \sum_{i=1}^{d} \ell_i(k_i) + \phi \, c \sum_{(i,j) \in E} I(k_i = k_j),
```

whose penalty rewards neighbouring locations sharing a regime; $K = 1$ is the
pooled model. An optional **ridge, lasso or elastic net** acts on the regression
coefficients of every regime.

Details: [the models](docs/models.md), [the penalty on the coefficients](docs/regularization.md).

## Estimation algorithms

Every model, the pooled one and each regime of SC-STEM, is estimated by an
EM-type algorithm on the state-space form: the **Kalman filter** gives the
log-likelihood, the **Kalman smoother** the conditional moments of the latent
states (the E-step), and the M-step updates $\beta$, $\sigma^2_\omega$, $G$,
$\Sigma_\eta$ and $m_0$ in closed form and the spatial parameters by
Newton-Raphson. Four algorithms are available through
`STEM_control(algorithm = )`; they reach the same maximum and differ in speed
and in how close to it a given stopping rule leaves them.

| `algorithm` | what it does | time | distance from the maximum at the stop |
|---|---|---|---|
| `"EM"` | the EM algorithm | 1 | 0.020 |
| `"ECME"` | as EM, with $\beta$ and $m_0$ updated on the observed likelihood (a GLS step computed by the Kalman filter) | 1.17 | 0.017 |
| **`"SQUAREM"`** (default) | EM iterations extrapolated along the direction of slowest convergence, with a monotonicity safeguard | **0.53** | 0.0027 |
| `"SQUAREM-ECME"` | the same on the ECME iterations | 0.62 | 0.0023 |

Time relative to EM on whole SC-STEM grids of the simulation study; distance:
median shortfall of the log-likelihood on single fits, at the default
tolerances. On the same data SQUAREM selected the same models as EM, and its
estimates matched the exact maximum against the true values.

The iterations stop when the largest relative change of a parameter or the
absolute change of the log-likelihood falls below `1e-3`. D-STEM uses `1e-4`
on the relative change of either: stricter on the parameters, but its criterion
on the log-likelihood amounts to 0.4 to 1.5 units on typical fits, so it stops
several units short of the maximum, where `Stem` stops about 0.003 units short.
No ridge is added to the matrices the iterations invert
(`regularization = 0`): with one, they converged to a point that is not the
maximum.

### The range at the boundary

On the few locations of a small regime the likelihood is flat in the range
beyond two limits: a range shorter than the distance between the two closest
locations, where the spatial field cannot be told from the nugget, and one much
longer than the largest distance, where it cannot be told from an effect common
to the locations. `Stem` keeps

```math
\frac{-\log 0.95}{h_{\max}} \;\le\; \theta \;\le\; \frac{-\log 0.05}{h_{\min}},
```

accepts a Newton step only if it does not worsen the objective, and reports a
range that ends at a limit (`theta_bound`). Where nothing binds the estimates
are unchanged; where something does, the fits were better, and the boundary
cases that used to dominate the computing time are gone.

Details, with formulas and the numerical experiments:
[computational aspects](docs/computational-aspects.md) and the
Computational Supplement that ships with the package (see
[Documentation](#documentation)).

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
| `SCSTEM_CV()` | optional validation: blocked cross-validation (locations, times, both, random cells) of one or more fits on the same folds, ranked by predictive error; it does not choose `(K, phi)` |

**The two engines**, which `STEM_Fit()` dispatches to. Call them directly only
if you want to bypass the dispatch; the arguments and the return values are the
same either way, and the historical API is preserved.

| Function | Purpose |
|---|---|
| `STEM_Estimation()` | the pooled fit: Kalman filter and smoother inside an EM-type algorithm (SQUAREM by default) |
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

## Documentation

Topic pages, in [`docs/`](docs/README.md):

| page | what it covers |
|---|---|
| [The models](docs/models.md) | STEM and SC-STEM, their parameters, and what the model class covers and does not |
| [Computational aspects](docs/computational-aspects.md) | the Kalman filter and smoother, the EM, ECME and SQUAREM algorithms, the stopping rule, the limits of the range, the numerical experiments |
| [Ridge, lasso and elastic net](docs/regularization.md) | the penalty on the regression coefficients |
| [Design notes](docs/sc-stem-design.md) | the implementation choices of SC-STEM |
| [How the functions fit together](docs/function-map.md) | the map of the package |
| [Theoretical references](docs/references.md) | the works behind the package |

Shipped with the package:

```r
vignette("getting-started", package = "Stem")      # the classical STEM workflow
vignette("SCSTEM", package = "Stem")               # spatially-clustered STEM
vignette("function-map", package = "Stem")         # map of the package
vignette("computational-notes", package = "Stem")  # how the filter and the M-step are computed

# the Computational Supplement: estimation algorithms, derivations, experiments
browseURL(system.file("extdata", "Stem-computational-supplement.pdf", package = "Stem"))
# a reading guide from each equation of the papers to the function that implements it
browseURL(system.file("extdata", "STEM_model_verification.html", package = "Stem"))
```

The development record - the defects fixed in this release and the reasoning
behind the algorithmic changes - is kept in the repository, outside the built
package, at [`dev/code-changes-report.html`](dev/code-changes-report.html). The
user-facing summary is [`NEWS.md`](NEWS.md), and the dated development log is
[`CHANGELOG.md`](CHANGELOG.md). The theory of the regularized estimator is at
[`dev/regularization/stem-elastic-net.tex`](dev/regularization/stem-elastic-net.tex).

## Authors

| Author | Role | ORCID | GitHub |
|---|---|---|---|
| Michela Cameletti | aut, cph | [0000-0002-6502-7779](https://orcid.org/0000-0002-6502-7779) | [@michelacameletti](https://github.com/michelacameletti) |
| Francesco Caccia | aut, cph | | |
| Paolo Maranzano | aut, cre, cph | [0000-0002-9228-2759](https://orcid.org/0000-0002-9228-2759) | [@PaoloMaranzano](https://github.com/PaoloMaranzano) |

Licensed under GPL (>= 2).
