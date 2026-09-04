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
across spatial regimes that are estimated from the data rather than imposed.

## The models

### STEM: one process, three sources of variation

A pollutant is observed at $d$ locations over $T$ days. The model separates what
is measured from what is real, and splits the real part into a regression
surface, a dynamic component shared by the whole network, and small-scale
spatial noise. In matrix form, with $Z_t$ the $d$-vector of observations at
time $t$:

```math
Z_t = X_t \beta + K y_t + e_t
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
| $K$ | $d \times p$ loading matrix | How strongly each location feels the shared process. Usually $K = 1$, meaning one common level; it can also carry EOF loadings. |
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
assigned to one of $k$ latent regimes, and conditionally on belonging to regime
$k$ it follows its own STEM model:

```math
z_{it} = x_{it}^{\top} \beta_k + K_i y_t^{(k)} + e_{it}
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

and the model estimates $k$ such sets together with the partition itself. Not
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
`?SCSTEM_Estim`.

The penalty is what makes the regimes *spatial*. The first term rewards fit and
would happily scatter the labels; the second rewards neighboring locations
sharing a label. The hyperparameter $\phi \ge 0$ arbitrates between them:
$\phi = 0$ gives ordinary clusterwise regression with no spatial structure,
large $\phi$ gives contiguous and rigid regimes. Setting $k = 1$ returns the
pooled model, which stays available as the reference against which any clustered
fit must justify itself.

## Installation

```r
# install.packages("remotes")
remotes::install_github("PaoloMaranzano/Stem")
```

## What is in the package

| Function | Purpose |
|---|---|
| `STEM_Model()` | build the model object from data and starting values |
| `STEM_Estimation()` | maximum likelihood via EM and the Kalman filter |
| `STEM_Simulation()` | simulate from a fitted or specified model |
| `STEM_Kriging()` | spatial prediction at unobserved locations |
| `STEM_Bootstrap()` | parametric bootstrap for the pooled model |
| `SCSTEM_Estim()` | fit a spatially-clustered STEM model |
| `SCSTEM_Infocrit()` | information criteria over a grid of `(k, phi)` |
| `SCSTEM_Select()` | two-step selection of the hyperparameters |
| `SCSTEM_Bootstrap()` | refit-with-clustering parametric bootstrap |
| `SCSTEM_BootInference()` | standard errors, intervals and between-regime tests |

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
                  phi = phi, K = matrix(1, d, 1))

# explore the grid and let the two-step rule choose k and phi
ic  <- SCSTEM_Infocrit(mod, k_grid = 1:4, phi_grid = seq(0, 1, by = 0.25),
                       distance = "geo")
sel <- SCSTEM_Select(ic, band = c(0.25, 1))
sel

# uncertainty, with the partition re-estimated at every draw
boot <- SCSTEM_Bootstrap(sel[["fit"]], B = 200, seed = 1)
SCSTEM_BootInference(boot)
```

## Design notes

Three points make SC-STEM behave sensibly in practice, and are documented in
detail in `?SCSTEM_Estim`.

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

## Documentation

```r
vignette("getting-started", package = "Stem") # the classical STEM workflow
vignette("SCSTEM", package = "Stem")          # spatially-clustered STEM
vignette("function-map", package = "Stem")    # map of the package
```

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

## Authors

| Author | Role | ORCID | GitHub |
|---|---|---|---|
| Michela Cameletti | aut, cph | [0000-0002-6502-7779](https://orcid.org/0000-0002-6502-7779) | [@michelacameletti](https://github.com/michelacameletti) |
| Francesco Caccia | aut, cph | | |
| Paolo Maranzano | aut, cre, cph | [0000-0002-9228-2759](https://orcid.org/0000-0002-9228-2759) | [@PaoloMaranzano](https://github.com/PaoloMaranzano) |

Licensed under GPL (>= 2).
