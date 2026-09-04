# Stem <img src="man/figures/logo.svg" align="right" height="132" alt="" />

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

For a network of `d` locations observed over `T` time points, the STEM model is
a hierarchical state-space model with a measurement equation and a transition
equation:

```math
z_t = X_t \beta + K y_t + e_t, \qquad e_t \sim N(0, \Sigma_e)
```

```math
y_t = G y_{t-1} + \eta_t, \qquad \eta_t \sim N(0, \Sigma_\eta), \qquad y_0 \sim N(m_0, C_0)
```

```math
\Sigma_e = \sigma^2_\varepsilon I + \sigma^2_\omega \exp(-\theta h)
```

so the network shares a latent temporal process `y_t`, loaded on the locations
through `K`, on top of a spatially correlated error with a nugget.

SC-STEM assigns each location to one of `k` latent regimes and fits a separate
STEM model within each of them, estimating labels and parameters jointly by
maximising a Potts-penalised log-likelihood:

```math
Q = \sum_{i=1}^{d} \ell_{i k_i} + \phi\, c \sum_{i<j} w_{ij}\, \mathbb{1}(k_i = k_j)
```

where `phi >= 0` tunes how strongly neighbouring locations are pushed into the
same regime, and `c` is a scale factor documented in `?SCSTEM_Estim`. Setting
`k = 1` returns the pooled model.

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
Tfull <- nrow(povalley$z); d <- ncol(povalley$z)
keep <- as.vector(outer(seq_len(Tn), (seq_len(d) - 1L) * Tfull, "+"))

phi <- list(beta = matrix(c(1.25, -0.00003, 0.64), 3, 1),
            sigma2eps = 18.66, sigma2omega = 1e-06, theta = 2e-06,
            G = matrix(0.59, 1, 1), Sigmaeta = matrix(4.25, 1, 1),
            m0 = as.matrix(0), C0 = as.matrix(1))

mod <- STEM_Model(z = povalley$z[seq_len(Tn), ],
                  covariates = povalley$covariates[keep, ],
                  coordinates = povalley$coords,
                  phi = phi, K = matrix(1, d, 1))

# explore the grid and let the two-step rule choose k and phi
ic  <- SCSTEM_Infocrit(mod, k_grid = 1:4, phi_grid = seq(0, 1, by = 0.25),
                       distance = "geo")
sel <- SCSTEM_Select(ic, band = c(0.25, 1))
sel

# uncertainty, with the partition re-estimated at every draw
boot <- SCSTEM_Bootstrap(sel$fit, B = 200, seed = 1)
SCSTEM_BootInference(boot)
```

## Design notes

Three points make SC-STEM behave sensibly in practice, and are documented in
detail in `?SCSTEM_Estim`.

**Labels are updated sequentially (ICM).** Each location maximises its own
penalised contribution given the current labels of all the others, so a sweep
cannot decrease the objective at fixed parameters and cannot cycle — unlike a
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

The implementation has been checked equation by equation against the three
reference papers. That verification ships with the package as a standalone
document:

```r
browseURL(system.file("extdata", "STEM_model_verification.html", package = "Stem"))
```

## Theoretical references

The STEM model implemented here comes from three companion works by Fassò and
Cameletti, which play different roles.

**Fassò, A., Cameletti, M. and Nicolis, O. (2007).** Air quality monitoring
using heterogeneous networks. *Environmetrics*, 18, 245–264.
<https://doi.org/10.1002/env.837>
> Introduces the *geostatistical dynamical calibration* (GDC) model, the most
> general of the three: it adds instrument calibration components — an additive
> bias `A(t)` and a multiplicative bias `B(t)` — so that a network of
> heterogeneous instruments (gravimetric and TEOM monitors) can be modelled
> jointly, with the loading matrix obtained by empirical orthogonal functions.
> The model estimated by this package is the special case with no calibration
> bias.

**Fassò, A. and Cameletti, M. (2007).** A general spatio-temporal model for
environmental data. *GRASPA Technical Report* n. 27.
> The direct theoretical reference for this package, which it announces by name.
> It states the three-stage hierarchy, the scaled spatial covariance
> `Gamma(h) = 1 + gamma` at `h = 0` and `C_theta(h)` otherwise with
> `gamma = sigma2eps/sigma2omega`, the choice of estimating `log(gamma)` rather
> than `sigma2eps` for positive-definiteness, the EM algorithm with closed-form
> M-steps for `beta`, `sigma2omega`, `G`, `Sigma_eta` and `m0` and
> Newton-Raphson for the spatial covariance parameters, and the spatio-temporal
> parametric bootstrap.

**Fassò, A. and Cameletti, M. (2010).** A unified statistical approach for
simulation, modeling, analysis and mapping of environmental data. *Simulation*,
86, 139–153. <https://doi.org/10.1177/0037549709102150>
> The most complete published statement: the same model and EM algorithm, plus
> the kriging predictor used by `STEM_Kriging()`, the bootstrap of
> `STEM_Bootstrap()`, and a sensitivity analysis of the model components. Its
> Equations (12)–(18) are what the estimation code implements.

The spatially-clustered layer adapts algorithmic devices from a separate
literature; the statistical model remains the STEM one above.

**Besag, J. (1986).** On the statistical analysis of dirty pictures.
*JRSS-B*, 48, 259–302.
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

## Authors

Michela Cameletti (aut),
[Paolo Maranzano](https://orcid.org/0000-0002-9228-2759) (aut, cre).

Licensed under GPL (>= 2).
