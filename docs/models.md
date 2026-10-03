# The models: STEM and SC-STEM

[Back to the README](../README.md) | [All pages](README.md)

The model estimated by the package, its spatially-clustered extension, and what the model class covers and does not.


## STEM: one process, three sources of variation

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
Newton-Raphson step. By default the iterations are accelerated by SQUAREM;
the plain EM algorithm, the ECME algorithm (which updates $\beta$ and $m_0$ on
the observed likelihood) and SQUAREM on ECME are the alternatives, all reaching
the same maximum (`STEM_control(algorithm = )`).

## SC-STEM: spatial regimes

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

## Statistical features and scope

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

---

[Back to the README](../README.md) | [All pages](README.md)
