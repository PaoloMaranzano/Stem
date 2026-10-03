# Ridge, lasso and elastic net on the coefficients

[Back to the README](../README.md) | [All pages](README.md)

The optional penalty on the regression coefficients, how it enters the EM algorithm, and how to choose it.


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
                latent = FALSE, spatial = FALSE)
```

`regularization` must stay at its default of 0: a ridge added for conditioning
would be the only thing separating the two, and on a collinear design it is not
negligible. The agreement with an elastic net
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
[`dev/regularization/stem-elastic-net.tex`](../dev/regularization/stem-elastic-net.tex).

---

[Back to the README](../README.md) | [All pages](README.md)
