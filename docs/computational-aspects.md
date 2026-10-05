# Computational aspects of Stem

How `Stem` estimates its models: the state-space form, the Kalman filter and
smoother every algorithm rests on, the EM algorithm and where it is slow, the
two accelerations the package offers (ECME and SQUAREM), the stopping rule and
how it relates to that of D-STEM, the limits within which the range is
estimated, and the numerical experiments behind these choices.

The full technical account, with the proofs and the complete tables, is
Section S1 of the **Computational Supplement**, which ships with the package:

```r
browseURL(system.file("extdata", "Stem-computational-supplement.pdf", package = "Stem"))
```

Contents

1. [The model in state-space form](#1-the-model-in-state-space-form)
2. [The Kalman filter and smoother](#2-the-kalman-filter-and-smoother)
3. [The EM algorithm](#3-the-em-algorithm)
4. [The EM algorithm within SC-STEM](#4-the-em-algorithm-within-sc-stem)
5. [Where the EM algorithm is slow](#5-where-the-em-algorithm-is-slow)
6. [ECME](#6-ecme)
7. [SQUAREM](#7-squarem)
8. [The four algorithms at a glance](#8-the-four-algorithms-at-a-glance)
9. [The stopping rule, and the rule of D-STEM](#9-the-stopping-rule-and-the-rule-of-d-stem)
10. [The range at the boundary of what the data identify](#10-the-range-at-the-boundary-of-what-the-data-identify)
11. [No ridge in the inversions](#11-no-ridge-in-the-inversions)
12. [Numerical experiments](#12-numerical-experiments)
13. [Where it is in the code](#13-where-it-is-in-the-code)

## 1. The model in state-space form

Within a regime, or in the pooled fit, STEM is the state-space model

```math
z_t = X_t \beta + A\, y_t + e_t, \qquad e_t \sim N(0, \Sigma_e), \qquad \Sigma_e = \sigma^2_\omega C_\theta + \sigma^2_\epsilon I_d,
```

```math
y_t = G\, y_{t-1} + \eta_t, \qquad \eta_t \sim N(0, \Sigma_\eta), \qquad y_0 \sim N(m_0, C_0),
```

for $t = 1, \dots, T$, with $z_t$ the $d$ observations at time $t$, $A$ the
$d \times p$ loading matrix and $C_0$ known. The package writes
$\Sigma_e = \sigma^2_\omega \Gamma$ with $\Gamma = C_\theta + \gamma I_d$ and
$\gamma = \sigma^2_\epsilon / \sigma^2_\omega$, and estimates
$\psi = (\beta, \sigma^2_\omega, \log\gamma, \log\theta, G, \Sigma_\eta, m_0)$
by the EM algorithm of Shumway and Stoffer (1982): the E-step is the Kalman
smoother, the M-step updates $G$, $\Sigma_\eta$, $m_0$, $\sigma^2_\omega$ and
$\beta$ in closed form and $(\log\theta, \log\gamma)$ by Newton-Raphson. The
log-likelihood $\ell(\psi)$ is the by-product of the Kalman filter (the
prediction-error decomposition).

## 2. The Kalman filter and smoother

Every algorithm of the package rests on three recursions, run at given
parameters by the internal function `kalman()`.

**The filter.** From $y_{0|0} = m_0$, $P_{0|0} = C_0$, for $t = 1, \dots, T$:

```math
y_{t|t-1} = G y_{t-1|t-1}, \quad P_{t|t-1} = G P_{t-1|t-1} G^\top + \Sigma_\eta, \quad
v_t = z_t - X_t \beta - A y_{t|t-1}, \quad F_t = A P_{t|t-1} A^\top + \Sigma_e,
```

```math
K_t = P_{t|t-1} A^\top F_t^{-1}, \quad y_{t|t} = y_{t|t-1} + K_t v_t, \quad P_{t|t} = P_{t|t-1} - K_t A P_{t|t-1},
```

```math
\ell(\psi) = -\tfrac12 \sum_{t=1}^T \left\{ d \log(2\pi) + \log|F_t| + v_t^\top F_t^{-1} v_t \right\}.
```

The $d \times d$ matrix $F_t$ is never formed: $\Sigma_e$ does not depend on $t$
and $A P_{t|t-1} A^\top$ has rank $p$, so the Woodbury identity and the matrix
determinant lemma reduce every step to $p \times p$ algebra against a Cholesky
factor of $\Sigma_e$ computed once per pass. A pass costs $O(d^3 + T d^2 p)$
instead of $O(T d^3)$, and the log-density is evaluated directly, since on a
few hundred locations the density itself underflows. Missing responses are
handled by dropping the unobserved rows of the measurement equation (Durbin and
Koopman 2012, Sect. 4.10); the vignette `computational-notes` gives the
details.

**The smoother.** With $J_t = P_{t|t} G^\top P_{t+1|t}^{-1}$, for
$t = T-1, \dots, 0$:

```math
y_{t|T} = y_{t|t} + J_t (y_{t+1|T} - y_{t+1|t}), \qquad P_{t|T} = P_{t|t} + J_t (P_{t+1|T} - P_{t+1|t}) J_t^\top .
```

**The lag-one covariance smoother** (Shumway and Stoffer 1982):

```math
P_{T,T-1|T} = (I - K_T A) G P_{T-1|T-1}, \qquad
P_{t-1,t-2|T} = P_{t-1|t-1} J_{t-2}^\top + J_{t-1} (P_{t,t-1|T} - G P_{t-1|t-1}) J_{t-2}^\top .
```

## 3. The EM algorithm

**E-step.** The expected complete-data log-likelihood depends on the states
only through the smoothed initial state and

```math
S_{11} = \sum_t (y_{t|T} y_{t|T}^\top + P_{t|T}), \quad
S_{00} = \sum_t (y_{t-1|T} y_{t-1|T}^\top + P_{t-1|T}), \quad
S_{10} = \sum_t (y_{t|T} y_{t-1|T}^\top + P_{t,t-1|T}),
```

```math
W = \sum_t (A P_{t|T} A^\top + r_t r_t^\top), \qquad r_t = z_t - X_t \beta - A y_{t|T},
```

completed by their conditional expectations where responses are missing.

**M-step.** The parameters are updated in turn, each given the updates that
precede it (an ECM algorithm): $m_0 \leftarrow y_{0|T}$; $\Sigma_\eta$ and $G$
from $S_{11}$, $S_{10}$, $S_{00}$ (element by element when diagonal, the
default; otherwise
$\Sigma_\eta \leftarrow T^{-1}(S_{11} - S_{10} S_{00}^{-1} S_{10}^\top)$ and
$G \leftarrow S_{10} S_{00}^{-1}$); $\sigma^2_\omega \leftarrow \mathrm{tr}(\Gamma^{-1} W)/(dT)$;
$\beta$ by generalized least squares given the smoothed states,

```math
\beta \leftarrow \Big(\sum_t X_t^\top \Sigma_e^{-1} X_t\Big)^{-1} \sum_t X_t^\top \Sigma_e^{-1} (\hat z_t - A y_{t|T}),
```

(or the corresponding penalized quadratic under a ridge or an elastic net);
and $(\log\theta, \log\gamma)$ by a safeguarded Newton-Raphson step on
$Q_e = T \log|\sigma^2_\omega \Gamma| + \mathrm{tr}\{(\sigma^2_\omega \Gamma)^{-1} W\}$
(Section 10). No step decreases the expected complete-data log-likelihood, so
the algorithm is a generalized EM algorithm and $\ell$ does not decrease.

## 4. The EM algorithm within SC-STEM

SC-STEM fits one STEM model per regime and alternates between the parameters
and the labels. The EM algorithm enters twice.

- **In the alternation**, where it only has to rank the regimes: a few plain
  EM (or ECME) iterations per regime, tolerances `1e-2` and `1`, at most 50
  iterations, never accelerated. A regime first starts from the least-squares
  coefficients of its own locations (an unpenalized fit) or from the final
  refit of the unpenalized solution (a penalized fit), and afterwards from its
  own estimates of the previous iteration.
- **In the final refit**, where the regimes are re-estimated once on the final
  partition, with the settings of the pooled fit and the algorithm chosen,
  starting from the last estimates of the alternation. The reported estimates
  and log-likelihoods come from this refit; fits of a grid that end at the same
  partition share one refit.

## 5. Where the EM algorithm is slow

The EM algorithm is a fixed-point iteration $\psi^{(j+1)} = M(\psi^{(j)})$
which, near the maximum, is linear:
$\psi^{(j+1)} - \hat\psi \approx J (\psi^{(j)} - \hat\psi)$ with
$J = I_c^{-1} I_m$, the fraction of the complete-data information that is
missing because the states are not observed (Dempster, Laird and Rubin 1977).
Convergence is linear with rate $\lambda$, the largest eigenvalue of $J$, and a
stopping rule that looks at the change between two iterations leaves a distance
from the maximum of about that change divided by $1 - \lambda$. On the regimes
of SC-STEM three situations make $\lambda$ close to one or the objective flat:

1. **The split of the common variation between the latent process and the
   field.** The latent process is common to the locations of a regime and the
   field, with a range comparable to its size, nearly so; only their dynamics
   tell them apart. Rate about 0.93 on the regimes of the simulation study:
   this direction sets the number of iterations.
2. **The intercept and the initial level of the latent process.** A change of
   $m_0$ moves the first periods, which the intercept and the smoothed states
   partly absorb, and the EM algorithm corrects only a small part of the error
   at every iteration. Rate about 0.988, little likelihood: this is where the
   error at the stop sits. With $G$ near zero, $G$ and $m_0$ trade along a
   ridge.
3. **The range at the boundary of what the distances identify** (Section 10).

## 6. ECME

ECME (Liu and Rubin 1994) lets some conditional maximizations of the M-step
maximize the observed log-likelihood instead of the expected complete-data one.
`Stem` updates $\delta = (\beta^\top, m_0^\top)^\top$ in this way, after the
other parameters, which keeps the algorithm monotone. Given the other
parameters the log-likelihood depends on $\delta$ only through the mean of the
data, $E(z_t) = X_t \beta + A G^t m_0$, so it is a quadratic in $\delta$ whose
maximizer is a generalized least-squares estimator. The augmented Kalman filter
(de Jong 1991) computes it without the $dT \times dT$ covariance: the filter of
the current parameters is run on every column of $[\,z_t,\ X_t,\ A G^t\,]$, the
columns share the gains and $F_t$, and with $[\,e_t,\ E_t\,]$ their innovations,

```math
\hat\delta = \Big(\sum_t E_t^\top F_t^{-1} E_t\Big)^{-1} \sum_t E_t^\top F_t^{-1} e_t .
```

The step removes the slow direction of the intercept and $m_0$ (and the ridge
between $G$ and $m_0$) in one go, at about 25% more time per iteration, and
leaves the split of the variance as it is.

## 7. SQUAREM

SQUAREM (Varadhan and Roland 2008; the steplength of the R package `SQUAREM`,
Du and Varadhan 2020) accelerates the EM, or the ECME, iterations by
extrapolating from two of them. From $\psi_0$, with $\psi_1 = M(\psi_0)$,
$\psi_2 = M(\psi_1)$ and $u$ the parameters on a scale on which every value is
admissible (logarithms of the variances and of $\theta$):

```math
r = u_1 - u_0, \quad w = u_2 - 2u_1 + u_0, \quad
\alpha = \min\{\alpha_{\max}, \max(1, \|r\|/\|w\|)\}, \quad u' = u_0 + 2\alpha r + \alpha^2 w .
```

For a linear iteration with rate $\lambda$, $\|r\|/\|w\| = 1/(1-\lambda)$ and
$u'$ is the maximum itself. One more iteration from $\psi'$ stabilizes the
extrapolation, which is kept only if $\ell(\psi') \ge \ell(\psi_0)$ (otherwise
the cycle ends at $\psi_2$); the safeguard costs nothing, since every iteration
evaluates $\ell$ at its input. $\alpha_{\max}$ starts at 1, grows by a factor 4
when it binds and shrinks after a rejection. SQUAREM is not used inside the
alternation of SC-STEM, whose fits stop after a few iterations.

With a ridge on $\beta$ the objective is the penalized log-likelihood
$\ell(\psi) - \tfrac{\lambda}{2}\sum_j w_j s_j(\psi)^2 \beta_j^2$, whose scale
$s_j^2 = [\sum_t X_t'\Sigma_e^{-1}X_t]_{jj}$ moves with $\Sigma_e$: there is no
fixed objective, and the safeguard compares this one with the scale measured
at each point. Along the plain iterations it decreases by at most $10^{-5}$,
and on dynamic designs with collinear lags (five seeds, three penalties) the
accelerated iterations reach the fixed point of the plain ones to $10^{-7}$ in
about half the time. With the lasso or the elastic net the fits run the
plain iterations.

## 8. The four algorithms at a glance

All four maximize the same likelihood and share the E-step, the Kalman filter
and smoother, which in a linear Gaussian state-space model is exact inference
on the latent states. In the language of machine learning, EM is coordinate
ascent on the evidence lower bound, made tight by the E-step; ECME replaces
some of the coordinate updates by "collapsed" ones, with the states integrated
out; SQUAREM is an extrapolation of the fixed-point map, in the family of
Aitken's and Anderson's accelerations, with a monotonicity safeguard.

| `STEM_control(algorithm = )` | M-step | iterations combined | time | iterations | shortfall |
|---|---|---|---|---|---|
| `"EM"` | conditional maximizations of $Q$ | plain | 1 | 1 | 0.020 |
| `"ECME"` | as EM, with $(\beta, m_0)$ on $\ell$ | plain | 1.17 | 0.96 | 0.017 |
| **`"SQUAREM"`** (default) | as EM | two, extrapolated | **0.53** | 0.41 | 0.0027 |
| `"SQUAREM-ECME"` | as ECME | two, extrapolated | 0.62 | 0.39 | 0.0023 |

Time and iterations relative to EM on whole SC-STEM grids of six cells of the
simulation study (two replications each); shortfall: median distance of the
log-likelihood from the maximum at the stop, on single fits; all at the default
tolerances.

## 9. The stopping rule, and the rule of D-STEM

Two criteria are computed at every iteration (at every cycle under SQUAREM):
the largest relative change of a free parameter, taken one at a time, with a
floor of $10^{-3}$ on the denominator, and the absolute change of the
log-likelihood. The algorithm stops when either falls below its tolerance
(`em_stop = "any"`), or after `em_maxit` iterations.

| | Stem | D-STEM v2 |
|---|---|---|
| parameters | $\max_l \lvert\Delta\psi_l\rvert / \max(\lvert\psi_l\rvert, 10^{-3}) < 10^{-3}$ | $\max_l \lvert\Delta\psi_l\rvert / \lvert\psi_l\rvert < 10^{-4}$ |
| log-likelihood | $\lvert\Delta\ell\rvert < 10^{-3}$ (absolute) | $\lvert\Delta\ell\rvert / \lvert\ell\rvert < 10^{-4}$ (relative) |
| combination | either criterion | either criterion |
| iterations | at most 500 (50 in the alternation) | at most 100 |

D-STEM is the stricter on the parameters and by far the looser on the
log-likelihood: with log-likelihoods between $-4000$ and $-15000$ its relative
tolerance is 0.4 to 1.5 units, 400 to 1500 times the absolute tolerance of
`Stem`. Since either criterion stops the iterations, D-STEM stops while several
units can still be gained; `Stem` with SQUAREM stops about $0.003$ units from
the maximum. The criterion on $\ell$ is absolute because the differences that
matter, between fits and in the information criteria, are absolute.

## 10. The range at the boundary of what the data identify

On a small set of locations the likelihood is flat in $\theta$ beyond two
limits: when the range is shorter than the distance between the two closest
locations the field is white and cannot be told from the nugget; when it is
much longer than the largest distance the field is constant over the locations
and cannot be told from an effect common to them. Left free, the Newton-Raphson
step took steps of any size ($\theta$ from 0.02 to $e^{-700}$ in one M-step was
seen), every later iteration paid up to thirty grid searches, and SQUAREM ran
to the limit on the iterations. `Stem` therefore keeps

```math
\frac{-\log 0.95}{h_{\max}} \;\le\; \theta \;\le\; \frac{-\log 0.05}{h_{\min}},
```

with $h_{\min}$ and $h_{\max}$ the smallest and largest distance between two
locations of the fit, in the Newton-Raphson step and in the extrapolation of
SQUAREM; a Newton step is taken only with a positive definite Hessian and
accepted only if it does not worsen the objective (halved otherwise); and one
grid search per step is allowed, after which a Hessian that is still not
positive definite ends the step. A fit whose range ends at a limit reports it
(`convergence.par$theta.bound` of `STEM_Estimation()`, `theta_bound` of
`SCSTEM_Estimation()`).

Where nothing binds the results are unchanged (identical to $10^{-9}$ on the
fits checked). On whole SC-STEM grids the fits at the true number of regimes
were unchanged in 40 of 48 cases and better in the other 8 (by 0.3 to 16.5
log-likelihood units), the selections unchanged, and the total time fell by
19% for EM and by 40% for SQUAREM; on the replication that had cost most, from
362 to 29 seconds for EM and from 1416 to 25 for SQUAREM.

## 11. No ridge in the inversions

Earlier versions added $0.01 I$ to every matrix the iterations inverted, while
the filter that computes $\ell$ added nothing. The iterations then converged to
a point that is not the maximum, displaced along the slow directions by about
$1/(1-\lambda)$ times the perturbation: the variance parameters were off by 3
to 13% on the regimes of the simulation study and $\Sigma_\eta$ of a pooled fit
by 84%, at any tolerance. Every matrix inverted is positive definite at
admissible parameters, and `regularization = 0` is the default.

## 12. Numerical experiments

All on the data-generating process of the simulation study of the paper.

**Correctness.** With the plain EM algorithm and the old ridge, the current
code reproduces the previous implementation to $10^{-13}$, iteration by
iteration, on complete and incomplete data, without latent process or spatial
correlation, under a ridge and an elastic net, and with $p = 2$. The ECME step
reproduces the GLS estimator on the explicit $dT \times dT$ covariance to
$10^{-10}$. Without the ridge, the iterations converge to the BFGS maximum of
the exact log-likelihood.

**Exactness.** Against computations on the full covariance of the data, with
gaps in the response: the log-likelihood returned is the Gaussian
log-likelihood of the observed values to $10^{-14}$, for EM, ECME and SQUAREM;
the estimates are the maximum that `optim()` finds on it (largest relative
difference $2 \cdot 10^{-5}$); the smoothed states and `STEM_Signal()` are the
conditional means to $3 \cdot 10^{-15}$. `STEM_Simulation()`, which the
bootstraps use, reproduces the mean and the covariance of the model within
Monte Carlo error over 20000 data sets. The first four checks are tests of the
package (`tests/testthat/test-exact.R`). Every version up to 2.0.0, 1.0
included, returned four times the log-likelihood in `estimates$loglik` (the
EM loop stored $-2\ell$ and multiplied it by $-2$ again), so the information
criteria of SC-STEM weighted the likelihood four times against the penalty;
the value returned is now $\ell$.

**Information criteria.** AIC, BIC and KIC use the exact total log-likelihood
of the final refit and $K_{\mathrm{eff}}(r + 3 + 3p)$ parameters; the sample
size of the BIC is the number of observed values of the response, $dT$ when
the panel is complete.

**The same maximum.** Run to a tolerance of $10^{-7}$, the four algorithms stop
at the same log-likelihood to $10^{-9}$, after 177 (EM), 178 (ECME), 39
(SQUAREM) and 30 (SQUAREM-ECME) iterations.

**Single fits** (three regimes and the pooled model of the reference cell, ten
replications), default tolerances:

| algorithm | seconds per replication | iterations per fit | median shortfall of $\ell$ | median largest relative error |
|---|---|---|---|---|
| EM | 12.7 | 89.0 | 0.020 | 3.8% |
| ECME | 15.2 | 83.5 | 0.017 | 3.2% |
| **SQUAREM** | **4.6** | 34.3 | 0.0027 | 1.4% |
| SQUAREM-ECME | 5.2 | 29.6 | 0.0023 | 1.1% |

**Accuracy against the true values** (relative RMSE of the regimes, %): at the
default tolerances SQUAREM matches the exact maximum on every parameter (for
instance 18.8% against 18.5% on $\theta$), while EM, ECME and SQUAREM-ECME,
stopping farther from it, lose accuracy on $\theta$ (24.5 to 24.9%) and on
$\sigma^2_\epsilon$. The regression coefficients are unaffected.

**Whole SC-STEM grids** (six cells from $n = 40$, $T = 60$ to $n = 400$ and
$n = 200$, $T = 365$, two replications; the four algorithms run in the same
process, time net of the position): SQUAREM 0.53 times the time of EM (0.58 on
the light cells, 0.46 on the heavy ones), SQUAREM-ECME 0.62, ECME 1.17; no fit
at the limit on the iterations; the same number of regimes selected in every
job.

## 13. Where it is in the code

| file | content |
|---|---|
| `R/kalman.R` | `kalman()`: the filter (`filtering()`, `filterstep()`), the smoother (`smoothing()`), the smoothed initial state and the lag-one covariances |
| `R/estep.R` | `stem_estep()`: the expected sufficient statistics |
| `R/mstep-updates.R` | the conditional updates: latent process, $\sigma^2_\omega$, $\beta$, the spatial parameters (safeguarded Newton-Raphson), the limits of $\theta$ |
| `R/mstep-em.R`, `R/mstep-ecme.R` | the M-steps of EM and ECME; the augmented filter `stem_gls_mean()` |
| `R/squarem.R` | SQUAREM and the unconstrained parametrization |
| `R/em-fit.R` | the wrapper `stem_em_fit()`, the plain iterations, the stopping rule |
| `R/STEM_control.R` | the settings: `algorithm`, `em_tol_par`, `em_tol_loglik`, `em_maxit`, `em_stop`, `alt_em_*`, `alt_maxit`, `nr_maxit` |

## References

- de Jong, P. (1991). The diffuse Kalman filter. *The Annals of Statistics*, 19, 1073-1083.
- Dempster, A. P., Laird, N. M. and Rubin, D. B. (1977). Maximum likelihood from incomplete data via the EM algorithm. *JRSS-B*, 39, 1-38.
- Du, Y. and Varadhan, R. (2020). SQUAREM: an R package for off-the-shelf acceleration of EM, MM and other EM-like monotone algorithms. *Journal of Statistical Software*, 92(7).
- Durbin, J. and Koopman, S. J. (2012). *Time Series Analysis by State Space Methods*, 2nd ed. Oxford University Press.
- Liu, C. and Rubin, D. B. (1994). The ECME algorithm. *Biometrika*, 81, 633-648.
- Meng, X.-L. and Rubin, D. B. (1993). Maximum likelihood estimation via the ECM algorithm. *Biometrika*, 80, 267-278.
- Shumway, R. H. and Stoffer, D. S. (1982). An approach to time series smoothing and forecasting using the EM algorithm. *Journal of Time Series Analysis*, 3, 253-264.
- Varadhan, R. and Roland, C. (2008). Simple and globally convergent methods for accelerating the convergence of any EM algorithm. *Scandinavian Journal of Statistics*, 35, 335-353.
- Wang, Y., Finazzi, F. and Fasso, A. (2021). D-STEM v2: a software for modeling functional spatio-temporal data. *Journal of Statistical Software*, 99(10).
