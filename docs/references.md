# Theoretical references

[Back to the README](../README.md) | [All pages](README.md)

The works the package implements, adapts or follows, with the role each plays.


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

### Estimation algorithms

**Dempster, A. P., Laird, N. M. and Rubin, D. B. (1977).** Maximum likelihood from incomplete data via the EM algorithm. *JRSS-B*, 39, 1-38.
> The EM algorithm, and the rate of convergence given by the fraction of missing information.

**Shumway, R. H. and Stoffer, D. S. (1982).** An approach to time series smoothing and forecasting using the EM algorithm. *Journal of Time Series Analysis*, 3, 253-264. <https://doi.org/10.1111/j.1467-9892.1982.tb00349.x>
> The EM algorithm for linear Gaussian state-space models, with the lag-one covariance smoother.

**Liu, C. and Rubin, D. B. (1994).** The ECME algorithm: a simple extension of EM and ECM with faster monotone convergence. *Biometrika*, 81, 633-648. <https://doi.org/10.1093/biomet/81.4.633>
> The ECME algorithm offered by `STEM_control(algorithm = "ECME")`.

**de Jong, P. (1991).** The diffuse Kalman filter. *The Annals of Statistics*, 19, 1073-1083. <https://doi.org/10.1214/aos/1176348139>
> The augmented Kalman filter that computes the generalized least-squares step of ECME.

**Varadhan, R. and Roland, C. (2008).** Simple and globally convergent methods for accelerating the convergence of any EM algorithm. *Scandinavian Journal of Statistics*, 35, 335-353. <https://doi.org/10.1111/j.1467-9469.2007.00585.x>
> SQUAREM, the default algorithm of the package.

**Du, Y. and Varadhan, R. (2020).** SQUAREM: an R package for off-the-shelf acceleration of EM, MM and other EM-like monotone algorithms. *Journal of Statistical Software*, 92(7). <https://doi.org/10.18637/jss.v092.i07>
> The steplength scheme the package follows.

**Wang, Y., Finazzi, F. and Fasso, A. (2021).** D-STEM v2: a software for modeling functional spatio-temporal data. *Journal of Statistical Software*, 99(10). <https://doi.org/10.18637/jss.v099.i10>
> The software whose stopping rule the package compares with its own.

---

[Back to the README](../README.md) | [All pages](README.md)
