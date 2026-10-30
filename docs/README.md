# Stem: documentation pages

[Back to the README](../README.md)

The README of the repository is the overview. These pages take one topic each.
They are the full account: the supplements of the paper report a subset of
what is here, and with the same numbers.

| page | what it covers |
|---|---|
| [The models: STEM and SC-STEM](models.md) | the model, its parameters and their reading; the spatially-clustered extension and its Potts-penalized objective; what the model class covers and does not, including missing values |
| [Computational aspects](computational-aspects.md) | how the models are estimated, following Supplement A of the paper section by section: the Kalman filter and smoother, the EM algorithm, its stopping rule and where it is slow, ECME and SQUAREM, the limits of the range, missing responses; for SC-STEM, the start of a fit and the scale of the Potts penalty, degeneracy and the swap pass, the label update and the monotonicity of the alternation, the information criteria, the refit-with-clustering bootstrap and its tests, the Potts penalty on point-referenced data; the numerical experiments and where everything is in the code |
| [Ridge, lasso and elastic net](regularization.md) | the optional penalty on the regression coefficients, how it enters the EM algorithm and the accelerations, and how to choose it |
| [How the functions fit together](function-map.md) | the map of the package |
| [Theoretical references](references.md) | the works the package implements, adapts or follows |

Two documents ship with the package:

```r
# Supplement A of the paper: the estimation algorithms and their derivations
browseURL(system.file("extdata", "Stem-computational-supplement.pdf", package = "Stem"))
# a reading guide from each equation of the papers to the function that implements it
browseURL(system.file("extdata", "STEM_model_verification.html", package = "Stem"))
```

and five vignettes:

```r
vignette("getting-started", package = "Stem")        # the classical STEM workflow
vignette("SCSTEM", package = "Stem")                 # spatially-clustered STEM
vignette("function-map", package = "Stem")           # map of the package
vignette("computational-notes", package = "Stem")    # how the filter and the M-step are computed
vignette("estimation-experiments", package = "Stem") # the numerical experiments behind the defaults
```
