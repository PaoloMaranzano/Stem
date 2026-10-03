# Stem: documentation pages

[Back to the README](../README.md)

The README of the repository is the overview. These pages take one topic each.

| page | what it covers |
|---|---|
| [The models: STEM and SC-STEM](models.md) | the model, its parameters and their reading; the spatially-clustered extension and its Potts-penalized objective; what the model class covers and does not, including missing values |
| [Computational aspects](computational-aspects.md) | how the models are estimated: the Kalman filter and smoother, the EM algorithm and where it is slow, ECME and SQUAREM, the stopping rule compared with D-STEM, the limits of the range, the numerical experiments |
| [Ridge, lasso and elastic net](regularization.md) | the optional penalty on the regression coefficients, how it enters the EM algorithm, and how to choose it |
| [Design notes](sc-stem-design.md) | the implementation choices of SC-STEM: starts, the label update, degeneracy, the bootstrap, the cost of the filter |
| [How the functions fit together](function-map.md) | the map of the package |
| [Theoretical references](references.md) | the works the package implements, adapts or follows |

Two documents ship with the package:

```r
# the Computational Supplement: estimation algorithms, derivations, numerical experiments
browseURL(system.file("extdata", "Stem-computational-supplement.pdf", package = "Stem"))
# a reading guide from each equation of the papers to the function that implements it
browseURL(system.file("extdata", "STEM_model_verification.html", package = "Stem"))
```

and four vignettes:

```r
vignette("getting-started", package = "Stem")      # the classical STEM workflow
vignette("SCSTEM", package = "Stem")               # spatially-clustered STEM
vignette("function-map", package = "Stem")         # map of the package
vignette("computational-notes", package = "Stem")  # how the filter and the M-step are computed
```
