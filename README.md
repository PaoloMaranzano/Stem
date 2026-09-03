
# Stem <img src="man/figures/logo.svg" align="right" height="140" alt="Stem logo" />

`Stem` provides tools for estimating hierarchical spatio-temporal models using the Expectation-Maximization (EM) algorithm in a likelihood-based framework. The package is designed for spatio-temporal panel data where observations are indexed over space and time and where spatial dependence, latent dynamics, and prediction uncertainty are relevant modelling components.

The current development version also includes spatially-clustered STEM (SC-STEM) routines. SC-STEM extends the classical STEM workflow by allowing regression parameters to vary across spatial clusters, so that groups of monitoring locations can be characterized by homogeneous local dynamics rather than being forced into a single global regression structure.

## Main features

- Construction of STEM model objects from observations, covariates, coordinates, and starting parameters.
- Maximum-likelihood estimation through an EM algorithm.
- Kalman filtering and smoothing utilities used internally by the estimation routines.
- Spatial prediction and kriging-style mapping utilities.
- Parametric bootstrap for uncertainty assessment.
- Spatially-clustered STEM estimation and information criteria for selecting the number of clusters.

## Installation

After creating the GitHub repository, the development version can be installed with:

```r
install.packages("remotes")
remotes::install_github("paolomaranzano/Stem")
```

## Basic workflow

```r
library(Stem)
data(pm10)

coordinates <- pm10$coords * 1000
covariates <- pm10$covariates
z <- pm10$z

phi <- list(
  beta = matrix(c(3.65, 0.046, -0.904), 3, 1),
  sigma2eps = 0.1,
  sigma2omega = 0.2,
  theta = 0.01,
  G = matrix(0.77, 1, 1),
  Sigmaeta = matrix(0.3, 1, 1),
  m0 = as.matrix(0),
  C0 = as.matrix(1)
)

K <- matrix(1, ncol(z), 1)

mod <- STEM_Model(
  z = z,
  covariates = covariates,
  coordinates = coordinates,
  phi = phi,
  K = K
)

fit <- STEM_Estimation(mod, distance = "euclidean")
```

## Spatially-clustered STEM

```r
sc_fit <- SCSTEM_Estim(
  StemModel = mod,
  distance = "euclidean",
  crs = 32632,
  k = 3,
  verbose = TRUE
)
```

The `verbose` argument controls progress messages. By default, estimation functions avoid unrequested console output.

## Function map

The package includes a supporting function map in `inst/extdata/`:

- `STEM_function_map.svg`: updated call graph generated from the current R source files;
- `STEM_function_map.dot`: Graphviz source used to render the updated map;
- `STEM_function_map_original.pdf`: original conceptual map supplied with the development materials.

## Development notes

This version implements the first round of CRAN-oriented fixes requested by maintainers: the package title no longer ends with the redundant phrase “in R”; references in `DESCRIPTION` are linked using DOI, ISBN, or URL syntax; `T`/`F` have been replaced by `TRUE`/`FALSE`; direct `cat()` diagnostics have been replaced by suppressible `message()`/`warning()` calls; and the top-level `.Random.seed` assignment has been removed.
