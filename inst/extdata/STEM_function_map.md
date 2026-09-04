# Map of the Stem package

Conceptual map of how the functions of the package call one another. The
diagram is in `STEM_function_map.svg`; the narrative version, with an
explanation of each layer, is the vignette:

```r
vignette("function-map", package = "Stem")
```

Exported functions are marked with `*`. The edge list below is obtained by
parsing the sources in `R/` and recording, for every function defined in the
package, which other package functions appear in its body. Two edges that a
purely static scan cannot see are added by hand and marked `[indirect]`: they go
through `lapply()` and `do.call()`.

## Layers

| Layer | Functions |
|---|---|
| Data | `pm10`, `povalley` |
| Model object | `STEM_Model`*, `STEM_Skeleton`, `STEM_Data`, `is.STEM_*` |
| Estimation engine | `STEM_Estimation`*, `kalman`, `filtering`, `filterstep`, `smoothing`, `smootherstep`, `smootherstep.uni`, `Q_function_addendo1/2/3`, `d1_Q`, `d2_Q`, `d12_Q`, `d1/d2_Sigmastar_logb.exp`, `d1/d2_Sigmastar_logtheta.exp`, `Sigmastar.exp`, `B_function`, `cov_lagone`, `sumMatrices`, `changedimension_covariates` |
| Simulation and prediction | `STEM_Simulation`*, `STEM_Kriging`*, `spatial.pred`, `STEM_Bootstrap`*, `STEM_Bootstrap.fn` |
| Spatially-clustered STEM | `SCSTEM_Estim`*, `SCSTEM_Infocrit`*, `SCSTEM_Select`*, `SCSTEM_Bootstrap`*, `SCSTEM_BootInference`*, and the `scstem_*` helpers |

## Edges

```
STEM_Model               -> STEM_Skeleton
STEM_Model               -> STEM_Data

STEM_Estimation          -> kalman
STEM_Estimation          -> changedimension_covariates
kalman                   -> filtering
kalman                   -> smoothing
kalman                   -> Q_function_addendo1
kalman                   -> Q_function_addendo2
kalman                   -> Q_function_addendo3
kalman                   -> d1_Q
kalman                   -> d2_Q
kalman                   -> d12_Q
kalman                   -> d1_Sigmastar_logb.exp
kalman                   -> d2_Sigmastar_logb.exp
kalman                   -> d1_Sigmastar_logtheta.exp
kalman                   -> d2_Sigmastar_logtheta.exp
kalman                   -> B_function
kalman                   -> cov_lagone
kalman                   -> sumMatrices
filtering                -> filterstep
smoothing                -> smootherstep
smoothing                -> smootherstep.uni

STEM_Simulation          -> changedimension_covariates
STEM_Simulation          -> Sigmastar.exp
STEM_Kriging             -> changedimension_covariates
STEM_Kriging             -> spatial.pred
STEM_Bootstrap           -> STEM_Bootstrap.fn          [indirect, via lapply]
STEM_Bootstrap.fn        -> STEM_Simulation
STEM_Bootstrap.fn        -> STEM_Estimation

SCSTEM_Estim             -> STEM_Model
SCSTEM_Estim             -> STEM_Estimation
SCSTEM_Estim             -> scstem_neighbors
SCSTEM_Estim             -> scstem_covariate_means
SCSTEM_Estim             -> scstem_init
SCSTEM_Estim             -> scstem_repair_partition
SCSTEM_Estim             -> scstem_loglike_i
SCSTEM_Estim             -> scstem_potts_pairs
SCSTEM_Estim             -> scstem_swap_pass
SCSTEM_Estim             -> scstem_npar
SCSTEM_Estim             -> scstem_rows
SCSTEM_Estim             -> scstem_with_seed
scstem_init              -> scstem_repair_partition
scstem_swap_pass         -> scstem_pen_local

SCSTEM_Infocrit          -> SCSTEM_Estim
SCSTEM_Select            -> scstem_ari
SCSTEM_Bootstrap         -> STEM_Model
SCSTEM_Bootstrap         -> STEM_Simulation
SCSTEM_Bootstrap         -> SCSTEM_Estim               [indirect, via do.call]
SCSTEM_Bootstrap         -> scstem_with_seed
SCSTEM_BootInference     -> scstem_align_labels
SCSTEM_BootInference     -> scstem_ari
```

## Reading the graph

Three things are worth noticing.

1. **The SC-STEM layer is a client of the STEM layer.** `SCSTEM_Estim()` does
   not reimplement the estimation: it builds one `STEM_Model` per regime and
   calls `STEM_Estimation()` on it. Anything that improves the engine improves
   the clustered models for free - and, conversely, any fragility of the engine
   is amplified, because the clustered algorithm fits the model on hundreds of
   different subsets of locations.

2. **`STEM_Simulation()` is the hinge of both bootstraps.** The classical
   `STEM_Bootstrap()` regenerates the whole network from a single pooled fit;
   `SCSTEM_Bootstrap()` regenerates it regime by regime and then re-runs the
   clustering, which is what propagates the uncertainty of the partition.

3. **The selection layer is stateless.** `SCSTEM_Select()` never refits: it
   consumes the grid produced by `SCSTEM_Infocrit()`, which already stores every
   fitted model, so the selected configuration is returned without any further
   computation.
