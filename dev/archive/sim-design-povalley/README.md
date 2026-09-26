# Archived: the Po Valley simulation design

The simulation design of the SC-STEM paper as it stood at commit `0bdf027`
(2026-09-26), before it was made generic. Kept for possible later use; it is
not part of the replication material and nothing here is run by the current
scripts.

## What distinguished it

The baseline regime was calibrated on the pooled STEM fit of the Po Valley
PM2.5 network of the application, and the abstract plane of the overlap design
was mapped onto a geographic box around the Po Valley:

| quantity | value |
|---|---|
| coefficients | beta = (2.34, 0.60), one standardised covariate |
| nugget, partial sill | sigma2eps = 21.4, sigma2omega = 4.19 (nugget share 0.84) |
| residual sd | 5.06 |
| range | 1/theta = 123 km, exponential covariance on geodesic distances |
| latent process | G = 0.90, stationary variance 6.0 |
| geography | 100 km per abstract unit, centred at lon 9.5, lat 45.5 |
| covariate | AR(1) with a = 0.7, spatial innovations exp(-h / 100 km) |

The signal-to-noise ratio was therefore that of the PM data,
(0.36 + 6) / 25.6, about 0.25, with a nugget share of 84%.

Contrasts of the extreme regime (the middle one halfway), by level 1 / 2 / 3:

| scenario | contrast |
|---|---|
| S1 | beta_1 + 0.25 / 0.50 / 1.00 residual sd |
| S2 | G = 0.80 / 0.60 / 0.30 against 0.90 |
| S3 | nugget share 0.70 / 0.50 / 0.30 against 0.84, range 60 / 30 / 15 km against 123 |
| S4 | all three at level 2 |
| S5 | S4 with rho = 0 and the error field by regime |

The design of the overlap (omega, total spatial variance held fixed), the
blocks, the factors and the reference cell were the same as in the generic
design that replaced it; only the parameter values, the geography and the
distance changed. The generic design adds the scenario S0b.

## Files

| file | content |
|---|---|
| `run-simulations.R` | the runner, as committed; standalone, installs Stem from GitHub |
| `analyse-simulations.R` | its analysis script; it sources the runner from this folder |
| `paper-main-section3.tex` | Section 3 of the paper (simulation experiments) for this design |
| `paper-supplement-S5.tex` | Section S5 of the supplement (the data-generating process) |
| `tab_dgp.tex`, `tab_overlap.tex` | the scenario and overlap tables for this design |

To run it, open `run-simulations.R` in RStudio and press Source, exactly as
for the current runner. Note that the Stem installed from GitHub is the current
one, so the estimator is today's even though the design is the old one.
