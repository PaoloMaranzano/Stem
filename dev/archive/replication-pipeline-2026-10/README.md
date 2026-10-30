# The original pipeline of the paper (archived 2026-10-30)

These are the scripts that produced the results of the paper "Latent Spatial
Regimes in Spatio-Temporal State-Space Models: Evidence on Local Price
Leadership in Italian Fuel Markets", as they were run on the two machines in
September and October 2026. They were condensed, with the code unchanged, into
the replication package of the paper, `papers/sc-stem-local-leader/code`, and
are kept here for the record; they are not maintained.

| archived | in the replication package |
|---|---|
| `run-simulations.R` | `Sims_00_Setup.R` (generator, Monte Carlo, design) and `Sims_01_Run.R` |
| `design-figures.R`, `analyse-simulations.R` | `Sims_02_Analysis.R` (parts B and A) |
| `sims-make-replication-results.R` | the compact results in `results/simulations` (it wrote them from the raw results on the Drive) |
| `fuels-data.R` | stage 1 (the pretreatment of the daily file of Osservaprezzi Carburanti) is upstream of the package: `fuels-make-replication-data.R` wrote its result as the weekly dataset `data/`; stages 2 to 5 (relative prices, completion, pairs, data of the settings) are `fu_setting_data()` of `AppFuels_00_Setup.R`, with the nearest X only; the exploratory analysis is `dev/replication/AppFuels_Internal_EDA.R` |
| `fuels-functions.R` | `AppFuels_00_Setup.R` |
| `fuels-main-setting1.R`, `-setting2.R`, `-setting3.R`, `fuels-run.R`, `fuels-wald-bootstrap.R`, `fuels-check-optima.R` | `AppFuels_01_RunSetting1.R` to `AppFuels_03_RunSetting3.R` |
| `fuels-make-replication-results.R` | the compact results in `results/application` (it converted the outputs of these scripts to the format of the run scripts) |
| `fuels-analysis.R` | `AppFuels_04_Analysis.R` |
| `run-application.R` | the earlier air-quality application, removed from the paper |
| `run-ridge-experiment.R` | an experiment on the ridge in SC-STEM, not part of the paper |

`README-pipeline.md` is the README these scripts had. The scripts expect to
sit together and read and write `fuels/`, `results/` and `output/` beside
them; those folders (1.4 GB, not versioned) were moved here with them.
`fuels-functions.R` and `fuels/` are still used by the private bootstrap of
the likelihood-ratio tests, `dev/lr-bootstrap/fuels-lr-bootstrap.R`.
