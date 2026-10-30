# Internal material of the replication package

The replication package of the paper "Latent Spatial Regimes in
Spatio-Temporal State-Space Models: Evidence on Local Price Leadership in
Italian Fuel Markets" is `papers/sc-stem-local-leader` (data, code, results,
manuscript). This folder holds what stays outside it.

- `AppFuels_Internal_EDA.R`: the exploratory analysis of the weekly data of
  the application (Sections 5.1 and 5.2 of the paper, Figure C1 of Appendix C,
  Sections C1 and C2 of Supplement C). It sources `AppFuels_00_Setup.R` of the
  package, builds the data of the three settings as the run scripts do, and
  writes its tables and figure into `results/application/eda` of the package,
  where `AppFuels_04_Analysis.R` finds them. Run it again only if the weekly
  data change: `Rscript AppFuels_Internal_EDA.R` from this folder.

The weekly dataset of the package (`data/`) was built from the daily file of
Osservaprezzi Carburanti by the original pipeline, archived with its other
scripts on the Drive, in `SC-STEM-replication/_archive/01_original-pipeline`
(see its `README-archive-repo.md`).

On the Drive, `dev/sync-gdrive.R` delivers the package to
`STEM_Cameletti/SC-STEM-replication`, with the structure it has in the
repository (`README.md`, `manuscript/`, `data/`, `code/`, `results/`, and
`output/` for what the scripts write when they run there), and this folder to
its subfolder `internal/`, beside the daily file of Osservaprezzi Carburanti
(`internal/App_FuelsITA/`). The scripts of the original pipeline, the outputs
of its 66 models (`fuels/`), the raw results of the simulation study
(`results/`, from which `results/simulations` of the package was compacted)
and the superseded outputs are in `_archive/01_original-pipeline`;
`_archive/README.md` indexes the rest of the archive.
