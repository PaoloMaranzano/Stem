# internal/: the exploratory analysis of the data

Nothing in this folder is needed to reproduce the models or the tables and
figures of the paper: `code/AppFuels_04_Analysis.R` reads the stored outputs
in `results/application/eda`. This folder holds the script that computed them.

- `AppFuels_Internal_EDA.R`: the exploratory analysis of the weekly data of
  the application (Sections 5.1 and 5.2 of the paper, Figure C1 of Appendix C,
  Sections C1 and C2 of Supplement C): the stations by city, the window and the
  missing weeks, unit roots, persistence and stickiness of the prices, the
  pairs (Y, X) and their distances, the memory of the relative prices (which
  led to four lags for every block), and the figure of the two steps of the
  centering. It sources `code/AppFuels_00_Setup.R`, builds the data of the
  three settings as the run scripts do, and writes its tables and figure into
  `results/application/eda`. Run it again only if the weekly data change:

  ```
  cd papers/sc-stem-local-leader/internal
  Rscript AppFuels_Internal_EDA.R
  ```

  It needs, besides the requirements of the application, the package `urca`.
  The completion of the missing weeks is cached in
  `output/runs/application/cache/completion`, which the run scripts share.

The weekly dataset (`data/`) was built from the daily file of Osservaprezzi
Carburanti, which is not distributed here, by a pre-treatment described in
Section 5.1 of the paper and Section C1 of Supplement C.
