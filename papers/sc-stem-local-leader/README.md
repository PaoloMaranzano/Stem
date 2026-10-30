# Latent spatial regimes and local price leadership: paper, data and code

This folder holds the paper

> P. Maranzano, *Latent Spatial Regimes in Spatio-Temporal State-Space Models:
> Evidence on Local Price Leadership in Italian Fuel Markets*

with its supplementary material, the dataset of its application and the code
that reproduces the application to Italian fuel prices (Section 5,
Appendix C and Supplement C) and the simulation study (Section 4, Appendix B
and Supplement B).

The methods are implemented in the R package **Stem**
(<https://github.com/PaoloMaranzano/Stem>). The scripts here use Stem as an
ordinary dependency; this folder is not part of the package and is excluded
from its builds.

## Contents

```
papers/sc-stem-local-leader/
  README.md                       this file
  manuscript/                     the paper and its supplements (PDF)
    sc-stem-local-leader.pdf        the manuscript
    supplement-A-estimation.pdf     Supplement A: estimation and computation
    supplement-B-simulations.pdf    Supplement B: the simulation study
    supplement-C-application.pdf    Supplement C: the application
    code-supplement.pdf             the implementation in Stem
  data/                           the dataset of the application
    fuel_stations.csv               the 3,136 stations of the eleven cities
    fuel_prices.csv.gz              their Sunday prices, January 2021 to June 2026
    fuel_national_means.csv         the national mean price of every Sunday
  code/
    AppFuels_00_Setup.R             application: the settings and the functions
                                    the other AppFuels scripts source
    AppFuels_01_RunSetting1.R       application: the 22 models of setting 1
    AppFuels_02_RunSetting2.R       application: the 22 models of setting 2
    AppFuels_03_RunSetting3.R       application: the 22 models of setting 3
    AppFuels_04_Analysis.R          application: every table, figure and number
    Sims_00_Setup.R                 simulations: the generator, the Monte Carlo,
                                    the runner and the design
    Sims_01_Run.R                   simulations: the run
    Sims_02_Analysis.R              simulations: every table, figure and number
  results/                        the results of the runs of the paper (read by
                                  the analysis scripts)
    application/                    the compact results of the 66 models, and
                                    (eda/) the outputs of the exploratory
                                    analysis of the data
    simulations/                    the results of the 66,400 replications
  output/                         what the scripts write (not versioned)
    application/                    the tables, figures and numbers of the
                                    application (AppFuels_04_Analysis.R)
    simulations/                    the tables, figures and numbers of the
                                    simulations (Sims_02_Analysis.R)
    runs/                           the results of a new run of the models or of
                                    the simulations, with their caches
```

`results/` and `output/` differ in who made them. `results/` holds what the
long runs of the paper produced (some 150 hours for the 66 models, 2,170
core-hours for the simulations), shipped so that the analysis scripts rebuild
every table and figure in minutes. `output/` is empty in the repository and
fills up as the scripts run here: the tables and figures of the two analysis
scripts, and, only if the models or the simulations are run again, their new
results in `output/runs/`, which the analysis scripts read with `--results=`.

Every script opens with the title of the paper, its author, what the script
does and the outputs it writes; in the analysis scripts every section names
the table or figure of the paper it produces.

## Requirements

- R 4.2 or later, with the packages `data.table` and `geodist` from CRAN (and
  `remotes`, to install Stem).
- For the application, Stem at commit `09e9cdb` of its GitHub repository, the
  version the results were computed with. The run scripts check the installed
  version and, when it differs, install this one with
  `remotes::install_github("PaoloMaranzano/Stem@09e9cdb2d13fdfc94d8bc74e62645ebbc0b16a54")`.
  Restart R before running a script in a session where another version of Stem
  is loaded. `AppFuels_04_Analysis.R` does not need Stem.
- For the simulations, Stem at commit `0f7b744`, and the packages `spdep` and
  `callr` (see "The simulation study").

Every script runs from a terminal (`Rscript <script>`, in `code/`) or from
RStudio (open it and press *Source*). Options are given on the command line as
`--name=value`; the header of every script lists them.

The scripts were checked on R 4.5.1 and R 4.6.1 (Windows): the two analysis
scripts, a model of the application run from the weekly data (completion of
the missing weeks included) and a replication of the simulation study give the
same tables, numbers and results on both, and the shipped results to within
$10^{-11}$. The figures differ between the two versions of R only in where the
text is placed, by a few pixels, since the two versions render fonts slightly
differently.

## The data

The prices come from Osservaprezzi Carburanti, the platform on which every
Italian fuel station posts its prices, published by the Italian Ministry of
Enterprises and Made in Italy as open data
(<https://www.mimit.gov.it/it/open-data/elenco-dataset/carburanti-prezzi-praticati-e-anagrafica-degli-impianti>,
license IODL 2.0). The dataset of the paper is the weekly one in `data/`: the
stations of the eleven metropolitan cities of the application, selected and
pre-treated from the daily file of the platform (codes of the same site
linked, highway stations excluded, regular reporting and active price setting
required; see Section 5.1 of the paper and Section C1 of Supplement C), and
their Sunday prices.

`fuel_stations.csv`: one row per station, in the order the models use.

| column | content |
|---|---|
| `city` | code of the metropolitan city: RM Rome, MI Milan, NA Naples, TO Turin, PA Palermo, BA Bari, CT Catania, BO Bologna, ME Messina, FI Florence, VE Venice |
| `site` | identifier of the station |
| `municipality`, `province`, `area` | its municipality, province and geographic area (North, Center, South and islands) |
| `brand` | its brand at the end of the window; independent stations have brand "Pompe Bianche" |
| `group` | the brand for the five major brands (Agip Eni, Api-Ip, Esso, Q8, Tamoil), "independent", or "other brands" |
| `type` | the type of station: "Stradale" (street station) for all, highway stations being excluded |
| `lat`, `lon` | its coordinates (degrees, WGS84) |

`fuel_prices.csv.gz`: one row per station and Sunday, from 10 January 2021 to
28 June 2026 (286 Sundays); columns `city`, `site`, `week` (the date of the
Sunday), `gasoline` and `diesel` (self-service prices, euro per liter; empty
when the station did not report that day). The 52 Sundays of 2021 are a
presample, which provides the lags of the first weeks of the window
(January 2022 to June 2026, 234 weeks).

`fuel_national_means.csv`: one row per Sunday; `gasoline` and `diesel`, the mean
price over all the reporting non-highway stations of Italy, and `stations`,
how many reported gasoline.

The empty field marks a missing value. Read the files with
`data.table::fread(..., na.strings = "")` or `read.csv(..., na.strings = "")`:
the code of Naples is the string `NA`.

## The application

### Quick path: the tables and figures from the shipped results

```
cd papers/sc-stem-local-leader/code
Rscript AppFuels_04_Analysis.R
```

It reads `results/application` and writes the tables, the figures and the
numbers of the paper into `output/application` in about a minute.

### Full path

```
cd papers/sc-stem-local-leader/code
Rscript AppFuels_01_RunSetting1.R --cores=7
Rscript AppFuels_02_RunSetting2.R --cores=7
Rscript AppFuels_03_RunSetting3.R --cores=7
Rscript AppFuels_04_Analysis.R --results=../output/runs/application
```

1. Each run script starts from the weekly data of `data/` and builds the data
   of its setting: the relative prices (Sunday price minus the national mean,
   minus the mean of the station over the window), their missing weeks
   completed for the lags by a pooled STEM model of every city and fuel, the
   pairs of every station Y with its station X, and the pulses of the excise
   duties. The completion takes about an hour on a first run, most of it for
   Rome; it is cached and the three settings share it.
2. Then it fits the 22 models of its setting (eleven cities, two fuels): the
   grid of (K, phi) and the two-step rule on the BIC, the
   refit-with-clustering bootstrap of the selected fit and of the pooled model
   with B = 100, the Wald tests, the station-by-station benchmark, the
   residual diagnostics and the check of the local maxima. This is the long
   step: in the run of the paper, settings 1 and 3 took about 73 hours each
   with the bootstraps on seven cores (Rome alone about a day), and setting 2
   about 5 hours. The largest models (Rome, Milan, Turin in settings 1 and 3)
   need some GB of memory. Every model is saved as soon as it is done and a
   stopped run resumes; a subset of the models can be run with `--cities=`
   and `--fuels=`, and the three settings can run on different machines,
   since the seed of a model depends only on its city and fuel.
3. `AppFuels_04_Analysis.R` as in the quick path, on the new results.

Each run reproduces the shipped results up to the numerical precision of the
platform.

### Where every result of the paper comes from

All the outputs below are written by `AppFuels_04_Analysis.R` into
`output/application`; the section of the script that writes each is
given in its header.

| output | place in the paper |
|---|---|
| `tab_app_main.tex` | Section 5.3, Table 1 |
| `fig_app_leadership.png` | Section 5.3, Figure 3 |
| `fig_app_maps.png` | Section 5.3, Figure 4 |
| `fig_app_passthrough.png` | Section 5.3, Figure 5 |
| `tab_app_models.tex` | Supplement C, Section C3, Table C1 |
| `fig_app_map_s1_RM_gasoline.png`, `fig_app_map_s3_NA_gasoline.png` | Supplement C, Section C4, Figure C1 |
| `tab_app_regimes_s1.tex`, `_s2.tex`, `_s3.tex` | Supplement C, Section C4, Tables C2 to C4 |
| `tab_app_joint.tex` | Supplement C, Section C4, Table C5 |
| `tab_app_distance.tex` | Supplement C, Section C5, Table C6 |
| `tab_app_pulses.tex` | Supplement C, Section C7, Table C7 |
| `app_numbers.csv` | every other number quoted in Sections 5.3 and 5.4 and in Sections C3 to C8 of Supplement C, one row each, with its place |
| `eda/eda_stations_by_city.csv` | Section 5.1; Supplement C, Section C1 (table of the stations) |
| `eda/eda_window.csv`, `eda/eda_missing_weeks.csv` | Section 5.1; Supplement C, Sections C1 and C2 |
| `eda/eda_unit_roots_and_stickiness.csv` | Section 5.1 (ADF test); Supplement C, Section C1 (ADF, autocorrelation, unchanged prices) |
| `eda/eda_pairs.csv` | Section 5.1; Supplement C, Section C1 (stations Y, distances to X) |
| `eda/eda_pacf_share_outside_band.csv`, `eda/eda_ar_order.csv` | Section 5.2; Supplement C, Section C2 (the choice of the lags) |
| `eda/fig_app_centring.png`, `eda/fig_app_centring_stations.csv` | Appendix C, Figure C1 and its text |

The `eda/` outputs describe the weekly data before any model is fitted. They
come from the exploratory analysis of the data, which is not part of this
folder; they are stored in `results/application/eda` and
`AppFuels_04_Analysis.R` copies them among its outputs.

### The compact results

`results/application/setting<s>/<city>_<fuel>/` holds, for each of the 66
models, what the run scripts write: `model.csv` (the selected number of
regimes and penalty, the sizes of the regimes, the usable bootstrap draws, the
stability of the partition, the minutes taken), `grid.csv` (the information
criteria over the grid), `regimes.csv` (the stations Y, their station X and
their regime), `estimates.csv` (the estimates of the selected and of the
pooled fit, with bootstrap standard errors and intervals), `draws.csv.gz` (the
aligned bootstrap draws of the regression coefficients of the selected fit),
`wald.csv` (the bootstrap Wald tests), `pump_tests.csv` and
`dumitrescu_hurlin.csv` (the station-by-station benchmark),
`diagnostics.csv` (the residual diagnostics) and `optima.csv`,
`optima_coefficients.csv` (the check of the local maxima).
`results/application/setting<s>/models.csv` stacks the rows of `model.csv`.

## The simulation study

```
code/
  Sims_00_Setup.R      the generator, the Monte Carlo of one replication, the
                       runner, and the design: scenario-variants, margins,
                       blocks, parameter values
  Sims_01_Run.R        the run: every cell and replication, on several cores
  Sims_02_Analysis.R   the tables, figures and numbers of Section 4, Appendix B
                       and Supplement B
```

The study ran on Stem at commit `0f7b744`
(`remotes::install_github("PaoloMaranzano/Stem@0f7b74479b03a7fe5508eafe68fc7510c21422ad")`),
which `Sims_01_Run.R` installs when another version is installed. 332 cells,
200 replications each, in two streams of 100 run on two machines:

```
cd papers/sc-stem-local-leader/code
Rscript Sims_00_Setup.R                               # prints the design
Rscript Sims_01_Run.R --stream=1 --nrep=100 --cores=7 # machine 1
Rscript Sims_01_Run.R --stream=2 --nrep=100 --cores=7 # machine 2
Rscript Sims_02_Analysis.R --results=../output/runs/simulations
```

The seed of a replication depends only on its stream and number, so the
results do not depend on the cores or on the machine; a stopped run resumes
where it stopped. The 66,400 replications took about 2,170 core-hours, that
is, about a week on two desktop computers.

`results/simulations` holds the results of the two streams, compressed:
`main3_s<stream>.csv.gz` (one row per replication: the selected number of
regimes and penalty, the recovery of the partition, the error of the fitted
signal, the times), `main3_s<stream>-params.csv.gz` (one row per replication,
regime and parameter: the truth and the estimate at the true number of
regimes) and `main3_s<stream>-grid.csv.gz` (one row per replication and
configuration of the grid). The per-location records and the full
observations that `Sims_01_Run.R` also writes are left out: the analysis does
not use them, and any replication can be drawn again from its seed. A few
lines broken by the synchronization client while the runs appended to the
files were dropped, and 17 replications failed for reasons unrelated to the
estimator (Section B6 of Supplement B).

`Sims_02_Analysis.R` without options reads `results/simulations` and writes
into `output/simulations` in about an hour and a half on a desktop
computer, most of it for re-applying the selection rule with seven criteria
to every replication (part R9 of the script); its header lists every output
with its place in the paper and the part of the script that writes it.
