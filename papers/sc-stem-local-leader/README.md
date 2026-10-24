# Latent spatial regimes and local price leadership: paper, data and code

This folder holds the paper

> P. Maranzano, *Latent Spatial Regimes in Spatio-Temporal State-Space Models,
> with an Application to Local Price Leadership in Fuel Markets*

with its supplementary material, the data of its application and the code
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
  README.md                    this file
  manuscript/                  the paper and its supplements (PDF)
    sc-stem-local-leader.pdf     the manuscript
    supplement-A-estimation.pdf  Supplement A: estimation and computation
    supplement-B-simulations.pdf Supplement B: the simulation study
    supplement-C-application.pdf Supplement C: the application
    code-supplement.pdf          the implementation in Stem
  data/                        the data of the application
    fuel_stations.csv            the 3,136 stations of the eleven cities
    fuel_prices.csv.gz           their Sunday prices, January 2021 to June 2026
    fuel_national_means.csv      the national mean price of every Sunday
  code/
    application/               the scripts of the application
      functions.R                auxiliary functions, sourced by the scripts
      01-data-and-eda.R          the data of the models and the exploratory analysis
      02-models.R                the 66 SC-STEM models, their bootstrap and the tests
      03-results.R               the tables, figures and numbers of the paper
    simulations/               the scripts of the simulation study
      functions.R                the generator, the Monte Carlo, the runner
      01-design.R                the design
      02-run.R                   the run
      03-results.R               the tables, figures and numbers of the paper
  results/
    application/               the compact results of the 66 models (02-models.R),
                               so that 03-results.R runs without refitting
    simulations/               the results of the simulation study (02-run.R)
  output/                      everything the scripts write (not versioned)
```

## Requirements

- R 4.2 or later, with the packages `data.table`, `geodist`, `urca` and
  `remotes` from CRAN.
- Stem at commit `09e9cdb` of its GitHub repository, the version the results
  were computed with. The scripts check the installed version and, when it
  differs, install this one with
  `remotes::install_github("PaoloMaranzano/Stem@09e9cdb2d13fdfc94d8bc74e62645ebbc0b16a54")`.
  Restart R before running a script in a session where another version of Stem
  is loaded. Only `01-data-and-eda.R` and `02-models.R` need Stem.
- For the simulations, Stem at commit `0f7b744`, and the packages `spdep` and
  `callr` (see "The simulation study").

Every script runs from a terminal (`Rscript <script>`, in its folder) or from
RStudio (open it and press *Source*). Options are given on the command line as
`--name=value`; the header of every script lists them.

## The data

The prices come from Osservaprezzi Carburanti, the platform on which every
Italian fuel station posts its prices, published by the Italian Ministry of
Enterprises and Made in Italy. The full file, one row per station and day
since 2016, is not included: the stations of the eleven metropolitan cities of
the application were selected and pre-treated from it (codes of the same site
linked, highway stations excluded, regular reporting and active price setting
required; see Section 5.1 of the paper and Section C1 of Supplement C), and
this folder contains the result.

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
cd papers/sc-stem-local-leader/code/application
Rscript 03-results.R
```

It reads `results/application` and writes the tables, the figures and the
numbers of the paper into `output/application/paper` in about a minute.

### Full path

```
cd papers/sc-stem-local-leader/code/application
Rscript 01-data-and-eda.R
Rscript 02-models.R --cores=7
Rscript 03-results.R --results=../../output/application/results
```

1. `01-data-and-eda.R` builds the relative prices, completes the missing weeks
   of the lagged series with a pooled STEM model of every city and fuel, pairs
   every station Y with its station X in the three settings, and runs the
   exploratory analysis. About one hour, most of it for the completion of
   Rome.
2. `02-models.R` fits the 66 models (three settings, eleven cities, two
   fuels): the grid of (K, phi) and the two-step rule on the BIC, the
   refit-with-clustering bootstrap of the selected fit and of the pooled model
   with B = 100, the Wald tests, the station-by-station benchmark, the
   residual diagnostics and the check of the local maxima. This is the long
   step: about 70 hours on a desktop computer with the bootstraps on seven
   cores. The largest models (Rome, Milan, Turin in settings 1 and 3) take
   several hours each and some GB of memory. Every model is saved as soon as it
   is done and a stopped run resumes; a subset of the models can be run with
   `--settings=`, `--cities=` and `--fuels=`, and several machines can share
   the models, since the seed of a model depends only on its city and fuel.
3. `03-results.R` as in the quick path, on the new results.

Each run of `02-models.R` reproduces the shipped results up to the numerical
precision of the platform.

### Where every result of the paper comes from

| output (in `output/application/`) | place in the paper |
|---|---|
| `eda/eda_stations_by_city.csv` | Section 5.1; Supplement C, Section C1 (table of the stations) |
| `eda/eda_window.csv`, `eda/eda_missing_weeks.csv` | Section 5.1; Supplement C, Sections C1 and C2 |
| `eda/eda_unit_roots_and_stickiness.csv` | Section 5.1 (ADF test); Supplement C, Section C1 (ADF, autocorrelation, unchanged prices) |
| `eda/eda_pairs.csv` | Section 5.1; Supplement C, Section C1 (stations Y, distances to X) |
| `eda/eda_pacf_share_outside_band.csv`, `eda/eda_ar_order.csv` | Section 5.2; Supplement C, Section C2 (the choice of the lags) |
| `eda/fig_app_centring.png`, `eda/fig_app_centring_stations.csv` | Appendix C, Figure C1 and its text |
| `paper/tab_app_main.tex` | Section 5.3, Table 1 |
| `paper/fig_app_leadership.png` | Section 5.3, Figure 3 |
| `paper/fig_app_maps.png` | Section 5.3, Figure 4 |
| `paper/fig_app_passthrough.png` | Section 5.3, Figure 5 |
| `paper/tab_app_models.tex` | Supplement C, Section C3, Table C1 |
| `paper/fig_app_map_s1_RM_gasoline.png`, `paper/fig_app_map_s3_NA_gasoline.png` | Supplement C, Section C4, Figure C1 |
| `paper/tab_app_regimes_s1.tex`, `_s2.tex`, `_s3.tex` | Supplement C, Section C4, Tables C2 to C4 |
| `paper/tab_app_joint.tex` | Supplement C, Section C4, Table C5 |
| `paper/tab_app_distance.tex` | Supplement C, Section C5, Table C6 |
| `paper/tab_app_pulses.tex` | Supplement C, Section C7, Table C7 |
| `paper/app_numbers.csv` | every other number quoted in Sections 5.3 and 5.4 and in Sections C3 to C8 of Supplement C, one row each, with its place |

### The compact results

`results/application/setting<s>/<city>_<fuel>/` holds, for each of the 66
models, what `02-models.R` writes: `model.csv` (the selected number of regimes
and penalty, the sizes of the regimes, the usable bootstrap draws, the
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
simulations/
  functions.R      the generator, the Monte Carlo of one replication, the runner
  01-design.R      the design: scenario-variants, margins, blocks, parameter values
  02-run.R         the run: every cell and replication, on several cores
  03-results.R     the tables, figures and numbers of Section 4, Appendix B and
                   Supplement B
```

The study ran on Stem at commit `0f7b744`
(`remotes::install_github("PaoloMaranzano/Stem@0f7b74479b03a7fe5508eafe68fc7510c21422ad")`),
which `02-run.R` installs when another version is installed. 332 cells, 200
replications each, in two streams of 100 run on two machines:

```
cd papers/sc-stem-local-leader/code/simulations
Rscript 01-design.R                                 # prints the design
Rscript 02-run.R --stream=1 --nrep=100 --cores=7    # machine 1
Rscript 02-run.R --stream=2 --nrep=100 --cores=7    # machine 2
Rscript 03-results.R --results=../../output/simulations/results
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
observations that `02-run.R` also writes are left out: the analysis does not
use them, and any replication can be drawn again from its seed. A few lines
broken by the synchronization client while the runs appended to the files
were dropped, and 17 replications failed for reasons unrelated to the
estimator (Section B9 of Supplement B).

`03-results.R` without options reads `results/simulations` and writes into
`output/simulations/paper` in about half an hour; its header lists every
output with its place in the paper.
