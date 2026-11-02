# Changelog

Development log for the `Stem` package: what changed, when, and why.

Commit messages in this repository are deliberately terse (`Update <date>`), so
that the file listing on GitHub stays readable. The reasoning behind each change
lives here instead. Entries are newest first.

For the user-facing summary of a release, see [`NEWS.md`](NEWS.md). For the
detailed account of the defects fixed in 2.0.0, see
[`dev/code-changes-report.html`](dev/code-changes-report.html). For the mapping
between the reference papers and the code, see
`inst/extdata/STEM_model_verification.html`.

---

## Unreleased

### 2026-11-02

**Diagrams in the documentation.** Flowcharts in Mermaid, which GitHub draws
from the text of the page, with no image files: the workflow of the functions
(README); how the computations nest and Algorithm 1 of the paper, one
iteration and the four algorithms, the cycle of SQUAREM, the grid and the
two-step rule (`docs/computational-aspects.md`, with the bootstrap of
2026-10-30); the layers of the estimation from `STEM_Fit()` to `kalman()`
(`docs/function-map.md`); the structure of the STEM model (`docs/models.md`);
the scripts, data and results of the replication package
(`papers/sc-stem-local-leader/README.md`). They are not shipped to CRAN:
`docs/` and `papers/` are outside the build.

**The paper.** At the suggestion of a colleague close to the journal, the
manuscript now places the application in the fuel industry and the data in
industry and data science (abstract, introduction, application, conclusion),
without changing its content; the clean PDFs are rebuilt.

### 2026-10-30 (night)

**ECME and the bootstrap explained step by step.** `docs/computational-aspects.md`
now builds ECME as it builds SQUAREM: why EM is slow on the intercept and on the
initial level of the latent process, why on the observed likelihood the update
is a regression, how the augmented Kalman filter computes the GLS step without
the $dT \times dT$ covariance, the update, the order that keeps it monotone,
and one iteration line by line with the lines that differ from EM marked. The
section on the refit-with-clustering bootstrap gains a flowchart (Mermaid,
drawn by GitHub), the pseudocode of the draws and of the alignment, and every
output with its formula: bootstrap mean, bias, standard error, the normal,
basic, percentile and bias-corrected intervals, the percentile tests between
regimes, the adjusted Rand index and the co-clustering matrix, the diagnostics
of the draws, and the software notes (seeds, cores, settings of the refits).
Supplement A has the same construction of ECME and the algorithm of the
bootstrap; the list of the outputs moved from Supplement A to this page, and
Supplement A no longer refers to the package (the model or the algorithm
instead), at the author's request.

**The exploratory analysis inside the replication package.**
`dev/replication` is now `papers/sc-stem-local-leader/internal/`, so the
package on GitHub has the structure of its copy on the Drive. The script finds
the package in the folder above its own; `AppFuels_04_Analysis.R` and the
README point to it; `dev/sync-gdrive.R` delivers it with the rest of the
package. The sentences on the version of Stem and on the seeds, which the
supplements no longer carry, are in the README of the package.

**The documentation aligned with the supplements of the paper.** The rule, at
the author's request: the supplements of the paper are a subset of what the
repository documents, with the same numbers. `docs/computational-aspects.md`
now follows Supplement A section by section (estimation, missing responses,
initialization and the scale of the penalty, degeneracy and the swap pass, the
label update and the monotonicity of the alternation, the information
criteria, the bootstrap and its tests, the Potts penalty on point-referenced
data, scope), with the implementation details, the numerical evidence and the
places in the code that the supplement leaves out; the design notes
(`docs/sc-stem-design.md`) are merged into it. As in the paper, the comparison
of the stopping rule with that of D-STEM is gone (from the page, the README,
the help of `STEM_control()` and NEWS), and the $0.01 I$ that earlier versions
added to the matrices they inverted is called a jitter, not a ridge; the
paragraph on SQUAREM under a ridge moved to `docs/regularization.md`. SQUAREM
is now explained step by step, for a reader who meets it for the first time
(the geometric error of EM, Aitken's extrapolation from three iterates, the
extrapolated point and why $\alpha = 1$ is EM, several directions, the scale,
the safeguards, one cycle line by line), and a synoptic table sets the four
algorithms side by side: what they share and where each departs from EM. The
same explanation and table are in Section A1 of Supplement A.
`inst/extdata/Stem-computational-supplement.pdf`, which still held the
supplement of 2026-10-04, is now the current Supplement A.

**`dev/` reduced to what the package and the paper need.** Kept: the sync
scripts, `make-logo.R`, `code-changes-report.html`, `regularization/` and
`replication/`. The archive of superseded material (`dev/archive/`), the paper
scripts of early October (`dev/paper/`), the internal design document of the
simulation study (`dev/simulation-design.tex`) and the private bootstrap of the
likelihood-ratio tests (`dev/lr-bootstrap/`) moved to the Drive, into
`SC-STEM-replication/_archive/` (indexed by its README), together with the
1.4 GB of outputs of the original pipeline that were never versioned. The
vignette `computational-notes` no longer points to `dev/paper/`.

**The replication package checked on R 4.6.1.** On R 4.6.1, with Stem
`09e9cdb` and `0f7b744` installed in temporary libraries:
`AppFuels_04_Analysis.R` and `Sims_02_Analysis.R` write the same tables and
csv files as on R 4.5.1, byte for byte; setting 2, Messina, gasoline, run from
the weekly data with its completion, matches the shipped results within 5e-12
(the completed series identical); the internal exploratory analysis gives the
same tables; one replication of the simulation study matches the stored one
exactly. The figures differ from those of R 4.5.1 only in the placement of the
text, by a few pixels (1 to 4% of the pixel values), since the two versions
render fonts slightly differently. The README of the package says so.

### 2026-10-30

**The replication package of the paper in eight numbered scripts.** At the
author's request, `papers/sc-stem-local-leader/code` is now one folder of
scripts named after the part of the paper they reproduce, with the code
unchanged: `AppFuels_00_Setup.R` (the settings and the functions of the
application), `AppFuels_01_RunSetting1.R` to `AppFuels_03_RunSetting3.R` (the
22 models of each setting), `AppFuels_04_Analysis.R` (every table, figure and
number of Section 5 and Supplement C), `Sims_00_Setup.R` (the generator, the
Monte Carlo, the runner and the design), `Sims_01_Run.R` and
`Sims_02_Analysis.R`. Every script opens with the title of the paper, its
author, what it does and what it writes, and every section of the two analysis
scripts names the table or figure it produces. The run scripts start from the
weekly dataset of `data/`, the dataset of the paper: `fu_setting_data()` builds
the data of a setting (relative prices, completion of the missing weeks by a
pooled STEM model, cached and shared by the settings, pairs, pulses), which the
old `01-data-and-eda.R` wrote to `cases.rds`. The exploratory analysis of the
data (Sections 5.1 and 5.2, Figure C1 of Appendix C, Sections C1 and C2 of
Supplement C) is no longer part of the package: it is
`dev/replication/AppFuels_Internal_EDA.R`, and its outputs are stored in
`results/application/eda`, from where `AppFuels_04_Analysis.R` copies them.
Checked against the old scripts: the data of the three settings are identical
(`identical()`); the stored exploratory outputs are identical byte for byte;
the analysis writes the same 28 files byte for byte; setting 2, Messina,
gasoline, run again from the weekly data matches the shipped results within
5e-12; the design printed by `Sims_00_Setup.R` is identical; one replication
of the simulation study run by `Sims_01_Run.R` on Stem `0f7b744` matches the
stored one exactly; and `Sims_02_Analysis.R` writes the same 38 tables and csv files byte for byte, its figures differing from those of the old script only in the anti-aliasing (at most 679 of 28 million pixel values), as two runs of the same script do. The scripts of the original pipeline that
were in `dev/replication`, with their outputs, are archived in
`dev/archive/replication-pipeline-2026-10` (its README maps them to the new
scripts); `dev/lr-bootstrap/fuels-lr-bootstrap.R` finds its functions there.
`results/` keeps the results of the runs of the paper, which the analysis
scripts read; `output/` (not versioned) receives the tables and figures of the
two analysis scripts directly in `output/application` and `output/simulations`,
and the results of a new run of the models or of the simulations in
`output/runs`. `dev/sync-gdrive.R` now delivers the package itself to
`SC-STEM-replication` on the Drive, with the structure it has here, and
`dev/replication` to its subfolder `internal/`; the original pipeline and the
raw outputs of the runs are in `_archive/` there.

### 2026-10-29

**A vignette for the numerical experiments on the estimation algorithms.** The
experiments that support the estimation algorithms and their defaults (the
checks of correctness and exactness, the four algorithms at three tolerances,
the accuracy against the true values, the time of whole grids, the effect of
the limits of the range) were Section A1.13 of the paper's Computational and
Algorithmic Supplement. At the author's request they now live in the package
documentation, as `vignettes/estimation-experiments.Rmd` ("5. Numerical
experiments on the estimation algorithms"), text and tables unchanged; the
supplement points to it. Two more parts of that section of the supplement
moved to the vignettes: the account of the jitter of $0.01 I$ that the earlier
versions added to the matrices they inverted, with its table (Section 8 of the
same vignette), and the organization of the iterations, with the functions and
the settings of `STEM_control()` that run them (vignette "3. Map of the
package", section "2. The estimation engine"). `dev/sync-paper.R` builds Supplement A from
`supplement-computational.tex`, its new name in the Overleaf project.

### 2026-10-27

**The regimes of the application against a purely spatial partition.** The
paper says that the regimes of the fuel application are groups of brands
rather than districts, on the strength of an adjusted Rand index of 0.23
(median) between the regimes and the brands. A reader asked for a benchmark,
and `papers/sc-stem-local-leader/code/application/03-results.R` now gives one:
the same index between the regimes and the k-means partition of the
coordinates of the stations into as many groups (`fu_ari_space()` in
`functions.R`). It is -0.02 to 0.10 (median 0.01) and below the index with the
brands in all 44 models of settings 1 and 3. The new column "ARI space" of
Table C1 and two rows of `app_numbers.csv` carry it; the k-means restores the
random stream, so every other table and figure is unchanged byte for byte.
The script also reads the compressed draws through `readLines(gzfile())`:
`fread` on a `.gz` file took about 17 seconds per file on Windows, and the
script now runs in seconds again. The manuscript PDFs were rebuilt after the
text revision that followed (predictive leadership, the limits of the
evidence, the data source, the titles of the supplements).

### 2026-10-23

**The simulation study is complete: 200 replications per cell.** The second
stream ended (100 replications of the 332 cells, as the first), and the
analysis runs on the 66,383 replications that ended normally (17 failed for
reasons unrelated to the estimator: R processes that could not start). The
tables and figures of Section 4 and Supplement B were rebuilt with
`papers/sc-stem-local-leader/code/simulations/03-results.R`, and the numbers
of the text updated; the conclusions are unchanged. The numbers were computed
by one script on the outputs of the analysis at 142 and at 200 replications:
it reproduces the 142-replication values of the text, so the update uses the
definitions of the text. One number of the text was not reproducible even at
142 replications (the regression coefficients "up to three times" more
dispersed when K is wrong: the ratio of the interquartile ranges reaches 3.8)
and is corrected. The design tables and figures differ from those of the
paper only in spelling (center, neighbors).

`papers/sc-stem-local-leader/results/simulations` holds the results of the two
streams compressed (63 MB): the replications, the parameters and the grids,
cleaned of the lines broken by the synchronization client by
`dev/replication/sims-make-replication-results.R`; the per-location records
and the full observations are left out. `03-results.R` now reads `.csv.gz`
files; run on them it gives the same tables and the same 19 tables of numbers
as on the original files, and the same figures (four PNG files differ in their
bytes only, from the rendering of the mathematical labels).

### 2026-10-22

**The paper gets its own folder: `papers/sc-stem-local-leader`.** The
replication package below moves there (`git mv`, history kept), with the
scripts under `code/application` and `code/simulations`, beside `data/`,
`results/` and a new `manuscript/` that holds the PDFs of the paper and of its
supplements, never their sources (user's decision: Zenodo or Mendeley only at
publication; until then the repository of the package hosts the paper).
`dev/sync-paper.R` compiles the PDFs from the Overleaf project on a temporary
copy in which the notes are dropped and the changes marked in orange are
printed in black, and reports any text still in red. `papers/` replaces
`replication/` in `.Rbuildignore`, and the default folders of the scripts and
of the two `dev/replication/fuels-make-replication-*.R` follow the move. The
data availability and code availability statements of the manuscript now
point to the folder.

**The replication package of the paper.** A new top-level folder,
`replication/`, reproduces the application and the simulation study of the
paper with few, condensed scripts. It is excluded from the package build
(`.Rbuildignore`), and what its scripts write (`replication/output/`) is not
versioned. The development scripts in `dev/replication` stay as they are: the
package is the reader-facing version of them.

- `replication/data`: the minimal data of the application, the Sunday prices of
  the 3,136 stations of the eleven cities (2021 to June 2026), their registry
  and the national mean prices. They are written from the output of the
  pre-treatment by `dev/replication/fuels-make-replication-data.R`; the full
  station-level file stays out. The coordinates and the national means are
  written with 17 significant digits: with 15, the distances changed in the
  14th digit and the completed series in the 11th, enough to break the exact
  reproduction of the fits.
- `replication/application`: `functions.R`, `01-data-and-eda.R` (relative
  prices, Kalman completion, pairs, exploratory analysis), `02-models.R` (the
  66 models: grid and two-step rule, bootstrap of the selected and pooled fits,
  Wald tests, station-by-station benchmark, diagnostics, local maxima) and
  `03-results.R` (every table, figure and number of Section 5 and Supplement C,
  each labeled with its place in the paper). The radius r*, the ridge and the
  likelihood-ratio test of the development scripts are left out, as in the
  paper.
- `replication/results/application`: the compact results of the 66 models, in
  the format of `02-models.R`, converted from the original runs by
  `dev/replication/fuels-make-replication-results.R` (the Wald tests are
  recomputed and checked against the original ones; the pooled estimates are
  read from the cached grids at full precision, since the bootstrap summary
  rounds them to five digits).
- `replication/simulations`: `run-simulations.R`, `analyse-simulations.R` and
  `design-figures.R` split into `functions.R`, `01-design.R`, `02-run.R` and
  `03-results.R`, the code moved as it is (seeds unchanged); the worker
  processes now load `02-run.R`, told its path, and Stem is checked against the
  pinned commit without asking GitHub.

Checks: the data rebuild the relative prices, the completed series and the
data of the three settings of the original pipeline exactly; `02-models.R` on
one model (setting 2, Messina, gasoline) gives the same partition and
estimates as the original run, the grid within 5e-12 and the bootstrap draws
within 2e-11; `03-results.R` on the shipped results reproduces the eight
tables and the five figures of the paper byte for byte; `02-run.R` reproduces a
stored replication of the simulation study exactly. Its numbers showed seven
statements of the text to be corrected: the median adjusted Rand index between
regimes and brands (0.23); the partition reproduced in the median draw in all
models but four; the sign changes at the local maxima; Naples among the cities
where half or fewer of the independents follow; and three descriptive figures
that came from an older exploratory script on daily prices carried forward,
now computed on the relative prices of the models (the ADF test rejects for 96
to 100% of the stations, the median first-order autocorrelation is 0.52 to
0.71, 12% of the major-brand stations have an independent within one
kilometer).

### 2026-10-20

**The likelihood-ratio test leaves the replication material.** The bootstrap
likelihood-ratio test of the fuel application added little to the Wald test
(they agreed on 93% of the regimes and 97% of the joint tests in the first
32 models) and lengthened a paper that is already long. The paper and the
replication material therefore use the bootstrap Wald test only, and the
results of the likelihood-ratio run are kept for a possible request of a
referee. `fuels-lr-bootstrap.R` moves to `dev/lr-bootstrap/`, where it still
runs: it finds `fuels-functions.R` and the fuels folder in `dev/replication`.
The README row, the comparison table in the header of
`fuels-wald-bootstrap.R` and the mentions in the comments of
`fuels-functions.R` and of the three main scripts are removed. The functions
of section 10 of `fuels-functions.R` are unchanged, since the private script
uses them. The Drive copies of `fuels-functions.R` and `fuels-lr-bootstrap.R`
are not updated while the run on the other machine may still be using them.

**American spelling in the replication material.** Comments, figure labels
and table captions of the fuel scripts, of `design-figures.R` and of
`analyse-simulations.R` now use American spelling (neighbor, liter,
kilometer, centered, color, behavior), as the paper does. Function names
(`dgp_centres()`), the area label "Centre" of the pre-treatment data and
`run-simulations.R` (runs in progress) are unchanged. `fuels-analysis.R` was
rerun: the tables differ from the previous ones only in their captions.

**The centering figure picks stations further apart.** The figure of the two
steps from posted to relative prices (Rome gasoline) took an Eni station at
the third quartile of the premiums of its group and an independent one at the
first. In Rome the premiums are small relative to the weekly fluctuations,
and the second step was hard to see. The figure now takes the ninth and the
first decile (premiums +0.8 and -2.7 cents).

### 2026-10-17

**The likelihood-ratio test of the fuel application, on converged fits.** The
run of the three settings (K at most 5, lambda = 0) ended without errors, all
127 bootstraps with 100 usable draws out of 100. In one of the 66 models,
Rome gasoline of setting 1, the likelihood-ratio test on the selected partition
(`fu_lr_fixed()`) gave negative statistics in regime 1 (127 pumps): -335, -351
and -261 on the three blocks, so that the full model had a lower
log-likelihood than the restricted ones. The refits of the test started from
least squares and stopped by the default rule (`em_stop = "any"`), which ends
the EM algorithm as soon as the parameters move little, even while the
log-likelihood is still rising: the full model of that regime stopped short.
The refits of the test now stop only when both criteria are met
(`em_stop = "all"`, at most 2000 iterations), and the full model of a regime
is the better of its refit and of the final refit of the selected fit, which
has the same partition and covariates and starts from the alternation. The
output reports the log-likelihoods, the source of the full one and the
convergence of every fit. The bootstrap Wald tests, the main inference, were
not affected. The settings of the cached test carry a marker, so that a rerun
of stage final recomputes the test; only Rome gasoline of setting 1 was rerun
for now, the other 65 models keep the test computed with the default rule.

The rerun did not remove the negative statistics: in regime 1 all the fits
now converge under the stricter rule, and the full model still ends at
-57742 (from least squares; the selected fit at -57836) against -57567 to
-57611 for the restricted ones. The cause is not an early stop but a
likelihood with several local maxima in the regime: the EM algorithm reaches
different ones from different starts.

**How common the local maxima are (`fuels-check-optima.R`).** For the 66
models, the selected partition was refitted regime by regime from least
squares (em_stop = "all") and set beside the selected fit, whose regimes start
from the alternation; output `fuels/checks/optima_by_regime.csv`, 25 minutes on
6 processes. Of 268 regimes (the 5 pooled models have one start only), 191
agree within 1 unit of log-likelihood; 52 differ by more than 10 and 17 by
more than 100, every fit converged. Least squares reaches the higher maximum
in 38 regimes by more than 10, the selected fit in 14. 35 of the 61 models
with K* > 1 have at least one regime differing by more than 10 (15 of 22 in
setting 1, 6 of 17 in setting 2, 14 of 22 in setting 3); the share grows with
the size of the regime (11% of the regimes up to 10 pumps, 47% above 100).
The two maxima of a regime usually differ in the range of the spatial field
and the persistence of the latent process (Rome gasoline, setting 3, regime 3:
theta 5.7e-6 and G = 0.93 in the selected fit, theta 3.7e-3 and G = 0.19 from
least squares, 264 units higher). The final refit of a regime starts from one
point (`SCSTEM_Estimation`, Final refit); whether to change it is open.

**The coefficients do not depend on the maximum.** The check was extended to
the coefficients at the two maxima (with the bootstrap standard errors of the
selected fit) and to the bootstrap Wald tests at both (on the covariance of the
same draws); outputs `optima_coefficients.csv` and `optima_wald.csv`, 28
minutes. In the 52 regimes whose maxima differ by more than 10 units, the
lag coefficients move by a median of 0.13 to 0.17 bootstrap standard errors;
of the 416 neighbour coefficients, one moves by more than 2; no significant
coefficient changes sign. The two maxima differ in the variance parameters
(range of the spatial field, persistence of the latent process), not in the
regression. At 5%, the Wald test of the neighbour (both fuels) changes
decision in 4 of 268 regimes: in two the selected fit is the higher maximum,
so its decision stands; in the other two the selected fit is the lower one,
Turin diesel of setting 3 (p = 0.048 against 0.050) and Catania gasoline of
setting 2 (a regime of 5 pumps, p = 0.39 against 0.027). On this evidence the
estimates and the bootstrap tests can stand without a rerun (the user's
preference), the check serving as a robustness note; the decision is the
user's. The likelihood-ratio test on the fixed partition is the output the
multiple maxima do affect.

**A bootstrap likelihood-ratio test (`fuels-lr-bootstrap.R`).** To put the
likelihood-ratio test on the same footing as the bootstrap Wald test, which
propagates the uncertainty of the partition, its null distribution is now
obtained by a parametric bootstrap under the null: samples drawn regime by
regime from the restricted fits (the generator of `SCSTEM_Bootstrap`), SC-STEM
re-estimated on each at (K*, phi*), the statistic recomputed on the regimes of
the draw and matched to the observed ones by the majority rule; p-value
(1 + #{LR* >= LR}) / (1 + B). The two maximized likelihoods are made
consistent by nested starts: the restricted model from two starts (the full
estimates without the tested lags, least squares), the full model refitted
from the restricted solution with the tested lags at zero, so that the
statistic is never negative. One bootstrap per block imposes the null in all
the regimes. Only exported functions of Stem are used; by the user's
decision the test stays in the replication material, and the package is
left as it is (a joint test on a subset of coefficients is computed by the
user from the aligned draws of `SCSTEM_BootInference`). Tested on Bari and
Venice diesel, setting 2 (4 to 6 draws), including the resumption of a
stopped run.

**The Wald test in a script of the same form (`fuels-wald-bootstrap.R`).** At
the user's request the two tests are written in parallel: the same header
(with the table of what each takes from the bootstrap), the same setup, the
same sections and the same output columns (`test_wald.csv`, `test_lr.csv`,
and their summaries per setting); the machinery they share is section 10 of
`fuels-functions.R`. The Wald script reads the bootstrap of stage `final`
(same seed, same B) and reproduces its statistics exactly (Bari and Venice
diesel, setting 2); it adds the test in all the regimes at once (the
coefficients of the K regimes stacked, with their joint bootstrap
covariance) and a bootstrap p-value from the draws centred at the estimate,
the counterpart of the p-value of the LR test.

**The old likelihood-ratio test leaves the main scripts.** With the user's
approval, `fu_lr_fixed()` (the test on the fixed partition with the chi-square
reference) is removed from `fuels-functions.R` and from stage `final` of the
three main scripts, and the change of the morning (EM rule and start from the
selected fit) is undone with it: the bootstrap LR test of
`fuels-lr-bootstrap.R` replaces it. The 66 files `lr_fixed_partition.csv`
already written stay where they are.

**The bootstrap LR test failed on the other machine.** Run with Source in an
open RStudio session, every model stopped at once ("the restricted model of a
regime could not be fitted"), while the same script ran here from a terminal
and through `source()`. The fits of the test call `STEM_Model()` and
`STEM_Estimation()` in the session itself, and the likely cause is another
version of Stem loaded in that session: `fu_require_stem()` checked the
version installed, not the one loaded. It now also stops, with a message to
restart R, when the loaded Stem is not the pinned one; and a failed fit of the
LR test reports its error instead of a bare NULL.

**The analysis of the application (`fuels-analysis.R`) and its text.** A new
script reads the outputs of the 66 models (the Wald tests regime by regime,
the bootstraps, the regimes, the pulses, the pump-by-pump tests, the
diagnostics, the check of the local maxima) and writes the tables and figures
of Section 5 of the paper and of a new Supplementary Material C, the
Application Supplement. The long-run pass-through from the neighbour,
lambda_k = sum_l d_kl / (1 - sum_l b_kl), and its percentile interval come
from the aligned draws of the bootstrap of stage `final`. Findings, Wald only
(the LR test is running): the neighbour leads in 83% of the regimes of setting
3 (major after another major; 89-93% of the pumps, long-run pass-through median
0.15-0.17), 68% in setting 1 and 56% in setting 2; the regimes are pricing
groups of the brands rather than areas (phi* = 0 in 29 models, ARI with the
brands median 0.24, one Eni-only regime in 25 of the 44 models of settings 1
and 3); Eni-dominated regimes follow less (72% against 86% in setting 3); the
pooled model misses the leadership in 9 city-fuel-settings where most pumps
follow, and attributes it to the whole city in 5 where a minority does; the
2022 excise cut and its end pass through very unevenly across the regimes of a
city (spread about 25 cents in the week), the small changes of 2025-2026
uniformly. Paper: Section 5 rewritten (questions, data, model, results,
implications; three figures and one table), the abstract and the outline
mention the application, Section 3.3 points to the new Section A7.

**Supplement A, Section A7: tests of several coefficients under the
bootstrap.** The Wald test (the bootstrap gives the covariance; with the draws
centred at the estimate, also the null distribution) and the likelihood-ratio
test (the bootstrap of the restricted model gives the null distribution), with
the classical F test as the exact case they generalize, the nested fits that
make the LR statistic immune to local maxima, and a table setting the two side
by side. Simulation experiments and code move to A8 and A9; the code section
lists the scripts of the application.

### 2026-10-13

**At most five regimes in the fuel application (night).** With the cap of
three, the rule selected K* = 3, the cap itself, in every model fitted so far
(Rome gasoline, 167/248/136 pumps; Venice gasoline and diesel): a selection on
the boundary of the grid cannot tell whether four or five regimes would be
preferred, and on networks of hundreds of pumps three regimes may be too few.
The user raised the cap to five (`k_cap = 5`, the default of the three main
scripts and an option of `fuels-run.R`): K_max = min(5, floor(n / 5)), the
pooled model below ten pumps. The grid grows from 15 to 29 fits per model;
its cache is keyed on the K of the grid, so the grids computed with the cap
of three are refitted. The run launched on the other machine with the cap of
three was stopped by the user before any grid was saved. The dry plan now
counts the bootstrap at lambda* only when the grid of lambda has a positive
value. The over-selection that the BIC showed on wide networks in the
simulation study (up to 27% at d = 200 with weak overlapping regimes) makes
the EBIC a check worth computing on the saved grids.

**The ridge leaves the fuel application and the paper (night).** The questions
of the application are inferential, and the ridge corrected for its bias is the
estimate at lambda = 0 (Proposition 1(c) of the earlier draft), so it changed
only point estimates (Venice pilot: lambda* = 0.003, standard errors
unchanged); pumps are point installations, nothing is mapped at unmonitored
sites, and in the scratch checks the ridge improved the interpolation of
missing cells by less than 0.1%. The user decided to take it out of the paper
and keep it for a second paper on interpolation and extrapolation (air quality
over Europe, ERA5 covariates in complex terrain), whose pieces are in the
Overleaf folder "SC-STEM Ridge". `fuels-run.R` gains `ridge = FALSE` (the
default): it passes `--lambdas=0`, so that stage lambda reduces to one fit at
lambda = 0 (its cache, keyed on the grid of lambda, is recomputed), lambda* is
0, the fit at lambda* and its bootstrap are skipped, and no ridge table is
written; checked end to end on Venice gasoline with two draws (lambda stage
0.15 min, bootstraps at lambda = 0 and of the pooled model only). The
bootstrap at lambda* was one of the three bootstraps of stage final.

**Boxplots everywhere, K_0, the paper restructured (afternoon).** The user
prefers boxplots to intervals throughout, for coherence. The core selection
figure now shows the selection error K_hat - K_0 of every replication as
boxplots (a boxplot of the 0/1 outcome "K right" collapses; this one also shows
whether the errors are too few or too many regimes; the shares stay in the
tables and in the summary heatmap); the robustness figure shows K_hat - K_0 and
the ARI as boxplots, S2 out of the ARI panel; the penalty figure is one panel
of boxplots per scenario-variant. The true number of regimes is written K_0 in
every figure and in the text, the reference cell is spelled out wherever it is
used, the parameters figure has three columns, the summary figure a common
title and the four geometries named. The fitted-signal figure moved to the
main text. The Wilson helper is gone. The paper (Overleaf, not in this
repository) moved to the Springer Nature template in referee mode for ASMBI,
with the computational section before the simulations and two supplements, the
Computational and Algorithmic Supplement (A) and the Simulation Supplement (B).
In the evening the manuscript lost every reference to the package (only the
code-availability declaration names it; the package appendix became a Code
Supplement, not submitted, for the GitHub page and CRAN), the EM appendix moved
into Supplement A, and the internal design notes of the simulation study, whose
content now lives in Supplement B, moved from the Overleaf project to
`dev/simulation-design.tex` (their development-extension section stays out of
the paper).

**The simulation results, first reading (replication material only).** Stream
1 of main3 ended on 2026-10-12 (332 cells x 100 replications, 1,344
core-hours, three replications failed because R could not enter its working
directory, no timeout); stream 2 was at replication 42 of 100. The large
results files of stream 1, written on Google Drive, had lost blocks of
appended lines (7 broken lines in the grids, 8 in the parameters, 28 in the
stations; the small summary file intact). The runner writes from one process
only, one replication at a time, so the synchronization client is the cause:
it loses an appended block while it uploads a file of hundreds of megabytes.
The losses do not depend on the results. `analyse-simulations.R` now reads the
files defensively: a line that does not start with the quoted cell or has the
wrong number of fields (commas inside quotes excluded) is dropped, a last line
without its end of line too, and a replication with fewer rows than the others
of its cell is dropped from that file (15 in the grids, 10 in the parameters).
A stream still running counts the replications that all its cells have reached
(`--balanced`, default TRUE), and `--upto=100,42` caps them per stream, so
that the numbers stay fixed while the paper is written.

**Figures, at the user's request.** All the figures of the analysis and of the
design are PNG (300 dpi). The core figures show d as ordered categories, the
ARI as boxplots of the replications of a cell, the selection of K as a share
with its 95% Wilson interval (a boxplot of a 0/1 outcome is degenerate), and
the legend below the panels in a larger font. The robustness figure says what
the arms are (one factor changed at a time, the others at the reference cell,
named in a subtitle) and carries 95% intervals. New: R10, the parameters at the
true K split by whether the rule selected the true K in the same replication
(`fig_sim_params`, `tab_sim_params_by_selection.csv`), and the fitted signal of
the selected model against the pooled fit (`fig_sim_signal`, `tab_sim_signal.csv`);
the runner records the parameters at the true K only, so the parameters of a
selected fit with the wrong K are not available without refitting. R11, the
one summary figure of Section 3 (`fig_sim_main`): the share selecting the true
K by scenario-variant and (T, d) over the four geometries, the ARI and the
relative error of the parameters at the reference cell. `design-figures.R`
adds `fig_sim_overview`, the regimes in space (four geometries, unbalanced
sizes) and the series of the nine scenario-variants, for the appendix of the
paper, in the place of the old `fig_dgp_examples` of the first design.

**What the results say** (142 replications per cell). S2 is never split
(1,704 replications). The well separated scenario-variants are identified in
93-100% of the replications from d = 100 on (over the geometries), the
failures being at d = 40 with T <= 120; S1w-shr and S3beta-shr in 25% and 24%
at the reference cell, up to 88% and 67% at d = 200, T = 365. Separation in
space makes K harder to identify (70% total overlap against 33% strong
separation at d = 40, T = 60), while the ARI at the true K hardly changes; the
over-selections occur with overlapping regimes (never with strongly separated
weak regimes), up to 27% of the
replications of the weak scenario-variants at d = 200 and 12-20% at d = 400.
At the true K and with K selected correctly the parameters are centred (median
relative error within 0.3% for the coefficients, 4% for the rest); with K
wrong the coefficients stay centred but are more dispersed and the variance
components are distorted. The AIC selects the true K more often than the BIC
in the weak scenario-variants (32-33% against 26-28%) with more
over-selection. Section 3 of the paper is rewritten short, in bullets, around
this figure, the extended results go to Supplementary Material B, Section B9.

### 2026-10-10

**`SCSTEM_Bootstrap()` in parallel.** The pilot of the fuel application on
Venice (100 pumps) took 32-46 minutes for one bootstrap of 100 draws, 19-27 s
a draw; on Bari (200 pumps) a draw takes 33 s with one fit per draw, and the
cost grows faster than the number of locations, so a bootstrap of Rome (551
pumps) would run for many hours in one process. The draws are independent, and
the user asked for an argument `cores` (default 1) so that the user chooses how
many processes to use. With `cores > 1` a socket cluster (package parallel,
now in Imports; it works on every platform) shares the draws in blocks and is
stopped on exit; the processes load the installed Stem. To make the result
independent of the number of processes, every draw now runs under its own
seed, drawn from `seed` before the draws start (before, the draws shared one
stream, so a draw depended on those computed before it): the draws of a given
seed differ from those of the previous code, and are the same whatever
`cores`. Test: a bootstrap on two processes returns the draws of the one in
the session. R CMD check: Status OK. Script-level chunking of the bootstrap
was written first and dropped in favour of this.

**The fuel scripts on the parallel bootstrap.** `FU_STEM_REF` and
`RX_STEM_REF` pinned to 09e9cdb; the main scripts get `boot_cores` (default 1,
`--boot_cores=`). Between 989790b and 09e9cdb only `R/SCSTEM_Bootstrap.R`
changed, so the cached grids, penalty paths and likelihood-ratio fits of the
pilot were copied to `cache/09e9cdb/`; the bootstraps were not (every draw now
has its own seed). Measured on setting 1, gasoline, one draw in one process:
Venice (100 pumps, phi* = 0.025) 19 s, Bari (200, phi* = 0) 33 s, Rome (551,
phi* = 0.025) 168 s; grids 2.6, 5.5 and 81 minutes. Eight draws of Rome on
eight processes took 306 s, a speed-up of 4.4 (with one draw per process the
start of the processes and the transfer of the fitted model weigh most).
Pilot results on Venice (setting 1): K* = 3 for both fuels, lambda* = 0.003;
the neighbour leads Y in every regime (bootstrap Wald p 0.00003-0.007), while
the pooled model does not see it for gasoline (p = 0.48).

**The fuel scripts on several machines.** The user will run the application on
several machines, choosing the cities of each. `--cities` and `--fuels`
already did that; added `--stage=all` (grid, lambda and final in one process)
and `--stage=summary` (the tables of the three stages rebuilt from the files
of every model of the case in the folder, nothing fitted), and the tables of
every stage are now rebuilt from all the models present rather than from
those of the process. Tested on a synthetic case of two cities run in two
folders, merged, then summarized. The test also hit the 260-character limit
of Windows paths, which made a cached bootstrap fail to save after it had run:
the lag key in the names of the cached files is shorter ("L1t4-1t4-1t4-1t4",
runs of consecutive lags written "atb"; the pilot caches were renamed), and
`fu_cached()` stops before running a stage whose file path would exceed 250
characters.

**The fuel application in four stages.** Decided with the user after checks in
scratch: the regimes are selected without the ridge (grid of (K, phi) at
lambda = 0, two-step rule on the BIC), then lambda is chosen on that partition
held fixed by the AIC, not the BIC (checks on the DGP of the ridge
experiment: the BIC shrank too much and made the coefficients worse than the
MLE in most replications at moderate collinearity; AIC = AICc = GCV came close
to the best lambda). All tests are made at lambda = 0 (the bias-corrected
ridge is the MLE); the ridge gives the lag coefficients and the completed
response. The interpolation gain of the ridge was negligible (< 0.1%): the
paper will say "no loss", not "improvement". Main scripts: stages lags, grid,
lambda, final, dry; `fu_select()` replaced by `fu_lambda()`; `fu_ridge_table()`;
`fu_cached(valid =)` reruns a stage whose settings or input changed. A defect
caught before any run: the fit at lambda* is made with max_iter = 0, which its
control keeps, so its bootstrap would not have re-estimated the partition; the
draws now get the control of the fit at lambda = 0. Smoke-tested end to end on
synthetic data. `run-ridge-experiment.R` is being adapted and is blocked from
running meanwhile.

**The ridge experiment and the criteria for K and phi.** `run-ridge-experiment.R`
now runs the procedure of the paper: the grid of (K, phi) at lambda = 0 and the
two-step rule, then a path in lambda on the selected partition and on the true
one, recording the log-likelihood and the effective degrees of freedom of
every fit so that the analysis lets seven criteria choose lambda (AIC, AICc,
GCV, HQ, KIC, BIC, EBIC); the error of the coefficients is measured location
by location (the regime a location was assigned to against its true one), so
that it is defined for any number of regimes selected; 5% of the cells and
every 20th week are blanked and their interpolation is scored; and the grid is
run a second time at the lambda of the AIC, to measure whether the order of
the steps changes (K, phi). Supplementary Material A had recommended choosing
lambda for every (K, phi) and called the reverse sequence wrong in the
direction of a cheaper regime under the ridge; that section now states the
sequence the paper uses, why its cost is the smaller one (both couplings are
of the order of the effective coefficients the penalty removes: 0.05-0.5 of 8
per regime at moderate collinearity, 0.9-2 of 16 at strong, in preliminary
runs) and that the experiment measures it. `analyse-simulations.R` gets a
part R9: the two-step rule applied again from the stored grids with each of
the seven criteria (N recovered as log N = (BIC - AIC)/df + 2); with the BIC it
reproduces the recorded K in all 3,320 replications of the first ten of
stream 1. The blocking stop() in the experiment is removed; smoke-tested.

**The lags of the fuel application, fixed by the pre-analysis.** Stage `lags`
of the three main scripts was run on the data of `fuels-data.R` (no model). The
PACF of the own price is beyond the band 2/sqrt(T) for 100% of the pumps at
lag 1, 56-72% at lag 2, 26-41% at lag 3, 14-22% at lag 4 and 0-17% at lags 5
to 8 (medians over the cities), and the order of the AR chosen by the AIC pump
by pump has median 3-4: the own dynamics reach lag 4, with a weak tail. Lag 52
is not needed (PACF beyond the band for at most 1% of the pumps; the ACF at 52
is the persistence of the series). The cross-correlations of the raw prices
exceed the band at almost every lag, an effect of the persistence, and were not
used. The user asked for one indication valid for every setting instead of a
BIC choice among candidate sets city by city: `LAG_CANDIDATES` is replaced by
one set, `LAGS`, the last four weeks of each of the four blocks (16 lags, as in
a VAR(4), the convention of the Granger tests, and the strong level of the
ridge experiment). Stage `select` now chooses the ridge penalty only, by the
BIC of the pooled model, for each setting, city and fuel, as the user decided:
5 pooled fits per city and fuel, 110 per setting. `FIT_FIXED` becomes
`FIT_LAMBDA`. The lags and the penalties now enter the names of the cached
fits, so that a change of `LAGS` never reuses a fit made with other lags
(before, the cache of stage `select` was keyed on the city and the fuel only).

### 2026-10-05

**SQUAREM with a ridge, and the scripts of the fuel application.** The user
decided the application: local leaders and followers among the pumps of the
eleven metropolitan cities, the three settings of the exploration, with a
multi-product model (four blocks of lags: own price on the fuel, own price on
the other fuel, the neighbour's price on the fuel, the neighbour's on the other
fuel), a ridge on the lags, the missing weeks handled by the Kalman filter
instead of carrying the last price forward, and both F and Wald tests. The ridge
made every fit run the plain EM iterations, so the user asked whether SQUAREM
could be extended to it safely, and to keep EM if not.

The obstacle: the ridge penalizes lambda/2 sum_j w_j s_j^2 beta_j^2 with
s_j^2 the GLS scale of column j, which moves with Sigma_e at every iteration,
so the safeguard of SQUAREM has no fixed objective. The candidate merit, the
penalized log-likelihood with the scale measured at the point itself, was
checked in scratch (`squarem-ridge-check.R`) on a dynamic design like the fuel
models (60 locations, 150 periods, six collinear lags, latent AR(1), spatial
error; 5 seeds x lambda 0.1, 1, 5): along the plain iterations the merit
decreases by at most 1e-5; SQUAREM with it reaches the fixed point of the plain
iterations to 1.1e-7 at a tight tolerance, never fails, and takes 0.38 of the
time at the default tolerances, as close to the fixed point. In the package
the stopping rule stays on the log-likelihood, as in the plain iterations, and
the same check on the implementation (`squarem-ridge-package.R`) gives the
fixed point to 1.4e-7 and 0.47 of the time. Implemented:
`stem_ridge_merit()` wraps the map when lambda > 0 and alpha = 0, and the
safeguard of `stem_iterate_squarem()` compares `merit` when the map returns it
(`stem_merit()`); without a penalty the code path is the previous one (the two
main3 replications reproduced to 4e-12). Lasso and elastic net stay plain. The
test "a penalized fit runs the plain iterations" becomes "a ridge fit reaches
the fixed point of the plain iterations" (and a lasso one stays plain).

The fuel application, `dev/replication/`: `fuels-functions.R` (auxiliary),
`fuels-data.R` (data management, five cached stages, one data.frame per case in
`fuels/<case>/data.RData`), `fuels-main-setting1.R`, `-2`, `-3` (stages lags,
select, fit). Not run on the data, as asked; the functions were exercised once
on synthetic data in scratch (`fuels-smoke.R`). The first scripts are archived
in `dev/archive/fuels-first-scripts/`.

### 2026-10-04

**A general check of the package, the replication scripts and the paper,
after the defect of the log-likelihood.** Comparing the first replications
of `main3` with `main2` showed log-likelihoods four times smaller on the same
data. The cause: `STEM_Estimation()`, in every version since Stem 1.0
(checked on the tarball of 1.0), stored -2 loglik in its iteration matrix
and returned that value multiplied by -2, that is 4 loglik. The rewrite of
2026-10-03 (ea73a9b) returns the log-likelihood of `kalman()` and so fixed it,
but the change went unnoticed and undocumented, because the tests checked
only that the value returned was that of the filter. Consequences: the AIC,
BIC and KIC of SC-STEM in `main`, `main2` and every earlier run weighted the
likelihood four times against the penalty. The two-step rule re-implemented
on the grid files reproduces the recorded selection in all 951 replications
of `main2` and all 48 of `main3` checked; on the 951 of `main2`, the BIC on
the true log-likelihood selects fewer regimes in 40 (correct K 81.2% ->
77.8%, mostly S1w-shr, S3beta-shr, S3theta-shr, S3error-shr). The partitions
were not affected: the label step and the automatic penalty scale use the
location-wise scores, which were always right. `main3` (0f7b744) runs on the
correct log-likelihood. The user asked for a check of everything before the
launch that has to be the last one. Checked, against independent computations
(scratch `audit-*.R`):

- the log-likelihood equals the Gaussian log-likelihood of the observed values
  on the full covariance of the data, gaps included (1e-14), for EM, ECME and
  SQUAREM; the estimates equal the maximum found by `optim()` on that
  likelihood (largest relative difference 2e-5); the smoothed states and
  `STEM_Signal()` equal the conditional means (3e-15);
- `STEM_Simulation()`, which both bootstraps use: mean and covariance of
  20000 simulated data sets against the exact ones, within Monte Carlo error;
- the whole replication of the simulation study, end to end: two
  replications of `main3` reproduced here from the pinned code (same data,
  same selection, same ARI and RMSE, grid log-likelihoods to 4e-12);
- `analyse-simulations.R` runs on the partial results of `main3`;
- the application scripts (design, Granger test on the bootstrap covariance,
  diagnostics), the bootstrap and the cross-validation predictor against
  what the paper states.

Fixed:

- the BIC of SC-STEM counted `d * T` observations also when the response
  has gaps; it now counts the observed values. Nothing changes on complete
  data, which is the case of the simulation study and of the fuel
  application (gaps carried forward), so `main3` and its pin are unaffected;
- `STEM_Kriging()`: argument checks that never fired, a sign in the
  documented formula, and what `se.pred` measures;
- `tests/testthat/test-exact.R`: the log-likelihood, the smoothed states, the
  maximum and the BIC against dense computations;
- the example of `STEM_Fit()`, inside `\donttest{}` and so never run by an
  ordinary check, used `povalley$altitude_std`, which does not exist; it now
  takes the covariates of `povalley` as the other examples do, and
  `R CMD check --run-donttest` passes;
- the paper: the minimum regime size (r + 2, as in the code), the functions
  of Appendix B (`STEM_Data()` and `STEM_Skeleton()` are internal), the
  sample size of the BIC, the settings of Supplement B (the stopping rule of
  1e-4 and the experiments on the old estimator are replaced by the current
  settings and Supplement A);
- this file: the entry of 2026-10-03 (first) gave "3-41" and "up to 33"
  log-likelihood units, measured on the old returned value; on the
  log-likelihood they are about 0.75-10 and 8.

Replication: `run-application.R` and `run-fuels-application.R` are pinned to
ddf90a8, the commit of the fixes above (their caches move to folders named
after it). `run-simulations.R` stays on 0f7b744, on which `main3` is running:
the two commits give identical results on complete data, and the simulated
data are complete. Nothing was launched.

### 2026-10-03 (fourth entry)

**Documentation split by topic, and the Computational Supplement shipped with
the package.** The user asked for the technical material on the estimation to
be in the Computational Supplement of the paper (Supplementary Material A),
attached to the package as documentation, and for the GitHub page to be
divided by topic. The code of the package did not change.

- `docs/` (outside the built package) holds the topic pages: the models,
  computational aspects (the Kalman filter and smoother, the EM algorithm and
  where it is slow, ECME, SQUAREM, the stopping rule against that of D-STEM,
  the limits of the range, the ridge, the numerical experiments), the penalty
  on the coefficients, the design notes of SC-STEM, the function map and the
  references.
- `README.md` keeps the short account of the models, of the estimation
  algorithms and of the range at the boundary, and points to the pages; the
  long sections moved to `docs/`.
- `inst/extdata/Stem-computational-supplement.pdf` is the compiled
  `supplement.tex` of the paper; it has to be refreshed from the Overleaf
  project before a release.
- The function map shows the new estimation engine (`stem_em_fit()`, the
  plain and SQUAREM iterations, `stem_estep()` on `kalman()`, the M-steps of
  EM and ECME and their conditional updates).
- The getting-started vignette states the new stopping rule and the
  algorithms; the computational notes point to the supplement.
- Replication (not run, as asked): `run-fuels-application.R` and
  `run-application.R` pinned to 0f7b744, the commit of the simulation study
  `main3` (the fuel application was still on 4325536, before the algorithms
  and the limits of the range); both pass `STEM_control()` explicitly, as the
  simulation runner does, and `run-application.R` no longer overrides the
  tolerances with those of the old rule (`precision = 0.1`, `max_iter = 8` on
  the grid, `max.iter = 60` on the pooled fit). Their caches live in a folder
  named after the commit, and the summary of the fuel application leaves out
  rows written by an earlier Stem: a cache of the old estimator would
  otherwise be read as a result of the new one.

### 2026-10-03 (third entry)

**The range kept within what the distances identify, and a safer
Newton-Raphson step.** The comparison of the algorithms had left SQUAREM 20%
slower than EM on the small cells. The cause was in the M-step, not in
SQUAREM: on the few locations of a small regime the likelihood is flat in
theta beyond two limits (the range shorter than the closest pair, where the
field cannot be told from the nugget; much longer than the farthest pair,
where it cannot be told from an effect common to the locations). There the
Newton-Raphson step took steps of any size (theta from 0.02 to e^-700 in one
M-step was seen), every later step paid up to 30 grid searches of 100
evaluations each, and SQUAREM, extrapolating along the flat direction, ran to
`em_maxit`. The user asked for the fix proposed (implement, then verify):

- theta is kept between -log(0.95)/h_max and -log(0.05)/h_min
  (`stem_theta_limits()`), in the Newton-Raphson step and in the
  extrapolation of SQUAREM; a range at a limit is reported in
  `convergence.par$theta.bound` and in `theta_bound` of SC-STEM, and the
  runner records it in the params file;
- a Newton step is accepted only if it does not increase the objective, and
  is halved otherwise (the M-step is now a generalized M-step, which it was
  not: a step with a Hessian that was not positive definite could worsen it);
- one grid search per step, within the limits, after which a second Hessian
  that is not positive definite ends the step; `nr_hess_maxit` is removed.

Verification (scratch `theta-v1.R`, `sq-diag-A.R`, `alg-compare-grid3.R`):
on the 24 fits of the reference cell, where nothing binds, the results are
those of the previous commit to 1e-9, iterations included; on the small cell
that cost most, EM went from 362 to 29 s and SQUAREM from 1416 to 25 s, no fit
reaching `em_maxit`; on the paired grid (6 cells, 2 replications) the fits at
the true K are unchanged in 40 of 48 and better in the 8 others (all in the
weak scenario at n = 40, by 0.3 to 16.5 log-likelihood units), the selections
unchanged, and the time relative to EM is SQUAREM 0.53 (0.58 on the small
cells, 0.46 on the large), SQUAREM-ECME 0.62, ECME 1.17. Over the grid,
SQUAREM takes 1057 s against 1747 s before the fix and EM 2102 s against
2585 s.

### 2026-10-03 (second entry)

**`kalman()` restored, and the simulation study re-pinned as `main3`.** The
user wanted `kalman()` back: the model is in state-space form, and the Kalman
filter and the smoother should be visible next to the E- and M-steps. It is
now the state-space core of an iteration (the filter, the fixed-interval
smoother, the smoothed initial state and the lag-one covariance smoother), the
E-step builds the expected sufficient statistics from its output, and
`STEM_Estimation()` calls it once more at the estimates. A pure move: the
equivalence test against the old code still agrees to 1e-13 and the tests
pass. The user asked to relaunch main2 with the new code: the runner is
pinned to 6808cfe and writes under the tag `main3`, so that its results never
mix with those of main2 (same design, EM algorithm, regularization 0.01).
main2 was still running when the runner was changed, so the runner is not
synced to Drive until the user has stopped it: every capped replication
re-sources the runner.

### 2026-10-03

**Four estimation algorithms, `regularization = 0`, tolerances 1e-3, and the
EM code split by role.** The second run of the simulation study (main2) was
2.8 times as costly as the first, and the user asked why and what Fasso and
Finazzi use. D-STEM v2 stops when the relative change of the parameters or of
the log-likelihood is below 1e-4 (at most 100 iterations); its criterion on the
log-likelihood is relative, about 0.4 to 1.5 units on our fits, hence much
looser than our absolute 1e-3. The user asked to apply it as well: it is kept
absolute for now and the question is put back to the user, because an earlier
check showed it stopping several units short and the differences that matter
(between fits, in the criteria) are absolute.

What makes the EM algorithm slow was measured by tracing it parameter by
parameter on a regime of the reference cell. Two directions are slow: the split
of the common variation between the latent process and the spatial field (rate
about 0.93, which sets the number of iterations), and the joint move of the
intercept and m0 (rate 0.988, little likelihood, where the stopping rule leaves
the error; with G near 0 it becomes a G-m0 ridge). An earlier statement that
the 100-200 iterations of the refits came from the intercept-latent ridge was
imprecise, and was corrected with the user.

The same comparison found that the iterations did not converge to the maximum
of the likelihood: `kalman()` added `regularization = 0.01` to every matrix it
inverted, the filter that computes the likelihood did not, and the fixed point
moved along the slow directions (variance parameters off by 3-13% on the
regimes, Sigmaeta of a pooled fit by 84%, the fit at the true K of the larger
cells about 0.75-10 log-likelihood units short; measured as 3-41 on the value
then returned, four times the log-likelihood, see 2026-10-04). With 0 the same
iterations reach the BFGS maximum of the exact likelihood. Every matrix
inverted is positive definite at admissible parameters, so the default is now
0.

Tolerance test (scratch `tol-test*.R`: 6 cells from n40/T60 to n400 and
n200/T365, omega 0 included, 2 replications, whole grid): 1e-3 on both criteria
left the estimates within 1-2% of those at the previous defaults and the
selection unchanged; 1e-2 lost up to about 8 log-likelihood units (33 on the
old returned value) at the true K and moved theta by up to seven times.
Defaults: 1e-3 and 1e-3.

ECME and SQUAREM were prototyped in scratch (`ecme-proto.R`, `em-accel.R`) and
written up for the user in an internal note (`em-algorithm-notes.tex` in the
Overleaf project). The user then asked for both in the package, with SQUAREM
as the default, the choice as a setting, and clean code: separate E-step and
M-step files per algorithm and a wrapper. So:

- `R/estep.R` (E-step, shared), `R/mstep-updates.R` (the conditional
  updates), `R/mstep-em.R` and `R/mstep-ecme.R` (the two M-steps, with the
  augmented Kalman filter of the ECME step), `R/squarem.R`, `R/em-fit.R` (the
  wrapper, the plain iterations, the stopping rule, the data of a fit).
  `kalman()`, `B_function()`, `cov_lagone()` and `Q_function_addendo2/3()` are
  gone; the Q values they computed were never used.
- Equivalence (scratch `equiv-test.R`): with `algorithm = "EM"` and the old
  regularization the new iterations reproduce `kalman()` to 1e-13, step by step
  and over whole fits, on complete and missing data, `latent`/`spatial` off,
  ridge and elastic net, p = 2. One deliberate difference: the smoother gain of
  the initial state uses t(G) (it used G, which matters only when G is not
  symmetric).
- SQUAREM: a first driver with repeated backtracking took 246 steps on the
  G-m0 ridge of a regime where the prototype took 82; the scheme of the R
  package SQUAREM (Du and Varadhan 2020: steplength capped by an adaptive
  bound, one extrapolation per cycle, fallback to the two plain steps) took
  66, and 45 instead of 258 on a pooled fit with its maximum on the boundary.
  That scheme is the one implemented.
- Penalized fits run the plain iterations: the scale of the penalty is
  re-measured at every iteration, so there is no fixed objective for the
  safeguard of SQUAREM.
- `STEM_Estimation()` now returns the log-likelihood and the smoothed states at
  the estimates (one filter and smoother pass after the iterations); they used
  to be those of the last-but-one iteration.
- The user had asked that fits start from estimates, never from distant
  values. The alternation and the final refit already did; the penalized fits
  did not (only the partition of the unpenalized fit was passed on). They now
  start their regimes from the refit of that partition, read from the shared
  refit cache, so no new argument was needed. The test of the start of a
  penalized fit passes the cache explicitly.

Comparison of the four algorithms (scratch `alg-compare-*.R`). Single fits of
the regimes and the pooled model of the reference cell, 10 replications, at
three tolerances: at 1e-3 SQUAREM took 4.5 s per replication against 11.9 s
for EM, 14.4 s for ECME and 5.1 s for SQUAREM-ECME. It was the only one whose
estimates had the statistical accuracy of the exact maximum (relative RMSE of
theta 18.8% against 18.5%; EM 24.5%). In one regime EM, ECME and SQUAREM-ECME
stopped on a plateau where the field takes the place of the latent process
(theta 0.21 against 2.45, 4 units short), which SQUAREM crossed.

The whole grid (6 cells, 2 replications) first showed SQUAREM changing the
partitions of a weak scenario at n = 40 (K = 1 selected where EM selected 2,
the fit at the true K 15 units lower): its extrapolations inside the
alternation, from fits stopped at a log-likelihood tolerance of 1, moved the
path of the labels. The alternation now runs the plain iterations of the map
(`stem_control_em()`), and with that the partitions and the selections are
those of EM in all 12 jobs. Wall times across processes were not comparable
on this machine (cores of different speed, 10 workers), so the final
comparison runs the four algorithms in the same process, order rotated, and
estimates the ratio to EM net of the position (the first algorithm of a
process runs slower): SQUAREM 0.75 (0.54 on the large cells, 1.20 on the
small), SQUAREM-ECME 0.86, ECME 1.09; at the true K SQUAREM is at most 0.10
log-likelihood units from the best fit, EM up to 1.29. The cost on the small
cells comes from regimes whose range runs to the boundary (theta to infinity,
field indistinguishable from the nugget): SQUAREM follows that direction until
`em_maxit`, and every step pays repeated grid searches in the Newton-Raphson
step. It is the boundary problem the user put aside for a separate study.

### 2026-10-02 (sixth entry)

**Names: the number of regimes `K`, the loading matrix `A`, the parameter
count `df`.** The user chose the notation after a preview of the paper and of
the code: `K` regimes with index `k`, as in the SCDA package; `G` the
transition matrix, as in STEM and D-STEM; `A` the loading matrix, which was
`K` in the code and Lambda in the paper (Lambda would have collided with the
lambda of the ridge). The parameter count was `k` in `info_crit` and `npar` in
the grid; next to `K` for the regimes that would be a trap, so it is `df`, as
in the paper. `mink`/`maxk`, already deprecated, are removed. No aliases, at
the user's request; but an old name stops with a message, because the end-to-
end check showed that `SCSTEM_Estimation(mod, k = 2)` ran silently with
`knn = 2` and three regimes, R completing `k` to `knn` by partial matching
(`stem_renamed_args()`). Renamed everywhere: the sources, the docs, the tests,
the vignettes (in the Woodbury identity of the computational notes the
generic matrix is now M, to leave A to the loadings), the README, the
replication scripts (the simulation runner writes `K_hat`, `K_correct`; the
analysis maps the names of the first pass), the fuel scripts, and dev/paper.
In the paper Lambda -> A in the manuscript and both supplements, and in
Section S1 of Supplement A the loading matrix, written Z there, is A too.

### 2026-10-02 (fifth entry)

**`run-simulations.R`: streams, to run the study on several machines at
once.** The user runs the study on two machines with different numbers of
replications (for instance 10 on 3 cores and 40 on 7) and wants the data
different and the results cumulable. A new setting `stream` (one per machine):
replication r of stream s is recorded as `rep = 1000 (s - 1) + r` and drawn
with seed `1000 rep + 1`, so no two streams share a seed; each stream writes
its own files, `<tag>_s<stream>*`, so the results folders of the machines are
put together by copying the files side by side, and `analyse-simulations.R`
stacks every stream of a tag (`sim_read_streams()`, which also reads the
`<tag>.csv` of the first pass). The seed of a redraw after a timeout is now
`1000 rep + 1 + 100 (attempt - 1)`, at most 9 attempts: the former offset of
100000 per attempt collided with the seeds of replications above 100. Stream
1 at its first attempt keeps the seeds of the first design. `rep_to` defaults
to `nrep`. The time limits and the cost are measured on every stream found in
the results folder. Checked: stream 1 with one replication and stream 2 with
two, written to one folder, give three different data sets (seeds 1001,
1001001, 1002001), and the analysis reads them as three replications of the
cell.

### 2026-10-02 (fourth entry)

**Shared final refits: `refit_cache` in `SCSTEM_Estimation()`, one cache per
grid in `SCSTEM_Infocrit()`.** With the EM algorithm run to convergence the
final refits take 85-92% of a grid (instrumented on replication 1 of the
reference cell: 398 of 462 s in S1-strong-ind, 342 of 401 s in S1-weak-shr,
224 of 244 s at d = 40, T = 60), and many of them refit the same partition:
the penalized fits at one k often end where the unpenalized one ended (13 of
22 configurations in S1-strong-ind, 7 of 22 in S1-weak-shr). The refit
depends on the data, the partition and the settings, not on phi, so it is now
computed once per (k, partition) and shared. A penalized fit run on its own
shares the refit of the unpenalized fit it starts from. The bootstrap and the
cross-validation never share refits across samples or folds. `refit_reused`
in the output says whether the refit was shared.

Checked on replication 1 of the reference cell against the grid without the
cache: the configurations whose refit is not shared have identical BICs; the
grid takes 107 s instead of 185 s in S1-strong-ind (about 2.9 times the old,
unconverged package, against 5 without the cache); the selection is unchanged
in S1-weak-shr and moves from phi = 0.025 to phi = 0 in S1-strong-ind, where
the two penalties had the same partition and the BICs differed by 0.005. The
shared refits are not always within tolerance of the ones they replace: at
k = 2, 3 the BICs moved by less than 0.01, at k = 4 of S1-strong-ind by 2.1
(the shared refit the better), at k = 3, 4 of S1-weak-shr by 30 (the shared
refit the worse). The 30 is a regime of 6 locations whose refit, warm-started
from the alternation of phi = 0, sits on the boundary theta = 0 (an infinite
range) and stays there after 3000 more iterations under em_stop = "all"
(log-likelihood -3323.8), while from the alternation of phi = 0.025 it
converges to theta = 0.24 (-3309.0). The trap exists without the cache (the
fit at phi = 0 reported it); the cache only makes the fits of one partition
agree. A guard on warm starts at the boundary is a candidate, to be verified
before it is adopted.

### 2026-10-02 (third entry)

**The simulation study, second design (`tag = "main2"`): new values, a fourth
geometry, a time limit per replication, and the figures of the design.**
Replication material only (`dev/replication`), nothing in the package.

- `run-simulations.R`. Regime 3 of S1-strong has nugget 0.40, partial sill
  0.20 and range 1 (it had 0.50, 0.10 and 0.5), regime 2 range 1.5: with a sill
  one sixth of the error variance and a range of three nearest-neighbour
  spacings the spatial parameters of regime 3 were not identified, and a
  long-range field absorbed part of its almost white latent process. The total
  error variances, all the assignment score sees, are unchanged, so the
  separation of the scenarios is the same. A fourth geometry, block `spread`:
  omega = 1 with the spread of a regime fixed at its value at omega = 0.7
  (medium-low overlap, 20% of the locations nearer another centre), beside
  omega = 1 with the spread of the network fixed (strong separation, 10%); 72
  cells more, 332 in all. Every replication runs in its own process (callr)
  under an adaptive time limit (10 times the median of its cell after three
  replications, 20 times the median of its (n, T) before, 3600 s with nothing
  recorded); a replication over the limit is written to
  `<tag>-timeouts.csv` and drawn again with seed 1000 rep + 1 + 100000
  (attempt - 1), up to three attempts. `setTimeLimit()` could not do this: it
  fires once and the `tryCatch()` blocks inside the estimator swallow it. The
  grid runs with `control = STEM_control()` passed explicitly; pinned to
  9f9c81e. Checked: the dry run, a replication stopped at a 5 s limit, drawn
  again and recorded as failed, and a replication completed inside the
  process.
- `analyse-simulations.R`: the core read by the four geometries (labels
  total overlap, medium-high overlap, medium-low overlap, strong separation);
  the cost also by geometry; new part R8, the timeouts and the replications
  drawn again, by scenario-variant, geometry, n and T.
- `design-figures.R`, new: the figures and tables that describe the design
  (the four geometries, one data set in time and space at each level of
  separation, the decay of the correlations, the variance of every regime by
  source, the parameter values, the separation measure, the blocks), computed
  from the definitions of the runner so that the description cannot drift
  from the design. They go into the new Supplementary Material B of the paper.

### 2026-10-02 (second entry)

**`em_stop`: the two criteria of the EM algorithm combined by "any"
(default) or "all".** The first version of the rule required both criteria at
the same iteration. The user disagreed: D-STEM stops at either, and requiring
both costs several times the iterations for a gain in the log-likelihood of
hundredths of a unit. `STEM_control(em_stop = "any")` is now the default,
`"all"` the option. Checked on the slow ridge (regime 1 of the reference cell,
G = 0.8, traced to 3000 iterations): "any" stopped at 86-110 iterations within
0.03 log-likelihood units of the maximum, "all" at 129-507 within 0.002. On the
fit at k = 3, phi = 0.05 of five replications of the reference cell, the
relative bias and RMSE of the regime parameters against the truth are the same
under the two rules (G: bias -15.5% against -17.0%, RMSE 32% against 35%; the
other parameters within about three points, "any" never the worse), with 85/35/79 EM iterations per regime
against 193/182/363, and 18 s per fit against 63 s. The grid of the simulation
study on one replication of the reference cell takes 185 s under "any", 669 s
under "all", against about 40 s with the old, unconverged rule.

Two things learned on the way, recorded for whoever repeats these checks:

- `pkgload::load_all()` sources the test helpers by default, and
  `tests/testthat/helper-models.R` sets the loose session option
  `Stem.control` that makes the test suite fast. A scratch run with
  `control = NULL` then runs with the test settings. Use
  `load_all(helpers = FALSE)`, or pass `control` explicitly.
- A grid with loose final refits followed by refitting only the reported fits
  at the defaults is NOT the same estimator: the departures that initialize the
  unpenalized fits are computed from the pooled fit of the grid, which then
  runs loose too. On replication 1 of the reference cell the initial
  partitions of S1-weak-shr changed (ARI 0.31 between the two grids) and
  S3-beta-shr selected k = 3 instead of 2. The simulation study runs at the
  package defaults.

### 2026-10-02

**The EM algorithm stopped too early, and the regime fits of SC-STEM started
from the pooled values: `STEM_control()`, a new stopping rule, data-based and
warm starts.** The first ten replications of the simulation study showed biased
regime parameters (the persistence G of the third regime 0.46 against 0.20, the
intercepts pulled towards the pooled one). Traced on single replications:

- the rule stopped when the relative change of the whole parameter vector AND
  the relative change of the log-likelihood fell below `precision`. The vector
  contained the d fixed loadings (and C0), whose norm grows with sqrt(d), so the
  rule loosened with the network; with log-likelihoods of 1e4 a relative
  tolerance of 1e-3 stops with several units still to gain. On a regime with
  G = 0.8 the old rule stopped at 32-39 iterations, 1 log-likelihood unit short,
  and at the pooled tolerance 0.01 at 5 iterations, 9-12 units short; the D-STEM
  rule (per-parameter relative change OR relative log-likelihood change, 1e-4)
  stopped at 12-18 iterations, 4-6 units short, because the relative
  log-likelihood fires first. The per-parameter relative change (1e-4) AND the
  absolute log-likelihood change (1e-3) stopped at 129-507 iterations within
  0.002 units of the maximum. That is the new rule.
- in SC-STEM every regime fit, in the alternation and in the final refit,
  started from the user's (pooled) starting values and ran at precision 0.1,
  while the pooled fit ran at 0.01, so the log-likelihoods compared across k had
  different accuracies. A regime with a persistent latent process stayed on the
  ridge between its intercept and the level of the latent process: a warm start
  from the alternation reached a lower log-likelihood than a start from the
  regime's own least squares in every replication tried. Now the first fit of a
  regime starts from its own least-squares coefficients, every later one from
  its previous estimates, and the final refit runs with the settings of the
  pooled fit. The pooled fit used only to initialize the partition from the
  departures runs with the loose settings of the alternation.
- the bootstrap refits follow the same starting rule as the original fit
  (least squares on the bootstrap data), not the estimates, which in the
  bootstrap world are the truth.

`STEM_control()` holds the settings (tolerances and iteration limits of the
final fits and of the EM inside the alternation, of the alternation, of the
Newton-Raphson step), all the estimation functions take `control`, and
`options(Stem.control = )` sets them for a session; the historical arguments
remain as overrides. The tests run under a loose session setting, since they
check behaviour and the new defaults cost ten to twenty times more per fit:
they are what convergence costs, the old fits were simply not converged. On the
reference cell of the simulation study a fit at k = 3 takes about 60 s instead
of 3, and the whole grid 11 minutes instead of 40 s.

### 2026-10-01 (third entry)

**`dev/replication/run-fuels-application.R`, the full fuel application, written
but not run: the user checks it and launches it.** 88 models on the eleven
cities: part A, gasoline and diesel at the same pump in both directions (22);
part B, leaders and followers, three settings by two fuels (66). Each model is
the weekly relative price (Sunday price, national mean and pump mean removed,
cents) on its own lags, the lags of x (the other fuel, or the nearest pump of
the X set; the mean within r* as an option), the fiscal pulses, the latent
process and the spatial error; the grid k = 1..K_max with the default phi grid,
K_max = min(3, floor(n / m)) with m = 5 fixed by the user without validation,
so that fewer than 10 pumps give the pooled STEM; the two-step rule; the
refit-with-clustering bootstrap of the selected model and of the pooled one;
the Granger test regime by regime as the Wald statistic on the bootstrap
covariance of the lags of x; Ljung-Box diagnostics of the residuals and of
their squares. `min_cluster_size` is set to m: the package default, the number
of covariates plus two, would be 17 here and would forbid the small regimes the
rule allows. Stages are cached, models can be split across processes
(`--job=i/N`), `--dry-run` prints the plan and `--B=20` runs a pilot whose
grids are reused by the full run; every stage records its minutes. It pins
Stem at 4325536, the commit of the simulation study (the package code has not
changed since).

A smoke test on two tiny models (Messina and Florence, setting 2, B = 2) ran
every stage. It also found that some bootstrap draws are pathological: one draw
of eight on the 14 independents of Messina took about 13 minutes against 3 to 4
seconds for the others. The profile puts all the time in the Newton-Raphson
step of the spatial parameters inside the EM (`Q_function_addendo1` and its
numerical Hessian, up to 50 iterations and 30 Hessian retries per EM
iteration), which runs to its limits when the likelihood is flat in the range,
as in a regime of five sparse pumps. The package is not touched: the cost and
the remedies go to the user first. `run-fuels.R` is marked as superseded.

### 2026-10-01 (second entry)

**The cities of the fuel application.** The frame is the 14 metropolitan
cities; a city enters if it keeps at least 10 independent sites after the
pre-treatment, because the application also asks about independents and major
brands, and with fewer independents the share of rejecting pairs in a city moves
by ten points with one pump. Nine pass (Naples, Rome, Catania, Turin, Bologna,
Milan, Palermo, Messina, Bari); the user added Florence and Venice by choice.
Genoa (no independent) leaves, Messina enters. The list does not depend on the
filters: allowing spells up to 180 days, or coverage down to 90%, the same nine
cities pass, so the spell stays at 90 days. `fuels-describe.R` and
`fuels-leader-follower.R` now take the cities from the pre-treatment output
instead of a fixed list.

**Both Granger tests in `fuels-leader-follower.R`.** The classical F test and
the HC1 Wald test are computed on the same regression and reported side by
side, at the user's request: the comparison shows the reader how much the
non-constant variance of sticky prices inflates the F test.

### 2026-10-01

**`dev/replication/fuels-pretreatment.R`, the data pre-treatment of the fuel
application, with a register of every step.** Two passes over the station-level
files. The first summarizes every code over its life, computes the national and
metropolitan daily means of all reporting non-highway pumps (the centring
references) and shards the window rows by province. The second, province by
province, links successive codes at the same place into sites (within 50 m,
starting no more than 7 days before and 120 days after the predecessor ends;
the newer code wins on overlapping days), then filters: not on a highway, valid
coordinates, coverage of both self prices >= 95%, no gap > 28 days, >= 12
changes a year, no unchanged spell > 90 days, one site per point. Each step is
counted by province, brand group, area and type, because the notes need a "Data
pre-treatment" section that explains what every step removes and how it shifts
the sample. The window, 2022-01-03 to 2026-06-28, avoids the 2018 and 2021 holes
and starts with the 2022 change in stickiness; it is a parameter. A part of the
files whose province column is missing altogether is labelled unknown rather
than Naples: the "NA" fix applies only to the code read as a missing value.
Step 4 also drops sites more than 150 km from the median location of their
province: three kept sites were geocoded hundreds of kilometres away (one of
Nuoro in Milan, one of Cosenza in Rome, one of Palermo near Messina), while the
genuine edges of provinces stay within 95 km.

**`dev/replication/fuels-describe.R`**: maps of the kept sites, residual
stillness, and pump-by-pump ADF and KPSS tests and autocorrelations of the
centred series (national and metropolitan centring) in the ten cities, with
daily and weekly series of representative sites.

**`dev/replication/fuels-leader-follower.R`**, an exploratory analysis asked for
by the user: leader and follower between major brands and independents (major
Y and independent X; independent Y and major X; major Y and other-major X), with
X the nearest pump or the mean within a radius. The radius r* is read from the
correlogram of weekly relative prices, per city, as the distance at which the
excess correlation halves; the fixed radii 1, 2 and 5 km are the sensitivity.
An exponential fit was tried for r* and kept only for reference, because one
scale does not fit every city (Palermo decays on two, and the fit returns
7.7 km). Pearson and Spearman correlations, cross-correlations and Granger
tests in both directions are a screening device only; the formal inference of
the application remains the package bootstrap. The tests are Wald tests with a
heteroskedasticity-consistent (HC1) covariance: the plain F test of a first run
rejected for 77-86% of the pairs whatever their distance, and a check on Rome
showed the cause was the non-constant variance of the sticky daily changes, not
the common factor (adding its lags changed little; the robust test brought the
rejection rate of pairs more than 10 km apart to 13%). Because 13% is still
above the nominal 5%, every Y also gets a reference X drawn at random beyond
10 km, and the rates are read against it.

### 2026-09-30

**`dev/replication/run-fuels.R`, a first runner for the fuel-price
application.** Province by province and in both directions, it builds the
weekly panel of the self-service prices (the price in force at the end of the
week by default, since averaging within the week blurs the timing a Granger test
reads; non-highway pumps with both fuels on at least 95% of the weeks, gaps
carried forward), and models the weekly changes in cents with the own lags, the
lags of the other fuel, the lagged gasoline-diesel spread demeaned by pump and
fiscal period (the error-correction term) and a pulse at each change of the
excise duties. It fits the grid k = 1..3 with the default phi grid, selects by
`SCSTEM_Select()`, and tests the lags of the other fuel by the likelihood ratio
against the fit without them on the same partition (`max_iter = 0`), for the
pooled model and regime by regime; the benchmark is the same test pump by pump
by least squares, combined as Dumitrescu and Hurlin (2012). The Naples code
"NA", read as missing in the data, is restored; pumps at the same point are
moved by about a metre. Tested end to end on Isernia (36 pumps, 74 weeks): 76
seconds. There the likelihood-ratio test of the model rejects in both
directions while the pump-level tests reject at the nominal rate in one of
them, a first sign that the model-based test may be anti-conservative, to be
checked by bootstrap or placebo before any reading.

### 2026-09-29 (third entry)

**`dev/replication/analyse-simulations.R`, the analysis of the new design.** It
reads the design from the runner beside it and the results of one or more tags
(a (cell, replication) recorded under two counts once), and writes to
`output/`: R1 the reference cell, every scenario-variant; R2 the core,
selection and recovery against n, T and omega, as wide tables and as a figure
with the scenario-variants by rows and T by columns; R3 how often S2 is split;
R4 the robustness blocks against the reference cell; R5 bias and RMSE of the
parameters at the reference, the reference case regime by regime; R6 the ARI at
the true k against phi, from the grid, and the phi selected; R7 the seconds per
replication by n and T. Each part is guarded, so it runs on partial results.
Tested on the check run of 2026-09-28, on the 83 runs of the stopped factorial,
and on a constructed set exercising R4. It replaces, on the Drive, the analysis
of the previous design, which stays archived in `dev/archive/sim-design-generic`.

### 2026-09-29 (second entry)

**The simulation design is organized in blocks.** The full factorial of the
previous entry was stopped after 83 of its 17,640 runs: at the measured mean of
about 130 seconds a run it would have taken some 650 core-hours, over three
days on the 7 cores of the machine that runs it. The SETUP now holds the
scenario-variants (`SIM_SCENARIOS`), the reference level of every margin
(`SIM_REFERENCE`) and the blocks (`SIM_BLOCKS`), each the full factorial of the
margins it lists around the reference; the design is their union, a shared cell
run once. The design adopted:

| block | margins | cells |
|---|---|---|
| core | n {40, 100, 200} x T {60, 120, 365} x omega {0, 0.70, 1} | 225 |
| balance | unbalanced | 8 |
| knn | 3, 10 | 18 |
| n400 | n = 400 | 9 |

260 cells, 2600 runs at `nrep = 10`, about 58 core-hours on the costs measured
one replication per (n, T), more where omega = 0, which cost about twice as
much in the stopped runs. The core crosses the three factors that govern
recovery; the balance, the graph and the national scale of n = 400 are
robustness margins, varied one at a time. `--blocks=` runs a subset of the
blocks; margins given on the command line replace the blocks with one
factorial block around the reference. Of the 83 runs already recorded, 8 are
cells of the new design.

### 2026-09-29

**The SETUP of the simulation runner is the full experiment, first pass.** The
nine scenario-variants of the paper (S2, S1s-ind, S1s-shr, S1w-shr and the five
S3x-shr) crossed with n {40, 100, 200, 400}, T {60, 120, 365}, omega
{0, 0.70, 1}, balanced and unbalanced regimes and knn {3, 5, 10}: 1728 cells
with three regimes and 36 for S2, where overlap and balance are vacuous, 1764
in all, with `nrep = 10`, 17,640 runs, about 51 hours on 13 cores on the costs
measured one replication per (n, T). The ten additional scenario-variants stay
listed in the SETUP comment, run at the reference levels only. No network below
n = 40: the grid reaches k = 4 regimes of at least six locations. `keep_obs`
drops from 5 to 1, since five full records per cell of the full factorial
would take several gigabytes; any replication can be regenerated from its seed.

### 2026-09-28 (fifth entry)

**The simulation runner is ready for the main run.** `run-simulations.R` pins
Stem to commit 4325536, records in `<tag>-grid.csv` the ARI of every partition
of the grid (so that the recovery at any penalty, or under another rule, can be
read without refitting), and checks the installed commit from the DESCRIPTION
file without loading the package: loading it before the update made the session
read the new files with the old index (`R_decompress1` warnings). One
replication of the 19 cells on 13 cores: no failure, 22 fits per replication,
14 minutes of core time per replication of the design (median 33 seconds per
cell, 118 for S3Seta-shr), so 100 replications take about 23 core-hours, 1.8
hours on 13 cores.

### 2026-09-28 (fourth entry)

**`SCSTEM_CV()`: blocked cross-validation as an optional validation.** The
function takes one `SCSTEM_Estimation` fit or a named list of them, refits each
on the folds of one or more schemes (`"LKLO"`, `"LKTO"`, `"LKLHTO"`,
`"random"`, as in Otto, Fasso and Maranzano, 2024), predicts the removed cells
and returns the root mean squared and mean absolute errors, fold by fold and
over the scored cells, with the rank of each model within each scheme. The
models must share the response, the locations and the time points; their
covariates, number of regimes and penalty may differ, so the pooled model, the
selected one and models with other covariates are compared on the same folds.
The folds depend only on the seed, the scheme and the dimensions of the data,
so separate calls with the same seed use the same folds (tested).

A removed cell of a retained location is predicted by `SCSTEM_Complete()`. A
removed location inherits the regime of its nearest retained location and is
predicted by the conditional mean given the retained locations of that regime
at the same time point (`scstem_predict_new()`, internal): `STEM_Complete()`
on the regime augmented by the new locations. The test checks the two against
each other to machine precision.

**Why, and why it does not select.** The two-step rule of `SCSTEM_Select()`
chooses `(k, phi)` in sample: which partition describes the monitored
locations best. An application often needs something else, prediction at
unmonitored sites or in missing periods, and the partition exists only at the
monitored locations. The cross-validation answers that question for a model
already chosen, and ranks alternative models by it; it is never used to choose
the hyperparameters, and it costs one refit per fold and model, which is why it
is optional and separate. It replaces, in the package, the development script
`dev/paper/04-cv.R`, whose four schemes it implements; the script predicted a
removed location with the full inverse covariance and zeros for the missing
residuals, while the function conditions exactly on what is observed at each
time point.

### 2026-09-28 (third entry)

**Step (S2) of `SCSTEM_Select()` chooses `phi` by the criterion.** At the `k`
chosen in (S1), the rule now takes the admissible fit with the smallest
criterion (BIC by default) over the whole grid, ties towards the smaller `phi`,
in place of the smallest `phi` on the stability plateau of the Adjusted Rand
Index. (S1) is unchanged: the modal BIC winner over the band, the pooled model
competing. The argument `delta` of the plateau is removed, and `step2` now
holds the criterion of every fit at the selected `k`. The default grid of
`SCSTEM_Infocrit()` gains the strong values 0.5 and 1, which lie outside the
default band `c(0.025, 0.2)` and so compete only in (S2).

**Why.** Under the warm start the plateau rule picked the bottom of the band in
359 of 360 replications: it could not find where the penalty helps. The rules
were compared on the same grids, `k = 1..6` by
`phi = c(0, 0.025, 0.05, 0.1, 0.2, 0.5, 1)`, in 190 replications of the 19
scenario-variants (18 with three regimes, S2 with one):

| rule | k correct | k = 2 | k > 3 | ARI | ARI, independent | ARI, shared | S2: k = 1 |
|---|---|---|---|---|---|---|---|
| (S1) + plateau (before) | 91.1% | 5.0% | 3.9% | 0.948 | 0.996 | 0.918 | 100% |
| (S1) + smallest BIC at that k (now) | 91.1% | 5.0% | 3.9% | 0.952 | 0.996 | 0.925 | 100% |
| smallest BIC over the whole grid | 90.0% | 5.6% | 4.4% | 0.953 | 0.996 | 0.926 | 100% |
| elbow per phi, then smallest BIC | 85.6% | 13.3% | 1.1% | 0.925 | 0.987 | 0.885 | 100% |
| phi of the smallest BIC, then elbow | 83.9% | 15.0% | 1.1% | 0.918 | 0.988 | 0.873 | 100% |
| (S1) by EBIC, gamma = 1, then smallest BIC | 91.7% | 5.0% | 3.3% | 0.953 | 0.996 | 0.925 | 100% |

The elbow rules under-select, in the shared variants; the extended BIC with
`gamma <= 0.5` changes no choice and with `gamma = 1` one in 180, so (S1) stays
as it was. On `phi` the two rules are equivalent: at the true `k` the ARI is
0.954 at `phi = 0`, 0.955 at 0.025, 0.952 with the smallest BIC, and 0.960 for
the best `phi` chosen knowing the truth, so the warm start and the departures
initialization leave the penalty little to correct. The smallest BIC gains in
the two cells of low separation (+0.03) and loses in S3Seta (-0.04 and -0.07),
where a penalized fit reaches a higher likelihood with a partition further from
the truth. The rule is kept because it is the one of the published
spatially-clustered regression papers (phi chosen by the criterion), because it
drops a step that did nothing, and because it lets strong penalties sit on the
grid harmlessly: at the true `k` the ARI is 0.865 at `phi = 0.5` and 0.629 at 1,
and the rule chose 0.5 once in 164 replications and 1 never. The two strong
values raise the fits of a grid over `k = 1..6` from 26 to 36; one grid took
106 seconds on average at `d = 100`, `T = 120`. The scripts and results are in
`dev/replication/results/init-verification/` (`rules-*.R`, `rules/`).

The application script and the paper figure of the selection follow the new
rule; the application grid and band, still on the scale that preceded the warm
start, are set to the package defaults.

### 2026-09-28 (second entry)

**The default initialization is the departures from the pooled fit.**
`init_method` gains the value `"departures"`, now the default: the pooled STEM
model is fitted, and each location is described by the mean of its residual
from the pooled signal, the slopes of that residual on its covariates, and the
lag-one autocorrelation and log-variance of what the slopes leave
(`scstem_departures()`); the starting partition is the k-means of these
features, through the same PCA compression, restarts and admissibility filter
as before. It initializes the unpenalized fits only: a penalized fit starts
from the unpenalized solution.

**Why.** The previous default clustered the time means of the covariates, which
depend on the regime only when the covariates do. In the simulation design the
covariate is exogenous, as a regressor usually is, and those starts had an ARI
of 0.04 with the true partition, against 0.78 for the departures. Under the
warm start of the previous entries the departures give an ARI at the true k of
0.953 and the right k in 91% of the replications, against 0.898 and 81% for the
covariate means (190 replications of the 19 scenario-variants). The
covariate means remain available as `init_method = "kmeans"`.

**Cost, and consistency with the grid.** A single fit at `k > 1` now fits the
pooled model first. `SCSTEM_Infocrit()` visits `k = 1` first, computes the
departures once from that fit and passes their k-means to the unpenalized fits
at every `k > 1`, so the grid fits the pooled model once and each of its
members still coincides with the single fit (tested). Should the pooled fit
fail, the covariate means are used, with a warning. The pooled fit uses the
`lambda` of the pooled model, so that it is the same fit as the member at
`k = 1` also under `lambda_by = "size"`. Algorithm 1 of the paper already said
"fit the pooled STEM model to initialize"; the code now does.

### 2026-09-28

**Verification of the warm start, complete.** On the 190
cells-by-replications of the simulation design (19 scenario-variants x 10
replications), grid `phi = c(0, 0.025, 0.05, 0.1, 0.2)`, band `c(0.025, 0.2)`:

| rule, start of the phi = 0 fits | ARI at the true k | ARI = 1 | k correct | k = 2 | k = 4 | seconds |
|---|---|---|---|---|---|---|
| before (first-sweep scale), default | 0.868 | 64% | 76% | 15% | 9% | 49 |
| before, M | 0.834 | 26% | 60% | 36% | 4% | 46 |
| warm start, default | 0.898 | 72% | 81% | 8% | 12% | 39 |
| warm start, M | 0.953 | 85% | 91% | 6% | 3% | 29 |

(The mean time of the warm start with the default start carries one
replication of S3Seta-shr that ran for about forty minutes; without it the two
warm-start rows take 27 and 24 seconds.) With the warm start and the M start
the independent variants are recovered almost exactly (ARI 0.999, k correct
96%) and the shared ones at 0.924 and 88%;
the two cells that stay low, S1w-shr (0.58) and S3beta-shr (0.64), are the two
the separation measure of the design notes marks as hard (4.5 and 6.6
log-likelihood units per location). S2 selects k = 1 in every replication. The
warm-started penalized fits converge in one to three iterations, so a
replication of the grid is also about half as long. The selection rule picks
phi = 0.025, the bottom of the band, in all replications but one. The default
initialization is still the k-means on the covariate means; whether M replaces
it is to be decided.

### 2026-09-27 (third entry)

**A penalized fit starts from the unpenalized solution.** With
`phi_penalty > 0`, `SCSTEM_Estimation()` now fits `phi = 0` at the same `k`
first, with the same arguments, and starts the penalized alternation from its
partition; `SCSTEM_Infocrit()`, which visits `phi = 0` first at every `k`,
passes that partition on through `init_partition`, so no fit is repeated and a
single fit coincides with the member of the grid. The alternating algorithm is
unchanged, and no argument is added: this is a rule of initialization.

**Why.** Verifying the initialization M (k-means on the departures of each
location from the pooled fit) on the new simulation design showed that a better
start gave a worse result: M started at an ARI of 0.78 against 0.04 for the
default, and ended at 0.83 against 0.87, choosing k = 2 twice as often. The
cause was the automatic scale of the penalty, fixed at the first sweep from the
spread of the location scores across the STARTING clusters: nearly identical
clusters under the default gave a multiplier of 1.3 to 3.8, distinct clusters
under M one of 10.6 to 28.3, so the same phi was a penalty five to ten times
stronger. Measuring the scale on the unpenalized fit alone made things worse
(ARI at the true k 0.41 with the default start): a penalty strong from the
first sweep freezes whatever partition it is given, 0.40 on average over twelve
checked cases against 0.81 when the same penalty starts from the unpenalized
solution, and 0.95 at phi = 0.05. Starting from the unpenalized solution fixes
both: the first sweep then scores the locations with the parameters of the
unpenalized fit, so the scale computed there, by the same code, is the one of
the unpenalized fit and does not depend on the start.

**Recalibrated defaults.** Measured on fitted clusters the scale is of the order
of the whole gain of the right regime over the wrong ones, several times larger
than the first-sweep scale on a poor start, and the old grid `phi` in `[0, 1]`
became far too strong: the grid of `SCSTEM_Infocrit()` is now
`c(0, 0.025, 0.05, 0.1, 0.2)`, the band of `SCSTEM_Select()` `c(0.025, 0.2)`,
and the default `phi_penalty` 0.05 (also in `STEM_Fit()`); examples, vignettes,
README and the simulation runner follow. A bootstrap refit with `phi > 0` now
fits `phi = 0` first on every draw, which doubles its cost: that is the
estimator reproduced.

**Status of the verification.** On the first 63 of 190 cells-by-replications
of the design (19 scenario-variants x 10 replications), the new rule gives an
ARI at the true k of 0.97 and the right k in 90% of cases with the M start,
0.87 and 74% with the default start (which chooses k = 4 in a quarter of
them), against 0.87/76% and 0.83/60% before; S2 selects k = 1 always. The rule
selects phi = 0.025, the bottom of the band, every time. The full results, and
whether M becomes the default initialization, are still to come.

### 2026-09-27 (second entry)

**The simulation study is redesigned from the model.** The scenarios of the
previous design separated the regimes one component at a time around a
baseline, with the latent paths coupled at rho = 1 by default; the question
they answered was what has to differ for the difference to be found, not
whether the procedure recovers the model it is built for. The new design starts
from SC-STEM itself, with every block of parameters regime-specific, and is
written up in `simulation-design.tex` in the Overleaf project of the paper:
S1 complete heterogeneity at two levels of separation; S2 the pooled STEM
model, for the false positives; S3 one block common to the regimes (beta, G,
Sigma_eta, theta, or the whole error), the others as S1-strong. Every scenario
comes in an independent variant (rho = 0, error fields by regime: exactly the
SC-STEM model) and a shared one (rho = 1, one field where theta is common),
plus the two mixed combinations where theta is common: 19 scenario-variants.
The runner is rewritten around them as the full factorial of margins set at its
top; the previous runner and its analysis script are kept in
`dev/archive/sim-design-generic`, with the paper text they supported.

**Why the independent variant is not the whole study.** With independent
regimes the latent paths separate the regimes on their own: the expected gap in
the assignment score between the true regime and the nearest wrong one is at
least (v_g + v_h) / (2 s^2) per period whatever the parameters, 45 to 132
log-likelihood units per location over 120 periods in this design. The
parameter contrasts then govern the estimation of the regimes and the choice of
k, hardly the partition, which is why the shared variant, where the same
measure is 4.5 (weak) and 24.6 (strong), is needed to see them act.

### 2026-09-27

**The simulation design is generic.** The baseline regime was the pooled fit of
the Po Valley PM2.5 network of the application, and the plane was mapped onto
a geographic box around the Po Valley; the study then described that network
rather than the method, and the application is going to change. The plane is
now used as it is, with Euclidean distances in its own units, and every value
has a reading that depends on no application: a standardized response with its
unit variance split into covariate 30%, common dynamics 30%, spatial field 20%
and nugget 20% (beta = (2, 0.55), a coefficient of variation of 0.5), a
persistence G = 0.8 (a shock halves in 3.1 periods), and practical ranges of 2
for the error field and 4 for the covariate against a network about 4 units
across. The contrasts read the same way: the covariate effect times 1.25, 1.5
or 2; G down to 0.7, 0.5 or 0.2; a local field, practical range down to 1, 0.5
or 0.25 with a nugget share up to 0.6, 0.7 or 0.8. The Po Valley design, its
scripts, its paper sections and its tables are kept in
`dev/archive/sim-design-povalley`.

**Scenario S0b.** A scratch check on the archived design showed that an error
field drawn regime by regime reveals the partition with no parameter differing:
at `n = 100`, `T = 120`, the true 3-regime partition beat the pooled model by
178 to 287 in log-likelihood, against a BIC cost of about 75, while with one
field over the network it lost by 570 to 663. Moving S1 and S2 to a field by
regime, as first proposed, would therefore have credited them with that
information; they keep the single field, and S0b -- S0 with the field by regime
-- measures what S3, S4 and S5, whose field cannot be global, receive for free.
The design has 255 cells, 252 of them feasible.

**The scenario table** of the paper shows every parameter in every regime at
`K = 3`, with the coupling of the latent paths and the kind of error field; the
analysis writes it from the generator, so the two cannot disagree. The paper
also gains a table of the baseline regime and a statement, per scenario, of
the question it answers.

### 2026-09-26 (sixth entry)

**`SCSTEM_Select()` can select the pooled model.** Step (S1) ranked only the
configurations with `k > 1`, because step (S2) has no meaning at `k = 1`; the
exclusion was applied to the whole rule rather than to (S2) alone, so the rule
could never say that a network has no regimes and left that comparison to the
user. The pooled fit now competes in (S1) at every penalty of the band (its
criterion does not depend on `phi`); when it wins, it is returned with
`step2 = NULL`, and it is also returned when no configuration with `k > 1` is
admissible, where the function used to stop. Checked in a scratch script before
adoption, on 8 replications of five cells: on homogeneous networks (`K = 1`,
`n = 50` and `100`, `T = 120`) the old rule selected `k = 2` in 16 of 16 cases,
the new one `k = 1` in 16 of 16, the pooled BIC being lower by 350 to 900; with
`K = 3` the two rules agree in 23 of 24 replications. The runner drops
`bic_beats_pooled` and `gain_over_pooled`, which only stood in for this; the
analysis reads the null block off `k_hat`; the application skips the bootstrap
when `k = 1` is selected and no longer caches the selection.

### 2026-09-26 (fifth entry)

**The replication scripts keep Stem in step with GitHub.** On a second machine
`install_github()` failed with HTTP 404 because the repository was private; it
is being made public, to be released on CRAN once the paper is done. The
feature check alone could not see a change of behaviour at an unchanged version
number, so `run-simulations.R` and `run-application.R` now also compare the
installed `RemoteSha` with the commit GitHub holds and reinstall when they
differ. The check runs only when a script is run, not when the analysis or the
worker processes read the runner's definitions; offline it falls back to the
feature check. `run-application.R` also finds its own folder when started with
Source in RStudio.

### 2026-09-26 (fourth entry)

**The simulation runner is launched from RStudio and uses every core by
itself.** `dev/replication/run-simulations.R` is opened in RStudio on any machine
that sees the Drive folder and started with Source; the SETUP block gains `mode`
(`"run"`, `"coverage"`, `"dry"`) and `cores` (all physical cores but one by
default). The replications are dispatched one at a time to a PSOCK cluster as
each worker frees up, so the slow cells (n = 400 with T = 60) hold one core and
not the run; only the master writes the files and prints progress with the time
left. Results are identical to a single-core run, since every seed derives from
the replication. Writes retry when a synchronization client holds a file. This
replaces the eight manual processes of the previous README. The whole design is
kept, all 237 cells; the prior cost table is replaced by times measured on this
machine, which put it at about 680 core-hours, some 85 hours on 8 cores.

### 2026-09-26 (third entry)

**The theory of the assignment score is withdrawn from the package and the
paper.** The `score` argument of `SCSTEM_Estimation()` (`"conditional"`,
`"corrected"`), the internal `scstem_cond_scores()`, the diagnostics
`score_last` and `objective_before`, their tests, `dev/replication/check-theory.R`
and the `--score` option of the runner are removed. The label step is exactly
what it was before. The Lemma and the propositions are removed from the
manuscript: they did not change anything a user of the model does, and the
correction that motivated them turned out not to be one. The entry below is
kept as the record of what was tried.

Kept, because they have practical consequences:

* the bootstrap fix: the refits carry `alpha`, `lambda`, `penalize`,
  `lambda_scale`, `lambda_by`, `latent` and `spatial` from the original fit,
  so a fit with a ridge is resampled with it (tested);
* the runner records `bic_beats_pooled` and `gain_over_pooled`, because
  `SCSTEM_Select()` ranks only configurations with `k > 1` and cannot by itself
  declare that a network has no regimes;
* in the paper, the argument for a Potts penalty on point-referenced data, in
  prose, and the supplementary sections on the implementation of the ridge and
  on the order of the hyperparameters.

### 2026-09-26 (second entry)

**The theory of the assignment score, checked before it is adopted.** The
paper's Lemma -- the expected gap between the exact within-regime log-density
and the pseudo-likelihood of the label step is `-(T/2) log|R_h|`, that is `T`
times the Kullback-Leibler divergence of the regime's error law from its
independence approximation -- is proved and verified by simulation to within
2.2 standard errors over twelve configurations, together with its variance.
Its per-location form, `Delta_i(h) = -(T/2) log(1 - q_ih)`, is the Schur
complement of the regime's correlation matrix.

What does NOT survive is the reading first proposed for it. Three statements
were wrong and are corrected in the paper:

* the marginal score does not over-credit a location entering a correlated
  regime: relative to the exact criterion it UNDER-credits the right regime by
  `Delta_i(h)` in expectation, and over-credits the wrong ones;
* the discarded term is not the same for every regime when the covariance
  parameters are shared: `q_ih` depends on where the members of the regime are;
* adding `Delta_i(h)` to the score is not a correction. It is the expected gain
  of the exact score under the hypothesis that the location belongs to the
  regime, granted whether or not the hypothesis holds: a reward for proximity.
  With the parameters known it helps where regimes are spatially compact and
  harms where they are mixed and at their boundaries.

The object that does pay the price is the CONDITIONAL score, the density of a
location's residuals given those of the other members of the regime, which is
the exact change of the block likelihood when the location joins. With the
parameters known it halves the misassignment when the regimes are independent
sub-networks, as the model assumes, and it sees covariance-only differences to
which the marginal score is blind (chance level there). But when one error
field spans the regimes it reads co-movement as membership and loses exactly
where proximity misleads -- mixed regimes and boundaries -- while the marginal
score is immune. And with the latent path smoothed rather than known, the
realized gap of the Lemma is close to zero: the latent process absorbs the
common part of the spatial error.

* `SCSTEM_Estimation()` gains `score = c("marginal", "conditional",
  "corrected")`. The default is unchanged. The new internal
  `scstem_cond_scores()` computes the conditional score and `Delta_i(h)` for
  every location against one regime, grouping the periods by the set of
  members observed, and is pinned by tests against a brute-force computation
  from the joint Gaussian density, with and without missing data.
* The fit returns two diagnostics of the label step: `objective_before` in
  `obj_trace`, the objective at the new parameters and the old labels, which
  isolates what the label step and the parameter step each did; and
  `score_last`, the score matrix of the last label step.
* The four propositions of the paper -- ascent and finite termination of the
  label step and non-monotonicity of the alternation; non-identifiability of
  the partition under S0; the two thresholds of the penalty; the gain of a
  spurious regime, `d/pi` in the canonical case -- are proved in the paper and
  checked by `dev/replication/check-theory.R`.

**A bug in the bootstrap.** `SCSTEM_Bootstrap()` rebuilt the arguments of its
refits from a list that did not carry `alpha`, `lambda`, `penalize`,
`lambda_scale`, `lambda_by`, `latent` or `spatial`, so a fit with a ridge on
the coefficients was resampled WITHOUT the ridge: the bootstrap described a
different estimator from the one reported. The fit now stores those settings,
and the score, and the bootstrap passes them on.

### 2026-09-26

**The simulation study and the application become standalone replication
material, in `dev/replication/`, disconnected from the package.** For the
paper there will be a replication package of its own, separate from the CRAN
package, so the scripts can no longer lean on the repository: they must run
from whatever folder they sit in, with nothing but an installed Stem.

| script | replaces |
|---|---|
| `run-simulations.R` | `07-simulation.R`, `14-bootstrap-coverage.R`, `design.R`, `06-dgp.R`, and the retired `02-simulation.R` |
| `analyse-simulations.R` | `09-overlap-figures.R`, `12-dgp-illustration.R`, `13-dgp-table.R`, `15-design-figures.R`, and the simulation half of `03-figures.R` |
| `run-application.R` | `01-application.R` |

* **One script runs the study, one analyses it.** `run-simulations.R` carries
  the generator, the design and the driver; `--coverage` runs the bootstrap
  experiment instead of the Monte Carlo, `--dry` prices a run without fitting
  anything. `analyse-simulations.R` draws the design and turns the results
  into tables and figures; it takes the generator and the design from the
  runner beside it, in a definitions-only mode, so the two cannot disagree.
* **The design is the SETUP block** near the top of the runner: levels,
  reference cell, blocks. Every factor can also be restricted from the command
  line with `--only_<factor>`.
* **Each script writes beside itself**, in `results/`, `output/` and
  `application/`, so a folder copied to another machine is a working unit.
* **Stem is checked on its features, not its version number.** The development
  builds all report 2.0.0, and the first standalone run proved the point: the
  Stem installed in this machine's library was an old build, without the
  fitted-signal functions and with the planar neighbour graph, and every
  replication failed. The runner now verifies what it needs and installs the
  current version from GitHub when it is missing or stale. `SIM_STEM_REF`
  should be pinned to a commit once the study is run.
* **A failed replication says why**, on the console and in a new `error` column.
  It used to leave a row of NA and nothing else, which on a long run on another
  machine is a day lost.

**The estimator is held to the generator's size floor.** Probing the first cell
showed that at `n = 20` a `k = 4` fit, with regimes of four and five locations,
took 364 s and 125 s where `k = 3` took one second, and one replication took 658
s instead of two. Those regimes are below `N_MIN = 6`, the floor the design
itself sets because a range is not identified on fewer locations. The runner
now passes `min_cluster_size = N_MIN` to the estimator, so the package refuses
such a `k` at once, the grid records it among the failed configurations and the
selection rule never sees it. The same replication now takes two seconds.

**The Drive.** `dev/sync-gdrive.R` delivers `dev/replication` to a folder of its
own beside the mirror, `STEM_Cameletti/SC-STEM-replication/`, with the opposite
rules: nothing is ever deleted there, and a script edited on the Drive is not
overwritten -- the repository version is written beside it as
`<name>.from-repo` and the conflict is reported. The mirror itself protects
`dev/replication/{results,output,application}/` and `dev/paper/cache/` from its
deletion pass.

`dev/paper` keeps only the development diagnostics (`03`, `04`, `05`, `08`,
`10`, `11`) and `00-setup.R`, which now diverts their outputs out of the mirror
when they are run from it.

### 2026-09-25 (fourth entry)

**The intermediate overlap is rounded to `omega = 0.70`, and the design is
written into the paper.** Since the separation is no longer proportional to
`omega`, a round value of `omega` is not a round value of anything that matters,
and the question is what a rounding costs on the scale that does. Measured at
`K = 3` by the Bayes error of the optimal assignment to the nearest centre, and
by the Adjusted Rand Index that rule attains -- what *any* purely spatial
procedure could achieve:

| `omega` | separation | Bayes error | attainable ARI | crossing edges |
|---|---|---|---|---|
| 0 | 0.00 | 0.667 | 0.000 | 0.67 |
| 0.50 | 1.05 | 0.443 | 0.117 | 0.55 |
| 0.686 | 1.58 | 0.333 | 0.252 | 0.44 |
| 0.70 | 1.63 | 0.324 | 0.265 | 0.41 |
| 1 | 3.16 | 0.100 | 0.722 | 0.15 |

`0.70` is within a rounding error of `0.686` and keeps the middle level where it
was; `0.50` would more than halve the attainable index, from 0.252 to 0.117, and
move the intermediate cell two thirds of the way back to the null case. The
design takes `omega = 0.70`. At `omega = 0` the Bayes error is `1 - 1/K`, the
error of guessing, which is the definition of the null case, and the crossing
share of the graph is the same number.

**The simulation study is now described in the paper.** Section 3.1.1 carries
the fixed-total-variance construction, the separation it implies, the table of
the three levels and the residual dependence of the distance distribution on the
overlap; a new Section 3.1.5 sets out the neighbourhood graph -- why it is a
modelling choice with point-referenced data and therefore a factor, the
symmetrization, the metric, the intractable partition function, and the fact
that spatial proximity enters the model twice; Section 3.1.6 carries the factors
with their levels, the five blocks with their cell counts, the evaluation
metrics, why the station is the unit of the error measures, and what the four
output records contain. The supplement's constants, pseudocode and regime-size
tables follow. Both documents compile, 25 and 11 pages.

### 2026-09-25 (third entry)

**The overlap parameter of the generator was doing two things at once.** In
`dgp_locations()` the within-regime variance was held at `nu_sp = 0.4` while the
centres were placed `2*omega` apart, so the physical extent of the network grew
with the separation. Measured over 20 draws at `n = 200`, `K = 3`: the spatial
standard deviation went from 63 km at `omega = 0` to 103 km at `omega = 1`, and
the median pairwise distance from 105 km to 180 km, against a true correlation
range of 123 km. So `omega` separated the regimes *and* enlarged the map, and an
effect attributed to the separation carried a share of the second.

The second effect is not cosmetic. At `omega = 0` the network is a blob smaller
than the range, the exponential decay is barely resolved over the observed
distances, and `theta` is close to unidentified. The estimator then crawls: on a
single *pooled* fit at `n = 20`, `T = 60` -- no clustering, no graph -- 0.4 s at
`omega = 1` against **more than ten minutes at `omega = 0` without converging**.
Every `omega = 0` cell of the study was effectively unaffordable, and nothing in
the timings said so, because they had been taken at one overlap.

* `dgp_locations()` now holds the TOTAL variance of a coordinate fixed,
  `nu_sp(K, omega) = NU_TOT - Var(mu | K, omega)`, with
  `NU_TOT = 0.4 + 2/3`: the total the design used to have at its most separated
  cell, so that cell is unchanged and the smaller overlaps are the ones that
  widen. The footprint is now 102-104 km of spatial standard deviation in every
  cell, at every `K` and every `omega`, and the pooled fit at `omega = 0` takes
  0.7 s.
* `K = 1` needs no special case any more. It has no between-centre variance, so
  it takes the whole of `NU_TOT` and covers the same area as everything else;
  the `k_ref` argument is gone.
* The separation is no longer proportional to `omega`, so the levels of the
  design are chosen with the new `dgp_omega_for()`: `omega = 0, 0.686, 1` give
  0, 1.58 and 3.16 within-regime standard deviations, the same three the design
  has always been read on. `dgp_nu_sp()` and `dgp_separation()` report what a
  cell actually is, and the scripts use them instead of recomputing
  `2*omega/sqrt(NU_SP)` by hand.
* `09-overlap-figures.R` keeps illustrating the design AS PUBLISHED, with the
  within-cluster variance fixed, which is recovered from the new signature by
  giving each cell the total its centres imply. Its third panel was already
  making this exact point about the absolute scale.

One residual dependence is worth stating in the paper rather than designing
away: at a constant footprint, `omega` still changes the *distribution* of
pairwise distances, from unimodal at `omega = 0` to short-within-and-long-between
at `omega = 1`, and a spread of lags identifies a range better than a single
scale. The estimated range at `n = 20` was 10, 42 and 75 km against a true
123 km. That is what clustered locations mean, not an artefact of the
parameterisation, but it should be said out loud.

### 2026-09-25 (second entry)

**The neighbourhood graph of the Potts penalty was built with a different
metric from the covariance.** `scstem_neighbors()` passed the coordinates to
`spdep::knearneigh()` without `longlat`, so on geographic coordinates the
nearest neighbours were those of a *planar* metric on raw degrees, in which one
degree of longitude and one of latitude count the same. They do not: at 45.5
degrees of latitude a degree of longitude is 78 km against 110.6 km, so a
north-south pair at a given separation in kilometres looked closer than the
east-west pair at the same separation, and the graph systematically preferred
north-south neighbours. Meanwhile `Sigma_e` was built from
`geodist::geodist(measure = "geodesic")`. The two spatial ingredients of the
model therefore disagreed with each other about which locations are close.

Measured on 300 points over a Po Valley box, at `knn = 5`: the planar graph
shares 78.5% of its edges with the geodesic one, and only 40.8% of its edges are
more east-west than north-south against 50.9% for the correct metric.

* `scstem_neighbors()` gains a `distance` argument and passes
  `longlat = TRUE` when it is `"geo"`. `SCSTEM_Estimation()` hands it the same
  `distance` it uses for the covariance, so the two can no longer diverge.
  Recovering 99.0% of the geodesic graph's edges, against 78.5% before.
* The coordinate columns are named `lon`/`lat` before the call, which is the
  convention the package follows throughout and which `spdep` otherwise has to
  assume out loud.

The default is unchanged in name only: `distance = "geo"` was already the
default of `SCSTEM_Estimation()`, so a fit that did not set it explicitly now
gets a different -- correct -- graph and may return a slightly different
partition. Nothing was released with the old behaviour.

**The design of the simulation study moves into `dev/paper/design.R`.** It was
spread between `dgp_dims()` in the generator and an `expand.grid()` in the
driver, with a third copy in the table script, which is three places to keep in
agreement and two too many. The new file holds the factors, their levels, the
reference cell and the blocks, and the driver and the table read it. It also
prices whatever is written in it from the timings measured on this machine, so
the budget of a change can be read before a run is launched:

```
Rscript dev/paper/design.R
```

The study is deliberately not a full factorial: crossing every factor with every
other is 2700 cells at `K = 3` alone, most of them answering no question. It is a
core factorial in the three factors that interact -- how many locations, how
separated the regimes, how long the series -- plus one-factor-at-a-time margins
around a reference cell for the factors that are there to show the results do not
turn on them. 240 cells, 237 of them feasible, 206 core-hours at M = 100.

`07-simulation.R` loses its per-factor command-line options, which redefined the
design, and gains `--blocks` and a `--only_<factor>` for each factor, which
select from it.

**Three figures for the design**, from `dev/paper/15-design-figures.R`:
`fig_design_omega.pdf` (the point configuration as the overlap varies, at every
network size), `fig_design_knn.pdf` (the graph the penalty lives on, as `knn` and
the overlap vary) and `fig_design_graph.pdf` (what the graph does in numbers).
The last one earns its place: it shows that `knn` multiplies the number of edges,
and so the weight the penalty carries at a given `phi`, but from n = 100 on it
leaves unchanged the share of edges that cross a true regime boundary -- which is
what decides whether the penalty helps. That share runs from 1 - 1/K, the value
of a graph carrying no information about the partition, at `omega = 0`, to about
0.13 at `omega = 1`.

### 2026-09-25

**`STEM_Fitted()` is split in two, because it was one function doing two jobs
under a name that belongs to only one of them.** In statistics and econometrics
*fitted values* means what the model predicts. What the function returned was
`E[z | observed]`: the imputation of the missing response, which wherever the
response *was* observed returns the observation itself, because the conditional
expectation of something already seen is that thing. On a complete record it was
therefore the data, and scoring it against anything measured nothing -- a trap
that caught the simulation driver, where the clustered and the pooled model came
out with root mean squared errors identical to fifteen digits.

* `STEM_Signal()` and `SCSTEM_Signal()` are new and carry the conventional
  meaning: the conditional mean given the latent path,
  `muhat_ti = x_ti' betahat + K_i yhat_t`, the regression surface plus the
  latent process, without the measurement error. This is what a fit is scored
  on.
* `STEM_Complete()` and `SCSTEM_Complete()` are the former `*_Fitted()`,
  renamed to say what they do. Use them to impute, not to fit.

No alias is kept for `*_Fitted()`: 2.0.0 is unreleased and the pair was added
during this development cycle, so the name never reached a user. The help pages
of each of the four point at the other, because the distinction is exactly the
one that is easy to get wrong.

**The simulation driver.** Every output now carries a primary key -- `cell`, a
string built from the design factors, together with `rep` -- so that the four
files join on two columns and on nothing else. The key is a string rather than
the tuple of factors because the overlap is a double: `2/3` does not survive a
round trip through a CSV exactly, and a join on a floating-point column is a
defect waiting to happen. The factors are kept beside it for filtering.

The driver also records the per-station error measures for every replication and
the full per-observation record for a configurable few, takes the neighbourhood
size `knn` as a factor of the design rather than a setting, and accepts
`--rep_from` / `--rep_to` so that a long study can be cut into blocks run on
different machines.

### 2026-09-07 (later)

**`SCSTEM_Estim()` is now `SCSTEM_Estimation()`.** The package exports two
estimation engines and they were called by names of different shapes. Both stay
exported -- the release brings the package back after its time off CRAN with a
wider interface, and deprecating the historical `STEM_Estimation()` would buy
nothing -- but they now read as a pair, with `STEM_Fit()` above them as the
entry point. No alias is kept: 2.0.0 has not been released, so `SCSTEM_Estim`
never existed for a user. The S3 class and its `print` method were renamed with
the function, and every reference in the sources, the tests, the vignettes, the
function map, the README and the development scripts was updated with it.

**`lambda` is now dimensionless, and only half the penalty is rescaled.** The
first version measured `lambda` against the scaled normal equations, which
leaves it carrying the units of the response for the L1 part: the same `lambda`
then meant different things at different error variances and, in the clustered
model, in different regimes. The obvious repair -- rescale the whole penalty by
the size of the gradient -- is wrong, and instructively so. After the columns of
the design are scaled to put the diagonal of `M` at one, the L2 part is
*already* unit-free: it multiplies a curvature that is one, so it shrinks by
`1/(1+lambda)` whatever the units. Rescaling it as well would make a ridge
depend on the error variance when it need not, to the tune of `1e-1` in the
coefficients on a collinear design.

So `lambda_scale = "relative"`, the new default, rescales the **L1 part only**,
by the largest partial gradient once the unpenalized coordinates are profiled
out -- the `lambda_max` of the lasso path. Every penalized coefficient is then
exactly zero once `lambda * alpha >= 1`, and `lambda` in `(0, 1]` traverses the
whole path at `alpha = 1`, which is the convention of `glmnet`.
`lambda_scale = "absolute"` keeps the previous meaning and is what to use when
comparing against an external implementation. The test suite asserts the
invariance directly: scaling the response by a constant scales the coefficients
by the same constant, and scaling the GLS information leaves them alone, for
every `alpha`.

**`latent = FALSE` and `spatial = FALSE`.** The first sets the loading matrix to
zero, so the state contributes nothing to the measurement equation, and holds
`G`, `Sigma_eta` and `m0` -- unidentified without it -- at their input values
rather than letting them chase a flat likelihood. The second replaces the
exponential correlation by the identity, so `Sigma_e` is a single variance, and
skips the Newton-Raphson step, which would otherwise be estimating a range that
has nothing to estimate. Together, and with `regularization = 0`, they reduce
the model **exactly** to penalized linear regression: the agreement with an
elastic net computed directly on `X'X` and `X'y` is between `1e-14` and `1e-11`
for the ridge, the lasso, the elastic net and the unpenalized case alike.

That last condition is not a detail. The small ridge the package adds to every
matrix it inverts is comparable, on a collinear design, with the smallest
eigenvalue of `X'X` in the offending direction, and it displaces the
coefficients by about `1e-2`. A tolerance that is invisible in a well-conditioned
problem is not invisible in the problem a penalty exists for.

**`lambda_by = "size"`.** With `k > 1` the penalty can be spread over the
regimes in two ways. `"common"`, the default, gives every regime the same
`lambda`, which under the relative scale already means the same *proportional*
shrinkage, since the reference is computed inside each regime. `"size"` sets
`lambda_g = lambda * nbar / n_g` with `nbar = d/k`, shrinking a regime of half
the average size twice as hard, on the argument that its coefficients are
noisier than proportionality alone accounts for. Either way `lambda` stays one
hyperparameter; genuinely cluster-specific `(alpha_g, lambda_g)` is a different
model and is deliberately not offered.

**Where the theory lives.** `dev/regularization/stem-elastic-net.tex` is
repositioned as the draft of a second, separate paper: cluster-specific elastic
net, with penalized linear regression as a special case. It gains a section
deriving why the four hyperparameters are not separable -- the optimal `lambda`
grows with `k` at rate `k`, the selected `k` grows with `lambda` because
shrinkage lowers the effective parameter count, and the automatic scale of the
Potts penalty moves with `lambda` because it is calibrated on a spread that
shrinkage narrows -- and a section on why the criterion must use the
*unpenalized* likelihood at the penalized estimate rather than the penalized
objective, which would charge for the penalty twice. What the manuscript for
Metron uses is a strict subset: a pooled ridge in closed form.

### 2026-09-07

**Ridge, lasso and elastic net on the regression coefficients.** `STEM_Fit()`
is a single entry point for the four estimators the package now provides, and
what runs is decided by two arguments and nothing else: `k = 1` or `k > 1`
chooses between the pooled and the clustered model, `lambda = 0` or `lambda > 0`
between the ordinary and the penalized one. The defaults `k = 1`, `alpha = 0`,
`lambda = 0` reproduce `STEM_Estimation()` bit for bit, which the test suite
asserts with `expect_identical()` on the whole parameter vector and on the
log-likelihood.

In the parameterization of `glmnet` the objective subtracts
`lambda { alpha ||D beta||_1 + (1-alpha)/2 beta' D beta }`, so `alpha = 0` is
ridge, `alpha = 1` the lasso and anything between the elastic net. Only the
regression coefficients are penalized: the variance components, the range, the
transition matrix and the initial state stay at their maximum likelihood values.

Three things made this cheap and one made it subtle.

*The E-step does not change.* The penalty is a function of `beta` alone and does
not involve the latent states, so it passes through the conditional expectation
unchanged and the usual Jensen argument applies verbatim to `Q - pen`. The EM
algorithm therefore still increases the PENALIZED observed-data log-likelihood
at every iteration, and the Kalman filter and smoother -- the expensive part --
are untouched. Only the point returned by the M-step differs.

*The M-step keeps a closed form under a ridge.* `beta = (M + lambda D)^{-1} v`
with `M` and `v` already assembled: one addition on a diagonal.

*The lasso and the elastic net are solved exactly.* The M-step objective is a
quadratic plus a separable penalty, which is the setting in which cyclic
coordinate descent converges to the global maximizer, so each step is soft
thresholding and the result is a genuine EM step rather than a generalized one.
The implementation is checked against the KKT conditions of the problem it
claims to solve, not against another implementation.

*The metric is the subtle part.* `M = sum_t X_t' Sigma_e^{-1} X_t` is a
GENERALIZED least squares cross-product: it carries an estimated spatial
covariance. Standardizing the columns of X in the usual Euclidean sense --- what
every lasso implementation does --- is therefore not the normalization that
makes `lambda` mean the same thing for every covariate here. The scaling that
does is the one putting the diagonal of `M` at one, and that is what the code
applies internally, returning the coefficients on the original scale. In the
clustered model this is also what makes a single `lambda` comparable across
regimes of different sizes and different error covariances.

The intercept is left unpenalized by default, and not merely by convention: the
model already carries a latent process whose initial mean absorbs the level, so
shrinking the intercept would not shrink "the level" but move it into `m0` at a
rate depending on `G`.

**The information criteria count effective parameters.** With a penalty in force
the nominal `r` overstates the flexibility of the fit. `SCSTEM_Estimation()` now sums
the effective count each regime reports: `tr(M (M + lambda D)^{-1})` for a
ridge, the number of active coefficients for a lasso, and the corresponding
trace on the active set for an elastic net. With `lambda = 0` every regime
returns `r` and the count reduces to the one the package has always used.

The derivation -- the penalized EM and its monotonicity, both forms of the
M-step, the GLS metric, the degrees of freedom, and what is still open -- is in
`dev/regularization/stem-elastic-net.tex`. It is deliberately outside the
package: the material is a study of its own and not documentation of software.

**A vignette on how the engine is computed.** `vignette("computational-notes")`
explains what a cost of `O(T d^3)` means, why inverting a `d x d` matrix costs
`O(d^3)`, how the Woodbury identity and the matrix determinant lemma bring a
forward pass down to `O(d^3 + T d^2)`, under exactly which assumptions, and why
evaluating the log-density directly is a correctness fix rather than a speed-up.
This material is deliberately kept out of the paper.

**The simulation design.** The overlap parameter is now called `omega`
throughout rather than `d`, which already means the number of locations; the
design gains a pooled cell `K = 1` and an unbalanced allocation of the regime
sizes; `T = 730` is dropped, which removes 42 per cent of the running time of
the whole design. Every factor of `dev/paper/07-simulation.R` is a command-line
option, so a run can be narrowed to the margins actually reported, and every
paper script now resolves the repository from its own location and reads its
output directories from environment variables, so they run from any working
directory and on a virtual machine.

### 2026-09-06 (later)

**The Newton-Raphson step no longer rebuilds what it already has.** With the
filter and the M-step down to `O(d^3 + T d^2)`, the inner Newton-Raphson loop
that updates `theta` and `log b` became the next thing worth looking at. Four
separate redundancies were in it.

*The exponential kernel.* `Sigmastar.exp()`, `d1_Sigmastar_logtheta.exp()` and
`d2_Sigmastar_logtheta.exp()` each evaluated `exp(-exp(logtheta) * dist)`
independently, and the two derivatives each formed `exp(logtheta) * dist`
twice. That is three `d x d` elementwise exponentials and four scalings per
inner iteration, of which one and one are needed. All three functions now take
the kernel as an optional argument and `kalman()` computes it once. A covariance
function that does not accept it still builds its own, so a user-supplied
`cov.spat` keeps working.

*The derivative with respect to `log b` is not a matrix.* It is `exp(logb)`
times the identity. It was being materialised as a dense `d x d` matrix --
160,000 doubles at `d = 400` to carry 400 of them -- and then multiplied into
`X^{-1}`, an `O(d^3)` matrix product standing in for an `O(d^2)` scaling.
`d1_Sigmastar_logb.exp()` and `d2_Sigmastar_logb.exp()` now return the
multiplier itself, and `d1_Q()`, `d2_Q()` and `d12_Q()` recognise a scalar
argument as a multiple of the identity through the two new helpers
`stem_xprod()` and `stem_xtrace()`.

*The shared products.* `X^{-1} dSigma/dlogtheta` entered three of the five
derivative evaluations and was formed three times; the same for `log b`. The
five functions now accept the products, and `kalman()` forms each once. Between
this and the previous point the `d x d` matrix products per inner iteration go
from twelve to five.

*Two factorisations where one does.* `Q_function_addendo1()` called `solve()`,
which factorises its argument, and then `determinant()`, which factorises it
again. One Cholesky now gives both: the log determinant is read off the diagonal
of the factor and `chol2inv()` inverts from it, with a fallback for a matrix
that has lost positive definiteness numerically -- which the grid search can
produce at extreme parameter values. This function is called a hundred times
whenever that grid search fires, so the saving is largest exactly where the
cost was worst.

**The smoother no longer computes a quantity nobody reads.** `smoothing()`
built, row by row, the `n x d` matrix of fitted means `F' m_t` and returned it
as `ss$mu`. Nothing in the package ever looked at it: `kalman()` forms its own
fitted values from the smoothed states. At `n = 730` and `d = 400` it allocated
and discarded 2.3 megabytes on every EM iteration.

Each changed expression was checked against the one it replaces at random
inputs: the covariance and both derivatives with respect to `log(theta)` agree
exactly, and the five derivatives of the objective agree to between 0 and
1.3e-15 in relative terms, on both the path `kalman()` takes and the path an
outside caller takes.

### 2026-09-06

**The Kalman filter no longer forms or inverts the `d x d` predictive
covariance.** The forward pass used to build

```
Q_t = Z R_t Z' + Sigma_e
```

at every time point and hand it to `solve()`. That is `O(d^3)` per step and
`O(T d^3)` over a pass, and on anything past a few dozen locations it was
essentially the whole running time of the package: SC-STEM fits one STEM model
per regime, on a grid of `(k, phi)`, on every bootstrap draw, so the cubic term
is paid thousands of times.

Two facts make it avoidable. `Sigma_e` does not depend on `t` -- it is rebuilt
once per EM iteration -- and `Z R_t Z'` has rank `p`, the dimension of the
latent state, which is one in the default specification. So the Woodbury
identity gives

```
Q_t^-1 = Sigma_e^-1 - U (R_t^-1 + Z' Sigma_e^-1 Z)^-1 U',   U = Sigma_e^-1 Z
```

and the matrix determinant lemma gives

```
log|Q_t| = log|Sigma_e| + log|R_t| + log|R_t^-1 + Z' Sigma_e^-1 Z| .
```

`U`, `Z' Sigma_e^-1 Z`, the Cholesky factor of `Sigma_e` and `log|Sigma_e|` are
constants of the pass, so `filtering()` computes them once and `filterstep()`
completes each step with `p x p` algebra plus one triangular solve. A pass costs
`O(d^3 + T d^2)`. The recursion never needs `Q_t` itself: the gain enters only
through `Z' Q^-1 Z` and `Z' Q^-1 e`, which are `p x p` and `p x 1`. With missing
values the constants depend on which rows are observed, so they are cached per
distinct missingness pattern; in a monitoring network a station is out of
service for a stretch of consecutive days, so the cache almost always hits.

**The log density is evaluated directly.** `log(mvtnorm::dmvnorm(...))` was not
only slower -- `dmvnorm` re-checks the symmetry of `Q` with `all.equal` at every
step, about a seventh of the running time on a 36-station network -- but wrong
on a large network. The Gaussian density on `d` observations is of order
`exp(-d)`, so past roughly 300 locations it falls below the smallest
representable double and its logarithm is `-Inf`. The previous code could not
fit a 400-station network at all; this is not a speed-up but a correctness fix.

**The M-step accumulations are matrix products rather than loops.** Three loops
built lists of `T` matrices and then summed them. Because the loading matrix
does not depend on time,

```
SUM_t { Z C_t Z' + r_t r_t' } = Z (SUM_t C_t) Z' + R'R
```

with `R` the `T x d` matrix of residuals, so one cross-product replaces `T`
outer products and the `T` matrices of size `d x d` that had to be held at once.
The three sums of outer products of the smoothed states are cross-products of
the `T x p` matrix of those states; the two accumulations entering the update of
`beta` are cross-products of the design blocks stacked by period. Every trace of
a matrix product is now evaluated as `tr(AB) = SUM_ij A_ij B_ji`, which is
`O(d^2)` instead of forming the product, and the derivative helpers `d1_Q`,
`d2_Q`, `d12_Q` and `Q_function_addendo1` receive the inverse of the scaled
covariance and its product with `B` from the caller instead of recomputing
`solve()` up to seven times each. `Q_function_addendo1` also takes the log
determinant through `determinant(., logarithm = TRUE)`, which does not overflow
at large `d`.

Measured against the state of the code at commit `7d0b094`, ten EM iterations
on a synthetic network generated from the model itself:

```
   d     T | before (s)  after (s)  speedup | max rel. diff
  36   365 |        2.1        0.9     2.3x | 1.18e-13
  60   365 |        2.1        0.6     3.3x | 7.96e-14
 100   365 |        5.9        1.2     5.0x | 1.17e-13
 200   365 |       15.5        1.6     9.8x | 2.85e-13
 400   200 |     FAILED       10.8        - | reference underflows to -Inf
```

Both changes are algebraic identities, so the point of the table is the last
column, not the speed-up: the estimates coincide to floating point. The
benchmark is `dev/paper/10-engine-benchmark.R`, which reconstructs the old
sources from git and runs the two side by side.

**Consequences elsewhere.** `sumMatrices()` was a two-line wrapper around
`Reduce("+", .)` called only from the three loops that are gone, so it has been
removed; it was internal and undocumented, so nothing user-facing changes.
`mvtnorm` is no longer used by any function in `R/` and has moved from `Imports`
to `Suggests`, where the test suite still uses it as an independent reference
implementation of the likelihood. The function map gains the missing-data
helpers and the two `*_Fitted()` functions, loses `sumMatrices`, and its
generator now refuses to emit a diagram whose edges run through boxes -- which
caught two long-standing routing defects in the diagram itself.

### 2026-09-04

**Missing values in the response are now supported.**
Following Durbin and Koopman (2012, 2nd ed.), Sections 2.7 and 4.10. At each
time point the measurement equation is restricted to the locations actually
observed, through the selection matrix whose rows are a subset of the rows of
the identity: `z`, the loading matrix and the measurement covariance are all
pre-multiplied by it, and the Kalman recursion proceeds on an observation vector
whose dimension varies over time. A time point at which nothing is observed
contributes no update and no likelihood term, which is their `Z_t = 0` device.
The likelihood therefore stays the exact likelihood of the observed data by the
prediction-error decomposition. The backward smoothing recursions needed no
change at all: they read the filtered moments and the transition, never the
data, and the same is true of `B_function()`, `cov_lagone()` and the M-step
updates of `G`, `Sigmaeta` and `m0`, which are functions of the smoothed states
alone.

The work is in the M-step, and this is where a naive implementation goes wrong.
EM maximizes the expected COMPLETE-data log-likelihood, so the sufficient
statistics have to be completed rather than truncated. A missing value enters
through its conditional expectation given everything observed, and its
conditional variance is added back as a correction term; the divisor of the
`sigma2omega` update stays the complete-data count `n*d`, not the number of
observed values.

The subtlety is what that conditional expectation is. Durbin and Koopman remark
that a missing element can be estimated by the corresponding element of
`Z_t yhat_t`, which is exact when the measurement covariance is diagonal - and
that is not this model. Here `Sigma_e` couples the locations, so the conditional
mean of a missing observation is its signal PLUS the part of the measurement
error predicted from the neighbors observed at the same instant, by the usual
Gaussian conditioning. It is the same algebra as kriging, applied at a fixed
time point. Ignoring that term would bias the variance components downwards, by
charging to noise a residual the model can explain spatially.

The rest of the chain follows. `STEM_Kriging()` conditions on the sub-vector
observed at the chosen time point instead of the whole of `z_t`, and falls back
on the unconditional mean where nothing was observed. `scstem_loglike_i()` scores
a location on the time points at which it was observed, with its own `T_i`;
since a location's missingness pattern does not depend on the cluster it is
being scored against, the scores stay comparable across clusters. Both
bootstraps reimpose the observed pattern of gaps on every replicate, since a
replicate with a complete response would understate the uncertainty of a fit
obtained from an incomplete one. `STEM_Model()` now rejects missing values only
in the covariates and the coordinates, and additionally rejects a location with
no observed value at all.

Verified on four counts. On complete data the estimates are **bit-for-bit
identical** to those of the previous commit - maximum absolute difference
exactly zero across all parameters and the log-likelihood, checked by running
the working tree and `HEAD` side by side on the `pm10` example. With 15 percent
of the response blanked at random the fit runs and the estimates move little.
A time point at which nothing is observed is handled. And on data simulated from
a known truth, a 20 percent MCAR gap recovers the true parameters as accurately
as the complete data do, which is what an unbiased treatment should look like.
`tests/testthat/test-missing.R` adds 35 assertions, including the algebra of the
conditional blocks against direct Gaussian conditioning.

One bug found and fixed while testing: pasting an empty observed-index gave the
same cache key as a complete time point, so a fully missing row was served the
complete-data shortcut and the accumulation failed.

**The conceptual map is now generated, in both formats.**
`inst/scripts/make-function-map.R` holds the map as tables of nodes and edges
and emits both `inst/extdata/STEM_function_map.svg` and the new
`inst/extdata/STEM_function_map.pdf`, so the two cannot drift apart. It ships
with the package, so a user can rerun it. The layout stays explicit rather than
computed by a graph-drawing algorithm, because the point of the map is the
reading imposed on the package - four bands, from the data up to the clustered
layer - which an automatic layout would not reproduce.

`STEM_function_map_original.pdf` has been removed. It was the map of version
1.0: dot-separated names (`Stem.Model`, `Stem.Estimation`) that no longer exist,
and no SC-STEM layer at all. It documented a package that is not this one.

The regeneration also caught a name the diagram had missed: it still said
`scstem_neighbours`, whereas the function was renamed `scstem_neighbors` in the
American-English pass.

**Documented the statistical features and the scope of the model.**
A new section of the README, and a matching block in `?SCSTEM_Estimation`, state what
the family covers and what it does not: Gaussian response only; univariate
response, with the multivariate part being the latent state; exponential spatial
covariance, hence isotropic and stationary within a regime, the partition being
the only source of non-stationarity across the domain; point-referenced
locations, fixed over time; a discrete and regularly spaced time index; a known
loading matrix, common across regimes; `C0` fixed while `m0` is estimated.

Writing that section is what surfaced the missing-value question. At the time it
recorded that missing values were not supported anywhere - `STEM_Model()`
rejected them and the Kalman recursion had no partial-observation branch - and
noted that the response case was the standard state-space treatment while the
covariate case would need a stochastic E-step. The response half was implemented
straight afterwards, in the entry above; the table now reflects that, and the
restriction stands only for the covariates and the coordinates.

**Package logo.**
The hexagon now carries the three ideas the package is about: a relief
silhouette for space, a node-and-edge network split by color into two regimes,
and one time series per regime along the base. The faint gray edges are the ones
the partition cuts, which is what the spatial penalty pays for.

The wordmark was the hard part, and the difficulty was structural rather than
typographic. In a pointy-top hexagon the width collapses below the widest band -
at the baseline where the word had been sitting the shape is only about 24 units
across - so the word could not grow without being clipped. Moving it up into the
full-width band took it from 15 units to 26. Two further adjustments came out of
looking at the render: the series were made angular rather than wavy, because
smooth curves under a mountain range read as water, and the ridge was broadened,
because sharp peaks and sharp series were competing for the same reading.

`man/figures/logo.svg` is the scalable master and `man/figures/logo.png` the
rendered artifact, both produced from one set of coordinates. The README points
at the PNG: the SVG carries live text, so its wordmark depends on whichever
serif the viewer has, and GitHub sanitizes SVG attributes before serving them,
possibly including the `textLength` that guarantees the word fits.
`dev/make-logo.R` redraws the PNG at any size. It uses base graphics, which can
only clip to a rectangle, so every element of the composition is laid out to
fall inside the hexagon and no clipping is needed.

**American English throughout, and ASCII only.**
Every source, vignette and document converted to American spelling: `-ize`
rather than `-ise`, `neighbor` rather than `neighbour`, `modeling`, `center`,
`meter`, `analyze`. Non-ASCII characters removed as well - em and en dashes
folded to plain hyphens, curly quotes to straight ones, accented letters in the
references to their ASCII forms. Beyond consistency this keeps `R CMD check`
quiet about non-ASCII characters in the sources.

**Commit messages simplified.**
The narrative that used to live in commit subjects moved here. GitHub shows the
last commit message as the caption of every file it touched, and long subjects
made that listing unreadable.

**Fixed the README maths error.**
GitHub processes the markdown before handing the maths to MathJax, and its HTML
sanitizer read `<j}` in `\sum_{i<j}` as the start of a tag. What reached MathJax
was truncated, with a brace opened and never closed - exactly the reported
*"Extra open brace or missing close brace"*. This is why the error survived
three rewrites: every version carried `\sum_{i<j}`. It also explains why a
MathJax probe cleared the same expression: the probe fed the LaTeX straight to
MathJax and skipped the very stage that breaks it.

The penalty is now written as a sum over the edge set of the neighborhood graph,
which carries no character hostile to HTML, is the more natural way to state a
Potts term, and renders in every engine. The MathJax macro `\lt` would also have
worked on GitHub but is not standard LaTeX and would break the PDF manual, so
the edge-set notation is used in the README, in the SC-STEM vignette (GitHub
renders `.Rmd` files too) and in the `\deqn` of `?SCSTEM_Estimation`.

**README: parameter roles, cluster-wise equations, more references.**
Equations returned to LaTeX, in fenced math blocks. A new table gives the role
and the reading of every parameter: `beta`, `y_t`, `K`, `G`, `Sigma_eta`,
`sigma2eps` as the nugget, `sigma2omega` as the partial sill, `theta` as the
range with `1/theta` its characteristic length, `Sigma_e`, and the initial
state. The SC-STEM part is now written cluster-wise - a measurement equation, a
state equation and a covariance per regime, then the full parameter set
`Psi_k` - making explicit that a regime differs from another not only in level
but in how it behaves in time and space. It is motivated by spatial
heterogeneity as the second law of geography, against Tobler's first law which
`Sigma_e` already encodes. Added Goodchild (2004), Zhu and Turner (2022), and
Maranzano, Mattera and Sugasawa (2026, arXiv:2608.13638).

**Authorship.**
Francesco Caccia restored among the package authors: he contributed to the
package, which is a separate matter from the paper. Michela Cameletti's ORCID
added, and the author list in the README turned into a table with ORCID and
GitHub handles.

**Documentation split in two.**
`inst/extdata/STEM_model_verification.html` now does one thing only: map the
equations of the three reference papers onto the functions that implement them.
It no longer mentions defects or changes, so it stays valid as documentation of
the package rather than of one release. The development record moved to
`dev/code-changes-report.html`, which `.Rbuildignore` excludes, so a
changelog-shaped document does not ship to CRAN.

**Verified the implementation against the three reference papers.**
Read Fasso, Cameletti and Nicolis (2007), the GRASPA technical report of Fasso
and Cameletti (2007) and Fasso and Cameletti (2010), and checked the estimation
code against them equation by equation. Verified as matching: the scaled
covariance and the reparameterization to `log(gamma)`; the likelihood by
prediction-error decomposition; `Qtilde` and `W`; the closed-form updates of
`beta`, `sigma2omega`, `G`, `Sigma_eta` and `m0`; the Newton-Raphson step,
whose gradient and Hessian were re-derived by hand, including the mixed term
that vanishes because the two covariance parameters enter additively; the
kriging predictor; the parametric bootstrap; and the two joint convergence
criteria.

One defect found: in `filterstep` the univariate branch computed the predictive
mean as `f - X beta` where Equation (9) gives `mu_t = X_t beta + K y_t`, so the
regression term must be added. The branch is reached only when a fit has a
single location - essentially never for a pooled model over a network, but
reachable cluster by cluster. It produced a silently wrong log-likelihood, hence
wrong information criteria and a wrong assignment, rather than an error.

The check also confirmed the SC-STEM assignment score: its per-location variance
`sigma2eps + sigma2omega` is exactly `Sigma_e[i,i]` under Equation (7) of the
technical report.

**Google Drive mirror.**
`dev/sync-gdrive.R` mirrors the git-tracked files onto the Google Drive working
folder and deletes what is no longer tracked, so the two copies cannot drift. It
refuses to run on a dirty tree, so the mirror always corresponds to a real
commit, and writes a `MIRROR.txt` recording which one.

**Restored the introductory vignette.**
It had been written as `vignettes/Stem.Rmd` and then deleted by an
`rm -f vignettes/STEM.Rmd` aimed at the older placeholder: the Windows
filesystem is case-insensitive. Rewritten as `vignettes/getting-started.Rmd`, a
name that cannot collide with the package name under any case folding.

---

## 2.0.0 - 2026-09-03

The full release notes are in [`NEWS.md`](NEWS.md). In summary:

* **New SC-STEM layer**: `SCSTEM_Estimation()`, `SCSTEM_Infocrit()`,
  `SCSTEM_Select()`, `SCSTEM_Bootstrap()` and `SCSTEM_BootInference()`, with
  their classes and `print()` methods.
* **Algorithmic rewrite**: ICM label update, explicit penalized objective with
  cycle detection, minimum-size constraint, size-preserving swap pass,
  symmetrized neighborhood graph, data-driven penalty scale, final refit.
* **Four corrections to the cluster assignment score**, none of which raised an
  error and all of which changed the estimated partition.
* **Numerical robustness** of `kalman()` and `STEM_Estimation()`: pre-existing
  defects, harmless for a pooled fit and fatal for a clusterwise algorithm. On
  the `pm10` example they took the bootstrap from 4 usable draws out of 12 to
  12 out of 12.
* **Information criteria** on the exact log-likelihood of the final refit, with
  the correct parameter count and sample size.
* **CRAN compliance**: all five review comments addressed, plus the `pm10`
  dataset converted to `.rda`, the non-CRAN `SCDA` dependency removed, and
  `dplyr` and `sf` dropped as unused.
* **New material**: the `povalley` dataset, three vignettes, a `testthat` suite
  of 43 assertions, and a conceptual map of the package.

`R CMD check --as-cran`: `Status: OK`, no errors, warnings or notes.

---

## 1.0 - 2010

Original CRAN release by Michela Cameletti, subsequently archived.
