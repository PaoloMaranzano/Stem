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

### 2026-09-28

**Verification of the warm start, complete.** On 189 of the 190
cells-by-replications of the simulation design (19 scenario-variants x 10
replications; one replication of S3Seta-shr still running), grid
`phi = c(0, 0.025, 0.05, 0.1, 0.2)`, band `c(0.025, 0.2)`:

| rule, start of the phi = 0 fits | ARI at the true k | ARI = 1 | k correct | k = 2 | k = 4 | seconds |
|---|---|---|---|---|---|---|
| before (first-sweep scale), default | 0.868 | 64% | 76% | 15% | 9% | 49 |
| before, M | 0.834 | 26% | 60% | 36% | 4% | 46 |
| warm start, default | 0.900 | 72% | 81% | 8% | 11% | 27 |
| warm start, M | 0.955 | 86% | 92% | 6% | 2% | 24 |

With the warm start and the M start the independent variants are recovered
almost exactly (ARI 0.999, k correct 96%) and the shared ones at 0.927 and 89%;
the two cells that stay low, S1w-shr (0.58) and S3beta-shr (0.64), are the two
the separation measure of the design notes marks as hard (4.5 and 6.6
log-likelihood units per location). S2 selects k = 1 in every replication. The
warm-started penalized fits converge in one to three iterations, so a
replication of the grid is also about half as long. The selection rule picks
phi = 0.025, the bottom of the band, in 357 of 358 replications. The default
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
