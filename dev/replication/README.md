# SC-STEM: replication material

The scripts that produce the simulation study and the application of the paper.
This folder is **not part of the Stem package** and is never shipped with it:
Stem is a dependency like any other, and the scripts are deliberately
disconnected from its sources.

| script | what it does |
|---|---|
| `run-simulations.R` | runs the simulation study: the Monte Carlo, or the bootstrap coverage experiment |
| `analyse-simulations.R` | turns the results into the tables and figures of the paper, and draws the design |
| `run-application.R` | runs the application, staged and cached |

## Requirements

R, and nothing else to set up. Each script installs Stem from GitHub
(`PaoloMaranzano/Stem`) when it is missing, and installs it again whenever
GitHub holds a newer commit than the installed one, so the study always runs on
the current code. Offline, an installed Stem that carries what the scripts use
is accepted as it is. The check is on the commit and on the features, not on
the version number, because the development builds all report 2.0.0.

For the paper the reference installed (`SIM_STEM_REF` and `APP_STEM_REF`) is
to be pinned to the commit the results were produced with, or replaced by the
CRAN release.

## Running the simulations from RStudio

1. Open `run-simulations.R` in RStudio, from wherever it sits (the Google Drive
   folder included).
2. If you want, change the SETUP block near the top. In "4. THE RUN":
   - `mode`: `"run"` (the Monte Carlo), `"coverage"` (the bootstrap coverage
     experiment) or `"dry"` (print the design and its cost, run nothing);
   - `cores`: all the physical cores but one by default, or a number;
   - `rep_from`, `rep_to`: the replications to run.
3. Press **Source** (Ctrl+Shift+S).

The replications are handed to the cores one at a time, as each core frees up,
and every one prints a line with its time and an estimate of the time left. The
results go to `results/` beside the script.

The run can be stopped at any moment with the red Stop button and resumed by
pressing Source again: everything is keyed on (cell, replication), so what was
recorded is skipped. Keep the RStudio session open while it runs, and keep the
machine from going to sleep.

Do not run the same `tag` on two machines at once: they would write to the same
files. To share the work between machines, give each its own `tag` and its own
range of replications, and analyse them together with
`analyse-simulations.R --tag=a,b`.

If R on that machine uses a multithreaded BLAS (MKL, OpenBLAS), set
`OMP_NUM_THREADS=1` in `~/.Renviron` before starting, otherwise every core
tries to use every core.

## Running from a terminal

Every SETUP setting can also be given on the command line, which overrides it:

```
Rscript run-simulations.R --mode=dry         # the design and its cost
Rscript run-simulations.R --cores=8          # the whole study   -> results/
Rscript run-simulations.R --mode=coverage    # the coverage      -> results/
Rscript analyse-simulations.R                # tables, figures   -> output/
Rscript run-application.R                    # the application   -> application/
```

Slices of the design compose: `--blocks=`, `--rep_from=`, `--rep_to=`,
`--tag=`, `--nrep=`, `--only_n=`, `--only_TN=`, `--only_K=`, `--only_omega=`,
`--only_knn=`, `--only_scenario=`, `--only_balance=` and `--out=`.

`analyse-simulations.R` takes the generator and the design from
`run-simulations.R`, so the two have to sit in the same folder; that is also
what guarantees they cannot disagree.

## The design

It is the SETUP block near the top of `run-simulations.R`: the levels of every
factor, the reference cell, the blocks. Edit the numbers there, and run with
`mode = "dry"` to read the number of cells and the cost before committing to it.

The study is not a full factorial. It is a factorial in the three factors that
interact -- number of locations, separation of the regimes, length of the
series -- plus one-factor-at-a-time margins around a reference cell for the
factors that are there to show the results do not turn on them.

The design is generic: no value is calibrated on the application. The plane is
used with Euclidean distances in its own units, the response is standardized,
and every parameter has a reading of its own (see `dgp_base()` and `dgp_psi()`
in the runner). The earlier design, calibrated on the Po Valley network of the
application, is kept in `dev/archive/sim-design-povalley` of the repository.

The cost table of the runner was measured on that earlier design and is
pessimistic for the generic one, which fits faster; trust the time left printed
on the console over the "dry" estimate until it is measured again.

A replication that fails says why, both on the console and in the `error`
column of the results.

## The outputs

`results/<tag>.csv`, `<tag>-params.csv`, `<tag>-stations.csv` and the folder
`<tag>-obs/` share the primary key (cell, rep) and join on it and on nothing
else. What each holds is described at the top of `run-simulations.R`.
