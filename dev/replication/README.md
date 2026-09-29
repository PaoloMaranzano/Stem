# SC-STEM: replication material

The scripts that produce the simulation study and the application of the paper.
This folder is **not part of the Stem package** and is never shipped with it:
Stem is a dependency like any other, and the scripts are deliberately
disconnected from its sources.

| script | what it does |
|---|---|
| `run-simulations.R` | runs the simulation study |
| `analyse-simulations.R` | turns its results into tables (`.tex`, `.csv`) and figures (`.pdf`) in `output/` |
| `run-application.R` | runs the application, staged and cached |
| `run-fuels.R` | first version of the fuel-price application: Granger causality between gasoline and diesel at the pump, by province, on weekly changes; pooled STEM against SC-STEM and a pump-by-pump benchmark. Reads `App_FuelsITA/station_level.zip` beside it, writes to `fuels/` |

The design of the simulation study is described in full in
`simulation-design.tex`, in the Overleaf project of the paper. The analysis
takes the design from `run-simulations.R`, which must sit beside it, and can be
run while the study is in progress: each part produces what the results
support so far. The analysis of the previous design is archived in
`dev/archive/sim-design-generic`.

```
Rscript analyse-simulations.R                  # the results of tag "main"
Rscript analyse-simulations.R --tag=main,ref   # several tags, stacked
```

## Requirements

R, and nothing else to set up. `run-simulations.R` installs Stem from GitHub
at the commit it is pinned to (`SIM_STEM_REF`, currently 4325536), and
reinstalls it whenever the installed Stem is a different commit; the commit is
read from the installed DESCRIPTION file, without loading the package. Restart
R before the first run on a machine where Stem is already loaded in the
session. Offline, an installed Stem that carries what the script uses is
accepted as it is. `run-application.R` installs the current GitHub commit
(`APP_STEM_REF`), to be pinned in the same way before the results go into the
paper.

## Running the simulations from RStudio

1. Open `run-simulations.R` in RStudio, from wherever it sits (the Google Drive
   folder included).
2. If you want, change the SETUP block near the top:
   - "1. THE DESIGN": the scenario-variants (`SIM_SCENARIOS`), the reference
     level of every margin (`SIM_REFERENCE`: `n`, `TN`, `omega`, `balance`,
     `knn`) and the blocks (`SIM_BLOCKS`). Each block is the full factorial of
     the margins it lists, every other margin at its reference level, crossed
     with the scenario-variants; the design is the union of the blocks, and a
     cell two blocks share is run once. The current design is a core crossing
     `n`, `TN` and `omega`, plus blocks varying the balance, `knn` and a large
     `n` one at a time around the reference: 260 cells;
   - "2. THE PARAMETER VALUES": the baseline regime and the scenarios;
   - "3. THE RUN": `mode` (`"run"`, or `"dry"` to print the cells and the cost
     and run nothing), `cores`, `nrep`, `rep_from`, `rep_to`.
3. Press **Source** (Ctrl+Shift+S).

The replications are handed to the cores one at a time, as each core frees up,
and every one prints a line with its time and an estimate of the time left. The
results go to `results/` beside the script.

The run can be stopped at any moment with the red Stop button and resumed by
pressing Source again: everything is keyed on (cell, replication), so what was
recorded is skipped. Keep the RStudio session open while it runs, and keep the
machine from going to sleep.

To know what a design costs before running it: run the first replication of
every cell (`rep_to = 1`), then `mode = "dry"`. The cost printed is measured on
those replications, for `nrep` replications on `cores` cores.

Do not run the same `tag` on two machines at once: they would write to the same
files. To share the work between machines, give each its own `tag` and its own
range of replications.

If R on that machine uses a multithreaded BLAS (MKL, OpenBLAS), set
`OMP_NUM_THREADS=1` in `~/.Renviron` before starting, otherwise every core
tries to use every core.

## Running from a terminal

Every SETUP setting of sections 1 and 3 can also be given on the command line,
which overrides it; lists are separated by commas:

```
Rscript run-simulations.R --mode=dry                  # the cells and the cost
Rscript run-simulations.R --cores=8                   # the whole study -> results/
Rscript run-simulations.R --blocks=core               # some blocks only
Rscript run-simulations.R --TN=60,120,365 --n=50,100  # one block of other margins
Rscript run-application.R                             # the application -> application/
```

Margins given on the command line replace the blocks with one factorial block
of those margins around the reference.

## The outputs

`results/<tag>.csv`, `<tag>-params.csv`, `<tag>-grid.csv`,
`<tag>-stations.csv` and the folder `<tag>-obs/` share the primary key
(cell, rep) and join on it and on nothing else. What each holds is described at
the top of `run-simulations.R`.
