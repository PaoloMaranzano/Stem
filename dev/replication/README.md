# SC-STEM: replication material

The scripts that produce the simulation study and the application of the paper.
This folder is **not part of the Stem package** and is never shipped with it:
Stem is a dependency like any other, and the scripts are deliberately
disconnected from its sources.

| script | what it does |
|---|---|
| `run-simulations.R` | runs the simulation study: the Monte Carlo, and with `--coverage` the bootstrap coverage experiment |
| `analyse-simulations.R` | turns the results into the tables and figures of the paper, and draws the design |
| `run-application.R` | runs the application, staged and cached |
| `check-theory.R` | verifies numerically the Lemma and the Propositions of the paper, with a verdict for each claim, and writes their tables (`theory/`) |

## Requirements

R, and nothing else to set up. Each script checks that an installed Stem
carries what it needs, and installs or updates it from GitHub when it does not.
The check is on the features, not on the version number, because the
development builds all report 2.0.0.

## Running

Each script runs from whatever folder it sits in and writes beside itself:

```
Rscript run-simulations.R --dry              # the design and its cost, nothing run
Rscript run-simulations.R                    # the whole study   -> results/
Rscript run-simulations.R --coverage         # the coverage      -> results/
Rscript analyse-simulations.R                # tables, figures   -> output/
Rscript run-application.R                    # the application   -> application/
```

`analyse-simulations.R` takes the generator and the design from
`run-simulations.R`, so the two have to sit in the same folder; that is also
what guarantees they cannot disagree.

## The design

It is the SETUP block near the top of `run-simulations.R`: the levels of every
factor, the reference cell, the blocks. Edit the numbers there. Run with
`--dry` to read the number of cells and the cost before committing to it.

The study is not a full factorial. It is a factorial in the three factors that
interact -- number of locations, separation of the regimes, length of the
series -- plus one-factor-at-a-time margins around a reference cell for the
factors that are there to show the results do not turn on them.

## Running in blocks

The whole design is about 200 core-hours at 100 replications, so it is meant to
be cut into slices that compose, on whichever machine is free:

```
Rscript run-simulations.R --blocks=core,null
Rscript run-simulations.R --blocks=scenarios --rep_from=1  --rep_to=50
Rscript run-simulations.R --blocks=scenarios --rep_from=51 --rep_to=100
```

Everything is appended and keyed on (cell, replication), so an interrupted run
keeps what it produced and a resumed one skips it. `--tag=` names the output
files, so several designs sit side by side. The other options are `--nrep`,
`--only_n`, `--only_TN`, `--only_K`, `--only_omega`, `--only_knn`,
`--only_scenario`, `--only_balance`, `--out`, and `--score=conditional` for the
assignment score that conditions on the other members of a regime (run it under
its own `--tag`).

A replication that fails says why, both on the console and in the `error`
column of the results.

## The outputs

`results/<tag>.csv`, `<tag>-params.csv`, `<tag>-stations.csv` and the folder
`<tag>-obs/` share the primary key (cell, rep) and join on it and on nothing
else. What each holds is described at the top of `run-simulations.R`.
