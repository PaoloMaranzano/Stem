# Archived: the generic simulation design in the paper

On 2026-09-27 the simulation study was removed from the SC-STEM paper, to be
designed again from scratch. The text as it stood just before the removal is
kept here in full:

| file | content |
|---|---|
| `main-before-sim-removal.tex` | the manuscript, with Section 3 "Simulation experiments": the overlap parameter, the generator, the baseline table, the coupling of the latent processes, the scenarios S0 to S5, the neighbourhood graph as a factor, the blocks and the metrics |
| `supplement-before-sim-removal.tex` | the supplement, with the data-generating process as an algorithm, the constants, the regime sizes, the fitted grid, and the placeholders for the per-scenario results and the bootstrap coverage |
| `Figures/` | the tables and figures Section 3 input: `tab_overlap`, `tab_dgp`, `fig_design_omega`, `fig_dgp_examples`, `fig_design_knn`, `fig_design_graph` |

The runner that implements this design is `dev/replication/run-simulations.R`
at commit f2a0ad8. The design before it, calibrated on the Po Valley network,
is in `dev/archive/sim-design-povalley`.
