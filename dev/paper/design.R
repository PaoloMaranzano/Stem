## ---------------------------------------------------------------------------
## THE DESIGN OF THE SIMULATION STUDY. This file is the setup: it is the only
## place where the factors and their levels are written down, and everything
## else -- the driver in 07-simulation.R, the table in 13-dgp-table.R -- reads
## them from here. Edit the numbers below and both follow.
##
## The study is NOT a full factorial. Crossing every factor with every other
## would be 2700 cells at K = 3 alone, most of them answering no question. It is
## built instead as a core factorial in the factors that interact -- how many
## locations, how separated the regimes are, how long the series is -- plus
## one-factor-at-a-time margins around a reference cell for the factors that are
## there to show the results do not depend on them.
##
## stem_design_cost() prices whatever is written here from the timings measured
## on this machine, so the budget can be read off before a run is launched:
##
##   Rscript dev/paper/design.R
## ---------------------------------------------------------------------------


## ---------------------------------------------------------------------------
## 1. THE LEVELS. Every factor, every level. Edit freely.
## ---------------------------------------------------------------------------
stem_design_levels <- function() list(

  ## locations. The cost of a cell is cubic in this, so 400 is the ceiling
  n = c(20L, 50L, 100L, 200L, 400L),

  ## periods. 60 is a season, 120 a half-year, 365 a year of daily data
  TN = c(60L, 120L, 365L),

  ## regimes. K = 1 is the null case: one regime, so nothing to recover
  K = c(1L, 3L),

  ## separation of the regime centres, in units of the within-regime dispersion:
  ## the centres sit 2*omega apart and each coordinate has variance nu_sp, so the
  ## standardised separation is 2*omega/sqrt(nu_sp) -- 0, 1.58 and 3.16 sd
  omega = c(0, 0.5, 1),

  ## relative sizes of the regimes
  balance = c("balanced", "unbalanced"),

  ## which parameters differ between regimes; the ten rows of dgp_scenarios()
  scenario = c("S0", "S1a", "S1b", "S1c", "S2a", "S2b", "S3a", "S3b", "S4", "S5"),

  ## The neighbourhood graph the Potts penalty lives on. With point-referenced
  ## data there is no canonical adjacency -- unlike areal data, where a shared
  ## boundary defines it -- so the graph is a modelling choice, and the results
  ## have to be shown not to turn on it.
  knn = c(3L, 5L, 10L)
)


## ---------------------------------------------------------------------------
## 2. THE REFERENCE CELL. The margins below vary one factor at a time around
##    this configuration, so it should be the cell the paper talks about most.
## ---------------------------------------------------------------------------
stem_design_reference <- function() list(
  n        = 100L,
  TN       = 120L,
  K        = 3L,
  omega    = 0.5,
  balance  = "balanced",
  scenario = "S4",
  knn      = 5L
)


## ---------------------------------------------------------------------------
## 3. THE BLOCKS. Each one answers a question and can be switched off on its
##    own with --blocks=core,scenarios on the command line.
##
##    core       recovery as a function of the three things that govern it:
##               how many locations, how separated the regimes, how long the
##               series. Fully crossed, at the reference scenario and graph.
##    null       the same grid at K = 1, where there is no regime to find and
##               the question is whether the procedure invents one.
##    scenarios  what has to differ between regimes for the difference to be
##               found. All ten scenarios, crossed with n and omega, at the
##               reference T.
##    graph      the robustness margin: knn crossed with n and omega, because a
##               graph that is dense relative to the network behaves differently
##               from one that is sparse.
##    balance    the robustness margin for unequal regime sizes.
## ---------------------------------------------------------------------------
stem_design_blocks <- function() c("core", "null", "scenarios", "graph", "balance")


## ---------------------------------------------------------------------------
## 4. Assembling the cells. Nothing below needs editing to change the design.
## ---------------------------------------------------------------------------
stem_design_cells <- function(lv = stem_design_levels(),
                              ref = stem_design_reference(),
                              blocks = stem_design_blocks()) {

  grid <- function(...) expand.grid(..., stringsAsFactors = FALSE,
                                    KEEP.OUT.ATTRS = FALSE)
  out <- list()

  ## core: n x omega x T, at the reference scenario, balance and graph
  if ("core" %in% blocks)
    out$core <- grid(n = lv$n, TN = lv$TN, K = setdiff(lv$K, 1L),
                     omega = lv$omega, balance = ref$balance,
                     id = ref$scenario, knn = ref$knn)

  ## null: K = 1, where omega, balance and the scenario are all vacuous
  if ("null" %in% blocks && 1L %in% lv$K)
    out$null <- grid(n = lv$n, TN = lv$TN, K = 1L, omega = 0,
                     balance = ref$balance, id = ref$scenario, knn = ref$knn)

  ## scenarios: what has to differ, crossed with n and omega at the reference T
  if ("scenarios" %in% blocks)
    out$scen <- grid(n = lv$n, TN = ref$TN, K = setdiff(lv$K, 1L),
                     omega = lv$omega, balance = ref$balance,
                     id = lv$scenario, knn = ref$knn)

  ## graph: knn crossed with n and omega, the two things it interacts with
  if ("graph" %in% blocks)
    out$graph <- grid(n = lv$n, TN = ref$TN, K = setdiff(lv$K, 1L),
                      omega = lv$omega, balance = ref$balance,
                      id = ref$scenario, knn = lv$knn)

  ## balance: unequal regime sizes, crossed with n and omega
  if ("balance" %in% blocks)
    out$bal <- grid(n = lv$n, TN = ref$TN, K = setdiff(lv$K, 1L),
                    omega = lv$omega, balance = setdiff(lv$balance, ref$balance),
                    id = ref$scenario, knn = ref$knn)

  cells <- do.call(rbind, out)
  if (is.null(cells)) return(cells)
  cells <- cells[, c("n", "TN", "K", "omega", "id", "balance", "knn")]
  ## the blocks overlap on the reference configuration: keep one copy
  cells <- unique(cells)
  ## the graph needs strictly fewer neighbours than locations
  cells <- cells[cells$knn < cells$n, ]
  ## cheapest first, so that an interrupted run still covers the design
  cells <- cells[order(cells$n * cells$TN), ]
  rownames(cells) <- NULL
  cells
}


## ---------------------------------------------------------------------------
## 5. What a cell costs.
##
## Seconds for ONE replication -- one (k, phi) grid through SCSTEM_Infocrit()
## and SCSTEM_Select() -- measured on the author's machine, single core, at
## k_grid = 1:4 and phi_grid of length 3, at omega = 0.5. The entries are affine
## in T at fixed n and interpolated in log n. The estimate is meant for
## budgeting, not for reporting.
##
## THE MODEL DOES NOT COVER omega, AND omega MATTERS. At omega = 0 the regime
## centres coincide, so the network is a tight blob: at n = 20 the median
## pairwise distance is 105 km against a true correlation range of 123 km, the
## exponential decay is barely resolved over the observed distances, and theta
## is close to unidentified. The EM then crawls. Measured at n = 20, T = 60,
## K = 3, on a single pooled fit: 0.6 s at omega = 0.5, 0.4 s at omega = 1.0,
## and more than 150 s at omega = 0 without converging. Until that is settled --
## by holding the total spatial variance fixed across omega, so that the
## footprint of the network stops depending on the separation -- the figures
## below understate the cost of every omega = 0 cell by orders of magnitude.
## ---------------------------------------------------------------------------
STEM_COST <- data.frame(
  n    = rep(c(20L, 50L, 100L, 200L, 400L), times = 6L),
  K    = rep(c(1L, 3L), each = 15L),
  TN   = rep(rep(c(60L, 120L, 365L), each = 5L), times = 2L),
  secs = c(  5.4,  7.4, 10.9, 27.5,  99.0,     # K = 1, T = 60
             7.4, 10.3, 15.2, 32.6, 110.1,     # K = 1, T = 120
            15.9, 22.2, 32.6, 53.6, 155.7,     # K = 1, T = 365
             4.1, 12.4, 10.0, 20.6,  68.9,     # K = 3, T = 60
             7.8, 17.2, 15.3, 27.7,  77.5,     # K = 3, T = 120
            23.1, 36.8, 37.2, 56.8, 112.4)     # K = 3, T = 365
)

## seconds for one replication of one cell, interpolated in log(n) and linear
## in T within the measured envelope, held flat outside it
stem_cell_secs <- function(n, TN, K) {
  mapply(function(nn, tt, kk) {
    tab <- STEM_COST[STEM_COST$K == (if (kk == 1L) 1L else 3L), ]
    Ts  <- sort(unique(tab$TN))
    ## interpolate in log n at each measured T, then in T
    at_T <- vapply(Ts, function(t0) {
      s <- tab[tab$TN == t0, ]
      s <- s[order(s$n), ]
      stats::approx(log(s$n), s$secs, xout = log(nn), rule = 2)$y
    }, numeric(1))
    stats::approx(Ts, at_T, xout = tt, rule = 2)$y
  }, n, TN, K)
}

stem_design_cost <- function(cells, nrep = 100L, cores = 12L) {
  secs <- stem_cell_secs(cells$n, cells$TN, cells$K)
  tot  <- sum(secs) * nrep
  list(cells = nrow(cells), nrep = nrep, cores = cores,
       core_hours = tot / 3600, wall_hours = tot / 3600 / cores,
       one_rep_hours = sum(secs) / 3600)
}


## ---------------------------------------------------------------------------
## 6. The marginal levels, for the scripts that describe the design
## ---------------------------------------------------------------------------
stem_design_dims <- function(lv = stem_design_levels())
  lv[c("TN", "n", "K", "omega", "balance")]


## ---------------------------------------------------------------------------
## Run this file directly to price the design that is written above.
## ---------------------------------------------------------------------------
if (sys.nframe() == 0L) {
  cl <- stem_design_cells()
  bl <- stem_design_blocks()
  cat("\nblocks:", paste(bl, collapse = ", "), "\n\n")
  for (b in bl) {
    one <- stem_design_cells(blocks = b)
    cat(sprintf("  %-10s %5d cells, %6.1f core-h at M = 100\n", b, nrow(one),
                stem_design_cost(one)$core_hours))
  }
  co <- stem_design_cost(cl)
  cat(sprintf("\n  %-10s %5d cells after removing the overlap between blocks\n",
              "TOTAL", co$cells))
  cat(sprintf("\n  one replication of the whole design : %6.2f h\n", co$one_rep_hours))
  cat(sprintf("  M = 100                            : %6.0f core-h\n", co$core_hours))
  cat(sprintf("  on %d cores                        : %6.0f h wall clock\n\n",
              co$cores, co$wall_hours))
  cat("cells by number of locations\n")
  print(table(n = cl$n, T = cl$TN))
  cat("\n")
}
