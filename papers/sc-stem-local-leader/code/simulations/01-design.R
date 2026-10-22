## ===========================================================================
## SIMULATIONS, SCRIPT 1 OF 3: THE DESIGN
##
## The general setting of the simulation study (Section 4.1 of the paper and
## Supplement B): the scenario-variants, the reference levels of the margins,
## the blocks of the design and the parameter values of the regimes. The
## generator that uses them is in functions.R.
##
## This file is loaded by 02-run.R and 03-results.R. Run on its own,
##
##     Rscript 01-design.R
##
## it prints the scenario-variants, the parameters of their regimes and the
## number of cells of every block (332 cells in all). The figures and tables
## that describe the design in the paper are drawn by 03-results.R.
## ===========================================================================

## Run on its own: the functions first (they need Stem, see functions.R)
SIM_DESIGN_MAIN <- !exists("sim_scenarios", mode = "function")
if (SIM_DESIGN_MAIN) {
  local({
    a <- commandArgs(trailingOnly = FALSE); m <- grep("^--file=", a, value = TRUE)
    f <- if (length(m)) sub("^--file=", "", m[1]) else NULL
    if (is.null(f)) for (i in rev(seq_len(sys.nframe()))) {
      of <- sys.frame(i)$ofile
      if (!is.null(of)) { f <- of; break }
    }
    here <- if (is.null(f)) getwd() else dirname(normalizePath(f, winslash = "/"))
    if (!exists("SIM_DEFINE_ONLY", envir = globalenv())) assign("SIM_DEFINE_ONLY", FALSE, envir = globalenv())
    source(file.path(here, "functions.R"), local = globalenv())
  })
}

## ---------------------------------------------------------------------------
## 1. THE DESIGN: the scenario-variants, the reference levels of the margins,
##    and the blocks. Each block is the full factorial of the margins it lists,
##    every other margin at its reference level, crossed with the
##    scenario-variants; the design is the union of the blocks, and a cell that
##    two blocks share is run once.
## ---------------------------------------------------------------------------

## The 19 scenario-variants (simulation-design.tex, Section 6.4). The name is
## <scenario>-<variant>:
##   S1w, S1s        complete heterogeneity, weak and strong separation
##   S2              complete homogeneity: the pooled STEM model
##   S3beta, S3G, S3Seta, S3theta, S3error
##                   one block common to the regimes, the others as S1s
##   ind  rho = 0, fields by regime: exactly the SC-STEM model
##   shr  rho = 1, one field where theta is common, otherwise by regime
##   lat  rho = 1, fields by regime  (only where theta is common)
##   fld  rho = 0, one field         (only where theta is common)
## The nine of the paper. The ten additional ones, to be added back when
## wanted, are "S1w-ind", "S3beta-ind", "S3G-ind", "S3Seta-ind",
## "S3theta-ind", "S3theta-lat", "S3theta-fld", "S3error-ind",
## "S3error-lat", "S3error-fld".
SIM_SCENARIOS <- c("S2", "S1s-ind", "S1s-shr", "S1w-shr",
                   "S3beta-shr", "S3G-shr", "S3Seta-shr", "S3theta-shr", "S3error-shr")

## The reference levels of the margins:
##   n        number of locations; not below 40, since the grid reaches K = 4
##            regimes of at least N_MIN = 6 locations each
##   TN       length of the series
##   omega    spatial overlap of the three regimes: the centres sit at the
##            vertices of an equilateral triangle of side 2*omega, with the
##            total variance of a coordinate held fixed. 0: the regimes
##            coincide in space; 0.70: contiguous areas interpenetrating along
##            their borders (1.63 within-regime standard deviations between
##            centres); 1: well apart. Vacuous in S2
##   spread   what is held fixed as omega moves the centres apart. "total": the
##            variance of a coordinate over the whole network, so the regimes
##            shrink as they separate (omega = 1: 3.2 standard deviations
##            between centres, 10% of the locations nearer another centre);
##            "regime": the variance within a regime, at its value at omega =
##            0.70, so the network widens instead (omega = 1: 2.3 standard
##            deviations, 20% nearer another centre). The two coincide at
##            omega = 0.70. Vacuous in S2
##   balance  relative sizes of the regimes, "balanced" or "unbalanced"
##            (1:2:3). Vacuous in S2
##   knn      neighbours of the graph of the Potts penalty
SIM_REFERENCE <- list(n = 100L, TN = 120L, omega = 0.70, spread = "total", balance = "balanced",
                      knn = 5L)

## The blocks (simulation-design.tex, Section 7). The core crosses the three
## factors that govern recovery; the others vary one factor at a time around
## the reference.
SIM_BLOCKS <- list(
  core    = list(n = c(40L, 100L, 200L), TN = c(60L, 120L, 365L), omega = c(0, 0.70, 1)),
  ## the fourth geometry of the core: regimes 2 omega apart with the spread
  ## they have at omega = 0.70 (medium-low overlap)
  spread  = list(n = c(40L, 100L, 200L), TN = c(60L, 120L, 365L), omega = 1, spread = "regime"),
  balance = list(balance = "unbalanced"),
  knn     = list(knn = c(3L, 10L)),
  n400    = list(n = 400L)
)

## ---------------------------------------------------------------------------
## 2. THE PARAMETER VALUES (simulation-design.tex, Sections 5 and 6).
##
## Each regime is described by
##   b0, b1  intercept and effect of the covariate
##   G       persistence of the latent process
##   v       stationary variance of the latent process; sigma2_eta = v (1 - G^2)
##   se, so  nugget and partial sill of the measurement error
##   R       practical range of the error field; theta = 3 / R
## ---------------------------------------------------------------------------

## The baseline: regime 1 of every scenario and the only regime of S2. The
## variance of z is then 0.55^2 + 0.30 + 0.20 + 0.20 = 1, so every contrast
## reads in standard deviations of the response.
SIM_BASE <- list(b0 = 2, b1 = 0.55, G = 0.80, v = 0.30, se = 0.20, so = 0.20, R = 2)

## S1, complete heterogeneity: regimes 1, 2, 3 at the two levels of separation
SIM_S1 <- list(
  weak = list(b0 = c(2.00, 2.10, 2.20), b1 = c(0.55, 0.65, 0.75),
              G  = c(0.80, 0.70, 0.60), v  = c(0.30, 0.35, 0.40),
              se = c(0.20, 0.25, 0.30), so = c(0.20, 0.175, 0.15),
              R  = c(2, 1.5, 1)),
  ## regime 3 has nugget 0.40, sill 0.20 and range 1 (it had 0.50, 0.10, 0.5):
  ## with a sill one sixth of the error variance and a range of three
  ## nearest-neighbour spacings its spatial parameters were not identified, and
  ## a long-range field absorbed part of its almost white latent process. The
  ## total error variance and the range, all the assignment score sees, are
  ## unchanged or invisible to it, so the separation of the scenarios is the same.
  ## Regime 2 has range 1.5 (it had 1), so that the ranges still decrease.
  strong = list(b0 = c(2.00, 2.30, 2.60), b1 = c(0.55, 0.85, 1.15),
                G  = c(0.80, 0.50, 0.20), v  = c(0.30, 0.40, 0.50),
                se = c(0.20, 0.35, 0.40), so = c(0.20, 0.15, 0.20),
                R  = c(2, 1.5, 1))
)

## S3 keeps the values of S1 at this level for the blocks that differ
SIM_S3_LEVEL <- "strong"


## The table of the scenario-variants, from the definitions above
SIM_SCEN <- sim_scenarios()

## Run on its own: print the design
if (SIM_DESIGN_MAIN) {
  if (!exists("SIM_CLI_MARGINS")) SIM_CLI_MARGINS <- character(0)
  cat("THE SCENARIO-VARIANTS\n"); print(SIM_SCEN[SIM_SCEN$id %in% SIM_SCENARIOS, ], row.names = FALSE)
  for (id in SIM_SCENARIOS) {
    cat("\n", id, ": the parameters of the regimes\n", sep = "")
    print(dgp_psi(SIM_SCEN[SIM_SCEN$id == id, ]), row.names = FALSE, digits = 3)
  }
  cells <- sim_cells(list(scenario = SIM_SCENARIOS, blocks = names(SIM_BLOCKS)), SIM_BLOCKS)
  cat(sprintf("\nTHE CELLS: %d\n", nrow(cells))); print(table(block = cells$block))
}
