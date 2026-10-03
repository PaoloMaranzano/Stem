## ===========================================================================
## SC-STEM: the simulation study. THE RUNNER.
##
## One script, self-contained. It needs nothing but R, runs from whatever folder
## it sits in -- the Google Drive folder included -- and writes its results
## beside itself, in results/.
##
## The design is described in full in simulation-design.tex, in the Overleaf
## project of the paper: the generator, the two variants of dependence between
## regimes, the baseline regime, the scenarios and their values, and the
## separation of the regimes. This file implements it and nothing else.
##
## HOW TO RUN IT
##
##   1. Open this file in RStudio.
##   2. If you want, change the SETUP below:
##        - "1. THE DESIGN": the scenario-variants, the reference level of every
##          margin, and the blocks. Each block is the full factorial of the
##          margins it lists, the others at the reference; the design is the
##          union of the blocks;
##        - "2. THE PARAMETER VALUES": the baseline regime and the scenarios;
##        - "3. THE RUN": mode ("run" or "dry"), cores, replications.
##   3. Press Source (Ctrl+Shift+S).
##
## The first time, Stem is installed or updated from GitHub. The replications
## are then spread over the cores, and every one prints a line with its time and
## the time left. The run can be stopped at any moment with the red Stop button
## and resumed by pressing Source again: what was recorded is skipped. Keep the
## RStudio session open while it runs.
##
## mode = "dry" prints the cells and the cost of the whole design, run nothing.
## Once some replications are recorded the cost is measured on them: run the
## first replication of every cell (rep_to = 1), then "dry", and the estimate
## printed is the one to plan on.
##
## SEVERAL MACHINES. Give each machine its own `stream` (1, 2, ...) and as many
## replications as it can carry, for instance
##
##     machine 1:  stream = 1, nrep = 10, cores = 3
##     machine 2:  stream = 2, nrep = 40, cores = 7
##
## Each stream draws its own data, with seeds that no other stream uses, and
## writes its own files, <tag>_s<stream>*. Copy the results/ files of one
## machine beside those of the other: the analysis stacks every stream of the
## tag, here 50 replications per cell. Replication r of stream s is recorded
## as rep = 1000 (s - 1) + r, its seed is 1000 rep + 1. Never run the same
## stream on two machines at once: they would draw the same data.
##
## THE TIME LIMIT. Every replication runs in a separate R process (package
## callr) under a time limit set from the times already recorded: cap_mult times
## the median of the cell once cap_after replications of it are done, before
## that cap_first_mult times the median of the cells with the same (n, T), and
## cap_first_min seconds when nothing is known yet. A replication over its limit
## is stopped, written to <tag>_s<stream>-timeouts.csv, and drawn again with
## the seed of the next attempt (1000 rep + 1 + 100 (attempt - 1)), up to
## cap_attempts; the last one is recorded as failed. cap_mult = 0 removes the
## limit.
##
## Every setting of the SETUP can also be given from a terminal, which overrides
## it; lists are separated by commas:
##
##     Rscript run-simulations.R --stream=2 --nrep=40 --cores=7
##     Rscript run-simulations.R --cores=8 --rep_from=11 --rep_to=20
##     Rscript run-simulations.R --blocks=core --mode=dry
##
## Margins given this way (--n, --TN, --omega, --balance, --knn) replace the
## blocks with one factorial block of those margins around the reference:
##
##     Rscript run-simulations.R --TN=60,120,365 --mode=dry
##
## This is part of the REPLICATION MATERIAL of the paper, not of the Stem
## package, and nothing here is shipped with it.
##
## SIX OUTPUTS per stream, keyed on (cell, rep) and joining on it and on
## nothing else; rep is unique across the streams, so the outputs of all the
## streams stack. Below, <tag> stands for <tag>_s<stream>:
##
##   <tag>.csv           one row per replication: the selected (K, phi), the
##                       recovery of the partition, the error of the fitted
##                       signal, the times
##   <tag>-params.csv    one row per replication, regime and parameter: the
##                       truth beside the estimate at the true number of
##                       regimes, after aligning the estimated regimes on the
##                       true ones
##   <tag>-stations.csv  one row per replication and location
##   <tag>-grid.csv      one row per replication and (K, phi) of the fitted
##                       grid: log-likelihood, parameters, criteria, and the
##                       ARI of its partition with the truth
##   <tag>-obs/*.rds     the full observations of the first keep_obs
##                       replications of every cell
##   <tag>-timeouts.csv  one row per replication stopped by the time limit
## ===========================================================================


## ---------------------------------------------------------------------------
## Where this script is, and therefore where its results go
## ---------------------------------------------------------------------------
sim_this_file <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[1]), winslash = "/"))
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(normalizePath(of, winslash = "/"))
  }
  NULL
}
SIM_FILE <- sim_this_file()
SIM_HERE <- if (is.null(SIM_FILE)) normalizePath(getwd(), winslash = "/") else dirname(SIM_FILE)

## The cores to use by default: all the physical ones but one, so that the
## machine stays usable while the study runs.
sim_default_cores <- function() {
  n <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (is.na(n)) n <- suppressWarnings(parallel::detectCores())
  if (is.na(n)) 1L else max(1L, as.integer(n) - 1L)
}

## Definitions only? The worker processes source this file with
## SIM_DEFINE_ONLY <- TRUE: they get the generator and the functions, and
## nothing is run, installed or read from the command line.
if (!exists("SIM_DEFINE_ONLY", inherits = FALSE)) SIM_DEFINE_ONLY <- FALSE


## ---------------------------------------------------------------------------
## The only prerequisite: Stem
##
## Stem is installed from GitHub (SIM_STEM_REF), and again whenever GitHub holds
## a newer commit than the installed one, so the study always runs on the
## current code. Offline, an installed Stem that carries what the study uses is
## accepted as it is. For the paper SIM_STEM_REF should be pinned to the commit
## the study was run with (e.g. "PaoloMaranzano/Stem@0bdf027").
## ---------------------------------------------------------------------------
sim_require <- function(pkgs) {
  for (p in pkgs) {
    if (requireNamespace(p, quietly = TRUE)) next
    message("installing ", p, " ...")
    utils::install.packages(p, repos = "https://cloud.r-project.org")
    if (!requireNamespace(p, quietly = TRUE))
      stop("could not install '", p, "'", call. = FALSE)
  }
}

## Pinned to the commit of 2026-10-03: the commit of 2026-10-02 (STEM_control,
## the stopping rule with em_stop = "any", the regimes started from their own
## least squares and then warm, the final refits shared across the grid, the
## names K and A) plus the four estimation algorithms with SQUAREM as the
## default, tolerances 1e-3 on both criteria, regularization 0, and the
## penalized fits started from the refits of the unpenalized solution
SIM_STEM_REF <- "PaoloMaranzano/Stem@6808cfeb5a1b3fff9258e38945b93d8e1b495e51"

## The commit GitHub holds for SIM_STEM_REF, or NA when it cannot be reached.
sim_github_sha <- function(ref = SIM_STEM_REF) {
  repo <- sub("@.*$", "", ref)
  at   <- if (grepl("@", ref, fixed = TRUE)) sub("^.*@", "", ref) else "HEAD"
  j <- tryCatch({
    con <- url(sprintf("https://api.github.com/repos/%s/commits/%s", repo, at))
    on.exit(close(con))
    paste(readLines(con, warn = FALSE), collapse = "")
  }, error = function(e) "", warning = function(w) "")
  m <- regmatches(j, regexpr("\"sha\" *: *\"[0-9a-f]{40}\"", j))
  if (length(m)) gsub("[^0-9a-f]", "", sub("^\"sha\" *: *", "", m)) else NA_character_
}

## The installed Stem is GitHub's commit when GitHub can be reached, and
## otherwise carries what the study uses. The version number cannot tell: the
## development builds all say 2.0.0. The commit is read from the DESCRIPTION
## file without loading the package: a namespace loaded before the update would
## read the new files with the old index.
sim_stem_ok <- function(latest = NA_character_) {
  if (!nzchar(system.file(package = "Stem"))) return(FALSE)
  if (!is.na(latest)) {
    here <- utils::packageDescription("Stem")$RemoteSha
    return(!is.null(here) && identical(here, latest))
  }
  if (!requireNamespace("Stem", quietly = TRUE)) return(FALSE)
  ns <- asNamespace("Stem")
  have <- c("STEM_control", "STEM_Signal", "SCSTEM_Signal", "SCSTEM_Infocrit", "SCSTEM_Select",
            "scstem_neighbors", "scstem_align_labels", "scstem_ari")
  all(vapply(have, exists, logical(1), envir = ns, inherits = FALSE)) &&
    "distance" %in% names(formals(get("scstem_neighbors", envir = ns))) &&
    "K_grid" %in% names(formals(get("SCSTEM_Infocrit", envir = ns))) &&  # the names of 2.0.0
    "algorithm" %in% names(formals(get("STEM_control", envir = ns)))     # the algorithms
}

SIM_STEM_SHA <- if (SIM_DEFINE_ONLY) NA_character_ else sim_github_sha()
if (!sim_stem_ok(SIM_STEM_SHA)) {
  message(if (nzchar(system.file(package = "Stem")))
            "the installed Stem is not the pinned commit: installing it" else
            "Stem is not installed: installing it from GitHub")
  sim_require("remotes")
  if ("Stem" %in% loadedNamespaces()) try(unloadNamespace("Stem"), silent = TRUE)
  remotes::install_github(SIM_STEM_REF, upgrade = "never", force = TRUE, quiet = TRUE)
  if (!sim_stem_ok(SIM_STEM_SHA))
    stop("Stem could not be installed or brought up to date. Restart R ",
         "(Session > Restart R in RStudio) and run again; if it still fails, ",
         "install it by hand with\n  remotes::install_github(\"", SIM_STEM_REF,
         "\")", call. = FALSE)
}
sim_require("spdep")
if (!SIM_DEFINE_ONLY) sim_require("callr")   # runs a replication under its time limit
suppressPackageStartupMessages(library("Stem"))


## ---------------------------------------------------------------------------
## Command-line overrides: --name=value, a value split on commas and coerced to
## the type of the default. --name= with nothing after the sign keeps the
## default.
## ---------------------------------------------------------------------------
sim_config <- function(defaults, args = commandArgs(trailingOnly = TRUE)) {
  kv <- grep("^--[^=]+=", args, value = TRUE)
  for (a in kv) {
    nm  <- sub("^--([^=]+)=.*$", "\\1", a)
    val <- sub("^--[^=]+=", "", a)
    if (!nzchar(val)) next
    if (!nm %in% names(defaults))
      stop("unknown option --", nm, "; the options are: ",
           paste(names(defaults), collapse = ", "), call. = FALSE)
    parts <- trimws(strsplit(val, ",", fixed = TRUE)[[1]])
    d <- defaults[[nm]]
    defaults[[nm]] <-
      if (is.logical(d)) as.logical(parts)
      else if (is.integer(d)) as.integer(parts)
      else if (is.numeric(d))
        vapply(parts, function(s) eval(parse(text = s)), numeric(1),
               USE.NAMES = FALSE)
      else parts
    if (anyNA(defaults[[nm]])) stop("cannot read --", nm, "=", val, call. = FALSE)
  }
  defaults
}

sim_print_config <- function(cfg) {
  cat("configuration\n")
  for (nm in names(cfg))
    cat(sprintf("  %-12s %s\n", nm, paste(format(cfg[[nm]]), collapse = ", ")))
  cat("\n")
}


## ===========================================================================
## THE SETUP. This is the part meant to be edited.
## ===========================================================================

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

## ---------------------------------------------------------------------------
## 3. THE RUN. What this invocation does.
## ---------------------------------------------------------------------------
CFG <- sim_config(args = if (SIM_DEFINE_ONLY) character(0) else commandArgs(TRUE), c(
  list(scenario = SIM_SCENARIOS),
  SIM_REFERENCE,
  list(
    ## the blocks this invocation runs, by name
    blocks   = names(SIM_BLOCKS),

    ## "run": the Monte Carlo; "dry": print the cells and the cost, run nothing
    mode     = "run",

    ## How many cores. Each core runs one replication at a time.
    cores    = sim_default_cores(),

    ## The stream: one per machine. Two machines that run the study at the
    ## same time give themselves different streams; each draws its own data
    ## (different seeds) and writes its own files, <tag>_s<stream>*, and the
    ## analysis stacks every stream of a tag. Copy the results/ files of one
    ## machine beside those of the other to put them together.
    stream   = 1L,

    ## replications per cell in this stream, and the slice this invocation
    ## covers (rep_to defaults to nrep). A stream can also be run in slices:
    ## the files are appended, so the slices compose. At most 999 per stream.
    nrep     = 10L,
    rep_from = 1L,
    rep_to   = NA_integer_,

    ## What the estimator searches over. phi_ref is the penalty at which the
    ## recovery at the true number of regimes is read, a point of phi_grid.
    K_grid   = 1:4,
    phi_grid = c(0, 0.025, 0.05, 0.1, 0.2, 0.5, 1),
    phi_ref  = 0.05,

    ## How many replications per cell keep their full per-observation record.
    ## One: every record is n x T rows, and five per cell would take gigabytes;
    ## any replication can be regenerated from its seed.
    keep_obs = 1L,

    ## The time limit of a replication, adaptive cell by cell. Once a cell has
    ## cap_after replications that ended normally, the limit is cap_mult times
    ## their median (at least cap_min seconds); before that, cap_first_mult
    ## times the median of the replications with the same (n, T) in any cell
    ## (at least cap_first_min seconds), or cap_first_min when there are none.
    ## A replication over its limit is stopped, written to <tag>-timeouts.csv
    ## and drawn again with the seed of its next attempt, up to cap_attempts
    ## attempts; after the last it is recorded as failed. Each replication then
    ## runs in a process of its own (package callr), which is what lets it be
    ## stopped. cap_mult = 0 switches the limit off.
    cap_mult       = 10,
    cap_min        = 600,
    cap_after      = 3L,
    cap_first_mult = 20,
    cap_first_min  = 3600,
    cap_attempts   = 3L,

    ## Bookkeeping. `out` defaults to a results/ folder beside this script.
    ## "main3": the design of October 2026 (new values of the strong level,
    ## the fourth geometry) with the estimator of the commit pinned above
    ## (SQUAREM, tolerances 1e-3, no regularization). "main2", the same design
    ## with the EM algorithm at regularization 0.01, and "main", the first pass
    ## of September, are not comparable with it.
    tag      = "main3",
    out      = file.path(SIM_HERE, "results")
  )))

## Margins given on the command line replace the blocks with one factorial
## block of those margins around the reference, as in
##   Rscript run-simulations.R --TN=60,365 --n=50
SIM_CLI_MARGINS <- if (SIM_DEFINE_ONLY) character(0) else
  intersect(names(SIM_REFERENCE),
            sub("^--([^=]+)=.*$", "\\1", grep("^--[^=]+=.", commandArgs(TRUE), value = TRUE)))


## ===========================================================================
## Helpers
## ===========================================================================

## Two internal routines of Stem score a partition against the truth: the
## relocation of labels by the majority rule and the Adjusted Rand Index.
scstem_align_labels <- utils::getFromNamespace("scstem_align_labels", "Stem")
scstem_ari          <- utils::getFromNamespace("scstem_ari", "Stem")

## The model object, with starting values from the pooled OLS fit, as a user
## would build it.
sim_model <- function(dat) {
  n   <- ncol(dat$z)
  ols <- stats::lm.fit(x = dat$covariates, y = as.vector(dat$z))
  s2  <- stats::var(ols$residuals)
  Stem::STEM_Model(
    z = dat$z, covariates = dat$covariates, coordinates = dat$coordinates,
    phi = list(beta = matrix(ols$coefficients, ncol = 1),
               sigma2eps = 0.7 * s2, sigma2omega = 0.3 * s2,
               theta = 1 / stats::median(stats::dist(dat$coordinates)),
               G = matrix(0.8, 1, 1),
               Sigmaeta = matrix(0.2 * s2, 1, 1),
               m0 = as.matrix(0), C0 = as.matrix(1)),
    A = matrix(1, n, 1))
}

## Write, and if the file is momentarily locked -- a synchronization client
## holds a file while it uploads it -- wait and try again.
sim_retry <- function(write, tries = 20L, wait = 3) {
  for (i in seq_len(tries)) {
    ok <- tryCatch({ write(); TRUE }, error = function(e) {
      if (i == tries) stop(e)
      FALSE
    })
    if (ok) return(invisible(NULL))
    Sys.sleep(wait)
  }
}

## Append to a CSV, writing the header only when the file is new.
sim_append <- function(df, path) {
  if (is.null(df) || !nrow(df)) return(invisible(NULL))
  sim_retry(function()
    utils::write.table(df, path, sep = ",", row.names = FALSE,
                       col.names = !file.exists(path), append = file.exists(path)))
}

sim_hours <- function(x) if (x < 1) sprintf("%.0f min", 60 * x) else sprintf("%.1f h", x)


## ===========================================================================
## THE GENERATOR (simulation-design.tex, Section 3). Machinery: nothing below
## needs editing to change the design.
## ===========================================================================

## ---------------------------------------------------------------------------
## Locations: the overlap design
##
##   s_i | g_i = g ~ N_2(mu_g, nu_sp I_2),   nu_sp = NU_TOT - Var(mu),
##
## with the three centres at the vertices of an equilateral triangle of side
## 2*omega. The total variance of a coordinate is held at NU_TOT, so the network
## covers the same area, about four units across, whatever omega is. In S2
## there are no regimes and the locations form one cloud of variance NU_TOT.
## ---------------------------------------------------------------------------
NU_TOT <- 0.4 + 2/3
N_MIN  <- 6L        # smallest regime the design allows
IMB_MIN <- 1.5      # smallest largest/smallest ratio that counts as unbalanced

dgp_centres <- function(K, omega) {
  if (K == 1L) return(cbind(0, 0))
  if (K == 3L) {
    r <- 2 * omega / sqrt(3); a <- c(90, 210, 330) * pi / 180
    return(cbind(r * cos(a), r * sin(a)))
  }
  stop("the design has K = 1 or K = 3 regimes")
}

dgp_centre_var <- function(K, omega) {
  if (K == 1L) return(0)
  mu <- dgp_centres(K, omega)
  mean(apply(mu, 2, function(v) mean((v - mean(v))^2)))
}

## regime sizes: equal up to the remainder, or proportional to 1:2:3 with no
## regime below N_MIN
dgp_sizes <- function(n, K, balance = "balanced", n_min = N_MIN) {
  if (K == 1L) return(n)
  if (balance == "balanced") {
    s <- rep(n %/% K, K)
    if (n %% K) s[seq_len(n %% K)] <- s[seq_len(n %% K)] + 1L
    return(as.integer(s))
  }
  w <- seq_len(K) / sum(seq_len(K))
  s <- pmax(1L, as.integer(round(n * w)))
  s[K] <- n - sum(s[-K])
  if (any(s < n_min)) { s <- pmax(s, n_min); s[K] <- n - sum(s[-K]) }
  as.integer(s)
}

dgp_feasible <- function(n, K, balance = "balanced") {
  s <- dgp_sizes(n, K, balance)
  ok <- all(s >= N_MIN) && sum(s) == n
  if (ok && K > 1L && balance == "unbalanced") ok <- max(s) / min(s) >= IMB_MIN
  ok
}

## spread = "total": the within-regime variance is what NU_TOT leaves to the
## centres; spread = "regime": it is fixed at its value at the reference omega.
dgp_locations <- function(n, K, omega, balance = "balanced", seed = 1, spread = "total") {
  set.seed(seed)
  g  <- rep(seq_len(K), times = dgp_sizes(n, K, balance))
  mu <- dgp_centres(K, omega)
  nu_sp <- NU_TOT - dgp_centre_var(K, if (identical(spread, "regime")) SIM_REFERENCE$omega else omega)
  if (nu_sp <= 0) stop("omega = ", omega, " spreads the centres beyond NU_TOT", call. = FALSE)
  xy <- cbind(sx = mu[g, 1] + stats::rnorm(n, sd = sqrt(nu_sp)),
              sy = mu[g, 2] + stats::rnorm(n, sd = sqrt(nu_sp)))
  list(coords = xy, labels = g)
}

## ---------------------------------------------------------------------------
## Covariate: an AR(1) in time with spatially correlated innovations,
##
##   x_t = a x_{t-1} + sqrt(1 - a^2) u_t,   u_t ~ N_n(0, C_x),
##   (C_x)_ij = exp(-3 h_ij / R_x),
##
## a = 0.7 and R_x = 4, then standardized over all n*T values.
## ---------------------------------------------------------------------------
dgp_covariate <- function(coords, TN, a_time = 0.7, range = 4, seed = 1) {
  d  <- nrow(coords)
  dm <- as.matrix(stats::dist(coords))
  L  <- chol(exp(-3 * dm / range) + diag(1e-8, d))
  set.seed(seed)
  x <- matrix(NA_real_, TN, d)
  x[1, ] <- as.numeric(crossprod(L, stats::rnorm(d)))
  for (tt in 2:TN) {
    x[tt, ] <- a_time * x[tt - 1, ] +
      sqrt(1 - a_time^2) * as.numeric(crossprod(L, stats::rnorm(d)))
  }
  (x - mean(x)) / stats::sd(x)
}

## ---------------------------------------------------------------------------
## The scenario-variants (simulation-design.tex, Section 6)
##
## Each row: the scenario, its level, the common block, the variant, rho and
## the construction of the error field. Built by rule, so that it cannot drift
## from the definition: a variant "shr" shares the field only where theta is
## common, and the mixed variants exist only there.
## ---------------------------------------------------------------------------
sim_scenarios <- function() {
  base <- data.frame(
    scen   = c("S1w", "S1s", "S3beta", "S3G", "S3Seta", "S3theta", "S3error"),
    family = c("S1", "S1", "S3", "S3", "S3", "S3", "S3"),
    level  = c("weak", "strong", rep(SIM_S3_LEVEL, 5)),
    common = c("", "", "beta", "G", "Seta", "theta", "error"),
    stringsAsFactors = FALSE)
  rows <- list()
  for (i in seq_len(nrow(base))) {
    b <- base[i, ]
    theta_common <- b$common %in% c("theta", "error")
    v <- rbind(
      data.frame(variant = "ind", rho = 0, field = "regime"),
      data.frame(variant = "shr", rho = 1, field = if (theta_common) "one" else "regime"))
    if (theta_common) v <- rbind(v,
      data.frame(variant = "lat", rho = 1, field = "regime"),
      data.frame(variant = "fld", rho = 0, field = "one"))
    rows[[i]] <- cbind(id = paste0(b$scen, "-", v$variant), b[rep(1, nrow(v)), ], v)
  }
  rows[[length(rows) + 1L]] <- data.frame(
    id = "S2", scen = "S2", family = "S2", level = "", common = "all",
    variant = "pooled", rho = 1, field = "one")
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
SIM_SCEN <- sim_scenarios()

## ---------------------------------------------------------------------------
## The parameters of the regimes of a scenario, one row per regime. In S3 the
## common block takes the baseline value; the dynamics are parametrized by
## (G, v), so a common G keeps the amplitudes v of S1, and a common sigma2_eta
## keeps the persistences G of S1 and recomputes v = sigma2_eta / (1 - G^2).
## ---------------------------------------------------------------------------
dgp_psi <- function(row) {
  if (row$family == "S2") {
    p <- as.data.frame(SIM_BASE)
  } else {
    p <- as.data.frame(SIM_S1[[row$level]])
    b <- SIM_BASE
    if (nzchar(row$common)) switch(row$common,
      beta  = { p$b0 <- b$b0; p$b1 <- b$b1 },
      G     = { p$G <- b$G },
      Seta  = { p$v <- b$v * (1 - b$G^2) / (1 - p$G^2) },
      theta = { p$R <- b$R },
      error = { p$se <- b$se; p$so <- b$so; p$R <- b$R },
      stop("unknown common block '", row$common, "'", call. = FALSE))
  }
  p$s2eta <- p$v * (1 - p$G^2)
  p$theta <- 3 / p$R
  p
}

## ---------------------------------------------------------------------------
## Latent processes (Section 3.3):
##
##   y^(g)_t = G_g y^(g)_{t-1} + eta^(g)_t,
##   eta^(g)_t = sigma_eta,g (sqrt(rho) xi_t + sqrt(1 - rho) zeta_g,t),
##
## xi_t the shock common to all regimes, zeta_g,t the specific one; started at
## zero and run through a burn-in that is discarded.
## ---------------------------------------------------------------------------
dgp_latent <- function(TN, psi, rho, burn = 200L) {
  K  <- nrow(psi)
  se <- sqrt(psi$s2eta)
  y  <- matrix(NA_real_, TN, K)
  cur <- rep(0, K)
  for (tt in seq_len(burn + TN)) {
    eta <- se * (sqrt(rho) * stats::rnorm(1) + sqrt(1 - rho) * stats::rnorm(K))
    cur <- psi$G * cur + eta
    if (tt > burn) y[tt - burn, ] <- cur
  }
  y
}

## ---------------------------------------------------------------------------
## Measurement error (Section 3.4), independent over t:
##
##   "regime"  e_{t,I_g} ~ N(0, se_g I + so_g exp(-theta_g h)), independent
##             across regimes
##   "one"     e_ti = sqrt(se_g) eps_ti + sqrt(so_g) w_t(s_i), with one field
##             w_t ~ N(0, exp(-theta h)) over the whole network; theta common
## ---------------------------------------------------------------------------
dgp_error <- function(labels, psi, coords, TN, field) {
  d  <- length(labels)
  dm <- as.matrix(stats::dist(coords))
  e  <- matrix(NA_real_, TN, d)
  if (field == "one") {
    th <- unique(round(psi$theta, 12))
    if (length(th) != 1L) stop("one error field needs a common theta", call. = FALSE)
    L <- chol(exp(-th * dm) + diag(1e-10, d))
    W <- matrix(stats::rnorm(TN * d), TN, d) %*% L
    E <- matrix(stats::rnorm(TN * d), TN, d)
    e <- sweep(E, 2, sqrt(psi$se[labels]), "*") + sweep(W, 2, sqrt(psi$so[labels]), "*")
  } else {
    for (g in seq_len(nrow(psi))) {
      idx <- which(labels == g)
      Sig <- psi$se[g] * diag(length(idx)) +
             psi$so[g] * exp(-psi$theta[g] * dm[idx, idx, drop = FALSE])
      e[, idx] <- matrix(stats::rnorm(TN * length(idx)), TN, length(idx)) %*% chol(Sig)
    }
  }
  e
}

## ---------------------------------------------------------------------------
## Replications and seeds. A STREAM is the run of one machine: its
## replications are numbered r = 1, ..., 999 within it, and the number recorded
## in the results, `rep`, is unique across streams,
##
##   rep = 1000 (stream - 1) + r ,
##
## so that two machines with different streams draw different data and their
## results stack. The seed of a replication follows from `rep` and from the
## attempt (a replication over its time limit is drawn again):
##
##   seed = 1000 rep + 1 + 100 (attempt - 1) ,
##
## and the data set uses seed, seed + 1 and seed + 2. With at most 9 attempts
## no two (rep, attempt) share a seed. Stream 1, attempt 1 gives 1000 r + 1, the
## seeds of the first design.
## ---------------------------------------------------------------------------
SIM_REP_MAX <- 999L
sim_rep_id    <- function(r, stream) as.integer((SIM_REP_MAX + 1L) * (as.integer(stream) - 1L) + r)
sim_stream_of <- function(rep) as.integer((as.integer(rep) - 1L) %/% (SIM_REP_MAX + 1L) + 1L)
sim_rep_local <- function(rep) as.integer((as.integer(rep) - 1L) %% (SIM_REP_MAX + 1L) + 1L)
sim_seed      <- function(rep, attempt = 1L)
  as.integer(1000L * as.integer(rep) + 1L + 100L * (as.integer(attempt) - 1L))

## ---------------------------------------------------------------------------
## One complete data set. Seeds are derived from the replication, so
## replication r uses the same locations and covariate in every
## scenario-variant (S2, with its single cloud, apart), and the same random
## numbers for the latent and the error, so that two variants differ by their
## definition and not by their draw.
## ---------------------------------------------------------------------------
dgp_draw <- function(cell, rep, attempt = 1L) {
  row  <- SIM_SCEN[SIM_SCEN$id == cell$scenario, , drop = FALSE]
  K    <- if (row$family == "S2") 1L else 3L
  seed <- sim_seed(rep, attempt)
  spread <- if (is.null(cell$spread) || K == 1L) "total" else cell$spread
  loc  <- dgp_locations(cell$n, K, if (K == 1L) 0 else cell$omega,
                        balance = if (K == 1L) "balanced" else cell$balance,
                        seed = seed, spread = spread)
  x    <- dgp_covariate(loc$coords, cell$TN, seed = seed + 1L)
  psi  <- dgp_psi(row)
  set.seed(seed + 2L)
  y    <- dgp_latent(cell$TN, psi, row$rho)
  e    <- dgp_error(loc$labels, psi, loc$coords, cell$TN, row$field)
  g    <- loc$labels
  mu   <- sweep(sweep(x, 2, psi$b1[g], "*"), 2, psi$b0[g], "+") + y[, g, drop = FALSE]
  list(z = mu + e, covariates = cbind(intercept = 1, xcov = as.vector(x)),
       coordinates = loc$coords, labels = g, psi = psi, K = K, row = row,
       latent = y, mu = mu, x = x)
}

## The true parameters of every regime, on the names the fitted object uses
dgp_truth <- function(psi) {
  do.call(rbind, lapply(seq_len(nrow(psi)), function(g) data.frame(
    regime    = g,
    parameter = c("beta1", "beta2", "sigma2eps", "sigma2omega", "theta", "G", "Sigmaeta"),
    truth     = c(psi$b0[g], psi$b1[g], psi$se[g], psi$so[g], psi$theta[g],
                  psi$G[g], psi$s2eta[g]),
    stringsAsFactors = FALSE)))
}


## ===========================================================================
## The cells: the union of the blocks, each the full factorial of its margins
## around the reference. In S2 the overlap and the balance are vacuous, so it
## enters once per (n, T, knn). A cell that two blocks share is kept once,
## under the first of them.
## ===========================================================================

## The blocks this invocation runs: those named in `blocks`, or, when margins
## are given on the command line, one factorial block of those margins.
sim_blocks <- function(cfg = CFG) {
  if (length(SIM_CLI_MARGINS)) return(list(command_line = cfg[SIM_CLI_MARGINS]))
  bad <- setdiff(cfg$blocks, names(SIM_BLOCKS))
  if (length(bad)) stop("unknown block: ", paste(bad, collapse = ", "),
                        "; the blocks are ", paste(names(SIM_BLOCKS), collapse = ", "),
                        call. = FALSE)
  SIM_BLOCKS[cfg$blocks]
}

sim_cells <- function(cfg = CFG, blocks = sim_blocks(cfg)) {
  bad <- setdiff(cfg$scenario, SIM_SCEN$id)
  if (length(bad)) stop("unknown scenario: ", paste(bad, collapse = ", "),
                        "; the scenarios are ", paste(SIM_SCEN$id, collapse = ", "),
                        call. = FALSE)
  cells <- do.call(rbind, lapply(names(blocks), function(b) {
    m <- utils::modifyList(SIM_REFERENCE, blocks[[b]])
    x <- expand.grid(scenario = cfg$scenario, n = m$n, TN = m$TN, omega = m$omega,
                     spread = m$spread, balance = m$balance, knn = m$knn,
                     stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
    x$block <- b
    x
  }))
  s2 <- cells$scenario == "S2"
  cells$omega[s2] <- NA_real_
  cells$spread[s2] <- "-"
  cells$balance[s2] <- "-"
  cells$cell <- sim_cell_id(cells)
  cells <- cells[!duplicated(cells$cell), , drop = FALSE]
  ok <- mapply(function(sc, n, bal) if (sc == "S2") TRUE else dgp_feasible(n, 3L, bal),
               cells$scenario, cells$n, cells$balance)
  cells <- cells[ok, , drop = FALSE]
  rownames(cells) <- NULL
  cells
}

## The primary key of a cell: a string, because omega is a double and a join on
## a floating-point column is a defect waiting to happen.
sim_cell_id <- function(cells) {
  ## omega carries a "b" when the regimes keep their spread (spread = "regime")
  sprintf("%s_n%d_T%d_w%s_%s_knn%d", cells$scenario, as.integer(cells$n),
          as.integer(cells$TN),
          ifelse(is.na(cells$omega), "NA", paste0(sprintf("%.2f", cells$omega),
                 ifelse(cells$spread %in% "regime", "b", ""))),
          substr(cells$balance, 1, 3), as.integer(cells$knn))
}

## The cost of the design: every cell priced at the mean time of its recorded
## replications, a cell not yet run at the mean over the cells that were.
sim_cost <- function(cells, measured, nrep, cores) {
  if (is.null(measured) || !nrow(measured)) return(NULL)
  measured <- measured[is.finite(measured$secs_total), ]
  if (!nrow(measured)) return(NULL)
  m <- tapply(measured$secs_total, measured$cell, mean)
  secs <- as.numeric(m[cells$cell])
  n_meas <- sum(!is.na(secs))
  secs[is.na(secs)] <- mean(m)
  list(per_cell = stats::setNames(secs, cells$cell), measured_cells = n_meas,
       one_rep_hours = sum(secs) / 3600,
       core_hours = sum(secs) * nrep / 3600,
       wall_hours = sum(secs) * nrep / 3600 / cores)
}


## ===========================================================================
## THE MONTE CARLO: one replication is the whole procedure a user would run
## ===========================================================================
sim_one <- function(cell, rep, attempt = 1L) {

  t_all <- proc.time()[["elapsed"]]
  dat <- dgp_draw(cell, rep, attempt)
  mod <- sim_model(dat)

  ## The computational settings are the defaults of the package, passed
  ## explicitly so that a session option Stem.control cannot change them. A
  ## grid with looser final refits would be another estimator: the departures
  ## that start the unpenalized fits come from the pooled fit of the grid.
  t0 <- proc.time()[["elapsed"]]
  ic <- Stem::SCSTEM_Infocrit(mod, K_grid = CFG$K_grid, phi_grid = CFG$phi_grid,
                              distance = "euclidean", verbose = FALSE,
                              knn = cell$knn, min_cluster_size = N_MIN,
                              control = Stem::STEM_control())
  sel <- Stem::SCSTEM_Select(ic)
  secs <- proc.time()[["elapsed"]] - t0

  ## SELECTION is read off the selected fit; RECOVERY off the fit at the TRUE
  ## number of regimes, so that a poor estimate is not confounded with a poor
  ## selection. Both, and the pooled benchmark, are in the grid.
  grab <- function(kk, pp) {
    j <- which(ic$table$K == kk & abs(ic$table$phi - pp) < 1e-8)
    if (!length(j)) NULL else ic$fits[[j[1]]]
  }
  K <- dat$K
  fit_sel  <- sel$fit
  fit_true <- grab(K, if (K == 1L) CFG$phi_grid[1] else CFG$phi_ref[1])
  fit_pool <- grab(1L, CFG$phi_grid[1])

  acc <- function(fit) {
    if (is.null(fit) || is.null(fit$group)) return(c(ari = NA_real_, share = NA_real_))
    g  <- fit$group
    kk <- max(max(g), K)
    map <- scstem_align_labels(reference = dat$labels, refit = g, K = kk)
    c(ari = scstem_ari(g, dat$labels), share = mean(map[g] == dat$labels, na.rm = TRUE))
  }
  a_sel  <- acc(fit_sel)
  a_true <- acc(fit_true)

  ## the fitted systematic part x'beta_g + y^(g), scored against the
  ## conditional mean mu, not against z, whose nugget is irreducible
  signal <- function(fit) {
    if (is.null(fit)) return(NULL)
    mh <- try(Stem::SCSTEM_Signal(fit), silent = TRUE)
    if (inherits(mh, "try-error")) NULL else mh
  }
  mh_true <- signal(fit_true); mh_pool <- signal(fit_pool); mh_sel <- signal(fit_sel)
  rmse <- function(mh) if (is.null(mh) || all(is.na(mh))) NA_real_ else
    sqrt(mean((mh - dat$mu)^2, na.rm = TRUE))

  row <- dat$row
  key <- data.frame(cell = cell$cell, stream = sim_stream_of(rep), rep = rep, attempt = as.integer(attempt), scenario = cell$scenario,
                    family = row$family, level = row$level, common = row$common,
                    variant = row$variant, rho = row$rho, field = row$field,
                    n = cell$n, TN = cell$TN, omega = cell$omega,
                    spread = if (is.null(cell$spread)) "total" else cell$spread,
                    balance = cell$balance, knn = cell$knn, K_true = K,
                    stringsAsFactors = FALSE)

  summ <- cbind(key, data.frame(
    K_hat = sel$K_selected, phi_hat = sel$phi_selected,
    K_correct = as.integer(sel$K_selected == K),
    ari_sel = a_sel[["ari"]],   share_sel = a_sel[["share"]],
    ari_true = a_true[["ari"]], share_true = a_true[["share"]],
    rmse_sel = rmse(mh_sel), rmse_true = rmse(mh_true), rmse_pooled = rmse(mh_pool),
    nconf = nrow(ic$table), nfail = if (is.null(ic$failed)) 0L else nrow(ic$failed),
    secs = secs, secs_total = NA_real_, error = "", stringsAsFactors = FALSE))

  ## parameters at the true number of regimes, the estimated regimes aligned
  ## on the true ones
  map <- if (is.null(fit_true) || is.null(fit_true$group)) NULL else
    scstem_align_labels(reference = dat$labels, refit = fit_true$group, K = K)
  tru <- dgp_truth(dat$psi)
  est <- rep(NA_real_, nrow(tru))
  if (!is.null(fit_true) && !is.null(fit_true$phi_hat) && !is.null(map)) {
    ph  <- fit_true$phi_hat
    est <- vapply(seq_len(nrow(tru)), function(i) {
      g_fit <- which(map == tru$regime[i])
      p <- tru$parameter[i]
      if (!length(g_fit) || !(p %in% colnames(ph))) NA_real_ else ph[g_fit[1], p]
    }, numeric(1))
  }
  params <- cbind(key[rep(1L, nrow(tru)), ], tru, estimate = est)
  rownames(params) <- NULL

  ## the fitted grid, one row per (K, phi), with the ARI of every partition, so
  ## that the recovery at any penalty, or under another selection rule, can be
  ## read afterwards without refitting
  grid <- cbind(key[rep(1L, nrow(ic$table)), ], ic$table)
  tags <- paste0("K=", ic$table$K, ", phi=", ic$table$phi)
  grid$ari <- if (K == 1L) NA_real_ else vapply(tags, function(tg)
    if (tg %in% colnames(ic$groups)) scstem_ari(ic$groups[, tg], dat$labels) else NA_real_,
    numeric(1), USE.NAMES = FALSE)
  rownames(grid) <- NULL

  ## per location
  g_true_hat <- if (is.null(map)) rep(NA_integer_, cell$n) else as.integer(map)[fit_true$group]
  g_sel_hat  <- if (is.null(fit_sel)) rep(NA_integer_, cell$n) else fit_sel$group
  colstat <- function(M) if (is.null(M)) rep(NA_real_, cell$n) else
    sqrt(colMeans((M - dat$mu)^2, na.rm = TRUE))
  station <- cbind(key[rep(1L, cell$n), ], data.frame(
    station = seq_len(cell$n),
    sx = dat$coordinates[, 1], sy = dat$coordinates[, 2],
    g_true = dat$labels, g_hat = g_true_hat, g_hat_sel = g_sel_hat,
    correct = as.integer(g_true_hat == dat$labels),
    rmse_true = colstat(mh_true), rmse_pool = colstat(mh_pool),
    stringsAsFactors = FALSE))
  rownames(station) <- NULL

  obs <- NULL
  if (sim_rep_local(rep) <= CFG$keep_obs[1]) {
    ii <- rep(seq_len(cell$n), each = cell$TN)
    tt <- rep(seq_len(cell$TN), times = cell$n)
    obs <- cbind(key[rep(1L, cell$n * cell$TN), ], data.frame(
      t = tt, station = ii, x = as.vector(dat$x), z = as.vector(dat$z),
      mu = as.vector(dat$mu), g_true = dat$labels[ii],
      g_hat = g_true_hat[ii], g_hat_sel = g_sel_hat[ii],
      mu_hat = if (is.null(mh_true)) NA_real_ else as.vector(mh_true),
      mu_hat_sel = if (is.null(mh_sel)) NA_real_ else as.vector(mh_sel),
      mu_hat_pool = if (is.null(mh_pool)) NA_real_ else as.vector(mh_pool),
      stringsAsFactors = FALSE))
    rownames(obs) <- NULL
  }

  summ$secs_total <- proc.time()[["elapsed"]] - t_all
  list(summary = summ, params = params, grid = grid, station = station, obs = obs)
}

## The row a replication leaves when it fails: the key, NA for every measure,
## and the reason.
sim_failed <- function(cell, rep, msg, attempt = 1L, secs = NA_real_) {
  row <- SIM_SCEN[SIM_SCEN$id == cell$scenario, , drop = FALSE]
  list(summary = data.frame(
    cell = cell$cell, stream = sim_stream_of(rep), rep = rep, attempt = as.integer(attempt), scenario = cell$scenario,
    family = row$family, level = row$level, common = row$common,
    variant = row$variant, rho = row$rho, field = row$field,
    n = cell$n, TN = cell$TN, omega = cell$omega,
    spread = if (is.null(cell$spread)) "total" else cell$spread, balance = cell$balance,
    knn = cell$knn, K_true = if (row$family == "S2") 1L else 3L,
    K_hat = NA_integer_, phi_hat = NA_real_, K_correct = NA_integer_,
    ari_sel = NA_real_, share_sel = NA_real_, ari_true = NA_real_, share_true = NA_real_,
    rmse_sel = NA_real_, rmse_true = NA_real_, rmse_pooled = NA_real_,
    nconf = NA_integer_, nfail = NA_integer_, secs = NA_real_, secs_total = secs,
    error = msg, stringsAsFactors = FALSE),
    params = NULL, grid = NULL, station = NULL, obs = NULL)
}

sim_one_safe <- function(cell, rep, attempt = 1L) {
  tryCatch(sim_one(cell, rep, attempt),
           error = function(e) sim_failed(cell, rep, conditionMessage(e), attempt))
}

## One replication under its time limit, in a process of its own that is killed
## when the limit is reached. Within one R process a limit cannot be enforced:
## setTimeLimit() raises an error once, which the estimation code catches as the
## failure of one configuration and goes on. The process loads this file in
## definitions-only mode and takes the settings of the session.
sim_one_capped <- function(task) {
  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch(
    callr::r(function(path, cfg, cell, rep, attempt) {
      assign("SIM_DEFINE_ONLY", TRUE, envir = globalenv())
      source(path, local = globalenv())
      assign("CFG", cfg, envir = globalenv())
      sim_one_safe(cell, rep, attempt)
    }, args = list(task$path, CFG, task$cell, task$rep, task$attempt),
    timeout = task$cap),
    error = function(e) e)
  secs <- proc.time()[["elapsed"]] - t0
  if (inherits(res, "callr_timeout_error"))
    return(sim_failed(task$cell, task$rep, sprintf("timeout after %.0f s", secs), task$attempt, secs))
  if (inherits(res, "error"))
    return(sim_failed(task$cell, task$rep, conditionMessage(res), task$attempt, secs))
  res
}


## ===========================================================================
## Running the tasks on several cores
## ===========================================================================
sim_task_fun <- function(task) {
  if (is.finite(task$cap) && !is.null(task$path) && requireNamespace("callr", quietly = TRUE))
    sim_one_capped(task) else sim_one_safe(task$cell, task$rep, task$attempt)
}

## What a core does once, when it starts: load this very file in
## definitions-only mode and take the settings of the session that launched it.
sim_worker_init <- function(path, cfg) {
  assign("SIM_DEFINE_ONLY", TRUE, envir = globalenv())
  source(path, local = globalenv())
  assign("CFG", cfg, envir = globalenv())
  invisible(TRUE)
}

## Runs `fun` on every task and hands each result to `on_result` as soon as it
## arrives; the tasks are handed out one at a time as the processes free up.
## `prepare` completes a task just before it is handed out (its time limit,
## computed on what has been recorded so far), and `on_result` may return a task
## to be appended to the queue (a replication drawn again after a timeout).
sim_run_tasks <- function(tasks, fun, on_result, cores, prepare = function(task) task) {
  n <- length(tasks)
  if (!n) return(invisible(0L))
  cores <- max(1L, min(as.integer(cores), n))
  if (cores > 1L && is.null(SIM_FILE)) {
    message("the path of this script is unknown, so it runs on one core; ",
            "open the file and press Source rather than pasting it")
    cores <- 1L
  }
  if (cores == 1L) {
    k <- 0L
    while (k < length(tasks)) {
      k <- k + 1L
      tasks[[k]] <- prepare(tasks[[k]])
      again <- on_result(fun(tasks[[k]]), tasks[[k]], k, length(tasks))
      if (!is.null(again)) tasks[[length(tasks) + 1L]] <- again
    }
    return(invisible(length(tasks)))
  }
  cat(sprintf("starting %d R processes ...\n", cores)); utils::flush.console()
  cl <- parallel::makePSOCKcluster(cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, sim_worker_init, SIM_FILE, CFG)
  send <- utils::getFromNamespace("sendCall", "parallel")
  recv <- utils::getFromNamespace("recvOneResult", "parallel")
  for (i in seq_len(cores)) {
    tasks[[i]] <- prepare(tasks[[i]])
    send(cl[[i]], fun, list(tasks[[i]]), tag = i)
  }
  nxt <- cores + 1L; inflight <- cores; k <- 0L
  while (inflight > 0L) {
    d <- recv(cl); inflight <- inflight - 1L; k <- k + 1L
    again <- on_result(d$value, tasks[[d$tag]], k, length(tasks))
    if (!is.null(again)) tasks[[length(tasks) + 1L]] <- again
    if (nxt <= length(tasks)) {
      tasks[[nxt]] <- prepare(tasks[[nxt]])
      send(cl[[d$node]], fun, list(tasks[[nxt]]), tag = nxt)
      nxt <- nxt + 1L; inflight <- inflight + 1L
    }
  }
  invisible(length(tasks))
}

## The results of every stream of a tag found in a folder, stacked: the files
## <tag>_s<stream><suffix>.csv, and <tag><suffix>.csv of the runs made before
## the streams existed. NULL when there are none.
sim_read_streams <- function(dir, tag, suffix = "") {
  pat <- sprintf("^%s(_s[0-9]+)?%s[.]csv$", gsub(".", "[.]", tag, fixed = TRUE), suffix)
  f <- list.files(dir, pattern = pat, full.names = TRUE)
  if (!length(f)) return(NULL)
  d <- lapply(f, utils::read.csv, stringsAsFactors = FALSE)
  cols <- Reduce(union, lapply(d, names))
  d <- lapply(d, function(x) { x[setdiff(cols, names(x))] <- NA; x[cols] })
  do.call(rbind, d)
}

## Replication by replication, so an interrupted run leaves whole replications
## of the design behind; within a replication in a shuffled but fixed order of
## the cells, so the running estimate of the time left is representative early.
sim_tasks <- function(cells, reps, done) {
  set.seed(20260927)
  ord <- sample(nrow(cells))
  tasks <- list()
  for (r in reps) for (i in ord) {
    cell <- cells[i, , drop = FALSE]
    if (paste(cell$cell, r, sep = "|") %in% done) next
    tasks[[length(tasks) + 1L]] <- list(cell = cell, rep = r)
  }
  tasks
}

sim_print_cost <- function(cells, prev, cfg) {
  cost <- sim_cost(cells, prev, cfg$nrep[1], cfg$cores[1])
  if (is.null(cost)) {
    cat("no replication recorded yet: run the first replication of every cell\n",
        "(rep_to = 1), then mode = \"dry\", to measure the cost\n\n", sep = "")
    return(invisible(NULL))
  }
  cat(sprintf(paste0("cost, measured on %d of %d cells: one replication of the ",
                     "design takes %s of core time;\n  %d replications take %.1f ",
                     "core-hours, %s on %d cores\n\n"),
              cost$measured_cells, nrow(cells), sim_hours(cost$one_rep_hours),
              cfg$nrep[1], cost$core_hours, sim_hours(cost$wall_hours), cfg$cores[1]))
  invisible(cost)
}


## ===========================================================================
## Main
## ===========================================================================
sim_main <- function() {

  OUT <- normalizePath(CFG$out[1], winslash = "/", mustWork = FALSE)
  sim_print_config(CFG)

  stream <- CFG$stream[1]
  rep_to <- if (is.na(CFG$rep_to[1])) CFG$nrep[1] else CFG$rep_to[1]
  if (!(stream >= 1L)) stop("stream must be a positive integer", call. = FALSE)
  if (CFG$rep_from[1] < 1L || rep_to < CFG$rep_from[1] || rep_to > SIM_REP_MAX)
    stop("replications run from rep_from to rep_to, within 1 to ", SIM_REP_MAX,
         " in a stream", call. = FALSE)
  if (CFG$cap_attempts[1] > 9L) stop("cap_attempts must be at most 9", call. = FALSE)

  ## the files of this stream; the analysis stacks every stream of the tag
  STEM_F  <- sprintf("%s_s%d", CFG$tag, stream)
  CSV     <- file.path(OUT, sprintf("%s.csv", STEM_F))
  CSV_PAR <- file.path(OUT, sprintf("%s-params.csv", STEM_F))
  CSV_GRD <- file.path(OUT, sprintf("%s-grid.csv", STEM_F))
  CSV_STA <- file.path(OUT, sprintf("%s-stations.csv", STEM_F))
  CSV_TO  <- file.path(OUT, sprintf("%s-timeouts.csv", STEM_F))
  DIR_OBS <- file.path(OUT, sprintf("%s-obs", STEM_F))

  cells <- sim_cells()
  prev  <- if (file.exists(CSV)) utils::read.csv(CSV, stringsAsFactors = FALSE) else NULL
  ## every stream of the tag found here, for the times: the cost and the time
  ## limits are measured on all of them
  all_prev <- sim_read_streams(OUT, CFG$tag)
  per_block <- table(factor(cells$block, levels = unique(cells$block)))
  cat(sprintf("%d cells (%s) x %d replications in stream %d, writing to\n  %s\n\n", nrow(cells),
              paste(names(per_block), per_block, sep = " ", collapse = ", "),
              CFG$nrep, stream, file.path(OUT, paste0(STEM_F, "*"))))

  if (identical(CFG$mode[1], "dry")) {
    print(cells[, c("block", "scenario", "n", "TN", "omega", "spread", "balance", "knn")],
          row.names = FALSE)
    cat("\n")
    sim_print_cost(cells, all_prev, CFG)
    return(invisible(cells))
  }

  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  dir.create(DIR_OBS, recursive = TRUE, showWarnings = FALSE)
  done <- if (is.null(prev)) character(0) else paste(prev$cell, prev$rep, sep = "|")

  ## The durations of the replications that ended normally, by cell and by
  ## (n, T), from which the time limits are computed; they start from what is
  ## already recorded.
  T_CELL <- new.env(); T_NT <- new.env()
  remember <- function(cell, n, TN, secs) {
    if (!is.finite(secs)) return(invisible(NULL))
    assign(cell, c(get0(cell, T_CELL, inherits = FALSE), secs), envir = T_CELL)
    key <- paste(n, TN); assign(key, c(get0(key, T_NT, inherits = FALSE), secs), envir = T_NT)
  }
  if (!is.null(all_prev)) {
    ok <- is.na(all_prev$error) | !nzchar(all_prev$error)
    for (i in which(ok)) remember(all_prev$cell[i], all_prev$n[i], all_prev$TN[i], all_prev$secs_total[i])
  }
  cap_of <- function(cell) {
    if (!(CFG$cap_mult[1] > 0)) return(Inf)
    own <- get0(cell$cell, T_CELL, inherits = FALSE)
    if (length(own) >= CFG$cap_after[1])
      return(max(CFG$cap_min[1], CFG$cap_mult[1] * stats::median(own)))
    ref <- get0(paste(cell$n, cell$TN), T_NT, inherits = FALSE)
    if (length(ref)) return(max(CFG$cap_first_min[1], CFG$cap_first_mult[1] * stats::median(ref)))
    CFG$cap_first_min[1]
  }
  prepare <- function(task) {
    if (is.null(task$attempt)) task$attempt <- 1L
    task$cap <- cap_of(task$cell)
    task$path <- SIM_FILE
    task
  }

  save_obs <- function(df, cell, rep) {
    if (is.null(df)) return(invisible(NULL))
    f <- sprintf("obs_%s_rep%04d.rds", gsub("[^A-Za-z0-9_.-]", "", cell$cell), rep)
    sim_retry(function() saveRDS(df, file.path(DIR_OBS, f), compress = "xz"))
  }

  reps  <- sim_rep_id(seq.int(CFG$rep_from[1], rep_to), stream)
  tasks <- sim_tasks(cells, reps, done)
  cat(sprintf("stream %d, replications %d to %d (recorded as rep %d to %d): %d tasks to run, %d already recorded\n\n",
              stream, CFG$rep_from[1], rep_to, min(reps), max(reps), length(tasks),
              nrow(cells) * length(reps) - length(tasks)))
  if (!length(tasks)) { cat("nothing to do\n"); return(invisible(NULL)) }

  cat(sprintf("%13s %-12s %4s %5s %4s | %8s %5s %6s | %s\n",
              "", "scenario", "n", "T", "rep", "time", "K_hat", "ARI", "elapsed, and left"))
  t_start <- proc.time()[["elapsed"]]
  on_result <- function(out, task, k, n) {
    if (!is.list(out) || is.null(out$summary))
      out <- sim_failed(task$cell, task$rep, paste(as.character(out), collapse = " "), task$attempt)
    s  <- out$summary
    el <- (proc.time()[["elapsed"]] - t_start) / 3600
    ## a replication over its time limit: recorded apart and, while attempts
    ## are left, drawn again with the seed of the next attempt (the last
    ## attempt is also written to the results, as a failed replication)
    if (grepl("^timeout", s$error))
      sim_append(data.frame(cell = task$cell$cell, stream = sim_stream_of(task$rep), rep = task$rep, attempt = task$attempt,
                            scenario = task$cell$scenario, n = task$cell$n, TN = task$cell$TN,
                            omega = task$cell$omega, spread = task$cell$spread,
                            cap = round(task$cap), secs = round(s$secs_total),
                            seed = sim_seed(task$rep, task$attempt)), CSV_TO)
    if (grepl("^timeout", s$error) && task$attempt < CFG$cap_attempts[1]) {
      cat(sprintf("[%5d/%5d] %-12s %4d %5d %4d | TIMEOUT after %.0f s (limit %.0f s), attempt %d of %d: drawn again\n",
                  k, n, s$scenario, s$n, s$TN, s$rep, s$secs_total, task$cap,
                  task$attempt, CFG$cap_attempts[1]))
      utils::flush.console()
      task$attempt <- task$attempt + 1L
      return(task)
    }
    sim_append(out$summary, CSV)
    sim_append(out$params,  CSV_PAR)
    sim_append(out$grid,    CSV_GRD)
    sim_append(out$station, CSV_STA)
    save_obs(out$obs, task$cell, task$rep)
    if (!nzchar(s$error)) remember(task$cell$cell, task$cell$n, task$cell$TN, s$secs_total)
    cat(sprintf("[%5d/%5d] %-12s %4d %5d %4d | %6.1f s %5s %6s | %s, ~%s left\n",
                k, n, s$scenario, s$n, s$TN, s$rep, s$secs_total, format(s$K_hat),
                format(round(s$ari_true, 3)), sim_hours(el), sim_hours(el / k * (n - k))))
    if (nzchar(s$error)) cat("              FAILED: ", s$error, "\n", sep = "")
    utils::flush.console()
    NULL
  }

  sim_run_tasks(tasks, sim_task_fun, on_result, CFG$cores[1], prepare = prepare)
  cat(sprintf("\ndone in %s\n\n", sim_hours((proc.time()[["elapsed"]] - t_start) / 3600)))
  sim_print_cost(cells, sim_read_streams(OUT, CFG$tag), CFG)
  invisible(NULL)
}

if (!SIM_DEFINE_ONLY) sim_main()
