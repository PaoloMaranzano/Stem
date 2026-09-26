## ===========================================================================
## SC-STEM: the simulation study. THE RUNNER.
##
## One script, self-contained. It needs nothing but R, runs from whatever folder
## it sits in -- the Google Drive folder included -- and writes its results
## beside itself, in results/. Its companion is analyse-simulations.R, which
## turns those results into the tables and figures of the paper.
##
## HOW TO RUN IT
##
##   1. Open this file in RStudio.
##   2. If you want, change the SETUP below: the levels of the design, and in
##      "4. THE RUN" the mode ("run", "coverage", "dry"), the number of cores,
##      the replications.
##   3. Press Source (Ctrl+Shift+S).
##
## The first time, Stem is installed or updated from GitHub. The replications
## are then spread over the cores, and every one prints a line with its time and
## the time left. The run can be stopped at any moment with the red Stop button
## and resumed by pressing Source again: what was recorded is skipped. Keep the
## RStudio session open while it runs.
##
## Do not run the same tag on two machines at once: they would write to the same
## files. To share the work between machines give each its own `tag` and its own
## range of replications, and analyse them together with --tag=a,b.
##
## The same settings can be given from a terminal, which overrides the SETUP:
##
##     Rscript run-simulations.R --cores=8 --rep_from=51 --rep_to=100
##     Rscript run-simulations.R --mode=dry
##
## This is part of the REPLICATION MATERIAL of the paper, not of the Stem
## package. It is deliberately disconnected from it: Stem is a dependency like
## any other, installed from CRAN or from GitHub, and nothing here is shipped
## with it.
##
## THE DESIGN IS THE SETUP BLOCK, a few screens below. Edit the levels there.
## Everything between here and it is machinery.
##
## FOUR OUTPUTS, because the questions have four shapes.
##
##   <tag>.csv           one row per replication: what the selection rule chose,
##                       how well the partition was recovered, the errors and
##                       their ratio to the pooled model, the wall time
##   <tag>-params.csv    one row per replication, regime and parameter: the
##                       truth beside the estimate, after the estimated regimes
##                       have been relocated onto the true ones. Bias and RMSE
##                       are Monte Carlo summaries of this file
##   <tag>-stations.csv  one row per replication and STATION: its coordinates,
##                       its true and estimated regime, its RMSE and MAE over
##                       time. This is the unit the tables are built on
##   <tag>-obs/*.rds     one file per kept replication, one row per station and
##                       period: x, z, the true conditional mean, the labels,
##                       and the fitted signal of the clustered and the pooled
##                       model
##
## They share the primary key (cell, rep) and join on it and on nothing else.
## Everything is appended, so an interrupted run keeps what it produced and a
## resumed run skips it.
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
## ---------------------------------------------------------------------------
## Definitions only?
##
## analyse-simulations.R sources this file to get the generator and the design,
## so that the two scripts cannot disagree about either. It sets
## SIM_DEFINE_ONLY <- TRUE first, and then nothing is run: the command line is
## not read and nothing is fitted. Stem is still required, as it is by both.
## ---------------------------------------------------------------------------
if (!exists("SIM_DEFINE_ONLY", inherits = FALSE)) SIM_DEFINE_ONLY <- FALSE



## ---------------------------------------------------------------------------
## The only prerequisite: Stem
##
## Stem is installed from GitHub (SIM_STEM_REF), and again whenever GitHub holds
## a newer commit than the installed one, so the study always runs on the
## current code; installing is the one thing the script does to the machine. It
## happens when the script is run, not when analyse-simulations.R or the worker
## processes read its definitions. Offline, an installed Stem that carries what
## the study uses is accepted as it is.
##
## For the paper SIM_STEM_REF should be pinned to the commit the study was run
## with (e.g. "PaoloMaranzano/Stem@d7cfe78"), so that a replication installs
## exactly that code.
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

SIM_STEM_REF <- "PaoloMaranzano/Stem"

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

## The installed Stem carries what the study uses and, when GitHub can be
## reached and the script is run rather than read, is GitHub's commit. The
## version number cannot tell: the development builds all say 2.0.0.
sim_stem_ok <- function(latest = NA_character_) {
  if (!requireNamespace("Stem", quietly = TRUE)) return(FALSE)
  ns <- asNamespace("Stem")
  have <- c("STEM_Signal", "SCSTEM_Signal", "SCSTEM_Infocrit", "SCSTEM_Select",
            "SCSTEM_Bootstrap", "scstem_neighbors", "scstem_align_labels",
            "scstem_ari")
  feats <- all(vapply(have, exists, logical(1), envir = ns, inherits = FALSE)) &&
    "distance" %in% names(formals(get("scstem_neighbors", envir = ns)))
  here <- utils::packageDescription("Stem")$RemoteSha
  feats && (is.na(latest) || (!is.null(here) && identical(here, latest)))
}

SIM_STEM_SHA <- if (SIM_DEFINE_ONLY) NA_character_ else sim_github_sha()
if (!sim_stem_ok(SIM_STEM_SHA)) {
  message(if (requireNamespace("Stem", quietly = TRUE))
            "a newer Stem is on GitHub: updating it" else
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
sim_require(c("geodist", "spdep"))
suppressPackageStartupMessages(library("Stem"))


## ---------------------------------------------------------------------------
## Command-line overrides
##
## Every element of the setup can be overridden with --name=value; a value is
## split on commas and coerced to the type of the default, so
##
##     Rscript run-simulations.R --nrep=25 --only_n=50,100
##
## runs 25 replications on two network sizes. Passing --name= with nothing
## after the equals sign keeps the default.
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
sim_flag <- function(x) any(commandArgs(trailingOnly = TRUE) == paste0("--", x))

sim_print_config <- function(cfg) {
  cat("configuration\n")
  for (nm in names(cfg))
    cat(sprintf("  %-14s %s\n", nm, paste(format(cfg[[nm]]), collapse = ", ")))
  cat("\n")
}


## ===========================================================================
## THE SETUP. This is the part meant to be edited.
## ===========================================================================

## ---------------------------------------------------------------------------
## 1. THE LEVELS. Every factor, every level.
## ---------------------------------------------------------------------------
SIM_LEVELS <- list(

  ## locations. The cost of a cell is cubic in this, so 400 is the ceiling
  n = c(20L, 50L, 100L, 200L, 400L),

  ## periods. 60 is a season, 120 a half-year, 365 a year of daily data
  TN = c(60L, 120L, 365L),

  ## regimes. K = 1 is the null case: one regime, so nothing to recover
  K = c(1L, 3L),

  ## Separation of the regime centres. The centres sit 2*omega apart and the
  ## dispersion within a regime is whatever is left of the fixed total variance
  ## NU_TOT, so the separation in within-regime standard deviations is
  ## 2*omega/sqrt(NU_TOT - Var(mu)) and is NOT proportional to omega. These
  ## levels are dgp_omega_for(c(0, 1.58, 3.16), K = 3) = 0, 0.686, 1 with the
  ## middle one rounded to 0.70, which costs 1.63 standard deviations against
  ## 1.58. Choose new ones with dgp_omega_for(), not by hand: rounding the
  ## middle level to 0.50 would drop it to 1.05 standard deviations, two thirds
  ## of the way back to the null case.
  omega = c(0, 0.70, 1),

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
## 2. THE REFERENCE CELL. The margins vary one factor at a time around this
##    configuration, so it should be the cell the paper talks about most.
## ---------------------------------------------------------------------------
SIM_REFERENCE <- list(
  n = 100L, TN = 120L, K = 3L, omega = 0.70,
  balance = "balanced", scenario = "S4", knn = 5L
)

## ---------------------------------------------------------------------------
## 3. THE BLOCKS. Each answers a question and can be switched off on its own
##    with --blocks=core,scenarios.
##
##    core       recovery against the three factors that govern it: how many
##               locations, how separated the regimes, how long the series.
##               Fully crossed, at the reference scenario and graph.
##    null       the same grid at K = 1, where there is no regime to find and
##               the question is whether the procedure invents one.
##    scenarios  what has to differ between regimes for the difference to be
##               found. All ten scenarios, crossed with n and omega.
##    graph      the robustness margin: knn crossed with n and omega, a dense
##               graph on a small network behaving unlike a sparse one on a
##               large network.
##    balance    the robustness margin for unequal regime sizes.
## ---------------------------------------------------------------------------
SIM_BLOCKS <- c("core", "null", "scenarios", "graph", "balance")

## ---------------------------------------------------------------------------
## 4. THE RUN. What this invocation does, all overridable on the command line.
## ---------------------------------------------------------------------------
CFG <- sim_config(args = if (SIM_DEFINE_ONLY) character(0) else commandArgs(TRUE), list(
  ## What to do when the file is sourced:
  ##   "run"       the Monte Carlo
  ##   "coverage"  the bootstrap coverage experiment
  ##   "dry"       print the design and an estimate of its cost, run nothing
  mode     = "run",

  ## How many cores to use. The default is all the physical cores of the machine
  ## but one; write a number to choose. Each core runs one replication at a time.
  cores    = sim_default_cores(),

  blocks   = SIM_BLOCKS,
  nrep     = 100L,

  ## Restrict the design to given levels without editing the block above, one
  ## option per factor. An empty option keeps every level of that factor.
  only_n        = integer(0),
  only_TN       = integer(0),
  only_K        = integer(0),
  only_omega    = numeric(0),
  only_knn      = integer(0),
  only_scenario = character(0),
  only_balance  = character(0),

  ## What the estimator searches over. phi_ref is the penalty at which parameter
  ## recovery is read off, and has to be a point of phi_grid.
  k_grid   = 1:4,
  phi_grid = c(0, 0.5, 1),
  phi_ref  = 0.5,

  ## Regularization of the regression coefficients; 0, 0 is the unpenalized
  ## estimator, which is what the paper reports.
  alpha    = 0,
  lambda   = 0,

  ## Which replications this invocation covers. A long study is executed in
  ## blocks on whatever machine is free: --rep_from=51 --rep_to=100 runs that
  ## slice and nothing else, and the files are appended, so the blocks compose.
  rep_from = 1L,
  rep_to   = 100L,

  ## How many replications per cell keep their full per-observation record.
  ## That record is n*T rows, so keeping it everywhere runs to gigabytes; the
  ## per-station summaries are kept for all replications and are what the
  ## tables are built from. Use --keep_obs=0 for none.
  keep_obs = 5L,

  ## The coverage experiment: how many bootstrap resamples, and on how many
  ## replications. Only used with --coverage.
  boot_B    = 200L,
  boot_reps = 100L,

  ## Bookkeeping. `out` defaults to a results/ folder beside this script.
  tag      = "full",
  seed0    = 1000L,
  out      = file.path(SIM_HERE, "results")
))



## ===========================================================================
## Helpers shared by the Monte Carlo and the coverage experiment
## ===========================================================================

## Two internal routines of Stem are used to score a partition against the
## truth: the relocation of labels by the majority rule and the Adjusted Rand
## Index. They are the package's own, so the study scores partitions exactly as
## the bootstrap does; being internal they are reached through the namespace.
scstem_align_labels <- utils::getFromNamespace("scstem_align_labels", "Stem")
scstem_ari          <- utils::getFromNamespace("scstem_ari", "Stem")

## The model object, with starting values from the pooled OLS fit, as a user
## would build it.
sim_model <- function(dat, n) {
  ols <- stats::lm.fit(x = dat$covariates, y = as.vector(dat$z))
  s2  <- stats::var(ols$residuals)
  Stem::STEM_Model(
    z = dat$z, covariates = dat$covariates, coordinates = dat$coordinates,
    phi = list(beta = matrix(ols$coefficients, ncol = 1),
               sigma2eps = 0.7 * s2, sigma2omega = 0.3 * s2,
               theta = 1 / 100000, G = matrix(0.8, 1, 1),
               Sigmaeta = matrix(0.2 * s2, 1, 1),
               m0 = as.matrix(0), C0 = as.matrix(1)),
    K = matrix(1, n, 1))
}

## Write, and if the file is momentarily locked -- a synchronization client such
## as Google Drive holds a file while it uploads it -- wait and try again rather
## than stop a run of several days.
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


## ===========================================================================
## THE GENERATOR. Machinery: nothing below needs editing to change the design.
##
##
## THE GEOMETRY: THE OVERLAP DESIGN
##
## The spatial configuration follows Morelli, Maranzano and Otto (2026), Spatial
## Statistics 73, 100960, Section 4. THE OVERLAP PARAMETER IS THEIR d, WHICH IS
## CALLED omega THROUGHOUT THIS DESIGN: in the package "d" already means the
## number of locations, and carrying two meanings for one letter through the
## simulation code and the paper was a defect waiting to happen.
##
## In the source design the K = 4 cluster centres sit at the corners of a square
## of half-side omega,
##
##   mu_sp = ((w,w), (-w,w), (w,-w), (-w,-w)),   Sigma_sp = nu_sp * I_2,
##
## so that omega alone controls how much the clusters overlap in space: at
## omega = 0 the Gaussians coincide and the partition has no spatial signature
## at all; as omega grows the clusters separate. Reducing the design to K = 2
## and K = 3 we keep the NEAREST-NEIGHBOUR centre distance at 2*omega rather
## than the radius, so that a given omega means the same degree of overlap
## whatever K:
##
##   K = 2   centres (-omega, 0) and (omega, 0)
##   K = 3   equilateral triangle of side 2*omega, circumradius 2*omega/sqrt(3)
##   K = 4   the square of the source paper, radius omega*sqrt(2)
##
## The dispersion within a cluster is NOT held fixed. What is held fixed is the
## TOTAL variance of a coordinate, between-centre plus within-cluster, so that
## the network covers the same area at every K and every omega:
##
##   nu_sp(K, omega) = NU_TOT - Var(mu | K, omega) ,
##
## and the standardised separation is 2*omega / sqrt(nu_sp(K, omega)), which is
## no longer proportional to omega. It still runs from 0 at omega = 0 to 3.16
## standard deviations at omega = 1, and the intermediate level of the design is
## chosen with dgp_omega_for() so that it lands on 1.58. See NU_TOT below for
## why the total rather than the within is the thing to fix.
##
## The abstract plane is mapped onto a geographic box centred on the Po Valley,
## one abstract unit being UNIT_KM kilometres. This is not cosmetic: it is what
## makes the covariance parameters interpretable. With UNIT_KM = 100 the network
## has a standard deviation of 103 km whatever the overlap, and at omega = 1 the
## centres are 200 km apart, against a baseline correlation range of 123 km --
## the regime in which a monitoring network of the Po Valley actually sits.
##
## WHY THE GENERATOR IS NOT STEM_Simulation() CALLED REGIME BY REGIME
##
## The first version of this file generated each regime by calling
## STEM_Simulation() on its own sub-network, which is what the SC-STEM model
## literally says: regime k has its own latent process y^(k), and two locations
## in different regimes are uncorrelated. That is faithful to the model and it
## is also useless as an experiment, because it makes the partition identifiable
## from the correlation structure alone, whatever the parameters. Generating
## with IDENTICAL parameters in the two regimes and one latent path per regime
## recovers the true partition with ARI = 1 in every replication; generating the
## same data with one shared latent path gives ARI = 0.14. The recovery measured
## in that design was the recovery of the latent split, not of any difference in
## Psi_k, and a design that separates the regimes only in the dynamics or only
## in the spatial covariance was testing nothing of the sort.
##
## THE FIX
##
## The coupling between the latent processes is an explicit factor of the
## design. The innovations of the K processes are drawn with cross-correlation
## rho:
##
##   rho = 1  the regimes share their dynamics, so the partition can be found
##            ONLY through the parameters -- the honest, hard case
##   rho = 0  independent processes: the model-consistent, easy case, kept as a
##            reference so that the size of the confound can be quantified
##
## With equal transition coefficients the correlation of the latent processes is
## exactly rho; with different ones it is
## rho * sqrt((1-G1^2)(1-G2^2)) / (1 - G1 G2), which is below one however large
## rho is. That is not a defect of the design but a fact about the object: two
## AR(1) processes with different persistence cannot be perfectly correlated, so
## a scenario that separates the regimes on their dynamics necessarily leaks a
## little information through the latent paths. The design measures that leak
## instead of hiding it.
##
## The measurement error is drawn as ONE global spatial field whenever the
## covariance parameters are common to the regimes, and block by regime only
## when they are not -- in which case the block structure is intrinsic to the
## scenario rather than an artefact of the generator.
##
## Effect sizes are in interpretable units: the coefficient contrast in residual
## standard deviations, the variance contrast as a nugget share, the range in
## kilometres.
## ---------------------------------------------------------------------------

## ---------------------------------------------------------------------------
## THE GENERATOR, IN FULL
##
## Write K for the number of regimes, n for the number of locations, T for the
## number of time points, g_i in {1,...,K} for the regime of location i and
## I_g = {i : g_i = g} for its index set.
##
## (1) Regime sizes.  n_1,...,n_K from dgp_sizes(): equal up to the remainder
##     under "balanced", proportional to 1:2:...:K under "unbalanced", in both
##     cases with every n_g >= N_MIN. Labels are assigned in blocks, so the
##     partition is fixed within a cell and only the coordinates are random.
##
## (2) Coordinates.  In the abstract plane of the overlap design,
##
##        s_i | g_i = g  ~  N_2( mu_g , nu_sp I_2 ) ,        independent over i,
##
##     with mu_g = dgp_centres(K, omega): centres at nearest-neighbour distance
##     2*omega, so 2*omega / sqrt(nu_sp) is the standardised separation and is
##     the same for every K. For K = 1, nu_sp is replaced by
##     nu_sp + Var(mu | K = 2, omega). The plane is mapped affinely to longitude
##     and latitude at UNIT_KM kilometres per unit, centred on (LON0, LAT0), and
##     h_ij is the geodesic distance between s_i and s_j.
##
## (3) Covariate.  One standardised, exogenous covariate, an AR(1) in time whose
##     innovations are a spatially correlated field:
##
##        x_1 = w_1 ,   x_t = a x_{t-1} + sqrt(1 - a^2) w_t ,
##        w_t ~ N_n(0, C_x) ,   (C_x)_ij = exp(-h_ij / rho_x) ,
##
##     with a = 0.7 and rho_x = 100 km, then centred and scaled over all nT
##     values. A covariate independent across stations would make a coefficient
##     contrast trivially visible; a perfectly common one would make it
##     indistinguishable from the latent process. This sits in between.
##
## (4) Latent processes.  One AR(1) per regime, with cross-correlated
##     innovations:
##
##        y_t = diag(G_1,...,G_K) y_{t-1} + eta_t ,   eta_t ~ N_K(0, Sigma_eta),
##        Sigma_eta = D R_rho D ,  D = diag(sigma_eta,1 , ... , sigma_eta,K) ,
##        (R_rho)_gh = rho for g != h and 1 for g = h ,
##        sigma^2_eta,g = vbar_y (1 - G_g^2) ,   y_1 ~ N_K(0, vbar_y R_rho) .
##
##     The innovation variance is tied to G so that every regime has the same
##     stationary variance vbar_y: a scenario that separates the regimes on G
##     then separates them on persistence alone, not on the amplitude of the
##     signal.
##
## (5) Measurement error.  When all regimes share (sigma^2_eps, sigma^2_omega,
##     theta) and force_block is FALSE, ONE global field:
##
##        e_t ~ N_n(0, Sigma) ,  Sigma = sigma^2_eps I_n + sigma^2_omega exp(-theta h),
##
##     independent over t. Otherwise one field per regime, independent across
##     regimes:
##
##        e_{t,I_g} ~ N_{n_g}(0, Sigma_g) ,
##        Sigma_g = sigma^2_eps,g I_{n_g} + sigma^2_omega,g exp(-theta_g h_{I_g,I_g}) .
##
## (6) Response.
##
##        z_ti = beta_0,g_i + beta_1,g_i x_ti + y_t,g_i + e_ti .
##
##     This is the STEM measurement equation with loading matrix K_g = 1_{n_g}
##     within each regime, which is what SCSTEM_Estimation() fits.
##
## PSEUDOCODE
##
##   input  n, T, K, omega, balance, scenario s, level l, rho, replication r
##   1  n_1..n_K  <- sizes(n, K, balance)          ; g <- labels(n_1..n_K)
##   2  mu        <- centres(K, omega)
##   3  for i in 1..n:  s_i <- mu_{g_i} + N_2(0, nu_sp I)      ; map to lon/lat
##   4  h         <- geodesic distances between the s_i
##   5  C_x       <- exp(-h / rho_x)  ;  L_x <- chol(C_x)
##      x_1 <- L_x' N(0, I) ;  for t in 2..T:  x_t <- a x_{t-1} + sqrt(1-a^2) L_x' N(0,I)
##      x   <- (x - mean(x)) / sd(x)
##   6  Psi_1..Psi_K <- parameters(s, l, K)        ; regime g at fraction (g-1)/(K-1)
##   7  Sigma_eta <- D R_rho D  ;  L_eta <- chol(Sigma_eta)
##      y_1 <- sqrt(vbar_y) chol(R_rho)' N(0,I)
##      for t in 2..T:  y_t <- diag(G) y_{t-1} + L_eta' N(0, I)
##   8  if the regimes share the covariance and not force_block:
##         L <- chol(Sigma) ;  for t in 1..T:  e_t <- L' N(0, I)
##      else for g in 1..K:
##         L_g <- chol(Sigma_g) ;  for t in 1..T:  e_{t,I_g} <- L_g' N(0, I)
##   9  for t in 1..T, i in 1..n:
##         z_ti <- beta_0,g_i + beta_1,g_i x_ti + y_t,g_i + e_ti
##   output z (T x n), X = [1, vec(x)], coordinates, g, Psi_1..Psi_K
##
## Seeds are derived from the replication index, so a cell is reproducible on
## its own and the same replication uses the same geometry across scenarios.
## ---------------------------------------------------------------------------

## The TOTAL variance of each coordinate: between-centre plus within-regime.
## It is held fixed across K and omega, and the within-regime dispersion is
## whatever is left over,
##
##     nu_sp(K, omega) = NU_TOT - Var(mu | K, omega) .
##
## WHY IT IS THE TOTAL, AND NOT THE WITHIN, THAT IS FIXED. With the within-regime
## dispersion held fixed instead -- centres 2*omega apart around clusters of
## constant spread -- the map grows with the separation, and omega then does two
## things at once: it separates the regimes AND it enlarges the network. The
## second is not innocuous here, because the covariance has a range: at omega = 0
## and n = 20 the median pairwise distance was 105 km against a true range of
## 123 km, so the exponential decay was barely resolved over the observed
## distances, theta was close to unidentified, and the EM crawled -- 0.4 s for a
## pooled fit at omega = 1 against more than ten minutes at omega = 0, on the
## same n and T. Any effect attributed to the separation would have carried a
## share of that. Holding the total fixed leaves the footprint of the network,
## and so the identifiability of theta, the same in every cell.
##
## The value is the total the design used to have at its most separated cell,
## K = 3 and omega = 1, so that cell is unchanged and the smaller overlaps are
## the ones that widen.
NU_TOT  <- 0.4 + 2/3    # = NU_SP + dgp_centre_var(3, 1)
NU_SP   <- 0.4          # kept for reference: the former within-cluster variance
UNIT_KM <- 100          # kilometres per abstract unit of the overlap design
LON0    <- 9.5          # centre of the geographic box, Po Valley
LAT0    <- 45.5
## The smallest regime the design allows. The package itself refuses fewer than
## r + 2 units, but the binding constraint is spatial rather than parametric: a
## regime estimates a range from its own pairwise distances, and six locations
## already give fifteen of them. Below that the range is not identified in any
## practical sense and a failure to recover the partition would be a failure to
## fit, not a failure to separate.
N_MIN   <- 6L
## The smallest ratio of largest to smallest regime that still counts as an
## unbalanced design. Below it the correction for N_MIN has flattened the
## allocation back to nearly equal sizes, and the cell would be a duplicate of
## the balanced one under a different name.
IMB_MIN <- 1.5

## ---------------------------------------------------------------------------
## Geometry
## ---------------------------------------------------------------------------

## centres of the K clusters at overlap omega, nearest-neighbour distance 2*omega
dgp_centres <- function(K, omega) {
  if (K == 1L) return(cbind(0, 0))
  if (K == 2L) return(cbind(c(-omega, omega), c(0, 0)))
  if (K == 3L) {
    r <- 2 * omega / sqrt(3); a <- c(90, 210, 330) * pi / 180
    return(cbind(r * cos(a), r * sin(a)))
  }
  if (K == 4L) return(cbind(c(omega, -omega, omega, -omega), c(omega, omega, -omega, -omega)))
  stop("K must be 1, 2, 3 or 4")
}

## Mean per-coordinate variance of the K centres. With K = 1 the design has no
## between-cluster spread, so the WITHIN-cluster dispersion is inflated by this
## amount: the pooled configuration then covers the same area as the clustered
## one at the same overlap, and the selection rule is not handed a free
## geometric cue for telling k = 1 from k > 1.
dgp_centre_var <- function(K, omega) {
  mu <- dgp_centres(K, omega)
  if (K == 1L) return(0)
  mean(apply(mu, 2, function(v) mean((v - mean(v))^2)))
}

## What is left of the total variance for the dispersion within a regime, and
## the separation of the centres it implies, in within-regime standard
## deviations. The second is the quantity the design is really indexed by:
## omega is the knob, this is what it means.
## dgp_centres() takes one omega at a time, so these vectorise over it
dgp_nu_sp <- function(K, omega, nu_tot = NU_TOT)
  nu_tot - vapply(omega, function(w) dgp_centre_var(K, w), numeric(1))

dgp_separation <- function(K, omega, nu_tot = NU_TOT) {
  if (K == 1L) return(rep(0, length(omega)))
  d <- dgp_nu_sp(K, omega, nu_tot)
  ifelse(d > 0, 2 * omega / sqrt(d), NA_real_)
}

## The overlap that delivers a wanted separation. Because the dispersion now
## depends on omega, the two are no longer proportional: with
## Var(mu) = c_K omega^2 the separation is 2 omega / sqrt(NU_TOT - c_K omega^2),
## which inverts to the expression below. This is what the levels of omega in
## design.R are chosen with.
dgp_omega_for <- function(sep, K, nu_tot = NU_TOT) {
  cK <- if (K == 1L) 0 else dgp_centre_var(K, 1)     # Var(mu) at omega = 1
  sep * sqrt(nu_tot / (4 + sep^2 * cK))
}

## Regime sizes.
##
##   "balanced"    exactly equal, the remainder spread over the first regimes
##   "unbalanced"  sizes proportional to 1 : 2 : ... : K, so the largest regime
##                 is K times the smallest, then corrected so that no regime
##                 falls below n_min, the excess being taken from the largest
##
## The correction is what makes the unbalanced case usable: SC-STEM fits a full
## STEM model inside every regime, so a regime with fewer units than the r + 6
## parameters it has to estimate is not a hard case, it is an infeasible one,
## and a design that produced those would confound the recovery of the partition
## with the feasibility of the fit. `dgp_feasible()` says whether a cell admits
## the requested imbalance at all.
dgp_sizes <- function(n, K, balance = c("balanced", "unbalanced"),
                      n_min = N_MIN) {
  balance <- match.arg(balance)
  if (K == 1L) return(n)
  if (balance == "balanced") {
    s <- rep(n %/% K, K)
    if (n %% K) s[seq_len(n %% K)] <- s[seq_len(n %% K)] + 1L
    return(as.integer(s))
  }
  w <- seq_len(K) / sum(seq_len(K))
  s <- pmax(1L, as.integer(round(n * w)))
  s[K] <- n - sum(s[-K])
  short <- pmax(0L, n_min - s)
  if (any(short)) {
    s <- pmax(s, n_min)
    s[K] <- n - sum(s[-K])
  }
  as.integer(s)
}

## TRUE when the cell can carry K regimes of at least n_min units each AND, for
## an unbalanced design, when the imbalance survives that constraint
dgp_feasible <- function(n, K, balance = "balanced", n_min = N_MIN,
                         imb_min = IMB_MIN) {
  s <- dgp_sizes(n, K, balance, n_min)
  ok <- all(s >= n_min) && sum(s) == n
  if (ok && K > 1L && balance == "unbalanced") ok <- max(s) / min(s) >= imb_min
  ok
}

dgp_labels <- function(n, K, balance = "balanced", n_min = N_MIN) {
  rep(seq_len(K), times = dgp_sizes(n, K, balance, n_min))
}

## Locations: the abstract cloud of the overlap design, mapped to longitude and
## latitude so that distances are kilometres and the covariance parameters keep
## their meaning.
##
##   s_i | g_i = g  ~  N_2( mu_g , nu_sp(K, omega) I_2 )
##
## with mu_g = dgp_centres(K, omega) and nu_sp(K, omega) = NU_TOT minus the
## between-centre variance, so that every cell of the design covers the same
## area whatever K and omega are. K = 1 needs no special case: it has no
## between-centre variance, so it takes the whole of NU_TOT.
dgp_locations <- function(n, K, omega, nu_tot = NU_TOT, balance = "balanced",
                          seed = 1) {
  set.seed(seed)
  mu <- dgp_centres(K, omega)
  g  <- dgp_labels(n, K, balance)
  nu_sp <- nu_tot - dgp_centre_var(K, omega)
  if (nu_sp <= 0) {
    stop("the centres of K = ", K, " at omega = ", omega, " already spread more ",
         "than the total variance NU_TOT = ", signif(nu_tot, 4),
         ": raise NU_TOT or lower omega", call. = FALSE)
  }
  sd_i <- sqrt(nu_sp)
  xy <- cbind(mu[g, 1] + stats::rnorm(n, sd = sd_i),
              mu[g, 2] + stats::rnorm(n, sd = sd_i))
  coords <- cbind(
    lon = LON0 + xy[, 1] * UNIT_KM / (111.320 * cos(LAT0 * pi / 180)),
    lat = LAT0 + xy[, 2] * UNIT_KM / 110.574)
  list(coords = coords, labels = g, xy = xy, mu = mu, sizes = tabulate(g, K))
}

## ---------------------------------------------------------------------------
## Baseline parameters, from the pooled fit on the real network
## ---------------------------------------------------------------------------
dgp_base <- function() {
  list(beta        = c(2.34, 0.60),   # intercept and one standardised covariate
       sigma2eps   = 21.4,            # nugget
       sigma2omega = 4.19,            # partial sill
       theta       = 8.1e-06,         # 1/theta = 123 km
       G           = 0.90,            # persistence of the latent process
       var_y       = 6.0)             # stationary variance of the latent process
}

## The innovation variance that holds the stationary variance of an AR(1) at
## var_y, so that a scenario separating the regimes on G separates them on
## persistence alone and not on the amplitude of the latent signal.
sigma_eta_of <- function(G, var_y) var_y * (1 - G^2)

## ---------------------------------------------------------------------------
## A standardised, exogenous covariate with realistic structure: an AR(1) in
## time whose innovations are a spatially correlated field. A covariate that was
## independent across stations would make a coefficient contrast trivially
## visible; one that was perfectly common would make it indistinguishable from
## the latent process. This sits in between, and the two knobs are explicit.
## ---------------------------------------------------------------------------
dgp_covariate <- function(coords, TN, a_time = 0.7, range_km = 100, seed = 1) {
  d <- nrow(coords)
  dm <- geodist::geodist(coords, measure = "geodesic")
  Cx <- exp(-dm / (range_km * 1000))
  L <- chol(Cx + diag(1e-8, d))

  set.seed(seed)
  x <- matrix(NA_real_, TN, d)
  w <- as.numeric(crossprod(L, stats::rnorm(d)))
  x[1, ] <- w
  for (tt in 2:TN) {
    w <- as.numeric(crossprod(L, stats::rnorm(d)))
    x[tt, ] <- a_time * x[tt - 1, ] + sqrt(1 - a_time^2) * w
  }
  (x - mean(x)) / stats::sd(x)
}

## the covariate stacked by station, as STEM_Model() expects
dgp_design <- function(x) {
  cbind(intercept = 1, xcov = as.vector(x))
}

## ---------------------------------------------------------------------------
## Parameter sets of the K regimes
##
## `scenario` says which component separates them, `level` how far apart the two
## EXTREME regimes are. Regime 1 always carries the baseline and regime K the
## full contrast; with K = 3 the middle regime sits halfway, so that K = 2
## reproduces exactly the earlier design and K = 3 extends it without changing
## the meaning of `level`.
## ---------------------------------------------------------------------------
dgp_psi <- function(scenario, level = 2L, K = 2L, base = dgp_base()) {

  res_sd <- sqrt(base$sigma2eps + base$sigma2omega)
  tot    <- base$sigma2eps + base$sigma2omega

  ## the contrast of the extreme regime
  delta_beta <- c(0.25, 0.5, 1.0)[level]   # in residual standard deviations
  G_last     <- c(0.80, 0.60, 0.30)[level] # against 0.90 in the baseline
  share_last <- c(0.70, 0.50, 0.30)[level] # nugget share, against 0.836
  range_last <- c(60, 30, 15)[level]       # km, against 123 km

  ## regime g sits at fraction w of the way from the baseline to the extreme
  w_of <- function(g) if (K == 1L) 0 else (g - 1) / (K - 1)

  lapply(seq_len(K), function(g) {
    w <- w_of(g)
    p <- base
    if (scenario %in% c("S1", "S4")) {
      p$beta <- c(base$beta[1], base$beta[2] + w * delta_beta * res_sd)
    }
    if (scenario %in% c("S2", "S4")) {
      p$G <- base$G + w * (G_last - base$G)
    }
    if (scenario %in% c("S3", "S4")) {
      share <- (base$sigma2eps / tot) + w * (share_last - base$sigma2eps / tot)
      rng   <- (1 / base$theta / 1000) + w * (range_last - 1 / base$theta / 1000)
      p$sigma2eps   <- tot * share
      p$sigma2omega <- tot * (1 - share)
      p$theta       <- 1 / (rng * 1000)
    }
    p
  })
}

## TRUE when all regimes share every parameter of the measurement covariance, in
## which case the error field can be drawn globally
dgp_common_field <- function(psi) {
  key <- function(p) unlist(p[c("sigma2eps", "sigma2omega", "theta")])
  all(vapply(psi[-1], function(p) isTRUE(all.equal(key(p), key(psi[[1]]))),
             logical(1)))
}

## ---------------------------------------------------------------------------
## The generator
##
##   labels  regime of each location
##   psi     list of parameter sets, one per regime
##   rho     cross-correlation of the innovations of the latent processes
##   force_block  draw the measurement error block by regime even when the
##                covariance parameters are common: this reproduces the first
##                design and is what makes the confound measurable
## ---------------------------------------------------------------------------
dgp_simulate <- function(labels, psi, x, coords, rho = 1, seed = 1,
                         force_block = FALSE) {

  TN <- nrow(x); d <- ncol(x)
  k  <- length(psi)
  dm <- geodist::geodist(coords, measure = "geodesic")

  set.seed(seed)

  ## ---- latent processes, with innovations correlated at rho ----------------
  Gs   <- vapply(psi, function(p) p$G, numeric(1))
  vy   <- psi[[1]]$var_y
  sds  <- sqrt(vapply(Gs, sigma_eta_of, numeric(1), var_y = vy))
  R    <- matrix(rho, k, k); diag(R) <- 1
  Seta <- diag(sds, k) %*% R %*% diag(sds, k)
  Leta <- chol(Seta + diag(1e-10, k))

  y <- matrix(NA_real_, TN, k)
  y[1, ] <- sqrt(vy) * as.numeric(crossprod(chol(R + diag(1e-10, k)),
                                            stats::rnorm(k)))
  for (tt in 2:TN) {
    eta <- as.numeric(crossprod(Leta, stats::rnorm(k)))
    y[tt, ] <- Gs * y[tt - 1, ] + eta
  }

  ## ---- measurement error ---------------------------------------------------
  e <- matrix(NA_real_, TN, d)
  if (dgp_common_field(psi) && !force_block) {
    ## one global spatial field: the partition leaves no trace in the residual
    ## covariance, which is what makes S0 a genuine null
    p1 <- psi[[1]]
    Sig <- p1$sigma2eps * diag(d) + p1$sigma2omega * exp(-p1$theta * dm)
    L <- chol(Sig)
    for (tt in seq_len(TN)) e[tt, ] <- as.numeric(crossprod(L, stats::rnorm(d)))
  } else {
    for (g in seq_len(k)) {
      idx <- which(labels == g)
      p <- psi[[g]]
      Sig <- p$sigma2eps * diag(length(idx)) +
             p$sigma2omega * exp(-p$theta * dm[idx, idx, drop = FALSE])
      L <- chol(Sig)
      for (tt in seq_len(TN)) {
        e[tt, idx] <- as.numeric(crossprod(L, stats::rnorm(length(idx))))
      }
    }
  }

  ## ---- assemble ------------------------------------------------------------
  ## `mu` is the conditional mean given the latent path, which is what a fitted
  ## value estimates: it is the target against which predictive accuracy is
  ## measured, and it is not recoverable from z once the noise is added.
  z  <- matrix(NA_real_, TN, d)
  mu <- matrix(NA_real_, TN, d)
  for (i in seq_len(d)) {
    g <- labels[i]
    mu[, i] <- psi[[g]]$beta[1] + psi[[g]]$beta[2] * x[, i] + y[, g]
    z[, i]  <- mu[, i] + e[, i]
  }
  attr(z, "latent") <- y
  attr(z, "mu")     <- mu
  z
}

## ---------------------------------------------------------------------------
## The scenarios: what separates the regimes, and by how much
##
## rho = 1 throughout except for the reference cell S5, which reproduces the
## model-consistent generator of the first design so that the contribution of
## the latent split can be read off directly.
## ---------------------------------------------------------------------------
dgp_scenarios <- function() {
  rbind(
    data.frame(id = "S0", scenario = "S0", level = 2L, rho = 1, force_block = FALSE,
               label = "null: no regime at all"),
    data.frame(id = "S1a", scenario = "S1", level = 1L, rho = 1, force_block = FALSE,
               label = "coefficients, 0.25 residual sd"),
    data.frame(id = "S1b", scenario = "S1", level = 2L, rho = 1, force_block = FALSE,
               label = "coefficients, 0.50 residual sd"),
    data.frame(id = "S1c", scenario = "S1", level = 3L, rho = 1, force_block = FALSE,
               label = "coefficients, 1.00 residual sd"),
    data.frame(id = "S2a", scenario = "S2", level = 1L, rho = 1, force_block = FALSE,
               label = "persistence, G = 0.90 vs 0.80"),
    data.frame(id = "S2b", scenario = "S2", level = 3L, rho = 1, force_block = FALSE,
               label = "persistence, G = 0.90 vs 0.30"),
    data.frame(id = "S3a", scenario = "S3", level = 1L, rho = 1, force_block = FALSE,
               label = "covariance, nugget share .84 vs .70, range 123 vs 60 km"),
    data.frame(id = "S3b", scenario = "S3", level = 3L, rho = 1, force_block = FALSE,
               label = "covariance, nugget share .84 vs .30, range 123 vs 15 km"),
    data.frame(id = "S4", scenario = "S4", level = 2L, rho = 1, force_block = FALSE,
               label = "all three contrasts, middle level"),
    data.frame(id = "S5", scenario = "S4", level = 2L, rho = 0, force_block = TRUE,
               label = "reference: the model-consistent generator"),
    stringsAsFactors = FALSE
  )
}

## ---------------------------------------------------------------------------
## The dimensions of the design are NOT here. They live in design.R, which is
## the single place where the factors and their levels are written down, so
## that the driver and the table of the paper cannot disagree about what was
## run. This file is the generator: given a cell, it produces the data.
##
## For the record of what the levels mean: T is read as a real observation
## window, 60 and 120 being five and ten years of monthly data and 365 one year
## of daily data; n spans the sizes a regional network actually takes, Northern
## Italy carrying about 260 air quality stations, so 20 and 50 are a small
## sub-network, 100 a regional one, 200 and 400 a national or multi-regional
## one.
## ---------------------------------------------------------------------------

## ---------------------------------------------------------------------------
## One complete data set of the design, ready for STEM_Model()
## ---------------------------------------------------------------------------
dgp_draw <- function(n, TN, K, omega, scenario_row, rep = 1L,
                     balance = "balanced") {
  seed <- 1000L * rep + 1L
  loc  <- dgp_locations(n, K, omega, balance = balance, seed = seed)
  x    <- dgp_covariate(loc$coords, TN, seed = seed + 1L)
  psi  <- dgp_psi(scenario_row$scenario, scenario_row$level, K = K)
  z    <- dgp_simulate(loc$labels, psi, x, loc$coords,
                       rho = scenario_row$rho, seed = seed + 2L,
                       force_block = scenario_row$force_block)
  list(z = z, covariates = dgp_design(x), coordinates = loc$coords,
       labels = loc$labels, psi = psi, sizes = loc$sizes,
       latent = attr(z, "latent"), mu = attr(z, "mu"), x = x)
}

## ---------------------------------------------------------------------------
## The true parameters of a regime, flattened onto the names the fitted object
## uses, so that estimate and truth can be compared element by element. The fit
## reports sigma2eps and sigma2omega separately, and theta rather than a range.
## ---------------------------------------------------------------------------
dgp_truth <- function(psi) {
  do.call(rbind, lapply(seq_along(psi), function(g) {
    p <- psi[[g]]
    data.frame(
      regime    = g,
      parameter = c("beta1", "beta2", "sigma2eps", "sigma2omega", "theta",
                    "G", "Sigmaeta", "m0"),
      truth     = c(p$beta[1], p$beta[2], p$sigma2eps, p$sigma2omega, p$theta,
                    p$G, sigma_eta_of(p$G, p$var_y), 0),
      stringsAsFactors = FALSE)
  }))
}


## ---------------------------------------------------------------------------
## ASSEMBLING THE DESIGN
##
## The LEVELS of the design are not here: they are the setup block of
## run-simulations.R, which is the file meant to be edited. What is here is the
## machinery that turns a set of levels into a list of cells, and the cost model
## that prices them, neither of which changes when the levels do.
##
## The study is deliberately not a full factorial. Crossing every factor with
## every other is 2700 cells at K = 3 alone, most of them answering no question.
## It is a core factorial in the factors that interact -- how many locations,
## how separated the regimes, how long the series -- plus one-factor-at-a-time
## margins around a reference cell for the factors that are there to show the
## results do not turn on them.
## ---------------------------------------------------------------------------
sim_cells <- function(lv = SIM_LEVELS,
                      ref = SIM_REFERENCE,
                      blocks = SIM_BLOCKS) {

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
## 5. What a cell costs.
##
## Seconds for ONE replication -- one (k, phi) grid through SCSTEM_Infocrit()
## and SCSTEM_Select() -- on one core, at k_grid = 1:4 and phi_grid of length 3,
## interpolated in log n and in T for sizes not in the table. The estimate is
## meant for budgeting, not for reporting; once replications are recorded, the
## "dry" mode recalibrates it on them.
## ---------------------------------------------------------------------------
## MEASURED on 2026-09-26, one replication of every cell of the core and null
## blocks, on the author's machine, averaged over the three overlaps. The times
## are heavy-tailed -- a replication in which the EM struggles can take ten times
## the typical one -- so a single replication per cell is a rough guide. The
## cells with n = 400 and T = 60 are the extreme case: there the pooled fit
## alone took more than twenty minutes and a whole replication more than an
## hour without finishing, so they are entered at one hour, a lower bound.
STEM_COST <- data.frame(
  n    = rep(c(20L, 50L, 100L, 200L, 400L), times = 6L),
  K    = rep(c(1L, 3L), each = 15L),
  TN   = rep(rep(c(60L, 120L, 365L), each = 5L), times = 2L),
  secs = c(  3.3, 17.2, 13.5,  24.9, 3600,     # K = 1, T = 60
             4.9, 17.6, 20.0,  27.4,  97.0,    # K = 1, T = 120
            15.7, 18.8, 28.9,  65.3, 159.0,    # K = 1, T = 365
             4.2, 43.6, 31.1,  26.5, 3600,     # K = 3, T = 60
             7.5, 32.8, 37.4,  37.1,  86.8,    # K = 3, T = 120
            19.0, 58.5, 42.7, 247.6, 122.6)    # K = 3, T = 365
)

## seconds for one replication of one cell, interpolated in log(n) and linear
## in T within the measured envelope, held flat outside it
sim_cell_secs <- function(n, TN, K) {
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

## THE TABLE ABOVE IS A PRIOR, AND A POOR ONE. The cost of a replication is
## heavy-tailed: it depends on how hard the EM works on that particular draw.
## Measured on one cell, n = 50 and T = 60, the pooled fit took 8.4 s on one
## replication and 0.2 s on the next, and one replication of the whole grid
## took 169 s against 25 s for its neighbours. A table built from a handful of
## replications cannot see a tail like that. So whenever results already exist
## the budget is recalibrated from them: a cell that has been run is priced at
## the mean of its own recorded times, and a cell that has not is priced by the
## table scaled by the ratio of measured to predicted over the cells that have.
## Run the first few replications of the whole design, then --dry, and the
## figure printed is one to plan on.
sim_cost <- function(cells, nrep = 100L, cores = 12L, measured = NULL) {
  prior <- sim_cell_secs(cells$n, cells$TN, cells$K)
  secs  <- prior
  calib <- NA_real_; n_meas <- 0L
  if (!is.null(measured) && nrow(measured)) {
    measured <- measured[is.finite(measured$secs), ]
    key <- function(d) paste(d$n, d$TN, d$K, sprintf("%.2f", d$omega), d$id,
                             substr(d$balance, 1, 3), d$knn)
    m <- tapply(measured$secs, key(measured), mean)
    hit <- key(cells) %in% names(m)
    if (any(hit)) {
      obs <- as.numeric(m[key(cells)[hit]])
      calib <- sum(obs) / sum(prior[hit])
      secs[hit]  <- obs
      secs[!hit] <- prior[!hit] * calib
      n_meas <- nrow(measured)
    }
  }
  tot <- sum(secs) * nrep
  list(cells = nrow(cells), nrep = nrep, cores = cores,
       core_hours = tot / 3600, wall_hours = tot / 3600 / cores,
       one_rep_hours = sum(secs) / 3600,
       calibration = calib, measured_reps = n_meas,
       measured_cells = if (is.na(calib)) 0L else sum(key(cells) %in% names(m)))
}


## ---------------------------------------------------------------------------


## ===========================================================================
## THE MONTE CARLO
##
## One replication is the whole procedure a user would run: build the model
## object, fit the (k, phi) grid with SCSTEM_Infocrit(), apply the two-step rule
## with SCSTEM_Select(). What is recorded is described at the top of the file.
##
## The work is a list of TASKS, one per (cell, replication). The tasks are
## handed to the cores one at a time, as each core frees up, so a replication
## that happens to be slow holds up one core and not the others; every result
## comes back to this R session, which alone writes the files and prints the
## progress. A slow cell -- n = 400 with T = 60 can take an hour for a single
## replication -- therefore costs its own time and nothing more.
## ===========================================================================
SIM_SCEN <- dgp_scenarios()

## ---------------------------------------------------------------------------
## One replication
## ---------------------------------------------------------------------------
sim_one <- function(cell, rep) {

  row <- SIM_SCEN[SIM_SCEN$id == cell$id, , drop = FALSE]
  dat <- dgp_draw(cell$n, cell$TN, cell$K, cell$omega, row, rep = rep,
                  balance = cell$balance)

  mod <- sim_model(dat, cell$n)

  t0 <- proc.time()[["elapsed"]]
  ## The estimator is held to the same floor as the generator. N_MIN is the
  ## smallest regime in which a range is identified in any practical sense, so
  ## a configuration below it is one the design itself calls unidentified --
  ## and fitting it is not merely uninformative but ruinously slow: at n = 20
  ## a k = 4 fit, with regimes of four and five locations, took 364 s where
  ## k = 3 took one. With the floor the package refuses such a k at once, the
  ## grid records it among the failed configurations, and the selection rule
  ## never sees it.
  ic <- Stem::SCSTEM_Infocrit(mod, k_grid = CFG$k_grid, phi_grid = CFG$phi_grid,
                        distance = "geo", verbose = FALSE, knn = cell$knn,
                        min_cluster_size = N_MIN,
                        alpha = CFG$alpha[1], lambda = CFG$lambda[1])
  sel <- Stem::SCSTEM_Select(ic)
  secs <- proc.time()[["elapsed"]] - t0

  ## ------------------------------------------------------------------------
  ## Two questions, two fits, deliberately separated.
  ##
  ## SELECTION asks whether the rule finds the truth, and is read off the fit
  ## the rule chose. RECOVERY asks whether the estimator gets the parameters
  ## right when it is told the truth, and is read off the fit at the TRUE
  ## number of regimes. Measuring recovery on the selected fit would confound
  ## the two: a poor estimate would be indistinguishable from a poor selection.
  ## Both fits are already in the grid, so neither costs an extra run.
  ## ------------------------------------------------------------------------
  grab <- function(kk, pp) {
    j <- which(ic$table$k == kk & abs(ic$table$phi - pp) < 1e-8)
    if (!length(j)) NULL else ic$fits[[j[1]]]
  }
  fit_sel  <- sel$fit
  ## with one regime the penalty is vacuous, and the grid holds k = 1 only at
  ## the first value of phi
  fit_true <- grab(cell$K, if (cell$K == 1L) CFG$phi_grid[1] else CFG$phi_ref[1])
  fit_pool <- grab(1L, CFG$phi_grid[1])

  ## Clustering accuracy. The ARI is invariant to label switching; the share of
  ## correctly assigned locations is not, so the labels are first relocated onto
  ## the truth by the majority rule, exactly as the bootstrap does.
  acc <- function(fit) {
    if (is.null(fit) || is.null(fit$group)) return(c(ari = NA_real_, share = NA_real_))
    g  <- fit$group
    kk <- max(max(g), cell$K)
    map <- scstem_align_labels(reference = dat$labels, refit = g, K = kk)
    c(ari   = scstem_ari(g, dat$labels),
      share = mean(map[g] == dat$labels, na.rm = TRUE))
  }
  a_sel  <- acc(fit_sel)
  a_true <- acc(fit_true)

  ## Predictive accuracy against the CONDITIONAL MEAN, not against z: the noise
  ## is irreducible, and scoring against z would compress every comparison
  ## towards one.
  ##
  ## NOTE. This is SCSTEM_Signal(), NOT SCSTEM_Complete(), and the distinction
  ## matters. SCSTEM_Complete() returns E[z | observed]: it fills the gaps and
  ## therefore returns z itself wherever z was observed, so on complete data it
  ## IS the data and scoring it against anything measures nothing. What is
  ## wanted here is the systematic part the model fits,
  ##
  ##     muhat_ti = x_ti' betahat_g + K_i yhat_t^(g) .
  signal <- function(fit) {
    if (is.null(fit)) return(NULL)
    mh <- try(Stem::SCSTEM_Signal(fit), silent = TRUE)
    if (inherits(mh, "try-error")) NULL else mh
  }
  rmse <- function(fit) {
    mh <- signal(fit)
    if (is.null(mh) || all(is.na(mh))) return(NA_real_)
    sqrt(mean((mh - dat$mu)^2, na.rm = TRUE))
  }
  r_true <- rmse(fit_true)
  r_pool <- rmse(fit_pool)

  ## THE PRIMARY KEY. Every one of the four outputs carries `cell` and `rep`,
  ## and the pair identifies a run: the four files join on it and on nothing
  ## else. `cell` is a string rather than the tuple of factors because omega is
  ## a double -- 2/3 does not survive a round trip through a CSV exactly, and a
  ## join on a floating-point column is a defect waiting to happen. The factor
  ## columns are kept beside it for filtering, not for joining.
  cell_id <- sprintf("n%d_T%d_K%d_w%.2f_%s_%s_knn%d_a%g_l%g",
                     cell$n, cell$TN, cell$K, cell$omega, cell$id,
                     substr(cell$balance, 1, 3), cell$knn,
                     CFG$alpha[1], CFG$lambda[1])

  key <- data.frame(cell = cell_id, rep = rep,
                    n = cell$n, TN = cell$TN, K = cell$K, omega = cell$omega,
                    id = cell$id, balance = cell$balance, knn = cell$knn,
                    alpha = CFG$alpha[1], lambda = CFG$lambda[1],
                    stringsAsFactors = FALSE)

  summ <- cbind(key, data.frame(
    k_hat = sel$k_selected, phi_hat = sel$phi_selected,
    k_correct = as.integer(sel$k_selected == cell$K),
    ari_sel = a_sel[["ari"]],   share_sel = a_sel[["share"]],
    ari_true = a_true[["ari"]], share_true = a_true[["share"]],
    rmse_true = r_true, rmse_pooled = r_pool, rmse_ratio = r_true / r_pool,
    nconf = nrow(ic$table), nfail = nrow(ic$failed),
    secs = secs, stringsAsFactors = FALSE))

  ## ------------------------------------------------------------------------
  ## Parameter recovery, regime by regime, at the true number of regimes and
  ## after relocating the estimated labels onto the true ones. What is stored is
  ## the estimate beside the truth, one row each: bias, RMSE and coverage are
  ## Monte Carlo summaries of this file, not quantities a single replication
  ## could compute.
  ## ------------------------------------------------------------------------
  ## the relocation of the estimated regimes onto the true ones, used both here
  ## and by the per-station record below
  map <- if (is.null(fit_true) || is.null(fit_true$group)) NULL else
    scstem_align_labels(reference = dat$labels, refit = fit_true$group,
                        K = cell$K)

  tru <- dgp_truth(dat$psi)
  est <- rep(NA_real_, nrow(tru))
  if (!is.null(fit_true) && !is.null(fit_true$phi_hat) && !is.null(map)) {
    ph  <- fit_true$phi_hat
    est <- vapply(seq_len(nrow(tru)), function(i) {
      g_fit <- which(map == tru$regime[i])   # the fitted regime playing that role
      p     <- tru$parameter[i]
      if (!length(g_fit) || !(p %in% colnames(ph))) NA_real_ else ph[g_fit[1], p]
    }, numeric(1))
  }
  params <- cbind(key[rep(1L, nrow(tru)), ], tru, estimate = est)
  rownames(params) <- NULL

  ## ------------------------------------------------------------------------
  ## Per-station and per-observation records.
  ##
  ## WHY BOTH, AND WHY THE STATION IS THE UNIT. With point-referenced data the
  ## unit that carries a regime is the STATION, so the error measure that
  ## answers "does the clustering help, and where" is the one computed within a
  ## station over time and then looked at across stations. Pooling over all n*T
  ## cells at once gives a single number that is, on a balanced panel, exactly
  ## the root mean of the per-station mean squared errors -- so it is a summary
  ## OF the per-station measure, not an alternative to it, and it destroys the
  ## distribution that shows the regimes at work. On an unbalanced panel it is
  ## not even that: it weights a station by how many periods it was observed.
  ## Both are recorded; the station file is what the tables are built from.
  ##
  ## The target is mu, the conditional mean given the latent path, not z. The
  ## nugget in z is irreducible, so scoring against z would add the same
  ## constant to every model and compress every comparison towards one. In the
  ## APPLICATION mu is not available and the honest measure is out-of-sample
  ## against z, with spatio-temporal blocking; that is a different script.
  ## ------------------------------------------------------------------------
  mh_true <- signal(fit_true)
  mh_pool <- signal(fit_pool)
  mh_sel  <- signal(fit_sel)
  g_true_hat <- if (is.null(map)) rep(NA_integer_, cell$n) else
    as.integer(map)[fit_true$group]
  g_sel_hat <- if (is.null(fit_sel)) rep(NA_integer_, cell$n) else fit_sel$group

  colstat <- function(M, f) if (is.null(M)) rep(NA_real_, cell$n) else
    apply(f(M), 2, mean, na.rm = TRUE)
  sq <- function(M) (M - dat$mu)^2
  ab <- function(M) abs(M - dat$mu)
  sqz <- function(M) (M - dat$z)^2

  station <- cbind(key[rep(1L, cell$n), ], data.frame(
    station  = seq_len(cell$n),
    lon = dat$coordinates[, 1], lat = dat$coordinates[, 2],
    g_true   = dat$labels,
    g_hat    = g_true_hat,
    g_hat_sel = g_sel_hat,
    correct  = as.integer(g_true_hat == dat$labels),
    rmse_true = sqrt(colstat(mh_true, sq)),
    mae_true  = colstat(mh_true, ab),
    rmse_pool = sqrt(colstat(mh_pool, sq)),
    mae_pool  = colstat(mh_pool, ab),
    rmse_z_true = sqrt(colstat(mh_true, sqz)),
    stringsAsFactors = FALSE))
  rownames(station) <- NULL

  obs <- NULL
  if (rep <= CFG$keep_obs[1]) {
    ii <- rep(seq_len(cell$n), each = cell$TN)
    tt <- rep(seq_len(cell$TN), times = cell$n)
    obs <- cbind(key[rep(1L, cell$n * cell$TN), ], data.frame(
      t = tt, station = ii,
      lon = dat$coordinates[ii, 1], lat = dat$coordinates[ii, 2],
      x = as.vector(dat$x), z = as.vector(dat$z), mu = as.vector(dat$mu),
      g_true = dat$labels[ii], g_hat = g_true_hat[ii], g_hat_sel = g_sel_hat[ii],
      mu_hat = if (is.null(mh_true)) NA_real_ else as.vector(mh_true),
      mu_hat_sel = if (is.null(mh_sel)) NA_real_ else as.vector(mh_sel),
      mu_hat_pool = if (is.null(mh_pool)) NA_real_ else as.vector(mh_pool),
      stringsAsFactors = FALSE))
    rownames(obs) <- NULL
  }

  list(summary = summ, params = params, station = station, obs = obs)
}

## The row a replication leaves when it fails: the key of the run, NA for every
## measure, and the reason. On a long run on another machine a silent NA is a day
## lost, so the reason is always recorded.
sim_failed <- function(cell, rep, msg) {
  k0 <- data.frame(
    cell = sprintf("n%d_T%d_K%d_w%.2f_%s_%s_knn%d_a%g_l%g",
                   cell$n, cell$TN, cell$K, cell$omega, cell$id,
                   substr(cell$balance, 1, 3), cell$knn,
                   CFG$alpha[1], CFG$lambda[1]),
    rep = rep, n = cell$n, TN = cell$TN, K = cell$K, omega = cell$omega,
    id = cell$id, balance = cell$balance, knn = cell$knn,
    alpha = CFG$alpha[1], lambda = CFG$lambda[1],
    stringsAsFactors = FALSE)
  list(summary = cbind(k0, data.frame(
         k_hat = NA_integer_, phi_hat = NA_real_, k_correct = NA_integer_,
         ari_sel = NA_real_, share_sel = NA_real_,
         ari_true = NA_real_, share_true = NA_real_,
         rmse_true = NA_real_, rmse_pooled = NA_real_, rmse_ratio = NA_real_,
         nconf = NA_integer_, nfail = NA_integer_, secs = NA_real_,
         error = msg, stringsAsFactors = FALSE)),
       params = NULL, station = NULL, obs = NULL)
}

sim_one_safe <- function(cell, rep) {
  out <- tryCatch(sim_one(cell, rep),
                  error = function(e) sim_failed(cell, rep, conditionMessage(e)))
  if (is.null(out$summary$error)) out$summary$error <- ""
  out
}

## ===========================================================================
## THE COVERAGE EXPERIMENT (mode = "coverage")
##
## The point estimator is one thing; the interval built around it is another,
## and the second does not follow from the first. This experiment asks whether
## the refit-with-clustering bootstrap of SCSTEM_Bootstrap() and
## SCSTEM_BootInference() is calibrated:
##
##   (a) COVERAGE. Does a nominal 95 per cent interval for a regime-specific
##       parameter contain the true value 95 per cent of the time? Reported for
##       the four interval families the package returns -- normal, basic,
##       percentile, bias-corrected -- at two nominal levels.
##   (b) STANDARD ERRORS. Is the average bootstrap standard error the same size
##       as the Monte Carlo standard deviation of the estimator it describes?
##       Their ratio separates an interval of the wrong width from one of the
##       wrong shape.
##   (c) THE PARTITION. How often does a bootstrap refit recover the partition
##       of the fit it was generated from? The bootstrap re-runs the clustering
##       on every draw, so its draws carry label uncertainty; if the refits
##       scatter, the intervals widen for a reason that is real.
##
## Each replication costs one fit plus boot_B refits of the WHOLE procedure, so
## this is the expensive experiment and runs on a deliberately small set of
## cells: the reference cell, at two overlaps and two scenarios, at the TRUE
## number of regimes and a fixed penalty. Coverage is asked of the estimator,
## not of the selection rule, and mixing the two would leave a failure
## unattributable.
##
## Two outputs, both keyed on (cell, rep):
##   <tag>-coverage.csv    one row per replication, level, regime, parameter
##   <tag>-stability.csv   one row per replication: how far the refits wander
## ===========================================================================
SIM_COVERAGE <- list(
  n        = SIM_REFERENCE$n,
  TN       = SIM_REFERENCE$TN,
  K        = SIM_REFERENCE$K,
  omega    = c(0, SIM_REFERENCE$omega),
  balance  = SIM_REFERENCE$balance,
  scenario = c("S1b", "S4"),
  phi      = 0.5,
  levels   = c(0.90, 0.95)
)

cov_one <- function(cell, rep) {
  row <- SIM_SCEN[SIM_SCEN$id == cell$id, , drop = FALSE]
  dat <- dgp_draw(cell$n, cell$TN, cell$K, cell$omega, row,
                  rep = CFG$seed0[1] + 5000L + rep, balance = cell$balance)
  mod <- sim_model(dat, cell$n)

  t0 <- proc.time()[["elapsed"]]
  fit <- Stem::SCSTEM_Estimation(mod, k = cell$K, phi_penalty = SIM_COVERAGE$phi,
                                 distance = "geo", verbose = FALSE,
                                 min_cluster_size = N_MIN,
                                 alpha = CFG$alpha[1], lambda = CFG$lambda[1])
  boot <- Stem::SCSTEM_Bootstrap(fit, B = CFG$boot_B[1],
                                 seed = CFG$seed0[1] + 5000L + rep,
                                 verbose = FALSE)
  secs <- proc.time()[["elapsed"]] - t0

  ## the bootstrap aligns every refit onto the ORIGINAL fit; aligning the
  ## original fit onto the TRUTH is this experiment's job
  map <- scstem_align_labels(reference = dat$labels, refit = fit$group,
                             K = cell$K)
  tru <- dgp_truth(dat$psi)
  cell_id <- sprintf("n%d_T%d_K%d_w%.2f_%s_%s_phi%g", cell$n, cell$TN,
                     cell$K, cell$omega, cell$id, substr(cell$balance, 1, 3),
                     SIM_COVERAGE$phi)
  key <- data.frame(cell = cell_id, rep = rep, n = cell$n, TN = cell$TN,
                    K = cell$K, omega = cell$omega, id = cell$id,
                    balance = cell$balance, stringsAsFactors = FALSE)

  out <- list(); inf <- NULL
  for (lev in SIM_COVERAGE$levels) {
    inf <- try(Stem::SCSTEM_BootInference(boot, level = lev, digits = 12),
               silent = TRUE)
    if (inherits(inf, "try-error")) { inf <- NULL; next }
    s <- inf$summary
    for (i in seq_len(nrow(tru))) {
      g_fit <- which(map == tru$regime[i])
      p <- tru$parameter[i]
      j <- if (length(g_fit)) which(s$cluster == g_fit[1] & s$parameter == p) else integer(0)
      if (!length(j)) next
      r <- s[j[1], ]
      out[[length(out) + 1L]] <- cbind(key, data.frame(
        level = lev, regime = tru$regime[i], parameter = p,
        truth = tru$truth[i], estimate = r$estimate, se = r$se,
        normal_lo = r$normal_lo, normal_up = r$normal_up,
        basic_lo = r$basic_lo, basic_up = r$basic_up,
        perc_lo = r$perc_lo, perc_up = r$perc_up,
        bc_lo = r$bc_lo, bc_up = r$bc_up,
        n_draws = r$n_draws, stringsAsFactors = FALSE))
    }
  }
  stab <- cbind(key, data.frame(
    B_used = if (is.null(boot$B_used)) CFG$boot_B[1] else boot$B_used,
    ari_mean = if (is.null(inf)) NA_real_ else mean(inf$stability$ARI, na.rm = TRUE),
    ari_min = if (is.null(inf)) NA_real_ else min(inf$stability$ARI, na.rm = TRUE),
    ari_share1 = if (is.null(inf)) NA_real_ else mean(inf$stability$ARI >= 0.999, na.rm = TRUE),
    secs = secs, stringsAsFactors = FALSE))
  list(coverage = if (length(out)) do.call(rbind, out) else NULL,
       stability = stab)
}

cov_one_safe <- function(cell, rep) {
  tryCatch(cov_one(cell, rep), error = function(e)
    list(coverage = NULL, stability = NULL, error = conditionMessage(e)))
}


## ===========================================================================
## Running the tasks on several cores
## ===========================================================================

## What a core does with one task. Top-level functions, so that only their
## name travels to the cores, not the environment they were created in.
sim_task_fun <- function(task) sim_one_safe(task$cell, task$rep)
cov_task_fun <- function(task) cov_one_safe(task$cell, task$rep)

## What a core does once, when it starts: load this very file in
## definitions-only mode -- Stem, the generator, the functions above -- and
## take the settings of the session that launched it.
sim_worker_init <- function(path, cfg) {
  assign("SIM_DEFINE_ONLY", TRUE, envir = globalenv())
  source(path, local = globalenv())
  assign("CFG", cfg, envir = globalenv())
  invisible(TRUE)
}

## Runs `fun` on every task and hands each result to `on_result` in this
## session as soon as it arrives. With one core the tasks simply run here; with
## more, a cluster of R processes is started and the tasks are handed out one at
## a time as the processes free up -- the scheduling of
## parallel::clusterApplyLB(), written out so that each result can be saved the
## moment it arrives instead of all at the end. The cluster is stopped on exit,
## also when the run is interrupted from RStudio.
sim_run_tasks <- function(tasks, fun, on_result, cores) {
  n <- length(tasks)
  if (!n) return(invisible(0L))
  cores <- max(1L, min(as.integer(cores), n))
  if (cores > 1L && is.null(SIM_FILE)) {
    message("the path of this script is unknown, so it runs on one core; ",
            "open the file and press Source rather than pasting it")
    cores <- 1L
  }
  if (cores == 1L) {
    for (k in seq_len(n)) on_result(fun(tasks[[k]]), tasks[[k]], k, n)
    return(invisible(n))
  }
  cat(sprintf("starting %d R processes ...\n", cores)); utils::flush.console()
  cl <- parallel::makePSOCKcluster(cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, sim_worker_init, SIM_FILE, CFG)
  send <- utils::getFromNamespace("sendCall", "parallel")
  recv <- utils::getFromNamespace("recvOneResult", "parallel")
  for (i in seq_len(cores)) send(cl[[i]], fun, list(tasks[[i]]), tag = i)
  nxt <- cores + 1L
  for (k in seq_len(n)) {
    d <- recv(cl)
    if (nxt <= n) {
      send(cl[[d$node]], fun, list(tasks[[nxt]]), tag = nxt)
      nxt <- nxt + 1L
    }
    on_result(d$value, tasks[[d$tag]], k, n)
  }
  invisible(n)
}

## The tasks in the order they are run: replication by replication, so that an
## interrupted run leaves whole replications of the design behind, and within a
## replication in a shuffled but fixed order of the cells, so that the running
## average of the times -- hence the estimate of the time left -- is
## representative from the first minutes rather than biased by cheap cells.
sim_tasks <- function(cells, reps, done_keys, key_of) {
  set.seed(20260926)
  ord <- sample(nrow(cells))
  tasks <- list()
  for (r in reps) for (i in ord) {
    cell <- cells[i, , drop = FALSE]
    if (key_of(cell, r) %in% done_keys) next
    tasks[[length(tasks) + 1L]] <- list(cell = cell, rep = r)
  }
  tasks
}

sim_hours <- function(x) if (x < 1) sprintf("%.0f min", 60 * x) else sprintf("%.1f h", x)


## ===========================================================================
## The Monte Carlo
## ===========================================================================
sim_main <- function() {

  OUT <- normalizePath(CFG$out[1], winslash = "/", mustWork = FALSE)
  sim_print_config(CFG)

  CSV     <- file.path(OUT, sprintf("%s.csv", CFG$tag))
  CSV_PAR <- file.path(OUT, sprintf("%s-params.csv", CFG$tag))
  CSV_STA <- file.path(OUT, sprintf("%s-stations.csv", CFG$tag))
  DIR_OBS <- file.path(OUT, sprintf("%s-obs", CFG$tag))

  ## ---------------------------------------------------------------------------
  ## The cells
  ##
  ## K = 1 has a single regime, so the overlap, the balance and the scenario are
  ## all vacuous there: it enters once per (n, T). It covers the same area as
  ## every other cell, because the generator holds the total spatial variance
  ## fixed and a single regime simply takes all of it -- see NU_TOT above --
  ## so the selection rule is not handed a free geometric cue for telling k = 1
  ## from k > 1. Cells whose imbalance cannot be realised with regimes of at least
  ## N_MIN units are dropped rather than silently rebalanced.
  ## ---------------------------------------------------------------------------
  cells <- sim_cells(blocks = CFG$blocks)
  ## only_<factor> keeps the named levels of that factor and nothing else.
  ## The option name is on the left, the column of `cells` it filters on the right.
  for (opt in list(c("only_n", "n"), c("only_TN", "TN"), c("only_K", "K"),
                   c("only_knn", "knn"), c("only_scenario", "id"),
                   c("only_balance", "balance"))) {
    lev <- CFG[[opt[1]]]
    if (length(lev)) cells <- cells[cells[[opt[2]]] %in% lev, ]
  }
  ## omega is a double, so it is matched with a tolerance rather than with %in%
  if (length(CFG$only_omega))
    cells <- cells[vapply(cells$omega, function(w)
      any(abs(w - CFG$only_omega) < 1e-8), logical(1)), ]
  cells <- cells[mapply(dgp_feasible, cells$n, cells$K, cells$balance), ]
  rownames(cells) <- NULL

  prev <- if (file.exists(CSV)) utils::read.csv(CSV, stringsAsFactors = FALSE) else NULL
  budget <- sim_cost(cells, nrep = CFG$nrep, cores = CFG$cores[1], measured = prev)
  cat(sprintf("%d cells x %d replications, writing to\n  %s\n", nrow(cells),
              CFG$nrep, OUT))
  if (is.na(budget$calibration)) {
    cat(sprintf("estimated %.0f core-hours (%s on %d cores), from the prior table:\n",
                budget$core_hours, sim_hours(budget$wall_hours), CFG$cores[1]))
    cat("  a rough guide only; the estimate printed as the run goes is the one to trust\n\n")
  } else {
    cat(sprintf("estimated %.0f core-hours (%s on %d cores), calibrated on the %d\n",
                budget$core_hours, sim_hours(budget$wall_hours), CFG$cores[1],
                budget$measured_reps))
    cat("  replications already recorded\n\n")
  }

  ## "dry" prices the run and stops: nothing is fitted and nothing is written
  if (identical(CFG$mode[1], "dry") || sim_flag("dry")) {
    cat("cells by number of locations and periods\n")
    print(table(n = cells$n, T = cells$TN))
    cat("\ncells by block\n")
    for (b in CFG$blocks)
      cat(sprintf("  %-10s %4d\n", b, nrow(sim_cells(blocks = b))))
    return(invisible(cells))
  }

  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  dir.create(DIR_OBS, recursive = TRUE, showWarnings = FALSE)

  ## the resume test is on the primary key of a run
  cell_key <- function(cl, r) {
    sprintf("n%d_T%d_K%d_w%.2f_%s_%s_knn%d_a%g_l%g|%d",
            cl$n, cl$TN, cl$K, cl$omega, cl$id, substr(cl$balance, 1, 3), cl$knn,
            CFG$alpha[1], CFG$lambda[1], r)
  }
  done <- if (is.null(prev)) character(0) else paste(prev$cell, prev$rep, sep = "|")

  ## The per-observation record goes to one compressed file per cell and
  ## replication rather than into a shared CSV: it is n*T rows, so appending it to
  ## a single file would produce something no editor opens and no resume could
  ## check cheaply.
  save_obs <- function(df, cell, rep) {
    if (is.null(df)) return(invisible(NULL))
    f <- sprintf("obs_n%d_T%d_K%d_w%s_%s_%s_knn%d_rep%04d.rds",
                 cell$n, cell$TN, cell$K, sub(".", "", sprintf("%.2f", cell$omega), fixed = TRUE),
                 cell$id, substr(cell$balance, 1, 3), cell$knn, rep)
    sim_retry(function() saveRDS(df, file.path(DIR_OBS, f), compress = "xz"))
  }

  reps  <- seq.int(CFG$rep_from[1], CFG$rep_to[1])
  tasks <- sim_tasks(cells, reps, done, cell_key)
  cat(sprintf("replications %d to %d: %d tasks to run, %d already recorded;",
              CFG$rep_from[1], CFG$rep_to[1], length(tasks),
              nrow(cells) * length(reps) - length(tasks)))
  cat(sprintf(" per-observation records kept for the first %d\n\n", CFG$keep_obs[1]))
  if (!length(tasks)) { cat("nothing to do\n"); return(invisible(NULL)) }

  cat(sprintf("%13s %4s %5s %2s %5s %-4s %-3s %3s %4s | %9s %5s %6s | %s\n",
              "", "n", "T", "K", "omega", "scen", "bal", "knn", "rep",
              "time", "k_hat", "ARI", "elapsed, and left"))
  t_start <- proc.time()[["elapsed"]]
  on_result <- function(out, task, k, n) {
    ## a result that is not a list is an R process that died on the task
    if (!is.list(out) || is.null(out$summary))
      out <- sim_failed(task$cell, task$rep, paste(as.character(out), collapse = " "))
    if (is.null(out$summary$error)) out$summary$error <- ""
    sim_append(out$summary, CSV)
    sim_append(out$params,  CSV_PAR)
    sim_append(out$station, CSV_STA)
    save_obs(out$obs, task$cell, task$rep)
    s  <- out$summary
    el <- (proc.time()[["elapsed"]] - t_start) / 3600
    cat(sprintf("[%5d/%5d] %4d %5d %2d %5.2f %-4s %-3s %3d %4d | %7.1f s %5s %6s | %s, ~%s left\n",
                k, n, s$n, s$TN, s$K, s$omega, s$id, substr(s$balance, 1, 3), s$knn,
                s$rep, s$secs, format(s$k_hat), format(round(s$ari_true, 3)),
                sim_hours(el), sim_hours(el / k * (n - k))))
    if (nzchar(s$error)) cat("              FAILED: ", s$error, "\n", sep = "")
    utils::flush.console()
  }

  sim_run_tasks(tasks, sim_task_fun, on_result, CFG$cores[1])
  cat(sprintf("\ndone in %s\n", sim_hours((proc.time()[["elapsed"]] - t_start) / 3600)))
  invisible(NULL)
}


## ===========================================================================
## The coverage experiment
## ===========================================================================
sim_coverage <- function() {

  OUT <- normalizePath(CFG$out[1], winslash = "/", mustWork = FALSE)
  cv  <- SIM_COVERAGE
  CSV <- file.path(OUT, sprintf("%s-coverage.csv", CFG$tag))
  CSA <- file.path(OUT, sprintf("%s-stability.csv", CFG$tag))

  cells <- expand.grid(n = cv$n, TN = cv$TN, K = cv$K, omega = cv$omega,
                       id = cv$scenario, balance = cv$balance,
                       stringsAsFactors = FALSE)
  cells <- cells[mapply(dgp_feasible, cells$n, cells$K, cells$balance), ]
  rownames(cells) <- NULL
  cat(sprintf("coverage: %d cells x %d replications x %d refits, writing to\n  %s\n\n",
              nrow(cells), CFG$boot_reps, CFG$boot_B, OUT))
  if (identical(CFG$mode[1], "dry") || sim_flag("dry")) return(invisible(cells))
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

  cov_key <- function(cell, r)
    paste(sprintf("n%d_T%d_K%d_w%.2f_%s_%s_phi%g", cell$n, cell$TN, cell$K,
                  cell$omega, cell$id, substr(cell$balance, 1, 3), cv$phi), r)
  prev <- if (file.exists(CSA)) utils::read.csv(CSA, stringsAsFactors = FALSE) else NULL
  done <- if (is.null(prev)) character(0) else paste(prev$cell, prev$rep)
  tasks <- sim_tasks(cells, seq_len(CFG$boot_reps[1]), done, cov_key)
  cat(sprintf("%d tasks to run\n\n", length(tasks)))
  if (!length(tasks)) { cat("nothing to do\n"); return(invisible(NULL)) }

  t_start <- proc.time()[["elapsed"]]
  on_result <- function(out, task, k, n) {
    el <- (proc.time()[["elapsed"]] - t_start) / 3600
    cell <- task$cell
    if (!is.list(out) || !is.null(out$error) || is.null(out$stability)) {
      msg <- if (is.list(out)) out$error else paste(as.character(out), collapse = " ")
      cat(sprintf("[%4d/%4d] %4d %5d %2d %5.2f %-4s %4d | FAILED: %s\n", k, n, cell$n,
                  cell$TN, cell$K, cell$omega, cell$id, task$rep, msg))
      return(invisible(NULL))
    }
    sim_append(out$coverage, CSV)
    sim_append(out$stability, CSA)
    s <- out$stability
    cat(sprintf("[%4d/%4d] %4d %5d %2d %5.2f %-4s %4d | %8.1f s  ARI %.3f | %s, ~%s left\n",
                k, n, s$n, s$TN, s$K, s$omega, s$id, s$rep, s$secs, s$ari_mean,
                sim_hours(el), sim_hours(el / k * (n - k))))
    utils::flush.console()
  }
  sim_run_tasks(tasks, cov_task_fun, on_result, CFG$cores[1])
  invisible(NULL)
}


## ===========================================================================
## Run
## ===========================================================================
if (!SIM_DEFINE_ONLY) {
  if (identical(CFG$mode[1], "coverage") || sim_flag("coverage")) sim_coverage() else sim_main()
}
