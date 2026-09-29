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
##        - "1. THE MARGINS": the levels of every factor. The design is the full
##          factorial of these lists, so adding a value to TN, say, adds that
##          series length to every scenario;
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
## Do not run the same tag on two machines at once: they would write to the same
## files. To share the work between machines give each its own `tag` and its own
## range of replications.
##
## Every setting of the SETUP can also be given from a terminal, which overrides
## it; lists are separated by commas:
##
##     Rscript run-simulations.R --cores=8 --rep_from=51 --rep_to=100
##     Rscript run-simulations.R --TN=60,120,365 --mode=dry
##
## This is part of the REPLICATION MATERIAL of the paper, not of the Stem
## package, and nothing here is shipped with it.
##
## FIVE OUTPUTS, keyed on (cell, rep) and joining on it and on nothing else:
##
##   <tag>.csv           one row per replication: the selected (k, phi), the
##                       recovery of the partition, the error of the fitted
##                       signal, the times
##   <tag>-params.csv    one row per replication, regime and parameter: the
##                       truth beside the estimate at the true number of
##                       regimes, after aligning the estimated regimes on the
##                       true ones
##   <tag>-stations.csv  one row per replication and location
##   <tag>-grid.csv      one row per replication and (k, phi) of the fitted
##                       grid: log-likelihood, parameters, criteria, and the
##                       ARI of its partition with the truth
##   <tag>-obs/*.rds     the full observations of the first keep_obs
##                       replications of every cell
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

## Pinned to the commit of 2026-09-28 (departures initialization, warm start,
## two-step rule with phi by the smallest BIC at the selected k)
SIM_STEM_REF <- "PaoloMaranzano/Stem@4325536872218e27a39f9fd21ba533196c9ace50"

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
  have <- c("STEM_Signal", "SCSTEM_Signal", "SCSTEM_Infocrit", "SCSTEM_Select",
            "scstem_neighbors", "scstem_align_labels", "scstem_ari")
  all(vapply(have, exists, logical(1), envir = ns, inherits = FALSE)) &&
    "distance" %in% names(formals(get("scstem_neighbors", envir = ns)))
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
## 1. THE MARGINS. The design is the full factorial of these lists. Remove a
##    value to drop it, add one to cross it with everything else.
## ---------------------------------------------------------------------------
SIM_MARGINS <- list(

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
  scenario = c("S2", "S1s-ind", "S1s-shr", "S1w-shr",
               "S3beta-shr", "S3G-shr", "S3Seta-shr", "S3theta-shr", "S3error-shr"),

  ## number of locations; not below 40, since the grid reaches k = 4 regimes
  ## of at least N_MIN = 6 locations each
  n = c(40L, 100L, 200L, 400L),

  ## length of the series
  TN = c(60L, 120L, 365L),

  ## Spatial overlap of the three regimes: the centres sit at the vertices of
  ## an equilateral triangle of side 2*omega, with the total variance of a
  ## coordinate held fixed. 0: the regimes coincide in space; 0.70: contiguous
  ## areas interpenetrating along their borders (1.63 within-regime standard
  ## deviations between centres); 1: well apart. Vacuous in S2.
  omega = c(0, 0.70, 1),

  ## relative sizes of the regimes: "balanced" or "unbalanced" (1:2:3).
  ## Vacuous in S2.
  balance = c("balanced", "unbalanced"),

  ## neighbours of the graph of the Potts penalty
  knn = c(3L, 5L, 10L)
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
  strong = list(b0 = c(2.00, 2.30, 2.60), b1 = c(0.55, 0.85, 1.15),
                G  = c(0.80, 0.50, 0.20), v  = c(0.30, 0.40, 0.50),
                se = c(0.20, 0.35, 0.50), so = c(0.20, 0.15, 0.10),
                R  = c(2, 1, 0.5))
)

## S3 keeps the values of S1 at this level for the blocks that differ
SIM_S3_LEVEL <- "strong"

## ---------------------------------------------------------------------------
## 3. THE RUN. What this invocation does.
## ---------------------------------------------------------------------------
CFG <- sim_config(args = if (SIM_DEFINE_ONLY) character(0) else commandArgs(TRUE), c(
  SIM_MARGINS,
  list(
    ## "run": the Monte Carlo; "dry": print the cells and the cost, run nothing
    mode     = "run",

    ## How many cores. Each core runs one replication at a time.
    cores    = sim_default_cores(),

    ## replications per cell, and the slice this invocation covers. A long
    ## study is executed in slices on whatever machine is free: the files are
    ## appended, so the slices compose.
    nrep     = 10L,
    rep_from = 1L,
    rep_to   = 10L,

    ## What the estimator searches over. phi_ref is the penalty at which the
    ## recovery at the true number of regimes is read, a point of phi_grid.
    k_grid   = 1:4,
    phi_grid = c(0, 0.025, 0.05, 0.1, 0.2, 0.5, 1),
    phi_ref  = 0.05,

    ## How many replications per cell keep their full per-observation record.
    ## One: across the full factorial every record is n x T rows, and five per
    ## cell would take several gigabytes; any replication can be regenerated
    ## from its seed.
    keep_obs = 1L,

    ## Bookkeeping. `out` defaults to a results/ folder beside this script.
    tag      = "main",
    out      = file.path(SIM_HERE, "results")
  )))


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
    K = matrix(1, n, 1))
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

dgp_locations <- function(n, K, omega, balance = "balanced", seed = 1) {
  set.seed(seed)
  g  <- rep(seq_len(K), times = dgp_sizes(n, K, balance))
  mu <- dgp_centres(K, omega)
  nu_sp <- NU_TOT - dgp_centre_var(K, omega)
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
## One complete data set. Seeds are derived from the replication index, so
## replication r uses the same locations and covariate in every
## scenario-variant (S2, with its single cloud, apart), and the same random
## numbers for the latent and the error, so that two variants differ by their
## definition and not by their draw.
## ---------------------------------------------------------------------------
dgp_draw <- function(cell, rep) {
  row  <- SIM_SCEN[SIM_SCEN$id == cell$scenario, , drop = FALSE]
  K    <- if (row$family == "S2") 1L else 3L
  seed <- 1000L * rep + 1L
  loc  <- dgp_locations(cell$n, K, if (K == 1L) 0 else cell$omega,
                        balance = if (K == 1L) "balanced" else cell$balance,
                        seed = seed)
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
## The cells: the full factorial of the margins. In S2 the overlap and the
## balance are vacuous, so it enters once per (n, T, knn).
## ===========================================================================
sim_cells <- function(cfg = CFG) {
  bad <- setdiff(cfg$scenario, SIM_SCEN$id)
  if (length(bad)) stop("unknown scenario: ", paste(bad, collapse = ", "),
                        "; the scenarios are ", paste(SIM_SCEN$id, collapse = ", "),
                        call. = FALSE)
  cells <- expand.grid(scenario = cfg$scenario, n = cfg$n, TN = cfg$TN,
                       omega = cfg$omega, balance = cfg$balance, knn = cfg$knn,
                       stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  s2 <- cells$scenario == "S2"
  cells$omega[s2] <- NA_real_
  cells$balance[s2] <- "-"
  cells <- unique(cells)
  ok <- mapply(function(sc, n, bal) if (sc == "S2") TRUE else dgp_feasible(n, 3L, bal),
               cells$scenario, cells$n, cells$balance)
  cells <- cells[ok, , drop = FALSE]
  cells$cell <- sim_cell_id(cells)
  rownames(cells) <- NULL
  cells
}

## The primary key of a cell: a string, because omega is a double and a join on
## a floating-point column is a defect waiting to happen.
sim_cell_id <- function(cells) {
  sprintf("%s_n%d_T%d_w%s_%s_knn%d", cells$scenario, as.integer(cells$n),
          as.integer(cells$TN),
          ifelse(is.na(cells$omega), "NA", sprintf("%.2f", cells$omega)),
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
sim_one <- function(cell, rep) {

  t_all <- proc.time()[["elapsed"]]
  dat <- dgp_draw(cell, rep)
  mod <- sim_model(dat)

  t0 <- proc.time()[["elapsed"]]
  ic <- Stem::SCSTEM_Infocrit(mod, k_grid = CFG$k_grid, phi_grid = CFG$phi_grid,
                              distance = "euclidean", verbose = FALSE,
                              knn = cell$knn, min_cluster_size = N_MIN)
  sel <- Stem::SCSTEM_Select(ic)
  secs <- proc.time()[["elapsed"]] - t0

  ## SELECTION is read off the selected fit; RECOVERY off the fit at the TRUE
  ## number of regimes, so that a poor estimate is not confounded with a poor
  ## selection. Both, and the pooled benchmark, are in the grid.
  grab <- function(kk, pp) {
    j <- which(ic$table$k == kk & abs(ic$table$phi - pp) < 1e-8)
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
  key <- data.frame(cell = cell$cell, rep = rep, scenario = cell$scenario,
                    family = row$family, level = row$level, common = row$common,
                    variant = row$variant, rho = row$rho, field = row$field,
                    n = cell$n, TN = cell$TN, omega = cell$omega,
                    balance = cell$balance, knn = cell$knn, K_true = K,
                    stringsAsFactors = FALSE)

  summ <- cbind(key, data.frame(
    k_hat = sel$k_selected, phi_hat = sel$phi_selected,
    k_correct = as.integer(sel$k_selected == K),
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

  ## the fitted grid, one row per (k, phi), with the ARI of every partition, so
  ## that the recovery at any penalty, or under another selection rule, can be
  ## read afterwards without refitting
  grid <- cbind(key[rep(1L, nrow(ic$table)), ], ic$table)
  tags <- paste0("k=", ic$table$k, ", phi=", ic$table$phi)
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
  if (rep <= CFG$keep_obs[1]) {
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
sim_failed <- function(cell, rep, msg) {
  row <- SIM_SCEN[SIM_SCEN$id == cell$scenario, , drop = FALSE]
  list(summary = data.frame(
    cell = cell$cell, rep = rep, scenario = cell$scenario,
    family = row$family, level = row$level, common = row$common,
    variant = row$variant, rho = row$rho, field = row$field,
    n = cell$n, TN = cell$TN, omega = cell$omega, balance = cell$balance,
    knn = cell$knn, K_true = if (row$family == "S2") 1L else 3L,
    k_hat = NA_integer_, phi_hat = NA_real_, k_correct = NA_integer_,
    ari_sel = NA_real_, share_sel = NA_real_, ari_true = NA_real_, share_true = NA_real_,
    rmse_sel = NA_real_, rmse_true = NA_real_, rmse_pooled = NA_real_,
    nconf = NA_integer_, nfail = NA_integer_, secs = NA_real_, secs_total = NA_real_,
    error = msg, stringsAsFactors = FALSE),
    params = NULL, grid = NULL, station = NULL, obs = NULL)
}

sim_one_safe <- function(cell, rep) {
  tryCatch(sim_one(cell, rep),
           error = function(e) sim_failed(cell, rep, conditionMessage(e)))
}


## ===========================================================================
## Running the tasks on several cores
## ===========================================================================
sim_task_fun <- function(task) sim_one_safe(task$cell, task$rep)

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

  CSV     <- file.path(OUT, sprintf("%s.csv", CFG$tag))
  CSV_PAR <- file.path(OUT, sprintf("%s-params.csv", CFG$tag))
  CSV_GRD <- file.path(OUT, sprintf("%s-grid.csv", CFG$tag))
  CSV_STA <- file.path(OUT, sprintf("%s-stations.csv", CFG$tag))
  DIR_OBS <- file.path(OUT, sprintf("%s-obs", CFG$tag))

  cells <- sim_cells()
  prev  <- if (file.exists(CSV)) utils::read.csv(CSV, stringsAsFactors = FALSE) else NULL
  cat(sprintf("%d cells x %d replications, writing to\n  %s\n\n", nrow(cells),
              CFG$nrep, OUT))

  if (identical(CFG$mode[1], "dry")) {
    print(cells[, c("scenario", "n", "TN", "omega", "balance", "knn")], row.names = FALSE)
    cat("\n")
    sim_print_cost(cells, prev, CFG)
    return(invisible(cells))
  }

  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  dir.create(DIR_OBS, recursive = TRUE, showWarnings = FALSE)
  done <- if (is.null(prev)) character(0) else paste(prev$cell, prev$rep, sep = "|")

  save_obs <- function(df, cell, rep) {
    if (is.null(df)) return(invisible(NULL))
    f <- sprintf("obs_%s_rep%04d.rds", gsub("[^A-Za-z0-9_.-]", "", cell$cell), rep)
    sim_retry(function() saveRDS(df, file.path(DIR_OBS, f), compress = "xz"))
  }

  reps  <- seq.int(CFG$rep_from[1], CFG$rep_to[1])
  tasks <- sim_tasks(cells, reps, done)
  cat(sprintf("replications %d to %d: %d tasks to run, %d already recorded\n\n",
              CFG$rep_from[1], CFG$rep_to[1], length(tasks),
              nrow(cells) * length(reps) - length(tasks)))
  if (!length(tasks)) { cat("nothing to do\n"); return(invisible(NULL)) }

  cat(sprintf("%13s %-12s %4s %5s %4s | %8s %5s %6s | %s\n",
              "", "scenario", "n", "T", "rep", "time", "k_hat", "ARI", "elapsed, and left"))
  t_start <- proc.time()[["elapsed"]]
  on_result <- function(out, task, k, n) {
    if (!is.list(out) || is.null(out$summary))
      out <- sim_failed(task$cell, task$rep, paste(as.character(out), collapse = " "))
    sim_append(out$summary, CSV)
    sim_append(out$params,  CSV_PAR)
    sim_append(out$grid,    CSV_GRD)
    sim_append(out$station, CSV_STA)
    save_obs(out$obs, task$cell, task$rep)
    s  <- out$summary
    el <- (proc.time()[["elapsed"]] - t_start) / 3600
    cat(sprintf("[%5d/%5d] %-12s %4d %5d %4d | %6.1f s %5s %6s | %s, ~%s left\n",
                k, n, s$scenario, s$n, s$TN, s$rep, s$secs_total, format(s$k_hat),
                format(round(s$ari_true, 3)), sim_hours(el), sim_hours(el / k * (n - k))))
    if (nzchar(s$error)) cat("              FAILED: ", s$error, "\n", sep = "")
    utils::flush.console()
  }

  sim_run_tasks(tasks, sim_task_fun, on_result, CFG$cores[1])
  cat(sprintf("\ndone in %s\n\n", sim_hours((proc.time()[["elapsed"]] - t_start) / 3600)))
  sim_print_cost(cells, utils::read.csv(CSV, stringsAsFactors = FALSE), CFG)
  invisible(NULL)
}

if (!SIM_DEFINE_ONLY) sim_main()
