## ===========================================================================
## SC-STEM: the application. One script, self-contained.
##
## It needs nothing but an installed Stem package, it runs from whatever folder
## it sits in, and it writes its results beside itself, in application/<commit>/. The
## script is staged and every stage caches its result, so an interrupted run
## resumes where it stopped: delete a cache file to force that stage to run
## again.
##
##     Rscript run-application.R
##     Rscript run-application.R --no-bootstrap     stop before the expensive stage
##     Rscript run-application.R --data=fuels.rds   another data set
##     Rscript run-application.R --lambda=0.1       with a ridge on the coefficients
##
## Stages
##   A  the data, the model object and its starting values
##   B  the pooled STEM fit, which is the k = 1 reference
##   C  the (k, phi) grid, through SCSTEM_Infocrit()
##   D  the two-step selection rule
##   E  the refit-with-clustering bootstrap at the selected configuration
##   F  the outputs
##
## THE DATA. By default the Po Valley data shipped with Stem. Any other data set
## is passed with --data=<file>.rds, holding a list with three elements in the
## shape STEM_Model() expects:
##
##   z           T x n matrix of the response, NA where missing
##   covariates  (n*T) x p matrix, stacked by location, intercept included
##   coords      n x 2 matrix of longitude and latitude, in that order
##
## This is part of the REPLICATION MATERIAL of the paper, not of the Stem
## package. It shares no code with the simulation scripts on purpose: each
## script of the replication material runs on its own.
## ===========================================================================


## ---------------------------------------------------------------------------
## Where this script is: from Rscript, or from Source in RStudio
## ---------------------------------------------------------------------------
app_here <- local({
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  f <- if (length(m)) normalizePath(sub("^--file=", "", m[1]), winslash = "/") else NULL
  if (is.null(f)) for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) { f <- normalizePath(of, winslash = "/"); break }
  }
  if (is.null(f)) normalizePath(getwd(), winslash = "/") else dirname(f)
})

## ---------------------------------------------------------------------------
## The only prerequisite: Stem, installed from GitHub at the commit of the
## simulation study (main3: the four estimation algorithms with SQUAREM as the
## default, tolerances 1e-3, regularization 0, the range within the limits the
## distances identify), and again whenever the installed copy is another
## commit. Offline, an installed Stem that carries what is used here is
## accepted as it is. The development builds all say 2.0.0, so the check is on
## the features and on the commit.
## ---------------------------------------------------------------------------
APP_STEM_REF <- "PaoloMaranzano/Stem@0f7b74479b03a7fe5508eafe68fc7510c21422ad"

app_github_sha <- function(ref = APP_STEM_REF) {
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

app_stem_ok <- function(latest = NA_character_) {
  if (!requireNamespace("Stem", quietly = TRUE)) return(FALSE)
  ns <- asNamespace("Stem")
  have <- c("STEM_Model", "STEM_Estimation", "STEM_Signal", "STEM_control",
            "SCSTEM_Infocrit", "SCSTEM_Select", "SCSTEM_Bootstrap",
            "SCSTEM_BootInference", "scstem_neighbors")
  feats <- all(vapply(have, exists, logical(1), envir = ns, inherits = FALSE)) &&
    "distance" %in% names(formals(get("scstem_neighbors", envir = ns))) &&
    "algorithm" %in% names(formals(get("STEM_control", envir = ns)))
  here <- utils::packageDescription("Stem")$RemoteSha
  feats && (is.na(latest) || (!is.null(here) && identical(here, latest)))
}

APP_STEM_SHA <- app_github_sha()
if (!app_stem_ok(APP_STEM_SHA)) {
  message(if (requireNamespace("Stem", quietly = TRUE))
            "a newer Stem is on GitHub: updating it" else
            "Stem is not installed: installing it from GitHub")
  if (!requireNamespace("remotes", quietly = TRUE))
    utils::install.packages("remotes", repos = "https://cloud.r-project.org")
  if ("Stem" %in% loadedNamespaces()) try(unloadNamespace("Stem"), silent = TRUE)
  remotes::install_github(APP_STEM_REF, upgrade = "never", force = TRUE, quiet = TRUE)
  if (!app_stem_ok(APP_STEM_SHA))
    stop("Stem could not be installed or brought up to date. Restart R ",
         "(Session > Restart R in RStudio) and run again; if it still fails, ",
         "install it by hand with\n  remotes::install_github(\"", APP_STEM_REF,
         "\")", call. = FALSE)
}
suppressPackageStartupMessages(library("Stem"))

## ---------------------------------------------------------------------------
## Settings, all overridable with --name=value
## ---------------------------------------------------------------------------
app_config <- function(defaults, args = commandArgs(trailingOnly = TRUE)) {
  for (a in grep("^--[^=]+=", args, value = TRUE)) {
    nm  <- sub("^--([^=]+)=.*$", "\\1", a)
    val <- sub("^--[^=]+=", "", a)
    if (!nzchar(val)) next
    if (!nm %in% names(defaults))
      stop("unknown option --", nm, "; the options are: ",
           paste(names(defaults), collapse = ", "), call. = FALSE)
    parts <- trimws(strsplit(val, ",", fixed = TRUE)[[1]])
    d <- defaults[[nm]]
    defaults[[nm]] <- if (is.integer(d)) as.integer(parts)
      else if (is.numeric(d)) as.numeric(parts) else parts
  }
  defaults
}

CFG <- app_config(list(
  data     = "",                                  # empty: Stem's povalley
  K_grid   = 1:4,
  phi_grid = c(0, 0.025, 0.05, 0.1, 0.2, 0.5, 1),
  band     = c(0.025, 0.2),  # the band in which the tuning rule chooses k
  knn      = 5L,
  alpha    = 0,              # 0 is ridge; with lambda = 0 no penalty at all
  lambda   = 0,
  B        = 200L,
  seed     = 20260904L,
  ## one folder per pinned commit, so that a cache of an earlier Stem is never reused
  out      = file.path(app_here, "application", substr(sub("^.*@", "", APP_STEM_REF), 1L, 7L))
))
OUT <- CFG$out[1]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
OUT <- normalizePath(OUT, winslash = "/")
message("application: results in ", OUT)

cached <- function(name, expr) {
  f <- file.path(OUT, paste0(name, ".rds"))
  if (file.exists(f)) {
    message("  [cache] ", name)
    return(readRDS(f))
  }
  message("  [run  ] ", name)
  t0 <- Sys.time()
  val <- force(expr)
  message("  [done ] ", name, " in ",
          round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
  saveRDS(val, f)
  val
}

## ---- A: data and starting values -------------------------------------------
message("A. data and model object")
if (nzchar(CFG$data[1])) {
  src <- CFG$data[1]
  if (!file.exists(src)) src <- file.path(app_here, src)
  dat <- readRDS(src)
  miss <- setdiff(c("z", "covariates", "coords"), names(dat))
  if (length(miss))
    stop("the data file lacks ", paste(miss, collapse = ", "), call. = FALSE)
  X <- dat$covariates; z <- dat$z; coords <- dat$coords
} else {
  povalley <- NULL
  utils::data("povalley", package = "Stem", envir = environment())
  X <- povalley$covariates; z <- povalley$z; coords <- povalley$coords
}
d <- ncol(z); Tn <- nrow(z)
message(sprintf("  %d locations, %d periods, %d covariates, %.1f%% missing",
                d, Tn, ncol(as.matrix(X)), 100 * mean(is.na(z))))

## Starting values from the pooled regression. The residual variance is split
## between a nugget and a spatially structured part; the range starts at 60 km
## and the latent process starts persistent.
ok  <- is.finite(as.vector(z))
ols <- stats::lm.fit(x = as.matrix(X)[ok, , drop = FALSE], y = as.vector(z)[ok])
s2  <- stats::var(ols$residuals)
phi0 <- list(beta = matrix(ols$coefficients, ncol = 1),
             sigma2eps = 0.5 * s2, sigma2omega = 0.4 * s2,
             theta = 1 / 60000,
             G = matrix(0.7, 1, 1), Sigmaeta = matrix(0.1 * s2, 1, 1),
             m0 = as.matrix(0), C0 = as.matrix(1))
mod <- Stem::STEM_Model(z = z, covariates = X, coordinates = coords,
                        phi = phi0, A = matrix(1, d, 1))

## ---- B: pooled reference ----------------------------------------------------
message("B. pooled STEM fit")
pooled <- cached("pooled", {
  Stem::STEM_Estimation(mod, distance = "geo", control = Stem::STEM_control(),
                        alpha = CFG$alpha[1], lambda = CFG$lambda[1])
})

## ---- C: the grid ------------------------------------------------------------
message("C. grid over k and phi")
grid <- cached("grid", {
  Stem::SCSTEM_Infocrit(mod, K_grid = CFG$K_grid, phi_grid = CFG$phi_grid,
                        knn = CFG$knn[1], distance = "geo",
                        seed = CFG$seed[1], verbose = TRUE,
                        alpha = CFG$alpha[1], lambda = CFG$lambda[1],
                        control = Stem::STEM_control())
})

## ---- D: selection -----------------------------------------------------------
message("D. two-step selection")
## Not cached: it takes no time, and it has to follow the rule of the installed
## Stem rather than of whichever version wrote a cache.
sel <- Stem::SCSTEM_Select(grid, band = CFG$band, criterion = "BIC")
best <- sel$fit
stopifnot(inherits(best, "SCSTEM_Estimation"))

if ("--no-bootstrap" %in% commandArgs(trailingOnly = TRUE)) {
  message("\nStopping before the bootstrap, as requested.")
  message("Selected configuration: K = ", sel$K_selected,
          ", phi = ", sel$phi_selected)
  print(table(best$group))
  print(grid$table)
  quit(save = "no", status = 0)
}

## ---- E: bootstrap -----------------------------------------------------------
## With k = 1 the rule finds no partition that improves on the pooled model,
## and there is no clustering to bootstrap.
boot <- inf <- NULL
if (sel$K_selected == 1L) {
  message("E. skipped: the selection rule chooses K = 1, no partition ",
          "improves on the pooled model")
} else {
  message("E. refit-with-clustering bootstrap")
  boot <- cached("bootstrap", {
    Stem::SCSTEM_Bootstrap(best, B = CFG$B[1], seed = CFG$seed[1], verbose = TRUE)
  })
  inf <- cached("bootinference", {
    Stem::SCSTEM_BootInference(boot, level = 0.95)
  })
}

## ---- F: outputs -------------------------------------------------------------
message("F. outputs")
saveRDS(list(pooled = pooled, grid = grid, sel = sel, best = best,
             boot = boot, inf = inf, mod = mod, settings = CFG),
        file.path(OUT, "all.rds"))
utils::write.csv(grid$table, file.path(OUT, "grid.csv"), row.names = FALSE)
utils::write.csv(data.frame(location = seq_len(d), lon = coords[, 1],
                            lat = coords[, 2], regime = best$group),
                 file.path(OUT, "regimes.csv"), row.names = FALSE)
if (!is.null(inf))
  utils::write.csv(inf$summary, file.path(OUT, "estimates.csv"), row.names = FALSE)

message("\nSelected configuration: K = ", sel$K_selected, ", phi = ", sel$phi_selected)
print(table(best$group))
message("\nGrid:")
print(grid$table)
