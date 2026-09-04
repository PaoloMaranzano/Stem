## ---------------------------------------------------------------------------
## Application of the SC-STEM paper: daily PM2.5 over the Po Valley, 2019-2023.
##
## Reproduces Section 5 of the manuscript. The script is staged and every stage
## caches its result, so an interrupted run resumes where it stopped: delete a
## cache file to force that stage to run again.
##
##   Rscript dev/paper/01-application.R
##
## Stages
##   A  the model object and its starting values
##   B  the pooled STEM fit, which is the k = 1 reference
##   C  the (k, phi) grid, through SCSTEM_Infocrit()
##   D  the two-step selection rule
##   E  the refit-with-clustering bootstrap at the selected configuration
##   F  figures and tables
## ---------------------------------------------------------------------------

suppressMessages(pkgload::load_all(
  "C:/Users/paulm/OneDrive/Documenti/GitHub/Stem", quiet = TRUE))

## ---- settings --------------------------------------------------------------
CACHE   <- file.path("dev", "paper", "cache")
FIGDIR  <- "C:/Users/paulm/Dropbox/Applicazioni/Overleaf/SC-STEM package paper/Figures"
K_GRID  <- 1:4
PHI_GRID <- c(0, 0.25, 0.5, 0.75, 1)
BAND    <- c(0.25, 1)          # the moderate-penalty band of the tuning rule
KNN     <- 5
B_BOOT  <- 200
SEED    <- 20260904

dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

cached <- function(name, expr) {
  f <- file.path(CACHE, paste0(name, ".rds"))
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
message("A. model object")
data(povalley)

X <- povalley$covariates
z <- povalley$z
coords <- povalley$coords
d <- ncol(z); Tn <- nrow(z)

## Starting values from the pooled regression. The residual variance is split
## between a nugget and a spatially structured part in the proportion the
## variogram of the residuals suggests; the range starts at 60 km, roughly half
## the median inter-station distance, and the latent process starts persistent.
ols <- stats::lm.fit(x = as.matrix(X), y = as.vector(z))
b0  <- matrix(ols$coefficients, ncol = 1)
s2  <- stats::var(ols$residuals)

phi0 <- list(beta = b0,
             sigma2eps = 0.5 * s2, sigma2omega = 0.4 * s2,
             theta = 1 / 60000,
             G = matrix(0.7, 1, 1), Sigmaeta = matrix(0.1 * s2, 1, 1),
             m0 = as.matrix(0), C0 = as.matrix(1))

mod <- STEM_Model(z = z, covariates = X, coordinates = coords,
                  phi = phi0, K = matrix(1, d, 1))

## ---- B: pooled reference ----------------------------------------------------
message("B. pooled STEM fit")
pooled <- cached("pooled", {
  STEM_Estimation(mod, precision = 0.001, max.iter = 60, distance = "geo")
})

## ---- C: the grid ------------------------------------------------------------
message("C. grid over k and phi")
grid <- cached("grid", {
  SCSTEM_Infocrit(mod, k_grid = K_GRID, phi_grid = PHI_GRID,
                  knn = KNN, distance = "geo",
                  precision = 0.1, precision_full_dataset = 0.01,
                  max_iter = 8, seed = SEED, verbose = TRUE)
})

## ---- D: selection -----------------------------------------------------------
message("D. two-step selection")
sel <- cached("selection", {
  SCSTEM_Select(grid, band = BAND, criterion = "BIC", delta = 0.05)
})

## SCSTEM_Select() carries the selected fit itself, taken from the grid without
## refitting, so there is nothing to look up.
best <- sel$fit
stopifnot(inherits(best, "SCSTEM_Estim"))

## The bootstrap is the expensive stage, so it can be held back until the
## selected configuration has been inspected:
##   Rscript dev/paper/01-application.R --no-bootstrap
if ("--no-bootstrap" %in% commandArgs(trailingOnly = TRUE)) {
  message("\nStopping before the bootstrap, as requested.")
  message("Selected configuration: k = ", sel$k_selected,
          ", phi = ", sel$phi_selected)
  print(table(best$group))
  print(grid$table)
  quit(save = "no", status = 0)
}

## ---- E: bootstrap -----------------------------------------------------------
message("E. refit-with-clustering bootstrap")
boot <- cached("bootstrap", {
  SCSTEM_Bootstrap(best, B = B_BOOT, seed = SEED, verbose = TRUE)
})
inf <- cached("bootinference", {
  SCSTEM_BootInference(boot, level = 0.95)
})

## ---- F: outputs -------------------------------------------------------------
message("F. figures and tables")
saveRDS(list(pooled = pooled, grid = grid, sel = sel, best = best,
             boot = boot, inf = inf, mod = mod),
        file.path(CACHE, "all.rds"))

message("\nSelected configuration: k = ", sel$k_selected, ", phi = ", sel$phi_selected)
print(table(best$group))
message("\nGrid:")
print(grid$table)
