## ===========================================================================
## SC-STEM fuel application: THE FULL APPLICATION, on the eleven metropolitan
## cities of the application (Naples, Rome, Catania, Turin, Bologna, Milan,
## Palermo, Messina, Bari, Florence, Venice).
##
## Part A  Granger causality between gasoline and diesel at the same pump. For
##         each city and direction, the weekly relative price y* of the fuel
##         explained on its own lags, the lags of the relative price x* of the
##         other fuel, the fiscal pulses, the latent AR(1) process of the regime
##         and the spatially correlated error (notes, eq:model). All the kept
##         pumps of the city.
## Part B  Leaders and followers. For each city, setting and fuel, the weekly
##         relative price of a pump Y on its own lags, the lags of the relative
##         price of its X (the nearest pump of the X set, or the mean of the X
##         set within r*), the fiscal pulses, the latent process and the error:
##           setting 1  Y the major-brand pumps,  X the independents
##           setting 2  Y the independent pumps,  X the major brands
##           setting 3  Y the major-brand pumps,  X the other major brands
##
## Every model:
##   1. the grid k = 1, ..., K_max with the default grid of phi
##      (SCSTEM_Infocrit), K_max = min(3, floor(n / m)) with m = 5 the smallest
##      regime, n the pumps of the model: with fewer than 10 pumps the model is
##      the pooled STEM;
##   2. the two-step selection rule (SCSTEM_Select);
##   3. the refit-with-clustering bootstrap at the selected configuration, and
##      of the pooled model for the comparison (SCSTEM_Bootstrap,
##      SCSTEM_BootInference);
##   4. the Granger test regime by regime: all the lags of x (or X) at zero,
##      judged by the Wald statistic on the bootstrap covariance of their
##      estimates;
##   5. residual diagnostics regime by regime: Ljung-Box on the residuals and
##      on their squares, pump by pump.
##
## Prices: the Sunday price (in force at the end of the week), in cents per
## litre, minus the national mean of the same Sunday, minus the mean of the
## pump over the window (notes, Section "Non-stationarity"). The fiscal events
## enter as pulses, the week of the event and the next.
##
## Inputs: the outputs of fuels-pretreatment.R and, for part B, the pair
## datasets of fuels-leader-follower.R. Every stage of every model is cached in
## <out>/cache/<commit>, so an interrupted run resumes where it stopped; the models can
## be split across R processes, and every process rewrites the summary with
## all the models finished so far.
##
##     Rscript run-fuels-application.R --dry-run       the plan, nothing fitted
##     Rscript run-fuels-application.R                 every model
##     Rscript run-fuels-application.R --job=1/3       one of three processes
##     Rscript run-fuels-application.R --parts=B --cities=NA,RM --settings=2
##     Rscript run-fuels-application.R --B=20          a pilot: the grids (kept in
##                                                     the cache for the full run)
##                                                     and a short bootstrap, with
##                                                     the minutes of every stage
##                                                     in summary.csv
##
## COST. The grid of a model is cheap; the bootstrap is not: every draw re-runs
## the whole SC-STEM fit, partition included. With the Stem of 4325536 some
## draws were pathological: when a regime was small and the spatial range not
## identified, the Newton-Raphson step of the spatial parameters ran to its
## limits at every iteration (on a model of 14 independents one draw of eight
## took about 13 minutes against 3 to 4 seconds for the others). The Stem pinned
## below keeps the range within the limits that the distances identify, which
## removes that cause, and estimates with SQUAREM by default. The cache is kept
## per pinned commit, so that fits of an earlier Stem are never reused.
##
## This is part of the REPLICATION MATERIAL, not of the Stem package.
## ===========================================================================
suppressPackageStartupMessages(library(data.table))
invisible(Sys.setlocale("LC_TIME", "C"))

## ---------------------------------------------------------------------------
## Where this script is, and the command-line overrides --name=value
## ---------------------------------------------------------------------------
fa_here <- local({
  a <- commandArgs(trailingOnly = FALSE); m <- grep("^--file=", a, value = TRUE)
  f <- if (length(m)) normalizePath(sub("^--file=", "", m[1]), winslash = "/") else NULL
  if (is.null(f)) for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) { f <- normalizePath(of, winslash = "/"); break }
  }
  if (is.null(f)) normalizePath(getwd(), winslash = "/") else dirname(f)
})
fa_config <- function(defaults, args = commandArgs(trailingOnly = TRUE)) {
  for (a in grep("^--[^=]+=.", args, value = TRUE)) {
    nm <- gsub("-", "_", sub("^--([^=]+)=.*$", "\\1", a))
    if (!nm %in% names(defaults)) stop("unknown option --", nm, "; the options are: ",
                                       paste(names(defaults), collapse = ", "), call. = FALSE)
    parts <- trimws(strsplit(sub("^--[^=]+=", "", a), ",", fixed = TRUE)[[1]])
    d <- defaults[[nm]]
    defaults[[nm]] <- if (is.logical(d)) as.logical(parts) else if (is.integer(d))
      as.integer(parts) else if (is.numeric(d)) as.numeric(parts) else parts
  }
  defaults
}


## ===========================================================================
## THE SETUP. This is the part meant to be edited.
## ===========================================================================
CFG <- fa_config(list(
  pre         = file.path(fa_here, "fuels", "pretreatment"),     # fuels-pretreatment.R
  pairs       = file.path(fa_here, "fuels", "leader-follower"),  # fuels-leader-follower.R
  out         = file.path(fa_here, "fuels", "application"),
  cities      = "",              # empty: every city of the pre-treatment
  parts       = c("A", "B"),
  settings    = 1:3,             # part B
  fuels       = c("g", "d"),     # part B: g gasoline, d diesel
  x           = "nearest",       # part B: "nearest", or "rstar", the mean within r*
                                 # (the nearest where r* holds no X)
  p           = 2L,              # lags of y and of x (weekly)
  m           = 5L,              # the smallest regime: K_max = min(k_cap, floor(n / m))
  k_cap       = 3L,
  knn         = 5L,
  phi_grid    = c(0, 0.025, 0.05, 0.1, 0.2, 0.5, 1),
  band        = c(0.025, 0.2),
  B           = 100L,            # bootstrap draws (100, decided 2026-10-01)
  boot_pooled = TRUE,            # bootstrap the pooled model too, when k > 1 is selected
  level       = 0.95,
  seed        = 20261002L,
  job         = "1/1"            # this process runs the models i, i + N, ... of "i/N"
))
DRY <- "--dry-run" %in% commandArgs(trailingOnly = TRUE)

## The changes of the excise duties in the window (to be checked against the
## decrees): pulses at the week of the event and the next.
FA_EVENTS <- data.frame(
  date  = as.Date(c("2022-03-22", "2022-12-01", "2023-01-01", "2025-05-15", "2026-01-01")),
  label = c("excise_cut", "cut_reduced", "cut_ends", "realignment", "equal_excise"),
  stringsAsFactors = FALSE)
FA_FUEL <- c(g = "gasoline", d = "diesel")
FA_MAJOR <- c("Agip Eni", "Api-Ip", "Esso", "Q8", "Tamoil")
FA_CITY <- c(RM = "Roma", MI = "Milano", "NA" = "Napoli", TO = "Torino", PA = "Palermo", BA = "Bari",
             CT = "Catania", BO = "Bologna", ME = "Messina", FI = "Firenze", VE = "Venezia",
             GE = "Genova", RC = "Reggio Calabria", CA = "Cagliari")
FA_SETTING <- c("1" = "Y major, X independent", "2" = "Y independent, X major",
                "3" = "Y major, X other major")


## ===========================================================================
## Stem, pinned to ddf90a8: the code of the simulation study (0f7b744, main3: the four
## estimation algorithms with SQUAREM as the default, tolerances 1e-3,
## regularization 0, the range within the limits the distances identify) with
## the BIC on the observed values, which differs only when the response has gaps;
## installed from GitHub when the installed copy is another commit
## ===========================================================================
FA_STEM_REF <- "PaoloMaranzano/Stem@ddf90a80d14b579aac805dd15fbf18c3c6f9f9a7"
fa_stem_ok <- function() {
  if (!nzchar(system.file(package = "Stem"))) return(FALSE)
  identical(utils::packageDescription("Stem")$RemoteSha, sub("^.*@", "", FA_STEM_REF))
}
if (!DRY && !fa_stem_ok()) {
  message("installing Stem at ", FA_STEM_REF)
  if (!requireNamespace("remotes", quietly = TRUE))
    utils::install.packages("remotes", repos = "https://cloud.r-project.org")
  if ("Stem" %in% loadedNamespaces()) try(unloadNamespace("Stem"), silent = TRUE)
  remotes::install_github(FA_STEM_REF, upgrade = "never", force = TRUE, quiet = TRUE)
  if (!fa_stem_ok()) stop("Stem could not be installed at ", FA_STEM_REF, ": restart R and ",
                          "run again, or install it by hand with remotes::install_github()", call. = FALSE)
}

OUT <- normalizePath(CFG$out[1], winslash = "/", mustWork = FALSE)
CACHE <- file.path(OUT, "cache", substr(sub("^.*@", "", FA_STEM_REF), 1L, 7L))
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
## a cached stage keeps the minutes it took, for the summary
cached <- function(name, expr) {
  f <- file.path(CACHE, paste0(name, ".rds"))
  if (file.exists(f)) { message("    [cache] ", name); return(readRDS(f)) }
  message("    [run  ] ", name); t0 <- Sys.time()
  val <- force(expr)
  mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
  message("    [done ] ", name, " in ", mins, " min")
  attr(val, "fa_minutes") <- mins
  saveRDS(val, f)
  val
}
fa_minutes <- function(x) if (is.null(x) || is.null(attr(x, "fa_minutes"))) NA_real_ else attr(x, "fa_minutes")
locf <- function(v) { ok <- !is.na(v); if (!any(ok)) return(v)
  idx <- cummax(ifelse(ok, seq_along(v), 0L)); idx[idx == 0L] <- which(ok)[1]; v[idx] }


## ===========================================================================
## A. The data: weekly relative prices of the kept pumps
## ===========================================================================
ser <- readRDS(file.path(CFG$pre[1], "metro_series.rds"))
mn  <- readRDS(file.path(CFG$pre[1], "means.rds"))
CITIES <- if (any(nzchar(CFG$cities))) CFG$cities else names(FA_CITY)[names(FA_CITY) %in% names(ser)]
miss <- setdiff(CITIES, names(ser))
if (length(miss)) stop("cities not in the pre-treatment output: ", paste(miss, collapse = ", "), call. = FALSE)
days <- ser[[1]]$days; sun <- which(format(days, "%u") == "7"); WEEKS <- days[sun]
nat <- mn$national; setkey(nat, date)

## sites x weeks, in cents per litre: minus the national mean of the Sunday,
## minus the mean of the pump
fa_relative <- function(city, fu) {
  s <- ser[[city]]
  R <- t(apply(if (fu == "g") s$G else s$D, 1, locf))
  C <- 100 * sweep(R, 2, locf(nat[J(days)][[fu]]))[, sun, drop = FALSE]
  C - rowMeans(C)
}
REL <- list()
fa_rel <- function(city, fu) {
  key <- paste(city, fu)
  if (is.null(REL[[key]])) REL[[key]] <<- fa_relative(city, fu)
  REL[[key]]
}

## the fiscal pulses, weeks x 2 events
fa_pulses <- function() {
  w <- vapply(FA_EVENTS$date, function(e) which(WEEKS >= e)[1], 1L)
  P <- matrix(0, length(WEEKS), 2 * nrow(FA_EVENTS),
              dimnames = list(NULL, paste0("pulse_", rep(FA_EVENTS$label, each = 2), c("", "_next"))))
  for (e in seq_along(w)) if (!is.na(w[e])) P[cbind(w[e] + 0:1, 2 * e - 1:0)] <- 1
  P[seq_along(WEEKS), , drop = FALSE]
}
PULSES <- fa_pulses()

## the design of a model: z (T x n) and the covariates stacked by location,
## intercept, the lags of y, the lags of every series in Xs, the pulses
fa_design <- function(Y, Xs) {
  p <- CFG$p[1]; tt <- (p + 1):ncol(Y)
  Pt <- PULSES[tt, , drop = FALSE]; Pt <- Pt[, colSums(Pt) > 0, drop = FALSE]
  lagm <- function(M, i, nm) {
    o <- vapply(seq_len(p), function(l) M[i, tt - l], numeric(length(tt)))
    colnames(o) <- paste0(nm, "_lag", seq_len(p)); o
  }
  covs <- lapply(seq_len(nrow(Y)), function(i)
    cbind(intercept = 1, lagm(Y, i, "y"), do.call(cbind, lapply(names(Xs), function(nm) lagm(Xs[[nm]], i, nm))), Pt))
  list(z = t(Y[, tt, drop = FALSE]), covariates = do.call(rbind, covs), names = colnames(covs[[1]]),
       Tn = length(tt), d = nrow(Y))
}

## the X series of part B, one row per Y: the nearest pump, or the mean of
## the X set within r* (the nearest where r* holds none)
fa_xseries <- function(city, fu, pr) {
  M <- fa_rel(city, fu); ids <- rownames(M)
  rows <- lapply(seq_len(nrow(pr)), function(k) {
    mem <- if (CFG$x[1] == "rstar" && nzchar(pr$members_rstar[k])) strsplit(pr$members_rstar[k], ";")[[1]] else pr$nn_site[k]
    colMeans(M[match(mem, ids), , drop = FALSE])
  })
  do.call(rbind, rows)
}


## ===========================================================================
## B. The list of models
## ===========================================================================
PAIRS <- if ("B" %in% CFG$parts) lapply(setNames(1:3, 1:3), function(s) {
  f <- file.path(CFG$pairs[1], sprintf("setting%d_pairs.rds", s))
  if (!file.exists(f)) stop("part B needs ", f, ": run fuels-leader-follower.R first", call. = FALSE)
  readRDS(f)
}) else NULL

MODELS <- list()
for (city in CITIES) {
  meta <- ser[[city]]$meta
  if ("A" %in% CFG$parts) for (dirn in c("d_to_g", "g_to_d")) {
    y <- substr(dirn, 6, 6); x <- substr(dirn, 1, 1)
    MODELS[[length(MODELS) + 1]] <- list(
      id = sprintf("A_%s_%s-to-%s", city, FA_FUEL[x], FA_FUEL[y]), part = "A", city = city,
      setting = NA_integer_, fuel = FA_FUEL[[y]], other = FA_FUEL[[x]], n = nrow(meta),
      build = local({ cc <- city; yy <- y; xx <- x; function() {
        Y <- fa_rel(cc, yy); list(Y = Y, Xs = list(x = fa_rel(cc, xx)), sites = rownames(Y)) } }))
  }
  if ("B" %in% CFG$parts) for (s in CFG$settings) for (fu in CFG$fuels) {
    pr <- PAIRS[[as.character(s)]][province == city]
    if (!nrow(pr)) next
    MODELS[[length(MODELS) + 1]] <- list(
      id = sprintf("B_%s_s%d_%s", city, s, FA_FUEL[fu]), part = "B", city = city, setting = s,
      fuel = FA_FUEL[[fu]], other = NA_character_, n = nrow(pr),
      build = local({ cc <- city; ff <- fu; pp <- pr; function() {
        M <- fa_rel(cc, ff); Y <- M[match(pp$y_site, rownames(M)), , drop = FALSE]
        list(Y = Y, Xs = list(x = fa_xseries(cc, ff, pp)), sites = pp$y_site) } }))
  }
}
for (k in seq_along(MODELS)) MODELS[[k]]$kmax <- max(1L, min(CFG$k_cap[1], MODELS[[k]]$n %/% CFG$m[1]))

## the share of this process: models ordered by size, dealt in turn
job <- as.integer(strsplit(CFG$job[1], "/", fixed = TRUE)[[1]])
ord <- order(-vapply(MODELS, `[[`, 1, "n"))
MINE <- ord[(seq_along(ord) - 1L) %% job[2] == job[1] - 1L]

plan <- rbindlist(lapply(seq_along(MODELS), function(k) with(MODELS[[k]], data.table(
  model = id, part, city = FA_CITY[city], setting = if (is.na(setting)) "" else FA_SETTING[as.character(setting)],
  fuel, pumps = n, weeks = length(WEEKS) - CFG$p[1], K_max = kmax,
  grid_fits = 1L + (kmax - 1L) * length(CFG$phi_grid),
  boot_refits = CFG$B[1] * (1L + (kmax > 1L && isTRUE(CFG$boot_pooled[1]))),
  this_job = k %in% MINE))))
options(width = 200)
cat(sprintf("\nTHE PLAN: %d models (%d in this job), weeks %s to %s, p = %d, m = %d, B = %d, x = %s\n",
            nrow(plan), length(MINE), format(WEEKS[1 + CFG$p[1]]), format(max(WEEKS)), CFG$p[1], CFG$m[1],
            CFG$B[1], CFG$x[1]))
print(plan[order(-pumps)])
cat(sprintf("\ntotal: %d grid fits, up to %d bootstrap refits\n", sum(plan$grid_fits), sum(plan$boot_refits)))
fwrite(plan, file.path(OUT, "plan.csv"))
if (DRY) quit(save = "no", status = 0)
suppressPackageStartupMessages(library(Stem))


## ===========================================================================
## C. Fitting a model
## ===========================================================================
fa_model <- function(des, coords) {
  ols <- stats::lm.fit(des$covariates, as.vector(des$z))       # z is T x d: stacked by location
  b <- ols$coefficients; b[!is.finite(b)] <- 0
  s2 <- stats::var(ols$residuals)
  dm <- stats::median(geodist::geodist(coords, measure = "geodesic"))
  Stem::STEM_Model(z = des$z, covariates = des$covariates, coordinates = coords,
                   phi = list(beta = matrix(b, ncol = 1), sigma2eps = 0.6 * s2, sigma2omega = 0.4 * s2,
                              theta = 3 / dm, G = matrix(0.3, 1, 1), Sigmaeta = matrix(0.2 * s2, 1, 1),
                              m0 = as.matrix(0), C0 = as.matrix(1)),
                   A = matrix(1, des$d, 1))
}

## the Granger test regime by regime: the Wald statistic of the lags of x on
## the covariance of their aligned bootstrap draws
fa_granger <- function(fit, inf, names, label) {
  j <- grep("^x_lag", names); par <- paste0("beta", j)
  A <- inf$aligned; ip <- match(par, dimnames(A)[[3]])
  k <- nrow(fit$phi_hat)
  rbindlist(lapply(seq_len(k), function(g) {
    est <- unname(fit$phi_hat[g, par])
    D <- matrix(A[, g, ip], ncol = length(ip)); D <- D[stats::complete.cases(D), , drop = FALSE]
    V <- if (nrow(D) > length(ip)) stats::cov(D) else matrix(NA, length(ip), length(ip))
    W <- tryCatch(drop(crossprod(est, solve(V, est))), error = function(e) NA_real_)
    data.table(model = label, regime = g, pumps = sum(fit$group == g),
               coef = paste(sprintf("%.4f", est), collapse = "; "),
               se = paste(sprintf("%.4f", sqrt(diag(V))), collapse = "; "),
               W = W, df = length(ip), p_value = stats::pchisq(W, length(ip), lower.tail = FALSE),
               B_used = nrow(D))
  }))
}

## residual diagnostics, pump by pump, summarized by regime
fa_diagnostics <- function(fit, des) {
  R <- des$z - Stem::SCSTEM_Signal(fit)
  pp <- rbindlist(lapply(seq_len(ncol(R)), function(i) {
    r <- R[, i]; r <- r[is.finite(r)]
    data.table(regime = fit$group[i],
               lb = stats::Box.test(r, lag = 4, type = "Ljung-Box")$p.value,
               lb_sq = stats::Box.test(r^2, lag = 4, type = "Ljung-Box")$p.value,
               ac1 = stats::acf(r, lag.max = 1, plot = FALSE)$acf[2],
               kurt = mean((r - mean(r))^4) / stats::var(r)^2)
  }))
  pp[, .(pumps = .N, `Ljung-Box rejects %` = round(100 * mean(lb < 0.05), 1),
         `Ljung-Box on squares rejects %` = round(100 * mean(lb_sq < 0.05), 1),
         `ACF 1 median` = round(stats::median(ac1), 3), `kurtosis median` = round(stats::median(kurt), 2)),
     by = regime][order(regime)]
}

fa_map <- function(meta, fit, file, title) {
  cols <- c("#1f6f8b", "#e0a458", "#a8516e")[fit$group]
  grDevices::cairo_pdf(file, width = 6.4, height = 6)
  graphics::par(mar = c(1, 1, 2.4, 1))
  graphics::plot(meta$lon, meta$lat, asp = 1 / cos(mean(meta$lat) * pi / 180), pch = 16, col = cols,
                 cex = 0.8, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  graphics::legend("topright", bty = "n", pch = 16, col = c("#1f6f8b", "#e0a458", "#a8516e")[seq_len(max(fit$group))],
                   legend = paste("regime", seq_len(max(fit$group))), cex = 0.8)
  graphics::title(title, cex.main = 0.85)
  invisible(grDevices::dev.off())
}

fa_run <- function(M, seed) {
  message(sprintf("\n== %s: %d pumps, K_max = %d", M$id, M$n, M$kmax))
  dir.create(file.path(OUT, M$id), showWarnings = FALSE)
  b <- M$build(); des <- fa_design(b$Y, b$Xs)
  meta <- ser[[M$city]]$meta[match(b$sites, ser[[M$city]]$meta$site)]
  coords <- cbind(lon = meta$lon, lat = meta$lat)
  dup <- duplicated(round(coords, 6)); coords[dup, ] <- coords[dup, ] + 1e-5 * seq_len(sum(dup))
  tag <- sprintf("%s_p%d_m%d_x%s", M$id, CFG$p[1], CFG$m[1], CFG$x[1])
  ## 1-2. the grid and the selection
  ic <- cached(paste0(tag, "_grid"), {
    mod <- fa_model(des, coords)
    Stem::SCSTEM_Infocrit(mod, K_grid = seq_len(M$kmax), phi_grid = CFG$phi_grid,
                          knn = min(CFG$knn[1], M$n - 1L), distance = "geo",
                          min_cluster_size = CFG$m[1], seed = seed, verbose = TRUE,
                          control = Stem::STEM_control())
  })
  sel <- Stem::SCSTEM_Select(ic, band = CFG$band, criterion = "BIC")
  pooled <- ic$fits[[which(ic$table$K == 1)[1]]]
  ## 3. the bootstrap of the selected model, and of the pooled one
  boot_tag <- sprintf("%s_B%d", tag, CFG$B[1])
  inf_sel <- cached(paste0(boot_tag, "_k", sel$K_selected, "_boot"), {
    bt <- Stem::SCSTEM_Bootstrap(sel$fit, B = CFG$B[1], seed = seed, verbose = TRUE)
    list(info = bt$info, B_valid = bt$B_valid, inf = Stem::SCSTEM_BootInference(bt, level = CFG$level[1]))
  })
  inf_pool <- if (sel$K_selected == 1L) inf_sel else if (isTRUE(CFG$boot_pooled[1]))
    cached(paste0(boot_tag, "_k1_boot"), {
      bt <- Stem::SCSTEM_Bootstrap(pooled, B = CFG$B[1], seed = seed, verbose = TRUE)
      list(info = bt$info, B_valid = bt$B_valid, inf = Stem::SCSTEM_BootInference(bt, level = CFG$level[1]))
    }) else NULL
  ## 4. the Granger tests
  gr <- rbind(fa_granger(sel$fit, inf_sel$inf, des$names, sprintf("selected, K = %d", sel$K_selected)),
              if (sel$K_selected > 1L && !is.null(inf_pool)) fa_granger(pooled, inf_pool$inf, des$names, "pooled"))
  ## 5. the residuals
  dg <- fa_diagnostics(sel$fit, des)
  ## outputs
  od <- file.path(OUT, M$id)
  fwrite(ic$table, file.path(od, "grid.csv"))
  fwrite(gr, file.path(od, "granger.csv"))
  fwrite(dg, file.path(od, "diagnostics.csv"))
  est <- as.data.table(inf_sel$inf$summary)
  est[, covariate := des$names[suppressWarnings(as.integer(sub("^beta", "", parameter)))]]
  fwrite(est, file.path(od, "estimates.csv"))
  stab <- inf_sel$inf$stability
  fwrite(data.table(site = b$sites, brand = meta$brand_end, group = meta$group, lon = meta$lon, lat = meta$lat,
                    regime = sel$fit$group), file.path(od, "regimes.csv"))
  fa_map(meta, sel$fit, file.path(od, "map.pdf"),
         sprintf("%s: K = %d, phi = %s", M$id, sel$K_selected, format(sel$phi_selected)))
  bic <- function(f) unname(f$info_crit[["BIC"]])
  row <- data.table(model = M$id, part = M$part, city = FA_CITY[M$city],
                    setting = if (is.na(M$setting)) "" else FA_SETTING[as.character(M$setting)],
                    fuel = M$fuel, pumps = M$n, K_max = M$kmax, K = sel$K_selected, phi = sel$phi_selected,
                    BIC_pooled = bic(pooled), BIC_selected = bic(sel$fit),
                    p_regimes = paste(format.pval(gr[model != "pooled"]$p_value, digits = 2), collapse = "; "),
                    p_pooled = if (sel$K_selected == 1L) gr$p_value[1] else
                      if (nrow(gr[model == "pooled"])) gr[model == "pooled"]$p_value[1] else NA_real_,
                    B_valid = inf_sel$B_valid,
                    ARI_boot_median = if (!is.null(stab) && nrow(stab)) stats::median(stab$ARI, na.rm = TRUE) else NA_real_,
                    min_grid = fa_minutes(ic), min_boot = fa_minutes(inf_sel),
                    min_boot_pooled = if (sel$K_selected > 1L) fa_minutes(inf_pool) else NA_real_,
                    stem = basename(CACHE))
  saveRDS(row, file.path(od, "summary_row.rds"))
  message(sprintf("   K = %d, phi = %s | Granger p by regime: %s", sel$K_selected,
                  format(sel$phi_selected), row$p_regimes))
  invisible(row)
}


## ===========================================================================
## D. Run this process's models, then the summary of all the models finished
## ===========================================================================
for (k in MINE) {
  res <- tryCatch(fa_run(MODELS[[k]], seed = CFG$seed[1] + k), error = function(e) e)
  if (inherits(res, "error")) message("   FAILED: ", conditionMessage(res))
}
## only the models fitted with the pinned Stem: a row of an earlier run is left out
rows <- lapply(MODELS, function(M) {
  f <- file.path(OUT, M$id, "summary_row.rds")
  if (file.exists(f)) { r <- readRDS(f); if (identical(r$stem, basename(CACHE))) r }
})
summ <- rbindlist(rows, fill = TRUE)
if (nrow(summ)) {
  fwrite(summ, file.path(OUT, "summary.csv"))
  cat(sprintf("\nSUMMARY: %d of %d models finished\n", nrow(summ), length(MODELS)))
  print(summ)
}
cat("\noutputs in", OUT, "\n")
