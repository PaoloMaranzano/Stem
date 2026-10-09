## ===========================================================================
## SC-STEM fuel application: AUXILIARY FUNCTIONS
##
## Sourced by the data-management script (fuels-data.R) and by the three main
## scripts of the cases (fuels-main-setting1.R, -setting2.R, -setting3.R).
## Nothing here runs on its own: every function is called from those scripts.
##
## The application in one paragraph. In the eleven metropolitan cities, for
## every pump Y of a set (major brands, or independents) and its neighbour X of
## another set, the weekly relative price of Y on one fuel is regressed on four
## blocks of lags -- its own price on that fuel, its own price on the other
## fuel, the neighbour's price on that fuel, the neighbour's price on the other
## fuel -- plus the pulses of the changes of the excise duties, within an
## SC-STEM model: a latent AR(1) process per regime and a spatially correlated
## error. The three settings:
##
##   setting 1  Y a major-brand pump,  X its neighbour among the independents
##   setting 2  Y an independent pump, X its neighbour among the major brands
##   setting 3  Y a major-brand pump,  X its neighbour among the OTHER major brands
##
## The questions: does the neighbour lead Y (the two neighbour blocks), on the
## same fuel or across fuels; does the relation differ across spatial regimes;
## and, from the same fits, the pass-through of the excise duties by regime
## (the pulse coefficients) and the behaviour on small samples (setting 2,
## 6 to 71 independents per city).
##
## Contents
##   0. settings shared by every script
##   1. weekly relative prices and fiscal pulses
##   2. Kalman completion of the missing weeks
##   3. the radius r* and the pairs (Y, X)
##   4. the data.frame of a case
##   5. the pre-analysis of the lags (ACF, PACF, cross-correlations, lag plots)
##   6. the design and the model of a city and a fuel
##   7. the selection of the lags and of the ridge penalty
##   8. tests, diagnostics and tables of the fitted models
##   9. bookkeeping
##
## Conventions. Prices are in cents per litre. A "relative price" is the
## Sunday price minus the national mean of that Sunday, minus the mean of the
## pump over the window. Matrices of series are sites x weeks. Distances are
## geodesic, in metres for the package and in kilometres for the pairs.
##
## This is part of the REPLICATION MATERIAL of the paper, not of the Stem
## package. The code uses data.table (attached by the scripts) and calls every
## other package with package::function().
## ===========================================================================


## ===========================================================================
## 0. Settings shared by every script
## ===========================================================================

## The commit of Stem the application runs on: installed from GitHub when the
## installed copy is another commit (fu_require_stem()). It is the code of the
## simulation study (0f7b744) with two changes that matter here and not there:
## the BIC counts the observed values (the response has gaps here), and SQUAREM
## also accelerates the fits with a ridge (all the fits of the application).
FU_STEM_REF <- "PaoloMaranzano/Stem@989790bf0d6d321de2cc1d3599ffb397e7db4dc8"

## The eleven cities: the nine with at least 10 kept independents, plus
## Florence and Venice by choice (notes, Section "Data pre-treatment").
FU_CITY <- c(RM = "Roma", MI = "Milano", "NA" = "Napoli", TO = "Torino", PA = "Palermo",
             BA = "Bari", CT = "Catania", BO = "Bologna", ME = "Messina", FI = "Firenze",
             VE = "Venezia")
FU_MAJOR <- c("Agip Eni", "Api-Ip", "Esso", "Q8", "Tamoil")
FU_FUEL <- c(g = "gasoline", d = "diesel")
FU_OTHER <- c(g = "d", d = "g")

## The changes of the excise duties in the window: a pulse at the week of the
## event and one at the next (dates to be checked against the decrees).
FU_EVENTS <- data.frame(
  date  = as.Date(c("2022-03-22", "2022-12-01", "2023-01-01", "2025-05-15", "2026-01-01")),
  label = c("excise_cut", "cut_reduced", "cut_ends", "realignment", "equal_excise"),
  stringsAsFactors = FALSE)

## The three cases: the folder of each under <root>/fuels, and who is Y and X.
FU_CASES <- data.frame(
  setting = 1:3,
  folder  = c("setting1-Ymajor-Xindependent", "setting2-Yindependent-Xmajor",
              "setting3-Ymajor-Xothermajor"),
  y_set   = c("major", "independent", "major"),
  x_set   = c("independent", "major", "other major"),
  stringsAsFactors = FALSE)

## The four blocks of lags of the multi-product model. Each block is built from
## a column of the case data.frame (see fu_case_frame()), for the fuel f of Y
## and the other fuel o:
##   own       Y's own price on f            (column y<f>_c, completed)
##   own_other Y's own price on o            (column y<o>_c, completed)
##   nb        the neighbour's price on f    (column x<f>_<xdef>)
##   nb_other  the neighbour's price on o    (column x<o>_<xdef>)
FU_BLOCKS <- c("own", "own_other", "nb", "nb_other")

## Where the script that sources this file sits: from Rscript, or from Source
## in RStudio.
fu_here <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  f <- if (length(m)) normalizePath(sub("^--file=", "", m[1]), winslash = "/") else NULL
  if (is.null(f)) for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) { f <- normalizePath(of, winslash = "/"); break }
  }
  if (is.null(f)) normalizePath(getwd(), winslash = "/") else dirname(f)
}

## The settings of a script, overridable on the command line with
## --name=value (a dash in the name is read as an underscore; lists are
## separated by commas). A setting that is a list (such as the lags) can only
## be edited in the script.
fu_config <- function(defaults, args = commandArgs(trailingOnly = TRUE)) {
  for (a in grep("^--[^=]+=.", args, value = TRUE)) {
    nm <- gsub("-", "_", sub("^--([^=]+)=.*$", "\\1", a))
    if (!nm %in% names(defaults))
      stop("unknown option --", nm, "; the options are: ", paste(names(defaults), collapse = ", "),
           call. = FALSE)
    parts <- trimws(strsplit(sub("^--[^=]+=", "", a), ",", fixed = TRUE)[[1]])
    d <- defaults[[nm]]
    defaults[[nm]] <- if (is.logical(d)) as.logical(parts) else if (is.integer(d))
      as.integer(parts) else if (is.numeric(d)) as.numeric(parts) else parts
  }
  defaults
}

## Stem at the pinned commit. The commit is read from the installed
## DESCRIPTION file without loading the package; Stem is installed from GitHub
## when it is missing or at another commit. Restart R before a run on a machine
## where another Stem is already loaded.
fu_require_stem <- function(ref = FU_STEM_REF) {
  sha <- sub("^.*@", "", ref)
  ok <- function() nzchar(system.file(package = "Stem")) &&
    identical(utils::packageDescription("Stem")$RemoteSha, sha)
  if (!ok()) {
    message("installing Stem at ", ref)
    if (!requireNamespace("remotes", quietly = TRUE))
      utils::install.packages("remotes", repos = "https://cloud.r-project.org")
    if ("Stem" %in% loadedNamespaces()) try(unloadNamespace("Stem"), silent = TRUE)
    remotes::install_github(ref, upgrade = "never", force = TRUE, quiet = TRUE)
    if (!ok()) stop("Stem could not be installed at ", ref, ": restart R and run again, ",
                    "or install it by hand with remotes::install_github()", call. = FALSE)
  }
  invisible(substr(sha, 1L, 7L))
}

## The folder of a case, created if needed.
fu_case_dir <- function(root, setting) {
  d <- file.path(root, FU_CASES$folder[FU_CASES$setting == setting])
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  normalizePath(d, winslash = "/")
}


## ===========================================================================
## 1. Weekly relative prices and fiscal pulses
## ===========================================================================

## The Sundays of a vector of days (the week runs from Monday to Sunday).
fu_sundays <- function(days) which(format(days, "%u") == "7")

## Weekly relative prices of the sites of a city, sites x weeks, in cents per
## litre. P is the daily price matrix of one fuel (sites x days, NA where the
## site did not report), ref the daily reference mean of the same fuel (the
## national mean of all reporting non-highway pumps, or the metropolitan one).
## The Sunday price is taken as it is: a Sunday without a report stays NA, and
## the Kalman filter of the package fills it (Section 2), instead of carrying
## the last reported price forward. The pump mean is the mean over the weeks
## of the WINDOW the pump reported (`in_window`: a logical over the Sundays;
## by default every Sunday), and it is subtracted from the presample weeks too.
fu_weekly_relative <- function(P, days, ref, in_window = NULL) {
  sun <- fu_sundays(days)
  C <- 100 * sweep(P[, sun, drop = FALSE], 2, ref[sun])
  if (is.null(in_window)) in_window <- rep(TRUE, length(sun))
  C - rowMeans(C[, in_window, drop = FALSE], na.rm = TRUE)
}

## The fiscal pulses, weeks x (2 x events): 1 in the week that contains the
## event (the first Sunday on or after it) and in the next one.
fu_pulses <- function(weeks, events = FU_EVENTS) {
  w <- vapply(events$date, function(e) which(weeks >= e)[1], 1L)
  P <- matrix(0, length(weeks), 2 * nrow(events),
              dimnames = list(NULL, paste0("pulse_", rep(events$label, each = 2), c("", "_next"))))
  for (e in seq_along(w)) if (!is.na(w[e])) {
    rows <- w[e] + 0:1
    rows <- rows[rows <= length(weeks)]
    P[cbind(rows, 2 * e - 1:0)[seq_along(rows), , drop = FALSE]] <- 1
  }
  P
}

## Coordinates of a set of sites as the package wants them (longitude,
## latitude); sites at identical coordinates are moved by a few centimetres, so
## that the spatial covariance stays non-singular.
fu_coords <- function(meta) {
  co <- cbind(lon = meta$lon, lat = meta$lat)
  dup <- duplicated(round(co, 6))
  if (any(dup)) co[dup, ] <- co[dup, ] + 1e-5 * seq_len(sum(dup))
  co
}


## ===========================================================================
## 2. Kalman completion of the missing weeks
## ===========================================================================
## The response of a model may have gaps, which the Kalman filter and smoother
## of the package handle exactly. The covariates may not, and the lags of a
## price are covariates. So the missing weeks are filled once, before any
## model: a pooled STEM model of the relative prices of every kept site of the
## city (intercept and fiscal pulses as covariates, the latent AR(1) common to
## the city and the spatially correlated error), and STEM_Complete(), which
## returns the conditional mean of every missing value given all the observed
## ones -- the smoothed latent state plus the part of the measurement error
## predicted from the sites observed in the same week. The completed series are
## used ONLY to build the lags; the response keeps its gaps.
##
## Y: relative prices, sites x weeks (NA where missing); meta: the sites (lon,
## lat); pulses: weeks x pulses. Returns the completed matrix (sites x weeks)
## and a summary of the fit.
fu_complete <- function(Y, meta, pulses) {
  d <- nrow(Y); Tn <- ncol(Y)
  pul <- pulses[, colSums(pulses) > 0, drop = FALSE]
  X1 <- cbind(intercept = 1, pul)
  X <- do.call(rbind, rep(list(X1), d))                       # stacked by site
  z <- t(Y)                                                   # weeks x sites
  ok <- !is.na(as.vector(z))
  ols <- stats::lm.fit(X[ok, , drop = FALSE], as.vector(z)[ok])
  b <- ols$coefficients; b[!is.finite(b)] <- 0
  s2 <- stats::var(ols$residuals)
  co <- fu_coords(meta)
  dm <- stats::median(geodist::geodist(co, measure = "geodesic"))   # metres
  mod <- Stem::STEM_Model(z = z, covariates = X, coordinates = co,
                          phi = list(beta = matrix(b, ncol = 1), sigma2eps = 0.6 * s2,
                                     sigma2omega = 0.4 * s2, theta = 3 / dm,
                                     G = matrix(0.5, 1, 1), Sigmaeta = matrix(0.2 * s2, 1, 1),
                                     m0 = as.matrix(0), C0 = as.matrix(1)),
                          A = matrix(1, d, 1))
  fit <- Stem::STEM_Estimation(mod, distance = "geo", control = Stem::STEM_control())
  comp <- Stem::STEM_Complete(fit, distance = "geo")          # weeks x sites
  ph <- fit$estimates$phi.hat
  list(completed = t(comp),
       summary = data.frame(sites = d, weeks = Tn, missing = sum(is.na(Y)),
                            missing_pct = round(100 * mean(is.na(Y)), 3),
                            loglik = fit$estimates$loglik,
                            converged = isTRUE(fit$estimates$convergence.par$converged),
                            G = ph$G[1, 1], sigma2eps = ph$sigma2eps, sigma2omega = ph$sigma2omega,
                            range_km = 3 / ph$theta / 1000))
}


## ===========================================================================
## 3. The radius r* and the pairs (Y, X)
## ===========================================================================

## The correlogram of a city: pairwise Pearson correlation of the weekly prices
## centred on the METROPOLITAN mean (the cleaner reference for local
## interaction, notes Section "Exploration"), against distance, binned. M:
## sites x weeks; Dkm: distances in km. Pairwise-complete weeks.
fu_correlogram <- function(M, Dkm, bins = c(0, .5, 1, 1.5, 2, 3, 4, 5, 7.5, 10, 15, 20, 30)) {
  rr <- suppressWarnings(stats::cor(t(M), use = "pairwise.complete.obs"))
  up <- upper.tri(Dkm)
  z <- data.table(d = Dkm[up], r = rr[up])[d <= max(bins) & is.finite(r)]
  z[, bin := cut(d, bins, include.lowest = TRUE)]
  z[, .(n = .N, d = stats::median(d), r = stats::median(r)), by = bin][order(d)]
}

## r*: the distance at which the excess of the binned median correlation over
## its plateau (pairs beyond 10 km) halves, interpolated between bins. cg: the
## correlogram of a city, both fuels stacked.
fu_rstar <- function(cg) {
  z <- cg[, .(d = mean(d), r = mean(r), n = sum(n)), by = bin][order(d)]
  a0 <- stats::median(z$r[z$d > 10]); b0 <- z$r[1] - a0; h <- a0 + b0 / 2
  k <- which(z$r < h)[1]
  rs <- if (is.na(k)) NA_real_ else if (k == 1) z$d[1] else
    z$d[k - 1] + (z$d[k] - z$d[k - 1]) * (z$r[k - 1] - h) / (z$r[k - 1] - z$r[k])
  list(r_star = if (is.finite(rs)) round(rs, 1) else 1, plateau = a0, excess_0 = b0)
}

## The pairs of a setting in a city: for every Y of the Y set, the nearest pump
## of the X set and every pump of the X set within r*. Pairs stay within the
## city. meta: the kept sites of the city; Dkm: their distances in km.
## Returns a list: `table`, one row per Y (rows of meta), and `rs`, for each Y
## the rows of its X within r* (the nearest X when r* holds none).
fu_pairs <- function(meta, Dkm, rstar, setting) {
  maj <- which(meta$brand_end %in% FU_MAJOR)
  ind <- which(meta$group == "independent")
  ys <- switch(as.character(setting), "1" = maj, "2" = ind, "3" = maj)
  xset <- function(i) switch(as.character(setting), "1" = ind, "2" = maj,
                             "3" = maj[meta$brand_end[maj] != meta$brand_end[i]])
  tab <- list(); rs <- list()
  for (i in ys) {
    xs <- xset(i)
    if (!length(xs)) next
    dx <- Dkm[i, xs]
    nn <- xs[which.min(dx)]
    mem <- xs[dx <= rstar]
    tab[[length(tab) + 1]] <- data.frame(y_row = i, nn_row = nn, nn_dist_km = round(min(dx), 3),
                                         r_star_km = rstar, n_rstar = length(mem),
                                         rs_fallback = !length(mem), stringsAsFactors = FALSE)
    rs[[length(rs) + 1]] <- if (length(mem)) mem else nn
  }
  list(table = do.call(rbind, tab), rs = rs)
}


## ===========================================================================
## 4. The data.frame of a case
## ===========================================================================
## One row per city, pump Y and week. Columns:
##   city, city_name, setting               the case
##   y_site, y_brand, y_group, y_type,      the pump Y and its attributes
##   y_municipality, y_lon, y_lat
##   week, t                                the Sunday of the week, its index 1..W
##   presample                              TRUE in the 52 weeks before the window,
##                                          which give lags only
##   yg, yd                                 Y's relative prices: the RESPONSE,
##                                          NA where the Sunday price is missing
##   yg_c, yd_c                             the same, completed (for the lags)
##   xg_nn, xd_nn                           the nearest X: its completed relative prices
##   xg_rs, xd_rs                           the mean of the X within r*, completed
##                                          (the nearest when r* holds none)
##   nn_site, nn_brand, nn_dist_km          the nearest X
##   r_star_km, n_rstar, rs_fallback        the radius, how many X it holds,
##                                          whether the nearest stands in
##   pulse_*                                the fiscal pulses
## Yrel and Ycomp: lists by fuel ("g", "d") of sites x weeks matrices of the
## city (relative and completed); meta: the kept sites; pairs: fu_pairs().
fu_case_frame <- function(city, setting, meta, Yrel, Ycomp, weeks, pairs, pulses, presample = rep(FALSE, length(weeks))) {
  W <- length(weeks)
  pt <- pairs$table
  one <- function(k) {
    i <- pt$y_row[k]; nn <- pt$nn_row[k]; rs <- pairs$rs[[k]]
    xr <- function(f) colMeans(Ycomp[[f]][rs, , drop = FALSE])
    data.frame(
      city = city, city_name = unname(FU_CITY[city]), setting = setting,
      y_site = meta$site[i], y_brand = meta$brand_end[i], y_group = meta$group[i],
      y_type = meta$type[i], y_municipality = meta$city[i], y_lon = meta$lon[i], y_lat = meta$lat[i],
      week = weeks, t = seq_len(W), presample = presample,
      yg = Yrel$g[i, ], yd = Yrel$d[i, ],
      yg_c = Ycomp$g[i, ], yd_c = Ycomp$d[i, ],
      xg_nn = Ycomp$g[nn, ], xd_nn = Ycomp$d[nn, ],
      xg_rs = xr("g"), xd_rs = xr("d"),
      nn_site = meta$site[nn], nn_brand = meta$brand_end[nn], nn_dist_km = pt$nn_dist_km[k],
      r_star_km = pt$r_star_km[k], n_rstar = pt$n_rstar[k], rs_fallback = pt$rs_fallback[k],
      pulses, stringsAsFactors = FALSE, row.names = NULL)
  }
  do.call(rbind, lapply(seq_len(nrow(pt)), one))
}


## ===========================================================================
## 5. The pre-analysis of the lags
## ===========================================================================
## Pump by pump, before any model: which lags of each block carry information.
## For every pump Y of a city and the fuel f of Y, on the completed relative
## prices:
##   ACF  of Y's own price on f                 (own block)
##   PACF of Y's own price on f                 (own block, net of the shorter lags)
##   CCF  cor(y_t, z_{t-l}) with z = Y's own price on the other fuel (own_other),
##        the neighbour's price on f (nb) and on the other fuel (nb_other)
## at the lags `lags` (by default 1 to 8 and the annual lag 52). For each lag:
## the median, the quartiles and the share of pumps whose coefficient exceeds
## the band 2 / sqrt(T) in absolute value, against about 5% by chance.
fu_lag_stats <- function(dfc, fuel, xdef, lags = c(1:8, 52)) {
  S <- fu_series(dfc, fuel, xdef)
  Tn <- nrow(S$own); band <- 2 / sqrt(Tn)
  ccf1 <- function(y, z, l) suppressWarnings(stats::cor(y[(l + 1):Tn], z[1:(Tn - l)]))
  per_pump <- lapply(seq_len(ncol(S$own)), function(i) {
    y <- S$own[, i]
    a <- stats::acf(y, lag.max = max(lags), plot = FALSE)$acf[lags + 1]
    pa <- stats::pacf(y, lag.max = max(lags), plot = FALSE)$acf[lags]
    cbind(acf_own = a, pacf_own = pa,
          ccf_own_other = vapply(lags, function(l) ccf1(y, S$own_other[, i], l), 1),
          ccf_nb = vapply(lags, function(l) ccf1(y, S$nb[, i], l), 1),
          ccf_nb_other = vapply(lags, function(l) ccf1(y, S$nb_other[, i], l), 1))
  })
  A <- simplify2array(per_pump)                              # lags x measures x pumps
  out <- lapply(dimnames(A)[[2]], function(ms) {
    V <- A[, ms, , drop = TRUE]
    if (is.null(dim(V))) V <- matrix(V, nrow = length(lags))
    data.frame(measure = ms, lag = lags,
               median = apply(V, 1, stats::median, na.rm = TRUE),
               q25 = apply(V, 1, stats::quantile, 0.25, na.rm = TRUE),
               q75 = apply(V, 1, stats::quantile, 0.75, na.rm = TRUE),
               share_sig = apply(abs(V) > band, 1, mean, na.rm = TRUE),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  out$city <- dfc$city[1]; out$fuel <- FU_FUEL[[fuel]]; out$pumps <- ncol(S$own); out$weeks <- Tn
  out[, c("city", "fuel", "pumps", "weeks", "measure", "lag", "median", "q25", "q75", "share_sig")]
}

## The four series of the blocks of a city and a fuel, weeks x pumps.
fu_series <- function(dfc, fuel, xdef) {
  fo <- FU_OTHER[[fuel]]
  dfc <- dfc[order(match(dfc$y_site, unique(dfc$y_site)), dfc$t), ]
  W <- max(dfc$t)
  m <- function(col) matrix(dfc[[col]], nrow = W)
  list(resp = m(paste0("y", fuel)), own = m(paste0("y", fuel, "_c")),
       own_other = m(paste0("y", fo, "_c")),
       nb = m(paste0("x", fuel, "_", xdef)), nb_other = m(paste0("x", fo, "_", xdef)),
       sites = unique(dfc$y_site))
}

## One page of plots per city and fuel: the median ACF, PACF and the three
## cross-correlations against the lag, with the interquartile band of the
## pumps and the band 2 / sqrt(T); and lag plots, pooled over the pumps of the
## city (smoothed scatter): y_t against y_{t-1}, y_{t-2}, y_{t-4}, y_{t-52},
## and against the neighbour's price on the same fuel at t-1.
fu_lag_plots <- function(dfc, fuel, xdef, st) {
  S <- fu_series(dfc, fuel, xdef)
  Tn <- nrow(S$own); band <- 2 / sqrt(Tn)
  graphics::par(mfrow = c(2, 5), mar = c(3.4, 3.4, 2.2, 0.6), mgp = c(2.1, 0.6, 0))
  for (ms in c("acf_own", "pacf_own", "ccf_own_other", "ccf_nb", "ccf_nb_other")) {
    z <- st[st$measure == ms, ]
    x <- seq_along(z$lag)
    graphics::plot(x, z$median, type = "n", ylim = range(c(z$q25, z$q75, band, -band), na.rm = TRUE),
                   xaxt = "n", las = 1, bty = "n", xlab = "lag (weeks)", ylab = "coefficient")
    graphics::axis(1, at = x, labels = z$lag)
    graphics::polygon(c(x, rev(x)), c(z$q25, rev(z$q75)), col = "#1f6f8b33", border = NA)
    graphics::lines(x, z$median, type = "b", pch = 16, col = "#1f6f8b")
    graphics::abline(h = c(-band, 0, band), lty = c(3, 1, 3), col = "#8a939f")
    graphics::title(sprintf("%s\nshare beyond band: %s", ms,
                            paste(round(100 * z$share_sig[z$lag %in% c(1, 2, 4, 52)]), collapse = "/")),
                    cex.main = 0.75)
  }
  lp <- function(l, z, lab) {
    y <- as.vector(S$own[(l + 1):Tn, ]); x <- as.vector(z[1:(Tn - l), ])
    ok <- is.finite(x) & is.finite(y)
    graphics::smoothScatter(x[ok], y[ok], xlab = lab, ylab = "y_t", nrpoints = 0, bty = "n")
    graphics::abline(stats::lm(y[ok] ~ x[ok]), col = "#a8516e")
    graphics::title(sprintf("cor %.2f", stats::cor(x[ok], y[ok])), cex.main = 0.8)
  }
  for (l in c(1, 2, 4)) lp(l, S$own, sprintf("y_{t-%d}", l))
  if (Tn > 60) lp(52, S$own, "y_{t-52}") else graphics::plot.new()
  lp(1, S$nb, "neighbour, same fuel, t-1")
  graphics::mtext(sprintf("%s, %s (Y), X %s: %d pumps, %d weeks", FU_CITY[dfc$city[1]],
                          FU_FUEL[[fuel]], xdef, ncol(S$own), Tn),
                  side = 3, outer = TRUE, line = -1.2, cex = 0.8)
}


## ===========================================================================
## 6. The design and the model of a city and a fuel
## ===========================================================================
## The design of one model: the response (weeks x pumps, with gaps), the
## covariates stacked by pump (the convention of STEM_Model()), their names,
## the columns of each block, and the columns the ridge acts on.
##
## lags: a named list with one integer vector per block of FU_BLOCKS, the lags
## of that block; lags may be non-contiguous (for instance c(1, 2, 4, 52)), and
## integer(0) leaves a block out. t_start: the first week of the response; by
## default the first week of the window (the weeks before it are the presample,
## which gives lags only), or the largest lag plus one if that is later. The
## selection of the lags passes a common t_start, so that the candidates are
## compared on the same weeks.
##
## NOTE on the annual lag. The presample of 52 weeks (2021) gives every lag up
## to 52 to the first week of the window: lag 52 costs no week of the window,
## and the pulses of the excise events of 2022 stay in. A lag beyond 52 would
## move the start of the response into the window.
fu_design <- function(dfc, fuel, xdef, lags, t_start = NULL) {
  S <- fu_series(dfc, fuel, xdef)
  W <- nrow(S$own); d <- ncol(S$own)
  lags <- lags[FU_BLOCKS]
  maxlag <- max(c(0L, unlist(lags)))
  first <- if (is.null(dfc$presample)) 1L else min(dfc$t[!dfc$presample])
  t0 <- if (is.null(t_start)) max(maxlag + 1L, first) else t_start
  if (t0 <= maxlag) stop("t_start must exceed the largest lag", call. = FALSE)
  tt <- t0:W
  pul <- as.matrix(dfc[dfc$y_site == S$sites[1], grep("^pulse_", names(dfc)), drop = FALSE])
  pul <- pul[order(dfc$t[dfc$y_site == S$sites[1]]), , drop = FALSE][tt, , drop = FALSE]
  pul <- pul[, colSums(pul) > 0, drop = FALSE]
  block_cols <- function(i) {
    out <- lapply(FU_BLOCKS, function(b) {
      L <- lags[[b]]
      if (!length(L)) return(NULL)
      m <- vapply(L, function(l) S[[b]][tt - l, i], numeric(length(tt)))
      if (is.null(dim(m))) m <- matrix(m, ncol = length(L))
      colnames(m) <- paste0(b, "_l", L)
      m
    })
    do.call(cbind, out)
  }
  covs <- lapply(seq_len(d), function(i) cbind(intercept = 1, block_cols(i), pul))
  nm <- colnames(covs[[1]])
  blocks <- lapply(stats::setNames(FU_BLOCKS, FU_BLOCKS), function(b) grep(paste0("^", b, "_l"), nm))
  blocks$pulses <- grep("^pulse_", nm)
  meta <- unique(dfc[, c("y_site", "y_brand", "y_group", "y_type", "y_municipality", "y_lon", "y_lat",
                         "nn_site", "nn_brand", "nn_dist_km", "n_rstar", "rs_fallback")])
  meta <- meta[match(S$sites, meta$y_site), ]
  list(z = S$resp[tt, , drop = FALSE], covariates = do.call(rbind, covs), names = nm,
       blocks = blocks, penalize = unlist(blocks[FU_BLOCKS], use.names = FALSE),
       coords = fu_coords(data.frame(lon = meta$y_lon, lat = meta$y_lat)),
       meta = meta, weeks = sort(unique(dfc$week))[tt], d = d, Tn = length(tt),
       city = dfc$city[1], fuel = fuel, xdef = xdef, lags = lags)
}

## The STEM_Model of a design, with starting values from least squares on the
## observed weeks; the range starts at the median distance between the pumps.
fu_model <- function(des) {
  zv <- as.vector(des$z)                        # weeks x pumps: stacked by pump
  ok <- !is.na(zv)
  ols <- stats::lm.fit(des$covariates[ok, , drop = FALSE], zv[ok])
  b <- ols$coefficients; b[!is.finite(b)] <- 0
  s2 <- stats::var(ols$residuals)
  dm <- if (des$d > 1) stats::median(geodist::geodist(des$coords, measure = "geodesic")) else 1000
  Stem::STEM_Model(z = des$z, covariates = des$covariates, coordinates = des$coords,
                   phi = list(beta = matrix(b, ncol = 1), sigma2eps = 0.6 * s2, sigma2omega = 0.4 * s2,
                              theta = 3 / dm, G = matrix(0.3, 1, 1), Sigmaeta = matrix(0.2 * s2, 1, 1),
                              m0 = as.matrix(0), C0 = as.matrix(1)),
                   A = matrix(1, des$d, 1))
}

## The largest number of regimes a model can have: at least m pumps each, at
## most k_cap; with fewer than 2m pumps, the pooled STEM model only.
fu_kmax <- function(n, m, k_cap) max(1L, min(as.integer(k_cap), n %/% m))


## ===========================================================================
## 7. The selection of the lags and of the ridge penalty
## ===========================================================================
## On the POOLED STEM model of a city and a fuel (cheap: no clustering), every
## candidate set of lags is crossed with every value of the ridge penalty, and
## the BIC of each is recorded; the BIC counts the effective number of
## coefficients the ridge leaves. All the candidates are fitted on the same
## weeks (the window, or from the largest lag of all the candidates plus one if later), so that their
## log-likelihoods are comparable. The ridge acts on the four blocks of lags
## only, never on the intercept and the fiscal pulses. lambda is on the
## relative scale of the package: on an orthogonal design a coefficient is
## multiplied by 1 / (1 + lambda). The main scripts pass one set of lags (LAGS),
## so that only the penalty is chosen, city by city and fuel by fuel.
fu_select <-function(dfc, fuel, xdef, candidates, lambdas, verbose = FALSE) {
  first <- if (is.null(dfc$presample)) 1L else min(dfc$t[!dfc$presample])
  t0 <- max(max(unlist(candidates)) + 1L, first)
  out <- list()
  for (cn in names(candidates)) {
    des <- fu_design(dfc, fuel, xdef, candidates[[cn]], t_start = t0)
    mod <- fu_model(des)
    for (lam in lambdas) {
      f <- tryCatch(suppressWarnings(Stem::SCSTEM_Estimation(
        mod, K = 1, distance = "geo", alpha = 0, lambda = lam, penalize = des$penalize,
        control = Stem::STEM_control(), verbose = verbose)), error = function(e) e)
      out[[length(out) + 1]] <- if (inherits(f, "error")) {
        data.frame(candidate = cn, lambda = lam, loglik = NA_real_, df = NA_real_, BIC = NA_real_,
                   error = conditionMessage(f), stringsAsFactors = FALSE)
      } else {
        data.frame(candidate = cn, lambda = lam, loglik = f$info_crit[["loglik"]],
                   df = f$info_crit[["df"]], BIC = f$info_crit[["BIC"]], error = "",
                   stringsAsFactors = FALSE)
      }
    }
  }
  res <- do.call(rbind, out)
  res$city <- dfc$city[1]; res$fuel <- FU_FUEL[[fuel]]; res$xdef <- xdef; res$t_start <- t0
  res$selected <- FALSE
  if (any(is.finite(res$BIC))) res$selected[which.min(res$BIC)] <- TRUE
  res
}


## ===========================================================================
## 8. Tests, diagnostics and tables of the fitted models
## ===========================================================================

## The hypotheses tested, as sets of columns of the design. Each is the null
## that the coefficients of those lags are all zero:
##   nb_both    the neighbour does not lead Y (its two fuels together): THE
##              leader-follower question
##   nb         the neighbour does not lead Y on the same fuel
##   nb_other   the neighbour does not lead Y across fuels
##   own_other  Y's own price on the other fuel does not lead its price on this
##              one (the multi-product link within the station)
fu_tests <- function(des) {
  b <- des$blocks
  out <- list(nb_both = c(b$nb, b$nb_other), nb = b$nb, nb_other = b$nb_other, own_other = b$own_other)
  out[vapply(out, length, 1L) > 0]
}

## The Wald test of a block of coefficients, regime by regime, on the
## covariance of their aligned bootstrap draws (the refit-with-clustering
## bootstrap of the package, which re-estimates the partition at every draw).
## cols: the columns of the block in the design. With the ridge the estimates
## are shrunk towards zero, so the test is conservative.
fu_wald_boot <- function(fit, inf, des, cols, label) {
  if (!length(cols)) return(NULL)
  par <- paste0("beta", cols)
  A <- inf$aligned
  ip <- match(par, dimnames(A)[[3]])
  K <- nrow(fit$phi_hat)
  do.call(rbind, lapply(seq_len(K), function(g) {
    est <- unname(fit$phi_hat[g, par])
    D <- matrix(A[, g, ip], ncol = length(ip)); D <- D[stats::complete.cases(D), , drop = FALSE]
    V <- if (nrow(D) > length(ip)) stats::cov(D) else matrix(NA_real_, length(ip), length(ip))
    W <- tryCatch(drop(crossprod(est, solve(V, est))), error = function(e) NA_real_)
    data.frame(block = label, regime = g, pumps = sum(fit$group == g),
               coef = paste(sprintf("%.4f", est), collapse = "; "),
               se_boot = paste(sprintf("%.4f", sqrt(diag(V))), collapse = "; "),
               W = W, df = length(ip), p_value = stats::pchisq(W, length(ip), lower.tail = FALSE),
               draws_used = nrow(D), stringsAsFactors = FALSE)
  }))
}

## The likelihood-ratio test of each hypothesis of fu_tests() on the SELECTED
## partition, without the ridge: the model with and without the lags tested,
## refitted regime by regime on
## the partition held fixed (no alternation), and 2 (l_full - l_restricted) for
## each regime and in total, against a chi-square with as many degrees of
## freedom as lags are tested. It is the likelihood counterpart of the F
## test; on the first test of the application (Isernia) it rejected more often
## than the pump-level tests, so it is read next to the bootstrap Wald test.
fu_lr_fixed <- function(des, partition, K, cfg, seed) {
  refit <- function(d) {
    mod <- fu_model(d)
    if (K == 1L) {
      Stem::SCSTEM_Estimation(mod, K = 1, distance = "geo", lambda = 0,
                              control = Stem::STEM_control(), verbose = FALSE)
    } else {
      Stem::SCSTEM_Estimation(mod, K = K, phi_penalty = 0, init_partition = partition, max_iter = 0,
                              distance = "geo", knn = min(cfg$knn[1], d$d - 1L),
                              min_cluster_size = cfg$m[1], lambda = 0, seed = seed,
                              control = Stem::STEM_control(), verbose = FALSE)
    }
  }
  full <- refit(des)
  tests <- fu_tests(des)
  rows <- lapply(names(tests), function(b) {
    cols <- tests[[b]]
    dr <- fu_subset_design(des, setdiff(seq_along(des$names), cols))
    rest <- tryCatch(refit(dr), error = function(e) NULL)
    if (is.null(rest)) return(data.frame(block = b, regime = NA, LR = NA, df = length(cols), p_value = NA))
    lr_g <- 2 * (full$loglik_g - rest$loglik_g)
    data.frame(block = b, regime = c(seq_len(K), 0L),
               LR = c(lr_g, sum(lr_g)), df = c(rep(length(cols), K), K * length(cols)),
               p_value = stats::pchisq(c(lr_g, sum(lr_g)), c(rep(length(cols), K), K * length(cols)),
                                       lower.tail = FALSE))
  })
  out <- do.call(rbind, rows)
  out$regime <- ifelse(out$regime == 0L, "all", as.character(out$regime))
  out
}

## A design restricted to some of its columns (keep: indices), with the blocks
## and the penalized columns re-indexed.
fu_subset_design <- function(des, keep) {
  des$covariates <- des$covariates[, keep, drop = FALSE]
  des$names <- des$names[keep]
  des$blocks <- lapply(stats::setNames(c(FU_BLOCKS, "pulses"), c(FU_BLOCKS, "pulses")), function(b)
    if (b == "pulses") grep("^pulse_", des$names) else grep(paste0("^", b, "_l"), des$names))
  des$penalize <- unlist(des$blocks[FU_BLOCKS], use.names = FALSE)
  des
}

## The pump-by-pump benchmark: for every pump Y, the same regression by least
## squares on the weeks with an observed response (own lags, own other fuel,
## neighbour, fiscal pulses; no latent process, no spatial correlation), and for
## every hypothesis of fu_tests() the classical F test and the Wald test with the
## heteroskedasticity-consistent covariance HC1 (White 1980; MacKinnon and
## White 1985). With sticky prices the error variance is far from constant, and
## the F test over-rejects (notes, Section "Exploration"): the two are reported
## side by side.
fu_pump_tests <- function(des) {
  rows <- lapply(seq_len(des$d), function(i) {
    Z <- des$covariates[(i - 1L) * des$Tn + seq_len(des$Tn), , drop = FALSE]
    y <- des$z[, i]; ok <- is.finite(y)
    Z <- Z[ok, , drop = FALSE]; y <- y[ok]
    q <- qr(Z); keep <- sort(q$pivot[seq_len(q$rank)])
    Zk <- Z[, keep, drop = FALSE]; n <- nrow(Zk); k <- ncol(Zk)
    A <- tryCatch(solve(crossprod(Zk)), error = function(e) NULL)
    if (is.null(A) || n <= k) return(NULL)
    b <- A %*% crossprod(Zk, y); e <- drop(y - Zk %*% b)
    V_hc1 <- A %*% crossprod(Zk * e) %*% A * n / (n - k)
    s2 <- sum(e^2) / (n - k)
    tests <- fu_tests(des)
    do.call(rbind, lapply(names(tests), function(bl) {
      j <- match(tests[[bl]], keep); j <- j[!is.na(j)]
      if (!length(j)) return(NULL)
      Fs <- drop(crossprod(b[j], solve(A[j, j, drop = FALSE], b[j]))) / (s2 * length(j))
      Wh <- drop(crossprod(b[j], solve(V_hc1[j, j, drop = FALSE], b[j])))
      data.frame(y_site = des$meta$y_site[i], block = bl, q = length(j), n = n, k = k,
                 F = Fs, p_F = stats::pf(Fs, length(j), n - k, lower.tail = FALSE),
                 W_hc1 = Wh, p_hc1 = stats::pchisq(Wh, length(j), lower.tail = FALSE),
                 stringsAsFactors = FALSE)
    }))
  })
  do.call(rbind, rows[!vapply(rows, is.null, TRUE)])
}

## The heterogeneous-panel test of Dumitrescu and Hurlin (2012) on the
## pump-level Wald statistics of a block: Zbar = sqrt(N / (2 q)) (mean(W) - q),
## asymptotically standard normal under non-causality for every pump. Computed
## on the classical Wald statistic (q times F) and on the HC1 one.
fu_dh <- function(pt) {
  do.call(rbind, lapply(split(pt, pt$block), function(z) {
    q <- z$q[1]; N <- nrow(z)
    zb <- function(W) sqrt(N / (2 * q)) * (mean(W, na.rm = TRUE) - q)
    zf <- zb(q * z$F); zh <- zb(z$W_hc1)
    data.frame(block = z$block[1], pumps = N, q = q,
               reject_F_pct = round(100 * mean(z$p_F < 0.05, na.rm = TRUE), 1),
               reject_hc1_pct = round(100 * mean(z$p_hc1 < 0.05, na.rm = TRUE), 1),
               Zbar_F = zf, p_F = 2 * stats::pnorm(-abs(zf)),
               Zbar_hc1 = zh, p_hc1 = 2 * stats::pnorm(-abs(zh)), stringsAsFactors = FALSE)
  }))
}

## Residual diagnostics, pump by pump, summarized by regime: the Ljung-Box test
## on the residuals and on their squares (4 lags), the first autocorrelation
## and the kurtosis. Residuals are taken on the weeks with an observed response.
fu_diagnostics <- function(fit, des) {
  R <- des$z - Stem::SCSTEM_Signal(fit)
  pp <- do.call(rbind, lapply(seq_len(ncol(R)), function(i) {
    r <- R[, i]; r <- r[is.finite(r)]
    data.frame(regime = fit$group[i],
               lb = stats::Box.test(r, lag = 4, type = "Ljung-Box")$p.value,
               lb_sq = stats::Box.test(r^2, lag = 4, type = "Ljung-Box")$p.value,
               ac1 = stats::acf(r, lag.max = 1, plot = FALSE)$acf[2],
               kurt = mean((r - mean(r))^4) / stats::var(r)^2)
  }))
  do.call(rbind, lapply(split(pp, pp$regime), function(z) data.frame(
    regime = z$regime[1], pumps = nrow(z),
    ljung_box_reject_pct = round(100 * mean(z$lb < 0.05), 1),
    ljung_box_sq_reject_pct = round(100 * mean(z$lb_sq < 0.05), 1),
    acf1_median = round(stats::median(z$ac1), 3), kurtosis_median = round(stats::median(z$kurt), 2))))
}

## The pass-through of the excise duties by regime: the coefficients of the
## pulses (unpenalized), with their bootstrap standard errors and percentile
## intervals. inf: SCSTEM_BootInference().
fu_pulse_table <- function(inf, des) {
  cols <- des$blocks$pulses
  if (!length(cols)) return(NULL)
  s <- as.data.frame(inf$summary)
  s <- s[s$parameter %in% paste0("beta", cols), , drop = FALSE]
  s$pulse <- des$names[as.integer(sub("^beta", "", s$parameter))]
  s
}

## The Adjusted Rand Index between two partitions (to compare the regimes with
## the brands).
fu_ari <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  c2 <- function(v) sum(v * (v - 1) / 2)
  s_ij <- c2(as.vector(tab)); s_i <- c2(rowSums(tab)); s_j <- c2(colSums(tab)); tot <- n * (n - 1) / 2
  e <- s_i * s_j / tot; mx <- (s_i + s_j) / 2
  if (isTRUE(all.equal(mx, e))) return(1)
  (s_ij - e) / (mx - e)
}

## Are the regimes spatial, or do they follow the brands? The regimes against
## the brands of the pumps Y, and their Adjusted Rand Index (about 0 when the
## regimes are unrelated to the brands). Meaningful where Y has several brands
## (settings 1 and 3).
fu_brand_table <- function(fit, des) {
  tab <- as.data.frame.matrix(table(regime = fit$group, brand = des$meta$y_brand))
  tab$regime <- rownames(tab)
  list(table = tab[, c("regime", setdiff(names(tab), "regime"))],
       ari = if (length(unique(des$meta$y_brand)) > 1 && length(unique(fit$group)) > 1)
               fu_ari(fit$group, des$meta$y_brand) else NA_real_)
}

## The map of the regimes of a fit: the pumps Y coloured by regime, their
## nearest X in grey.
fu_map <- function(des, fit, file, title) {
  pal <- c("#1f6f8b", "#e0a458", "#a8516e")
  grDevices::cairo_pdf(file, width = 6.4, height = 6)
  graphics::par(mar = c(1, 1, 2.4, 1))
  m <- des$meta
  graphics::plot(m$y_lon, m$y_lat, asp = 1 / cos(mean(m$y_lat) * pi / 180), pch = 16,
                 col = pal[fit$group], cex = 0.8, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  graphics::legend("topright", bty = "n", pch = 16, col = pal[seq_len(max(fit$group))],
                   legend = paste("regime", seq_len(max(fit$group))), cex = 0.8)
  graphics::title(title, cex.main = 0.85)
  invisible(grDevices::dev.off())
}


## ===========================================================================
## 9. Bookkeeping
## ===========================================================================

## A cached stage: its value is saved in <dir>/<name>.rds and read back on a
## later run, with the minutes it took; delete the file to run the stage again.
fu_cached <- function(dir, name, expr) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(dir, paste0(name, ".rds"))
  if (file.exists(f)) { message("    [cache] ", name); return(readRDS(f)) }
  message("    [run  ] ", name); t0 <- Sys.time()
  val <- force(expr)
  mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
  message("    [done ] ", name, " in ", mins, " min")
  attr(val, "fu_minutes") <- mins
  saveRDS(val, f)
  val
}
fu_minutes <- function(x) if (is.null(x) || is.null(attr(x, "fu_minutes"))) NA_real_ else attr(x, "fu_minutes")

## The models of this process when the work is split over several processes
## ("i/N": this one runs the i-th of every N, larger models first).
fu_my_share <- function(sizes, job) {
  j <- as.integer(strsplit(job, "/", fixed = TRUE)[[1]])
  ord <- order(-sizes)
  ord[(seq_along(ord) - 1L) %% j[2] == j[1] - 1L]
}
