## ===========================================================================
## SC-STEM fuel application, exploratory: leader and follower between the
## major brands and the independent pumps (pompe bianche) in the ten
## metropolitan cities. Reads the outputs of fuels-pretreatment.R.
##
##   setting 1  Y = a major-brand pump,   X = independents
##   setting 2  Y = an independent pump,  X = major-brand pumps
##   setting 3  Y = a major-brand pump,   X = pumps of the other major brands
##
## X is the nearest pump of the X set, or the mean of the X set within a
## radius: the data-driven radius r* (half range of the correlogram of weekly
## relative prices, per city) and the fixed radii 1, 2 and 5 km.
## For each pair: Pearson and Spearman correlations (lag 0 and cross-
## correlations) and Granger tests in both directions (Wald, heteroskedasticity-
## consistent), gasoline and diesel, daily (first differences) and weekly
## (levels), centred on the national or on the metropolitan mean. The tests are
## a SCREENING device, read against a reference pair (X drawn beyond 10 km):
## the formal inference of the application is the package bootstrap.
##
## This is part of the REPLICATION MATERIAL, not of the Stem package.
## ===========================================================================
suppressPackageStartupMessages(library(data.table))
invisible(Sys.setlocale("LC_TIME", "C"))
fu_here <- local({
  a <- commandArgs(trailingOnly = FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]), winslash = "/")) else
    normalizePath(getwd(), winslash = "/")
})
args <- commandArgs(trailingOnly = TRUE)
opt <- function(nm, def) { a <- grep(paste0("^--", nm, "="), args, value = TRUE)
  if (length(a)) sub("^--[^=]+=", "", a[1]) else def }
CFG <- list(
  pre      = opt("pre", file.path(fu_here, "fuels", "pretreatment")),
  out      = opt("out", file.path(fu_here, "fuels", "leader-follower")),
  majors   = c("Agip Eni", "Api-Ip", "Esso", "Q8", "Tamoil"),
  radii    = c(1, 2, 5),                                  # fixed radii, km
  p_daily  = 7L, p_weekly = 2L,                           # Granger lags
  l_daily  = 7L, l_weekly = 4L,                           # cross-correlation lags
  events   = as.Date(c("2022-03-22", "2022-12-01", "2023-01-01", "2025-05-15", "2026-01-01")),
  bins     = c(0, .5, 1, 1.5, 2, 3, 4, 5, 7.5, 10, 15, 20, 30),
  far_km   = 10,                                          # the reference pair: X beyond this distance
  seed     = 20261001L
)
OUT <- CFG$out; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
options(width = 170)
METRO <- c(RM = "Roma", MI = "Milano", "NA" = "Napoli", TO = "Torino", PA = "Palermo", GE = "Genova",
           BO = "Bologna", FI = "Firenze", BA = "Bari", CT = "Catania")
GCOL <- c("Agip Eni" = "#f2c500", "Api-Ip" = "#1f6f8b", "Esso" = "#c0392b", "Q8" = "#5b8c5a",
          "Tamoil" = "#e67e22", "independent" = "#4d4d4d", "other brands" = "#c9b3d9")
ser <- readRDS(file.path(CFG$pre, "metro_series.rds")); mn <- readRDS(file.path(CFG$pre, "means.rds"))
days <- ser[[1]]$days; ND <- length(days); sun <- which(format(days, "%u") == "7")
locf <- function(v) { ok <- !is.na(v); if (!any(ok)) return(v)
  idx <- cummax(ifelse(ok, seq_along(v), 0L)); idx[idx == 0L] <- which(ok)[1]; v[idx] }

## ---- the transformed series ------------------------------------------------
## weekly: the Sunday price (in force at the end of the week) centred on the
## national or on the metropolitan mean, then on the pump mean; daily: first
## differences of the raw, the national-centred and the metro-centred prices.
## Missing days carry the last reported price (the price in force).
TRANS <- c("weekly, national", "weekly, metro", "daily diff, raw", "daily diff, national", "daily diff, metro")
nat <- mn$national; setkey(nat, date)
series <- function(p, fu) {
  s <- ser[[p]]; mm <- mn$metro[prov == p]; setkey(mm, date)
  R <- t(apply(if (fu == "g") s$G else s$D, 1, locf))
  nref <- locf(nat[J(days)][[fu]]); mref <- locf(mm[J(days)][[fu]])
  Cn <- sweep(R, 2, nref); Cm <- sweep(R, 2, mref)
  wn <- Cn[, sun]; wm <- Cm[, sun]
  list("weekly, national" = wn - rowMeans(wn), "weekly, metro" = wm - rowMeans(wm),
       "daily diff, raw" = t(apply(R, 1, diff)), "daily diff, national" = t(apply(Cn, 1, diff)),
       "daily diff, metro" = t(apply(Cm, 1, diff)))
}
## fiscal pulses: the event week and the next (weekly); the event day and the
## six following, one dummy each (daily differences, which start on day 2)
pulses <- function(tr) {
  if (startsWith(tr, "weekly")) {
    w <- vapply(CFG$events, function(e) which(days[sun] >= e)[1], 1L)
    P <- matrix(0, length(sun), 2 * length(w)); for (k in seq_along(w)) P[cbind(w[k] + 0:1, 2 * k - 1:0)] <- 1
  } else {
    d <- vapply(CFG$events, function(e) match(e, days[-1]), 1L)
    P <- matrix(0, ND - 1, 7 * length(d)); for (k in seq_along(d)) P[cbind(d[k] + 0:6, 7 * (k - 1) + 1:7)] <- 1
  }
  P[, colSums(P) > 0, drop = FALSE]
}
PUL <- stats::setNames(lapply(TRANS, pulses), TRANS)

## ---- 1. the radius: correlogram of weekly relative prices ------------------
## pairwise Pearson correlation of the weekly metro-centred prices of all kept
## pumps against distance, binned; r* is the distance at which the excess of the
## binned median correlation over its plateau (pairs beyond 10 km) halves.
## An exponential fit r(d) = a + b exp(-d / phi) is reported for reference: one
## scale does not fit every city (Palermo decays on two)
cg <- list(); SER <- list(); DIST <- list()
for (p in names(METRO)) {
  s <- ser[[p]]; if (is.null(s)) next
  DIST[[p]] <- geodist::geodist(s$meta[, .(lon, lat)], measure = "geodesic") / 1000
  SER[[p]] <- list(g = series(p, "g"), d = series(p, "d"))
  up <- upper.tri(DIST[[p]])
  for (fu in c("g", "d")) {
    rr <- stats::cor(t(SER[[p]][[fu]][["weekly, metro"]]))
    z <- data.table(d = DIST[[p]][up], r = rr[up])[d <= max(CFG$bins)]
    z[, bin := cut(d, CFG$bins, include.lowest = TRUE)]
    cg[[length(cg) + 1]] <- z[, .(province = p, fuel = fu, n = .N, d = stats::median(d), r = stats::median(r)), by = bin]
  }
  message("series and correlogram: ", p)
}
cg <- rbindlist(cg)
rstar <- cg[, .(d = mean(d), r = mean(r), n = sum(n)), by = .(province, bin)][order(province, d)][, {
  a0 <- stats::median(r[d > 10]); b0 <- r[1] - a0; h <- a0 + b0 / 2
  ## r*: the first crossing of the half-excess level, interpolated between bins
  k <- which(r < h)[1]
  rs <- if (is.na(k)) NA_real_ else if (k == 1) d[1] else d[k - 1] + (d[k] - d[k - 1]) * (r[k - 1] - h) / (r[k - 1] - r[k])
  f <- tryCatch(stats::nls(r ~ a + b * exp(-d / phi), start = list(a = a0, b = b0, phi = 1), weights = sqrt(n),
                           control = stats::nls.control(warnOnly = TRUE)), error = function(e) NULL)
  cf <- if (is.null(f)) c(a = NA, b = NA, phi = NA) else stats::coef(f)
  .(plateau = a0, excess_0 = b0, r_star = round(rs, 1), exp_a = cf[["a"]], exp_b = cf[["b"]], exp_phi = cf[["phi"]])
}, by = province]
rstar[is.na(r_star), r_star := 1]
rstar[, city := METRO[province]]
cat("\nRADIUS: correlogram of weekly metro-centred prices (fuels averaged); r* = distance at which the excess\n",
    "correlation over the plateau (pairs beyond 10 km) halves; exponential fit for reference (km)\n")
print(rstar[, .(city, plateau = round(plateau, 3), excess_0 = round(excess_0, 3), r_star,
                exp_half_range = round(exp_phi * log(2), 1), exp_practical_range = round(3 * exp_phi, 1))])
fwrite(cg, file.path(OUT, "correlogram.csv")); fwrite(rstar, file.path(OUT, "radius.csv"))
grDevices::cairo_pdf(file.path(OUT, "correlogram.pdf"), width = 12, height = 5.5)
graphics::par(mfrow = c(2, 5), mar = c(3.4, 3.4, 2, 0.6), mgp = c(2.1, 0.6, 0))
for (p in names(METRO)) {
  z <- cg[province == p]; k <- rstar[province == p]
  graphics::plot(z$d, z$r, type = "n", log = "x", las = 1, bty = "n", ylim = c(-0.05, 0.55),
                 xlab = "distance (km, log scale)", ylab = "median correlation")
  graphics::points(z[fuel == "g"]$d, z[fuel == "g"]$r, pch = 16, col = "#1f6f8b")
  graphics::points(z[fuel == "d"]$d, z[fuel == "d"]$r, pch = 17, col = "#a8516e")
  dd <- exp(seq(log(0.2), log(30), length.out = 100))
  if (is.finite(k$exp_phi)) graphics::lines(dd, k$exp_a + k$exp_b * exp(-dd / k$exp_phi), col = "#4d4d4d")
  graphics::abline(v = k$r_star, lty = 2, col = "#8a939f"); graphics::abline(h = k$plateau, lty = 3, col = "#8a939f")
  graphics::title(sprintf("%s: r* = %.1f km", METRO[p], k$r_star), cex.main = 0.95)
}
invisible(grDevices::dev.off())

## ---- 2. the pairs: the three georeferenced datasets -------------------------
KINDS <- c("nearest", "r*", paste0(CFG$radii, " km"), "far (reference)")
pairs <- list(); sets <- list(); set.seed(CFG$seed)
for (p in names(METRO)) {
  m <- ser[[p]]$meta; D <- DIST[[p]]; rs <- rstar[province == p]$r_star
  maj <- which(m$brand_end %in% CFG$majors); ind <- which(m$group == "independent")
  defs <- list(`1` = list(Y = maj, X = function(i) ind),
               `2` = list(Y = ind, X = function(i) maj),
               `3` = list(Y = maj, X = function(i) maj[m$brand_end[maj] != m$brand_end[i]]))
  for (st in names(defs)) for (i in defs[[st]]$Y) {
    xs <- defs[[st]]$X(i); if (!length(xs)) next
    dx <- D[i, xs]; nn <- xs[which.min(dx)]
    rad <- c(rs, CFG$radii)
    ## the reference: one X of the same set drawn at random beyond far_km
    fx <- xs[dx > CFG$far_km]; far <- if (length(fx)) fx[sample.int(length(fx), 1)] else integer(0)
    mem <- lapply(rad, function(r) xs[dx <= r])
    sets[[length(sets) + 1]] <- data.table(setting = as.integer(st), province = p, row = i,
      kind = KINDS, radius_km = c(NA, rad, NA), members = c(list(nn), mem, list(far)))
    pairs[[length(pairs) + 1]] <- data.table(setting = as.integer(st), province = p, city_metro = METRO[p],
      y_site = m$site[i], y_brand = m$brand_end[i], y_group = m$group[i], y_type = m$type[i], y_municipality = m$city[i],
      y_lon = m$lon[i], y_lat = m$lat[i],
      nn_site = m$site[nn], nn_brand = m$brand_end[nn], nn_group = m$group[nn], nn_lon = m$lon[nn], nn_lat = m$lat[nn],
      nn_dist_km = round(min(dx), 3), r_star_km = rs,
      n_rstar = length(mem[[1]]), n_1km = length(mem[[2]]), n_2km = length(mem[[3]]), n_5km = length(mem[[4]]),
      members_rstar = paste(m$site[mem[[1]]], collapse = ";"), members_1km = paste(m$site[mem[[2]]], collapse = ";"),
      members_2km = paste(m$site[mem[[3]]], collapse = ";"), members_5km = paste(m$site[mem[[4]]], collapse = ";"),
      far_site = if (length(far)) m$site[far] else NA_character_, far_dist_km = if (length(far)) round(D[i, far], 3) else NA_real_)
  }
}
pairs <- rbindlist(pairs); sets <- rbindlist(sets)
sets[, n_x := lengths(members)]
for (st in 1:3) {
  fwrite(pairs[setting == st], file.path(OUT, sprintf("setting%d_pairs.csv", st)))
  saveRDS(pairs[setting == st], file.path(OUT, sprintf("setting%d_pairs.rds", st)))
}
cat("\nPAIRS: Y pumps per setting and city; nearest-X distance (median km); X pumps in the radius (median; share of Y with none)\n")
print(pairs[, .(Y = .N, nn_km = round(stats::median(nn_dist_km), 2), r_star = r_star_km[1],
                n_rstar = as.numeric(stats::median(n_rstar)), none_rstar = round(mean(n_rstar == 0), 2),
                n_1km = as.numeric(stats::median(n_1km)), none_1km = round(mean(n_1km == 0), 2),
                n_2km = as.numeric(stats::median(n_2km)), n_5km = as.numeric(stats::median(n_5km))), by = .(setting, city_metro)][order(setting, -Y)])

## ---- 3. correlations and Granger tests -------------------------------------
rowcor <- function(A, B) { A <- A - rowMeans(A); B <- B - rowMeans(B)
  rowSums(A * B) / sqrt(rowSums(A^2) * rowSums(B^2)) }
rowrank <- function(A) t(apply(A, 1, rank))
## cor(Y_t, X_{t-l}): l > 0, X leads; l < 0, Y leads
ccf_rows <- function(A, B, L) { T <- ncol(A)
  sapply(-L:L, function(l) if (l >= 0) rowcor(A[, (l + 1):T, drop = FALSE], B[, 1:(T - l), drop = FALSE]) else
    rowcor(A[, 1:(T + l), drop = FALSE], B[, (1 - l):T, drop = FALSE])) }
lagm <- function(v, p, ix) vapply(seq_len(p), function(l) v[ix - l], numeric(length(ix)))
## Wald test, heteroskedasticity-consistent (HC1), that the lags of the second
## series add nothing, both directions. The plain F test assumes a constant
## error variance, which the sticky daily changes violate: on Rome it rejected
## for 77% of the pairs more than 10 km apart, the robust test for 13%.
granger2 <- function(y, x, p, P) {
  ix <- (p + 1):length(y); Ly <- lagm(y, p, ix); Lx <- lagm(x, p, ix); Pm <- P[ix, , drop = FALSE]
  Pm <- Pm[, colSums(Pm != 0) > 0, drop = FALSE]
  one <- function(resp, own, other) {
    Z <- cbind(1, own, Pm, other); q <- qr(Z); keep <- sort(q$pivot[seq_len(q$rank)])
    j <- which(keep > ncol(Z) - ncol(other)); if (!length(j)) return(c(NA, NA))
    Z <- Z[, keep, drop = FALSE]; A <- solve(crossprod(Z)); b <- A %*% crossprod(Z, resp); e <- drop(resp - Z %*% b)
    V <- A %*% crossprod(Z * e) %*% A * nrow(Z) / (nrow(Z) - ncol(Z))
    W <- drop(crossprod(b[j], solve(V[j, j, drop = FALSE], b[j])))
    c(W, stats::pchisq(W, length(j), lower.tail = FALSE))
  }
  c(one(y[ix], Ly, Lx), one(x[ix], Lx, Ly))
}
## --reuse=1 reads the results of an earlier run and redoes only the summaries and maps
if (opt("reuse", "0") == "1" && file.exists(file.path(OUT, "results.rds"))) res <- readRDS(file.path(OUT, "results.rds")) else {
res <- list()
for (p in names(METRO)) {
  sp <- sets[province == p & n_x > 0]; if (!nrow(sp)) next
  nS <- nrow(ser[[p]]$meta)
  W <- matrix(0, nrow(sp), nS)
  W[cbind(rep(seq_len(nrow(sp)), sp$n_x), unlist(sp$members))] <- rep(1 / sp$n_x, sp$n_x)
  for (fu in c("g", "d")) for (tr in TRANS) {
    M <- SER[[p]][[fu]][[tr]]
    Y <- M[sp$row, , drop = FALSE]; X <- W %*% M
    L <- if (startsWith(tr, "weekly")) CFG$l_weekly else CFG$l_daily
    pl <- if (startsWith(tr, "weekly")) CFG$p_weekly else CFG$p_daily
    cc <- ccf_rows(Y, X, L); colnames(cc) <- paste0("ccf_", -L:L)
    gt <- t(vapply(seq_len(nrow(sp)), function(k) granger2(Y[k, ], X[k, ], pl, PUL[[tr]]), numeric(4)))
    res[[length(res) + 1]] <- data.table(sp[, .(setting, province, row, kind, radius_km, n_x)], fuel = fu, trans = tr,
      pearson = rowcor(Y, X), spearman = rowcor(rowrank(Y), rowrank(X)),
      W_xy = gt[, 1], p_xy = gt[, 2], W_yx = gt[, 3], p_yx = gt[, 4], cc)
  }
  message("tests: ", p)
}
res <- rbindlist(res, fill = TRUE)
res[, y_site := NA_character_]
for (p in unique(res$province)) res[province == p, y_site := ser[[p]]$meta$site[row]]
res[pairs, on = c("setting", "province", "y_site"), `:=`(y_brand = i.y_brand, nn_brand = i.nn_brand, nn_dist_km = i.nn_dist_km)]
res[, row := NULL]
saveRDS(res, file.path(OUT, "results.rds"))
fwrite(res, file.path(OUT, "results.csv"))
}

## ---- 4. summaries ----------------------------------------------------------
SET <- c("1" = "1: Y major, X independent", "2" = "2: Y independent, X major", "3" = "3: Y major, X other major")
FUEL <- c(g = "gasoline", d = "diesel")
summ <- function(z) z[, .(pairs = .N, pearson = round(stats::median(pearson, na.rm = TRUE), 3),
  spearman = round(stats::median(spearman, na.rm = TRUE), 3),
  `X->Y %` = round(100 * mean(p_xy < 0.05, na.rm = TRUE), 1), `Y->X %` = round(100 * mean(p_yx < 0.05, na.rm = TRUE), 1),
  `ccf X leads 1` = round(stats::median(ccf_1, na.rm = TRUE), 3), `ccf Y leads 1` = round(stats::median(`ccf_-1`, na.rm = TRUE), 3))]
cat("\nSUMMARY (medians over pairs; Granger: share of pairs rejecting at 5%, both directions)\n")
for (tr in TRANS) {
  cat(sprintf("\n== %s ==\n", tr))
  z <- res[trans == tr, summ(.SD), by = .(setting, kind, fuel)]
  z[, setting := SET[as.character(setting)]]; z[, fuel := FUEL[fuel]]
  print(z[order(fuel, setting, match(kind, KINDS))])
}
main <- res[kind %in% c("nearest", "r*", "far (reference)") & trans %in% c("weekly, metro", "daily diff, metro")]
cat("\nBY BRAND OF Y (settings 1 and 3) and BY BRAND OF THE NEAREST X (setting 2): metro-centred, both fuels\n")
by_y <- main[setting %in% c(1, 3), summ(.SD), by = .(setting, trans, kind, brand = y_brand)]
by_x <- main[setting == 2 & kind == "nearest", summ(.SD), by = .(setting, trans, kind, brand = nn_brand)]
print(rbind(by_y, by_x)[order(setting, trans, kind, brand)])
fwrite(rbind(by_y, by_x), file.path(OUT, "summary_by_brand.csv"))
cat("\nSETTING 3, nearest X: share of pairs where X Granger-causes Y (rows Y brand, columns X brand), weekly metro-centred, both fuels\n")
print(dcast(res[setting == 3 & kind == "nearest" & trans == "weekly, metro",
                .(v = round(100 * mean(p_xy < 0.05, na.rm = TRUE), 1)), by = .(y_brand, nn_brand)], y_brand ~ nn_brand, value.var = "v"))
cat("\nBY DISTANCE of the nearest X (metro-centred, both fuels)\n")
dz <- res[kind %in% c("nearest", "far (reference)") & trans %in% c("weekly, metro", "daily diff, metro")]
dz[, dist := fifelse(kind == "nearest", as.character(cut(nn_dist_km, c(0, .5, 1, 2, 5, 10, Inf), include.lowest = TRUE)),
                     "reference, > 10 km")]
print(dz[, summ(.SD), by = .(setting, trans, dist)][order(setting, trans, dist)])
cat("\nBY CITY: share rejecting X->Y / Y->X, nearest X and the far reference, weekly metro-centred, both fuels\n")
print(dcast(res[kind %in% c("nearest", "far (reference)") & trans == "weekly, metro",
  .(v = sprintf("%.0f/%.0f", 100 * mean(p_xy < .05, na.rm = TRUE), 100 * mean(p_yx < .05, na.rm = TRUE))),
  by = .(col = paste0("s", setting, ifelse(kind == "nearest", " nearest", " far")), city = METRO[province])],
  city ~ col, value.var = "v"))
all_s <- res[, summ(.SD), by = .(setting, kind, fuel, trans)]
fwrite(all_s, file.path(OUT, "summary.csv"))

## ---- 5. maps: the nearest X and the radius selection ------------------------
circ <- function(lon, lat, r) { th <- seq(0, 2 * pi, length.out = 90)
  cbind(lon + r / (111.32 * cos(lat * pi / 180)) * cos(th), lat + r / 110.57 * sin(th)) }
grDevices::cairo_pdf(file.path(OUT, "maps_pairs.pdf"), width = 15, height = 10.5, onefile = TRUE)
for (p in names(METRO)) {
  m <- ser[[p]]$meta; rs <- rstar[province == p]$r_star
  asp <- 1 / cos(mean(m$lat) * pi / 180)
  graphics::par(mfrow = c(2, 3), mar = c(0.6, 0.6, 2.4, 0.6), oma = c(0, 0, 2.2, 0))
  ## top row: every Y joined to its nearest X
  for (st in 1:3) {
    z <- pairs[setting == st & province == p]
    graphics::plot(m$lon, m$lat, pch = 16, cex = 0.45, col = "#b3bac3", asp = asp, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    if (st == 1) graphics::legend("bottomleft", bty = "n", pch = 16, col = c(GCOL[1:6], "#b3bac3"), cex = 0.8,
                                 legend = c(names(GCOL)[1:6], "other kept pumps"))
    if (nrow(z)) {
      graphics::segments(z$y_lon, z$y_lat, z$nn_lon, z$nn_lat, col = grDevices::adjustcolor("#4d4d4d", 0.6), lwd = 0.8)
      graphics::points(z$nn_lon, z$nn_lat, pch = 1, cex = 0.8, col = GCOL[z$nn_group])
      graphics::points(z$y_lon, z$y_lat, pch = 16, cex = 0.7, col = GCOL[ifelse(z$y_group == "independent", "independent", z$y_brand)])
    }
    graphics::title(sprintf("%s: Y (dot) to nearest X (circle), %d pairs, median %.2f km", SET[st], nrow(z),
                            if (nrow(z)) stats::median(z$nn_dist_km) else NA), cex.main = 0.85)
  }
  ## bottom row: zoom on the centre, four Y with their r* circle and members
  cx <- stats::median(m$lon); cy <- stats::median(m$lat); hw <- 4
  xl <- cx + c(-1, 1) * hw / (111.32 * cos(cy * pi / 180)); yl <- cy + c(-1, 1) * hw / 110.57
  for (st in 1:3) {
    z <- pairs[setting == st & province == p & n_rstar > 0 & y_lon > xl[1] & y_lon < xl[2] & y_lat > yl[1] & y_lat < yl[2]]
    graphics::plot(m$lon, m$lat, pch = 16, cex = 0.7, col = "#b3bac3", asp = asp, xlim = xl, ylim = yl,
                   xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    if (nrow(z)) {
      d0 <- (z$y_lon - cx)^2 + ((z$y_lat - cy) / asp)^2; pk <- which.min(d0)
      while (length(pk) < min(4, nrow(z))) {
        dm <- sapply(seq_len(nrow(z)), function(k) min((z$y_lon[k] - z$y_lon[pk])^2 + ((z$y_lat[k] - z$y_lat[pk]) / asp)^2))
        pk <- c(pk, which.max(dm))
      }
      for (k in pk) {
        ids <- strsplit(z$members_rstar[k], ";")[[1]]; mm <- m[match(ids, m$site)]
        graphics::polygon(circ(z$y_lon[k], z$y_lat[k], rs), border = "#8a939f", lty = 2)
        graphics::segments(z$y_lon[k], z$y_lat[k], mm$lon, mm$lat, col = grDevices::adjustcolor("#4d4d4d", 0.7))
        graphics::points(mm$lon, mm$lat, pch = 1, cex = 1.1, lwd = 1.3, col = GCOL[mm$group])
        graphics::points(z$y_lon[k], z$y_lat[k], pch = 16, cex = 1.3,
                         col = GCOL[ifelse(z$y_group[k] == "independent", "independent", z$y_brand[k])])
      }
    }
    graphics::title(sprintf("centre, 8 x 8 km: four Y, their r* = %.1f km circle and the X inside", rs), cex.main = 0.85)
  }
  graphics::mtext(sprintf("%s: Y filled, X open; grey, the other kept pumps", METRO[p]), outer = TRUE, font = 2)
}
invisible(grDevices::dev.off())
cat("\noutputs in", OUT, "\n")
