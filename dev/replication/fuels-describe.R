## ===========================================================================
## SC-STEM fuel application: the kept pumps. Where they are, how much they
## still stay still, their daily and weekly series and the properties of the
## centred series (ACF, unit-root and stationarity tests). Reads the outputs of
## fuels-pretreatment.R.
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
PRE <- opt("pre", file.path(fu_here, "fuels", "pretreatment"))
OUT <- opt("out", file.path(fu_here, "fuels", "describe")); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
pt <- readRDS(file.path(PRE, "pretreatment.rds")); ser <- readRDS(file.path(PRE, "metro_series.rds"))
mn <- readRDS(file.path(PRE, "means.rds"))
sites <- pt$sites[kept == TRUE]
days <- ser[[1]]$days; ND <- length(days); sun <- which(format(days, "%u") == "7")
METRO <- c(RM = "Roma", MI = "Milano", "NA" = "Napoli", TO = "Torino", PA = "Palermo", GE = "Genova",
           BO = "Bologna", FI = "Firenze", BA = "Bari", CT = "Catania")
GCOL <- c("Agip Eni" = "#f2c500", "Api-Ip" = "#1f6f8b", "Esso" = "#c0392b", "Q8" = "#5b8c5a",
          "Tamoil" = "#e67e22", "independent" = "#7f8c8d", "other brands" = "#8e44ad")
options(width = 160)
locf <- function(v) { ok <- !is.na(v); if (!any(ok)) return(v)
  idx <- cummax(ifelse(ok, seq_along(v), 0L)); idx[idx == 0L] <- which(ok)[1]; v[idx] }

## ---- 1. how many, and where ------------------------------------------------
cat("KEPT SITES, Italy:", nrow(sites), "\n")
print(sites[, .N, by = area][order(-N)])
tabm <- dcast(sites[province %in% names(METRO), .N, by = .(province, group)], province ~ group, value.var = "N", fill = 0)
tabm[, total := rowSums(.SD), .SDcols = -1]; tabm[, city := METRO[province]]
cat("\nKEPT SITES in the ten metropolitan cities, by brand group\n"); print(tabm[order(-total)])
fwrite(tabm, file.path(OUT, "kept_metros_by_group.csv"))
ita <- rnaturalearth::ne_countries(country = "Italy", scale = "medium", returnclass = "sf")
draw_map <- function() {
  graphics::layout(cbind(1, matrix(2:11, 5, byrow = TRUE)), widths = c(1.35, 0.6, 0.6))
  graphics::par(mar = c(0.5, 0.5, 2, 0.5))
  plot(sf::st_geometry(ita), col = "#f4f5f7", border = "#b8bec6", lwd = 0.6, xlim = c(6.6, 18.5), ylim = c(36.6, 47.1))
  graphics::points(sites$lon, sites$lat, pch = 16, cex = 0.2, col = GCOL[sites$group])
  graphics::title(sprintf("Kept sites, Italy (%d)", nrow(sites)), cex.main = 1)
  graphics::legend("bottomleft", bty = "n", pch = 16, col = GCOL, legend = names(GCOL), cex = 0.8, pt.cex = 1.3)
  for (p in names(METRO)) {
    z <- sites[province == p]
    graphics::plot(z$lon, z$lat, pch = 16, cex = 0.55, col = GCOL[z$group], xaxt = "n", yaxt = "n",
                   xlab = "", ylab = "", asp = 1 / cos(mean(z$lat) * pi / 180))
    graphics::title(sprintf("%s (%d)", METRO[p], nrow(z)), cex.main = 0.9, line = 0.3)
  }
}
grDevices::png(file.path(OUT, "map_kept.png"), width = 2200, height = 2000, res = 170); draw_map(); invisible(grDevices::dev.off())
grDevices::cairo_pdf(file.path(OUT, "map_kept.pdf"), width = 13, height = 12); draw_map(); invisible(grDevices::dev.off())

## ---- 2. do the kept sites still stay still? --------------------------------
sites[, `:=`(chg_min = pmin(chg_g, chg_d), spell_max = pmax(spell_g, spell_d))]
cat("\nRESIDUAL STILLNESS among kept sites\n")
cat("changes a year, less active fuel (quantiles 5-25-50-75-95):",
    round(stats::quantile(sites$chg_min, c(.05, .25, .5, .75, .95)), 1), "\n")
cat("longest unchanged spell, worse fuel (days):", stats::quantile(sites$spell_max, c(.05, .25, .5, .75, .95)), "\n")
still <- sites[, .(sites = .N, `< 24 changes/yr` = sum(chg_min < 24), `spell > 60 days` = sum(spell_max > 60),
                   both = sum(chg_min < 24 & spell_max > 60)), by = group]
print(still)

## ---- 3. the properties of the centred series, metropolitan sites ------------
nat <- mn$national; setkey(nat, date)
feat <- list()
for (p in names(METRO)) {
  s <- ser[[p]]; if (is.null(s)) next
  mm <- mn$metro[prov == p]; setkey(mm, date)
  for (fu in c("g", "d")) {
    X <- if (fu == "g") s$G else s$D
    nref <- nat[J(days)][[fu]]; mref <- mm[J(days)][[fu]]
    for (i in seq_len(nrow(X))) {
      v <- locf(X[i, ])
      for (ctr in c("national", "metro")) {
        cv <- v - (if (ctr == "national") nref else mref)
        cw <- cv[sun]; cw <- cw - mean(cw, na.rm = TRUE); cw <- cw[is.finite(cw)]
        dd <- diff(cv); dd <- dd[is.finite(dd)]
        adf <- urca::ur.df(cw, type = "drift", lags = 8, selectlags = "BIC")
        kp <- urca::ur.kpss(cw, type = "mu", lags = "short"); kl <- urca::ur.kpss(cw, type = "mu", lags = "long")
        a <- stats::acf(cw, lag.max = 13, plot = FALSE)$acf
        ad <- stats::acf(dd, lag.max = 7, plot = FALSE)$acf
        feat[[length(feat) + 1]] <- data.table(province = p, site = rownames(X)[i], fuel = fu, centre = ctr,
          adf_reject = adf@teststat[1] < adf@cval[1, "5pct"], kpss_reject = kp@teststat > kp@cval[, "5pct"],
          kpss_long_reject = kl@teststat > kl@cval[, "5pct"],
          ar1 = a[2], acf4 = a[5], acf13 = a[14], dacf1 = ad[2], dacf7 = ad[8],
          still_weeks = mean(diff(v[sun]) == 0, na.rm = TRUE))
      }
    }
  }
  message("features: ", p)
}
feat <- rbindlist(feat)
fwrite(feat, file.path(OUT, "features_metro_sites.csv"))
sm <- feat[, .(sites = .N, `ADF rejects unit root %` = round(100 * mean(adf_reject), 1),
               `KPSS rejects %` = round(100 * mean(kpss_reject), 1),
               `KPSS long lags rejects %` = round(100 * mean(kpss_long_reject), 1),
               `AR(1) weekly` = round(stats::median(ar1), 2), `ACF 4 wk` = round(stats::median(acf4), 2),
               `ACF 13 wk` = round(stats::median(acf13), 2), `daily diff ACF 1` = round(stats::median(dacf1), 2),
               `daily diff ACF 7` = round(stats::median(dacf7), 2)),
           by = .(centre, fuel, province)]
sm[, city := METRO[province]]
cat("\nPROPERTIES of the centred series: tests and ACF on weekly levels, ACF of daily differences; medians over sites\n")
for (ce in c("national", "metro")) for (fu in c("g", "d")) {
  cat(sprintf("\n-- centred on the %s mean, %s\n", ce, c(g = "gasoline", d = "diesel")[fu]))
  print(sm[centre == ce & fuel == fu, -(1:2)][order(province)])
}
cat("\nweeks with an unchanged Sunday price (median over sites, %)\n")
print(dcast(feat[centre == "metro", .(w = round(100 * stats::median(still_weeks), 1)), by = .(province, fuel)],
            province ~ fuel, value.var = "w"))
fwrite(sm, file.path(OUT, "features_summary.csv"))

## ---- 4. the series of two kept sites per city: a major brand and an independent
pick <- function(z, grp) {
  z <- copy(z)[, chg_min := pmin(chg_g, chg_d)]
  z <- z[if (grp == "major") group %in% names(GCOL)[1:5] else group == "independent"]
  if (!nrow(z)) return(NULL)
  z[order(abs(chg_min - stats::median(chg_min)))][1]
}
zoom <- which(days >= as.Date("2025-04-01") & days <= as.Date("2025-06-29"))
grDevices::cairo_pdf(file.path(OUT, "series_metros.pdf"), width = 12, height = 9, onefile = TRUE)
for (p in names(METRO)) {
  s <- ser[[p]]; if (is.null(s)) next
  mm <- mn$metro[prov == p]; setkey(mm, date)
  graphics::par(mfrow = c(2, 3), mar = c(3, 4, 2.6, 0.8), oma = c(0, 0, 2, 0), mgp = c(2.2, 0.6, 0))
  for (grp in c("major", "independent")) {
    z <- pick(s$meta, grp)
    if (is.null(z)) { for (k in 1:3) graphics::plot.new(); next }
    i <- match(z$site, rownames(s$G))
    g <- locf(s$G[i, ]); d <- locf(s$D[i, ])
    graphics::plot(days[zoom], g[zoom], type = "s", col = "#1f6f8b", las = 1, bty = "n", xlab = "", ylab = "euro per litre",
                   ylim = range(c(g[zoom], d[zoom], mm[J(days[zoom])]$g, mm[J(days[zoom])]$d), na.rm = TRUE))
    graphics::lines(days[zoom], d[zoom], type = "s", col = "#a8516e")
    graphics::lines(days[zoom], mm[J(days[zoom])]$g, col = "#1f6f8b", lty = 3)
    graphics::lines(days[zoom], mm[J(days[zoom])]$d, col = "#a8516e", lty = 3)
    graphics::title(sprintf("%s, %s: daily, Apr-Jun 2025", z$brand_end, tools::toTitleCase(tolower(z$city))), cex.main = 0.85)
    graphics::plot(days[sun], g[sun], type = "s", col = "#1f6f8b", las = 1, bty = "n", xlab = "", ylab = "euro per litre",
                   ylim = range(c(g[sun], d[sun]), na.rm = TRUE))
    graphics::lines(days[sun], d[sun], type = "s", col = "#a8516e")
    graphics::title("weekly, Sunday price", cex.main = 0.85)
    cg <- g[sun] - mm[J(days[sun])]$g; cd <- d[sun] - mm[J(days[sun])]$d
    cn <- g[sun] - nat[J(days[sun])]$g
    graphics::plot(days[sun], cg - mean(cg, na.rm = TRUE), type = "l", col = "#1f6f8b", las = 1, bty = "n",
                   xlab = "", ylab = "centred", ylim = range(c(cg - mean(cg, na.rm = TRUE), cd - mean(cd, na.rm = TRUE)), na.rm = TRUE))
    graphics::lines(days[sun], cd - mean(cd, na.rm = TRUE), col = "#a8516e")
    graphics::lines(days[sun], cn - mean(cn, na.rm = TRUE), col = "#1f6f8b", lty = 3)
    graphics::abline(h = 0, col = "#8a939f", lty = 2)
    graphics::title("weekly centred on the metro mean (dotted gasoline: national)", cex.main = 0.85)
  }
  graphics::mtext(sprintf("%s: a major-brand site (top) and an independent (bottom); blue gasoline, red diesel, dotted metro mean",
                          METRO[p]), outer = TRUE, font = 2, cex = 0.9)
}
invisible(grDevices::dev.off())
cat("\noutputs in", OUT, "\n")
