## ===========================================================================
## APPLICATION, SCRIPT 1 OF 3: THE DATA OF THE MODELS AND THE EXPLORATORY
## ANALYSIS
##
## From the three files of replication/data to the data of the 66 models, and
## the exploratory analysis that describes the data and motivates the choice
## of the lags. The steps:
##
##   1. the data: Sunday prices of the stations of the eleven cities, the
##      national mean of every Sunday, and the pulses of the excise duties;
##   2. the relative prices: Sunday price minus the national mean, minus the
##      mean of the station over the window;
##   3. the completion of the missing weeks of the lagged series, by a pooled
##      STEM model of each city and fuel (the slow step: see below);
##   4. the pairs (Y, X) of the three settings and the data of every setting,
##      saved for 02-models.R;
##   5. the exploratory analysis: the stations, the missing weeks, unit roots,
##      persistence and stickiness of the prices, the distances between the
##      pairs, the memory of the relative prices (partial autocorrelations and
##      autoregressive orders), and the figure of the two steps of the
##      centering.
##
## Run from a terminal, in the folder of this script,
##
##     Rscript 01-data-and-eda.R
##
## or open it in RStudio and press Source. Options (command line):
##   --data=<folder>    the data (default: ../data)
##   --out=<folder>     the outputs (default: ../output/application)
##   --cities=RM,MI     a subset of the cities (default: all eleven); the
##                      exploratory tables then describe that subset only
##
## Outputs, in <out>:
##   data/cases.rds         the data of the three settings, read by 02-models.R
##   data/completed_*.rds   the completed series of every city and fuel (a cache:
##                          a later run reads them instead of fitting again)
##   eda/*.csv, eda/*.png   the tables and the figure listed in step 5, each with
##                          the place of the paper it reproduces
##
## Requirements: R (>= 4.2) with data.table, geodist, urca and Stem at the
## commit pinned in functions.R (installed from GitHub when missing). Time:
## step 3 fits 22 pooled STEM models, from a few minutes for the smaller cities
## to about an hour for Rome on one core; the rest takes a few minutes.
## ===========================================================================

suppressPackageStartupMessages(library(data.table))
invisible(Sys.setlocale("LC_TIME", "C"))
HERE <- local({
  a <- commandArgs(trailingOnly = FALSE); m <- grep("^--file=", a, value = TRUE)
  f <- if (length(m)) normalizePath(sub("^--file=", "", m[1]), winslash = "/") else NULL
  if (is.null(f)) for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) { f <- normalizePath(of, winslash = "/"); break }
  }
  if (is.null(f)) normalizePath(getwd(), winslash = "/") else dirname(f)
})
source(file.path(HERE, "functions.R"))

OPT <- fu_options(list(
  data   = file.path(HERE, "..", "data"),
  out    = file.path(HERE, "..", "output", "application"),
  cities = names(FU_CITY)
))
OUT_DATA <- file.path(OPT$out, "data"); OUT_EDA <- file.path(OPT$out, "eda")
dir.create(OUT_DATA, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_EDA, recursive = TRUE, showWarnings = FALSE)
CITIES <- intersect(names(FU_CITY), OPT$cities)
options(width = 160)


## ===========================================================================
## 1. The data
## ===========================================================================
## The Sundays run from January 2021 to June 2026: the 52 Sundays before
## FU_WINDOW_FROM are the presample, which provides the lags of the first weeks
## of the window only; every model is estimated on the weeks of the window.
DAT <- fu_read_data(OPT$data)
WEEKS <- DAT$national$week
IN_WINDOW <- WEEKS >= FU_WINDOW_FROM
PULSES <- fu_pulses(WEEKS)
NAT <- list(g = DAT$national$gasoline, d = DAT$national$diesel)
message(sprintf("%d stations in %d cities; %d Sundays, %d of them in the window (%s to %s)",
                sum(DAT$stations$city %in% CITIES), length(CITIES), length(WEEKS), sum(IN_WINDOW),
                min(WEEKS[IN_WINDOW]), max(WEEKS)))


## ===========================================================================
## 2. The relative prices
## ===========================================================================
## For every city: the stations (in the order of the station file), and, for
## each fuel, the Sunday prices and the relative prices, stations x weeks.
CITY <- lapply(stats::setNames(CITIES, CITIES), function(cc) {
  meta <- DAT$stations[DAT$stations$city == cc, ]
  P <- lapply(c(g = "g", d = "d"), function(fu) fu_price_matrix(DAT, cc, fu, WEEKS))
  list(meta = meta, price = P,
       rel = lapply(c(g = "g", d = "d"), function(fu) fu_relative(P[[fu]], NAT[[fu]], IN_WINDOW)))
})


## ===========================================================================
## 3. Completion of the missing weeks of the lagged series
## ===========================================================================
## A pooled STEM model of every city and fuel fills the missing weeks of the
## relative prices (fu_complete() in functions.R). The completed series give
## the lags of the models; the response keeps its missing weeks. Each fit is
## saved in <out>/data and read back on a later run.
SHA7 <- fu_require_stem()
COMP <- list(); comp_summary <- list()
for (cc in CITIES) {
  COMP[[cc]] <- list()
  for (fu in c("g", "d")) {
    message(sprintf("completion: %s, %s (%d stations)", FU_CITY[[cc]], FU_FUEL[[fu]], nrow(CITY[[cc]]$meta)))
    res <- fu_cached(OUT_DATA, sprintf("completed_%s_%s", cc, FU_FUEL[[fu]]),
                     fu_complete(CITY[[cc]]$rel[[fu]], CITY[[cc]]$meta, PULSES))
    COMP[[cc]][[fu]] <- res$completed
    comp_summary[[length(comp_summary) + 1]] <- cbind(city = cc, fuel = FU_FUEL[[fu]], res$summary,
                                                      minutes = fu_minutes(res))
  }
}
fu_write(do.call(rbind, comp_summary), file.path(OUT_EDA, "completion_fits.csv"))


## ===========================================================================
## 4. The pairs and the data of the three settings
## ===========================================================================
## For every setting and city, the nearest station X of every station Y
## (fu_pairs()), and one row per station Y and week with the response, the
## completed series and the pulses (fu_case_frame()). The three data.frames
## are saved together in <out>/data/cases.rds.
CASES <- lapply(stats::setNames(FU_SETTINGS$setting, FU_SETTINGS$setting), function(s) {
  d <- do.call(rbind, lapply(CITIES, function(cc) {
    pr <- fu_pairs(CITY[[cc]]$meta, s)
    if (is.null(pr)) return(NULL)
    fu_case_frame(cc, s, CITY[[cc]]$meta, CITY[[cc]]$rel, COMP[[cc]], WEEKS, pr, PULSES, !IN_WINDOW)
  }))
  rownames(d) <- NULL
  d
})
saveRDS(list(cases = CASES, weeks = WEEKS, in_window = IN_WINDOW, stem = SHA7),
        file.path(OUT_DATA, "cases.rds"))
message("wrote ", file.path(OUT_DATA, "cases.rds"))


## ===========================================================================
## 5. The exploratory analysis
## ===========================================================================

## --- 5.1 The stations ------------------------------------------------------
## Reproduces: the table of the stations by city in Section C1 of Supplement C
## ("The eleven cities"), and the totals of Section 5.1 of the paper (3,136
## stations, 264 independent; Rome the largest city, Messina the smallest).
st <- DAT$stations[DAT$stations$city %in% CITIES, ]
tab_st <- data.frame(city = FU_CITY_EN[FU_CITY[CITIES]], stations = as.integer(table(st$city)[CITIES]),
                     independent = as.integer(table(factor(st$city[st$group == "independent"], CITIES))))
tab_st <- tab_st[order(-tab_st$stations), ]
tab_st <- rbind(tab_st, data.frame(city = "all", stations = sum(tab_st$stations),
                                   independent = sum(tab_st$independent)))
fu_write(tab_st, file.path(OUT_EDA, "eda_stations_by_city.csv"))
print(tab_st, row.names = FALSE)

## --- 5.2 The window and the missing weeks ---------------------------------
## The share of the station-weeks without a Sunday price, in the window (the
## response) and over presample and window (the lagged series, which the
## completion fills).
## Reproduces: Section 5.1 of the paper and Section C1 of Supplement C (234
## weeks in the window); Section C2 of Supplement C (1.4 to 2.1% of the
## station-weeks of the lagged series are missing); Section C1 of Supplement C
## (the number of stations behind the national mean of a Sunday, some 20,000).
miss <- do.call(rbind, lapply(CITIES, function(cc) do.call(rbind, lapply(c("g", "d"), function(fu) data.frame(
  city = FU_CITY_EN[[FU_CITY[[cc]]]], fuel = FU_FUEL[[fu]],
  missing_window_pct = round(100 * mean(is.na(CITY[[cc]]$price[[fu]][, IN_WINDOW])), 3),
  missing_lagged_series_pct = round(100 * mean(is.na(CITY[[cc]]$price[[fu]])), 3))))))
fu_write(miss, file.path(OUT_EDA, "eda_missing_weeks.csv"))
window <- data.frame(first_sunday = min(WEEKS[IN_WINDOW]), last_sunday = max(WEEKS), weeks = sum(IN_WINDOW),
                     presample_weeks = sum(!IN_WINDOW),
                     national_mean_stations_median = stats::median(DAT$national$stations[IN_WINDOW]),
                     missing_lagged_min_pct = min(miss$missing_lagged_series_pct),
                     missing_lagged_max_pct = max(miss$missing_lagged_series_pct))
fu_write(window, file.path(OUT_EDA, "eda_window.csv"))
print(window, row.names = FALSE)

## --- 5.3 Unit roots, persistence and stickiness ---------------------------
## Station by station over the window (fu_station_stats()): the augmented
## Dickey-Fuller test on the relative price, its first-order autocorrelation,
## and the share of the weeks whose Sunday price is unchanged from the previous
## one. Summarized by city and fuel: share of the stations for which the test
## rejects the unit root at 5%, median autocorrelation, median share of
## unchanged weeks.
## Reproduces: Section 5.1 of the paper (the ADF test rejects for 96 to 100% of
## the stations of every city) and Section C1 of Supplement C (the same, the
## median first-order autocorrelation, and the stickiness of the Sunday price).
ss <- do.call(rbind, lapply(CITIES, function(cc) do.call(rbind, lapply(c("g", "d"), function(fu) {
  R <- CITY[[cc]]$rel[[fu]][, IN_WINDOW, drop = FALSE]; P <- CITY[[cc]]$price[[fu]][, IN_WINDOW, drop = FALSE]
  x <- do.call(rbind, lapply(seq_len(nrow(R)), function(i) fu_station_stats(R[i, ], P[i, ])))
  data.frame(city = FU_CITY_EN[[FU_CITY[[cc]]]], fuel = FU_FUEL[[fu]], stations = nrow(R),
             adf_reject_pct = round(100 * mean(x$adf_reject), 1), ar1_median = round(stats::median(x$ar1), 2),
             unchanged_weeks_median_pct = round(100 * stats::median(x$unchanged), 1))
}))))
fu_write(ss, file.path(OUT_EDA, "eda_unit_roots_and_stickiness.csv"))
print(ss, row.names = FALSE)

## --- 5.4 The pairs ---------------------------------------------------------
## For every setting and city: the stations Y, and the median distance to
## their station X; in setting 1, the share of the major-brand stations with an
## independent station within one kilometer.
## Reproduces: Section 5.1 of the paper and Section C1 of Supplement C (2,427
## stations Y in settings 1 and 3, 264 in setting 2, from 6 in Florence to 71
## in Naples; median distance to X: 0.6 to 1.6 km in setting 3, 1.4 km in
## Naples and 2.8 to 7.9 km elsewhere in setting 1; 12% of the major-brand
## stations have an independent within about one kilometer).
pairs <- do.call(rbind, lapply(FU_SETTINGS$setting, function(s) {
  d <- unique(CASES[[as.character(s)]][, c("city", "y_site", "x_dist_km")])
  do.call(rbind, lapply(split(d, d$city), function(z) data.frame(
    setting = s, city = FU_CITY_EN[[FU_CITY[[z$city[1]]]]], stations_Y = nrow(z),
    x_dist_km_median = round(stats::median(z$x_dist_km), 2),
    x_within_1km_pct = round(100 * mean(z$x_dist_km <= 1), 1))))
}))
pairs <- pairs[order(pairs$setting, match(pairs$city, FU_CITY_EN[FU_CITY_ORDER])), ]
fu_write(pairs, file.path(OUT_EDA, "eda_pairs.csv"))
print(stats::aggregate(cbind(stations_Y) ~ setting, pairs, sum), row.names = FALSE)
s1 <- unique(CASES[["1"]][, c("y_site", "x_dist_km")])
message(sprintf("setting 1: %.0f%% of the major-brand stations have an independent within 1 km",
                100 * mean(s1$x_dist_km <= 1)))

## --- 5.5 The memory of the relative prices: the choice of the lags --------
## For every setting, city and fuel, on the completed series of the stations Y
## (fu_memory()): the share of the stations whose partial autocorrelation lies
## outside the band +-2/sqrt(T) at lags 1 to 8 and 52, and the order of the
## autoregression chosen by the AIC. The shares are summarized by their median
## over the cities, for every setting and fuel; the orders by their quartiles
## over the stations Y of every setting and fuel.
## Reproduces: Section 5.2 of the paper and Section C2 of Supplement C (the
## partial autocorrelation exceeds the band for all the stations at lag 1, for
## 56-72% at lag 2, 26-41% at lag 3 and 14-22% at lag 4, at most 17% from lag
## 5 on; it is within the band at lag 52 for 99-100% of the stations; the
## autoregressive order chosen by the AIC has median 3 to 4 and quartiles 2
## and 5 to 6). These shares led to four lags for every block.
mem <- lapply(FU_SETTINGS$setting, function(s) {
  d <- CASES[[as.character(s)]]
  lapply(unique(d$city), function(cc) lapply(c("g", "d"), function(fu) {
    m <- fu_memory(d[d$city == cc, ], fu)
    list(pacf = cbind(setting = s, m$pacf), ar = cbind(setting = s, m$ar_order))
  }))
})
mem <- unlist(unlist(mem, recursive = FALSE), recursive = FALSE)
pacf <- do.call(rbind, lapply(mem, `[[`, "pacf"))
ar_order <- do.call(rbind, lapply(mem, `[[`, "ar"))
pacf_med <- stats::aggregate(share_outside ~ setting + fuel + lag, pacf, stats::median)
pacf_med$share_outside <- round(100 * pacf_med$share_outside, 1)
pacf_wide <- stats::reshape(pacf_med, idvar = c("setting", "fuel"), timevar = "lag", direction = "wide")
names(pacf_wide) <- sub("^share_outside[.]", "lag_", names(pacf_wide))
fu_write(pacf_wide, file.path(OUT_EDA, "eda_pacf_share_outside_band.csv"))
cat("\nPARTIAL AUTOCORRELATION OUTSIDE THE BAND, % of the stations Y (median over the cities)\n")
print(pacf_wide, row.names = FALSE)
ar_q <- stats::aggregate(order ~ setting + fuel, ar_order, function(v) stats::quantile(v, c(0.25, 0.5, 0.75)))
ar_q <- do.call(data.frame, ar_q)
names(ar_q) <- c("setting", "fuel", "order_q25", "order_median", "order_q75")
fu_write(ar_q, file.path(OUT_EDA, "eda_ar_order.csv"))
cat("\nORDER OF THE AUTOREGRESSION CHOSEN BY THE AIC, quartiles over the stations Y\n")
print(ar_q, row.names = FALSE)

## --- 5.6 Figure C1 of Appendix C: from posted prices to relative prices ---
## Two gasoline stations of Rome, an Eni station and an independent one, with
## ordinary dynamics (largest absolute relative price below the median of
## their group) and premiums apart (the Eni station at the ninth decile of the
## premiums of its group, the independent at the first). Three panels: the
## Sunday price with the national mean; the price net of the national mean,
## with the mean of each station (its premium); the relative price.
## Reproduces: Figure C1 of Appendix C (fig:app:centring), and the premiums of
## the two stations quoted in its text (0.8 cents above the national price and
## 2.7 cents below it).
if ("RM" %in% CITIES) {
  P <- CITY$RM$price$g; meta <- CITY$RM$meta
  C <- 100 * sweep(P, 2, NAT$g)                                # cents above the national mean
  Z <- CITY$RM$rel$g                                           # the relative price
  full <- rowSums(is.na(C[, IN_WINDOW])) == 0
  pick <- function(brand, q) {
    i <- which(meta$brand == brand & full)
    dev <- apply(abs(Z[i, IN_WINDOW]), 1, max); i <- i[dev <= stats::median(dev)]
    prem <- rowMeans(C[i, IN_WINDOW]); i[which.min(abs(prem - stats::quantile(prem, q)))]
  }
  two <- c(pick("Agip Eni", 0.9), pick("Pompe Bianche", 0.1))
  col2 <- c("#1f6f8b", "#a8516e"); lab <- c("Eni station", "independent station")
  w <- which(IN_WINDOW); days <- WEEKS
  f <- file.path(OUT_EDA, "fig_app_centring.png")
  grDevices::png(f, width = 9, height = 7.6, units = "in", res = 300,
                 type = if (isTRUE(capabilities("cairo"))) "cairo" else "windows")
  graphics::par(mfrow = c(3, 1), mar = c(2.4, 4.6, 2.4, 1), mgp = c(2.6, 0.7, 0))
  panel <- function(Y, ylab, title, extra = NULL) {
    plot(days[w], Y[1, w], type = "n", ylim = range(Y[, w], extra, na.rm = TRUE), xlab = "", ylab = ylab,
         main = title, cex.main = 1.05, font.main = 1, las = 1)
    graphics::abline(v = FU_EVENTS$date, lty = 3, col = "grey60")
    for (j in 1:2) graphics::lines(days[w], Y[j, w], col = col2[j], lwd = 1.3)
  }
  panel(P[two, , drop = FALSE], "euro per liter", "(a) Sunday price of gasoline", NAT$g[w])
  graphics::lines(days[w], NAT$g[w], col = "grey20", lty = 2, lwd = 1.2)
  graphics::legend("topright", bty = "n", lty = c(1, 1, 2), col = c(col2, "grey20"), lwd = 1.3,
                   legend = c(lab, "national mean"), cex = 1)
  panel(C[two, , drop = FALSE], "cents per liter", "(b) Net of the national mean of the week")
  graphics::abline(h = rowMeans(C[two, IN_WINDOW]), col = col2, lty = 2)
  panel(Z[two, , drop = FALSE], "cents per liter", "(c) Net of the mean of the station: the relative price of the models")
  graphics::abline(h = 0, lty = 2)
  invisible(grDevices::dev.off()); message("wrote ", f)
  fu_write(data.frame(site = meta$site[two], brand = meta$brand[two], premium_cents = round(rowMeans(C[two, IN_WINDOW]), 2)),
           file.path(OUT_EDA, "fig_app_centring_stations.csv"))
}
message("done: the data of the models in ", OUT_DATA, ", the exploratory analysis in ", OUT_EDA)
