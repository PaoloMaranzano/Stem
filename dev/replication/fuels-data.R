## ===========================================================================
## SC-STEM fuel application: DATA MANAGEMENT
##
## From the station-level file as received to one data.frame per case, saved
## as data.RData in the folder of the case:
##
##   <root>/setting1-Ymajor-Xindependent/data.RData
##   <root>/setting2-Yindependent-Xmajor/data.RData
##   <root>/setting3-Ymajor-Xothermajor/data.RData
##
## each holding the object `fuels_data` (one row per city, pump Y and week; the
## columns are listed in fuels-functions.R, fu_case_frame()). The main script
## of each case reads it and writes its outputs in the same folder.
##
## Five stages, each saving its output so that a later run skips it (delete the
## file, or pass --redo=<stage>, to run it again):
##
##   1. pretreatment  the panel of pumps that set their prices actively, with
##                    the register of every step (unchanged from
##                    fuels-pretreatment.R: the numbers of the notes)
##                    -> <root>/pretreatment/
##   2. weekly        weekly relative prices: the Sunday price minus the
##                    national mean of the Sunday, minus the pump mean; a
##                    Sunday without a report stays MISSING (no last price
##                    carried forward)                   -> <root>/common/weekly.rds
##   3. completion    the missing weeks filled by the Kalman filter and
##                    smoother of a pooled STEM model of each city and fuel,
##                    used only to build the lags        -> <root>/common/completed.rds
##   4. pairs         the radius r* of each city (correlogram of the prices
##                    centred on the metropolitan mean) and, for each setting,
##                    the nearest X and the X within r* of every Y
##                                                       -> <root>/common/pairs.rds
##   5. cases         the data.frame of each case        -> <root>/<case>/data.RData
##
##     Rscript fuels-data.R
##     Rscript fuels-data.R --redo=completion,pairs,cases
##
## Stage 1 needs the station-level file (App_FuelsITA/station_level.zip beside
## this script, or --data=); stage 3 needs Stem, installed at the commit pinned
## in fuels-functions.R.
##
## This is part of the REPLICATION MATERIAL of the paper, not of the Stem
## package.
## ===========================================================================
suppressPackageStartupMessages(library(data.table))
invisible(Sys.setlocale("LC_TIME", "C"))

## This script's folder, and the auxiliary functions beside it.
HERE <- local({
  a <- commandArgs(trailingOnly = FALSE); m <- grep("^--file=", a, value = TRUE)
  f <- if (length(m)) normalizePath(sub("^--file=", "", m[1]), winslash = "/") else NULL
  if (is.null(f)) for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) { f <- normalizePath(of, winslash = "/"); break }
  }
  if (is.null(f)) normalizePath(getwd(), winslash = "/") else dirname(f)
})
source(file.path(HERE, "fuels-functions.R"))


## ===========================================================================
## THE SETUP
## ===========================================================================
CFG <- fu_config(list(
  data        = file.path(HERE, "App_FuelsITA", "station_level.zip"),
  root        = file.path(HERE, "fuels"),
  redo        = "",              # stages to run again even if their output exists
  ## the pre-treatment (decided; see the notes, Section "Data pre-treatment")
  from        = "2022-01-03",    # a Monday: the start of the window
  to          = "2026-06-28",    # a Sunday: its end
  ## the PRESAMPLE: the 52 weeks before the window, read only to give the lags
  ## of the first weeks of the window (lag 52 included). The filters, the
  ## register and the weeks the models are estimated on stay those of the
  ## window; the presample prices of the kept sites enter only as lags, and
  ## where a site did not report in 2021 they are filled by the Kalman
  ## completion of stage 3 (stage 2 counts how many).
  presample_from = "2021-01-04", # a Monday, 52 weeks before `from`
  metros      = names(FU_CITY),  # the eleven cities
  d_link      = 50,              # metres: successive codes closer than this are one site
  overlap     = 7L,              # days two successive codes may overlap
  link_gap    = 120L,            # days that may separate successive codes
  min_cover   = 0.95,            # both prices on at least 95% of the days
  max_gap     = 28L,             # no gap longer than 28 days
  min_chg_yr  = 12,              # at least 12 changes a year of each price
  max_spell   = 90L,             # no unchanged spell longer than 90 days
  max_prov_km = 150              # km: distance from the median location of the province
))
ROOT <- normalizePath(CFG$root[1], winslash = "/", mustWork = FALSE)
PRE <- file.path(ROOT, "pretreatment"); COM <- file.path(ROOT, "common")
dir.create(PRE, recursive = TRUE, showWarnings = FALSE)
dir.create(COM, recursive = TRUE, showWarnings = FALSE)
REDO <- trimws(CFG$redo)
run_stage <- function(stage, files) stage %in% REDO || !all(file.exists(files))
options(width = 160)


## ===========================================================================
## Stage 1. The pre-treatment
## ===========================================================================
## The code of fuels-pretreatment.R, unchanged: the steps and the register are
## those of the notes. Its outputs are the daily price matrices of the kept
## sites of the metropolitan cities (metro_series.rds, NA where a site did not
## report), the national and metropolitan daily means of all reporting
## non-highway pumps (means.rds) and the register.
if (run_stage("pretreatment", file.path(PRE, c("metro_series.rds", "means.rds", "pretreatment.rds")))) {
  message("\n== stage 1: pre-treatment")
  FROM <- as.Date(CFG$from[1]); TO <- as.Date(CFG$to[1]); PFROM <- as.Date(CFG$presample_from[1])
  ## DAYS: the window, on which every statistic and filter is computed;
  ## DAYS_ALL: presample and window, for the prices of the kept sites and the
  ## national and metropolitan means
  DAYS <- seq(FROM, TO, by = "day"); ND <- length(DAYS); YRS <- ND / 365.25
  DAYS_ALL <- seq(PFROM, TO, by = "day"); NDA <- length(DAYS_ALL)
  SH <- file.path(PRE, "shards")
  dir.create(SH, recursive = TRUE, showWarnings = FALSE)
  bgroup <- function(b) { b <- as.character(b); fifelse(b %in% FU_MAJOR, b, fifelse(b %in% "Pompe Bianche", "independent", "other brands")) }
  NORTH <- c("AO","TO","VC","NO","CN","AT","AL","BI","VB","GE","IM","SP","SV","MI","BG","BS","CO","CR","LC","LO",
             "MN","MB","PV","SO","VA","BZ","TN","VR","VI","BL","TV","VE","PD","RO","UD","GO","TS","PN","BO","FE",
             "FC","MO","PR","PC","RA","RE","RN")
  CENTRE <- c("FI","AR","GR","LI","LU","MS","PI","PT","PO","SI","PG","TR","AN","AP","FM","MC","PU","RM","FR","LT","RI","VT")
  area <- function(p) fifelse(p %in% NORTH, "North", fifelse(p %in% CENTRE, "Centre", fifelse(p == "unknown", "unknown", "South and islands")))
  mode_chr <- function(v) { v <- v[!is.na(v) & nzchar(v)]; if (!length(v)) NA_character_ else names(which.max(table(v))) }

  ## Pass 1: the codes over their whole life; the window rows, sharded by
  ## province; the national daily means of all reporting non-highway pumps
  parts <- if (dir.exists(CFG$data[1])) CFG$data[1] else {
    dz <- file.path(tempdir(), "fuels_station_level")
    if (!dir.exists(dz)) utils::unzip(CFG$data[1], exdir = dz)
    dz
  }
  parts <- list.files(parts, "^part_[0-9]+[.]RDS$", recursive = TRUE, full.names = TRUE)
  if (!length(parts)) stop("no station-level parts found in ", CFG$data[1], call. = FALSE)
  ids <- list(); nat <- list()
  for (k in seq_along(parts)) {
    x <- as.data.table(readRDS(parts[k]))[, .(id_pump, g = price_gasoline_self, d = price_diesel_self, date,
                                              brand, type = station_type, city, province,
                                              lat = latitude, lon = longitude)]
    ## the code of Naples is read as a missing value; a part whose province
    ## column is missing altogether is left unknown, and its codes leave at step 4
    x[, province := if (is.character(province)) fifelse(is.na(province), "NA", province) else rep("unknown", .N)]
    x[, type := fifelse(is.na(type), "", as.character(type))]
    x[, `:=`(brand = as.character(brand), city = as.character(city), lat = as.numeric(lat), lon = as.numeric(lon))]
    x[!is.finite(g) | g < 0.5 | g > 3.5, g := NA_real_]
    x[!is.finite(d) | d < 0.5 | d > 3.5, d := NA_real_]
    setorder(x, id_pump, date)
    s <- x[, .(first = min(date), last = max(date), lat = stats::median(lat, na.rm = TRUE),
               lon = stats::median(lon, na.rm = TRUE), brand_first = brand[1], brand_last = brand[.N],
               type = mode_chr(type), city = mode_chr(city), province = mode_chr(province)), by = id_pump]
    ids[[k]] <- s
    w <- x[date >= PFROM & date <= TO]                         # presample and window
    w[s, on = "id_pump", prov := i.province]
    nat[[k]] <- w[type != "Autostradale", .(sg = sum(g, na.rm = TRUE), ng = sum(!is.na(g)),
                                             sd = sum(d, na.rm = TRUE), nd = sum(!is.na(d))), by = .(prov, date)]
    for (pv in unique(w$prov)) {
      dir.create(file.path(SH, pv), showWarnings = FALSE)
      saveRDS(w[prov == pv, .(id_pump, date, g, d, brand, type)], file.path(SH, pv, sprintf("part_%02d.rds", k)))
    }
    message("pass 1: ", basename(parts[k]))
  }
  ids <- rbindlist(ids)
  nat <- rbindlist(nat)[, .(sg = sum(sg), ng = sum(ng), sd = sum(sd), nd = sum(nd)), by = .(prov, date)]
  means <- list(
    national = nat[, .(g = sum(sg) / sum(ng), d = sum(sd) / sum(nd), n = sum(ng)), by = date][order(date)],
    metro = nat[prov %in% CFG$metros, .(g = sg / ng, d = sd / nd, n = ng), by = .(prov, date)][order(prov, date)])
  saveRDS(means, file.path(PRE, "means.rds"))

  ## Pass 2, province by province: linking, statistics, filters, register
  link_pairs <- list(); reg <- list(); kept <- list(); series <- list()
  runmax <- function(v) { if (!any(v)) return(0L); r <- rle(v); max(r$lengths[r$values]) }
  for (pv in sort(unique(ids$province))) {
    fl <- list.files(file.path(SH, pv), full.names = TRUE)
    if (!length(fl)) next
    w <- rbindlist(lapply(fl, readRDS))
    ## the codes with data in the WINDOW (the presample does not enter the register)
    L <- ids[province == pv & id_pump %in% unique(w[date >= FROM]$id_pump)]
    A <- ids[province == pv & is.finite(lat) & is.finite(lon)]
    setorder(A, first)
    ## Step 2: link successive codes at the same point
    succ <- setNames(rep(NA_real_, nrow(A)), A$id_pump)
    if (nrow(A) > 1) {
      D <- geodist::geodist(A[, .(lon, lat)], measure = "geodesic")
      cand <- which(D <= 200, arr.ind = TRUE)
      cand <- cand[cand[, 1] != cand[, 2], , drop = FALSE]
      a <- cand[, 1]; b <- cand[, 2]
      ok <- A$first[b] > A$first[a] & A$first[b] >= A$last[a] - CFG$overlap[1] &
            A$first[b] <= A$last[a] + CFG$link_gap[1]
      cp <- data.table(a = A$id_pump[a[ok]], b = A$id_pump[b[ok]], dist = D[cbind(a[ok], b[ok])],
                       gap = as.integer(A$first[b[ok]] - A$last[a[ok]]), province = pv)
      link_pairs[[pv]] <- cp
      cp <- cp[dist <= CFG$d_link[1]][order(dist, abs(gap))]
      used_a <- used_b <- character(0)
      for (r in seq_len(nrow(cp))) {
        ka <- as.character(cp$a[r]); kb <- as.character(cp$b[r])
        if (ka %in% used_a || kb %in% used_b) next
        succ[ka] <- cp$b[r]; used_a <- c(used_a, ka); used_b <- c(used_b, kb)
      }
    }
    pred <- setNames(rep(NA_real_, length(succ)), names(succ))
    pred[as.character(stats::na.omit(succ))] <- as.numeric(names(succ)[!is.na(succ)])
    root <- vapply(names(succ), function(k) { while (!is.na(pred[k])) k <- as.character(pred[k]); k }, "")
    nocoord <- setdiff(as.character(unique(w$id_pump)), names(root))
    root <- c(root, stats::setNames(nocoord, nocoord))
    w[, site := root[as.character(id_pump)]]
    w[ids, on = "id_pump", idfirst := i.first]
    setorder(w, site, date, -idfirst)
    w <- w[, .SD[1], by = .(site, date)]                         # overlapping days: the newer code
    ## every statistic and filter on the window only: wa keeps the presample
    ## for the prices of the kept sites
    wa <- w; w <- wa[date >= FROM]
    last_code <- w[, .(last_id = id_pump[which.max(idfirst)], n_codes = uniqueN(id_pump),
                       brand_start = brand[1], brand_end = brand[.N], n_brands = uniqueN(stats::na.omit(brand)),
                       type = type[.N]), by = site]
    last_code[ids, on = c(last_id = "id_pump"), `:=`(lat = i.lat, lon = i.lon, city = i.city)]
    st <- w[, {
      di <- as.integer(date - FROM) + 1L
      gv <- dv <- rep(NA_real_, ND); gv[di] <- g; dv[di] <- d
      both <- !is.na(gv) & !is.na(dv)
      cg <- diff(gv); cd <- diff(dv)
      .(cover = mean(both), max_gap = runmax(!both),
        chg_g = sum(cg != 0, na.rm = TRUE) / YRS, chg_d = sum(cd != 0, na.rm = TRUE) / YRS,
        spell_g = runmax(c(FALSE, cg == 0) %in% TRUE) + 1L, spell_d = runmax(c(FALSE, cd == 0) %in% TRUE) + 1L)
    }, by = site]
    st <- last_code[st, on = "site"]
    st[, `:=`(province = pv, area = area(pv), group = bgroup(brand_end), metro = pv %in% CFG$metros)]
    st[, d_prov := NA_real_]
    okc <- is.finite(st$lat) & is.finite(st$lon)
    if (any(okc)) st[okc, d_prov := geodist::geodist(cbind(lon = lon, lat = lat), cbind(lon = stats::median(lon), lat = stats::median(lat)),
                                                     measure = "geodesic")[, 1] / 1000]
    steps <- list(
      "1 codes with data in the window" = NULL,
      "2 sites after linking codes"     = rep(TRUE, nrow(st)),
      "3 not on a highway"              = st$type != "Autostradale",
      "4 valid coordinates"             = is.finite(st$lat) & is.finite(st$lon) & st$lat > 35.4 & st$lat < 47.2 &
                                          st$lon > 6.5 & st$lon < 18.6 & st$d_prov <= CFG$max_prov_km[1],
      "5 coverage >= 95%"               = st$cover >= CFG$min_cover[1],
      "6 no gap > 28 days"              = st$max_gap <= CFG$max_gap[1],
      "7 >= 12 changes a year"          = st$chg_g >= CFG$min_chg_yr[1] & st$chg_d >= CFG$min_chg_yr[1],
      "8 no unchanged spell > 90 days"  = st$spell_g <= CFG$max_spell[1] & st$spell_d <= CFG$max_spell[1])
    reg[[paste(pv, 1)]] <- L[, .(step = names(steps)[1], province = pv, group = bgroup(brand_last),
                                 type = type, area = area(pv), n = 1L)]
    keep <- rep(TRUE, nrow(st))
    for (s in names(steps)[-1]) {
      keep <- keep & (steps[[s]] %in% TRUE)
      reg[[paste(pv, s)]] <- st[keep, .(step = s, province = pv, group, type, area, n = 1L)]
    }
    ## Step 9: one site per point
    k9 <- st[keep][order(-cover)][!duplicated(round(cbind(lat, lon), 6))]$site
    keep <- st$site %in% k9
    reg[[paste(pv, 9)]] <- st[keep, .(step = "9 one site per point", province = pv, group, type, area, n = 1L)]
    st[, kept := keep]
    kept[[pv]] <- st
    if (pv %in% CFG$metros) {
      ## the daily prices of the kept sites, presample and window
      ks <- st[kept == TRUE]$site
      wk <- wa[site %in% ks]
      di <- as.integer(wk$date - PFROM) + 1L
      G <- Dm <- matrix(NA_real_, length(ks), NDA, dimnames = list(ks, NULL))
      G[cbind(match(wk$site, ks), di)] <- wk$g; Dm[cbind(match(wk$site, ks), di)] <- wk$d
      series[[pv]] <- list(G = G, D = Dm, meta = st[kept == TRUE], days = DAYS_ALL, window_from = FROM)
    }
    message("pass 2: ", pv, " ", sum(keep), " kept of ", nrow(st), " sites")
  }
  reg <- rbindlist(reg)
  sites <- rbindlist(kept, fill = TRUE)
  links <- rbindlist(link_pairs)
  saveRDS(list(register = reg, sites = sites, links = links, cfg = CFG), file.path(PRE, "pretreatment.rds"))
  saveRDS(series, file.path(PRE, "metro_series.rds"))
  ## the register, written
  tot <- reg[, .(pumps = sum(n)), by = step][order(step)]
  tot[, removed := shift(pumps) - pumps]
  cat("\nREGISTER, Italy\n"); print(tot)
  fwrite(tot, file.path(PRE, "register_italy.csv"))
  comp <- function(by) {
    z <- reg[step %in% c("2 sites after linking codes", "9 one site per point"), .(n = sum(n)), by = c("step", by)]
    z[, share := round(100 * n / sum(n), 1), by = step]
    dcast(z, as.formula(paste(by, "~ step")), value.var = c("n", "share"), fill = 0)
  }
  fwrite(comp("group"), file.path(PRE, "composition_group.csv"))
  fwrite(comp("area"), file.path(PRE, "composition_area.csv"))
  mt <- dcast(reg[province %in% CFG$metros, .(pumps = sum(n)), by = .(province, step)], province ~ step, value.var = "pumps")
  cat("\nREGISTER, the metropolitan cities of the application\n"); print(mt)
  fwrite(mt, file.path(PRE, "register_metros.csv"))
  unlink(SH, recursive = TRUE)
}
ser <- readRDS(file.path(PRE, "metro_series.rds"))
mn  <- readRDS(file.path(PRE, "means.rds"))
CITIES <- intersect(names(FU_CITY), names(ser))
days <- ser[[1]]$days                                   # presample and window
if (is.null(ser[[1]]$window_from) || min(days) > as.Date(CFG$presample_from[1]))
  stop("the pre-treatment output has no presample: run again with --redo=pretreatment", call. = FALSE)
WEEKS <- days[fu_sundays(days)]
IN_WINDOW <- WEEKS >= as.Date(CFG$from[1])              # FALSE in the presample
PULSES <- fu_pulses(WEEKS)


## ===========================================================================
## Stage 2. Weekly relative prices
## ===========================================================================
## For each city and fuel, sites x weeks in cents per litre, presample and
## window: the Sunday price minus the NATIONAL mean of that Sunday, minus the
## pump mean over the window (decided: the national centring leaves to the
## latent process the movement common to the city). A Sunday without a report
## stays NA. The same prices centred on the METROPOLITAN mean are kept for the
## radius r* only.
F_WEEKLY <- file.path(COM, "weekly.rds")
if (run_stage("weekly", F_WEEKLY)) {
  message("\n== stage 2: weekly relative prices")
  nat <- mn$national[match(days, mn$national$date)]
  weekly <- lapply(stats::setNames(CITIES, CITIES), function(cc) {
    s <- ser[[cc]]
    mm <- mn$metro[prov == cc]; mm <- mm[match(days, mm$date)]
    list(meta = as.data.frame(s$meta),
         rel = list(g = fu_weekly_relative(s$G, days, nat$g, IN_WINDOW),
                    d = fu_weekly_relative(s$D, days, nat$d, IN_WINDOW)),
         rel_metro = list(g = fu_weekly_relative(s$G, days, mm$g, IN_WINDOW),
                          d = fu_weekly_relative(s$D, days, mm$d, IN_WINDOW)))
  })
  saveRDS(list(weekly = weekly, weeks = WEEKS, in_window = IN_WINDOW, pulses = PULSES), F_WEEKLY)
  ## the gaps of the window (the response) and of the presample (lags only):
  ## a site with little of 2021 gets most of its presample lags from the
  ## Kalman completion
  miss <- rbindlist(lapply(CITIES, function(cc) {
    R <- weekly[[cc]]$rel
    pre_obs <- rowMeans(!is.na(R$g[, !IN_WINDOW, drop = FALSE]) & !is.na(R$d[, !IN_WINDOW, drop = FALSE]))
    data.table(city = cc, sites = nrow(weekly[[cc]]$meta),
               window_missing_g_pct = round(100 * mean(is.na(R$g[, IN_WINDOW])), 3),
               window_missing_d_pct = round(100 * mean(is.na(R$d[, IN_WINDOW])), 3),
               presample_observed_median_pct = round(100 * stats::median(pre_obs), 1),
               sites_presample_below_50pct = sum(pre_obs < 0.5),
               sites_presample_none = sum(pre_obs == 0))
  }))
  cat("\nMISSING SUNDAY PRICES: % of pump-weeks in the window; presample coverage of the sites\n"); print(miss)
  fwrite(miss, file.path(COM, "weekly_missing.csv"))
}
WK <- readRDS(F_WEEKLY)


## ===========================================================================
## Stage 3. Kalman completion of the missing weeks
## ===========================================================================
## For each city and fuel, a pooled STEM model of the relative prices of all
## the kept sites (intercept and fiscal pulses, the latent AR(1) common to the
## city, the spatially correlated error), and STEM_Complete(): the conditional
## mean of every missing week given all the observed ones. The completed series
## build the lags of the models; the response keeps its gaps, which the Kalman
## filter of each model handles exactly.
F_COMP <- file.path(COM, "completed.rds")
if (run_stage("completion", F_COMP)) {
  message("\n== stage 3: Kalman completion")
  sha7 <- fu_require_stem()
  completed <- list(); csum <- list()
  for (cc in CITIES) {
    completed[[cc]] <- list()
    for (fu in c("g", "d")) {
      message(sprintf("  %s, %s: %d sites", FU_CITY[cc], FU_FUEL[[fu]], nrow(WK$weekly[[cc]]$meta)))
      res <- fu_complete(WK$weekly[[cc]]$rel[[fu]], WK$weekly[[cc]]$meta, WK$pulses)
      completed[[cc]][[fu]] <- res$completed
      csum[[length(csum) + 1]] <- cbind(city = cc, fuel = FU_FUEL[[fu]], res$summary)
    }
  }
  csum <- do.call(rbind, csum)
  cat("\nCOMPLETION: the pooled STEM fits that fill the missing weeks\n"); print(csum)
  saveRDS(list(completed = completed, summary = csum, stem = sha7), F_COMP)
  utils::write.csv(csum, file.path(COM, "completion_summary.csv"), row.names = FALSE)
}
CP <- readRDS(F_COMP)


## ===========================================================================
## Stage 4. The radius r* and the pairs
## ===========================================================================
## r* of each city: the correlogram of the weekly prices centred on the
## metropolitan mean (pairwise-complete weeks), both fuels; r* is the distance
## at which the excess correlation over the plateau (pairs beyond 10 km)
## halves. For each setting and each Y: the nearest X, and every X within r*.
F_PAIRS <- file.path(COM, "pairs.rds")
if (run_stage("pairs", F_PAIRS)) {
  message("\n== stage 4: radius and pairs")
  radius <- list(); pairs <- list(); cgs <- list()
  for (cc in CITIES) {
    meta <- WK$weekly[[cc]]$meta
    Dkm <- geodist::geodist(meta[, c("lon", "lat")], measure = "geodesic") / 1000
    ## on the weeks of the window, as decided
    cg <- rbind(cbind(fuel = "g", fu_correlogram(WK$weekly[[cc]]$rel_metro$g[, WK$in_window], Dkm)),
                cbind(fuel = "d", fu_correlogram(WK$weekly[[cc]]$rel_metro$d[, WK$in_window], Dkm)))
    rs <- fu_rstar(as.data.table(cg))
    radius[[cc]] <- data.frame(city = cc, city_name = FU_CITY[[cc]], r_star_km = rs$r_star,
                               plateau = rs$plateau, excess_0 = rs$excess_0)
    cgs[[cc]] <- cbind(city = cc, cg)
    pairs[[cc]] <- lapply(stats::setNames(1:3, 1:3), function(s) fu_pairs(meta, Dkm, rs$r_star, s))
  }
  radius <- do.call(rbind, radius)
  cat("\nRADIUS r* (km)\n"); print(radius, row.names = FALSE)
  saveRDS(list(radius = radius, pairs = pairs, correlogram = do.call(rbind, cgs)), F_PAIRS)
  utils::write.csv(radius, file.path(COM, "radius.csv"), row.names = FALSE)
  grDevices::cairo_pdf(file.path(COM, "correlogram.pdf"), width = 12, height = 2.75 * ceiling(length(CITIES) / 4))
  graphics::par(mfrow = c(ceiling(length(CITIES) / 4), 4), mar = c(3.4, 3.4, 2, 0.6), mgp = c(2.1, 0.6, 0))
  for (cc in CITIES) {
    z <- cgs[[cc]]
    graphics::plot(z$d, z$r, type = "n", log = "x", las = 1, bty = "n", ylim = c(-0.05, 0.55),
                   xlab = "distance (km, log scale)", ylab = "median correlation")
    graphics::points(z$d[z$fuel == "g"], z$r[z$fuel == "g"], pch = 16, col = "#1f6f8b")
    graphics::points(z$d[z$fuel == "d"], z$r[z$fuel == "d"], pch = 17, col = "#a8516e")
    graphics::abline(v = radius$r_star_km[radius$city == cc], lty = 2, col = "#8a939f")
    graphics::title(sprintf("%s: r* = %.1f km", FU_CITY[[cc]], radius$r_star_km[radius$city == cc]), cex.main = 0.95)
  }
  invisible(grDevices::dev.off())
}
PR <- readRDS(F_PAIRS)


## ===========================================================================
## Stage 5. The data.frame of each case
## ===========================================================================
## For each setting, one data.frame with every city: one row per city, pump Y
## and week (columns: fu_case_frame() in fuels-functions.R). Saved as
## `fuels_data` in <root>/<case>/data.RData, with a summary of the pairs.
for (s in 1:3) {
  cdir <- fu_case_dir(ROOT, s)
  F_CASE <- file.path(cdir, "data.RData")
  if (!run_stage("cases", F_CASE)) next
  message("\n== stage 5: case ", FU_CASES$folder[s])
  fuels_data <- do.call(rbind, lapply(CITIES, function(cc) {
    pr <- PR$pairs[[cc]][[as.character(s)]]
    if (is.null(pr$table) || !nrow(pr$table)) return(NULL)
    fu_case_frame(cc, s, WK$weekly[[cc]]$meta, WK$weekly[[cc]]$rel, CP$completed[[cc]],
                  WK$weeks, pr, WK$pulses, presample = !WK$in_window)
  }))
  rownames(fuels_data) <- NULL
  attr(fuels_data, "setting") <- FU_CASES[s, ]
  attr(fuels_data, "stem") <- CP$stem
  save(fuels_data, file = F_CASE)
  ## the pairs of the case, one row per Y
  site_level <- unique(fuels_data[, c("city", "city_name", "y_site", "y_brand", "y_group", "nn_site", "nn_brand",
                                      "nn_dist_km", "r_star_km", "n_rstar", "rs_fallback")])
  utils::write.csv(site_level, file.path(cdir, "pairs.csv"), row.names = FALSE)
  sm <- do.call(rbind, lapply(split(site_level, site_level$city), function(z) data.frame(
    city = z$city[1], city_name = z$city_name[1], Y_pumps = nrow(z),
    nn_km_median = round(stats::median(z$nn_dist_km), 2), r_star_km = z$r_star_km[1],
    n_rstar_median = stats::median(z$n_rstar), rstar_fallback_pct = round(100 * mean(z$rs_fallback), 1))))
  cat(sprintf("\nCASE %s: %d rows, %d pumps Y\n", FU_CASES$folder[s], nrow(fuels_data), nrow(site_level)))
  print(sm, row.names = FALSE)
  utils::write.csv(sm, file.path(cdir, "pairs_summary.csv"), row.names = FALSE)
}
cat("\ndone\n")
