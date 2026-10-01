## ===========================================================================
## SC-STEM fuel application: DATA PRE-TREATMENT.
##
## From the station-level file as received to the panel of pumps that set
## their prices actively, with a register of every step: how many pumps each
## step removes, and the composition of the sample before and after, for the
## section "Data pre-treatment" of the paper.
##
##   Step 1  the window
##   Step 2  linking: a station that changes operator is registered under a new
##           code; successive codes at the same point are joined into one site
##   Step 3  sites not on highways
##   Step 4  valid coordinates, within 150 km of the province median location
##   Step 5  coverage: both self-service prices on at least 95% of the days
##   Step 6  no gap longer than 28 days
##   Step 7  at least 12 price changes a year on each fuel
##   Step 8  no unchanged spell longer than 90 days on either fuel
##   Step 9  one site per point
##
## Outputs in <out>: the register (national, by province, by composition), the
## link diagnostics, the list of kept sites with their statistics, the national
## and metropolitan daily means of all reporting pumps (for centring), and the
## daily series of the kept sites of the metropolitan cities of the application.
##
##     Rscript fuels-pretreatment.R
##     Rscript fuels-pretreatment.R --from=2021-01-04 --max_spell=120
##
## This is part of the REPLICATION MATERIAL, not of the Stem package.
## ===========================================================================
suppressPackageStartupMessages(library(data.table))

fu_here <- local({
  a <- commandArgs(trailingOnly = FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]), winslash = "/")) else
    normalizePath(getwd(), winslash = "/")
})
fu_config <- function(defaults, args = commandArgs(trailingOnly = TRUE)) {
  for (a in grep("^--[^=]+=.", args, value = TRUE)) {
    nm <- sub("^--([^=]+)=.*$", "\\1", a)
    if (!nm %in% names(defaults)) stop("unknown option --", nm, call. = FALSE)
    parts <- trimws(strsplit(sub("^--[^=]+=", "", a), ",", fixed = TRUE)[[1]])
    d <- defaults[[nm]]
    defaults[[nm]] <- if (is.integer(d)) as.integer(parts) else if (is.numeric(d)) as.numeric(parts) else parts
  }
  defaults
}

## ===========================================================================
## THE SETUP
## ===========================================================================
CFG <- fu_config(list(
  data       = file.path(fu_here, "App_FuelsITA", "station_level.zip"),
  from       = "2022-01-03",     # a Monday: after the change of pricing behaviour of 2022
  to         = "2026-06-28",     # a Sunday
  ## the metropolitan cities of the application: those with at least 10 kept
  ## independents (Naples, Rome, Catania, Bologna, Turin, Milan, Palermo, Messina,
  ## Bari), plus Florence and Venice by choice
  metros     = c("RM", "MI", "NA", "TO", "PA", "BA", "CT", "BO", "ME", "FI", "VE"),
  d_link     = 50,               # metres: successive codes closer than this are one site
  overlap    = 7L,               # days two successive codes may overlap
  link_gap   = 120L,             # days that may separate successive codes
  min_cover  = 0.95,
  max_gap    = 28L,
  min_chg_yr = 12,
  max_spell  = 90L,
  max_prov_km = 150,             # km: step 4, distance from the median location of the province
  out        = file.path(fu_here, "fuels", "pretreatment")
))
FROM <- as.Date(CFG$from[1]); TO <- as.Date(CFG$to[1])
DAYS <- seq(FROM, TO, by = "day"); ND <- length(DAYS); YRS <- ND / 365.25
OUT <- CFG$out[1]; SH <- file.path(OUT, "shards")
dir.create(SH, recursive = TRUE, showWarnings = FALSE)
MAJOR <- c("Agip Eni", "Api-Ip", "Esso", "Q8", "Tamoil")
bgroup <- function(b) { b <- as.character(b); fifelse(b %in% MAJOR, b, fifelse(b %in% "Pompe Bianche", "independent", "other brands")) }
NORTH <- c("AO","TO","VC","NO","CN","AT","AL","BI","VB","GE","IM","SP","SV","MI","BG","BS","CO","CR","LC","LO",
           "MN","MB","PV","SO","VA","BZ","TN","VR","VI","BL","TV","VE","PD","RO","UD","GO","TS","PN","BO","FE",
           "FC","MO","PR","PC","RA","RE","RN")
CENTRE <- c("FI","AR","GR","LI","LU","MS","PI","PT","PO","SI","PG","TR","AN","AP","FM","MC","PU","RM","FR","LT","RI","VT")
area <- function(p) fifelse(p %in% NORTH, "North", fifelse(p %in% CENTRE, "Centre", fifelse(p == "unknown", "unknown", "South and islands")))
mode_chr <- function(v) { v <- v[!is.na(v) & nzchar(v)]; if (!length(v)) NA_character_ else names(which.max(table(v))) }

## ===========================================================================
## Pass 1: the codes over their whole life; the window rows, sharded by
## province; the national daily means of all reporting non-highway pumps
## ===========================================================================
parts <- if (dir.exists(CFG$data[1])) CFG$data[1] else {
  dz <- file.path(tempdir(), "fuels_station_level")
  if (!dir.exists(dz)) utils::unzip(CFG$data[1], exdir = dz)
  dz
}
parts <- list.files(parts, "^part_[0-9]+[.]RDS$", recursive = TRUE, full.names = TRUE)
ids <- list(); nat <- list()
for (k in seq_along(parts)) {
  x <- as.data.table(readRDS(parts[k]))[, .(id_pump, g = price_gasoline_self, d = price_diesel_self, date,
                                            brand, type = station_type, city, province,
                                            lat = latitude, lon = longitude)]
  ## the code of Naples is read as a missing value; a part whose province column
  ## is missing altogether (codes without province, type or coordinates) is
  ## left unknown, and its codes leave at step 4
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
  w <- x[date >= FROM & date <= TO]
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
saveRDS(means, file.path(OUT, "means.rds"))

## ===========================================================================
## Pass 2, province by province: linking, statistics, filters, register
## ===========================================================================
link_pairs <- list(); reg <- list(); kept <- list(); series <- list()
runmax <- function(v) { if (!any(v)) return(0L); r <- rle(v); max(r$lengths[r$values]) }
for (pv in sort(unique(ids$province))) {
  fl <- list.files(file.path(SH, pv), full.names = TRUE)
  if (!length(fl)) next
  w <- rbindlist(lapply(fl, readRDS))
  L <- ids[province == pv & id_pump %in% unique(w$id_pump)]
  ## all codes of the province alive at any time, for the linking
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
  ## codes without coordinates cannot be linked: each is a site of its own,
  ## removed at step 4
  nocoord <- setdiff(as.character(unique(w$id_pump)), names(root))
  root <- c(root, stats::setNames(nocoord, nocoord))
  w[, site := root[as.character(id_pump)]]
  w[ids, on = "id_pump", idfirst := i.first]
  setorder(w, site, date, -idfirst)
  w <- w[, .SD[1], by = .(site, date)]                         # overlapping days: the newer code
  ## site attributes: those of its latest code
  last_code <- w[, .(last_id = id_pump[which.max(idfirst)], n_codes = uniqueN(id_pump),
                     brand_start = brand[1], brand_end = brand[.N], n_brands = uniqueN(stats::na.omit(brand)),
                     type = type[.N]), by = site]
  last_code[ids, on = c(last_id = "id_pump"), `:=`(lat = i.lat, lon = i.lon, city = i.city)]
  ## statistics over the window
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
  ## distance from the median location of the province: a few sites are geocoded
  ## hundreds of kilometres away from their province
  st[, d_prov := NA_real_]
  okc <- is.finite(st$lat) & is.finite(st$lon)
  if (any(okc)) st[okc, d_prov := geodist::geodist(cbind(lon = lon, lat = lat), cbind(lon = stats::median(lon), lat = stats::median(lat)),
                                                   measure = "geodesic")[, 1] / 1000]
  ## the register: each step and what it leaves
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
    ks <- st[kept == TRUE]$site
    wk <- w[site %in% ks]
    di <- as.integer(wk$date - FROM) + 1L
    G <- Dm <- matrix(NA_real_, length(ks), ND, dimnames = list(ks, NULL))
    G[cbind(match(wk$site, ks), di)] <- wk$g; Dm[cbind(match(wk$site, ks), di)] <- wk$d
    series[[pv]] <- list(G = G, D = Dm, meta = st[kept == TRUE], days = DAYS)
  }
  message("pass 2: ", pv, " ", sum(keep), " kept of ", nrow(st), " sites")
}
reg <- rbindlist(reg)
sites <- rbindlist(kept, fill = TRUE)
links <- rbindlist(link_pairs)
saveRDS(list(register = reg, sites = sites, links = links, cfg = CFG), file.path(OUT, "pretreatment.rds"))
saveRDS(series, file.path(OUT, "metro_series.rds"))

## ===========================================================================
## The register, printed and written
## ===========================================================================
options(width = 160)
tot <- reg[, .(pumps = sum(n)), by = step][order(step)]
tot[, removed := shift(pumps) - pumps]
cat("\nREGISTER, Italy\n"); print(tot)
fwrite(tot, file.path(OUT, "register_italy.csv"))
comp <- function(by) {
  z <- reg[step %in% c("2 sites after linking codes", "9 one site per point"), .(n = sum(n)), by = c("step", by)]
  z[, share := round(100 * n / sum(n), 1), by = step]
  dcast(z, as.formula(paste(by, "~ step")), value.var = c("n", "share"), fill = 0)
}
cat("\nCOMPOSITION before (sites after linking) and after (kept), by brand group\n"); print(comp("group"))
cat("\nby area\n"); print(comp("area"))
cat("\nby station type\n"); print(comp("type"))
fwrite(comp("group"), file.path(OUT, "composition_group.csv"))
fwrite(comp("area"), file.path(OUT, "composition_area.csv"))
mt <- reg[province %in% CFG$metros, .(pumps = sum(n)), by = .(province, step)]
mt <- dcast(mt, province ~ step, value.var = "pumps")
cat("\nREGISTER, the metropolitan cities of the application\n"); print(mt)
fwrite(mt, file.path(OUT, "register_metros.csv"))
cat("\nLINKING: successive codes within 200 m, distance quantiles (m):\n")
print(round(stats::quantile(links$dist, c(.1, .25, .5, .75, .9)), 1))
cat("pairs within", CFG$d_link[1], "m:", sum(links$dist <= CFG$d_link[1]), "of", nrow(links), "\n")
unlink(SH, recursive = TRUE)
cat("\ndone\n")
