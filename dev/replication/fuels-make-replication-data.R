## ===========================================================================
## SC-STEM fuel application: THE DATA OF THE REPLICATION PACKAGE
##
## Writes the minimal data set of the replication package (replication/data/
## in the repository) from the output of the pre-treatment (stage 1 of
## fuels-data.R), which needs the full station-level file and is not part of
## the package. Three files:
##
##   fuel_stations.csv        one row per kept station of the eleven cities, in
##                            the order of the pre-treatment (the order matters:
##                            it is the order of the stations in every model)
##   fuel_prices.csv.gz       the Sunday price of gasoline and diesel of every
##                            station, euro per liter, presample and window
##                            (missing where the station did not report)
##   fuel_national_means.csv  the national mean price of every Sunday, over all
##                            the reporting non-highway stations of Italy
##
##     Rscript fuels-make-replication-data.R --root=<the fuels folder> --out=<replication/data>
##
## This is part of the development material, not of the replication package.
## ===========================================================================
suppressPackageStartupMessages(library(data.table))
HERE <- local({
  a <- commandArgs(trailingOnly = FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]), winslash = "/")) else getwd()
})
source(file.path(HERE, "fuels-functions.R"))
args <- commandArgs(trailingOnly = TRUE)
opt <- function(k, d) { a <- grep(paste0("^--", k, "="), args, value = TRUE); if (length(a)) sub("^--[^=]+=", "", a[1]) else d }
ROOT <- normalizePath(opt("root", file.path(HERE, "fuels")), winslash = "/", mustWork = TRUE)
OUT <- normalizePath(opt("out", file.path(HERE, "..", "..", "replication", "data")), winslash = "/", mustWork = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

ser <- readRDS(file.path(ROOT, "pretreatment", "metro_series.rds"))
mn <- readRDS(file.path(ROOT, "pretreatment", "means.rds"))
CITIES <- intersect(names(FU_CITY), names(ser))
days <- ser[[1]]$days
sun <- fu_sundays(days); weeks <- days[sun]

## the stations, city by city in the order of FU_CITY, each city in the order
## of the pre-treatment
## (the coordinates with 17 significant digits, so that they are read back
## exactly: the distances and the spatial covariances depend on them)
AREA <- c(North = "North", Centre = "Center", "South and islands" = "South and islands")
full <- function(v) ifelse(is.na(v), "", sprintf("%.17g", v))
st <- rbindlist(lapply(CITIES, function(cc) {
  m <- as.data.table(ser[[cc]]$meta)
  m[, .(city = cc, site = site, municipality = city, province = province, area = AREA[area],
        brand = brand_end, group = group, type = type, lat = full(lat), lon = full(lon))]
}))
fwrite(st, file.path(OUT, "fuel_stations.csv"))

## the Sunday prices, long: one row per station and Sunday
pr <- rbindlist(lapply(CITIES, function(cc) {
  s <- ser[[cc]]
  data.table(city = cc, site = rep(s$meta$site, times = length(sun)), week = rep(weeks, each = nrow(s$meta)),
             gasoline = as.vector(s$G[, sun]), diesel = as.vector(s$D[, sun]))
}))
setorder(pr, city, site, week)
fwrite(pr, file.path(OUT, "fuel_prices.csv.gz"))

## the national means of the Sundays, written with 17 significant digits so
## that they are read back exactly (they are averages, with long decimals)
nat <- mn$national[match(weeks, mn$national$date)]
fwrite(data.table(week = weeks, gasoline = full(nat$g), diesel = full(nat$d), stations = nat$n),
       file.path(OUT, "fuel_national_means.csv"), quote = FALSE)
cat(sprintf("stations: %d; Sundays: %d (%s to %s); files in %s\n", nrow(st), length(weeks),
            min(weeks), max(weeks), OUT))
