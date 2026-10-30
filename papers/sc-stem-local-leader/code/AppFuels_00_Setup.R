## ===========================================================================
## Latent Spatial Regimes in Spatio-Temporal State-Space Models: Evidence on
## Local Price Leadership in Italian Fuel Markets
##
## Author:  Paolo Maranzano (University of Milano-Bicocca; Fondazione Eni Enrico
##          Mattei)
## Script:  AppFuels_00_Setup.R, the setup of the application to Italian fuel
##          prices (Section 5 of the paper, Appendix C, Supplement C)
##
## What it does: defines the settings shared by the scripts of the application
## and the functions they call. It is sourced by AppFuels_01_RunSetting1.R,
## AppFuels_02_RunSetting2.R, AppFuels_03_RunSetting3.R and
## AppFuels_04_Analysis.R, and nothing in it runs on its own.
##
## Outputs: none.
## ===========================================================================
##
## The application in brief. In each of eleven Italian cities, the weekly
## relative price of a fuel at a station Y is regressed on four blocks of lags
## (the station's own price of the same fuel and of the other fuel, and the
## prices of the two fuels at the nearest station X of a given group) and on
## pulses for the changes of the excise duties. The regression is the
## observation equation of an SC-STEM model: the stations of a city are
## partitioned into latent regimes, each with its own coefficients, its own
## latent AR(1) process and its own spatially correlated error. The three
## settings define who is Y and who is X:
##
##   setting 1  Y a station of a major brand, X the nearest independent station
##   setting 2  Y an independent station,    X the nearest major-brand station
##   setting 3  Y a station of a major brand, X the nearest station of another
##              major brand
##
## A model is one setting, one city and one fuel of Y (3 x 11 x 2 = 66 models).
## Every model starts from the weekly data of data/ (the Sunday prices of the
## 3,136 stations, their registry and the national mean prices), from which
## fu_setting_data() builds the data of a setting: the relative prices, their
## missing weeks completed for the lags, the pairs (Y, X) and the pulses.
##
## Conventions. Prices are in cents per liter. A relative price is the Sunday
## price minus the national mean of that Sunday, minus the mean of the station
## over the estimation window. Matrices of prices are stations x weeks.
## Distances are geodesic: in meters inside the Stem package, in kilometers in
## the tables. Every function of a package other than base R is called as
## package::function().
##
## Contents
##   0. Settings shared by the scripts
##   1. Utilities: options, the Stem version, caching
##   2. The data: reading, relative prices, pulses, coordinates
##   3. Completion of the missing weeks of the lagged series
##   4. The pairs (Y, X) and the data of a setting
##   5. The series and the design of one city and fuel
##   6. Inference: Wald tests, benchmarks, diagnostics
##   7. The check of the local maxima of a selected fit
##   8. One model: the steps of a run
## ===========================================================================





## ===========================================================================
## 0. Settings shared by the scripts
## ===========================================================================

## The version of the Stem package the results were computed with (a commit of
## its GitHub repository). fu_require_stem() installs it when another version
## is installed.
FU_STEM_REF <- "PaoloMaranzano/Stem@09e9cdb2d13fdfc94d8bc74e62645ebbc0b16a54"

## The eleven cities (codes of the provinces) and their names.
FU_CITY <- c(RM = "Roma", MI = "Milano", "NA" = "Napoli", TO = "Torino", PA = "Palermo",
             BA = "Bari", CT = "Catania", BO = "Bologna", ME = "Messina", FI = "Firenze",
             VE = "Venezia")
## The English names, used in the tables and figures, and the order of the
## cities there (from north to south).
FU_CITY_EN <- c(Torino = "Turin", Milano = "Milan", Venezia = "Venice", Bologna = "Bologna",
                Firenze = "Florence", Roma = "Rome", Napoli = "Naples", Bari = "Bari",
                Palermo = "Palermo", Catania = "Catania", Messina = "Messina")
FU_CITY_ORDER <- names(FU_CITY_EN)
## The three geographic areas of Italy (ISTAT), by city.
FU_AREA <- c(Torino = "North", Milano = "North", Venezia = "North", Bologna = "North",
             Firenze = "Center", Roma = "Center", Napoli = "South and islands",
             Bari = "South and islands", Palermo = "South and islands",
             Catania = "South and islands", Messina = "South and islands")

## The major brands; the independent stations have group "independent" in the
## station file.
FU_MAJOR <- c("Agip Eni", "Api-Ip", "Esso", "Q8", "Tamoil")

## The fuels: the code of each and the code of the other one.
FU_FUEL <- c(g = "gasoline", d = "diesel")
FU_OTHER <- c(g = "d", d = "g")

## The estimation window starts on this Monday; the 52 weeks before it are a
## presample, which provides the lags of the first weeks of the window only.
FU_WINDOW_FROM <- as.Date("2022-01-03")

## The changes of the excise duties in the window. Each gives two pulses: one
## in the week that contains the change and one in the following week.
FU_EVENTS <- data.frame(
  date  = as.Date(c("2022-03-22", "2022-12-01", "2023-01-01", "2025-05-15", "2026-01-01")),
  label = c("excise_cut", "cut_reduced", "cut_ends", "realignment", "equal_excise"),
  stringsAsFactors = FALSE)

## The three settings.
FU_SETTINGS <- data.frame(
  setting = 1:3,
  y_set   = c("major", "independent", "major"),
  x_set   = c("independent", "major", "other major"),
  label   = c("Y major, X independent", "Y independent, X major", "Y major, X other major"),
  stringsAsFactors = FALSE)

## The four blocks of lags and the lags of each block (the last four weeks).
FU_BLOCKS <- c("own", "own_other", "nb", "nb_other")
FU_LAGS <- list(own = 1:4, own_other = 1:4, nb = 1:4, nb_other = 1:4)

## The settings of the SC-STEM models.
FU_M     <- 5L                                  # smallest size of a regime
FU_KCAP  <- 5L                                  # largest number of regimes
FU_KNN   <- 5L                                  # neighbors of the Potts graph
FU_PHI   <- c(0, 0.025, 0.05, 0.1, 0.2, 0.5, 1) # grid of the spatial penalty
FU_BAND  <- c(0.025, 0.2)                       # band of the two-step rule
FU_B     <- 100L                                # bootstrap draws
FU_LEVEL <- 0.95                                # level of the bootstrap intervals
FU_SEED  <- 20261005L                           # base of the seeds

## The seed of a model depends on its city and fuel only, so that a model
## gives the same results whichever subset of models is run.
fu_seed <- function(city, fuel)
  FU_SEED + 10L * match(city, names(FU_CITY)) + match(fuel, names(FU_FUEL))

## The largest number of regimes of a model with n stations: every regime must
## have at least FU_M stations, and there are at most FU_KCAP regimes. Below
## 2 FU_M stations only the pooled model (one regime) is fitted.
fu_kmax <- function(n) max(1L, min(FU_KCAP, n %/% FU_M))





## ===========================================================================
## 1. Utilities: options, the Stem version, caching
## ===========================================================================

## The options of a script: the defaults, overridden by the command line in
## the form --name=value (lists separated by commas). The type of each option
## is that of its default.
fu_options <- function(defaults, args = commandArgs(trailingOnly = TRUE)) {
  for (a in grep("^--[^=]+=.", args, value = TRUE)) {
    nm <- sub("^--([^=]+)=.*$", "\\1", a)
    if (!nm %in% names(defaults))
      stop("unknown option --", nm, "; the options are: ", paste(names(defaults), collapse = ", "),
           call. = FALSE)
    v <- trimws(strsplit(sub("^--[^=]+=", "", a), ",", fixed = TRUE)[[1]])
    d <- defaults[[nm]]
    defaults[[nm]] <- if (is.logical(d)) as.logical(v) else if (is.integer(d)) as.integer(v) else
      if (is.numeric(d)) as.numeric(v) else v
  }
  defaults
}

## Stem at the version of the results. The installed version is read from its
## DESCRIPTION file; when it is another one, Stem is installed from GitHub.
## Returns the first seven characters of the commit.
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
    if (!ok()) stop("Stem could not be installed at ", ref, ": restart R and run again",
                    call. = FALSE)
  }
  ## a session that already has another version of Stem loaded would use it
  if ("Stem" %in% loadedNamespaces()) {
    p <- getNamespaceInfo("Stem", "path")
    s <- tryCatch(unname(read.dcf(file.path(p, "DESCRIPTION"), fields = "RemoteSha")[1, 1]),
                  error = function(e) NA_character_)
    if (!identical(s, sha))
      stop("this R session has another version of Stem loaded: restart R and run again",
           call. = FALSE)
  }
  invisible(substr(sha, 1L, 7L))
}

## A step whose value is saved to <dir>/<name>.rds and read back on a later
## run instead of being computed again; delete the file to compute it again.
## The minutes the computation took are kept as an attribute of the value.
fu_cached <- function(dir, name, expr) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(dir, paste0(name, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  t0 <- Sys.time()
  val <- force(expr)
  attr(val, "minutes") <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
  saveRDS(val, f)
  val
}
fu_minutes <- function(x) if (is.null(attr(x, "minutes"))) NA_real_ else attr(x, "minutes")

## Writing a table: a .csv file, and a message.
fu_write <- function(d, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(d, file, row.names = FALSE)
  message("wrote ", file)
}





## ===========================================================================
## 2. Data: reading, relative prices, pulses, coordinates
## ===========================================================================

## The three files of the data (see the README of the folder of the paper):
## the stations, their Sunday prices and the national mean prices. The empty
## field marks a missing value, so that the code of Naples, "NA", is read as a
## string. Returns a list: stations (a data.frame, in the order of the file,
## which is the order of the stations in every model), prices (long: one row
## per station and Sunday) and national (one row per Sunday).
fu_read_data <- function(dir) {
  rd <- function(f, ...) as.data.frame(data.table::fread(file.path(dir, f), na.strings = "", ...))
  st <- rd("fuel_stations.csv", colClasses = list(character = c("city", "site", "province")))
  pr <- rd("fuel_prices.csv.gz", colClasses = list(character = c("city", "site")))
  nat <- rd("fuel_national_means.csv")
  pr$week <- as.Date(pr$week); nat$week <- as.Date(nat$week)
  list(stations = st, prices = pr, national = nat)
}

## The prices of one fuel of the stations of one city, stations x weeks
## (euro per liter; NA where the station did not report on that Sunday). The
## rows follow the station file, the columns the vector of Sundays `weeks`.
fu_price_matrix <- function(dat, city, fuel, weeks) {
  s <- dat$stations$site[dat$stations$city == city]
  x <- dat$prices[dat$prices$city == city, ]
  P <- matrix(NA_real_, length(s), length(weeks), dimnames = list(s, NULL))
  P[cbind(match(x$site, s), match(x$week, weeks))] <- x[[FU_FUEL[[fuel]]]]
  P
}

## The relative prices, stations x weeks, in cents per liter: the Sunday price
## minus the national mean of that Sunday (this removes the stochastic trend of
## costs and taxes common to all the stations of Italy), minus the mean of the
## station over the weeks of the window it reported (this removes its permanent
## premium). The mean of the window is also subtracted from the presample. A
## Sunday without a report stays missing.
fu_relative <- function(P, national, in_window) {
  C <- 100 * sweep(P, 2, national)
  C - rowMeans(C[, in_window, drop = FALSE], na.rm = TRUE)
}

## The pulses of the changes of the excise duties, weeks x (2 x events): 1 in
## the week that contains the change (the first Sunday on or after it) and in
## the following week, 0 elsewhere.
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

## The coordinates of a set of stations as the Stem package wants them
## (longitude, latitude). Stations at identical coordinates are moved apart by
## a few centimeters, so that the spatial covariance matrix stays non-singular.
fu_coords <- function(meta) {
  co <- cbind(lon = meta$lon, lat = meta$lat)
  dup <- duplicated(round(co, 6))
  if (any(dup)) co[dup, ] <- co[dup, ] + 1e-5 * seq_len(sum(dup))
  co
}





## ===========================================================================
## 3. Completion of the missing weeks of the lagged series
## ===========================================================================
## The response of a model may have missing weeks, which the Kalman filter of
## the model handles exactly. Its covariates may not, and the lagged prices are
## covariates. The missing weeks are therefore filled once, before any model is
## fitted, by a pooled STEM model of the relative prices of all the stations of
## the city (an intercept and the pulses as covariates, one latent AR(1)
## process common to the city, and a spatially correlated error), through
## Stem::STEM_Complete(), which returns the conditional mean of every missing
## value given all the observed ones. The completed series are used only to
## build the lags: the response keeps its missing weeks.
##
## Y: relative prices, stations x weeks; meta: the stations (lon, lat);
## pulses: weeks x pulses. Returns the completed matrix (stations x weeks) and
## a one-row summary of the pooled fit.
fu_complete <- function(Y, meta, pulses) {
  d <- nrow(Y); Tn <- ncol(Y)
  pul <- pulses[, colSums(pulses) > 0, drop = FALSE]
  X1 <- cbind(intercept = 1, pul)
  X <- do.call(rbind, rep(list(X1), d))                        # stacked by station
  z <- t(Y)                                                    # weeks x stations
  ok <- !is.na(as.vector(z))
  ## starting values: least squares for the coefficients, a split of the
  ## residual variance between the two error components, and a range equal to
  ## the median distance between the stations
  ols <- stats::lm.fit(X[ok, , drop = FALSE], as.vector(z)[ok])
  b <- ols$coefficients; b[!is.finite(b)] <- 0
  s2 <- stats::var(ols$residuals)
  co <- fu_coords(meta)
  dm <- stats::median(geodist::geodist(co, measure = "geodesic"))
  mod <- Stem::STEM_Model(z = z, covariates = X, coordinates = co,
                          phi = list(beta = matrix(b, ncol = 1), sigma2eps = 0.6 * s2,
                                     sigma2omega = 0.4 * s2, theta = 3 / dm,
                                     G = matrix(0.5, 1, 1), Sigmaeta = matrix(0.2 * s2, 1, 1),
                                     m0 = as.matrix(0), C0 = as.matrix(1)),
                          A = matrix(1, d, 1))
  fit <- Stem::STEM_Estimation(mod, distance = "geo", control = Stem::STEM_control())
  comp <- Stem::STEM_Complete(fit, distance = "geo")           # weeks x stations
  ph <- fit$estimates$phi.hat
  list(completed = t(comp),
       summary = data.frame(stations = d, weeks = Tn, missing = sum(is.na(Y)),
                            missing_pct = round(100 * mean(is.na(Y)), 3),
                            loglik = fit$estimates$loglik,
                            converged = isTRUE(fit$estimates$convergence.par$converged),
                            G = ph$G[1, 1], sigma2eps = ph$sigma2eps, sigma2omega = ph$sigma2omega,
                            range_km = 3 / ph$theta / 1000))
}





## ===========================================================================
## 4. The pairs (Y, X) and the data of a setting
## ===========================================================================

## The pairs of a setting in a city: for every station Y of the Y group, its
## nearest station X of the X group (in setting 3, the nearest major-brand
## station of a brand other than that of Y). meta: the stations of the city.
## Returns one row per Y: the rows of Y and X in meta and their distance (km).
fu_pairs <- function(meta, setting) {
  Dkm <- geodist::geodist(meta[, c("lon", "lat")], measure = "geodesic") / 1000
  maj <- which(meta$brand %in% FU_MAJOR)
  ind <- which(meta$group == "independent")
  ys <- switch(as.character(setting), "1" = maj, "2" = ind, "3" = maj)
  xset <- function(i) switch(as.character(setting), "1" = ind, "2" = maj,
                             "3" = maj[meta$brand[maj] != meta$brand[i]])
  rows <- lapply(ys, function(i) {
    xs <- xset(i)
    if (!length(xs)) return(NULL)
    dx <- Dkm[i, xs]
    data.frame(y_row = i, x_row = xs[which.min(dx)], x_dist_km = round(min(dx), 3))
  })
  do.call(rbind, rows)
}

## The data of a setting in a city: one row per station Y and week, with
##   city, setting                          the model
##   y_site, y_brand, y_group, y_type,      the station Y and its attributes
##   y_municipality, y_lon, y_lat
##   week, t, presample                     the Sunday, its index, and whether it
##                                          belongs to the presample (lags only)
##   yg, yd                                 the relative prices of Y: the response,
##                                          with its missing weeks
##   yg_c, yd_c                             the same, completed (for the lags)
##   xg, xd                                 the completed relative prices of X
##   x_site, x_brand, x_dist_km             the station X and its distance to Y
##   pulse_*                                the pulses of the excise duties
## meta: the stations of the city; Yrel and Ycomp: lists by fuel ("g", "d") of
## stations x weeks matrices (relative and completed); pairs: fu_pairs().
fu_case_frame <- function(city, setting, meta, Yrel, Ycomp, weeks, pairs, pulses, presample) {
  one <- function(k) {
    i <- pairs$y_row[k]; x <- pairs$x_row[k]
    data.frame(
      city = city, setting = setting,
      y_site = meta$site[i], y_brand = meta$brand[i], y_group = meta$group[i],
      y_type = meta$type[i], y_municipality = meta$municipality[i],
      y_lon = meta$lon[i], y_lat = meta$lat[i],
      week = weeks, t = seq_along(weeks), presample = presample,
      yg = Yrel$g[i, ], yd = Yrel$d[i, ], yg_c = Ycomp$g[i, ], yd_c = Ycomp$d[i, ],
      xg = Ycomp$g[x, ], xd = Ycomp$d[x, ],
      x_site = meta$site[x], x_brand = meta$brand[x], x_dist_km = pairs$x_dist_km[k],
      pulses, stringsAsFactors = FALSE, row.names = NULL)
  }
  do.call(rbind, lapply(seq_len(nrow(pairs)), one))
}

## The data of the models of setting s, from the weekly data in data_dir: for
## every city of `cities`, the relative prices of its stations (fu_relative()),
## their missing weeks completed by a pooled STEM model of each fuel
## (fu_complete(); each completion is saved in cache_dir and read back on a
## later run, and the three settings share it), the pairs (Y, X) of the setting
## (fu_pairs()), and one row per station Y and week (fu_case_frame()). The 52
## Sundays before FU_WINDOW_FROM are the presample, which provides the lags of
## the first weeks of the window only. Returns the data.frame of the setting,
## the cities in the order of FU_CITY.
fu_setting_data <- function(s, data_dir, cache_dir, cities = names(FU_CITY)) {
  dat <- fu_read_data(data_dir)
  weeks <- dat$national$week
  in_window <- weeks >= FU_WINDOW_FROM
  pulses <- fu_pulses(weeks)
  nat <- list(g = dat$national$gasoline, d = dat$national$diesel)
  d <- do.call(rbind, lapply(intersect(names(FU_CITY), cities), function(cc) {
    meta <- dat$stations[dat$stations$city == cc, ]
    pairs <- fu_pairs(meta, s)
    if (is.null(pairs)) return(NULL)
    rel <- lapply(c(g = "g", d = "d"), function(fu)
      fu_relative(fu_price_matrix(dat, cc, fu, weeks), nat[[fu]], in_window))
    comp <- lapply(c(g = "g", d = "d"), function(fu) {
      message(sprintf("completion: %s, %s (%d stations)", FU_CITY[[cc]], FU_FUEL[[fu]], nrow(meta)))
      fu_cached(cache_dir, sprintf("completed_%s_%s", cc, FU_FUEL[[fu]]),
                fu_complete(rel[[fu]], meta, pulses))$completed
    })
    fu_case_frame(cc, s, meta, rel, comp, weeks, pairs, pulses, !in_window)
  }))
  rownames(d) <- NULL
  d
}





## ===========================================================================
## 5. The series and the design of one city and fuel
## ===========================================================================

## The series of the four blocks of one city and fuel of Y, weeks x stations Y
## (presample and window): the response and the completed series of the own
## price, of the own price of the other fuel and of the two prices of X.
fu_series <- function(dfc, fuel) {
  fo <- FU_OTHER[[fuel]]
  dfc <- dfc[order(match(dfc$y_site, unique(dfc$y_site)), dfc$t), ]
  W <- max(dfc$t)
  m <- function(col) matrix(dfc[[col]], nrow = W)
  list(resp = m(paste0("y", fuel)), own = m(paste0("y", fuel, "_c")),
       own_other = m(paste0("y", fo, "_c")), nb = m(paste0("x", fuel)),
       nb_other = m(paste0("x", fo)), sites = unique(dfc$y_site))
}

## The design of one model: the response (weeks x stations Y, with its missing
## weeks), the covariates stacked by station (the layout of Stem::STEM_Model()),
## their names, and the columns of each block. Column names: intercept, the
## lags of the four blocks (own_l1, ..., nb_other_l4) and the pulses of the
## events that fall in the window. The response starts at the first week of the
## window: the presample provides the lags of its first weeks.
fu_design <- function(dfc, fuel, lags = FU_LAGS) {
  S <- fu_series(dfc, fuel)
  W <- nrow(S$own); d <- ncol(S$own)
  lags <- lags[FU_BLOCKS]
  maxlag <- max(c(0L, unlist(lags)))
  t0 <- max(maxlag + 1L, min(dfc$t[!dfc$presample]))
  tt <- t0:W
  first <- dfc$y_site == S$sites[1]
  pul <- as.matrix(dfc[first, grep("^pulse_", names(dfc)), drop = FALSE])
  pul <- pul[order(dfc$t[first]), , drop = FALSE][tt, , drop = FALSE]
  pul <- pul[, colSums(pul) > 0, drop = FALSE]
  ## the lags of the four blocks of station i, one column per lag
  block_cols <- function(i) do.call(cbind, lapply(FU_BLOCKS, function(b) {
    L <- lags[[b]]
    if (!length(L)) return(NULL)
    m <- vapply(L, function(l) S[[b]][tt - l, i], numeric(length(tt)))
    if (is.null(dim(m))) m <- matrix(m, ncol = length(L))
    colnames(m) <- paste0(b, "_l", L)
    m
  }))
  covs <- lapply(seq_len(d), function(i) cbind(intercept = 1, block_cols(i), pul))
  nm <- colnames(covs[[1]])
  blocks <- lapply(stats::setNames(FU_BLOCKS, FU_BLOCKS), function(b) grep(paste0("^", b, "_l"), nm))
  blocks$pulses <- grep("^pulse_", nm)
  meta <- unique(dfc[, c("y_site", "y_brand", "y_group", "y_type", "y_municipality", "y_lon", "y_lat",
                         "x_site", "x_brand", "x_dist_km")])
  meta <- meta[match(S$sites, meta$y_site), ]
  list(z = S$resp[tt, , drop = FALSE], covariates = do.call(rbind, covs), names = nm,
       blocks = blocks, coords = fu_coords(data.frame(lon = meta$y_lon, lat = meta$y_lat)),
       meta = meta, weeks = sort(unique(dfc$week))[tt], d = d, Tn = length(tt),
       city = dfc$city[1], fuel = fuel)
}

## The STEM model of a design, with starting values: least squares on the
## observed weeks for the coefficients, a split of the residual variance
## between the two error components, and a range equal to the median distance
## between the stations. SC-STEM starts every regime from these values.
fu_model <- function(des) {
  zv <- as.vector(des$z)                                       # stacked by station
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





## ===========================================================================
## 6. Inference: Wald tests, benchmarks, diagnostics
## ===========================================================================

## The hypotheses, as sets of columns of the design; each is the null that the
## coefficients of those lags are zero:
##   nb_both    the neighbor X does not lead Y (its eight lags, both fuels):
##              the test of leadership of the paper
##   nb         X does not lead Y on the same fuel (four lags)
##   nb_other   X does not lead Y across fuels (four lags)
##   own_other  the other fuel of Y does not lead the fuel of Y (four lags)
fu_tests <- function(des) {
  b <- des$blocks
  list(nb_both = c(b$nb, b$nb_other), nb = b$nb, nb_other = b$nb_other, own_other = b$own_other)
}

## The bootstrap Wald test of every hypothesis of fu_tests(), regime by regime
## and in all the regimes at once.
##
## est: the estimated coefficients of the fit, regimes x parameters, with
## columns named beta1, beta2, ... (the columns of the design); group: the
## regime of every station; inf: Stem::SCSTEM_BootInference() of the
## refit-with-clustering bootstrap of the fit, whose `aligned` array holds the
## draws of every parameter with the regimes of each draw matched to the
## estimated ones. For a block of q coefficients b_k of regime k:
##   V_k  the covariance of the draws of b_k;
##   W_k  = b_k' V_k^-1 b_k, against the chi-square on q degrees of freedom;
##   p_boot: W_k against its bootstrap distribution under the null, the draws
##        (b*_k - b_k)' V_k^-1 (b*_k - b_k) centered at the estimate.
## Row "all" stacks the coefficients of the K regimes (q K degrees of freedom)
## with their joint bootstrap covariance. Draws with a missing value are left
## out.
fu_wald <- function(est, group, inf, des, label) {
  K <- nrow(est); A <- inf$aligned
  one <- function(e, D) {
    D <- D[stats::complete.cases(D), , drop = FALSE]
    if (nrow(D) <= length(e)) return(c(NA_real_, NA_real_, nrow(D)))
    V <- stats::cov(D)
    W <- tryCatch(drop(crossprod(e, solve(V, e))), error = function(err) NA_real_)
    Vi <- tryCatch(solve(V), error = function(err) NULL)
    if (is.null(Vi) || !is.finite(W)) return(c(W, NA_real_, nrow(D)))
    C <- sweep(D, 2, e)
    null <- rowSums((C %*% Vi) * C)
    c(W, (1 + sum(null >= W)) / (1 + length(null)), nrow(D))
  }
  tests <- fu_tests(des)
  do.call(rbind, lapply(names(tests), function(h) {
    par <- paste0("beta", tests[[h]]); ip <- match(par, dimnames(A)[[3]]); q <- length(ip)
    e <- lapply(seq_len(K), function(g) unname(est[g, par]))
    D <- lapply(seq_len(K), function(g) matrix(A[, g, ip], ncol = q))
    r <- rbind(do.call(rbind, lapply(seq_len(K), function(g) one(e[[g]], D[[g]]))),
               one(unlist(e), do.call(cbind, D)))
    data.frame(fit = label, block = h, regime = c(as.character(seq_len(K)), "all"),
               pumps = c(tabulate(group, K), length(group)), statistic = r[, 1],
               df = c(rep(q, K), K * q), p_chisq = stats::pchisq(r[, 1], c(rep(q, K), K * q), lower.tail = FALSE),
               p_boot = r[, 2], draws_used = as.integer(r[, 3]), stringsAsFactors = FALSE)
  }))
}

## The estimates of a fit with their bootstrap standard errors and percentile
## intervals (the summary of Stem::SCSTEM_BootInference()), with the name of
## the covariate of every regression coefficient.
fu_estimates <- function(inf, names, label) {
  e <- as.data.frame(inf$summary)
  num <- suppressWarnings(as.integer(sub("^beta", "", e$parameter)))
  e$covariate <- ifelse(is.na(num), e$parameter, names[num])
  cbind(fit = label, e, stringsAsFactors = FALSE)
}

## The aligned bootstrap draws of the regression coefficients of a fit, long:
## one row per draw, regime and covariate.
fu_draws <- function(inf, names) {
  A <- inf$aligned
  par <- grep("^beta", dimnames(A)[[3]], value = TRUE)
  g <- expand.grid(draw = seq_len(dim(A)[1]), regime = seq_len(dim(A)[2]), parameter = par,
                   KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  data.frame(draw = g$draw, regime = g$regime,
             covariate = names[as.integer(sub("^beta", "", g$parameter))],
             value = as.vector(A[, , par, drop = FALSE]), stringsAsFactors = FALSE)
}

## The station-by-station benchmark: for every station Y, the regression of
## the model by least squares on the weeks with an observed response, without
## latent process, spatial error or regimes; for every hypothesis of
## fu_tests(), the classical F test and the Wald test with the
## heteroskedasticity-consistent covariance HC1.
fu_pump_tests <- function(des) {
  tests <- fu_tests(des)
  rows <- lapply(seq_len(des$d), function(i) {
    Z <- des$covariates[(i - 1L) * des$Tn + seq_len(des$Tn), , drop = FALSE]
    y <- des$z[, i]; ok <- is.finite(y)
    Z <- Z[ok, , drop = FALSE]; y <- y[ok]
    ## the columns that are not collinear (a pulse may fall in a missing week)
    q <- qr(Z); keep <- sort(q$pivot[seq_len(q$rank)])
    Zk <- Z[, keep, drop = FALSE]; n <- nrow(Zk); k <- ncol(Zk)
    A <- tryCatch(solve(crossprod(Zk)), error = function(e) NULL)
    if (is.null(A) || n <= k) return(NULL)
    b <- A %*% crossprod(Zk, y); e <- drop(y - Zk %*% b)
    V_hc1 <- A %*% crossprod(Zk * e) %*% A * n / (n - k)
    s2 <- sum(e^2) / (n - k)
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

## The panel test of Dumitrescu and Hurlin (2012) on the station-level Wald
## statistics of a block: Zbar = sqrt(N / (2 q)) (mean(W) - q), standard normal
## under non-causality for every station. Computed on the classical statistic
## (q times F) and on the HC1 one, with the shares of the stations whose own
## test rejects at 5%.
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

## Residual diagnostics, station by station on the weeks with an observed
## response, summarized by regime: the shares of stations for which the
## Ljung-Box test at lag 4 rejects at 5% on the residuals and on their squares,
## and the medians of the first-order autocorrelation and of the kurtosis. The
## residuals are the response minus the fitted signal of the regime.
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

## The adjusted Rand index between two partitions (Hubert and Arabie 1985):
## 1 when they coincide, about 0 when they are unrelated.
fu_ari <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  c2 <- function(v) sum(v * (v - 1) / 2)
  s_ij <- c2(as.vector(tab)); s_i <- c2(rowSums(tab)); s_j <- c2(colSums(tab)); tot <- n * (n - 1) / 2
  e <- s_i * s_j / tot; mx <- (s_i + s_j) / 2
  if (isTRUE(all.equal(mx, e))) return(1)
  (s_ij - e) / (mx - e)
}

## The adjusted Rand index between a partition of the stations and a purely
## spatial partition with the same number of groups: k-means on the
## coordinates of the stations, projected on a plane in kilometers (longitude
## scaled by the cosine of the mean latitude), with many random starts. It is
## the benchmark of the index between the regimes and the brands: if the
## regimes were districts, they would agree with the spatial partition more
## than with the brands. The seed makes the random starts reproducible; the
## random stream of the caller is restored on exit, so that the random jitter
## of the figures drawn afterwards does not change.
fu_ari_space <- function(regime, lon, lat, seed, nstart = 50L) {
  K <- length(unique(regime))
  if (K < 2L) return(NA_real_)
  xy <- 111.2 * cbind(lon * cos(mean(lat) * pi / 180), lat)
  if (exists(".Random.seed", envir = globalenv())) {
    old <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old, envir = globalenv()))
  } else on.exit(rm(".Random.seed", envir = globalenv()))
  set.seed(seed)
  fu_ari(regime, stats::kmeans(xy, K, nstart = nstart, iter.max = 100L)$cluster)
}





## ===========================================================================
## 7. The check of the local maxima of a selected fit
## ===========================================================================
## The likelihood of a regime may have more than one local maximum, and the EM
## algorithm reaches the one its start leads to. The selected partition is
## refitted regime by regime from the least-squares start of fu_model(), with
## the alternation switched off (max_iter = 0) and the EM algorithm stopping
## only when both the parameters and the log-likelihood have settled. Returns,
## by regime, the log-likelihoods and the main parameters at the two maxima
## with the Wald test of leadership at both (on the same bootstrap draws), and,
## by regime and covariate, the coefficients at both with their bootstrap
## standard errors. NULL for a pooled model, which has one start only.
fu_check_optima <- function(fit, des, inf, seed) {
  K <- nrow(fit$phi_hat)
  if (K == 1L) return(NULL)
  r <- Stem::SCSTEM_Estimation(
    fu_model(des), K = K, phi_penalty = 0, init_partition = fit$group, max_iter = 0,
    distance = "geo", knn = min(FU_KNN, des$d - 1L), min_cluster_size = FU_M, seed = seed,
    control = Stem::STEM_control(em_stop = "all", em_maxit = 2000L), verbose = FALSE)
  w_sel <- fu_wald(fit$phi_hat, fit$group, inf, des, "selected")
  w_ls <- fu_wald(r$phi_hat, fit$group, inf, des, "least squares")
  pick <- function(w, col) w[[col]][w$block == "nb_both" & w$regime != "all"]
  regimes <- data.frame(
    regime = seq_len(K), size = tabulate(fit$group, K),
    loglik_selected = fit$loglik_g, loglik_ls = r$loglik_g, diff = r$loglik_g - fit$loglik_g,
    converged_selected = fit$em_converged, converged_ls = r$em_converged,
    theta_selected = unname(fit$phi_hat[, "theta"]), theta_ls = unname(r$phi_hat[, "theta"]),
    G_selected = unname(fit$phi_hat[, "G"]), G_ls = unname(r$phi_hat[, "G"]),
    sigma2eps_selected = unname(fit$phi_hat[, "sigma2eps"]), sigma2eps_ls = unname(r$phi_hat[, "sigma2eps"]),
    W_nb_both_selected = pick(w_sel, "statistic"), p_nb_both_selected = pick(w_sel, "p_chisq"),
    W_nb_both_ls = pick(w_ls, "statistic"), p_nb_both_ls = pick(w_ls, "p_chisq"))
  sm <- as.data.frame(inf$summary)
  block <- rep("intercept", length(des$names))
  for (b in names(des$blocks)) block[des$blocks[[b]]] <- b
  p <- paste0("beta", seq_along(des$names))
  coefficients <- do.call(rbind, lapply(seq_len(K), function(g) data.frame(
    regime = g, covariate = des$names, block = block,
    beta_selected = unname(fit$phi_hat[g, p]), beta_ls = unname(r$phi_hat[g, p]),
    se_boot = sm$se[sm$cluster == g][match(p, sm$parameter[sm$cluster == g])],
    stringsAsFactors = FALSE)))
  list(regimes = regimes, coefficients = coefficients)
}





## ===========================================================================
## 8. One model: the steps of a run
## ===========================================================================
## Fits one model (setting s, city cc, fuel fu) and writes its compact results
## into <opt$out>/setting<s>/<city>_<fuel>/. cases: the data of the setting
## (fu_setting_data()); opt: the options of the run script (out,
## cache, cores); stem: the commit of Stem, recorded in model.csv. A model whose
## model.csv is already there is skipped, and the grid, the bootstraps and the
## check of the local maxima are cached in <opt$cache>/setting<s>, so a model
## stopped halfway resumes from its last completed step. The seed depends on
## the city and the fuel only (fu_seed()).
##
## The steps:
##   1. the design: the response and the covariates (fu_design());
##   2. the regimes: the grid of (K, phi) and the two-step rule on the BIC;
##   3. the refit-with-clustering bootstrap of the selected fit and of the
##      pooled model;
##   4. the tests: the bootstrap Wald tests and the station-by-station
##      benchmark;
##   5. the residual diagnostics of the selected fit;
##   6. the check of the local maxima of the selected fit;
##   7. the compact results of the model, read by AppFuels_04_Analysis.R.
fu_run_model <- function(s, cc, fu, cases, opt, stem) {
  id <- sprintf("%s_%s", cc, FU_FUEL[[fu]])
  od <- file.path(opt$out, sprintf("setting%d", s), id)
  if (file.exists(file.path(od, "model.csv"))) { message("done already: setting ", s, ", ", id); return(invisible()) }
  cache <- file.path(opt$cache, sprintf("setting%d", s))
  seed <- fu_seed(cc, fu)

  ## 1. The design
  des <- fu_design(cases[cases$city == cc, , drop = FALSE], fu)
  kmax <- fu_kmax(des$d)
  message(sprintf("\n== setting %d, %s, %s: %d stations Y, %d weeks, K_max = %d", s, FU_CITY[[cc]],
                  FU_FUEL[[fu]], des$d, des$Tn, kmax))

  ## 2. The regimes: the grid of (K, phi) and the two-step rule on the BIC.
  ## Only the table of the grid, the selected fit and the pooled fit are kept
  ## (the grid of a large city takes some GB).
  g <- fu_cached(cache, paste0(id, "_grid"), {
    ic <- Stem::SCSTEM_Infocrit(fu_model(des), K_grid = seq_len(kmax), phi_grid = FU_PHI,
                                knn = min(FU_KNN, des$d - 1L), distance = "geo", min_cluster_size = FU_M,
                                seed = seed, verbose = TRUE, control = Stem::STEM_control())
    sel <- Stem::SCSTEM_Select(ic, band = FU_BAND, criterion = "BIC")
    list(table = ic$table, K = sel$K_selected, phi = sel$phi_selected, fit = sel$fit,
         pooled = ic$fits[[which(ic$table$K == 1)[1]]])
  })
  K <- g$K
  message(sprintf("   K = %d, phi = %s, regimes of %s stations", K, format(g$phi),
                  paste(tabulate(g$fit$group, K), collapse = "/")))

  ## 3. The refit-with-clustering bootstrap: each draw simulates the data from
  ## the fitted model, regime by regime, and re-estimates SC-STEM at the
  ## selected (K, phi), partition included; the regimes of a draw are matched
  ## to the estimated ones (Stem::SCSTEM_BootInference()). The pooled model is
  ## bootstrapped in the same way, with the same seed.
  boot <- function(fit, name) fu_cached(cache, paste0(id, "_", name), {
    bt <- Stem::SCSTEM_Bootstrap(fit, B = FU_B, seed = seed, verbose = TRUE, cores = opt$cores)
    list(B_valid = bt$B_valid, inf = Stem::SCSTEM_BootInference(bt, level = FU_LEVEL))
  })
  b0 <- boot(g$fit, "boot")
  bP <- if (K == 1L) b0 else boot(g$pooled, "boot_pooled")

  ## 4. The tests
  wald <- rbind(fu_wald(g$fit$phi_hat, g$fit$group, b0$inf, des, "selected"),
                if (K > 1L) fu_wald(g$pooled$phi_hat, g$pooled$group, bP$inf, des, "pooled"))
  pt <- fu_pump_tests(des)
  dh <- fu_dh(pt)

  ## 5. The residual diagnostics
  dg <- fu_diagnostics(g$fit, des)

  ## 6. The check of the local maxima (NULL for a pooled model)
  opt_check <- fu_cached(cache, paste0(id, "_optima"), fu_check_optima(g$fit, des, b0$inf, seed))

  ## 7. The compact results
  dir.create(od, recursive = TRUE, showWarnings = FALSE)
  w <- function(d, f) if (!is.null(d)) utils::write.csv(cbind(setting = s, model = id, d), file.path(od, f),
                                                         row.names = FALSE)
  w(g$table[, !vapply(g$table, is.list, TRUE)], "grid.csv")
  w(data.frame(des$meta, regime = g$fit$group, theta_bound = g$fit$theta_bound[g$fit$group]), "regimes.csv")
  w(rbind(fu_estimates(b0$inf, des$names, "selected"),
          if (K > 1L) fu_estimates(bP$inf, des$names, "pooled")), "estimates.csv")
  data.table::fwrite(cbind(setting = s, model = id, fu_draws(b0$inf, des$names)), file.path(od, "draws.csv.gz"))
  w(wald, "wald.csv"); w(pt, "pump_tests.csv"); w(dh, "dumitrescu_hurlin.csv"); w(dg, "diagnostics.csv")
  if (!is.null(opt_check)) { w(opt_check$regimes, "optima.csv"); w(opt_check$coefficients, "optima_coefficients.csv") }
  brands <- des$meta$y_brand
  row <- data.frame(
    setting = s, model = id, city = cc, fuel = FU_FUEL[[fu]], pumps = des$d, weeks = des$Tn,
    K_max = kmax, K = K, phi = g$phi, BIC = g$table$BIC[g$table$K == K & g$table$phi == g$phi][1],
    sizes = paste(tabulate(g$fit$group, K), collapse = "/"),
    B_valid = b0$B_valid, B_valid_pooled = bP$B_valid,
    ari_regimes_brands = if (length(unique(brands)) > 1 && K > 1) fu_ari(g$fit$group, brands) else NA_real_,
    ari_draws_median = stats::median(b0$inf$stability$ARI),
    minutes_grid = fu_minutes(g), minutes_boot = fu_minutes(b0),
    minutes_boot_pooled = if (K > 1L) fu_minutes(bP) else 0, minutes_optima = fu_minutes(opt_check),
    stem = stem, stringsAsFactors = FALSE)
  utils::write.csv(row, file.path(od, "model.csv"), row.names = FALSE)   # written last: the model is done
  message(sprintf("   done: X leads Y (8 lags, 5%%) in %d of %d regimes", sum(wald$fit == "selected" &
    wald$block == "nb_both" & wald$regime != "all" & wald$p_chisq < 0.05), K))
  invisible()
}

## Every model of setting s in the cities and fuels of opt, the smaller first,
## then the table of the models of the setting (<opt$out>/setting<s>/models.csv,
## from the files of all its models). A model that fails is reported and the
## run goes on.
fu_run_setting <- function(s, cases, opt, stem) {
  n <- tapply(cases$y_site, cases$city, function(v) length(unique(v)))
  todo <- expand.grid(city = intersect(names(FU_CITY), intersect(opt$cities, names(n))),
                      fuel = opt$fuels, stringsAsFactors = FALSE)
  todo$n <- as.integer(n[todo$city])
  todo <- todo[order(todo$n), ]
  message(sprintf("setting %d: %d models; Stem %s; bootstrap on %d process(es)", s, nrow(todo), stem, opt$cores))
  for (i in seq_len(nrow(todo))) {
    r <- tryCatch(fu_run_model(s, todo$city[i], todo$fuel[i], cases, opt, stem), error = function(e) e)
    if (inherits(r, "error"))
      message("   FAILED: setting ", s, ", ", todo$city[i], ", ", todo$fuel[i], ": ", conditionMessage(r))
  }
  f <- Sys.glob(file.path(opt$out, sprintf("setting%d", s), "*", "model.csv"))
  if (length(f)) fu_write(do.call(rbind, lapply(f, utils::read.csv, stringsAsFactors = FALSE)),
                          file.path(opt$out, sprintf("setting%d", s), "models.csv"))
  message("done: results of setting ", s, " in ", opt$out)
  invisible()
}
