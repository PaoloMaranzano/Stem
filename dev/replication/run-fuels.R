## ===========================================================================
## SC-STEM application: Granger causality between gasoline and diesel at the
## pump, province by province, on weekly data. THE RUNNER (first version).
##
## The idea, the data and the choices behind this script are described in
## fuels-application-notes.tex, in the Overleaf project of the paper.
##
## For each province and each direction (diesel -> gasoline, gasoline ->
## diesel) the script
##   A. builds the weekly panel of the self-service prices of the pumps of the
##      province that sell both fuels over the whole window;
##   B. writes the model in weekly changes, in cents per litre:
##        dy_it = a + sum_l b_l dy_i,t-l + sum_l c_l dx_i,t-l
##                + g (y - x)*_i,t-1 + fiscal pulses + K y_t + e_it
##      with the own lags of y, the lags of the other fuel x, the lagged
##      gasoline-diesel spread demeaned by pump and fiscal period (the error
##      correction term), a pulse at each change of the excise duties, the
##      latent AR(1) common to the regime and the spatially correlated error;
##   C. fits the grid of (k, phi) and selects by the two-step rule of
##      SCSTEM_Select(), the pooled STEM model competing; x does not
##      Granger-cause y in a regime when all c_l are zero there, tested by the
##      likelihood ratio against the fit without the lags of x on the same
##      partition, for the pooled model and regime by regime;
##   D. runs the same Granger test pump by pump by least squares, and combines
##      it across pumps as Dumitrescu and Hurlin (2012), the fully
##      heterogeneous benchmark;
##   E. writes the tables and a map per province and direction, and a summary.
##
## The stages are cached in <out>/cache, so a run can be stopped and resumed.
##
##     Rscript run-fuels.R                          the provinces of the SETUP
##     Rscript run-fuels.R --provinces=BO,PR --p=3  other provinces, other lags
##
## This is part of the REPLICATION MATERIAL, not of the Stem package.
## ===========================================================================


## ---------------------------------------------------------------------------
## Where this script is
## ---------------------------------------------------------------------------
fu_here <- local({
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  f <- if (length(m)) normalizePath(sub("^--file=", "", m[1]), winslash = "/") else {
    of <- NULL
    for (i in rev(seq_len(sys.nframe()))) {
      of <- sys.frame(i)$ofile
      if (!is.null(of)) break
    }
    if (is.null(of)) NULL else normalizePath(of, winslash = "/")
  }
  if (is.null(f)) normalizePath(getwd(), winslash = "/") else dirname(f)
})

## command-line overrides, --name=value, lists separated by commas
fu_config <- function(defaults, args = commandArgs(trailingOnly = TRUE)) {
  for (a in grep("^--[^=]+=.", args, value = TRUE)) {
    nm <- sub("^--([^=]+)=.*$", "\\1", a)
    if (!nm %in% names(defaults)) stop("unknown option --", nm, call. = FALSE)
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
CFG <- fu_config(list(
  ## the station-level data, as received
  data       = file.path(fu_here, "App_FuelsITA", "station_level.zip"),

  ## provinces (NUTS-3 codes) and window. Two to four provinces of different
  ## structure: a metropolitan one, one of the Po Valley, a rural one of the
  ## South. The window starts on a Monday and covers the fiscal breaks below.
  provinces  = c("BO", "PR", "CS"),
  from       = "2021-01-04",
  to         = "2026-06-28",

  ## weekly price of a pump: "last", the price in force on the last day of
  ## the week (point sampling keeps the timing that a Granger test reads), or
  ## "mean", the mean of the week (averaging blurs it)
  weekly     = "last",

  ## pumps kept: not on highways, both self-service prices present in at
  ## least this share of the weeks of the window; gaps carried forward
  min_cover  = 0.95,

  ## lags of y and of x, and the error-correction term
  p          = 2L,
  ecm        = TRUE,

  ## the grid of the SC-STEM fit: few regimes
  k_max      = 3L,
  knn        = 5L,
  ## penalty on the coefficients: 0, since the likelihood-ratio test assumes
  ## the unpenalized estimator
  lambda     = 0,

  ## the station-by-station benchmark
  stations   = TRUE,

  out        = file.path(fu_here, "fuels")
))

## The changes of the excise duties in the window (to be checked against the
## decrees). Each enters as a pulse in the weekly changes: one week with
## weekly = "last", the week and the next with weekly = "mean".
FU_EVENTS <- data.frame(
  date  = as.Date(c("2022-03-22", "2022-12-01", "2023-01-01", "2025-05-15", "2026-01-01")),
  label = c("excise cut", "cut reduced", "cut ends", "realignment", "equal excise"),
  stringsAsFactors = FALSE)

## the two directions: y the fuel explained, x the fuel whose lags are tested
FU_DIRECTIONS <- list(
  diesel_to_gasoline = c(y = "price_gasoline_self", x = "price_diesel_self"),
  gasoline_to_diesel = c(y = "price_diesel_self",   x = "price_gasoline_self"))


## ===========================================================================
## Helpers
## ===========================================================================
if (!requireNamespace("Stem", quietly = TRUE))
  stop("Stem is not installed: remotes::install_github(\"PaoloMaranzano/Stem\")", call. = FALSE)
OUT   <- normalizePath(CFG$out[1], winslash = "/", mustWork = FALSE)
CACHE <- file.path(OUT, "cache")
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
FROM  <- as.Date(CFG$from[1]); TO <- as.Date(CFG$to[1])
if (format(FROM, "%u") != "1") stop("'from' must be a Monday", call. = FALSE)

mode_chr <- function(v) { v <- v[!is.na(v)]; if (!length(v)) NA_character_ else names(which.max(table(v))) }
locf <- function(v) {                       # carry forward, then backward at the start
  ok <- !is.na(v)
  if (!any(ok)) return(v)
  idx <- cummax(ifelse(ok, seq_along(v), 0L))
  idx[idx == 0L] <- which(ok)[1]
  v[idx]
}
cached <- function(name, expr) {
  f <- file.path(CACHE, paste0(name, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  val <- force(expr)
  saveRDS(val, f)
  val
}


## ===========================================================================
## A. The weekly panel of a province
## ===========================================================================
## the parts of the data: `data` is the zip as received, unzipped once per
## session into a temporary folder, or a folder where it was already unzipped
fu_parts <- function() {
  if (!file.exists(CFG$data[1])) stop("data not found: ", CFG$data[1], call. = FALSE)
  dir <- if (dir.exists(CFG$data[1])) CFG$data[1] else {
    dz <- file.path(tempdir(), "fuels_station_level")
    if (!dir.exists(dz)) {
      message("unzipping the station-level data ...")
      utils::unzip(CFG$data[1], exdir = dz)
    }
    dz
  }
  f <- list.files(dir, "^part_[0-9]+[.]RDS$", recursive = TRUE, full.names = TRUE)
  f[order(as.integer(gsub("\\D", "", basename(f))))]
}

## the daily rows of the requested provinces, read part by part
fu_daily <- function(provinces) {
  cols <- c("id_pump", "price_gasoline_self", "price_diesel_self", "date", "brand",
            "station_type", "city", "province", "latitude", "longitude")
  out <- lapply(fu_parts(), function(p) {
    x <- as.data.frame(readRDS(p))[, cols]
    ## the code of Naples, "NA", was read as a missing value
    x$province[is.na(x$province)] <- "NA"
    x[x$province %in% provinces & x$date >= FROM & x$date <= TO, , drop = FALSE]
  })
  do.call(rbind, out)
}

fu_panel <- function(x) {
  W  <- as.integer(floor(as.numeric(TO - FROM) / 7)) + 1L
  x$station_type[is.na(x$station_type)] <- ""
  x  <- x[x$station_type != "Autostradale", , drop = FALSE]
  x  <- x[order(x$id_pump, x$date), ]        # in date order within a pump, for "last"
  wk <- factor(as.integer(floor(as.numeric(x$date - FROM) / 7)) + 1L, levels = seq_len(W))
  pick <- if (CFG$weekly[1] == "mean") function(v) mean(v) else function(v) v[length(v)]
  weekly <- function(col) {
    v <- x[[col]]; ok <- is.finite(v) & v > 0.5 & v < 3.5
    m <- tapply(v[ok], list(factor(x$id_pump[ok]), wk[ok]), pick)
    100 * m                                  # cents per litre
  }
  G <- weekly("price_gasoline_self"); D <- weekly("price_diesel_self")
  ids <- intersect(rownames(G), rownames(D))
  G <- G[ids, , drop = FALSE]; D <- D[ids, , drop = FALSE]
  keep <- rowMeans(!is.na(G)) >= CFG$min_cover[1] & rowMeans(!is.na(D)) >= CFG$min_cover[1]
  G <- t(apply(G[keep, , drop = FALSE], 1, locf))
  D <- t(apply(D[keep, , drop = FALSE], 1, locf))
  xk <- x[as.character(x$id_pump) %in% rownames(G), , drop = FALSE]
  meta <- do.call(rbind, lapply(split(xk, xk$id_pump),
    function(z) data.frame(id_pump = z$id_pump[1], brand = mode_chr(z$brand),
                           type = mode_chr(z$station_type), city = mode_chr(z$city),
                           lon = stats::median(z$longitude, na.rm = TRUE),
                           lat = stats::median(z$latitude, na.rm = TRUE),
                           stringsAsFactors = FALSE)))
  meta <- meta[match(rownames(G), as.character(meta$id_pump)), ]
  ok <- is.finite(meta$lon) & is.finite(meta$lat)
  list(G = G[ok, , drop = FALSE], D = D[ok, , drop = FALSE], meta = meta[ok, ],
       weeks = FROM + 7 * (seq_len(W) - 1L))
}


## ===========================================================================
## B. The model in weekly changes
## ===========================================================================
## pulses at the week of each fiscal event (and the next, with weekly means);
## `weeks` are the weeks the changes arrive in, and the change into the week of
## an event is the one that carries it
fu_pulses <- function(weeks) {
  w <- as.integer(floor(as.numeric(FU_EVENTS$date - weeks[1]) / 7)) + 1L
  shift <- if (CFG$weekly[1] == "mean") 0:1 else 0L
  P <- do.call(cbind, lapply(seq_along(w), function(e) vapply(shift, function(s)
    as.numeric(seq_along(weeks) == w[e] + s), numeric(length(weeks)))))
  colnames(P) <- paste0("pulse_", rep(gsub(" ", "_", FU_EVENTS$label), each = length(shift)),
                        if (length(shift) > 1) rep(c("", "_next"), nrow(FU_EVENTS)) else "")
  P[, colSums(P) > 0, drop = FALSE]          # events outside the window drop out
}

## the fiscal period of each week, for the error-correction term
fu_period <- function(weeks) findInterval(weeks, sort(FU_EVENTS$date)) + 1L

fu_design <- function(panel, dir, with_x = TRUE) {
  Y <- if (dir[["y"]] == "price_gasoline_self") panel$G else panel$D
  X <- if (dir[["x"]] == "price_gasoline_self") panel$G else panel$D
  p <- CFG$p[1]
  dY <- Y[, -1, drop = FALSE] - Y[, -ncol(Y), drop = FALSE]   # week t = 2..W
  dX <- X[, -1, drop = FALSE] - X[, -ncol(X), drop = FALSE]
  wk <- panel$weeks[-1]
  tt <- (p + 1):ncol(dY)                                       # usable weeks
  d <- nrow(Y); Tn <- length(tt)
  ## the spread of the previous week, demeaned by pump and fiscal period
  S <- Y - X
  per <- fu_period(panel$weeks)
  Sd <- S
  for (q in unique(per)) Sd[, per == q] <- S[, per == q] - rowMeans(S[, per == q, drop = FALSE])
  pul <- fu_pulses(wk)[tt, , drop = FALSE]
  pul <- pul[, colSums(pul) > 0, drop = FALSE]                 # events in the first p weeks
  cov_i <- function(i) {
    m <- cbind(intercept = 1,
               sapply(seq_len(p), function(l) dY[i, tt - l]),
               if (with_x) sapply(seq_len(p), function(l) dX[i, tt - l]),
               if (CFG$ecm[1]) Sd[i, tt],          # spread at week t-1: S has one more column
               pul)
    colnames(m) <- c("intercept", paste0("dy_lag", seq_len(p)),
                     if (with_x) paste0("dx_lag", seq_len(p)),
                     if (CFG$ecm[1]) "ecm", colnames(pul))
    m
  }
  covs <- lapply(seq_len(d), cov_i)
  list(z = t(dY[, tt, drop = FALSE]),                          # T x d
       covariates = do.call(rbind, covs),                      # stacked by location
       names = colnames(covs[[1]]), weeks = wk[tt], d = d, Tn = Tn)
}

fu_model <- function(des, coords) {
  X <- des$covariates; y <- as.vector(des$z)                   # z is T x d: stacked by location
  ols <- stats::lm.fit(X, y)
  s2 <- stats::var(ols$residuals)
  dm <- stats::median(as.matrix(geodist::geodist(coords, measure = "geodesic")))
  Stem::STEM_Model(z = des$z, covariates = X, coordinates = coords,
                   phi = list(beta = matrix(ols$coefficients, ncol = 1),
                              sigma2eps = 0.6 * s2, sigma2omega = 0.4 * s2,
                              theta = 3 / dm, G = matrix(0.3, 1, 1),
                              Sigmaeta = matrix(0.2 * s2, 1, 1),
                              m0 = as.matrix(0), C0 = as.matrix(1)),
                   K = matrix(1, des$d, 1))
}


## ===========================================================================
## C. SC-STEM against the pooled model, and the Granger tests
## ===========================================================================
fu_fit <- function(panel, dir) {
  coords <- cbind(lon = panel$meta$lon, lat = panel$meta$lat)
  ## two pumps at the same point (one site under two codes) would make the
  ## neighbour graph and the covariance degenerate: move the later one by
  ## about a metre
  dup <- duplicated(round(coords, 6))
  coords[dup, ] <- coords[dup, ] + 1e-5 * seq_len(sum(dup))
  des  <- fu_design(panel, dir, with_x = TRUE)
  desR <- fu_design(panel, dir, with_x = FALSE)
  mod  <- fu_model(des, coords); modR <- fu_model(desR, coords)
  r <- length(des$names)
  kmax <- min(CFG$k_max[1], floor(des$d / (r + 2)))            # regimes of at least r + 2 pumps
  ic <- Stem::SCSTEM_Infocrit(mod, k_grid = seq_len(kmax), distance = "geo",
                              knn = CFG$knn[1], lambda = CFG$lambda[1], verbose = FALSE)
  sel <- Stem::SCSTEM_Select(ic)
  pooled <- ic$fits[[paste0("k=1, phi=", ic$table$phi[ic$table$k == 1][1])]]
  ## the restricted fits, without the lags of x, on the same partitions
  refit <- function(fit, k, phi) {
    if (k == 1L) Stem::SCSTEM_Estimation(modR, k = 1, distance = "geo", knn = CFG$knn[1],
                                         lambda = CFG$lambda[1], verbose = FALSE)
    else Stem::SCSTEM_Estimation(modR, k = k, phi_penalty = phi, distance = "geo",
                                 knn = CFG$knn[1], lambda = CFG$lambda[1],
                                 init_partition = fit$group, max_iter = 0, verbose = FALSE)
  }
  pooledR <- refit(pooled, 1L, 0)
  selR <- if (sel$k_selected == 1L) pooledR else refit(sel$fit, sel$k_selected, sel$phi_selected)
  list(des = des, ic = ic, sel = sel, pooled = pooled, pooledR = pooledR, selR = selR)
}

## the likelihood-ratio test of no Granger causality, regime by regime and in
## total, and the coefficients of the lags of x
fu_granger <- function(fit, fitR, names, model) {
  p <- CFG$p[1]
  jx <- grep("^dx_lag", names)
  k <- nrow(fit$phi_hat)
  lr <- 2 * (fit$loglik_g - fitR$loglik_g)
  b <- fit$phi_hat[, paste0("beta", jx), drop = FALSE]
  colnames(b) <- names[jx]
  out <- data.frame(model = model, regime = seq_len(k),
                    pumps = as.integer(table(factor(fit$group, levels = seq_len(k)))),
                    b, LR = lr, df = p,
                    p_value = stats::pchisq(pmax(lr, 0), df = p, lower.tail = FALSE),
                    stringsAsFactors = FALSE)
  rbind(out, data.frame(model = model, regime = 0L, pumps = sum(out$pumps),
                        b[1, , drop = FALSE] * NA, LR = sum(lr), df = p * k,
                        p_value = stats::pchisq(pmax(sum(lr), 0), df = p * k, lower.tail = FALSE)))
}


## ===========================================================================
## D. The station-by-station benchmark
## ===========================================================================
fu_stations <- function(des) {
  p <- CFG$p[1]
  jx <- grep("^dx_lag", des$names)
  Tn <- des$Tn
  res <- t(vapply(seq_len(des$d), function(i) {
    rows <- (i - 1L) * Tn + seq_len(Tn)
    X <- des$covariates[rows, , drop = FALSE]; y <- des$z[, i]
    X <- X[, colSums(abs(X)) > 0, drop = FALSE]                # pulses with no variation
    u <- stats::lm.fit(X, y); r0 <- stats::lm.fit(X[, setdiff(colnames(X), des$names[jx]), drop = FALSE], y)
    rss1 <- sum(u$residuals^2); rss0 <- sum(r0$residuals^2)
    dfr <- Tn - ncol(X)
    ## a pump whose price never moved has no test
    if (rss1 <= 1e-10 || u$rank < ncol(X)) return(c(F = NA, p_value = NA, W = NA))
    F <- ((rss0 - rss1) / p) / (rss1 / dfr)
    c(F = F, p_value = stats::pf(F, p, dfr, lower.tail = FALSE), W = p * F)
  }, numeric(3)))
  ok <- is.finite(res[, "W"])
  N <- sum(ok)
  zbar <- sqrt(N / (2 * p)) * (mean(res[ok, "W"]) - p)
  list(pumps = as.data.frame(res),
       dh = c(N = N, Wbar = mean(res[ok, "W"]), Zbar = zbar,
              p_value = stats::pnorm(zbar, lower.tail = FALSE),
              share_rejecting = mean(res[ok, "p_value"] < 0.05)))
}


## ===========================================================================
## E. Run and write
## ===========================================================================
fu_map <- function(panel, fit, st, file, title) {
  grDevices::cairo_pdf(file, width = 6.4, height = 6)
  graphics::par(mar = c(3, 3, 2.4, 1))
  cols <- c("#1f6f8b", "#e0a458", "#a8516e", "#5b8c5a")[fit$group]
  sig <- if (is.null(st)) rep(FALSE, length(fit$group)) else st$pumps$p_value < 0.05
  graphics::plot(panel$meta$lon, panel$meta$lat, asp = 1 / cos(mean(panel$meta$lat) * pi / 180),
                 pch = ifelse(sig, 19, 21), bg = "white", col = cols, cex = 0.9,
                 xlab = "", ylab = "", las = 1, bty = "n")
  graphics::legend("topright", bty = "n", cex = 0.75, pch = c(rep(19, max(fit$group)), 21),
                   col = c(unique(cols[order(fit$group)]), "grey40"),
                   legend = c(paste("regime", seq_len(max(fit$group))),
                              "open: pump-level test not rejected"))
  graphics::title(title, cex.main = 0.9)
  invisible(grDevices::dev.off())
}

summary_rows <- list()
for (prov in CFG$provinces) {
  panel <- cached(paste0("panel_", prov, "_", CFG$from, "_", CFG$to, "_", CFG$weekly),
                  fu_panel(fu_daily(prov)))
  cat(sprintf("%s: %d pumps, %d weeks\n", prov, nrow(panel$G), ncol(panel$G)))
  if (nrow(panel$G) < 20) { cat("  too few pumps, skipped\n"); next }
  for (dn in names(FU_DIRECTIONS)) {
    tag <- paste(prov, dn, sep = "_")
    fit <- cached(paste0("fit_", tag, "_p", CFG$p, "_ecm", CFG$ecm, "_k", CFG$k_max),
                  fu_fit(panel, FU_DIRECTIONS[[dn]]))
    g <- rbind(fu_granger(fit$pooled, fit$pooledR, fit$des$names, "pooled"),
               if (fit$sel$k_selected > 1)
                 fu_granger(fit$sel$fit, fit$selR, fit$des$names, "SC-STEM"))
    utils::write.csv(fit$ic$table, file.path(OUT, paste0(tag, "_grid.csv")), row.names = FALSE)
    utils::write.csv(g, file.path(OUT, paste0(tag, "_granger.csv")), row.names = FALSE)
    st <- if (isTRUE(CFG$stations[1])) fu_stations(fit$des) else NULL
    if (!is.null(st))
      utils::write.csv(cbind(panel$meta, regime = fit$sel$fit$group, st$pumps),
                       file.path(OUT, paste0(tag, "_stations.csv")), row.names = FALSE)
    fu_map(panel, fit$sel$fit, st, file.path(OUT, paste0(tag, "_map.pdf")),
           sprintf("%s, %s: k = %d, phi = %s", prov, gsub("_", " ", dn),
                   fit$sel$k_selected, format(fit$sel$phi_selected)))
    bic <- function(f) unname(f$info_crit[["BIC"]])
    summary_rows[[tag]] <- data.frame(
      province = prov, direction = dn, pumps = fit$des$d, weeks = fit$des$Tn,
      k = fit$sel$k_selected, phi = fit$sel$phi_selected,
      BIC_pooled = bic(fit$pooled), BIC_selected = bic(fit$sel$fit),
      p_pooled = g$p_value[g$model == "pooled" & g$regime == 0],
      p_regimes = if (fit$sel$k_selected > 1)
        paste(format.pval(g$p_value[g$model == "SC-STEM" & g$regime > 0], digits = 2), collapse = "; ")
        else NA_character_,
      DH_Zbar = if (is.null(st)) NA_real_ else st$dh[["Zbar"]],
      DH_p = if (is.null(st)) NA_real_ else st$dh[["p_value"]],
      pumps_rejecting = if (is.null(st)) NA_real_ else st$dh[["share_rejecting"]],
      stringsAsFactors = FALSE)
    cat(sprintf("  %-20s k = %d, phi = %s | pooled p = %s | DH p = %s\n", dn,
                fit$sel$k_selected, format(fit$sel$phi_selected),
                format.pval(summary_rows[[tag]]$p_pooled, digits = 2),
                format.pval(summary_rows[[tag]]$DH_p, digits = 2)))
  }
}
if (length(summary_rows)) {
  s <- do.call(rbind, summary_rows); rownames(s) <- NULL
  utils::write.csv(s, file.path(OUT, "fuels_summary.csv"), row.names = FALSE)
  cat("\nsummary written to", file.path(OUT, "fuels_summary.csv"), "\n")
}
