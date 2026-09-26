## ===========================================================================
## SC-STEM: the simulation study. THE ANALYSIS.
##
## One script, self-contained. It needs nothing but an installed Stem package
## and its companion run-simulations.R in the same folder, from which it takes
## the generator and the design -- so the two cannot disagree about either. It
## runs from whatever folder it sits in.
##
##     Rscript analyse-simulations.R                    everything it can
##     Rscript analyse-simulations.R --parts=design     the design only
##     Rscript analyse-simulations.R --tag=pilot        another set of results
##     Rscript analyse-simulations.R --tag=p1,p2,p3     several, stacked: the
##                                                      processes of a parallel run
##     Rscript analyse-simulations.R --out=D:/paper/Figures
##
## Two parts.
##
##   DESIGN    what the study is, before anything is run: the configurations,
##             the neighbourhood graph, one replication of each kind, and the
##             tables of the scenarios and the overlap levels. Needs no results.
##   RESULTS   what the study found, read from the outputs of run-simulations.R.
##             Every block is guarded: it produces what the results currently
##             support and says what is missing, so it can be run while the
##             study is still in progress.
##
## Inputs  : <here>/results/<tag>*.csv        (change with --results=)
## Outputs : <here>/output/                   (change with --out=)
##
## This is part of the REPLICATION MATERIAL of the paper, not of the Stem
## package.
## ===========================================================================


## ---------------------------------------------------------------------------
## Where this script is
## ---------------------------------------------------------------------------
ana_here <- local({
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

## ---------------------------------------------------------------------------
## The generator and the design, from the runner beside this script
## ---------------------------------------------------------------------------
runner <- file.path(ana_here, "run-simulations.R")
if (!file.exists(runner))
  stop("run-simulations.R must sit in the same folder as this script:\n  ",
       ana_here, call. = FALSE)
SIM_DEFINE_ONLY <- TRUE
source(runner, local = globalenv())
## The neighbour graph as the estimator builds it. On the smallest networks of
## the design spdep warns that k exceeds a third of the locations or that the
## graph splits into components; both are true and both are the point of those
## panels, so the warnings are not repeated forty times over.
neighbors <- function(...)
  suppressWarnings(utils::getFromNamespace("scstem_neighbors", "Stem")(...))

## ---------------------------------------------------------------------------
## Options
## ---------------------------------------------------------------------------
AN <- sim_config(list(
  parts   = c("design", "results"),
  tag     = "full",
  results = file.path(ana_here, "results"),
  out     = file.path(ana_here, "output")
))
OUT <- AN$out[1]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
OUT <- normalizePath(OUT, winslash = "/")
RES <- normalizePath(AN$results[1], winslash = "/", mustWork = FALSE)
cat("analysis\n  results from ", RES, "\n  outputs to   ", OUT, "\n\n", sep = "")

## ---------------------------------------------------------------------------
## Style, shared by every figure
## ---------------------------------------------------------------------------
COL   <- c("#1f6f8b", "#e0a458", "#5b8c5a", "#a8516e")
GREY  <- "#8a939f"
INK   <- "#1b2430"
EDGE  <- "#9aa4b0"
CROSS <- "#c0392b"
PAL3  <- c("#1f6f8b", "#e0a458", "#a8516e")

KTRUE   <- SIM_REFERENCE$K
N_GRID  <- SIM_LEVELS$n
OM_GRID <- SIM_LEVELS$omega
KNN     <- SIM_LEVELS$knn
om_lab  <- function(w) sprintf("omega = %s", format(w))

pdf_open <- function(name, width, height) {
  grDevices::cairo_pdf(file.path(OUT, name), width = width, height = height)
}
pdf_close <- function(name) {
  invisible(grDevices::dev.off())
  message("wrote ", name)
}

## a row label in the outer margin, right-aligned against the panels
row_label <- function(top, bottom, main, sub) {
  y <- (top + bottom) / 2
  graphics::mtext(main, side = 2, outer = TRUE, at = y, line = 2.4,
                  las = 1, adj = 1, cex = 1.0, font = 2, col = INK)
  graphics::mtext(sub, side = 2, outer = TRUE, at = y - 0.035, line = 2.4,
                  las = 1, adj = 1, cex = 0.72, col = GREY)
}

draw_loc <- function(n, K, omega, seed = 1, balance = "balanced") {
  loc <- dgp_locations(n, K, omega, balance = balance, seed = seed)
  list(xy = loc$xy, g = loc$labels, mu = loc$mu, coords = loc$coords)
}

cloud_panel <- function(s, lim, cex_pt = 0.9, centres = TRUE) {
  graphics::plot(NA, xlim = lim, ylim = lim, asp = 1, xlab = "", ylab = "",
                 xaxt = "n", yaxt = "n", bty = "n")
  graphics::rect(lim[1], lim[1], lim[2], lim[2], border = "#dfe4ea", lwd = 0.8)
  graphics::abline(h = 0, v = 0, col = "#eef1f4", lwd = 0.8)
  if (centres) graphics::points(s$mu[, 1], s$mu[, 2], pch = 3, col = "#5a626d",
                                cex = 1.3, lwd = 1.6)
  graphics::points(s$xy[, 1], s$xy[, 2], pch = 21, cex = cex_pt, lwd = 0.4,
                   bg = grDevices::adjustcolor(COL[s$g], 0.9), col = "white")
}

tex_write <- function(lines, name) {
  writeLines(lines, file.path(OUT, name))
  message("wrote ", name)
}
fmt <- function(x, d = 2) ifelse(is.na(x), "---", formatC(x, format = "f", digits = d))


## ===========================================================================
## PART A. THE DESIGN
## ===========================================================================
if ("design" %in% AN$parts) {

  cat("DESIGN\n")
  lim <- c(-2.7, 2.7)

  ## -------------------------------------------------------------------------
  ## A1. The overlap, across the network sizes of the design
  ## -------------------------------------------------------------------------
  pdf_open("fig_design_omega.pdf", 11.0, 7.4)
  graphics::par(mfrow = c(length(OM_GRID), length(N_GRID)),
                mar = c(0.4, 0.4, 0.4, 0.4), oma = c(3.2, 12.5, 4.2, 0.8))
  nr <- length(OM_GRID)
  for (wi in seq_along(OM_GRID)) {
    w <- OM_GRID[wi]
    for (n in N_GRID) {
      s <- draw_loc(n, KTRUE, w, seed = 100 + n)
      cloud_panel(s, lim, cex_pt = if (n > 200) 0.55 else if (n > 100) 0.75 else
                                   if (n > 50) 0.95 else 1.25)
      if (wi == 1L)
        graphics::mtext(paste0("n = ", n), side = 3, line = 0.8, cex = 0.95, font = 2)
      if (n == N_GRID[1])
        row_label(1 - (wi - 1) / nr, 1 - wi / nr, om_lab(w),
                  sprintf("%.2f sd apart", dgp_separation(KTRUE, w)))
    }
  }
  graphics::mtext(sprintf("K = %d regimes: the overlap parameter across the network sizes of the design", KTRUE),
                  outer = TRUE, side = 3, line = 2.4, cex = 1.15, font = 2, col = INK)
  graphics::mtext(paste0("centres 2*omega apart, total spatial variance held fixed at ",
                         round(NU_TOT, 3), "; crosses are the true centres, colour the true regime"),
                  outer = TRUE, side = 1, line = 1.2, cex = 0.78, col = GREY)
  pdf_close("fig_design_omega.pdf")

  ## -------------------------------------------------------------------------
  ## A2. The graph the penalty lives on
  ##
  ## knn does not change the data -- coordinates, covariate and response are
  ## the same whatever the graph -- it changes which pairs the penalty rewards
  ## for agreeing. An edge that joins two different true regimes is one the
  ## penalty rewards when it should not, so the edges are coloured by that.
  ## -------------------------------------------------------------------------
  graph_panel <- function(s, knn, lim, cex_pt = 0.8) {
    W  <- neighbors(s$coords, knn = knn, distance = "geo")$W
    graphics::plot(NA, xlim = lim, ylim = lim, asp = 1, xlab = "", ylab = "",
                   xaxt = "n", yaxt = "n", bty = "n")
    graphics::rect(lim[1], lim[1], lim[2], lim[2], border = "#dfe4ea", lwd = 0.8)
    ij <- which(upper.tri(W) & W > 0, arr.ind = TRUE)
    cross <- s$g[ij[, 1]] != s$g[ij[, 2]]
    graphics::segments(s$xy[ij[!cross, 1], 1], s$xy[ij[!cross, 1], 2],
                       s$xy[ij[!cross, 2], 1], s$xy[ij[!cross, 2], 2],
                       col = EDGE, lwd = 0.5)
    graphics::segments(s$xy[ij[cross, 1], 1], s$xy[ij[cross, 1], 2],
                       s$xy[ij[cross, 2], 1], s$xy[ij[cross, 2], 2],
                       col = grDevices::adjustcolor(CROSS, 0.75), lwd = 0.8)
    graphics::points(s$xy[, 1], s$xy[, 2], pch = 21, cex = cex_pt, lwd = 0.4,
                     bg = grDevices::adjustcolor(COL[s$g], 0.9), col = "white")
    c(deg = mean(rowSums(W > 0)), cross = mean(cross))
  }

  NG <- SIM_REFERENCE$n
  pdf_open("fig_design_knn.pdf", 9.8, 9.8)
  graphics::par(mfrow = c(length(KNN), length(OM_GRID)),
                mar = c(1.4, 0.4, 0.4, 0.4), oma = c(4.0, 12.0, 4.2, 0.8))
  nr <- length(KNN)
  for (ki in seq_along(KNN)) {
    for (w in OM_GRID) {
      s  <- draw_loc(NG, KTRUE, w, seed = 100 + NG)
      st <- graph_panel(s, KNN[ki], lim)
      if (ki == 1L)
        graphics::mtext(om_lab(w), side = 3, line = 0.8, cex = 1.0, font = 2)
      if (identical(w, OM_GRID[1]))
        row_label(1 - (ki - 1) / nr, 1 - ki / nr, paste0("knn = ", KNN[ki]),
                  sprintf("mean degree %.1f", st[["deg"]]))
      graphics::mtext(sprintf("%.0f%% of edges cross", 100 * st[["cross"]]),
                      side = 1, line = 0.3, cex = 0.8, col = CROSS)
    }
  }
  graphics::mtext(sprintf("The neighbourhood graph of the Potts penalty, n = %d, K = %d",
                          NG, KTRUE),
                  outer = TRUE, side = 3, line = 2.4, cex = 1.15, font = 2, col = INK)
  graphics::mtext("grey edges join two locations of the same true regime; red edges cross a regime boundary,",
                  outer = TRUE, side = 1, line = 1.3, cex = 0.76, col = GREY)
  graphics::mtext("and are the pairs the penalty rewards for agreeing when they should not",
                  outer = TRUE, side = 1, line = 2.2, cex = 0.76, col = GREY)
  pdf_close("fig_design_knn.pdf")

  ## -------------------------------------------------------------------------
  ## A3. The graph in numbers, over the whole design
  ## -------------------------------------------------------------------------
  gtab <- do.call(rbind, lapply(N_GRID, function(n)
    do.call(rbind, lapply(KNN, function(kk)
      do.call(rbind, lapply(OM_GRID, function(w) {
        s  <- draw_loc(n, KTRUE, w, seed = 100 + n)
        W  <- neighbors(s$coords, knn = kk, distance = "geo")$W
        ij <- which(upper.tri(W) & W > 0, arr.ind = TRUE)
        dm <- as.matrix(geodist::geodist(s$coords, measure = "geodesic")) / 1000
        data.frame(n = n, knn = kk, omega = w,
                   degree = mean(rowSums(W > 0)),
                   cross  = mean(s$g[ij[, 1]] != s$g[ij[, 2]]),
                   len    = stats::median(dm[ij]))
      }))))))

  pdf_open("fig_design_graph.pdf", 11.4, 4.2)
  graphics::layout(matrix(1:4, nrow = 1), widths = c(1, 1, 1, 0.48))
  graphics::par(mar = c(4.2, 4.6, 2.8, 0.8), oma = c(1.0, 0, 2.4, 0))
  for (v in c("degree", "cross", "len")) {
    ttl <- switch(v, degree = "mean degree",
                  cross = "share of edges crossing a regime",
                  len = "median edge length (km)")
    ylim <- range(gtab[[v]]); if (v == "cross") ylim <- c(0, max(ylim))
    graphics::plot(NA, xlim = range(N_GRID), ylim = ylim, log = "x",
                   xlab = "", ylab = "", xaxt = "n", bty = "n",
                   col.axis = "#3d4753", las = 1)
    graphics::axis(1, at = N_GRID, labels = N_GRID, col.axis = "#3d4753")
    graphics::mtext("number of locations", side = 1, line = 2.5, cex = 0.8,
                    col = "#3d4753")
    ## for K balanced regimes a graph that knows nothing of the partition
    ## has a crossing share of 1 - 1/K
    if (v == "cross")
      graphics::abline(h = 1 - 1 / KTRUE, col = "#c9ced6", lty = "13", lwd = 1.2)
    for (ki in seq_along(KNN)) for (wi in seq_along(OM_GRID)) {
      s <- gtab[gtab$knn == KNN[ki] & gtab$omega == OM_GRID[wi], ]
      s <- s[order(s$n), ]
      graphics::lines(s$n, s[[v]], col = PAL3[ki], lwd = 1.7,
                      lty = c("solid", "22", "44")[wi])
      graphics::points(s$n, s[[v]], col = PAL3[ki], pch = 19, cex = 0.7)
    }
    if (v == "cross")
      graphics::text(N_GRID[length(N_GRID)], 1 - 1 / KTRUE, adj = c(1, 1.6),
                     cex = 0.72, col = GREY, labels = "no information: 1 - 1/K")
    graphics::mtext(ttl, side = 3, line = 0.8, cex = 0.88, font = 2, col = INK)
  }
  graphics::par(mar = c(4.2, 0.2, 2.8, 0.2))
  graphics::plot.new()
  graphics::legend("left", bty = "n", cex = 0.86,
                   legend = c(paste0("knn = ", KNN), "", paste0("omega = ", OM_GRID)),
                   col = c(PAL3, NA, rep("#3d4753", 3)),
                   lty = c(1, 1, 1, NA, 1, 2, 3), lwd = 1.7, seg.len = 1.8)
  graphics::mtext(sprintf("What the neighbourhood graph does, K = %d", KTRUE),
                  outer = TRUE, side = 3, line = 0.9, cex = 1.1, font = 2, col = INK)
  graphics::mtext(paste("knn multiplies the number of edges, but from n = 100 on it leaves unchanged the information",
                        "the graph holds about the partition; under the default penalty scale it also leaves its weight unchanged"),
                  outer = TRUE, side = 1, line = -0.2, cex = 0.76, col = GREY)
  pdf_close("fig_design_graph.pdf")
  cat("\nthe graph in numbers\n")
  print(gtab[order(gtab$knn, gtab$omega, gtab$n), ], row.names = FALSE, digits = 3)

  ## -------------------------------------------------------------------------
  ## A4. One replication of each kind
  ## -------------------------------------------------------------------------
  scen <- dgp_scenarios()
  wmid <- OM_GRID[2]; wmax <- max(OM_GRID)
  ex_rows <- list(
    list(K = 1L, omega = 0, id = "S0", bal = "balanced",
         lab = "K = 1 (pooled)", sub = "no regime at all"),
    list(K = KTRUE, omega = 0, id = "S4", bal = "balanced",
         lab = sprintf("K = %d, omega = 0", KTRUE),
         sub = "all three contrasts, no spatial separation"),
    list(K = KTRUE, omega = wmid, id = "S4", bal = "balanced",
         lab = sprintf("K = %d, omega = %s", KTRUE, format(wmid)),
         sub = sprintf("all three contrasts, %.2f sd apart", dgp_separation(KTRUE, wmid))),
    list(K = KTRUE, omega = wmax, id = "S4", bal = "balanced",
         lab = sprintf("K = %d, omega = %s", KTRUE, format(wmax)),
         sub = sprintf("all three contrasts, %.2f sd apart", dgp_separation(KTRUE, wmax))),
    list(K = KTRUE, omega = wmax, id = "S1a", bal = "balanced",
         lab = sprintf("K = %d, omega = %s", KTRUE, format(wmax)),
         sub = "separated in space, coefficients 0.25 sd apart"),
    list(K = KTRUE, omega = wmax, id = "S4", bal = "unbalanced",
         lab = sprintf("K = %d, omega = %s", KTRUE, format(wmax)),
         sub = "unbalanced regimes")
  )
  EX_N <- SIM_REFERENCE$n; EX_T <- 365L; NS <- 3L; WIN <- 120L

  ex <- lapply(ex_rows, function(r) {
    row <- scen[scen$id == r$id, , drop = FALSE]
    d   <- dgp_draw(EX_N, EX_T, r$K, r$omega, row, rep = 1L, balance = r$bal)
    loc <- dgp_locations(EX_N, r$K, r$omega, balance = r$bal, seed = 1000L + 1L)
    set.seed(11)
    pick <- unlist(lapply(seq_len(r$K), function(k)
      sample(which(d$labels == k), min(NS, sum(d$labels == k)))))
    list(z = d$z, y = d$latent, g = d$labels, xy = loc$xy, mu = loc$mu,
         pick = pick, sizes = d$sizes, lab = r$lab, sub = r$sub)
  })
  elim <- max(abs(unlist(lapply(ex, function(d) range(d$xy))))) * 1.05
  elim <- c(-elim, elim)
  ylim <- range(unlist(lapply(ex, function(d) d$z[seq_len(WIN), d$pick])))

  pdf_open("fig_dgp_examples.pdf", 9.2, 12.6)
  graphics::par(mfcol = c(length(ex), 2), mar = c(1.6, 1.0, 1.4, 0.6),
                oma = c(3.2, 6.4, 4.0, 1.0))
  for (d in ex) cloud_panel(list(xy = d$xy, g = d$g, mu = d$mu), elim, cex_pt = 0.7)
  for (d in ex) {
    tt <- seq_len(WIN)
    graphics::plot(NA, xlim = range(tt), ylim = ylim, xlab = "", ylab = "",
                   xaxt = "n", yaxt = "n", bty = "n")
    graphics::abline(h = 0, col = "#eef1f4", lwd = 0.8)
    for (i in d$pick)
      graphics::lines(tt, d$z[tt, i],
                      col = grDevices::adjustcolor(COL[d$g[i]], 0.22), lwd = 0.4)
    for (k in seq_len(ncol(d$y))) {
      idx <- which(d$g == k)
      graphics::lines(tt, d$y[tt, k] + mean(d$z[, idx]), col = COL[k], lwd = 2.0)
    }
    graphics::axis(1, at = pretty(tt, 4), cex.axis = 0.7, col = GREY,
                   col.axis = GREY, tck = -0.04, mgp = c(2, 0.35, 0))
    graphics::axis(2, at = pretty(ylim, 3), cex.axis = 0.7, col = GREY,
                   col.axis = GREY, las = 1, tck = -0.03, mgp = c(2, 0.45, 0))
  }
  h <- 1 / length(ex)
  for (i in seq_along(ex)) {
    yy <- 1 - (i - 0.5) * h
    graphics::mtext(ex[[i]]$lab, side = 2, outer = TRUE, at = yy, line = 4.4,
                    las = 1, cex = 0.92, font = 2, col = INK, adj = 0)
    graphics::mtext(ex[[i]]$sub, side = 2, outer = TRUE, at = yy - 0.030,
                    line = 4.4, las = 1, cex = 0.72, col = GREY, adj = 0)
    graphics::mtext(paste0("sizes ", paste(ex[[i]]$sizes, collapse = "/")),
                    side = 2, outer = TRUE, at = yy - 0.055, line = 4.4,
                    las = 1, cex = 0.68, col = GREY, adj = 0)
  }
  graphics::mtext("locations", side = 3, outer = TRUE, at = 0.30, line = 1.0,
                  cex = 1.0, font = 2, col = INK)
  graphics::mtext("the series they carry", side = 3, outer = TRUE, at = 0.76,
                  line = 1.0, cex = 1.0, font = 2, col = INK)
  graphics::mtext(sprintf("One replication of the design at n = %d, T = %d", EX_N, EX_T),
                  outer = TRUE, side = 3, line = 2.5, cex = 1.15, font = 2, col = INK)
  graphics::mtext("crosses mark the true centres; thin lines are three stations per regime, over the first 120 periods",
                  outer = TRUE, side = 1, line = 0.4, cex = 0.74, col = GREY)
  graphics::mtext("thick lines are the latent process of each regime, shifted to the regime mean",
                  outer = TRUE, side = 1, line = 1.4, cex = 0.74, col = GREY)
  pdf_close("fig_dgp_examples.pdf")

  ## -------------------------------------------------------------------------
  ## A5. The scenarios: what each sets in the last regime
  ## -------------------------------------------------------------------------
  b   <- dgp_base()
  tot <- b$sigma2eps + b$sigma2omega
  srows <- do.call(rbind, lapply(seq_len(nrow(scen)), function(i) {
    s <- scen[i, ]
    p <- dgp_psi(s$scenario, s$level, K = 2L)[[2]]
    data.frame(id = s$id, label = s$label,
               beta1 = fmt(p$beta[2]),
               d_beta = fmt((p$beta[2] - b$beta[2]) / sqrt(tot)),
               G = fmt(p$G), nugget = fmt(p$sigma2eps / tot),
               range = fmt(1 / p$theta / 1000, 0), rho = fmt(s$rho, 0),
               stringsAsFactors = FALSE)
  }))
  tex_write(c(
    "\\begin{tabular}{llrrrrrr}",
    "\\toprule",
    "& & \\multicolumn{6}{c}{parameters of the last regime} \\\\",
    "\\cmidrule(l){3-8}",
    "id & separation & $\\beta_1$ & $\\Delta\\beta_1/\\sigma$ & $G$ & nugget share & range (km) & $\\rho$ \\\\",
    "\\midrule",
    sprintf("%s & %s & %s & %s & %s & %s & %s & %s \\\\", srows$id, srows$label,
            srows$beta1, srows$d_beta, srows$G, srows$nugget, srows$range, srows$rho),
    "\\bottomrule",
    "\\end{tabular}"), "tab_dgp.tex")

  ## -------------------------------------------------------------------------
  ## A6. The overlap levels, read on the scale that matters
  ##
  ## The Bayes error of the optimal rule that assigns a location to the nearest
  ## centre, and the Adjusted Rand Index that rule attains: what ANY purely
  ## spatial procedure could extract from the configuration.
  ## -------------------------------------------------------------------------
  ari <- utils::getFromNamespace("scstem_ari", "Stem")
  otab <- do.call(rbind, lapply(OM_GRID, function(w) {
    nu <- dgp_nu_sp(KTRUE, w); mu <- dgp_centres(KTRUE, w)
    set.seed(7); M <- 200000L
    g  <- sample.int(KTRUE, M, TRUE)
    xy <- cbind(mu[g, 1] + stats::rnorm(M, sd = sqrt(nu)),
                mu[g, 2] + stats::rnorm(M, sd = sqrt(nu)))
    d2 <- sapply(seq_len(KTRUE), function(j) (xy[, 1] - mu[j, 1])^2 + (xy[, 2] - mu[j, 2])^2)
    gh <- max.col(-d2)
    data.frame(omega = w, nu_sp = nu, sd_km = sqrt(nu) * UNIT_KM,
               centre_km = 2 * w * UNIT_KM, sep = dgp_separation(KTRUE, w),
               bayes = mean(gh != g), ari = ari(gh[1:20000], g[1:20000]))
  }))
  cat("\nthe overlap levels\n"); print(otab, row.names = FALSE, digits = 3)
  tex_write(c(
    "\\begin{tabular}{rrrrrrr}",
    "\\toprule",
    "$\\omega$ & $\\nu_{sp}$ & cluster sd & centre distance & $\\delta$ & Bayes error & attainable ARI \\\\",
    "\\midrule",
    sprintf("$%s$ & $%s$ & $%s$ km & $%s$ km & $%s$ & $%s$ & $%s$ \\\\",
            fmt(otab$omega, 2), fmt(otab$nu_sp, 3), fmt(otab$sd_km, 0),
            fmt(otab$centre_km, 0), fmt(otab$sep, 2), fmt(otab$bayes, 3),
            fmt(otab$ari, 3)),
    "\\bottomrule",
    "\\end{tabular}"), "tab_overlap.tex")
  cat("\n")
}


## ===========================================================================
## PART B. THE RESULTS
## ===========================================================================
if ("results" %in% AN$parts) {

  cat("RESULTS\n")
  f_sum <- file.path(RES, sprintf("%s.csv", AN$tag))
  f_par <- file.path(RES, sprintf("%s-params.csv", AN$tag))
  f_sta <- file.path(RES, sprintf("%s-stations.csv", AN$tag))
  f_cov <- file.path(RES, sprintf("%s-coverage.csv", AN$tag))
  f_stb <- file.path(RES, sprintf("%s-stability.csv", AN$tag))
  ## --tag may name several sets of results -- one per process when the study
  ## was run in parallel -- and they are read and stacked
  rd <- function(f) {
    f <- f[file.exists(f)]
    if (!length(f)) return(NULL)
    do.call(rbind, lapply(f, utils::read.csv, stringsAsFactors = FALSE))
  }
  S  <- rd(f_sum); P <- rd(f_par); ST <- rd(f_sta); CV <- rd(f_cov); SB <- rd(f_stb)
  have <- c(summary = !is.null(S), params = !is.null(P), stations = !is.null(ST),
            coverage = !is.null(CV), stability = !is.null(SB))
  for (nm in names(have))
    cat(sprintf("  %-10s %s\n", nm, if (have[[nm]]) "found" else "missing"))
  cat("\n")

  ## Block membership is a property of the DESIGN, not of the results: a cell
  ## is in a block if the design puts it there. The cells were deduplicated
  ## across blocks when they were run, so a row can belong to several.
  fkey <- function(d) paste(d$n, d$TN, d$K, sprintf("%.2f", d$omega), d$id,
                            substr(d$balance, 1, 3), d$knn)
  in_block <- function(d, b) fkey(d) %in% fkey(sim_cells(blocks = b))

  ## a Monte Carlo summary with its standard error
  mc <- function(x) {
    x <- x[is.finite(x)]
    c(mean = if (length(x)) mean(x) else NA_real_,
      se = if (length(x) > 1) stats::sd(x) / sqrt(length(x)) else NA_real_,
      M = length(x))
  }
  agg <- function(d, by, f) {
    sp <- split(d, d[by], drop = TRUE)
    do.call(rbind, lapply(sp, function(s) cbind(s[1, by, drop = FALSE], f(s))))
  }

  if (!is.null(S)) {
    S$block_core  <- in_block(S, "core")
    S$block_null  <- in_block(S, "null")
    S$block_scen  <- in_block(S, "scenarios")
    S$block_graph <- in_block(S, "graph")
    S$block_bal   <- in_block(S, "balance")

    ## -----------------------------------------------------------------------
    ## B1. The core: selection, recovery and prediction against n, omega, T
    ## -----------------------------------------------------------------------
    core <- S[S$block_core, ]
    if (nrow(core)) {
      tab <- agg(core, c("TN", "omega", "n"), function(s) {
        a <- mc(s$ari_true); k <- mc(s$k_correct); r <- mc(s$rmse_ratio)
        data.frame(M = a[["M"]], sel = k[["mean"]], sel_se = k[["se"]],
                   ari = a[["mean"]], ari_se = a[["se"]],
                   ari_sel = mc(s$ari_sel)[["mean"]],
                   ratio = stats::median(s$rmse_ratio, na.rm = TRUE),
                   secs = stats::median(s$secs, na.rm = TRUE))
      })
      tab <- tab[order(tab$TN, tab$omega, tab$n), ]
      cat("CORE: selection, recovery and prediction\n")
      print(tab, row.names = FALSE, digits = 3)

      tex_write(c(
        "\\begin{tabular}{rrrrrrrr}",
        "\\toprule",
        "$T$ & $\\omega$ & $n$ & $M$ & $\\Pr(\\hat k = K)$ & ARI at $K$ & ARI at $\\hat k$ & RMSE ratio \\\\",
        "\\midrule",
        sprintf("%d & %s & %d & %d & %s & %s & %s & %s \\\\", tab$TN,
                format(tab$omega), tab$n, as.integer(tab$M),
                paste0(fmt(tab$sel), " (", fmt(tab$sel_se), ")"),
                paste0(fmt(tab$ari), " (", fmt(tab$ari_se), ")"),
                fmt(tab$ari_sel), fmt(tab$ratio)),
        "\\bottomrule",
        "\\end{tabular}"), "tab_sim_core.tex")

      Ts <- sort(unique(tab$TN))
      pdf_open("fig_sim_core.pdf", 3.6 * length(Ts) + 1.6, 7.2)
      graphics::layout(matrix(c(seq_len(2 * length(Ts)), rep(2 * length(Ts) + 1, 2)),
                              nrow = 2), widths = c(rep(1, length(Ts)), 0.42))
      graphics::par(mar = c(4.0, 4.4, 2.4, 0.6), oma = c(0, 0, 2.2, 0))
      for (tt in Ts) for (v in c("ari", "sel")) {
        s <- tab[tab$TN == tt, ]
        graphics::plot(NA, xlim = range(N_GRID), ylim = c(0, 1), log = "x",
                       xlab = "number of locations",
                       ylab = if (v == "ari") "ARI at the true K" else "Pr(k selected = K)",
                       xaxt = "n", bty = "n", las = 1)
        graphics::axis(1, at = N_GRID, labels = N_GRID)
        for (wi in seq_along(OM_GRID)) {
          q <- s[abs(s$omega - OM_GRID[wi]) < 1e-8, ]
          q <- q[order(q$n), ]
          if (!nrow(q)) next
          se <- if (v == "ari") q$ari_se else q$sel_se
          graphics::polygon(c(q$n, rev(q$n)),
                            c(q[[v]] - 2 * se, rev(q[[v]] + 2 * se)),
                            col = grDevices::adjustcolor(PAL3[wi], 0.15), border = NA)
          graphics::lines(q$n, q[[v]], col = PAL3[wi], lwd = 1.8)
          graphics::points(q$n, q[[v]], col = PAL3[wi], pch = 19, cex = 0.7)
        }
        graphics::mtext(sprintf("T = %d", tt), side = 3, line = 0.6, cex = 0.86,
                        font = 2, col = INK)
      }
      graphics::par(mar = c(4, 0.2, 2.4, 0.2)); graphics::plot.new()
      graphics::legend("left", bty = "n", cex = 0.86,
                       legend = paste0("omega = ", OM_GRID, " (",
                                       sprintf("%.2f", dgp_separation(KTRUE, OM_GRID)),
                                       " sd)"),
                       col = PAL3, lwd = 1.8, pch = 19)
      graphics::mtext("Recovery and selection against the network size, the overlap and the series length",
                      outer = TRUE, side = 3, line = 0.6, cex = 1.05, font = 2, col = INK)
      pdf_close("fig_sim_core.pdf")
    }

    ## -----------------------------------------------------------------------
    ## B2. The null: how often a homogeneous network is split
    ## -----------------------------------------------------------------------
    ## A homogeneous network is split when the rule selects k > 1 instead of
    ## the pooled model.
    nul <- S[S$block_null, ]
    if (nrow(nul)) {
      tn <- agg(nul, c("TN", "n"), function(s) {
        k <- mc(as.integer(s$k_hat > 1))
        data.frame(M = k[["M"]], false_split = k[["mean"]], se = k[["se"]])
      })
      cat("\nNULL (K = 1): how often the procedure splits a homogeneous network\n")
      print(tn[order(tn$TN, tn$n), ], row.names = FALSE, digits = 3)
      tex_write(c(
        "\\begin{tabular}{rrrr}", "\\toprule",
        "$T$ & $n$ & $M$ & $\\Pr(\\hat k > 1)$ \\\\", "\\midrule",
        sprintf("%d & %d & %d & %s (%s) \\\\", tn$TN, tn$n, as.integer(tn$M),
                fmt(tn$false_split), fmt(tn$se)),
        "\\bottomrule", "\\end{tabular}"), "tab_sim_null.tex")
    }

    ## -----------------------------------------------------------------------
    ## B3. The scenarios: what has to differ for the difference to be found
    ## -----------------------------------------------------------------------
    sc <- S[S$block_scen, ]
    if (nrow(sc)) {
      ts <- agg(sc, c("id", "omega"), function(s) {
        a <- mc(s$ari_true); k <- mc(s$k_correct)
        data.frame(M = a[["M"]], ari = a[["mean"]], ari_se = a[["se"]],
                   sel = k[["mean"]],
                   ratio = stats::median(s$rmse_ratio, na.rm = TRUE))
      })
      ids <- SIM_LEVELS$scenario[SIM_LEVELS$scenario %in% ts$id]
      cat("\nSCENARIOS: recovery by what separates the regimes\n")
      print(ts[order(match(ts$id, ids), ts$omega), ], row.names = FALSE, digits = 3)

      pdf_open("fig_sim_scenarios.pdf", 10.0, 4.4)
      graphics::par(mfrow = c(1, 2), mar = c(4.6, 4.4, 2.4, 0.6), oma = c(0, 0, 1.8, 0))
      for (v in c("ari", "sel")) {
        m <- sapply(OM_GRID, function(w) sapply(ids, function(i) {
          q <- ts[ts$id == i & abs(ts$omega - w) < 1e-8, v]
          if (length(q)) q else NA_real_ }))
        m <- matrix(m, nrow = length(ids), dimnames = list(ids, format(OM_GRID)))
        graphics::barplot(t(m), beside = TRUE, col = PAL3, border = NA, las = 2,
                          ylim = c(0, 1), cex.names = 0.8,
                          ylab = if (v == "ari") "ARI at the true K" else "Pr(k selected = K)")
        graphics::mtext(if (v == "ari") "recovery of the partition" else "selection of K",
                        side = 3, line = 0.5, cex = 0.86, font = 2, col = INK)
      }
      graphics::legend("topright", bty = "n", cex = 0.78, fill = PAL3, border = NA,
                       legend = paste0("omega = ", OM_GRID))
      graphics::mtext("What has to differ between regimes for the difference to be found",
                      outer = TRUE, side = 3, line = 0.2, cex = 1.05, font = 2, col = INK)
      pdf_close("fig_sim_scenarios.pdf")
    }

    ## -----------------------------------------------------------------------
    ## B4. The graph: does the answer turn on knn?
    ## -----------------------------------------------------------------------
    gr <- S[S$block_graph, ]
    if (nrow(gr)) {
      tg <- agg(gr, c("knn", "omega", "n"), function(s) {
        a <- mc(s$ari_true)
        data.frame(M = a[["M"]], ari = a[["mean"]], ari_se = a[["se"]],
                   sel = mc(s$k_correct)[["mean"]])
      })
      cat("\nGRAPH: recovery against the neighbourhood size\n")
      print(tg[order(tg$omega, tg$n, tg$knn), ], row.names = FALSE, digits = 3)

      pdf_open("fig_sim_graph.pdf", 3.4 * length(OM_GRID) + 1.4, 4.0)
      graphics::layout(matrix(seq_len(length(OM_GRID) + 1), nrow = 1),
                       widths = c(rep(1, length(OM_GRID)), 0.42))
      graphics::par(mar = c(4.2, 4.4, 2.4, 0.6), oma = c(0, 0, 1.8, 0))
      for (w in OM_GRID) {
        graphics::plot(NA, xlim = range(N_GRID), ylim = c(0, 1), log = "x",
                       xlab = "number of locations", ylab = "ARI at the true K",
                       xaxt = "n", bty = "n", las = 1)
        graphics::axis(1, at = N_GRID, labels = N_GRID)
        for (ki in seq_along(KNN)) {
          q <- tg[tg$knn == KNN[ki] & abs(tg$omega - w) < 1e-8, ]
          q <- q[order(q$n), ]
          if (!nrow(q)) next
          graphics::lines(q$n, q$ari, col = PAL3[ki], lwd = 1.8)
          graphics::points(q$n, q$ari, col = PAL3[ki], pch = 19, cex = 0.7)
        }
        graphics::mtext(om_lab(w), side = 3, line = 0.6, cex = 0.86, font = 2, col = INK)
      }
      graphics::par(mar = c(4.2, 0.2, 2.4, 0.2)); graphics::plot.new()
      graphics::legend("left", bty = "n", cex = 0.86, legend = paste0("knn = ", KNN),
                       col = PAL3, lwd = 1.8, pch = 19)
      graphics::mtext("Does the recovery turn on the neighbourhood graph?",
                      outer = TRUE, side = 3, line = 0.2, cex = 1.05, font = 2, col = INK)
      pdf_close("fig_sim_graph.pdf")
    }

    ## -----------------------------------------------------------------------
    ## B5. Balance
    ## -----------------------------------------------------------------------
    bl <- S[S$block_bal | (S$block_core & S$TN == SIM_REFERENCE$TN), ]
    if (nrow(bl) && length(unique(bl$balance)) > 1) {
      tb <- agg(bl, c("balance", "omega", "n"), function(s) {
        a <- mc(s$ari_true)
        data.frame(M = a[["M"]], ari = a[["mean"]], sel = mc(s$k_correct)[["mean"]])
      })
      cat("\nBALANCE: equal against unequal regime sizes\n")
      print(tb[order(tb$omega, tb$n, tb$balance), ], row.names = FALSE, digits = 3)
    }
  }

  ## -------------------------------------------------------------------------
  ## B6. The parameters: bias and RMSE, regime by regime
  ## -------------------------------------------------------------------------
  if (!is.null(P)) {
    P$err <- P$estimate - P$truth
    ref <- P[P$TN == SIM_REFERENCE$TN & P$K == SIM_REFERENCE$K &
             P$id == SIM_REFERENCE$scenario & P$knn == SIM_REFERENCE$knn &
             P$balance == SIM_REFERENCE$balance, ]
    if (nrow(ref)) {
      tp <- agg(ref, c("parameter", "omega", "n"), function(s) {
        e <- s$err[is.finite(s$err)]
        data.frame(M = length(e),
                   bias = if (length(e)) mean(e) else NA_real_,
                   rmse = if (length(e)) sqrt(mean(e^2)) else NA_real_,
                   rel = if (length(e)) sqrt(mean(e^2)) / mean(abs(s$truth)) else NA_real_)
      })
      tp <- tp[order(tp$parameter, tp$omega, tp$n), ]
      cat("\nPARAMETERS at the reference scenario, pooled over regimes\n")
      print(tp, row.names = FALSE, digits = 3)
      tex_write(c(
        "\\begin{tabular}{lrrrrrr}", "\\toprule",
        "parameter & $\\omega$ & $n$ & $M$ & bias & RMSE & relative RMSE \\\\",
        "\\midrule",
        sprintf("%s & %s & %d & %d & %s & %s & %s \\\\",
                gsub("_", "\\\\_", tp$parameter), format(tp$omega), tp$n,
                as.integer(tp$M), fmt(tp$bias, 3), fmt(tp$rmse, 3), fmt(tp$rel, 2)),
        "\\bottomrule", "\\end{tabular}"), "tab_sim_params.tex")
    } else {
      cat("\nPARAMETERS: no replication at the reference scenario, T and graph yet\n")
    }
  }

  ## -------------------------------------------------------------------------
  ## B7. The stations: where the clustering helps
  ##
  ## The station is the unit of the error measures. On a balanced panel the
  ## pooled RMSE is exactly the root mean of the per-station MSEs, so it is a
  ## summary of this distribution and not an alternative to it.
  ## -------------------------------------------------------------------------
  if (!is.null(ST)) {
    ST$ratio <- ST$rmse_true / ST$rmse_pool
    sref <- ST[ST$TN == SIM_REFERENCE$TN & ST$K == SIM_REFERENCE$K &
               ST$id == SIM_REFERENCE$scenario & ST$knn == SIM_REFERENCE$knn &
               ST$balance == SIM_REFERENCE$balance, ]
    if (nrow(sref)) {
      tsn <- agg(sref, c("omega", "n"), function(s) {
        r <- s$ratio[is.finite(s$ratio)]
        data.frame(stations = length(r),
                   improved = if (length(r)) mean(r < 1) else NA_real_,
                   q10 = if (length(r)) stats::quantile(r, 0.10) else NA_real_,
                   median = if (length(r)) stats::median(r) else NA_real_,
                   q90 = if (length(r)) stats::quantile(r, 0.90) else NA_real_,
                   misassigned_ratio = {
                     x <- s$ratio[s$correct == 0 & is.finite(s$ratio)]
                     if (length(x)) stats::median(x) else NA_real_ })
      })
      cat("\nSTATIONS: per-station RMSE, clustered over pooled\n")
      print(tsn[order(tsn$omega, tsn$n), ], row.names = FALSE, digits = 3)

      pdf_open("fig_sim_stations.pdf", 3.4 * length(OM_GRID), 3.8)
      graphics::par(mfrow = c(1, length(OM_GRID)), mar = c(4.2, 4.4, 2.4, 0.6),
                    oma = c(0, 0, 1.8, 0))
      for (w in OM_GRID) {
        q <- sref[abs(sref$omega - w) < 1e-8 & is.finite(sref$ratio), ]
        ns <- sort(unique(q$n))
        if (!length(ns)) { graphics::plot.new(); next }
        graphics::boxplot(ratio ~ n, data = q, log = "y", outline = FALSE,
                          col = grDevices::adjustcolor(PAL3[match(w, OM_GRID)], 0.35),
                          border = PAL3[match(w, OM_GRID)], las = 1,
                          xlab = "number of locations",
                          ylab = "station RMSE, clustered / pooled")
        graphics::abline(h = 1, col = GREY, lty = 2)
        graphics::mtext(om_lab(w), side = 3, line = 0.6, cex = 0.86, font = 2, col = INK)
      }
      graphics::mtext("Where the clustering helps: the per-station error, relative to the pooled model",
                      outer = TRUE, side = 3, line = 0.2, cex = 1.05, font = 2, col = INK)
      pdf_close("fig_sim_stations.pdf")
    } else {
      cat("\nSTATIONS: no replication at the reference scenario, T and graph yet\n")
    }
  }

  ## -------------------------------------------------------------------------
  ## B8. Coverage of the bootstrap
  ## -------------------------------------------------------------------------
  if (!is.null(CV)) {
    fam <- c(normal = "normal", basic = "basic", perc = "percentile", bc = "bias-corrected")
    CV$se_mc <- stats::ave(CV$estimate, CV$cell, CV$level, CV$regime, CV$parameter,
                           FUN = function(x) stats::sd(x, na.rm = TRUE))
    for (f in names(fam))
      CV[[paste0("in_", f)]] <- CV[[paste0(f, "_lo")]] <= CV$truth &
                                CV$truth <= CV[[paste0(f, "_up")]]
    tc <- agg(CV, c("level", "parameter"), function(s) {
      out <- data.frame(M = nrow(s))
      for (f in names(fam)) out[[f]] <- mean(s[[paste0("in_", f)]], na.rm = TRUE)
      out$se_ratio <- mean(s$se, na.rm = TRUE) / mean(s$se_mc, na.rm = TRUE)
      out
    })
    tc <- tc[order(tc$level, tc$parameter), ]
    cat("\nCOVERAGE of the bootstrap intervals, and bootstrap se over Monte Carlo sd\n")
    print(tc, row.names = FALSE, digits = 3)
    tex_write(c(
      "\\begin{tabular}{rlrrrrrr}", "\\toprule",
      "level & parameter & $M$ & normal & basic & percentile & bias-corrected & se ratio \\\\",
      "\\midrule",
      sprintf("%s & %s & %d & %s & %s & %s & %s & %s \\\\", fmt(tc$level),
              gsub("_", "\\\\_", tc$parameter), as.integer(tc$M), fmt(tc$normal, 3),
              fmt(tc$basic, 3), fmt(tc$perc, 3), fmt(tc$bc, 3), fmt(tc$se_ratio)),
      "\\bottomrule", "\\end{tabular}"), "tab_sim_coverage.tex")
  }
  if (!is.null(SB)) {
    cat("\nSTABILITY of the partition across bootstrap refits\n")
    print(agg(SB, c("omega", "id"), function(s)
      data.frame(M = nrow(s), ari_mean = mean(s$ari_mean, na.rm = TRUE),
                 share_identical = mean(s$ari_share1, na.rm = TRUE))),
      row.names = FALSE, digits = 3)
  }

  if (!any(have))
    cat("no results yet: run run-simulations.R first, or point --results= at them\n")
}

cat("\ndone\n")
