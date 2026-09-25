## ---------------------------------------------------------------------------
## The two factors that act on the geometry, drawn.
##
##   fig_design_omega.pdf   the point configuration as the overlap varies, at
##                          every network size of the design
##   fig_design_knn.pdf     the neighbourhood graph the Potts penalty lives on,
##                          as knn and the overlap vary
##   fig_design_graph.pdf   what the graph does in numbers: degree, edge length,
##                          and the share of edges that cross a regime boundary
##
## WHY THE THIRD FIGURE IS THE IMPORTANT ONE. knn does not change the data --
## the coordinates, the covariate and the response are identical whatever the
## graph -- it changes which pairs the penalty rewards for agreeing. The
## quantity that decides whether the penalty helps is therefore the share of
## edges that join two DIFFERENT true regimes: those are the pairs the penalty
## pushes towards a common label when they should not have one. At omega = 0
## the clusters coincide and almost every edge crosses; at omega = 1 few do. A
## larger knn reaches further and crosses more. That trade-off is the whole
## argument for treating knn as a factor rather than a setting.
##
##   Rscript dev/paper/15-design-figures.R
## ---------------------------------------------------------------------------

local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  h <- if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
  source(file.path(h, "00-setup.R"), chdir = TRUE)
}, envir = globalenv())
source(file.path(stem_paper_dir(), "design.R"), chdir = TRUE)
source(file.path(stem_paper_dir(), "06-dgp.R"), chdir = TRUE)
stem_load()
OUT <- stem_fig_dir()

## the levels of the design
N_GRID  <- stem_design_levels()$n
OM_GRID <- stem_design_levels()$omega
KNN     <- stem_design_levels()$knn
KTRUE   <- 3L

COL  <- c("#1f6f8b", "#e0a458", "#5b8c5a", "#a8516e")
GREY <- "#8a939f"
INK  <- "#1b2430"
EDGE <- "#9aa4b0"
CROSS <- "#c0392b"

om_lab <- function(w) sprintf("omega = %s", format(w))

## scstem_neighbors() is internal: reach it through the namespace so that the
## script runs both under pkgload::load_all() and against an installed package
neighbors <- utils::getFromNamespace("scstem_neighbors", "Stem")

draw <- function(n, K, omega, seed = 1) {
  loc <- dgp_locations(n, K, omega, balance = "balanced", seed = seed)
  list(xy = loc$xy, g = loc$labels, mu = loc$mu, coords = loc$coords)
}

panel <- function(s, lim, cex_pt = 0.9, centres = TRUE) {
  graphics::plot(NA, xlim = lim, ylim = lim, asp = 1, xlab = "", ylab = "",
                 xaxt = "n", yaxt = "n", bty = "n")
  graphics::rect(lim[1], lim[1], lim[2], lim[2], border = "#dfe4ea", lwd = 0.8)
  graphics::abline(h = 0, v = 0, col = "#eef1f4", lwd = 0.8)
  if (centres) graphics::points(s$mu[, 1], s$mu[, 2], pch = 3, col = "#5a626d",
                                cex = 1.3, lwd = 1.6)
  graphics::points(s$xy[, 1], s$xy[, 2], pch = 21, cex = cex_pt, lwd = 0.4,
                   bg = grDevices::adjustcolor(COL[s$g], 0.9), col = "white")
}

## a row label written into the outer margin, right-aligned against the panels
row_label <- function(top, bottom, main, sub) {
  y <- (top + bottom) / 2
  graphics::mtext(main, side = 2, outer = TRUE, at = y, line = 2.4,
                  las = 1, adj = 1, cex = 1.0, font = 2, col = INK)
  graphics::mtext(sub, side = 2, outer = TRUE, at = y - 0.035, line = 2.4,
                  las = 1, adj = 1, cex = 0.72, col = GREY)
}

## ---------------------------------------------------------------------------
## Figure 1: the overlap, across the network sizes of the design
## ---------------------------------------------------------------------------
lim <- c(-2.7, 2.7)
grDevices::cairo_pdf(file.path(OUT, "fig_design_omega.pdf"),
                     width = 11.0, height = 7.4)
graphics::par(mfrow = c(length(OM_GRID), length(N_GRID)),
              mar = c(0.4, 0.4, 0.4, 0.4), oma = c(3.2, 12.5, 4.2, 0.8))
nr <- length(OM_GRID)
for (wi in seq_along(OM_GRID)) {
  w <- OM_GRID[wi]
  for (n in N_GRID) {
    s <- draw(n, KTRUE, w, seed = 100 + n)
    panel(s, lim, cex_pt = if (n > 200) 0.55 else if (n > 100) 0.75 else
                           if (n > 50) 0.95 else 1.25)
    if (identical(w, OM_GRID[1]))
      graphics::mtext(paste0("n = ", n), side = 3, line = 0.8, cex = 0.95, font = 2)
    if (n == N_GRID[1])
      row_label(1 - (wi - 1) / nr, 1 - wi / nr, om_lab(w),
                sprintf("%.2f sd apart", dgp_separation(KTRUE, w)))
  }
}
graphics::mtext(sprintf("K = %d regimes: the overlap parameter across the network sizes of the design", KTRUE),
                outer = TRUE, side = 3, line = 2.4, cex = 1.15, font = 2, col = INK)
graphics::mtext(paste0("centres 2*omega apart, total spatial variance held fixed at ", round(NU_TOT, 3),
                       "; crosses are the true centres, colour the true regime"),
                outer = TRUE, side = 1, line = 1.2, cex = 0.78, col = GREY)
invisible(grDevices::dev.off())
message("wrote fig_design_omega.pdf")

## ---------------------------------------------------------------------------
## Figure 2: the graph the penalty lives on
## ---------------------------------------------------------------------------
graph_panel <- function(s, knn, lim, cex_pt = 0.8) {
  nb <- neighbors(s$coords, knn = knn, distance = "geo")
  W  <- nb$W
  graphics::plot(NA, xlim = lim, ylim = lim, asp = 1, xlab = "", ylab = "",
                 xaxt = "n", yaxt = "n", bty = "n")
  graphics::rect(lim[1], lim[1], lim[2], lim[2], border = "#dfe4ea", lwd = 0.8)
  ij <- which(upper.tri(W) & W > 0, arr.ind = TRUE)
  ## an edge that joins two different true regimes is the one the penalty
  ## rewards for agreeing when it should not
  cross <- s$g[ij[, 1]] != s$g[ij[, 2]]
  graphics::segments(s$xy[ij[!cross, 1], 1], s$xy[ij[!cross, 1], 2],
                     s$xy[ij[!cross, 2], 1], s$xy[ij[!cross, 2], 2],
                     col = EDGE, lwd = 0.5)
  graphics::segments(s$xy[ij[cross, 1], 1], s$xy[ij[cross, 1], 2],
                     s$xy[ij[cross, 2], 1], s$xy[ij[cross, 2], 2],
                     col = grDevices::adjustcolor(CROSS, 0.75), lwd = 0.8)
  graphics::points(s$xy[, 1], s$xy[, 2], pch = 21, cex = cex_pt, lwd = 0.4,
                   bg = grDevices::adjustcolor(COL[s$g], 0.9), col = "white")
  c(deg = mean(rowSums(W > 0)), cross = mean(cross), nedge = nrow(ij))
}

NG <- 100L
grDevices::cairo_pdf(file.path(OUT, "fig_design_knn.pdf"),
                     width = 9.8, height = 9.8)
graphics::par(mfrow = c(length(KNN), length(OM_GRID)),
              mar = c(1.4, 0.4, 0.4, 0.4), oma = c(4.0, 12.0, 4.2, 0.8))
nr <- length(KNN)
for (ki in seq_along(KNN)) {
  kk <- KNN[ki]
  for (w in OM_GRID) {
    s <- draw(NG, KTRUE, w, seed = 100 + NG)
    st <- graph_panel(s, kk, lim)
    if (ki == 1L)
      graphics::mtext(om_lab(w), side = 3, line = 0.8, cex = 1.0, font = 2)
    if (identical(w, OM_GRID[1]))
      row_label(1 - (ki - 1) / nr, 1 - ki / nr, paste0("knn = ", kk),
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
invisible(grDevices::dev.off())
message("wrote fig_design_knn.pdf")

## ---------------------------------------------------------------------------
## Figure 3: the graph in numbers, over the whole design
## ---------------------------------------------------------------------------
tab <- do.call(rbind, lapply(N_GRID, function(n)
  do.call(rbind, lapply(KNN, function(kk)
    do.call(rbind, lapply(OM_GRID, function(w) {
      s  <- draw(n, KTRUE, w, seed = 100 + n)
      nb <- neighbors(s$coords, knn = kk, distance = "geo")
      ij <- which(upper.tri(nb$W) & nb$W > 0, arr.ind = TRUE)
      dm <- as.matrix(geodist::geodist(s$coords, measure = "geodesic")) / 1000
      data.frame(n = n, knn = kk, omega = w,
                 degree = mean(rowSums(nb$W > 0)),
                 cross  = mean(s$g[ij[, 1]] != s$g[ij[, 2]]),
                 len    = stats::median(dm[ij]),
                 stringsAsFactors = FALSE)
    }))))))

grDevices::cairo_pdf(file.path(OUT, "fig_design_graph.pdf"),
                     width = 11.4, height = 4.2)
graphics::layout(matrix(c(1, 2, 3, 4), nrow = 1), widths = c(1, 1, 1, 0.48))
graphics::par(mar = c(4.2, 4.6, 2.8, 0.8), oma = c(1.0, 0, 2.4, 0))
pc <- c("#1f6f8b", "#e0a458", "#a8516e")
for (v in c("degree", "cross", "len")) {
  ylab <- switch(v, degree = "mean degree",
                 cross = "share of edges crossing a regime",
                 len = "median edge length (km)")
  ylim <- range(tab[[v]]); if (v == "cross") ylim <- c(0, max(ylim))
  graphics::plot(NA, xlim = range(N_GRID), ylim = ylim, log = "x",
                 xlab = "", ylab = "", xaxt = "n", bty = "n",
                 col.axis = "#3d4753", las = 1)
  graphics::axis(1, at = N_GRID, labels = N_GRID, col.axis = "#3d4753")
  graphics::mtext("number of locations", side = 1, line = 2.5, cex = 0.8,
                  col = "#3d4753")
  ## for K balanced regimes a graph carrying no information about the
  ## partition has a crossing share of 1 - 1/K
  if (v == "cross")
    graphics::abline(h = 1 - 1 / KTRUE, col = "#c9ced6", lty = "13", lwd = 1.2)
  for (ki in seq_along(KNN)) for (wi in seq_along(OM_GRID)) {
    s <- tab[tab$knn == KNN[ki] & tab$omega == OM_GRID[wi], ]
    s <- s[order(s$n), ]
    graphics::lines(s$n, s[[v]], col = pc[ki], lwd = 1.7,
                    lty = c("solid", "22", "44")[wi])
    graphics::points(s$n, s[[v]], col = pc[ki], pch = 19, cex = 0.7)
  }
  if (v == "cross")
    graphics::text(N_GRID[length(N_GRID)], 1 - 1 / KTRUE, adj = c(1, 1.6), cex = 0.72,
                   col = "#8a939f", labels = "no information: 1 - 1/K")
  graphics::mtext(ylab, side = 3, line = 0.8, cex = 0.88, font = 2, col = INK)
}
graphics::par(mar = c(4.2, 0.2, 2.8, 0.2))
graphics::plot.new()
graphics::legend("left", bty = "n", cex = 0.86,
                 legend = c(paste0("knn = ", KNN), "", paste0("omega = ", OM_GRID)),
                 col = c(pc, NA, rep("#3d4753", 3)),
                 lty = c(1, 1, 1, NA, 1, 2, 3), lwd = 1.7, seg.len = 1.8)
graphics::mtext(sprintf("What the neighbourhood graph does, K = %d", KTRUE),
                outer = TRUE, side = 3, line = 0.9, cex = 1.1, font = 2, col = INK)
graphics::mtext(paste("knn multiplies the number of edges, and so the weight the penalty carries at a given phi,",
                      "but from n = 100 on it leaves the information the graph holds about the partition unchanged"),
                outer = TRUE, side = 1, line = -0.2, cex = 0.76, col = GREY)
invisible(grDevices::dev.off())
message("wrote fig_design_graph.pdf")

cat("\nthe graph in numbers\n")
print(tab[order(tab$knn, tab$omega, tab$n), ], row.names = FALSE, digits = 3)
