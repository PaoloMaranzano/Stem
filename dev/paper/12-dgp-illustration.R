## ---------------------------------------------------------------------------
## What a simulated data set looks like: the geometry and the series it carries.
##
## The overlap figures of 09-overlap-figures.R show only where the stations are.
## This one puts the two halves of the design side by side, because the point of
## the experiment is that the two need not agree: at omega = 0 the regimes are
## spatially indistinguishable and can only be told apart from the series, and a
## scenario that separates them weakly in the parameters produces series that
## look alike even when the map is perfectly separated.
##
##   Rscript dev/paper/12-dgp-illustration.R
## ---------------------------------------------------------------------------

local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  h <- if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
  source(file.path(h, "00-setup.R"), chdir = TRUE)
  source(file.path(h, "06-dgp.R"), chdir = TRUE)
}, envir = globalenv())
OUT <- stem_fig_dir()

COL  <- c("#1f6f8b", "#e0a458", "#5b8c5a", "#a8516e")
GREY <- "#8a939f"
INK  <- "#1b2430"

scen <- dgp_scenarios()

## the rows of the figure: what varies is the number of regimes, the overlap,
## the balance and how far apart the parameters are
## Rows 2 to 4 walk the overlap up at a fixed scenario; row 5 holds the overlap
## at its largest and weakens the parameters instead, which is the case the
## design exists to separate; row 6 adds a third regime and an unbalanced
## allocation. omega is the d of the source paper -- see 06-dgp.R.
rows <- list(
  list(K = 1L, omega = 2/3, id = "S0",  bal = "balanced",
       lab = "K = 1 (pooled)",     sub = "no regime at all"),
  list(K = 2L, omega = 0,   id = "S4",  bal = "balanced",
       lab = "K = 2, omega = 0",   sub = "all three contrasts, no spatial separation"),
  list(K = 2L, omega = 2/3, id = "S4",  bal = "balanced",
       lab = "K = 2, omega = 2/3", sub = "all three contrasts, 2.11 sd apart"),
  list(K = 2L, omega = 1,   id = "S4",  bal = "balanced",
       lab = "K = 2, omega = 1",   sub = "all three contrasts, 3.16 sd apart"),
  list(K = 2L, omega = 1,   id = "S1a", bal = "balanced",
       lab = "K = 2, omega = 1",   sub = "separated in space, coefficients 0.25 sd apart"),
  list(K = 3L, omega = 1,   id = "S4",  bal = "unbalanced",
       lab = "K = 3, omega = 1",   sub = "unbalanced regimes")
)

N  <- 60L
TN <- 365L
NS <- 3L        # stations drawn per regime in the series panel

## ---------------------------------------------------------------------------
map_panel <- function(dat, lim) {
  graphics::plot(NA, xlim = lim, ylim = lim, asp = 1, xlab = "", ylab = "",
                 xaxt = "n", yaxt = "n", bty = "n")
  graphics::rect(lim[1], lim[1], lim[2], lim[2], border = "#dfe4ea", lwd = 0.8)
  graphics::abline(h = 0, v = 0, col = "#eef1f4", lwd = 0.8)
  graphics::points(dat$mu[, 1], dat$mu[, 2], pch = 3, col = GREY, cex = 1.0,
                   lwd = 1.2)
  graphics::points(dat$xy[, 1], dat$xy[, 2], pch = 21, cex = 0.85, lwd = 0.5,
                   bg = grDevices::adjustcolor(COL[dat$g], 0.85), col = "white")
}

## The measurement error has twice the standard deviation of the latent process,
## which is what a real network looks like and what makes the problem hard -- but
## it also means that a panel of raw series is a band of noise. The panel
## therefore shows a WINDOW of the record rather than all of it, the individual
## stations very light, and the signal each regime carries in full strength.
WIN <- 120L

series_panel <- function(z, y, g, pick, ylim) {
  tt <- seq_len(WIN)
  graphics::plot(NA, xlim = range(tt), ylim = ylim, xlab = "", ylab = "",
                 xaxt = "n", yaxt = "n", bty = "n")
  graphics::abline(h = 0, col = "#eef1f4", lwd = 0.8)
  for (i in pick) {
    graphics::lines(tt, z[tt, i], col = grDevices::adjustcolor(COL[g[i]], 0.22),
                    lwd = 0.4)
  }
  ## the signal of each regime: its latent process on the scale of the response
  for (k in seq_len(ncol(y))) {
    idx <- which(g == k)
    graphics::lines(tt, y[tt, k] + mean(z[, idx]), col = COL[k], lwd = 2.0)
  }
  graphics::axis(1, at = pretty(tt, 4), cex.axis = 0.7, col = GREY,
                 col.axis = GREY, tck = -0.04, mgp = c(2, 0.35, 0))
  graphics::axis(2, at = pretty(ylim, 3), cex.axis = 0.7, col = GREY,
                 col.axis = GREY, las = 1, tck = -0.03, mgp = c(2, 0.45, 0))
}

## ---------------------------------------------------------------------------
draws <- lapply(rows, function(r) {
  row <- scen[scen$id == r$id, , drop = FALSE]
  d   <- dgp_draw(N, TN, r$K, r$omega, row, rep = 1L, balance = r$bal)
  loc <- dgp_locations(N, r$K, r$omega, balance = r$bal, seed = 1000L + 1L)
  set.seed(11)
  pick <- unlist(lapply(seq_len(r$K), function(k)
    sample(which(d$labels == k), min(NS, sum(d$labels == k)))))
  list(z = d$z, y = d$latent, g = d$labels, xy = loc$xy, mu = loc$mu,
       pick = pick, sizes = d$sizes, lab = r$lab, sub = r$sub)
})

lim  <- max(abs(unlist(lapply(draws, function(d) range(d$xy))))) * 1.05
lim  <- c(-lim, lim)
ylim <- range(unlist(lapply(draws, function(d) d$z[seq_len(WIN), d$pick])))

grDevices::cairo_pdf(file.path(OUT, "fig_dgp_examples.pdf"),
                     width = 9.2, height = 12.6)
graphics::par(mfcol = c(length(rows), 2), mar = c(1.6, 1.0, 1.4, 0.6),
              oma = c(3.2, 6.4, 4.0, 1.0))

for (d in draws) map_panel(d, lim)
for (d in draws) series_panel(d$z, d$y, d$g, d$pick, ylim)

## row labels down the left margin, one per row of the left column
h <- 1 / length(rows)
for (i in seq_along(draws)) {
  d <- draws[[i]]
  yy <- 1 - (i - 0.5) * h
  graphics::mtext(d$lab, side = 2, outer = TRUE, at = yy, line = 4.4,
                  las = 1, cex = 0.92, font = 2, col = INK, adj = 0)
  graphics::mtext(d$sub, side = 2, outer = TRUE, at = yy - 0.030, line = 4.4,
                  las = 1, cex = 0.72, col = GREY, adj = 0)
  graphics::mtext(paste0("sizes ", paste(d$sizes, collapse = "/")),
                  side = 2, outer = TRUE, at = yy - 0.055, line = 4.4,
                  las = 1, cex = 0.68, col = GREY, adj = 0)
}

graphics::mtext("locations", side = 3, outer = TRUE, at = 0.30, line = 1.0,
                cex = 1.0, font = 2, col = INK)
graphics::mtext("the series they carry", side = 3, outer = TRUE, at = 0.76,
                line = 1.0, cex = 1.0, font = 2, col = INK)
graphics::mtext(sprintf("One replication of the design at n = %d, T = %d", N, TN),
                outer = TRUE, side = 3, line = 2.5, cex = 1.15, font = 2,
                col = INK)
graphics::mtext("crosses mark the true cluster centres; thin lines are three stations per regime, over the first 120 periods",
                outer = TRUE, side = 1, line = 0.4, cex = 0.74, col = GREY)
graphics::mtext("thick lines are the latent process of each regime, shifted to the regime mean",
                outer = TRUE, side = 1, line = 1.4, cex = 0.74, col = GREY)
invisible(grDevices::dev.off())
message("wrote ", file.path(OUT, "fig_dgp_examples.pdf"))
