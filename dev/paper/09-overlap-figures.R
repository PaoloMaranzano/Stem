## Illustration of the overlapping design of Morelli, Maranzano and Otto (2026),
## Spatial Statistics 73, 100960, Section 4, adapted to K = 2 and K = 3 and to
## the sample sizes we can afford to fit.
##
## Their design places the K = 4 cluster centres at the corners of a square of
## half-side d, mu_sp = ((d,d), (-d,d), (d,-d), (-d,-d)), with Sigma_sp =
## nu_sp * I_2, so that adjacent centres are 2d apart. Reducing to K = 2 and
## K = 3 we keep the NEAREST-NEIGHBOUR centre distance at 2d rather than the
## radius, so that a given d means the same degree of overlap whatever K:
##
##   K = 2   centres (-d, 0) and (d, 0)
##   K = 3   equilateral triangle of side 2d, i.e. radius 2d/sqrt(3)
##   K = 4   the square of the paper, radius d*sqrt(2)
##
## The standardised separation is then 2*omega / sqrt(nu_sp), the same for every K:
## with nu_sp = 0.4 it runs from 0 at omega = 0 to 3.16 sd at omega = 1. The
## overlap parameter, their d, is called omega here; see 06-dgp.R.

local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  h <- if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
  source(file.path(h, "00-setup.R"), chdir = TRUE)
}, envir = globalenv())
source(file.path(stem_paper_dir(), "06-dgp.R"), chdir = TRUE)
OUT <- stem_fig_dir()

D_GRID <- c(0, 1/3, 2/3, 1)
N_GRID <- c(20, 40, 60, 80, 100)
COL    <- c("#1f6f8b", "#e0a458", "#5b8c5a", "#a8516e")
GREY   <- "#8a939f"

## One draw, in the abstract plane of the design: dgp_locations() also maps it
## to longitude and latitude, which the figures do not need.
##
## These figures illustrate the design AS PUBLISHED, in which the within-cluster
## variance nu_sp is what is held fixed. Our own design holds the TOTAL variance
## fixed instead -- see NU_TOT in 06-dgp.R and fig_design_omega.pdf -- because
## fixing the within-cluster variance lets the network grow with the separation,
## and the third figure below is precisely the one that shows why that matters.
## The published parameterisation is recovered from the new one by giving each
## cell the total its centres and its nu_sp imply.
draw <- function(n, K, d, nu_sp = NU_SP, seed = 1) {
  loc <- dgp_locations(n, K, d, nu_tot = nu_sp + dgp_centre_var(K, d),
                       balance = "balanced", seed = seed)
  list(xy = loc$xy, g = loc$labels, mu = loc$mu)
}

panel <- function(s, lim, cex_pt = 0.9, show_centres = TRUE) {
  graphics::plot(NA, xlim = lim, ylim = lim, asp = 1, xlab = "", ylab = "",
                 xaxt = "n", yaxt = "n", bty = "n")
  graphics::rect(lim[1], lim[1], lim[2], lim[2], border = "#dfe4ea", lwd = 0.8)
  graphics::abline(h = 0, v = 0, col = "#eef1f4", lwd = 0.8)
  if (show_centres) {
    graphics::points(s$mu[, 1], s$mu[, 2], pch = 3, col = GREY, cex = 1.1, lwd = 1.2)
  }
  graphics::points(s$xy[, 1], s$xy[, 2], pch = 21, cex = cex_pt, lwd = 0.6,
                   bg = grDevices::adjustcolor(COL[s$g], 0.85), col = "white")
}

## ---------------------------------------------------------------------------
## Figures 1 and 2: the d by n grid, one figure per K
## ---------------------------------------------------------------------------
overlap_grid <- function(K, file) {
  lim <- c(-3.2, 3.2)
  grDevices::cairo_pdf(file, width = 10.5, height = 9.2)
  graphics::par(mfrow = c(length(D_GRID), length(N_GRID)),
                mar = c(0.3, 0.3, 0.3, 0.3), oma = c(3.0, 7.0, 4.2, 0.8))
  for (i in seq_along(D_GRID)) {
    for (j in seq_along(N_GRID)) {
      s <- draw(N_GRID[j], K, D_GRID[i], seed = 100 * i + j)
      panel(s, lim, cex_pt = if (N_GRID[j] > 60) 0.8 else 1.0)
      if (i == 1) graphics::mtext(paste0("n = ", N_GRID[j]), side = 3, line = 0.8,
                                  cex = 0.85, font = 2)
      if (j == 1) {
        lab <- c("0", "1/3", "2/3", "1")[i]
        sdl <- sprintf("%.2f sd apart", 2 * D_GRID[i] / sqrt(NU_SP))
        ## one stacked expression rather than two calls, which cannot collide
        graphics::mtext(bquote(atop(bold(.(paste0("ω = ", lab))),
                                    scriptstyle(.(sdl)))),
                        side = 2, line = 2.4, cex = 0.95, las = 1)
      }
    }
  }
  graphics::mtext(sprintf("K = %d clusters: spatial overlap and sample size", K),
                  outer = TRUE, side = 3, line = 2.5, cex = 1.15, font = 2)
  graphics::mtext(paste0("cluster centres 2ω apart, isotropic dispersion ν_sp = ",
                         NU_SP, "; crosses mark the true centres, colour the true cluster"),
                  outer = TRUE, side = 1, line = 1.2, cex = 0.78, col = GREY)
  invisible(grDevices::dev.off())
  message("wrote ", file)
}

overlap_grid(2L, file.path(OUT, "fig_overlap_K2.pdf"))
overlap_grid(3L, file.path(OUT, "fig_overlap_K3.pdf"))

## ---------------------------------------------------------------------------
## Figure 3: the role of the dispersion, at fixed n
## ---------------------------------------------------------------------------
NU_GRID <- c(0.1, 0.4, 1.0)
grDevices::cairo_pdf(file.path(OUT, "fig_overlap_dispersion.pdf"),
               width = 9.0, height = 6.2)
graphics::par(mfrow = c(length(NU_GRID), length(D_GRID)),
              mar = c(0.4, 0.4, 0.4, 0.4), oma = c(3.4, 5.6, 4.2, 1.0))
for (i in seq_along(NU_GRID)) {
  for (j in seq_along(D_GRID)) {
    s <- draw(60L, 2L, D_GRID[j], nu_sp = NU_GRID[i], seed = 700 + 10 * i + j)
    panel(s, c(-3.6, 3.6))
    if (i == 1) {
      lab <- c("0", "1/3", "2/3", "1")[j]
      graphics::mtext(paste0("ω = ", lab), side = 3, line = 0.8,
                      cex = 0.9, font = 2)
    }
    if (j == 1) {
      sdl <- sprintf("sd = %.2f", sqrt(NU_GRID[i]))
      graphics::mtext(bquote(atop(bold(.(paste0("ν_sp = ", NU_GRID[i]))),
                                  scriptstyle(.(sdl)))),
                      side = 2, line = 2.2, cex = 0.95, las = 1)
    }
  }
}
graphics::mtext("K = 2, n = 60: dispersion against overlap",
                outer = TRUE, side = 3, line = 2.3, cex = 1.15, font = 2)
graphics::mtext("only the ratio 2ω / sqrt(ν_sp) governs the separation; the absolute scale governs how the domain compares with the spatial range of the model",
                outer = TRUE, side = 1, line = 1.4, cex = 0.72, col = GREY)
invisible(grDevices::dev.off())
message("wrote fig_overlap_dispersion.pdf")

## ---------------------------------------------------------------------------
## What the design implies, in numbers
## ---------------------------------------------------------------------------
cat("\nStandardised separation 2*omega/sqrt(nu_sp), by omega:\n")
print(round(setNames(2 * D_GRID / sqrt(NU_SP), c("0", "1/3", "2/3", "1")), 2))

cat("\nUnits per cluster, balanced allocation:\n")
tab <- outer(N_GRID, c(2, 3), function(n, K) floor(n / K))
dimnames(tab) <- list(paste0("n=", N_GRID), c("K=2", "K=3"))
print(tab)
cat("\nMinimum admissible regime size in the package is r+2; with r = 2\n",
    "covariates that is 4, with r = 6 it is 8.\n", sep = "")
