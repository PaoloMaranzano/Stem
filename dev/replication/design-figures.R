## ===========================================================================
## SC-STEM: the simulation study. THE FIGURES AND TABLES OF THE DESIGN.
##
## One script, self-contained. It needs its companion run-simulations.R in the
## same folder, from which it takes the generator, the parameter values and the
## blocks, so that the description of the design cannot disagree with the
## design that is run. It fits nothing and takes a few seconds.
##
##     Rscript design-figures.R
##     Rscript design-figures.R --out=D:/paper/Figures
##
##   fig_sim_geometry.pdf     the four geometries of the regimes in space
##   fig_sim_data.pdf         one data set of S1-strong-ind and one of
##                            S1-weak-shr, in time and in space
##   fig_sim_correlation.pdf  the correlation of the measurement error against
##                            distance, and of the latent process against lag
##   fig_sim_variance.pdf     the variance of the response of every regime of
##                            every scenario, by source
##   tab_sim_scenarios.tex    the parameter values of every regime
##   tab_sim_geometry.tex     the four geometries in numbers
##   tab_sim_separation.tex   the separation of the regimes seen by the
##                            assignment step
##   tab_sim_blocks.tex       the blocks of the design and their cells
##
## Outputs : <here>/output/  (--out=)
##
## This is part of the REPLICATION MATERIAL of the paper, not of the Stem
## package.
## ===========================================================================


## ---------------------------------------------------------------------------
## Where this script is, and the design from the runner beside it
## ---------------------------------------------------------------------------
fig_here <- local({
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
runner <- file.path(fig_here, "run-simulations.R")
if (!file.exists(runner))
  stop("run-simulations.R must sit in the same folder as this script:\n  ",
       fig_here, call. = FALSE)
SIM_DEFINE_ONLY <- TRUE
source(runner, local = globalenv())

FG <- sim_config(list(out = file.path(fig_here, "output")))
OUT <- FG$out[1]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
OUT <- normalizePath(OUT, winslash = "/")
cat("design figures and tables to ", OUT, "\n\n", sep = "")

## ---------------------------------------------------------------------------
## Style and small helpers
## ---------------------------------------------------------------------------
RC   <- c("#1f6f8b", "#e0a458", "#a8516e")       # regimes 1, 2, 3
GREY <- "#8a939f"
INK  <- "#1b2430"
pdf_open <- function(name, width, height)
  grDevices::cairo_pdf(file.path(OUT, name), width = width, height = height)
pdf_close <- function(name) {
  invisible(grDevices::dev.off())
  message("wrote ", name)
}
tex_write <- function(lines, name) {
  writeLines(lines, file.path(OUT, name))
  message("wrote ", name)
}
fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)
par_base <- function(...) graphics::par(mgp = c(2.1, 0.6, 0), las = 1, tcl = -0.3, ...)

## The four geometries of the core, from the most to the least overlapping
GEO <- data.frame(code   = c("0", "0.7", "1b", "1"),
                  omega  = c(0, 0.7, 1, 1),
                  spread = c("total", "total", "regime", "total"),
                  label  = c("total overlap", "medium-high overlap",
                             "medium-low overlap", "strong separation"),
                  stringsAsFactors = FALSE)
GEO$nu <- NU_TOT - mapply(function(om, sp)
  dgp_centre_var(3L, if (sp == "regime") SIM_REFERENCE$omega else om), GEO$omega, GEO$spread)
## share of the locations of a regime nearer the centre of another regime, by
## Monte Carlo on a fixed seed (a draw from N(mu_g, nu I) falls outside the
## 120-degree wedge of its own centre)
GEO$share_other <- vapply(seq_len(nrow(GEO)), function(i) {
  if (GEO$omega[i] == 0) return(2 / 3)
  mu <- dgp_centres(3L, GEO$omega[i])
  set.seed(20261002)
  M <- 4e5
  xy <- cbind(mu[1, 1] + stats::rnorm(M, sd = sqrt(GEO$nu[i])),
              mu[1, 2] + stats::rnorm(M, sd = sqrt(GEO$nu[i])))
  d2 <- sapply(1:3, function(g) (xy[, 1] - mu[g, 1])^2 + (xy[, 2] - mu[g, 2])^2)
  mean(max.col(-d2, ties.method = "first") != 1L)
}, numeric(1))
GEO$network_sd <- sqrt(GEO$nu + mapply(function(om) dgp_centre_var(3L, om), GEO$omega))

## The scenarios with three regimes, in the order of the design, and their
## parameters (the variant does not change them)
SC <- data.frame(id   = c("S1w-ind", "S1s-ind", "S3beta-ind", "S3G-ind", "S3Seta-ind",
                          "S3theta-ind", "S3error-ind"),
                 name = c("S1-weak", "S1-strong", "S3-beta", "S3-G", "S3-Seta",
                          "S3-theta", "S3-error"),
                 stringsAsFactors = FALSE)
PSI <- lapply(SC$id, function(i) dgp_psi(SIM_SCEN[SIM_SCEN$id == i, , drop = FALSE]))
names(PSI) <- SC$name
PSI[["S2"]] <- dgp_psi(SIM_SCEN[SIM_SCEN$id == "S2", , drop = FALSE])

## The reference cell of a scenario-variant
ref_cell <- function(scen) {
  data.frame(scenario = scen, n = SIM_REFERENCE$n, TN = SIM_REFERENCE$TN,
             omega = if (scen == "S2") NA_real_ else SIM_REFERENCE$omega,
             spread = SIM_REFERENCE$spread, balance = SIM_REFERENCE$balance,
             knn = SIM_REFERENCE$knn, stringsAsFactors = FALSE)
}


## ===========================================================================
## Figure: the four geometries
## ===========================================================================
pdf_open("fig_sim_geometry.pdf", 10, 10)
par_base(mfrow = c(2, 2), mar = c(3, 3, 3.8, 1))
th <- seq(0, 2 * pi, length.out = 120)
for (i in seq_len(nrow(GEO))) {
  loc <- dgp_locations(SIM_REFERENCE$n, 3L, GEO$omega[i], seed = 1001L, spread = GEO$spread[i])
  mu <- dgp_centres(3L, GEO$omega[i])
  graphics::plot(loc$coords, col = RC[loc$labels], pch = 16, cex = 0.9, asp = 1,
                 xlim = c(-3.3, 3.3), ylim = c(-3.3, 3.3), xlab = "", ylab = "", bty = "n")
  for (g in 1:3) {
    graphics::lines(mu[g, 1] + 2 * sqrt(GEO$nu[i]) * cos(th),
                    mu[g, 2] + 2 * sqrt(GEO$nu[i]) * sin(th), col = RC[g], lty = 2)
    graphics::points(mu[g, 1], mu[g, 2], pch = 3, col = INK, cex = 1.4, lwd = 2)
  }
  graphics::title(sprintf("(%s) %s: omega = %s, spread fixed %s\ncentres %.1f apart, regime sd %.2f, %.0f%% of the locations nearer another centre",
                          letters[i], GEO$label[i], format(GEO$omega[i]),
                          if (GEO$spread[i] == "total") "over the network" else "within a regime",
                          2 * GEO$omega[i], sqrt(GEO$nu[i]), 100 * GEO$share_other[i]),
                  cex.main = 0.9, font.main = 1)
}
pdf_close("fig_sim_geometry.pdf")


## ===========================================================================
## Figure: one data set in time and in space, at the two levels of separation
## ===========================================================================
draws <- list(`S1-strong-ind` = dgp_draw(ref_cell("S1s-ind"), 1L),
              `S1-weak-shr`   = dgp_draw(ref_cell("S1w-shr"), 1L))
TT <- 60L                                          # the period shown in space
map_col <- function(v, lim) {
  cl <- grDevices::hcl.colors(21, "Blue-Red 2")
  cl[findInterval(pmax(pmin(v, lim - 1e-9), -lim + 1e-9), seq(-lim, lim, length.out = 22))]
}
## one colour scale per kind of map, common to the two data sets
## (the 90th percentile of the absolute values, beyond which the colour saturates)
lim_mean <- stats::quantile(unlist(lapply(draws, function(dd) { m <- colMeans(dd$z); abs(m - mean(m)) })), 0.9)
lim_err  <- stats::quantile(unlist(lapply(draws, function(dd) abs(dd$z[TT, ] - dd$mu[TT, ]))), 0.9)
pdf_open("fig_sim_data.pdf", 12, 11)
## by column a data set: (a) and (b) span the column, (c) and (d) share a row
graphics::layout(matrix(c(1, 1, 5, 5,
                          2, 2, 6, 6,
                          3, 4, 7, 8), 3, byrow = TRUE), heights = c(1, 1, 1.1))
par_base(mar = c(3.4, 3.6, 3, 1))
for (nm in names(draws)) {
  dd <- draws[[nm]]; psi <- dd$psi; TN <- nrow(dd$z)
  ## (a) the latent paths
  graphics::matplot(seq_len(TN), dd$latent, type = "l", lty = 1, lwd = 1.3, col = RC,
                    xlab = "t", ylab = "latent process", bty = "n")
  graphics::legend("topleft", bty = "n", lwd = 2, col = RC, cex = 0.8, ncol = 3,
                   legend = sprintf("G = %.1f, v = %.2f", psi$G, psi$v))
  graphics::title(sprintf("%s: (a) the latent process of each regime", nm), cex.main = 0.95, font.main = 1)
  ## (b) the response at three locations of every regime
  graphics::plot(NA, xlim = c(1, TN), ylim = range(dd$z), xlab = "t", ylab = "response", bty = "n")
  for (g in 1:3) for (i in which(dd$labels == g)[1:3])
    graphics::lines(seq_len(TN), dd$z[, i], col = grDevices::adjustcolor(RC[g], 0.8), lwd = 0.9)
  graphics::title("(b) the response at three locations of every regime", cex.main = 0.95, font.main = 1)
  ## (c) the mean of the response over time, by location
  zm <- colMeans(dd$z)
  graphics::par(mar = c(2, 1.5, 3.4, 0.5))
  graphics::plot(dd$coordinates, pch = c(21, 22, 24)[dd$labels], bg = map_col(zm - mean(zm), lim_mean),
                 col = INK, cex = 1.25, asp = 1, xlab = "", ylab = "", axes = FALSE)
  graphics::title("(c) the mean of the response over time,\nas a deviation from the network mean", cex.main = 0.9, font.main = 1)
  ## (d) the measurement error at one period
  e <- dd$z[TT, ] - dd$mu[TT, ]
  graphics::plot(dd$coordinates, pch = c(21, 22, 24)[dd$labels], bg = map_col(e, lim_err),
                 col = INK, cex = 1.25, asp = 1, xlab = "", ylab = "", axes = FALSE)
  graphics::title(sprintf("(d) the measurement error at t = %d\n", TT), cex.main = 0.9, font.main = 1)
  graphics::par(mar = c(3.4, 3.6, 3, 1))
}
pdf_close("fig_sim_data.pdf")


## ===========================================================================
## Figure: the correlations and their decay
## ===========================================================================
## nearest-neighbour spacing of a network of n locations at the reference
nn_spacing <- vapply(SIM_BLOCKS$core$n, function(n) {
  xy <- dgp_locations(n, 3L, SIM_REFERENCE$omega, seed = 1001L)$coords
  D <- as.matrix(stats::dist(xy)); diag(D) <- Inf
  stats::median(apply(D, 1, min))
}, numeric(1))
pdf_open("fig_sim_correlation.pdf", 11, 8.5)
par_base(mfrow = c(2, 2), mar = c(3.6, 3.8, 3, 1))
hh <- seq(0, 4, length.out = 300)
for (lv in c("S1-weak", "S1-strong")) {
  psi <- PSI[[lv]]
  graphics::plot(NA, xlim = c(0, 4), ylim = c(0, 0.55), xlab = "distance between two locations",
                 ylab = "correlation of the measurement errors", bty = "n")
  for (g in 1:3)
    graphics::lines(hh, psi$so[g] / (psi$se[g] + psi$so[g]) * exp(-psi$theta[g] * hh), col = RC[g], lwd = 2)
  graphics::abline(v = nn_spacing, col = GREY, lty = 3)
  graphics::text(max(nn_spacing), 0.02, sprintf("nearest-neighbour spacing\nat n = %s",
                                                 paste(rev(SIM_BLOCKS$core$n), collapse = ", ")),
                 pos = 4, cex = 0.7, col = GREY)
  graphics::legend("topright", bty = "n", lwd = 2, col = RC, cex = 0.8,
                   legend = sprintf("regime %d: sill share %.2f, range %.1f", 1:3,
                                    psi$so / (psi$se + psi$so), psi$R))
  graphics::title(sprintf("%s: the error, against distance", lv), cex.main = 0.95, font.main = 1)
}
lag <- 0:12
for (lv in c("S1-weak", "S1-strong")) {
  psi <- PSI[[lv]]
  graphics::plot(NA, xlim = range(lag), ylim = c(0, 1), xlab = "lag (periods)",
                 ylab = "autocorrelation of the latent process", bty = "n")
  for (g in 1:3) graphics::lines(lag, psi$G[g]^lag, col = RC[g], lwd = 2, type = "o", pch = 16, cex = 0.6)
  graphics::abline(h = 0.5, col = GREY, lty = 3)
  graphics::legend("topright", bty = "n", lwd = 2, col = RC, cex = 0.8,
                   legend = sprintf("regime %d: G = %.1f, half-life %.1f", 1:3, psi$G, log(0.5) / log(psi$G)))
  graphics::title(sprintf("%s: the latent process, against lag", lv), cex.main = 0.95, font.main = 1)
}
pdf_close("fig_sim_correlation.pdf")


## ===========================================================================
## Figure: the variance of the response by source, regime by regime
## ===========================================================================
VS <- do.call(rbind, lapply(names(PSI), function(nm) {
  p <- PSI[[nm]]
  data.frame(scenario = nm, regime = seq_len(nrow(p)), covariate = p$b1^2, latent = p$v,
             sill = p$so, nugget = p$se, stringsAsFactors = FALSE)
}))
VS$scenario <- factor(VS$scenario, levels = c("S2", SC$name))
VS <- VS[order(VS$scenario, VS$regime), ]
src <- c("covariate", "latent", "sill", "nugget")
src_col <- c("#5b8c5a", "#1f6f8b", "#e0a458", "#a8516e")
pdf_open("fig_sim_variance.pdf", 11, 5.6)
par_base(mar = c(6.5, 4, 2.5, 9), xpd = NA)
M <- t(as.matrix(VS[, src]))
gap <- unlist(lapply(levels(VS$scenario), function(s) c(0.9, rep(0.15, sum(VS$scenario == s) - 1))))
bp <- graphics::barplot(M, col = src_col, border = NA, space = gap, ylab = "variance of the response",
                        names.arg = VS$regime, cex.names = 0.75, ylim = c(0, 2.6))
mids <- tapply(bp, VS$scenario, mean)
graphics::text(mids, -0.32, names(mids), cex = 0.8)
graphics::legend(max(bp) + 1.4, 2.4, bty = "n", fill = rev(src_col), border = NA, cex = 0.8,
                 legend = rev(expression(paste("covariate, ", beta[1]^2), paste("latent, ", v),
                                         paste("partial sill, ", sigma[omega]^2),
                                         paste("nugget, ", sigma[epsilon]^2))))
graphics::title("The variance of the response of every regime (1, 2, 3) by source", cex.main = 0.95, font.main = 1)
pdf_close("fig_sim_variance.pdf")


## ===========================================================================
## Table: the parameter values of every regime
## ===========================================================================
row_tex <- function(nm, p, g) {
  s2 <- p$se[g] + p$so[g]
  sprintf("%s & %d & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
          if (g == 1L) nm else "", g, fmt(p$b0[g]), fmt(p$b1[g]), fmt(p$G[g], 1), fmt(p$v[g], 3),
          fmt(p$s2eta[g], 3), fmt(p$se[g]), fmt(p$so[g], 3), fmt(p$R[g], 1), fmt(p$theta[g], 2),
          fmt(s2), fmt(p$se[g] / s2))
}
lines <- c("\\begin{tabular}{lc rrrrr rrrr rr}", "\\toprule",
           "scenario & $g$ & $\\beta_0$ & $\\beta_1$ & $G$ & $v$ & $\\sigma^2_\\eta$ & $\\sigma^2_\\epsilon$ & $\\sigma^2_\\omega$ & $R$ & $\\theta$ & $s^2$ & nugget share \\\\",
           "\\midrule")
for (nm in c("S2", SC$name)) {
  p <- PSI[[nm]]
  lines <- c(lines, vapply(seq_len(nrow(p)), function(g) row_tex(nm, p, g), character(1)))
  if (nm != utils::tail(SC$name, 1)) lines <- c(lines, "\\addlinespace")
}
tex_write(c(lines, "\\bottomrule", "\\end{tabular}"), "tab_sim_scenarios.tex")


## ===========================================================================
## Table: the four geometries
## ===========================================================================
tex_write(c("\\begin{tabular}{lcl rrrrr}", "\\toprule",
            "geometry & $\\omega$ & fixed spread & centre distance & regime sd & distance in sd & nearer another centre & network sd \\\\",
            "\\midrule",
            sprintf("%s & %s & %s & %s & %s & %s & %s\\%% & %s \\\\", GEO$label, fmt(GEO$omega, 1),
                    ifelse(GEO$spread == "total", "network", "regime"), fmt(2 * GEO$omega), fmt(sqrt(GEO$nu)),
                    fmt(2 * GEO$omega / sqrt(GEO$nu)), fmt(100 * GEO$share_other, 0), fmt(GEO$network_sd)),
            "\\bottomrule", "\\end{tabular}"), "tab_sim_geometry.tex")


## ===========================================================================
## Table: the separation of the regimes seen by the assignment step
## ===========================================================================
## D_gh = 1/2 [ (s2_g + Delta2_gh) / s2_h - 1 - log(s2_g / s2_h) ],
## Delta2_gh = (b0_g - b0_h)^2 + (b1_g - b1_h)^2 + E b^2,
## E b^2 = v_g + v_h - 2 rho kappa_gh sqrt(v_g v_h),
## kappa_gh = sqrt((1 - G_g^2)(1 - G_h^2)) / (1 - G_g G_h)
sep_D <- function(p, rho, g, h) {
  s2 <- p$se + p$so
  kap <- sqrt((1 - p$G[g]^2) * (1 - p$G[h]^2)) / (1 - p$G[g] * p$G[h])
  lat <- p$v[g] + p$v[h] - 2 * rho * kap * sqrt(p$v[g] * p$v[h])
  mn  <- (p$b0[g] - p$b0[h])^2 + (p$b1[g] - p$b1[h])^2
  c(D = 0.5 * ((s2[g] + mn + lat) / s2[h] - 1 - log(s2[g] / s2[h])), mean = mn, latent = lat,
    rk = rho * kap)
}
SEP <- do.call(rbind, lapply(setdiff(SIM_SCEN$id, "S2"), function(id) {
  row <- SIM_SCEN[SIM_SCEN$id == id, , drop = FALSE]
  p <- dgp_psi(row)
  pr <- list(c(1, 2), c(2, 1), c(2, 3), c(3, 2))
  v <- lapply(pr, function(q) sep_D(p, row$rho, q[1], q[2]))
  Ds <- vapply(v, `[[`, 1, "D")
  j <- which.min(Ds)
  data.frame(id = id, rho = row$rho, D12 = Ds[1], D21 = Ds[2], D23 = Ds[3], D32 = Ds[4],
             minD = Ds[j], mean = v[[j]][["mean"]], latent = v[[j]][["latent"]],
             T60 = 60 * Ds[j], T120 = 120 * Ds[j], T365 = 365 * Ds[j], paper = id %in% SIM_SCENARIOS,
             stringsAsFactors = FALSE)
}))
SEP <- SEP[order(!SEP$paper, match(SEP$id, SIM_SCENARIOS)), ]
tex_write(c("\\begin{tabular}{lc rrrr r rr rrr}", "\\toprule",
            "scenario-variant & $\\rho$ & $D_{12}$ & $D_{21}$ & $D_{23}$ & $D_{32}$ & $\\min D$ & mean & latent & $60\\min D$ & $120\\min D$ & $365\\min D$ \\\\",
            "\\midrule",
            sprintf("%s & %d & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\%s",
                    sprintf("\\texttt{%s}", SEP$id), as.integer(SEP$rho), fmt(SEP$D12, 3), fmt(SEP$D21, 3),
                    fmt(SEP$D23, 3), fmt(SEP$D32, 3), fmt(SEP$minD, 3), fmt(SEP$mean), fmt(SEP$latent, 3),
                    fmt(SEP$T60, 1), fmt(SEP$T120, 1), fmt(SEP$T365, 1),
                    ifelse(seq_len(nrow(SEP)) == sum(SEP$paper), " \\midrule", "")),
            "\\bottomrule", "\\end{tabular}"), "tab_sim_separation.tex")


## ===========================================================================
## Table: the blocks of the design and their cells
## ===========================================================================
cfg_all <- CFG
all_cells <- sim_cells(cfg_all, SIM_BLOCKS)
blk_desc <- c(core    = "$n \\in \\{40, 100, 200\\}$, $T \\in \\{60, 120, 365\\}$, $\\omega \\in \\{0, 0.7, 1\\}$",
              spread  = "$n$ and $T$ as the core, $\\omega = 1$ with the spread of a regime fixed",
              balance = "regime sizes $1:2:3$",
              knn     = "$m \\in \\{3, 10\\}$ neighbours",
              n400    = "$n = 400$")
cnt <- table(factor(all_cells$block, levels = names(SIM_BLOCKS)))
tex_write(c("\\begin{tabular}{ll r}", "\\toprule",
            "block & margins (the others at the reference) & new cells \\\\", "\\midrule",
            sprintf("%s & %s & %d \\\\", names(SIM_BLOCKS), blk_desc[names(SIM_BLOCKS)], as.integer(cnt)),
            "\\midrule", sprintf("total & & %d \\\\", nrow(all_cells)),
            "\\bottomrule", "\\end{tabular}"), "tab_sim_blocks.tex")
cat("\ncells by block (each cell counted under the first block that has it):\n")
print(cnt)
cat("\nthe four geometries:\n")
print(GEO[, c("label", "omega", "spread", "nu", "share_other", "network_sd")], digits = 3, row.names = FALSE)
cat("\ndone\n")
