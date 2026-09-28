## ---------------------------------------------------------------------------
## Figures and tables of the SC-STEM paper, built from the caches written by
## 01-application.R and 02-simulation.R. Every block is guarded, so the script
## produces whatever the caches currently support and reports what is missing.
##
##   Rscript dev/paper/03-figures.R
## ---------------------------------------------------------------------------

PKG <- "C:/Users/paulm/OneDrive/Documenti/GitHub/Stem"
suppressMessages(pkgload::load_all(PKG, quiet = TRUE))

CACHE  <- file.path("dev", "paper", "cache")
FIGDIR <- "C:/Users/paulm/Dropbox/Applicazioni/Overleaf/SC-STEM package paper/Figures"
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

have <- function(x) file.exists(file.path(CACHE, paste0(x, ".rds")))
grab <- function(x) readRDS(file.path(CACHE, paste0(x, ".rds")))
fig  <- function(name, w, h) grDevices::pdf(file.path(FIGDIR, name),
                                            width = w, height = h)

## a colour-blind-safe pair, matching the package logo
COL <- c("#1f6f8b", "#e0a458", "#5b8c5a", "#a8516e")
GREY <- "#8a939f"

data(povalley)
coords <- povalley$coords
alt    <- povalley$altitude

## ===========================================================================
## Application
## ===========================================================================

if (have("grid") && have("selection")) {

  grid <- grab("grid")
  sel  <- grab("selection")
  best <- sel$fit
  g    <- best$group
  k    <- sel$k_selected

  ## ---- Figure: the estimated regimes on the network ------------------------
  nbg <- Stem:::scstem_neighbors(coords, knn = best$input_args$knn)
  W   <- nbg$W

  fig("app_regimes.pdf", 7.2, 4.2)
  graphics::par(mar = c(3.6, 3.6, 1.2, 0.6), mgp = c(2.2, 0.7, 0), las = 1)
  graphics::plot(coords, type = "n", xlab = "Longitude", ylab = "Latitude")
  ## edges: within a regime in its colour, cut edges faint and dashed
  for (i in seq_len(nrow(W) - 1L)) for (j in (i + 1L):ncol(W)) {
    if (W[i, j] == 0) next
    if (g[i] == g[j]) {
      graphics::segments(coords[i, 1], coords[i, 2], coords[j, 1], coords[j, 2],
                         col = grDevices::adjustcolor(COL[g[i]], 0.55), lwd = 1.6)
    } else {
      graphics::segments(coords[i, 1], coords[i, 2], coords[j, 1], coords[j, 2],
                         col = grDevices::adjustcolor(GREY, 0.8), lwd = 1, lty = 2)
    }
  }
  graphics::points(coords, pch = 21, bg = COL[g], col = "white", cex = 1.9, lwd = 1.2)
  graphics::legend("topleft", bty = "n",
                   legend = c(paste0("regime ", seq_len(k),
                                     " (", as.integer(table(g)), " stations)"),
                              "edge cut by the partition"),
                   pch = c(rep(21, k), NA), pt.bg = c(COL[seq_len(k)], NA),
                   col = c(rep("white", k), GREY),
                   lty = c(rep(NA, k), 2), lwd = c(rep(NA, k), 1), pt.cex = 1.6)
  grDevices::dev.off()
  message("wrote app_regimes.pdf")

  ## ---- Figure: model selection ---------------------------------------------
  tb <- grid$table
  fig("app_selection.pdf", 7.6, 3.4)
  graphics::par(mfrow = c(1, 2), mar = c(3.8, 4.0, 1.6, 0.8),
                mgp = c(2.4, 0.7, 0), las = 1)

  ok <- tb[tb$admissible, , drop = FALSE]
  phis <- sort(unique(ok$phi))
  graphics::plot(range(ok$k), range(ok$BIC), type = "n",
                 xlab = "number of regimes k", ylab = "BIC",
                 main = "(a) BIC over the grid", cex.main = 1)
  for (i in seq_along(phis)) {
    s <- ok[ok$phi == phis[i], ]
    s <- s[order(s$k), ]
    graphics::lines(s$k, s$BIC, col = COL[(i - 1) %% length(COL) + 1], lwd = 1.8)
    graphics::points(s$k, s$BIC, col = COL[(i - 1) %% length(COL) + 1], pch = 19, cex = 0.8)
  }
  bad <- tb[!tb$admissible, , drop = FALSE]
  if (nrow(bad)) graphics::points(bad$k, bad$BIC, pch = 4, col = GREY, cex = 1.1)
  graphics::legend("topright", bty = "n", cex = 0.8,
                   legend = c(paste0("phi = ", phis), "inadmissible"),
                   col = c(COL[seq_along(phis)], GREY), lwd = c(rep(1.8, length(phis)), NA),
                   pch = c(rep(19, length(phis)), 4))

  st <- sel$step2
  if (!is.null(st) && nrow(st)) {
    ## the grid is roughly logarithmic, so the penalties are equally spaced
    ## and labelled with their values
    x <- seq_len(nrow(st))
    graphics::plot(x, st$value, type = "b", pch = 19, lwd = 1.8, xaxt = "n",
                   col = COL[1], xlab = expression(phi), ylab = "BIC",
                   main = bquote("(b) BIC over" ~ phi ~ "at" ~ hat(k) == .(sel$k_selected)),
                   cex.main = 1)
    graphics::axis(1, at = x, labels = format(st$phi, drop0trailing = TRUE))
    graphics::abline(v = x[st$selected], lty = 2, col = COL[2], lwd = 1.6)
    graphics::legend("topleft", bty = "n", cex = 0.8,
                     legend = bquote(hat(phi) == .(sel$phi_selected)),
                     lty = 2, col = COL[2])
  }
  grDevices::dev.off()
  message("wrote app_selection.pdf")

  ## ---- Figure: the latent process of each regime ---------------------------
  ok_fits <- which(best$final_refit)
  if (length(ok_fits)) {
    fig("app_latent.pdf", 7.2, 3.2)
    graphics::par(mar = c(3.6, 3.8, 1.2, 0.6), mgp = c(2.3, 0.7, 0), las = 1)
    ys <- lapply(ok_fits, function(gg) as.numeric(best$fit_list[[gg]]$estimates$y.smoothed))
    n1 <- min(365L, length(ys[[1]]))
    dts <- povalley$dates[seq_len(n1)]
    graphics::plot(dts, ys[[1]][seq_len(n1)], type = "n",
                   ylim = range(unlist(lapply(ys, function(v) v[seq_len(n1)]))),
                   xlab = "", ylab = "smoothed latent state")
    for (i in seq_along(ok_fits)) {
      graphics::lines(dts, ys[[i]][seq_len(n1)], col = COL[ok_fits[i]], lwd = 1.5)
    }
    Gs <- vapply(ok_fits, function(gg) best$phi_hat[gg, "G"], numeric(1))
    graphics::legend("topright", bty = "n", cex = 0.85,
                     legend = sprintf("regime %d  (G = %.3f)", ok_fits, Gs),
                     col = COL[ok_fits], lwd = 1.8)
    grDevices::dev.off()
    message("wrote app_latent.pdf")
  }

  ## ---- Table: the grid ------------------------------------------------------
  tt <- tb[order(tb$k, tb$phi), c("k", "phi", "k_eff", "admissible",
                                  "loglik", "npar", "BIC")]
  con <- file(file.path(FIGDIR, "tab_grid.tex"), "w")
  writeLines(c("\\begin{tabular}{rrrlrrr}", "\\toprule",
               "$k$ & $\\phi$ & $k_{\\mathrm{eff}}$ & admissible & $\\ell$ & $\\mathcal{K}$ & BIC \\\\",
               "\\midrule"), con)
  for (i in seq_len(nrow(tt))) {
    writeLines(sprintf("%d & %.2f & %d & %s & %.1f & %d & %.1f \\\\",
                       tt$k[i], tt$phi[i], tt$k_eff[i],
                       ifelse(tt$admissible[i], "yes", "no"),
                       tt$loglik[i], tt$npar[i], tt$BIC[i]), con)
  }
  writeLines(c("\\bottomrule", "\\end{tabular}"), con)
  close(con)
  message("wrote tab_grid.tex")

} else {
  message("application caches not ready: skipping the application figures")
}

## ---- bootstrap-dependent output ---------------------------------------------
if (have("bootinference")) {
  inf <- grab("bootinference")
  sel <- grab("selection")
  best <- sel$fit
  pooled <- if (have("pooled")) grab("pooled") else NULL

  keep <- c("beta1", "beta2", "beta3", "sigma2eps", "sigma2omega",
            "theta", "G", "Sigmaeta")
  s <- inf$summary[inf$summary$parameter %in% keep, , drop = FALSE]

  lab <- c(beta1 = "intercept", beta2 = "altitude", beta3 = "PM10",
           sigma2eps = "nugget", sigma2omega = "partial sill",
           theta = "range parameter", G = "persistence G",
           Sigmaeta = "innovation variance")

  fig("app_params.pdf", 7.6, 4.6)
  graphics::par(mfrow = c(2, 4), mar = c(2.6, 3.4, 2.0, 0.6),
                mgp = c(2.2, 0.6, 0), las = 1)
  for (pnm in keep) {
    ss <- s[s$parameter == pnm, , drop = FALSE]
    if (!nrow(ss)) next
    lo <- ss[["perc_lower"]]; hi <- ss[["perc_upper"]]
    if (is.null(lo)) { lo <- ss$estimate - 2 * ss$se; hi <- ss$estimate + 2 * ss$se }
    yl <- range(c(lo, hi, ss$estimate), na.rm = TRUE)
    graphics::plot(seq_len(nrow(ss)), ss$estimate, ylim = yl, xaxt = "n",
                   xlim = c(0.5, nrow(ss) + 0.5), pch = 19, cex = 1.3,
                   col = COL[seq_len(nrow(ss))], xlab = "",
                   ylab = "", main = lab[[pnm]], cex.main = 1)
    graphics::arrows(seq_len(nrow(ss)), lo, seq_len(nrow(ss)), hi,
                     angle = 90, code = 3, length = 0.04,
                     col = COL[seq_len(nrow(ss))], lwd = 1.8)
    graphics::axis(1, at = seq_len(nrow(ss)), labels = paste0("R", seq_len(nrow(ss))))
    if (!is.null(pooled)) {
      pv <- unlist(pooled$estimates$phi.hat)[pnm]
      if (!is.na(pv)) graphics::abline(h = pv, lty = 2, col = GREY)
    }
  }
  grDevices::dev.off()
  message("wrote app_params.pdf")

  ## co-clustering heat map
  cc <- inf$coclustering
  ord <- order(best$group)
  fig("app_coclust.pdf", 4.4, 4.2)
  graphics::par(mar = c(3.2, 3.2, 1.4, 1.0), mgp = c(2, 0.7, 0), las = 1)
  graphics::image(seq_len(nrow(cc)), seq_len(ncol(cc)), cc[ord, rev(ord)],
                  col = grDevices::hcl.colors(64, "Blues", rev = TRUE),
                  xlab = "stations, ordered by regime", ylab = "",
                  main = "bootstrap co-clustering", cex.main = 1)
  grDevices::dev.off()
  message("wrote app_coclust.pdf")
} else {
  message("bootstrap cache not ready: skipping the parameter figure")
}

## ===========================================================================
## Simulation
## ===========================================================================
sims <- list.files(CACHE, pattern = "^sim_.*[.]rds$", full.names = TRUE)
if (length(sims)) {
  rows <- list(); selrows <- list()
  for (f in sims) {
    x <- readRDS(f)
    des <- x$scenario$design; lev <- x$scenario$level
    for (r in x$res) {
      if (!is.null(r$grid)) {
        rows[[length(rows) + 1L]] <- data.frame(design = des, level = lev,
                                                phi = r$grid$phi, ari = r$grid$ari)
      }
      if (!is.null(r$k_selected)) {
        selrows[[length(selrows) + 1L]] <- data.frame(
          design = des, level = lev, k_selected = r$k_selected,
          phi_selected = r$phi_selected, ari = r$ari_selected)
      }
    }
  }
  ari <- do.call(rbind, rows)
  saveRDS(ari, file.path(CACHE, "sim_ari_long.rds"))

  fig("sim_ari.pdf", 7.6, 4.4)
  designs <- sort(unique(ari$design))
  graphics::par(mfrow = c(2, length(designs)), mar = c(3.2, 3.6, 2.0, 0.6),
                mgp = c(2.2, 0.7, 0), las = 1)
  for (lev in c("clear", "poor")) for (des in designs) {
    sub <- ari[ari$design == des & ari$level == lev, ]
    if (!nrow(sub)) { graphics::plot.new(); next }
    graphics::boxplot(ari ~ phi, data = sub, ylim = c(0, 1),
                      col = grDevices::adjustcolor(COL[1], 0.35),
                      border = COL[1], outcex = 0.4,
                      xlab = expression(phi), ylab = "Adjusted Rand Index",
                      main = paste0(des, " - ", lev), cex.main = 1)
  }
  grDevices::dev.off()
  message("wrote sim_ari.pdf")

  if (length(selrows)) {
    sl <- do.call(rbind, selrows)
    tabk <- table(paste(sl$design, sl$level), sl$k_selected)
    con <- file(file.path(FIGDIR, "tab_selection.tex"), "w")
    ks <- colnames(tabk)
    writeLines(c(paste0("\\begin{tabular}{l", paste(rep("r", length(ks)), collapse = ""), "}"),
                 "\\toprule",
                 paste0("design & ", paste0("$k=", ks, "$", collapse = " & "), " \\\\"),
                 "\\midrule"), con)
    for (i in seq_len(nrow(tabk))) {
      pr <- round(100 * tabk[i, ] / sum(tabk[i, ]), 1)
      writeLines(paste0(rownames(tabk)[i], " & ",
                        paste0(sprintf("%.1f", pr), collapse = " & "), " \\\\"), con)
    }
    writeLines(c("\\bottomrule", "\\end{tabular}"), con)
    close(con)
    message("wrote tab_selection.tex")

    cat("\nSelection frequencies (%):\n")
    print(round(100 * prop.table(tabk, 1), 1))
  }

  cat("\nMedian ARI by design, level and penalty:\n")
  print(round(tapply(ari$ari, list(paste(ari$design, ari$level), ari$phi),
                     stats::median, na.rm = TRUE), 3))
} else {
  message("no simulation caches: skipping the simulation figures")
}

message("\nDone.")
