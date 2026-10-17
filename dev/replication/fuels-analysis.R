## ===========================================================================
## SC-STEM fuel application: THE ANALYSIS OF THE RESULTS
##
## Reads the outputs of the main scripts (stage "final" of the three settings),
## of fuels-wald-bootstrap.R and of fuels-check-optima.R, and writes the tables
## (.tex, .csv) and the figures (.png) of the application section of the paper
## and of its Application Supplement into <root>/analysis/. Nothing is fitted.
##
##     Rscript fuels-analysis.R --root=<the fuels folder>
##
## The questions the outputs answer:
##   Q1 local leadership: in which regimes do the lagged relative prices of the
##      neighbouring pump X predict those of pump Y (bootstrap Wald test of the
##      eight lags of X), how strongly (the long-run pass-through of a
##      permanent deviation of X), and what the pooled model would have said;
##   Q2 what the regimes are: their number, spatial penalty, brands, distance
##      to the neighbour, stability under the bootstrap;
##   Q3 the pass-through of the excise duties, regime by regime (the pulses);
##   Q4 the small samples of setting 2 (the independents, 6 to 71 per city);
## and the checks: the pump-by-pump tests, the residual diagnostics, the
## local maxima of the selected fits.
##
## This is part of the REPLICATION MATERIAL of the paper, not of the Stem
## package.
## ===========================================================================

suppressPackageStartupMessages(library(data.table))
HERE <- local({
  a <- commandArgs(trailingOnly = FALSE); m <- grep("^--file=", a, value = TRUE)
  f <- if (length(m)) normalizePath(sub("^--file=", "", m[1]), winslash = "/") else NULL
  if (is.null(f)) for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) { f <- normalizePath(of, winslash = "/"); break }
  }
  if (is.null(f)) normalizePath(getwd(), winslash = "/") else dirname(f)
})
source(file.path(HERE, "fuels-functions.R"))
ROOT <- local({
  a <- grep("^--root=", commandArgs(trailingOnly = TRUE), value = TRUE)
  normalizePath(if (length(a)) sub("^--root=", "", a[1]) else file.path(HERE, "fuels"), winslash = "/", mustWork = TRUE)
})
OUT <- file.path(ROOT, "analysis"); dir.create(OUT, showWarnings = FALSE)
SHA7 <- substr(sub("^.*@", "", FU_STEM_REF), 1L, 7L)
KEY <- FU_BT$lag_key
SETTING_LAB <- c("1" = "Y major, X independent", "2" = "Y independent, X major", "3" = "Y major, X other major")
CITY_ORDER <- c("Torino", "Milano", "Venezia", "Bologna", "Firenze", "Roma", "Napoli", "Bari", "Palermo",
                "Catania", "Messina")
## the names printed in the tables and figures
CITY_EN <- c(Torino = "Turin", Milano = "Milan", Venezia = "Venice", Bologna = "Bologna", Firenze = "Florence",
             Roma = "Rome", Napoli = "Naples", Bari = "Bari", Palermo = "Palermo", Catania = "Catania",
             Messina = "Messina")
en <- function(cc) unname(CITY_EN[as.character(cc)])
## the caption, numbered by longtable, and the header of a longtable, repeated
## on every page
lt_head <- function(cap, lab, header)
  c(sprintf("\\caption{%s}\\label{%s}\\\\", cap, lab), "\\toprule", header, "\\midrule", "\\endfirsthead",
    "\\toprule", header, "\\midrule", "\\endhead")
PAL3 <- c("#1f6f8b", "#e0a458", "#a8516e")
PAL5 <- c("#1b2430", "#1f6f8b", "#e0a458", "#5b8c5a", "#a8516e")
EVENTS <- c(pulse_excise_cut = "excise cut\n22 Mar 2022", pulse_cut_reduced = "cut reduced\n1 Dec 2022",
            pulse_cut_ends = "cut ends\n1 Jan 2023", pulse_realignment = "realignment\n15 May 2025",
            pulse_equal_excise = "equal excise\n1 Jan 2026")

fig_open <- function(name, width, height)
  grDevices::png(file.path(OUT, name), width = width, height = height, units = "in", res = 300,
                 type = if (isTRUE(capabilities("cairo"))) "cairo" else "windows")
fig_close <- function(name) { invisible(grDevices::dev.off()); message("wrote ", name) }
## minus signs in math mode, and no negative zero
tex_write <- function(lines, name) {
  lines <- gsub("(?<![0-9A-Za-z.$])-(0\\.0+)(?![0-9])", "\\1", lines, perl = TRUE)
  lines <- gsub("(?<![0-9A-Za-z.$-])-(?=[0-9])", "$-$", lines, perl = TRUE)
  writeLines(lines, file.path(OUT, name)); message("wrote ", name)
}
csv_write <- function(d, name) { utils::write.csv(d, file.path(OUT, name), row.names = FALSE); message("wrote ", name) }
rd <- function(f) if (file.exists(f)) utils::read.csv(f, stringsAsFactors = FALSE) else NULL
pct <- function(x) round(100 * mean(x, na.rm = TRUE))
fmt_p <- function(p) ifelse(is.na(p), "--", ifelse(p < 0.001, "$<$0.001", sprintf("%.3f", p)))


## ===========================================================================
## 1. The outputs of every model
## ===========================================================================
CASE <- vapply(1:3, function(s) fu_case_dir(ROOT, s), "")
SM <- do.call(rbind, lapply(1:3, function(s) cbind(setting = s, rd(file.path(CASE[s], "nn/final/summary.csv")))))
SM$city <- factor(SM$city, CITY_ORDER)
fdir <- function(s, m) file.path(CASE[s], "nn", "final", m)

## the regimes, one row per regime of the selected fit
REG <- do.call(rbind, lapply(seq_len(nrow(SM)), function(i) {
  s <- SM$setting[i]; m <- SM$model[i]
  g <- rd(file.path(fdir(s, m), "regimes.csv"))
  w <- rd(file.path(fdir(s, m), "wald_bootstrap.csv"))
  w <- w[startsWith(w$model, "selected"), ]
  pw <- function(b) w$p_value[w$block == b][order(w$regime[w$block == b])]
  b0 <- readRDS(file.path(CASE[s], "nn", "cache", SHA7, sprintf("%s_%s_boot_lambda0.rds", m, KEY)))
  A <- b0$inf$aligned; sm <- as.data.frame(b0$inf$summary)
  K <- SM$K[i]
  ## the long-run pass-through of a permanent deviation of X (same fuel, and
  ## both fuels), sum of its lags over one minus the sum of the own lags, with
  ## the percentile interval of its bootstrap draws
  lrm <- t(vapply(seq_len(K), function(g) {
    e <- function(ix) sum(sm$estimate[sm$cluster == g & sm$parameter %in% paste0("beta", ix)])
    d <- function(ix) rowSums(matrix(A[, g, paste0("beta", ix)], ncol = length(ix)))
    own <- 2:5; nb <- 10:13; nbo <- 14:17
    l1 <- d(nb) / (1 - d(own)); l2 <- (d(nb) + d(nbo)) / (1 - d(own))
    c(e(nb) / (1 - e(own)), stats::quantile(l1, c(.025, .975), na.rm = TRUE),
      (e(nb) + e(nbo)) / (1 - e(own)), stats::quantile(l2, c(.025, .975), na.rm = TRUE), e(own))
  }, numeric(7)))
  br <- table(factor(g$regime, seq_len(K)), g$y_brand)
  top <- apply(br, 1, function(r) names(r)[which.max(r)])
  data.frame(setting = s, model = m, city = SM$city[i], fuel = SM$fuel[i], K = K, phi = SM$phi[i],
             regime = seq_len(K), pumps = as.integer(table(factor(g$regime, seq_len(K)))),
             p_nb_both = pw("nb_both"), p_nb = pw("nb"), p_nb_other = pw("nb_other"), p_own_other = pw("own_other"),
             lrm_nb = lrm[, 1], lrm_nb_lo = lrm[, 2], lrm_nb_hi = lrm[, 3],
             lrm_both = lrm[, 4], lrm_both_lo = lrm[, 5], lrm_both_hi = lrm[, 6], sum_own = lrm[, 7],
             dist_median = as.numeric(tapply(g$nn_dist_km, factor(g$regime, seq_len(K)), stats::median)),
             top_brand = top, top_share = as.numeric(apply(br, 1, max) / rowSums(br)),
             eni_share = if ("Agip Eni" %in% colnames(br)) as.numeric(br[, "Agip Eni"] / rowSums(br)) else 0,
             theta_bound = as.character(tapply(g$theta_bound, factor(g$regime, seq_len(K)), `[`, 1)),
             ari_boot = stats::median(b0$inf$stability$ARI), B_valid = b0$B_valid,
             stringsAsFactors = FALSE)
}))
REG$lead <- REG$p_nb_both < 0.05

## the pooled model of every model with K* > 1 (and the selected one when K* = 1)
POOL <- do.call(rbind, lapply(seq_len(nrow(SM)), function(i) {
  s <- SM$setting[i]; m <- SM$model[i]
  w <- rd(file.path(fdir(s, m), "wald_bootstrap.csv"))
  w <- if (SM$K[i] == 1) w[startsWith(w$model, "selected"), ] else w[w$model == "pooled", ]
  data.frame(setting = s, model = m, city = SM$city[i], fuel = SM$fuel[i], block = w$block, p = w$p_value)
}))

## the joint tests of fuels-wald-bootstrap.R (row "all"), when available
TW <- do.call(rbind, lapply(seq_len(nrow(SM)), function(i) {
  x <- rd(file.path(fdir(SM$setting[i], SM$model[i]), "test_wald.csv"))
  if (is.null(x)) NULL else cbind(city = SM$city[i], fuel = SM$fuel[i], x)
}))

## the pulses of the excise duties, one row per regime and pulse
PUL <- do.call(rbind, lapply(seq_len(nrow(SM)), function(i) {
  s <- SM$setting[i]; m <- SM$model[i]
  e <- rd(file.path(fdir(s, m), "excise_pulses.csv")); e <- e[startsWith(e$model, "selected"), ]
  data.frame(setting = s, model = m, city = SM$city[i], fuel = SM$fuel[i], regime = e$cluster, pulse = e$pulse,
             est = e$estimate, se = e$se, lo = e$perc_lower_95_CI, hi = e$perc_upper_95_CI)
}))
PUL$sig <- PUL$lo > 0 | PUL$hi < 0

## the pump-by-pump tests and their panel combination
PT <- do.call(rbind, lapply(seq_len(nrow(SM)), function(i) {
  s <- SM$setting[i]; m <- SM$model[i]
  p <- rd(file.path(fdir(s, m), "pump_tests.csv")); g <- rd(file.path(fdir(s, m), "regimes.csv"))
  data.frame(setting = s, model = m, city = SM$city[i], fuel = SM$fuel[i], block = p$block, p_F = p$p_F,
             p_hc1 = p$p_hc1, dist = g$nn_dist_km[match(p$y_site, g$y_site)],
             regime = g$regime[match(p$y_site, g$y_site)])
}))
DH <- do.call(rbind, lapply(seq_len(nrow(SM)), function(i) {
  d <- rd(file.path(fdir(SM$setting[i], SM$model[i]), "dumitrescu_hurlin.csv"))
  cbind(setting = SM$setting[i], model = SM$model[i], city = SM$city[i], fuel = SM$fuel[i], d)
}))
DG <- do.call(rbind, lapply(seq_len(nrow(SM)), function(i) {
  d <- rd(file.path(fdir(SM$setting[i], SM$model[i]), "diagnostics.csv"))
  cbind(setting = SM$setting[i], model = SM$model[i], d)
}))
OPT <- rd(file.path(ROOT, "checks", "optima_by_regime.csv"))
OPTC <- rd(file.path(ROOT, "checks", "optima_coefficients.csv"))
OPTW <- rd(file.path(ROOT, "checks", "optima_wald.csv"))
csv_write(REG, "app_regimes.csv"); csv_write(PUL, "app_pulses.csv")


## ===========================================================================
## 2. The main table: one row per setting and fuel
## ===========================================================================
main_rows <- do.call(rbind, lapply(1:3, function(s) do.call(rbind, lapply(c("gasoline", "diesel"), function(f) {
  r <- REG[REG$setting == s & REG$fuel == f, ]; sm <- SM[SM$setting == s & SM$fuel == f, ]
  po <- POOL[POOL$setting == s & POOL$fuel == f & POOL$block == "nb_both", ]
  dh <- DH[DH$setting == s & DH$fuel == f & DH$block == "nb_both", ]
  data.frame(setting = s, fuel = f, pumps = sum(sm$pumps), regimes = nrow(r),
             K = paste(range(sm$K), collapse = "--"),
             lead_regimes = pct(r$lead), lead_pumps = round(100 * sum(r$pumps[r$lead]) / sum(r$pumps)),
             pooled_lead = sum(po$p < 0.05), cities = nrow(po),
             lrm_med = stats::median(r$lrm_nb), lrm_q1 = stats::quantile(r$lrm_nb, .25), lrm_q3 = stats::quantile(r$lrm_nb, .75),
             lrm_pos = pct(r$lrm_nb_lo > 0),
             hc1 = stats::median(dh$reject_hc1_pct), F = stats::median(dh$reject_F_pct))
}))))
csv_write(main_rows, "tab_app_main.csv")
tex_write(c(
  "\\begin{tabular}{llrrcrrcccr}",
  "\\toprule",
  " & & & & & \\multicolumn{2}{c}{X leads Y (\\%)} & pooled & \\multicolumn{2}{c}{long-run pass-through} & per pump \\\\",
  "\\cmidrule{6-7}\\cmidrule{9-10}",
  "setting & fuel & pumps & regimes & $\\hat K$ & regimes & pumps & model & median (IQR) & $>0$ (\\%) & HC1 (\\%) \\\\",
  "\\midrule",
  unlist(lapply(seq_len(nrow(main_rows)), function(j) { x <- main_rows[j, ]
    sprintf("%s & %s & %d & %d & %s & %d & %d & %d/%d & %.2f (%.2f, %.2f) & %d & %.0f \\\\%s",
            if (x$fuel == "gasoline") as.character(x$setting) else "",
            x$fuel, x$pumps, x$regimes, x$K, x$lead_regimes, x$lead_pumps, x$pooled_lead, x$cities,
            x$lrm_med, x$lrm_q1, x$lrm_q3, x$lrm_pos, x$hc1, if (x$fuel == "diesel" && x$setting < 3) "[2pt]" else "") })),
  "\\bottomrule", "\\end{tabular}"), "tab_app_main.tex")


## ===========================================================================
## 3. Figure: leadership city by city
## ===========================================================================
city_tab <- do.call(rbind, lapply(1:3, function(s) do.call(rbind, lapply(c("gasoline", "diesel"), function(f)
  do.call(rbind, lapply(CITY_ORDER, function(cc) {
    r <- REG[REG$setting == s & REG$fuel == f & REG$city == cc, ]
    po <- POOL[POOL$setting == s & POOL$fuel == f & POOL$city == cc & POOL$block == "nb_both", ]
    if (!nrow(r)) return(NULL)
    data.frame(setting = s, fuel = f, city = cc, pumps = sum(r$pumps), K = nrow(r),
               lead_pumps = 100 * sum(r$pumps[r$lead]) / sum(r$pumps), pooled = po$p[1] < 0.05)
  }))))))
csv_write(city_tab, "app_leadership_by_city.csv")
fig_open("fig_app_leadership.png", 9, 4.6)
graphics::par(mfrow = c(1, 3), mar = c(3.6, 0.8, 3.2, 0.6), oma = c(2.6, 5.6, 1.6, 0), mgp = c(2.1, 0.6, 0))
for (s in 1:3) {
  x <- city_tab[city_tab$setting == s, ]
  yy <- match(x$city, rev(CITY_ORDER))
  plot(NA, xlim = c(0, 100), ylim = c(0.5, length(CITY_ORDER) + 0.5), yaxt = "n", xlab = "% of the pumps Y",
       ylab = "", main = sprintf("Setting %d\n%s", s, SETTING_LAB[as.character(s)]), cex.main = 1.05, font.main = 1)
  graphics::abline(h = seq_along(CITY_ORDER), col = "grey92")
  if (s == 1) graphics::axis(2, at = seq_along(CITY_ORDER), labels = en(rev(CITY_ORDER)), las = 1, cex.axis = 1.05)
  off <- ifelse(x$fuel == "gasoline", 0.15, -0.15)
  graphics::points(x$lead_pumps, yy + off, pch = ifelse(x$pooled, 19, 1), cex = 1.5,
                   col = ifelse(x$fuel == "gasoline", PAL3[1], PAL3[3]), lwd = 1.6)
}
graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
graphics::legend("bottom", horiz = TRUE, bty = "n", cex = 1.05, pt.cex = 1.4, pt.lwd = 1.6,
                 legend = c("gasoline", "diesel", "pooled model: X leads", "pooled model: no"),
                 pch = c(19, 19, 19, 1), col = c(PAL3[1], PAL3[3], "grey30", "grey30"))
graphics::mtext("Share of the pumps Y in the regimes where the neighbour X leads (bootstrap Wald test, 5%)",
                side = 3, line = -1.2, cex = 0.95)
fig_close("fig_app_leadership.png")


## ===========================================================================
## 4. Figure: how much passes through -- from the neighbour, from the excise
## ===========================================================================
fig_open("fig_app_passthrough.png", 9.6, 4.8)
graphics::layout(matrix(1:2, 1), widths = c(1, 1.4))
graphics::par(mar = c(5.2, 4.2, 2.6, 0.8), mgp = c(2.4, 0.7, 0))
grp <- interaction(REG$fuel, REG$setting, sep = " ", lex.order = TRUE)
grp <- factor(grp, c(outer(c("gasoline", "diesel"), 1:3, paste)))
graphics::boxplot(REG$lrm_nb ~ grp, outline = FALSE, col = rep(c(PAL3[1], PAL3[3]), 3), border = "grey25",
                  xaxt = "n", ylab = "cents per cent", xlab = "", ylim = range(REG$lrm_nb, na.rm = TRUE),
                  main = "(a) Long-run pass-through from X", cex.main = 1, font.main = 1)
set.seed(1); graphics::points(as.integer(grp) + stats::runif(length(grp), -0.18, 0.18), REG$lrm_nb, pch = 16, cex = 0.5,
                col = ifelse(REG$lrm_nb_lo > 0 | REG$lrm_nb_hi < 0, "grey10", "grey65"))
graphics::abline(h = 0, lty = 2)
graphics::mtext(rep(c("G", "D"), 3), side = 1, at = 1:6, line = 0.5, cex = 0.8)
graphics::mtext(paste("setting", 1:3), side = 1, at = c(1.5, 3.5, 5.5), line = 1.8, cex = 0.9)
ev <- PUL[PUL$pulse %in% names(EVENTS), ]
ev$pulse <- factor(ev$pulse, names(EVENTS))
eg <- interaction(ev$fuel, ev$pulse, sep = " ", lex.order = TRUE)
eg <- factor(eg, c(t(outer(names(EVENTS), c("gasoline", "diesel"), function(a, b) paste(b, a)))))
bs <- graphics::boxplot(ev$est ~ eg, plot = FALSE)
graphics::boxplot(ev$est ~ eg, outline = FALSE, col = rep(c(PAL3[1], PAL3[3]), length(EVENTS)), border = "grey25",
                  xaxt = "n", ylab = "cents per litre", xlab = "", ylim = range(bs$stats),
                  main = "(b) Relative price in the week of an excise change", cex.main = 1, font.main = 1)
graphics::abline(h = 0, lty = 2)
graphics::abline(v = seq(2.5, 2 * length(EVENTS) - 1, 2), col = "grey85")
for (j in seq_along(EVENTS)) {
  lab <- strsplit(EVENTS[[j]], "\n")[[1]]
  graphics::mtext(lab[1], side = 1, at = 2 * j - 0.5, line = 0.6, cex = 0.75)
  graphics::mtext(lab[2], side = 1, at = 2 * j - 0.5, line = 1.5, cex = 0.75)
}
graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
graphics::legend("bottom", horiz = TRUE, bty = "n", cex = 0.9, fill = c(PAL3[1], PAL3[3], NA, NA), border = c("grey25", "grey25", NA, NA),
                 pch = c(NA, NA, 16, 16), col = c(NA, NA, "grey10", "grey65"),
                 legend = c("gasoline", "diesel", "regime: interval excludes 0", "regime: interval includes 0"))
fig_close("fig_app_passthrough.png")


## ===========================================================================
## 5. Figure: the regimes of one city, in space and by brand
## ===========================================================================
PALM <- c("#1f6f8b", "#e0a458", "#a8516e", "#5b8c5a", "#6c5b7b")
draw_map <- function(s, m, title) {
  g <- rd(file.path(fdir(s, m), "regimes.csv")); r <- REG[REG$setting == s & REG$model == m, ]
  col <- PALM[seq_len(nrow(r))]
  graphics::par(mar = c(5.4, 0.8, 2.4, 0.8))
  plot(g$y_lon, g$y_lat, type = "n", asp = 1 / cos(mean(g$y_lat) * pi / 180), axes = FALSE, xlab = "", ylab = "",
       main = title, cex.main = 1, font.main = 1)
  graphics::box(col = "grey70")
  graphics::points(g$y_lon, g$y_lat, pch = ifelse(r$lead[g$regime], 19, 1), col = col[g$regime], cex = 0.7, lwd = 1.2)
  graphics::legend("bottom", inset = c(0, -0.27), xpd = NA, ncol = 2, bty = "n", cex = 0.78, pt.cex = 1.1,
                   pch = ifelse(r$lead, 19, 1), col = col,
                   legend = sprintf("%d: %d pumps, %s %.0f%%", r$regime, r$pumps,
                                    sub("Agip Eni", "Eni", r$top_brand), 100 * r$top_share))
}
for (cs in list(c(3, "RM_gasoline"), c(3, "MI_diesel"), c(1, "RM_gasoline"), c(3, "NA_gasoline"))) {
  s <- as.integer(cs[1]); m <- cs[2]
  fig_open(sprintf("fig_app_map_s%d_%s.png", s, m), 5, 5.6)
  draw_map(s, m, sprintf("%s, %s, setting %d", en(FU_CITY[[sub("_.*$", "", m)]]), sub("^.*_", "", m), s))
  fig_close(sprintf("fig_app_map_s%d_%s.png", s, m))
}
fig_open("fig_app_maps.png", 9, 5.4)
graphics::par(mfrow = c(1, 2))
draw_map(3, "RM_gasoline", "(a) Rome, gasoline, setting 3")
draw_map(3, "MI_diesel", "(b) Milan, diesel, setting 3")
fig_close("fig_app_maps.png")


## ===========================================================================
## 6. Supplement: the models, the regimes, the joint tests
## ===========================================================================
## the selected models
mod <- SM[order(SM$setting, SM$city, SM$fuel), ]
mod$sizes <- vapply(seq_len(nrow(mod)), function(i)
  paste(REG$pumps[REG$setting == mod$setting[i] & REG$model == mod$model[i]], collapse = "/"), "")
mod$ari_boot <- vapply(seq_len(nrow(mod)), function(i)
  REG$ari_boot[REG$setting == mod$setting[i] & REG$model == mod$model[i]][1], 1)
csv_write(mod, "tab_app_models.csv")
tex_write(c("\\begin{longtable}{llrrrlrrr}",
  lt_head(paste("The selected models. $K_{\\max}$, the largest number of regimes of the grid; $\\hat K$, $\\hat\\phi$",
                "and the sizes of the regimes of the selected fit; ARI brands, the adjusted Rand index between the regimes",
                "and the brands of the stations $Y$ (not defined in setting 2, where they are all independent); usable",
                "draws, of $B = 100$; ARI draws, the median adjusted Rand index between the partitions of the bootstrap",
                "draws and the estimated one."), "tab:sc:models",
          "setting & city & fuel & pumps & $K_{\\max}$ & $\\hat K$, $\\hat\\phi$, sizes & ARI brands & usable draws & ARI draws \\\\"),
  sprintf("%d & %s & %s & %d & %d & %d, %s, %s & %s & %d & %.2f \\\\", mod$setting, en(mod$city), mod$fuel, mod$pumps, mod$K_max,
          mod$K, format(mod$phi), mod$sizes, ifelse(is.na(mod$ari_regimes_brands), "--", sprintf("%.2f", mod$ari_regimes_brands)),
          mod$B_valid, mod$ari_boot),
  "\\bottomrule", "\\end{longtable}"), "tab_app_models.tex")

## the regimes, setting by setting
for (s in 1:3) {
  r <- REG[REG$setting == s, ]; r <- r[order(r$city, r$fuel, r$regime), ]
  tex_write(c("\\begin{longtable}{llrrlrrrrrr}",
    lt_head(paste(sprintf("The regimes of setting %d.", s),
                  "Brand: the most frequent brand of the stations $Y$ of the regime, with its share; X: median distance",
                  "to the neighbour, in kilometres; $p$-values of the bootstrap Wald test of the eight lags of $X$, of its",
                  "four lags on the same fuel, of its four lags on the other fuel, and of the four lags of the other fuel",
                  "of $Y$; long-run pass-through from $X$ on the same fuel, with its bootstrap 95\\% interval."),
            sprintf("tab:sc:reg%d", s),
            c("city & fuel & regime & pumps & brand (\\%) & X (km) & \\multicolumn{4}{c}{bootstrap Wald, $p$-value} & long-run pass-through \\\\",
              "\\cmidrule(lr){7-10}", " & & & & & & X & X same & X other & Y other & from X, same fuel \\\\")),
    sprintf("%s & %s & %d & %d & %s %.0f & %.1f & %s & %s & %s & %s & %.2f (%.2f, %.2f) \\\\", en(r$city), r$fuel, r$regime, r$pumps,
            sub("Agip Eni", "Eni", gsub("Pompe Bianche", "indep.", r$top_brand)), 100 * r$top_share, r$dist_median, fmt_p(r$p_nb_both), fmt_p(r$p_nb),
            fmt_p(r$p_nb_other), fmt_p(r$p_own_other), r$lrm_nb, r$lrm_nb_lo, r$lrm_nb_hi),
    "\\bottomrule", "\\end{longtable}"), sprintf("tab_app_regimes_s%d.tex", s))
}

## leadership by the brand of the regime and by the distance to X
brand_tab <- do.call(rbind, lapply(c(1, 3), function(s) {
  r <- REG[REG$setting == s, ]
  data.frame(setting = s, group = c("Eni > 50%", "other regimes"), regimes = c(sum(r$eni_share > .5), sum(r$eni_share <= .5)),
             lead = c(pct(r$lead[r$eni_share > .5]), pct(r$lead[r$eni_share <= .5])),
             lrm = c(stats::median(r$lrm_nb[r$eni_share > .5]), stats::median(r$lrm_nb[r$eni_share <= .5])))
}))
csv_write(brand_tab, "tab_app_brand.csv")
PT$dbin <- cut(PT$dist, c(0, 0.5, 1, 2, 5, 10, Inf), labels = c("0--0.5", "0.5--1", "1--2", "2--5", "5--10", "$>$10"))
dist_tab <- do.call(rbind, lapply(1:3, function(s) {
  x <- PT[PT$setting == s & PT$block == "nb_both", ]
  data.frame(setting = s, dist = levels(PT$dbin), pumps = as.integer(table(x$dbin)),
             F = as.numeric(tapply(x$p_F < 0.05, x$dbin, function(v) round(100 * mean(v)))),
             hc1 = as.numeric(tapply(x$p_hc1 < 0.05, x$dbin, function(v) round(100 * mean(v)))))
}))
csv_write(dist_tab, "tab_app_distance.csv")
tex_write(c("\\begin{tabular}{lrrrrrrrrr}", "\\toprule",
  " & \\multicolumn{3}{c}{setting 1} & \\multicolumn{3}{c}{setting 2} & \\multicolumn{3}{c}{setting 3} \\\\",
  "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}",
  "distance to X (km) & pumps & $F$ & HC1 & pumps & $F$ & HC1 & pumps & $F$ & HC1 \\\\", "\\midrule",
  vapply(levels(PT$dbin), function(b) { z <- dist_tab[dist_tab$dist == b, ]
    paste(b, paste(sprintf("%d & %s & %s", z$pumps, ifelse(z$pumps > 0, z$F, "--"), ifelse(z$pumps > 0, z$hc1, "--")), collapse = " & "), sep = " & ") |>
      paste("\\\\") }, ""),
  "\\bottomrule", "\\end{tabular}"), "tab_app_distance.tex")

## the joint tests, all regimes at once, and the pooled model
if (!is.null(TW)) {
  j <- TW[TW$regime == "all" & TW$block == "nb_both", ]
  j <- merge(j, POOL[POOL$block == "nb_both", c("setting", "model", "p")], by = c("setting", "model"))
  j <- j[order(j$setting, factor(j$city, CITY_ORDER), j$fuel), ]
  csv_write(j, "tab_app_joint.csv")
  tex_write(c("\\begin{longtable}{llrrrrrr}",
    lt_head(paste("The Wald test of the eight lags of $X$ in all the regimes of a model at once: statistic and degrees",
                  "of freedom; $p$-value against the chi-square and against the bootstrap distribution of the statistic",
                  "centred at the estimate (its smallest value is $1/101$); $p$-value of the same test on the pooled",
                  "STEM model."), "tab:sc:joint",
            "setting & city & fuel & regimes & $W$ (df) & $p_{\\chi^2}$ & $p_{\\mathrm{boot}}$ & pooled $p$ \\\\"),
    sprintf("%d & %s & %s & %d & %.1f (%d) & %s & %s & %s \\\\", j$setting, en(j$city), j$fuel, j$df %/% 8, j$statistic, j$df,
            fmt_p(j$p_chisq), fmt_p(j$p_boot), fmt_p(j$p)), "\\bottomrule", "\\end{longtable}"), "tab_app_joint.tex")
}

## the pulses, event by event
pul_tab <- do.call(rbind, lapply(c("gasoline", "diesel"), function(f) do.call(rbind, lapply(unique(PUL$pulse), function(p) {
  x <- PUL[PUL$fuel == f & PUL$pulse == p, ]
  sp <- tapply(x$est, x$model, function(v) diff(range(v)))
  data.frame(fuel = f, pulse = p, regimes = nrow(x), sig = pct(x$sig), q10 = stats::quantile(x$est, .1),
             med = stats::median(x$est), q90 = stats::quantile(x$est, .9), spread_med = stats::median(sp[is.finite(sp)]))
}))))
csv_write(pul_tab, "tab_app_pulses.csv")
tex_write(c("\\begin{tabular}{llrrrrr}", "\\toprule",
  "fuel & pulse & $\\ne 0$ (\\%) & 10\\% & median & 90\\% & spread in a city \\\\", "\\midrule",
  sprintf("%s & %s & %d & %.1f & %.1f & %.1f & %.1f \\\\", pul_tab$fuel, gsub("_", " ", sub("^pulse_", "", pul_tab$pulse)), pul_tab$sig,
          pul_tab$q10, pul_tab$med, pul_tab$q90, pul_tab$spread_med), "\\bottomrule", "\\end{tabular}"), "tab_app_pulses.tex")

## residual diagnostics and the identification of the range
dg <- merge(DG, unique(REG[, c("setting", "model", "fuel")]))
diag_tab <- do.call(rbind, lapply(1:3, function(s) { x <- dg[dg$setting == s, ]; r <- REG[REG$setting == s, ]
  data.frame(setting = s, regimes = nrow(x), lb = stats::median(x$ljung_box_reject_pct), lb_sq = stats::median(x$ljung_box_sq_reject_pct),
             acf1 = stats::median(x$acf1_median), kurt = stats::median(x$kurtosis_median),
             theta_upper = sum(r$theta_bound == "upper"), theta_lower = sum(r$theta_bound == "lower")) }))
csv_write(diag_tab, "tab_app_diagnostics.csv")

## the local maxima (fuels-check-optima.R)
if (!is.null(OPT) && !is.null(OPTC) && !is.null(OPTW)) {
  d <- OPT[is.finite(OPT$diff), ]
  big <- OPTC[abs(OPTC$diff_loglik) > 10 & OPTC$block %in% c("nb", "nb_other"), ]
  w <- OPTW[OPTW$block == "nb_both", ]
  opt <- data.frame(regimes = nrow(d), within1 = sum(abs(d$diff) <= 1), over10 = sum(abs(d$diff) > 10),
                    over100 = sum(abs(d$diff) > 100), ls_higher10 = sum(d$diff > 10),
                    models_over10 = length(unique(paste(d$setting, d$model)[abs(d$diff) > 10])),
                    z_med = stats::median(abs(big$z_diff)), z_over2 = sum(abs(big$z_diff) > 2), coef = nrow(big),
                    wald_changes = sum((w$p_selected < 0.05) != (w$p_ls < 0.05)))
  csv_write(opt, "tab_app_optima.csv")
}

## the small samples of setting 2
small <- SM[SM$setting == 2, c("city", "fuel", "pumps", "K_max", "K", "phi", "B_valid")]
csv_write(small, "tab_app_small.csv")

## pump by pump against the regimes: the Dumitrescu-Hurlin test
dh <- DH[DH$block == "nb_both", c("setting", "city", "fuel", "pumps", "reject_F_pct", "reject_hc1_pct", "p_F", "p_hc1")]
csv_write(dh, "tab_app_dh.csv")

## the numbers quoted in the text
key_numbers <- list(
  regimes = nrow(REG), models = nrow(SM),
  lead_s3 = pct(REG$lead[REG$setting == 3]), lead_s1 = pct(REG$lead[REG$setting == 1]), lead_s2 = pct(REG$lead[REG$setting == 2]),
  pooled_misses = with(merge(city_tab, unique(REG[, c("setting", "fuel")])), sum(!pooled & lead_pumps >= 50)),
  own_other_lead = pct(REG$p_own_other < 0.05),
  lrm_pos_share = pct(REG$lrm_nb_lo > 0), lrm_neg_share = pct(REG$lrm_nb_hi < 0),
  phi_zero = sum(SM$phi == 0), phi_band = sum(SM$phi > 0 & SM$phi <= 0.2),
  ari_brand_median = stats::median(SM$ari_regimes_brands, na.rm = TRUE),
  boot_ari_min = min(REG$ari_boot), B_valid_min = min(REG$B_valid))
str(key_numbers)
saveRDS(list(REG = REG, POOL = POOL, PUL = PUL, PT = PT, DH = DH, DG = DG, TW = TW, city = city_tab,
             main = main_rows, brand = brand_tab, dist = dist_tab, pulses = pul_tab, key = key_numbers),
        file.path(OUT, "app_results.rds"))
message("done: outputs in ", OUT)
