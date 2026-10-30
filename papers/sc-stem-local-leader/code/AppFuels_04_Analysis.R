## ===========================================================================
## Latent Spatial Regimes in Spatio-Temporal State-Space Models: Evidence on
## Local Price Leadership in Italian Fuel Markets
##
## Author:  Paolo Maranzano (University of Milano-Bicocca; Fondazione Eni Enrico
##          Mattei)
## Script:  AppFuels_04_Analysis.R, the analysis of the application (Section 5
##          of the paper, Appendix C, Supplement C)
##
## What it does: reads the compact results of the 66 models (written by
## AppFuels_01_RunSetting1.R to AppFuels_03_RunSetting3.R; by default those of
## the run of the paper, shipped in results/application) and writes every
## table, figure and number of the application: those of Section 5 of the
## paper and of Supplement C, computed here, and those of Appendix C and of the
## description of the data (Sections 5.1, 5.2, C1 and C2), which come from the
## exploratory analysis of the weekly data, kept outside this package, and are
## stored in results/application/eda. Nothing is fitted: it runs in a minute.
##
## Outputs, in <out> (output/application by default), with the place of
## the paper each reproduces ("Section", "Table", "Figure" of the paper;
## "Section C.., Table C.., Figure C.." of Supplement C) and the section of
## this script that writes it:
##
##   tab_app_main.tex/.csv          Table 1 (Section 5.3)                     [2]
##   fig_app_leadership.png         Figure 3 (Section 5.3)                    [3]
##   fig_app_passthrough.png        Figure 5 (Section 5.3)                    [4]
##   fig_app_maps.png               Figure 4 (Section 5.3)                    [5]
##   fig_app_map_s1_RM_gasoline.png Supplement C, Figure C1 (Section C4), left [5]
##   fig_app_map_s3_NA_gasoline.png Supplement C, Figure C1 (Section C4), right [5]
##   tab_app_models.tex/.csv        Supplement C, Table C1 (Section C3)       [6]
##   tab_app_regimes_s1/2/3.tex     Supplement C, Tables C2 to C4 (Section C4) [6]
##   tab_app_joint.tex/.csv         Supplement C, Table C5 (Section C4)       [6]
##   tab_app_distance.tex/.csv      Supplement C, Table C6 (Section C5)       [6]
##   tab_app_pulses.tex/.csv        Supplement C, Table C7 (Section C7)       [6]
##   app_numbers.csv                every other number quoted in Sections 5.3
##                                  and 5.4 and in Sections C3 to C8, one row
##                                  each, with its place in the paper         [7]
##   app_*.csv                      the tables behind them (one row per regime,
##                                  per city, per pulse...)                [1-7]
##   eda/fig_app_centring.png       Figure C1 of Appendix C                   [8]
##   eda/fig_app_centring_stations.csv  the premiums quoted in its text      [8]
##   eda/eda_*.csv                  the numbers of Sections 5.1 and 5.2 of the
##                                  paper and of Sections C1 and C2 of
##                                  Supplement C (the map is in section 8)    [8]
## ===========================================================================
##
## Run from a terminal, in the folder of this script,
##
##     Rscript AppFuels_04_Analysis.R
##
## or open it in RStudio and press Source. Options (command line):
##   --results=<folder>  the compact results (default: ../results/application;
##                       ../output/runs/application after a run of the
##                       AppFuels_0<s>_RunSetting<s>.R scripts)
##   --eda=<folder>      the stored outputs of the exploratory analysis
##                       (default: ../results/application/eda)
##   --out=<folder>      the outputs (default: ../output/application)
##
## Requirements: R (>= 4.2) with data.table; Stem is not needed.
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
source(file.path(HERE, "AppFuels_00_Setup.R"))

OPT <- fu_options(list(
  results = file.path(HERE, "..", "results", "application"),
  eda     = file.path(HERE, "..", "results", "application", "eda"),
  out     = file.path(HERE, "..", "output", "application")
))
OUT <- OPT$out; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
options(width = 160)

## --- helpers -----------------------------------------------------------------
en <- function(cc) unname(FU_CITY_EN[as.character(cc)])      # English name of a city
pct <- function(x) round(100 * mean(x, na.rm = TRUE))
fmt_p <- function(p) ifelse(is.na(p), "--", ifelse(p < 0.001, "$<$0.001", sprintf("%.3f", p)))
csv_write <- function(d, name) fu_write(d, file.path(OUT, name))
## a .tex file: minus signs in math mode, and no negative zero
tex_write <- function(lines, name) {
  lines <- gsub("(?<![0-9A-Za-z.$])-(0\\.0+)(?![0-9])", "\\1", lines, perl = TRUE)
  lines <- gsub("(?<![0-9A-Za-z.$-])-(?=[0-9])", "$-$", lines, perl = TRUE)
  writeLines(lines, file.path(OUT, name)); message("wrote ", file.path(OUT, name))
}
## the caption and the header of a longtable, repeated on every page
lt_head <- function(cap, lab, header)
  c(sprintf("\\caption{%s}\\label{%s}\\\\", cap, lab), "\\toprule", header, "\\midrule", "\\endfirsthead",
    "\\toprule", header, "\\midrule", "\\endhead")
fig_open <- function(name, width, height)
  grDevices::png(file.path(OUT, name), width = width, height = height, units = "in", res = 300,
                 type = if (isTRUE(capabilities("cairo"))) "cairo" else "windows")
fig_close <- function(name) { invisible(grDevices::dev.off()); message("wrote ", file.path(OUT, name)) }
PAL3 <- c("#1f6f8b", "#e0a458", "#a8516e")
EVENTS <- c(pulse_excise_cut = "excise cut\n22 Mar 2022", pulse_cut_reduced = "cut reduced\n1 Dec 2022",
            pulse_cut_ends = "cut ends\n1 Jan 2023", pulse_realignment = "realignment\n15 May 2025",
            pulse_equal_excise = "equal excise\n1 Jan 2026")
## the numbers quoted in the text: one row each, with its place in the paper
NUM <- list()
num <- function(place, quantity, value)
  NUM[[length(NUM) + 1]] <<- data.frame(place = place, quantity = quantity,
                                        value = paste(format(value, digits = 3), collapse = "; "))





## ===========================================================================
## 1. The results of every model
## ===========================================================================
## A compressed file is decompressed by R and handed to fread as text: on some
## systems fread decompresses through a temporary file, which is far slower.
rd <- function(s, m, f, ...) {
  p <- file.path(OPT$results, sprintf("setting%d", s), m, f)
  if (!file.exists(p)) return(NULL)
  if (grepl("[.]gz$", p)) as.data.frame(data.table::fread(text = readLines(gzfile(p)), na.strings = c("", "NA"), ...))
  else as.data.frame(data.table::fread(p, na.strings = c("", "NA"), ...))
}
SM <- do.call(rbind, lapply(1:3, function(s)
  utils::read.csv(file.path(OPT$results, sprintf("setting%d", s), "models.csv"), stringsAsFactors = FALSE)))
## the city from the name of the model (the code of Naples, "NA", would be read
## as a missing value), as its Italian name, ordered from north to south
SM$city <- factor(FU_CITY[sub("_.*$", "", SM$model)], FU_CITY_ORDER)
## the models in the order of the runs (setting, city as in FU_CITY, gasoline
## first), which the random jitter of Figure 5 follows
SM <- SM[order(SM$setting, match(sub("_.*$", "", SM$model), names(FU_CITY)), match(SM$fuel, FU_FUEL)), ]
## The benchmark of the adjusted Rand index between the regimes and the brands
## (column ari_regimes_brands, written by the run scripts): the same index between
## the regimes and a purely spatial partition of the stations Y with as many
## groups (k-means on their coordinates; see fu_ari_space). The regimes are
## groups of brands rather than districts if the first index exceeds the
## second. The seed of a model is the one of its fit.
SM$ari_regimes_space <- vapply(seq_len(nrow(SM)), function(i) {
  cc <- sub("_.*$", "", SM$model[i])
  g <- rd(SM$setting[i], SM$model[i], "regimes.csv")
  fu_ari_space(g$regime, g$y_lon, g$y_lat, seed = fu_seed(cc, names(FU_FUEL)[match(SM$fuel[i], FU_FUEL)]))
}, numeric(1))

## One row per regime of the selected fits: size, brands, distance to X, the
## p-values of the four bootstrap Wald tests, and the long-run pass-through
## from X (Equation 10 of the paper): the sum of the coefficients of the lags of
## X over one minus the sum of the own lags, with the percentile interval of
## its aligned bootstrap draws, on the same fuel and on both fuels.
LAG <- function(b) paste0(b, "_l", 1:4)
REG <- do.call(rbind, lapply(seq_len(nrow(SM)), function(i) {
  s <- SM$setting[i]; m <- SM$model[i]; K <- SM$K[i]
  g <- rd(s, m, "regimes.csv", colClasses = list(character = c("y_site", "x_site")))
  w <- rd(s, m, "wald.csv"); w <- w[w$fit == "selected" & w$regime != "all", ]
  pw <- function(b) w$p_chisq[w$block == b][order(as.integer(w$regime[w$block == b]))]
  e <- rd(s, m, "estimates.csv"); e <- e[e$fit == "selected", ]
  dr <- rd(s, m, "draws.csv.gz")
  cov <- unique(dr$covariate)
  A <- array(NA_real_, c(max(dr$draw), K, length(cov)), dimnames = list(NULL, NULL, cov))
  A[cbind(dr$draw, dr$regime, match(dr$covariate, cov))] <- dr$value
  lrm <- t(vapply(seq_len(K), function(k) {
    est <- function(v) sum(e$estimate[e$cluster == k & e$covariate %in% v])
    drw <- function(v) rowSums(matrix(A[, k, v], ncol = length(v)))
    l1 <- drw(LAG("nb")) / (1 - drw(LAG("own"))); l2 <- (drw(LAG("nb")) + drw(LAG("nb_other"))) / (1 - drw(LAG("own")))
    c(est(LAG("nb")) / (1 - est(LAG("own"))), stats::quantile(l1, c(.025, .975), na.rm = TRUE),
      (est(LAG("nb")) + est(LAG("nb_other"))) / (1 - est(LAG("own"))), stats::quantile(l2, c(.025, .975), na.rm = TRUE),
      est(LAG("own")))
  }, numeric(7)))
  br <- table(factor(g$regime, seq_len(K)), g$y_brand)
  data.frame(setting = s, model = m, city = SM$city[i], fuel = SM$fuel[i], K = K, phi = SM$phi[i],
             regime = seq_len(K), pumps = as.integer(table(factor(g$regime, seq_len(K)))),
             p_nb_both = pw("nb_both"), p_nb = pw("nb"), p_nb_other = pw("nb_other"), p_own_other = pw("own_other"),
             lrm_nb = lrm[, 1], lrm_nb_lo = lrm[, 2], lrm_nb_hi = lrm[, 3],
             lrm_both = lrm[, 4], lrm_both_lo = lrm[, 5], lrm_both_hi = lrm[, 6], sum_own = lrm[, 7],
             dist_median = as.numeric(tapply(g$x_dist_km, factor(g$regime, seq_len(K)), stats::median)),
             top_brand = apply(br, 1, function(r) names(r)[which.max(r)]),
             top_share = as.numeric(apply(br, 1, max) / rowSums(br)),
             eni_share = if ("Agip Eni" %in% colnames(br)) as.numeric(br[, "Agip Eni"] / rowSums(br)) else 0,
             theta_bound = as.character(tapply(g$theta_bound, factor(g$regime, seq_len(K)), `[`, 1)),
             ari_boot = SM$ari_draws_median[i], B_valid = SM$B_valid[i], stringsAsFactors = FALSE)
}))
REG$lead <- REG$p_nb_both < 0.05                                # X leads Y in the regime
REG$area <- unname(FU_AREA[as.character(REG$city)])

## the test of leadership on the pooled model (the selected fit when K = 1)
POOL <- do.call(rbind, lapply(seq_len(nrow(SM)), function(i) {
  s <- SM$setting[i]; m <- SM$model[i]
  w <- rd(s, m, "wald.csv"); w <- w[w$fit == (if (SM$K[i] == 1) "selected" else "pooled") & w$regime != "all", ]
  data.frame(setting = s, model = m, city = SM$city[i], fuel = SM$fuel[i], block = w$block, p = w$p_chisq)
}))
## the tests in all the regimes of a model at once
TW <- do.call(rbind, lapply(seq_len(nrow(SM)), function(i) {
  w <- rd(SM$setting[i], SM$model[i], "wald.csv"); w <- w[w$fit == "selected" & w$regime == "all", ]
  cbind(city = SM$city[i], fuel = SM$fuel[i], w)
}))
## the pulses of the excise duties, one row per regime and pulse
PUL <- do.call(rbind, lapply(seq_len(nrow(SM)), function(i) {
  s <- SM$setting[i]; m <- SM$model[i]
  e <- rd(s, m, "estimates.csv"); e <- e[e$fit == "selected" & startsWith(e$covariate, "pulse_"), ]
  data.frame(setting = s, model = m, city = SM$city[i], fuel = SM$fuel[i], regime = e$cluster, pulse = e$covariate,
             est = e$estimate, se = e$se, lo = e$perc_lower_95_CI, hi = e$perc_upper_95_CI)
}))
PUL$sig <- PUL$lo > 0 | PUL$hi < 0                              # the interval excludes zero
PUL$area <- unname(FU_AREA[as.character(PUL$city)])
## the station-by-station tests, with the distance to X and the regime
PT <- do.call(rbind, lapply(seq_len(nrow(SM)), function(i) {
  s <- SM$setting[i]; m <- SM$model[i]
  p <- rd(s, m, "pump_tests.csv", colClasses = list(character = "y_site"))
  g <- rd(s, m, "regimes.csv", colClasses = list(character = c("y_site", "x_site")))
  data.frame(setting = s, model = m, city = SM$city[i], fuel = SM$fuel[i], block = p$block, p_F = p$p_F,
             p_hc1 = p$p_hc1, dist = g$x_dist_km[match(p$y_site, g$y_site)])
}))
DH <- do.call(rbind, lapply(seq_len(nrow(SM)), function(i)
  cbind(city = SM$city[i], fuel = SM$fuel[i], rd(SM$setting[i], SM$model[i], "dumitrescu_hurlin.csv")[, -(1:2)],
        setting = SM$setting[i], model = SM$model[i])))
DG <- do.call(rbind, lapply(seq_len(nrow(SM)), function(i) rd(SM$setting[i], SM$model[i], "diagnostics.csv")))
OPTR <- do.call(rbind, lapply(seq_len(nrow(SM)), function(i) rd(SM$setting[i], SM$model[i], "optima.csv")))
OPTC <- do.call(rbind, lapply(seq_len(nrow(SM)), function(i) rd(SM$setting[i], SM$model[i], "optima_coefficients.csv")))
csv_write(REG, "app_regimes.csv"); csv_write(PUL, "app_pulses.csv")





## ===========================================================================
## 2. Table 1 (Section 5.3): one row per setting and fuel
## ===========================================================================
## Regimes and their range over the cities; the share of the regimes, and of
## the stations Y in them, where X leads; the cities where the pooled model
## finds leadership; the median and quartiles of the long-run pass-through
## over the regimes and the share whose interval lies above zero; the median
## over the cities of the share of stations whose HC1 (and F) test rejects.
main_rows <- do.call(rbind, lapply(1:3, function(s) do.call(rbind, lapply(c("gasoline", "diesel"), function(f) {
  r <- REG[REG$setting == s & REG$fuel == f, ]; sm <- SM[SM$setting == s & SM$fuel == f, ]
  po <- POOL[POOL$setting == s & POOL$fuel == f & POOL$block == "nb_both", ]
  dh <- DH[DH$setting == s & DH$fuel == f & DH$block == "nb_both", ]
  data.frame(setting = s, fuel = f, pumps = sum(sm$pumps), regimes = nrow(r),
             K = paste(range(sm$K), collapse = "--"),
             lead_regimes = pct(r$lead), lead_pumps = round(100 * sum(r$pumps[r$lead]) / sum(r$pumps)),
             pooled_lead = sum(po$p < 0.05), cities = nrow(po),
             lrm_med = stats::median(r$lrm_nb), lrm_q1 = stats::quantile(r$lrm_nb, .25),
             lrm_q3 = stats::quantile(r$lrm_nb, .75), lrm_pos = pct(r$lrm_nb_lo > 0),
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
print(main_rows, digits = 3, row.names = FALSE)
num("Section 5.3, Q1, per pump", "F test: share of stations rejecting, median over the cities, range over settings and fuels (%)",
    range(main_rows$F))





## ===========================================================================
## 3. Figure 3 (Section 5.3): leadership city by city
## ===========================================================================
## For every setting, city and fuel: the share of the stations Y in the regimes
## where X leads, and whether the pooled model finds leadership.
city_tab <- do.call(rbind, lapply(1:3, function(s) do.call(rbind, lapply(c("gasoline", "diesel"), function(f)
  do.call(rbind, lapply(FU_CITY_ORDER, function(cc) {
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
  yy <- match(x$city, rev(FU_CITY_ORDER))
  plot(NA, xlim = c(0, 100), ylim = c(0.5, length(FU_CITY_ORDER) + 0.5), yaxt = "n", xlab = "% of the pumps Y",
       ylab = "", main = sprintf("Setting %d\n%s", s, FU_SETTINGS$label[s]), cex.main = 1.05, font.main = 1)
  graphics::abline(h = seq_along(FU_CITY_ORDER), col = "grey92")
  if (s == 1) graphics::axis(2, at = seq_along(FU_CITY_ORDER), labels = en(rev(FU_CITY_ORDER)), las = 1, cex.axis = 1.05)
  off <- ifelse(x$fuel == "gasoline", 0.15, -0.15)
  graphics::points(x$lead_pumps, yy + off, pch = ifelse(x$pooled, 19, 1), cex = 1.5,
                   col = ifelse(x$fuel == "gasoline", PAL3[1], PAL3[3]), lwd = 1.6)
}
graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
graphics::legend("bottom", horiz = TRUE, bty = "n", cex = 1.05, pt.cex = 1.4, pt.lwd = 1.6,
                 legend = c("gasoline", "diesel", "pooled model: X leads", "pooled model: no"),
                 pch = c(19, 19, 19, 1), col = c(PAL3[1], PAL3[3], "grey30", "grey30"))
graphics::mtext("Share of the pumps Y in the regimes where the neighbor X leads (bootstrap Wald test, 5%)",
                side = 3, line = -1.2, cex = 0.95)
fig_close("fig_app_leadership.png")





## ===========================================================================
## 4. Figure 5 (Section 5.3): the pass-through from X and from the excise
## ===========================================================================
## (a) the long-run pass-through from X by setting and fuel, one value per
## regime; (b) the coefficient of the pulse of the week of each change of the
## excise duties, one value per regime.
fig_open("fig_app_passthrough.png", 9.6, 4.8)
graphics::layout(matrix(1:2, 1), widths = c(1, 1.4))
graphics::par(mar = c(5.2, 4.2, 2.6, 0.8), mgp = c(2.4, 0.7, 0))
grp <- interaction(REG$fuel, REG$setting, sep = " ", lex.order = TRUE)
grp <- factor(grp, c(outer(c("gasoline", "diesel"), 1:3, paste)))
graphics::boxplot(REG$lrm_nb ~ grp, outline = FALSE, col = rep(c(PAL3[1], PAL3[3]), 3), border = "grey25",
                  xaxt = "n", ylab = "cents per cent", xlab = "", ylim = range(REG$lrm_nb, na.rm = TRUE),
                  main = "(a) Long-run pass-through from X", cex.main = 1, font.main = 1)
set.seed(1)                                                     # the horizontal jitter of the points
graphics::points(as.integer(grp) + stats::runif(length(grp), -0.18, 0.18), REG$lrm_nb, pch = 16, cex = 0.5,
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
                  xaxt = "n", ylab = "cents per liter", xlab = "", ylim = range(bs$stats),
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
graphics::legend("bottom", horiz = TRUE, bty = "n", cex = 0.9, fill = c(PAL3[1], PAL3[3], NA, NA),
                 border = c("grey25", "grey25", NA, NA), pch = c(NA, NA, 16, 16), col = c(NA, NA, "grey10", "grey65"),
                 legend = c("gasoline", "diesel", "regime: interval excludes 0", "regime: interval includes 0"))
fig_close("fig_app_passthrough.png")





## ===========================================================================
## 5. Figure 4 (Section 5.3) and Figure C1 of Supplement C: maps of regimes
## ===========================================================================
## Every station Y at its coordinates, colored by its regime; filled marks
## for the regimes where X leads. The legend gives the size of every regime and
## its most frequent brand with its share.
PALM <- c("#1f6f8b", "#e0a458", "#a8516e", "#5b8c5a", "#6c5b7b")
draw_map <- function(s, m, title) {
  g <- rd(s, m, "regimes.csv", colClasses = list(character = c("y_site", "x_site")))
  r <- REG[REG$setting == s & REG$model == m, ]
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
## Figure C1 of Supplement C: Rome, gasoline, setting 1, and Naples, gasoline,
## setting 3 (one file each)
for (cs in list(c(1, "RM_gasoline"), c(3, "NA_gasoline"))) {
  s <- as.integer(cs[1]); m <- cs[2]; nm <- sprintf("fig_app_map_s%d_%s.png", s, m)
  fig_open(nm, 5, 5.6)
  draw_map(s, m, sprintf("%s, %s, setting %d", en(FU_CITY[[sub("_.*$", "", m)]]), sub("^.*_", "", m), s))
  fig_close(nm)
}
## Figure 4: Rome, gasoline, and Milan, diesel, setting 3
fig_open("fig_app_maps.png", 9, 5.4)
graphics::par(mfrow = c(1, 2))
draw_map(3, "RM_gasoline", "(a) Rome, gasoline, setting 3")
draw_map(3, "MI_diesel", "(b) Milan, diesel, setting 3")
fig_close("fig_app_maps.png")
r <- REG[REG$setting == 3 & REG$model == "RM_gasoline", ]
num("Section 5.3, Figure 4 caption", "Rome, gasoline, setting 3: regimes; regimes where X leads", c(nrow(r), sum(r$lead)))
r <- REG[REG$setting == 3 & REG$model == "MI_diesel", ]
num("Section 5.3, Figure 4 caption", "Milan, diesel, setting 3: size, Eni share of regime 3", c(r$pumps[3], r$eni_share[3]))





## ===========================================================================
## 6. Supplement C: the tables
## ===========================================================================

## --- Table C1 (Section C3): the selected models -----------------------------
mod <- SM[order(SM$setting, SM$city, SM$fuel), ]
csv_write(mod, "tab_app_models.csv")
tex_write(c("\\begin{longtable}{llrrrlrrrr}",
  lt_head(paste("The selected models. $K_{\\max}$, the largest number of regimes of the grid; $\\hat K$, $\\hat\\phi$",
                "and the sizes of the regimes of the selected fit; ARI brands, the adjusted Rand index between the regimes",
                "and the brands of the stations $Y$ (not defined in setting 2, where they are all independent); ARI space,",
                "the same index between the regimes and the $k$-means partition of the coordinates of the stations $Y$",
                "into $\\hat K$ groups; usable draws, of $B = 100$; ARI draws, the median adjusted Rand index between the",
                "partitions of the bootstrap draws and the estimated one."), "tab:sc:models",
          "setting & city & fuel & pumps & $K_{\\max}$ & $\\hat K$, $\\hat\\phi$, sizes & ARI brands & ARI space & usable draws & ARI draws \\\\"),
  sprintf("%d & %s & %s & %d & %d & %d, %s, %s & %s & %s & %d & %.2f \\\\", mod$setting, en(mod$city), mod$fuel, mod$pumps,
          mod$K_max, mod$K, format(mod$phi), mod$sizes,
          ifelse(is.na(mod$ari_regimes_brands), "--", sprintf("%.2f", mod$ari_regimes_brands)),
          ifelse(is.na(mod$ari_regimes_space), "--", sprintf("%.2f", mod$ari_regimes_space)),
          mod$B_valid, mod$ari_draws_median),
  "\\bottomrule", "\\end{longtable}"), "tab_app_models.tex")

## --- Tables C2 to C4 (Section C4): the regimes, setting by setting ----------
for (s in 1:3) {
  r <- REG[REG$setting == s, ]; r <- r[order(r$city, r$fuel, r$regime), ]
  tex_write(c("\\begin{longtable}{llrrlrrrrrr}",
    lt_head(paste(sprintf("The regimes of setting %d.", s),
                  "Brand: the most frequent brand of the stations $Y$ of the regime, with its share; X: median distance",
                  "to the neighbor, in kilometers; $p$-values of the bootstrap Wald test of the eight lags of $X$, of its",
                  "four lags on the same fuel, of its four lags on the other fuel, and of the four lags of the other fuel",
                  "of $Y$; long-run pass-through from $X$ on the same fuel, with its bootstrap 95\\% interval."),
            sprintf("tab:sc:reg%d", s),
            c("city & fuel & regime & pumps & brand (\\%) & X (km) & \\multicolumn{4}{c}{bootstrap Wald, $p$-value} & long-run pass-through \\\\",
              "\\cmidrule(lr){7-10}", " & & & & & & X & X same & X other & Y other & from X, same fuel \\\\")),
    sprintf("%s & %s & %d & %d & %s %.0f & %.1f & %s & %s & %s & %s & %.2f (%.2f, %.2f) \\\\", en(r$city), r$fuel,
            r$regime, r$pumps, sub("Agip Eni", "Eni", gsub("Pompe Bianche", "indep.", r$top_brand)), 100 * r$top_share,
            r$dist_median, fmt_p(r$p_nb_both), fmt_p(r$p_nb), fmt_p(r$p_nb_other), fmt_p(r$p_own_other),
            r$lrm_nb, r$lrm_nb_lo, r$lrm_nb_hi),
    "\\bottomrule", "\\end{longtable}"), sprintf("tab_app_regimes_s%d.tex", s))
}

## --- Table C5 (Section C4): all the regimes at once, and the pooled model --
j <- TW[TW$block == "nb_both", ]
j <- merge(j, POOL[POOL$block == "nb_both", c("setting", "model", "p")], by = c("setting", "model"))
j <- j[order(j$setting, factor(j$city, FU_CITY_ORDER), j$fuel), ]
csv_write(j, "tab_app_joint.csv")
tex_write(c("\\begin{longtable}{llrrrrrr}",
  lt_head(paste("The Wald test of the eight lags of $X$ in all the regimes of a model at once: statistic and degrees",
                "of freedom; $p$-value against the chi-square and against the bootstrap distribution of the statistic",
                "centered at the estimate (its smallest value is $1/101$); $p$-value of the same test on the pooled",
                "STEM model."), "tab:sc:joint",
          "setting & city & fuel & regimes & $W$ (df) & $p_{\\chi^2}$ & $p_{\\mathrm{boot}}$ & pooled $p$ \\\\"),
  sprintf("%d & %s & %s & %d & %.1f (%d) & %s & %s & %s \\\\", j$setting, en(j$city), j$fuel, j$df %/% 8, j$statistic,
          j$df, fmt_p(j$p_chisq), fmt_p(j$p_boot), fmt_p(j$p)), "\\bottomrule", "\\end{longtable}"), "tab_app_joint.tex")

## --- Table C6 (Section C5): the station-by-station tests by distance to X --
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
    paste(paste(b, paste(sprintf("%d & %s & %s", z$pumps, ifelse(z$pumps > 0, z$F, "--"),
                                 ifelse(z$pumps > 0, z$hc1, "--")), collapse = " & "), sep = " & "), "\\\\") }, ""),
  "\\bottomrule", "\\end{tabular}"), "tab_app_distance.tex")

## --- Table C7 (Section C7): the pulses of the excise duties ----------------
## Over the regimes: the share whose interval excludes zero, the 10th, 50th and
## 90th percentiles of the coefficients, and the median over the models of the
## spread (largest minus smallest coefficient) between the regimes of a city.
pul_tab <- do.call(rbind, lapply(c("gasoline", "diesel"), function(f) do.call(rbind, lapply(unique(PUL$pulse), function(p) {
  x <- PUL[PUL$fuel == f & PUL$pulse == p, ]
  sp <- tapply(x$est, x$model, function(v) diff(range(v)))
  data.frame(fuel = f, pulse = p, regimes = nrow(x), sig = pct(x$sig), q10 = stats::quantile(x$est, .1),
             med = stats::median(x$est), q90 = stats::quantile(x$est, .9), spread_med = stats::median(sp[is.finite(sp)]))
}))))
csv_write(pul_tab, "tab_app_pulses.csv")
tex_write(c("\\begin{tabular}{llrrrrr}", "\\toprule",
  "fuel & pulse & $\\ne 0$ (\\%) & 10\\% & median & 90\\% & spread in a city \\\\", "\\midrule",
  sprintf("%s & %s & %d & %.1f & %.1f & %.1f & %.1f \\\\", pul_tab$fuel, gsub("_", " ", sub("^pulse_", "", pul_tab$pulse)),
          pul_tab$sig, pul_tab$q10, pul_tab$med, pul_tab$q90, pul_tab$spread_med), "\\bottomrule", "\\end{tabular}"),
  "tab_app_pulses.tex")





## ===========================================================================
## 7. The numbers quoted in the text
## ===========================================================================
SL <- function(s) sprintf("setting %d", s)

## --- Section 5.3, Q1: local leadership --------------------------------------
for (s in c(3, 1, 2)) {
  x <- main_rows[main_rows$setting == s, ]
  num("Section 5.3, Q1", paste0(SL(s), ": regimes where X leads, gasoline; diesel (%)"), x$lead_regimes)
  num("Section 5.3, Q1", paste0(SL(s), ": stations Y in those regimes, gasoline; diesel (%)"), x$lead_pumps)
  num("Section 5.3, Q1", paste0(SL(s), ": median long-run pass-through, gasoline; diesel"), round(x$lrm_med, 2))
  num("Section 5.3, Q1", paste0(SL(s), ": regimes whose pass-through interval lies above zero, gasoline; diesel (%)"), x$lrm_pos)
}
num("Section 5.3, Q1", "setting 3: cities with a regime where X leads, gasoline; diesel",
    sapply(c("gasoline", "diesel"), function(f) length(unique(REG$city[REG$setting == 3 & REG$fuel == f & REG$lead]))))
x <- city_tab[city_tab$setting == 1 & city_tab$city == "Bari", ]
num("Section 5.3, Q1", "setting 1, Bari: stations in the regimes where X leads, gasoline; diesel (%)", round(x$lead_pumps))
num("Section 5.3, Q1", "setting 1, Bari: median distance to X over the regimes, range (km)",
    range(REG$dist_median[REG$setting == 1 & REG$city == "Bari"]))
x <- city_tab[city_tab$setting == 2, ]
num("Section 5.3, Q1", "setting 2: cities where all the stations Y follow X (both fuels)",
    en(names(which(tapply(x$lead_pumps, x$city, min) == 100))))
num("Section 5.3, Q1", "setting 2: cities where at most half of the stations Y follow X (both fuels, %)",
    en(names(which(tapply(x$lead_pumps, x$city, max) <= 50))))
brand_tab <- do.call(rbind, lapply(c(1, 3), function(s) {
  r <- REG[REG$setting == s, ]; e <- r$eni_share > 0.5
  data.frame(setting = s, regimes = c("Eni > 50%", "other"), number = c(sum(e), sum(!e)),
             lead_pct = c(pct(r$lead[e]), pct(r$lead[!e])),
             lrm_median = round(c(stats::median(r$lrm_nb[e]), stats::median(r$lrm_nb[!e])), 2))
}))
csv_write(brand_tab, "tab_app_brand.csv")
num("Section 5.3, Q1, and Section C4 (table Brands)", "regimes with Eni > 50%: number, X leads (%), median pass-through; settings 1 and 3",
    paste(brand_tab$setting, brand_tab$regimes, brand_tab$number, brand_tab$lead_pct, brand_tab$lrm_median))
x <- city_tab[!city_tab$pooled & city_tab$lead_pumps >= 50, ]
num("Section 5.3, Q1", "the pooled model misses leadership where most stations follow X (setting, city, fuel)",
    paste(x$setting, en(x$city), x$fuel))
x <- city_tab[city_tab$pooled & city_tab$lead_pumps < 50, ]
num("Section 5.3, Q1", "the pooled model finds leadership in regimes covering less than half of the stations: models",
    nrow(x))
num("Section 5.3, Q1", "setting 3: models where the pooled model finds leadership, of 22",
    sum(city_tab$pooled[city_tab$setting == 3]))
num("Section 5.3, Q1, per pump", "HC1 test: share of the stations rejecting, median over the cities, range (%)",
    range(main_rows$hc1))
x <- dist_tab[dist_tab$setting == 3, ]
num("Section 5.3, Q1, per pump, and Section C5", "setting 3: HC1 rejections within 0.5 km; 2 to 10 km (%)",
    c(x$hc1[x$dist == "0--0.5"], round(100 * mean(PT$p_hc1[PT$setting == 3 & PT$block == "nb_both" &
                                                       PT$dist > 2 & PT$dist <= 10] < 0.05))))

## --- Section 5.3, Q2: the regimes -------------------------------------------
s13 <- SM[SM$setting != 2, ]
num("Section 5.3, Q2, and Section C3", "settings 1 and 3: range of K; models at the cap of five, of 44",
    c(range(s13$K), sum(s13$K == 5)))
num("Section 5.3, Q2", "setting 2: range of K", range(SM$K[SM$setting == 2]))
num("Section 5.3, Q2, and Section C3", "models with phi = 0; with phi <= 0.2 (of 66)", c(sum(SM$phi == 0), sum(SM$phi <= 0.2)))
num("Section 5.3, Q2, and Section C3", "ARI between regimes and brands, settings 1 and 3: min, max, median",
    round(c(range(s13$ari_regimes_brands), stats::median(s13$ari_regimes_brands)), 3))
num("Section 5.3, Q2, and Section C3", "ARI between regimes and the k-means partition of the coordinates, settings 1 and 3: min, max, median",
    round(c(range(s13$ari_regimes_space), stats::median(s13$ari_regimes_space)), 3))
num("Section 5.3, Q2, and Section C3", "settings 1 and 3: models where the ARI with the brands exceeds the ARI with space, of 44",
    sum(s13$ari_regimes_brands > s13$ari_regimes_space))
eni1 <- REG[REG$setting != 2 & REG$top_share == 1 & REG$top_brand == "Agip Eni", ]
num("Section 5.3, Q2, and Section C3", "settings 1 and 3: models with a regime of Eni stations only; their sizes",
    c(length(unique(paste(eni1$setting, eni1$model))), range(eni1$pumps)))
q8 <- REG[REG$setting != 2 & REG$top_share == 1 & REG$top_brand == "Q8", ]
num("Section C3", "settings 1 and 3: models with a regime of Q8 stations only", length(unique(paste(q8$setting, q8$model))))

## --- Section 5.3, Q3: the excise duties -------------------------------------
x <- pul_tab[pul_tab$pulse == "pulse_excise_cut", ]
num("Section 5.3, Q3", "excise cut: regimes whose interval excludes 0, gasoline; diesel (%)", x$sig)
num("Section 5.3, Q3", "excise cut: 10th percentile, gasoline; diesel (cents)", round(x$q10, 1))
num("Section 5.3, Q3", "excise cut: 90th percentile, gasoline; diesel (cents)", round(x$q90, 1))
num("Section 5.3, Q3", "excise cut: median spread between the regimes of a city, gasoline; diesel (cents)", round(x$spread_med, 1))
num("Section 5.3, Q3", "realignment 2025: regimes whose interval excludes 0, both fuels (%)",
    pct(PUL$sig[PUL$pulse == "pulse_realignment"]))

## --- Section 5.4: the geography of leadership -------------------------------
## Pooling the two fuels, the share of the stations Y in the regimes where X
## leads, by area; the median distance to X, the median pass-through and the
## spread of the excise-cut pulse by area; the cities with a regime of Eni
## stations only.
a <- stats::aggregate(cbind(pl = pumps * lead, pumps) ~ setting + area, REG, sum)
a$pct <- round(100 * a$pl / a$pumps)
csv_write(a, "app_leadership_by_area.csv")
for (s in c(1, 3)) num("Section 5.4", paste0(SL(s), ": stations Y in the regimes where X leads, North; Center; South and islands (%)"),
                       a$pct[a$setting == s][match(c("North", "Center", "South and islands"), a$area[a$setting == s])])
for (s in c(1, 3)) { x <- city_tab[city_tab$setting == s, ]; x <- x[order(x$lead_pumps), ][1:5, ]
  num("Section 5.4", paste0(SL(s), ": the five weakest city-fuels (stations following X, %)"),
      paste(en(x$city), x$fuel, round(x$lead_pumps))) }
dd <- stats::aggregate(dist_median ~ setting + area, REG, function(v) round(stats::median(v), 2))
for (s in c(1, 3)) num("Section 5.4", paste0(SL(s), ": median distance to X by area, North; Center; South (km)"),
                       dd$dist_median[dd$setting == s][match(c("North", "Center", "South and islands"), dd$area[dd$setting == s])])
ll <- stats::aggregate(lrm_nb ~ setting + area, REG, function(v) round(stats::median(v), 2))
num("Section 5.4", "setting 3: median pass-through over the regimes by area, North; Center; South",
    ll$lrm_nb[ll$setting == 3][match(c("North", "Center", "South and islands"), ll$area[ll$setting == 3])])
e <- PUL[PUL$pulse == "pulse_excise_cut" & PUL$setting != 2, ]
sp <- stats::aggregate(est ~ setting + model + city + area, e, function(v) diff(range(v)))
spa <- stats::aggregate(est ~ area, sp, function(v) round(stats::median(v)))
num("Section 5.4", "excise cut, settings 1 and 3: median spread between the regimes of a city, North; Center; South (cents)",
    spa$est[match(c("North", "Center", "South and islands"), spa$area)])
spc <- stats::aggregate(est ~ city, sp, function(v) round(stats::median(v)))
spc <- spc[order(-spc$est), ]
csv_write(spc, "app_excise_spread_by_city.csv")
num("Section 5.4", "excise cut, settings 1 and 3: median spread by city, largest to smallest (cents)",
    paste(en(spc$city), spc$est))
num("Section 5.4", "cities with a regime of Eni stations only (settings 1 and 3)",
    en(intersect(FU_CITY_ORDER, as.character(unique(eni1$city)))))
if (file.exists(file.path(OPT$eda, "eda_unit_roots_and_stickiness.csv"))) {
  st <- utils::read.csv(file.path(OPT$eda, "eda_unit_roots_and_stickiness.csv"), stringsAsFactors = FALSE)
  st <- st[order(-st$unchanged_weeks_median_pct), ]
  num("Section 5.4", "the stickiest cities: median share of Sunday prices unchanged from the previous week (%)",
      paste(st$city, st$fuel, st$unchanged_weeks_median_pct)[1:4])
}

## --- Section C3: the selected models ----------------------------------------
s2 <- SM[SM$setting == 2, ]
num("Section C3 and Section C8", "setting 2: models at the cap floor(n/5), of 22", sum(s2$K == s2$K_max))
num("Section C3", "models with phi in the band [0.025, 0.2]; with phi in {0.5, 1}",
    c(sum(SM$phi >= FU_BAND[1] & SM$phi <= FU_BAND[2]), sum(SM$phi >= 0.5)))
num("Section C3 and Section C8", "usable bootstrap draws: min, max", range(SM$B_valid))
x <- SM[SM$ari_draws_median < 1, ]
num("Section C3 and Section C8", "models whose partition is not reproduced in the median draw (ARI)",
    paste(SL(x$setting), en(x$city), x$fuel, x$ari_draws_median))

## --- Section C4: leadership by regime ---------------------------------------
for (b in c("nb", "nb_other"))
  num("Section C4", paste0("the four lags of X on the ", if (b == "nb") "same" else "other", " fuel are significant: % of the regimes, settings 1; 2; 3"),
      sapply(1:3, function(s) pct(REG[[paste0("p_", b)]][REG$setting == s] < 0.05)))
num("Section C4", "the other fuel of Y predicts its price: % of the regimes", pct(REG$p_own_other < 0.05))
jt <- j
num("Section C4", "all the regimes at once: models rejecting, chi-square; bootstrap; pooled model (of 66)",
    c(sum(jt$p_chisq < 0.05), sum(jt$p_boot < 0.05), sum(jt$p < 0.05)))

## --- Section C5: the station-by-station benchmark ---------------------------
dh <- DH[DH$block == "nb_both", ]
csv_write(dh, "tab_app_dh.csv")
num("Section C5", "F test: share of the stations rejecting, range over the models (%)", range(dh$reject_F_pct))
num("Section C5", "HC1 test: share of the stations rejecting, range over the models (%)", range(dh$reject_hc1_pct))
num("Section C5", "Dumitrescu-Hurlin test rejects at 5%: models, F; HC1 (of 66)", c(sum(dh$p_F < 0.05), sum(dh$p_hc1 < 0.05)))
x <- dist_tab[dist_tab$setting == 1, ]
num("Section C5", "setting 1: HC1 rejections by distance to X, 0-0.5; 0.5-1; 1-2; 2-5; 5-10 km (%)", x$hc1[1:5])

## --- Section C6: diagnostics and robustness ---------------------------------
dg <- merge(DG, unique(REG[, c("setting", "model")]))
diag_tab <- do.call(rbind, lapply(1:3, function(s) { x <- dg[dg$setting == s, ]; r <- REG[REG$setting == s, ]
  data.frame(setting = s, regimes = nrow(x), lb = stats::median(x$ljung_box_reject_pct),
             lb_sq = stats::median(x$ljung_box_sq_reject_pct), acf1 = stats::median(x$acf1_median),
             kurt = stats::median(x$kurtosis_median), theta_upper = sum(r$theta_bound == "upper"),
             theta_lower = sum(r$theta_bound == "lower")) }))
csv_write(diag_tab, "tab_app_diagnostics.csv")
num("Section C6, Residuals", "Ljung-Box rejections, median over the regimes, range over the settings (%)", range(diag_tab$lb))
num("Section C6, Residuals", "Ljung-Box on the squares, median over the regimes, range over the settings (%)", range(diag_tab$lb_sq))
num("Section C6, Residuals", "median kurtosis, range over the settings", range(diag_tab$kurt))
num("Section C6, Residuals", "median first-order autocorrelation, range over the settings", range(diag_tab$acf1))
tb <- REG$theta_bound %in% c("upper", "lower")
num("Section C6, The range", "regimes with the range at the upper bound; at the lower bound; of all the regimes",
    c(sum(REG$theta_bound == "upper"), sum(REG$theta_bound == "lower"), nrow(REG)))
num("Section C6, The range", "of these, regimes of at most ten stations", sum(tb & REG$pumps <= 10))
## the local maxima
OPTC <- merge(OPTC, OPTR[, c("setting", "model", "regime", "diff")], by = c("setting", "model", "regime"))
OPTC$z_diff <- (OPTC$beta_ls - OPTC$beta_selected) / OPTC$se_boot
d <- OPTR[is.finite(OPTR$diff), ]
num("Section C6, Local maxima", "regimes refitted; within 1; differing by more than 10; by more than 100",
    c(nrow(d), sum(abs(d$diff) <= 1), sum(abs(d$diff) > 10), sum(abs(d$diff) > 100)))
num("Section C6, Local maxima", "higher maximum by more than 10: least squares; selected fit", c(sum(d$diff > 10), sum(d$diff < -10)))
num("Section C6, Local maxima", "share of regimes differing by more than 10: at most 10 stations; more than 100 (%)",
    c(pct(abs(d$diff[d$size <= 10]) > 10), pct(abs(d$diff[d$size > 100]) > 10)))
big <- OPTC[abs(OPTC$diff) > 10, ]
zb <- sapply(FU_BLOCKS, function(b) round(stats::median(abs(big$z_diff[big$block == b]), na.rm = TRUE), 2))
num("Section C6, Local maxima", "regimes differing by more than 10: median |change|/se of the lag coefficients, by block", zb)
bx <- big[big$block %in% c("nb", "nb_other"), ]
num("Section C6, Local maxima", "coefficients of the lags of X in those regimes; moving by more than two se",
    c(nrow(bx), sum(abs(bx$z_diff) > 2, na.rm = TRUE)))
sig_s <- abs(big$beta_selected / big$se_boot) > 1.96; sig_l <- abs(big$beta_ls / big$se_boot) > 1.96
num("Section C6, Local maxima", "coefficients that change sign: significant at both maxima; at one of them",
    c(sum(sig_s & sig_l & sign(big$beta_selected) != sign(big$beta_ls), na.rm = TRUE),
      sum(xor(sig_s, sig_l) & sign(big$beta_selected) != sign(big$beta_ls), na.rm = TRUE)))
ch <- d[(d$p_nb_both_selected < 0.05) != (d$p_nb_both_ls < 0.05), ]
num("Section C6, Local maxima", "regimes whose 5% decision of the test of leadership changes", nrow(ch))
num("Section C6, Local maxima", "those regimes: setting, model, regime, size, log-lik difference, p selected, p least squares",
    paste(ch$setting, ch$model, ch$regime, ch$size, round(ch$diff, 1), round(ch$p_nb_both_selected, 3),
          round(ch$p_nb_both_ls, 3)))

## --- Section C7 and the caption of Table C7 ---------------------------------
num("Section C7, Table C7 caption", "regimes of the 66 models", nrow(REG))

## --- Section C8: the small samples of setting 2 -----------------------------
csv_write(s2[, c("city", "fuel", "pumps", "K_max", "K", "phi", "sizes", "B_valid", "ari_draws_median")], "tab_app_small.csv")
num("Section C8", "setting 2: pooled models (K = 1)", paste(en(s2$city[s2$K == 1]), s2$fuel[s2$K == 1]))
num("Section C8", "setting 2, K > 1: range of the regime sizes", range(REG$pumps[REG$setting == 2 & REG$K > 1]))
num("Section C8", "setting 2: models with phi = 0; models with phi >= 0.5", c(sum(s2$phi == 0),
    paste(en(s2$city[s2$phi >= 0.5]), s2$fuel[s2$phi >= 0.5], s2$pumps[s2$phi >= 0.5])))
r2 <- REG[REG$setting == 2, ]
num("Section C8", "setting 2: regimes with the range at a bound; regimes", c(sum(r2$theta_bound %in% c("upper", "lower")), nrow(r2)))

NUM <- do.call(rbind, NUM)
csv_write(NUM, "app_numbers.csv")
print(NUM, right = FALSE, row.names = FALSE)





## ===========================================================================
## 8. Figure C1 of Appendix C and the numbers of Sections 5.1, 5.2, C1 and C2:
##    the stored outputs of the exploratory analysis
## ===========================================================================
## These outputs describe the weekly data of data/ before any model is fitted.
## They were computed by the exploratory analysis of the data, kept outside
## this package, and are stored in results/application/eda. They are copied
## here into <out>/eda, so that every output of the application is in one
## place:
##   eda_stations_by_city.csv        Section 5.1 (3,136 stations, 264
##                                   independent); Supplement C, Section C1
##                                   (the table of the stations by city)
##   eda_window.csv                  Section 5.1 and Supplement C, Section C1
##                                   (234 weeks in the window; the stations
##                                   behind the national mean of a Sunday)
##   eda_missing_weeks.csv           Supplement C, Section C2 (the missing
##                                   station-weeks of the lagged series)
##   eda_unit_roots_and_stickiness.csv  Section 5.1 (ADF test); Supplement C,
##                                   Section C1 (ADF test, autocorrelation,
##                                   unchanged prices); Section 5.4 (read in
##                                   section 7 above)
##   eda_pairs.csv                   Section 5.1 and Supplement C, Section C1
##                                   (the stations Y, their distances to X)
##   eda_pacf_share_outside_band.csv, eda_ar_order.csv
##                                   Section 5.2 and Supplement C, Section C2
##                                   (the choice of the lags)
##   completion_fits.csv             the pooled STEM models that completed the
##                                   lagged series (not quoted in the paper)
##   fig_app_centring.png            Figure C1 of Appendix C
##   fig_app_centring_stations.csv   the premiums of its two stations, quoted
##                                   in its text
EDA_FILES <- c("eda_stations_by_city.csv", "eda_window.csv", "eda_missing_weeks.csv",
               "eda_unit_roots_and_stickiness.csv", "eda_pairs.csv", "eda_pacf_share_outside_band.csv",
               "eda_ar_order.csv", "completion_fits.csv", "fig_app_centring.png",
               "fig_app_centring_stations.csv")
if (all(file.exists(file.path(OPT$eda, EDA_FILES)))) {
  dir.create(file.path(OUT, "eda"), showWarnings = FALSE)
  file.copy(file.path(OPT$eda, EDA_FILES), file.path(OUT, "eda"), overwrite = TRUE, copy.date = TRUE)
  message("copied the stored outputs of the exploratory analysis into ", file.path(OUT, "eda"))
} else message("the stored outputs of the exploratory analysis are not in ", OPT$eda, ": skipped")
message("done: outputs in ", OUT)
