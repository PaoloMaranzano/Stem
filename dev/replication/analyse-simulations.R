## ===========================================================================
## SC-STEM: the simulation study. THE ANALYSIS.
##
## One script, self-contained. It needs an installed Stem and its companion
## run-simulations.R in the same folder, from which it takes the design (the
## scenario-variants, the reference levels and the blocks), so that the two
## cannot disagree about it. It runs from whatever folder it sits in.
##
##     Rscript analyse-simulations.R                  the results of tag "main"
##     Rscript analyse-simulations.R --tag=main,ref   several tags, stacked
##     Rscript analyse-simulations.R --out=D:/paper/Figures
##
## Every part is guarded: it produces what the results currently support and
## says what is missing, so it can be run while the study is in progress.
##
##   R1  the reference cell: every scenario-variant at the reference levels
##   R2  the core: selection and recovery against n, T and the geometry
##   R3  S2: how often a homogeneous network is split
##   R4  robustness: balance, graph and n = 400 against the reference
##   R5  the parameters: bias and RMSE at the reference
##   R6  the penalty: recovery against phi, and the phi selected
##   R7  the cost of a replication, by n and T and by geometry
##   R8  the replications stopped by the time limit and drawn again
##
## Inputs  : <here>/results/<tag>_s<stream>.csv, -params.csv, -grid.csv, -timeouts.csv,
##           every stream of the tag stacked; <tag>.csv of the first pass  (--results=)
## Outputs : <here>/output/  tables (.tex, .csv) and figures (.pdf)  (--out=)
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
## The design, from the runner beside this script
## ---------------------------------------------------------------------------
runner <- file.path(ana_here, "run-simulations.R")
if (!file.exists(runner))
  stop("run-simulations.R must sit in the same folder as this script:\n  ",
       ana_here, call. = FALSE)
SIM_DEFINE_ONLY <- TRUE
source(runner, local = globalenv())

## ---------------------------------------------------------------------------
## Options
## ---------------------------------------------------------------------------
AN <- sim_config(list(
  tag     = "main2",
  results = file.path(ana_here, "results"),
  out     = file.path(ana_here, "output")
))
OUT <- AN$out[1]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
OUT <- normalizePath(OUT, winslash = "/")
RES <- normalizePath(AN$results[1], winslash = "/", mustWork = FALSE)
cat("analysis\n  results from ", RES, " (tag ", paste(AN$tag, collapse = ", "),
    ")\n  outputs to   ", OUT, "\n\n", sep = "")

## ---------------------------------------------------------------------------
## Style and small helpers
## ---------------------------------------------------------------------------
PAL3 <- c("#1f6f8b", "#e0a458", "#a8516e")
PAL4 <- c("#1b2430", "#1f6f8b", "#e0a458", "#a8516e")
PAL5 <- c("#1b2430", "#1f6f8b", "#e0a458", "#5b8c5a", "#a8516e")
GREY <- "#8a939f"
INK  <- "#1b2430"

CORE_N  <- SIM_BLOCKS$core$n
CORE_T  <- SIM_BLOCKS$core$TN
CORE_OM <- SIM_BLOCKS$core$omega

## The four geometries of the core, from the most to the least overlapping: the
## overlap omega and what is held fixed as it grows (spread = "total": the
## variance of the whole network; "regime": the variance within a regime).
## Share of the locations nearer the centre of another regime: 67, 32, 20, 10%.
GEO_LEV <- c("0", "0.7", "1b", "1")
GEO_LAB <- c("0" = "total overlap", "0.7" = "medium-high overlap",
             "1b" = "medium-low overlap", "1" = "strong separation")
geo_code <- function(omega, spread) {
  out <- rep(NA_character_, length(omega))
  ok <- !is.na(omega)
  out[ok & abs(omega) < 1e-8] <- "0"
  out[ok & abs(omega - 0.7) < 1e-8] <- "0.7"
  one <- ok & abs(omega - 1) < 1e-8
  out[one] <- ifelse(spread[one] %in% "regime", "1b", "1")
  out
}

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
csv_write <- function(d, name) {
  utils::write.csv(d, file.path(OUT, name), row.names = FALSE)
  message("wrote ", name)
}
fmt <- function(x, d = 2) ifelse(is.na(x), "---", formatC(x, format = "f", digits = d))
tt  <- function(x) sprintf("\\texttt{%s}", x)

## a Monte Carlo mean with its standard error
mc <- function(x) {
  x <- x[is.finite(x)]
  c(mean = if (length(x)) mean(x) else NA_real_,
    se = if (length(x) > 1) stats::sd(x) / sqrt(length(x)) else NA_real_,
    M = length(x))
}
## one row per group of `by`, built by f
agg <- function(d, by, f) {
  sp <- split(d, lapply(d[by], function(v) ifelse(is.na(v), "NA", as.character(v))),
              drop = TRUE)
  out <- do.call(rbind, lapply(sp, function(s) cbind(s[1, by, drop = FALSE], f(s))))
  rownames(out) <- NULL
  out
}
mode_of <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  t <- table(x)
  as.numeric(names(t)[which.max(t)])
}

## the scenario-variants in the order of the design, the additional ones after
scen_order <- function(ids) {
  lev <- c(SIM_SCENARIOS, setdiff(SIM_SCEN$id, SIM_SCENARIOS))
  lev[lev %in% ids]
}

## the cells at the reference levels, one margin possibly moved; in S2 the
## overlap and the balance are vacuous
at_ref <- function(d, ...) {
  m <- utils::modifyList(SIM_REFERENCE, list(...))
  s2 <- d$scenario == "S2"
  d$n == m$n & d$TN == m$TN & d$knn == m$knn &
    (s2 | (abs(d$omega - m$omega) < 1e-8 & d$spread %in% m$spread & d$balance == m$balance))
}

## the metrics of a group of replications. With three regimes: how often k is
## right, and the ARI at the true k and at the selected one; in S2, how often
## the network is split.
metrics <- function(s) {
  k <- mc(s$k_correct); a <- mc(s$ari_true)
  data.frame(
    M       = nrow(s),
    k_ok    = k[["mean"]], k_ok_se = k[["se"]],
    k_under = mean(s$k_hat < s$K_true, na.rm = TRUE),
    k_over  = mean(s$k_hat > s$K_true, na.rm = TRUE),
    ari     = if (s$K_true[1] > 1) a[["mean"]] else NA_real_,
    ari_se  = if (s$K_true[1] > 1) a[["se"]] else NA_real_,
    ari_sel = if (s$K_true[1] > 1) mc(s$ari_sel)[["mean"]] else NA_real_,
    phi_mode = mode_of(s$phi_hat[s$k_hat == s$K_true]),
    ratio   = stats::median(s$rmse_true / s$rmse_pooled, na.rm = TRUE),
    secs    = stats::median(s$secs_total, na.rm = TRUE))
}

## ---------------------------------------------------------------------------
## The results. Several tags are stacked; a (cell, replication) recorded under
## two tags counts once.
## ---------------------------------------------------------------------------
## every stream of every tag, stacked (rep is unique across the streams)
rd <- function(suffix, key) {
  d <- lapply(AN$tag, function(tg) sim_read_streams(RES, tg, suffix))
  d <- d[!vapply(d, is.null, logical(1))]
  if (!length(d)) return(NULL)
  cols <- Reduce(union, lapply(d, names))
  d <- do.call(rbind, lapply(d, function(x) { x[setdiff(cols, names(x))] <- NA; x[cols] }))
  d[!duplicated(d[, key, drop = FALSE]), , drop = FALSE]
}
S <- rd("", c("cell", "rep"))
P <- rd("-params", c("cell", "rep", "regime", "parameter"))
G <- rd("-grid", c("cell", "rep", "k", "phi"))
for (nm in c("S", "P", "G"))
  cat(sprintf("  %-8s %s\n", c(S = "summary", P = "params", G = "grid")[[nm]],
              if (is.null(get(nm))) "missing" else sprintf("%d rows", nrow(get(nm)))))
cat("\n")
if (is.null(S)) stop("no results for tag ", paste(AN$tag, collapse = ", "), " in ", RES,
                     call. = FALSE)

S$error[is.na(S$error)] <- ""
## results written before the fourth geometry and the time limit existed
if (is.null(S$spread)) S$spread <- ifelse(is.na(S$omega), "-", "total")
if (is.null(S$attempt)) S$attempt <- 1L
S$geo <- geo_code(S$omega, S$spread)
if (!is.null(P) && is.null(P$spread)) P$spread <- ifelse(is.na(P$omega), "-", "total")
bad <- S[nzchar(S$error), , drop = FALSE]
if (nrow(bad)) {
  cat("FAILED replications:", nrow(bad), "\n")
  print(utils::head(table(bad$scenario, bad$error), 10))
  cat("\n")
}
S <- S[!nzchar(S$error), , drop = FALSE]

## block membership is a property of the design: a cell is in a block if the
## design puts it there
blk <- lapply(names(SIM_BLOCKS), function(b) sim_cells(blocks = SIM_BLOCKS[b])$cell)
names(blk) <- names(SIM_BLOCKS)
S$in_core <- S$cell %in% unlist(blk[intersect(c("core", "spread"), names(blk))])
cat("replications:", nrow(S), "in", length(unique(S$cell)), "cells;",
    sum(S$in_core), "in the core\n\n")


## ===========================================================================
## R1. The reference cell: every scenario-variant at the reference levels
## ===========================================================================
ref <- S[at_ref(S), , drop = FALSE]
if (nrow(ref)) {
  tr <- agg(ref, "scenario", metrics)
  tr <- tr[match(scen_order(tr$scenario), tr$scenario), ]
  cat("R1. THE REFERENCE CELL (n = ", SIM_REFERENCE$n, ", T = ", SIM_REFERENCE$TN,
      ", omega = ", SIM_REFERENCE$omega, ", balanced, knn = ", SIM_REFERENCE$knn, ")\n", sep = "")
  print(tr, row.names = FALSE, digits = 3)
  csv_write(tr, "tab_sim_reference.csv")
  tex_write(c(
    "\\begin{tabular}{lrrrrrr}", "\\toprule",
    paste("scenario-variant & $M$ & $\\Pr(\\hat k = K)$ & ARI at $K$ & ARI at $\\hat k$",
          "& modal $\\hat\\phi$ & RMSE ratio \\\\"),
    "\\midrule",
    sprintf("%s & %d & %s & %s & %s & %s & %s \\\\", tt(tr$scenario), as.integer(tr$M),
            paste0(fmt(tr$k_ok), " (", fmt(tr$k_ok_se), ")"),
            ifelse(is.na(tr$ari), "---", paste0(fmt(tr$ari), " (", fmt(tr$ari_se), ")")),
            fmt(tr$ari_sel), fmt(tr$phi_mode, 3), fmt(tr$ratio)),
    "\\bottomrule", "\\end{tabular}"), "tab_sim_reference.tex")
  cat("\n")
} else cat("R1: no replication at the reference cell yet\n\n")


## ===========================================================================
## R2. The core: selection and recovery against n, T and the geometry
## ===========================================================================
core <- S[S$in_core & S$K_true > 1, , drop = FALSE]
if (nrow(core)) {
  tc <- agg(core, c("scenario", "geo", "TN", "n"), metrics)
  tc <- tc[order(match(tc$scenario, scen_order(tc$scenario)), match(tc$geo, GEO_LEV), tc$TN, tc$n), ]
  csv_write(tc, "tab_sim_core.csv")
  cat("R2. THE CORE:", nrow(tc), "cells with three regimes\n")
  print(utils::head(tc[, c("scenario", "geo", "TN", "n", "M", "k_ok", "ari", "ari_sel")], 12),
        row.names = FALSE, digits = 3)
  if (nrow(tc) > 12) cat("  ... (", nrow(tc) - 12, " more rows in tab_sim_core.csv)\n", sep = "")

  ## a wide table: scenario and geometry by rows, T and n by columns
  wide_tex <- function(v, name, caption_head) {
    ids <- scen_order(tc$scenario)
    head1 <- paste0(" & & ", paste(sprintf("\\multicolumn{%d}{c}{$T = %d$}", length(CORE_N), CORE_T),
                                   collapse = " & "), " \\\\")
    rules <- paste(sprintf("\\cmidrule(lr){%d-%d}", 3 + (seq_along(CORE_T) - 1) * length(CORE_N),
                           2 + seq_along(CORE_T) * length(CORE_N)), collapse = " ")
    head2 <- paste0("scenario-variant & geometry & ",
                    paste(rep(CORE_N, length(CORE_T)), collapse = " & "), " \\\\")
    body <- unlist(lapply(ids, function(i) vapply(GEO_LEV, function(w) {
      vals <- unlist(lapply(CORE_T, function(Tn) vapply(CORE_N, function(nn) {
        q <- tc[tc$scenario == i & tc$geo %in% w & tc$TN == Tn & tc$n == nn, v]
        if (length(q)) fmt(q) else "---"
      }, character(1))))
      paste0(if (w == GEO_LEV[1]) tt(i) else "", " & ", GEO_LAB[[w]], " & ",
             paste(vals, collapse = " & "), " \\\\")
    }, character(1))))
    tex_write(c(sprintf("%% %s", caption_head),
                sprintf("\\begin{tabular}{ll%s}", strrep("r", length(CORE_T) * length(CORE_N))),
                "\\toprule", head1, rules, head2, "\\midrule", body,
                "\\bottomrule", "\\end{tabular}"), name)
  }
  wide_tex("ari", "tab_sim_core_ari.tex", "ARI at the true number of regimes")
  wide_tex("k_ok", "tab_sim_core_select.tex", "share of replications selecting the true k")

  ## figures: rows the scenario-variants, columns T, n on the x axis, one line
  ## per overlap
  core_figure <- function(v, se_v, ylab, name, title) {
    ids <- scen_order(tc$scenario)
    pdf_open(name, 2.5 * length(CORE_T) + 1.9, 1.75 * length(ids) + 1.0)
    graphics::layout(cbind(matrix(seq_len(length(ids) * length(CORE_T)), nrow = length(ids),
                                  byrow = TRUE), length(ids) * length(CORE_T) + 1),
                     widths = c(rep(1, length(CORE_T)), 0.55))
    graphics::par(mar = c(2.6, 3.0, 1.2, 0.4), oma = c(1.4, 6.2, 2.6, 0), mgp = c(1.7, 0.5, 0))
    for (i in ids) for (Tn in CORE_T) {
      s <- tc[tc$scenario == i & tc$TN == Tn, ]
      graphics::plot(NA, xlim = range(CORE_N), ylim = c(0, 1), log = "x", xaxt = "n",
                     xlab = "", ylab = "", bty = "n", las = 1, cex.axis = 0.8)
      graphics::axis(1, at = CORE_N, labels = CORE_N, cex.axis = 0.8)
      for (wi in seq_along(GEO_LEV)) {
        q <- s[s$geo %in% GEO_LEV[wi], ]
        q <- q[order(q$n), ]
        if (!nrow(q)) next
        se <- q[[se_v]]; se[!is.finite(se)] <- 0
        graphics::polygon(c(q$n, rev(q$n)), pmin(1, pmax(0, c(q[[v]] - 2 * se, rev(q[[v]] + 2 * se)))),
                          col = grDevices::adjustcolor(PAL4[wi], 0.15), border = NA)
        graphics::lines(q$n, q[[v]], col = PAL4[wi], lwd = 1.6)
        graphics::points(q$n, q[[v]], col = PAL4[wi], pch = 19, cex = 0.6)
      }
      if (i == ids[1]) graphics::mtext(sprintf("T = %d", Tn), side = 3, line = 0.4,
                                       cex = 0.8, font = 2, col = INK)
      if (Tn == CORE_T[1]) graphics::mtext(i, side = 2, line = 3.2, las = 1, adj = 1,
                                           cex = 0.72, font = 2, col = INK)
    }
    graphics::par(mar = c(2.6, 0.2, 1.2, 0.2)); graphics::plot.new()
    graphics::legend("center", bty = "n", cex = 0.75, legend = GEO_LAB[GEO_LEV],
                     col = PAL4, lwd = 1.6, pch = 19)
    graphics::mtext("number of locations", side = 1, outer = TRUE, line = 0.2, cex = 0.8)
    graphics::mtext(paste(title, "-", ylab), side = 3, outer = TRUE, line = 1.0, cex = 0.95,
                    font = 2, col = INK)
    pdf_close(name)
  }
  core_figure("ari", "ari_se", "ARI at the true K", "fig_sim_core_ari.pdf",
              "Recovery of the partition")
  core_figure("k_ok", "k_ok_se", "Pr(k selected = K)", "fig_sim_core_select.pdf",
              "Selection of the number of regimes")
  cat("\n")
} else cat("R2: no replication of the core with three regimes yet\n\n")


## ===========================================================================
## R3. S2: how often a homogeneous network is split
## ===========================================================================
nul <- S[S$scenario == "S2", , drop = FALSE]
if (nrow(nul)) {
  tn <- agg(nul, c("n", "TN", "knn"), function(s) {
    f <- mc(as.integer(s$k_hat > 1))
    data.frame(M = f[["M"]], split = f[["mean"]], se = f[["se"]],
               k_hat_max = max(s$k_hat, na.rm = TRUE))
  })
  tn <- tn[order(tn$knn, tn$TN, tn$n), ]
  cat("R3. S2: share of replications with k > 1\n")
  print(tn, row.names = FALSE, digits = 3)
  csv_write(tn, "tab_sim_null.csv")
  tex_write(c(
    "\\begin{tabular}{rrrrr}", "\\toprule",
    "$m$ & $T$ & $n$ & $M$ & $\\Pr(\\hat k > 1)$ \\\\", "\\midrule",
    sprintf("%d & %d & %d & %d & %s (%s) \\\\", as.integer(tn$knn), as.integer(tn$TN),
            as.integer(tn$n), as.integer(tn$M), fmt(tn$split), fmt(tn$se)),
    "\\bottomrule", "\\end{tabular}"), "tab_sim_null.tex")
  cat("\n")
} else cat("R3: no replication of S2 yet\n\n")


## ===========================================================================
## R4. Robustness: balance, graph and n = 400 against the reference
## ===========================================================================
arms <- list(reference     = at_ref(S),
             unbalanced    = at_ref(S, balance = "unbalanced") & S$scenario != "S2",
             `knn = 3`     = at_ref(S, knn = 3L),
             `knn = 10`    = at_ref(S, knn = 10L),
             `n = 400`     = at_ref(S, n = 400L))
rob <- do.call(rbind, lapply(names(arms), function(a) {
  d <- S[arms[[a]], , drop = FALSE]
  if (!nrow(d)) return(NULL)
  cbind(arm = a, agg(d, "scenario", metrics))
}))
if (!is.null(rob) && length(unique(rob$arm)) > 1) {
  rob$arm <- factor(rob$arm, levels = names(arms))
  rob <- rob[order(match(rob$scenario, scen_order(rob$scenario)), rob$arm), ]
  cat("R4. ROBUSTNESS against the reference\n")
  print(rob[, c("scenario", "arm", "M", "k_ok", "ari", "ari_sel")], row.names = FALSE, digits = 3)
  csv_write(rob, "tab_sim_robustness.csv")
  ids <- scen_order(rob$scenario)
  rob_tex <- function(v, name) {
    body <- vapply(ids, function(i) {
      vals <- vapply(names(arms), function(a) {
        q <- rob[rob$scenario == i & rob$arm == a, v]
        if (length(q)) fmt(q) else "---"
      }, character(1))
      paste0(tt(i), " & ", paste(vals, collapse = " & "), " \\\\")
    }, character(1))
    tex_write(c("\\begin{tabular}{lrrrrr}", "\\toprule",
                paste0("scenario-variant & ", paste(names(arms), collapse = " & "), " \\\\"),
                "\\midrule", body, "\\bottomrule", "\\end{tabular}"), name)
  }
  rob_tex("ari", "tab_sim_robustness_ari.tex")
  rob_tex("k_ok", "tab_sim_robustness_select.tex")

  pdf_open("fig_sim_robustness.pdf", 9.0, 0.42 * length(ids) + 2.0)
  graphics::par(mfrow = c(1, 2), mar = c(3.6, 7.2, 2.2, 0.6), oma = c(2.0, 0, 1.6, 0),
                mgp = c(2.0, 0.6, 0))
  for (v in c("ari", "k_ok")) {
    graphics::plot(NA, xlim = c(0, 1), ylim = c(0.5, length(ids) + 0.5), yaxt = "n",
                   xlab = if (v == "ari") "ARI at the true K" else "Pr(k selected = K)",
                   ylab = "", bty = "n", las = 1)
    graphics::axis(2, at = rev(seq_along(ids)), labels = ids, las = 1, cex.axis = 0.8,
                   tick = FALSE)
    graphics::abline(h = seq_along(ids), col = "#eef1f4")
    off <- seq(-0.28, 0.28, length.out = length(arms))
    for (ai in seq_along(arms)) {
      q <- rob[rob$arm == names(arms)[ai], ]
      y <- length(ids) + 1 - match(q$scenario, ids) + off[ai]
      graphics::points(q[[v]], y, pch = 19, cex = 0.8, col = PAL5[ai])
    }
    graphics::mtext(if (v == "ari") "recovery of the partition" else "selection of K",
                    side = 3, line = 0.4, cex = 0.85, font = 2, col = INK)
  }
  graphics::mtext("One margin at a time, against the reference cell", outer = TRUE,
                  side = 3, line = 0.2, cex = 0.95, font = 2, col = INK)
  ## the legend in the outer margin, below both panels
  graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  graphics::plot.new()
  graphics::legend("bottom", bty = "n", cex = 0.8, pch = 19, col = PAL5, horiz = TRUE,
                   legend = names(arms))
  pdf_close("fig_sim_robustness.pdf")
  cat("\n")
} else cat("R4: the robustness blocks have no replication yet\n\n")


## ===========================================================================
## R5. The parameters: bias and RMSE at the reference, the estimated regimes
## aligned on the true ones
## ===========================================================================
if (!is.null(P)) {
  P <- P[P$cell %in% S$cell[at_ref(S)] & P$K_true > 1, , drop = FALSE]
  P$err <- P$estimate - P$truth
}
if (!is.null(P) && nrow(P)) {
  tp <- agg(P, c("scenario", "parameter", "regime"), function(s) {
    e <- s$err[is.finite(s$err)]
    data.frame(M = length(e), truth = s$truth[1],
               mean = if (length(e)) mean(s$estimate[is.finite(s$err)]) else NA_real_,
               bias = if (length(e)) mean(e) else NA_real_,
               rmse = if (length(e)) sqrt(mean(e^2)) else NA_real_)
  })
  PAR <- c("beta1", "beta2", "G", "Sigmaeta", "sigma2eps", "sigma2omega", "theta")
  tp <- tp[order(match(tp$scenario, scen_order(tp$scenario)), match(tp$parameter, PAR),
                 tp$regime), ]
  csv_write(tp, "tab_sim_params.csv")
  cat("R5. PARAMETERS at the reference: see tab_sim_params.csv (",
      nrow(tp), " rows)\n", sep = "")

  ## the reference case in full, regime by regime
  r1 <- tp[tp$scenario == "S1s-ind", ]
  if (nrow(r1)) {
    lab <- c(beta1 = "$\\beta_0$", beta2 = "$\\beta_1$", G = "$G$", Sigmaeta = "$\\sigma^2_\\eta$",
             sigma2eps = "$\\sigma^2_\\epsilon$", sigma2omega = "$\\sigma^2_\\omega$",
             theta = "$\\theta$")
    tex_write(c(
      "\\begin{tabular}{lrrrrrr}", "\\toprule",
      "parameter & regime & $M$ & truth & mean & bias & RMSE \\\\", "\\midrule",
      sprintf("%s & %d & %d & %s & %s & %s & %s \\\\", lab[r1$parameter], as.integer(r1$regime),
              as.integer(r1$M), fmt(r1$truth, 3), fmt(r1$mean, 3), fmt(r1$bias, 3),
              fmt(r1$rmse, 3)),
      "\\bottomrule", "\\end{tabular}"), "tab_sim_params_reference.tex")
  }
  ## every scenario-variant: RMSE relative to the true value, averaged over
  ## the regimes
  rel <- agg(tp, c("scenario", "parameter"), function(s)
    data.frame(rel = mean(s$rmse / abs(s$truth), na.rm = TRUE)))
  ids <- scen_order(rel$scenario)
  body <- vapply(ids, function(i) paste0(tt(i), " & ", paste(vapply(PAR, function(p) {
    q <- rel$rel[rel$scenario == i & rel$parameter == p]
    if (length(q)) fmt(q) else "---"
  }, character(1)), collapse = " & "), " \\\\"), character(1))
  tex_write(c("\\begin{tabular}{lrrrrrrr}", "\\toprule",
              "scenario-variant & $\\beta_0$ & $\\beta_1$ & $G$ & $\\sigma^2_\\eta$ & $\\sigma^2_\\epsilon$ & $\\sigma^2_\\omega$ & $\\theta$ \\\\",
              "\\midrule", body, "\\bottomrule", "\\end{tabular}"), "tab_sim_params_relrmse.tex")
  cat("\n")
} else cat("R5: no parameter estimates at the reference yet\n\n")


## ===========================================================================
## R6. The penalty: recovery against phi at the true k, and the phi selected
## ===========================================================================
if (!is.null(G) && "ari" %in% names(G)) {
  g <- G[G$K_true > 1 & G$k == G$K_true & G$admissible & G$cell %in% S$cell[at_ref(S)], ,
         drop = FALSE]
  if (nrow(g)) {
    tg <- agg(g, c("scenario", "phi"), function(s) {
      a <- mc(s$ari)
      data.frame(M = a[["M"]], ari = a[["mean"]], ari_se = a[["se"]])
    })
    tg <- tg[order(match(tg$scenario, scen_order(tg$scenario)), tg$phi), ]
    csv_write(tg, "tab_sim_penalty.csv")
    cat("R6. ARI at the true k against phi, at the reference\n")
    phis <- sort(unique(tg$phi))
    ids <- scen_order(tg$scenario)
    m <- sapply(phis, function(p) sapply(ids, function(i) {
      q <- tg$ari[tg$scenario == i & abs(tg$phi - p) < 1e-9]
      if (length(q)) q else NA_real_ }))
    m <- matrix(m, nrow = length(ids), dimnames = list(ids, format(phis)))
    print(round(m, 3))

    pdf_open("fig_sim_penalty.pdf", 6.4, 4.4)
    graphics::par(mar = c(4.0, 4.2, 2.2, 7.6), xpd = FALSE)
    x <- seq_along(phis)
    graphics::plot(NA, xlim = range(x), ylim = c(0, 1), xaxt = "n", las = 1, bty = "n",
                   xlab = "spatial penalty phi", ylab = "ARI at the true K")
    graphics::axis(1, at = x, labels = format(phis, drop0trailing = TRUE))
    cols <- grDevices::hcl.colors(length(ids), "Dark 3")
    for (i in seq_along(ids)) {
      graphics::lines(x, m[i, ], col = cols[i], lwd = 1.5)
      graphics::points(x, m[i, ], col = cols[i], pch = 19, cex = 0.6)
    }
    graphics::par(xpd = NA)
    graphics::legend(max(x) + 0.3, 1, bty = "n", cex = 0.72, legend = ids, col = cols,
                     lwd = 1.5, pch = 19)
    graphics::mtext("Recovery against the penalty, at the reference cell", side = 3,
                    line = 0.6, cex = 0.95, font = 2, col = INK)
    pdf_close("fig_sim_penalty.pdf")
  } else cat("R6: no grid at the reference yet\n")
  ## the phi selected, over every cell with three regimes where k is right
  sel <- S[S$K_true > 1 & S$k_hat == S$K_true, , drop = FALSE]
  if (nrow(sel)) {
    ts <- as.data.frame.matrix(table(factor(sel$scenario, levels = scen_order(sel$scenario)),
                                     sel$phi_hat))
    ts <- cbind(scenario = rownames(ts), ts)
    csv_write(ts, "tab_sim_phi_selected.csv")
    cat("\nphi selected when k is right, all cells\n")
    print(ts, row.names = FALSE)
  }
  cat("\n")
} else cat("R6: no grid results, or a grid without the ARI column\n\n")


## ===========================================================================
## R7. The cost of a replication
## ===========================================================================
tcost <- agg(S, c("n", "TN"), function(s)
  data.frame(runs = nrow(s), median = stats::median(s$secs_total, na.rm = TRUE),
             mean = mean(s$secs_total, na.rm = TRUE), max = max(s$secs_total, na.rm = TRUE)))
tcost <- tcost[order(tcost$n, tcost$TN), ]
cat("R7. SECONDS per replication, by n and T\n")
print(tcost, row.names = FALSE, digits = 3)
csv_write(tcost, "tab_sim_cost.csv")
tex_write(c("\\begin{tabular}{rrrrrr}", "\\toprule",
            "$n$ & $T$ & runs & median & mean & max \\\\", "\\midrule",
            sprintf("%d & %d & %d & %s & %s & %s \\\\", as.integer(tcost$n), as.integer(tcost$TN),
                    as.integer(tcost$runs), fmt(tcost$median, 0), fmt(tcost$mean, 0),
                    fmt(tcost$max, 0)),
            "\\bottomrule", "\\end{tabular}"), "tab_sim_cost.tex")

## by geometry, in the core only: the overlap of the regimes changes the
## number of iterations of the alternation and of the EM algorithm
core_c <- S[S$in_core, , drop = FALSE]
if (nrow(core_c)) {
  tgeo <- agg(core_c, c("geo", "n", "TN"), function(s)
    data.frame(runs = nrow(s), median = stats::median(s$secs_total, na.rm = TRUE),
               max = max(s$secs_total, na.rm = TRUE)))
  tgeo <- tgeo[order(match(tgeo$geo, GEO_LEV), tgeo$n, tgeo$TN), ]
  cat("\nR7. SECONDS per replication in the core, by geometry, n and T\n")
  print(tgeo, row.names = FALSE, digits = 3)
  csv_write(tgeo, "tab_sim_cost_geometry.csv")
}
cat("\n")


## ===========================================================================
## R8. The replications stopped by the time limit and drawn again
## ===========================================================================
## Every timeout is logged by the runner in <tag>-timeouts.csv; the replication
## is then drawn again with another seed (attempt 2, 3, ...). A cell with many
## timeouts is one where the estimator is systematically slow, not a few bad
## draws.
TO <- rd("-timeouts", c("cell", "rep", "attempt"))
if (!is.null(TO) && nrow(TO)) {
  if (is.null(TO$spread)) TO$spread <- ifelse(is.na(TO$omega), "-", "total")
  TO$geo <- geo_code(TO$omega, TO$spread)
  done_key <- paste(S$cell, S$rep)          # S holds the completed replications only
  tto <- agg(TO, c("scenario", "geo", "n", "TN"), function(s) {
    reps <- unique(paste(s$cell, s$rep))
    data.frame(timeouts = nrow(s), reps = length(reps), recovered = sum(reps %in% done_key),
               cap_median = stats::median(s$cap, na.rm = TRUE))
  })
  redrawn <- S[S$attempt > 1, , drop = FALSE]
  tto <- tto[order(match(tto$scenario, scen_order(tto$scenario)), match(tto$geo, GEO_LEV),
                   tto$n, tto$TN), ]
  cat("R8. TIMEOUTS by scenario-variant, geometry, n and T (recovered = replications\n",
      "    that succeeded at a later attempt)\n", sep = "")
  print(tto, row.names = FALSE, digits = 3)
  cat(sprintf("\n  %d timeouts in total; %d replications drawn again and completed\n",
              nrow(TO), nrow(redrawn)))
  csv_write(tto, "tab_sim_timeouts.csv")
} else cat("R8: no timeouts logged\n")
cat("\ndone\n")
