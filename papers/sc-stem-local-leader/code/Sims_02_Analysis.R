## ===========================================================================
## Latent Spatial Regimes in Spatio-Temporal State-Space Models: Evidence on
## Local Price Leadership in Italian Fuel Markets
##
## Author:  Paolo Maranzano (University of Milano-Bicocca; Fondazione Eni Enrico
##          Mattei)
## Script:  Sims_02_Analysis.R, the analysis of the simulation study (Section 4
##          of the paper, Appendix B, Supplement B)
##
## What it does: reads the results of the run (Sims_01_Run.R; by default the
## compact results of the study of the paper, shipped in results/simulations)
## and writes every table (.tex, .csv), figure (.png) and number of Section 4
## and Appendix B of the paper and of Supplement B. Part A analyses the
## results; part B draws the figures and tables that describe the design, from
## its definitions (Sims_00_Setup.R), and needs no result. Nothing is fitted.
##
## Outputs, in <out> (output/simulations by default), with the place of
## the paper each reproduces and the part of this script that writes it:
##
##   Part A, the results
##     fig_sim_main.png               Figure 1 (Section 4.2)                   R11
##     fig_sim_signal.png             Figure 2 (Section 4.2)                   R10
##     tab_sim_reference.tex          Supplement B, Table B7 (Section B6.1)     R1
##     fig_sim_core_select.png        Supplement B, Figure B6 (Section B6.2)    R2
##     tab_sim_core_select.tex        Supplement B, Table B8                    R2
##     tab_sim_null.tex               Supplement B, Table B9                    R3
##     tab_sim_core_ari.tex           Supplement B, Table B10 (Section B6.3)    R2
##     fig_sim_core_ari.png           Supplement B, Figure B7                   R2
##     fig_sim_params.png             Supplement B, Figure B8 (Section B6.4)   R10
##     tab_sim_params_reference.tex   Supplement B, Table B11                   R5
##     tab_sim_params_relrmse.tex     Supplement B, Table B12                   R5
##     fig_sim_robustness.png         Supplement B, Figure B9 (Section B6.5)    R4
##     tab_sim_robustness_select.tex  Supplement B, Table B13                   R4
##     tab_sim_robustness_ari.tex     Supplement B, Table B14                   R4
##     fig_sim_penalty.png            Supplement B, Figure B10 (Section B6.6)   R6
##     tab_sim_criteria_correct.tex   Supplement B, Table B15 (Section B6.7)    R9
##     tab_sim_cost.tex               the time of a replication by d and T; not
##                                    in the paper (Section B5 gives the total
##                                    core-hours)                               R7
##   Part B, the design
##     tab_sim_geometry.tex           Supplement B, Table B1 (Section B2.2)
##     fig_sim_geometry.png           Supplement B, Figure B1 (Section B2.2)
##     fig_sim_variance.png           Supplement B, Figure B2 (Section B3.2)
##     fig_sim_correlation.png        Supplement B, Figure B3 (Section B3.2)
##     fig_sim_data.png               Supplement B, Figure B4 (Section B3.2)
##     tab_sim_separation.tex         Supplement B, Table B5 (Section B4)
##     tab_sim_blocks.tex             Supplement B, Table B6 (Section B5)
##     fig_sim_overview.png           Figure B1 of Appendix B; Supplement B,
##                                    Figure B5 (Section B5)
##     tab_sim_scenarios.tex          every parameter value with the derived
##                                    quantities; not in the paper (the values
##                                    are in Table B2 of Supplement B)
## The other outputs (.csv) are the tables behind them and the numbers quoted
## in the text of Section 4 and Supplement B; every part below names what it
## writes.
## ===========================================================================
##
##     Rscript Sims_02_Analysis.R
##     Rscript Sims_02_Analysis.R --results=../output/runs/simulations
##
## or open it in RStudio and press Source. Options (command line): --results
## (the folder of the results), --out (the outputs; default:
## ../output/simulations), --tag (default "main3", the study of the
## paper), --upto (at most this many replications of each stream, e.g.
## --upto=100,100), --balanced.
##
## Part A, by part:
##   R1  the reference cell: every scenario-variant at the reference levels
##   R2  the core: selection and recovery against n, T and the geometry
##   R3  S2: how often a homogeneous network is split
##   R4  robustness: balance, graph and n = 400 against the reference
##   R5  the parameters: bias and RMSE at the reference
##   R6  the penalty: recovery against phi, and the phi selected
##   R7  the cost of a replication, by n and T and by geometry
##   R8  the replications stopped by the time limit and drawn again
##   R9  the information criteria in the two-step rule: (K, phi) selected
##       again with AIC, AICc, GCV, HQ, KIC, BIC and EBIC, from the grids
##   R10 the estimates against the truth, when K is selected correctly and
##       when it is not
##   R11 the summary figure of the paper
##
## The results are read defensively: a line broken by a synchronization client
## while the run appended to a file is dropped, and so is every replication
## left with fewer rows than the others of its cell. A stream still running
## counts only the replications that every one of its cells has reached
## (--balanced=FALSE keeps them all).
##
## Requirements: R (>= 4.2) with spdep and Stem; nothing is fitted, so any
## version of Stem with the function names of 2.0.0 will do (Sims_00_Setup.R
## checks them). Time: about an hour and a half on a desktop computer, most of
## it for R9.
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
## The design and the functions: Sims_01_Run.R in definitions-only mode, which
## loads Sims_00_Setup.R
## ---------------------------------------------------------------------------
SIM_DEFINE_ONLY <- TRUE
SIM_FILE <- file.path(ana_here, "Sims_01_Run.R")
source(SIM_FILE, local = globalenv())

## ---------------------------------------------------------------------------
## Options
## ---------------------------------------------------------------------------
AN <- sim_config(list(
  tag      = "main3",
  results  = file.path(ana_here, "..", "results", "simulations"),
  out      = file.path(ana_here, "..", "output", "simulations"),
  balanced = TRUE,
  ## at most this many replications of each stream (one value per stream, in
  ## order; NA: all those complete), e.g. --upto=100,42
  upto     = NA_integer_
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
## Share of the locations nearer the center of another regime: 67, 32, 20, 10%.
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

## figures in PNG, 300 dpi, sizes in inches
fig_open <- function(name, width, height)
  grDevices::png(file.path(OUT, name), width = width, height = height, units = "in", res = 300,
                 type = if (isTRUE(capabilities("cairo"))) "cairo" else "windows")
fig_close <- function(name) {
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
  k <- mc(s$K_correct); a <- mc(s$ari_true)
  data.frame(
    M       = nrow(s),
    k_ok    = k[["mean"]], k_ok_se = k[["se"]],
    k_under = mean(s$K_hat < s$K_true, na.rm = TRUE),
    k_over  = mean(s$K_hat > s$K_true, na.rm = TRUE),
    ari     = if (s$K_true[1] > 1) a[["mean"]] else NA_real_,
    ari_se  = if (s$K_true[1] > 1) a[["se"]] else NA_real_,
    ari_sel = if (s$K_true[1] > 1) mc(s$ari_sel)[["mean"]] else NA_real_,
    phi_mode = mode_of(s$phi_hat[s$K_hat == s$K_true]),
    ratio   = stats::median(s$rmse_true / s$rmse_pooled, na.rm = TRUE),
    secs    = stats::median(s$secs_total, na.rm = TRUE))
}

## ---------------------------------------------------------------------------
## The results. Several tags are stacked; a (cell, replication) recorded under
## two tags counts once.
## ---------------------------------------------------------------------------
## the column names of the results written before the renaming of the package
## (number of regimes k -> K, parameter count npar -> df)
legacy_names <- function(d) {
  if (is.null(d)) return(d)
  map <- c(k_hat = "K_hat", k_correct = "K_correct", k = "K", k_eff = "K_eff", npar = "df")
  for (o in names(map)) if (o %in% names(d) && !map[[o]] %in% names(d)) names(d)[names(d) == o] <- map[[o]]
  d
}
## One results file, read defensively. A row starts with the quoted cell and
## has as many fields as the header (commas inside quotes do not count); a line
## that does not, the remainder of a row whose beginning was lost, is dropped,
## and so is a last line without its end of line (a file still being written).
an_read_file <- function(f) {
  L <- readLines(f, warn = FALSE)                  # a .csv.gz file is read as it is
  if (length(L) < 2L) return(NULL)
  if (!grepl("[.]gz$", f)) {
    sz <- file.size(f)
    con <- file(f, "rb"); seek(con, sz - 1); last <- readBin(con, "raw", 1L); close(con)
    if (last != as.raw(10L)) L <- L[-length(L)]
  }
  n_quote <- nchar(L) - nchar(gsub("\"", "", L, fixed = TRUE))
  bare <- gsub("\"[^\"]*\"", "", L)
  n_field <- nchar(bare) - nchar(gsub(",", "", bare, fixed = TRUE)) + 1L
  ok <- n_quote %% 2L == 0L & n_field == n_field[1] & startsWith(L, "\"")
  ok[1] <- TRUE
  if (any(!ok)) message(sprintf("  %s: %d broken lines dropped", basename(f), sum(!ok)))
  tf <- tempfile(fileext = ".csv")
  on.exit(unlink(tf))
  writeLines(L[ok], tf)
  utils::read.csv(tf, stringsAsFactors = FALSE)
}
## Every stream of a tag, stacked; in a file with several rows per replication,
## a replication with fewer rows than the others of its cell lost some of them
## and is dropped.
an_read_streams <- function(dir, tag, suffix = "") {
  pat <- sprintf("^%s(_s[0-9]+)?%s[.]csv([.]gz)?$", gsub(".", "[.]", tag, fixed = TRUE), suffix)
  f <- list.files(dir, pattern = pat, full.names = TRUE)
  if (!length(f)) return(NULL)
  d <- lapply(f, an_read_file)
  d <- d[!vapply(d, is.null, logical(1))]
  if (!length(d)) return(NULL)
  cols <- Reduce(union, lapply(d, names))
  d <- do.call(rbind, lapply(d, function(x) { x[setdiff(cols, names(x))] <- NA; x[cols] }))
  if (suffix %in% c("-grid", "-params", "-stations") && all(c("cell", "rep") %in% names(d))) {
    key <- paste(d$cell, d$rep)
    cnt <- table(key)
    rows <- as.integer(cnt[key])
    full <- stats::ave(rows, d$cell, FUN = function(v) as.integer(names(which.max(table(v)))))
    short <- rows < full
    if (any(short)) {
      message(sprintf("  %s%s: %d replications with missing rows dropped", tag, suffix,
                      length(unique(key[short]))))
      d <- d[!short, , drop = FALSE]
    }
  }
  d
}
## every stream of every tag, stacked (rep is unique across the streams)
rd <- function(suffix, key) {
  d <- lapply(AN$tag, function(tg) legacy_names(an_read_streams(RES, tg, suffix)))
  d <- d[!vapply(d, is.null, logical(1))]
  if (!length(d)) return(NULL)
  cols <- Reduce(union, lapply(d, names))
  d <- do.call(rbind, lapply(d, function(x) { x[setdiff(cols, names(x))] <- NA; x[cols] }))
  d[!duplicated(d[, key, drop = FALSE]), , drop = FALSE]
}
S <- rd("", c("cell", "rep"))
P <- rd("-params", c("cell", "rep", "regime", "parameter"))
G <- rd("-grid", c("cell", "rep", "K", "phi"))
for (nm in c("S", "P", "G"))
  cat(sprintf("  %-8s %s\n", c(S = "summary", P = "params", G = "grid")[[nm]],
              if (is.null(get(nm))) "missing" else sprintf("%d rows", nrow(get(nm)))))
cat("\n")
if (is.null(S)) stop("no results for tag ", paste(AN$tag, collapse = ", "), " in ", RES,
                     call. = FALSE)

## A stream still running has reached a replication in some cells and not in
## others (the runner goes replication by replication). Each stream counts the
## replications that all its cells have reached, failed ones included, so that
## every cell has the same replications.
if (isTRUE(AN$balanced[1])) {
  st <- sim_stream_of(S$rep)
  reach <- tapply(S$rep, list(st, S$cell), max)
  upto <- apply(reach, 1, min, na.rm = TRUE)
  for (s in names(upto)) {
    cap <- AN$upto[as.integer(s)]
    if (!is.na(cap)) upto[[s]] <- min(upto[[s]], (as.integer(s) - 1L) * (SIM_REP_MAX + 1L) + cap)
  }
  for (s in names(upto)) {
    k <- st == as.integer(s)
    cat(sprintf("  stream %s: replications %d to %d, complete in all %d cells%s\n", s,
                min(S$rep[k]), upto[[s]], length(unique(S$cell[k])),
                if (any(S$rep[k] > upto[[s]]))
                  sprintf(" (%d later ones left out)", sum(S$rep[k] > upto[[s]])) else ""))
  }
  keep_rep <- function(d) if (is.null(d)) d else
    d[d$rep <= upto[as.character(sim_stream_of(d$rep))], , drop = FALSE]
  S <- keep_rep(S); P <- keep_rep(P); G <- keep_rep(G)
  cat("\n")
}

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
## Writes: tab_sim_reference.tex, Supplement B, Table B7 (Section B6.1)
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
    paste("scenario-variant & $M$ & $\\Pr(\\hat K = K_0)$ & ARI at $K_0$ & ARI at $\\hat K$",
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
## Writes: tab_sim_core_select.tex and fig_sim_core_select.png, Supplement B,
## Table B8 and Figure B6 (Section B6.2); tab_sim_core_ari.tex and
## fig_sim_core_ari.png, Supplement B, Table B10 and Figure B7 (Section B6.3)
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
                    paste(paste0("$d = ", rep(CORE_N, length(CORE_T)), "$"), collapse = " & "), " \\\\")
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

  ## figures: rows the scenario-variants, columns T; in every panel the numbers
  ## of locations are ordered categories, each with the four geometries side by
  ## side, and the legend of the geometries runs below the panels
  ids <- scen_order(tc$scenario)
  ng <- length(GEO_LEV)
  xpos <- function(j, wi) (j - 1) * (ng + 1.2) + wi
  xmid <- vapply(seq_along(CORE_N), function(j) mean(xpos(j, seq_len(ng))), numeric(1))
  core_panels <- function(name, title, draw, ylim = c(0, 1), yat = seq(0, 1, 0.25),
                          ylab = c("0", "0.25", "0.5", "0.75", "1"), href = NULL) {
    nr <- length(ids); nc <- length(CORE_T)
    fig_open(name, 8.6, 1.35 * nr + 1.7)
    graphics::layout(rbind(matrix(seq_len(nr * nc), nr, byrow = TRUE), nr * nc + 1),
                     heights = c(rep(1, nr), 0.5))
    graphics::par(mar = c(1.8, 2.6, 0.9, 0.4), oma = c(0, 7.4, 4.4, 0), mgp = c(1.7, 0.5, 0))
    for (i in ids) for (Tn in CORE_T) {
      graphics::plot(NA, xlim = c(0.4, xpos(length(CORE_N), ng) + 0.6), ylim = ylim, axes = FALSE,
                     xlab = "", ylab = "")
      graphics::abline(h = yat, col = "#e6eaee")
      if (!is.null(href)) graphics::abline(h = href, col = GREY, lty = 2)
      graphics::axis(2, at = yat, labels = ylab, las = 1, cex.axis = 0.85)
      graphics::axis(1, at = xmid, labels = paste0("d = ", CORE_N), tick = FALSE,
                     cex.axis = 0.85, line = -0.7)
      draw(i, Tn)
      if (i == ids[1]) graphics::mtext(sprintf("T = %d", Tn), side = 3, line = 0.4,
                                       cex = 0.9, font = 2, col = INK)
      if (Tn == CORE_T[1]) graphics::mtext(i, side = 2, line = 2.6, las = 1, adj = 1,
                                           cex = 0.85, font = 2, col = INK)
    }
    graphics::par(mar = c(0, 0, 0, 0)); graphics::plot.new()
    graphics::legend("center", ncol = ng, bty = "n", cex = 1.15, legend = GEO_LAB[GEO_LEV],
                     fill = PAL4, border = NA, x.intersp = 0.6)
    graphics::mtext(title, side = 3, outer = TRUE, line = 2.4, cex = 0.95, font = 2, col = INK,
                    padj = 0)
    fig_close(name)
  }
  ## the ARI of every replication, as boxplots
  ari_of <- split(core$ari_true, list(core$scenario, core$TN, core$n, core$geo), drop = TRUE)
  core_box <- function(of) function(i, Tn) for (j in seq_along(CORE_N)) for (wi in seq_len(ng)) {
    v <- of[[paste(i, Tn, CORE_N[j], GEO_LEV[wi], sep = ".")]]
    v <- v[is.finite(v)]
    if (!length(v)) next
    graphics::boxplot(v, at = xpos(j, wi), add = TRUE, axes = FALSE, boxwex = 0.9,
                      col = grDevices::adjustcolor(PAL4[wi], 0.45), border = PAL4[wi],
                      lwd = 0.7, staplewex = 0.4, outpch = 20, outcex = 0.4, outcol = PAL4[wi])
  }
  core_panels("fig_sim_core_ari.png",
              expression(bold("Recovery of the partition: ARI at the true number of regimes") ~ bolditalic(K)[0]),
              core_box(ari_of))
  ## the selection error of every replication, K selected minus K_0, as
  ## boxplots: a boxplot of the 0/1 outcome "K right" would collapse, while
  ## this one shows the direction of the errors (the shares are in the tables)
  core$dK <- core$K_hat - core$K_true
  dk_of <- split(core$dK, list(core$scenario, core$TN, core$n, core$geo), drop = TRUE)
  core_panels("fig_sim_core_select.png",
              expression(bold("Selection of the number of regimes: error") ~ hat(K) - K[0] ~
                           bold("of every replication")),
              core_box(dk_of), ylim = c(-2.2, 1.2), yat = -2:1, ylab = c("-2", "-1", "0", "+1"),
              href = 0)
  cat("\n")
} else cat("R2: no replication of the core with three regimes yet\n\n")


## ===========================================================================
## R3. S2: how often a homogeneous network is split
## Writes: tab_sim_null.tex, Supplement B, Table B9 (Section B6.2)
## ===========================================================================
nul <- S[S$scenario == "S2", , drop = FALSE]
if (nrow(nul)) {
  tn <- agg(nul, c("n", "TN", "knn"), function(s) {
    f <- mc(as.integer(s$K_hat > 1))
    data.frame(M = f[["M"]], split = f[["mean"]], se = f[["se"]],
               k_hat_max = max(s$K_hat, na.rm = TRUE))
  })
  tn <- tn[order(tn$knn, tn$TN, tn$n), ]
  cat("R3. S2: share of replications with k > 1\n")
  print(tn, row.names = FALSE, digits = 3)
  csv_write(tn, "tab_sim_null.csv")
  tex_write(c(
    "\\begin{tabular}{rrrrr}", "\\toprule",
    "$m$ & $T$ & $d$ & $M$ & $\\Pr(\\hat K > 1)$ \\\\", "\\midrule",
    sprintf("%d & %d & %d & %d & %s (%s) \\\\", as.integer(tn$knn), as.integer(tn$TN),
            as.integer(tn$n), as.integer(tn$M), fmt(tn$split), fmt(tn$se)),
    "\\bottomrule", "\\end{tabular}"), "tab_sim_null.tex")
  cat("\n")
} else cat("R3: no replication of S2 yet\n\n")


## ===========================================================================
## R4. Robustness: balance, graph and n = 400 against the reference
## Writes: tab_sim_robustness_select.tex, tab_sim_robustness_ari.tex and
## fig_sim_robustness.png, Supplement B, Tables B13 and B14 and Figure B9
## (Section B6.5)
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

  ## each arm changes one factor of the reference cell and keeps the others;
  ## boxplots over the replications: the selection error K selected - K_0, and
  ## the ARI at K_0
  arm_lab <- c(reference = "reference cell", unbalanced = "unbalanced regimes (1:2:3)",
               `knn = 3` = "sparser graph (m = 3)", `knn = 10` = "denser graph (m = 10)",
               `n = 400` = "larger network (d = 400)")
  S$dK <- S$K_hat - S$K_true
  fig_open("fig_sim_robustness.png", 9.6, 0.62 * length(ids) + 3.0)
  graphics::layout(rbind(c(1, 2), c(3, 3)), heights = c(1, 0.2))
  graphics::par(mar = c(3.8, 7.8, 2.0, 0.8), oma = c(0, 0, 3.6, 0), mgp = c(2.2, 0.6, 0))
  off <- seq(0.32, -0.32, length.out = length(arms))
  for (v in c("dK", "ari_true")) {
    graphics::plot(NA, xlim = if (v == "dK") c(-2.2, 1.2) else c(0, 1),
                   ylim = c(0.5, length(ids) + 0.5), yaxt = "n", xaxt = "n",
                   xlab = if (v == "dK") expression(hat(K) - K[0]) else expression("ARI at" ~ K[0]),
                   ylab = "", bty = "n")
    if (v == "dK") {
      graphics::axis(1, at = -2:1, labels = c("-2", "-1", "0", "+1"))
      graphics::abline(v = 0, col = GREY, lty = 2)
    } else graphics::axis(1)
    graphics::axis(2, at = rev(seq_along(ids)), labels = ids, las = 1, cex.axis = 0.9,
                   tick = FALSE)
    graphics::abline(h = seq_along(ids) + 0.5, col = "#e6eaee")
    for (ai in seq_along(arms)) {
      d <- S[arms[[ai]], , drop = FALSE]
      for (k in seq_along(ids)) {
        x <- d[[v]][d$scenario == ids[k] & (v == "dK" | d$K_true > 1)]
        x <- x[is.finite(x)]
        if (!length(x)) next
        graphics::boxplot(x, at = length(ids) + 1 - k + off[ai], add = TRUE, horizontal = TRUE,
                          axes = FALSE, boxwex = 0.14, col = grDevices::adjustcolor(PAL5[ai], 0.45),
                          border = PAL5[ai], lwd = 0.7, staplewex = 0.4, outpch = 20,
                          outcex = 0.45, outcol = PAL5[ai])
      }
    }
    graphics::mtext(if (v == "ari_true") "recovery of the partition" else "selection of the number of regimes",
                    side = 3, line = 0.5, cex = 0.9, font = 2, col = INK)
  }
  graphics::par(mar = c(0, 0, 0, 0)); graphics::plot.new()
  graphics::legend("center", bty = "n", cex = 1.05, fill = grDevices::adjustcolor(PAL5, 0.45),
                   border = PAL5, ncol = 3, x.intersp = 0.8,
                   text.width = max(graphics::strwidth(arm_lab, cex = 1.05)),
                   legend = arm_lab[names(arms)])
  graphics::mtext("Robustness: one factor of the design changed at a time, the others kept at the reference cell",
                  outer = TRUE, side = 3, line = 1.9, cex = 1.0, font = 2, col = INK)
  graphics::mtext(paste("reference cell: d = 100 locations, T = 120, medium-high overlap, balanced regimes,",
                        "m = 5 neighbors; boxplots over the replications; S2 has no partition"),
                  outer = TRUE, side = 3, line = 0.6, cex = 0.85, col = GREY)
  fig_close("fig_sim_robustness.png")
  cat("\n")
} else cat("R4: the robustness blocks have no replication yet\n\n")


## ===========================================================================
## R5. The parameters: bias and RMSE at the reference, the estimated regimes
## aligned on the true ones
## Writes: tab_sim_params_reference.tex and tab_sim_params_relrmse.tex,
## Supplement B, Tables B11 and B12 (Section B6.4)
## ===========================================================================
## how often the range of a regime ends at a limit of what its distances
## identify, by number of locations and length of the series (all cells)
if (!is.null(P) && "theta_bound" %in% names(P)) {
  tb <- P[P$parameter == "theta" & !is.na(P$theta_bound), , drop = FALSE]
  if (nrow(tb)) {
    tbs <- agg(tb, c("n", "TN"), function(s) data.frame(
      regimes = nrow(s), lower = mean(s$theta_bound == "lower"), upper = mean(s$theta_bound == "upper")))
    csv_write(tbs, "tab_sim_theta_bound.csv")
    cat("R5. RANGES AT A LIMIT (share of the regimes, by n and T): see tab_sim_theta_bound.csv\n")
    print(tbs, row.names = FALSE, digits = 3)
  }
}
P_all <- P          # every cell, for R10 and R11
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
## Writes: fig_sim_penalty.png, Supplement B, Figure B10 (Section B6.6)
## ===========================================================================
if (!is.null(G) && "ari" %in% names(G)) {
  g <- G[G$K_true > 1 & G$K == G$K_true & G$admissible & G$cell %in% S$cell[at_ref(S)], ,
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

    ## one panel per scenario-variant, one boxplot per penalty over the
    ## replications of the reference cell
    fig_open("fig_sim_penalty.png", 9.6, 5.8)
    nc4 <- 4L; nr4 <- ceiling(length(ids) / nc4)
    graphics::par(mfrow = c(nr4, nc4), mar = c(3.4, 3.2, 2.0, 0.6), oma = c(0, 0, 3.4, 0),
                  mgp = c(2.0, 0.55, 0))
    x <- seq_along(phis)
    for (i in ids) {
      graphics::plot(NA, xlim = c(0.4, length(phis) + 0.6), ylim = c(0, 1), xaxt = "n", las = 1,
                     bty = "n", xlab = expression("spatial penalty" ~ phi), ylab = expression("ARI at" ~ K[0]),
                     cex.axis = 0.8)
      graphics::axis(1, at = x, labels = format(phis, drop0trailing = TRUE), cex.axis = 0.75, las = 2)
      for (k in x) {
        v <- g$ari[g$scenario == i & abs(g$phi - phis[k]) < 1e-9]
        v <- v[is.finite(v)]
        if (!length(v)) next
        graphics::boxplot(v, at = k, add = TRUE, axes = FALSE, boxwex = 0.7,
                          col = grDevices::adjustcolor("#1f6f8b", 0.45), border = "#1f6f8b",
                          lwd = 0.7, staplewex = 0.4, outpch = 20, outcex = 0.35, outcol = "#1f6f8b")
      }
      graphics::mtext(i, side = 3, line = 0.4, cex = 0.8, font = 2, col = INK)
    }
    graphics::mtext("Recovery of the partition against the spatial penalty",
                    outer = TRUE, side = 3, line = 1.8, cex = 1.0, font = 2, col = INK)
    graphics::mtext(paste("reference cell: d = 100 locations, T = 120, medium-high overlap,",
                          "balanced regimes, m = 5 neighbors; boxplots over the replications"),
                    outer = TRUE, side = 3, line = 0.5, cex = 0.82, col = GREY)
    fig_close("fig_sim_penalty.png")
  } else cat("R6: no grid at the reference yet\n")
  ## the phi selected, over every cell with three regimes where k is right
  sel <- S[S$K_true > 1 & S$K_hat == S$K_true, , drop = FALSE]
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
## Writes: tab_sim_cost.tex, not in the paper (Section B5 gives the total
## core-hours)
## ===========================================================================
tcost <- agg(S, c("n", "TN"), function(s)
  data.frame(runs = nrow(s), median = stats::median(s$secs_total, na.rm = TRUE),
             mean = mean(s$secs_total, na.rm = TRUE), max = max(s$secs_total, na.rm = TRUE)))
tcost <- tcost[order(tcost$n, tcost$TN), ]
cat("R7. SECONDS per replication, by n and T\n")
print(tcost, row.names = FALSE, digits = 3)
csv_write(tcost, "tab_sim_cost.csv")
tex_write(c("\\begin{tabular}{rrrrrr}", "\\toprule",
            "$d$ & $T$ & runs & median & mean & max \\\\", "\\midrule",
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
## Writes: csv files only, behind the numbers of Section B6 on the
## replications drawn again
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


## ===========================================================================
## R9. The information criteria in the two-step rule
## Writes: tab_sim_criteria_correct.tex, Supplement B, Table B15 (Section
## B6.7)
## ===========================================================================
## The study selects (K, phi) by the two-step rule of SCSTEM_Select() on the
## BIC. The grid of every replication records, for every configuration, the
## log-likelihood l, the number of parameters df and the AIC and BIC, from
## which the number of observed values N follows, log N = (BIC - AIC) / df + 2;
## the rule is applied again here with every criterion:
##   AIC  = -2l + 2 df                 AICc = AIC + 2 df (df + 1) / (N - df - 1)
##   GCV  = -2l - 2 N log(1 - df / N)  HQ   = -2l + 2 log(log(N)) df
##   KIC  = -2l + 3 df                 BIC  = -2l + log(N) df
##   EBIC = BIC + 2 log S(d, K), S the Stirling number of the second kind: the
##          number of partitions of the d locations into K regimes (Chen and
##          Chen 2008, gamma = 1; the check of the paper for wide networks).
## With the BIC the rule must return what the study recorded (agreement below).
AN_CRIT <- c("AIC", "AICc", "GCV", "HQ", "KIC", "BIC", "EBIC")
## log S(d, K) = d log K - log K! + log sum_j (-1)^j C(K, j) (1 - j / K)^d
log_stirling2 <- function(d, K) {
  if (K == 1) return(0)
  j <- 0:(K - 1)
  d * log(K) - lfactorial(K) + log(sum((-1)^j * choose(K, j) * (1 - j / K)^d))
}
## the two-step rule of SCSTEM_Select() on the column `crit` of one grid
an_select <- function(tab, crit, band = c(0.025, 0.2)) {
  ok <- tab$admissible %in% TRUE & is.finite(tab[[crit]])
  adm <- tab[ok & tab$K > 1, , drop = FALSE]
  pooled <- tab[ok & tab$K == 1, , drop = FALSE][seq_len(min(1L, sum(ok & tab$K == 1))), , drop = FALSE]
  phis <- sort(unique(adm$phi))
  in_band <- phis[phis >= band[1] & phis <= band[2]]
  if (length(phis) && length(in_band) < 2L) in_band <- phis
  winners <- vapply(in_band, function(p) {
    sub <- rbind(adm[adm$phi == p, , drop = FALSE], pooled)
    sub$K[order(sub[[crit]], sub$K)][1]
  }, numeric(1))
  if (!length(winners)) winners <- 1
  votes <- table(winners)
  K_sel <- min(as.numeric(names(votes)[votes == max(votes)]))
  if (K_sel == 1) return(data.frame(K = 1, phi = if (nrow(pooled)) pooled$phi else NA_real_,
                                    ari = if (nrow(pooled)) pooled$ari else NA_real_))
  sk <- adm[adm$K == K_sel, , drop = FALSE]
  b <- order(sk[[crit]], sk$phi)[1]
  data.frame(K = K_sel, phi = sk$phi[b], ari = sk$ari[b])
}
if (!is.null(G) && all(c("loglik", "df", "AIC", "BIC", "admissible") %in% names(G))) {
  g9 <- G[G$cell %in% S$cell & is.finite(G$loglik) & G$df > 0, , drop = FALSE]
  g9$N <- exp((g9$BIC - g9$AIC) / g9$df + 2)
  g9$AICc <- g9$AIC + 2 * g9$df * (g9$df + 1) / (g9$N - g9$df - 1)
  g9$GCV  <- -2 * g9$loglik - 2 * g9$N * log(1 - g9$df / g9$N)
  g9$HQ   <- -2 * g9$loglik + 2 * log(log(g9$N)) * g9$df
  g9$KIC  <- -2 * g9$loglik + 3 * g9$df
  g9$EBIC <- g9$BIC + 2 * mapply(log_stirling2, g9$n, g9$K)
  if (!"ari" %in% names(g9)) g9$ari <- NA_real_
  sp <- split(g9, paste(g9$cell, g9$rep), drop = TRUE)
  sel9 <- do.call(rbind, lapply(sp, function(tab) {
    do.call(rbind, lapply(AN_CRIT, function(cr) cbind(tab[1, c("cell", "rep", "scenario", "n", "TN", "K_true")],
                                                       criterion = cr, an_select(tab, cr))))
  }))
  rownames(sel9) <- NULL
  sel9$criterion <- factor(sel9$criterion, levels = AN_CRIT)
  ## the check: with the BIC the rule returns the K the study recorded
  chk <- merge(sel9[sel9$criterion == "BIC", c("cell", "rep", "K")], S[, c("cell", "rep", "K_hat")])
  cat(sprintf("R9. THE CRITERIA IN THE TWO-STEP RULE (%d replications; with the BIC the rule returns the recorded K in %.1f%% of them)\n",
              length(sp), 100 * mean(chk$K == chk$K_hat)))
  sel9$ok <- sel9$K == sel9$K_true
  sel9$over <- sel9$K > sel9$K_true; sel9$under <- sel9$K < sel9$K_true
  sel9$ari[sel9$K_true == 1] <- NA_real_
  wide9 <- function(d, v, by, digits = 2) {
    a <- stats::aggregate(stats::as.formula(paste(v, "~", paste(c(by, "criterion"), collapse = "+"))),
                          data = d, FUN = mean, na.action = stats::na.pass)
    a[[v]] <- round(a[[v]], digits)
    w <- stats::reshape(a, idvar = by, timevar = "criterion", direction = "wide")
    names(w) <- sub(paste0("^", v, "[.]"), "", names(w)); w
  }
  ## by scenario-variant, over every cell: right, too many, too few regimes
  t9ok <- wide9(sel9, "ok", "scenario"); t9ok <- t9ok[match(scen_order(t9ok$scenario), t9ok$scenario), ]
  t9ov <- wide9(sel9, "over", "scenario"); t9ov <- t9ov[match(scen_order(t9ov$scenario), t9ov$scenario), ]
  t9un <- wide9(sel9, "under", "scenario"); t9un <- t9un[match(scen_order(t9un$scenario), t9un$scenario), ]
  t9ari <- wide9(sel9[sel9$K_true > 1, ], "ari", "scenario")
  t9ari <- t9ari[match(scen_order(t9ari$scenario), t9ari$scenario), ]
  ## the over-selection against the size of the network (AIC: penalty 2 per
  ## parameter whatever d; BIC: log(dT); EBIC: adds the partitions)
  t9n <- wide9(sel9, "over", "n")
  cat("\nShare of replications with K selected = K true, by scenario-variant\n"); print(t9ok, row.names = FALSE)
  cat("\nToo many regimes\n"); print(t9ov, row.names = FALSE)
  cat("\nToo few regimes\n"); print(t9un, row.names = FALSE)
  cat("\nARI of the selected partition (scenarios with regimes)\n"); print(t9ari, row.names = FALSE)
  cat("\nToo many regimes, by n\n"); print(t9n, row.names = FALSE)
  csv_write(sel9, "tab_sim_criteria_selections.csv")
  csv_write(t9ok, "tab_sim_criteria_correct.csv"); csv_write(t9ov, "tab_sim_criteria_over.csv")
  csv_write(t9un, "tab_sim_criteria_under.csv"); csv_write(t9ari, "tab_sim_criteria_ari.csv")
  csv_write(t9n, "tab_sim_criteria_over_by_n.csv")
  tex_write(c(paste0("\\begin{tabular}{l", strrep("r", length(AN_CRIT)), "}"), "\\toprule",
              paste("scenario-variant &", paste(AN_CRIT, collapse = " & "), "\\\\"), "\\midrule",
              apply(t9ok, 1, function(r) paste(tt(r[["scenario"]]), "&",
                                               paste(fmt(as.numeric(r[AN_CRIT])), collapse = " & "), "\\\\")),
              "\\bottomrule", "\\end{tabular}"), "tab_sim_criteria_correct.tex")
  cat("\n")
} else cat("R9: no grid results with the log-likelihood and the criteria\n")


## ===========================================================================
## R10. The estimates against the truth, when K is selected correctly and when
## it is not
## Writes: fig_sim_params.png, Supplement B, Figure B8 (Section B6.4), and
## fig_sim_signal.png, Figure 2 of the paper (Section 4.2)
## ===========================================================================
## The runner records, in every replication, the parameters of the fit at the
## true K (phi = 0.05), the estimated regimes aligned on the true ones: they are
## split here by whether the two-step rule selected the true K in the same
## replication. The selected fit, whatever its K, is scored by its fitted signal
## x'beta_g + y^(g) against the conditional mean mu, as a ratio to the RMSE of
## the pooled STEM fit (below 1: better than ignoring the regimes).
PAR10 <- c("beta1", "beta2", "G", "Sigmaeta", "sigma2eps", "sigma2omega", "theta")
PAR_EXPR <- list(beta1 = quote(beta[0]), beta2 = quote(beta[1]), G = quote(G),
                 Sigmaeta = quote(sigma[eta]^2), sigma2eps = quote(sigma[epsilon]^2),
                 sigma2omega = quote(sigma[omega]^2), theta = quote(theta))
SEL_COL <- c(true = GREY, right = "#1f6f8b", wrong = "#e0a458")
SEL_LAB <- list(true = quote("fit at" ~ K[0] * ", every replication"),
                right = quote("selected" ~ hat(K) == K[0]), wrong = quote("selected" ~ hat(K) != K[0]))
## a box with whiskers to 1.5 IQR, the outliers not drawn; at least 10 values
box_at <- function(v, at, col, wex = 0.8, horizontal = FALSE) {
  v <- v[is.finite(v)]
  if (length(v) < 10L) return(invisible(NULL))
  graphics::boxplot(v, at = at, add = TRUE, axes = FALSE, outline = FALSE, boxwex = wex,
                    horizontal = horizontal, col = grDevices::adjustcolor(col, 0.5), border = col,
                    lwd = 0.8, staplewex = 0.4)
}
## the range of the whiskers of a set of groups
whisk <- function(vs) {
  r <- unlist(lapply(vs, function(v) {
    v <- v[is.finite(v)]
    if (length(v) >= 10L) range(grDevices::boxplot.stats(v)$stats) else NULL
  }))
  if (length(r)) range(r) else c(-1, 1)
}

if (!is.null(P_all) && nrow(P_all)) {
  pe <- merge(P_all[P_all$K_true > 1, c("cell", "rep", "scenario", "parameter", "regime",
                                         "truth", "estimate", "K_true")],
              S[S$K_true > 1, c("cell", "rep", "K_hat")], by = c("cell", "rep"))
  pe$rel <- (pe$estimate - pe$truth) / abs(pe$truth)
  pe$sel <- ifelse(pe$K_hat == pe$K_true, "right", "wrong")
  t10 <- agg(pe, c("scenario", "parameter", "sel"), function(s) {
    v <- s$rel[is.finite(s$rel)]
    q <- if (length(v)) stats::quantile(v, c(0.25, 0.5, 0.75)) else rep(NA_real_, 3)
    data.frame(replications = length(unique(paste(s$cell, s$rep))), values = length(v),
               q25 = q[[1]], median = q[[2]], q75 = q[[3]],
               mean_abs = if (length(v)) mean(abs(v)) else NA_real_)
  })
  t10 <- t10[order(match(t10$scenario, scen_order(t10$scenario)), match(t10$parameter, PAR10),
                   t10$sel), ]
  csv_write(t10, "tab_sim_params_by_selection.csv")
  cat("R10. PARAMETERS at the true K, all cells with three regimes: median relative error\n")
  m10 <- stats::xtabs(median ~ scenario + interaction(parameter, sel, sep = ":"), t10)
  print(round(m10[scen_order(rownames(m10)), , drop = FALSE], 3))

  ids10 <- scen_order(unique(pe$scenario))
  sp10 <- split(pe$rel, list(pe$parameter, pe$scenario, pe$sel), drop = TRUE)
  fig_open("fig_sim_params.png", 11.0, 9.6)
  graphics::layout(matrix(1:9, 3, byrow = TRUE))
  graphics::par(mar = c(5.6, 3.8, 2.2, 0.6), oma = c(0, 0, 3.4, 0), mgp = c(2.4, 0.6, 0))
  for (p in PAR10) {
    vs <- sp10[startsWith(names(sp10), paste0(p, "."))]
    graphics::plot(NA, xlim = c(0.5, length(ids10) + 0.5), ylim = whisk(vs), xaxt = "n",
                   xlab = "", ylab = "relative error", bty = "n", las = 1, cex.axis = 0.85)
    graphics::abline(h = 0, col = INK, lty = 2)
    graphics::axis(1, at = seq_along(ids10), labels = ids10, las = 2, cex.axis = 0.85, tick = FALSE)
    for (k in seq_along(ids10)) {
      box_at(sp10[[paste(p, ids10[k], "right", sep = ".")]], k - 0.19, SEL_COL[["right"]], wex = 0.34)
      box_at(sp10[[paste(p, ids10[k], "wrong", sep = ".")]], k + 0.19, SEL_COL[["wrong"]], wex = 0.34)
    }
    graphics::mtext(as.expression(PAR_EXPR[[p]]), side = 3, line = 0.4, cex = 1.1)
  }
  graphics::par(mar = c(0, 0, 0, 0)); graphics::plot.new()
  graphics::legend("center", bty = "n", cex = 1.2, border = SEL_COL[c("right", "wrong")],
                   fill = grDevices::adjustcolor(SEL_COL[c("right", "wrong")], 0.5),
                   legend = as.expression(SEL_LAB[c("right", "wrong")]),
                   title = "in the same replication, the two-step rule has")
  graphics::mtext(expression(bold("The parameters at the true number of regimes") ~ bolditalic(K)[0] *
                               bold(": relative error (estimate - truth) / |truth|")),
                  outer = TRUE, side = 3, line = 1.8, cex = 1.0, col = INK)
  graphics::mtext(paste("every cell with three regimes, the three regimes together;",
                        "boxes from the first to the third quartile, whiskers to 1.5 IQR"),
                  outer = TRUE, side = 3, line = 0.5, cex = 0.85, col = GREY)
  fig_close("fig_sim_params.png")
  cat("\n")
}

sg <- S[S$K_true > 1 & is.finite(S$rmse_pooled) & S$rmse_pooled > 0, , drop = FALSE]
if (nrow(sg)) {
  sg$r_sel  <- sg$rmse_sel / sg$rmse_pooled
  sg$r_true <- sg$rmse_true / sg$rmse_pooled
  sg$sel <- ifelse(sg$K_hat == sg$K_true, "right", ifelse(sg$K_hat < sg$K_true, "too few", "too many"))
  t10s <- agg(sg, c("scenario", "sel"), function(s) {
    q <- stats::quantile(s$r_sel, c(0.25, 0.5, 0.75), na.rm = TRUE)
    data.frame(replications = nrow(s), q25 = q[[1]], median = q[[2]], q75 = q[[3]],
               median_true = stats::median(s$r_true, na.rm = TRUE))
  })
  t10s <- t10s[order(match(t10s$scenario, scen_order(t10s$scenario)), t10s$sel), ]
  csv_write(t10s, "tab_sim_signal.csv")
  cat("R10. FITTED SIGNAL of the selected model, RMSE relative to the pooled STEM fit\n")
  print(t10s, row.names = FALSE, digits = 3)

  ids10 <- scen_order(unique(sg$scenario))
  wrong <- sg$sel != "right"
  fig_open("fig_sim_signal.png", 9.0, 5.0)
  graphics::layout(rbind(1, 2), heights = c(1, 0.16))
  graphics::par(mar = c(6.2, 4.2, 3.6, 0.8), mgp = c(2.6, 0.6, 0))
  vs <- c(split(sg$r_true, sg$scenario), split(sg$r_sel[!wrong], sg$scenario[!wrong]),
          split(sg$r_sel[wrong], sg$scenario[wrong]))
  graphics::plot(NA, xlim = c(0.5, length(ids10) + 0.5), ylim = range(c(whisk(vs), 1)), xaxt = "n",
                 xlab = "", ylab = "RMSE relative to the pooled STEM fit", bty = "n", las = 1)
  graphics::abline(h = 1, col = INK, lty = 2)
  graphics::axis(1, at = seq_along(ids10), labels = ids10, las = 2, cex.axis = 0.9, tick = FALSE)
  for (k in seq_along(ids10)) {
    i <- ids10[k]
    box_at(sg$r_true[sg$scenario == i], k - 0.27, SEL_COL[["true"]], wex = 0.25)
    box_at(sg$r_sel[sg$scenario == i & !wrong], k, SEL_COL[["right"]], wex = 0.25)
    box_at(sg$r_sel[sg$scenario == i & wrong], k + 0.27, SEL_COL[["wrong"]], wex = 0.25)
  }
  graphics::mtext("The fitted signal against the true conditional mean, every cell with three regimes",
                  side = 3, line = 1.9, cex = 1.0, font = 2, col = INK)
  graphics::mtext("RMSE of x'beta + y against mu, as a ratio to the pooled STEM fit (dashed line); below 1, better than ignoring the regimes",
                  side = 3, line = 0.7, cex = 0.82, col = GREY)
  graphics::par(mar = c(0, 0, 0, 0)); graphics::plot.new()
  graphics::legend("center", ncol = 3, bty = "n", cex = 1.0, border = SEL_COL,
                   fill = grDevices::adjustcolor(SEL_COL, 0.5),
                   legend = expression("fit at" ~ K[0] * ", every replication",
                                       "selected fit," ~ hat(K) == K[0], "selected fit," ~ hat(K) != K[0]))
  fig_close("fig_sim_signal.png")
  cat("\n")
}


## ===========================================================================
## R11. The summary figure of the manuscript
## Writes: fig_sim_main.png, Figure 1 of the paper (Section 4.2)
## ===========================================================================
## (a) the share of replications selecting the true K, by scenario-variant and
##     (T, d), over the four geometries of the core (S2: the share with K = 1);
## (b) the ARI at the true K at the reference cell;
## (c) the relative error of the parameters at the true K at the reference
##     cell, well and weakly separated scenario-variants apart.
## The scenario-variants run from the most to the least separated (the
## separation measure of Supplementary Material B, Section B4).
SEP_ORDER <- c("S2", "S1s-ind", "S3error-shr", "S3theta-shr", "S1s-shr", "S3Seta-shr",
               "S3G-shr", "S3beta-shr", "S1w-shr")
WEAK <- c("S3beta-shr", "S1w-shr")
hm_d <- S[(S$in_core & S$K_true > 1) |
            (S$scenario == "S2" & S$knn == SIM_REFERENCE$knn & S$n %in% CORE_N & S$TN %in% CORE_T), ,
          drop = FALSE]
ref3 <- S[at_ref(S) & S$K_true > 1, , drop = FALSE]
if (nrow(hm_d) && nrow(ref3) && !is.null(P) && nrow(P)) {
  rows <- SEP_ORDER[SEP_ORDER %in% hm_d$scenario]
  cols <- expand.grid(n = CORE_N, TN = CORE_T)
  hm <- tapply(hm_d$K_correct, list(hm_d$scenario, paste(hm_d$TN, hm_d$n)), mean)
  hm <- hm[rows, paste(cols$TN, cols$n), drop = FALSE]
  csv_write(data.frame(scenario = rows, hm, check.names = FALSE), "tab_sim_main_select.csv")
  pal <- grDevices::hcl.colors(101, "YlGnBu", rev = TRUE)

  fig_open("fig_sim_main.png", 7.6, 8.3)
  graphics::layout(matrix(c(1, 1, 2, 3), 2, byrow = TRUE), heights = c(1.12, 1))
  graphics::par(oma = c(0, 0, 2.4, 0))
  ## (a)
  graphics::par(mar = c(1.0, 8.4, 7.2, 1.0))
  nr <- length(rows); nc <- nrow(cols)
  graphics::plot(NA, xlim = c(0, nc), ylim = c(0, nr), axes = FALSE, xlab = "", ylab = "",
                 xaxs = "i", yaxs = "i")
  for (r in seq_len(nr)) for (cc in seq_len(nc)) {
    v <- hm[r, cc]
    if (!is.finite(v)) next
    graphics::rect(cc - 1, nr - r, cc, nr - r + 1, col = pal[round(100 * v) + 1], border = "white")
    graphics::text(cc - 0.5, nr - r + 0.5, sprintf("%.0f", 100 * v), cex = 0.9,
                   col = if (v > 0.55) "white" else INK)
  }
  graphics::segments(c(3, 6), 0, c(3, 6), nr, col = "white", lwd = 4)
  graphics::segments(0, nr - 1, nc, nr - 1, col = "white", lwd = 4)
  graphics::axis(2, at = nr - seq_len(nr) + 0.5, labels = rows, las = 1, tick = FALSE,
                 cex.axis = 0.9, line = -0.6)
  graphics::text(seq_len(nc) - 0.5, nr + 0.28, paste0("d = ", cols$n), cex = 0.85, xpd = NA)
  graphics::text(c(1.5, 4.5, 7.5), nr + 0.75, paste0("T = ", CORE_T), cex = 0.95, font = 2, xpd = NA)
  graphics::mtext(expression(bold("(a) Selection of the number of regimes")),
                  side = 3, line = 5.5, adj = 0, cex = 0.85, col = INK)
  graphics::mtext(expression("% of replications with" ~ hat(K) == K[0] ~ "(S2:" ~ hat(K) == 1 * "),"
                             ~ "pooled over the four geometries of the regimes in space:"),
                  side = 3, line = 4.3, adj = 0, cex = 0.75, col = GREY)
  graphics::mtext("total overlap, medium-high overlap, medium-low overlap, strong separation",
                  side = 3, line = 3.3, adj = 0, cex = 0.75, col = GREY)
  ## (b)
  graphics::par(mar = c(4.2, 8.4, 4.0, 1.0), mgp = c(2.3, 0.6, 0))
  r3 <- rows[rows != "S2"]
  graphics::plot(NA, xlim = c(0, 1), ylim = c(0.5, length(r3) + 0.5), yaxt = "n", bty = "n",
                 xlab = expression("ARI at" ~ K[0]), ylab = "", las = 1, cex.axis = 0.85)
  graphics::axis(2, at = rev(seq_along(r3)), labels = r3, las = 1, tick = FALSE, cex.axis = 0.85)
  for (k in seq_along(r3))
    box_at(ref3$ari_true[ref3$scenario == r3[k]], length(r3) + 1 - k,
           if (r3[k] %in% WEAK) "#e0a458" else "#1f6f8b", wex = 0.6, horizontal = TRUE)
  graphics::mtext(expression(bold("(b) Recovery of the partition")), side = 3, line = 2.6,
                  adj = 0, cex = 0.85, col = INK)
  graphics::mtext("reference cell: d = 100, T = 120, medium-high", side = 3, line = 1.5, adj = 0,
                  cex = 0.75, col = GREY)
  graphics::mtext("overlap, balanced regimes, m = 5", side = 3, line = 0.5, adj = 0, cex = 0.75,
                  col = GREY)
  ## (c)
  P$rel <- P$err / abs(P$truth)
  P$weak <- P$scenario %in% WEAK
  vs <- c(split(P$rel[!P$weak], P$parameter[!P$weak]), split(P$rel[P$weak], P$parameter[P$weak]))
  graphics::par(mar = c(4.2, 3.6, 4.0, 1.0))
  graphics::plot(NA, xlim = whisk(vs), ylim = c(0.5, length(PAR10) + 0.5), yaxt = "n", bty = "n",
                 xlab = expression("relative error at" ~ K[0]), ylab = "", las = 1, cex.axis = 0.85)
  graphics::abline(v = 0, col = INK, lty = 2)
  graphics::axis(2, at = rev(seq_along(PAR10)), labels = as.expression(PAR_EXPR[PAR10]), las = 1,
                 tick = FALSE, cex.axis = 0.95)
  for (k in seq_along(PAR10)) {
    y <- length(PAR10) + 1 - k
    box_at(P$rel[!P$weak & P$parameter == PAR10[k]], y + 0.17, "#1f6f8b", wex = 0.3, horizontal = TRUE)
    box_at(P$rel[P$weak & P$parameter == PAR10[k]], y - 0.17, "#e0a458", wex = 0.3, horizontal = TRUE)
  }
  graphics::mtext(expression(bold("(c) Parameters")), side = 3, line = 2.6, adj = 0, cex = 0.85,
                  col = INK)
  graphics::mtext("reference cell, as in (b); blue: well separated", side = 3,
                  line = 1.5, adj = 0, cex = 0.75, col = GREY)
  graphics::mtext("scenario-variants, orange: S1w-shr and S3beta-shr", side = 3,
                  line = 0.5, adj = 0, cex = 0.75, col = GREY)
  graphics::mtext("Simulation study: number of regimes, partition and parameters",
                  outer = TRUE, side = 3, line = 0.8, cex = 1.0, font = 2, col = INK)
  fig_close("fig_sim_main.png")
}
cat("\ndone\n")

## ===========================================================================
## PART B: THE FIGURES AND TABLES OF THE DESIGN
## ===========================================================================
## From the definitions and the generator of Sims_00_Setup.R, so
## that the description of the design cannot disagree with the design that is
## run. One data set is drawn for some of the figures; nothing is fitted.
FG <- list(out = OUT)
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
## figures in PNG, 300 dpi, sizes in inches
fig_open <- function(name, width, height)
  grDevices::png(file.path(OUT, name), width = width, height = height, units = "in", res = 300,
                 type = if (isTRUE(capabilities("cairo"))) "cairo" else "windows")
fig_close <- function(name) {
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
## share of the locations of a regime nearer the center of another regime, by
## Monte Carlo on a fixed seed (a draw from N(mu_g, nu I) falls outside the
## 120-degree wedge of its own center)
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
## Figure B1 of Supplement B (Section B2.2), fig_sim_geometry.png: the four
## geometries
## ===========================================================================
fig_open("fig_sim_geometry.png", 10, 10)
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
  graphics::title(sprintf("(%s) %s: omega = %s, spread fixed %s\ncenters %.1f apart, regime sd %.2f, %.0f%% of the locations nearer another center",
                          letters[i], GEO$label[i], format(GEO$omega[i]),
                          if (GEO$spread[i] == "total") "over the network" else "within a regime",
                          2 * GEO$omega[i], sqrt(GEO$nu[i]), 100 * GEO$share_other[i]),
                  cex.main = 0.9, font.main = 1)
}
fig_close("fig_sim_geometry.png")


## ===========================================================================
## Figure B1 of Appendix B and Figure B5 of Supplement B (Section B5),
## fig_sim_overview.png: an overview of the design. Above, where
## the regimes lie: the four geometries and the unbalanced sizes, at d = 100.
## Below, the series of every scenario-variant at the reference cell, from the
## first replication: three locations of every regime (thin) and the mean of
## the response over the locations of each regime (thick).
## ===========================================================================
OV_DESC <- c(`S2` = "no regimes", `S1s-ind` = "every block differs, strongly; independent regimes",
             `S1s-shr` = "every block differs, strongly; shared shocks",
             `S1w-shr` = "every block differs, weakly",
             `S3beta-shr` = "common regression coefficients", `S3G-shr` = "common persistence G",
             `S3Seta-shr` = "common innovation variance", `S3theta-shr` = "common range",
             `S3error-shr` = "common measurement error")
ov <- lapply(stats::setNames(SIM_SCENARIOS, SIM_SCENARIOS), function(id) dgp_draw(ref_cell(id), 1L))
ov_lines <- function(dd) {
  out <- list()
  for (g in sort(unique(dd$labels))) {
    ii <- which(dd$labels == g)
    out[[length(out) + 1L]] <- list(g = g, thin = dd$z[, ii[1:3], drop = FALSE],
                                    thick = rowMeans(dd$z[, ii, drop = FALSE]))
  }
  out
}
ov_l <- lapply(ov, ov_lines)
ov_ylim <- range(unlist(lapply(ov_l, function(l) lapply(l, function(q) c(q$thin, q$thick)))))
fig_open("fig_sim_overview.png", 11, 10.2)
graphics::layout(rbind(rep(1:5, each = 3), rep(6:8, each = 5), rep(9:11, each = 5),
                       rep(12:14, each = 5)), heights = c(0.82, 1, 1, 1))
graphics::par(oma = c(0.5, 0, 0.5, 0), mgp = c(2.0, 0.55, 0), las = 1, tcl = -0.3)
maps <- rbind(GEO[, c("omega", "spread", "label")],
              data.frame(omega = SIM_REFERENCE$omega, spread = SIM_REFERENCE$spread,
                         label = "unbalanced sizes 1:2:3", stringsAsFactors = FALSE))
maps$balance <- c(rep("balanced", nrow(GEO)), "unbalanced")
for (i in seq_len(nrow(maps))) {
  graphics::par(mar = c(1.2, 1.2, 2.6, 1.2))
  loc <- dgp_locations(SIM_REFERENCE$n, 3L, maps$omega[i], balance = maps$balance[i],
                       seed = 1001L, spread = maps$spread[i])
  mu <- dgp_centres(3L, maps$omega[i])
  graphics::plot(loc$coords, col = RC[loc$labels], pch = 16, cex = 0.75, asp = 1,
                 xlim = c(-3, 3), ylim = c(-3, 3), xlab = "", ylab = "", axes = FALSE)
  graphics::box(col = "#d5dbe1")
  graphics::points(mu[, 1], mu[, 2], pch = 3, col = INK, cex = 1.3, lwd = 2)
  graphics::title(sprintf("(%s) %s", letters[i], maps$label[i]), cex.main = 1.0, font.main = 1,
                  line = 0.6)
}
for (k in seq_along(SIM_SCENARIOS)) {
  id <- SIM_SCENARIOS[k]
  graphics::par(mar = c(3.0, 3.4, 2.6, 0.8))
  graphics::plot(NA, xlim = c(1, SIM_REFERENCE$TN), ylim = ov_ylim, xlab = "t", ylab = "response",
                 bty = "n", cex.axis = 0.85)
  for (q in ov_l[[id]]) {
    col <- if (id == "S2") RC[1] else RC[q$g]
    for (j in seq_len(ncol(q$thin)))
      graphics::lines(q$thin[, j], col = grDevices::adjustcolor(col, 0.35), lwd = 0.7)
  }
  for (q in ov_l[[id]])
    graphics::lines(q$thick, col = if (id == "S2") RC[1] else RC[q$g], lwd = 2)
  graphics::title(sprintf("(%s) %s: %s", letters[nrow(maps) + k], id, OV_DESC[[id]]),
                  cex.main = 0.95, font.main = 1, line = 0.8, adj = 0)
}
fig_close("fig_sim_overview.png")


## ===========================================================================
## Figure B4 of Supplement B (Section B3.2), fig_sim_data.png: one data set in
## time and in space, at the two levels of separation
## ===========================================================================
draws <- list(`S1-strong-ind` = dgp_draw(ref_cell("S1s-ind"), 1L),
              `S1-weak-shr`   = dgp_draw(ref_cell("S1w-shr"), 1L))
TT <- 60L                                          # the period shown in space
map_col <- function(v, lim) {
  cl <- grDevices::hcl.colors(21, "Blue-Red 2")
  cl[findInterval(pmax(pmin(v, lim - 1e-9), -lim + 1e-9), seq(-lim, lim, length.out = 22))]
}
## one color scale per kind of map, common to the two data sets
## (the 90th percentile of the absolute values, beyond which the color saturates)
lim_mean <- stats::quantile(unlist(lapply(draws, function(dd) { m <- colMeans(dd$z); abs(m - mean(m)) })), 0.9)
lim_err  <- stats::quantile(unlist(lapply(draws, function(dd) abs(dd$z[TT, ] - dd$mu[TT, ]))), 0.9)
fig_open("fig_sim_data.png", 12, 11)
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
fig_close("fig_sim_data.png")


## ===========================================================================
## Figure B3 of Supplement B (Section B3.2), fig_sim_correlation.png: the
## correlations and their decay
## ===========================================================================
## nearest-neighbor spacing of a network of n locations at the reference
nn_spacing <- vapply(SIM_BLOCKS$core$n, function(n) {
  xy <- dgp_locations(n, 3L, SIM_REFERENCE$omega, seed = 1001L)$coords
  D <- as.matrix(stats::dist(xy)); diag(D) <- Inf
  stats::median(apply(D, 1, min))
}, numeric(1))
fig_open("fig_sim_correlation.png", 11, 8.5)
par_base(mfrow = c(2, 2), mar = c(3.6, 3.8, 3, 1))
hh <- seq(0, 4, length.out = 300)
for (lv in c("S1-weak", "S1-strong")) {
  psi <- PSI[[lv]]
  graphics::plot(NA, xlim = c(0, 4), ylim = c(0, 0.55), xlab = "distance between two locations",
                 ylab = "correlation of the measurement errors", bty = "n")
  for (g in 1:3)
    graphics::lines(hh, psi$so[g] / (psi$se[g] + psi$so[g]) * exp(-psi$theta[g] * hh), col = RC[g], lwd = 2)
  graphics::abline(v = nn_spacing, col = GREY, lty = 3)
  graphics::text(max(nn_spacing), 0.02, sprintf("nearest-neighbor spacing\nat n = %s",
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
fig_close("fig_sim_correlation.png")


## ===========================================================================
## Figure B2 of Supplement B (Section B3.2), fig_sim_variance.png: the variance
## of the response by source, regime by regime
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
fig_open("fig_sim_variance.png", 11, 5.6)
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
fig_close("fig_sim_variance.png")


## ===========================================================================
## tab_sim_scenarios.tex, not in the paper (the values are in Table B2 of
## Supplement B): the parameter values of every regime
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
## Table B1 of Supplement B (Section B2.2), tab_sim_geometry.tex: the four
## geometries
## ===========================================================================
tex_write(c("\\begin{tabular}{lcl rrrrr}", "\\toprule",
            "geometry & $\\omega$ & fixed spread & center distance & regime sd & distance in sd & nearer another center & network sd \\\\",
            "\\midrule",
            sprintf("%s & %s & %s & %s & %s & %s & %s\\%% & %s \\\\", GEO$label, fmt(GEO$omega, 1),
                    ifelse(GEO$spread == "total", "network", "regime"), fmt(2 * GEO$omega), fmt(sqrt(GEO$nu)),
                    fmt(2 * GEO$omega / sqrt(GEO$nu)), fmt(100 * GEO$share_other, 0), fmt(GEO$network_sd)),
            "\\bottomrule", "\\end{tabular}"), "tab_sim_geometry.tex")


## ===========================================================================
## Table B5 of Supplement B (Section B4), tab_sim_separation.tex: the
## separation of the regimes seen by the assignment step
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
## the nine scenario-variants of the study; the additional ones defined in
## SIM_SCEN are not run, and not described
SEP <- do.call(rbind, lapply(setdiff(SIM_SCENARIOS, "S2"), function(id) {
  row <- SIM_SCEN[SIM_SCEN$id == id, , drop = FALSE]
  p <- dgp_psi(row)
  pr <- list(c(1, 2), c(2, 1), c(2, 3), c(3, 2))
  v <- lapply(pr, function(q) sep_D(p, row$rho, q[1], q[2]))
  Ds <- vapply(v, `[[`, 1, "D")
  j <- which.min(Ds)
  data.frame(id = id, rho = row$rho, D12 = Ds[1], D21 = Ds[2], D23 = Ds[3], D32 = Ds[4],
             minD = Ds[j], mean = v[[j]][["mean"]], latent = v[[j]][["latent"]],
             T60 = 60 * Ds[j], T120 = 120 * Ds[j], T365 = 365 * Ds[j],
             stringsAsFactors = FALSE)
}))
tex_write(c("\\begin{tabular}{lc rrrr r rr rrr}", "\\toprule",
            "scenario-variant & $\\rho$ & $D_{12}$ & $D_{21}$ & $D_{23}$ & $D_{32}$ & $\\min D$ & mean & latent & $60\\min D$ & $120\\min D$ & $365\\min D$ \\\\",
            "\\midrule",
            sprintf("%s & %d & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
                    sprintf("\\texttt{%s}", SEP$id), as.integer(SEP$rho), fmt(SEP$D12, 3), fmt(SEP$D21, 3),
                    fmt(SEP$D23, 3), fmt(SEP$D32, 3), fmt(SEP$minD, 3), fmt(SEP$mean), fmt(SEP$latent, 3),
                    fmt(SEP$T60, 1), fmt(SEP$T120, 1), fmt(SEP$T365, 1)),
            "\\bottomrule", "\\end{tabular}"), "tab_sim_separation.tex")


## ===========================================================================
## Table B6 of Supplement B (Section B5), tab_sim_blocks.tex: the blocks of the
## design and their cells
## ===========================================================================
cfg_all <- CFG
all_cells <- sim_cells(cfg_all, SIM_BLOCKS)
blk_desc <- c(core    = "$n \\in \\{40, 100, 200\\}$, $T \\in \\{60, 120, 365\\}$, $\\omega \\in \\{0, 0.7, 1\\}$",
              spread  = "$n$ and $T$ as the core, $\\omega = 1$ with the spread of a regime fixed",
              balance = "regime sizes $1:2:3$",
              knn     = "$m \\in \\{3, 10\\}$ neighbors",
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
