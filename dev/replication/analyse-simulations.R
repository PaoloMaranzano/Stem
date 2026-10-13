## ===========================================================================
## SC-STEM: the simulation study. THE ANALYSIS.
##
## One script, self-contained. It needs an installed Stem and its companion
## run-simulations.R in the same folder, from which it takes the design (the
## scenario-variants, the reference levels and the blocks), so that the two
## cannot disagree about it. It runs from whatever folder it sits in.
##
##     Rscript analyse-simulations.R                  the results of tag "main3"
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
##   R9  the information criteria in the two-step rule: (K, phi) selected
##       again with AIC, AICc, GCV, HQ, KIC, BIC and EBIC, from the grids
##   R10 the estimates against the truth, when K is selected correctly and
##       when it is not: the parameters at the true K, the fitted signal of
##       the selected model
##   R11 the summary figure of the manuscript
##
## Inputs  : <here>/results/<tag>_s<stream>.csv, -params.csv, -grid.csv, -timeouts.csv,
##           every stream of the tag stacked; <tag>.csv of the first pass  (--results=)
## Outputs : <here>/output/  tables (.tex, .csv) and figures (.png)  (--out=)
##
## The results files are read defensively. A synchronization client (Google
## Drive) can lose a block of the lines that the runner appends while it is
## uploading a large file: the line after the gap then starts in the middle of
## a row. Such lines are dropped, and so is every replication left with fewer
## rows than the others of its cell. A stream still running counts only the
## replications that every one of its cells has reached, so that every cell
## has the same replications (--balanced=FALSE keeps them all).
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
  tag      = "main3",
  results  = file.path(ana_here, "results"),
  out      = file.path(ana_here, "output"),
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
## the Wilson 95% interval of a proportion p out of M
wilson <- function(p, M, z = 1.96) {
  ctr <- (p + z^2 / (2 * M)) / (1 + z^2 / M)
  hw  <- z * sqrt(p * (1 - p) / M + z^2 / (4 * M^2)) / (1 + z^2 / M)
  cbind(lo = pmax(0, ctr - hw), hi = pmin(1, ctr + hw))
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
  L <- readLines(f, warn = FALSE)
  if (length(L) < 2L) return(NULL)
  sz <- file.size(f)
  con <- file(f, "rb"); seek(con, sz - 1); last <- readBin(con, "raw", 1L); close(con)
  if (last != as.raw(10L)) L <- L[-length(L)]
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
  pat <- sprintf("^%s(_s[0-9]+)?%s[.]csv$", gsub(".", "[.]", tag, fixed = TRUE), suffix)
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
    paste("scenario-variant & $M$ & $\\Pr(\\hat K = K)$ & ARI at $K$ & ARI at $\\hat K$",
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
  core_panels <- function(name, title, draw) {
    nr <- length(ids); nc <- length(CORE_T)
    fig_open(name, 8.6, 1.35 * nr + 1.7)
    graphics::layout(rbind(matrix(seq_len(nr * nc), nr, byrow = TRUE), nr * nc + 1),
                     heights = c(rep(1, nr), 0.5))
    graphics::par(mar = c(1.8, 2.6, 0.9, 0.4), oma = c(0, 7.4, 4.4, 0), mgp = c(1.7, 0.5, 0))
    for (i in ids) for (Tn in CORE_T) {
      graphics::plot(NA, xlim = c(0.4, xpos(length(CORE_N), ng) + 0.6), ylim = c(0, 1), axes = FALSE,
                     xlab = "", ylab = "")
      graphics::abline(h = c(0.25, 0.5, 0.75, 1), col = "#e6eaee")
      graphics::axis(2, at = seq(0, 1, 0.25), labels = c("0", "0.25", "0.5", "0.75", "1"),
                     las = 1, cex.axis = 0.85)
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
  core_panels("fig_sim_core_ari.png",
              "Recovery of the partition: ARI at the true number of regimes, one box per cell",
              function(i, Tn) for (j in seq_along(CORE_N)) for (wi in seq_len(ng)) {
                v <- ari_of[[paste(i, Tn, CORE_N[j], GEO_LEV[wi], sep = ".")]]
                v <- v[is.finite(v)]
                if (!length(v)) next
                graphics::boxplot(v, at = xpos(j, wi), add = TRUE, axes = FALSE, boxwex = 0.9,
                                  col = grDevices::adjustcolor(PAL4[wi], 0.45), border = PAL4[wi],
                                  lwd = 0.7, staplewex = 0.4, outpch = 20, outcex = 0.25)
              })
  ## the share selecting the true K, with its 95% Wilson interval
  core_panels("fig_sim_core_select.png",
              "Selection of the number of regimes: share of replications with K selected = K, with 95% intervals",
              function(i, Tn) for (j in seq_along(CORE_N)) for (wi in seq_len(ng)) {
                q <- tc[tc$scenario == i & tc$TN == Tn & tc$n == CORE_N[j] & tc$geo %in% GEO_LEV[wi], ]
                if (!nrow(q)) next
                ci <- wilson(q$k_ok, q$M)
                graphics::segments(xpos(j, wi), ci[, "lo"], xpos(j, wi), ci[, "hi"],
                                   col = PAL4[wi], lwd = 1.8)
                graphics::points(xpos(j, wi), q$k_ok, pch = 19, col = PAL4[wi], cex = 0.8)
              })
  cat("\n")
} else cat("R2: no replication of the core with three regimes yet\n\n")


## ===========================================================================
## R3. S2: how often a homogeneous network is split
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
  ## the bars are 95% intervals: Wilson for the share selecting the true K,
  ## mean +- 1.96 standard errors for the ARI
  arm_lab <- c(reference = "reference cell", unbalanced = "unbalanced regimes (1:2:3)",
               `knn = 3` = "sparser graph (m = 3)", `knn = 10` = "denser graph (m = 10)",
               `n = 400` = "larger network (d = 400)")
  fig_open("fig_sim_robustness.png", 9.6, 0.5 * length(ids) + 3.0)
  graphics::layout(rbind(c(1, 2), c(3, 3)), heights = c(1, 0.2))
  graphics::par(mar = c(3.8, 7.8, 2.0, 0.8), oma = c(0, 0, 3.6, 0), mgp = c(2.2, 0.6, 0))
  off <- seq(0.3, -0.3, length.out = length(arms))
  for (v in c("k_ok", "ari")) {
    graphics::plot(NA, xlim = c(0, 1), ylim = c(0.5, length(ids) + 0.5), yaxt = "n",
                   xlab = if (v == "ari") "ARI at the true K" else "share of replications with K selected = K",
                   ylab = "", bty = "n", las = 1)
    graphics::axis(2, at = rev(seq_along(ids)), labels = ids, las = 1, cex.axis = 0.9,
                   tick = FALSE)
    graphics::abline(h = seq_along(ids) + 0.5, col = "#e6eaee")
    for (ai in seq_along(arms)) {
      q <- rob[rob$arm == names(arms)[ai] & is.finite(rob[[v]]), ]
      if (!nrow(q)) next
      y <- length(ids) + 1 - match(q$scenario, ids) + off[ai]
      ci <- if (v == "k_ok") wilson(q$k_ok, q$M) else
        cbind(lo = pmax(0, q$ari - 1.96 * q$ari_se), hi = pmin(1, q$ari + 1.96 * q$ari_se))
      graphics::segments(ci[, "lo"], y, ci[, "hi"], y, col = PAL5[ai], lwd = 1.6)
      graphics::points(q[[v]], y, pch = 19, cex = 0.8, col = PAL5[ai])
    }
    graphics::mtext(if (v == "ari") "recovery of the partition" else "selection of the number of regimes",
                    side = 3, line = 0.5, cex = 0.9, font = 2, col = INK)
  }
  graphics::par(mar = c(0, 0, 0, 0)); graphics::plot.new()
  graphics::legend("center", bty = "n", cex = 1.05, pch = 19, lwd = 1.6, col = PAL5, ncol = 3, x.intersp = 0.8,
                   text.width = max(graphics::strwidth(arm_lab, cex = 1.05)),
                   legend = arm_lab[names(arms)])
  graphics::mtext("Robustness: one factor of the design changed at a time, the others kept at the reference cell",
                  outer = TRUE, side = 3, line = 1.9, cex = 1.0, font = 2, col = INK)
  graphics::mtext(paste("reference cell: d = 100 locations, T = 120, medium-high overlap, balanced regimes,",
                        "m = 5 neighbours; bars: 95% intervals; S2 has no partition"),
                  outer = TRUE, side = 3, line = 0.6, cex = 0.85, col = GREY)
  fig_close("fig_sim_robustness.png")
  cat("\n")
} else cat("R4: the robustness blocks have no replication yet\n\n")


## ===========================================================================
## R5. The parameters: bias and RMSE at the reference, the estimated regimes
## aligned on the true ones
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

    fig_open("fig_sim_penalty.png", 6.4, 4.4)
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
SEL_LAB <- c(true = "fit at the true K, every replication", right = "K selected correctly",
             wrong = "K selected wrongly")
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
  fig_open("fig_sim_params.png", 9.0, 11.0)
  graphics::layout(matrix(1:8, 4, byrow = TRUE))
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
  graphics::legend("center", bty = "n", cex = 1.15, border = SEL_COL[c("right", "wrong")],
                   fill = grDevices::adjustcolor(SEL_COL[c("right", "wrong")], 0.5),
                   legend = SEL_LAB[c("right", "wrong")],
                   title = "in the same replication, the two-step rule has")
  graphics::mtext("The parameters at the true number of regimes: relative error (estimate - truth) / |truth|",
                  outer = TRUE, side = 3, line = 1.8, cex = 1.0, font = 2, col = INK)
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
                   legend = c(SEL_LAB[["true"]], "selected fit, K right", "selected fit, K wrong"))
  fig_close("fig_sim_signal.png")
  cat("\n")
}


## ===========================================================================
## R11. The summary figure of the manuscript
## ===========================================================================
## (a) the share of replications selecting the true K, by scenario-variant and
##     (T, d), over the four geometries of the core (S2: the share with K = 1);
## (b) the ARI at the true K at the reference cell;
## (c) the relative error of the parameters at the true K at the reference
##     cell, well and weakly separated scenario-variants apart.
## The scenario-variants run from the most to the least separated (the
## separation measure of Supplementary Material B, Section B6).
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

  fig_open("fig_sim_main.png", 7.6, 7.8)
  graphics::layout(matrix(c(1, 1, 2, 3), 2, byrow = TRUE), heights = c(1.08, 1))
  ## (a)
  graphics::par(mar = c(1.0, 8.4, 6.0, 1.0))
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
                  side = 3, line = 4.3, adj = 0, cex = 0.85, col = INK)
  graphics::mtext(expression("% of replications with" ~ hat(K) == K ~ "(S2:" ~ hat(K) == 1 * "),"
                             ~ "over the four geometries"),
                  side = 3, line = 3.1, adj = 0, cex = 0.75, col = GREY)
  ## (b)
  graphics::par(mar = c(4.2, 8.4, 3.0, 1.0), mgp = c(2.3, 0.6, 0))
  r3 <- rows[rows != "S2"]
  graphics::plot(NA, xlim = c(0, 1), ylim = c(0.5, length(r3) + 0.5), yaxt = "n", bty = "n",
                 xlab = "ARI at the true K", ylab = "", las = 1, cex.axis = 0.85)
  graphics::axis(2, at = rev(seq_along(r3)), labels = r3, las = 1, tick = FALSE, cex.axis = 0.85)
  for (k in seq_along(r3))
    box_at(ref3$ari_true[ref3$scenario == r3[k]], length(r3) + 1 - k,
           if (r3[k] %in% WEAK) "#e0a458" else "#1f6f8b", wex = 0.6, horizontal = TRUE)
  graphics::mtext(expression(bold("(b) Recovery of the partition")), side = 3, line = 1.6,
                  adj = 0, cex = 0.85, col = INK)
  graphics::mtext("reference cell", side = 3, line = 0.5, adj = 0, cex = 0.75, col = GREY)
  ## (c)
  P$rel <- P$err / abs(P$truth)
  P$weak <- P$scenario %in% WEAK
  vs <- c(split(P$rel[!P$weak], P$parameter[!P$weak]), split(P$rel[P$weak], P$parameter[P$weak]))
  graphics::par(mar = c(4.2, 3.6, 3.0, 1.0))
  graphics::plot(NA, xlim = whisk(vs), ylim = c(0.5, length(PAR10) + 0.5), yaxt = "n", bty = "n",
                 xlab = "relative error at the true K", ylab = "", las = 1, cex.axis = 0.85)
  graphics::abline(v = 0, col = INK, lty = 2)
  graphics::axis(2, at = rev(seq_along(PAR10)), labels = as.expression(PAR_EXPR[PAR10]), las = 1,
                 tick = FALSE, cex.axis = 0.95)
  for (k in seq_along(PAR10)) {
    y <- length(PAR10) + 1 - k
    box_at(P$rel[!P$weak & P$parameter == PAR10[k]], y + 0.17, "#1f6f8b", wex = 0.3, horizontal = TRUE)
    box_at(P$rel[P$weak & P$parameter == PAR10[k]], y - 0.17, "#e0a458", wex = 0.3, horizontal = TRUE)
  }
  graphics::mtext(expression(bold("(c) Parameters")), side = 3, line = 1.6, adj = 0, cex = 0.85,
                  col = INK)
  graphics::mtext("reference cell; orange: S1w-shr and S3beta-shr", side = 3,
                  line = 0.5, adj = 0, cex = 0.75, col = GREY)
  fig_close("fig_sim_main.png")
}
cat("\ndone\n")
