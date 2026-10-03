## Generates the conceptual map of the package, in both formats:
##
##   inst/extdata/STEM_function_map.svg   for the README and the vignette
##   inst/extdata/STEM_function_map.pdf   for print and for slides
##
## Both come from the tables below, so the two files cannot drift apart. The
## layout is deliberately explicit rather than computed by a graph-drawing
## algorithm: the point of the map is the reading imposed on the package -- four
## horizontal bands, from the data upwards to the spatially-clustered layer --
## and an automatic layout would not reproduce it.
##
## The coordinate system is the SVG one throughout: 960 by 972 user units, y
## growing downwards. The PDF device is opened so that one user unit is one
## point, and the font sizes below are converted to cex, so that text has the
## same size in both outputs.
##
## Run it from the root of the package:
##
##   Rscript inst/scripts/make-function-map.R
##
## The edge list is documented, function by function, in
## inst/extdata/STEM_function_map.md, which is the place to look when adding a
## function: update that file, then the tables here, then rerun the script.

DY  <- 80
W   <- 960
H   <- 892 + DY
out <- file.path("inst", "extdata")

## ---------------------------------------------------------------------------
## Palette, shared by the two emitters
## ---------------------------------------------------------------------------
pal <- list(
  paper     = "#ffffff",
  band_fill = "#f4f6f8", band_line = "#dfe4ea", band_lab = "#8a939f",
  exp_fill  = "#e8f2f6", exp_line  = "#1f6f8b",
  int_fill  = "#ffffff", int_line  = "#b9c1cb",
  dat_fill  = "#f5efe4", dat_line  = "#b08d47",
  ink       = "#1b2430", ink_soft  = "#3d4753", note = "#6b7480",
  edge      = "#5b6472", accent    = "#1f6f8b"
)

## ---------------------------------------------------------------------------
## Bands: the four layers of the package
## ---------------------------------------------------------------------------
bands <- data.frame(
  x     = c(24, 24, 24, 24),
  y     = c(70, 180, 438 + DY, 574 + DY),
  w     = rep(912, 4),
  h     = c(92, 240 + DY, 118, 262),
  label = c("DATA AND MODEL OBJECT",
            "ESTIMATION ENGINE (KALMAN + EM)",
            "SIMULATION, PREDICTION AND UNCERTAINTY",
            "SPATIALLY-CLUSTERED STEM (SC-STEM)"),
  stringsAsFactors = FALSE
)

## ---------------------------------------------------------------------------
## Boxes: one row per node. `kind` drives the styling and says what the node is:
## "exported" for the user-facing API, "internal" for the machinery, "dataset"
## for the two shipped datasets.
## ---------------------------------------------------------------------------
box <- function(x, y, w, h, label, kind, size = NA) {
  data.frame(x = x, y = y, w = w, h = h, label = label, kind = kind,
             size = size, stringsAsFactors = FALSE)
}

boxes_1 <- rbind(
  box( 44, 100,  96, 34, "pm10",                          "dataset"),
  box(152, 100, 104, 34, "povalley",                      "dataset"),
  box(330, 100, 140, 34, "STEM_Model",                    "exported"),
  box(530, 100, 130, 34, "STEM_Skeleton",                 "internal"),
  box(686, 100, 130, 34, "STEM_Data",                     "internal"),
  box(826, 100, 100, 34, "STEM_Fit",                      "exported")
)
boxes_2 <- rbind(
  box( 44, 212, 150, 36, "STEM_Estimation",               "exported"),
  box(214, 212, 112, 36, "stem_em_fit",                   "internal"),
  box(346, 198, 165, 28, "stem_iterate_plain",            "internal"),
  box(346, 234, 165, 28, "stem_iterate_squarem",          "internal"),
  box(346, 274, 150, 24, "stem_par_vec / _unvec",         "internal", 10.5),
  box(531, 212, 150, 36, "stem_map_em / _ecme",           "internal"),
  box(531, 290, 150, 30, "stem_estep",                    "internal"),
  box(701, 290, 100, 30, "kalman",                        "internal"),
  box(821, 270, 105, 26, "filtering",                     "internal"),
  box(821, 310, 105, 26, "smoothing",                     "internal"),
  box(346, 330, 165, 30, "stem_mstep_em / _ecme",         "internal"),
  box(346, 382, 165, 28, "stem_update_*",                 "internal"),
  box(531, 382, 150, 28, "stem_gls_mean",                 "internal"),
  box(346, 428, 210, 26, "Q_function_addendo1, d1_Q, d2_Q, d12_Q", "internal", 10),
  box(576, 428, 190, 26, "Sigmastar.exp, d1/d2_Sigmastar_*", "internal", 10),
  box( 44, 456, 150, 24, "changedimension_covariates",    "internal", 9.5),
  box(576, 462, 350, 24, "stem_obs_index, stem_blocks_cache, stem_missing_blocks", "internal", 10)
)
boxes_34 <- rbind(
  ## band 3
  box( 44, 472, 150, 34, "STEM_Simulation",               "exported"),
  box(230, 472, 140, 34, "STEM_Kriging",                  "exported"),
  box(230, 515, 140, 30, "spatial.pred",                  "internal"),
  box(420, 472, 150, 34, "STEM_Bootstrap",                "exported"),
  box(420, 515, 150, 30, "STEM_Bootstrap.fn",             "internal"),
  box(596, 472, 128, 34, "STEM_Signal",                    "exported"),
  box(596, 512, 128, 30, "STEM_Complete",                  "exported", 10.5),
  box(742, 472, 140, 34, "SCSTEM_Signal",                  "exported", 11.5),
  box(742, 512, 140, 30, "SCSTEM_Complete",                "exported", 10.5),
  ## band 4
  box( 44, 612, 160, 38, "SCSTEM_Estimation",             "exported"),
  box( 44, 676, 160, 34, "SCSTEM_Infocrit",               "exported"),
  box( 44, 732, 160, 34, "SCSTEM_Select",                 "exported"),
  box( 44, 788, 160, 34, "SCSTEM_Bootstrap",              "exported", 11.5),
  box(250, 788, 185, 34, "SCSTEM_BootInference",          "exported", 11.5),
  box(250, 606, 200, 28, "scstem_neighbors",              "internal"),
  box(250, 640, 200, 28, "scstem_init / _repair_partition",   "internal", 10.5),
  box(470, 606, 200, 28, "scstem_loglike_i",              "internal"),
  box(470, 640, 200, 28, "scstem_swap_pass",              "internal"),
  box(690, 606, 200, 28, "scstem_potts_pairs / _pen_local",   "internal", 10.5),
  box(690, 640, 200, 28, "scstem_npar / _covariate_means",    "internal", 10.5),
  box(470, 700, 200, 28, "scstem_ari",                    "internal"),
  box(690, 700, 200, 28, "scstem_align_labels",           "internal")
)
boxes_34$y <- boxes_34$y + DY
boxes <- rbind(boxes_1, boxes_2, boxes_34)

## ---------------------------------------------------------------------------
## Edges. A straight edge is given by its two endpoints; a curved one by the
## four control points of a cubic Bezier, in the order start, c1, c2, end.
## `accent` marks the spine of the package -- the calls that carry the actual
## estimation -- and `dashed` the two indirect calls, which a static scan of the
## sources cannot see because they go through lapply() and do.call().
## ---------------------------------------------------------------------------
seg <- function(x1, y1, x2, y2, accent = FALSE, dashed = FALSE) {
  list(type = "L", pts = c(x1, y1, x2, y2), accent = accent, dashed = dashed)
}
cur <- function(x1, y1, cx1, cy1, cx2, cy2, x2, y2,
                accent = FALSE, dashed = FALSE) {
  list(type = "C", pts = c(x1, y1, cx1, cy1, cx2, cy2, x2, y2),
       accent = accent, dashed = dashed)
}
## shift the y coordinates of an edge that lie in the two lower bands
shift <- function(e, from = 430) {
  i <- seq(2, length(e$pts), by = 2)
  e$pts[i] <- ifelse(e$pts[i] >= from, e$pts[i] + DY, e$pts[i])
  e
}

edges_1 <- list(
  seg(256, 117, 324, 117),
  seg(470, 117, 524, 117),
  seg(660, 117, 680, 117),
  ## STEM_Fit -> STEM_Estimation, above the boxes of the engine
  cur(876, 138, 860, 178, 200, 172, 150, 210, accent = TRUE)
)
edges_2 <- list(
  ## the spine: STEM_Estimation -> stem_em_fit -> the iterations -> the map
  seg(194, 230, 208, 230, accent = TRUE),
  seg(326, 224, 340, 213, accent = TRUE),
  seg(326, 236, 340, 247, accent = TRUE),
  seg(511, 212, 525, 224, accent = TRUE),
  seg(511, 248, 525, 236, accent = TRUE),
  seg(428, 262, 428, 268),
  ## the map: E-step and M-step
  seg(606, 248, 606, 284, accent = TRUE),
  seg(681, 305, 695, 305, accent = TRUE),
  seg(801, 298, 815, 284, accent = TRUE),
  seg(801, 312, 815, 322, accent = TRUE),
  seg(531, 240, 504, 324, accent = TRUE),
  seg(428, 360, 428, 376),
  seg(511, 350, 525, 390),
  seg(428, 410, 428, 422),
  seg(511, 404, 570, 436),
  ## the missing-data bookkeeping, from the E-step and from the filter
  cur(640, 320, 700, 345, 800, 400, 812, 456),
  cur(926, 283, 934, 320, 934, 430, 910, 456),
  ## STEM_Estimation -> changedimension_covariates (through stem_em_data)
  seg(119, 248, 119, 450)
)
edges_34 <- lapply(list(
  ## band 3
  seg(300, 506, 300, 509),
  seg(495, 506, 495, 509),
  seg(738, 489, 728, 489),
  seg(738, 527, 728, 527),
  seg(419, 530, 200, 500, dashed = TRUE),
  ## band 4: the pipeline
  seg(124, 710, 124, 726, accent = TRUE),
  seg(124, 766, 124, 782, accent = TRUE),
  seg(204, 805, 244, 805, accent = TRUE),
  seg(124, 676, 124, 656, accent = TRUE),
  cur( 30, 693,  14, 693,  14, 631,  38, 631, accent = TRUE),
  ## band 4: the helpers
  seg(204, 626, 244, 620),
  seg(204, 636, 244, 650),
  cur(204, 748, 380, 748, 400, 726, 464, 716),
  cur(435, 800, 560, 800, 640, 746, 700, 728),
  cur(196, 786, 235, 700, 235, 540, 150, 510, dashed = TRUE)
), shift)
edges_x <- list(
  ## SCSTEM_Estimation -> STEM_Estimation: one STEM fit per regime
  cur(60, 612 + DY, 18, 560 + DY, 18, 300, 38, 232, accent = TRUE),
  ## STEM_Bootstrap.fn -> STEM_Estimation
  seg(445, 515 + DY, 160, 254, dashed = TRUE)
)
edges <- c(edges_1, edges_2, edges_34, edges_x)

## ---------------------------------------------------------------------------
## Free text
## ---------------------------------------------------------------------------
titles <- list(
  list(x = 24, y = 30, size = 17, bold = TRUE,  col = pal$ink,
       text = "Stem: how the functions fit together"),
  list(x = 24, y = 50, size = 11, bold = FALSE, col = pal$note,
       text = "Blue boxes are exported; white boxes are internal. Arrows point from caller to callee.")
)

footnotes <- list(
  list(x = 24, y = 856 + DY, text = "SCSTEM_Bootstrap regenerates the data regime by regime and re-runs"),
  list(x = 24, y = 872 + DY, text = "the whole SC-STEM procedure, clustering included, on every draw."),
  list(x = 470, y = 856 + DY, text = "STEM_Fit is the single entry point: it dispatches to STEM_Estimation"),
  list(x = 470, y = 872 + DY, text = "when K = 1 and to SCSTEM_Estimation when K > 1, with or without a penalty on beta.")
)

## Labels set along a vertical edge, rotated a quarter turn counterclockwise.
rotated <- list(
  list(x =  8, y = 662 + DY, text = "grid"),
  list(x = 14, y = 560, text = "one STEM fit per regime")
)

style_of <- function(kind) {
  switch(kind,
         exported = list(fill = pal$exp_fill, line = pal$exp_line, lwd = 1.6,
                         col = pal$ink,      size = 12.5, bold = TRUE),
         dataset  = list(fill = pal$dat_fill, line = pal$dat_line, lwd = 1.4,
                         col = pal$ink,      size = 12.5, bold = TRUE),
         list(fill = pal$int_fill, line = pal$int_line, lwd = 1.1,
              col = pal$ink_soft, size = 11.5, bold = FALSE))
}

## ---------------------------------------------------------------------------
## Bezier sampling, shared by the two emitters: the PDF device draws polylines,
## and the arrowhead needs the direction of the final segment in both formats.
## ---------------------------------------------------------------------------
bezier_pts <- function(p, n = 60) {
  tt <- seq(0, 1, length.out = n)
  b0 <- (1 - tt)^3; b1 <- 3 * (1 - tt)^2 * tt
  b2 <- 3 * (1 - tt) * tt^2; b3 <- tt^3
  cbind(x = b0 * p[1] + b1 * p[3] + b2 * p[5] + b3 * p[7],
        y = b0 * p[2] + b1 * p[4] + b2 * p[6] + b3 * p[8])
}

edge_path <- function(e, n = 60) {
  if (e$type == "L") {
    cbind(x = e$pts[c(1, 3)], y = e$pts[c(2, 4)])
  } else {
    bezier_pts(e$pts, n)
  }
}

## ---------------------------------------------------------------------------
## SVG emitter
## ---------------------------------------------------------------------------
emit_svg <- function(path) {
  esc <- function(s) {
    s <- gsub("&", "&amp;", s, fixed = TRUE)
    s <- gsub("<", "&lt;",  s, fixed = TRUE)
    gsub(">", "&gt;", s, fixed = TRUE)
  }
  L <- character(0)
  add <- function(...) L <<- c(L, paste0(...))

  add('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ', W, ' ', H,
      '" font-family="ui-sans-serif, system-ui, -apple-system, ',
      "'Segoe UI', Helvetica, Arial, sans-serif\">")
  add('  <!-- Generated by inst/scripts/make-function-map.R. Do not edit by hand. -->')
  add('  <defs>')
  for (m in c("arrow", "arrowAccent")) {
    fill <- if (m == "arrow") pal$edge else pal$accent
    add('    <marker id="', m, '" viewBox="0 0 10 10" refX="9.5" refY="5" ',
        'markerWidth="7" markerHeight="7" orient="auto-start-reverse">')
    add('      <path d="M 0 0 L 10 5 L 0 10 z" fill="', fill, '"/>')
    add('    </marker>')
  }
  add('    <style>')
  add('      .band     { fill:', pal$band_fill, '; stroke:', pal$band_line, '; stroke-width:1; rx:10; }')
  add('      .bandlab  { font-size:12px; font-weight:600; fill:', pal$band_lab, '; letter-spacing:.08em; }')
  add('      .exported { fill:', pal$exp_fill, '; stroke:', pal$exp_line, '; stroke-width:1.6; rx:7; }')
  add('      .internal { fill:', pal$int_fill, '; stroke:', pal$int_line, '; stroke-width:1.1; rx:6; }')
  add('      .dataset  { fill:', pal$dat_fill, '; stroke:', pal$dat_line, '; stroke-width:1.4; rx:7; }')
  add('      .edge     { stroke:', pal$edge, '; stroke-width:1.3; fill:none; marker-end:url(#arrow); }')
  add('      .edgeA    { stroke:', pal$accent, '; stroke-width:1.8; fill:none; marker-end:url(#arrowAccent); }')
  add('      .note     { font-size:11px; fill:', pal$note, '; }')
  add('    </style>')
  add('  </defs>')
  add('')
  add('  <rect width="', W, '" height="', H, '" fill="', pal$paper, '"/>')
  add('')

  for (t in titles) {
    add('  <text x="', t$x, '" y="', t$y, '" font-size="', t$size, '"',
        if (t$bold) ' font-weight="700"' else '',
        ' fill="', t$col, '">', esc(t$text), '</text>')
  }
  add('')

  for (i in seq_len(nrow(bands))) {
    b <- bands[i, ]
    add('  <rect class="band" x="', b$x, '" y="', b$y, '" width="', b$w,
        '" height="', b$h, '"/>')
    add('  <text x="', b$x + 12, '" y="', b$y + 20, '" class="bandlab">',
        esc(b$label), '</text>')
  }
  add('')

  for (i in seq_len(nrow(boxes))) {
    b  <- boxes[i, ]
    st <- style_of(b$kind)
    sz <- if (is.na(b$size)) st$size else b$size
    add('  <rect class="', b$kind, '" x="', b$x, '" y="', b$y, '" width="',
        b$w, '" height="', b$h, '"/>')
    add('  <text x="', round(b$x + b$w / 2, 2), '" y="',
        round(b$y + b$h / 2 + sz / 3, 2),
        '" text-anchor="middle" font-size="', sz, '"',
        if (st$bold) ' font-weight="600"' else '',
        ' fill="', st$col, '">', esc(b$label), '</text>')
  }
  add('')

  for (e in edges) {
    cls <- if (e$accent) "edgeA" else "edge"
    d <- if (e$type == "L") {
      paste("M", e$pts[1], e$pts[2], "L", e$pts[3], e$pts[4])
    } else {
      paste("M", e$pts[1], e$pts[2], "C", e$pts[3], e$pts[4],
            e$pts[5], e$pts[6], e$pts[7], e$pts[8])
    }
    add('  <path class="', cls, '" d="', d, '"',
        if (e$dashed) ' stroke-dasharray="4 3"' else '', '/>')
  }
  add('')

  for (r in rotated) {
    add('  <text x="', r$x, '" y="', r$y, '" class="note" transform="rotate(-90 ',
        r$x, ' ', r$y, ')">', esc(r$text), '</text>')
  }
  for (f in footnotes) {
    add('  <text x="', f$x, '" y="', f$y, '" class="note">', esc(f$text), '</text>')
  }
  add('</svg>')

  writeLines(L, path, useBytes = TRUE)
  path
}

## ---------------------------------------------------------------------------
## PDF emitter. The device is opened at 1 user unit = 1 point, so the font
## sizes and stroke widths carry over from the SVG unchanged.
## ---------------------------------------------------------------------------
roundrect <- function(x, y, w, h, r, fill, line, lwd) {
  r <- min(r, w / 2, h / 2)
  ## The outline is built corner by corner. Because y grows downwards, an
  ## increasing angle runs clockwise on screen -- right, then down, then left,
  ## then up -- so each corner sweeps a quarter turn forwards from `from`.
  co <- function(cx, cy, from) {
    a <- seq(from, from + pi / 2, length.out = 12)
    cbind(cx + r * cos(a), cy + r * sin(a))
  }
  p <- rbind(co(x + w - r, y + r,     -pi / 2),     # top-right corner
             co(x + w - r, y + h - r,  0),          # bottom-right
             co(x + r,     y + h - r,  pi / 2),     # bottom-left
             co(x + r,     y + r,      pi))         # top-left
  graphics::polygon(p[, 1], p[, 2], col = fill, border = line, lwd = lwd)
}

arrowhead <- function(x1, y1, x2, y2, col, size = 7) {
  ang <- atan2(y2 - y1, x2 - x1)
  wing <- 0.42
  px <- c(x2,
          x2 - size * cos(ang - wing),
          x2 - size * cos(ang + wing))
  py <- c(y2,
          y2 - size * sin(ang - wing),
          y2 - size * sin(ang + wing))
  graphics::polygon(px, py, col = col, border = NA)
}

emit_pdf <- function(path) {
  grDevices::pdf(path, width = W / 72, height = H / 72, pointsize = 12,
                 bg = pal$paper)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i", family = "sans")
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, W), ylim = c(H, 0))

  ## The device is opened at its default 12 point size, so a font size given in
  ## SVG user units becomes the cex that reproduces it; and lwd is measured in
  ## 1/96 inch against the device's 1/72, hence the 96/72 factor.
  cx0 <- function(size) size / 12
  lw  <- function(x) x * 96 / 72

  for (t in titles) {
    graphics::text(t$x, t$y, t$text, adj = c(0, 0.5), cex = cx0(t$size),
                   col = t$col, font = if (t$bold) 2 else 1)
  }

  for (i in seq_len(nrow(bands))) {
    b <- bands[i, ]
    roundrect(b$x, b$y, b$w, b$h, 10, pal$band_fill, pal$band_line, lw(1))
    graphics::text(b$x + 12, b$y + 20, b$label, adj = c(0, 0.5),
                   cex = cx0(12), col = pal$band_lab, font = 2)
  }

  for (i in seq_len(nrow(boxes))) {
    b  <- boxes[i, ]
    st <- style_of(b$kind)
    sz <- if (is.na(b$size)) st$size else b$size
    roundrect(b$x, b$y, b$w, b$h, 7, st$fill, st$line, lw(st$lwd))
    ## shrink the label if the box is too narrow for it, which keeps the PDF
    ## faithful to the SVG, where the browser would simply let it overflow
    cx <- cx0(sz)
    repeat {
      tw <- graphics::strwidth(b$label, cex = cx, font = if (st$bold) 2 else 1)
      if (tw <= b$w - 10 || cx < cx0(5)) break
      cx <- cx * 0.96
    }
    graphics::text(b$x + b$w / 2, b$y + b$h / 2, b$label, adj = c(0.5, 0.5),
                   cex = cx, col = st$col, font = if (st$bold) 2 else 1)
  }

  for (e in edges) {
    col <- if (e$accent) pal$accent else pal$edge
    wd  <- if (e$accent) lw(1.8) else lw(1.3)
    p   <- edge_path(e)
    n   <- nrow(p)
    graphics::lines(p[, 1], p[, 2], col = col, lwd = wd,
                    lty = if (e$dashed) "42" else "solid")
    arrowhead(p[n - 1, 1], p[n - 1, 2], p[n, 1], p[n, 2], col)
  }

  for (r in rotated) {
    graphics::text(r$x, r$y, r$text, adj = c(0, 0.5), srt = 90,
                   cex = cx0(11), col = pal$note)
  }
  for (f in footnotes) {
    graphics::text(f$x, f$y, f$text, adj = c(0, 0.5), cex = cx0(11), col = pal$note)
  }
  invisible(path)
}

## ---------------------------------------------------------------------------
## Layout check. The coordinates above are written by hand, so a box added to a
## band can silently land under an edge that was routed through what used to be
## empty space. This refuses to emit anything in that case: no edge may run
## through the interior of a box other than the one it leaves or the one it
## enters, no two boxes may overlap, and every box must sit inside a band.
## ---------------------------------------------------------------------------
check_layout <- function() {
  pad <- 2
  bad <- character(0)
  for (k in seq_along(edges)) {
    p <- edge_path(edges[[k]], n = 200)
    for (i in seq_len(nrow(boxes))) {
      b <- boxes[i, ]
      inside <- p[, 1] > b$x + pad & p[, 1] < b$x + b$w - pad &
                p[, 2] > b$y + pad & p[, 2] < b$y + b$h - pad
      j <- which(inside)
      ## an edge is allowed to touch its own endpoints, so only the middle of
      ## the path counts as a crossing
      if (any(j > 0.10 * nrow(p) & j < 0.90 * nrow(p)))
        bad <- c(bad, sprintf("edge %d crosses box '%s'", k, b$label))
    }
  }
  for (i in seq_len(nrow(boxes))) for (j in seq_len(nrow(boxes))) if (i < j) {
    a <- boxes[i, ]; b <- boxes[j, ]
    if (a$x < b$x + b$w && b$x < a$x + a$w && a$y < b$y + b$h && b$y < a$y + a$h)
      bad <- c(bad, sprintf("box '%s' overlaps box '%s'", a$label, b$label))
  }
  for (i in seq_len(nrow(boxes))) {
    b <- boxes[i, ]
    if (!any(bands$x <= b$x & b$x + b$w <= bands$x + bands$w &
             bands$y <= b$y & b$y + b$h <= bands$y + bands$h))
      bad <- c(bad, sprintf("box '%s' is outside every band", b$label))
  }
  if (length(bad))
    stop("the layout is not clean:\n  ", paste(bad, collapse = "\n  "),
         call. = FALSE)
  invisible(TRUE)
}
check_layout()

dir.create(out, recursive = TRUE, showWarnings = FALSE)
svg_path <- emit_svg(file.path(out, "STEM_function_map.svg"))
pdf_path <- emit_pdf(file.path(out, "STEM_function_map.pdf"))

## The vignette carries its own copy of the SVG, so that it renders on a
## machine where the package is not installed.
if (dir.exists("vignettes")) {
  invisible(file.copy(svg_path, file.path("vignettes", "function-map.svg"), overwrite = TRUE))
}

cat("wrote ", svg_path, "\n", "wrote ", pdf_path, "\n",
    nrow(boxes), " boxes, ", length(edges), " edges\n", sep = "")
