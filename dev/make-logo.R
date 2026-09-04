## Renders man/figures/logo.png from the same geometry as man/figures/logo.svg.
##
## The SVG is the source of truth for the design; this script exists because the
## SVG carries live text, so the wordmark depends on whichever serif the viewer
## happens to have. The PNG is rendered here once, with the intended face, and is
## the artifact to use wherever the rendering must be fixed (slides, stickers,
## anything printed).
##
## The whole composition is laid out in the SVG user space, 100 wide by 115 tall,
## and every element is designed to fall inside the hexagon, so no clipping is
## needed - which matters, because base graphics can only clip to a rectangle.
##
##   Rscript dev/make-logo.R

out    <- file.path("man", "figures", "logo.png")
width  <- 480L                      # 2x the nominal 240 px, for crisp edges
height <- as.integer(round(width * 115 / 100))
res    <- 144L

S       <- width / 100              # device pixels per user unit
lwd_of  <- function(w) w * S * 96 / res
ps_of   <- function(size) size * S * 72 / res

col_bg    <- "#12262e"
col_land  <- "#1b3b46"
col_ridge <- "#2f5d6b"
col_a     <- "#5fb0cc"
col_b     <- "#e0a458"
col_cut   <- grDevices::adjustcolor("#4a6b78", alpha.f = 0.55)
col_word  <- "#f2f7f9"

hex_x <- c(50, 100, 100, 50,  0,   0)
hex_y <- c( 0, 28.75, 86.25, 115, 86.25, 28.75)

ridge_x <- c(0, 12, 22, 33, 44, 55, 66, 78, 89, 100)
ridge_y <- c(31, 25, 29, 19, 26, 17, 25, 20, 27,  30)

node_a <- data.frame(x = c(26, 38, 22, 30), y = c(32, 36, 40, 44))
node_b <- data.frame(x = c(58, 72, 52, 68), y = c(42, 36, 48, 52))

## edges given as index pairs into the node frames
edge_a <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 2), c(1, 3))
edge_b <- rbind(c(1, 2), c(1, 3), c(3, 4), c(4, 2), c(1, 4))

## the two edges the partition cuts, one endpoint in each regime
cut_from <- rbind(c(38, 36), c(30, 44))
cut_to   <- rbind(c(58, 42), c(52, 48))

series_x  <- c(14, 23, 32, 41, 50, 59, 68, 77, 86)
series_ya <- c(64, 59, 62, 58, 63, 59, 65, 60, 62)
series_yb <- c(70, 73, 68, 72, 67, 71, 69, 73, 69)

word          <- "Stem"
word_size     <- 26                 # user units, matching font-size in the SVG
word_baseline <- 95
word_target   <- 60                 # units, matching textLength in the SVG

disc <- function(x, y, r, col, n = 64) {
  a <- seq(0, 2 * pi, length.out = n + 1)
  graphics::polygon(x + r * cos(a), y + r * sin(a), col = col, border = NA)
}

segs <- function(nodes, edges, col, lwd) {
  graphics::segments(nodes$x[edges[, 1]], nodes$y[edges[, 1]],
                     nodes$x[edges[, 2]], nodes$y[edges[, 2]],
                     col = col, lwd = lwd, lend = "round")
}

## Georgia where it exists, a generic serif otherwise.
fam <- "serif"
if (.Platform$OS.type == "windows") {
  ok <- tryCatch({
    grDevices::windowsFonts(stemlogo = grDevices::windowsFont("Georgia"))
    TRUE
  }, error = function(e) FALSE)
  if (ok) fam <- "stemlogo"
}

dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
grDevices::png(out, width = width, height = height, res = res,
               bg = "transparent", type = "windows", antialias = "cleartype",
               pointsize = ps_of(word_size))

graphics::par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
graphics::plot.new()
graphics::plot.window(xlim = c(0, 100), ylim = c(115, 0))

graphics::polygon(hex_x, hex_y, col = col_bg, border = NA)

## space: the valley floor, framed by the reliefs
graphics::polygon(c(ridge_x, 100, 50, 0), c(ridge_y, 86.25, 115, 86.25),
                  col = col_land, border = NA)
graphics::lines(ridge_x, ridge_y, col = col_ridge, lwd = lwd_of(1),
                ljoin = "round")

graphics::segments(cut_from[, 1], cut_from[, 2], cut_to[, 1], cut_to[, 2],
                   col = col_cut, lwd = lwd_of(1))

segs(node_a, edge_a, col_a, lwd_of(1.7))
segs(node_b, edge_b, col_b, lwd_of(1.7))

for (i in seq_len(nrow(node_a))) disc(node_a$x[i], node_a$y[i], 3.2, col_a)
for (i in seq_len(nrow(node_b))) disc(node_b$x[i], node_b$y[i], 3.2, col_b)

## time: one series per regime
graphics::lines(series_x, series_ya, col = col_a, lwd = lwd_of(1.8),
                ljoin = "round", lend = "round")
graphics::lines(series_x, series_yb, col = col_b, lwd = lwd_of(1.8),
                ljoin = "round", lend = "round")

## The SVG forces the wordmark to word_target units through textLength, so that
## it fits the hexagon whatever the font. Here the equivalent is to measure the
## string and shrink it if the face is wider than the one it was laid out for.
cex <- 1
w   <- graphics::strwidth(word, cex = cex, font = 2, family = fam)
if (w > word_target) cex <- cex * word_target / w
graphics::text(50, word_baseline, word, adj = c(0.5, 0), cex = cex,
               col = col_word, font = 2, family = fam)

graphics::polygon(hex_x, hex_y, border = col_a, lwd = lwd_of(2.5),
                  ljoin = "round")

invisible(grDevices::dev.off())
cat("wrote ", out, " (", width, "x", height, ", wordmark ",
    round(min(w, word_target), 1), " units)\n", sep = "")
