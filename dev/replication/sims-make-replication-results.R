## ===========================================================================
## SC-STEM simulation study: THE RESULTS OF THE REPLICATION PACKAGE
##
## Writes the results of the simulation study shipped with the paper
## (papers/sc-stem-local-leader/results/simulations) from the files written by
## the runs: for every stream, the replications (<tag>_s<stream>.csv), the
## parameters (-params) and the grids (-grid), compressed (.csv.gz), which is
## what the analysis (code/simulations/03-results.R) reads. The per-location
## records (-stations) and the full observations (-obs) are left out: the
## analysis does not use them, and any replication can be drawn again from its
## seed.
##
## Lines broken by the synchronization client while the runs appended to the
## files are dropped here, with the same rule as the analysis (a line must
## start with the quoted cell and have as many fields as the header); the
## replications left with missing rows are dropped by the analysis.
##
##     Rscript sims-make-replication-results.R --results=<the results folder> --out=<...>
##
## This is part of the development material, not of the replication package.
## ===========================================================================
HERE <- local({
  a <- commandArgs(trailingOnly = FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]), winslash = "/")) else getwd()
})
args <- commandArgs(trailingOnly = TRUE)
opt <- function(k, d) { a <- grep(paste0("^--", k, "="), args, value = TRUE); if (length(a)) sub("^--[^=]+=", "", a[1]) else d }
RES <- normalizePath(opt("results", file.path(HERE, "results")), winslash = "/", mustWork = TRUE)
OUT <- normalizePath(opt("out", file.path(HERE, "..", "..", "papers", "sc-stem-local-leader", "results", "simulations")),
                     winslash = "/", mustWork = FALSE)
TAG <- opt("tag", "main3")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

files <- list.files(RES, pattern = sprintf("^%s_s[0-9]+(-params|-grid|-timeouts)?[.]csv$", TAG))
for (f in files) {
  L <- readLines(file.path(RES, f), warn = FALSE)
  n_quote <- nchar(L) - nchar(gsub("\"", "", L, fixed = TRUE))
  bare <- gsub("\"[^\"]*\"", "", L)
  n_field <- nchar(bare) - nchar(gsub(",", "", bare, fixed = TRUE)) + 1L
  ok <- n_quote %% 2L == 0L & n_field == n_field[1] & startsWith(L, "\"")
  ok[1] <- TRUE
  con <- gzfile(file.path(OUT, paste0(f, ".gz")), "w", compression = 9)
  writeLines(L[ok], con); close(con)
  cat(sprintf("%-26s %8d lines, %d broken dropped -> %.1f MB\n", f, length(L), sum(!ok),
              file.size(file.path(OUT, paste0(f, ".gz"))) / 1e6))
}
