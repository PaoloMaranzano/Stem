## ===========================================================================
## The PDFs of the SC-STEM paper in the repository
##
## Compiles the manuscript and its supplements from the Overleaf project
## (synchronized in Dropbox) and copies the PDFs into
## papers/sc-stem-local-leader/manuscript. Only the PDFs go into the
## repository, never the sources.
##
## The PDFs are clean: the compilation runs on a copy of the project in a
## temporary folder, where the notes (\nota{}) are dropped and the changes
## marked in orange (\rev{} and the environment "modified") are printed in
## black. The project itself is not touched. Text still marked in red is
## reported, since it would appear in the PDFs.
##
##     Rscript dev/sync-paper.R
##     Rscript dev/sync-paper.R --src=<the Overleaf folder>
##
## Requires a TeX distribution with pdflatex and bibtex on the path (MiKTeX
## or TeX Live). Run it at the milestones of the paper, not at every edit.
## ===========================================================================

HERE <- local({
  a <- commandArgs(trailingOnly = FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]), winslash = "/")) else getwd()
})
args <- commandArgs(trailingOnly = TRUE)
opt <- function(k, d) { a <- grep(paste0("^--", k, "="), args, value = TRUE); if (length(a)) sub("^--[^=]+=", "", a[1]) else d }
SRC <- normalizePath(opt("src", "C:/Users/paulm/Dropbox/Applicazioni/Overleaf/SC-STEM package paper"),
                     winslash = "/", mustWork = TRUE)
OUT <- file.path(HERE, "..", "papers", "sc-stem-local-leader", "manuscript")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## the documents, and the name of each PDF in the repository
DOCS <- c("main"                   = "sc-stem-local-leader.pdf",
          "supplement"             = "supplement-A-estimation.pdf",
          "supplement-simulations" = "supplement-B-simulations.pdf",
          "supplement-application" = "supplement-C-application.pdf",
          "code-supplement"        = "code-supplement.pdf")

## a copy of the project, with the draft marks switched off just before
## \begin{document} of every document
TMP <- file.path(tempdir(), "sc-stem-paper")
unlink(TMP, recursive = TRUE); dir.create(TMP)
invisible(file.copy(list.files(SRC, full.names = TRUE, all.files = FALSE), TMP, recursive = TRUE))
CLEAN <- c("\\providecommand{\\rev}[1]{}\\renewcommand{\\rev}[1]{#1}",
           "\\providecommand{\\nota}[1]{}\\renewcommand{\\nota}[1]{}",
           "\\makeatletter\\@ifundefined{modified}{\\newenvironment{modified}{}{}}{\\renewenvironment{modified}{}{}}\\makeatother")
for (d in names(DOCS)) {
  f <- file.path(TMP, paste0(d, ".tex")); x <- readLines(f, warn = FALSE)
  i <- grep("^\\\\begin\\{document\\}", x)[1]
  writeLines(append(x, CLEAN, after = i - 1L), f)
  ## text still in red would appear in the PDF: report it
  red <- grep("\\\\textcolor\\{red\\}", x, value = TRUE)
  red <- red[!grepl("^\\s*(%|\\\\newcommand|\\\\providecommand|\\\\renewcommand)", red)]
  if (length(red)) message("red text in ", d, ".tex:\n  ", paste(substr(trimws(red), 1, 120), collapse = "\n  "))
}

## pdflatex, bibtex, pdflatex twice; MiKTeX installs nothing on the fly
tex <- function(cmd, doc) {
  a <- if (cmd == "pdflatex") c("-interaction=nonstopmode", "-halt-on-error",
                                if (.Platform$OS.type == "windows") "--disable-installer", doc) else doc
  system2(cmd, a, stdout = FALSE, stderr = FALSE)
}
owd <- setwd(TMP); on.exit(setwd(owd), add = TRUE)
for (d in names(DOCS)) {
  tex("pdflatex", d); tex("bibtex", d); tex("pdflatex", d); s <- tex("pdflatex", d)
  pdf <- paste0(d, ".pdf")
  if (s != 0 || !file.exists(pdf)) { message("FAILED: ", d, " (see ", file.path(TMP, paste0(d, ".log")), ")"); next }
  file.copy(pdf, file.path(OUT, DOCS[[d]]), overwrite = TRUE)
  message("wrote ", DOCS[[d]])
}
