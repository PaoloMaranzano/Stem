## ---------------------------------------------------------------------------
## Shared setup for the paper scripts.
##
## Everything here exists so that the scripts run from ANY working directory and
## on any machine: the repository root is resolved from the location of the file
## being executed, not from getwd(), and every output directory is derived from
## it or read from an environment variable. Nothing below hard-codes a path.
##
## Source it at the top of a script:
##
##     source(file.path(dirname(sys.frame(1)$ofile), "00-setup.R"))
##
## or, more simply and equivalently from a script run with Rscript,
##
##     source("00-setup.R")            # when the working directory is dev/paper
##     Rscript dev/paper/04-cv.R           # from anywhere: the script finds itself
##
## Environment variables, all optional:
##
##   STEM_ROOT     the repository root, when the automatic search fails
##   STEM_CACHE    where the result CSVs go       (default <root>/dev/paper/cache)
##   STEM_FIGURES  where the figures and tables go (default the Overleaf folder
##                 if it exists, otherwise <root>/dev/paper/figures)
## ---------------------------------------------------------------------------

## The path of the file currently being executed, whether by Rscript, by
## source(), or inside RStudio.
stem_this_file <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[1]), winslash = "/"))
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(normalizePath(of, winslash = "/"))
  }
  NULL
}

## The root of the repository: the nearest ancestor holding a DESCRIPTION file.
stem_root <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    env <- Sys.getenv("STEM_ROOT", "")
    if (nzchar(env) && file.exists(file.path(env, "DESCRIPTION"))) {
      cached <<- normalizePath(env, winslash = "/"); return(cached)
    }
    start <- stem_this_file()
    p <- if (is.null(start)) normalizePath(getwd(), winslash = "/") else dirname(start)
    repeat {
      if (file.exists(file.path(p, "DESCRIPTION"))) { cached <<- p; return(p) }
      up <- dirname(p)
      if (identical(up, p)) break
      p <- up
    }
    stop("cannot locate the repository root: set STEM_ROOT", call. = FALSE)
  }
})

stem_paper_dir <- function() file.path(stem_root(), "dev", "paper")

## ---------------------------------------------------------------------------
## Running from the Google Drive mirror
##
## The mirror is regenerated from the repository by dev/sync-gdrive.R, which
## DELETES whatever it finds there and does not track: that is what makes it a
## mirror rather than a copy. Results written inside it would therefore survive
## only until the next sync. So when the scripts notice they are running from a
## mirror -- a root carrying a MIRROR.txt -- the outputs are diverted to a
## sibling folder that the sync never touches:
##
##     .../STEM_Cameletti/Stem/         the mirror, code only, rewritten
##     .../STEM_Cameletti/Stem-runs/    the results, never touched
##
## STEM_CACHE and STEM_FIGURES always win, so a run can still be pointed
## anywhere. Nothing changes on a machine that works in the repository itself.
## ---------------------------------------------------------------------------
stem_is_mirror <- function() file.exists(file.path(stem_root(), "MIRROR.txt"))

stem_runs_dir <- function() {
  file.path(dirname(stem_root()), paste0(basename(stem_root()), "-runs"))
}

stem_out_dir <- local({
  warned <- FALSE
  function(env, inside, leaf) {
    p <- Sys.getenv(env, "")
    if (!nzchar(p)) {
      if (stem_is_mirror()) {
        p <- file.path(stem_runs_dir(), leaf)
        if (!warned) {
          warned <<- TRUE
          message("running from the mirror: results go to ", stem_runs_dir(),
                  "\n  (the mirror itself is rewritten by dev/sync-gdrive.R)")
        }
      } else {
        p <- inside
      }
    }
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
    normalizePath(p, winslash = "/")
  }
})

stem_cache_dir <- function()
  stem_out_dir("STEM_CACHE", file.path(stem_paper_dir(), "cache"), "cache")

## The figure directory. The Overleaf folder when it is there, so that a run on
## the author's machine writes straight into the paper; a folder inside the
## repository otherwise, so that a run on a virtual machine still produces
## something.
stem_fig_dir <- function() {
  ov <- file.path(Sys.getenv("USERPROFILE", Sys.getenv("HOME")), "Dropbox",
                  "Applicazioni", "Overleaf", "SC-STEM package paper", "Figures")
  inside <- if (dir.exists(ov)) ov else file.path(stem_paper_dir(), "figures")
  if (dir.exists(ov) && !nzchar(Sys.getenv("STEM_FIGURES", ""))) {
    ## Overleaf is on this machine: write straight into the paper, mirror or not
    dir.create(ov, recursive = TRUE, showWarnings = FALSE)
    return(normalizePath(ov, winslash = "/"))
  }
  stem_out_dir("STEM_FIGURES", inside, "figures")
}

## Load the package from source when pkgload is available, from the library
## otherwise. A virtual machine that has the package installed does not need the
## development toolchain.
stem_load <- function() {
  if (requireNamespace("pkgload", quietly = TRUE)) {
    suppressMessages(pkgload::load_all(stem_root(), quiet = TRUE))
  } else {
    library("Stem", character.only = TRUE)
  }
  invisible(TRUE)
}

## ---------------------------------------------------------------------------
## Command-line overrides
##
## Every element of a configuration list can be overridden with --name=value on
## the command line; a value is split on commas and coerced to the type of the
## default, so
##
##     Rscript some-script.R --n=50,100 --nrep=25
##
## sets those two options, whatever the defaults were, and leaves the rest.
## Passing --name= with nothing after the equals sign keeps the default.
##
## The simulation study and the application no longer use this file: they are
## the standalone scripts of dev/replication, which need nothing but an
## installed Stem. What remains here are the development diagnostics.
## ---------------------------------------------------------------------------
stem_config <- function(defaults, args = commandArgs(trailingOnly = TRUE)) {
  kv <- grep("^--[^=]+=", args, value = TRUE)
  for (a in kv) {
    nm  <- sub("^--([^=]+)=.*$", "\\1", a)
    val <- sub("^--[^=]+=", "", a)
    if (!nzchar(val)) next
    if (!nm %in% names(defaults)) {
      stop("unknown option --", nm, "; the options are: ",
           paste(names(defaults), collapse = ", "), call. = FALSE)
    }
    parts <- trimws(strsplit(val, ",", fixed = TRUE)[[1]])
    d <- defaults[[nm]]
    defaults[[nm]] <-
      if (is.logical(d))   as.logical(parts)
      else if (is.integer(d)) as.integer(parts)
      else if (is.numeric(d)) {
        ## allow fractions such as 2/3 on the command line
        vapply(parts, function(s) eval(parse(text = s)), numeric(1),
               USE.NAMES = FALSE)
      }
      else parts
    if (anyNA(defaults[[nm]])) {
      stop("cannot read --", nm, "=", val, call. = FALSE)
    }
  }
  defaults
}

stem_print_config <- function(cfg) {
  cat("configuration\n")
  for (nm in names(cfg)) {
    v <- cfg[[nm]]
    cat(sprintf("  %-10s %s\n", nm,
                paste(if (is.numeric(v)) signif(v, 4) else v, collapse = ", ")))
  }
  cat("\n")
}
