## ---------------------------------------------------------------------------
## Mirror the package sources onto the Google Drive working folder.
##
## The git repository is the single source of truth; the Google Drive copy is a
## mirror kept for offline access and for continuity with the historical
## material in that folder. Run this after every change that is committed, so
## that the two never diverge:
##
##     Rscript dev/sync-gdrive.R
##
## Only files TRACKED BY GIT are copied, so build artefacts, .Rcheck
## directories, tarballs and RStudio state never leak into the mirror. Files
## that exist in the mirror but no longer in the repository are deleted, which
## is what makes this a mirror rather than a copy.
##
## This script is excluded from the built package through .Rbuildignore.
## ---------------------------------------------------------------------------

repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!dir.exists(file.path(repo, ".git"))) {
  stop("Run this script from the root of the Stem repository:\n",
       "    Rscript dev/sync-gdrive.R", call. = FALSE)
}

## Destination can be overridden, e.g. when the drive letter changes:
##     Rscript dev/sync-gdrive.R "D:/some/other/folder/Stem"
args <- commandArgs(trailingOnly = TRUE)
dest <- if (length(args)) args[1] else
  "K:/Il mio Drive/Spatiotemporal_Modeling/STEM_Cameletti/Stem"

message("source : ", repo)
message("mirror : ", dest)

## --- files tracked by git ---------------------------------------------------
tracked <- system2("git", c("-C", shQuote(repo), "ls-files"),
                   stdout = TRUE, stderr = FALSE)
if (!length(tracked)) stop("git ls-files returned nothing.", call. = FALSE)

## --- refuse to run on a dirty tree ------------------------------------------
## Mirroring uncommitted work would put the two copies out of step with the
## history, which is exactly what this script exists to prevent.
dirty <- system2("git", c("-C", shQuote(repo), "status", "--porcelain"),
                 stdout = TRUE, stderr = FALSE)
if (length(dirty)) {
  message("\nUncommitted changes:\n", paste(" ", dirty, collapse = "\n"))
  stop("Commit (or stash) before syncing, so the mirror matches a real commit.",
       call. = FALSE)
}

## --- copy -------------------------------------------------------------------
if (!dir.exists(dest)) dir.create(dest, recursive = TRUE)

copied <- 0L
for (f in tracked) {
  src <- file.path(repo, f)
  tgt <- file.path(dest, f)
  if (!file.exists(src)) next
  d <- dirname(tgt)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  ## copy only when the content actually differs, to spare the sync client
  same <- file.exists(tgt) &&
    identical(tools::md5sum(src)[[1]], tools::md5sum(tgt)[[1]])
  if (!same) {
    file.copy(src, tgt, overwrite = TRUE, copy.date = TRUE)
    copied <- copied + 1L
  }
}

## --- remove what is no longer tracked ---------------------------------------
existing <- list.files(dest, recursive = TRUE, all.files = TRUE,
                       no.. = TRUE, full.names = FALSE)
## never touch the RStudio project state or anything the user keeps there by hand
protect <- grepl("^\\.Rproj\\.user/|^\\.Rhistory$|^\\.RData$", existing)
stale <- setdiff(existing[!protect], tracked)
if (length(stale)) {
  file.remove(file.path(dest, stale))
  message("removed ", length(stale), " stale file(s)")
}

## --- record which commit the mirror corresponds to ---------------------------
sha <- system2("git", c("-C", shQuote(repo), "rev-parse", "--short", "HEAD"),
               stdout = TRUE, stderr = FALSE)
writeLines(c(
  "This folder is a MIRROR of the git repository and must not be edited.",
  "",
  paste0("Repository : ", repo),
  paste0("Remote     : https://github.com/PaoloMaranzano/Stem"),
  paste0("Commit     : ", sha),
  paste0("Synced     : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "Regenerate with:  Rscript dev/sync-gdrive.R"
), file.path(dest, "MIRROR.txt"))

message("copied ", copied, " changed file(s) of ", length(tracked), " tracked")
message("mirror now at commit ", sha)
