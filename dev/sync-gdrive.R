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
## Never touch the RStudio project state, anything the user keeps there by hand,
## or the results of a run launched from the mirror. 00-setup.R diverts those to
## a sibling folder precisely so they cannot land here, but a run started with
## STEM_CACHE pointing inside the mirror would still put them here, and deleting
## somebody's simulation output to keep a mirror tidy is not a trade worth making.
protect <- grepl(paste0("^\\.Rproj\\.user/|^\\.Rhistory$|^\\.RData$|^dev/paper/cache/|",
                        "^dev/replication/(results|output|application)/"),
                 existing)
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


## ===========================================================================
## THE REPLICATION MATERIAL
##
## dev/replication is also delivered to a folder of its own on the Drive,
## BESIDE the mirror rather than inside it, because it is used differently:
## the scripts there are run, from the other machine too, and they write their
## results beside themselves. So this second copy follows opposite rules.
##
##   * nothing is ever deleted there: results, logs and whatever else a run
##     leaves are the point of the folder;
##   * a script is updated from the repository only if nobody has edited it
##     there since the last sync. The setup block of run-simulations.R is meant
##     to be edited, and an edit made on the other machine must not be silently
##     overwritten. When both sides changed, the repository version is written
##     beside it as <name>.from-repo and the conflict is reported.
##
## What was delivered last time is recorded in .synced, one md5 per script,
## which is how an edit made on the Drive is told apart from a stale copy.
## ===========================================================================
repl_src <- file.path(repo, "dev", "replication")
repl_dst <- if (length(args) >= 2) args[2] else
  file.path(dirname(dest), "SC-STEM-replication")
if (dir.exists(repl_src)) {
  if (!dir.exists(repl_dst)) dir.create(repl_dst, recursive = TRUE)
  manifest <- file.path(repl_dst, ".synced")
  last <- if (file.exists(manifest)) {
    m <- utils::read.csv(manifest, stringsAsFactors = FALSE)
    stats::setNames(m$md5, m$file)
  } else character(0)

  files <- tracked[startsWith(tracked, "dev/replication/")]
  files <- sub("^dev/replication/", "", files)
  n_new <- 0L; conflicts <- character(0); kept <- character(0)
  now <- character(0)
  for (f in files) {
    src <- file.path(repl_src, f)
    tgt <- file.path(repl_dst, f)
    if (!dir.exists(dirname(tgt))) dir.create(dirname(tgt), recursive = TRUE)
    md5_src <- unname(tools::md5sum(src))
    if (!file.exists(tgt)) {
      file.copy(src, tgt, copy.date = TRUE); n_new <- n_new + 1L
      now[f] <- md5_src; next
    }
    md5_tgt  <- unname(tools::md5sum(tgt))
    md5_last <- if (is.na(last[f])) NA_character_ else unname(last[f])
    if (identical(md5_tgt, md5_src)) { now[f] <- md5_src; next }
    drive_edited <- is.na(md5_last) || !identical(md5_tgt, md5_last)
    repo_changed <- is.na(md5_last) || !identical(md5_src, md5_last)
    if (!drive_edited) {
      ## untouched on the Drive: the repository version simply replaces it
      file.copy(src, tgt, overwrite = TRUE, copy.date = TRUE)
      n_new <- n_new + 1L; now[f] <- md5_src
    } else if (!repo_changed) {
      ## edited on the Drive, nothing new in the repository: leave it alone
      kept <- c(kept, f); now[f] <- md5_last
    } else {
      ## both sides changed, or a file of unknown origin: never overwrite
      file.copy(src, paste0(tgt, ".from-repo"), overwrite = TRUE, copy.date = TRUE)
      conflicts <- c(conflicts, f)
      now[f] <- if (is.na(md5_last)) md5_src else md5_last
    }
  }
  utils::write.csv(data.frame(file = names(now), md5 = unname(now)),
                   manifest, row.names = FALSE)
  message("\nreplication material: ", repl_dst)
  message("updated ", n_new, " of ", length(files), " script(s); nothing deleted")
  if (length(kept))
    message("edited on the Drive and unchanged in the repository, kept: ",
            paste(kept, collapse = ", "))
  if (length(conflicts))
    message("CONFLICT -- changed both on the Drive and in the repository, left as ",
            "they are on the Drive: ", paste(conflicts, collapse = ", "),
            "\n  the repository version is beside each as <name>.from-repo")
}
