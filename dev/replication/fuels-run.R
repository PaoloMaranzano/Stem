## ===========================================================================
## SC-STEM fuel application: THE LAUNCHER
##
## Runs the main scripts of the application (fuels-main-setting1.R, -2.R and
## -3.R) on the settings, the cities and the fuels chosen in THE SETUP below,
## every run in an R process of its own, so that the main scripts behave
## exactly as from a terminal. Edit THE SETUP, then press Source in RStudio,
## or from a terminal
##
##     Rscript fuels-run.R
##
## With processes = 1 the runs follow one another and their output appears in
## the console. With processes > 1 they run at the same time: the settings are
## shared among the processes, and a setting given more than one process is
## split by model (the models of its cities and fuels are dealt out to the
## processes, the largest first); each process writes its log in
## <root>/logs/, and the console shows, every few minutes, the last line of
## every log. The seed of a model depends on its city and fuel only, so neither
## the processes nor the machines change the results.
##
## When the runs are over, the launcher rebuilds the tables of every setting
## from the files of all its models (stage "summary"). With several machines,
## give each machine its own cities, never the same city and fuel to two
## machines; once all have finished, copy the folders <root>/<case>/<xdef>/
## of each machine into one, and run this launcher there with stage =
## "summary".
##
## Requirements: the data of the cases, <root>/<case>/data.RData, made by
## fuels-data.R; the package callr (installed when missing); Stem at the
## commit the scripts are pinned to (installed when missing, before any run
## starts, so that the processes do not install it at the same time).
##
## This is part of the REPLICATION MATERIAL of the paper, not of the Stem
## package.
## ===========================================================================


## ===========================================================================
## THE SETUP. This is the part meant to be edited.
## ===========================================================================
RUN <- list(
  ## the settings: 1 = Y a major, X its independent neighbour; 2 = Y
  ## independent, X its major neighbour; 3 = Y a major, X its neighbour of
  ## another major brand. One, several, or 1:3 for all
  settings   = 1:3,
  ## the cities: "all", or codes among
  ##   RM Roma, MI Milano, NA Napoli, TO Torino, PA Palermo, BA Bari,
  ##   CT Catania, BO Bologna, ME Messina, FI Firenze, VE Venezia
  ## for instance c("RM", "MI")
  cities     = "all",
  ## the fuel of Y: "g" gasoline, "d" diesel
  fuels      = c("g", "d"),
  ## the stage: "all" (grid, lambda and final in a row), "grid", "lambda",
  ## "final", "summary" (the tables only, nothing fitted) or "dry" (the plan)
  stage      = "all",
  ## R processes at the same time on this machine
  processes  = 1L,
  ## processes of every bootstrap (keep processes * boot_cores within the cores
  ## of the machine)
  boot_cores = 7L,
  ## the neighbour X: "nn" the nearest pump of the X set (the main analysis)
  xdef       = "nn",
  ## the folder with the data and the outputs: "" for fuels/ beside the
  ## scripts. A local folder is safer than a synchronized one (Google Drive):
  ## copy the data there first, and the outputs back at the end
  root       = "",
  ## minutes between two reports of the progress, with processes > 1
  report_min = 10
)


## ===========================================================================
## Machinery: nothing below needs editing
## ===========================================================================
fr_here <- local({
  a <- commandArgs(trailingOnly = FALSE); m <- grep("^--file=", a, value = TRUE)
  f <- if (length(m)) normalizePath(sub("^--file=", "", m[1]), winslash = "/") else NULL
  if (is.null(f)) for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) { f <- normalizePath(of, winslash = "/"); break }
  }
  if (is.null(f)) normalizePath(getwd(), winslash = "/") else dirname(f)
})
if (!requireNamespace("callr", quietly = TRUE)) {
  message("installing callr ...")
  utils::install.packages("callr", repos = "https://cloud.r-project.org")
}

fr_run <- function(run = RUN, here = fr_here) {
  cities_all <- c("RM", "MI", "NA", "TO", "PA", "BA", "CT", "BO", "ME", "FI", "VE")
  settings <- sort(unique(as.integer(run$settings)))
  if (!length(settings) || any(!settings %in% 1:3))
    stop("settings must be among 1, 2 and 3", call. = FALSE)
  cities <- if (identical(run$cities, "all")) character(0) else toupper(run$cities)
  bad <- setdiff(cities, cities_all)
  if (length(bad)) stop("unknown cities: ", paste(bad, collapse = ", "), call. = FALSE)
  fuels <- intersect(c("g", "d"), run$fuels)
  if (!length(fuels)) stop('fuels must be among "g" and "d"', call. = FALSE)
  stage <- match.arg(run$stage, c("all", "grid", "lambda", "final", "summary", "dry"))
  processes <- max(1L, as.integer(run$processes))
  root <- if (nzchar(run$root)) normalizePath(run$root, winslash = "/", mustWork = FALSE) else
    file.path(here, "fuels")
  scripts <- file.path(here, sprintf("fuels-main-setting%d.R", settings))
  if (any(!file.exists(scripts)))
    stop("missing beside the launcher: ", paste(basename(scripts[!file.exists(scripts)]), collapse = ", "),
         call. = FALSE)
  cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (is.finite(cores) && processes * run$boot_cores > cores)
    message(sprintf("note: %d processes x %d bootstrap processes exceed the %d cores of this machine",
                    processes, run$boot_cores, cores))

  ## the runs: every setting gets an equal share of the processes, and a
  ## setting with more than one is split into jobs "i/N"
  share <- rep(processes %/% length(settings), length(settings))
  extra <- processes %% length(settings)
  if (extra) share[seq_len(extra)] <- share[seq_len(extra)] + 1L
  share <- pmax(1L, share)
  if (stage %in% c("summary", "dry")) share[] <- 1L
  tasks <- do.call(rbind, lapply(seq_along(settings), function(i)
    data.frame(setting = settings[i], job = sprintf("%d/%d", seq_len(share[i]), share[i]),
               stringsAsFactors = FALSE)))
  args_of <- function(tk, stage_tk = stage) c(
    paste0("--stage=", stage_tk), paste0("--xdef=", run$xdef), paste0("--fuels=", paste(fuels, collapse = ",")),
    if (length(cities)) paste0("--cities=", paste(cities, collapse = ",")),
    paste0("--boot_cores=", as.integer(run$boot_cores)), paste0("--job=", tk$job),
    if (nzchar(run$root)) paste0("--root=", root))

  cat(sprintf("fuel application: settings %s; cities %s; fuels %s; stage %s\n",
              paste(settings, collapse = ", "), if (length(cities)) paste(cities, collapse = ", ") else "all",
              paste(fuels, collapse = ", "), stage))
  cat(sprintf("  %d run(s) on %d process(es) at a time, %d processes per bootstrap; outputs in %s\n\n",
              nrow(tasks), processes, as.integer(run$boot_cores), root))

  ## Stem at the pinned commit, once, in a separate process
  if (!stage %in% c("summary", "dry")) {
    sha <- callr::r(function(f) { source(f); fu_require_stem() },
                    args = list(file.path(here, "fuels-functions.R")), show = TRUE)
    cat("Stem at commit", sha, "\n\n")
  }

  logs <- file.path(root, "logs"); dir.create(logs, recursive = TRUE, showWarnings = FALSE)
  stamp <- format(Sys.time(), "%Y%m%d-%H%M")
  log_of <- function(tk) file.path(logs, sprintf("setting%d_%s_job%s_%s.log", tk$setting, stage,
                                                 sub("/", "of", tk$job), stamp))
  script_of <- function(tk) file.path(here, sprintf("fuels-main-setting%d.R", tk$setting))
  t0 <- Sys.time()

  if (processes == 1L || nrow(tasks) == 1L) {
    ## one run after the other, with the output in the console (and in the log)
    for (k in seq_len(nrow(tasks))) {
      tk <- tasks[k, ]
      cat(sprintf("=== setting %d, job %s, started %s ===\n", tk$setting, tk$job, format(Sys.time(), "%H:%M")))
      res <- callr::rscript(script_of(tk), cmdargs = args_of(tk), wd = here, show = TRUE,
                            stdout = log_of(tk), stderr = "2>&1", fail_on_status = FALSE)
      cat(sprintf("=== setting %d, job %s: %s after %.0f min ===\n\n", tk$setting, tk$job,
                  if (identical(res$status, 0L)) "done" else paste("FAILED, status", res$status),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    }
  } else {
    ## several runs at the same time, each with its log
    queue <- seq_len(nrow(tasks)); running <- list(); done <- character(0)
    last_report <- Sys.time() - 3600
    repeat {
      while (length(queue) && length(running) < processes) {
        k <- queue[1]; queue <- queue[-1]; tk <- tasks[k, ]
        opts <- callr::rscript_process_options(script = script_of(tk), cmdargs = args_of(tk), wd = here,
                                               stdout = log_of(tk), stderr = "2>&1")
        running[[as.character(k)]] <- callr::rscript_process$new(opts)
        cat(sprintf("%s  started setting %d, job %s -> %s\n", format(Sys.time(), "%H:%M"), tk$setting,
                    tk$job, basename(log_of(tk))))
      }
      for (key in names(running)) {
        p <- running[[key]]
        if (!p$is_alive()) {
          tk <- tasks[as.integer(key), ]
          st <- p$get_exit_status()
          cat(sprintf("%s  setting %d, job %s: %s\n", format(Sys.time(), "%H:%M"), tk$setting, tk$job,
                      if (identical(st, 0L)) "done" else paste("FAILED, status", st, "(see the log)")))
          running[[key]] <- NULL; done <- c(done, key)
        }
      }
      if (!length(running) && !length(queue)) break
      if (as.numeric(difftime(Sys.time(), last_report, units = "mins")) >= run$report_min) {
        cat(sprintf("--- %s, %.0f min ---\n", format(Sys.time(), "%H:%M"),
                    as.numeric(difftime(Sys.time(), t0, units = "mins"))))
        for (key in names(running)) {
          tk <- tasks[as.integer(key), ]
          l <- tryCatch(readLines(log_of(tk), warn = FALSE), error = function(e) character(0))
          l <- l[nzchar(trimws(l))]
          cat(sprintf("  setting %d, job %s: %s\n", tk$setting, tk$job,
                      if (length(l)) substr(l[length(l)], 1, 150) else "(no output yet)"))
        }
        last_report <- Sys.time()
      }
      Sys.sleep(20)
    }
  }

  ## the tables of every setting, from the files of all its models
  if (stage %in% c("all", "grid", "lambda", "final") && (processes > 1L || length(settings) > 0L)) {
    cat("\nthe tables of every setting (stage summary)\n")
    for (s in settings) {
      tk <- data.frame(setting = s, job = "1/1", stringsAsFactors = FALSE)
      callr::rscript(script_of(tk), cmdargs = args_of(tk, "summary"), wd = here, show = TRUE,
                     stderr = "2>&1", fail_on_status = FALSE)
    }
  }
  cat(sprintf("\nall done in %.1f hours; logs in %s\n",
              as.numeric(difftime(Sys.time(), t0, units = "hours")), logs))
  invisible(tasks)
}

fr_run()
