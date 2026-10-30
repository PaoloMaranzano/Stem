## ===========================================================================
## Latent Spatial Regimes in Spatio-Temporal State-Space Models: Evidence on
## Local Price Leadership in Italian Fuel Markets
##
## Author:  Paolo Maranzano (University of Milano-Bicocca; Fondazione Eni Enrico
##          Mattei)
## Script:  Sims_01_Run.R, the run of the simulation study (Section 4 of the
##          paper, Supplement B)
##
## What it does: runs the Monte Carlo study: for every cell of the design
## (Sims_00_Setup.R) and every replication, one data set is drawn and the
## whole procedure is run on it (sim_one() in Sims_00_Setup.R). No table or
## figure of the paper is written here: Sims_02_Analysis.R draws them all from
## the results.
##
## Outputs, per stream, in <out> (output/runs/simulations by default),
## keyed on (cell, rep); rep is unique across the streams, so the outputs of
## all the streams stack. <tag> stands for <tag>_s<stream>:
##
##   <tag>.csv           one row per replication: the selected (K, phi), the
##                       recovery of the partition, the error of the fitted
##                       signal, the times
##   <tag>-params.csv    one row per replication, regime and parameter: the
##                       truth beside the estimate at the true number of
##                       regimes, after aligning the estimated regimes on the
##                       true ones, and whether the range of the regime ended
##                       at a limit of what its distances identify
##   <tag>-stations.csv  one row per replication and location
##   <tag>-grid.csv      one row per replication and (K, phi) of the fitted
##                       grid: log-likelihood, parameters, criteria, and the
##                       ARI of its partition with the truth
##   <tag>-obs/*.rds     the full observations of the first keep_obs
##                       replications of every cell
##   <tag>-timeouts.csv  one row per replication stopped by the time limit
## The results of the study of the paper (the first, second and fourth files,
## compressed) are in results/simulations.
## ===========================================================================
##
##     Rscript Sims_01_Run.R --stream=1 --nrep=100 --cores=7
##     Rscript Sims_01_Run.R --mode=dry                  the cells and the cost
##     Rscript Sims_01_Run.R --nrep=1 --cores=7          one replication of every cell
##
## or open it in RStudio, edit "THE RUN" below and press Source. The run can
## be stopped at any moment and resumed: what is recorded is skipped.
##
## The study of the paper: 332 cells, 200 replications per cell, run as two
## streams of 100 on two machines (--stream=1 and --stream=2, --nrep=100), on
## Stem at the commit pinned in Sims_00_Setup.R (SIM_STEM_REF). A stream draws
## its own data, with seeds no other stream uses: replication r of stream s is
## recorded as rep = 1000 (s - 1) + r, and its seed is 1000 rep + 1. Never run
## the same stream on two machines at once. The 66,400 replications took
## about 2,170 core-hours.
##
## THE TIME LIMIT. Every replication runs in a separate R process (package
## callr) under a time limit set from the times already recorded: cap_mult
## times the median of the cell once cap_after replications of it are done,
## before that cap_first_mult times the median of the cells with the same
## (n, T), and cap_first_min seconds when nothing is known yet. A replication
## over its limit is stopped, written to <tag>_s<stream>-timeouts.csv, and
## drawn again with the seed of the next attempt (1000 rep + 1 + 100 (attempt
## - 1)), up to cap_attempts; the last one is recorded as failed. cap_mult = 0
## removes the limit.
##
## Margins given on the command line (--n, --TN, --omega, --balance, --knn)
## replace the blocks with one factorial block of those margins around the
## reference, as in  Rscript Sims_01_Run.R --TN=60,120,365 --mode=dry
##
## Requirements: R (>= 4.2) with spdep and callr, and Stem at the commit pinned
## in Sims_00_Setup.R (installed from GitHub when another version is
## installed).
## ===========================================================================

## Definitions only? The worker processes, and Sims_02_Analysis.R, load this file
## with SIM_DEFINE_ONLY <- TRUE and SIM_FILE set to its path: they get the
## design, the generator and the functions, and nothing is run.
if (!exists("SIM_DEFINE_ONLY", inherits = FALSE)) SIM_DEFINE_ONLY <- FALSE
if (!exists("SIM_FILE", inherits = FALSE)) SIM_FILE <- local({
  a <- commandArgs(trailingOnly = FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[1]), winslash = "/"))
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(normalizePath(of, winslash = "/"))
  }
  NULL
})
SIM_HERE <- if (is.null(SIM_FILE)) normalizePath(getwd(), winslash = "/") else dirname(SIM_FILE)
source(file.path(SIM_HERE, "Sims_00_Setup.R"), local = globalenv())

## ---------------------------------------------------------------------------
## THE RUN. What it does.
## ---------------------------------------------------------------------------
CFG <- sim_config(args = if (SIM_DEFINE_ONLY) character(0) else commandArgs(TRUE), c(
  list(scenario = SIM_SCENARIOS),
  SIM_REFERENCE,
  list(
    ## the blocks this invocation runs, by name
    blocks   = names(SIM_BLOCKS),

    ## "run": the Monte Carlo; "dry": print the cells and the cost, run nothing
    mode     = "run",

    ## How many cores. Each core runs one replication at a time.
    cores    = sim_default_cores(),

    ## The stream: one per machine. Two machines that run the study at the
    ## same time give themselves different streams; each draws its own data
    ## (different seeds) and writes its own files, <tag>_s<stream>*, and the
    ## analysis stacks every stream of a tag. Copy the results/ files of one
    ## machine beside those of the other to put them together.
    stream   = 1L,

    ## replications per cell in this stream, and the slice this invocation
    ## covers (rep_to defaults to nrep). A stream can also be run in slices:
    ## the files are appended, so the slices compose. At most 999 per stream.
    nrep     = 10L,
    rep_from = 1L,
    rep_to   = NA_integer_,

    ## What the estimator searches over. phi_ref is the penalty at which the
    ## recovery at the true number of regimes is read, a point of phi_grid.
    K_grid   = 1:4,
    phi_grid = c(0, 0.025, 0.05, 0.1, 0.2, 0.5, 1),
    phi_ref  = 0.05,

    ## How many replications per cell keep their full per-observation record.
    ## One: every record is n x T rows, and five per cell would take gigabytes;
    ## any replication can be regenerated from its seed.
    keep_obs = 1L,

    ## The time limit of a replication, adaptive cell by cell. Once a cell has
    ## cap_after replications that ended normally, the limit is cap_mult times
    ## their median (at least cap_min seconds); before that, cap_first_mult
    ## times the median of the replications with the same (n, T) in any cell
    ## (at least cap_first_min seconds), or cap_first_min when there are none.
    ## A replication over its limit is stopped, written to <tag>-timeouts.csv
    ## and drawn again with the seed of its next attempt, up to cap_attempts
    ## attempts; after the last it is recorded as failed. Each replication then
    ## runs in a process of its own (package callr), which is what lets it be
    ## stopped. cap_mult = 0 switches the limit off.
    cap_mult       = 10,
    cap_min        = 600,
    cap_after      = 3L,
    cap_first_mult = 20,
    cap_first_min  = 3600,
    cap_attempts   = 3L,

    ## Bookkeeping. `out` defaults to output/runs/simulations in the
    ## folder of the paper.
    ## "main3": the design of October 2026 (new values of the strong level,
    ## the fourth geometry) with the estimator of the commit pinned above
    ## (SQUAREM, tolerances 1e-3, no regularization). "main2", the same design
    ## with the EM algorithm at regularization 0.01, and "main", the first pass
    ## of September, are not comparable with it.
    tag      = "main3",
    out      = file.path(SIM_HERE, "..", "output", "runs", "simulations")
  )))

## Margins given on the command line replace the blocks with one factorial
## block of those margins around the reference, as in
##   Rscript Sims_01_Run.R --TN=60,365 --n=50
SIM_CLI_MARGINS <- if (SIM_DEFINE_ONLY) character(0) else
  intersect(names(SIM_REFERENCE),
            sub("^--([^=]+)=.*$", "\\1", grep("^--[^=]+=.", commandArgs(TRUE), value = TRUE)))


## ===========================================================================
## Main
## ===========================================================================
sim_main <- function() {

  OUT <- normalizePath(CFG$out[1], winslash = "/", mustWork = FALSE)
  sim_print_config(CFG)

  stream <- CFG$stream[1]
  rep_to <- if (is.na(CFG$rep_to[1])) CFG$nrep[1] else CFG$rep_to[1]
  if (!(stream >= 1L)) stop("stream must be a positive integer", call. = FALSE)
  if (CFG$rep_from[1] < 1L || rep_to < CFG$rep_from[1] || rep_to > SIM_REP_MAX)
    stop("replications run from rep_from to rep_to, within 1 to ", SIM_REP_MAX,
         " in a stream", call. = FALSE)
  if (CFG$cap_attempts[1] > 9L) stop("cap_attempts must be at most 9", call. = FALSE)

  ## the files of this stream; the analysis stacks every stream of the tag
  STEM_F  <- sprintf("%s_s%d", CFG$tag, stream)
  CSV     <- file.path(OUT, sprintf("%s.csv", STEM_F))
  CSV_PAR <- file.path(OUT, sprintf("%s-params.csv", STEM_F))
  CSV_GRD <- file.path(OUT, sprintf("%s-grid.csv", STEM_F))
  CSV_STA <- file.path(OUT, sprintf("%s-stations.csv", STEM_F))
  CSV_TO  <- file.path(OUT, sprintf("%s-timeouts.csv", STEM_F))
  DIR_OBS <- file.path(OUT, sprintf("%s-obs", STEM_F))

  cells <- sim_cells()
  prev  <- if (file.exists(CSV)) utils::read.csv(CSV, stringsAsFactors = FALSE) else NULL
  ## every stream of the tag found here, for the times: the cost and the time
  ## limits are measured on all of them
  all_prev <- sim_read_streams(OUT, CFG$tag)
  per_block <- table(factor(cells$block, levels = unique(cells$block)))
  cat(sprintf("%d cells (%s) x %d replications in stream %d, writing to\n  %s\n\n", nrow(cells),
              paste(names(per_block), per_block, sep = " ", collapse = ", "),
              CFG$nrep, stream, file.path(OUT, paste0(STEM_F, "*"))))

  if (identical(CFG$mode[1], "dry")) {
    print(cells[, c("block", "scenario", "n", "TN", "omega", "spread", "balance", "knn")],
          row.names = FALSE)
    cat("\n")
    sim_print_cost(cells, all_prev, CFG)
    return(invisible(cells))
  }

  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  dir.create(DIR_OBS, recursive = TRUE, showWarnings = FALSE)
  done <- if (is.null(prev)) character(0) else paste(prev$cell, prev$rep, sep = "|")

  ## The durations of the replications that ended normally, by cell and by
  ## (n, T), from which the time limits are computed; they start from what is
  ## already recorded.
  T_CELL <- new.env(); T_NT <- new.env()
  remember <- function(cell, n, TN, secs) {
    if (!is.finite(secs)) return(invisible(NULL))
    assign(cell, c(get0(cell, T_CELL, inherits = FALSE), secs), envir = T_CELL)
    key <- paste(n, TN); assign(key, c(get0(key, T_NT, inherits = FALSE), secs), envir = T_NT)
  }
  if (!is.null(all_prev)) {
    ok <- is.na(all_prev$error) | !nzchar(all_prev$error)
    for (i in which(ok)) remember(all_prev$cell[i], all_prev$n[i], all_prev$TN[i], all_prev$secs_total[i])
  }
  cap_of <- function(cell) {
    if (!(CFG$cap_mult[1] > 0)) return(Inf)
    own <- get0(cell$cell, T_CELL, inherits = FALSE)
    if (length(own) >= CFG$cap_after[1])
      return(max(CFG$cap_min[1], CFG$cap_mult[1] * stats::median(own)))
    ref <- get0(paste(cell$n, cell$TN), T_NT, inherits = FALSE)
    if (length(ref)) return(max(CFG$cap_first_min[1], CFG$cap_first_mult[1] * stats::median(ref)))
    CFG$cap_first_min[1]
  }
  prepare <- function(task) {
    if (is.null(task$attempt)) task$attempt <- 1L
    task$cap <- cap_of(task$cell)
    task$path <- SIM_FILE
    task
  }

  save_obs <- function(df, cell, rep) {
    if (is.null(df)) return(invisible(NULL))
    f <- sprintf("obs_%s_rep%04d.rds", gsub("[^A-Za-z0-9_.-]", "", cell$cell), rep)
    sim_retry(function() saveRDS(df, file.path(DIR_OBS, f), compress = "xz"))
  }

  reps  <- sim_rep_id(seq.int(CFG$rep_from[1], rep_to), stream)
  tasks <- sim_tasks(cells, reps, done)
  cat(sprintf("stream %d, replications %d to %d (recorded as rep %d to %d): %d tasks to run, %d already recorded\n\n",
              stream, CFG$rep_from[1], rep_to, min(reps), max(reps), length(tasks),
              nrow(cells) * length(reps) - length(tasks)))
  if (!length(tasks)) { cat("nothing to do\n"); return(invisible(NULL)) }

  cat(sprintf("%13s %-12s %4s %5s %4s | %8s %5s %6s | %s\n",
              "", "scenario", "n", "T", "rep", "time", "K_hat", "ARI", "elapsed, and left"))
  t_start <- proc.time()[["elapsed"]]
  on_result <- function(out, task, k, n) {
    if (!is.list(out) || is.null(out$summary))
      out <- sim_failed(task$cell, task$rep, paste(as.character(out), collapse = " "), task$attempt)
    s  <- out$summary
    el <- (proc.time()[["elapsed"]] - t_start) / 3600
    ## a replication over its time limit: recorded apart and, while attempts
    ## are left, drawn again with the seed of the next attempt (the last
    ## attempt is also written to the results, as a failed replication)
    if (grepl("^timeout", s$error))
      sim_append(data.frame(cell = task$cell$cell, stream = sim_stream_of(task$rep), rep = task$rep, attempt = task$attempt,
                            scenario = task$cell$scenario, n = task$cell$n, TN = task$cell$TN,
                            omega = task$cell$omega, spread = task$cell$spread,
                            cap = round(task$cap), secs = round(s$secs_total),
                            seed = sim_seed(task$rep, task$attempt)), CSV_TO)
    if (grepl("^timeout", s$error) && task$attempt < CFG$cap_attempts[1]) {
      cat(sprintf("[%5d/%5d] %-12s %4d %5d %4d | TIMEOUT after %.0f s (limit %.0f s), attempt %d of %d: drawn again\n",
                  k, n, s$scenario, s$n, s$TN, s$rep, s$secs_total, task$cap,
                  task$attempt, CFG$cap_attempts[1]))
      utils::flush.console()
      task$attempt <- task$attempt + 1L
      return(task)
    }
    sim_append(out$summary, CSV)
    sim_append(out$params,  CSV_PAR)
    sim_append(out$grid,    CSV_GRD)
    sim_append(out$station, CSV_STA)
    save_obs(out$obs, task$cell, task$rep)
    if (!nzchar(s$error)) remember(task$cell$cell, task$cell$n, task$cell$TN, s$secs_total)
    cat(sprintf("[%5d/%5d] %-12s %4d %5d %4d | %6.1f s %5s %6s | %s, ~%s left\n",
                k, n, s$scenario, s$n, s$TN, s$rep, s$secs_total, format(s$K_hat),
                format(round(s$ari_true, 3)), sim_hours(el), sim_hours(el / k * (n - k))))
    if (nzchar(s$error)) cat("              FAILED: ", s$error, "\n", sep = "")
    utils::flush.console()
    NULL
  }

  sim_run_tasks(tasks, sim_task_fun, on_result, CFG$cores[1], prepare = prepare)
  cat(sprintf("\ndone in %s\n\n", sim_hours((proc.time()[["elapsed"]] - t_start) / 3600)))
  sim_print_cost(cells, sim_read_streams(OUT, CFG$tag), CFG)
  invisible(NULL)
}

if (!SIM_DEFINE_ONLY) sim_main()
