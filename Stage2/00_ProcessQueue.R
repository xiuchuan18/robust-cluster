# Shared fail-fast subprocess queue for the experiment drivers.

run_process_queue <- function(jobs, script, args_fun, done_fun, label_fun,
                              n_parallel = 8L, timeout_s = 7200,
                              poll_s = 1) {
  if (!requireNamespace("processx", quietly = TRUE))
    stop("package 'processx' is required")
  n_parallel <- as.integer(n_parallel)
  if (!is.finite(n_parallel) || n_parallel < 1L)
    stop("n_parallel must be a positive integer")
  if (!nrow(jobs)) return(invisible(NULL))

  logs <- file.path(RESULTS_DIR, "_logs")
  dir.create(logs, recursive = TRUE, showWarnings = FALSE)
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows")
    "Rscript.exe" else "Rscript")
  running <- list()
  next_job <- 1L
  completed <- 0L
  started <- Sys.time()

  stop_children <- function() {
    for (x in running) if (x$process$is_alive())
      try(x$process$kill_tree(), silent = TRUE)
  }
  on.exit(stop_children(), add = TRUE)

  while (next_job <= nrow(jobs) || length(running)) {
    while (length(running) < n_parallel && next_job <= nrow(jobs)) {
      job <- jobs[next_job, , drop = FALSE]
      label <- gsub("[^A-Za-z0-9_.-]", "_", label_fun(job))
      log <- file.path(logs, paste0(label, ".log"))
      p <- processx::process$new(rscript, c(script, args_fun(job)),
        stdout = log, stderr = "2>&1", cleanup = TRUE,
        windows_hide_window = TRUE)
      running[[label]] <- list(process = p, job = job, start = Sys.time(),
                               log = log)
      next_job <- next_job + 1L
    }

    Sys.sleep(poll_s)
    finished <- character()
    for (label in names(running)) {
      x <- running[[label]]
      elapsed <- as.numeric(difftime(Sys.time(), x$start, units = "secs"))
      if (x$process$is_alive() && elapsed > timeout_s) {
        try(x$process$kill_tree(), silent = TRUE)
        stop("job timed out after ", timeout_s, " seconds: ", label,
             "; log: ", x$log)
      }
      if (!x$process$is_alive()) {
        status <- x$process$get_exit_status()
        if (!identical(status, 0L) || !isTRUE(done_fun(x$job)))
          stop("job failed or produced no valid cache: ", label,
               " (status ", status, "); log: ", x$log)
        completed <- completed + 1L
        finished <- c(finished, label)
      }
    }
    running[finished] <- NULL
    cat(sprintf("\r  %d queued, %d running, %d done, %.1f min elapsed   ",
      nrow(jobs) - next_job + 1L, length(running), completed,
      as.numeric(difftime(Sys.time(), started, units = "mins"))))
    flush.console()
  }
  cat("\n")
  invisible(NULL)
}
