# Section4_Run.R
#
#   Rscript Section4_Run.R [n_parallel]
#
# Resumable: a replicate with a cached _rep_ file is skipped, so an interrupted
# run continues where it stopped.

setwd(normalizePath(".", winslash = "/"))
suppressWarnings(suppressMessages({
  source("00_CONFIG.R"); source("Section4_Lib.R") }))

args    <- commandArgs(trailingOnly = TRUE)
N_PAR   <- if (length(args) >= 1) as.integer(args[1]) else 8
CELLS   <- names(CONFIG_SIM$scenarios)
NSIM    <- CONFIG_SIM$n_sim

rep_path  <- function(cell, sid) file.path("results", sprintf("_rep_%s_%d.rds", cell, sid))
wbar_path <- function(cell) file.path("results", paste0("Section4_Wbar_", cell, ".rds"))

grid <- expand.grid(sid = seq_len(NSIM), cell = CELLS, stringsAsFactors = FALSE)
grid <- grid[order(grid$sid), ]                     # breadth first, so a partial
                                                    # run still covers every cell
done <- mapply(function(cell, sid) section4_rep_cache_valid(rep_path(cell, sid), cell),
               grid$cell, grid$sid)
todo <- grid[!done, , drop = FALSE]
cat(sprintf("Section 4: %d of %d replicates outstanding, %d at a time\n",
            nrow(todo), nrow(grid), N_PAR))

run_process_queue(todo, "Section4_RunReplicate.R",
  args_fun = function(j) c(j$sid, j$cell),
  done_fun = function(j) section4_rep_cache_valid(rep_path(j$cell, j$sid), j$cell),
  label_fun = function(j) sprintf("section4_%s_%03d", j$cell, j$sid),
  n_parallel = N_PAR, timeout_s = 7200)

# assembly
all_results <- list()
for (cell in CELLS) {
  fs <- mapply(rep_path, cell, seq_len(NSIM))
  valid <- mapply(function(f) section4_rep_cache_valid(f, cell), fs)
  if (!all(valid))
    stop("invalid or missing Section 4 cache(s) for ", cell)
  zs <- lapply(fs, readRDS)
  rows <- do.call(rbind, lapply(zs, `[[`, "row"))
  all_results[[cell]] <- rows[order(rows$sim_id), ]

  wb <- Filter(Negate(is.null), lapply(zs, `[[`, "wbar"))
  if (length(wb)) {
    names(wb) <- vapply(wb, function(x) as.character(x$sim_id), character(1))
    save_result(wb, wbar_path(cell), compress = "xz")
  }
  cat(sprintf("  [%s] %d replicates%s\n", cell, nrow(all_results[[cell]]),
              if (length(wb)) sprintf(", %d latent-scale records", length(wb)) else ""))
}
save_result(all_results, file.path(RESULTS_DIR, "Section4_Simulation_Results.rds"))
cat(sprintf("\nassembled %d cells into Section4_Simulation_Results.rds\n", length(all_results)))
