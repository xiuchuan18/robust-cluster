# Section4_RunReplicate.R
#
# One replicate, one operating-system process.
#
#   Rscript Section4_RunReplicate.R <sim_id> <cell>
#
# Writes results/_rep_<cell>_<sim_id>.rds containing list(row = <one-row data frame>, wbar = <latent scales or NULL>) and exits.


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: Rscript Section4_RunReplicate.R <sim_id> <cell>")
sid  <- as.integer(args[1])
cell <- args[2]

setwd(normalizePath(".", winslash = "/"))
out_file <- file.path("results", sprintf("_rep_%s_%d.rds", cell, sid))

# package startup chatter must not reach stdout: the caller polls on files, but a noisy child still costs I/O on a machine running several of these.
suppressWarnings(suppressMessages({
  source("00_CONFIG.R"); source_stage1(); source("Section4_Lib.R")
}))
if (section4_rep_cache_valid(out_file, cell)) quit(save = "no", status = 0)

sc <- CONFIG_SIM$scenarios[[cell]]
if (is.null(sc)) stop("unknown cell: ", cell)

row <- tryCatch(
  run_one_replicate(sid, cell, sc$N, sc$d, CONFIG_SIM$K,
                    CONFIG_SIM$seed, CONFIG_SIM$max_iter),
  error = function(e) NULL)

if (is.null(row)) quit(save = "no", status = 2)

wb <- attr(row, "wbar")
attr(row, "wbar") <- NULL
save_result(list(row = row, wbar = wb), out_file, compress = "xz")
quit(save = "no", status = 0)
