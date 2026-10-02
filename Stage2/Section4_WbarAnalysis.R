# Section4_WbarAnalysis.R
#
# What the latent scales are doing on the cell-wise cells.
#
# Output: results/Section4_Wbar_Summary.csv and the manuscript density figure.

source("00_CONFIG.R")
suppressWarnings(suppressMessages(source("Section4_Lib.R")))
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })
source("PlotStyle.R")

CELLS <- c("S7_CellWise10")
VARIANTS <- c("MixSIT_None", "MixSIT_U", "MixSIT_Full")

auc_neg <- function(w, y) {                 # AUC of -w as a detector of y == 1
  y <- as.logical(y); n1 <- sum(y); n0 <- sum(!y)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(-as.numeric(w))
  (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

rows <- list(); dens <- list()
for (cell in CELLS) {
  f <- file.path(RESULTS_DIR, paste0("Section4_Wbar_", cell, ".rds"))
  if (!file.exists(f)) { cat("missing:", f, "\n"); next }
  store <- readRDS(f)
  cat(sprintf("%-16s %d replicates\n", cell, length(store)))
  if (!all(vapply(store, function(rec)
      identical(rec$axis_match_version, WBAR_AXIS_MATCH_VERSION), logical(1))))
    stop("legacy or unmatched Wbar records found for ", cell,
         "; rerun the targeted Wbar recomputation")

  for (v in VARIANTS) {
    st <- list()
    for (rec in store) {
      WB <- rec$wbar[[v]]; C <- rec$contam
      if (is.null(WB) || is.null(C) ||
          !identical(WB$matching, WBAR_AXIS_MATCH_VERSION)) next
      W <- WB$weights
      Cl <- C == 1
      nbad <- rowSums(Cl)

      # within-row contrast, on rows with exactly one corrupted cell
      one <- which(nbad == 1)
      wr <- if (length(one)) {
        j   <- max.col(Cl[one, , drop = FALSE], "first")
        bad <- W[cbind(one, j)]
        sib <- vapply(seq_along(one), function(t)
                 median(W[one[t], !Cl[one[t], ]]), numeric(1))
        c(bad = mean(bad), sib = mean(sib), ratio = mean(bad / sib))
      } else c(bad = NA, sib = NA, ratio = NA)

      st[[length(st) + 1]] <- data.frame(
        cell = cell, variant = v, sim_id = rec$sim_id, d = ncol(W),
        w_contam = mean(W[Cl]), w_clean = mean(W[!Cl]),
        auc = auc_neg(W, Cl),
        rows_hit = mean(nbad > 0),          # a row-wise method must give up here
        mean_abs_cosine = mean(WB$axis_abs_cosine),
        min_abs_cosine = min(WB$axis_abs_cosine),
        w_bad_1 = wr[["bad"]], w_sib_1 = wr[["sib"]], w_ratio_1 = wr[["ratio"]])
    }
    if (!length(st)) next
    S <- do.call(rbind, st)
    rows[[length(rows) + 1]] <- data.frame(
      cell = cell, variant = v, d = S$d[1], n_rep = nrow(S),
      w_contam = mean(S$w_contam), w_clean = mean(S$w_clean),
      auc = mean(S$auc), auc_sd = sd(S$auc),
      rows_hit = mean(S$rows_hit),
      mean_abs_cosine = mean(S$mean_abs_cosine),
      min_abs_cosine = min(S$min_abs_cosine),
      w_bad_1 = mean(S$w_bad_1, na.rm = TRUE),
      w_sib_1 = mean(S$w_sib_1, na.rm = TRUE),
      w_ratio_1 = mean(S$w_ratio_1, na.rm = TRUE))
  }

  # Pool all 100 replicates for the density figure. Each replicate has the same
  # number of cells, so this weights experimental replicates equally and avoids
  # presenting a visually convenient but potentially idiosyncratic run.
  for (rec in store) for (v in VARIANTS) {
    WB <- rec$wbar[[v]]; C <- rec$contam
    if (is.null(WB) || is.null(C) ||
        !identical(WB$matching, WBAR_AXIS_MATCH_VERSION)) next
    W <- WB$weights
    dens[[length(dens) + 1]] <- data.frame(
      cell = "S5: Cell-wise contamination (m = 10)",
      variant = METHOD_CODE_TO_LABEL[[v]],
      status = ifelse(as.vector(C) == 1, "Contaminated cells", "Clean cells"),
      wbar = as.vector(W))
  }
}

if (!length(rows)) stop("no Wbar files found; run Section4_Run.R first")

SUM <- do.call(rbind, rows)
write.csv(SUM, file.path(RESULTS_DIR, "Section4_Wbar_Summary.csv"), row.names = FALSE)

cat("\n=============== latent scales on the cell-wise cells ===============\n")
cat(sprintf("%-16s %-11s %3s %8s %8s %8s %9s\n",
            "cell", "variant", "d", "W|bad", "W|clean", "AUC", "rows hit"))
for (i in seq_len(nrow(SUM)))
  cat(sprintf("%-16s %-11s %3d %8.3f %8.3f %8.3f %8.1f%%\n",
      SUM$cell[i], METHOD_CODE_TO_LABEL[[SUM$variant[i]]], SUM$d[i],
      SUM$w_contam[i], SUM$w_clean[i], SUM$auc[i], 100 * SUM$rows_hit[i]))

cat("\n===== within-row contrast (rows with exactly one corrupted cell) =====\n")
cat("A row-wise robust method has ONE weight for the whole row and cannot\n")
cat("produce a gap here. MixSIT can.\n\n")
cat(sprintf("%-16s %-11s %10s %10s %8s\n",
            "cell", "variant", "W corrupted", "W siblings", "ratio"))
for (i in seq_len(nrow(SUM)))
  cat(sprintf("%-16s %-11s %10.3f %10.3f %8.3f\n",
      SUM$cell[i], METHOD_CODE_TO_LABEL[[SUM$variant[i]]],
      SUM$w_bad_1[i], SUM$w_sib_1[i], SUM$w_ratio_1[i]))

# ---- figures ---------------------------------------------------------------
D <- do.call(rbind, dens)
D$cell <- factor(D$cell, levels = "S5: Cell-wise contamination (m = 10)")

p1 <- ggplot(D, aes(x = wbar, fill = status)) +
  geom_density(alpha = 0.58, colour = NA, adjust = 1.05) +
  facet_wrap(~ variant, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c("Clean cells" = "#4E79A7",
                               "Contaminated cells" = "#E15759")) +
  coord_cartesian(xlim = c(0, 1.65)) +
  labs(x = expression(bar(w)[ij]), y = "Density", fill = NULL) +
  theme_journal() +
  theme(legend.position = "top", panel.grid.major.x = element_blank())

journal_save(file.path(FIGURES_DIR, "Section4_Wbar_Density.pdf"), p1,
             width = 7.2, height = 2.55)

cat("\nWritten:\n")
cat("  ", file.path(RESULTS_DIR, "Section4_Wbar_Summary.csv"), "\n")
cat("  ", file.path(FIGURES_DIR, "Section4_Wbar_Density.pdf"), "\n")
