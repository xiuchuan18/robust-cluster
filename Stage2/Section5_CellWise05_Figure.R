setwd(normalizePath(".", winslash = "/"))
suppressWarnings(suppressMessages({
  source("00_CONFIG.R"); source("Section5_ContaminationLib.R"); library(ggplot2)
}))
source("PlotStyle.R")

datasets <- c("Banknote", "Biopsy", "AIS", "Wine")
methods <- c("GMM", "tMM", "SNMM", "STMM", "MSGHD", "CGMM",
             "MixSIT_None", "MixSIT_U", "MixSIT_Full")

# The first three datasets come from the final Section 5 run. Banknote follows
# exactly the same alpha=.05, eta=5, 50-replicate protocol and is stored by the
# calibrated-flexibility runner.
D_main <- readRDS(section5_contam_result("rds"))
C_main <- readRDS(file.path(RESULTS_DIR, "Section5_RealData.rds"))
D_bank <- readRDS(file.path(RESULTS_DIR, "calibrated-flexibility",
                            "CalFlex_Banknote_CellWise05.rds"))
C_bank <- readRDS(file.path(RESULTS_DIR, "calibrated-flexibility",
                            "CalFlex_Banknote_Clean.rds"))

keep_contam <- c("dataset", "rep", "method", "ARI", "success")
keep_clean <- c("dataset", "method", "ARI", "success")
D <- rbind(D_main[keep_contam], D_bank[keep_contam])
C <- rbind(C_main[keep_clean], C_bank[keep_clean])
D <- D[D$dataset %in% datasets & D$method %in% methods, ]
C <- C[C$dataset %in% datasets & C$method %in% methods, ]

labs <- c(Banknote = "Banknote (m = 6)", Biopsy = "Biopsy (m = 9)",
          AIS = "AIS (m = 11)", Wine = "Wine (m = 13)")
D$Method <- factor(unname(METHOD_CODE_TO_LABEL[D$method]),
                   levels = rev(METHOD_LABEL_ORDER))
C$Method <- factor(unname(METHOD_CODE_TO_LABEL[C$method]),
                   levels = rev(METHOD_LABEL_ORDER))
D$Dataset <- factor(labs[D$dataset], levels = unname(labs[datasets]))
C$Dataset <- factor(labs[C$dataset], levels = unname(labs[datasets]))
D$status <- ifelse(is.finite(D$ARI), "contaminated fit", "fit not returned")
D$plotted <- ifelse(is.finite(D$ARI), D$ARI, -0.04)
C$status <- "original data"

p <- ggplot(D, aes(plotted, Method, colour = Method)) +
  geom_vline(xintercept = 0, colour = "grey72", linewidth = 0.3) +
  geom_boxplot(data = subset(D, is.finite(ARI)), width = 0.62,
               outlier.shape = NA, linewidth = 0.45, na.rm = TRUE) +
  geom_point(aes(shape = status), alpha = 0.34, size = 0.65, stroke = 0.45,
             position = position_jitter(height = 0.11, width = 0, seed = 2026),
             na.rm = TRUE) +
  geom_point(data = C, aes(ARI, Method, shape = status), inherit.aes = FALSE,
             fill = "white", colour = "black", size = 1.9,
             stroke = 0.7, na.rm = TRUE) +
  facet_wrap(~Dataset, nrow = 1) +
  scale_colour_methods(drop = FALSE, guide = "none") +
  scale_shape_manual(values = c("original data" = 23,
                                "contaminated fit" = 16,
                                "fit not returned" = 4),
                     breaks = c("original data", "contaminated fit",
                                "fit not returned")) +
  scale_x_continuous(limits = c(-0.06, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Adjusted Rand index", y = NULL, shape = NULL) +
  theme_journal() +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank())

out <- file.path(FIGURES_DIR, "Section5_CellWise05_ARI.pdf")
journal_save(out, p, width = 7.2, height = 3.65)
cat("figure written to", out, "\n")
