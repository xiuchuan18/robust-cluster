# Formal median-time curves, generated solely from the aggregate RDS.
source("00_CONFIG.R")
library(ggplot2)
MTH <- c("GMM", "MSGHD", "MixSIT_None", "MixSIT_U", "MixSIT_Full")
LAB <- c(GMM="GMM", MSGHD="MSGHD", MixSIT_None="MixSIT-VV",
         MixSIT_U="MixSIT-EV", MixSIT_Full="MixSIT-EE")
SC <- readRDS("results/Section6_Scale.rds")
SC <- SC[SC$ok %in% c(TRUE, NA) & is.finite(SC$t_tot), ]
anchor <- list(N = 800, m = 10, K = 3)
sweeps <- list(
  N = SC[SC$m == anchor$m & SC$K == anchor$K, ],
  m = SC[SC$N == anchor$N & SC$K == anchor$K, ],
  K = SC[SC$N == anchor$N & SC$m == anchor$m, ])

# figure
source("PlotStyle.R")
PL <- do.call(rbind, lapply(names(sweeps), function(s) {
  g <- sweeps[[s]]
  a <- aggregate(list(secs = g$t_tot), list(x = g[[s]], method = g$method),
                 median, na.rm = TRUE)
  a$sweep <- switch(s, N = "sample size N", m = "dimension m", K = "components K")
  a }))
PL <- PL[PL$method %in% MTH, ]
PL$Method <- factor(LAB[PL$method], levels = unname(LAB[MTH]))
PL$sweep <- factor(PL$sweep,
                   levels = c("sample size N", "dimension m", "components K"))
p <- ggplot(PL, aes(x = x, y = secs, colour = Method, shape = Method,
                    linetype = Method, group = Method)) +
  geom_line(linewidth = 0.65) + geom_point(size = 1.75, stroke = 0.65) +
  facet_wrap(~ sweep, scales = "free_x") +
  # the default log breaks land between the design points (300, 500, 1000 for a sweep run at 200, 400, 800, 1600), so the axis is set to the configurations
  # actually fitted; each panel shows only those inside its own range
  scale_x_log10(breaks = sort(unique(PL$x))) +
  scale_y_log10(breaks = c(0.2, 0.5, 1, 2, 5, 10, 20, 50, 100)) +
  scale_colour_methods() + scale_shape_methods() + scale_linetype_methods() +
  labs(x = NULL, y = "median seconds per fit (log scale)", colour = NULL) +
  guides(shape = "none", linetype = "none") +
  theme_journal() +
  theme(legend.position = "top")
journal_save(file.path(FIGURES_DIR, "Section6_Scaling.pdf"), p,
             width = 7.2, height = 2.75)
cat(sprintf("\nwritten: %s\n\n", file.path(FIGURES_DIR, "Section6_Scaling.pdf")))
write.csv(PL, file.path(RESULTS_DIR, "Section6_Scaling_Points.csv"), row.names = FALSE)

