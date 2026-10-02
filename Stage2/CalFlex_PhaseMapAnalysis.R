setwd(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")))
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(tidyr) })
source("00_CONFIG.R"); source("Section4_Lib.R"); source("Section5_Datasets.R")
source("CalFlex_FormalLib.R"); source("PlotStyle.R")
OUT <- file.path(RESULTS_DIR, "calibrated-flexibility")
FIG <- file.path(FIGURES_DIR, "calibrated-flexibility")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
D <- readRDS(file.path(OUT, "CalFlex_PhaseMap.rds"))

paths <- D %>% group_by(mechanism, severity, separation, method, rep, K_true) %>%
  summarise(
    complete = all(ok & is.finite(BIC)),
    Khat = if (complete) K[which.min(BIC)] else NA_integer_,
    selected_ARI = if (complete) ARI[which.min(BIC)] else NA_real_,
    true_K_ARI = ARI[K == first(K_true)][1],
    Khat_package = if (all(ok & is.finite(BIC_package)))
      K[which.min(BIC_package)] else NA_integer_,
    Khat_over = if (all(ok & is.finite(BIC_over)))
      K[which.min(BIC_over)] else NA_integer_,
    seconds_path = sum(seconds), .groups = "drop") %>%
  mutate(Method = unname(calflex_method_label[method]),
         outcome = case_when(!complete ~ "incomplete", Khat > K_true ~ "split",
                             Khat < K_true ~ "merge", TRUE ~ "correct"),
         package_disagrees = Khat != Khat_package,
         over_disagrees = Khat != Khat_over)

summary <- paths %>% group_by(mechanism, severity, separation, Method) %>%
  summarise(attempted = n(), complete_rate = mean(complete),
    correct_rate = mean(outcome == "correct"), split_rate = mean(outcome == "split"),
    merge_rate = mean(outcome == "merge"), mean_Khat = mean(Khat, na.rm = TRUE),
    mean_selected_ARI = mean(selected_ARI, na.rm = TRUE),
    mean_true_K_ARI = mean(true_K_ARI, na.rm = TRUE),
    se_true_K_ARI = sd(true_K_ARI, na.rm = TRUE) /
      sqrt(sum(is.finite(true_K_ARI))),
    mean_seconds_path = mean(seconds_path), .groups = "drop")
calflex_write_csv(paths, file.path(OUT, "CalFlex_PhaseMap_Paths.csv"))
calflex_write_csv(summary, file.path(OUT, "CalFlex_PhaseMap_Summary.csv"))

condition_levels <- c(
  "Cell scale; sep. 2.5; severity 0",
  "Cell scale; sep. 2.5; severity 2/3",
  "Cell scale; sep. 2.5; severity 1",
  "Cell scale; sep. 4; severity 0",
  "Cell scale; sep. 4; severity 2/3",
  "Cell scale; sep. 4; severity 1",
  "Skew mismatch; sep. 2.5; severity 1",
  "Tail mismatch; sep. 2.5; severity 1")
condition_label <- function(mechanism, severity, separation) {
  sev <- ifelse(abs(severity) < 1e-10, "0",
    ifelse(abs(severity - 2/3) < 1e-10, "2/3", "1"))
  sep <- ifelse(abs(separation - 4) < 1e-10, "4", "2.5")
  family <- c(cellSD = "Cell scale", skewMM = "Skew mismatch",
              tailMM = "Tail mismatch")[mechanism]
  paste0(family, "; sep. ", sep, "; severity ", sev)
}
summary$Condition <- factor(condition_label(summary$mechanism, summary$severity,
                                             summary$separation),
                            levels = rev(condition_levels))
summary$Method <- factor(summary$Method,
  levels = c("GMM", "MSGHD", "MixSIT-VV", "MixSIT-EV", "MixSIT-EE"))

H <- bind_rows(
  summary %>% transmute(Condition, Method, Metric = "Recovery of true K",
                        value = correct_rate,
                        label = scales::percent(correct_rate, accuracy = 1)),
  summary %>% transmute(Condition, Method, Metric = "ARI at the true K",
                        value = mean_true_K_ARI,
                        label = sprintf("%.2f", mean_true_K_ARI))) %>%
  mutate(Metric = factor(Metric,
    levels = c("Recovery of true K", "ARI at the true K")),
    text_colour = ifelse(value >= 0.63, "white", "#202020"))

p <- ggplot(H, aes(Method, Condition, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_text(aes(label = label, colour = text_colour), size = 2.35,
            show.legend = FALSE) +
  facet_grid(. ~ Metric) +
  scale_fill_gradientn(colours = c("#F7FBFF", "#9ECAE1", "#3182BD", "#08306B"),
                       limits = c(0, 1), breaks = c(0, .25, .5, .75, 1),
                       labels = scales::percent, name = NULL) +
  scale_colour_identity() +
  labs(x = NULL, y = NULL) +
  theme_journal() +
  theme(axis.text.x = element_text(angle = 38, hjust = 1, size = 7.2),
        axis.text.y = element_text(size = 7.3),
        panel.grid = element_blank(), legend.position = "bottom",
        legend.key.width = grid::unit(28, "pt"))
journal_save(file.path(FIG, "CalFlex_PhaseMap.pdf"), p,
             width = 7.2, height = 4.55)
cat("phase-map analysis complete\n")
