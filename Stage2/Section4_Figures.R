# Section4_Figures.R
#
# Draws Figure 1 (ARI boxplot) and Figure 2 (BIC boxplot) from Section4_Simulation_Results.rds. Uses Visualization_0220.Rmd palette and
# the manuscript-style method/scenario labels.

source("00_CONFIG.R")
source("PlotStyle.R")

res_file <- file.path(RESULTS_DIR, "Section4_Simulation_Results.rds")
if (!file.exists(res_file))
  stop("Run Section4_Run.R first; results not found at ", res_file)

res_list <- readRDS(res_file)

scenario_map <- tibble::tibble(
  ScenarioKey = c("S1_Gaussian", "S2_Imbalanced", "S3_GaussCopula",
                  "S4_Skewed", "S7_CellWise10", "S6_AlignedWide"),
  Scenario    = c("S1: Gaussian",
                  "S2: Imbalanced",
                  "S3: Heavy tailed",
                  "S4: Strongly skewed",
                  "S5: Cell-wise contamination",
                  "S6: Aligned, m = 10")
)

build_long <- function(metric = c("ARI", "BIC")) {
  metric <- match.arg(metric)
  out <- purrr::imap_dfr(res_list, function(df, scen_key) {
    if (!(scen_key %in% scenario_map$ScenarioKey)) return(NULL)
    metric_cols <- grep(paste0("_", metric, "$"), names(df), value = TRUE)
    raw_prefixes <- sub(paste0("_", metric, "$"), "", metric_cols)
    raw_prefixes <- raw_prefixes[raw_prefixes %in% names(METHOD_CODE_TO_LABEL)]
    if (length(raw_prefixes) == 0) return(NULL)

    purrr::map_dfr(raw_prefixes, function(pref) {
      vals    <- df[[paste0(pref, "_", metric)]]
      success <- df[[paste0(pref, "_Success")]]
      ok      <- (success == 1) & !is.na(vals)
      tibble::tibble(
        ScenarioKey = scen_key,
        Model       = unname(METHOD_CODE_TO_LABEL[pref]),
        Value       = as.numeric(vals[ok])
      )
    })
  })

  out %>%
    dplyr::left_join(scenario_map, by = "ScenarioKey") %>%
    dplyr::mutate(
      Scenario = factor(Scenario, levels = scenario_map$Scenario),
      Model    = factor(Model,    levels = METHOD_LABEL_ORDER)
    ) %>%
    dplyr::select(Scenario, Model, Value)
}

plot_box <- function(dat, ylab) {
  ggplot(dat, aes(x = Model, y = Value, fill = Model)) +
    geom_boxplot(outlier.size = 0.45, linewidth = 0.35, width = 0.72) +
    facet_wrap(~ Scenario, ncol = 3, scales = "free_y") +
    scale_fill_methods(drop = FALSE) +
    labs(x = NULL, y = ylab) +
    theme_journal() +
    theme(
      axis.text.x        = element_text(size = 7.2, angle = 40, hjust = 1),
      panel.grid.major.x = element_blank(),
      legend.position    = "none"
    )
}

ari_long <- build_long("ARI")
bic_long <- build_long("BIC")

p_ari <- plot_box(ari_long, "Adjusted Rand Index")
p_bic <- plot_box(bic_long, "BIC")

journal_save(file.path(FIGURES_DIR, "Figure1_ARI_Boxplots.pdf"),
             p_ari, width = 7.2, height = 5.15)
journal_save(file.path(FIGURES_DIR, "Figure2_BIC_Boxplots.pdf"),
             p_bic, width = 7.2, height = 5.15)

cat("Figures written:\n")
cat("  ", file.path(FIGURES_DIR, "Figure1_ARI_Boxplots.pdf"), "\n")
cat("  ", file.path(FIGURES_DIR, "Figure2_BIC_Boxplots.pdf"), "\n")
