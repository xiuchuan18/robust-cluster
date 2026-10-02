setwd(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")))
suppressPackageStartupMessages({ library(dplyr) })
source("00_CONFIG.R"); source("Section4_Lib.R"); source("Section5_Datasets.R")
source("CalFlex_FormalLib.R")

OUT <- file.path(RESULTS_DIR, "calibrated-flexibility")
C <- readRDS(file.path(OUT, "CalFlex_Banknote_Clean.rds"))
D <- readRDS(file.path(OUT, "CalFlex_Banknote_CellWise05.rds"))
S <- D %>% group_by(method) %>% summarise(
  attempted = n(), returned = sum(success), return_rate = mean(success),
  convergence_rate = mean(converged, na.rm = TRUE),
  mean_ARI = mean(ARI, na.rm = TRUE), sd_ARI = sd(ARI, na.rm = TRUE),
  se_ARI = sd_ARI / sqrt(returned), mean_seconds = mean(seconds),
  .groups = "drop") %>%
  left_join(C %>% dplyr::select(method, clean_ARI = ARI, clean_success = success),
            by = "method") %>%
  mutate(Method = unname(calflex_method_label[method]),
         ARI_drop = clean_ARI - mean_ARI) %>%
  rename(method_code = method)
calflex_write_csv(S, file.path(OUT, "CalFlex_Banknote_Summary.csv"))
cat("Banknote analysis complete; visualization is integrated into ",
    "Section5_CellWise05_ARI.pdf\n", sep = "")
