source("R/00_config.R")
for (s in c("01_data", "02_scores", "03_cluster", "04_mechanism",
            "05_metabolic", "06_validation", "07_figures", "08_fig9",
            "09_tcga_survival", "10_tide_supp")) {
  message("
================ ", s, " ================")
  source(file.path("R", paste0(s, ".R")))
}
message("
All done.")
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

