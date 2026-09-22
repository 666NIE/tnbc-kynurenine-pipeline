source("R/00_config.R")
for (s in c("01_data", "02_scores", "03_cluster", "04_mechanism",
            "05_metabolic", "06_validation", "07_figures", "08_fig9",
            "09_tcga_survival", "10_tide_supp", "11_extra_analyses",
            "12_k2_gap", "13_metabric_cluster", "14_sensitivity",
            "15_cluster_abs_sensitivity", "16_impress_validation",
            "17_scrna_kp_heatmap", "18_reviewer_round9", "19_final_robustness", "20_immune_correction", "21_spatial_kp", "22_fig7_combined", "23_fig56_merge", "24_cibersortx_check", "25_random_geneset_control", "26_expression_matched_null", "27_multivariable_arbitration", "28_metabric_replication", "29_metabric_m1_survival",
            "30_icb_validation", "31_icb_validation_figure")) {
  message("\n================ ", s, " ================")
  source(file.path("R", paste0(s, ".R")))
}
message("\nAll done.")
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
