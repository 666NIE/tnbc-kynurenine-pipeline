# 36_metabolic_only_clustering.R — 终版 v3：if-else 加大括号修复
source("R/00_config.R")
suppressPackageStartupMessages({ library(tidyverse); library(ConsensusClusterPlus) })

# ---------- 1. TCGA：加载 + 标签 ----------
x <- readRDS("results/02_scores.rds")
if (is.matrix(x) || is.data.frame(x)) {
  pw_tcg <- as.data.frame(x)
} else {
  el <- names(x)[sapply(x, function(e) !is.null(dim(e)) && ncol(e) >= 10)][1]
  pw_tcg <- as.data.frame(x[[el]])
}
cat("TCGA 通路矩阵:", nrow(pw_tcg), "×", ncol(pw_tcg), "\n")

lab_df <- read.csv("results/03_cluster_assignments.csv", stringsAsFactors = FALSE)
lab_col <- names(lab_df)[sapply(lab_df, function(v) any(str_detect(as.character(v), "^C[123]$")))][1]
lab_tcg <- as.character(lab_df[[lab_col]])
stopifnot(length(lab_tcg) == nrow(pw_tcg))

# ---------- 2. 显式选取 8 条非免疫通路 ----------
keep <- c("Kynurenine", "Glycolysis", "Glutamine", "Adenosine",
          "Tryptophan_metabolism", "Hypoxia", "EMT", "Stemness")
message("保留 ", length(keep), " 条: ", paste(keep, collapse = ", "))

pw8 <- pw_tcg[, keep] %>% mutate(across(everything(), as.numeric))
cc  <- complete.cases(pw8)
cat("TCGA complete.cases:", sum(cc), "/", nrow(pw8), "\n")
pw8     <- pw8[cc, ]
lab_tcg <- lab_tcg[cc]

# ---------- 3. TCGA：非免疫重聚类（参数与原文一致，k=3） ----------
set.seed(42)
cc8 <- ConsensusClusterPlus(t(as.matrix(pw8)), maxK = 6, reps = 100, pItem = 0.8,
                            clusterAlg = "km", distance = "euclidean", seed = 42)
km3 <- cc8[[3]]$consensusClass
tab_tcg <- table(recluster = km3, original = lab_tcg)
print(tab_tcg)
newc1 <- rownames(tab_tcg)[which.max(tab_tcg[, "C1"])]
cat(sprintf("新簇 %s 中 C1 占比: %.1f%%\n", newc1,
            100 * tab_tcg[newc1, "C1"] / sum(tab_tcg[newc1, ])))
kyncol <- names(pw_tcg)[str_detect(names(pw_tcg), regex("kyn", ignore_case = TRUE))][1]
cat("各新簇平均", kyncol, ":\n")
print(tapply(as.numeric(pw_tcg[[kyncol]][cc]), km3, mean))
saveRDS(list(tab = tab_tcg, labels = km3), "results/36_metabolic_only_tcg.rds")

# ---------- 4. METABRIC：同流程 ----------
xm <- readRDS("results/13_metabric_cluster.rds")
if (is.matrix(xm) || is.data.frame(xm)) {
  pw_met <- as.data.frame(xm)
} else {
  elm <- names(xm)[sapply(xm, function(e) !is.null(dim(e)) && ncol(e) >= 10)][1]
  pw_met <- as.data.frame(xm[[elm]])
}
cat("METABRIC 通路矩阵:", nrow(pw_met), "×", ncol(pw_met), "\n")
print(names(pw_met))

keepm <- intersect(keep, names(pw_met))
labm <- read.csv("results/13_metabric_assignments.csv", stringsAsFactors = FALSE)
labm_col <- names(labm)[sapply(labm, function(v) any(str_detect(as.character(v), "^M[123]$")))][1]
lab_met <- as.character(labm[[labm_col]])

pw8m <- pw_met[, keepm] %>% mutate(across(everything(), as.numeric))
ccm  <- complete.cases(pw8m)
cat("METABRIC complete.cases:", sum(ccm), "/", nrow(pw8m), "\n")
pw8m    <- pw8m[ccm, ]
lab_met <- lab_met[ccm]

set.seed(42)
cc8m <- ConsensusClusterPlus(t(as.matrix(pw8m)), maxK = 6, reps = 100, pItem = 0.8,
                             clusterAlg = "km", distance = "euclidean", seed = 42)
km3m <- cc8m[[3]]$consensusClass
tab_met <- table(recluster = km3m, original = lab_met)
print(tab_met)
saveRDS(list(tab = tab_met, labels = km3m), "results/36_metabolic_only_met.rds")
message("36 done")
