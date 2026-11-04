#!/usr/bin/env Rscript
# 43_purity_rank.R (v2.7 终版) — scores行号为纯数字，统一改用 coh$patients 对齐
source("R/00_config.R")
source("R/01_data.R"); source("R/02_scores.R"); source("R/03_cluster.R")
suppressPackageStartupMessages({library(tidyverse); library(data.table)})

dir.create("results", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)

## ---------- 1. 纯度代理 ----------
immune_cols <- intersect(c("T_cell","NK_cell","Treg","Macrophage"), colnames(scores))
stopifnot(length(immune_cols) == 4)
imm_z <- scores %>% as.data.frame() %>%
  dplyr::select(all_of(c(immune_cols, "EMT"))) %>%
  dplyr::mutate(dplyr::across(everything(), ~ as.numeric(scale(.x))))
purity_df <- imm_z %>%
  dplyr::mutate(immuno_stromal = rowMeans(dplyr::select(., all_of(immune_cols), EMT)),
                tumorPurity     = as.numeric(scale(-immuno_stromal))) %>%
  dplyr::mutate(sample = coh$patients) %>%          # 行号 → 真实条码
  dplyr::select(sample, tumorPurity)
cat("纯度代理（负向免疫-基质评分），范围:",
    round(range(purity_df$tumorPurity), 2), "\n")

## ---------- 2. 纯度校正后的 Kyn 差异 ----------
df <- scores %>% as.data.frame() %>%
  dplyr::mutate(
    sample  = coh$patients,
    cluster = as.character(cluster)[match(coh$patients, names(cluster))],
    grp     = ifelse(cluster == "C1", "C1", "non-C1")) %>%
  dplyr::inner_join(purity_df, by = "sample")
stopifnot(all(!is.na(df$grp)), nrow(df) == length(coh$patients))

cat("\n===== 1) 纯度代理与 Kyn 的相关 =====\n")
print(suppressWarnings(cor.test(df$Kynurenine, df$tumorPurity, method = "spearman")))
print(df %>% dplyr::group_by(grp) %>%
  dplyr::summarise(purity = round(mean(tumorPurity), 3), n = dplyr::n(), .groups = "drop"))

cat("\n===== 2) 校正后 C1 vs non-C1 的 Kyn 差异 =====\n")
df$Kyn_resid <- resid(lm(Kynurenine ~ tumorPurity, data = df))
cat("--- 残差法 ---\n");  print(wilcox.test(Kyn_resid ~ grp, data = df))
cat("--- 多变量模型 ---\n"); print(summary(lm(Kynurenine ~ grp + tumorPurity, data = df))$coefficients)
cat("--- 对照：未校正 ---\n"); print(wilcox.test(Kynurenine ~ grp, data = df))
write.csv(df %>% dplyr::select(sample, grp, Kynurenine, tumorPurity, Kyn_resid),
          "results/38_purity_correction.csv", row.names = FALSE)

## ---------- 3. KP 基因 rank 贡献分解 ----------
expr_full <- fread("~/afmid_project/HiSeqV2")
g0 <- colnames(expr_full)[1]
mat <- expr_full %>% column_to_rownames(g0) %>% as.matrix()
colnames(mat) <- sub("-01[A-Z]?$", "", colnames(mat))
mat <- mat[, intersect(coh$patients, colnames(mat)), drop = FALSE]
cat("\n表达矩阵:", nrow(mat), "基因 x", ncol(mat), "样本\n")

rk  <- apply(mat, 2, function(x) rank(x, ties.method = "average"))
pct <- 1 - rk / nrow(rk)

kp <- intersect(params$gene_sets$Kynurenine, rownames(pct))
contrib <- data.frame(gene = rownames(pct)) %>%
  dplyr::mutate(
    C1_pct    = rowMeans(pct[, df$sample[df$grp == "C1"],    drop = FALSE], na.rm = TRUE),
    nonC1_pct = rowMeans(pct[, df$sample[df$grp != "C1"],    drop = FALSE], na.rm = TRUE)) %>%
  dplyr::filter(gene %in% kp) %>%
  dplyr::mutate(delta = C1_pct - nonC1_pct,
                drives_C1_score = C1_pct > 0.5) %>%
  dplyr::arrange(desc(C1_pct))
print(contrib)
write.csv(contrib, "results/38_kp_rank_contribution.csv", row.names = FALSE)

cat("\n===== 3) 贡献解读 =====\n")
cat("C1 中排名靠前（贡献高分）:", paste(head(contrib$gene, 5), collapse = ", "), "\n")
cat("C1 中排名靠后（拉低或不贡献）:", paste(tail(contrib$gene, 5), collapse = ", "), "\n")
cat("done -> results/43_*.csv\n")
