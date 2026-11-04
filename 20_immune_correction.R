# 20_immune_correction.R -- 标志基因面板免疫浸润校正（导师评审稳健性验证）
source("R/00_config.R")
suppressPackageStartupMessages(library(tidyverse))
coh <- readRDS(file.path(dirs$results, "01_cohort.rds"))
scr <- readRDS(file.path(dirs$results, "02_scores.rds"))
cl  <- readRDS(file.path(dirs$results, "03_cluster.rds"))
scores <- scr$scores; cluster <- cl$cluster
markers <- list(
  B_cells    = c("CD19","MS4A1","CD79A","CD79B"),
  CD8_T      = c("CD8A","CD8B"),
  CD4_T      = c("CD4","IL7R","CCR7"),
  NK         = c("NKG7","GNLY","KLRD1"),
  Macrophage = c("CD68","CD14","CSF1R","CD163"),
  Neutrophil = c("FCGR3B","CSF3R","CEACAM8"),
  DC         = c("CD1C","XCR1","CLEC9A","LILRB4"),
  Treg       = c("FOXP3","IL2RA","CTLA4")
)
set_score <- function(genes) {
  g <- intersect(genes, rownames(coh$expr))
  stopifnot(length(g) >= 2)
  colMeans(t(scale(t(coh$expr[g, , drop = FALSE]))))
}
S <- sapply(markers, set_score)
imm <- rowSums(S); names(imm) <- coh$patients
df <- scores %>% mutate(cluster = cluster, imm = imm[patient])
for (nm in colnames(S)) df[[nm]] <- S[, nm]
print(df %>% group_by(cluster) %>%
  summarise(across(all_of(colnames(S)), function(x) round(mean(x), 2)), n = n(), .groups = "drop"))
cat(sprintf("免疫评分 KW p = %s\n", format.pval(kruskal.test(imm ~ cluster, data = df)$p.value)))
pcor_resid <- function(x, y, covs) {
  rx <- resid(lm(x ~ ., data = as.data.frame(covs)))
  ry <- resid(lm(y ~ ., data = as.data.frame(covs)))
  suppressWarnings(cor.test(rx, ry, method = "spearman"))
}
r2 <- pcor_resid(df$Kynurenine, df$Cytotoxicity, data.frame(T = df$T_cell, imm = df$imm))
cat(sprintf("校正 T_cell+免疫评分: rho = %.3f, p = %s\n", unname(r2$estimate), format.pval(r2$p.value)))
df$resid_kyn <- resid(lm(Kynurenine ~ imm, data = df))
print(df %>% group_by(cluster) %>% summarise(resid_Kyn = round(mean(resid_kyn), 2), n = n(), .groups = "drop"))
cat(sprintf("校正免疫浸润后 Kyn 评分 KW p = %s\n",
    format.pval(kruskal.test(resid_kyn ~ cluster, data = df)$p.value)))
message("20 done")
