# 18_reviewer_round9.R -- 第9轮评审稳健性补充包
source("R/00_config.R")
suppressPackageStartupMessages({library(tidyverse); library(ConsensusClusterPlus)})
coh <- readRDS(file.path(dirs$results, "01_cohort.rds"))
scr <- readRDS(file.path(dirs$results, "02_scores.rds"))
cl  <- readRDS(file.path(dirs$results, "03_cluster.rds"))
scores <- scr$scores; cluster <- cl$cluster

cat("\n===== (a) TIDE 效应量 =====\n")
f9 <- read.csv(file.path(dirs$results, "08_tide.csv"))
kw <- kruskal.test(TIDE ~ cluster, data = f9)
H <- unname(kw$statistic); n <- nrow(f9); k <- 3
d1 <- f9$TIDE[f9$cluster == "C1"]; d23 <- f9$TIDE[f9$cluster != "C1"]
sp <- sqrt(((length(d1)-1)*var(d1) + (length(d23)-1)*var(d23)) / (length(d1)+length(d23)-2))
cat(sprintf("KW H = %.2f, epsilon^2 = %.3f\n", H, (H-k+1)/(n-k)))
cat(sprintf("C1 vs C2+C3: mean %.2f vs %.2f, Cohen's d = %.2f\n", mean(d1), mean(d23), (mean(d1)-mean(d23))/sp))

cat("\n===== (b) 全局转录背景偏相关（Kyn vs 细胞毒性）=====\n")
bg_med <- apply(coh$expr, 2, median)
hk <- intersect(c("GAPDH","ACTB","B2M","HPRT1","RPLP0","PPIA"), rownames(coh$expr))
pcor_resid <- function(x, y, covs) {
  rx <- resid(lm(x ~ ., data = as.data.frame(covs)))
  ry <- resid(lm(y ~ ., data = as.data.frame(covs)))
  suppressWarnings(cor.test(rx, ry, method = "spearman"))
}
r1 <- pcor_resid(scores$Kynurenine, scores$Cytotoxicity, data.frame(T = scores$T_cell))
r2 <- pcor_resid(scores$Kynurenine, scores$Cytotoxicity, data.frame(T = scores$T_cell, bg = bg_med))
cat(sprintf("校正 T_cell:          rho = %.3f, p = %s\n", unname(r1$estimate), format.pval(r1$p.value)))
cat(sprintf("校正 T_cell+全局中位: rho = %.3f, p = %s\n", unname(r2$estimate), format.pval(r2$p.value)))
if (length(hk) >= 3) {
  r3 <- pcor_resid(scores$Kynurenine, scores$Cytotoxicity,
                   data.frame(T = scores$T_cell, bg = colMeans(coh$expr[hk, , drop = FALSE])))
  cat(sprintf("校正 T_cell+管家基因(%d个): rho = %.3f, p = %s\n", length(hk), unname(r3$estimate), format.pval(r3$p.value)))
}

cat("\n===== (c) 单细胞细胞数范围与剔除稳健性（fig7 评分文件）=====\n")
mf <- path.expand("~/Downloads/fig7_output/Fig7_metadata_with_scores.csv")
if (!file.exists(mf)) stop("找不到 ", mf)
meta7 <- read.csv(mf, check.names = FALSE, row.names = 1)
imm <- meta7[meta7$Cell_class %in% c("Myeloid", "T_cell"), ]
cnt <- imm %>% count(orig.ident, Cell_class) %>% pivot_wider(names_from = Cell_class, values_from = n)
paired <- cnt %>% filter(!is.na(Myeloid), !is.na(T_cell))
cat(sprintf("配对患者 %d; 髓系 %d-%d (中位 %.0f); T细胞 %d-%d (中位 %.0f)\n",
    nrow(paired), min(paired$Myeloid), max(paired$Myeloid), median(paired$Myeloid),
    min(paired$T_cell), max(paired$T_cell), median(paired$T_cell)))
agg <- imm %>% group_by(orig.ident, Cell_class) %>%
  summarise(Kyn = mean(Kyn_score, na.rm = TRUE), Exh = mean(Exhaustion_score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Cell_class, values_from = c(Kyn, Exh)) %>%
  filter(!is.na(Kyn_Myeloid), !is.na(Exh_T_cell))
tot <- paired$Myeloid + paired$T_cell; names(tot) <- paired$orig.ident
agg$tot <- tot[agg$orig.ident]
rf <- suppressWarnings(cor.test(agg$Kyn_Myeloid, agg$Exh_T_cell, method = "spearman"))
q <- quantile(agg$tot, c(0.1, 0.9)); keep <- agg$tot > q[1] & agg$tot < q[2]
rt <- suppressWarnings(cor.test(agg$Kyn_Myeloid[keep], agg$Exh_T_cell[keep], method = "spearman"))
cat(sprintf("全样本 rho = %.3f (n=%d); 剔除细胞数最高/最低 10%% 后 rho = %.3f (n=%d)\n",
    unname(rf$estimate), nrow(agg), unname(rt$estimate), sum(keep)))

cat("\n===== (d) 全通路统一 ssGSEA 聚类 =====\n")
allss <- bind_cols(tibble(patient = coh$patients),
  as_tibble(sapply(names(params$gene_sets), function(nm) simple_ssgsea(coh$expr, params$gene_sets[[nm]])))) %>%
  mutate(across(-patient, ~ as.numeric(scale(.x))))
pm2 <- as.matrix(t(allss %>% select(-patient)))
colnames(pm2) <- coh$patients
set.seed(params$seed)
cc3 <- ConsensusClusterPlus(pm2, maxK = params$maxK, reps = params$reps, pItem = params$pItem,
  clusterAlg = "km", distance = "euclidean", plot = "pdf", title = "ConsensusClustering_allssGSEA", seed = params$seed)
cl3 <- cc3[[params$k_final]]$consensusClass; names(cl3) <- colnames(pm2)
orig <- as.character(cluster[colnames(pm2)])
print(table(origC1 = orig == "C1", allss_cluster = cl3))
top3 <- names(which.max(tapply(allss$Kynurenine, cl3, mean)))
cat(sprintf("全ssGSEA 高Kyn簇 = %s (n=%d); 原C1 落入比例 %.1f%%\n",
    top3, sum(cl3 == top3), 100*sum(orig=="C1" & cl3==top3)/sum(orig=="C1")))
print(data.frame(cl = cl3, Kyn = allss$Kynurenine, T_cell = scores$T_cell, Cyto = scores$Cytotoxicity) %>%
  group_by(cl) %>% summarise(across(everything(), mean), n = n(), .groups = "drop"))

cat("\n===== (e) C1 vs C3 下游节点差异倍数 =====\n")
c1 <- names(cluster)[cluster == "C1"]; c3 <- names(cluster)[cluster == "C3"]
for (g in c("AFMID","AADAT","GPT2")) {
  tt <- t.test(coh$expr[g, c1], coh$expr[g, c3])
  cat(sprintf("%s: C1 %.2f vs C3 %.2f, log2FC = %+.2f, p = %s\n",
      g, mean(coh$expr[g, c1]), mean(coh$expr[g, c3]),
      mean(coh$expr[g, c1])-mean(coh$expr[g, c3]), format.pval(tt$p.value)))
}
message("18 done")
