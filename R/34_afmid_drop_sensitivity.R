# 34_afmid_drop_sensitivity.R —— 终极修正版（去一基因检验，ConsensusClusterPlus 原参数重跑）
# 回应评审：① 去除 AFMID 后分型是否仍恢复 ② 完整 BH 框架 ③ M2 分子画像
source("R/00_config.R")
suppressPackageStartupMessages(library(ConsensusClusterPlus))

mexpr <- readRDS("results/06_metabric_expr.rds")
asg   <- read.csv("results/13_metabric_assignments.csv", check.names = FALSE)
obj   <- readRDS("results/13_metabric_cluster.rds")     # list(scores, cluster, cc)

common <- intersect(colnames(mexpr), asg$patient)
e      <- mexpr[, common]
grp    <- asg$cluster[match(common, asg$patient)]

kyn_full <- params$gene_sets$Kynurenine
kyn_noAF <- setdiff(kyn_full, "AFMID")
present  <- intersect(kyn_full, rownames(e))

# 基准矩阵：与原聚类输入逐元素一致（仅取 16 通路列，样本对齐）
mss <- obj$scores
sc0 <- as.matrix(mss[match(common, mss$patient), names(params$gene_sets)])
rownames(sc0) <- common
stopifnot(all(colnames(sc0) == names(params$gene_sets)))

# ============ 工具函数 ============
ari_fun <- function(a, b){
  tab <- table(a, b); n <- sum(tab)
  sij <- sum(tab * (tab - 1) / 2)
  ca  <- sum(rowSums(tab) * (rowSums(tab) - 1) / 2)
  cb  <- sum(colSums(tab) * (colSums(tab) - 1) / 2)
  exp <- ca * cb / (n * (n - 1) / 2)
  (sij - exp) / ((ca + cb) / 2 - exp)
}
auc_fun <- function(score, pos){
  r  <- rank(score); n1 <- sum(pos); n0 <- sum(!pos)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# 关键修复：simple_ssgsea 返回无名字向量，必须先命名再按 common 对齐
kyn_col <- function(kyn_set){
  v <- simple_ssgsea(e, intersect(kyn_set, rownames(e)))
  names(v) <- colnames(e)
  as.numeric(scale(v[common]))
}

swap <- function(kyn_set, tag){
  sc <- sc0
  sc[, "Kynurenine"] <- kyn_col(kyn_set)
  stopifnot(!any(is.na(sc)))                    # 防御：有 NA 直接报错，不进聚类
  set.seed(42)
  mcc <- ConsensusClusterPlus(t(sc), maxK = 3, reps = 100, pItem = 0.8,
                              clusterAlg = "km", distance = "euclidean",
                              plot = "pdf", title = paste0("34_swap_", tag), seed = 42)
  cl <- mcc[[3]]$consensusClass
  ord <- order(tapply(sc[, "Kynurenine"], cl, mean), decreasing = TRUE)
  setNames(c("M1", "M2", "M3")[match(cl, ord)], names(cl))
}

lab_full <- swap(kyn_full, "full")              # 阳性对照：重算全基因列，应还原原分型
lab_noAF <- swap(kyn_noAF, "noAFMID")

ari_full  <- ari_fun(grp, lab_full[common])
ari_noAF  <- ari_fun(grp, lab_noAF[common])
m1_recall <- mean(lab_noAF[common][grp == "M1"] == "M1")
auc_full  <- auc_fun(sc0[, "Kynurenine"], grp == "M1")
auc_noAF  <- auc_fun(kyn_col(kyn_noAF), grp == "M1")

cat(sprintf("阳性对照 ARI(重算全基因列重聚类 vs 原分型) = %.3f   <- 应为 1.000\n", ari_full))
cat(sprintf("ARI(去AFMID重聚类 vs 原分型) = %.3f\n", ari_noAF))
cat(sprintf("M1 召回率 = %.1f%%\n", 100 * m1_recall))
cat(sprintf("M1 vs 其余 AUC：全基因 %.3f → 去AFMID %.3f\n", auc_full, auc_noAF))

write.csv(data.frame(patient = common, orig = grp,
                     lab_full = lab_full[common], lab_noAF = lab_noAF[common]),
          "results/34_dropAFMID_recluster.csv", row.names = FALSE)

# ============ 检验 1：全亚型 × 两两比较 × 全基因 BH ============
pw <- function(g){
  v <- as.numeric(e[g, ])
  c(M1vM2 = wilcox.test(v[grp == "M1"], v[grp == "M2"])$p.value,
    M1vM3 = wilcox.test(v[grp == "M1"], v[grp == "M3"])$p.value,
    M2vM3 = wilcox.test(v[grp == "M2"], v[grp == "M3"])$p.value,
    KW    = kruskal.test(v ~ grp)$p.value)
}
res <- as.data.frame(t(sapply(present, pw)))
res$gene <- rownames(res)
res <- res[, c("gene", "M1vM2", "M1vM3", "M2vM3", "KW")]
res$q_M1vM2 <- p.adjust(res$M1vM2, "BH")
qmat <- matrix(p.adjust(as.matrix(res[, 2:4]), "BH"), nrow = nrow(res))
res$q_M1vM2_g <- qmat[, 1]; res$q_M1vM3_g <- qmat[, 2]; res$q_M2vM3_g <- qmat[, 3]
res$q_global  <- apply(qmat, 1, min)
res <- res[order(res$KW), ]
print(res[, c("gene", "M1vM3", "q_M1vM2", "q_M1vM3_g", "q_global")])
write.csv(res, "results/34_kpnodes_metabric_fullBH.csv", row.names = FALSE)

# ============ 检验 3：M2 分子画像 ============
panel <- c(present, "CD8A", "CD8B", "CD3D", "GZMA", "GZMB", "PRF1", "IFNG")
prof <- do.call(rbind, lapply(intersect(panel, rownames(e)), function(g){
  v <- as.numeric(e[g, ])
  data.frame(gene = g,
             M1 = round(mean(v[grp == "M1"]), 2),
             M2 = round(mean(v[grp == "M2"]), 2),
             M3 = round(mean(v[grp == "M3"]), 2),
             M2_vs_M1 = round(mean(v[grp == "M2"]) - mean(v[grp == "M1"]), 2),
             M2_vs_M3 = round(mean(v[grp == "M3"]) - mean(v[grp == "M2"]), 2))
}))
print(prof)
write.csv(prof, "results/34_M2_portrait.csv", row.names = FALSE)

message("done: results/34_dropAFMID_recluster.csv | 34_kpnodes_metabric_fullBH.csv | 34_M2_portrait.csv")
