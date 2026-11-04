# 32_afmid_metabric_check.R —— 独立版（无需 source config，任何工作目录可运行）
# KP 通路节点基因在 METABRIC M1 vs M3 中的表达检验
# 验证"AFMID 是 M1 中唯一维持表达的 KP 节点"这一课题假说的跨队列稳健性

# 00_config.R 中的 Kynurenine 基因集（照抄内置）
KYN_GENES <- c("IDO1","IDO2","TDO2","AFMID","KYNU","KMO","KYAT1","AADAT",
               "CCBL2","GPT2","HAAO","QPRT","NADSYN1","IL4I1")

# ---- 自动定位文件（多个候选目录里找，找不到就弹窗手动选）----
find_file <- function(fname){
  for (d in c("results", "02_results",
              "data-raw/manual_downloads/results", "data-raw/manual_downloads/02_results",
              "data-raw/manual_downloads")) {
    if (dir.exists(d)) {
      f <- list.files(d, pattern = paste0("^", fname, "$"),
                      full.names = TRUE, ignore.case = TRUE)
      if (length(f)) { message("found: ", f[1]); return(f[1]) }
    }
  }
  message("未找到 ", fname, "，请在弹窗中选择该文件：")
  file.choose()
}

f_expr <- find_file("06_metabric_expr.rds")
f_asn  <- find_file("13_metabric_assignments.csv")

mexpr <- readRDS(f_expr)
asg   <- read.csv(f_asn, check.names = FALSE)

# ---- 对齐样本（分型文件列名：patient / cluster；标签：M1–M4）----
common <- intersect(colnames(mexpr), asg$patient)
e   <- mexpr[, common, drop = FALSE]
grp <- asg$cluster[match(common, asg$patient)]
message("matched samples: ", length(common))

# ---- 正式检验：M1 vs M3 ----
genes <- intersect(KYN_GENES, rownames(e))
res <- do.call(rbind, lapply(genes, function(g){
  v  <- as.numeric(e[g, ])
  m1 <- grp == "M1"; m3 <- grp == "M3"
  data.frame(gene = g,
             M1_mean = round(mean(v[m1]), 3),
             M3_mean = round(mean(v[m3]), 3),
             logFC_M1_vs_M3 = round(mean(v[m1]) - mean(v[m3]), 3),
             p_wilcox = signif(wilcox.test(v[m1], v[m3])$p.value, 3))
}))
res$p_adj <- signif(p.adjust(res$p_wilcox, "BH"), 3)
res <- res[order(res$p_wilcox), ]
print(res)

write.csv(res, file.path(dirname(f_expr), "32_kpnodes_metabric_M1vsM3.csv"),
          row.names = FALSE)
message("saved: ", file.path(dirname(f_expr), "32_kpnodes_metabric_M1vsM3.csv"))
