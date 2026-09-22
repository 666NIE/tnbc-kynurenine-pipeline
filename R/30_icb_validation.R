# 30_icb_validation.R — I-SPY2 (GSE194040) ICB 队列验证：CD8/IDO1(+AFMID) panel vs pCR
# 设计：主检验 = pembro 臂 TNBC 内 CD8-low/IDO1-low(trap) vs pCR；交互 = trap x 臂
# 数据流向：series matrix(注释, 自解析) + gene-level 表达文件 -> 打分(00_config 同口径) -> 分析
source("R/00_config.R")
suppressPackageStartupMessages({ library(tidyverse); library(data.table) })

ispy2_dir <- "~/Downloads"   # 三个下载文件所在目录，按需修改

# ---- 1. 注释表（绕过 GEOquery，直接解析 series matrix） ----
parse_geo_meta <- function(f) {
  lines <- readLines(gzfile(f))
  grab <- function(prefix) {
    hit <- lines[startsWith(lines, prefix)]
    if (length(hit) == 0) return(NULL)
    lapply(hit, function(h) {
      v <- strsplit(sub(paste0("^", prefix, "\t"), "", h), "\t")[[1]]
      gsub("^\"|\"$", "", v)
    })
  }
  chars <- grab("!Sample_characteristics_ch1")
  keys  <- vapply(chars, function(x) sub(":.*$", "", x[1]), character(1))
  vals  <- do.call(rbind, lapply(chars, function(x) sub("^[^:]*:\\s*", "", x)))
  meta  <- as.data.frame(t(vals), stringsAsFactors = FALSE)
  colnames(meta) <- make.unique(keys)
  meta <- cbind(geo_accession = grab("!Sample_geo_accession")[[1]],
                title = grab("!Sample_title")[[1]], meta)
  rownames(meta) <- NULL
  meta
}
meta <- rbind(
  parse_geo_meta(file.path(ispy2_dir, "GSE194040-GPL20078_series_matrix.txt.gz")),
  parse_geo_meta(file.path(ispy2_dir, "GSE194040-GPL30493_series_matrix.txt.gz"))
)
meta$pid <- as.character(meta$`patient id`)
cat("样本总数:", nrow(meta), " 重复患者:", sum(duplicated(meta$pid)), "\n")
meta <- meta %>% filter(!duplicated(pid))   # 跨平台重复者保留首条

# ---- 2. 表达矩阵 ----
expr <- fread(file.path(ispy2_dir,
  "GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988.txt.gz"))
gids <- as.character(expr[[1]])
emat <- as.matrix(expr[, -1])
rownames(emat) <- gids
colnames(emat) <- as.character(colnames(emat))

common <- intersect(meta$pid, colnames(emat))
cat("匹配样本:", length(common), "\n")
meta <- meta %>% filter(pid %in% common)
emat <- emat[, meta$pid]

# ---- 3. 打分（与研究队列同口径） ----
ss_methods <- c("Kynurenine", "Glycolysis", "Glutamine", "Adenosine")
sc <- bind_cols(tibble(pid = colnames(emat)),
  as_tibble(sapply(names(params$gene_sets), function(nm) {
    v <- if (nm %in% ss_methods) simple_ssgsea(emat, params$gene_sets[[nm]])
         else calc_z(params$gene_sets[[nm]], emat)
    as.numeric(scale(v))
  })))
cd8_genes <- intersect(c("CD8A","CD8B","CD8G","CD8"), rownames(emat))
sc$CD8_score <- if (length(cd8_genes) > 0) as.numeric(scale(colMeans(emat[cd8_genes, , drop = FALSE]))) else NA
for (g in c("IDO1","AFMID","SLC7A5")) {
  if (g %in% rownames(emat)) sc[[g]] <- as.numeric(scale(emat[g, ]))
}

# ---- 4. 合并临床 + TNBC 定义 + trap 定义 ----
d <- meta %>%
  select(pid, hr, her2, mp, pcr, arm, `adjustment method`) %>%
  left_join(sc, by = "pid") %>%
  mutate(tnbc = hr == 0 & her2 == 0,
         pcr  = as.integer(pcr),
         trap = CD8_score < 0 & IDO1 < 0)   # CD8-low/IDO1-low，队列内 z=0 切分（同稿件口径）

pembro <- d %>% filter(tnbc, arm == "Paclitaxel + Pembrolizumab")
ctrl   <- d %>% filter(tnbc, arm == "Paclitaxel")
cat("TNBC n — pembro:", nrow(pembro), " paclitaxel:", nrow(ctrl), "\n")
cat("pembro 臂 pCR 率:", round(mean(pembro$pcr), 3), " trap 占比:", round(mean(pembro$trap), 3), "\n")

# ---- 分析 1：主检验（pembro 臂内 trap vs pCR） ----
tab <- table(pembro$trap, pembro$pcr)
print(tab)
ft  <- fisher.test(tab)
cat(sprintf("Fisher p = %.4g\n", ft$p.value))
m1 <- glm(pcr ~ trap, data = pembro, family = binomial)
cat(sprintf("OR(trap) = %.2f (%.2f-%.2f)\n",
    exp(coef(m1))[2], exp(confint(m1))[2,1], exp(confint(m1))[2,2]))
m1b <- glm(pcr ~ trap + mp, data = pembro, family = binomial)  # MammaPrint 校正
print(summary(m1b)$coefficients)

# ---- 分析 2：交互（trap x 治疗臂） ----
d2 <- d %>% filter(tnbc, arm %in% c("Paclitaxel + Pembrolizumab", "Paclitaxel")) %>%
  mutate(pembro = arm != "Paclitaxel")
m2 <- glm(pcr ~ trap * pembro, data = d2, family = binomial)
cat("\n=== interaction ===\n"); print(summary(m2)$coefficients)

# ---- 分析 3：结构复现（trap 相关谱是否在 ICB 队列重现） ----
ct <- cor.test(d$T_cell[d$tnbc], d$Kynurenine[d$tnbc], method = "spearman")
cat(sprintf("\nTNBC 内 Kyn~T_cell: rho=%.3f p=%.3g\n", unname(ct$estimate), ct$p.value))
ct2 <- cor.test(d$Cytotoxicity[d$tnbc], d$Kynurenine[d$tnbc], method = "spearman")
cat(sprintf("TNBC 内 Kyn~Cytotoxicity: rho=%.3f p=%.3g\n", unname(ct2$estimate), ct2$p.value))

# ---- 5. 落盘 ----
write.csv(d, file.path(dirs$results, "30_ispY2_annot_scores.csv"), row.names = FALSE)
write.csv(data.frame(test = c("fisher_trap_pcr","OR_trap","OR_lo","OR_hi",
                              "interaction_p","rho_kyn_tcell","rho_kyn_cyto"),
                     value = c(ft$p.value, exp(coef(m1))[2],
                               exp(confint(m1))[2,1], exp(confint(m1))[2,2],
                               summary(m2)$coefficients[4,4],
                               unname(ct$estimate), unname(ct2$estimate))),
          file.path(dirs$results, "30_icb_validation_summary.csv"), row.names = FALSE)
message("30 done")
