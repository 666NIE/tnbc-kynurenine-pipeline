# 45_cptac_stratified.R — CPTAC 按 KP-high/IDO1-low 分层 + 连续模型的 AFMID 蛋白分析
# 数据：Krug et al., Cell 2020 补充表（NIHMS1687926，Downloads 内自动定位）
# 产出：① 二分分层 Wilcoxon ② TNBC 描述 ③ Spearman 相关 ④ PAM50 调整线性模型
source("R/00_config.R")
suppressPackageStartupMessages({ library(tidyverse); library(readxl) })

# ---------- 自动定位：认列名不认文件名 ----------
find_sheet <- function(colname) {
  cands <- list.files("data-raw/manual_downloads", pattern = "xlsx$", full.names = TRUE)
  for (f in cands) {
    for (s in excel_sheets(f)) {
      hdr <- tryCatch(names(read_excel(f, sheet = s, n_max = 1)),
                      error = function(e) NULL)
      if (!is.null(hdr) && colname %in% hdr) {
        cat("定位成功:", basename(f), "| sheet:", s, "\n")
        return(list(file = f, sheet = s))
      }
    }
  }
  stop("Downloads 里没找到含 '", colname, "' 列的表")
}
pp   <- find_sheet("geneSymbol")
prot <- read_excel(pp$file, sheet = pp$sheet)
cc   <- find_sheet("Sample.ID")
clin <- read_excel(cc$file, sheet = cc$sheet)

# ---------- 数据整理 ----------
kp13 <- c("IDO1","IL4I1","KYNU","KMO","HAAO","AFMID","AADAT","GPT2","QPRT","CCBL2","NADSYN1")
samples <- intersect(colnames(prot), clin$Sample.ID)
m <- prot %>%
  filter(geneSymbol %in% kp13) %>%
  select(geneSymbol, all_of(samples)) %>%
  mutate(across(all_of(samples), as.numeric)) %>%
  group_by(geneSymbol) %>% slice(1) %>% ungroup() %>%
  column_to_rownames("geneSymbol") %>% as.matrix()
cat("矩阵维度:", paste(dim(m), collapse = " × "), "（10 = 可定量 KP 蛋白，同 S4 panel k 口径）\n")
mz <- t(scale(t(m)))

df <- tibble(
  sample   = samples,
  PAM50    = clin$PAM50[match(samples, clin$Sample.ID)],
  KP_score = colMeans(mz, na.rm = TRUE),
  IDO1_z   = mz["IDO1", ],
  AFMID_z  = mz["AFMID", ]
) %>%
  mutate(C1_like = KP_score > 0 & IDO1_z < 0,
         TNBC    = sample %in% (
           clin %>% filter(ER.Updated.Clinical.Status == "negative",
                           PR.Clinical.Status == "negative",
                           ERBB2.Updated.Clinical.Status == "negative") %>%
             pull(Sample.ID)))
df$PAM50 <- factor(df$PAM50)

# ---------- ① 二分分层 ----------
cat("\n===== ① 全体样本：C1 样 vs 非 C1 样 =====\n")
print(df %>% group_by(C1_like) %>%
  summarise(n = n(), AFMID_z_中位 = round(median(AFMID_z, na.rm = TRUE), 3), .groups = "drop"))
print(wilcox.test(AFMID_z ~ C1_like, data = df))

# ---------- ② TNBC 子集（描述性） ----------
cat("\n===== ② TNBC 子集（描述性）=====\n")
print(df %>% filter(TNBC) %>%
  group_by(C1_like) %>%
  summarise(n = n(), AFMID_z_中位 = round(median(AFMID_z, na.rm = TRUE), 3), .groups = "drop"))

# ---------- ③ Spearman 相关 ----------
cat("\n===== ③ AFMID 蛋白 vs KP 蛋白分数 =====\n")
print(cor.test(df$AFMID_z, df$KP_score, method = "spearman"))
print(df %>% filter(TNBC) %>%
  summarise(rho = cor(AFMID_z, KP_score, method = "spearman"), n = n()))

# ---------- ④ PAM50 调整线性模型（主分析，n=122） ----------
cat("\n===== ④ 线性模型：AFMID_z ~ KP_score + IDO1_z + PAM50 =====\n")
fit <- lm(AFMID_z ~ KP_score + IDO1_z + PAM50, data = df)
print(summary(fit)$coefficients[c("KP_score", "IDO1_z"), ])
cat("模型 R² =", round(summary(fit)$r.squared, 3), "\n")

write.csv(df, file.path(dirs$results, "45_cptac_stratified.csv"), row.names = FALSE)
message("45 done")
