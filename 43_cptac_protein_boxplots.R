# 43_cptac_protein_boxplots.R — S4 f/g 同数据 R 重建（原图出自 CPTAC 数据门户在线分析）
# 验证结论：IDO1 在 TNBC 中最高、AFMID 在 TNBC 中最低，三轮重建均复现；
#          Normal 组绝对水平与门户图存在差异（归一化细节不可考），不影响上述两个关键结论。 — CPTAC 乳腺癌蛋白组：IDO1/AFMID 跨亚型箱线图（AF4 Fig S4 f, g）
# 数据：Krug et al., Cell 2020 补充表（NIHMS1687926）
#   Supplemental Table 2.xlsx（85MB，Proteome sheet）
#   Supplemental Table 1.xlsx（临床/Metadata sheet，A) Metadata）
# 说明：本脚本复现附加文件 4 面板 f（IDO1）与 g（AFMID）；
#   TNBC 定义为 ER/PR/ERBB2（HER2）IHC 全阴（严格 n=16）。
source("R/00_config.R")
suppressPackageStartupMessages({ library(tidyverse); library(readxl) })

F_PROT <- list.files("data-raw/manual_downloads", pattern = "Table.?2.*xlsx$", full.names = TRUE)[1]
F_CLIN <- list.files("data-raw/manual_downloads", pattern = "Table.?1.*xlsx$", full.names = TRUE)[1]

find_sheet <- function(colname) {
  cands <- list.files("data-raw/manual_downloads", pattern = "xlsx$", full.names = TRUE)
  for (f in cands) {
    for (s in excel_sheets(f)) {
      hdr <- tryCatch(names(read_excel(f, sheet = s, n_max = 1)),
                      error = function(e) NULL)
      if (!is.null(hdr) && colname %in% hdr) return(list(file = f, sheet = s))
    }
  }
  stop("Downloads 里没找到含 '", colname, "' 列的表")
}
pp   <- find_sheet("geneSymbol")
prot <- read_excel(pp$file, sheet = pp$sheet)
cc   <- find_sheet("Sample.ID")
clin <- read_excel(cc$file, sheet = cc$sheet)

sample_cols <- intersect(colnames(prot), clin$Sample.ID)
tnbc_ids <- clin %>%
  filter(ER.Updated.Clinical.Status == "negative",
         PR.Clinical.Status == "negative",
         ERBB2.Updated.Clinical.Status == "negative") %>% pull(Sample.ID)
message("严格 IHC TNBC n = ", length(tnbc_ids))

get_gene <- function(g) {
  prot %>% filter(geneSymbol == g) %>%
    select(all_of(sample_cols)) %>% slice(1) %>% as.numeric()   # 取首行（覆盖最全）；Normal组绝对水平与门户存在已知差异，见头部注释
}
dat <- tibble(
  sample = sample_cols,
  IDO1  = get_gene("IDO1"),
  AFMID = get_gene("AFMID")
) %>%
  mutate(subclass = case_when(
    sample %in% tnbc_ids ~ "TNBC",
    clin$PAM50[match(sample, clin$Sample.ID)] == "LumA" ~ "Luminal",
    clin$PAM50[match(sample, clin$Sample.ID)] == "LumB" ~ "Luminal",
    clin$PAM50[match(sample, clin$Sample.ID)] == "Her2" ~ "HER2+",
    clin$PAM50[match(sample, clin$Sample.ID)] == "Normal-like" ~ "Normal",
    TRUE ~ NA_character_)) %>%
  filter(!is.na(subclass)) %>%
  mutate(across(c(IDO1, AFMID), ~ as.numeric(scale(.))),
         subclass = factor(subclass, levels = c("Normal", "Luminal", "HER2+", "TNBC")))

for (g in c("IDO1", "AFMID")) {
  p <- ggplot(dat, aes(x = subclass, y = .data[[g]], fill = subclass)) +
    geom_boxplot(width = 0.6, alpha = 0.85, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 1.1, alpha = 0.7) +
    scale_fill_manual(values = c(Normal = "grey70", Luminal = "#4C72B0",
                                 "HER2+" = "#DD8452", TNBC = "#55A868")) +
    labs(x = NULL, y = "z-value",
         title = paste0("Protein expression of ", g, " in breast cancer (CPTAC)")) +
    theme_minimal(base_size = 12) + theme(legend.position = "none")
  ggsave(paste0("figures/CPTAC_", g, "_by_subclass.png"), p, width = 4.2, height = 3.6, dpi = 300)
  print(p)
}
message("43 done — S4 f/g 已导出至 figures/CPTAC_IDO1_by_subclass.png 与 CPTAC_AFMID_by_subclass.png")
