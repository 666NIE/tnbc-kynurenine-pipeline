# ============================================================
# S4 新增面板：within-TNBC KP 家族蛋白丰度排序图（终版）
# 数据：Krug et al. Cell 2020 补充表（NIHMS1687926）
# 口径：严格 IHC 三阴性定义（ER/PR/ERBB2 全 negative，n=16，与论文 S4f 一致）
# ============================================================
library(tidyverse)
library(readxl)

# ---------- 0. 文件位置（自动匹配 Downloads 里的文件名，不用手抄） ----------
F_PROT <- list.files("data-raw/manual_downloads", pattern = "Table.?2.*xlsx$", full.names = TRUE)[1]
F_CLIN <- list.files("data-raw/manual_downloads", pattern = "Table.?1.*xlsx$", full.names = TRUE)[1]
message("蛋白表 = ", F_PROT)
message("临床表 = ", F_CLIN)

# ---------- 1. 读入：自动定位 Proteome sheet ----------
sheets_prot <- excel_sheets(F_PROT)
sheet_prot  <- sheets_prot[str_detect(sheets_prot, regex("proteome", ignore_case = TRUE))]
stopifnot(length(sheet_prot) == 1)
message("Proteome sheet = ", sheet_prot)

prot <- read_excel(F_PROT, sheet = sheet_prot)
clin <- read_excel(F_CLIN, sheet = "A) Metadata")

# ---------- 2. 严格 IHC 定义 TNBC（与论文口径一致，n=16） ----------
tnbc_ids <- clin %>%
  filter(ER.Updated.Clinical.Status == "negative",
         PR.Clinical.Status == "negative",
         ERBB2.Updated.Clinical.Status == "negative") %>%
  pull(Sample.ID)
message("n(TNBC, strict) = ", length(tnbc_ids))
clin %>% filter(Sample.ID %in% tnbc_ids) %>% count(PAM50) %>% print()

# ---------- 3. 样本列自动识别 ----------
sample_cols <- intersect(colnames(prot), clin$Sample.ID)
message("匹配到样本列数 = ", length(sample_cols))

# ---------- 4. KP 家族：缺失基因提示 + 中位丰度排序 ----------
kp_genes <- c("IDO1","IDO2","TDO2","IL4I1","KYNU","KMO","HAAO",
              "AFMID","AADAT","GPT2","QPRT","CCBL2","NADSYN1")
found <- prot %>% filter(geneSymbol %in% kp_genes) %>% distinct(geneSymbol) %>% pull()
message("未定量到的 KP 基因: ", paste(setdiff(kp_genes, found), collapse = ", "))

rank_df <- prot %>%
  filter(geneSymbol %in% kp_genes) %>%
  group_by(geneSymbol) %>%
  slice_max(numSpectraProteinObserved, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(geneSymbol, all_of(sample_cols)) %>%
  mutate(across(all_of(sample_cols), as.numeric)) %>%
  pivot_longer(-geneSymbol, names_to = "sample", values_to = "abund") %>%
  filter(sample %in% tnbc_ids) %>%
  group_by(geneSymbol) %>%
  summarise(median_log2 = median(abund, na.rm = TRUE),
            n_detected  = sum(!is.na(abund)),
            .groups = "drop") %>%
  arrange(desc(median_log2)) %>%
  mutate(geneSymbol = factor(geneSymbol, levels = rev(geneSymbol)),
         highlight = case_when(geneSymbol == "AFMID" ~ "AFMID",
                               geneSymbol == "IDO1"  ~ "IDO1",
                               TRUE                  ~ "other"))

print(rank_df)

# ---------- 5. 作图 ----------
p <- ggplot(rank_df, aes(x = geneSymbol, y = median_log2, fill = highlight)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = sprintf("%.2f", median_log2)), hjust = -0.15, size = 3.4) +
  coord_flip() +
  scale_fill_manual(values = c("AFMID" = "#B06A3B",
                               "IDO1"  = "#3A5A78",
                               "other" = "#BDBDBD")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL,
       y = "Median protein abundance in TNBC (log2 ratio)",
       title = "KP protein ranking within TNBC (CPTAC, Krug et al. 2020, n=16)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 13),
        panel.grid.major.y = element_blank())

ggsave("figures/FigS4_new_panel_TNBC_KP_protein_ranking.png",
       p, width = 6.2, height = 4.2, dpi = 300)
ggsave("figures/FigS4_new_panel_TNBC_KP_protein_ranking.pdf",
       p, width = 6.2, height = 4.2)

# ---------- 6. 图注用关键数值 ----------
cat(sprintf("n(TNBC) = %d | 可定量 KP 蛋白 = %d 个\n",
            length(tnbc_ids), nrow(rank_df)))
cat(sprintf("第1名 %s = %.2f | 最后一名 %s = %.2f | IDO1 = %.2f (第%d位)\n",
            as.character(rank_df$geneSymbol[1]), rank_df$median_log2[1],
            as.character(rank_df$geneSymbol[nrow(rank_df)]), rank_df$median_log2[nrow(rank_df)],
            rank_df$median_log2[rank_df$geneSymbol == "IDO1"],
            which(as.character(rank_df$geneSymbol) == "IDO1")))
