# 46_afmid_dual_origin.R — AFMID 细胞类型来源图（图3E 素材，终版）
# 数据：results/17_kp_by_major.csv（17 号脚本产物；缺失时自动先跑 17）
# 读图要点：AFMID 以肿瘤上皮为主要表达区室（32.8%），兼髓系（10.5%）；
#           IDO1 集中于髓系（8.9%）——代谢臂=肿瘤内在，炎症臂=免疫细胞
# 注：标题与坐标轴一律英文（Mac PDF 设备不支持中文渲染，标题由拼版大标题承担）
source("R/00_config.R")
suppressPackageStartupMessages({ library(tidyverse) })

f <- file.path(dirs$results, "17_kp_by_major.csv")
if (!file.exists(f)) {
  message("未找到 ", f, "，先运行 17 号脚本……")
  source(file.path("R", "17_scrna_kp_heatmap.R"))
}
long <- read.csv(f, stringsAsFactors = FALSE)

dat <- long %>%
  filter(gene %in% c("AFMID", "IDO1")) %>%
  mutate(gene = factor(gene, levels = c("AFMID", "IDO1")))

# ★ 硬验证（答辩报数用）：AFMID 应以 Cancer Epithelial 居首、Myeloid≈10.5%；IDO1 应 Myeloid≈8.9%
print(dat %>% filter(gene == "AFMID") %>%
        select(celltype_major, pct, avg) %>% arrange(desc(pct)))
print(dat %>% filter(gene == "IDO1") %>%
        select(celltype_major, pct, avg) %>% arrange(desc(pct)))

p <- ggplot(dat, aes(x = pct, y = reorder(celltype_major, pct))) +
  geom_col(width = 0.68, fill = "#B06A3B") +
  facet_wrap(~ gene, ncol = 2, scales = "free_x") +
  geom_text(aes(label = sprintf("%.1f%%", pct)), hjust = -0.12, size = 3.2) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "% of expressing cells", y = NULL, title = NULL) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 12))

ggsave(file.path(dirs$figures, "Fig3E_AFMID_cell_origin.pdf"), p, width = 7.4, height = 3.4)
ggsave(file.path(dirs$figures, "Fig3E_AFMID_cell_origin.png"), p, width = 7.4, height = 3.4, dpi = 300)
write.csv(dat, file.path(dirs$results, "46_afmid_dual_origin.csv"), row.names = FALSE)
message("46 done — figures/Fig3E_AFMID_cell_origin.pdf/png")
