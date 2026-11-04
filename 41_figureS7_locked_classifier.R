# 41_figureS7_locked_classifier.R — 独立 Figure S7（METABRIC 锁定阈值分类器 KM + 风险表）
source("R/00_config.R")
suppressPackageStartupMessages({ library(tidyverse); library(survival); library(survminer) })

# ---------- 1. 数据（与 3.6 节同口径） ----------
surv <- read.csv("results/06_metabric_surv.csv", check.names = FALSE)
df <- surv %>%
  mutate(kyn_z = as.numeric(scale(Kyn_score)),
         cls   = factor(ifelse(kyn_z > 0.0063, "KP-high", "KP-low"),
                        levels = c("KP-high", "KP-low")),
         time  = os_time, event = os_event,
         age   = as.numeric(AGE_AT_DIAGNOSIS)) %>%
  filter(!is.na(time), !is.na(event), !is.na(kyn_z), !is.na(age))
table(df$cls)

# ---------- 2. 作图 ----------
p <- ggsurvplot(
  survfit(Surv(time, event) ~ cls, data = df),
  data       = df,
  palette    = c("#B06A3B", "#3A5A78"),     # 顺序对齐 levels：KP-high 橙 / KP-low 蓝
  legend.title = "Classifier",
  legend.labs  = c("KP-high", "KP-low"),
  legend     = c(0.85, 0.85),
  pval       = TRUE,                        # 自动计算 log-rank p（显示为 0.12）
  pval.size  = 4.5,
  risk.table = TRUE,                        # 保留风险表（独立附图的信息量优势）
  risk.table.y.text.col = TRUE,
  xlab = "Months", ylab = "Overall survival probability",
  title = "METABRIC: locked-threshold classifier",   # 短标题防裁断；阈值放图注
  ggtheme = theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))
)

# ---------- 3. 导出（显式设备，ggsurvplot 对象走 print） ----------
pdf("figures/FigureS7_METABRIC_locked_classifier.pdf", width = 6, height = 5.6)
print(p); dev.off()
png("figures/FigureS7_METABRIC_locked_classifier.png", width = 6, height = 5.6,
    units = "in", res = 300)
print(p); dev.off()
message("Figure S7 已导出（图注里写：TCGA-derived cutoff z = 0.0063; log-rank p = 0.12）")
