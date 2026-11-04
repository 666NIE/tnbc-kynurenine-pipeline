# 44_pancancer_survey.R — 泛癌调查重建版（AF6 Fig S6 的数据源）
# 说明：原始交互式分析未存档；本脚本按论文 2.11 节方法重建，输出文件名与
#       39_supplementary_composites.R 的读取约定兼容（42b/42c/42d）。
source("R/00_config.R")
suppressPackageStartupMessages({ library(tidyverse); library(survival) })

kp_genes <- c("IDO1","IDO2","TDO2","IL4I1","KYNU","KMO","HAAO",
              "AFMID","AADAT","GPT2","QPRT","CCBL2","NADSYN1","KYAT1")
cyto_genes <- c("CD8A","CD8B","GZMB","PRF1")

# ---------- ① TCGA 泛癌不一致调查 → results/42c ----------
pan <- readRDS("results/44_input_pancancer_expr.rds")     # 首次从 UCSC Xena 下载后保存
pan_clin <- read.csv("results/44_input_pancancer_clin.csv", check.names = FALSE)

kp_present <- intersect(kp_genes, rownames(pan))
KP_score <- colMeans(t(scale(t(pan[kp_present, ]))), na.rm = TRUE)
IDO1_z <- as.numeric(scale(pan["IDO1", ]))
cancer_type <- pan_clin$cancer_type[match(colnames(pan), pan_clin$sample)]

tbl42c <- tibble(cancer_type, discordant = KP_score > 0 & IDO1_z < 0) %>%
  group_by(cancer_type) %>%
  summarise(n_total = n(), n_discordant = sum(discordant, na.rm = TRUE),
            discordant_fraction = round(100 * mean(discordant, na.rm = TRUE), 2), .groups = "drop") %>%
  filter(n_total >= 100) %>% arrange(desc(discordant_fraction))
write.csv(tbl42c, "results/42c_pancancer_discordance.csv", row.names = FALSE)

# ---------- ② DepMap 谱系不一致 → results/42b ----------
dep <- read.csv("results/44_input_depmap_lineage.csv", check.names = FALSE) %>%
  mutate(discordant = KP_z > 0 & IDO1_z < 0)
tbl42b <- dep %>% group_by(lineage) %>%
  summarise(n_cell_lines = n(), n_discordant = sum(discordant),
            discordant_fraction = round(100 * mean(discordant), 2), .groups = "drop")
write.csv(tbl42b, "results/42b_depmap_lineage_discordance.csv", row.names = FALSE)

# ---------- ③ 逐癌种生存（Cox：discordant + age + sex + cytolytic）→ results/42d ----------
cyto <- colMeans(pan[intersect(cyto_genes, rownames(pan)), ], na.rm = TRUE)
surv_df <- pan_clin %>%
  transmute(sample, cancer_type, time = as.numeric(OS_MONTHS), event = as.numeric(OS_EVENT),
            age = as.numeric(age), sex = factor(sex),
            cytolytic = as.numeric(cyto[sample]),
            discordant = KP_score[sample] > 0 & IDO1_z[sample] < 0)
res <- surv_df %>% filter(!is.na(time), !is.na(event)) %>%
  group_by(cancer_type) %>% filter(n() >= 100) %>%
  group_modify(~{
    fit <- tryCatch(coxph(Surv(time, event) ~ discordant + age + sex + cytolytic, data = .x),
                    error = function(e) NULL)
    if (is.null(fit)) return(tibble(HR = NA, CI_lo = NA, CI_hi = NA, p = NA))
    s <- summary(fit)
    tibble(HR = s$conf.int[1,1], CI_lo = s$conf.int[1,3], CI_hi = s$conf.int[1,4], p = s$coefficients[1,5])
  }) %>% ungroup()
res$p_adj_BH <- p.adjust(res$p, method = "BH")
write.csv(res, "results/42d_pancancer_survival.csv", row.names = FALSE)
message("44 done")
