# 04b_fig2_withE.R -- Figure 2 五面板版（A-D 原样 + E 换入聚类检验）
source("R/00_config.R")
suppressPackageStartupMessages({ library(tidyverse); library(patchwork); library(ConsensusClusterPlus) })

coh <- readRDS(file.path(dirs$results, "01_cohort.rds"))
scr <- readRDS(file.path(dirs$results, "02_scores.rds"))
cl  <- readRDS(file.path(dirs$results, "03_cluster.rds"))
scores <- scr$scores; cluster <- cl$cluster
sub_cols <- c(C1 = "#F8766D", C2 = "#619CFF", C3 = "#00BA38")
stopifnot(all(coh$patients == scores$patient), all(coh$patients == names(cluster)))
df <- scores %>% mutate(cluster = cluster)
cor_lab <- function(x, y) {
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", use = "complete.obs"))
  sprintf("Spearman rho = %.3f\n%s", unname(ct$estimate), fmt_p(ct$p.value))
}

p3a <- ggplot(df, aes(cluster, Kynurenine, fill = cluster)) +
  geom_boxplot(alpha = .8, outlier.shape = NA) +
  geom_jitter(width = .15, alpha = .6, size = 1.8) +
  scale_fill_manual(values = sub_cols) + theme_paper() + theme(legend.position = "none") +
  labs(title = "Kynurenine score", x = NULL,
       y = "Kynurenine ssGSEA (z)", caption = fmt_p(kw_p(Kynurenine ~ cluster, df)))
p3b <- ggplot(df, aes(cluster, IDO1, fill = cluster)) +
  geom_boxplot(alpha = .8, outlier.shape = NA) +
  geom_jitter(width = .15, alpha = .6, size = 1.8) +
  scale_fill_manual(values = sub_cols) + theme_paper() + theme(legend.position = "none") +
  labs(title = "IDO1 expression", x = NULL, y = "log2(FPKM+1)",
       caption = fmt_p(kw_p(IDO1 ~ cluster, df)))
p3c <- ggplot(df, aes(Kynurenine, Cytotoxicity, color = cluster)) +
  geom_point(size = 2.5, alpha = .8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = sub_cols) + theme_paper() +
  labs(title = "Kynurenine vs cytotoxicity", x = "Kynurenine score",
       y = "Cytotoxicity score", color = "Subtype",
       caption = cor_lab(df$Kynurenine, df$Cytotoxicity))
p3d <- ggplot(df, aes(Kynurenine, T_cell, color = cluster)) +
  geom_point(size = 2.5, alpha = .8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = sub_cols) + theme_paper() +
  labs(title = "Kynurenine vs T-cell infiltration", x = "Kynurenine score",
       y = "T-cell score", color = "Subtype",
       caption = cor_lab(df$Kynurenine, df$T_cell))

set.seed(params$seed)
kyn_abs <- calc_z(params$gene_sets$Kynurenine, coh$expr)
scores_swap <- scores %>% mutate(Kynurenine = kyn_abs)
pm <- as.matrix(t(scores_swap %>% select(-patient, -IDO1, -TDO2, -Down_score)))
colnames(pm) <- coh$patients
cc2 <- ConsensusClusterPlus(pm, maxK = params$maxK, reps = params$reps, pItem = params$pItem,
  clusterAlg = "km", distance = "euclidean",
  plot = "pdf", title = "ConsensusClustering_absKyn", seed = params$seed)
cl2 <- cc2[[params$k_final]]$consensusClass
names(cl2) <- colnames(pm)
orig <- as.character(cl$cluster[colnames(pm)])
top_cl <- names(which.max(tapply(kyn_abs, cl2, mean)))
pct_c1 <- round(100 * sum(orig == "C1" & cl2 == top_cl) / sum(orig == "C1"), 1)
e_dat <- data.frame(swap_cluster = factor(paste0("swap-", cl2), levels = paste0("swap-", sort(unique(cl2)))),
                    Kyn_abs = kyn_abs, T_cell = scores$T_cell, Cytotoxicity = scores$Cytotoxicity) %>%
  pivot_longer(c(Kyn_abs, T_cell, Cytotoxicity), names_to = "metric", values_to = "z") %>%
  group_by(swap_cluster, metric) %>% summarise(z = mean(z), .groups = "drop") %>%
  mutate(metric = factor(metric, levels = c("Kyn_abs", "T_cell", "Cytotoxicity"),
                         labels = c("Kyn absolute z", "T-cell score", "Cytotoxicity")))
e_cols <- c("Kyn absolute z" = "#8B1A1A", "T-cell score" = "#3B6FB5", "Cytotoxicity" = "#2E8B57")
p3e <- ggplot(e_dat, aes(swap_cluster, z, fill = metric)) +
  geom_col(position = position_dodge(width = .75), width = .68, alpha = .9) +
  geom_hline(yintercept = 0, linewidth = .4) +
  scale_fill_manual(values = e_cols) + theme_paper() +
  labs(title = "Swap-in clustering (absolute z replaces ssGSEA)",
       x = "Swap-in cluster", y = "mean z", fill = NULL,
       caption = sprintf("Original C1 in top absolute-Kyn cluster: %s%%", pct_c1)) +
  theme(legend.position = "bottom")

fig2 <- (p3a | p3b) / (p3c | p3d) / p3e + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))
ggsave(file.path(dirs$figures, "Figure2_trap.pdf"), fig2, width = 12, height = 15.5)
message("Figure2 五面板版已出: figures/Figure2_trap.pdf")
