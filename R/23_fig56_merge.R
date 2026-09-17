# 23_fig56_merge.R -- 图5合并版：A 偶联谱 + B 16通路 Spearman 相关矩阵
source("R/00_config.R")
suppressPackageStartupMessages({library(tidyverse); library(gridExtra)})

scr <- readRDS(file.path(dirs$results, "02_scores.rds"))
cl  <- readRDS(file.path(dirs$results, "03_cluster.rds"))
scores <- scr$scores %>% mutate(cluster = as.character(cl$cluster))
pw <- setdiff(colnames(scores), c("patient", "cluster", "IDO1", "TDO2", "Down_score"))

cor_df <- tibble(pathway = pw) %>%
  mutate(rho = sapply(pathway, function(p) suppressWarnings(cor.test(scores$Kynurenine, scores[[p]], method = "spearman")$estimate)),
         p   = sapply(pathway, function(p) suppressWarnings(cor.test(scores$Kynurenine, scores[[p]], method = "spearman")$p.value)),
         sig = p < 0.05)
cor_df$pathway <- factor(gsub("_", " ", cor_df$pathway), levels = gsub("_", " ", cor_df$pathway[order(cor_df$rho)]))
pA <- ggplot(cor_df, aes(pathway, rho, fill = sig)) +
  geom_col(width = .7, color = "black", linewidth = .3) +
  scale_fill_manual(values = c("TRUE" = "#C0392B", "FALSE" = "grey75"), guide = "none") +
  coord_flip() + theme_bw(base_size = 11) +
  labs(x = NULL, y = "Spearman rho with kynurenine score", title = "A | Coupling spectrum")

M <- cor(scores %>% select(all_of(pw)) %>% rename_with(~ gsub("_", " ", .)), method = "spearman")
Mdf <- as.data.frame(as.table(M))
Mdf$label <- ifelse(Mdf$Var1 != Mdf$Var2, sprintf("%.2f", Mdf$Freq), "")
pB <- ggplot(Mdf, aes(Var1, Var2, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 2.0) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1)) +
  theme_bw(base_size = 8.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                    axis.text.y = element_text(size = 7)) +
  labs(x = NULL, y = NULL, title = "B | Pathway correlation matrix", fill = "rho")

ggsave(file.path(dirs$figures, "Figure5_merged.pdf"),
       arrangeGrob(pA, pB, ncol = 2, widths = c(1, 1.15)), width = 13, height = 6.5)
message("23 done")
