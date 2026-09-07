# 04_mechanism.R — Figures 3 / 6 / 8
source("R/00_config.R")
suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })
coh <- readRDS(file.path(dirs$results, "01_cohort.rds"))
scr <- readRDS(file.path(dirs$results, "02_scores.rds"))
cl  <- readRDS(file.path(dirs$results, "03_cluster.rds"))
scores <- scr$scores; cluster <- cl$cluster
sub_cols <- c(C1 = "#F8766D", C2 = "#619CFF", C3 = "#00BA38")
stopifnot(all(coh$patients == scores$patient), all(coh$patients == names(cluster)))
df <- scores %>% mutate(cluster = cluster)
cor_lab <- function(x, y) {
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", use = "complete.obs"))
  sprintf("Spearman rho = %.3f
%s", unname(ct$estimate), fmt_p(ct$p.value))
}

p3a <- ggplot(df, aes(cluster, Kynurenine, fill = cluster)) +
  geom_boxplot(alpha = .8, outlier.shape = NA) +
  geom_jitter(width = .15, alpha = .6, size = 1.8) +
  scale_fill_manual(values = sub_cols) + theme_paper() + theme(legend.position = "none") +
  labs(title = "Kynurenine score across subtypes (cluster-defining)", x = NULL,
       y = "Kynurenine ssGSEA (z)", caption = fmt_p(kw_p(Kynurenine ~ cluster, df)))
p3b <- ggplot(df, aes(cluster, IDO1, fill = cluster)) +
  geom_boxplot(alpha = .8, outlier.shape = NA) +
  geom_jitter(width = .15, alpha = .6, size = 1.8) +
  scale_fill_manual(values = sub_cols) + theme_paper() + theme(legend.position = "none") +
  labs(title = "IDO1 expression across subtypes", x = NULL, y = "log2(FPKM+1)",
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
fig3 <- (p3a | p3b) / (p3c | p3d) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))
ggsave(file.path(dirs$figures, "Figure3_trap.pdf"), fig3, width = 11, height = 9)

feats <- setdiff(names(params$gene_sets), "Kynurenine")
coupling <- bind_rows(lapply(feats, function(f) {
  ct <- suppressWarnings(cor.test(df$Kynurenine, df[[f]], method = "spearman"))
  data.frame(feature = f, rho = unname(ct$estimate), p = ct$p.value)
})) %>%
  mutate(sig = p < 0.05, label = sprintf("rho=%.2f
%s", rho, fmt_p(p)))
p6 <- ggplot(coupling, aes(reorder(feature, rho), rho, fill = sig)) +
  geom_col(width = .7, color = "black", linewidth = .3) +
  geom_text(aes(label = label), size = 2.6,
            hjust = ifelse(coupling$rho > 0, -0.05, 1.05)) +
  scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#B0B0B0")) +
  coord_flip() + theme_paper() + theme(legend.position = "none") +
  ylim(min(coupling$rho) - .25, max(coupling$rho) + .25) +
  labs(title = "Kynurenine-immune coupling across all samples", x = NULL, y = "Spearman rho")
ggsave(file.path(dirs$figures, "Figure6_coupling.pdf"), p6, width = 7.5, height = 6)
write.csv(coupling, file.path(dirs$results, "04_coupling.csv"), row.names = FALSE)

mat <- df %>% select(all_of(names(params$gene_sets)))
cor_m <- cor(mat, method = "spearman", use = "pairwise.complete.obs")
cor_df <- as.data.frame(as.table(cor_m)) %>% filter(as.integer(Var1) < as.integer(Var2))
p8 <- ggplot(cor_df, aes(Var1, Var2, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       limits = c(-1, 1), name = "rho") +
  geom_text(aes(label = sprintf("%.2f", Freq)), size = 2.3) +
  theme_paper() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlation matrix of metabolic-immune features", x = NULL, y = NULL)
ggsave(file.path(dirs$figures, "Figure8_corrmatrix.pdf"), p8, width = 9, height = 7.5)
saveRDS(cor_m, file.path(dirs$results, "04_cormatrix.rds"))
message("04 done")
