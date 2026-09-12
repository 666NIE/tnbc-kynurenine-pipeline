# 05_metabolic.R — Figure 4 (PGC-1a / OXPHOS)
source("R/00_config.R")
suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })
coh <- readRDS(file.path(dirs$results, "01_cohort.rds"))
scr <- readRDS(file.path(dirs$results, "02_scores.rds"))
cl  <- readRDS(file.path(dirs$results, "03_cluster.rds"))
scores <- scr$scores; cluster <- cl$cluster
sub_cols <- c(C1 = "#F8766D", C2 = "#619CFF", C3 = "#00BA38")
df <- scores %>% mutate(cluster = cluster,
                        PPARGC1A = as.numeric(coh$expr["PPARGC1A", ]),
                        OXPHOS = calc_z(params$etc_genes, coh$expr))
cor_lab <- function(x, y) {
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", use = "complete.obs"))
  sprintf("Spearman rho = %.3f
%s", unname(ct$estimate), fmt_p(ct$p.value))
}
p4A <- ggplot(df, aes(cluster, PPARGC1A, fill = cluster)) +
  geom_boxplot(alpha = .8, outlier.shape = NA) +
  geom_jitter(width = .15, alpha = .6, size = 1.8) +
  scale_fill_manual(values = sub_cols) + theme_paper() + theme(legend.position = "none") +
  labs(title = "PPARGC1A expression across subtypes", x = NULL, y = "log2(FPKM+1)",
       caption = fmt_p(kw_p(PPARGC1A ~ cluster, df)))
p4B <- ggplot(df, aes(cluster, OXPHOS, fill = cluster)) +
  geom_boxplot(alpha = .8, outlier.shape = NA) +
  geom_jitter(width = .15, alpha = .6, size = 1.8) +
  scale_fill_manual(values = sub_cols) + theme_paper() + theme(legend.position = "none") +
  labs(title = "OXPHOS score (19-gene ETC panel) across subtypes", x = NULL,
       y = "OXPHOS score", caption = fmt_p(kw_p(OXPHOS ~ cluster, df)))
p4C <- ggplot(df, aes(PPARGC1A, OXPHOS, color = cluster)) +
  geom_point(size = 2.5, alpha = .8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = sub_cols) + theme_paper() +
  labs(title = "PPARGC1A expression does not predict OXPHOS activity",
       x = "PPARGC1A expression", y = "OXPHOS score",
       caption = cor_lab(df$PPARGC1A, df$OXPHOS))
p4D <- ggplot(df, aes(OXPHOS, Glycolysis, color = cluster)) +
  geom_point(size = 2.5, alpha = .8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = sub_cols) + theme_paper() +
  labs(title = "No compensatory glycolytic switch",
       x = "OXPHOS score", y = "Glycolysis score (cluster-defining)",
       caption = cor_lab(df$OXPHOS, df$Glycolysis))
fig4 <- (p4A | p4B) / (p4C | p4D) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))
ggsave(file.path(dirs$figures, "Figure3_PGC1a_OXPHOS.pdf"), fig4, width = 11, height = 9)
message("05 done")
