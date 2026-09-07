# 07_figures.R — Figures 5 / 7 (Fig 9 已迁移至 08_fig9.R)
source("R/00_config.R")
suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })
coh     <- readRDS(file.path(dirs$results, "01_cohort.rds"))
scr     <- readRDS(file.path(dirs$results, "02_scores.rds"))
cl      <- readRDS(file.path(dirs$results, "03_cluster.rds"))
scores  <- scr$scores; cluster <- cl$cluster
sub_cols <- c(C1 = "#F8766D", C2 = "#619CFF", C3 = "#00BA38")
df <- scores %>% mutate(cluster = cluster)
cor_lab <- function(x, y) {
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", use = "complete.obs"))
  sprintf("Spearman rho = %.3f, %s", unname(ct$estimate), fmt_p(ct$p.value))
}

# ============ Figure 5 ============
df$Immunosuppressive <- rowMeans(as.matrix(df[, c("Treg","Macrophage","Checkpoint","Exhaustion")]))
p5a <- ggplot(df %>% pivot_longer(c(Treg, Macrophage, Checkpoint, Exhaustion),
                                  names_to = "marker", values_to = "score"),
              aes(cluster, score, fill = cluster)) +
  geom_boxplot(alpha = .8, outlier.shape = NA) +
  geom_jitter(width = .12, alpha = .4, size = .9) +
  facet_wrap(~ marker, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = sub_cols) + theme_paper() + theme(legend.position = "none") +
  labs(title = "Immunosuppressive markers across subtypes (cluster-defining)", x = NULL, y = "score (z)")
p5b <- ggplot(df, aes(Kynurenine, Immunosuppressive, color = cluster)) +
  geom_point(size = 2.5, alpha = .8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = sub_cols) + theme_paper() +
  labs(title = "Kynurenine vs combined immunosuppressive score",
       x = "Kynurenine score", y = "Immunosuppressive score", color = "Subtype",
       caption = cor_lab(df$Kynurenine, df$Immunosuppressive))
p5c <- ggplot(df, aes(cluster, Immunosuppressive, fill = cluster)) +
  geom_boxplot(alpha = .8, outlier.shape = NA) +
  geom_jitter(width = .15, alpha = .6, size = 1.8) +
  scale_fill_manual(values = sub_cols) + theme_paper() + theme(legend.position = "none") +
  labs(title = "Combined immunosuppressive score", x = NULL, y = "score (z)",
       caption = fmt_p(kw_p(Immunosuppressive ~ cluster, df)))
fig5 <- (p5a / (p5b | p5c)) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))
ggsave(file.path(dirs$figures, "Figure5_immunosuppression.pdf"), fig5, width = 11, height = 9)
message("Figure 5 done")

# ============ Figure 7 (metadata 版; A/B/H 面板需计数矩阵, README TODO) ============
meta_file <- path.expand(paths$scrna_meta)
if (file.exists(meta_file)) {
  m <- read.csv(meta_file, check.names = FALSE, row.names = 1)
  need_cols <- c("orig.ident", "Kyn_score", "Exhaustion_score", "Cytotoxic_score", "Cell_class")
  lack <- setdiff(need_cols, colnames(m))
  if (length(lack) > 0) {
    message("Fig7 缺列 ", paste(lack, collapse = ", "), ", 跳过")
  } else {
    has_ido <- "IDO1" %in% colnames(m)
    message("Fig7: IDO1表达列=", has_ido, " (FALSE 则 A/B/H 面板需计数矩阵)")
    agg <- m %>% filter(Cell_class %in% c("Myeloid", "T_cell")) %>%
      group_by(orig.ident, Cell_class) %>%
      summarise(Kyn = mean(Kyn_score, na.rm = TRUE),
                Exh = mean(Exhaustion_score, na.rm = TRUE),
                Cyto = mean(Cytotoxic_score, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = Cell_class, values_from = c(Kyn, Exh, Cyto))
    p7a <- ggplot(agg, aes(Kyn_Myeloid, Exh_T_cell)) +
      geom_point(size = 3, alpha = .8, color = "#C55A11") +
      geom_smooth(method = "lm", se = TRUE, color = "black") + theme_paper() +
      labs(title = "Myeloid kynurenine score vs T-cell exhaustion",
           x = "Myeloid Kynurenine score", y = "T-cell exhaustion score",
           caption = if (sum(complete.cases(agg$Kyn_Myeloid, agg$Exh_T_cell)) > 5)
             cor_lab(agg$Kyn_Myeloid, agg$Exh_T_cell) else "insufficient paired samples")
    p7b <- ggplot(agg, aes(Kyn_Myeloid, Cyto_T_cell)) +
      geom_point(size = 3, alpha = .8, color = "#70AD47") +
      geom_smooth(method = "lm", se = TRUE, color = "black") + theme_paper() +
      labs(title = "Myeloid kynurenine score vs T-cell cytotoxicity",
           x = "Myeloid Kynurenine score", y = "T-cell cytotoxicity score",
           caption = cor_lab(agg$Kyn_Myeloid, agg$Cyto_T_cell))
    fig7 <- (p7a | p7b) + plot_annotation(tag_levels = "A") &
      theme(plot.tag = element_text(size = 14, face = "bold"))
    ggsave(file.path(dirs$figures, "Figure7_singlecell_metadata.pdf"), fig7, width = 11, height = 5)
    message("Figure 7 done")
  }
} else message("Fig7 metadata 文件不存在, 跳过")
message("07 all done")
