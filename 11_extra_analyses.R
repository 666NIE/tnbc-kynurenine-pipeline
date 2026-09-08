R/11_extra_analyses.R
# 11_extra_analyses.R — SLC7A5 分析 + C1 vs C3 差异程序分析
source("R/00_config.R")
suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })

coh     <- readRDS(file.path(dirs$results, "01_cohort.rds"))
scr     <- readRDS(file.path(dirs$results, "02_scores.rds"))
cl      <- readRDS(file.path(dirs$results, "03_cluster.rds"))
expr    <- coh$expr
scores  <- scr$scores
cluster <- cl$cluster
gs      <- params$gene_sets
sub_cols <- c(C1 = "#F8766D", C2 = "#619CFF", C3 = "#00BA38")
df <- scores %>% mutate(cluster = cluster)
cor_lab <- function(x, y) {
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", use = "complete.obs"))
  sprintf("Spearman rho = %.3f, %s", unname(ct$estimate), fmt_p(ct$p.value))
}
message("=== Part A: SLC7A5 ===")

# ---- A1: SLC7A5 表达 by 亚型 ----
if ("SLC7A5" %in% rownames(expr)) {
  df$SLC7A5 <- as.numeric(expr["SLC7A5", ])
  pA1 <- ggplot(df, aes(cluster, SLC7A5, fill = cluster)) +
    geom_boxplot(alpha = .8, outlier.shape = NA) +
    geom_jitter(width = .15, alpha = .6, size = 1.8) +
    scale_fill_manual(values = sub_cols) + theme_paper() + theme(legend.position = "none") +
    labs(title = "SLC7A5 expression across subtypes", x = NULL, y = "log2(FPKM+1)",
         caption = fmt_p(kw_p(SLC7A5 ~ cluster, df)))
  ct <- suppressWarnings(cor.test(df$Kynurenine, df$SLC7A5, method = "spearman"))
  message(sprintf("SLC7A5 by subtype: %s | Kyn vs SLC7A5: rho = %.3f, %s",
                  fmt_p(kw_p(SLC7A5 ~ cluster, df)),
                  unname(ct$estimate), fmt_p(ct$p.value)))
} else {
  pA1 <- ggplot() + theme_void() + labs(title = "SLC7A5 not found in matrix")
  message("警告: 表达矩阵中无 SLC7A5 基因")
}

# ---- A2: Kyn 评分 vs SLC7A5 ----
if ("SLC7A5" %in% colnames(df)) {
  pA2 <- ggplot(df, aes(Kynurenine, SLC7A5, color = cluster)) +
    geom_point(size = 2.5, alpha = .8) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    scale_color_manual(values = sub_cols) + theme_paper() +
    labs(title = "Kynurenine score vs SLC7A5 expression",
         x = "Kynurenine score", y = "SLC7A5 expression", color = "Subtype",
         caption = cor_lab(df$Kynurenine, df$SLC7A5))
  figA <- (pA1 | pA2) + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 14, face = "bold"))
  ggsave(file.path(dirs$figures, "Figure11_SLC7A5.pdf"), figA, width = 11, height = 5)
  message("Figure11_SLC7A5.pdf saved")
}

message("=== Part B: C1 vs C3 差异程序 ===")
c1 <- names(cluster)[cluster == "C1"]
c3 <- names(cluster)[cluster == "C3"]
pid_expr11 <- patient_id(colnames(expr))
e1 <- expr[, pid_expr11 %in% c1, drop = FALSE]
e3 <- expr[, pid_expr11 %in% c3, drop = FALSE]
n1 <- ncol(e1); n3 <- ncol(e3)
m1 <- rowMeans(e1); m3 <- rowMeans(e3)
v1 <- apply(e1, 1, var); v3 <- apply(e3, 1, var)
se <- sqrt(v1 / n1 + v3 / n3)
tstat <- (m1 - m3) / se
deg <- data.frame(gene = rownames(expr), logFC = m1 - m3, t = tstat) %>%
  arrange(desc(t))
write.csv(deg, file.path(dirs$results, "11_C1vsC3_all_genes.csv"), row.names = FALSE)
message("差异基因表: ", nrow(deg), " genes, top5 up: ",
        paste(head(deg$gene[deg$t > 0], 5), collapse = ", "))

# ---- B1: 火山图（标注 top 基因）----
top_up   <- deg %>% filter(t > 0) %>% slice_max(t, n = 8)
top_down <- deg %>% filter(t < 0) %>% slice_min(t, n = 8)
label_df <- bind_rows(top_up, top_down)
pB1 <- ggplot(deg, aes(t, logFC)) +
  geom_point(alpha = .25, size = 1) +
  geom_point(data = label_df, aes(color = t > 0), size = 2.2) +
  geom_text(data = label_df, aes(label = gene), size = 2.8, vjust = -0.8) +
  scale_color_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#4DBBD5")) +
  theme_paper() + theme(legend.position = "none") +
  labs(title = sprintf("C1 vs C3 differential expression (Welch t, n=%d vs %d)", n1, n3),
       x = "t-statistic (C1 higher > 0)", y = "log2 fold change")

# ---- B2: 16 通路 C1 vs C3（用现有评分, Wilcoxon）----
pw <- bind_rows(lapply(names(gs), function(nm) {
  w <- suppressWarnings(wilcox.test(df[[nm]][df$cluster == "C1"],
                                    df[[nm]][df$cluster == "C3"]))
  data.frame(pathway = nm,
             median_C1 = median(df[[nm]][df$cluster == "C1"], na.rm = TRUE),
             median_C3 = median(df[[nm]][df$cluster == "C3"], na.rm = TRUE),
             p = w$p.value)
})) %>% mutate(dir = ifelse(median_C1 > median_C3, "C1-high", "C3-high"),
              label = sprintf("%s
p=%s", pathway, fmt_p(p)))
write.csv(pw, file.path(dirs$results, "11_C1vsC3_pathways.csv"), row.names = FALSE)
pB2 <- ggplot(pw, aes(reorder(pathway, median_C1 - median_C3), median_C1 - median_C3, fill = dir)) +
  geom_col(width = .7, color = "black", linewidth = .3) +
  scale_fill_manual(values = c("C1-high" = "#E64B35", "C3-high" = "#4DBBD5")) +
  coord_flip() + theme_paper() + theme(legend.position = "none") +
  labs(title = "Pathway score differences: C1 vs C3",
       subtitle = "median z-score difference (Wilcoxon rank-sum)",
       x = NULL, y = "C1 minus C3")
figB <- (pB1 | pB2) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))
ggsave(file.path(dirs$figures, "Figure11_C1vsC3.pdf"), figB, width = 13, height = 6)
message("Figure11_C1vsC3.pdf saved")
message("11 all done")
