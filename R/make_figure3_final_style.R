# make_figure3_final_style.R
# Correct current-manuscript Figure 3: immune-desert / immunosuppression panels.
# A: Treg, Macrophage, Checkpoint, Exhaustion scores, C1 vs non-C1
# B: Kynurenine vs combined immunosuppression score
# C: Combined immunosuppression score, C1 vs non-C1
# Style: theme_bw(base_size=9), small points, unified C1 red / non-C1 cyan.
source("R/00_config.R")
suppressPackageStartupMessages({ library(tidyverse); library(ggpubr); library(patchwork) })

coh <- readRDS(file.path(dirs$results, "01_cohort.rds"))
scr <- readRDS(file.path(dirs$results, "02_scores.rds"))
cl  <- readRDS(file.path(dirs$results, "03_cluster.rds"))
scores <- scr$scores; cluster <- cl$cluster

df <- scores %>% mutate(pid = patient, cluster = as.character(cluster[patient]))
df$grp <- factor(ifelse(df$cluster == "C1", "C1", "non-C1"), levels = c("C1", "non-C1"))
df$Immunosuppressive <- rowMeans(as.matrix(df[, c("Treg", "Macrophage", "Checkpoint", "Exhaustion")]))

sub_cols <- c("C1" = "#E64B35", "C2" = "#4DBBD5", "C3" = "#00A087")
pal2 <- c("C1" = "#E64B35", "non-C1" = "#4DBBD5")

markers <- c("Checkpoint", "Exhaustion", "Macrophage", "Treg")
long_df <- df %>% pivot_longer(all_of(markers), names_to = "marker", values_to = "score")
pvals <- long_df %>% group_by(marker) %>%
  summarise(p = suppressWarnings(wilcox.test(score ~ grp)$p.value), .groups = "drop") %>%
  mutate(label = sprintf("%.1e", p))

p3a <- ggplot(long_df, aes(grp, score, fill = grp)) +
  geom_boxplot(width = .55, alpha = .8, outlier.shape = NA, linewidth = .5) +
  geom_jitter(aes(color = grp), width = .12, size = .55, alpha = .7) +
  facet_wrap(~ factor(marker, levels = markers), nrow = 1, scales = "free_y") +
  scale_fill_manual(values = pal2) +
  scale_color_manual(values = pal2) +
  geom_text(data = pvals, aes(x = 1.5, y = Inf, label = label), inherit.aes = FALSE,
            vjust = 1.4, size = 2.4) +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        strip.text = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 9)) +
  labs(title = "Immunosuppressive signatures", x = NULL, y = "Pathway z-score")

cor_lab <- function(x, y) {
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", use = "complete.obs"))
  sprintf("Spearman rho = %.3f\n%s", unname(ct$estimate), format.pval(ct$p.value, digits = 2, eps = 0.001))
}
p3b <- ggplot(df, aes(Kynurenine, Immunosuppressive, color = cluster)) +
  geom_point(size = .7, alpha = .85) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = sub_cols) +
  theme_bw(base_size = 9) +
  labs(title = "Kynurenine vs immunosuppression", x = "Kynurenine score",
       y = "Combined immunosuppression score", color = "Subtype")

p3c <- ggplot(df, aes(grp, Immunosuppressive, fill = grp)) +
  geom_boxplot(width = .55, alpha = .8, outlier.shape = NA, linewidth = .6) +
  geom_jitter(aes(color = grp), width = .13, size = .6, alpha = .75) +
  scale_fill_manual(values = pal2) + scale_color_manual(values = pal2) +
  ggpubr::stat_compare_means(comparisons = list(c("C1", "non-C1")), label = "p.format", size = 3,
                              label.y = max(df$Immunosuppressive, na.rm = TRUE) + 0.35) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  theme_bw(base_size = 9) +
  labs(title = "Combined immunosuppression score", x = NULL, y = "score (z)") +
  theme(legend.position = "none")

fig3 <- p3a / (p3b | p3c) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 12))
dir.create("figures", showWarnings = FALSE)
ggsave("figures/Figure3.png", fig3, width = 170, height = 140, units = "mm", dpi = 300, bg = "white")
message("Final Figure3 written: figures/Figure3.png")
