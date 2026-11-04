# make_figure2_final_style.R
# Build final Figure 2 in the same visual style as Figure 6:
# - theme_bw(base_size = 9)
# - small solid points
# - unified colors: C1 red, C2 cyan, C3 green; two-group non-C1 uses cyan
# - keeps manuscript statistics: A/B are C1 vs non-C1 (p ~2.3e-11 and ~1.2e-12)
# - E panel is built directly, not cropped from a PDF
source("R/00_config.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(patchwork)
  library(ConsensusClusterPlus)
})

coh <- readRDS(file.path(dirs$results, "01_cohort.rds"))
scr <- readRDS(file.path(dirs$results, "02_scores.rds"))
cl  <- readRDS(file.path(dirs$results, "03_cluster.rds"))
scores <- scr$scores
cluster <- cl$cluster

df <- scores %>% mutate(pid = patient, cluster = as.character(cluster[patient]))
df$grp <- factor(ifelse(df$cluster == "C1", "C1", "non-C1"), levels = c("C1", "non-C1"))

sub_cols <- c("C1" = "#E64B35", "C2" = "#4DBBD5", "C3" = "#00A087")
pal2 <- c("C1" = "#E64B35", "non-C1" = "#4DBBD5")

box2 <- function(y, ylab) {
  ggplot(df, aes(grp, .data[[y]])) +
    geom_boxplot(aes(fill = grp), width = .55, alpha = .8, outlier.shape = NA, linewidth = .6) +
    geom_jitter(aes(color = grp), width = .13, size = .6, alpha = .75) +
    scale_fill_manual(values = pal2) +
    scale_color_manual(values = pal2) +
    ggpubr::stat_compare_means(comparisons = list(c("C1", "non-C1")), label = "p.format", size = 3) +
    coord_cartesian(ylim = c(NA, max(df[[y]], na.rm = TRUE) + 1), clip = "off") +
    theme_bw(base_size = 9) +
    labs(x = NULL, y = ylab) +
    theme(legend.position = "none")
}

corr_p <- function(y, ylab) {
  ggplot(df, aes(Kynurenine, .data[[y]])) +
    geom_point(aes(color = cluster), size = .7, alpha = .85) +
    geom_smooth(method = "lm", color = "black", se = TRUE) +
    scale_color_manual(values = sub_cols) +
    theme_bw(base_size = 9) +
    labs(x = "Kynurenine score", y = ylab, color = "Subtype")
}

p2a <- box2("Kynurenine", "Kynurenine ssGSEA score (z)")
p2b <- box2("IDO1", "IDO1 expression (log2)")
p2c <- corr_p("Cytotoxicity", "Cytotoxicity score")
p2d <- corr_p("T_cell", "T-cell infiltration score")

# E panel: direct swap-in clustering
set.seed(params$seed)
if (!exists("calc_z")) {
  calc_z <- function(genes, expr) {
    genes <- intersect(genes, rownames(expr))
    rowMeans(t(scale(t(expr[genes, , drop = FALSE]))))
  }
}
kyn_abs <- calc_z(params$gene_sets$Kynurenine, coh$expr)
scores_swap <- scores %>% mutate(Kynurenine = kyn_abs)
pm <- as.matrix(t(scores_swap %>% select(-patient) %>% select(where(is.numeric))))
colnames(pm) <- coh$patients
cc2 <- ConsensusClusterPlus(pm, maxK = params$maxK, reps = params$reps, pItem = params$pItem,
                            clusterAlg = "km", distance = "euclidean",
                            plot = "pdf", title = "ConsensusClustering_absKyn_fig2style", seed = params$seed)
cl2 <- cc2[[params$k_final]]$consensusClass
names(cl2) <- colnames(pm)
orig <- as.character(cluster[colnames(pm)])
top_cl <- names(which.max(tapply(kyn_abs, cl2, mean)))
pct_c1 <- round(100 * sum(orig == "C1" & cl2 == top_cl) / sum(orig == "C1"), 1)
e_dat <- data.frame(swap_cluster = factor(paste0("swap-", cl2), levels = paste0("swap-", sort(unique(cl2)))),
                    Kyn_abs = kyn_abs, T_cell = scores$T_cell, Cytotoxicity = scores$Cytotoxicity) %>%
  pivot_longer(c(Kyn_abs, T_cell, Cytotoxicity), names_to = "metric", values_to = "z") %>%
  group_by(swap_cluster, metric) %>% summarise(z = mean(z), .groups = "drop") %>%
  mutate(metric = factor(metric, levels = c("Kyn_abs", "T_cell", "Cytotoxicity"),
                         labels = c("Kyn absolute z", "T-cell score", "Cytotoxicity")))
e_cols <- c("Kyn absolute z" = "#8B1A1A", "T-cell score" = "#3B6FB5", "Cytotoxicity" = "#2E8B57")
p2e <- ggplot(e_dat, aes(swap_cluster, z, fill = metric)) +
  geom_col(position = position_dodge(width = .75), width = .68, alpha = .9) +
  geom_text(aes(label = sprintf("%.2f", z)), position = position_dodge(width = .75), vjust = -0.4, size = 2) +
  geom_hline(yintercept = 0, linewidth = .4) +
  scale_fill_manual(values = e_cols) +
  theme_bw(base_size = 9) +
  labs(title = "Swap-in clustering (absolute z replaces ssGSEA)",
       x = "Swap-in cluster", y = "mean z", fill = NULL,
       caption = sprintf("Original C1 in top absolute-Kyn cluster: %s%%", pct_c1)) +
  theme(legend.position = "bottom")

fig2 <- (p2a | p2b) / (p2c | p2d) / p2e +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

dir.create("figures", showWarnings = FALSE)
ggsave("figures/Figure2.png", fig2, width = 170, height = 200, units = "mm", dpi = 300, bg = "white")
message("Final Figure2 written: figures/Figure2.png")
