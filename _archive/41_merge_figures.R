# 41_merge_figures.R — 原生重建 Figure 2 / 3 / 6（两组口径，终版，大写 panel 标号）
# 运行：在 pipeline 根目录 source("R/41_merge_figures.R")
source("R/00_config.R")
source("R/01_data.R"); source("R/02_scores.R"); source("R/03_cluster.R")
suppressPackageStartupMessages({
  library(tidyverse); library(ggpubr); library(patchwork); library(ggsci)
  library(pdftools); library(grid); library(png)
})

dir.create("figures", showWarnings = FALSE)

## ---------- 数据准备 ----------
df <- as.data.frame(scores)
df$pid     <- coh$patients
df$cluster <- as.character(cluster[df$pid])
df$grp     <- factor(ifelse(df$cluster == "C1", "C1", "non-C1"),
                     levels = c("C1", "non-C1"))
immuno4 <- c("Treg", "Macrophage", "Checkpoint", "Exhaustion")
df$Immunosuppressive_composite <- rowMeans(df[, immuno4], na.rm = TRUE)

pal2 <- c("C1" = "#E64B35", "non-C1" = "#4DBBD5")

mk2 <- function(y, ylab)
  ggboxplot(df, x = "grp", y = y, color = "grp", add = "jitter", palette = pal2,
            ylab = ylab, xlab = NULL) +
    stat_compare_means(comparisons = list(c("C1", "non-C1")), label = "p.format") +
    theme_classic(base_size = 12) + theme(legend.position = "none")

corr_p <- function(y, ylab)
  ggplot(df, aes(Kynurenine, .data[[y]])) +
    geom_point(aes(color = cluster), size = 1.5, alpha = .7) +
    geom_smooth(method = "lm", color = "black", se = TRUE) +
    scale_color_npg() + theme_classic(base_size = 12) +
    labs(x = "Kynurenine score", y = ylab, color = "Subtype")

## ---------- Figure 2 ----------
p2a <- mk2("Kynurenine", "Kynurenine ssGSEA score (z)")
p2b <- mk2("IDO1", "IDO1 expression (log2)")
p2c <- corr_p("Cytotoxicity", "Cytotoxicity score")
p2d <- corr_p("T_cell", "T-cell infiltration score")

## 2e: swap-in panel —— 渲染旧图→写临时PNG→读回→裁底部（起点 0.68，可调）
pg  <- pdf_render_page("figures/Figure2_trap.pdf", dpi = 150, numeric = TRUE)
tmp <- tempfile(fileext = ".png")
png::writePNG(pg / 255, tmp)
im  <- png::readPNG(tmp)
h2  <- dim(im)[1]
img <- im[as.integer(h2 * 0.68):h2, , , drop = FALSE]
p2e <- wrap_elements(rasterGrob(as.raster(img), interpolate = TRUE))

fig2 <- (p2a | p2b) / (p2c | p2d) / p2e +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 14))
ggsave("figures/Figure2.png", fig2, width = 170, height = 200, units = "mm", dpi = 300, bg = "white")
ggsave("figures/Figure2_twoGroup_native.pdf", fig2, width = 11, height = 13)

## ---------- Figure 3 ----------
p3a <- df %>%
  pivot_longer(all_of(immuno4), names_to = "signature", values_to = "score") %>%
  ggboxplot(x = "grp", y = "score", color = "grp", add = "jitter", palette = pal2,
            facet.by = "signature", nrow = 1, ylab = "Pathway z-score", xlab = NULL) +
  stat_compare_means(comparisons = list(c("C1", "non-C1")), label = "p.format") +
  theme_classic(base_size = 12) + theme(legend.position = "none")
p3b <- corr_p("Immunosuppressive_composite", "Combined immunosuppression score")
p3c <- mk2("Immunosuppressive_composite", "Composite immunosuppression z-score")

fig3 <- p3a / (p3b | p3c) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 14))
ggsave("figures/Figure3_twoGroup_native.pdf", fig3, width = 13, height = 8.5)

## ---------- Figure 6 ----------
tide <- read.csv("results/08_tide.csv", check.names = FALSE)
samp_col <- names(tide)[1]
tid  <- tide %>% mutate(pid = sub("-01$", "", sub("\\.", "-", as.character(.data[[samp_col]])))) %>%
  filter(pid %in% df$pid)
tid$cluster <- df$cluster[match(tid$pid, df$pid)]
tid$grp <- factor(ifelse(tid$cluster == "C1", "C1", "non-C1"),
                  levels = c("C1", "non-C1"))
tcol  <- grep("TIDE", names(tid), value = TRUE)[1]
comps <- intersect(c("Dysfunction", "Exclusion", "MDSC", "CAF",
                     "TAM M2", "TAM_M2", "TAM.M2"), names(tid))

rr <- tid %>% group_by(cluster) %>%
  summarise(resp = sum(.data[[tcol]] < 0), n = n(), .groups = "drop") %>%
  rowwise() %>% mutate(ci = list(binom.test(resp, n)$conf.int)) %>%
  unnest_wider(ci, names_sep = "_") %>% rename(lo = ci_1, hi = ci_2)
p6a <- ggplot(rr, aes(cluster, 100 * resp / n, fill = cluster)) +
  geom_col(width = .6) +
  geom_errorbar(aes(ymin = 100 * lo, ymax = 100 * hi), width = .15) +
  geom_text(aes(label = paste0(round(100 * resp / n, 1), "%\n(", resp, "/", n, ")")),
            vjust = -0.4, size = 3.2) +
  scale_fill_npg() + ylim(0, 60) + theme_classic(base_size = 12) +
  labs(x = NULL, y = "Predicted response rate (%)", fill = "Subtype") +
  theme(legend.position = "none")

tid$TIDE_num <- as.numeric(tid[[tcol]])
p6b <- ggboxplot(tid, x = "grp", y = "TIDE_num", color = "grp", add = "jitter",
                 palette = pal2, ylab = "TIDE score", xlab = NULL) +
  stat_compare_means(comparisons = list(c("C1", "non-C1")), label = "p.format") +
  theme_classic(base_size = 12) + theme(legend.position = "none")

kx  <- df$Kynurenine[match(tid$pid, df$pid)]
p6c <- ggplot(tid, aes(kx, TIDE_num)) +
  geom_point(aes(color = cluster), size = 1.5, alpha = .7) +
  geom_smooth(method = "lm", color = "black") +
  scale_color_npg() + theme_classic(base_size = 12) +
  labs(x = "Kynurenine score", y = "TIDE score", color = "Subtype")

if (length(comps) >= 3) {
  cmp <- map_dfr(comps, function(cp) {
    cc2 <- suppressWarnings(cor.test(kx, as.numeric(tid[[cp]]), method = "spearman"))
    tibble(component = cp, rho = unname(cc2$estimate), p = cc2$p.value)
  })
} else {
  cmp <- tibble(
    component = c("TAM M2", "Exclusion", "CAF", "MDSC", "Dysfunction"),
    rho       = c(0.545, 0.408, 0.234, 0.205, -0.223),
    p         = c(1.2e-11, 3.3e-06, 0.009, 0.023, 0.013))
}
p6d <- ggplot(cmp, aes(reorder(component, rho), rho, fill = rho > 0)) +
  geom_col(width = .6) + coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.28))) +
  scale_fill_manual(values = c("#4DBBD5", "#E64B35")) +
  geom_text(aes(label = sprintf("rho=%.2f\np=%.2g", rho, p)), hjust = -0.1, size = 3) +
  theme_classic(base_size = 12) + labs(x = NULL, y = "Spearman rho", fill = NULL) +
  theme(legend.position = "none")

fig6 <- (p6a | p6b) / (p6c | p6d) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 14))
ggsave("figures/Figure6_twoGroup_native.pdf", fig6, width = 11, height = 9)

cat("done -> figures/Figure2/3/6_twoGroup_native.pdf\n")
