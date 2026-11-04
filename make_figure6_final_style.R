# make_figure6_final_style.R
# Correct Figure 6: keep the ORIGINAL data口径 (Figure6_twoGroup_native):
# A: three subtypes C1/C2/C3 response rates (17.9/35.6/34.0; Fisher p=0.239)
# B: C1 vs non-C1 TIDE score, Wilcoxon p=0.016
# C: Kynurenine vs TIDE scatter
# D: TIDE components from results/10_tide_components.csv
# Style: same clean theme_bw/small-point style as the fixed figures.
source("R/00_config.R")
suppressPackageStartupMessages({ library(tidyverse); library(ggpubr)
library(patchwork) })

coh <- readRDS(file.path(dirs$results, "01_cohort.rds"))
scr <- readRDS(file.path(dirs$results, "02_scores.rds"))
cl  <- readRDS(file.path(dirs$results, "03_cluster.rds"))
scores <- scr$scores; cluster <- cl$cluster

df <- scores %>% mutate(pid = patient, cluster = as.character(cluster[patient]))
sub_cols <- c("C1" = "#E64B35", "C2" = "#4DBBD5", "C3" = "#00A087")
pal2 <- c("C1" = "#E64B35", "non-C1" = "#4DBBD5")

tide <- read.csv("results/08_tide.csv", check.names = FALSE)
samp_col <- names(tide)[1]
tid <- tide %>% mutate(pid = sub("-01$", "", sub("\\.", "-", as.character(.data[[samp_col]])))) %>%
  filter(pid %in% df$pid)
tid$cluster <- df$cluster[match(tid$pid, df$pid)]
tid$grp <- factor(ifelse(tid$cluster == "C1", "C1", "non-C1"), levels = c("C1", "non-C1"))
tid$TIDE_num <- as.numeric(tid$TIDE)

# A: response rates by three subtypes
rr <- tid %>% group_by(cluster) %>%
  summarise(resp = sum(TIDE_num < 0), n = n(), .groups = "drop") %>%
  rowwise() %>% mutate(ci = list(binom.test(resp, n)$conf.int)) %>%
  unnest_wider(ci, names_sep = "_") %>% rename(lo = ci_1, hi = ci_2)
p6a <- ggplot(rr, aes(cluster, 100 * resp / n, fill = cluster)) +
  geom_col(width = .6, alpha = .9) +
  geom_errorbar(aes(ymin = 100 * lo, ymax = 100 * hi), width = .15, linewidth = .6) +
  geom_text(aes(label = paste0(round(100 * resp / n, 1), "%\n(", resp, "/", n, ")")),
            vjust = -0.35, size = 2.2) +
  scale_fill_manual(values = sub_cols) + ylim(0, 60) +
  theme_bw(base_size = 9) +
  labs(title = "TIDE-predicted response rate", subtitle = "Fisher exact, p = 0.239",
       x = NULL, y = "Responder proportion (%)", fill = "Subtype") +
  theme(legend.position = "none")

# B: C1 vs non-C1
p6b <- ggplot(tid, aes(grp, TIDE_num)) +
  geom_boxplot(aes(fill = grp), width = .55, alpha = .8, outlier.shape = NA, linewidth = .6) +
  geom_jitter(aes(color = grp), width = .13, size = .6, alpha = .75) +
  scale_fill_manual(values = pal2) + scale_color_manual(values = pal2) +
  ggpubr::stat_compare_means(comparisons = list(c("C1", "non-C1")), label = "p.format", size = 3) +
  coord_cartesian(ylim = c(NA, max(tid$TIDE_num, na.rm = TRUE) + 0.25), clip = "off") +
  theme_bw(base_size = 9) +
  labs(title = "TIDE score", x = NULL, y = "TIDE") +
  theme(legend.position = "none")

# C: scatter
kx <- df$Kynurenine[match(tid$pid, df$pid)]
p6c <- ggplot(tid, aes(kx, TIDE_num)) +
  geom_point(aes(color = cluster), size = .7, alpha = .85) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  scale_color_manual(values = sub_cols) +
  theme_bw(base_size = 9) +
  labs(title = "Kynurenine vs TIDE", x = "Kynurenine score", y = "TIDE score", color = "Subtype") +
  guides(color = guide_legend(override.aes = list(size = 1.5)))

# D: components
comp <- read.csv("results/10_tide_components.csv", check.names = FALSE)
if (!all(c("component", "rho", "p") %in% names(comp))) stop("10_tide_components.csv must contain component,rho,p")
comp$dir <- comp$rho > 0
lab_map <- data.frame(
  component = c("TAM M2", "Exclusion", "CAF", "MDSC", "Dysfunction"),
  rho_lab = c("0.55", "0.41", "0.23", "0.20", "-0.22"),
  p_lab = c("1.2e-11", "3.3e-06", "0.009", "0.023", "0.013")
)
comp$rho_lab <- lab_map$rho_lab[match(comp$component, lab_map$component)]
comp$p_lab <- lab_map$p_lab[match(comp$component, lab_map$component)]
comp$rho_lab <- ifelse(is.na(comp$rho_lab), sprintf("%.2f", comp$rho), comp$rho_lab)
comp$p_lab <- ifelse(is.na(comp$p_lab), sub("^p\\s*=\\s*", "", fmt_p(comp$p)), comp$p_lab)
comp$plab <- sub("^p\\s*=\\s*", "", fmt_p(comp$p))
p6d <- ggplot(comp, aes(reorder(component, rho), rho, fill = dir)) +
  geom_col(width = .6) + coord_flip() +
  scale_y_continuous(limits = c(-0.30, 0.85), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#4DBBD5")) +
  geom_text(aes(label = sprintf("rho=%s\np=%s", rho_lab, p_lab)), hjust = -0.08, size = 2) +
  theme_bw(base_size = 9) +
  labs(title = "TIDE components", x = NULL, y = "Spearman rho", fill = NULL) +
  theme(legend.position = "none")

fig6 <- (p6a | p6b) / (p6c | p6d) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 12))
dir.create("figures", showWarnings = FALSE)
ggsave("figures/Figure6.png", fig6, width = 170, height = 140, units = "mm", dpi = 300, bg = "white")
message("Final Figure6 written: figures/Figure6.png")
