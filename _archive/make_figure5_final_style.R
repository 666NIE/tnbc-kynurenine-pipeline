# make_figure5_final_style.R
# Final Figure 5: patient-level myeloid KP vs T-cell exhaustion/cytotoxicity + spatial panels.
# Keeps the original colored boxplot style for compartments and uses short titles.
source("R/00_config.R")
suppressPackageStartupMessages({ library(tidyverse); library(gridExtra) })

meta_file <- "data-raw/manual_downloads/fig7_output/Fig7_metadata_with_scores.csv"
if (!file.exists(meta_file)) stop("Missing: ", meta_file)
m <- read.csv(meta_file, check.names = FALSE, row.names = 1)
agg <- m %>% filter(Cell_class %in% c("Myeloid", "T_cell")) %>%
  group_by(orig.ident, Cell_class) %>%
  summarise(Kyn = mean(Kyn_score, na.rm = TRUE), Exh = mean(Exhaustion_score, na.rm = TRUE),
            Cyto = mean(Cytotoxic_score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Cell_class, values_from = c(Kyn, Exh, Cyto)) %>%
  filter(!is.na(Kyn_Myeloid), !is.na(Exh_T_cell))
cor_lab <- function(x, y) {
  ct <- suppressWarnings(cor.test(x, y, method = "spearman"))
  sprintf("Spearman rho = %.3f\np = %s", unname(ct$estimate), format.pval(ct$p.value, digits = 2, eps = 0.001))
}
th <- theme_bw(base_size = 9)
pA <- ggplot(agg, aes(Kyn_Myeloid, Exh_T_cell)) +
  geom_point(size = .8, alpha = .8, color = "grey30") +
  geom_smooth(method = "lm", se = TRUE, color = "black") + th +
  labs(x = "Myeloid kynurenine score (patient mean)", y = "T-cell exhaustion score",
       title = "Myeloid KP vs T-cell exhaustion", caption = cor_lab(agg$Kyn_Myeloid, agg$Exh_T_cell))
pB <- ggplot(agg, aes(Kyn_Myeloid, Cyto_T_cell)) +
  geom_point(size = .8, alpha = .8, color = "grey30") +
  geom_smooth(method = "lm", se = TRUE, color = "black") + th +
  labs(x = "Myeloid kynurenine score (patient mean)", y = "T-cell cytotoxicity score",
       title = "Myeloid KP vs T-cell cytotoxicity", caption = cor_lab(agg$Kyn_Myeloid, agg$Cyto_T_cell))

spots <- readRDS(file.path(dirs$results, "21_spots.rds"))
tnbc <- c("1142243F", "1160920F", "CID4465", "CID44971")
sp <- spots %>% filter(sample %in% tnbc, !is.na(Classification), Classification != "", Classification != "Artefact")
short_lab <- c("Invasive cancer + stroma + lymphocytes" = "Cancer+Stroma+Lymph",
               "Invasive cancer + lymphocytes" = "Cancer+Lymph",
               "Invasive cancer + stroma" = "Cancer+Stroma",
               "Lymphocytes" = "Lymphocytes", "TLS" = "TLS", "DCIS" = "DCIS",
               "Stroma" = "Stroma", "Stroma + adipose tissue" = "Stroma+Adipose",
               "Adipose tissue" = "Adipose", "Necrosis" = "Necrosis",
               "Normal + stroma + lymphocytes" = "Normal+Stroma+Lymph",
               "Normal glands + lymphocytes" = "Normal glands+Lymph",
               "Normal duct" = "Normal duct",
               "Cancer trapped in lymphocyte aggregation" = "Cancer trapped in lymphoid agg.")
sp$class_short <- short_lab[sp$Classification]
grad <- function(hi) scale_color_gradient2(low = "grey85", mid = "white", high = hi, midpoint = 0)
smap <- function(smp, score, hi, ttl) {
  dd <- sp %>% filter(sample == smp)
  hi_dd <- dd[dd[[score]] > 2, , drop = FALSE]
  ggplot(dd, aes(x, y, color = .data[[score]])) +
    geom_point(size = .18, alpha = .55) +
    geom_point(data = hi_dd, aes(x, y), color = hi, size = .45, alpha = .95) +
    grad(hi) + scale_y_reverse() + coord_fixed() +
    theme_void(base_size = 9) + labs(title = ttl, color = NULL)
}
pC1 <- smap("1142243F", "KP", "#B22222", "1142243F — KP score")
pC2 <- smap("1142243F", "Tcell", "#1F4E8C", "1142243F — T-cell score")
pD <- ggplot(sp, aes(reorder(class_short, KP, median), KP, fill = sample)) +
  geom_boxplot(outlier.shape = NA, alpha = .85, width = .55, linewidth = .4) +
  coord_flip() + theme_bw(base_size = 8) +
  labs(title = "KP score across compartments", x = NULL, y = "KP spot score", fill = "Section") +
  theme(plot.title.position = "plot")
ra <- suppressWarnings(cor.test(sp$KP, sp$Tcell, method = "spearman"))
pE <- ggplot(sp, aes(KP, Tcell, color = sample)) +
  geom_point(size = .18, alpha = .45) +
  geom_smooth(method = "lm", se = TRUE, color = "black") + theme_bw(base_size = 9) +
  labs(title = "Spot-level KP vs T-cell association",
       subtitle = sprintf("rho = %.3f, p = %s", unname(ra$estimate), format.pval(ra$p.value, digits = 2, eps = 0.001)),
       x = "KP score", y = "T-cell score", color = "Section") +
  theme(plot.title.position = "plot")

fig5 <- arrangeGrob(pA, pB, pC1, pC2, pD, pE, ncol = 2,
                    layout_matrix = rbind(c(1,2), c(3,4), c(5,6)))
dir.create("figures", showWarnings = FALSE)
ggsave("figures/Figure5.png", fig5, width = 170, height = 225, units = "mm", dpi = 300, bg = "white")
message("Final Figure5 written: figures/Figure5.png")
