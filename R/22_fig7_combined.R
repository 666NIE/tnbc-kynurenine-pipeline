# 22_fig7_combined.R -- 图7五面板合并版（单细胞患者级 + 空间转录组）+ 精简版图S7
source("R/00_config.R")
suppressPackageStartupMessages({library(tidyverse); library(gridExtra)})

m <- read.csv(path.expand(paths$scrna_meta), check.names = FALSE, row.names = 1)
agg <- m %>% filter(Cell_class %in% c("Myeloid", "T_cell")) %>%
  group_by(orig.ident, Cell_class) %>%
  summarise(Kyn = mean(Kyn_score, na.rm = TRUE), Exh = mean(Exhaustion_score, na.rm = TRUE),
            Cyto = mean(Cytotoxic_score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Cell_class, values_from = c(Kyn, Exh, Cyto)) %>%
  filter(!is.na(Kyn_Myeloid), !is.na(Exh_T_cell))
cor_lab <- function(x, y) {
  ct <- suppressWarnings(cor.test(x, y, method = "spearman"))
  sprintf("Spearman rho = %.3f\n%s", unname(ct$estimate), fmt_p(ct$p.value)) }
th <- theme_paper()
pA <- ggplot(agg, aes(Kyn_Myeloid, Exh_T_cell)) + geom_point(size = 2.2, alpha = .8, color = "#555555") +
  geom_smooth(method = "lm", se = TRUE, color = "black") + th +
  labs(x = "Myeloid kynurenine score (patient mean)", y = "T-cell exhaustion score",
       title = "Myeloid KP vs T-cell exhaustion", caption = cor_lab(agg$Kyn_Myeloid, agg$Exh_T_cell))
pB <- ggplot(agg, aes(Kyn_Myeloid, Cyto_T_cell)) + geom_point(size = 2.2, alpha = .8, color = "#555555") +
  geom_smooth(method = "lm", se = TRUE, color = "black") + th +
  labs(x = "Myeloid kynurenine score (patient mean)", y = "T-cell cytotoxicity score",
       title = "Myeloid KP vs T-cell cytotoxicity", caption = cor_lab(agg$Kyn_Myeloid, agg$Cyto_T_cell))

spots <- readRDS(file.path(dirs$results, "21_spots.rds"))
tnbc <- c("1142243F", "1160920F", "CID4465", "CID44971")
sp <- spots %>% filter(sample %in% tnbc, !is.na(Classification), Classification != "", Classification != "Artefact")
grad <- function(hi) scale_color_gradient2(low = "grey85", mid = "white", high = hi, midpoint = 0)
smap <- function(smp, score, hi, ttl) ggplot(sp %>% filter(sample == smp), aes(x, y, color = .data[[score]])) +
  geom_point(size = 1.0) + grad(hi) + scale_y_reverse() + coord_fixed() +
  theme_void(base_size = 11) + labs(title = ttl, color = NULL)

pC1 <- smap("1142243F", "KP",    "#B22222", "1142243F — KP score")
pC2 <- smap("1142243F", "Tcell", "#1F4E8C", "1142243F — T-cell score")

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
pD <- ggplot(sp, aes(reorder(class_short, KP, median), KP, fill = sample)) +
  geom_boxplot(outlier.shape = NA, alpha = .7) + coord_flip() + theme_bw(base_size = 10) +
  labs(x = NULL, y = "KP spot score", title = "KP score across pathological compartments", fill = "Section")

ra <- suppressWarnings(cor.test(sp$KP, sp$Tcell, method = "spearman"))
pE <- ggplot(sp, aes(KP, Tcell, color = sample)) + geom_point(size = .6, alpha = .5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") + theme_bw(base_size = 11) +
  labs(title = sprintf("Spot-level KP vs T-cell score (pooled rho = %.3f, %s)",
                       unname(ra$estimate), fmt_p(ra$p.value)),
       x = "KP score", y = "T-cell score", color = "Section")

ggsave(file.path(dirs$figures, "Figure7_combined.pdf"),
       arrangeGrob(pA, pB, pC1, pC2, pD, pE, ncol = 2,
                   layout_matrix = rbind(c(1,2), c(3,4), c(5,6))), width = 13, height = 13)

pF1 <- smap("CID44971", "KP",    "#B22222", "CID44971 — KP score")
pF2 <- smap("CID44971", "Tcell", "#1F4E8C", "CID44971 — T-cell score")
ggsave(file.path(dirs$figures, "FigureS7_spatial_v2.pdf"),
       arrangeGrob(pF1, pF2, pD, ncol = 2, layout_matrix = rbind(c(1,2), c(3,3))),
       width = 13, height = 9)
message("22 done")
