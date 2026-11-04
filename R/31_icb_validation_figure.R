# 31_icb_validation_figure.R — I-SPY2 ICB validation figure (Supplementary Figure S9)
source("R/00_config.R")
suppressPackageStartupMessages({ library(tidyverse) })

d <- read_csv(file.path(dirs$results, "30_ispY2_annot_scores.csv"),
              col_types = cols(tnbc = col_logical(), trap = col_logical(), pcr = col_integer()))
ARMS <- c("Paclitaxel + Pembrolizumab", "Paclitaxel")

cp_ci <- function(k, n, conf = 0.95) {
  a <- (1 - conf) / 2
  tibble(lo = if (k == 0) 0 else qbeta(a, k, n - k + 1),
         hi = if (k == n) 1 else qbeta(1 - a, k + 1, n - k))
}

pa <- d %>% filter(tnbc, arm %in% ARMS) %>%
  group_by(arm, trap) %>% summarise(k = sum(pcr), n = n(), .groups = "drop") %>%
  mutate(rate = 100 * k / n,
         ci = purrr::map2(k, n, cp_ci)) %>%
  tidyr::unnest_wider(ci) %>%
  mutate(lo = 100 * lo, hi = 100 * hi,
         arm  = factor(arm, levels = ARMS),
         trap = factor(trap, levels = c(FALSE, TRUE),
                       labels = c("other tumors", "CD8-low/IDO1-low (trap-like)")))

arm_n <- pa %>% group_by(arm) %>% summarise(N = sum(n), .groups = "drop")
arm_lab <- c("Paclitaxel + Pembrolizumab" =
               paste0("Pembrolizumab\n(TNBC, n=", arm_n$N[arm_n$arm == ARMS[1]], ")"),
             "Paclitaxel" =
               paste0("Paclitaxel\n(TNBC, n=", arm_n$N[arm_n$arm == ARMS[2]], ")"))

pembro <- d %>% filter(tnbc, arm == ARMS[1])
ft <- fisher.test(table(pembro$trap, pembro$pcr))

tn <- d %>% filter(tnbc, arm %in% ARMS) %>%
  mutate(arm2 = ifelse(arm == ARMS[1], "pembrolizumab", "paclitaxel"))
ct <- suppressWarnings(cor.test(tn$Kynurenine, tn$T_cell, method = "spearman"))
labB <- sprintf("Spearman rho = %.3f\np = %.1e", unname(ct$estimate), ct$p.value)

th <- theme_bw(base_size = 8, base_family = "Arial") +
  theme(plot.title   = element_text(face = "bold", size = 8.5),
        legend.position = "top", legend.text = element_text(size = 7),
        legend.margin = margin(b = -2), legend.box.spacing = unit(2, "pt"),
        panel.grid   = element_line(colour = "grey90", linewidth = 0.3),
        axis.text    = element_text(size = 7.5))

pA <- ggplot(pa, aes(arm, rate, fill = trap)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.58,
           colour = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.72),
                width = 0.15, linewidth = 0.4) +
  geom_text(aes(y = hi + 4, label = paste0(k, "/", n)),
            position = position_dodge(width = 0.72), size = 2.6) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_fill_manual(values = c("other tumors" = "#619CFF",
                               "CD8-low/IDO1-low (trap-like)" = "#F8766D")) +
  scale_x_discrete(labels = arm_lab) +
  scale_y_continuous(limits = c(0, 112), expand = c(0, 0)) +
  labs(title = "A | pCR by CD8-low/IDO1-low profile and arm",
       subtitle = sprintf("pembrolizumab arm, Fisher's exact p = %.3f", ft$p.value),
       x = NULL, y = "pCR rate (%)", fill = NULL) + th +
  theme(plot.subtitle = element_text(size = 7.5))

pB <- ggplot(tn, aes(Kynurenine, T_cell)) +
  geom_point(aes(colour = arm2), size = 1.2, alpha = 0.85, stroke = 0) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = 0.7) +
  annotate("text", x = -Inf, y = -Inf, label = labB, size = 2.6,
           hjust = -0.05, vjust = -0.6) +
  scale_colour_manual(values = c(pembrolizumab = "#2C7FB8", paclitaxel = "#BDBDBD")) +
  labs(title = "B | Kyn-immune coupling in TNBC",
       x = "Kynurenine score (z)", y = "T-cell score (z)", colour = NULL) + th

g <- patchwork::wrap_plots(pA, pB, nrow = 1) + patchwork::plot_layout(widths = c(1.15, 1))
dir.create("figures", showWarnings = FALSE)
ggsave("figures/Figure_ICB_validation.pdf",  g, width = 17.6, height = 7.4, units = "cm")
ggsave("figures/Figure_ICB_validation.tiff", g, width = 17.6, height = 7.4, units = "cm",
       dpi = 300, compression = "lzw")
message("31 done")
