# 08_fig9.R — Figure 9: TIDE (4 panels, A-D merged)
source("R/00_config.R")
suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })
coh     <- readRDS(file.path(dirs$results, "01_cohort.rds"))
scr     <- readRDS(file.path(dirs$results, "02_scores.rds"))
cl      <- readRDS(file.path(dirs$results, "03_cluster.rds"))
scores  <- scr$scores; cluster <- cl$cluster
sub_cols <- c(C1 = "#F8766D", C2 = "#619CFF", C3 = "#00BA38")
cor_lab <- function(x, y) {
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", use = "complete.obs"))
  sprintf("Spearman rho = %.3f, %s", unname(ct$estimate), fmt_p(ct$p.value))
}

rdf <- read.csv(path.expand(paths$response), check.names = FALSE)
cn <- colnames(rdf)
id_c   <- grep("patient|sample|^ID$", cn, ignore.case = TRUE, value = TRUE)[1]
tide_c <- grep("TIDE", cn, ignore.case = TRUE, value = TRUE)[1]
resp_c <- intersect("Responder", cn)[1]
rdf$patient <- patient_id(as.character(rdf[[id_c]]))
rdf$TIDE <- suppressWarnings(as.numeric(rdf[[tide_c]]))
if (!is.na(resp_c)) {
  rv <- as.character(rdf[[resp_c]])
  rdf$Responder <- ifelse(grepl("true|responder|benefit|yes", rv, ignore.case = TRUE),
                          "Responder", "Non-responder")
} else {
  rdf$Responder <- ifelse(rdf$TIDE < median(rdf$TIDE, na.rm = TRUE), "Responder", "Non-responder")
}
rdf$Responder <- factor(rdf$Responder, levels = c("Non-responder", "Responder"))

f9 <- rdf %>% filter(patient %in% coh$patients, !is.na(TIDE)) %>%
  mutate(cluster = as.character(cluster[match(patient, names(cluster))])) %>%
  filter(!is.na(cluster)) %>%
  left_join(scores %>% select(patient, Kynurenine), by = "patient")
tab <- table(f9$cluster, f9$Responder)
message("Fig9 应答表:"); print(tab)
message("Fisher: ", fmt_p(fisher.test(tab)$p.value))

prop_df <- f9 %>% group_by(cluster, Responder) %>% summarise(n = n(), .groups = "drop") %>%
  group_by(cluster) %>% mutate(total = sum(n), rate = n / total,
    lo = binom.test(n, total)$conf.int[1], hi = binom.test(n, total)$conf.int[2]) %>%
  filter(Responder == "Responder")
pA <- ggplot(prop_df, aes(cluster, rate, fill = cluster)) +
  geom_col(width = .6, alpha = .85) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .15) +
  geom_text(aes(label = sprintf("%.1f%%
(n=%d/%d)", 100 * rate, n, total)),
            vjust = -0.3, size = 3.2) +
  scale_fill_manual(values = sub_cols) + theme_paper() + theme(legend.position = "none") +
  ylim(0, 1.2) +
  labs(title = "TIDE-predicted response rate",
       subtitle = paste0("Fisher exact, ", fmt_p(fisher.test(tab)$p.value)),
       x = NULL, y = "Responder proportion")
pB <- ggplot(f9, aes(cluster, TIDE, fill = cluster)) +
  geom_boxplot(alpha = .8, outlier.shape = NA) +
  geom_jitter(width = .15, alpha = .5, size = 1.5) +
  scale_fill_manual(values = sub_cols) + theme_paper() + theme(legend.position = "none") +
  labs(title = "TIDE score", x = NULL, y = "TIDE",
       caption = fmt_p(kw_p(TIDE ~ cluster, f9)))
pC <- ggplot(f9, aes(Kynurenine, TIDE, color = cluster)) +
  geom_point(size = 2.5, alpha = .8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = sub_cols) + theme_paper() +
  labs(title = "Kynurenine vs TIDE", x = "Kynurenine score", y = "TIDE score",
       color = "Subtype", caption = cor_lab(f9$Kynurenine, f9$TIDE))

# ---- Panel D: TIDE components (merged from former Fig S1) ----
comp <- read.csv(file.path(dirs$results, "10_tide_components.csv"))
comp$dir <- comp$rho > 0
pD <- ggplot(comp, aes(reorder(component, rho), rho, fill = dir)) +
  geom_col(width = .7, color = "black", linewidth = .3) +
  geom_text(aes(label = sprintf("rho=%.2f
%s", rho, fmt_p(p))),
            size = 2.8, hjust = ifelse(comp$rho > 0, -0.05, 1.08)) +
  scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#4DBBD5")) +
  coord_flip() + theme_paper() + theme(legend.position = "none") +
  ylim(-0.6, 0.78) +
  labs(title = "TIDE components", x = NULL, y = "Spearman rho")

fig9 <- (pA | pB | pC | pD) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))
ggsave(file.path(dirs$figures, "Figure9_TIDE.pdf"), fig9, width = 17, height = 5)
write.csv(f9 %>% select(patient, cluster, TIDE, Responder, Kynurenine),
          file.path(dirs$results, "08_tide.csv"), row.names = FALSE)
message("08 done: Figure9_TIDE.pdf 四联图（A-D 合并）已出")
