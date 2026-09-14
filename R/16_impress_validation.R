# 16_impress_validation.R -- IMPRES consistency check (reviewer Major 6)
# IMPRES: Auslander et al., Nat Med 2018;24:1545-1548. doi:10.1038/s41591-018-0157-9
source("R/00_config.R")
suppressPackageStartupMessages({ library(tidyverse); library(ggplot2); library(gridExtra) })
coh <- readRDS(file.path(dirs$results, "01_cohort.rds"))
scr <- readRDS(file.path(dirs$results, "02_scores.rds"))
cl  <- readRDS(file.path(dirs$results, "03_cluster.rds"))
scores <- scr$scores
if (!"VSIR" %in% rownames(coh$expr) && "C10orf54" %in% rownames(coh$expr)) {
  rownames(coh$expr)[rownames(coh$expr) == "C10orf54"] <- "VSIR"
  message("alias applied: C10orf54 -> VSIR")
}
pairs <- list(
  c("PDCD1","TNFRSF4"), c("CD27","PDCD1"), c("CTLA4","TNFRSF4"),
  c("CD40","CD28"), c("CD86","TNFRSF4"), c("CD28","CD86"),
  c("CD80","TNFSF9"), c("CD274","VSIR"), c("CD86","HAVCR2"),
  c("CD40","PDCD1"), c("CD86","CD200"), c("CD40","CD80"),
  c("CD28","CD276"), c("CD40","CD274"), c("TNFRSF14","CD86"))
genes_needed <- unique(unlist(pairs))
missing <- setdiff(genes_needed, rownames(coh$expr))
if (length(missing) > 0) message("genes not found, excluded: ", paste(missing, collapse = ", "))
present <- intersect(genes_needed, rownames(coh$expr))
qn <- apply(coh$expr[present, , drop = FALSE], 1, function(x) {
  r <- rank(x, ties.method = "average"); qnorm((r - 0.5) / length(x)) })
score <- vapply(pairs, function(p) {
  if (!all(p %in% colnames(qn))) return(rep(NA_integer_, nrow(qn)))
  as.integer(qn[, p[1]] > qn[, p[2]]) }, integer(nrow(qn)))
res <- data.frame(patient = coh$patients,
                  cluster = as.character(cl$cluster[coh$patients]),
                  IMPRES = rowSums(score, na.rm = TRUE))
res$Kyn <- scores$Kynurenine[match(res$patient, scores$patient)]
cat("=== IMPRES by subtype ===\n")
print(res %>% group_by(cluster) %>%
  summarise(mean = round(mean(IMPRES), 2), sd = round(sd(IMPRES), 2), n = n(), .groups = "drop"))
kw <- kruskal.test(IMPRES ~ cluster, data = res)
cat("Kruskal-Wallis p =", format.pval(kw$p.value), "\n")
rho <- cor.test(res$Kyn, res$IMPRES, method = "spearman")
cat("Kyn vs IMPRES: rho =", round(unname(rho$estimate), 3), " p =", format.pval(rho$p.value), "\n")
res$IMPRES_resp <- res$IMPRES >= 10
print(res %>% group_by(cluster) %>%
  summarise(resp_rate = paste0(sum(IMPRES_resp), "/", n(), " (", round(100 * mean(IMPRES_resp), 1), "%)"), .groups = "drop"))
write.csv(res, file.path(dirs$results, "16_impress.csv"), row.names = FALSE)
pal <- c(C1 = "#F8766D", C2 = "#7CAE00", C3 = "#00BFC4")
pA <- ggplot(res, aes(cluster, IMPRES, fill = cluster)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(aes(color = cluster), width = 0.15, size = 1.4, alpha = 0.6) +
  scale_fill_manual(values = pal) + scale_color_manual(values = pal) +
  theme_bw(base_size = 12) +
  labs(x = NULL, y = "IMPRES score (0-15)",
       title = paste0("IMPRES by subtype (KW p = ", signif(kw$p.value, 3), ")")) +
  theme(legend.position = "none")
pB <- ggplot(res, aes(Kyn, IMPRES, color = cluster)) +
  geom_point(size = 1.6, alpha = 0.75) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = pal) + theme_bw(base_size = 12) +
  labs(x = "Kynurenine score", y = "IMPRES score",
       title = paste0("rho = ", round(unname(rho$estimate), 3), ", p = ", signif(rho$p.value, 3)))
ggsave(file.path(dirs$figures, "FigureS5_IMPRES.pdf"),
       gridExtra::arrangeGrob(pA, pB, ncol = 2), width = 11, height = 5)
message("16 done")
