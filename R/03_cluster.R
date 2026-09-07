# 03_cluster.R — consensus clustering + Figure 1
source("R/00_config.R")
suppressPackageStartupMessages({
  library(tidyverse); library(ConsensusClusterPlus)
  library(patchwork); library(cluster)
})
set.seed(params$seed)
coh <- readRDS(file.path(dirs$results, "01_cohort.rds"))
scr <- readRDS(file.path(dirs$results, "02_scores.rds"))
scores <- scr$scores

pm <- as.matrix(t(scores %>% select(-patient, -IDO1, -TDO2, -Down_score)))
colnames(pm) <- coh$patients
message("clustering input: ", nrow(pm), " pathways x ", ncol(pm), " samples")

cc <- ConsensusClusterPlus(
  pm, maxK = params$maxK, reps = params$reps, pItem = params$pItem,
  clusterAlg = "km", distance = "euclidean",
  plot = "pdf", title = "ConsensusClustering", seed = params$seed)

stab <- t(sapply(2:params$maxK, function(k) {
  m <- cc[[k]]$consensusMatrix
  off <- m[upper.tri(m)]
  c(PAC = mean(off > 0.1 & off < 0.9),
    Silhouette = mean(silhouette(cc[[k]]$consensusClass, as.dist(1 - m))[, 3]))
}))
rownames(stab) <- paste0("k", 2:params$maxK)
message("stability:"); print(round(stab, 3))

k <- params$k_final
message("final k = ", k)
cl <- cc[[k]]$consensusClass
kyn_mean <- tapply(scores$Kynurenine, cl, mean)
ord <- names(sort(kyn_mean, decreasing = TRUE))
lab <- setNames(paste0("C", seq_along(ord)), ord)
cluster <- factor(lab[as.character(cl)], levels = paste0("C", seq_along(ord)))
names(cluster) <- coh$patients
message("Kynurenine mean by subtype (must decrease):"); print(round(tapply(scores$Kynurenine, cluster, mean), 3))
message("subtype sizes:"); print(table(cluster))

cl_ord <- order(cluster)
cm <- cc[[k]]$consensusMatrix[cl_ord, cl_ord]
hm <- expand.grid(x = seq_len(nrow(cm)), y = seq_len(ncol(cm))); hm$v <- as.vector(cm)
p1a <- ggplot(hm, aes(x, y, fill = v)) + geom_raster() +
  scale_fill_gradient(low = "white", high = "#2166AC", limits = c(0, 1)) +
  theme_void() + labs(title = paste0("Consensus matrix (k=", k, ")"), fill = "consensus") +
  theme(plot.title = element_text(hjust = 0.5))
pca <- prcomp(t(pm), scale. = FALSE)
pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], cluster = cluster[rownames(pca$x)])
ve <- round(100 * summary(pca)$importance[2, 1:2], 1)
p1b <- ggplot(pca_df, aes(PC1, PC2, color = cluster)) +
  geom_point(size = 3, alpha = .85) + stat_ellipse(level = .8, linetype = 2, show.legend = FALSE) +
  scale_color_manual(values = c(C1 = "#F8766D", C2 = "#619CFF", C3 = "#00BA38")) + theme_paper() +
  labs(title = "PCA of metabolic-immune pathways",
       x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)"), color = "Subtype")
p1c <- ggplot(data.frame(cluster), aes(cluster, fill = cluster)) +
  geom_bar(width = .65, color = "black") +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  scale_fill_manual(values = c(C1 = "#F8766D", C2 = "#619CFF", C3 = "#00BA38")) +
  theme_paper() + theme(legend.position = "none") +
  labs(title = "Sample distribution", x = "Subtype", y = "Count")
fig1 <- (p1a | p1b) / p1c + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))
ggsave(file.path(dirs$figures, "Figure1_subtyping.pdf"), fig1, width = 12, height = 9)

write.csv(data.frame(patient = names(cluster), cluster = as.character(cluster)),
          file.path(dirs$results, "03_cluster_assignments.csv"), row.names = FALSE)
saveRDS(list(cluster = cluster, stab = stab, k = k, cc = cc),
        file.path(dirs$results, "03_cluster.rds"))
message("saved 03_cluster.rds")
