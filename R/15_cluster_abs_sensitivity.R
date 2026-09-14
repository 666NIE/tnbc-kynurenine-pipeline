# 15_cluster_abs_sensitivity.R -- swap-in clustering: absolute-expression z replaces Kyn ssGSEA
source("R/00_config.R")
suppressPackageStartupMessages({ library(tidyverse); library(ConsensusClusterPlus) })
set.seed(params$seed)
coh <- readRDS(file.path(dirs$results, "01_cohort.rds"))
scr <- readRDS(file.path(dirs$results, "02_scores.rds"))
cl  <- readRDS(file.path(dirs$results, "03_cluster.rds"))
scores <- scr$scores
kyn_abs <- calc_z(params$gene_sets$Kynurenine, coh$expr)
scores_swap <- scores %>% mutate(Kynurenine = kyn_abs)
pm <- as.matrix(t(scores_swap %>% select(-patient, -IDO1, -TDO2, -Down_score)))
colnames(pm) <- coh$patients
message("swap-in input: ", nrow(pm), " pathways x ", ncol(pm), " samples (Kyn = absolute z)")
cc2 <- ConsensusClusterPlus(pm, maxK = params$maxK, reps = params$reps, pItem = params$pItem,
  clusterAlg = "km", distance = "euclidean",
  plot = "pdf", title = "ConsensusClustering_absKyn", seed = params$seed)
k <- params$k_final
cl2 <- cc2[[k]]$consensusClass
names(cl2) <- colnames(pm)
prof <- data.frame(cluster = paste0("abs_", cl2), Kyn_abs = kyn_abs,
                   T_cell = scores$T_cell, Cytotoxicity = scores$Cytotoxicity)
cat("=== swap-in k =", k, " cluster profiles ===\n")
print(prof %>% group_by(cluster) %>%
  summarise(across(everything(), function(x) mean(x, na.rm = TRUE)), n = n()) %>%
  mutate(across(where(is.numeric), function(x) round(x, 3))))
orig <- as.character(cl$cluster[colnames(pm)])
top_cl <- names(which.max(tapply(kyn_abs, cl2, mean)))
cat("\nabsolute-Kyn top cluster =", top_cl, " (n =", sum(cl2 == top_cl), ")\n")
print(table(origC1 = orig == "C1", absTop = cl2 == top_cl))
pct <- round(100 * sum(orig == "C1" & cl2 == top_cl) / sum(orig == "C1"), 1)
cat("share of original C1 in absolute-top cluster:", pct, "%\n")
write.csv(data.frame(patient = colnames(pm), abs_cluster = as.character(cl2), orig_cluster = orig),
          file.path(dirs$results, "15_cluster_abs_sensitivity.csv"), row.names = FALSE)
message("15 done")
