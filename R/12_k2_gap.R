# 12_k2_gap.R — k=2 敏感性 + gap statistic
source("R/00_config.R")
suppressPackageStartupMessages({ library(tidyverse); library(cluster) })
scr <- readRDS(file.path(dirs$results, "02_scores.rds"))
cl  <- readRDS(file.path(dirs$results, "03_cluster.rds"))
scores <- scr$scores
cc <- cl$cc

cl2 <- cc[[2]]$consensusClass
prof <- data.frame(cluster = paste0("k2_", cl2),
                   Kyn = scores$Kynurenine,
                   T_cell = scores$T_cell,
                   Cytotoxicity = scores$Cytotoxicity)
cat("=== k=2 簇大小 ===
"); print(table(cl2))
cat("=== k=2 各簇评分均值 ===
")
print(prof %>% group_by(cluster) %>%
  summarise(across(everything(), function(x) mean(x, na.rm = TRUE)), n = n()) %>%
  mutate(across(where(is.numeric), function(x) round(x, 3))))

pm <- as.matrix(t(scores %>% select(-patient, -IDO1, -TDO2, -Down_score)))
set.seed(42)
g <- clusGap(pm, FUN = kmeans, K.max = 6, B = 100)
cat("=== Gap statistic ===
")
print(round(g$Tab[, c("logW", "E.logW", "gap", "SE.sim")], 3))
write.csv(g$Tab, file.path(dirs$results, "12_gap_statistic.csv"), row.names = FALSE)
message("12 done")
