# 26_expression_matched_null.R -- n=1000 expression-matched random sets: empirical null for KP C1 z-score
source("R/00_config.R")
suppressPackageStartupMessages(library(tidyverse))
coh <- readRDS(file.path(dirs$results, "01_cohort.rds"))
scr <- readRDS(file.path(dirs$results, "02_scores.rds"))
cl  <- readRDS(file.path(dirs$results, "03_cluster.rds"))
scores <- scr$scores; cluster <- cl$cluster
kp <- params$gene_sets$Kynurenine
kp_mean <- rowMeans(coh$expr[kp, , drop = FALSE])
allg <- rownames(coh$expr); gmean <- rowMeans(coh$expr)
qs <- quantile(gmean, seq(0, 1, 0.1))
kp_bin <- cut(kp_mean, qs, include.lowest = TRUE)
bin_pool <- split(allg, cut(gmean, qs, include.lowest = TRUE))
set.seed(99); n_iter <- 1000; z_c1 <- numeric(n_iter)
for (i in seq_len(n_iter)) {
  gs <- vapply(seq_along(kp), function(j) sample(bin_pool[[as.character(kp_bin[j])]], 1), character(1))
  rss <- as.numeric(scale(as.numeric(simple_ssgsea(coh$expr, gs))))
  z_c1[i] <- mean(rss[cluster == "C1"])
  if (i %% 100 == 0) cat(i, "done\n")
}
kp_real <- as.numeric(scale(as.numeric(simple_ssgsea(coh$expr, kp))))
z_real <- mean(kp_real[cluster == "C1"])
cat(sprintf("KP set C1 z = %.3f | null mean %.3f (SD %.3f) | empirical two-sided p = %.4f\n",
            z_real, mean(z_c1), sd(z_c1), mean(abs(z_c1) >= abs(z_real))))
write.csv(data.frame(z_C1 = z_c1), file.path(dirs$results, "26_null_z_distribution.csv"), row.names = FALSE)
message("26 done")
