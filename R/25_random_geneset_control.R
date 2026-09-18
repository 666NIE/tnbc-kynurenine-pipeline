# 25_random_geneset_control.R -- random 14-gene set negative control for relative-score background sensitivity
source("R/00_config.R")
suppressPackageStartupMessages(library(tidyverse))
coh <- readRDS(file.path(dirs$results, "01_cohort.rds"))
scr <- readRDS(file.path(dirs$results, "02_scores.rds"))
cl  <- readRDS(file.path(dirs$results, "03_cluster.rds"))
scores <- scr$scores; cluster <- cl$cluster
allg <- rownames(coh$expr)
set.seed(2026)
for (i in 1:3) {
  gs <- sample(allg, 14)
  rss <- as.numeric(simple_ssgsea(coh$expr, gs))
  df <- scores %>% mutate(cluster = cluster, rss = as.numeric(scale(rss)))
  print(df %>% group_by(cluster) %>% summarise(rand_ssGSEA = round(mean(rss), 2), n = n()))
  kw <- kruskal.test(rss ~ cluster, data = df)
  rk <- suppressWarnings(cor.test(df$rss, df$Kynurenine, method = "spearman"))
  cat(sprintf("[random %d] KW across subtypes p = %s | corr with real Kyn score: rho = %.3f\n",
              i, format.pval(kw$p.value), unname(rk$estimate)))
}
message("25 done")
