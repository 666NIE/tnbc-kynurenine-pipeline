# 14_sensitivity.R — 评分体系敏感性分析
source("R/00_config.R")
suppressPackageStartupMessages(library(tidyverse))

coh <- readRDS(file.path(dirs$results, "01_cohort.rds"))
scr <- readRDS(file.path(dirs$results, "02_scores.rds"))
cl  <- readRDS(file.path(dirs$results, "03_cluster.rds"))
expr <- coh$expr; scores <- scr$scores; cluster <- cl$cluster

kyn_abs <- calc_z(params$gene_sets$Kynurenine, expr)
df <- data.frame(cluster = cluster, abs = kyn_abs,
                 ssgsea = scores$Kynurenine,
                 cyto = scores$Cytotoxicity, tcell = scores$T_cell)

cat("A1 绝对版Kyn by 亚型, KW p =", fmt_p(kruskal.test(abs ~ cluster, df)$p.value), "
")
print(round(tapply(df$abs, df$cluster, mean, na.rm = TRUE), 3))
for (v in c("abs", "ssgsea")) {
  ct <- suppressWarnings(cor.test(df[[v]], df$cyto, method = "spearman"))
  cat(sprintf("A  %s版 vs 细胞毒性: rho = %+.3f, %s
", v, unname(ct$estimate), fmt_p(ct$p.value)))
}
part <- function(x, y, z) {
  suppressWarnings(cor.test(resid(lm(x ~ z)), resid(lm(y ~ z)), method = "spearman"))
}
for (v in c("ssgsea", "abs")) {
  pt <- part(df[[v]], df$cyto, df$tcell)
  cat(sprintf("B  %s版 ~ 细胞毒性 | T细胞: rho = %+.3f, %s
", v, unname(pt$estimate), fmt_p(pt$p.value)))
}
m <- lm(cyto ~ cluster + tcell, df)
cat("C  细胞毒性 ~ 亚型 + T细胞: 亚型项 p =", fmt_p(anova(m)$"Pr(>F)"[1]), "
")
write.csv(df, file.path(dirs$results, "14_sensitivity.csv"), row.names = FALSE)
message("14 done")
