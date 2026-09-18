# 24_cibersortx_check.R -- CIBERSORTx LM22 deconvolution: independent immune correction (pre-submission defense)
source("R/00_config.R")
suppressPackageStartupMessages(library(tidyverse))
coh <- readRDS(file.path(dirs$results, "01_cohort.rds"))
scr <- readRDS(file.path(dirs$results, "02_scores.rds"))
cl  <- readRDS(file.path(dirs$results, "03_cluster.rds"))
scores <- scr$scores; cluster <- cl$cluster
out <- file.path(path.expand("~/Downloads"), "CIBERSORTx_input.txt")
df <- data.frame(GeneSymbol = rownames(coh$expr), coh$expr, check.names = FALSE)
write.table(df, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Exported:", out, "(", nrow(df), "genes x", ncol(df) - 1, "samples)\n")
res_f <- file.path(path.expand("~/Downloads"), "CIBERSORTx_Results.txt")
if (!file.exists(res_f)) {
  stop("Run CIBERSORTx web (LM22, disable Quantile normalization), download CIBERSORTx_Results.txt to ~/Downloads/, then re-source.")
}
cs <- read.delim(res_f, check.names = FALSE)
cd8  <- cs[["T cells CD8"]]
lymph_cols <- grep("^B cells|^Plasma cells|^T cells|^NK cells", colnames(cs), value = TRUE)
lymph <- rowSums(cs[, lymph_cols, drop = FALSE])
stopifnot(nrow(scores) == nrow(cs))
df2 <- scores %>% mutate(cluster = cluster, CD8frac = cd8, LymphFrac = lymph)
cat("\n=== CD8/Lymph fraction by subtype ===\n")
print(df2 %>% group_by(cluster) %>% summarise(CD8 = round(mean(CD8frac), 3),
      Lymph = round(mean(LymphFrac), 3), n = n()))
cat(sprintf("CD8 fraction KW p = %s\n",
    format.pval(kruskal.test(CD8frac ~ cluster, data = df2)$p.value)))
pcor <- function(x, y, covs) {
  rx <- resid(lm(x ~ ., data = as.data.frame(covs)))
  ry <- resid(lm(y ~ ., data = as.data.frame(covs)))
  suppressWarnings(cor.test(rx, ry, method = "spearman"))
}
r0 <- pcor(df2$Kynurenine, df2$Cytotoxicity, data.frame(T = df2$T_cell))
r1 <- pcor(df2$Kynurenine, df2$Cytotoxicity, data.frame(T = df2$T_cell, CD8 = df2$CD8frac))
r2 <- pcor(df2$Kynurenine, df2$Cytotoxicity, data.frame(T = df2$T_cell, Lymph = df2$LymphFrac))
cat(sprintf("adj T_cell:             rho = %.3f, p = %s\n", unname(r0$estimate), format.pval(r0$p.value)))
cat(sprintf("adj T_cell + CD8frac:   rho = %.3f, p = %s\n", unname(r1$estimate), format.pval(r1$p.value)))
cat(sprintf("adj T_cell + LymphFrac: rho = %.3f, p = %s\n", unname(r2$estimate), format.pval(r2$p.value)))
df2$resid_kyn_cd8 <- resid(lm(Kynurenine ~ CD8frac, data = df2))
df2$resid_kyn_lym <- resid(lm(Kynurenine ~ LymphFrac, data = df2))
print(df2 %>% group_by(cluster) %>%
  summarise(resid_Kyn_CD8 = round(mean(resid_kyn_cd8), 2),
            resid_Kyn_Lymph = round(mean(resid_kyn_lym), 2), n = n()))
cat(sprintf("Kyn score KW after CD8 correction:   p = %s\n",
    format.pval(kruskal.test(resid_kyn_cd8 ~ cluster, data = df2)$p.value)))
cat(sprintf("Kyn score KW after Lymph correction: p = %s\n",
    format.pval(kruskal.test(resid_kyn_lym ~ cluster, data = df2)$p.value)))
write.csv(df2, file.path(dirs$results, "24_cibersortx_correction.csv"), row.names = FALSE)
message("24 done")
