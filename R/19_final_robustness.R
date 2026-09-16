# 19_final_robustness.R -- 第10轮：14基因全表(q值) + IDO1免疫校正Cox
source("R/00_config.R")
suppressPackageStartupMessages({library(tidyverse); library(survival)})
coh <- readRDS(file.path(dirs$results, "01_cohort.rds"))
cl  <- readRDS(file.path(dirs$results, "03_cluster.rds"))
cluster <- cl$cluster
i1 <- which(as.character(cluster) == "C1"); i3 <- which(as.character(cluster) == "C3")
res <- t(sapply(params$gene_sets$Kynurenine, function(g) {
  if (!g %in% rownames(coh$expr)) return(c(mean1 = NA, mean3 = NA, log2FC = NA, p = NA))
  v1 <- coh$expr[g, i1]; v3 <- coh$expr[g, i3]
  tt <- t.test(v1, v3)
  c(mean1 = mean(v1), mean3 = mean(v3), log2FC = mean(v1) - mean(v3), p = tt$p.value)
}))
res <- as.data.frame(res); res$q <- p.adjust(res$p, "BH")
cat("=== 19a: 14基因 C1 vs C3 ===\n"); print(round(res[order(res$log2FC, decreasing = TRUE), ], 4))
write.csv(res, file.path(dirs$results, "19_kp_gene_stats.csv"), row.names = TRUE)
x <- readRDS(file.path(dirs$results, "06_metabric.rds"))
msurv <- x$msurv
age_c <- grep("AGE", colnames(msurv), value = TRUE)[1]
stg_c <- grep("STAGE", colnames(msurv), value = TRUE)[1]
msurv$age <- as.numeric(msurv[[age_c]])
msurv$stage <- as.numeric(factor(msurv[[stg_c]]))
f0  <- coxph(Surv(os_time, os_event) ~ scale(IDO1), data = msurv)
f1  <- coxph(Surv(os_time, os_event) ~ scale(IDO1) + scale(age) + stage, data = msurv)
f2a <- coxph(Surv(os_time, os_event) ~ scale(IDO1) + scale(age) + stage + scale(CD8_score), data = msurv)
f2b <- coxph(Surv(os_time, os_event) ~ scale(IDO1) + scale(age) + stage + scale(Cytotoxic), data = msurv)
cat("\n=== 19b: METABRIC IDO1 Cox ===\n")
for (nm in c("f0","f1","f2a","f2b")) {
  s <- summary(get(nm))
  cat(sprintf("%s (n=%d): IDO1 HR/SD = %.3f (%.3f-%.3f), p = %.4g\n", nm, s$n,
      s$conf.int[1,"exp(coef)"], s$conf.int[1,"lower .95"], s$conf.int[1,"upper .95"],
      s$coefficients[1,"Pr(>|z|)"]))
}
message("19 done")
