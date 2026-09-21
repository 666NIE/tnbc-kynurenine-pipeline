# 20_metabric_m1_survival.R — M1 vs non-M1 OS + Cox PH assumption tests
# v2：修正年龄提取（as.numeric(as.character(AGE_AT_DIAGNOSIS))，原 grep 方案抓到 0-3 编码列）
# 并加入全分层敏感性模型（stage x 年龄三分位）
source("R/00_config.R")
suppressPackageStartupMessages({ library(tidyverse); library(survival) })

assign <- read.csv("results/13_metabric_assignments.csv", stringsAsFactors = FALSE)
x      <- readRDS("results/06_metabric.rds")
msurv  <- x$msurv
msurv$M1 <- factor(ifelse(msurv$PATIENT_ID_EXPR %in% assign$patient[assign$cluster == "M1"],
                          "M1", "nonM1"), levels = c("nonM1", "M1"))

msurv$age   <- as.numeric(as.character(msurv$AGE_AT_DIAGNOSIS))
msurv$stage <- as.numeric(factor(msurv$TUMOR_STAGE))

cat("=== age check ===\n"); print(summary(msurv$age))
cat("=== M1 vs non-M1 ===\n"); print(table(msurv$M1))
sd <- survdiff(Surv(os_time, os_event) ~ M1, data = msurv)
cat(sprintf("log-rank p = %.4g\n", 1 - pchisq(sd$chisq, length(sd$n) - 1)))

f_m1  <- coxph(Surv(os_time, os_event) ~ M1, data = msurv)
f_adj <- coxph(Surv(os_time, os_event) ~ M1 + scale(age) + stage, data = msurv)
for (nm in c("f_m1", "f_adj")) {
  s <- summary(get(nm))
  cat(sprintf("%s (n=%d): M1 HR = %.2f (%.2f-%.2f), p = %.4g\n", nm, s$n,
      s$conf.int[1,"exp(coef)"], s$conf.int[1,"lower .95"], s$conf.int[1,"upper .95"],
      s$coefficients[1,"Pr(>|z|)"]))
}

# 全分层敏感性模型（协变量 PH 不满足时的最稳健设定）
q <- quantile(msurv$age, c(0, 1/3, 2/3, 1), na.rm = TRUE)
msurv$age_grp3 <- cut(msurv$age, breaks = unique(q), include.lowest = TRUE)
f_strat <- coxph(Surv(os_time, os_event) ~ M1 + strata(TUMOR_STAGE, age_grp3), data = msurv)
s5 <- summary(f_strat)
cat(sprintf("f_strat (n=%d): M1 HR = %.2f (%.2f-%.2f), p = %.4g\n", s5$n,
    s5$conf.int[1,"exp(coef)"], s5$conf.int[1,"lower .95"], s5$conf.int[1,"upper .95"],
    s5$coefficients[1,"Pr(>|z|)"]))

cat("\n=== cox.zph: PH assumption ===\n")
print(cox.zph(f_m1)); print(cox.zph(f_adj)); print(cox.zph(f_strat))

saveRDS(list(f_m1 = f_m1, f_adj = f_adj, f_strat = f_strat),
        file.path(dirs$results, "20_metabric_m1_survival.rds"))
message("20 done")
