# 29b_metabric_m1_survival_npi.R — M1 OS with age + NPI adjustment
# 目的：与锁定分类器（locked classifier）分析统一协变量口径（年龄 + NPI，n=319 无分期缺失）
# 原 29 号脚本（age + stage 版）保持不变，本脚本为并行口径
source("R/00_config.R")
suppressPackageStartupMessages({ library(tidyverse); library(survival) })

assign <- read.csv("results/13_metabric_assignments.csv", stringsAsFactors = FALSE)
x      <- readRDS("results/06_metabric.rds")
msurv  <- x$msurv
names(msurv)   # 确认 NPI 列名；若叫别的，改下一行
msurv$M1 <- factor(ifelse(msurv$PATIENT_ID_EXPR %in% assign$patient[assign$cluster == "M1"],
                          "M1", "nonM1"), levels = c("nonM1", "M1"))

msurv$age <- as.numeric(as.character(msurv$AGE_AT_DIAGNOSIS))
msurv$NPI <- as.numeric(as.character(msurv$NPI))

cat("=== NPI check ===\n"); print(summary(msurv$NPI))
cat("=== M1 vs non-M1 ===\n"); print(table(msurv$M1))

# 主模型：年龄 + NPI（与锁定分类器分析同一套协变量）
f_adj_npi <- coxph(Surv(os_time, os_event) ~ M1 + scale(age) + scale(NPI), data = msurv)
s <- summary(f_adj_npi)
cat(sprintf("f_adj_npi (n=%d): M1 HR = %.2f (%.2f-%.2f), p = %.4g\n", s$n,
    s$conf.int[1,"exp(coef)"], s$conf.int[1,"lower .95"], s$conf.int[1,"upper .95"],
    s$coefficients[1,"Pr(>|z|)"]))

# 分层稳健性：年龄三分位分层 + NPI（避免分期 86 例缺失，n=319 全覆盖）
q <- quantile(msurv$age, c(0, 1/3, 2/3, 1), na.rm = TRUE)
msurv$age_grp3 <- cut(msurv$age, breaks = unique(q), include.lowest = TRUE)
f_strat_npi <- coxph(Surv(os_time, os_event) ~ M1 + scale(NPI) + strata(age_grp3), data = msurv)
s2 <- summary(f_strat_npi)
cat(sprintf("f_strat_npi (n=%d): M1 HR = %.2f (%.2f-%.2f), p = %.4g\n", s2$n,
    s2$conf.int[1,"exp(coef)"], s2$conf.int[1,"lower .95"], s2$conf.int[1,"upper .95"],
    s2$coefficients[1,"Pr(>|z|)"]))

cat("\n=== cox.zph: PH assumption ===\n")
print(cox.zph(f_adj_npi)); print(cox.zph(f_strat_npi))

saveRDS(list(f_adj_npi = f_adj_npi, f_strat_npi = f_strat_npi),
        file.path(dirs$results, "29b_metabric_m1_survival_npi.rds"))
message("29b done")
