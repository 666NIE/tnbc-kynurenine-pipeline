# 09_tcga_survival.R — exploratory TCGA OS (underpowered; reported as one sentence)
source("R/00_config.R")
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(survminer))
coh <- readRDS(file.path(dirs$results, "01_cohort.rds"))
cl  <- readRDS(file.path(dirs$results, "03_cluster.rds"))
sdf <- data.frame(patient = coh$patients, cluster = cl$cluster) %>%
  left_join(coh$surv, by = "patient") %>%
  filter(!is.na(os_time), os_time > 0)
message("TCGA survival: n = ", nrow(sdf), ", events = ", sum(sdf$os_event))
fit <- survfit(Surv(os_time, os_event) ~ cluster, data = sdf)
message("log-rank: ", fmt_p(surv_pvalue(fit)$pval))
s <- summary(coxph(Surv(os_time, os_event) ~ cluster, data = sdf))
message(sprintf("Cox C2 HR = %.2f | C3 HR = %.2f",
                s$conf.int[1, "exp(coef)"], s$conf.int[2, "exp(coef)"]))
write.csv(data.frame(n = nrow(sdf), events = sum(sdf$os_event),
                     logrank_p = surv_pvalue(fit)$pval),
          file.path(dirs$results, "09_tcga_survival.csv"), row.names = FALSE)
message("09 done")
