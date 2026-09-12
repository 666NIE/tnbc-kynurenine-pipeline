# 06_validation.R — METABRIC external validation (Fig 10)
source("R/00_config.R")
suppressPackageStartupMessages({ library(tidyverse); library(survival); library(survminer) })
scr <- readRDS(file.path(dirs$results, "02_scores.rds"))
kyn_cutoff <- scr$kyn_cutoff
read_clin <- function(file) read.delim(file, comment.char = "#",
                                        check.names = FALSE, stringsAsFactors = FALSE)
pick_col <- function(pattern, dat) {
  hits <- grep(pattern, colnames(dat), value = TRUE, ignore.case = TRUE)
  if (length(hits) == 0) NA_character_ else hits[1]
}
meta_files <- list.files(path.expand(paths$metabric), recursive = TRUE, full.names = TRUE)
me <- meta_files[grepl("data_mrna.*microarray", meta_files, ignore.case = TRUE)]
me <- me[!grepl("zscore", me, ignore.case = TRUE)][1]
ms <- meta_files[grepl("data_clinical_sample",  meta_files, ignore.case = TRUE)][1]
mp <- meta_files[grepl("data_clinical_patient", meta_files, ignore.case = TRUE)][1]
stopifnot(!any(is.na(c(me, ms, mp))))

mraw <- read.delim(me, check.names = FALSE, stringsAsFactors = FALSE)
mexpr <- mraw %>%
  mutate(median_expr = rowMeans(across(where(is.numeric)), na.rm = TRUE)) %>%
  arrange(desc(median_expr)) %>% distinct(Hugo_Symbol, .keep_all = TRUE) %>%
  select(-median_expr) %>% column_to_rownames("Hugo_Symbol")
if ("Entrez_Gene_Id" %in% colnames(mexpr)) mexpr <- mexpr %>% select(-Entrez_Gene_Id)
mexpr <- as.matrix(mexpr)

mclin <- read_clin(ms) %>% left_join(read_clin(mp), by = "PATIENT_ID")
yc <- grep("\\.y$", colnames(mclin), value = TRUE)
if (length(yc) > 0) {
  mclin <- mclin %>% select(-all_of(yc))
  xc <- grep("\\.x$", colnames(mclin), value = TRUE)
  colnames(mclin)[match(xc, colnames(mclin))] <- sub("\\.x$", "", xc)
}
mclin$PATIENT_ID_EXPR <- mclin$PATIENT_ID
er_col  <- pick_col("ER[_\\.]STATUS|IHC[_\\.]ER", mclin)
pr_col  <- pick_col("PR[_\\.]STATUS|IHC[_\\.]PR", mclin)
her_col <- pick_col("HER2[_\\.]STATUS|IHC[_\\.]HER2", mclin)
mtnbc <- mclin %>% filter(.data[[er_col]] == "Negative",
                          .data[[pr_col]] == "Negative",
                          .data[[her_col]] == "Negative")
mids <- intersect(mtnbc$PATIENT_ID_EXPR, colnames(mexpr))
message("METABRIC TNBC: ", length(mids))
stopifnot(length(mids) >= 10)
mexpr <- mexpr[, mids, drop = FALSE]

kyn_m <- as.numeric(scale(simple_ssgsea(mexpr, params$gene_sets$Kynurenine)))
mdf <- mtnbc %>% filter(PATIENT_ID_EXPR %in% mids) %>%
  mutate(Kyn_score = kyn_m[match(PATIENT_ID_EXPR, colnames(mexpr))],
         Cytotoxic = calc_z(params$cyto_genes, mexpr)[match(PATIENT_ID_EXPR, colnames(mexpr))],
         CD8_score = calc_z(params$cd8_genes, mexpr)[match(PATIENT_ID_EXPR, colnames(mexpr))],
         IDO1 = as.numeric(mexpr["IDO1", ])[match(PATIENT_ID_EXPR, colnames(mexpr))])

osc <- pick_col("OS[_\\.]STATUS", mdf); otc <- pick_col("OS[_\\.]MONTHS", mdf)
stopifnot(!is.na(osc), !is.na(otc))
msurv <- mdf %>%
  mutate(os_event = ifelse(grepl("DECEASED|1:", .data[[osc]]), 1, 0),
         os_time  = as.numeric(.data[[otc]]),
         Kyn_group = factor(ifelse(Kyn_score > kyn_cutoff, "High", "Low"),
                            levels = c("Low", "High"))) %>%
  filter(!is.na(os_time), os_time > 0)
message("METABRIC groups (cutoff ", round(kyn_cutoff, 4), "):")
print(msurv %>% count(Kyn_group))

s1 <- summary(coxph(Surv(os_time, os_event) ~ Kyn_group, data = msurv))
message(sprintf("OS Cox Kyn(High vs Low): HR = %.2f (%.2f-%.2f), p = %.4g",
                s1$conf.int[1,"exp(coef)"], s1$conf.int[1,"lower .95"],
                s1$conf.int[1,"upper .95"], s1$coefficients[1,"Pr(>|z|)"]))
s2 <- summary(coxph(Surv(os_time, os_event) ~ scale(IDO1), data = msurv))
message(sprintf("OS Cox per SD IDO1: HR = %.3f (%.3f-%.3f), p = %.4g",
                s2$conf.int[1,"exp(coef)"], s2$conf.int[1,"lower .95"],
                s2$conf.int[1,"upper .95"], s2$coefficients[1,"Pr(>|z|)"]))
ct <- suppressWarnings(cor.test(msurv$IDO1, msurv$CD8_score, method = "spearman"))
message(sprintf("IDO1 vs CD8: rho = %.3f, %s", unname(ct$estimate), fmt_p(ct$p.value)))

msurv$IDO1_grp <- factor(ifelse(msurv$IDO1 > median(msurv$IDO1, na.rm = TRUE), "High", "Low"),
                         levels = c("Low", "High"))
fkm <- survfit(Surv(os_time, os_event) ~ IDO1_grp, data = msurv)
p10 <- ggsurvplot(fkm, data = msurv, pval = TRUE, risk.table = TRUE,
                  palette = c("#3498DB", "#E74C3C"),
                  title = "METABRIC TNBC: OS by IDO1 expression",
                  xlab = "Months", ylab = "Overall Survival",
                  legend.title = "IDO1", legend.labs = c("Low", "High"),
                  ggtheme = theme_paper())
pdf(file.path(dirs$figures, "Figure9_METABRIC_IDO1_OS.pdf"), width = 7, height = 7)
print(p10); dev.off()
write.csv(msurv, file.path(dirs$results, "06_metabric_surv.csv"), row.names = FALSE)
saveRDS(list(metabric = mdf, msurv = msurv), file.path(dirs$results, "06_metabric.rds"))
message("06 done")
