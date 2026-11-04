# 01_data.R
source("R/00_config.R")
suppressPackageStartupMessages(library(tidyverse))
stopifnot(file.exists(path.expand(paths$tcga_expr)),
          file.exists(path.expand(paths$tcga_surv)))

expr_raw <- readRDS(path.expand(paths$tcga_expr))
med <- rowMeans(expr_raw, na.rm = TRUE)
keep <- tibble(g = rownames(expr_raw), m = med) %>%
  arrange(desc(m)) %>% distinct(g) %>% pull(g)
expr <- expr_raw[keep, , drop = FALSE]
message("expression: ", nrow(expr), " genes x ", ncol(expr), " samples")

sv <- read.delim(path.expand(paths$tcga_surv), check.names = FALSE)
fr <- sapply(sv, function(v)
  suppressWarnings(mean(grepl("^TCGA-", as.character(v)), na.rm = TRUE)))
surv <- data.frame(
  patient  = patient_id(as.character(sv[[names(which.max(fr))]])),
  os_event = as.integer(sv[["OS"]]),
  os_time  = as.numeric(sv[["OS.time"]]),
  stringsAsFactors = FALSE) %>%
  distinct(patient, .keep_all = TRUE)

pid_expr <- patient_id(colnames(expr))
if (any(duplicated(pid_expr))) {
  expr <- expr[, !duplicated(pid_expr), drop = FALSE]
  pid_expr <- patient_id(colnames(expr))
}
patients <- intersect(pid_expr, surv$patient)
message("analysis cohort: ", length(patients), " patients")
expr <- expr[, match(patients, pid_expr), drop = FALSE]
surv <- surv[match(patients, surv$patient), ]

saveRDS(list(expr = expr, surv = surv, patients = patients),
        file.path(dirs$results, "01_cohort.rds"))
message("saved 01_cohort.rds")
