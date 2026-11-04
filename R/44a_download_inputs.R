# 44a_download_inputs.R
# Recreate the three inputs consumed by R/44_pancancer_survey.R:
#   results/44_input_pancancer_expr.rds   (genes x samples log-expression matrix; rownames=gene symbol)
#   results/44_input_pancancer_clin.csv   (columns: sample, cancer_type, OS_MONTHS, OS_EVENT, age, sex)
#   results/44_input_depmap_lineage.csv   (columns: model_id, lineage, KP_z, IDO1_z)
# Run from repository root:  Rscript R/44a_download_inputs.R
# NOTE: network/Xena dataset names can change. If use_UCSCXenaTools <- TRUE and it stops,
# inspect printed candidates, then either add XenaDownload calls or download CSVs manually.

source("R/00_config.R")
suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(stringr)
})

kp_genes <- c("IDO1","IDO2","TDO2","IL4I1","KYNU","KMO","HAAO",
              "AFMID","AADAT","GPT2","QPRT","CCBL2","NADSYN1","KYAT1")  # must match 44_pancancer_survey.R
cyto_genes <- c("CD8A","CD8B","GZMB","PRF1")                         # must match 44_pancancer_survey.R
need_genes <- c(kp_genes, cyto_genes)

results_dir <- "results"
data_raw <- "data-raw"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(data_raw, "xena"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(data_raw, "depmap"), showWarnings = FALSE, recursive = TRUE)
writeLines(capture.output(sessionInfo()), file.path(results_dir, "44a_sessionInfo.txt"))

first_existing <- function(nms, candidates) {
  hit <- nms[tolower(nms) %in% tolower(candidates)]
  if (length(hit)) hit[1] else NA_character_
}

# ---------- Helper: Xena matrix csv -> genes x samples matrix ----------
xena_csv_to_gene_by_sample <- function(dt, need_genes) {
  sample_col <- names(dt)[1]
  genes <- names(dt)[-1]
  hit <- intersect(need_genes, genes)
  missing <- setdiff(need_genes, hit)
  message("Matched genes: ", length(hit), "/", length(need_genes),
          ifelse(length(missing), paste0("; missing: ", paste(missing, collapse=",")), ""))
  if (length(intersect(kp_genes, hit)) < length(kp_genes)) stop("Too few KP genes matched. Check gene-name convention (symbol vs Ensembl).")
  m <- as.matrix(dt[, c(sample_col, hit), with = FALSE][, -1, with = FALSE])
  storage.mode(m) <- "numeric"
  rownames(m) <- dt[[sample_col]]
  t(m)  # genes x samples
}

# ========================= 1) TCGA Pan-Cancer from UCSC Xena =========================
use_UCSCXenaTools <- FALSE
xena_expr_file <- file.path(data_raw, "xena", "tcga_pancan_htseq_fpkm.csv")              # samples x genes
xena_surv_file <- file.path(data_raw, "xena", "tcga_pancan_survival_phenotype.csv")      # samples x clinical columns

if (use_UCSCXenaTools) {
  if (!requireNamespace("UCSCXenaTools", quietly = TRUE)) stop("Install UCSCXenaTools or set use_UCSCXenaTools <- FALSE and download CSVs manually.")
  library(UCSCXenaTools)
  # UCSCXenaTools::XenaDataUpdate()  # uncomment once if your local dataset table is stale
  xd <- UCSCXenaTools::XenaData
  print(xd[grep("TCGA.*Pan-Cancer.*HTSeq.*FPKM", paste(xd$XenaDatasets, xd$XenaCohorts), ignore.case = TRUE), ])
  print(xd[grep("TCGA.*Pan-Cancer.*(survival|phenotype)", paste(xd$XenaDatasets, xd$XenaCohorts), ignore.case = TRUE), ])
  stop("Inspect printed candidates, then add UCSCXenaTools::XenaDownload(...) or download via browser and set file paths above.")
}

stopifnot(file.exists(xena_expr_file), file.exists(xena_surv_file))
expr_dt <- fread(xena_expr_file, check.names = FALSE)
pan <- xena_csv_to_gene_by_sample(expr_dt, need_genes)

clin_raw <- fread(xena_surv_file, check.names = FALSE) %>% as_tibble()
colmap <- c(
  sample = first_existing(names(clin_raw), c("sample","_sample","sample_id","_sample_id")),
  cancer_type = first_existing(names(clin_raw), c("cancer_type","cancer type abbreviation","_cancer_type","oncotree_code")),
  OS_MONTHS = first_existing(names(clin_raw), c("os_months","_os_months","overall_survival_months")),
  OS_EVENT = first_existing(names(clin_raw), c("os_event","_os_ind","_os_event","overall_survival_event")),
  age = first_existing(names(clin_raw), c("age","_age","age_at_initial_pathologic_diagnosis")),
  sex = first_existing(names(clin_raw), c("sex","_sex","gender"))
)
if (any(is.na(colmap))) stop("Missing clinical columns after auto-mapping. Found: ", paste(names(clin_raw), collapse=", "))
clin <- clin_raw %>%
  rename(sample = !!sym(colmap[["sample"]]),
         cancer_type = !!sym(colmap[["cancer_type"]]),
         OS_MONTHS = !!sym(colmap[["OS_MONTHS"]]),
         OS_EVENT = !!sym(colmap[["OS_EVENT"]]),
         age = !!sym(colmap[["age"]]),
         sex = !!sym(colmap[["sex"]])) %>%
  transmute(sample = as.character(sample),
            cancer_type = as.character(cancer_type),
            OS_MONTHS = as.numeric(OS_MONTHS),
            OS_EVENT = as.numeric(OS_EVENT),
            age = as.numeric(age),
            sex = as.character(sex))

ov <- length(intersect(colnames(pan), clin$sample))
message("Pan-cancer sample overlap expr vs clin: ", ov)
if (ov < 10000) stop("Sample overlap too small. Check barcode conventions between expression and clinical files.")

saveRDS(pan, file.path(results_dir, "44_input_pancancer_expr.rds"))
write.csv(clin, file.path(results_dir, "44_input_pancancer_clin.csv"), row.names = FALSE)

# ========================= 2) DepMap Public 24Q4 lineage summary =========================
# Put DepMap Public 24Q4 download files here:
#   Model.csv
#   OmicsExpressionProteinCodingGenesTPMLogp1.csv   (or current protein-coding expression TPM file)
dep_model_file <- file.path(data_raw, "depmap", "Model.csv")
dep_expr_file  <- file.path(data_raw, "depmap", "OmicsExpressionProteinCodingGenesTPMLogp1.csv")
stopifnot(file.exists(dep_model_file), file.exists(dep_expr_file))

model <- fread(dep_model_file, check.names = FALSE) %>% as_tibble()
lineage_col <- first_existing(names(model), c("OncotreeLineage","lineage","PrimaryTumorOrigin","OncotreePrimaryDisease"))
if (is.na(lineage_col)) stop("Could not find a lineage column in Model.csv")
model <- model %>% transmute(model_id = ModelID, lineage = .data[[lineage_col]])

header <- names(fread(dep_expr_file, nrows = 0))
gene_pat <- paste0("\\((?:", paste(need_genes, collapse = "|"), ")\\)")
keep_cols <- header[str_detect(header, regex(gene_pat, ignore_case = TRUE))]
if (length(keep_cols) < length(need_genes)) {
  message("Only matched expression columns: ", paste(keep_cols, collapse = "; "))
  stop("Too few gene columns matched. DepMap expression headers must contain symbols like 'ENSG000001...(IDO1)'.")
}
sel <- c(1, which(header %in% keep_cols))
dep_dt <- fread(dep_expr_file, select = sel, check.names = FALSE)
mat <- as.matrix(dep_dt[, -1, with = FALSE]); storage.mode(mat) <- "numeric"
rownames(mat) <- dep_dt[[1]]
genesym <- str_match(keep_cols, "\\(([^)]+)\\)")[, 2]
colnames(mat) <- genesym
gdep <- t(mat)  # genes x cell lines

kp_present <- intersect(kp_genes, rownames(gdep))
if (length(kp_present) < length(kp_genes)) stop("Missing KP genes in DepMap expression: ", paste(setdiff(kp_genes, kp_present), collapse=", "))
KP_score <- colMeans(scale(t(gdep[kp_present, , drop = FALSE])), na.rm = TRUE)
IDO1_z <- as.numeric(scale(as.numeric(gdep["IDO1", ])))

dep <- tibble(model_id = colnames(gdep), KP_z = as.numeric(KP_score[colnames(gdep)]), IDO1_z = IDO1_z) %>%
  inner_join(model, by = "model_id") %>%
  select(model_id, lineage, KP_z, IDO1_z) %>%
  filter(!is.na(lineage), !is.na(KP_z), !is.na(IDO1_z))
write.csv(dep, file.path(results_dir, "44_input_depmap_lineage.csv"), row.names = FALSE)

message("44a done: wrote results/44_input_pancancer_expr.rds, results/44_input_pancancer_clin.csv, results/44_input_depmap_lineage.csv")
message("Next: Rscript R/44_pancancer_survey.R should regenerate results/42b/42c/42d.")
