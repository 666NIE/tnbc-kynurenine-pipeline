# 13_metabric_cluster.R — METABRIC 独立分型
source("R/00_config.R")
suppressPackageStartupMessages({ library(tidyverse); library(ConsensusClusterPlus) })

read_clin <- function(file) read.delim(file, comment.char = "#",
                                        check.names = FALSE, stringsAsFactors = FALSE)
mf <- list.files(path.expand(paths$metabric), recursive = TRUE, full.names = TRUE)
me <- mf[grepl("data_mrna.*microarray", mf, ignore.case = TRUE)]
me <- me[!grepl("zscore", me, ignore.case = TRUE)][1]
ms <- mf[grepl("data_clinical_sample",  mf, ignore.case = TRUE)][1]
mp <- mf[grepl("data_clinical_patient", mf, ignore.case = TRUE)][1]
stopifnot(!any(is.na(c(me, ms, mp))))

mraw <- read.delim(me, check.names = FALSE, stringsAsFactors = FALSE)
mexpr <- mraw %>%
  mutate(median_expr = rowMeans(across(where(is.numeric)), na.rm = TRUE)) %>%
  arrange(desc(median_expr)) %>% distinct(Hugo_Symbol, .keep_all = TRUE) %>%
  select(-median_expr) %>% column_to_rownames("Hugo_Symbol")
if ("Entrez_Gene_Id" %in% colnames(mexpr)) mexpr <- mexpr %>% select(-Entrez_Gene_Id)
mexpr <- as.matrix(mexpr)

mclin <- read_clin(ms) %>% left_join(read_clin(mp), by = "PATIENT_ID")
yc <- grep(".y", colnames(mclin), value = TRUE, fixed = TRUE)
if (length(yc) > 0) {
  mclin <- mclin %>% select(-all_of(yc))
  xc <- grep(".x", colnames(mclin), value = TRUE, fixed = TRUE)
  colnames(mclin)[match(xc, colnames(mclin))] <- sub(".x", "", xc, fixed = TRUE)
}
mclin$PATIENT_ID_EXPR <- mclin$PATIENT_ID
er <- grep("ER_STATUS|ER_IHC|IHC_ER", colnames(mclin), value = TRUE, ignore.case = TRUE)[1]
pr <- grep("PR_STATUS|PR_IHC|IHC_PR", colnames(mclin), value = TRUE, ignore.case = TRUE)[1]
he <- grep("HER2_STATUS|HER2_IHC|IHC_HER2", colnames(mclin), value = TRUE, ignore.case = TRUE)[1]
mtnbc <- mclin %>% filter(.data[[er]] == "Negative", .data[[pr]] == "Negative", .data[[he]] == "Negative")
mids <- intersect(mtnbc$PATIENT_ID_EXPR, colnames(mexpr))
mexpr <- mexpr[, mids, drop = FALSE]
message("METABRIC TNBC: ", ncol(mexpr))

ss_methods <- c("Kynurenine", "Glycolysis", "Glutamine", "Adenosine")
mss <- bind_cols(
  tibble(patient = mids),
  as_tibble(sapply(names(params$gene_sets), function(nm) {
    v <- if (nm %in% ss_methods) simple_ssgsea(mexpr, params$gene_sets[[nm]])
         else calc_z(params$gene_sets[[nm]], mexpr)
    as.numeric(scale(v))
  })))
mpm <- as.matrix(t(mss %>% select(-patient)))
colnames(mpm) <- mids
message("METABRIC 聚类输入: ", nrow(mpm), " x ", ncol(mpm))

set.seed(42)
mcc <- ConsensusClusterPlus(mpm, maxK = 6, reps = 100, pItem = 0.8,
                            clusterAlg = "km", distance = "euclidean",
                            plot = "pdf", title = "ConsensusClustering_METABRIC",
                            seed = 42)
mk <- 3
mcl <- mcc[[mk]]$consensusClass
kyn_mean <- tapply(mss$Kynurenine, mcl, mean)
ord <- names(sort(kyn_mean, decreasing = TRUE))
lab <- setNames(paste0("M", seq_along(ord)), ord)
mcluster <- factor(lab[as.character(mcl)], levels = paste0("M", seq_along(ord)))
cat("=== METABRIC 亚型频数 ===
"); print(table(mcluster))
cat("=== Kyn 均值（须递减）===
"); print(round(tapply(mss$Kynurenine, mcluster, mean), 3))
cat("=== T细胞/细胞毒性均值 ===
")
print(round(sapply(c("T_cell","Cytotoxicity"), function(f) tapply(mss[[f]], mcluster, mean)), 3))
saveRDS(list(scores = mss, cluster = mcluster, cc = mcc), "results/13_metabric_cluster.rds")
write.csv(data.frame(patient = mids, cluster = as.character(mcluster)),
          "results/13_metabric_assignments.csv", row.names = FALSE)
message("13 done")
