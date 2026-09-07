# 02_scores.R — 16-pathway ssGSEA + Kynurenine cutoff
source("R/00_config.R")
suppressPackageStartupMessages(library(tidyverse))
coh  <- readRDS(file.path(dirs$results, "01_cohort.rds"))
expr <- coh$expr
ss_methods <- c("Kynurenine", "Glycolysis", "Glutamine", "Adenosine")
ss <- bind_cols(
  tibble(patient = coh$patients),
  as_tibble(sapply(names(params$gene_sets), function(nm) {
    v <- if (nm %in% ss_methods) simple_ssgsea(expr, params$gene_sets[[nm]])
         else calc_z(params$gene_sets[[nm]], expr)
    as.numeric(scale(v))
  })))


scores <- ss
scores$IDO1 <- as.numeric(expr["IDO1", ])
scores$TDO2 <- as.numeric(expr["TDO2", ])
scores$Down_score <- calc_z(params$down_genes, expr)

ct <- suppressWarnings(cor.test(scores$Kynurenine, scores$Cytotoxicity, method = "spearman"))
message(sprintf("QC | Kyn vs Cytotoxicity: rho = %.3f, %s",
                unname(ct$estimate), fmt_p(ct$p.value)))
kyn_cutoff <- median(scores$Kynurenine, na.rm = TRUE)
message("Kynurenine cutoff = ", round(kyn_cutoff, 4))

saveRDS(list(scores = scores, kyn_cutoff = kyn_cutoff),
        file.path(dirs$results, "02_scores.rds"))
message("saved 02_scores.rds")
