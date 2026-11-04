# 17_scrna_kp_heatmap.R -- per-cell KP gene expression across GSE176078 (reviewer Major 7)
# Run from repository root:  Rscript R/17_scrna_kp_heatmap.R
# Expected inputs are repository-relative, NOT ~/Downloads:
#   data-raw/GSE176078_scRNASeq/count_matrix_sparse.mtx
#   data-raw/GSE176078_scRNASeq/count_matrix_barcodes.tsv
#   data-raw/GSE176078_scRNASeq/count_matrix_genes.tsv
#   data-raw/GSE176078_scRNASeq/metadata.csv   (must contain celltype_major, celltype_minor, orig.ident)
# Optional override: GSE176078_DIR=/path/to/folder Rscript R/17_scrna_kp_heatmap.R

source("R/00_config.R")
suppressPackageStartupMessages({ library(tidyverse); library(Matrix); library(ggplot2); library(gridExtra) })

sc_dir <- Sys.getenv("GSE176078_DIR", unset = file.path("data-raw", "GSE176078_scRNASeq"))
mtx_path   <- file.path(sc_dir, "count_matrix_sparse.mtx")
bar_path   <- file.path(sc_dir, "count_matrix_barcodes.tsv")
genes_path <- file.path(sc_dir, "count_matrix_genes.tsv")
meta_path  <- file.path(sc_dir, "metadata.csv")
need_files <- c(mtx_path, bar_path, genes_path, meta_path)
missing_files <- need_files[!file.exists(need_files)]
if (length(missing_files)) {
  stop("Missing GSE176078 processed inputs under: ", sc_dir,
       "\nMissing:\n  ", paste(missing_files, collapse = "\n  "),
       "\nSee data-raw/GSE176078_scRNASeq/README.md for the expected provenance.")
}

barcodes <- readLines(bar_path)
message("矩阵: ", mtx_path, " (", length(barcodes), " barcodes)")

meta <- read.csv(meta_path, row.names = 1, check.names = FALSE)
need_meta_cols <- c("celltype_major", "celltype_minor", "orig.ident")
miss_meta <- setdiff(need_meta_cols, colnames(meta))
if (length(miss_meta)) stop("metadata.csv missing columns: ", paste(miss_meta, collapse=", "))
best_n <- length(intersect(rownames(meta), barcodes))
cat(sprintf("metadata: %s -> 重叠 %d\n", basename(meta_path), best_n))
if (best_n < 10000) stop("metadata 与矩阵 barcode 重叠过少: ", best_n)

kp <- c("IDO1","IDO2","TDO2","AFMID","KMO","KYNU","KYAT1","AADAT","CCBL2","GPT2","HAAO","QPRT","NADSYN1","IL4I1")
genes_df <- read.delim(genes_path, header = FALSE, stringsAsFactors = FALSE)
sym_col <- if (ncol(genes_df) >= 2) 2 else 1
mtx <- readMM(mtx_path)
colnames(mtx) <- barcodes
rownames(mtx) <- genes_df[[sym_col]]
cells <- intersect(colnames(mtx), rownames(meta))
if (length(cells) < 10000) stop("barcode 重叠过少: ", length(cells))
message("cells in both: ", length(cells))
mtx <- mtx[, cells]; meta <- meta[cells, ]
libsize <- colSums(mtx)
kp_found <- intersect(kp, rownames(mtx))
message("KP genes found: ", length(kp_found), " / ", length(kp))
cell_ids <- colnames(mtx)
cnt <- t(mtx[kp_found, , drop = FALSE])
rownames(cnt) <- cell_ids
cnt <- Diagonal(x = 10000 / libsize[cell_ids]) %*% cnt
cnt <- as(cnt, "dgCMatrix"); cnt@x <- log1p(cnt@x)
rownames(cnt) <- cell_ids
expr <- as.data.frame(as.matrix(cnt)); expr$barcode <- rownames(expr)
idx <- match(expr$barcode, rownames(meta))
stopifnot(!any(is.na(idx)))
df <- cbind(expr, meta[idx, need_meta_cols])
stopifnot(nrow(df) == length(cell_ids))
kp_order <- c("IDO1","IDO2","TDO2","AFMID","KMO","KYNU","KYAT1","AADAT","CCBL2","GPT2","HAAO","QPRT","NADSYN1","IL4I1")
genes_plot <- intersect(kp_order, kp_found)
long <- df %>%
  pivot_longer(all_of(genes_plot), names_to = "gene", values_to = "expr") %>%
  group_by(celltype_major, gene) %>%
  summarise(avg = mean(expr), pct = 100 * mean(expr > 0), n = n(), .groups = "drop")
long$gene <- factor(long$gene, levels = genes_plot)
long_z <- long %>% group_by(gene) %>% mutate(avg_z = as.numeric(scale(avg))) %>% ungroup()
pA <- ggplot(long_z, aes(celltype_major, gene)) +
  geom_point(aes(size = pct, color = avg_z)) +
  scale_color_gradient2(low = "grey85", mid = "white", high = "#B22222", midpoint = 0, name = "Scaled avg") +
  scale_size_continuous(range = c(0.4, 7), name = "% expr") +
  theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = NULL, y = NULL, title = "KP genes across major cell types (GSE176078)")
imm <- df[df$celltype_major %in% c("Myeloid","T-cells"), ]
heat <- imm %>%
  pivot_longer(all_of(genes_plot), names_to = "gene", values_to = "expr") %>%
  group_by(celltype_minor, gene) %>%
  summarise(avg = mean(expr), .groups = "drop")
heat$gene <- factor(heat$gene, levels = genes_plot)
heat_z <- heat %>% group_by(gene) %>% mutate(avg_z = as.numeric(scale(avg))) %>% ungroup()
pB <- ggplot(heat_z, aes(celltype_minor, gene, fill = avg_z)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "grey85", mid = "white", high = "#B22222", midpoint = 0, name = "Scaled avg") +
  theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = NULL, y = NULL, title = "KP genes in myeloid / T-cell subsets")
ggsave(file.path(dirs$figures, "FigureS4_KP_scRNA.pdf"),
       gridExtra::arrangeGrob(pA, pB, ncol = 2, widths = c(1.15, 1)), width = 15, height = 7.5)
write.csv(long, file.path(dirs$results, "17_kp_by_major.csv"), row.names = FALSE)
write.csv(heat, file.path(dirs$results, "17_kp_by_minor.csv"), row.names = FALSE)
message("17 done")
